!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2023 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

module dynamics_module
  use crest_parameters
  use crest_calculator
  use strucrd
  use atmasses
  use shake_module
  use metadynamics_module
!$ use omp_lib
  implicit none

  !======================================================================================!
  !--- private module variables and parameters
  private

  !>--- some constants and name mappings
  real(wp),parameter :: amutoau = amutokg*kgtome
  !real(wp),parameter :: fstoau = 41.3413733365614_wp

!========================================================================================!
!>--- LOCAL EXTENSION: process-wide state of the "stop"-file poll.
!>
!> dynamics() looks for a user-created file (default name: 'stop') in the working
!> directory while the simulation is running.  If that file shows up, the MD is
!> ended early and everything sampled so far is handed on to the next stage, so a
!> long run can be cut short instead of being killed and losing everything.
!>
!> stopfile_seen0 remembers whether the file was ALREADY there when the very
!> first dynamics() of this process started.  If it was, it is a leftover from an
!> earlier run and must never trigger - otherwise the next crest run in the same
!> directory would stop after a few minutes for no visible reason, and the user
!> would have to guess why.  Both variables are only touched inside an
!> !$omp critical, because dynamics() runs concurrently in the parallel MTD.
!========================================================================================!
  logical :: stopfile_initialized = .false.
  logical :: stopfile_seen0 = .false.
!>--- the directory crest was STARTED from.  Only differs from the working
!> directory when "-scratch" is in play: scratch.f90:54 chdir()s into the scratch
!> directory, and from that moment on a relative "stop" would be looked up in a
!> directory the user is not watching.  Set once from parseflags(), where
!> env%homedir is filled by getcwd() before any chdir() can have happened.
  character(len=512) :: stopfile_launchdir = ''
!>--- LOCAL CHANGE (2026-10-05): the name as configured (CREST_STOP_FILE or
!> 'stop'), and the round-level flag "somebody really did stop on it".
!>
!> md_stop_cleanup() uses both to DELETE the file once the whole parallel round
!> is over, so it cannot cut the NEXT MTD round short as well.  Deleting it at
!> the moment of detection would be wrong: the jobs poll on their own 300 s
!> wall-clock schedule and a step can take minutes, so the second job's poll
!> typically lands well after the first one's - it would simply not find the
!> file any more and run to its full length, leaving a half-stopped batch.
  character(len=256) :: stopfile_name = ''
  logical :: stopfile_hit = .false.

  !>-- filetypes as integers
  integer,parameter,public :: type_md = 1
  integer,parameter,public :: type_mtd = 2

  !>-- REEXPORTS from metadynamics_module
  public :: mtdpot,mtd_ini,cv_dump,calc_mtd
  public :: cv_std_mtd,cv_rmsd,cv_rmsd_static

  !>-- REEXPORTS from shake_module
  public :: shakedata,init_shake
  public :: crest_hard_bonds_to_shake

  public :: mddata
  !======================================================================================!
  !data object that contains settings for a molecular dynamics simulation.
  type :: mddata

    logical :: requested = .false.

    integer :: md_index = 0      !> some index for parallelization
    integer :: simtype = type_md !> type of the molecular dynamics simulation
    logical :: restart = .false.
    character(len=:),allocatable :: restartfile
    character(len=:),allocatable :: trajectoryfile
    !>--- data
    real(wp) :: length_ps = 0.0_wp !20.0_wp  !> total simulation length in ps
    integer :: length_steps = 0 !20000       !> total simulation length in steps
    real(wp) :: tstep = 0.0_wp !1.0_wp       !> timestep in fs
    real(wp) :: dumpstep = 0.0_wp !1000.0_wp !> snapshot dump step in fs
    integer :: sdump = 0 !1000 !> snapshot dump to trajectory every x steps
    integer :: dumped = 0      !> count how many structures have been written
    integer :: printstep = 200 !> control how often (in steps) to print to stdout

    !>--- LOCAL EXTENSION: coarse progress marks.
    !> The MD printout above is either full per-step thermodynamics (pr=.true.)
    !> or completely silent (pr=.false.).  For a *silent* run - the trial MTD is
    !> one - there is no way to tell whether the job is progressing at all.
    !> Setting %progress prints exactly three progress lines, at the steps
    !>   progress_steps, 2*progress_steps, 3*progress_steps.
    !> %progress_steps <= 0 is derived automatically as length_steps/4, i.e.
    !> marks at 25 / 50 / 75 % - deliberately NOT at 100 % (see below).
    logical :: progress = .false.      !> print coarse progress marks?
    integer :: progress_steps = 0      !> every N steps (<=0 -> length_steps/4)
    character(len=32) :: progress_tag = 'MD'  !> label in front of the marks

    !>--- LOCAL EXTENSION: poll for the user-created "stop" file.
    !> The interval itself is a process-wide setting (CREST_STOP_CHECK_SEC, see
    !> md_stopcheck_setup below); this switch only lets a CALLER opt out.
    !> The one caller that must: setuptest.f90 (the trial MTD) reads io /= 0 as
    !> "did not converge -> shorten the time step and retry", and after six tries
    !> it aborts the whole run.  A stop file must never be able to trigger that.
    logical :: stopcheck = .true.      !> poll for the "stop" file?

    real(wp) :: md_hmass = 0.0_wp !1.00794075_wp !> hydrogen mass

    logical :: shake = .true. !> use SHAKE algorithm
    integer :: nshake = 0  !> number of bonds considered in SHAKE
    type(shakedata) :: shk !> SHAKE bond information

    real(wp) :: tsoll = 0.0_wp !298.15_wp  !> wanted temperature
    logical :: thermostat = .true. !> apply thermostat?
    character(len=64) :: thermotype = 'berendsen'
    real(wp) :: thermo_damp = 500.0_wp !> thermostat damping parameter
    logical :: samerand = .false.

    integer :: blockl = 0    !> block length in MD steps
    integer :: iblock = 0    !> block counter
    integer :: nblock = 0    !> continuous block counter
    integer :: blocknreg = 0 !> block regression points
    integer :: maxblock
    real(wp),allocatable :: blockrege(:)
    real(wp),allocatable :: blocke(:)
    real(wp),allocatable :: blockt(:)

    !>--- the collection of MTD potentials
    integer :: npot = 0
    type(mtdpot),allocatable :: mtd(:)
    integer,allocatable :: cvtype(:)

    !>--- on-the-fly multi-scale level choosing
    integer,allocatable :: active_potentials(:)

  contains
    generic,public :: add => md_add_mtd
    procedure,private :: md_add_mtd
    procedure :: defaults => md_defaults_fallback
  end type mddata

  public :: dynamics
  public :: mdautoset
!>--- called once from parseflags() with env%homedir, i.e. BEFORE any chdir()
  public :: md_stopcheck_setdir
!>--- LOCAL CHANGE (2026-10-05): round-level stop-file housekeeping.  The caller
!> asks md_stop_seen() AFTER its parallel region and, if it is .true., calls
!> md_stop_cleanup() - which deletes the file and resets the flag.
  public :: md_stop_seen,md_stop_cleanup

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine dynamics(mol,dat,calc,pr,term)
!*************************************************************
!* subroutine dynamics
!* perform a molecular dynamics simulation
!* the coordinate propagation is made with an 
!* Leap-Frog algorithm (Velert-type algo)
!*************************************************************
    implicit none

    type(coord) :: mol          !> molecule data (should be in Bohr)
    type(mddata) :: dat         !> MD data
    type(calcdata),target :: calc      !> calculation control
    logical,intent(in) :: pr    !> printout control
    integer,intent(out) :: term !> termination status

    integer :: t,nfreedom
    real(wp) :: tstep_au
    real(wp) :: epot,ekin,edum
    real(wp) :: temp,thermoscal
!>--- averages & errors
    real(wp) :: Tav,Epav,Ekav,Eerror

    real(wp),allocatable :: grd(:,:)
!>--- LOCAL EXTENSION (2026-09-27): scratch for the MTD bias gradient.
!> It used to be allocated/deallocated inside md_calc_mtd() on EVERY MD step,
!> and both the allocate and the deallocate were wrapped in an (unnamed, i.e.
!> process-wide) !$omp critical - which serialised every MD job against all the
!> others twice per step for no reason at all: the array is private to the
!> calling thread either way.  Allocating it once per MD run removes both the
!> heap round trip and the two global locks from the step loop.
    real(wp),allocatable :: grdmtd(:,:)
    real(wp),allocatable :: velo(:,:)
    real(wp),allocatable :: vel(:,:)
    real(wp),allocatable :: veln(:,:)
    real(wp),allocatable :: acc(:,:)
    real(wp),allocatable :: mass(:)
    real(wp),allocatable :: xyz_angstrom(:,:)
    real(wp),allocatable :: backupweights(:)
    type(coord) :: molo
    real(wp) :: f,rt,rtshift
    real(wp) :: molmass,tmass
    character(len=:),allocatable :: trajectory
    integer :: trj
    character(len=256) :: commentline
    integer :: i,j,k,l,ich,och,io
    integer :: dcount,printcount
    logical :: ex,fail,bdump
    !>--- LOCAL EXTENSION: wall clock for the coarse progress marks
    integer(int64) :: pcount0,pcount1,prate
!>--- LOCAL EXTENSION: wall clock + state of the "stop"-file poll.
!> sinterval <= 0 disables the poll; stoppre says the file was already there when
!> this process started (see stopfile_seen0 above), stopgo is the trigger.
    integer(int64) :: scount0
    real(wp) :: sinterval
    character(len=:),allocatable :: stopfile
    logical :: stoppre,stopgo

    call initsignal()

!>--- pre-settings and calculations
    !$omp critical
    call dat%defaults() !> check for unset parameters
!>--- derive the progress interval from the total length (marks at 25/50/75 %)
    if (dat%progress .and. dat%progress_steps <= 0) then
      dat%progress_steps = max(1,dat%length_steps/4)
    end if
    term = 0
    tstep_au = dat%tstep*fstoau
    nfreedom = 3*mol%nat
    if(calc%nfreeze > 0)then
      nfreedom = nfreedom - 3*calc%nfreeze
    endif
    if (dat%shake) then
      if(calc%nfreeze > 0)then
        dat%shk%freezeptr => calc%freezelist
      else
        nullify(dat%shk%freezeptr)
      endif
      call init_shake(mol%nat,mol%at,mol%xyz,dat%shk,pr)
      dat%nshake = dat%shk%ncons
      nfreedom = nfreedom-dat%nshake
    end if
!>--- averages
    tav = 0.0_wp
    eerror = 0.0_wp
    ekav = 0.0_wp
    epav = 0.0_wp
    temp = 0.0_wp

!>--- on-the-fly multiscale definition
    if(allocated(dat%active_potentials))then
      call calc%active(dat%active_potentials)
    endif

!>--- allocate data fields
    allocate (xyz_angstrom(3,mol%nat))
    allocate (molo%at(mol%nat),molo%xyz(3,mol%nat))
    allocate (grd(3,mol%nat),vel(3,mol%nat),velo(3,mol%nat),source=0.0_wp)
    allocate (grdmtd(3,mol%nat),source=0.0_wp)
    allocate (veln(3,mol%nat),acc(3,mol%nat),mass(mol%nat),source=0.0_wp)
    allocate (dat%blocke(dat%blockl),dat%blockt(dat%blockl))
    allocate (dat%blockrege(dat%maxblock))
    !$omp end critical

!>--- settings printout
    if (pr) then
      write (*,*)
      write (*,'("> ",a)') 'Molecular dynamics settings'
      write (*,'('' MD time /ps        :'',f10.2)') dat%length_ps
      write (*,'('' dt /fs             :'',f10.2)') dat%tstep
      write (*,'('' temperature /K     :'',f10.2)') dat%tsoll
      write (*,'('' max steps          :'',i10  )') dat%length_steps
      write (*,'('' block length (av.) :'',i10  )') dat%blockl
      write (*,'('' dumpstep(trj) /fs  :'',f10.2,i6)') dat%dumpstep,dat%sdump
      write (*,'('' # deg. of freedom  :'',i10  )') nfreedom
      if(calc%nfreeze > 0)then
      write (*,'('' # frozen atoms     :'',i10  )') calc%nfreeze
      endif
      call thermostatprint(dat,pr)
      write (*,'('' SHAKE constraint   :'',8x,l)') dat%shake
      if (dat%shake) then
        if (dat%shk%shake_mode == 2) then
          write (*,'('' # SHAKE bonds      :'',i10,a)') dat%nshake,' (all bonds)'
        elseif (dat%shk%shake_mode == 1) then
          write (*,'('' # SHAKE bonds      :'',i10,a)') dat%nshake,' (H only)'
        end if
      end if
      write (*,'('' hydrogen mass /u   :'',f10.5 )') dat%md_hmass
      if(allocated(dat%active_potentials))then
       write (*,'('' active potentials  :'',i10)') size(dat%active_potentials,1)
      endif
    end if

!>--- set atom masses
    !$omp critical
    molmass = 0.0_wp
    do i = 1,mol%nat
      molmass = molmass+ams(mol%at(i))
      mass(i) = ams(mol%at(i))*amutoau !>-- ams from module atmasses
    end do
    tmass = molmass*amutoau
    do i = 1,mol%nat
      if (mol%at(i) .eq. 1.and.dat%md_hmass .gt. 0.0_wp) then
        mass(i) = dat%md_hmass*amutoau
      end if
    end do
    molmass = molmass*amutokg
    !$omp end critical

!>--- initialize velocities (or read from restart file)
    if (dat%thermostat) then
      f = 1.0_wp
    else
      f = 2.0_wp
    end if
    edum = f*dat%tsoll*0.5_wp*kB*float(nfreedom)
    rtshift = 0.0_wp
    if(.not.dat%restart .or. .not.allocated(dat%restartfile))then 
     call mdinitu(mol,dat,velo,mass,edum,pr)
    else
     call rdmdrestart(mol,dat,velo,fail,rtshift)
     if(fail)then
        call mdinitu(mol,dat,velo,mass,edum,pr)
     else
       call ekinet(mol%nat,velo,mass,ekin)
       temp = 2.0_wp*ekin/float(nfreedom)/kB
       tav = temp
     endif
    endif 
    call ekinet(mol%nat,velo,mass,ekin)
    if(calc%nfreeze > 0)then
       do i = 1,mol%nat
         if(calc%freezelist(i))then
           acc(:,i) = 0.0_wp
           grd(:,i) = 0.0_wp
           velo(:,i) = 0.0_wp
         endif
       end do
     endif

!>--- initialize MTDs (if required)
    !$omp critical
    if (dat%simtype == type_mtd) then
      call md_init_mtd(mol,dat,pr)
    end if
    !$omp end critical

    !>--- initialize trajectory file
    if (allocated(dat%trajectoryfile)) then
      trajectory = dat%trajectoryfile
    else
      write (commentline,'(a,i0,a)') 'crest_',dat%md_index,'.xyz'
      trajectory = trim(commentline)
    end if
    !$omp critical
    open (newunit=trj,file=trajectory)
    !$omp end critical

!>--- begin printout
    if (pr) then
      write (*,'(/,"> ",a)') 'Starting simulation'
      if (.not.dat%thermostat) then
        write (*,'(/,11x,"time (ps)",7x,"<Epot>",8x,"Ekin",5x,"<T>",7x,"T",12x, &
           &         "Etot",7x,"error")')
      else
        write (*,'(/,11x,"time (ps)",7x,"<Epot>",8x,"Ekin",5x,"<T>",7x,"T",12x, &
           &         "Etot")')
      end if
    end if

    dcount = 0
    printcount = 1
!>--- wall clock reference for the (optional) progress marks
    call system_clock(count=pcount0,count_rate=prate)
!>--- LOCAL EXTENSION: arm the "stop"-file poll.
!> md_stopcheck_setup() returns sinterval <= 0 when the poll is switched off
!> (either CREST_STOP_CHECK_SEC <= 0, or this caller opted out via
!> dat%stopcheck) - then the two lines inside the MD loop below are dead code
!> and nothing at all is stat()ed.
    sinterval = 0.0_wp
    stoppre = .false.
    stopgo = .false.
    if (dat%stopcheck) then
      call md_stopcheck_setup(sinterval,stopfile,stoppre)
      if (sinterval > 0.0_wp) call system_clock(count=scount0)
    end if
!===============================================================!
!===============================================================!
!>--- begin MD loop
    MD: do t = 1,dat%length_steps
      call initsignal()

!>--- LOCAL EXTENSION: "stop"-file poll.
!> Deliberately at the TOP of the step and not at the end: the whole point is to
!> skip the next engrad(), which for an external calculator (ORCA) is by far the
!> most expensive thing in this loop.  The poll itself is one system_clock() per
!> step (a few nanoseconds) plus one stat() per interval, so it is invisible
!> against the gradient call it replaces.
      if (sinterval > 0.0_wp) then
        call md_stopcheck(stopfile,sinterval,scount0,prate,stoppre,stopgo)
        if (stopgo) then
!> term = 3 -> "ended on request, take as converged".  It has to be a value of
!> its own: the callers cannot tell 2 ("something broke halfway through") apart
!> from a deliberate stop, and crest_moleculardynamics() would report the latter
!> as "terminated ABNORMALLY" + status_failed.
          term = 3
          call md_stopreport(stopfile,t,dat%length_steps,dat%dumped)
          exit MD
        end if
      end if

      !>>-- STEP 1: calculate energy and forces
      !>--- singlepoint calculation
      epot = 0.0_wp
      grd = 0.0_wp
      call engrad(mol,calc,epot,grd,io)

      if (io /= 0) then
        if (dat%dumped > 0) then
          term = 2  !> termination during MD
        else
          term = 1  !> termination upon first engrad call
        end if
!>--- LOCAL EXTENSION: say which of the two things in this loop actually failed.
!> Upstream reported every failure through this one integer and nothing else, so
!> a broken calculator and a broken constraint set looked identical.
        call md_failreport('SINGLEPOINT (energy/gradient)',t,dat%length_steps,term)
        exit MD
      end if
      if (t == 1) then
        edum = epot+ekin
      end if

      !>>-- STEP 1.5: calculate metadynamics bias
      if (dat%simtype == type_mtd) then
        !> MTD energy and gradient are added to epot and grd, respectively.
        !> grdmtd is the scratch allocated once above; see its declaration.
        call md_calc_mtd(mol,dat,epot,grd,grdmtd,pr)
        !$omp critical
        call md_update_mtd(mol,dat,calc,pr)
        !$omp end critical
      end if

      !>--- block data printouts
      call u_block(mol,dat,epot,temp,pr,bdump)
      !>--- MD restart files written for each block rather than at each timestep to reduce I/O
      if(bdump)then
        rt = float(t)*dat%tstep + rtshift
        !$omp critical
        call wrmdrestart(mol,dat,velo,rt) 
        !$omp end critical
      endif 

      !===========================================!
      !>>-- write to trajectory and printout
      if (dcount == dat%sdump) then
        dcount = 0
        dat%dumped = dat%dumped+1
        !$omp critical
        xyz_angstrom = mol%xyz*bohr
        write (commentline,'(a,f22.12,1x,a)') 'Epot =',epot,''
        call wrxyz(trj,mol%nat,mol%at,xyz_angstrom,commentline)
!>--- LOCAL EXTENSION: push the snapshot out of the runtime I/O buffer.
!> The trajectory file is opened once (above) and only closed when the MD is
!> completely over. Everything written in between just sits in the runtime
!> buffer and only reaches the disk when that buffer happens to run full.
!> That has two visible consequences:
!>   * a trajectory that is inspected WHILE the MD is still running ends with
!>     a HALF-WRITTEN frame - the buffer boundary cuts a frame in the middle,
!>     which looks exactly like "the file was written but never completed";
!>   * every frame that is still buffered is LOST when the process dies
!>     (Ctrl+C, walltime, killed job, or creststop() -> exit()).
!> One flush per frame costs a single write(2) per snapshot, which is
!> irrelevant next to the energy/gradient call that produced that snapshot.
!> stdout is flushed on the same cadence: a redirected crest.out is otherwise
!> frozen for hours during a silent (pr=.false.) parallel MTD run, because
!> %progress defaults to .false. and therefore nothing flushes unit 6 at all.
        flush (trj)
        flush (stdout)
        !$omp end critical
      end if
      if ((printcount == dat%printstep).or.(t == 1)) then
        if (t > 1) printcount = 0
        if (pr) then
          rt = float(t)*dat%tstep + rtshift
          if (.not.dat%thermostat) then
            write (*,'(i7,f10.2,F16.5,F12.4,2F8.1,F16.5,4F10.4)') &
               &   t,0.001_wp*rt, (Epav+Epot)/float(t), &
               &   Ekin,Tav/float(t),temp,Epot+Ekin, &
               &   Edum/float(t)-Epot-Ekin
          else
            write (*,'(i7,f10.2,F16.5,F12.4,2F8.1,F16.5)') &
               &   t,0.001_wp*rt, (Epav+epot)/float(t), &
               &   Ekin,Tav/float(t),temp,Epot+Ekin
          end if
!>--- LOCAL EXTENSION: same reasoning as for the trajectory file above.
!> Without this the MD table only shows up in blocks of one buffer size.
          flush (stdout)
        end if
      end if

!>--- LOCAL EXTENSION: coarse progress marks, independent of "pr".
!> This is the only feedback a silent (pr=.false.) run gives, e.g. the trial
!> MTD.
!> EXACTLY THREE marks, at 25/50/75 %.  Do NOT use "mod(t,progress_steps)==0":
!> progress_steps is floor(length_steps/4), so whenever length_steps is not a
!> multiple of four a FOURTH mark is emitted at t = 4*floor(length_steps/4),
!> which is still < length_steps, and its rounded percentage comes out as 99
!> or 100, e.g.
!>      length_steps = 4k+1  ->  t = 4k  ->  pct = 400k/(4k+1)  ->  "100%"
!> That line is both useless (the MD is essentially over) and confusing: it
!> reports a remaining time of a couple of minutes and is immediately followed
!> by the job's own "... completed successfully ... <total runtime>" line.
!> Comparing t against the three target steps directly also makes the 100 %
!> mark impossible by construction, instead of relying on the displayed
!> percentage (which was the earlier - insufficient - attempt).
      if (dat%progress.and.dat%progress_steps > 0) then
        if (t < dat%length_steps.and. &
        &  (t == dat%progress_steps .or. t == 2*dat%progress_steps .or. &
        &   t == 3*dat%progress_steps)) then
          !$omp critical
          call md_progress(stdout,dat,t,pcount0,prate)
          !$omp end critical
        end if
      end if

      !>--- compute the acceleration at t
      do i = 1,mol%nat
        acc(:,i) = -grd(:,i)/mass(i)
      end do

      !>--- special setup for frozen atoms
      if(calc%nfreeze > 0)then
       do i = 1,mol%nat
         if(calc%freezelist(i))then
           acc(:,i) = 0.0_wp
           grd(:,i) = 0.0_wp
         endif
       end do
      endif

      !>--- store positions (at t); velocities are at t-1/2dt
!>--- LOCAL CHANGE (2026-09-27): no !$omp critical here.  `molo` is a local of
!> dynamics(), i.e. private to this thread's invocation, and the only thing that
!> was inside the critical was a copy of thread-private data - it merely
!> serialised every MD job against all the others once per step.  The
!> molo%nat / molo%at copies are gone as well: nothing in this module ever READS
!> them (only molo%xyz is used, at the position updates and in do_shake).
      molo%xyz = mol%xyz

      !>>-- STEP 2: temperature and pressure/density control
      !>--- estimate(!) velocities at t
      veln = velo+0.5_wp*tstep_au*acc

      !>--- compute kinetic energy and termperature
      call ekinet(mol%nat,veln,mass,ekin)
      temp = 2.0_wp*ekin/float(nfreedom)/kB

      !>--- THERMOSTATING (determine factor thermoscal)
      call thermostating(mol,dat,temp,thermoscal)

      !>>-- STEP 3: velocity and position update
      !>--- update velocities to t
      vel = thermoscal*(velo + acc*tstep_au)

      !>--- update positions to t+dt, except for frozen atoms
      if(calc%nfreeze > 0)then
       do i = 1,mol%nat
         if(.not.calc%freezelist(i))then
           mol%xyz(:,i) = molo%xyz(:,i)+vel(:,i)*tstep_au
         else
           vel(:,i) = 0.0_wp
         endif
       end do
      else  
        mol%xyz = molo%xyz+vel*tstep_au
      endif

      !>--- estimate new velocities at t
      veln = 0.5_wp*(velo+vel)

      !>--- compute kinetic energy and temperature for average tracking
      call ekinet(mol%nat,veln,mass,ekin)
      temp = 2.0_wp*ekin/float(nfreedom)/kB

      !>--- apply SHAKE at t+dt?
      if (dat%shake.and.dat%nshake > 0) then
        call do_shake(mol%nat,molo%xyz,mol%xyz,vel,acc,mass,tstep_au,dat%shk,pr,io)
        if (io /= 0) then
          if (dat%dumped > 0) then
            term = 2  !> termination during MD
          else
            term = 1  !> termination upon first engrad call
          end if
!>--- LOCAL EXTENSION: same reason as above - and here it matters even more,
!> because a SHAKE failure is entirely insensitive to the time step while the
!> standard reaction to a failed trial MTD is exactly that: shorten it.
!> The detailed per-bond report comes from shake_module:shake_failreport().
          call md_failreport('SHAKE (rigid constraints)',t,dat%length_steps,term)
          exit MD
        end if
      end if

      !>--- update velocities
      velo = vel

      !>--- remove translational and rotational componetnts of the velocity
      call rmrottr(mol%nat,mass,velo,mol%xyz)

      !>>-- Update averages and counter
      edum = edum+epot+ekin
      eerror = edum/float(t)-epot-ekin
      tav = tav+temp
      epav = epav+epot
      ekav = ekav+ekin
      dcount = dcount+1
      printcount = printcount+1

    end do MD
!>--- finish MD loop
!===============================================================!
!===============================================================!
!>--- close trajectory file
    !$omp critical
    close (trj)
    !$omp end critical

!>--- averages printout
    if (pr) then
      write (*,*)
      write (*,*) 'average properties '
      write (*,*) '----------------------'
      write (*,*) '<Epot> / Eh          :',Epav/float(t)
      write (*,*) '<Ekin> / Eh          :',Ekav/float(t)
      write (*,*) '<Etot> / Eh          :', (Ekav+Epav)/float(t)
      write (*,*) '<T> / K              :',Tav/float(t)
    end if

!>--- write restart file
    rt = float(dat%length_steps)*dat%tstep + rtshift
    call wrmdrestart(mol,dat,velo,rt)

!>--- termination printout
    if (pr) then
      select case (term)
      case (0)
        write (*,*) 'normal MD termination'
      case (1)
        write (stderr,*) 'error in MD calculation'
      case (2)
        write (*,*) 'MD terminated, but still taking as converged.'
      case (3)
        write (*,*) 'MD stopped on request, but still taking as converged.'
      end select
    end if

!>--- deallocate data
    deallocate (dat%blockrege,dat%blockt,dat%blocke)
    deallocate (mass,acc,veln)
    deallocate (vel,velo,grd)
    deallocate (grdmtd)
    deallocate (molo%xyz,molo%at)
    deallocate (xyz_angstrom)

!>--- restore weights if necessary
    if(allocated(dat%active_potentials))then
      call calc%active_restore()
    endif

    return
!========================================================================================!
!========================================================================================!
  contains
!========================================================================================!
!>--- LOCAL EXTENSION: tell apart the two failure modes of the MD loop.
!>
!> There are exactly two ways this loop can end early, and both of them used to
!> surface as the same integer `term` with no further comment:
!>   1. engrad() failed            -> the calculator (xtb/ORCA subprocess) is at fault
!>   2. do_shake() did not converge -> the CONSTRAINT SET is at fault
!> (1) can depend on the geometry and hence indirectly on the time step, (2) does
!> not.  Printing which one it is, plus the step at which it happened, is what
!> turns "Trial MTD n did not converge!" into something actionable.
!>
!> Printing unconditionally is safe: both paths exit the loop immediately, so a
!> given MD emits at most one report.
!========================================================================================!
    subroutine md_failreport(origin,t,nsteps,termstatus)
      implicit none
      character(len=*),intent(in) :: origin
      integer,intent(in) :: t,nsteps,termstatus
      !$omp critical
      write (stdout,'(/,1x,a)') '> MD terminated early'
      write (stdout,'(3x,a,a)') 'failed in                  : ',trim(origin)
      write (stdout,'(3x,a,i0,a,i0)') 'MD step                    : ',t,' of ',nsteps
      write (stdout,'(3x,a,i0,a)') 'termination status         : ',termstatus, &
      &   '  (1 = failed in step 1, 2 = failed later on)'
!>--- LOCAL FIX (2026-10-24): the old text claimed a first-step failure means the
!> starting geometry already violates the constraints.  That is WRONG for the
!> SHAKE origin: the r0 targets are taken from the starting/'coord' structure,
!> so it satisfies them by construction (measured 0.00000 Angstrom on a real
!> system, 0.0018 Angstrom after the pre-optimisation).  What is large at step 1
!> is how far ONE raw step moved the constrained distances - a timestep /
!> force-scale matter, not a geometry mismatch.
!> For the SINGLEPOINT origin there are no constraints at all, so the SHAKE
!> wording would be doubly misleading; that branch keeps a neutral message.
      if (termstatus == 1) then
        if (index(origin,'SHAKE') > 0) then
          write (stdout,'(3x,a)') 'This is NOT a geometry mismatch: the r0 targets come from the'
          write (stdout,'(3x,a)') 'starting structure, so it satisfies them by construction.'
          write (stdout,'(3x,a)') 'A step-1 failure means one single step moved the constrained'
          write (stdout,'(3x,a)') 'distances too far (see "deviation AFTER one raw step" above).'
          write (stdout,'(3x,a)') 'Check, in this order: (1) is the constraint group a nearly'
          write (stdout,'(3x,a)') 'coplanar/over-determined cluster? (2) is tstep too large?'
          write (stdout,'(3x,a)') '(3) is maxcyc (now 6000) still too small for rho ~ 0.997?'
        else
          write (stdout,'(3x,a)') 'The FIRST MD step already failed.  This is not a constraint'
          write (stdout,'(3x,a)') 'problem - see the energy/gradient report above for the cause.'
        end if
      end if
      write (stdout,*)
      flush (stdout)
      !$omp end critical
      return
    end subroutine md_failreport
!========================================================================================!
!>--- LOCAL EXTENSION: the "stop"-file poll itself.
!>
!> Called once per MD step, but does real work only every `interval` seconds of
!> WALL CLOCK - a step-count gate would be useless here, because how much wall
!> time a step costs differs by orders of magnitude between GFN2 (milliseconds)
!> and an external ORCA call (seconds to minutes).
!>
!> The trigger is `ex .and. .not. pre`: the file has to APPEAR while crest is
!> running.  `pre` carries stopfile_seen0, i.e. "it was there before we started".
!========================================================================================!
    subroutine md_stopcheck(fname,interval,last,rate,pre,hit)
      implicit none
      character(len=*),intent(in) :: fname
      real(wp),intent(in) :: interval
      integer(int64),intent(inout) :: last
      integer(int64),intent(in) :: rate
      logical,intent(in) :: pre
      logical,intent(out) :: hit
      integer(int64) :: now
      real(wp) :: el
      logical :: ex
      hit = .false.
      call system_clock(count=now)
      el = real(now-last,wp)/real(max(rate,1_int64),wp)
      if (el < interval) return
      !$omp critical
      last = now
      ex = md_stopfile_exists(fname)
      !$omp end critical
      hit = (ex.and.(.not.pre))
      return
    end subroutine md_stopcheck
!========================================================================================!
!>--- LOCAL EXTENSION: what the user needs to see when the poll fires.
!> Printed unconditionally (not behind `pr`): a parallel MTD runs with
!> pr=.false., and "why did my run stop after 5 minutes?" must be answerable
!> from the log alone.  At most one report per MD, same as md_failreport.
!========================================================================================!
    subroutine md_stopreport(fname,t,nsteps,dumped)
      implicit none
!>--- NOTE: `fname` is no longer printed (the headline names the event, and the
!> file it belongs to is the one the caller deletes) but the dummy is kept so
!> the single call site does not have to change.
      character(len=*),intent(in) :: fname
      integer,intent(in) :: t,nsteps,dumped
      !$omp critical
!>--- LOCAL CHANGE (2026-10-05): short report.  The "Remove the file before the
!> next run." paragraph is gone because crest now removes it itself
!> (md_stop_cleanup(), called by the caller once the whole round is over), and
!> the trailing blank line is gone with it.
      stopfile_hit = .true.
      write (stdout,'(/,1x,a)') '> stop file found, MD stopped on user request.'
      write (stdout,'(3x,a,i0,a,i0)') 'MD step                    : ',t,' of ',nsteps
      write (stdout,'(3x,a,i0)') 'snapshots written so far   : ',dumped
      flush (stdout)
      !$omp end critical
      return
    end subroutine md_stopreport
!========================================================================================!
!>--- LOCAL EXTENSION: one coarse progress line, e.g.
!> ("*MTD   2  25% (  1.50 /  6.00 ps, elapsed 3h32m21s, remaining 10h37m03s)").
!> Kept inside dynamics() so that no public interface changes.
!========================================================================================!
    subroutine md_progress(iounit,dat,t,t0,rate)
      implicit none
      integer,intent(in) :: iounit
      type(mddata),intent(in) :: dat
      integer,intent(in) :: t
      integer(int64),intent(in) :: t0,rate
      integer(int64) :: t1
      real(wp) :: pct,ela,eta,sim
      character(len=32) :: sela,seta
      call system_clock(count=t1)
      pct = 100.0_wp*real(t)/real(max(1,dat%length_steps))
      !> guard against a nonsensical rate (0 on exotic platforms)
      ela = real(t1-t0,wp)/real(max(rate,1_int64),wp)
      sim = real(t)*dat%tstep*0.001_wp
      if (pct > 1.0e-6_wp) then
        eta = ela*(100.0_wp-pct)/pct
      else
        eta = ela
      end if
      call md_sec2str(ela,sela)
      call md_sec2str(eta,seta)
!> LOCAL EXTENSION (2026-09-24): the line used to be
!>     '*MTD   2  25%  (   1.50 /   6.00 ps,  elapsed 3 h 32 min 21 sec, expected remaining time 10 h 37 min 3 sec)'
!> i.e. 108 columns - it wrapped on a narrow terminal and, because four jobs
!> interleave their marks inside !$omp critical, a wrapped line is unreadable.
!> Squeezed to 73 columns by (a) removing the double spaces, (b) f7.2 -> f6.2
!> (still good up to 999.99 ps), and (c) the compact h/m/s form below.  The
!> label is 'remaining', not the old 'expected remaining time' - inside the
!> 'elapsed X, remaining Y' pair it is unambiguous, and it is NOT the cryptic
!> 'ETA' that was rejected earlier.  Everything that is still padded
!> (the 'a' tag, the i3 percentage) stays padded ON PURPOSE: those are what
!> align the columns between jobs 1..N.
      write (iounit,'(1x,a,1x,i3,"% (",f6.2," /",f6.2," ps, elapsed ",a,", remaining ",a,")")') &
      &  trim(dat%progress_tag),nint(pct),sim,dat%length_ps,trim(sela),trim(seta)
      !> without this the line sits in the buffer and the progress is invisible
      !> in a redirected (SLURM) log until the run ends - which defeats the point
      flush (iounit)
      return
    end subroutine md_progress

    subroutine md_sec2str(sec,str)
      implicit none
      real(wp),intent(in) :: sec
      character(len=*),intent(out) :: str
      real(wp) :: hours,minutes,seconds
      hours = aint(sec/3600.0_wp)
      minutes = aint((sec-3600.0_wp*hours)/60.0_wp)
      seconds = anint(sec-(3600.0_wp*hours+60.0_wp*minutes))
!> anint() rounds, so 59.6 s becomes 60 s - and after the carry the minutes can
!> reach 60 as well.  Normalise BEFORE picking the format, otherwise the
!> compact string reads '3h32m60s'.  (The old spaced form had the same defect,
!> it was simply less obvious: '3 h 32 min 60 sec'.)
      if (nint(seconds) >= 60) then
        seconds = 0.0_wp
        minutes = minutes+1.0_wp
      end if
      if (nint(minutes) >= 60) then
        minutes = 0.0_wp
        hours = hours+1.0_wp
      end if
      if (hours > 0.0_wp) then
        write (str,'(i0,"h",i2.2,"m",i2.2,"s")') nint(hours),nint(minutes),nint(seconds)
      else if (minutes > 0.0_wp) then
        write (str,'(i2.2,"m",i2.2,"s")') nint(minutes),nint(seconds)
      else
        write (str,'(i0,"s")') nint(seconds)
      end if
      return
    end subroutine md_sec2str
!========================================================================================!
  end subroutine dynamics
!========================================================================================!
!>--- LOCAL FIX: round-level stop-file housekeeping.
!>
!> WHY THESE LIVE AT MODULE LEVEL AND NOT INSIDE dynamics(): md_stop_seen() and
!> md_stop_cleanup() are declared `public` in the specification part above and
!> are called from parallel.f90 AFTER its parallel region.  A `public ::
!> name` statement can only refer to MODULE entities, so the original layout -
!> the three routines nested inside subroutine dynamics() - did not compile:
!>     Symbol 'md_stop_seen' at (1) has no IMPLICIT type
!>     Symbol 'md_stop_cleanup' at (1) has no IMPLICIT type
!> They touch only module variables (stopfile_hit / stopfile_name /
!> stopfile_launchdir) and each other, so relocating them changes nothing else.
!========================================================================================!
!>--- did any MD of the round that just ended stop on the stop file?
!> (Ask AFTER the parallel region, never inside it.)
!========================================================================================!
  logical function md_stop_seen() result(hit)
    implicit none
    hit = stopfile_hit
    return
  end function md_stop_seen
!========================================================================================!
!>--- delete the honoured stop file and re-arm the flag, so the NEXT dynamic
!> round is not cut short by the same file.
!>
!> WHY THE CALLER HAS TO DO IT AND NOT THE DETECTING JOB: the jobs poll on
!> their own 300 s wall-clock schedule and a single MD step can take minutes,
!> so the second job's poll regularly lands well after the first one's.  A file
!> deleted at the moment of detection would simply not be there any more for
!> it, and that job would run to its full length - a half-stopped batch.  By
!> deleting after the parallel region, the file stays visible to every job for
!> the whole round and is gone before the next one starts.
!>
!> Both probe locations of md_stopfile_exists() are cleaned: the working
!> directory (the one the user watches) and, when it differs, the launch
!> directory.  An absolute name is only ever touched as given.
!========================================================================================!
  subroutine md_stop_cleanup()
    implicit none
    character(len=1024) :: alt
    logical :: ex
    if (.not.stopfile_hit) return
    if (len_trim(stopfile_name) < 1) return
    inquire (file=trim(stopfile_name),exist=ex)
    if (ex) call md_stopfile_delete(stopfile_name)
    if (len_trim(stopfile_launchdir) > 0) then
      if (stopfile_name(1:1) /= '/') then
        alt = trim(stopfile_launchdir)//'/'//trim(stopfile_name)
        inquire (file=trim(alt),exist=ex)
        if (ex) call md_stopfile_delete(trim(alt))
      end if
    end if
    stopfile_hit = .false.
    return
  end subroutine md_stop_cleanup
!========================================================================================!
!>--- delete one file.  open(status='old') + close(status='delete') is plain
!> Fortran 2003 and saves pulling iomod's remove() into this module.
!========================================================================================!
  subroutine md_stopfile_delete(fname)
    implicit none
    character(len=*),intent(in) :: fname
    integer :: u,ios
    open (newunit=u,file=trim(fname),status='old',iostat=ios)
    if (ios /= 0) return
    close (u,status='delete',iostat=ios)
    return
  end subroutine md_stopfile_delete
!========================================================================================!
!========================================================================================!
!>--- LOCAL EXTENSION: arm the "stop"-file poll of one dynamics() call.
!>
!> Two environment variables, mirroring the CREST_* knobs already used elsewhere
!> in this tree (and read BY VALUE, not by mere presence):
!>   CREST_STOP_CHECK_SEC  poll interval in SECONDS.  Default 300 (= 5 min).
!>                         0 or negative switches the whole feature off.
!>   CREST_STOP_FILE       name of the file to look for.  Default 'stop'.
!>
!> The interval is in seconds and not in minutes on purpose: "5" would be
!> ambiguous, and the unit has to be guessable from the name alone.
!>
!> The reference inquire() at the end ("was the file there when this process
!> started?") is what keeps a LEFTOVER stop file harmless.  Without it, a user
!> who stopped a run last week and starts a new one in the same directory today
!> would see that run die after five minutes with no explanation.
!========================================================================================!
!========================================================================================!
!>--- LOCAL EXTENSION: remember the directory crest was started from.
!========================================================================================!
  subroutine md_stopcheck_setdir(dir)
    implicit none
    character(len=*),intent(in) :: dir
    stopfile_launchdir = trim(dir)
    return
  end subroutine md_stopcheck_setdir
!========================================================================================!
!>--- LOCAL EXTENSION: one "is the stop file there?" probe, checking BOTH places.
!>
!> 1. the current working directory - which is where crest runs, and which under
!>    "-scratch" is the scratch directory, not the one the user is looking at;
!> 2. the launch directory (stopfile_launchdir), when it is known and different.
!>
!> Probing both is what makes "touch stop" next to the toml/xyz work in every
!> mode instead of only in the default one.  An ABSOLUTE name is only ever
!> probed as given - prepending the launch directory to "/home/me/stop" would
!> turn it into nonsense.
!========================================================================================!
  logical function md_stopfile_exists(fname) result(ex)
    implicit none
    character(len=*),intent(in) :: fname
    character(len=1024) :: alt
    logical :: e1,e2
    ex = .false.
    if (len_trim(fname) < 1) return
    inquire (file=fname,exist=e1)
    if (e1) then
      ex = .true.
      return
    end if
    if (len_trim(stopfile_launchdir) > 0) then
      if (fname(1:1) /= '/') then
        alt = trim(stopfile_launchdir)//'/'//trim(fname)
        inquire (file=trim(alt),exist=e2)
        ex = e2
      end if
    end if
    return
  end function md_stopfile_exists
!========================================================================================!
  subroutine md_stopcheck_setup(interval,fname,pre)
    implicit none
    real(wp),intent(out) :: interval
    character(len=:),allocatable,intent(out) :: fname
    logical,intent(out) :: pre
    character(len=256) :: envstr
    integer :: stat,rd
    real(wp) :: val
    logical :: ex

    interval = 300.0_wp
    fname = 'stop'
    pre = .false.

    call get_environment_variable('CREST_STOP_FILE',value=envstr,status=stat)
    if (stat == 0) then
      if (len_trim(envstr) > 0) fname = trim(envstr)
    end if
!>--- LOCAL CHANGE (2026-10-05): remember the name so md_stop_cleanup() can
!> delete exactly the file that was honoured (same value in every job).
    stopfile_name = trim(fname)

    envstr = ''
    call get_environment_variable('CREST_STOP_CHECK_SEC',value=envstr,status=stat)
    if (stat == 0) then
      if (len_trim(envstr) > 0) then
        read (envstr,*,iostat=rd) val
!> an unreadable value keeps the default rather than silently disabling a
!> feature the user clearly asked for
        if (rd == 0) interval = val
      end if
    end if

    if (interval <= 0.0_wp) then
      interval = 0.0_wp
      return
    end if
!> below one second the poll degenerates into one stat() per MD step
    if (interval < 1.0_wp) interval = 1.0_wp

    !$omp critical
    ex = md_stopfile_exists(fname)
    if (.not.stopfile_initialized) then
      stopfile_initialized = .true.
      stopfile_seen0 = ex
    end if
    pre = stopfile_seen0
    !$omp end critical
    return
  end subroutine md_stopcheck_setup
!========================================================================================!
!========================================================================================!

!========================================================================================!
  subroutine mdautoset(dat,iostatus)
!************************************************
!* subroutine mdautoset
!* convert real-time settings (ps,fs) to steps
!* Intended to restore settings once reset
!************************************************
    implicit none
    type(mddata) :: dat
    integer,intent(out) :: iostatus
    real(wp) :: dum
    integer :: idum
    integer :: i,j,k,l,ich,och,io
    logical :: ex

    iostatus = 0

    if (dat%length_ps .le. 0.0_wp.or. &
    &  dat%tstep .le. 0.0_wp) then
      write (stderr,*) 'need valid simulation length and time step!'
      write (stderr,*) 'abort MD.'
      iostatus = -1
      return
    end if

    !>--- MD length to steps
    dum = (dat%length_ps*1000.0_wp)/dat%tstep
    dat%length_steps = nint(dum)

    !>--- adjust structure dump to trajectory file
    dum = max(1.0_wp, (dat%dumpstep/dat%tstep))
    dat%sdump = nint(dum)

    return
  end subroutine mdautoset

!========================================================================================!
! subroutine ekinet
! calculate kinetic energy from velocities and masses
  subroutine ekinet(n,velo,mass,E)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: velo(3,n),mass(n)
    real(wp),intent(out) :: e
    integer :: i
    e = 0.0_wp
    do i = 1,n
      e = e+mass(i)*(velo(1,i)**2+velo(2,i)**2+velo(3,i)**2)
    end do
    e = e*0.5_wp
    return
  end subroutine ekinet

!========================================================================================!
! subroutine u_block
! update block data and printout
  subroutine u_block(mol,dat,epot,temp,pr,bdump)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(in) :: epot
    real(wp),intent(in) :: temp
    logical,intent(in) :: pr
    logical,intent(out) :: bdump
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    integer :: nreg
    real(wp) :: bave,bavt,slope

    bdump = .false.
    if (dat%iblock == dat%blockl) then
      dat%nblock = dat%nblock+1
      dat%iblock = 0
      bdump = .true.
      call blocksd(mol%nat,dat%blockl,dat%blocke,dat%blockt,bave,bavt)
      dat%blocknreg = dat%blocknreg+1
      nreg = dat%blocknreg
      dat%blockrege(nreg) = bave
      if (nreg .ge. 4) then
        call regress(nreg-3,nreg,dat%blockrege,slope)
      else
        slope = 99.0_wp
      end if
      if (pr) then
        write (*,'(''block <Epot> / <T> :'',f14.5,f7.1,4x, &
           &             ''drift:'',d10.2,3x,''Tbath :'',f6.1)')  &
           &             bave,bavt,slope,dat%tsoll
      end if
    else
      dat%iblock = dat%iblock+1
      dat%blocke(dat%iblock) = epot
      dat%blockt(dat%iblock) = temp
    end if

    return

  contains
    subroutine regress(n1,n2,rege,slope)
      implicit none
      real(wp) :: rege(*),slope
      integer :: n1,n2,n
      real(wp) :: sx,sy,sxx,sxy,x
      integer :: i,j,k,l,ich,och,io

      n = n2-n1+1
      sx = 0.0_wp
      sy = 0.0_wp
      sxy = 0.0_wp
      sxx = 0.0_wp
      x = 0.0_wp
      do i = n1,n2
        x = x+1.0_wp
        sx = sx+x
        sy = sy+rege(i)
        sxx = sxx+x**2
        sxy = sxy+x*rege(i)
      end do

      slope = (dble(n)*sxy-sx*sy)/(dble(n)*sxx-sx*sx)
      return
    end subroutine regress

    subroutine blocksd(n,nbl,ebl,tbl,esd,tsd)
      implicit none
      integer :: n,nbl
      real(wp) :: ebl(nbl),tbl(nbl),dum,av,esd,tsd
      integer :: i

      dum = 0.0_wp
      do i = 1,nbl
        dum = dum+ebl(i)
      end do
      av = dum/dble(nbl)
      esd = av

      dum = 0.0_wp
      do i = 1,nbl
        dum = dum+tbl(i)
      end do
      av = dum/dble(nbl)
      tsd = av

      return
    end subroutine blocksd
  end subroutine u_block

!========================================================================================!
! subroutines wrmdrestart & rdmdrestart
! write a file containing coordinates and velocities to restart the simulation
  subroutine wrmdrestart(mol,dat,velo,realtime_fs)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(in) :: velo(3,mol%nat)
    real(wp),intent(in) :: realtime_fs
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    character(len=256) :: atmp
    if (.not.allocated(dat%restartfile) .or. dat%restart) then
    !>--- we must not overwrite the user-provided restart file!
      write (atmp,'(a,i0,a)') 'crest_',dat%md_index,'.mdrestart'
    else
      atmp = dat%restartfile
    end if
    open (newunit=ich,file=trim(atmp)) !dat%restartfile)
    write (ich,*) realtime_fs
    do i = 1,mol%nat
      write (ich,'(6D22.14)') mol%xyz(1:3,i),velo(1:3,i)
    end do
    close (ich)
    return
  end subroutine wrmdrestart

  subroutine rdmdrestart(mol,dat,velo,fail,rtshift)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(inout) :: velo(3,mol%nat)
    logical,intent(out) :: fail
    real(wp),intent(out) :: rtshift
    real(wp) :: dum
    character(len=256) :: atmp
    integer :: i,j,k,l,ich,och,io
    logical :: ex

    fail = .false.
    rtshift = 0.0_wp
    if (allocated(dat%restartfile)) then
      inquire (file=dat%restartfile,exist=ex)
    end if
    if (dat%restart.and.ex) then
      open (newunit=ich,file=dat%restartfile)
      do
        read (ich,*,iostat=io) dum
        if (io < 0) exit
        if(dum > 0.0_wp) rtshift = dum 
        do i = 1,mol%nat
          read (ich,'(a)',iostat=io) atmp
          if (io < 0) exit
          read (atmp,'(6D22.14)',iostat=io) mol%xyz(1:3,i),velo(1:3,i)
          if (io /= 0) exit
        end do
        exit
      end do
      if (io /= 0) then
        write (0,*) 'failed to read MD restart file.'
        fail = .true.
      end if
      close (ich)
    else
      fail = .true.
    end if
    if(.not.fail)then
      write (*,'(1x,a,8x,l)') 'read restart file  :',.not.fail
    endif

    return
  end subroutine rdmdrestart

!========================================================================================!
! subroutine mdinitu
! initialize MD velocities, either by reading a restart file
! or by setting them randomly
  subroutine mdinitu(mol,dat,velo,mass,Ekin,pr)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(inout) :: velo(3,mol%nat)
    real(wp),intent(in) :: mass(mol%nat)
    real(wp),intent(in) :: Ekin
    logical,intent(in) :: pr
    real :: x(3),ranf
    integer :: n
    real(wp) :: eperat,v,f,t,edum,f2
    integer,allocatable :: iseed(:)
    logical :: newvelos
    integer :: i

    newvelos = .true.

    !>--- newly initialized
    if (newvelos) then
      if (dat%samerand) then
        call random_seed(size=n)
        allocate (iseed(n),source=1)
        call random_seed(put=iseed)
      else
        call random_seed()
      end if
      eperat = Ekin/(3.0_wp*float(mol%nat))
      do i = 1,mol%nat
        call random_number(x)
        f2 = 1.0_wp
        if (mol%at(i) .eq. 1) f2 = 2.0_wp
        v = sqrt(2*eperat/mass(i))
        f = 1.0_wp
        if (x(1) .gt. 0.5_wp) f = -1.0_wp
        velo(1,i) = v*f*f2
        f = 1.0_wp
        if (x(2) .gt. 0.5_wp) f = -1.0_wp
        velo(2,i) = v*f*f2
        f = 1.0d0
        if (x(3) .gt. 0.5_wp) f = -1.0_wp
        velo(3,i) = v*f*f2
      end do
    end if
    call ekinet(mol%nat,velo,mass,edum)
    t = edum/(0.5_wp*3.0_wp*float(mol%nat)*0.316681534524639E-05)
    return
  end subroutine mdinitu

!========================================================================================!
! subroutine thermostating
! helper routine to re-scale velocities,
! i.e., thermostating

  subroutine thermostating(mol,dat,t,scal)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(in) :: t
    real(wp),intent(out) :: scal
    integer :: i,j,k,l,ich,och,io

    scal = 1.0_wp

    if (.not.dat%thermostat) return

    select case (trim(dat%thermotype))
    case ('berendsen')
      scal = dsqrt(1.0d0+(dat%tstep/dat%thermo_damp) &
                  &     *(dat%tsoll/t-1.0_wp))
    case default 
    !>-- (no scaling, other thermostats require special implementation)
      scal = 1.0_wp
    end select

    return
  end subroutine thermostating

  subroutine thermostatprint(dat,pr)
    implicit none
    type(mddata) :: dat
    logical,intent(in) :: pr
    integer :: i,j,k,l,ich,och,io

    if (.not.pr) return
    if (dat%thermostat) then
      select case (trim(dat%thermotype))
      case ('berendsen')
        write (*,'('' thermostat         :'',1x,a  )') trim(dat%thermotype)
      case default !>-- (also berendsen thermostat)
        write (*,'('' thermostat         :'',1x,a  )') 'berendsen'
      end select
    else
      write (*,'('' thermostat         :'',1x,a  )') 'OFF'
    end if

    return
  end subroutine thermostatprint

!========================================================================================!
! subroutine zeroz
! remove z-directional acceleration of selected atoms
  subroutine zeroz(nat,acc,apply)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(inout) :: acc(3,nat)
    logical,intent(in) :: apply(nat)
    integer :: i
    do i = 1,nat
      if (apply(i)) acc(3,i) = 0.0_wp
    end do
    return
  end subroutine zeroz

!========================================================================================!
! subroutine rmrottr
! some MATHs to remove translational and rotational velocities
  subroutine rmrottr(natoms,mass,vel_atom,c)
    implicit none
    integer :: natoms
    real(wp) :: vel_atom(3,natoms),c(3,natoms),mass(natoms)
    real(wp) :: rlm(3),ram(3),omega(3)
    real(wp) :: ixx,iyy,izz,ixy,ixz,iyz,dummy
    real(wp) :: fixx,fiyy,fizz,fixy,fixz,fiyz,COM(3)
    real(wp) :: inertia(3,3),angmom(3)
    real(wp) :: tmass
    integer :: i

    rlm = 0.0
    ram = 0.0
    call centerofmass(natoms,c,mass,tmass,COM)
    angmom = 0.0
    do i = 1,natoms
      c(1,i) = c(1,i)-COM(1)
      c(2,i) = c(2,i)-COM(2)
      c(3,i) = c(3,i)-COM(3)
      angmom(1) = angmom(1)+mass(i)*(c(2,i)*vel_atom(3,i)-&
         &                                      c(3,i)*vel_atom(2,i))
      angmom(2) = angmom(2)+mass(i)*(c(3,i)*vel_atom(1,i)-&
         &                                      c(1,i)*vel_atom(3,i))
      angmom(3) = angmom(3)+mass(i)*(c(1,i)*vel_atom(2,i)-&
         &                                      c(2,i)*vel_atom(1,i))
    end do
    ixx = 0.0
    iyy = 0.0
    izz = 0.0
    ixy = 0.0
    ixz = 0.0
    iyz = 0.0
    do i = 1,natoms
      ixx = ixx+mass(i)*(c(2,i)*c(2,i)+c(3,i)*c(3,i))
      iyy = iyy+mass(i)*(c(3,i)*c(3,i)+c(1,i)*c(1,i))
      izz = izz+mass(i)*(c(1,i)*c(1,i)+c(2,i)*c(2,i))
      ixy = ixy-mass(i)*c(1,i)*c(2,i)
      ixz = ixz-mass(i)*c(1,i)*c(3,i)
      iyz = iyz-mass(i)*c(2,i)*c(3,i)
    end do
    inertia(1,1) = ixx
    inertia(2,2) = iyy
    inertia(3,3) = izz
    inertia(1,2) = ixy
    inertia(2,1) = ixy
    inertia(1,3) = ixz
    inertia(3,1) = ixz
    inertia(2,3) = iyz
    inertia(3,2) = iyz
    call dmatinv(inertia,3,3,dummy)
    omega = matmul(inertia,angmom)
    do i = 1,natoms
      rlm(1) = rlm(1)+mass(i)*vel_atom(1,i)
      rlm(2) = rlm(2)+mass(i)*vel_atom(2,i)
      rlm(3) = rlm(3)+mass(i)*vel_atom(3,i)
    end do
    do i = 1,natoms
      ram(1) = (omega(2)*c(3,i)-omega(3)*c(2,i))
      ram(2) = (omega(3)*c(1,i)-omega(1)*c(3,i))
      ram(3) = (omega(1)*c(2,i)-omega(2)*c(1,i))

      vel_atom(1,i) = vel_atom(1,i)-rlm(1)/tmass-ram(1)
      vel_atom(2,i) = vel_atom(2,i)-rlm(2)/tmass-ram(2)
      vel_atom(3,i) = vel_atom(3,i)-rlm(3)/tmass-ram(3)

      c(1,i) = c(1,i)+COM(1)
      c(2,i) = c(2,i)+COM(2)
      c(3,i) = c(3,i)+COM(3)
    end do

  contains
    subroutine centerofmass(natoms,c,mass,totmass,COM)
      implicit none
      integer natoms
      real(wp) :: c(3,natoms),totmass,mass(natoms),COM(3)
      integer :: i,j

      COM(1) = 0.0
      COM(2) = 0.0
      COM(3) = 0.0

      totmass = 0
      do i = 1,natoms
        totmass = totmass+mass(i)
        COM(1) = COM(1)+mass(i)*c(1,i)
        COM(2) = COM(2)+mass(i)*c(2,i)
        COM(3) = COM(3)+mass(i)*c(3,i)
      end do

      COM(1) = COM(1)/totmass
      COM(2) = COM(2)/totmass
      COM(3) = COM(3)/totmass
      return
    end subroutine centerofmass

    subroutine dmatinv(a,ldm,n,d)
      implicit none
      integer,intent(in) :: ldm,n
      real(wp),intent(out) :: d
      real(wp),intent(inout) :: a(ldm,*)
      integer :: i,j,k,l(n),m(n)
      real(wp) :: biga,temp
      real(wp),parameter :: tol = 1.0d-12
      !
      d = 1.0_wp
      !
      do k = 1,n
        l(k) = k
        m(k) = k
        biga = a(k,k)
        do j = k,n
          do i = k,n
            if (abs(biga) .lt. abs(a(j,i))) then
              biga = a(j,i)
              l(k) = i
              m(k) = j
            end if
          end do
        end do
        j = l(k)
        if (j .gt. k) then
          do i = 1,n
            temp = -a(i,k)
            a(i,k) = a(i,j)
            a(i,j) = temp
          end do
        end if
        i = m(k)
        if (i .gt. k) then
          do j = 1,n
            temp = -a(k,j)
            a(k,j) = a(i,j)
            a(i,j) = temp
          end do
        end if
        if (abs(biga) .lt. tol) then
          d = 0.0_wp
          return
        end if
        do i = 1,n
          if (i .ne. k) a(k,i) = a(k,i)/(-biga)
        end do
        do i = 1,n
          do j = 1,n
            if (i .ne. k) then
              if (j .ne. k) a(j,i) = a(k,i)*a(j,k)+a(j,i)
            end if
          end do
        end do
        do j = 1,n
          if (j .ne. k) a(j,k) = a(j,k)/biga
        end do
        d = max(-1.0d25,min(1.0d25,d))
        d = d*biga
        a(k,k) = 1.0_wp/biga
      end do
      !
      k = n
      do
        !
        k = k-1
        if (k .le. 0) exit
        i = l(k)
        if (i .gt. k) then
          do j = 1,n
            temp = a(k,j)
            a(k,j) = -a(i,j)
            a(i,j) = temp
          end do
        end if
        j = m(k)
        if (j .gt. k) then
          do i = 1,n
            temp = a(i,k)
            a(i,k) = -a(i,j)
            a(i,j) = temp
          end do
        end if
      end do
      !
      return
    end subroutine dmatinv
  end subroutine rmrottr

!========================================================================================!
  subroutine md_init_mtd(mol,dat,pr)
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    logical :: pr
    integer :: i

    if (dat%simtype .ne. type_mtd) return
    if (dat%npot < 1) return

    do i = 1,dat%npot
      call mtd_ini(mol,dat%mtd(i),dat%tstep,dat%length_ps,pr)
    end do

    return
  end subroutine md_init_mtd

!========================================================================================!
  subroutine md_add_mtd(self,mtd)
    implicit none
    class(mddata) :: self
    type(mtdpot) :: mtd
    type(mtdpot),allocatable :: mtdtmp(:)
    integer,allocatable :: cvtmp(:)
    integer :: i,j,k

    k = self%npot+1
    allocate (mtdtmp(k))
    allocate (cvtmp(k))
    if (k > 1) then
      do i = 1,k-1
        mtdtmp(i) = self%mtd(i)
        cvtmp(i) = self%cvtype(i)
      end do
    end if
    mtdtmp(k) = mtd
    cvtmp(k) = mtd%mtdtype
    call move_alloc(mtdtmp,self%mtd)
    call move_alloc(cvtmp,self%cvtype)
    self%npot = k

    !> if a metadynamics potential was added for the first time
    !> change MD runtype to MTD accordingly
    if (self%simtype == type_md) then
      self%simtype = type_mtd
    end if

    return
  end subroutine md_add_mtd

!========================================================================================!
  subroutine md_update_mtd(mol,dat,calc,pr)
!*********************************************
!* Update the collective variables for each
!* metadynamics potential saved in mddata
!*********************************************
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    type(calcdata) :: calc
    logical :: pr
    integer :: i

    if (dat%simtype .ne. type_mtd) return
    if (dat%npot < 1) return

    do i = 1,dat%npot
      select case (dat%cvtype(i))
      case (cv_rmsd,cv_rmsd_static)
        call cv_dump(mol,dat%mtd(i),0.0_wp,pr)
      case default
        cycle
      end select
    end do

    return
  end subroutine md_update_mtd

!========================================================================================!
  subroutine md_calc_mtd(mol,dat,epot,grd,grdmtd,pr)
!***********************************************
!* Calculate energy and gradient contributions
!* from metadynamics potentials and add them to
!* the current total energy and gradient
!***********************************************
!$  use omp_lib
    implicit none
    type(coord) :: mol
    type(mddata) :: dat
    real(wp),intent(inout) :: epot
    real(wp),intent(inout) :: grd(3,mol%nat)
!>--- LOCAL EXTENSION: scratch buffer, allocated ONCE per MD run by dynamics().
!> It used to be a local allocatable here, allocated and freed on every MD step
!> and both wrapped in an !$omp critical (which protects nothing: the array is
!> private to this call either way).
    real(wp),intent(inout) :: grdmtd(3,mol%nat)
    logical :: pr
    integer :: i
    real(wp) :: emtd

    if (dat%simtype .ne. type_mtd) return
    if (dat%npot < 1) return

    do i = 1,dat%npot
      call calc_mtd(mol,dat%mtd(i),emtd,grdmtd)
      epot = epot+emtd
      grd = grd+grdmtd
    end do

    return
  end subroutine md_calc_mtd

!========================================================================================!

  subroutine md_defaults_fallback(self)
!********************************************
!* Check if selected parameters that are
!* necessary to be able to run an MD
!* have been set, and restore default
!* values if not so
!********************************************
    implicit none
    class(mddata) :: self
    real(wp) :: dum

    if (self%length_ps <= 0.0_wp) then
      !> total runtime in ps
      self%length_ps = 20.0_wp
    end if
    if (self%tstep <= 0.0_wp) then
      !> time step in fs
      self%tstep = 1.0_wp
    end if
    if (self%length_steps <= 0) then
      !> simulation steps
      self%length_steps = nint(self%length_ps*1000.0_wp/self%tstep)
    end if
    if (self%tsoll <= 0.0_wp) then
      !> target temperature
      self%tsoll = 298.15_wp
    end if

    if (self%dumpstep <= 0.0_wp) then
      !> dump frequency in fs
      self%dumpstep = 1000.0_wp
    end if
    if (self%sdump <= 0) then
      !> trajectory structure dump every x steps
      dum = max(1.0_wp, (self%dumpstep/self%tstep))
      self%sdump = nint(dum)
    end if

    if (self%md_hmass <= 0.0_wp) then
      !> hydrogen mass
      self%md_hmass = 1.00794075_wp
    end if

    if (self%shake) then
      !> SHAKE, if turned on but no mode selected
      if (self%shk%shake_mode == 0) then
        self%shk%shake_mode = 2 !> all bonds
      end if
    end if

    !> block length (for average analysis and restart dump)
    if(self%blockl <= 0 )then
      self%blockl = min(5000,idint(5000.0_wp/self%tstep))
    endif
    self%maxblock = nint(self%length_steps/float(self%blockl))

  end subroutine md_defaults_fallback
!========================================================================================!
!========================================================================================!
!> LOCAL EXTENSION: turn the bond constraints that were read from the
!> constraint input file (--cinp / TOML "constraints") into RIGID SHAKE bonds.
!>
!> WHY THIS EXISTS
!> A bond constraint normally enters CREST as a HARMONIC RESTRAINT
!>   V(r) = 1/2 k (r - r_ref)^2        (calculator.F90 -> calc_constraint ->
!>                                      constraints.f90:bond_constraint)
!> whose force constant k is given in Eh/Bohr^2.  That is a SPRING, not a
!> rigid constraint: it introduces a vibrational mode with
!>   T = 2*pi*sqrt(mu/k).
!> The default force constants in this tree are 0.01 (fcdefault) and 0.05
!> (xtb-style "force constant="), i.e. a k of 3.5 is 70-350x the default and
!> about 12x stiffer than a real C-C bond.  For C-C that gives T ~ 8.5 fs, for
!> C-H ~ 3.3 fs, while the MD runs with a default time step of 5 fs.  Stable
!> integration needs dt <~ T/10, so the trajectory diverges and the only
!> visible symptom is "SHAKE did not converge" / "Trial MTD did not converge".
!>
!> A rigid constraint has no force constant at all, hence no frequency and no
!> extra restriction on the time step.  That is what a user asking to "fix the
!> reacting bond of a transition state" actually wants.
!>
!> The harmonic restraint is deliberately LEFT IN PLACE: SHAKE pins r at
!> r_ref exactly, so (r - r_ref) is zero to within tolshake and the spring
!> contributes nothing during an MD.  It does however remain active for the
!> geometry OPTIMIZATIONS and CREGEN, where SHAKE is not applied at all and
!> the spring is the only thing holding the bond.
!========================================================================================!
  subroutine crest_hard_bonds_to_shake(shk,n,atms,dist,pr)
    implicit none
    type(shakedata),intent(inout) :: shk
    integer,intent(in) :: n
    integer,intent(in) :: atms(2,n)   !> 1-based atom pairs
    real(wp),intent(in) :: dist(n)    !> target distances in BOHR
    logical,intent(in) :: pr
    integer :: i,j,m
    integer :: k,nu,ndof,ip
    integer,allocatable :: uat(:)
    logical :: dup,found

    if (n < 1) return
!>--- idempotent: if the list is already there, do not touch it
    if (shk%nusr > 0) return

    if (allocated(shk%conslistu)) deallocate (shk%conslistu)
    if (allocated(shk%distconsu)) deallocate (shk%distconsu)
!> NOTE: these two must NOT share one allocate statement.  A single
!> "allocate(a,b,source=0)" makes the source-expr conformable-check apply to
!> BOTH objects, and conslistu is integer while distconsu is real(wp) - gfortran
!> rejects that with "Type of entity ... incompatible with source-expr".
    allocate (shk%conslistu(2,n),source=0)
    allocate (shk%distconsu(n),source=0.0_wp)

    m = 0
    do i = 1,n
      if (atms(1,i) < 1.or.atms(2,i) < 1) cycle
      if (atms(1,i) == atms(2,i)) cycle
!>--- the same bond twice would be solved twice by SHAKE - harmless but
!> pointless, and it inflates nshake (and hence the reported DOF count).
      dup = .false.
      do j = 1,m
        if ((shk%conslistu(1,j) == atms(1,i).and.shk%conslistu(2,j) == atms(2,i)).or. &
        &   (shk%conslistu(1,j) == atms(2,i).and.shk%conslistu(2,j) == atms(1,i))) then
          dup = .true.
          exit
        end if
      end do
      if (dup) cycle
      m = m+1
      shk%conslistu(1,m) = atms(1,i)
      shk%conslistu(2,m) = atms(2,i)
      shk%distconsu(m) = dist(i)
    end do

    if (m < 1) then
      deallocate (shk%conslistu,shk%distconsu)
      return
    end if
    shk%nusr = m

    if (pr) then
!>--- LOCAL CHANGE (2026-10-04): no leading blank line - the caller already
!> ends the coordinate block, and two empty lines in a row are just noise.
      write (stdout,'(1x,a,i0,a)') 'rigid (SHAKE) bonds taken from the constraint input: ',m,' bond(s)'
      do j = 1,m
        write (stdout,'(5x,a,i0,a,i0,3x,a,f9.5,a)') 'bond ',shk%conslistu(1,j),' - ', &
        &  shk%conslistu(2,j),'fixed at ',shk%distconsu(j)*bohr,' Angstroem'
      end do
!>--- redundancy diagnostics.  When the list holds *every* pair of a set of nu
!> atoms (which is what "atoms: a,b,c,..." in the constraint file produces), that
!> subset can only supply 3*nu-6 internal degrees of freedom.  Past that the
!> constraint set is CONSISTENT (every distance came from a real geometry) but
!> LINEARLY DEPENDENT, and SHAKE then has to project onto a dependent set - the
!> classic cause of "SHAKE did not converge!".  Say so up front instead of
!> letting the run die 20 ps later.
!> LOCAL CHANGE (2026-10-04): only the "note:" line and the actionable hint are
!> printed.  The two lines that spelled out the DOF arithmetic in between were
!> removed - ndof is still computed, it is what the test below triggers on.
      allocate (uat(2*m))
      nu = 0
      do i = 1,m
        do k = 1,2
          ip = shk%conslistu(k,i)
          found = .false.
          do j = 1,nu
            if (uat(j) == ip) then
              found = .true.
              exit
            end if
          end do
          if (.not.found) then
            nu = nu+1
            uat(nu) = ip
          end if
        end do
      end do
      deallocate (uat)
      if (nu >= 3) then
        ndof = 3*nu-6
        if (m > ndof) then
          write (stdout,'(5x,a,i0,a,i0,a)') &
          &  'note: these bonds tie ',nu,' atoms together with ',m,' distances'
          write (stdout,'(5x,a)') &
          & '      If SHAKE stalls, narrow the atom list or fall back to -softcons.'
        end if
      end if
      write (stdout,*)
    end if
    return
  end subroutine crest_hard_bonds_to_shake
!========================================================================================!
!========================================================================================!
end module dynamics_module
