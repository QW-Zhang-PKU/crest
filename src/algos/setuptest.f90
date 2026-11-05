!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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
!================================================================================!

!> The following contains routines that are run in the automatic
!> setup of CREST to check for consistent settings.
!> First is the trialMD calculator that can automatically adjust the timestep
!> Second is the initial geometry optimization with topology check

!========================================================================================!
!========================================================================================!
!> Routines for checking metadynamics settings (timestep & shake, mostly)
!========================================================================================!
!========================================================================================!
subroutine trialMD_calculator(env)
!*****************************************************
!* subroutine trialMD_calculator
!* Runs a short metadynamics simulation and, in case
!* the simulation terminates with error, tries to
!* adjust the time step and tries again. Up to a
!* maximum of 6 tries.
!*****************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use dynamics_module
  use strucrd
  use iomod
  implicit none

  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA

  integer :: io,io2,Vdumpfreq,counter
  real(wp) :: prefac,alpha,length
  real(wp) :: rtime,tstep
  real(wp) :: trial_md_ps
!>--- Vbias (RMSD-CV) dump interval of the trial MTD, in fs (see the note where
!> it is set - it no longer scales with the trial length).
  real(wp) :: trial_vbdump_fs
  integer :: shakemode,maxiter
  integer :: nsteps   !> trial MTD length in steps (for the settings printout)
!>--- LOCAL EXTENSION: rigid (SHAKE) bonds coming from --cinp must survive the
!> automatic restart chain - turning SHAKE off would silently drop them.
  logical :: hardcons

  character(len=512) :: thispath
  character(len=80)  :: trajectory,atmp
  character(len=64)  :: cenv
!>--- where the trial-MTD length came from (ORCA default / internal default /
!> CREST_TRIAL_MD_PS) - printed in the settings block so the log is
!> self-explanatory when the number differs from the documented default.
  character(len=64)  :: trial_src
!> the "length set by" value is printed into a fixed 17-column field of the
!> settings box below, so it must fit into 17 characters
  character(len=17)  :: srcval

  logical :: shakefail,mdfail,ex,pr
  type(coord) :: mol,molstart
  type(mddata) :: MD,MDSTART
  type(mtdpot) :: MTD
  type(timer) :: profiler

  type(calcdata) :: tmpcalc
  type(calcdata) :: calcstart
  real(wp) :: energy
  integer :: T,Tn
  logical,external :: crest_calc_is_subprocess
!>--- LOCAL EXTENSION: the trial-MTD length now depends on the calculator, so
!> this predicate is needed here as well (defined next to its sibling in
!> parallel.f90).
  logical,external :: crest_calc_is_orca
!> LOCAL EXTENSION: same reasoning as MDFILES in parallel.f90 - the trial MTD
!> appends to this trajectory every sdump steps, and a symlink would not survive
!> the rmrf()+makedir() below.  Relocated to CREST_SCRATCH_DIR when set.
  character(len=:),allocatable :: dirnam

!>--- OMP settings.
!> Upstream set this to 'min' (one core) "to simulate max parallelization",
!> i.e. the trial MTD always ran fully serial.  For an external MPI program
!> (ORCA) that meant: a serial gradient at every MD step, for a step that
!> brings no scientific benefit in the first place - it only calibrates the
!> time step and extrapolates the runtime.
!> Instead hand the trial EXACTLY the budget one MTD gets in the production
!> run ('share', i.e. T/nmetadyn cores).  The trial is then representative
!> and the estimator below no longer has to *assume* a speed-up factor.
  if (crest_calc_is_subprocess(env%calc)) then
    call new_ompautoset(env,'share',env%nmetadyn,T,Tn)
  else
    call new_ompautoset(env,'min',0,T,Tn)
  end if

  call getcwd(thispath)

!>--- Make a dir
  call crest_hot_dir(dirnam,'TRIALMD')
  call rmrf(dirnam)
  io = makedir(trim(dirnam))

!>--- set up mol and MD calculator
  prefac = 0.003d0*env%rednat !> Vbias k
  alpha = 0.5d0         !> Vbias alpha

!>--- Length of the trial MTD, in ps.
!> The trial only (a) checks that the MD survives with the chosen time step /
!> SHAKE settings and (b) provides the wall time that is extrapolated to the
!> production run, so it does not have to be as long as a real MTD.
!> Upstream used 1.0 ps; we cut it to a quarter of that.
!>
!> LOCAL EXTENSION: when CREST drives ORCA the default is halved again, to
!> 0.125 ps.  An ORCA gradient is a full SCF - tens of seconds to minutes - so
!> the trial is nothing but overhead; its only jobs are "does the MD survive
!> the time step" and "give me one wall time sample", and neither needs a long
!> trajectory.  At the default 4 fs step, 0.125 ps is still 31 MD steps, and the
!> failures that actually matter show up almost immediately: the real one that
!> motivated this work died in STEP 1.  Everything else (the internal
!> GFN0/GFN1/GFN2/GFN-FF methods) keeps 0.25 ps - there the trial runs in
!> seconds, so there is nothing to save and the longer sample gives a better
!> runtime estimate.
!>
!> TRADE-OFF, read before shortening further: a shorter trial samples less wall
!> time, so the extrapolated production estimate gets noisier, and an
!> instability that only develops after ~50 MD steps would no longer be caught.
!> 31 steps is ample for a step-1 or first-few-step failure but it is a
!> calibration, not a proof - if you ever see a trial pass and the production
!> MTD blow up, raise this back rather than shortening it more.
!>
!> Override without recompiling:   export CREST_TRIAL_MD_PS=0.5
!> The env variable wins over BOTH defaults, ORCA or not.
  if (crest_calc_is_orca(env%calc)) then
    trial_md_ps = 0.125_wp
    trial_src = 'ORCA default'
  else
    trial_md_ps = 0.25_wp
    trial_src = 'internal default'
  end if
  cenv = ''
  call getenv('CREST_TRIAL_MD_PS',cenv)
  if (len_trim(cenv) > 0) then
    read (cenv,*,iostat=io2) length
    if (io2 == 0 .and. length > 0.0_wp) then
      trial_md_ps = length
      trial_src = 'CREST_TRIAL_MD_PS'
    else
      write (stdout,'(1x,a)') 'Warning: could not read CREST_TRIAL_MD_PS, using the default'
    end if
  end if

  call env%ref%to(molstart)
  MDSTART = env%mddat         !> env%mddat should already be set up at this point
  hardcons = (MDSTART%shk%nusr > 0)  !> rigid bonds from --cinp (choose_settings)
  MDSTART%length_ps = trial_md_ps !> we want only a short simulaton
  MDSTART%length_steps = 0    !> to zero so it will be calculated automatically
  MDSTART%dumpstep = 20.0_wp  !> fs dump step to trajectory file (the trj is deleted anyway)
  MDSTART%sdump = 0           !> also zero to reset
  if (allocated(env%ref%wbo)) then  !> should be allocated from main program
    MDSTART%shk%wbo = env%ref%wbo
  else !> otherwise, obtain from scratch
    tmpcalc = env%calc
    mol = molstart
    tmpcalc%calcs(1)%rdwbo = .true. !> obtain WBOs
!>--- LOCAL CHANGE (2026-09-26): the trial MTD needs the WBOs for SHAKE, nothing
!> else.  rd_ORCA_wbo() reads them from ORCA.out, so requesting a gradient made
!> the child solve a CP-SCF/Z-vector (and the CPCM gradient term) whose result
!> was discarded - `grd` only appears in its declaration, the allocate, this
!> call and the deallocate.  singlepoint() returns the same WBOs and energy.
    call singlepoint(mol,tmpcalc,energy,io)
    call move_alloc(tmpcalc%calcs(1)%wbo, env%ref%wbo)
    call tmpcalc%reset()
    MDSTART%shk%wbo = env%ref%wbo
  end if
  !MDSTART%printstep = 10

  MTD%kpush = prefac
  MTD%alpha = alpha
  MTD%mtdtype = cv_rmsd
!>--- Vbias (RMSD-CV) dump interval of the TRIAL MTD, in fs.
!>
!> History: upstream inherited 550 fs, which is the interval that made sense for
!> the 1 ps trial it used to run.  This tree shortened the trial to 0.25 ps (and
!> to 0.125 ps for ORCA) and scaled the interval with it, purely so that a short
!> trial would not end up with ZERO bias updates.
!>
!> LOCAL EXTENSION (2026-09-23): it is now a plain 50 fs.  Scaling was the wrong
!> lever - the point of the interval is "how many times does the bias adapt",
!> and the answer should not depend on how long the trial happens to be.  At the
!> ORCA default (0.125 ps) 50 fs gives 2 bias updates; at 0.25 ps it gives 4.
!> Compare the old behaviour: 68.75 fs over 0.125 ps gave exactly 1, and 137.5 fs
!> over 0.25 ps also gave exactly 1 - i.e. the trial was barely exercising the
!> bias at all.
!>
!> The min() below keeps the original safety property without the scaling: the
!> dump interval may never exceed half the trial, so ANY trial still contains at
!> least two bias updates.  That matters because CREST_TRIAL_MD_PS can make the
!> trial arbitrarily short, and without this guard a 0.02 ps trial would run as
!> plain MD while still being reported as a metadynamics test.
!> NOTE: the value in fs is converted to whole MD steps downstream
!> (metadynamics_module.f90:126, dum2 = max(1, cvdump_fs/tstep)), so at 4 fs the
!> actual interval is 13 steps = 52 fs.  The banner reports the requested value.
!>
!> Override without recompiling:   export CREST_TRIAL_VBDUMP_FS=100
  trial_vbdump_fs = 50.0_wp
  cenv = ''
  call getenv('CREST_TRIAL_VBDUMP_FS',cenv)
  if (len_trim(cenv) > 0) then
    read (cenv,*,iostat=io2) length
    if (io2 == 0 .and. length > 0.0_wp) then
      trial_vbdump_fs = length
    else
      write (stdout,'(1x,a)') 'Warning: could not read CREST_TRIAL_VBDUMP_FS, using the default'
    end if
  end if
  MTD%cvdump_fs = min(trial_vbdump_fs,trial_md_ps*500.0_wp)

!>--- LOCAL EXTENSION: the trial MTD runs with pr=.false. (silent), so without
!> this the log gives no indication of progress at all - for a slow external
!> calculator the run looks hung.  Report at 25/50/75 %; there is deliberately
!> no 100 % mark (it would duplicate the runtime line printed right after).
!> Only MDSTART carries the flag; env%mddat (production MTDs) is untouched.
  MDSTART%progress = .true.
  MDSTART%progress_tag = 'Trial MTD'

  call MDSTART%add(MTD)
  calcstart = env%calc !> Save clean state before loop

  pr = .false. !> supress stdout printout of MD

!>--- Header
  write (stdout,*)
  call smallhead('Starting trial MTD to test settings')

!>--- LOCAL EXTENSION: print the settings that are actually going to be used.
!> Upstream runs the trial with pr=.false., i.e. with NO printout whatsoever, so
!> the log never showed the time step, the length or the Vbias parameters - yet
!> those are exactly the numbers one needs to judge the runtime estimate below.
  nsteps = nint(MDSTART%length_ps*1000.0_wp/max(MDSTART%tstep,1.0e-6_wp))
!> NOTE on the layout: this block is deliberately IDENTICAL to the production
!> banner  "::::: starting MTD n :::::"  printed by
!> parallel_md_block_printout() in src/algos/parallel.f90, so that the trial and
!> the production settings line up character for character in the log:
!>   2 leading blanks | '|' | 3 blanks | label left-justified in 21 columns |
!>   ':' | 18-column value field | '|'
!> i.e. the ':' sits in column 26 and the closing '|' in column 47.  Every
!> *numeric* value is written into a fixed 8-character field (f8.x / i8) - do
!> NOT use i0, it emits the integer directly behind the colon and breaks the
!> column of the float lines.  Text values get a 17-column field with one
!> leading blank.
  write (stdout,'(2x,a,f8.3,a)') &
  &  '|   MD simulation time   :',MDSTART%length_ps,' ps       |'
!>--- why it is that length.  Without this line a log showing 0.125 ps looks
!> like a mis-set default; with it, the origin is unambiguous.
  srcval = adjustl(trial_src)
  write (stdout,'(2x,a,a,a)') '|   length set by        : ',srcval,'|'
  write (stdout,'(2x,a,f8.1,a)') '|   timestep dt          :',MDSTART%tstep,' fs       |'
  write (stdout,'(2x,a,i8,a)')   '|   steps                :',nsteps,'          |'
  write (stdout,'(2x,a,f8.1,a)') '|   target T             :',MDSTART%tsoll,' K        |'
  if (MDSTART%shake) then
    if (MDSTART%shk%shake_mode == 2) then
      write (stdout,'(2x,a,a)') '|   SHAKE algorithm      :',' true (all bonds) |'
    else if (MDSTART%shk%shake_mode == 1) then
      write (stdout,'(2x,a,a)') '|   SHAKE algorithm      :',' true (H only)    |'
    else
      write (stdout,'(2x,a,a)') '|   SHAKE algorithm      :',' true (user bonds)|'
    end if
  else
    write (stdout,'(2x,a,a)') '|   SHAKE algorithm      :','false             |'
  end if
  write (stdout,'(2x,a,f8.4,a)') '|   Vbias prefactor (k)  :',MTD%kpush,' Eh       |'
!>--- the trial CV is cv_rmsd (MTD%mtdtype above), i.e. exactly the case in which
!> the production banner prints the bohr^-2 unit, so print it here too.
  write (stdout,'(2x,a,f8.4,a)') '|   Vbias exponent (α)   :',MTD%alpha,' bohr⁻²   |'
  write (stdout,'(2x,a,f8.1,a)') '|   dump interval(Vbias) :',MTD%cvdump_fs,' fs       |'
  write (stdout,'(2x,a,i8,a)')   '|   cores per job        :',Tn,'          |'
!>--- LOCAL EXTENSION: closing rule of the box, same width as the
!> "::::: starting MTD n :::::" headline of the production banner (45 columns).
  write (stdout,'(2x,a)') repeat(':',45)
  write (stdout,*)

!>--- Iterative loop, since it is also tested if the MD runs at all
  counter = 1
  maxiter = 6
  tstep = MDSTART%tstep
  shakemode = MDSTART%shk%shake_mode
  call profiler%init(maxiter)
  iterativ: do

!>--- Restore initial starting geometry
    mol = molstart
!>--- Restore clean calculation state
    env%calc = calcstart
!>--- Modify MD output trajectory
    MD = MDSTART
    MD%tstep = tstep
    MD%shk%shake_mode = shakemode
!>--- LOCAL EXTENSION: rigid bonds from --cinp need SHAKE even when the user
!> asked for -shake 0 (a typical reaction-path setup: no X-H constraints, only
!> the reacting bond).  init_shake() with shake_mode==0 then constrains exactly
!> the user bonds and nothing else.
    MD%shake = (shakemode > 0).or.hardcons
!>--- LOCAL EXTENSION: the trial MTD must NOT react to a "stop" file.
!> This loop below reads io /= 0 as "the MD blew up -> halve the time step and
!> try again", and after maxiter tries it gives up with status_safety, i.e. the
!> WHOLE crest run aborts.  A user who created a stop file would therefore get
!> "Trial MTD n did not converge!" six times and no results at all - the exact
!> opposite of what the stop file is for.  The trial MTD is short anyway.
    MD%stopcheck = .false.
    write (trajectory,'("crest_trial_md_",i0,".xyz")') counter
    MD%trajectoryfile = trim(dirnam)//sep//trim(trajectory)

!>--- And run it
    io = 1
    !================================!
    call profiler%start(counter)
    call dynamics(mol,MD,env%calc,pr,io)
    call profiler%stop(counter)
    !================================!

    if (io == 0) then
      write (atmp,'(1x,"Trial MTD ",i0," runtime (",f4.2," ps)")') counter,MD%length_ps
      call profiler%write_timing(stdout,counter,trim(atmp))
      exit iterativ

    else

      write (stdout,'(1x,"Trial MTD ",i0," did not converge!")') counter
      if (counter >= maxiter) then
        write (stdout,'(1x,"Automatic MD restart failed ",i0," times!")') counter
        write (stdout,'(1x,"Please try other settings manually.")')
        write (stdout,*)
        env%iostatus_meta = status_safety
        return
      end if
      counter = counter+1

!>--- LOCAL EXTENSION: with rigid constraints from --cinp SHAKE may NOT be
!> switched off, so there is nothing left to trade a shorter time step against.
      if (tstep <= 1.0d0.and.hardcons) then
        write (stdout,'(1x,"Time step is already at 1.0 fs and the rigid")')
        write (stdout,'(1x,"constraints from the input file must stay on.")')
        write (stdout,'(1x,"Please check the constrained bonds manually.")')
        write (stdout,*)
        env%iostatus_meta = status_safety
        return

      else if (tstep <= 1.0d0.and.shakemode > 0) then
        write (stdout,'(1x,"Time step can´t be reduced further!")')
        write (stdout,'(1x,"Trying to turn SHAKE off ...")')
        write (stdout,*)
        shakemode = 0
        cycle iterativ

      else if (tstep <= 1.0d0.and.shakemode == 0) then
        write (stdout,'(1x,"Automatic MTD settings check failed!")')
        write (stdout,'(1x,"Please try other settings manually.")')
        write (stdout,*)
        env%iostatus_meta = status_safety
        return
      end if

      !> don't reduce the timestep below 1 fs automatically
      tstep = max(1.0d0, (tstep-1.0d0))
      write (stdout,'(1x,"Reducing the time step to ",f4.1," fs and trying again ...")') tstep
      write (stdout,*)
    end if
!>--- End loop
  end do iterativ

!>--- transfer final settings to global settings
  env%mdstep = MD%tstep
  env%mddat%tstep = MD%tstep
  env%shake = shakemode
  env%mddat%shk%shake_mode = shakemode
  env%mddat%shake = (shakemode > 0).or.hardcons

!>--- If we succeded and exited the above loop, estimate runtime for multiple MTDs
!> NOTE: rtime was measured with Tn cores (see the 'share' OMP setting above),
!> so it already IS a per-job runtime - hand Tn to the estimator.
  rtime = profiler%get(counter)
  call timeestimate(rtime,env%mdtime,env%nmetadyn,env%threads,Tn,trial_md_ps)

  return
!=========================================================================================!
contains
!=========================================================================================!
  subroutine timeestimate(rtime,mdtime,nmetadyn,threads,perjob,trialps)
!******************************************
!* calculate and print the time estimate
!*
!* rtime  : measured wall time of the trial MTD, which ran for "trialps"
!*          picoseconds with "perjob" cores - the very same budget one MTD
!*          gets in the production run ('share' mode) - so it must NOT be
!*          scaled down by a fictional speed-up factor the way the old
!*          single-thread version was.
!* trialps: length of the trial MTD in ps.  Must be divided out, because the
!*          estimate below wants "time per ps", not "time per trial".
!******************************************
    implicit none
    real(wp)  :: mdtime
    integer :: nmetadyn
    integer :: threads
    integer :: perjob
    real(wp)  :: rtime
    real(wp)  :: trialps
    integer   :: parallel,waves
    character(len=80) :: stime,stime2

!>--- rtime is the time for "trialps" ps; scale it to the production MTD length
    rtime = rtime*(mdtime/max(trialps,1.0e-6_wp))
    call crest_time2string(rtime,stime)
    write (stime2,'(f14.1)') mdtime
    stime2 = adjustl(stime2)
    write (stdout,'(1x,a,a,a,i0,a,a)') 'Estimated runtime for one MTD (', &
    &                      trim(stime2),' ps) on ',perjob,' cores: ',trim(stime)

    !>--- how many MTDs fit into the budget at the same time when each of them
    !> gets "perjob" cores; anything beyond that runs in a later wave.
    !> The min() matters: the production phase never runs more than "nmetadyn"
    !> MTDs even if the budget would allow more by the per-job accounting.
    parallel = max(1,min(nmetadyn,threads/max(1,perjob)))
    waves = max(1,int(ceiling(real(nmetadyn)/real(parallel))))

    rtime = rtime*real(waves)
    call crest_time2string(rtime,stime)
    if (waves > 1) then
      write (stdout,'(1x,a,i0,a,i0,a,i0,a,i0,a,a)') &
      &   'Estimated runtime for a batch of ',nmetadyn,' MTDs (',parallel, &
      &   ' at a time, ',waves,' waves) on ',threads,' threads: ',trim(stime)
    else
      write (stdout,'(1x,a,i0,a,i0,a,a)') &
      &   'Estimated runtime for a batch of ',nmetadyn, &
      &   ' MTDs (all parallel) on ',threads,' threads: ',trim(stime)
    end if

  end subroutine timeestimate

!========================================================================================!
end subroutine trialMD_calculator

!=========================================================================================!
!> Convert a number of seconds into a short human readable string.
!>
!> This used to be a contained subroutine of trialMD_calculator, but the initial
!> geometry optimization reports its runtime as well, so it lives at file scope
!> now.  It must NOT be called "time2string": legacy_algos/trialmd_legacy.f90
!> already defines that name at file scope and the linker would then see a
!> duplicate symbol.
!=========================================================================================!
subroutine crest_time2string(rtime,str)
  use crest_parameters,only: wp
  implicit none
  real(wp) :: rtime
  real(wp) :: hours,minutes,seconds
  character(len=*) :: str

  hours = aint(rtime/3600d0)
  minutes = aint((rtime-(3600d0*hours))/60d0)
  seconds = anint(rtime-(3600d0*hours+minutes*60d0))

  if (hours > 0.0d0) then
    write (str,'(i0,1x,a,1x,i0,1x,a,1x,i0,1x,a)') &
    & nint(hours),'h',nint(minutes),'min',nint(seconds),'sec'
  else if (minutes > 0.0d0) then
    write (str,'(i0,1x,a,1x,i0,1x,a)') &
    & nint(minutes),'min',nint(seconds),'sec'
  else
    write (str,'(i0,1x,a)') nint(seconds),'sec'
  end if

  return
end subroutine crest_time2string

!========================================================================================!
!========================================================================================!
!> Routines to check topology behaviour upon optimization of input
!========================================================================================!
!========================================================================================!
subroutine trialOPT_calculator(env)
!**********************************************************
!* subroutine trialOPT_calculator
!* A new calculator-based implementation of xtbopt_legacy
!* Performs a geometry optimization of the structure
!* saved to env%ref and checks for changes in the topology
!*
!* requires to have env%ref and env%calc initialized
!**********************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use optimize_module
  use strucrd
  use iomod
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  !> LOCAL
  type(coord) :: mol,molopt
  type(calcdata) :: tmpcalc
  integer :: io,T,Tn
  real(wp) :: energy
  real(wp),allocatable :: grd(:,:)
  logical :: success,pr,wr
!>--- LOCAL EXTENSION: runtime / cycle count of the initial optimization
  type(timer) :: optprof
  real(wp) :: optrtime
  integer :: noptcycles
  character(len=80) :: stime

!>--- get all available threads
  call new_ompautoset(env,'max',0,T,Tn)

!>--- small header
  write (stdout,*)
  call smallhead('Initial Geometry Optimization')

!>--- setup
  call env%ref%to(mol)
  call env%ref%to(molopt)
  allocate(grd(3,mol%nat), source=0.0_wp)
  tmpcalc = env%calc  !> create copy of calculator
  tmpcalc%optlev = -1 !> set loose convergence thresholds 

!>--- perform geometry optimization
  pr = .false. !> stdout printout
  wr = .true.  !> write crestopt.xyz
  if(wr)then
    call remove('crestopt.xyz')
  endif 
!>--- LOCAL EXTENSION: time it.  Upstream prints nothing at all between the
!> "Initial Geometry Optimization" header and "Geometry successfully optimized.",
!> so with a slow external calculator the run looks hung and there is no number
!> to compare against when something is off.
  call optprof%init(1)
  call optprof%start(1)
  call optimize_geometry(mol,molopt,tmpcalc,energy,grd,pr,wr,io)
  call optprof%stop(1)
  optrtime = optprof%get(1)
  call crest_time2string(optrtime,stime)
!> one frame per optimization cycle (gradient evaluation) is appended to
!> crestopt.log, so counting them gives the number of cycles
  noptcycles = count_xyz_frames('crestopt.xyz',molopt%nat)
  if (noptcycles > 0) then
    write (stdout,'(1x,a,i0,a,a,a)') 'Initial geometry optimization: ', &
    &   noptcycles,' cycles, ',trim(stime),'   (crestopt.xyz)'
  else
    write (stdout,'(1x,a,a)') 'Initial geometry optimization: ',trim(stime)
  end if

!>--- check success
  success = (io == 0)
  call trialOPT_warning(env,molopt,success)
!>--- if the checks were successfull, env%ref is overwritten
  env%ref%nat = molopt%nat
  env%ref%at = molopt%at
  env%ref%xyz = molopt%xyz
  env%ref%etot = energy  

  deallocate(grd)
  return
!=========================================================================================!
!=========================================================================================!
  contains
!=========================================================================================!
!> Count the structures (frames) in an XYZ trajectory file.
!>
!> All optimization engines append exactly one frame per cycle to crestopt.xyz,
!> so the frame count is the number of optimization cycles that were needed.
!>   >0 : number of frames
!>   -1 : file missing / unreadable (caller then just prints the runtime)
!=========================================================================================!
  integer function count_xyz_frames(fname,nat) result(nfr)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    character(len=512) :: line
    integer :: ich,io,i
    logical :: ex

    nfr = -1
    if (nat <= 0) return
    inquire (file=trim(fname),exist=ex)
    if (.not.ex) return
    nfr = 0
    open (newunit=ich,file=trim(fname),status='old',action='read',iostat=io)
    if (io /= 0) then
      nfr = -1
      return
    end if
    do
      read (ich,'(a)',iostat=io) line   !> atom count
      if (io /= 0) exit
      read (ich,'(a)',iostat=io) line   !> comment line
      if (io /= 0) exit
      do i = 1,nat
        read (ich,'(a)',iostat=io) line
        if (io /= 0) exit
      end do
      if (io /= 0) exit
      nfr = nfr+1
    end do
    close (ich)
    return
  end function count_xyz_frames
!=========================================================================================!
end subroutine trialOPT_calculator

!========================================================================================!
subroutine trialOPT_warning(env,mol,success)
!*******************************************************
!* subroutine trialOPT_warning
!* Processes the trialOPT status and, if successfull,
!* overwrites env%ref with the opt. structure
!* If the checks fail, CREST is stopped by this routine
!*******************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use utilities
  implicit none
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol  !> optimized geometry to be checked
  logical,intent(in) :: success  !> status of geometry optimization
  integer :: ntopo
  integer,allocatable :: topo(:)
  logical,allocatable :: changed(:)
  integer(dp) :: i  !> revlin expects int64 input
  integer :: j,k,l
  logical :: tchange

  if (.not.success) then
    write (stdout,*)
    write (stdout,*) ' Initial geometry optimization failed!'
    write (stdout,*) ' Please check your input and, if present, crestopt.xyz.'
    call creststop(status_failed)
  end if
  write (stdout,*) 'Geometry successfully optimized.'
!---- if necessary, check if the topology has changed!
  if (allocated(env%ref%topo)) then
    ntopo = mol%nat*(mol%nat+1)/2
    allocate (topo(ntopo))
    allocate (changed(mol%nat),source=.false.)
    call quicktopo(mol%nat,mol%at,mol%xyz,ntopo,topo)
    tchange = .false.
    do i = 1,ntopo
      if (topo(i) .ne. env%ref%topo(i)) then
        tchange = .true.
        call revlin(i,k,l)
        changed(k) = .true.
        changed(l) = .true.
      end if
    end do
    if (tchange) then
      write (stdout,*)
      write (stdout,'(1x,a)') '*WARNING* Change in topology detected!'
      write (stdout,'(1x,a)') 'Topology change compared to the input affects atoms:'
      do i = 1,mol%nat
        if (changed(i)) then
          write (stdout,'(1x,i0,"(",a,")")',advance='no') i,trim(i2e(mol%at(i),'nc'))
        end if
      end do
      write (stdout,*)

!>--- either update the topology (see option B below)
      if (.not.env%reftopo) then
        write (stdout,'(1x,a)') 'Taking new topology as reference and continue ...'
        env%ref%topo = topo

!>--- or abort the run
      else
        write (stdout,*)
        call smallhead('* READ THE FOLLOWING CAREFULLY *')
        write (stdout,'(1x,a)') 'A topology change was seen in the initial geometry optimization.'
        write (stdout,'(1x,a,a,a)') 'This could be an artifact of the chosen theory level (e.g. xTB).'
        if (env%legacy) then
          write (stdout,'(1x,a)') 'You can check the optimization trajectory in the "xtbopt.log" file.'
        else
          write (stdout,'(1x,a)') 'You can check the optimization trajectory in the "crestopt.xyz" file.'
        end if
        write (stdout,'(1x,a)') 'Try either of these options:'
        write (stdout,'(/,4x,a)') 'A) Pre-optimize your input seperately and use the optimized'
        write (stdout,'(4x,a)') '   structure as input for CREST. (Only recommended if structure is intact)'
        write (stdout,'(/,4x,a)') 'B) Restart the same CREST call as before, but ignore the topology change'
        write (stdout,'(4x,a)') '   by using the "--noreftopo" keyword. (May produce artifacts)'
        write (stdout,'(/,4x,a)') 'C) Fix the initial input geometry by introducing bond length constraints'
        write (stdout,'(4x,a)') '   or by using a method with fixed topology (e.g. GFN-FF).'
        write (stdout,*)
        call creststop(status_safety)
      end if
    end if
  end if

  return
end subroutine trialOPT_warning

!========================================================================================!
!========================================================================================!
