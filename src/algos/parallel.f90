!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2023  Philipp Pracht
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

!> A collection of routines to set up OMP-parallel runs of MDs and optimizations.

!========================================================================================!
!========================================================================================!
!> Interfaces to handle optional arguments
!========================================================================================!
!========================================================================================!
module parallel_interface
!*******************************************************
!* module to load an interface to the parallel routines
!* mandatory to handle any optional input arguments
!*******************************************************
  implicit none
  interface
    subroutine crest_sploop(env,nat,nall,at,xyz,eread)
      use crest_parameters,only:wp,stdout,sep
      use crest_calculator
      use omp_lib
      use crest_data
      use strucrd
      use optimize_module
      use iomod,only:makedir,directory_exist,remove
      implicit none
      type(systemdata),intent(inout) :: env
      real(wp),intent(inout) :: xyz(3,nat,nall)
      integer,intent(in)  :: at(nat)
      real(wp),intent(inout) :: eread(nall)
      integer,intent(in) :: nat,nall
    end subroutine crest_sploop
  end interface

  interface
    subroutine crest_oloop(env,nat,nall,at,xyz,eread,dump,customcalc)
      use crest_parameters,only:wp,stdout,sep
      use crest_calculator
      use omp_lib
      use crest_data
      use strucrd
      use optimize_module
      use iomod,only:makedir,directory_exist,remove
      use crest_restartlog,only:trackrestart,restart_write_dummy
      implicit none
      type(systemdata),target,intent(inout) :: env
      real(wp),intent(inout) :: xyz(3,nat,nall)
      integer,intent(in)  :: at(nat)
      real(wp),intent(inout) :: eread(nall)
      integer,intent(in) :: nat,nall
      logical,intent(in) :: dump
      type(calcdata),intent(in),target,optional :: customcalc
    end subroutine crest_oloop
  end interface
end module parallel_interface

!========================================================================================!
!========================================================================================!
!> LOCAL EXTENSION: dynamic re-distribution of the core budget among the MD/MTD jobs
!> that are STILL running ("share the whole node").
!========================================================================================!
!========================================================================================!
!>
!> WHY IT IS NEEDED
!> ----------------
!> crest_search_multimd*() gives every job a FIXED slice of the -T budget:
!> new_ompautoset(env,'auto',nsim,T,Tn) computes Tn = T/nsim once, exports it as
!> OMP_NUM_THREADS (ompenvset), and ORCA turns that into "%pal nprocs Tn".
!> The slice is never revisited.  A job that finishes early therefore leaves its
!> slice idle for the REST of the round - with 4 MTDs on 48 cores, one early
!> finisher means 12 cores doing nothing, two means 24.  Since all MTDs of a
!> round run the same number of steps, the tail of the round is exactly the part
!> where this hurts.
!>
!> WHAT IT DOES
!> ------------
!> Whenever a job finishes, the jobs that are still running split the WHOLE
!> budget again:  4 -> 3 -> 2 -> 1 job  =>  12 -> 16 -> 24 -> 48 cores per job.
!>
!> HOW
!> ---
!> The one channel that reaches a job WHILE IT IS RUNNING is the
!> OMP_NUM_THREADS environment variable:
!>   * ORCA: orca_sc.f90 re-writes ORCA.in for EVERY singlepoint, and
!>           write_orca_input() -> orca_nprocs() -> OMP_NUM_THREADS decides the
!>           "%pal nprocs N" of that step.
!>   * xtb / turbomole / generic: the child inherits OMP_NUM_THREADS directly.
!> So the new share is published with ompenvset() and takes effect at the NEXT
!> external call of the still-running jobs.  A job that is currently inside its
!> external call keeps the ranks it already has - which is exactly right: the
!> change is a hand-over, never a resize of a live calculation.
!>
!> SCOPE
!> -----
!> ONLY for SUBPROCESS calculators.  An in-process calculator (tblite, gfnff,
!> ...) consumes its budget as a nested OpenMP team that is created INSIDE
!> dynamics(); its size cannot be changed from out here, so the flag stays off
!> for those and nothing at all changes for them.
!>
!> THREAD SAFETY
!> -------------
!> All state lives here and is touched only inside
!>     !$omp critical (crest_dynshare_lock)
!> so the counter can never be lost and the printout can never interleave.
!> The remaining concern is setenv() being called while OTHER threads fork/exec:
!> glibc's setenv() re-uses the storage of an existing entry when the new value
!> fits into it (only a longer value makes it replace the pointer, which is the
!> case that can race with a concurrent getenv()).  crest_dynshare_init()
!> therefore first publishes the WIDEST value the round can ever produce
!> (ds_total) and only then the actual one - every later value is <= ds_total,
!> so every later setenv() writes in place and never swaps a pointer.
!>
!> KNOBS
!> -----
!>   CREST_DYNAMIC_SHARE=0        switch the whole thing off (default: ON)
!>   CREST_MAX_CORES_PER_JOB=<n>  cap the share a single job may get.
!>                                0 (default) = no cap.
!>                                With ORCA/MPI this is worth setting: 48 MPI
!>                                ranks on a MINIX/PBE system can easily LOSE
!>                                against 12, because every SCF iteration then
!>                                pays an Allreduce of the Fock matrix.  16-24
!>                                is a sane starting point.
!>   CREST_ORCA_NPROCS=<n>        (existing knob) pins the ORCA rank count, so
!>                                it ALSO caps this mechanism - a warning is
!>                                printed when both are active.
!========================================================================================!
module crest_dynshare
  use crest_parameters,only:stdout
  implicit none
  private
  public :: crest_dynshare_init,crest_dynshare_release,crest_dynshare_reset
  public :: crest_dynshare_barmark
!>--- state shared by the OpenMP threads (only ever touched under the critical)
  logical :: ds_on = .false.       !> is the dynamic sharing active for this round?
  integer :: ds_njobs = 0          !> total number of jobs of this round (nsim)
  integer :: ds_nslots = 0         !> jobs running concurrently (= new_ompautoset's T)
  integer :: ds_total = 0          !> the whole core budget (= env%Threads, i.e. -T)
  integer :: ds_cap = 0            !> CREST_MAX_CORES_PER_JOB, 0 = unlimited
  integer :: ds_share = 0          !> cores per job right now
  integer :: ds_share0 = 0         !> cores per job at the start (restored by reset)
  integer :: ds_nfin = 0           !> jobs finished so far
!>--- LOCAL CHANGE (2026-10-04): the reallocation is reported ONCE per round, and
!> the report must not land in the middle of a progress bar.
!>  ds_announced - has the first reallocation been reported already?
!>  ds_baropen   - is stdout currently in the middle of a "|>xx.x%" line?
!>                 (set by crest_dynshare_barmark(), see crest_oloop_pr_progress)
  logical :: ds_announced = .false.
  logical :: ds_baropen = .false.
contains
!========================================================================================!
!> Call ONCE per round, right after new_ompautoset() and BEFORE the OMP region.
!> njobs  = nsim              (how many jobs this round has in total)
!> nslots = T                 (how many of them run at the same time)
!> total  = env%Threads       (the -T budget)
!> share0 = Tn                (what new_ompautoset() gave a single job)
!========================================================================================!
  subroutine crest_dynshare_init(njobs,nslots,total,share0)
    implicit none
    integer,intent(in) :: njobs,nslots,total,share0
    integer :: r,ival
    character(len=64) :: val
    ds_on = .false.
    ds_njobs = max(0,njobs)
    ds_nslots = max(1,nslots)
    ds_total = max(1,total)
    ds_nfin = 0
    ds_cap = 0
    ds_share0 = max(1,share0)
    ds_share = ds_share0
    ds_announced = .false.
    ds_baropen = .false.
!>--- only meaningful with more than one concurrent job, and only if the budget
!>    is actually split (a single job already owns everything)
    if (ds_nslots <= 1) return
!>--- off switch
    call getenv('CREST_DYNAMIC_SHARE',val)
!> a failed internal READ leaves the variable UNCHANGED, and Fortran's .and.
!> does not guarantee short-circuit evaluation - so ival must be defined before
!> the test below looks at it, otherwise an unset variable can read as 0 and
!> silently switch dynamic sharing off.
    ival = 0
    read (val,*,iostat=r) ival
    if (r == 0.and.ival == 0) return
!>--- optional cap
    call getenv('CREST_MAX_CORES_PER_JOB',val)
    read (val,*,iostat=r) ival
    if (r == 0) ds_cap = max(0,ival)
    if (ds_cap > 0) ds_share = min(ds_share,ds_cap)
    ds_on = .true.
!>--- pre-size the environment entries to the widest value this round can ever
!>    produce, so that every later setenv() re-writes them in place (see the
!>    THREAD SAFETY note at the top of this module).
    call ompenvset(ds_total)
    call ompenvset(ds_share)
    write (stdout,'(1x,a,i0,a,i0,a,a)') '[cores] dynamic sharing: ', &
    &  ds_njobs,' jobs, ',ds_share, &
    &  ' cores each - a job that finishes hands its', &
    &  ' cores to the rest'
    if (ds_cap > 0) then
      write (stdout,'(1x,a,i0,a)') '[cores] capped at ',ds_cap,' cores per job ($CREST_MAX_CORES_PER_JOB)'
    end if
!>--- CREST_ORCA_NPROCS wins over OMP_NUM_THREADS in orca_nprocs(), i.e. it pins
!>    the rank count and silently disables this mechanism.  Say so.
    call getenv('CREST_ORCA_NPROCS',val)
    if (len_trim(val) > 0) then
      write (stdout,'(1x,a)') '[cores] **WARNING** $CREST_ORCA_NPROCS is set - it pins the ORCA rank count'
      write (stdout,'(1x,a)') '        and therefore overrides the dynamic sharing above.'
    end if
    flush (stdout)
  end subroutine crest_dynshare_init

!========================================================================================!
!> Call ONCE per job, right after dynamics() has returned (whatever the reason).
!========================================================================================!
  subroutine crest_dynshare_release()
    implicit none
    integer :: nrunning,newshare,oldshare
    if (.not.ds_on) return
    !$omp critical (crest_dynshare_lock)
    ds_nfin = ds_nfin+1
!>--- how many jobs are ACTUALLY still running: with waves (nsim > nslots) a
!>    finishing job is immediately replaced by a queued one, so the number only
!>    drops once the queue is drained.
    nrunning = max(1,min(ds_nslots,ds_njobs-ds_nfin))
    newshare = max(1,ds_total/nrunning)
    if (ds_cap > 0) newshare = min(newshare,ds_cap)
    if (newshare /= ds_share) then
      oldshare = ds_share
      ds_share = newshare
!>--- the new value is always applied; it is only the *report* that is limited
!> to the first change (LOCAL CHANGE 2026-10-04).  A 172-job batch used to emit
!> seven consecutive lines saying the same sentence with a larger number each
!> time, all of them appended to the progress bar line.
      call ompenvset(ds_share)
      if (.not.ds_announced) then
        ds_announced = .true.
!>--- break the progress bar first: it is written with advance='no', so without
!> this the message would land on the same line as "|>90.1%" and the remaining
!> bar tokens would be appended after it.
        if (ds_baropen) then
          write (stdout,'()')
          ds_baropen = .false.
        end if
        write (stdout,'(1x,a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
        &  '[cores] ',ds_nfin,' of ',ds_njobs,' jobs done -> ',nrunning, &
        &  ' still running, ',ds_share,' cores/job (was ',oldshare,')'
        flush (stdout)
      end if
    end if
    !$omp end critical (crest_dynshare_lock)
  end subroutine crest_dynshare_release

!========================================================================================!
!> Call this from any progress printer that writes with advance='no', i.e. every
!> time stdout is deliberately left in the middle of a line.  Only then may
!> crest_dynshare_release() write a newline of its own before its message.
!> The flag is a plain logical: the write side runs on the master thread and the
!> read side sits in the critical section, so a torn read cannot happen.
!========================================================================================!
  subroutine crest_dynshare_barmark()
    implicit none
    if (.not.ds_on) return
    ds_baropen = .true.
  end subroutine crest_dynshare_barmark

!========================================================================================!
!> Call once after the OMP region: restore the share the round started with.
!> Every later stage calls new_ompautoset() again and would reset it anyway, but
!> the code between the two must not inherit a "one job owns everything" value.
!========================================================================================!
  subroutine crest_dynshare_reset()
    implicit none
    if (.not.ds_on) return
    call ompenvset(ds_share0)
    ds_on = .false.
  end subroutine crest_dynshare_reset
end module crest_dynshare

!========================================================================================!
!========================================================================================!
!> Helper to decide whether an additional (nested) OpenMP team within CREST makes sense
!========================================================================================!
!========================================================================================!
logical function crest_calc_is_subprocess(calc)
!*******************************************************************
!* logical function crest_calc_is_subprocess
!* Returns .true. if any of the calculations contained in the
!* calcdata object is performed by an external program, i.e. by a
!* system call (xtb, orca, turbomole, generic script, ...).
!*
!* In that case the CREST thread is blocked (waitpid) for the entire
!* duration of the calculation, while the number of cores used by the
!* subprocess is governed by OMP_NUM_THREADS (= cores_per_job).
!* Opening an additional nested OMP team on the CREST side therefore
!* adds no computing power at all: it only doubles the amount of
!* threads that is booked for one single job.
!*******************************************************************
  use crest_calculator
  implicit none
  type(calcdata),intent(in) :: calc
  integer :: i

  crest_calc_is_subprocess = .false.
  do i = 1,calc%ncalculations
    select case (calc%calcs(i)%id)
    case (jobtype%xtbsys,jobtype%generic,jobtype%turbomole, &
    &     jobtype%orca,jobtype%terachem)
      crest_calc_is_subprocess = .true.
      return
    end select
  end do
  return
end function crest_calc_is_subprocess

!========================================================================================!
!========================================================================================!
!> Helper to decide whether CREST is driving ORCA
!========================================================================================!
!========================================================================================!
logical function crest_calc_is_orca(calc)
!*******************************************************************
!* logical function crest_calc_is_orca
!* Returns .true. if any of the calculations contained in the
!* calcdata object is performed by the ORCA program package.
!*
!* Why it exists: ORCA is a full SCF.  One gradient costs tens of
!* seconds to minutes, against milliseconds for the internal GFN
!* methods, so any run that exists only to calibrate the MD time
!* step and to sample the wall time is pure overhead and must be
!* kept as short as it is still meaningful.  That default now
!* depends on the calculator, and this is the predicate it uses
!* (see the trial-MTD length in setuptest.f90).
!*
!* The test is "any", not "the first one": a composite calculator
!* (ONIOM) can carry ORCA in a non-leading slot, and in that case
!* it is still the ORCA subprocess that dominates the cost.
!*******************************************************************
  use crest_calculator
  implicit none
  type(calcdata),intent(in) :: calc
  integer :: i

  crest_calc_is_orca = .false.
!>--- calcs() is allocated lazily, so an unset-up calculator must not crash here
  if (.not.allocated(calc%calcs)) return
!>--- size(), not ncalculations: this is a pure predicate and must not depend on
!> the counter being in sync with the array it describes.
  do i = 1,size(calc%calcs)
    if (calc%calcs(i)%id .eq. jobtype%orca) then
      crest_calc_is_orca = .true.
      return
    end if
  end do
  return
end function crest_calc_is_orca

!========================================================================================!
!========================================================================================!
!> Routines for concurrent singlepoint evaluations
!========================================================================================!
!========================================================================================!
subroutine crest_sploop(env,nat,nall,at,xyz,eread)
!***************************************************************
!* subroutine crest_sploop
!* This subroutine performs concurrent singlepoint evaluations
!* for the given ensemble. Input eread is overwritten
!***************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_calculator
  use omp_lib
  use crest_data
  use strucrd
  use optimize_module
  use iomod,only:makedir,directory_exist,remove
  use crest_dynshare,only:crest_dynshare_init,crest_dynshare_release, &
  &                        crest_dynshare_reset
  implicit none
  type(systemdata),intent(inout) :: env
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(in)  :: at(nat)
  real(wp),intent(inout) :: eread(nall)
  integer,intent(in) :: nat,nall

  type(coord),allocatable :: mols(:)
  integer :: i,j,k,l,io,ich,ich2,c,z,job_id,zcopy
  logical :: pr,wr,ex
  type(calcdata),allocatable :: calculations(:)
  real(wp) :: energy
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  real(wp) :: percent,runtime

  type(timer) :: profiler
  integer :: T,Tn  !> threads and threads per core
  logical :: nested
  logical :: subproc
  logical,external :: crest_calc_is_subprocess

!>--- check if we have any calculation settings allocated
  if (env%calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

!>--- prepare calculation objects for parallelization (one per thread)
!>--- for external (system call) calculators the CREST thread is blocked while
!>--- the subprocess runs, hence no nested team must be opened on our side
  subproc = crest_calc_is_subprocess(env%calc)
  ! NOTE: new_ompautoset is an external procedure (no explicit interface at the
  ! call sites), so the "no nesting" switch is encoded in the modus string and
  ! NOT in an optional argument (that would be undefined behaviour).
  ! For subprocess calculators cores_per_job is handed to the external program
  ! via OMP_NUM_THREADS -> CREST must not open a nested team on top of it.
  ! 'fixed' (instead of 'auto') reserves CREST_CORES_PER_JOB cores (default 4)
  ! per job and runs T/cores_per_job jobs at a time; the OpenMP task queue
  ! keeps the remaining structures waiting.  With plain 'auto' an ensemble that
  ! is larger than the core count (the normal case) would hand out one core per
  ! job, i.e. hundreds of serial ORCA gradients.
  if (subproc) then
    call new_ompautoset(env,'fixed',nall,T,Tn)
  else
    call new_ompautoset(env,'auto_nested',nall,T,Tn)
  end if
  nested = env%omp_allow_nested.and..not.subproc
!>--- LOCAL EXTENSION: dynamic core sharing (see crest_dynshare), subprocess
!> only, exactly as in crest_search_multimd*().  The batch singlepoint loop is
!> where the idle-core problem is at its worst: T/cores_per_job jobs run at
!> once and the rest queue, so the tail of the ensemble (the last few
!> structures) is executed with a fraction of the machine while the rest sits
!> idle.  With this, every job that finishes hands its cores to the ones still
!> running.  Note that new_ompautoset('fixed') already spreads the whole budget
!> when nall < T/cores_per_job, so there is nothing to hand over in that case -
!> crest_dynshare_init() then stays off because T == 1.
  if (subproc) call crest_dynshare_init(nall,T,env%Threads,Tn)

!>--- prepare objects for parallelization
  allocate (calculations(T),source=env%calc)
  allocate (mols(T))
  do i = 1,T
    do j = 1,env%calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories and io preparation
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if(allocated(calculations(i)%calcs(j)%calcfile)) deallocate(calculations(i)%calcs(j)%calcfile)
      if(allocated(calculations(i)%calcs(j)%systemcall)) deallocate(calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
!>--- LOCAL EXTENSION: this loop ranks structures by ENERGY ONLY.  It used to
!> call engrad() and fill grads(), which was then deallocated without a single
!> read.  calc%nogr makes potential_core() dispatch to the gradient-free variant
!> (ORCA_sp for ORCA: no EnGrad keyword in the template, energy parsed from
!> ORCA.out).  singlepoint() raises the flag itself, too - setting it here as
!> well documents the intent where the ensemble is prepared, costs nothing, and
!> is done once per thread COPY outside the parallel region, so there is no race
!> and env%calc itself stays untouched.
    calculations(i)%nogr = .true.
    allocate (mols(i)%at(nat),mols(i)%xyz(3,nat))
  end do

!>--- printout directions and timer initialization
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log
  call profiler%init(1)
  call profiler%start(1)

!>--- first progress printout (initializes progress variables)
  call crest_oloop_pr_progress(env,nall,0)

!>--- shared variables
  c = 0  !> counter of successfull optimizations
  k = 0  !> counter of total optimization (fail+success)
  z = 0  !> counter to perform optimization in right order (1...nall)
  eread(:) = 0.0_wp
!>--- loop over ensemble
  !$omp parallel &
  !$omp shared(env,calculations,nat,nall,at,xyz,eread,c,k,z,pr,wr) &
  !$omp shared(ich,ich2,mols, nested,Tn)
  !$omp single
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private(i,j,job,energy,io,thread_id,zcopy)
    call initsignal()

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !>--- modify calculation spaces
    !$omp critical
    z = z+1
    zcopy = z
    mols(job)%nat = nat
    mols(job)%at(:) = at(:)
    mols(job)%xyz(:,:) = xyz(:,:,z)
    !$omp end critical

    !>-- energy call (NO gradient: crest_sploop only ever ranks by energy,
    !>   see the calc%nogr note in the setup loop above)
    call singlepoint(mols(job),calculations(job),energy,io)

    !$omp critical
    if (io == 0) then
      !>--- successful optimization (io==0)
      c = c+1
      eread(zcopy) = energy
    else
      eread(zcopy) = 0.0_wp
    end if
    k = k+1
    !>--- print progress
    call crest_oloop_pr_progress(env,nall,k)
    !$omp end critical
!>--- LOCAL EXTENSION: this job is done -> hand its cores to the still running
!> ones.  Must be OUTSIDE the critical above: crest_dynshare_release() takes its
!> own NAMED critical, and nesting two regions of the same (unnamed) lock would
!> deadlock.
    call crest_dynshare_release()
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel
!>--- LOCAL EXTENSION: restore the share this round started with
  call crest_dynshare_reset()

!>--- finalize progress printout
  call crest_oloop_pr_progress(env,nall,-1)

!>--- stop timer
  call profiler%stop(1)

!>--- prepare some summary printout
  percent = float(c)/float(nall)*100.0_wp
  write (atmp,'(f5.1,a)') percent,'% success)'
  write (stdout,'(">",1x,i0,a,i0,a,a)') c,' of ',nall,' structures successfully evaluated (', &
  &     trim(adjustl(atmp))
  write (atmp,'(">",1x,a,i0,a)') 'Total runtime for ',nall,' singlepoint calculations:'
!>--- LOCAL EXTENSION: for an external calculator (xtb/ORCA/turbomole/...) the
!> cpu-time measured by CREST is meaningless - see the note in timer.f90.
  call profiler%write_timing(stdout,1,trim(atmp),.true.,.not.subproc)
  runtime = profiler%get(1)
  write (atmp,'(f16.3,a)') runtime/real(nall,wp),' sec'
  write (stdout,'(a,a,a)') '> Corresponding to approximately ',trim(adjustl(atmp)), &
  &                       ' per processed structure'

  call profiler%clear()
  deallocate (calculations)
  if (allocated(mols)) deallocate (mols)
  return
end subroutine crest_sploop

!========================================================================================!
!========================================================================================!
!> Routines for concurrent geometry optimization
!========================================================================================!
!========================================================================================!
subroutine crest_oloop(env,nat,nall,at,xyz,eread,dump,customcalc)
!*******************************************************************************
!* subroutine crest_oloop
!* This subroutine performs concurrent geometry optimizations
!* for the given ensemble. Inputs xyz and eread are overwritten
!* env        - contains parallelization and other program settings
!* dump       - decides on whether to dump an ensemble file
!*              WARNING: the ensemble file will NOT be in the same order
!*              as the input xyz array. However, the overwritten xyz will be! 
!* customcalc - customized (optional) calculation level data
!*
!* IMPORTANT: xyz should be in Bohr(!) for this routine
!******************************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_calculator
  use omp_lib
  use crest_data
  use strucrd
  use optimize_module
  use iomod,only:makedir,directory_exist,remove
  use crest_restartlog,only:trackrestart,restart_write_dummy
  use crest_dynshare,only:crest_dynshare_init,crest_dynshare_release, &
  &                        crest_dynshare_reset
  use minregistry_module,only:minreg_reset,minreg_active,minreg_stats
  implicit none
  type(systemdata),target,intent(inout) :: env
  real(wp),intent(inout) :: xyz(3,nat,nall)
  integer,intent(in)  :: at(nat)
  real(wp),intent(inout) :: eread(nall)
  integer,intent(in) :: nat,nall
  logical,intent(in) :: dump
  type(calcdata),intent(in),target,optional :: customcalc

  type(coord),allocatable :: mols(:)
  type(coord),allocatable :: molsnew(:)
  integer :: i,j,k,l,io,ich,ich2,c,z,job_id,zcopy
  logical :: pr,wr,ex
  type(calcdata),allocatable :: calculations(:)
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grads(:,:,:)
  integer :: thread_id,vz,job
  character(len=80) :: atmp
  real(wp) :: percent,runtime
  type(calcdata),pointer :: mycalc
  type(timer) :: profiler
  integer :: T,Tn  !> threads and threads per core
  logical :: nested
  logical :: subproc
  logical,external :: crest_calc_is_subprocess

!>--- decide wether to skip this call
  if (trackrestart(env)) then
    call restart_write_dummy(ensemblefile)
    return
  end if

!>--- check which calc to use
  if(present(customcalc))then
    mycalc => customcalc
  else
    mycalc => env%calc
  endif

!>--- check if we have any calculation settings allocated
  if (mycalc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  end if

!>--- prepare calculation objects for parallelization (one per thread)
!>--- for external (system call) calculators the CREST thread is blocked while
!>--- the subprocess runs, hence no nested team must be opened on our side
  subproc = crest_calc_is_subprocess(mycalc)
  ! NOTE: new_ompautoset is an external procedure (no explicit interface at the
  ! call sites), so the "no nesting" switch is encoded in the modus string and
  ! NOT in an optional argument (that would be undefined behaviour).
  ! For subprocess calculators cores_per_job is handed to the external program
  ! via OMP_NUM_THREADS -> CREST must not open a nested team on top of it.
  ! 'fixed' (instead of 'auto') reserves CREST_CORES_PER_JOB cores (default 4)
  ! per job and runs T/cores_per_job jobs at a time; the OpenMP task queue
  ! keeps the remaining structures waiting.  With plain 'auto' an ensemble that
  ! is larger than the core count (the normal case) would hand out one core per
  ! job, i.e. hundreds of serial ORCA gradients.
  if (subproc) then
    call new_ompautoset(env,'fixed',nall,T,Tn)
  else
    call new_ompautoset(env,'auto_nested',nall,T,Tn)
  end if
  nested = env%omp_allow_nested.and..not.subproc
!>--- LOCAL EXTENSION: dynamic core sharing (see crest_dynshare), subprocess
!> only, exactly as in crest_search_multimd*().
!> THIS is the stage the idle-core problem hurts most: the ensemble is drained
!> in waves of T = Threads/CREST_CORES_PER_JOB jobs, so once the queue is
!> shorter than one wave the remaining structures each keep their original
!> CREST_CORES_PER_JOB cores while the rest of the machine idles - a 200
!> structure ensemble on 48 cores with 4 cores/job ends with 1-2 optimizations
!> running on 4 cores and 40 doing nothing.  Every job that finishes now hands
!> its share to the ones still running; the new share takes effect at the next
!> external call of those jobs (see the module header for how).
!> new_ompautoset('fixed') already spreads the whole budget when there are fewer
!> structures than slots, in which case T == 1 and crest_dynshare_init() stays
!> off on its own - there is nothing left to hand over.
  if (subproc) call crest_dynshare_init(nall,T,env%Threads,Tn)

!>--- prepare objects for parallelization
  allocate (calculations(T),source=mycalc)
  allocate (mols(T),molsnew(T))
  do i = 1,T
    do j = 1,mycalc%ncalculations
      calculations(i)%calcs(j) = mycalc%calcs(j)
      !>--- directories and io preparation
      ex = directory_exist(mycalc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(mycalc%calcs(j)%calcspace))
      end if
      if(calculations(i)%calcs(j)%id == jobtype%tblite)then
         calculations(i)%optnewinit=.true.
      endif
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = mycalc%calcs(j)%calcspace//trim(atmp)
      if(allocated(calculations(i)%calcs(j)%calcfile)) deallocate(calculations(i)%calcs(j)%calcfile)
      if(allocated(calculations(i)%calcs(j)%systemcall)) deallocate(calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
    allocate (mols(i)%at(nat),mols(i)%xyz(3,nat))
    allocate (molsnew(i)%at(nat),molsnew(i)%xyz(3,nat))
  end do

!>--- LOCAL EXTENSION (anchor pruning): a fresh registry for THIS optimisation
!> stage.  Resetting per stage is essential - otherwise a structure re-optimised
!> at tighter thresholds in the NEXT stage would recognise its own result from the
!> previous stage and stop at once, skipping the tighter optimisation entirely.
  call minreg_reset(nat)
  if (minreg_active()) then
    write (stdout,'(1x,a)') &
    & 'Anchor pruning ON: a structure descending into an already located minimum is stopped early.'
  end if

!>--- printout directions and timer initialization
  pr = .false. !> stdout printout
  wr = .false. !> write crestopt.log
  if (dump) then
    open (newunit=ich,file=ensemblefile)
    open (newunit=ich2,file=ensembleelog)
  end if
  call profiler%init(1)
  call profiler%start(1)

!>--- first progress printout (initializes progress variables)
  call crest_oloop_pr_progress(env,nall,0)

!>--- shared variables
  allocate (grads(3,nat,T),source=0.0_wp)
  c = 0  !> counter of successfull optimizations
  k = 0  !> counter of total optimization (fail+success)
  z = 0  !> counter to perform optimization in right order (1...nall)
  eread(:) = 0.0_wp
  grads(:,:,:) = 0.0_wp
!>--- loop over ensemble
  !$omp parallel &
  !$omp shared(env,calculations,nat,nall,at,xyz,eread,grads,c,k,z,pr,wr,dump) &
  !$omp shared(ich,ich2,mols,molsnew, nested,Tn)
  !$omp single
  do i = 1,nall

    call initsignal()
    vz = i
    !$omp task firstprivate( vz ) private(j,job,energy,io,atmp,gnorm,thread_id,zcopy)
    call initsignal()

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !>--- modify calculation spaces
    !$omp critical
    z = z+1
    zcopy = z
    mols(job)%nat = nat
    mols(job)%at(:) = at(:)
    mols(job)%xyz(:,:) = xyz(:,:,z)

    molsnew(job)%nat = nat
    molsnew(job)%at(:) = at(:)
    molsnew(job)%xyz(:,:) = xyz(:,:,z)
    !$omp end critical

    !>-- geometry optimization
    call optimize_geometry(mols(job),molsnew(job),calculations(job),energy,grads(:,:,job),pr,wr,io)

    !$omp critical
    if (io == 0) then
      !>--- successful optimization (io==0)
      c = c+1
      if (dump) then
        gnorm = norm2(grads(:,:,job))
        write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
        molsnew(job)%comment = trim(atmp)
        call molsnew(job)%append(ich)
        call calc_eprint(calculations(job),energy,calculations(job)%etmp,gnorm,ich2)
      end if
      eread(zcopy) = energy
      xyz(:,:,zcopy) = molsnew(job)%xyz(:,:)
    else if(io==calculations(job)%maxcycle .and. calculations(job)%anopt) then
      !>--- allow partial optimization?
      c = c+1
      eread(zcopy) = energy
      xyz(:,:,zcopy) = molsnew(job)%xyz(:,:)
!>--- LOCAL FIX (2026-10-28): DUMP the partially optimized structure too.
!> Without this the counter `c` above reports it as a success (the log then
!> says "N of N structures successfully optimized (100.0% success)") while
!> the structure was NEVER written to crest_ensemble.xyz - the only dump
!> call sat inside the io==0 branch.  If ALL structures of an ensemble end
!> up on this path, crest_ensemble.xyz is created EMPTY, it gets renamed to
!> crest_rotamers_<n>.xyz, and the next multilevel step dies with
!>
!>     **ERROR** empty ensemble file crest_rotamers_<n>.xyz
!>
!> Observed on a constrained TS search where all 6 structures hit the
!> 120-step limit with anopt enabled.  Same failure mode as the
!> "empty ensemble file crest_rotamers_2.xyz" abort seen on fw19.
      if (dump) then
        gnorm = norm2(grads(:,:,job))
        write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
        molsnew(job)%comment = trim(atmp)
        call molsnew(job)%append(ich)
        call calc_eprint(calculations(job),energy,calculations(job)%etmp,gnorm,ich2)
      end if
    else
      eread(zcopy) = 1.0_wp
    end if
    k = k+1
    !>--- print progress
    call crest_oloop_pr_progress(env,nall,k)
    !$omp end critical
!>--- LOCAL EXTENSION: this job is done -> hand its cores to the still running
!> ones.  Outside the critical above on purpose (crest_dynshare_release() takes
!> its own NAMED critical; nesting two regions of the same unnamed lock would
!> deadlock).  Called for failed optimizations as well - the job is finished
!> either way, which is exactly what the hand-over counts.
    call crest_dynshare_release()
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel
!>--- LOCAL EXTENSION: restore the share this round started with
  call crest_dynshare_reset()

!>--- finalize progress printout
  call crest_oloop_pr_progress(env,nall,-1)

!>--- stop timer
  call profiler%stop(1)

!>--- prepare some summary printout
  percent = float(c)/float(nall)*100.0_wp
  write (atmp,'(f5.1,a)') percent,'% success)'
  write (stdout,'(">",1x,i0,a,i0,a,a)') c,' of ',nall,' structures successfully optimized (', &
  &     trim(adjustl(atmp))
  write (atmp,'(">",1x,a,i0,a)') 'Total runtime for ',nall,' optimizations:'
!>--- LOCAL EXTENSION: for an external calculator (xtb/ORCA/turbomole/...) the
!> cpu-time measured by CREST is meaningless - see the note in timer.f90.
  call profiler%write_timing(stdout,1,trim(atmp),.true.,.not.subproc)
  runtime = profiler%get(1)
  write (atmp,'(f16.3,a)') runtime/real(nall,wp),' sec'
  write (stdout,'(a,a,a)') '> Corresponding to approximately ',trim(adjustl(atmp)), &
  &                       ' per processed structure'

!>--- close files (if they are open)
  if (dump) then
    close (ich)
    close (ich2)
  end if

!>--- LOCAL EXTENSION: how much did anchor pruning actually save in this stage?
  call minreg_stats(stdout)

  deallocate (grads)
  call profiler%clear()
  deallocate (calculations)
  if (allocated(mols)) deallocate (mols)
  if (allocated(molsnew)) deallocate (molsnew)
  return
end subroutine crest_oloop

!========================================================================================!
subroutine crest_oloop_pr_progress(env,total,current)
!*********************************************
!* subroutine crest_oloop_pr_progress
!* A subroutine to print and track progress of
!* concurrent geometry optimizations
!*********************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use iomod,only:to_str
!>--- LOCAL CHANGE (2026-10-04): tell the dynamic sharding that a bar is open
  use crest_dynshare,only:crest_dynshare_barmark
  implicit none
  type(systemdata),intent(inout) :: env
  integer,intent(in) :: total,current
  real(wp) :: percent
  character(len=5) :: atmp
  real(wp),save :: increment
  real(wp),save :: progressbarrier

  percent = float(current)/float(total)*100.0_wp
  if (current == 0) then !> as a wrapper to start the printout
    progressbarrier = 0.0_wp
    if (env%niceprint) then
      percent = 0.0_wp
      call printprogbar(percent)
      call crest_dynshare_barmark()
    end if
    increment = 10.0_wp
    if (total > 1000) increment = 7.5_wp
    if (total > 5000) increment = 5.0_wp
    if (total > 10000) increment = 2.5_wp
    if (total > 20000) increment = 1.0_wp

  else if (current <= total.and.current > 0) then !> the regular printout case
    if (env%niceprint) then
      call printprogbar(percent)
      call crest_dynshare_barmark()

    else if (.not.env%legacy) then
      if (percent >= progressbarrier) then
        write (atmp,'(f5.1)') percent
        write (stdout,'(1x,a)',advance='no') '|>'//trim(adjustl(atmp))//'%'
        call crest_dynshare_barmark()
        progressbarrier = progressbarrier+increment
        progressbarrier = min(progressbarrier,100.0_wp)
        flush (stdout)
      end if
    else
      write (stdout,'(1x,i0)',advance='no') current
      call crest_dynshare_barmark()
      flush (stdout)
    end if

  else !> as a wrapper to finalize the printout
    if (.not.env%niceprint) then
      write (stdout,'(/,1x,a)') 'done.'
    else
      write (stdout,*)
    end if
  end if

end subroutine crest_oloop_pr_progress

!========================================================================================!
!========================================================================================!
!> Routines for parallel MDs
!========================================================================================!
!========================================================================================!

subroutine crest_search_multimd(env,mol,mddats,nsim)
!*****************************************************
!* subroutine crest_search_multimd
!* this runs #nsim MDs on the same structure (mol)
!*****************************************************
  use crest_parameters,only:wp,int64,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  use crest_restartlog,only:trackrestart,restart_write_dummy
  use crest_dynshare,only:crest_dynshare_init,crest_dynshare_release,crest_dynshare_reset
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mol
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich
  logical :: pr,ex,nested
  integer :: T,Tn
  real(wp) :: percent
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES'

  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
  real(wp) :: etmp
  real(wp),allocatable :: grdtmp(:,:)
  type(timer) :: profiler
  logical :: subproc
  logical,external :: crest_calc_is_subprocess
  !>--- LOCAL EXTENSION: wall clock and wave count for the phase summary.
  !> The per-job "*MTD n completed successfully ..." lines only carry the wall
  !> time of ONE dynamics() call.  Everything else this stage does - above all
  !> the sequential per-slot engrad() of the worker containers, which for an
  !> external calculator is one full external call PER SLOT - is not covered by
  !> any timer and therefore never shows up.  tbmid splits the two.
  integer(int64) :: tb0,tb1,tbmid,brate
  real(wp) :: setwall,runwall,smins
  integer :: nwaves
!>--- LOCAL EXTENSION (2026-09-29): explicit interface for
!> parallel_md_block_printout().  It is a module-less subroutine, and since
!> 2026-09-29 it takes an OPTIONAL third dummy (ntot), so Fortran now requires
!> an explicit interface at EVERY call site:
!>
!>   Error: Explicit interface required for 'parallel_md_block_printout'
!>          at (1): optional argument
!>
!> An INTERNAL procedure is not an option here: the calls below sit inside
!> orphaned `!$omp do` / `!$omp task` constructs, i.e. whatever thread happens
!> to execute the body would be the host, and a host-associated internal
!> procedure would then be called from a different thread.  Keep this block and
!> the declaration at parallel.f90:1665-1667 in sync.
  interface
    subroutine parallel_md_block_printout(MD,vz,ntot)
      use crest_parameters,only:wp
      use dynamics_module,only:mddata
      implicit none
      type(mddata),intent(in) :: MD
      integer,intent(in) :: vz
      integer,intent(in),optional :: ntot
    end subroutine parallel_md_block_printout
  end interface
!===========================================================!
!>--- decide wether to skip this call
  if (trackrestart(env)) then
    call restart_write_dummy('crest_dynamics.xyz')
    return
  end if

!>--- check if we have any MD & calculation settings allocated
  if (.not.env%mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

!>--- LOCAL EXTENSION: start of the phase clock (see the summary further down).
!> It deliberately starts BEFORE the worker-container setup below, because that
!> setup is part of what this stage costs.
  call system_clock(count=tb0,count_rate=brate)

!>--- prepare calculation containers for parallelization (one per thread)
!>--- NOTE: for external (system call) calculators, e.g. the default xtb binary,
!>--- the CREST thread is blocked while the subprocess runs. cores_per_job is
!>--- then exclusively used by the subprocess (OMP_NUM_THREADS), and CREST must
!>--- NOT open an additional nested team, otherwise every MD would book
!>--- 2*cores_per_job threads instead of cores_per_job.
  subproc = crest_calc_is_subprocess(env%calc)
  ! NOTE: new_ompautoset is an external procedure (no explicit interface at the
  ! call sites), so the "no nesting" switch is encoded in the modus string and
  ! NOT in an optional argument (that would be undefined behaviour).
  ! For subprocess calculators cores_per_job is handed to the external program
  ! via OMP_NUM_THREADS -> CREST must not open a nested team on top of it.
  if (subproc) then
    call new_ompautoset(env,'auto',nsim,T,Tn)
  else
    call new_ompautoset(env,'auto_nested',nsim,T,Tn)
  end if
  nested = env%omp_allow_nested.and..not.subproc
!>--- LOCAL EXTENSION: dynamic core sharing.  subprocess calculators only: an
!>    in-process calculator spends its budget on a nested OMP team created inside
!>    dynamics(), whose size cannot be changed from out here (see crest_dynshare).
  if (subproc) call crest_dynshare_init(nsim,T,env%Threads,Tn)

  allocate (calculations(T),source=env%calc)
  allocate (moltmps(T),source=mol)
  allocate (grdtmp(3,mol%nat),source=0.0_wp)
  do i = 1,T
    moltmps(i)%nat = mol%nat
    moltmps(i)%at = mol%at
    moltmps(i)%xyz = mol%xyz
    do j = 1,env%calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories and io preparation
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if(allocated(calculations(i)%calcs(j)%calcfile)) deallocate(calculations(i)%calcs(j)%calcfile)
      if(allocated(calculations(i)%calcs(j)%systemcall)) deallocate(calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
!>--- LOCAL CHANGE (2026-09-26): the per-slot "warm-up" below only makes sense
!> for an IN-PROCESS calculator (tblite/gfn0/gfnff build their API state on the
!> first call).  For a SUBPROCESS calculator (ORCA/xtb/turbomole) there is no
!> in-memory state to warm up - ORCA_setup() rewrites ORCA.in on every call -
!> and the result (etmp/grdtmp) is NEVER READ: `etmp` appears only in its own
!> declaration and here, `grdtmp` only in its declaration, the allocate below
!> and here.  Meanwhile the MD itself evaluates this very geometry again in its
!> first step (dynamics_module.f90:410).  So this was `T` SERIAL external calls
!> (for ORCA: T whole SCFs) per MTD batch, thrown away.
!> Evidence that it is not needed: the sister routine crest_search_multimd2()
!> has the same setup loop (parallel.f90:1503) without this call, and its MDs
!> run fine.
    if (.not.subproc) then
      !>--- initialize the calculations
      call engrad(moltmps(i),calculations(i),etmp,grdtmp,io)
    end if
  end do

  !>--- LOCAL EXTENSION: the sequential worker setup is done -> mark it, so the
  !> summary below can separate "setup" from "actual MD sampling".
  call system_clock(count=tbmid)

  !>--- other settings
  pr = .false.
  call profiler%init(nsim)

!>--- LOCAL CHANGE (2026-10-04): the per-simulation header boxes are printed
!> HERE - sequentially, before the OpenMP region - instead of from inside it.
!> Inside the region the order is the order in which the threads happen to
!> reach their !$omp critical, which for a batch of 4 came out as 4, 2, 3, 1:
!> a log that cannot be read linearly.  Everything a box contains (MD length,
!> k, alpha, SHAKE mode, active potentials) was already set up by
!> crest_search_multimd_init2(), which ran before this point.
  do i = 1,nsim
    call parallel_md_block_printout(mddats(i),i,nsim)
  end do

  !>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mol,pr,percent,ich, nsim, moltmps, nested,Tn) &
  !!$omp single
  !$omp private(vz,i,job,thread_id,io,ex)
  !$omp do
  do i = 1,nsim

    call initsignal()
    vz = i

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    !!$omp task firstprivate( vz ) private( job,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !$omp critical
    moltmps(job)%nat = mol%nat
    moltmps(job)%at = mol%at
    moltmps(job)%xyz = mol%xyz
    !$omp end critical

    !>--- the acutal MD call with timing
    call profiler%start(vz)
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)
    call profiler%stop(vz)

    !>--- finish printout (thread safe)
    call parallel_md_finish_printout(mddats(vz),vz,io,profiler)
!>--- LOCAL EXTENSION: this job is done -> hand its cores to the jobs still running
    call crest_dynshare_release()
    !!$omp end task
  end do
  !!$omp taskwait
  !$omp end parallel
!>--- LOCAL EXTENSION: restore the share this round started with
  call crest_dynshare_reset()

!>--- LOCAL CHANGE (2026-10-05): a stop file that was honoured by this round is
!> deleted HERE, not by the job that saw it.  Every job of the round has
!> finished by now, so the file was visible to all of them for the whole round
!> (their polls are 300 s apart and are not synchronised), and it is gone
!> before the NEXT round starts - which is the point of removing it at all.
  if (md_stop_seen()) then
    call md_stop_cleanup()
    write (stdout,'(1x,a)') &
    & 'Structures so far are handed on to the next stage. Stop file is automatically removed.'
  end if

!>--- collect trajectories into one
  call collect(nsim,mddats)

!>--- LOCAL EXTENSION: wall time of the WHOLE MD/MTD phase.
!> The "*MTD n completed successfully ..." lines printed above only ever carry
!> the wall time of ONE simulation, i.e. of that single dynamics() call.  They
!> cannot add up to the wall time of the stage, because this routine additionally
!>  (a) runs the worker-container setup sequentially (one engrad() per slot -
!>      for an external calculator that is one full external call per slot), and
!>  (b) spreads nsim jobs over ceil(nsim/T) "waves" when nsim exceeds the number
!>      T of parallel slots, so the per-job timers overlap rather than add.
!> The summary below is the number to compare against the wall clock (or the
!> batch scheduler).  It reports the two halves separately, so it is immediately
!> visible whether a short per-job line is explained by waves or by setup cost.
  call system_clock(count=tb1)
  setwall = real(tbmid-tb0,wp)/real(max(brate,1_int64),wp)
  runwall = real(tb1-tbmid,wp)/real(max(brate,1_int64),wp)
  nwaves = max(1,ceiling(real(nsim)/real(max(1,T))))
  write (stdout,'(1x,a)') repeat('-',60)
  if (nwaves > 1) then
    write (stdout,'(1x,a,i0,a,i0,a,i0,a)') 'MD/MTD phase: ',nsim,' jobs on ',T, &
    &   ' slots (',nwaves,' waves)'
  else
    write (stdout,'(1x,a,i0,a,i0,a)') 'MD/MTD phase: ',nsim,' jobs on ',T,' slots'
  end if
  smins = aint(setwall/60._wp)
  write (stdout,'(1x,a,1x,"...",i9," min, ",f6.3," sec")') '  worker setup (sequential)', &
  &   nint(smins),setwall-60._wp*smins
  smins = aint(runwall/60._wp)
  write (stdout,'(1x,a,1x,"...",i9," min, ",f6.3," sec")') '  MD sampling (parallel)  ', &
  &   nint(smins),runwall-60._wp*smins
  smins = aint((setwall+runwall)/60._wp)
  write (stdout,'(1x,a,1x,"...",i9," min, ",f6.3," sec")') '  TOTAL of this MTD round ', &
  &   nint(smins),setwall+runwall-60._wp*smins
  write (stdout,'(1x,a)') repeat('-',60)

  call profiler%clear()
  deallocate (calculations)
  if (allocated(moltmps)) deallocate (moltmps)
  return
contains
  subroutine collect(n,mddats)
    implicit none
    integer :: n
    type(mddata) :: mddats(n)
    logical :: ex
    integer :: i,io,ich,ich2
    character(len=:),allocatable :: atmp
    character(len=256) :: btmp
    open (newunit=ich,file='crest_dynamics.xyz')
    do i = 1,n
      atmp = mddats(i)%trajectoryfile
      inquire (file=atmp,exist=ex)
      if (ex) then
        open (newunit=ich2,file=atmp)
        io = 0
        do while (io == 0)
          read (ich2,'(a)',iostat=io) btmp
          if (io == 0) then
            write (ich,'(a)') trim(btmp)
          end if
        end do
        close (ich2)
!>--- LOCAL CHANGE (2026-09-27): drop the per-worker trajectory as soon as its
!> content is inside crest_dynamics.xyz.  Keeping it doubled the on-disk
!> footprint of the whole sampling (worker file + merged copy, byte for byte),
!> and did so for the entire length of the run on a shared/NFS filesystem.
!> SAFE because nothing reads anything under MDFILES/ after this merge: the only other
!> reference to mddats(i)%trajectoryfile is the identical collect() in
!> crest_search_multimd2, and cleanup.f90 does nothing but rmrf the directory.
!> The next dynamics() call re-creates its file (dynamics_module.f90:350 opens
!> without status=, i.e. it writes from the start).
        call remove(atmp)
      end if
    end do
    close (ich)
    return
  end subroutine collect
end subroutine crest_search_multimd

!========================================================================================!
subroutine crest_search_multimd_init(env,mol,mddat,nsim)
!*******************************************************
!* subroutine crest_search_multimd_init
!* This routine will initialize a copy of env%mddat
!* and save it to the local mddat. If we are about to
!* run RMSD metadynamics, the required number of
!* simulations (#nsim) is returned
!*******************************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddat
  type(coord) :: mol
  integer,intent(inout) :: nsim
  integer :: io
  logical :: pr
!=======================================================!
  type(calcdata),target :: calc
  type(shakedata) :: shk

  real(wp) :: energy
  character(len=*),parameter :: mdir = 'MDFILES'
!======================================================!

  !>--- check if we have any MD & calculation settings allocated
  mddat = env%mddat
  if (.not.mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

  !>--- init SHAKE?
  if (mddat%shake) then
    if (allocated(env%ref%wbo)) then
      shk%wbo = env%ref%wbo
    else
      calc = env%calc
      calc%calcs(1)%rdwbo = .true.
!>--- LOCAL CHANGE (2026-09-26): only the WBOs are wanted here, and
!> rd_ORCA_wbo() parses ORCA.out (orca_sc.f90), which the child writes whether
!> or not the EnGrad keyword is present.  engrad() therefore made ORCA solve a
!> CP-SCF / Z-vector (plus the CPCM gradient term) whose gradient landed in
!> `grad` and was thrown away - `grad` appeared only in its declaration, the
!> allocate, that call and the deallocate.  singlepoint() yields the same WBOs
!> and the same energy without paying for the gradient.
      call singlepoint(mol,calc,energy,io)
      calc%calcs(1)%rdwbo = .false.

      shk%shake_mode = env%mddat%shk%shake_mode
      call move_alloc(calc%calcs(1)%wbo,shk%wbo)
    end if

    if (calc%nfreeze > 0) then
      shk%freezeptr => calc%freezelist
    else
      nullify (shk%freezeptr)
    end if

    shk%shake_mode = env%shake
!>--- LOCAL EXTENSION: mddat%shk = shk REPLACES the shk that came out of
!> "mddat = env%mddat", so the rigid bonds filled by choose_settings() would be
!> lost here.  Carry them over into the local shakedata explicitly.
    if (env%mddat%shk%nusr > 0) then
      shk%nusr = env%mddat%shk%nusr
      if (allocated(env%mddat%shk%conslistu)) shk%conslistu = env%mddat%shk%conslistu
      if (allocated(env%mddat%shk%distconsu)) shk%distconsu = env%mddat%shk%distconsu
    end if
    mddat%shk = shk
    call init_shake(mol%nat,mol%at,mol%xyz,mddat%shk,pr)
    mddat%nshake = mddat%shk%ncons
  end if
  !>--- complete real-time settings to steps
  call mdautoset(mddat,io)

  !>--- (optional)  MTD initialization
  if (nsim < 0) then
    mddat%simtype = type_mtd  !>-- set runtype to MTD

    call defaultGF(env)
!>--- LOCAL CHANGE (2026-10-04): the "$metadyn <k> <alpha>" dump that used to be
!> printed here is gone.  It listed exactly the same k/alpha pairs that the
!> per-simulation header boxes carry below, one box per MTD, so it was pure
!> duplication - and for a large batch a long wall of numbers nobody reads.

    !>--- how many simulations
    nsim = env%nmetadyn
  end if

  return
end subroutine crest_search_multimd_init

!========================================================================================!
subroutine crest_search_multimd_init2(env,mddats,nsim)
  use crest_parameters,only:wp,stdout,sep,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use iomod,only:makedir,directory_exist,remove,crest_hot_dir
  use omp_lib
  implicit none
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  integer :: i,io,j,k
  logical :: ex
!========================================================!
  type(mtdpot),allocatable :: mtds(:)

  character(len=80) :: atmp
!>--- LOCAL EXTENSION (-iru): static bias built from the conformers found so far
  logical :: iruon
  integer :: irunat,irunall,iruio
  real(wp) :: irukscal
  integer,allocatable :: iruat(:)
  real(wp),allocatable :: iruxyz(:,:,:)
  real(wp),allocatable :: irueread(:)
  character(len=128) :: envstr
  logical,save :: irufirst = .true.
!>--- LOCAL EXTENSION (--loadconformer): which file the static bias is read from.
!>    -iru            -> conformerfile      ('crest_conformers.xyz', written by
!>                       the multilevel optimization, i.e. only exists from the
!>                       2nd MTD iteration on)
!>    --loadconformer -> loadconformerfile  ('readconformers.xyz', supplied by
!>                       the user in the CWD, i.e. available from the 1st MTD on)
!>    If both are given, --loadconformer wins: the user file is the more
!>    specific request and it is the one that is guaranteed to be complete.
  character(len=:),allocatable :: irufile
  logical :: iruuser
!> LOCAL EXTENSION: the trajectories are appended to every sdump steps, i.e.
!> this directory is the hottest thing CREST writes outside the calculation
!> directories.  A symlink would not survive the rmrf()+makedir() below, so the
!> path itself is relocated here (mdir is only used in this routine).
  character(len=:),allocatable :: mdir
  character(len=:),allocatable :: mdirtmp

  !>--- parallel MD setup
  call crest_hot_dir(mdirtmp,'MDFILES')
  mdir = mdirtmp
  ex = directory_exist(mdir)
  if (ex) then
    call rmrf(mdir)
  end if
  !> the parent exists: crest_scratch_dir() created the root with mkdir -p
  io = makedir(mdir)
  do i = 1,nsim
    mddats(i)%md_index = i
    write (atmp,'(a,i0,a)') 'crest_',i,'.xyz'
    mddats(i)%trajectoryfile = mdir//sep//trim(atmp)
    write (atmp,'(a,i0,a)') 'crest_',i,'.mdrestart'
    mddats(i)%restartfile = mdir//sep//trim(atmp)
  end do

!========================================================================================!
!> LOCAL EXTENSION: wire up "-iru" (re-use the conformers found so far as an
!> additional STATIC RMSD bias) and "--loadconformer" (the same mechanism, but
!> the ensemble comes from the user-supplied readconformers.xyz).
!>
!> Upstream declares env%iru (classes.f90:560), defaults it (confparse.f90:247)
!> and parses the flag (confparse.f90:917-918) - but NEVER reads it, so the
!> switch was dead code.  Here it adds a SECOND potential, i.e. npot becomes 2:
!>    mtd(1) = dynamic RMSD MTD   -> unchanged behaviour, still deposits hills
!>    mtd(2) = static RMSD bias   -> built from the conformers found so far
!> Both contribute to the energy/gradient (md_calc_mtd loops over npot).
!>
!> The reference structures are read ONCE here and handed to every job as a
!> pre-filled cvxyz (in Bohr).  mtd_ini() then takes its
!> "if (allocated(pot%cvxyz))" branch, so no job ever re-opens the file: that
!> avoids nsim concurrent readers and any race with cregen() rewriting
!> crest_conformers.xyz later in the same iteration.
!>
!> conformerfile ('crest_conformers.xyz') is only produced by the multilevel
!> optimization, so -iru automatically kicks in from the 2nd MTD iteration on.
!>
!> LOCAL EXTENSION, "--loadconformer": the SAME mechanism, but the ensemble is
!> read from loadconformerfile ('readconformers.xyz'), which the user places in
!> the CWD before the run.  Because it already exists when the first MTD batch
!> is set up, the static bias is active from the very FIRST MTD - that is the
!> whole point of the option (push the sampling away from structures you have
!> already characterised, e.g. from a previous crest run on the same molecule).
!> Everything below is shared; only the file name and two printouts differ.
!>
!> KNOB: CREST_IRU_KSCALE scales k of the static bias.  Default 1.0 = "normal",
!> i.e. exactly the kpush/alpha this MTD job already uses for its dynamic bias
!> (env%metadfac/metadexp, the "$metadyn k alpha" list).  NOTE that the static
!> bias sums over ALL conformers, so its magnitude grows with the ensemble size
!> - dial this down (e.g. 0.1) if the structures get blown apart.
!> For --loadconformer, CREST_LOADCONF_KSCALE is probed FIRST and overrides
!> CREST_IRU_KSCALE, so a large user ensemble can be scaled separately from the
!> (smaller) auto-generated one.
!========================================================================================!
  iruon = .false.
  iruuser = .false.
  irunat = 0
  irunall = 0
  irukscal = 1.0_wp
  irufile = conformerfile
!>--- LOCAL EXTENSION (2026-10-26): prefer the PREVIOUS iteration's FULL rotamer
!> ensemble over the few conformers that survived the ewin/rthr filtering.
!> search_conformers.f90 writes 'crest_irubias.xyz' at the end of every
!> MTD iteration; it is absent on the first one, so round 1 falls back to
!> crest_conformers.xyz exactly as before.  CREST_IRU_ROTAMERS=0 restores the
!> old behaviour.
  if (env%iru) then
    block
      character(len=32) :: rstr
      integer :: rstat
      logical :: rex
      rstr = ''
      call get_environment_variable('CREST_IRU_ROTAMERS',value=rstr,status=rstat)
      if (.not.(rstat == 0 .and. trim(rstr) == '0')) then
        inquire (file='crest_irubias.xyz',exist=rex)
        if (rex) irufile = 'crest_irubias.xyz'
      end if
    end block
  end if
  if (env%loadconformer) then
    irufile = loadconformerfile
    iruuser = .true.
  end if
!>--- env%metadfac is only allocated when the MTD setup ran (defaultGF).
!>    rotamermds() also calls this routine, but with plain MDs and no GFs.
  if ((env%iru.or.env%loadconformer).and.allocated(env%metadfac)) then
    envstr = ''
    if (iruuser) call get_environment_variable('CREST_LOADCONF_KSCALE',envstr)
!>--- fall back to the -iru knob if the --loadconformer-specific one is unset
    if (len_trim(envstr) == 0) call get_environment_variable('CREST_IRU_KSCALE',envstr)
    if (len_trim(envstr) > 0) then
      read (envstr,*,iostat=iruio) irukscal
      if (iruio /= 0) irukscal = 1.0_wp
    end if
    inquire (file=irufile,exist=ex)
    if (ex) then
      call rdensembleparam(irufile,irunat,irunall)
      if (irunall > 0.and.irunat == env%ref%nat) then
        allocate (iruxyz(3,irunat,irunall),source=0.0_wp)
        allocate (iruat(irunat),source=0)
        allocate (irueread(irunall),source=0.0_wp)
        call rdensemble(irufile,irunat,irunall,iruat,iruxyz,irueread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: cvxyz must be in Bohrs (mol%xyz is, too)
        iruxyz = iruxyz/bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
        iruon = .true.
        if (iruuser) then
          write (stdout,'(1x,a)') 'using the user-supplied conformers as an additional static MTD bias (--loadconformer)'
        else
          write (stdout,'(1x,a)') 're-using the conformers found so far as an additional static MTD bias (-iru)'
        end if
        write (stdout,'(3x,a,a)')    'bias ensemble   : ',irufile
        write (stdout,'(3x,a,i0)')   'bias structures : ',irunall
        write (stdout,'(3x,a,f10.5,a)') 'k (static, max) : ',maxval(env%metadfac)*irukscal,' Eh'
        write (stdout,'(3x,a,f8.4)') 'k-scaling       : ',irukscal
        write (stdout,*)
      else
        if (iruuser) then
          write (stdout,'(1x,a)') '**WARNING** --loadconformer: ensemble does not match the system, static bias disabled'
        else
          write (stdout,'(1x,a)') '**WARNING** -iru: conformer ensemble does not match the system, static bias disabled'
        end if
        write (stdout,'(3x,a,a)') 'ensemble file: ',irufile
        write (stdout,'(3x,a,i0,a,i0)') 'ensemble: ',irunall,' structures with ',irunat
        write (stdout,'(3x,a,i0)') 'expected nat: ',env%ref%nat
      end if
    else
      if (iruuser) then
        write (stdout,'(1x,a)') '**WARNING** --loadconformer: "'//irufile//'" not found, static bias disabled'
      else
        if (irufirst) then
          irufirst = .false.
          write (stdout,'(1x,a)') '-iru: no conformer ensemble yet, the static bias starts with the 2nd MTD iteration'
        end if
      end if
    end if
  end if

  allocate (mtds(nsim))
  do i = 1,nsim
    if (mddats(i)%simtype == type_mtd) then
!>--- a TOML [dynamics.meta] block may have left a potential here; drop it so
!>    the allocate() below cannot fail (mirrors search_entropy.f90:471-475)
      if (allocated(mddats(i)%mtd)) deallocate (mddats(i)%mtd)
      if (allocated(mddats(i)%cvtype)) deallocate (mddats(i)%cvtype)
      mtds(i)%kpush = env%metadfac(i)
      mtds(i)%alpha = env%metadexp(i)
      mtds(i)%cvdump_fs = float(env%mddump)
      mtds(i)%mtdtype = cv_rmsd

      if (iruon) then
        mddats(i)%npot = 2
        allocate (mddats(i)%mtd(2),source=mtds(i))
        allocate (mddats(i)%cvtype(2),source=cv_rmsd)
        mddats(i)%cvtype(2) = cv_rmsd_static
        mddats(i)%mtd(2)%mtdtype = cv_rmsd_static
        mddats(i)%mtd(2)%kpush = mtds(i)%kpush*irukscal
!>--- a static potential never dumps new CVs
        mddats(i)%mtd(2)%cvdump_fs = 0.0_wp
        mddats(i)%mtd(2)%cvdumpstep = 0
        allocate (mddats(i)%mtd(2)%cvxyz(3,irunat,irunall),source=iruxyz)
      else
        mddats(i)%npot = 1
        allocate (mddats(i)%mtd(1),source=mtds(i))
        allocate (mddats(i)%cvtype(1),source=cv_rmsd)
      end if
      !> if necessary exclude atoms from RMSD bias (all potentials: the static
      !> -iru bias must use the same atom mask as the dynamic one)
      if (sum(env%includeRMSD) /= env%ref%nat) then
        do k = 1,mddats(i)%npot
          if (.not.allocated(mddats(i)%mtd(k)%atinclude)) &
          & allocate (mddats(i)%mtd(k)%atinclude(env%ref%nat),source=.true.)
          do j = 1,env%ref%nat
            if (env%includeRMSD(j) .ne. 1) mddats(i)%mtd(k)%atinclude(j) = .false.
          end do
        end do
      end if
    end if

    !>--- LOCAL EXTENSION: coarse progress marks for the production MD/MTD.
    !> These runs use pr=.false. - upstream keeps the log completely silent for
    !> the whole simulation, which for a slow external calculator looks like a
    !> hang.  progress_steps<=0 lets dynamics() derive length_steps/4, i.e. one
    !> line at 25/50/75 %.  (There is deliberately NO 100 % mark: it duplicates
    !> the completion line that follows immediately.)
    !>
    !> The job index is part of the tag on purpose: nsim jobs share the core
    !> budget and their lines interleave, so "*MTD  3  25% ..." could otherwise
    !> not be attributed to a simulation.
    mddats(i)%progress = .true.
    mddats(i)%progress_steps = 0
    !> The prefix deliberately mirrors parallel_md_finish_printout's
    !> write(...,'(a,1x,i3,a)') '*MTD',vz,' completed successfully' so that the
    !> marks and the completion line of one job line up on the same column.
    if (mddats(i)%simtype == type_mtd) then
      if (allocated(mddats(i)%cvtype)) then
        if (mddats(i)%cvtype(1) == cv_rmsd_static) then
          write (atmp,'(a,1x,i3)') '*sMTD',i
        else
          write (atmp,'(a,1x,i3)') '*MTD',i
        end if
      else
        write (atmp,'(a,1x,i3)') '*MTD',i
      end if
    else
      write (atmp,'(a,1x,i3)') '*MD',i
    end if
    mddats(i)%progress_tag = atmp
  end do
  if (allocated(mtds)) deallocate (mtds)

  return
end subroutine crest_search_multimd_init2

!========================================================================================!
subroutine crest_search_multimd2(env,mols,mddats,nsim)
!*******************************************************************
!* subroutine crest_search_multimd2
!* this runs #nsim MDs on #nsim selected different structures (mols)
!*******************************************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod,only:makedir,directory_exist,remove
  use omp_lib
  use crest_restartlog,only:trackrestart,restart_write_dummy
  use crest_dynshare,only:crest_dynshare_init,crest_dynshare_release,crest_dynshare_reset
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(mddata) :: mddats(nsim)
  integer :: nsim
  type(coord) :: mols(nsim)
  type(coord),allocatable :: moltmps(:)
  integer :: i,j,io,ich
  logical :: pr,ex,nested
  integer :: T,Tn
  real(wp) :: percent
  character(len=80) :: atmp
  character(len=*),parameter :: mdir = 'MDFILES'

  type(calcdata),allocatable :: calculations(:)
  integer :: vz,job,thread_id
  type(timer) :: profiler
  logical :: subproc
  logical,external :: crest_calc_is_subprocess
!>--- LOCAL EXTENSION (2026-09-29): explicit interface for
!> parallel_md_block_printout() - see the identical block in
!> crest_search_multimd above for the full explanation (optional third dummy
!> ntot; internal procedure impossible inside an orphaned OpenMP construct).
  interface
    subroutine parallel_md_block_printout(MD,vz,ntot)
      use crest_parameters,only:wp
      use dynamics_module,only:mddata
      implicit none
      type(mddata),intent(in) :: MD
      integer,intent(in) :: vz
      integer,intent(in),optional :: ntot
    end subroutine parallel_md_block_printout
  end interface
!===========================================================!
!>--- decide wether to skip this call
  if (trackrestart(env)) then
    call restart_write_dummy('crest_dynamics.xyz')
    return
  end if

!>--- check if we have any MD & calculation settings allocated
  if (.not.env%mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    return
  else if (env%calc%ncalculations < 1) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    return
  end if

!>--- prepare calculation objects for parallelization (one per thread)
  subproc = crest_calc_is_subprocess(env%calc)
  ! NOTE: new_ompautoset is an external procedure (no explicit interface at the
  ! call sites), so the "no nesting" switch is encoded in the modus string and
  ! NOT in an optional argument (that would be undefined behaviour).
  ! For subprocess calculators cores_per_job is handed to the external program
  ! via OMP_NUM_THREADS -> CREST must not open a nested team on top of it.
  if (subproc) then
    call new_ompautoset(env,'auto',nsim,T,Tn)
  else
    call new_ompautoset(env,'auto_nested',nsim,T,Tn)
  end if
  nested = env%omp_allow_nested.and..not.subproc
!>--- LOCAL EXTENSION: dynamic core sharing (see crest_dynshare), subprocess only
  if (subproc) call crest_dynshare_init(nsim,T,env%Threads,Tn)

  allocate (calculations(T),source=env%calc)
  allocate (moltmps(T),source=mols(1))
  do i = 1,T
    do j = 1,env%calc%ncalculations
      calculations(i)%calcs(j) = env%calc%calcs(j)
      !>--- directories and io preparation
      ex = directory_exist(env%calc%calcs(j)%calcspace)
      if (.not.ex) then
        io = makedir(trim(env%calc%calcs(j)%calcspace))
      end if
      write (atmp,'(a,"_",i0)') sep,i
      calculations(i)%calcs(j)%calcspace = env%calc%calcs(j)%calcspace//trim(atmp)
      if(allocated(calculations(i)%calcs(j)%calcfile)) deallocate(calculations(i)%calcs(j)%calcfile)
      if(allocated(calculations(i)%calcs(j)%systemcall)) deallocate(calculations(i)%calcs(j)%systemcall)
      call calculations(i)%calcs(j)%printid(i,j)
    end do
    calculations(i)%pr_energies = .false.
  end do

!>--- other settings
  pr = .false.
  call profiler%init(nsim)

!>--- LOCAL CHANGE (2026-10-04): startup boxes in ascending order, outside the
!> OpenMP region - same reasoning as in crest_search_multimd() above.
  do i = 1,nsim
    call parallel_md_block_printout(mddats(i),i,nsim)
  end do

!>--- run the MDs
  !$omp parallel &
  !$omp shared(env,calculations,mddats,mols,pr,percent,ich, moltmps,profiler, nested,Tn)
  !$omp single
  do i = 1,nsim

    call initsignal()
    vz = i

    !>--- OpenMP nested region threads
    if (nested) call ompmklset(Tn)

    !$omp task firstprivate( vz ) private( job,thread_id,io,ex )
    call initsignal()

    thread_id = OMP_GET_THREAD_NUM()
    job = thread_id+1
    !$omp critical
    moltmps(job)%nat = mols(vz)%nat
    moltmps(job)%at = mols(vz)%at
    moltmps(job)%xyz = mols(vz)%xyz
    !$omp end critical

    !>--- the acutal MD call with timing
    call profiler%start(vz)
    call dynamics(moltmps(job),mddats(vz),calculations(job),pr,io)
    call profiler%stop(vz)

    !>--- finish printout (thread safe)
    call parallel_md_finish_printout(mddats(vz),vz,io,profiler)
!>--- LOCAL EXTENSION: this job is done -> hand its cores to the jobs still running
    call crest_dynshare_release()
    !$omp end task
  end do
  !$omp taskwait
  !$omp end single
  !$omp end parallel
!>--- LOCAL EXTENSION: restore the share this round started with
  call crest_dynshare_reset()

!>--- LOCAL CHANGE (2026-10-05): same as in crest_search_multimd() above - the
!> honoured stop file is removed once the whole round is over, so it cannot cut
!> the next one short as well.
  if (md_stop_seen()) then
    call md_stop_cleanup()
    write (stdout,'(1x,a)') &
    & 'Structures so far are handed on to the next stage. Stop file is automatically removed.'
  end if

!>--- collect trajectories into one
  call collect(nsim,mddats)

  call profiler%clear()
  deallocate (calculations)
  if (allocated(moltmps)) deallocate (moltmps)
  return
contains
  subroutine collect(n,mddats)
    implicit none
    integer :: n
    type(mddata) :: mddats(n)
    logical :: ex
    integer :: i,io,ich,ich2
    character(len=:),allocatable :: atmp
    character(len=256) :: btmp
    open (newunit=ich,file='crest_dynamics.xyz')
    do i = 1,n
      atmp = mddats(i)%trajectoryfile
      inquire (file=atmp,exist=ex)
      if (ex) then
        open (newunit=ich2,file=atmp)
        io = 0
        do while (io == 0)
          read (ich2,'(a)',iostat=io) btmp
          if (io == 0) then
            write (ich,'(a)') trim(btmp)
          end if
        end do
        close (ich2)
!>--- LOCAL CHANGE (2026-09-27): drop the per-worker trajectory as soon as its
!> content is inside crest_dynamics.xyz.  Keeping it doubled the on-disk
!> footprint of the whole sampling (worker file + merged copy, byte for byte),
!> and did so for the entire length of the run on a shared/NFS filesystem.
!> SAFE because nothing reads anything under MDFILES/ after this merge: the only other
!> reference to mddats(i)%trajectoryfile is the identical collect() in
!> crest_search_multimd2, and cleanup.f90 does nothing but rmrf the directory.
!> The next dynamics() call re-creates its file (dynamics_module.f90:350 opens
!> without status=, i.e. it writes from the start).
        call remove(atmp)
      end if
    end do
    close (ich)
    return
  end subroutine collect
end subroutine crest_search_multimd2

!========================================================================================!
subroutine parallel_md_block_printout(MD,vz,ntot)
!***********************************************
!* subroutine parallel_md_block_printout
!* This will print information about the MD/MTD
!* simulation. The execution is omp threadsave
!*
!* LOCAL EXTENSION: ntot is the number of simulations of this batch.  The
!* closing rule (45 colons, same width as the headline) is only printed for
!* vz == ntot, i.e. once, underneath the LAST box of the batch - the
!* intermediate boxes are left open, which keeps N consecutive boxes from
!* turning into 2N rules.  If ntot is absent (or <= 0) every box is closed,
!* which is the safe fallback for a caller that does not know the total.
!***********************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod,only:to_str
  implicit none
  type(mddata),intent(in) :: MD
  integer,intent(in) :: vz
  integer,intent(in),optional :: ntot   !> size of the batch, see above
  character(len=40) :: atmp
  integer :: il
  logical :: closebox
  !$omp critical

  if (MD%simtype == type_md) then
    write (atmp,'(a,1x,i3)') 'starting MD',vz
  else if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd_static) then
      write (atmp,'(a,1x,i3)') 'starting static MTD',vz
    else
      write (atmp,'(a,1x,i4)') 'starting MTD',vz
    end if
  end if
  il = (44-len_trim(atmp))/2
  write (stdout,'(2x,a,1x,a,1x,a)') repeat(':',il),trim(atmp),repeat(':',il)

  write (stdout,'(2x,"|   MD simulation time   :",f8.1," ps       |")') MD%length_ps
  write (stdout,'(2x,"|   target T             :",f8.1," K        |")') MD%tsoll
  write (stdout,'(2x,"|   timestep dt          :",f8.1," fs       |")') MD%tstep
  write (stdout,'(2x,"|   dump interval(trj)   :",f8.1," fs       |")') MD%dumpstep
  if (MD%shake.and.MD%shk%shake_mode > 0) then
    if (MD%shk%shake_mode == 2) then
      write (stdout,'(2x,"|   SHAKE algorithm      :",a5," (all bonds) |")') to_str(MD%shake)
    else
      write (stdout,'(2x,"|   SHAKE algorithm      :",a5," (H only) |")') to_str(MD%shake)
    end if
  end if
  if (allocated(MD%active_potentials)) then
    write (stdout,'(2x,"|   active potentials    :",i4," potential    |")') size(MD%active_potentials,1)
  end if
  if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd) then
      write (stdout,'(2x,"|   dump interval(Vbias) :",f8.2," ps       |")') &
          & MD%mtd(1)%cvdump_fs/1000.0_wp
    end if
    write (stdout,'(2x,"|   Vbias prefactor (k)  :",f8.4," Eh       |")') &
      &  MD%mtd(1)%kpush
    if (MD%cvtype(1) == cv_rmsd.or.MD%cvtype(1) == cv_rmsd_static) then
      write (stdout,'(2x,"|   Vbias exponent (α)   :",f8.4," bohr⁻²   |")') MD%mtd(1)%alpha
    else
      write (stdout,'(2x,"|   Vbias exponent (α)   :",f8.4,"          |")') MD%mtd(1)%alpha
    end if
!>--- LOCAL EXTENSION (-iru): report the additional static bias, if present
    if (allocated(MD%cvtype)) then
      if (size(MD%cvtype) >= 2) then
        if (MD%cvtype(2) == cv_rmsd_static.and.allocated(MD%mtd(2)%cvxyz)) then
          write (stdout,'(2x,"|   static bias (-iru)   :",i8," structs  |")') &
          &  size(MD%mtd(2)%cvxyz,3)
          write (stdout,'(2x,"|   static bias k        :",f8.4," Eh       |")') &
          &  MD%mtd(2)%kpush
        end if
      end if
    end if
  end if

!>--- LOCAL EXTENSION: closing rule.  Same width as the headline above
!> (il+1+len_trim(atmp)+1+il = 45 columns for the usual "starting MTD n").
  closebox = .true.
  if (present(ntot)) closebox = (ntot .le. 0).or.(vz .ge. ntot)
  if (closebox) write (stdout,'(2x,a)') repeat(':',45)

  !$omp end critical

end subroutine parallel_md_block_printout

subroutine parallel_md_finish_printout(MD,vz,io,profiler)
!*******************************************
!* subroutine parallel_md_finish_printout
!* This will print information termination
!* info about the MD/MTD simulation
!*******************************************
  use crest_parameters,only:wp,stdout,sep
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  implicit none
  type(mddata),intent(in) :: MD
  integer,intent(in) :: vz,io
  type(timer),intent(inout) :: profiler
  character(len=40) :: atmp
  character(len=80) :: btmp

  !$omp critical

  if (MD%simtype == type_mtd) then
    if (MD%cvtype(1) == cv_rmsd_static) then
      write (atmp,'(a)') '*sMTD'
    else
      write (atmp,'(a)') '*MTD'
    end if
  else
    write (atmp,'(a)') '*MD'
  end if
  select case (io)
  case (0)
    write (btmp,'(a,1x,i3,a)') trim(atmp),vz,' completed successfully'
!>--- LOCAL EXTENSION: io = 3 is the deliberate "stop" file.  It is still an
!> EARLY end, but saying so is what makes a log that ends after five minutes
!> self-explanatory instead of looking like a crash.
  case (3)
    write (btmp,'(a,1x,i3,a)') trim(atmp),vz,' stopped ("stop" file)'
  case default
    write (btmp,'(a,1x,i3,a)') trim(atmp),vz,' terminated EARLY'
  end select
  call profiler%write_timing(stdout,vz,trim(btmp))

  !$omp end critical

end subroutine parallel_md_finish_printout
!========================================================================================!
!========================================================================================!
