!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c OMP and MKL parallelization settings
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine ompmklset(threads)
  use omp_lib
  implicit none
  integer,intent(in) :: threads

  call OMP_Set_Num_Threads(threads)
#ifdef WITH_MKL
  call MKL_Set_Num_Threads(threads)
  call mkl_set_dynamic(0)
#endif
! call openblasset(threads)
end subroutine ompmklset

subroutine openblasset(threads)
  implicit none
  integer,intent(in) :: threads
#ifdef WITH_OPENBLAS
  call openblas_set_num_threads(threads)
#endif
  return
end subroutine openblasset

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c OMP and MKL parallelization settings (short routine)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine ompenvset(omp)
  use iomod
  implicit none
  integer,intent(in) :: omp
  integer :: io

  io = setenv('OMP_NUM_THREADS',omp)
  io = setenv('MKL_NUM_THREADS',omp)
  io = setenv('OPENBLAS_NUM_THREADS',omp)
  !>--- NOTE: this routine prepares the environment for an EXTERNAL program
  !> (xtb, orca, ...).  Those programs parallelise in more than one way: xtb
  !> runs OpenMP loops *and* calls a threaded BLAS from them, i.e. it literally
  !> spawns "threads within threads".  With an outer budget of N and an inner
  !> budget of N the child can end up with N*N runnable threads per job, which
  !> over-subscribes the machine and shows up as low per-process CPU together
  !> with a large amount of system (red) time.
  !> Setting OMP_MAX_ACTIVE_LEVELS=1 deactivates the nested OpenMP level, so
  !> the outer N threads remain the only OpenMP layer.  This is exactly what
  !> the xtb manual recommends for parallel runs:
  !>   export OMP_NUM_THREADS=...   OMP_MAX_ACTIVE_LEVELS=1
  !> It only affects child processes - CREST configures its own runtime with
  !> omp_set_max_active_levels()/OMP_Set_Num_Threads() in new_ompautoset().
  !> Remove this one line to revert to the plain behaviour.
  io = setenv('OMP_MAX_ACTIVE_LEVELS',1)

end subroutine ompenvset

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c OMP and MKL autoset switchcase routine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine new_ompautoset(env,modus,maxjobs,parallel_jobs,cores_per_job)
  use omp_lib
  use crest_data
  use crest_parameters,only:wp
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in)    :: modus
  integer,intent(in)  :: maxjobs
  integer,intent(out) :: parallel_jobs
  integer,intent(out) :: cores_per_job
  !> NOTE: this is an EXTERNAL procedure (no module, no explicit interface at
  !> any of the ~50 call sites).  An OPTIONAL dummy argument is therefore
  !> undefined behaviour here: callers that pass only five arguments leave the
  !> sixth slot filled with garbage, present() can then report .true. for an
  !> argument that was never passed, and the dereference segfaults.
  !> Any behavioural switch must be encoded in the "modus" string instead:
  !>   'auto'        -> distribute the threads, but do NOT open a nested team
  !>   'auto_nested' -> distribute the threads and allow a nested team
  !>   'fixed'       -> CREST_CORES_PER_JOB cores per job (default 4), run
  !>                    T/cores_per_job jobs at a time, the rest queue
  !> 'auto_nested' may only be used when cores_per_job is consumed inside CREST.
  !> If cores_per_job is handed to an external program via OMP_NUM_THREADS
  !> (xtb, orca, ... subprocesses) 'auto' has to be used, because the CREST
  !> thread blocks anyway and a nested team would book the budget twice.
  !> 'fixed' is for batch loops over an ensemble that is larger than the core
  !> count, where 'auto' would degenerate to one core per job.
  integer :: T,Tdiff
  integer :: r
  character(len=64) :: val
  real(wp) :: Tfrac,Tfloor

  !> The default, all threads allocated to CREST
  T = env%threads
  parallel_jobs = max(1,T)
  cores_per_job = 1
  !> More settings, nested parallelization reset
  call omp_set_max_active_levels(1)

  select case (modus)
  case ('auto','auto_nested')
    !> distribute jobs automatically:
    !> if more cores are available than maxjobs, try to distribute remaining
    !> threads EVENLY for each job
    if (maxjobs > 0.and.T > maxjobs) then
      parallel_jobs = maxjobs
      Tfrac = real(T)/real(maxjobs)
      Tfloor = floor(Tfrac)
      cores_per_job = max(nint(Tfloor),1)
    end if
    !> never ask for more cores in total than there are threads available
    if (parallel_jobs*cores_per_job > T) then
      cores_per_job = max(1,T/parallel_jobs)
    end if
    if (index(modus,'_nested') .ne. 0 .and. cores_per_job > 1) then
      if (env%omp_allow_nested) then
        !> We should never need more than two active nested layers
        call omp_set_max_active_levels(2)
      endif
    else
#ifdef WITH_MKL
      !call mkl_free_buffers()
      call mkl_set_dynamic(0)
#endif
    end if
    call openblasset(1)

  case ('fixed')
    !>---------------------------------------------------------------------
    !> Give every job a FIXED number of cores and run as many jobs at once as
    !> fit into the budget; the OpenMP task queue keeps the remaining jobs
    !> waiting, so nothing oversubscribes the machine:
    !>     cores_per_job = CREST_CORES_PER_JOB   (default 4)
    !>     parallel_jobs = T / cores_per_job
    !>
    !> This is what an MPI-parallel external program (ORCA) needs.  The generic
    !> 'auto' scheme distributes T over "maxjobs" instead, which for a real
    !> conformer ensemble (nall >> T) means T jobs with ONE core each - an ORCA
    !> gradient on one core, hundreds of times over.  Reserving 4 cores per job
    !> and running T/4 of them is far faster and is exactly what the user asked
    !> for: "N jobs in parallel, the rest queue".
    !>
    !> If there are FEWER structures than T/cores_per_job, the leftover cores
    !> are handed to the jobs that do exist (so a 2-structure ensemble on 40
    !> cores still gets 20 cores per job instead of 4).
    !>---------------------------------------------------------------------
    call getenv('CREST_CORES_PER_JOB',val)
    read (val,*,iostat=r) cores_per_job
    if (r /= 0) cores_per_job = 4
    cores_per_job = max(1,min(cores_per_job,T))
    if (maxjobs > 0.and.T > maxjobs) then
      cores_per_job = max(cores_per_job,T/maxjobs)
    end if
    cores_per_job = max(1,min(cores_per_job,T))
    parallel_jobs = max(1,T/cores_per_job)
    call openblasset(1)

  case ('share')
    !>---------------------------------------------------------------------
    !> Run ONE job, but hand it exactly the core budget a single job gets in
    !> the parallel phase (T/maxjobs).  Used by the trial MTD: measuring with
    !> the real per-job budget makes both the trial AND the runtime estimate
    !> derived from it representative.
    !>
    !> Upstream used 'min' here (one core) "to simulate max parallelization".
    !> That made every trial MTD needlessly slow and forced the estimator to
    !> *assume* a perfect speed-up of floor(T/nmetadyn) afterwards, which is
    !> where wildly optimistic predictions came from.
    !>---------------------------------------------------------------------
    parallel_jobs = 1
    if (maxjobs > 0) then
      cores_per_job = max(1,T/maxjobs)
    else
      cores_per_job = max(1,T)
    end if
    call openblasset(1)

  case ('max')
    !> Both intern and environment variable threads to max
    parallel_jobs = T
    cores_per_job = T

  case ('min','serial')
    !> Both intern and environment variable threads to one (like a serial program)
    parallel_jobs = 1
    cores_per_job = 1

  case ('subprocess','la-focus')
    !> CREST itself uses one thread, and but the environment variable is set to max
    !> which is useful when driving a single subprocess/systemcall
    parallel_jobs = 1
    cores_per_job = T 

    !> the setting may also be used for linear-algebra focused runs, in which case 
    !> nested parallelism should be active
    if (env%omp_allow_nested) then
      !> We should never need more than two active nested layers
      call omp_set_max_active_levels(2)
    endif

  end select

  !> apply the calculated settings
  call ompmklset(parallel_jobs)
  call ompenvset(cores_per_job)
#ifdef WITH_OPENBLAS
!  if(modus.ne.'auto'.and.modus.ne.'auto_nested')then
!      call openblasset(cores_per_job)
!  endif
   call openblasset(1)
#endif
end subroutine new_ompautoset

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c get omp/mkl automatically from the global variables
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

subroutine ompgetauto(threads,omp,maxrun)
  use omp_lib
  use iomod
  implicit none
  integer,intent(inout) :: threads,omp,maxrun
  integer :: nproc
  integer :: r
  character(len=256) :: val

  call getenv('OMP_NUM_THREADS',val)
  read (val,*,iostat=r) nproc
  if (r .ne. 0) then
    nproc = 1
  end if
  threads = nproc
  maxrun = 1
  omp = nproc

end subroutine ompgetauto

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!c print omp/mkl threads that are used at the moment
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
subroutine ompprint_intern()
  use omp_lib
  implicit none
  integer :: nproc,TID

!$OMP PARALLEL PRIVATE(TID)
  TID = OMP_GET_THREAD_NUM()
  IF (TID .EQ. 0) THEN
    nproc = OMP_GET_NUM_THREADS()
    write (*,*) '============================='
    write (*,*) ' # threads =',nproc
    write (*,*) '============================='
  END IF
!$OMP END PARALLEL
end subroutine ompprint_intern

