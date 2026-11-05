!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023-2024 Philipp Pracht
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

subroutine creststop(io)
  use crest_parameters
  use crest_data
  implicit none
  integer,intent(in) :: io

  select case(io)
  case (status_normal)
    write (stdout,*) 'CREST terminated normally.'
  case default
    write (stdout,*) 'CREST terminated abnormally.'
  case ( status_error )
    write (stdout,*) 'CREST terminated with errors.'
  case ( status_ioerr )
    write (stdout,*) 'CREST terminated with I/O errors.'
  case ( status_args )
    write (stdout,*) 'CREST terminated due to invalid parameters.'
  case ( status_input )
    write (stdout,*) 'CREST terminated due to failed input file read.'
  case ( status_config )
    write (stdout,*) 'CREST terminated due to invalid configuration.'
  case ( status_failed )
    write (stdout,*) 'CREST terminated with failures.'
  case ( status_safety )
    write (stdout,*) 'Safety termination of CREST.'
  end select
!>--- LOCAL EXTENSION: exit(3) does NOT reliably flush the Fortran runtime
!> I/O buffers. Without this the terminating message above (and everything
!> printed shortly before it) can be missing from a redirected crest.out,
!> which is what makes the log look "written but never finished".
  flush (stdout)
  call exit(io)

end subroutine creststop

!================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!================================================================================!

!> NOTE: the handlers for SIGINT/SIGQUIT/SIGTERM/SIGHUP/SIGXCPU live in
!> src/signal.c (crest_sighandler / crest_install_signals), NOT here.
!> They used to be the Fortran routines wsigint/wsigquit/wsigterm below, but
!> those are NOT async-signal-safe: write(*,*), cs_shutdown() and exit() all
!> take locks (the Fortran I/O lock, the allocator lock) that some OTHER
!> thread may already hold when the signal is delivered. With OpenMP the
!> signal goes to an arbitrary thread, so this deadlocked regularly - that
!> was half of the "CREST cannot be stopped with Ctrl+C" problem.
!> The C handler only does write(2) -> kill(child process groups) -> _exit(),
!> all three of which are guaranteed async-signal-safe.
subroutine initsignal()
  implicit none
  logical,save :: signalset = .false.
  interface
    subroutine crest_install_signals() bind(C,name="crest_install_signals")
    end subroutine crest_install_signals
  end interface

  !> signal dispositions are process-wide and (BSD/Linux semantics) stay
  !> installed once they are set. This routine used to be called for every
  !> single energy/gradient evaluation, i.e. several times per MD step and
  !> thread. Every call is a syscall that takes the process-wide sighand
  !> lock, which shows up as pure kernel (system) time once many threads
  !> hammer it concurrently. Install them exactly once instead.
  if (signalset) return

  !> SIGKILL (9) can never be caught and 69 is not a Linux signal at all -
  !> both old call signal(9,...) / signal(69,...) lines were no-ops.
  call crest_install_signals()
  signalset = .true.
end subroutine initsignal

