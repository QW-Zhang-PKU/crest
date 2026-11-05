!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

subroutine crest_moleculardynamics(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  real(wp) :: energy,gnorm

  character(len=80) :: atmp
  character(len=*),parameter :: trjf='crest_dynamics.xyz'
!========================================================================================!
  write(stdout,*)
  !call system('figlet dynamics')
  write(stdout,*) "      _                             _           " 
  write(stdout,*) "   __| |_   _ _ __   __ _ _ __ ___ (_) ___ ___  "
  write(stdout,*) "  / _` | | | | '_ \ / _` | '_ ` _ \| |/ __/ __| "
  write(stdout,*) " | (_| | |_| | | | | (_| | | | | | | | (__\__ \ "
  write(stdout,*) "  \__,_|\__, |_| |_|\__,_|_| |_| |_|_|\___|___/ "
  write(stdout,*) "        |___/                                   "
  write(stdout,*)
!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()
  call tim%start(14,'Molecular dynamics (MD)')
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!

  pr = .true.
  !>--- default settings from env
  call env_to_mddat(env)
  mddat = env%mddat
  calc = env%calc
  !>--- check if we have any MD & calculation settings allocated
  if (.not. mddat%requested) then
    write (stdout,*) 'MD requested, but no MD settings present.'
    env%iostatus_meta = status_config
    return
  else if (calc%ncalculations < 0) then
    write (stdout,*) 'MD requested, but no calculation settings present.'
    env%iostatus_meta = status_config
    return
  end if

  !>--- print calculation info
  call calc%info( stdout )

  !>--- init SHAKE? --> we need connectivity info
  if (mddat%shake) then
    calc%calcs(1)%rdwbo = .true.
    if(.not.calc%calcs(1)%active) calc%calcs(1)%active=.true.
!>--- LOCAL CHANGE (2026-09-26): only the WBOs are needed, and the bonds come
!> from ORCA.out (rd_ORCA_wbo), not from the gradient file.  engrad() filled
!> `grad`, which was allocated, passed and deallocated without ever being read
!> - for ORCA that is a whole CP-SCF/Z-vector solve plus the CPCM gradient term
!> per standalone MD run.  singlepoint() gives identical WBOs and energy.
    call singlepoint(mol,calc,energy,io)
    calc%calcs(1)%rdwbo = .false.
    call move_alloc(calc%calcs(1)%wbo,mddat%shk%wbo)
  end if

  !>--- complete real-time settings to steps
  mddat%trajectoryfile = trjf

  !>--- run the MD
  call dynamics(mol,mddat,calc,pr,io)

  select case (io)
  case (0)
    write (stdout,*) 'MD run completed successfully'
    write (stdout,*) 'Trajectory written to ',trjf
!>--- LOCAL EXTENSION: term = 3 is the deliberate "stop" file, NOT a failure.
!> It used to fall into the else branch below and was reported as
!> "terminated ABNORMALLY" with status_failed, i.e. the trajectory the user
!> asked for was written but crest claimed the run had failed.
  case (3)
    write (stdout,*) 'MD run stopped early on request (stop file found)'
    write (stdout,*) 'Trajectory written to ',trjf
  case default
    write (stdout,*) 'WARNING: MD run terminated ABNORMALLY'
    env%iostatus_meta = status_failed
  end select
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_moleculardynamics
