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

!> module subprocess_engrad
!> RE-EXPORTS of subprocess engrad routines

!=========================================================================================!

module subprocess_engrad
  use generic_sc
  use xtb_sc
  use turbom_sc
  use orca_sc
  implicit none
  !>--- private module variables and parameters
  private

  !>--- generic subrpocess (run.sh)
  public :: generic_engrad

  !>--- xtb subprocess
  public :: xtb_engrad

  !>--- Turbomole-style subprocesses
  public :: turbom_engrad

  !>--- ORCA subprocesses
  public :: ORCA_engrad
  !>--- LOCAL EXTENSION: gradient-free ORCA single point, see orca_sc.f90.
  !> crest_sploop() (the "Singlepoint re-ranking" stage) never reads a gradient,
  !> so it uses this instead of ORCA_engrad and saves the CP-SCF / Z-vector
  !> solve plus the CPCM gradient contribution on every structure.
  public :: ORCA_sp

!=========================================================================================!
!=========================================================================================!
contains    !> MODULE PROCEDURES START HERE
!=========================================================================================!
!=========================================================================================!

end module subprocess_engrad
