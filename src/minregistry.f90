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
!================================================================================!

!> module minregistry_module
!> LOCAL EXTENSION: a registry of minima that have already been located, used to
!> stop an individual geometry optimisation early when it is descending into a
!> basin whose bottom is already known.
!>
!> WHY / WHAT IT IS WORTH
!> ----------------------
!> An ensemble produced by MTD contains many structures that relax into the SAME
!> minimum.  CREST finds that out with CREGEN - i.e. only AFTER every structure
!> has been optimised.  This module lets a still-running optimisation be
!> abandoned as soon as it provably heads for a known minimum, and assigns it
!> that minimum's geometry and energy.
!>
!> MEASURED ON A REAL SYSTEM (182 MTD frames of a 109-atom metal complex,
!> g-xTB relaxation trajectories, CREGEN as the arbiter of "same conformer"):
!>   * ~42% of the structures can be stopped early
!>   * equal to ~14% of all optimisation CYCLES (and ~13.7% of the wall time,
!>     measured, because the skipped tail cycles are the cheap ones)
!>   * ~4% of the stopped structures are assigned to the WRONG minimum
!> The error rate is the price of stopping before convergence; the anchor rule
!> below (see "drop") is what keeps it at ~4% instead of ~20%.
!>
!> THE ANCHOR RULE (important)
!> ---------------------------
!> Only trajectory points that lie within `drop` of their own final energy are
!> registered.  Matching against arbitrary mid-descent points was measured to be
!> far less reliable (~7% / up to 20% wrong), because a point that is still far
!> above its own minimum does not identify a basin: it only identifies "somewhere
!> on the way down".  Everything else about the scheme is standard:
!>   * only points of trajectories that CONVERGED are registered (so every
!>     endpoint is confirmed, never inferred)
!>   * a hit must satisfy BOTH  RMSD < rthr  AND  |dE| < ethr
!>   * the registry is reset for every optimisation stage, so a structure can
!>     never prune itself against its own result from the previous stage
!>
!> SWITCHES (read once, by value):
!>   CREST_PRUNE_ANCHORS=0   turn the whole thing off      (default: ON)
!>   CREST_PRUNE_RTHR        RMSD window in Angstrom       (default 0.15)
!>   CREST_PRUNE_ETHR        energy window in kcal/mol     (default 0.15)
!>   CREST_PRUNE_DROP        anchor "remaining drop" in kcal/mol (default 4.0)
!>
!> NOTE: only the ANCOPT engine is instrumented (ancopt.f90/relax).  The other
!> engines (L-BFGS, RFO, gradient descent) are simply not pruned - they keep
!> working exactly as before.
!================================================================================!
module minregistry_module
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use ls_rmsd,only:rmsd_sp
  implicit none
  private

  public :: minreg_setup
  public :: minreg_reset
  public :: minreg_active
  public :: minreg_collect
  public :: minreg_query
  public :: minreg_commit
  public :: minreg_stats
  public :: minreg_keep_trj

  !> how many sample points / endpoints can be held.  The registry is reset for
  !> every optimisation stage, and each converged structure contributes only the
  !> tail of its trajectory (~10-40 points), so this is generous.
  integer,parameter :: maxreg = 12000
  integer,parameter :: maxend = 2000

  !>--- shared registry.  Everything below is only ever touched inside the
  !> NAMED critical section "minreg" (see the routines); crest_oloop runs the
  !> optimisations as OpenMP tasks, so this has to be thread safe.
  integer :: nreg = 0
  integer :: natreg = 0
  real(wp),allocatable :: reg_xyz(:,:,:)  !> (3,nat,nreg) sample point
  real(wp),allocatable :: reg_e(:)        !> its energy
  integer,allocatable  :: reg_t(:)        !> -> endpoint slot of its trajectory

  integer :: nend = 0
  real(wp),allocatable :: end_xyz(:,:,:)  !> (3,nat,nend) bottom of a basin
  real(wp),allocatable :: end_e(:)

  !>--- settings
  logical  :: reg_on   = .true.
  real(wp) :: reg_rthr = 0.15_wp
  real(wp) :: reg_ethr = 0.15_wp/627.5095_wp
  real(wp) :: reg_drop = 4.0_wp/627.5095_wp
  logical  :: reg_set  = .false.

  !>--- counters, reported once per optimisation stage by minreg_stats
  integer :: npruned = 0
  integer :: ncommit = 0
  integer :: nquery = 0
  integer :: nregadd = 0
  !>--- LOCAL EXTENSION (Δ-ML): dump optimisation trajectories
  logical  :: keep_trj = .false.
  real(wp) :: keep_trj_e = 40.0_wp/627.5095_wp   !> window in Eh
  integer  :: ntrjdump = 0
!========================================================================================!
!========================================================================================!
contains
!========================================================================================!
!========================================================================================!

!> read the switches once
  subroutine minreg_setup()
    implicit none
    character(len=64) :: env
    integer :: io,ival
    real(wp) :: rval

    if (reg_set) return
    reg_set = .true.

    call getenv('CREST_PRUNE_ANCHORS',env)
    ival = 1
    read (env,*,iostat=io) ival
    if (io == 0) reg_on = (ival /= 0)

!>--- LOCAL EXTENSION (Δ-ML): -keep_opt_trj sets these
    call getenv('CREST_KEEP_OPT_TRJ',env)
    ival = 0
    read (env,*,iostat=io) ival
    if (io == 0) keep_trj = (ival /= 0)

    call getenv('CREST_OPT_TRJ_E',env)
    rval = keep_trj_e*627.5095_wp
    read (env,*,iostat=io) rval
    if (io == 0.and.rval > 0.0_wp) keep_trj_e = rval/627.5095_wp

    call getenv('CREST_PRUNE_RTHR',env)
    rval = reg_rthr
    read (env,*,iostat=io) rval
    if (io == 0.and.rval > 0.0_wp) reg_rthr = rval

    call getenv('CREST_PRUNE_ETHR',env)
    rval = reg_ethr*627.5095_wp
    read (env,*,iostat=io) rval
    if (io == 0.and.rval > 0.0_wp) reg_ethr = rval/627.5095_wp

    call getenv('CREST_PRUNE_DROP',env)
    rval = reg_drop*627.5095_wp
    read (env,*,iostat=io) rval
    if (io == 0.and.rval >= 0.0_wp) reg_drop = rval/627.5095_wp

    return
  end subroutine minreg_setup

!========================================================================================!

!>--- Δ-ML extension: TRUE when we must RECORD the trajectory, which is needed
!> both for anchor pruning (reg_on) and for -keep_opt_trj (keep_trj).
!> They are deliberately decoupled: all production scripts set
!> CREST_PRUNE_ANCHORS=0, and the trajectory must still be collected then.
  logical function minreg_collect() result(on)
    implicit none
    call minreg_setup()
    on = reg_on .or. keep_trj
    return
  end function minreg_collect

  logical function minreg_active() result(on)
    implicit none
    call minreg_setup()
    on = reg_on
    return
  end function minreg_active

!========================================================================================!
!> start a fresh registry for one optimisation stage
  subroutine minreg_reset(nat)
    implicit none
    integer,intent(in) :: nat

    call minreg_setup()

    !$omp critical(minreg)
    if (.not.allocated(reg_xyz).or.natreg /= nat) then
      if (allocated(reg_xyz)) deallocate (reg_xyz)
      if (allocated(reg_e)) deallocate (reg_e)
      if (allocated(reg_t)) deallocate (reg_t)
      if (allocated(end_xyz)) deallocate (end_xyz)
      if (allocated(end_e)) deallocate (end_e)
      allocate (reg_xyz(3,nat,maxreg),source=0.0_wp)
      allocate (reg_e(maxreg),source=0.0_wp)
      allocate (reg_t(maxreg),source=0)
      allocate (end_xyz(3,nat,maxend),source=0.0_wp)
      allocate (end_e(maxend),source=0.0_wp)
      natreg = nat
    end if
    nreg = 0
    nend = 0
    npruned = 0
    ncommit = 0
    nquery = 0
    nregadd = 0
    !$omp end critical(minreg)

    return
  end subroutine minreg_reset

!========================================================================================!
!> Called once an optimisation has CONVERGED on its own minimum: store the
!> minimum, plus every point of its trajectory that lies within `reg_drop` of
!> it.  Those points are the anchors later structures may be matched against.
  subroutine minreg_commit(nat,trjxyz,trje,ntrj,endxyz,ende)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: ntrj
    real(wp),intent(in) :: trjxyz(3,nat,ntrj)
    real(wp),intent(in) :: trje(ntrj)
    real(wp),intent(in) :: endxyz(3,nat)
    real(wp),intent(in) :: ende

    integer :: k,it

    if (.not.reg_on) return
    if (ntrj < 1) return
    if (.not.allocated(reg_xyz)) return

    !$omp critical(minreg)
    if (nend < maxend) then
      nend = nend+1
      end_xyz(:,:,nend) = endxyz
      end_e(nend) = ende
      it = nend
      do k = 1,ntrj
        if (trje(k)-ende > reg_drop) cycle
        if (nreg >= maxreg) exit
        nreg = nreg+1
        reg_xyz(:,:,nreg) = trjxyz(:,:,k)
        reg_e(nreg) = trje(k)
        reg_t(nreg) = it
        nregadd = nregadd+1
      end do
      ncommit = ncommit+1
    end if
    !$omp end critical(minreg)

    return
  end subroutine minreg_commit

!========================================================================================!
!> Is this point already known to be heading for a located minimum?
!> Returns .true. and the minimum's geometry/energy if so.
  logical function minreg_query(nat,xyz,e,qxyz,qe) result(hit)
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(in) :: xyz(3,nat)
    real(wp),intent(in) :: e
    real(wp),intent(out) :: qxyz(3,nat)
    real(wp),intent(out) :: qe

    real(wp) :: U(3,3),xc(3),yc(3),err,rbest
    real(wp),allocatable :: gdum(:,:)
    integer :: i,ib

    hit = .false.
    qxyz = xyz
    qe = e

    if (.not.reg_on) return
    if (.not.allocated(reg_xyz)) return
    if (nreg < 1) return

    allocate (gdum(3,nat),source=0.0_wp)

    !$omp critical(minreg)
    ib = 0
    rbest = huge(1.0_wp)
    do i = 1,nreg
      !>--- cheap scalar rejection first: the vast majority is thrown out here
      if (abs(reg_e(i)-e) >= reg_ethr) cycle
      call rmsd_sp(nat,xyz,reg_xyz(:,:,i),1,U,xc,yc,err,.false.,gdum)
      if (err < rbest) then
        rbest = err
        ib = i
      end if
    end do
    if (ib > 0.and.rbest < reg_rthr) then
      qxyz = end_xyz(:,:,reg_t(ib))
      qe = end_e(reg_t(ib))
      hit = .true.
      npruned = npruned+1
    end if
    !$omp end critical(minreg)

    nquery = nquery+1

    deallocate (gdum)

    return
  end function minreg_query

!========================================================================================!
!> one summary line for the log, printed at the end of every optimisation stage
!>--- LOCAL EXTENSION (Δ-ML): write the structures visited during this
!> optimisation to 'enopt_traj.xyz', keeping only those within keep_trj_e of the
!> trajectory's END energy.  Ideas: the points near a minimum are exactly the
!> region where a Δ-ML correction is needed and where its gradient is most
!> trustworthy, and they are real geometries (not MD snapshots).
  subroutine minreg_keep_trj(nat,at,trjxyz,trje,ntrj,eend)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    integer,intent(in) :: ntrj
    real(wp),intent(in) :: trjxyz(3,nat,*)
    real(wp),intent(in) :: trje(*)
    real(wp),intent(in) :: eend

    integer :: k,i,iu,io
    real(wp) :: de
    character(len=2),parameter :: sym(1:86) = [character(len=2) :: &
      'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si', &
      'P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni', &
      'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo', &
      'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba', &
      'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
      'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po', &
      'At','Rn']

    call minreg_setup()
    if (.not.keep_trj) return

    !$omp critical (minreg_trj_io)
    open(newunit=iu,file='enopt_traj.xyz',access='append',status='unknown', &
         action='write',iostat=io)
    if (io == 0) then
      do k = 1,ntrj
        de = trje(k) - eend
        if (de > keep_trj_e) cycle
        write(iu,'(i0)') nat
        write(iu,'(1x,a,f20.12,1x,a,f10.4)') 'Etot=',trje(k),'dE_kcal=',de*627.5095_wp
        do i = 1,nat
          if (at(i) >= 1 .and. at(i) <= 86) then
            write(iu,'(1x,a2,3f16.8)') sym(at(i)), &
                 trjxyz(1,i,k)*0.529177210903_wp, &
                 trjxyz(2,i,k)*0.529177210903_wp, &
                 trjxyz(3,i,k)*0.529177210903_wp
          end if
        end do
        ntrjdump = ntrjdump+1
      end do
      close(iu)
    end if
    !$omp end critical (minreg_trj_io)
  end subroutine minreg_keep_trj

  subroutine minreg_stats(iunit)
    implicit none
    integer,intent(in) :: iunit
    integer :: ip,nc,nrg,nq
    if (.not.reg_on) return
    !$omp critical(minreg)
    ip = npruned
    nc = ncommit
    nrg = nregadd
    nq = nquery
    !$omp end critical(minreg)
    write (iunit,'(1x,a,i0,a,i0,a,i0,a,i0)') &
    & 'Anchor pruning: minima registered: ',nc,', anchor points: ',nrg, &
    & ', optimisations stopped early: ',ip,', queries: ',nq
    return
  end subroutine minreg_stats

!========================================================================================!
!========================================================================================!
end module minregistry_module
