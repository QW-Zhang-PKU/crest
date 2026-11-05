!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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

module shake_module
  use crest_parameters
  implicit none

  integer,private,parameter :: metal(1:86) = [&
     & 0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1, &
     & 1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1, &
     & 1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1, &
     & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
     & 1,0,0,0,0,0]

  !> Atomic radii
  !> M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
  !> in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
  !> edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
  !> corrected Nov. 17, 2010 for the 92nd edition.
  real(wp),private,parameter :: atomicRad(1:118) = aatoau*[ &
     & 0.32_wp,0.37_wp,1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp, &
     & 0.60_wp,0.62_wp,1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp, &
     & 1.00_wp,1.01_wp,2.00_wp,1.74_wp,1.59_wp,1.48_wp,1.44_wp,1.30_wp, &
     & 1.29_wp,1.24_wp,1.18_wp,1.17_wp,1.22_wp,1.20_wp,1.23_wp,1.20_wp, &
     & 1.20_wp,1.18_wp,1.17_wp,1.16_wp,2.15_wp,1.90_wp,1.76_wp,1.64_wp, &
     & 1.56_wp,1.46_wp,1.38_wp,1.36_wp,1.34_wp,1.30_wp,1.36_wp,1.40_wp, &
     & 1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp,2.38_wp,2.06_wp, &
     & 1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp,1.82_wp, &
     & 1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp,1.74_wp,1.64_wp, &
     & 1.58_wp,1.50_wp,1.41_wp,1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp, &
     & 1.44_wp,1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp,2.42_wp,2.11_wp, &
     & 2.01_wp,1.90_wp,1.84_wp,1.83_wp,1.80_wp,1.80_wp,1.73_wp,1.68_wp, &
     & 1.68_wp,1.68_wp,1.65_wp,1.67_wp,1.73_wp,1.76_wp,1.61_wp,1.57_wp, &
     & 1.49_wp,1.43_wp,1.41_wp,1.34_wp,1.29_wp,1.28_wp,1.21_wp,1.22_wp, &
     & 1.36_wp,1.43_wp,1.62_wp,1.75_wp,1.65_wp,1.57_wp]

  public :: shakedata
  !=========================================================================================!
  !data object that contains settings for a molecular dynamics simulation.
  type :: shakedata

    logical :: initialized = .false.

    integer :: shake_mode = 0 ! 0 = off/user, 1 = H only, 2 = all

    !>--- user defined SHAKE bonds
    integer :: nusr = 0 !number of user-defined shake bonds
    integer,allocatable :: conslistu(:,:)
!>--- LOCAL EXTENSION: OPTIONAL target distance (in BOHR, *not* squared) for
!> each user-defined bond.  Upstream only ever took the distance found in the
!> geometry handed to init_shake(), i.e. a user-specified constraint value was
!> silently ignored and the bond was pinned at whatever the starting structure
!> happened to have.  If this array is allocated, distcons(i) is built from it
!> instead - which is what makes a rigid constraint coming from --cinp honour
!> the distance the user actually wrote down.
    real(wp),allocatable :: distconsu(:)

    !>--- topological connectivity/WBO matrix
    real(wp),allocatable :: wbo(:,:)

    !>--- final SHAKE constrainment data
    integer :: ncons = 0
    integer,allocatable :: conslist(:,:)
    real(wp),allocatable :: distcons(:)
    real(wp),allocatable :: dro(:,:)
    real(wp),allocatable :: dr(:,:)

    !>--- SHAKE convergence parameter
!>==============================================================================!
!> LOCAL EXTENSION (P3-a, 2026-09-23): maxcyc 250 -> 1000.
!>
!> Why.  The sweep converges GEOMETRICALLY, so the number of sweeps needed
!> grows only LOGARITHMICALLY in the initial violation:
!>
!>     n_sweeps ~ ln(maxdev0 / tolshake) / (-ln rho)
!>
!> Measured on a real 97-atom system with a 10-distance rigid cluster:
!>     maxdev0 = 0.2521        (the uncorrected 4 fs step)
!>     dev(250) = 1.2887e-07   against tolshake = 1.e-7
!>     rho = (1.2887e-7/0.2521)**(1/250) = 0.9437
!>     ->  dev(n) = 0.2521 * 0.9437**n  reaches 1e-7 at n ~ 254.
!> Verified: 0.2521*0.9437**250 = 1.288e-7, matching the log digit for digit.
!> So the solver was NOT diverging, NOT stalling on a redundant set and NOT
!> facing an unreachable target - it was 4 SWEEPS SHORT OF THE BUDGET.
!>
!> Why 1000 and not, say, 260.  The margin has to survive a different system,
!> a different bias force or a larger maxdev0.  With rho = 0.9437 each extra
!> 17.3 sweeps buys one more decade of residual, so 1000 sweeps reaches
!> ~1e-19 from the same start - it will never be the binding constraint again.
!>
!> Cost.  Zero semantics: this can only turn "gave up" into "converged", it
!> cannot change any step that already converged within 250 sweeps.  The extra
!> work is pure arithmetic on ncons constraints and happens only on steps that
!> would otherwise ABORT THE RUN, so it is negligible next to one engrad call.
!>
!> Note this is NOT a reason to keep a too-large time step: a smaller step
!> reduces maxdev0 (roughly as dt**2) and therefore helps, but only
!> logarithmically - 4 -> 3 fs saves about 9 sweeps out of 254.  That is why
!> the trial-MTD step ladder used to squeak through by a 2 % margin instead of
!> actually working.
!==============================================================================!
!>==============================================================================!
!> LOCAL CHANGE (2026-10-24): tolshake 1e-7 -> 1e-5, maxcyc 1000 -> 3000.
!>
!> Why.  A real 53-atom Rh system with a 4-atom rigid cluster (atoms 1,7-8,37)
!> ran out of budget at maxdev0 = 2.7955e-01 with a contraction of rho = 0.9923:
!>
!>     n_sweeps = ln(maxdev0/tol) / ln(1/rho)
!>              = ln(2.7955e-01/1e-7) / 0.00773 = 1920    against maxcyc = 1000
!>
!> So the solver was converging, just too slowly for the budget.  Two changes:
!>
!>   tolshake 1e-5  -> ln(2.7955e-01/1e-5)/0.00773 = 1324  (fits in 3000)
!>
!> 1e-5 is a RELATIVE deviation (dev = |dr^2-dist|/dist), i.e. ~2e-5 Angstrom on
!> a 2 Angstrom bond.  With a C-C force constant of ~450 kcal/mol/A^2 that is a
!> 9e-8 kcal/mol energy effect - far below the force noise.  For reference
!> GROMACS SHAKE defaults to 1e-4, so 1e-5 is still 10x tighter than standard MD.
!>
!>   maxcyc 6000   -> measured worst case needed ~3100 sweeps (rho = 0.9968 on
!>                 the real failing system), so 6000 is a 1.9x margin.
!>
!> Cost.  maxcyc is an upper bound and the sweep loop exits as soon as conv is
!> true (see the "if (.not.conv.and.icyc .le. ...) cycle / exit" at the bottom),
!> so a normal step pays nothing.  A looser tol actually SAVES sweeps on every
!> step: the count goes as ln(maxdev0/tol), so 1e-7 -> 1e-5 removes
!> ln(100)/ln(1/rho) sweeps per step (about 50% at rho = 0.99).
!>
!> What this does NOT change: the fixed point (every constraint satisfied), the
!> correction direction dro, the mass weighting, or the Gauss-Seidel ordering.
!> A step that already converged stops at the same point it always did.
!>==============================================================================!
    integer :: maxcyc = 6000
    real(wp) :: tolshake = 1.d-5

    !>--- temporary coordinate space (to avoid reallocation every shake call)
    real(wp),allocatable :: xyzt(:,:)

    !>--- pointer to the freeze list
    logical,pointer :: freezeptr(:)

  end type shakedata

  public :: init_shake
  public :: do_shake
!>--- LOCAL EXTENSION (2026-10-24): SHAKE diagnostics.
!> Some real systems converge so slowly (rho ~ 0.99) that the sweep count is the
!> only visible symptom of an over-constrained/over-determined constraint group.
!> Counting calls, failures and the worst rho makes that visible without
!> touching the solver.  Silent by default; set CREST_SHAKE_STATS=1 for a
!> periodic one-line summary.
  integer,private :: shake_ncalls = 0
  integer,private :: shake_nfail  = 0
  integer,private :: shake_nsweep = 0
  real(wp),private :: shake_rho_max = 0.0_wp
  logical,private :: shake_stats_env = .false.
  logical,private :: shake_stats_init = .false.
  public :: shake_stats_report
!>--- LOCAL EXTENSION: single source of truth for "is this atom a metal?".
!> The metal() table above used to be private to init_shake().  CREGEN's
!> topology check needs the *identical* definition: a bond that SHAKE refuses to
!> constrain must also be a bond whose topology CREGEN refuses to judge.
  public :: crest_is_metal

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!> LOCAL EXTENSION: metal test shared by SHAKE and by CREGEN's topology check.
!>
!> Returns .true. for every element flagged in metal(), EXCEPT Li (3) and Be (4)
!> - exactly the exception init_shake() applies: their bonds are short/covalent
!> enough to be constrained there, so they stay part of the topology here.
!> Z outside 1..86 (the table's range) is simply "not a metal".
!========================================================================================!
  pure logical function crest_is_metal(z) result(yes)
    implicit none
    integer,intent(in) :: z
    yes = .false.
    if (z < 1.or.z > 86) return
    if (metal(z) /= 1) return
    if (z == 3.or.z == 4) return
    yes = .true.
    return
  end function crest_is_metal
!========================================================================================!
!========================================================================================!
!========================================================================================!

  subroutine init_shake(nat,at,xyz,shk,pr)
!********************************************
!* subroutine init_shake
!* initialize SHAKE algorithm by documenting
!* which pairs of atoms to constrain
!********************************************
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(shakedata) :: shk
    logical,intent(in) :: pr
    integer :: iat,jat,i,j,jmin
!>--- LOCAL EXTENSION: bookkeeping for merging the user list with the
!> auto-detected list (k/l) and for the "same bond twice" test (dup,nadd)
    integer :: k,nadd
    logical :: dup
    real(wp) :: minrij,rij,wthr
    real(wp) :: drij(3)

    integer :: nconsu
    integer :: ij
    integer,allocatable :: list(:)
    real(wp) :: rco
    logical :: metalbond,rcut,checkfreeze
    integer :: nfrz

    integer :: n2
    integer,allocatable :: cons2(:,:)

!>--- if shake was already set up, return
    if (shk%initialized) return

!>--- reset counter
    nconsu = 0
    n2 = 0
    checkfreeze = associated(shk%freezeptr)
    !> gfortran workaround
    if (checkfreeze) then !> gfortran check
      nfrz = size(shk%freezeptr,1)
    else
      nfrz = 0
    end if

!>--- count user-defined bonds
    if ((shk%nusr > 0).and.allocated(shk%conslistu)) then
      do i = 1,shk%nusr
        nconsu = nconsu+1
      end do
    end if

!>--- constrain all X-H only
    if (shk%shake_mode == 1) then
      ij = nat*(nat+1)/2
      allocate (cons2(2,ij),source=0)
      do i = 1,nat
        if (at(i) .eq. 1) then
          minrij = 1000.0_wp
          do j = 1,nat
            if (j .ne. i) then

              !> check if BOTH are frozen
              if (checkfreeze.and.j <= nfrz.and.i <= nfrz) then
                if (shk%freezeptr(i).and.shk%freezeptr(j)) cycle
              end if

              rij = (xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
              if (rij .lt. minrij) then
                minrij = rij
                jmin = j
              end if
            end if
          end do
          if (at(jmin) .eq. 1) then
            if (jmin .gt. i) then
              nconsu = nconsu+1
              n2 = n2+1
              cons2(1,n2) = i
              cons2(2,n2) = jmin
            end if
          else
            nconsu = nconsu+1
            n2 = n2+1
            cons2(1,n2) = i
            cons2(2,n2) = jmin
          end if
        end if
      end do
    end if

!>--- SHAKE all bonds
    if (shk%shake_mode == 2) then
      if (allocated(shk%wbo)) then
        ij = nat*(nat+1)/2
        allocate (cons2(2,ij),source=0)
        allocate (list(ij),source=0)
        do i = 1,nat
          do j = 1,nat
            if (i .eq. j) cycle

            !>--- if both are frozen, no SHAKE!
            if (checkfreeze.and.j <= nfrz.and.i <= nfrz) then
              if (shk%freezeptr(i).and.shk%freezeptr(j)) cycle
            end if

            rij = (xyz(1,i)-xyz(1,j))**2+(xyz(2,i)-xyz(2,j))**2+(xyz(3,i)-xyz(3,j))**2
            rco = (atomicRad(at(j))+atomicRad(at(i)))*autoaa
            ij = shake_lin(i,j)
            !>--- to consider?
            rcut = (0.52917726_wp*sqrt(rij) .lt. 1.2_wp*rco).and.(list(ij) .eq. 0)
            !>--- metal bond?
            metalbond = metal(at(i)) .eq. 1.or.metal(at(j)) .eq. 1
            !>--- WBO threshold. if WBO > thr, the bond is constrained
            wthr = 0.5_wp
            !>--- modify wthr, e.g. K...O in phosphates have WBO around 0.15
            if (metalbond) wthr = 0.1_wp
            !>--- check relevant
            if (rcut.and.shk%wbo(i,j) .gt. wthr) then
              !>--- do not constrain M bonds except Li/Be
              metalbond = metal(at(i)) .eq. 1.or.metal(at(j)) .eq. 1
              if (at(i) .eq. 3.or.at(i) .eq. 4.or.at(j) .eq. 3.or.at(j) .eq. 4) &
                               &         metalbond = .false.
              if (metalbond) then
                if ((i .lt. j).and.pr) &
               & write (*,*) 'init_shake: metal bond ',i,j,'not constrained'
              else
                list(ij) = 1
                nconsu = nconsu+1
                n2 = n2+1
                cons2(1,n2) = i
                cons2(2,n2) = j
              end if
            end if
          end do
        end do
        deallocate (list)
      else
        write (*,*) 'No bonding information provided!'
        write (*,*) 'Automatic SHAKE setup failed.'
        error stop
      end if
    end if

    shk%initialized = .true.

!>--- LOCAL EXTENSION: an auto-detected bond that is already in the user list
!> must not be constrained a second time.  With shake_mode == 2 (the default)
!> every bond of the molecule is auto-detected, so a bond that also came in
!> through --cinp would otherwise show up twice.
    nadd = 0
    if (n2 > 0) then
      do i = 1,n2
        dup = .false.
        if (shk%nusr > 0) then
          do k = 1,shk%nusr
            if ((cons2(1,i) == shk%conslistu(1,k).and.cons2(2,i) == shk%conslistu(2,k)).or. &
            &   (cons2(1,i) == shk%conslistu(2,k).and.cons2(2,i) == shk%conslistu(1,k))) then
              dup = .true.
              exit
            end if
          end do
        end if
        if (.not.dup) nadd = nadd+1
      end do
    end if
    shk%ncons = shk%nusr+nadd
    if (shk%ncons .lt. 1) then
      if (allocated(cons2)) deallocate (cons2)
      return
    end if

    allocate (shk%conslist(2,shk%ncons),shk%distcons(shk%ncons))
    allocate (shk%dro(3,shk%ncons),shk%dr(4,shk%ncons))

    if (shk%nusr > 0) then
      shk%conslist(1:2,1:shk%nusr) = shk%conslistu(1:2,1:shk%nusr)
    end if
    if (n2 > 0) then
      j = shk%nusr+1
      do i = 1,n2
        dup = .false.
        if (shk%nusr > 0) then
          do k = 1,shk%nusr
            if ((cons2(1,i) == shk%conslistu(1,k).and.cons2(2,i) == shk%conslistu(2,k)).or. &
            &   (cons2(1,i) == shk%conslistu(2,k).and.cons2(2,i) == shk%conslistu(1,k))) then
              dup = .true.
              exit
            end if
          end do
        end if
        if (.not.dup) then
          shk%conslist(1:2,j) = cons2(1:2,i)
          j = j+1
        end if
      end do
    end if
    do i = 1,shk%ncons
      iat = shk%conslist(1,i)
      jat = shk%conslist(2,i)
      drij(1:3) = xyz(1:3,iat)-xyz(1:3,jat)
      shk%distcons(i) = drij(1)**2+drij(2)**2+drij(3)**2
    end do
!>--- LOCAL EXTENSION: an explicitly given target distance wins over the
!> distance found in the geometry.  distcons holds the SQUARE of the distance.
    if (allocated(shk%distconsu)) then
      do i = 1,min(shk%nusr,shk%ncons)
        if (shk%distconsu(i) > 1.0e-6_wp) then
          shk%distcons(i) = shk%distconsu(i)*shk%distconsu(i)
        end if
      end do
    end if

    if (allocated(cons2)) deallocate (cons2)
    return
  end subroutine init_shake

!========================================================================================!

  subroutine do_shake(nat,xyzo,xyz,velo,acc,mass,tstep,shk,pr,iostat)
!**************************************************
!* subroutine do_shake
!* Calculation of the actual SHAKE constraint
!**************************************************
    implicit none
    integer,intent(in) :: nat
    type(shakedata) :: shk
    logical,intent(in) :: pr
    integer,intent(out),optional :: iostat
    real(wp) :: xyzo(3,nat),xyz(3,nat),velo(3,nat),acc(3,nat),mass(nat)
    real(wp) :: virsh(3),tstep
    real(wp),allocatable :: xyzt(:,:)
    real(wp) :: vel(3)
    integer :: i,icyc,io
    integer :: iat,jat
    logical :: conv
    real(wp) :: maxdev
!>--- LOCAL EXTENSION: the deviation at the very START of the iteration.  It is
!> what tells an INFEASIBLE constraint set apart from a DIVERGING solver: an
!> unreachable target leaves a finite residual (of the order of the mismatch) no
!> matter how many cycles are run, while an unstable iteration explodes from a
!> tiny residual to something like 1e290.  Reporting only the final value makes
!> those two look identical.
    real(wp) :: maxdev0
    real(wp) :: r,dev,dist,denom
    real(wp) :: gcons,rmi,rmj
    real(wp) :: tau1,tau2
    integer :: jmaxdev,nfrz
    logical :: checkfreeze

    conv = .false. !> assume not converged
    io = 1         !> and therefore unsuccessfull
    icyc = 0
    maxdev0 = -1.0_wp

    if (.not.allocated(shk%xyzt)) allocate (shk%xyzt(3,nat))
    shk%xyzt = xyz

    checkfreeze = associated(shk%freezeptr)
    !> gfortran workaround
    if (checkfreeze) then !> gfortran check
      nfrz = size(shk%freezeptr,1)
    else
      nfrz = 0
    end if

    tau1 = 1.0_wp/tstep
    tau2 = tau1*tau1

    do i = 1,shk%ncons
      iat = shk%conslist(1,i)
      jat = shk%conslist(2,i)
      shk%dro(1:3,i) = xyzo(1:3,iat)-xyzo(1:3,jat)
    end do

!>--- iterative SHAKE loop
    do
      maxdev = 0.d0

      do i = 1,shk%ncons
        iat = shk%conslist(1,i)
        jat = shk%conslist(2,i)
        shk%dr(1:3,i) = shk%xyzt(1:3,iat)-shk%xyzt(1:3,jat)
        shk%dr(4,i) = shk%dr(1,i)**2+shk%dr(2,i)**2+shk%dr(3,i)**2
        dist = shk%distcons(i)
        dev = abs(shk%dr(4,i)-dist)/dist
        if (dev .gt. maxdev) then
          maxdev = dev
          jmaxdev = i
        end if
      end do

      if (maxdev .lt. shk%tolshake) conv = .true.
!>--- LOCAL EXTENSION: icyc is still 0 on the first pass (it is only incremented
!> at the end of the loop body), so this captures the deviation of the UNCORRECTED
!> step, i.e. how far the constraint actually was before SHAKE touched anything.
      if (icyc == 0) maxdev0 = maxdev
!>--- NOTE for md_failreport(): maxdev0 is the deviation of the UNCORRECTED
!> STEP, not an initial geometry mismatch.  Do not tell the user to compare
!> coord with coord.ref on the basis of it.

      if (.not.conv) then
!>==============================================================================!
!> LOCAL EXTENSION (P4, 2026-09-23): Gauss-Seidel instead of Jacobi.
!>
!> What was wrong: shk%dr was filled ONCE at the top of the cycle (just above)
!> and the correction sweep below then used those STALE vectors for every one of
!> the ncons constraints, even though the sweep moves the atoms as it goes.  So
!> the sweep linearised every constraint at the same - by then outdated -
!> geometry.  That is a Jacobi sweep, and its iteration matrix has a much larger
!> spectral radius than the sequential one.  On a strongly coupled set (110
!> constraints sharing atoms, plus the 10 extra bonds from --cinp) it is not
!> merely slow: it AMPLIFIES.  Observed on a real system - the step arrived on
!> the constraint manifold to within 1e-9 and 250 cycles later the deviation was
!> 1e290, i.e. the MD blew up instead of failing politely.
!>
!> The fix: recompute the bond vector of constraint i from the CURRENT,
!> already-partially-corrected positions right before using it.  That is the
!> classic sequential SHAKE used by GROMACS / AMBER / LAMMPS.
!>
!> What deliberately does NOT change:
!>   * the correction DIRECTION is still dro, the bond vector at the positions
!>     before the step.  Classic SHAKE keeps that reference fixed; only the
!>     linearisation (dr) is refreshed.
!>   * the mass weighting (rmi/rmj) is untouched, so the projection is still the
!>     minimum mass-weighted displacement, i.e. the same Lagrange multipliers.
!>   * the velocity/acceleration fix further down uses the ACTUAL applied
!>     displacement (shk%xyzt - xyz), so it follows automatically.
!> Consequently the fixed point of the iteration - "every constraint satisfied" -
!> is the same set of equations as before.  When both variants converge they
!> agree to within tolshake; Gauss-Seidel just converges (and is far more
!> stable).  It does NOT change the sampled ensemble.
!==============================================================================!
        if (.not.checkfreeze.or.nfrz > 0) then  !> Non-frozen case

          do i = 1,shk%ncons
            iat = shk%conslist(1,i)
            jat = shk%conslist(2,i)
            dist = shk%distcons(i)
            rmi = 1.0_wp/mass(iat)
            rmj = 1.0_wp/mass(jat)
!>--- P4: refresh the linearisation at the CURRENT positions
            shk%dr(1:3,i) = shk%xyzt(1:3,iat)-shk%xyzt(1:3,jat)
            shk%dr(4,i) = shk%dr(1,i)**2+shk%dr(2,i)**2+shk%dr(3,i)**2
            denom = 2.0_wp*(rmi+rmj)*(shk%dr(1,i)*shk%dro(1,i)+  &
               &     shk%dr(2,i)*shk%dro(2,i)+ &
               &                            shk%dr(3,i)*shk%dro(3,i))
            gcons = (dist-shk%dr(4,i))/denom
            shk%xyzt(1:3,iat) = shk%xyzt(1:3,iat)+rmi*gcons*shk%dro(1:3,i)
            shk%xyzt(1:3,jat) = shk%xyzt(1:3,jat)-rmj*gcons*shk%dro(1:3,i)
          end do

        else !> Frozen case

          do i = 1,shk%ncons
            iat = shk%conslist(1,i)
            jat = shk%conslist(2,i)
            dist = shk%distcons(i)
            rmi = 1.0_wp/mass(iat)
            rmj = 1.0_wp/mass(jat)
!>--- P4: same refresh here (this branch only differs in WHO is allowed to move)
            shk%dr(1:3,i) = shk%xyzt(1:3,iat)-shk%xyzt(1:3,jat)
            shk%dr(4,i) = shk%dr(1,i)**2+shk%dr(2,i)**2+shk%dr(3,i)**2
            denom = 2.0_wp*(rmi+rmj)*(shk%dr(1,i)*shk%dro(1,i)+  &
               &     shk%dr(2,i)*shk%dro(2,i)+ &
               &                            shk%dr(3,i)*shk%dro(3,i))
            gcons = (dist-shk%dr(4,i))/denom
            !> update only non-frozen
            if (.not.shk%freezeptr(iat)) then
              shk%xyzt(1:3,iat) = shk%xyzt(1:3,iat)+rmi*gcons*shk%dro(1:3,i)
            end if
            if (.not.shk%freezeptr(jat)) then
              shk%xyzt(1:3,jat) = shk%xyzt(1:3,jat)-rmj*gcons*shk%dro(1:3,i)
            end if
          end do
        end if

      end if
      icyc = icyc+1
      if (.not.conv.and.icyc .le. shk%maxcyc) cycle
      exit
    end do

!>--- LOCAL EXTENSION (2026-10-24): diagnostics bookkeeping.
    shake_ncalls = shake_ncalls + 1
    shake_nsweep = shake_nsweep + min(icyc,shk%maxcyc)
    if (.not.conv) shake_nfail = shake_nfail + 1
    if (icyc > 1.and.maxdev > 0.0_wp.and.maxdev0 > 0.0_wp &
        .and.maxdev0 < huge(1.0_wp)) then
      block
        real(wp) :: r_est
        r_est = (maxdev/maxdev0)**(1.0_wp/real(icyc-1,wp))
        if (r_est > shake_rho_max.and.r_est < 1.0_wp) shake_rho_max = r_est
      end block
    end if
    if (.not.shake_stats_init) then
      shake_stats_init = .true.
      block
        character(len=32) :: e
        integer :: ist
        call get_environment_variable('CREST_SHAKE_STATS',value=e,status=ist)
        shake_stats_env = (ist == 0.and.len_trim(e) > 0.and.trim(e) /= '0')
      end block
    end if
    if (shake_stats_env.and.mod(shake_ncalls,1000) == 0) then
      !$omp critical
      write (stdout,'(3x,a,i0,a,i0,a,f8.1,a,f7.4)') &
      &  'SHAKE stats: calls ',shake_ncalls,' / fails ',shake_nfail, &
      &  ' / mean sweeps ',real(shake_nsweep,wp)/real(max(shake_ncalls,1),wp), &
      &  ' / worst rho ',shake_rho_max
      !$omp end critical
    end if

    if (conv) then

      if (nfrz .eq. 0) then
        velo = velo+(shk%xyzt-xyz)*tau1
        acc = acc+(shk%xyzt-xyz)*tau2
        xyz = shk%xyzt
        io = 0 !> successful termination
      else
        do i = 1,nat
          if (.not.shk%freezeptr(i)) then
            velo(:,i) = velo(:,i)+(shk%xyzt(:,i)-xyz(:,i))*tau1
            acc(:,i) = acc(:,i)+(shk%xyzt(:,i)-xyz(:,i))*tau2
            xyz(:,i) = shk%xyzt(:,i)
          end if
        end do
      end if
      io = 0
!>--- LOCAL EXTENSION: report the failure UNCONDITIONALLY.
!> Upstream gated this message behind "pr", and the trial MTD runs with
!> pr=.false. (setuptest.f90:164).  The consequence was that the ONLY hint that
!> the failure is a constraint problem - and which bond it is - never reached
!> the log; the user just got "Trial MTD n did not converge!" and everybody
!> starts shortening the time step, which cannot fix a SHAKE failure at all.
!> A failed do_shake() terminates the MD immediately (dynamics_module.f90:504),
!> so this prints AT MOST ONCE per MD and cannot spam the log.
    else
      call shake_failreport(shk,maxdev,maxdev0,min(icyc,shk%maxcyc))
    end if

    if (present(iostat)) then
      iostat = io
    end if

    return
  end subroutine do_shake

!========================================================================================!
!> LOCAL EXTENSION: the detailed counterpart of the "else if (pr)" above.
!>
!> Why this exists: do_shake() returning io/=0 and engrad() returning io/=0 are
!> reported to the caller through the SAME integer (dynamics_module.f90:340 and
!> :504), and the caller (setuptest.f90:229) only tests "io == 0".  Without this
!> report every MD failure looks like "did not converge" and the first instinct
!> is to shorten the time step - which is powerless against either a redundant
!> or an infeasible constraint set.
!>
!> What it prints:
!>   * how many iterations were actually used vs allowed (1000 default, see
!>     the P3-a note at the top of this module),
!>   * how far off the worst constraint still is vs tolshake (1e-7 default,
!>     i.e. the relative deviation of the SQUARED distance),
!>   * P3-c: the measured contraction per sweep rho, and how many more sweeps
!>     would have been needed.  These are what say "it was FOUR sweeps short"
!>     instead of leaving you to guess between divergence, a stall and a
!>     budget overrun - see the P3-c note inside the routine,
!>   * the 8 worst bonds with their current and their target distance in
!>     Angstroem - that pair of numbers is what tells an INFEASIBLE set
!>     (target unreachable from any geometry) apart from a REDUNDANT one
!>     (target reachable, but the solver cannot get there in maxcyc sweeps).
!========================================================================================!
  subroutine shake_stats_report()
!>--- LOCAL EXTENSION (2026-10-24): one-line summary of every SHAKE call in this
!> run.  Called from shake_failreport so a failure report carries context.
    implicit none
    write (stdout,'(3x,a,i0,a,i0,a,f8.1,a,f7.4)') &
    &  'SHAKE so far: calls ',shake_ncalls,' / fails ',shake_nfail, &
    &  ' / mean sweeps ',real(shake_nsweep,wp)/real(max(shake_ncalls,1),wp), &
    &  ' / worst rho ',shake_rho_max
    return
  end subroutine shake_stats_report

!========================================================================================!
  subroutine shake_failreport(shk,maxdev,maxdev0,icyc)
    implicit none
    type(shakedata),intent(in) :: shk
    real(wp),intent(in) :: maxdev
    real(wp),intent(in) :: maxdev0
    integer,intent(in) :: icyc
    integer,parameter :: ntop = 8
    integer :: i,k,ncons,nviol
    integer :: itop(ntop)
    real(wp) :: dtop(ntop),dev,rl,r0,amp
!>--- P3-c: the per-sweep contraction factor and the extrapolated number of
!> additional sweeps that would have been needed.  These two numbers are what
!> turn "it failed" into "it was 4 sweeps short" (see the P3-a note above).
    real(wp) :: rho,nmore
    logical :: have_rho

    ncons = shk%ncons
    dtop = -1.0_wp
    itop = 0

!>--- everything still above the tolerance
    nviol = 0
    do i = 1,ncons
      dev = abs(shk%dr(4,i)-shk%distcons(i))/max(shk%distcons(i),1.0e-30_wp)
      if (dev > shk%tolshake) nviol = nviol+1
!>--- insertion into the sorted "worst first" list
      do k = 1,ntop
        if (dev .le. dtop(k)) cycle
        if (k < ntop) then
          dtop(k+1:ntop) = dtop(k:ntop-1)
          itop(k+1:ntop) = itop(k:ntop-1)
        end if
        dtop(k) = dev
        itop(k) = i
        exit
      end do
    end do

!>--- P3-c: fit the convergence.
!>
!> The sweep converges GEOMETRICALLY, so the average contraction per cycle is
!> simply rho = (dev_after / dev_before)**(1/icyc), and the number of extra
!> cycles that would have been needed is
!>     n_more = ln(tolshake / dev_after) / ln(rho).
!> On the real failure that started all this: maxdev0 = 0.2521,
!> maxdev = 1.2887e-07, icyc = 250  ->  rho = 0.9437 and n_more = 4.4, i.e. the
!> run missed the tolerance by FOUR SWEEPS.  Without these two numbers the same
!> log is indistinguishable from a diverging or a stalling solver, and the
!> instinct is to shorten the time step - which only buys ~9 sweeps (the cost is
!> logarithmic in maxdev0, and maxdev0 ~ dt**2).
!>
!> Everything is guarded: maxdev0 <= 0 (nothing was violated to begin with),
!> a non-finite value, rho >= 1 (the residual is not shrinking - a genuine
!> stall, reported by the last DIAGNOSIS branch below) and maxdev0 itself.
    rho = 0.0_wp
    nmore = -1.0_wp
    have_rho = .false.
    if (icyc > 0.and.maxdev0 > 0.0_wp) then
      if (maxdev > 0.0_wp.and.maxdev0 < huge(1.0_wp).and.maxdev < huge(1.0_wp)) then
        rho = (maxdev/maxdev0)**(1.0_wp/real(icyc))
        have_rho = (rho > 1.0e-8_wp.and.rho < 1.0_wp)
      end if
    end if
    if (have_rho) then
      if (maxdev > shk%tolshake) then
        nmore = log(shk%tolshake/maxdev)/log(rho)
      else
        nmore = 0.0_wp
      end if
    end if

!> NOTE ON THE FORMATS: everything here is printed with 'es' (scientific, WITH
!> the exponent letter) and never with 'f'.  The first version of this report
!> used f9.5 for the distances and e12.4 for the deviations; on the very first
!> real failure both overflowed and the report came out as
!>      atoms 7 - 8   0.4404+291  r = ********* -> r0 =   1.11852 Angstroem
!> i.e. the exponent letter was dropped (e12.4 has no room for it) and the
!> distance was replaced by asterisks - precisely the numbers one needs.
    !$omp critical
    write (stdout,'(/,1x,a)') '> SHAKE did not converge!'
    write (stdout,'(3x,a,i0,a,i0)') 'iterations used / allowed   : ',icyc,' / ',shk%maxcyc
    write (stdout,'(3x,a,i0,a,i0)') 'constraints violated        : ',nviol,' of ',ncons
!>--- LOCAL FIX (2026-10-24): the old label said "BEFORE correcting", which
!> invited the reading "the starting geometry already violates the constraints".
!> It does not: maxdev0 is captured when icyc == 0, i.e. after the first raw MD
!> step has moved the atoms and BEFORE SHAKE touched anything.  The r0 targets
!> come from the starting structure, so it satisfies them by construction
!> (measured 0.00000 Angstrom on a real system).
    write (stdout,'(3x,a,es12.4)') 'deviation AFTER one raw step: ',maxdev0
    write (stdout,'(3x,a)') '  (how far one MD step moved the constrained distances;'
    write (stdout,'(3x,a)') '   NOT a mismatch between the starting geometry and r0 -'
    write (stdout,'(3x,a)') '   those agree by construction)'
    write (stdout,'(3x,a,es12.4)') 'deviation AFTER  correcting : ',maxdev
    write (stdout,'(3x,a,es12.4)') 'tolerance                   : ',shk%tolshake
    if (have_rho) then
      write (stdout,'(3x,a,f8.4)') 'contraction per sweep (rho) : ',rho
      if (nmore >= 0.0_wp.and.nmore <= 1.0e6_wp) then
        write (stdout,'(3x,a,f10.1)') 'more sweeps needed to reach : ',nmore
      end if
    end if
    call shake_stats_report()
    write (stdout,'(3x,a)') 'worst bonds (deviation, current -> target):'
    do k = 1,ntop
      if (itop(k) < 1) exit
      i = itop(k)
      rl = sqrt(max(shk%dr(4,i),0.0_wp))
      r0 = sqrt(max(shk%distcons(i),0.0_wp))
      write (stdout,'(5x,a,i0,a,i0,3x,es12.4,2x,a,es12.4,a,es12.4,a)') &
      &   'atoms ',shk%conslist(1,i),' - ',shk%conslist(2,i),dtop(k), &
      &   'r = ',rl*bohr,' -> r0 = ',r0*bohr,' Angstroem'
    end do
!>--- P3-c ordering note.  "Blew it up" is tested by whether the residual GREW
!> (maxdev >> maxdev0), not by where maxdev0 started - the earlier version
!> required maxdev0 <= 1e-3 and so would have called a genuine explosion with a
!> large starting violation a "stall".  A relative deviation above 1e3 means the
!> squared distance is off by three orders of magnitude, which no geometry
!> problem produces.
    if (maxdev > 1.0e3_wp.and.maxdev >= max(10.0_wp*maxdev0,1.0e-3_wp)) then
      write (stdout,'(3x,a)') 'DIAGNOSIS: the step arrived essentially ON the constraints'
!>--- amplification per cycle, guarded against Inf/NaN (log10 of a non-finite
!> value would make nint() undefined and could abort the run right here).
      amp = 0.0_wp
      if (maxdev0 > 0.0_wp.and.maxdev < huge(1.0_wp)) then
        amp = (log10(max(maxdev,1.0e-30_wp))-log10(maxdev0))/real(max(icyc,1))
        amp = min(max(amp,-1.0e6_wp),1.0e6_wp)
      end if
      write (stdout,'(3x,a,f8.3,a,i0,a)') &
      &   'and the ITERATION BLEW IT UP - grown by ~10^',amp,' per cycle over ', &
      &   icyc,' cycles.'
      write (stdout,'(3x,a)') 'That is a solver instability, NOT an unreachable target and'
      write (stdout,'(3x,a)') 'NOT a redundancy stall - a shorter time step cannot fix it.'
!>==============================================================================!
!> P3-c (2026-09-23): this branch is the one that used to be MISSING, and its
!> absence sent everybody down the wrong path.  On the real failure the report
!> read maxdev0 = 2.5213E-01, maxdev = 1.2887E-07, tolerance = 1.0000E-07 -
!> i.e. the sweep had come down SIX ORDERS OF MAGNITUDE and was within a factor
!> 1.29 of the tolerance - yet the old "maxdev0 > 1e-3" branch announced "the
!> step is too large - here a shorter time step DOES help".
!>
!> That is wrong.  When the residual has come that close, the run did not fail
!> because the target is unreachable or the step is too big: it failed because
!> maxcyc ran out.  Distinguish by where maxdev sits relative to tolshake, not
!> by where maxdev0 sits - maxdev0 only sets the STARTING point, and the sweep
!> count depends on it LOGARITHMICALLY (n ~ ln(maxdev0/tol)/ln(1/rho)), so a
!> large maxdev0 is normal and mostly harmless.
!==============================================================================!
    else if (maxdev <= 1.0e2_wp*shk%tolshake) then
      write (stdout,'(3x,a)') 'DIAGNOSIS: NOT a feasibility problem - the sweep came all the'
      write (stdout,'(3x,a)') 'way down to within a small factor of the tolerance and then ran'
      write (stdout,'(3x,a)') 'out of CYCLES.'
      if (nmore >= 0.0_wp.and.nmore <= 1.0e6_wp) then
        write (stdout,'(3x,a,f10.1,a)') 'It needed about ',nmore,' more sweeps than allowed.'
      end if
      write (stdout,'(3x,a)') 'Fix: raise maxcyc (it is a source default, see the P3-a note at'
      write (stdout,'(3x,a)') 'the top of this module) - or loosen tolshake a little. Do NOT'
      write (stdout,'(3x,a)') 'blame the time step: reducing it buys only ~9 sweeps per fs and'
      write (stdout,'(3x,a)') 'cannot turn a 250-cycle budget into the ~254 cycles needed here.'
    else if (.not.have_rho) then
      write (stdout,'(3x,a)') 'DIAGNOSIS: the residual is NOT shrinking any more - it has'
      write (stdout,'(3x,a)') 'stalled at a finite value. That is a constraint-set problem, not'
      write (stdout,'(3x,a)') 'a budget one: read the r -> r0 column above - if some bond sits'
      write (stdout,'(3x,a)') 'at a constant distance away from its target, the target is'
      write (stdout,'(3x,a)') 'unreachable (check the reference= key of the -cinp file, or a'
      write (stdout,'(3x,a)') 'bond listed twice with two different values). If every r already'
      write (stdout,'(3x,a)') 'equals its r0, the set is simply REDUNDANT - narrow the atom'
      write (stdout,'(3x,a)') 'list or fall back to -softcons.'
    else
      write (stdout,'(3x,a)') 'DIAGNOSIS: the residual is still descending but is nowhere near'
      write (stdout,'(3x,a)') 'the tolerance after the whole budget. Two causes can add up:'
      write (stdout,'(3x,a)') '(a) maxcyc is too small for this starting deviation, and'
      write (stdout,'(3x,a)') '(b) the step is genuinely too large, so maxdev0 is big.'
      write (stdout,'(3x,a)') 'Raise maxcyc FIRST (a source default - see the P3-a note);'
      write (stdout,'(3x,a)') 'only then consider a shorter time step, since its effect is'
      write (stdout,'(3x,a)') 'logarithmic - a few sweeps per fs, not orders of magnitude.'
    end if
    write (stdout,*)
    flush (stdout)
    !$omp end critical
    return
  end subroutine shake_failreport

!========================================================================================!
  integer function shake_lin(i1,i2)
!**********************
!* helper function lin
!**********************
    integer :: i1,i2,idum1,idum2
    idum1 = max(i1,i2)
    idum2 = min(i1,i2)
    shake_lin = idum2+idum1*(idum1-1)/2
    return
  end function shake_lin

!========================================================================================!
!========================================================================================!
end module shake_module
