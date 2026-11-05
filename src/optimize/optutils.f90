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

!> this module implements several utility routines for geometry optimization
!> including the threshold selection

module optimize_utils
  use crest_parameters
  use crest_calculator
  use axis_module
  use strucrd
  use ls_rmsd

  use optimize_type
  use optimize_maths
  use modelhessian_module
  use hessupdate_module
  implicit none
  public !> all public

  integer,parameter :: olev_crude = -3
  integer,parameter :: olev_sloppy = -2
  integer,parameter :: olev_loose = -1
  integer,parameter :: olev_lax = -4
  integer,parameter :: olev_normal = 0
  integer,parameter :: olev_tight = 1
  integer,parameter :: olev_vtight = 2
  integer,parameter :: olev_extreme = 3

!>--- LOCAL EXTENSION (2026-09-27): "tight" level used ONLY by the conformer-search
!> multilevel chain ("optimization with tight thresholds" in set_multilevel_options,
!> = search_conformers.f90).
!>
!> WHY A SEPARATE LEVEL: that stage only has to polish MTD snapshots so that
!> CREGEN can deduplicate them and so that the geometry is a sane starting point
!> for a later, better-level re-optimisation.  It is NOT the final answer, so
!> olev_tight (1e-6 Eh / 8e-4 Eh/a0) over-converges for it and burns cycles.
!>
!> WHY NOT loosen olev_tight ITSELF: the ladder must stay monotonic, and
!> olev_tight is also the target of the user-facing "tight" optlevel keyword
!> (classes.f90 optlevnum) and of -tight in other run types (MECP/TSCoord work
!> does care about tight geometry).  Loosening it in place would make
!> "tight" LOOSER than olev_normal (5e-6 / 1e-3) in the energy criterion, i.e.
!> it would invert the meaning of the ladder.  Hence a separate, additive level.
!>
!> VALUES: 1e-4 Eh / 4e-3 Eh/a0, i.e. 100x / 5x looser than olev_tight (tuned by
!> the maintainer on 2026-09-27; see the note on case (olev_tightcs) in
!> get_optthr below for the full history and for the consequences).
!> 1e-4 Eh = 0.0628 kcal/mol - two orders of magnitude below the error of the
!> MTD-level method itself, so invisible for conformer ranking on the ENERGY
!> side.  The GRADIENT criterion is the operative one in practice, so it (not
!> this number) is what controls the cost of the stage.
  integer,parameter :: olev_tightcs = 4

!>--- LOCAL EXTENSION (2026-09-29): "very tight" companion of olev_tightcs,
!> used by the LAST stage (6) of the conformer-search multilevel chain.
!>
!> WHY: stage 6 must sit exactly one factor of two below stage 5.  Pointing it
!> at the absolute olev_vtight level (1e-7 / 2e-4) would make the jump from
!> stage 5 (1e-4 / 4e-3) a factor of 1000 in energy and 20 in gradient, which
!> is not a ladder - it is a wall, and in practice stage 6 then exhausts
!> maxcycle on most structures and reports them as failures.
!>
!> VALUES: half of olev_tightcs in BOTH criteria:
!>   ethr = 5e-5 Eh      (olev_tightcs: 1e-4)   = 0.0314 kcal/mol
!>   gthr = 2e-3 Eh/a0   (olev_tightcs: 4e-3)
!> The energy criterion is shared with olev_loose in value only; the pair is
!> what identifies a level, and this pair is unique in the table.
!>
!> NOTE: this is a SEPARATE level from olev_vtight on purpose.  olev_vtight is
!> the target of the user-facing "vtight" optlevel keyword (classes.f90
!> optlevnum) and must keep its absolute meaning; stage 6 is a conformer-search
!> internal stage whose only job is to polish the final ensemble, so it is
!> defined relative to its predecessor instead.
  integer,parameter :: olev_vtightcs = 5

!>--- LOCAL EXTENSION (2026-09-27): the GLOBAL optimisation cycle budget.
!>
!> WHY: the per-level `maxcycle = n*K` values in get_optthr below were never
!> really in charge.  print_opt_data() calls get_optthr(0,...), i.e. n = 0, and
!> the old clamp `max(maxcycle,0*K,200)` turned that into a flat 200 which was
!> then written into calc%maxcycle and kept for every following stage (the
!> assignment is guarded by `calc%maxcycle <= 0`).  So the multilevel chain
!> always ran with 200 cycles regardless of the level, while a standalone --opt
!> run (n = nat) got nat*K, i.e. something entirely different.  Rather than
!> leave that coincidence in place, the budget is now an explicit constant.
!>
!> 120 cycles: ANCOPT reaches the (deliberately loose) criteria of these stages
!> in a few dozen cycles for well-behaved systems, so 120 is still generous as
!> a cap while cutting the tail.  A structure that exhausts it is no longer
!> discarded - see the anopt note in calculator/calc_type.f90.
!>
!> An explicit user value still wins: calc%maxcycle is only overwritten when it
!> is <= 0, so a TOML `maxcycle` (parse_calcdata.f90:541) or a deliberate
!> setting such as protonate.f90's 10000 is respected.
  integer,parameter :: opt_maxcycle_global = 120

  public :: get_optthr

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine get_optthr(n,olev,calc,ethr,gthr)
!************************************************
!* subroutine get_optthr
!* routine to set some optimization thresholds
!************************************************
    implicit none
    integer,intent(in) :: n
    integer,intent(inout) :: olev
    type(calcdata) :: calc
    real(wp),intent(out) :: ethr
    real(wp),intent(out) :: gthr
    integer :: maxcycle
!>--- NOTE (found 2026-09-27), two latent quirks in this routine:
!>  (a) `acc` is assigned in every case below but NEVER read - the optimiser
!>      takes its accuracy from calc%acc_opt instead (ancopt.f90:367
!>      `acc_in = calc%acc_opt`), a separate field with its own default (1.0).
!>      The per-level acc tuning that this table implies has therefore never
!>      taken effect.  The assignments are kept: removing them would hide the
!>      intent, and the values are right if anybody ever wires acc -> acc_opt.
!>  (b) HISTORICAL, now superseded (2026-09-27): print_opt_data() (called by
!>      set_multilevel_options BEFORE the batch loop) calls get_optthr(0,...),
!>      i.e. n = 0.  With the old clamps that made maxcycle collapse to the
!>      floor of 200 and - because it is only written when calc%maxcycle <= 0 -
!>      it STAYED 200 for every following stage, so the per-level n*3 / n*5
!>      budgets never governed the multilevel chain.  That coincidence is now
!>      replaced by the explicit global budget opt_maxcycle_global (120), so the
!>      flat value is deliberate.  The n*K lines in the cases below are kept as
!>      documentation of CREST's original per-level intent only.
    real(wp) :: acc
    select case (olev)
!> very approximate = crude
    case (olev_crude)
      ethr = 5.d-4
      gthr = 1.d-2
      maxcycle = n
      acc = 3.00d0
!> approximate = sloopy
    case (olev_sloppy)
      ethr = 1.d-4
      gthr = 6.d-3
      maxcycle = n
      acc = 3.00d0
!> loose
    case (olev_loose)
      ethr = 5.d-5
      gthr = 4.d-3
      maxcycle = n*2
      acc = 2.00d0
!>  for DCOSMO-RS opts with TM i.e. between loose and normal, keyword "lax"
    case (olev_lax)
      ethr = 2.d-5
      gthr = 2.5d-3
      maxcycle = n*2
      acc = 2.00d0
!> normal
    case default
      ethr = 5.d-6
      gthr = 1.d-3
      maxcycle = n*3
      acc = 1.0d0
!> tight
    case (olev_tight)
      ethr = 1.d-6
      gthr = 8.d-4
      maxcycle = n*5
      acc = 0.20d0
!> very tight
    case (olev_vtight)
      ethr = 1.d-7
      gthr = 2.d-4
      maxcycle = n*20
      acc = 0.05d0
!> extreme
    case (olev_extreme)
      ethr = 5.d-8
      gthr = 5.d-5
      maxcycle = n*20
      acc = 0.01d0
!> conformer-search "tight" (olev_tightcs): deliberately loose.
!> Maintainer tuning history (2026-09-27), starting from olev_tight:
!>   1e-6 / 8e-4 -> 1e-5 / 8e-3 -> 5e-5 / 1.6e-2 -> 5e-5 / 4e-3 -> 1e-4 / 4e-3
!>   (current: 100x looser than olev_tight in ethr, 5x in gthr)
!> Kept as its own level on purpose: tuning this stage must not move the
!> user-facing "loose" optlevel, and olev_tight has to keep its meaning for
!> -tight / MECP-style work.  (Note ethr = 1e-4 now coincides with olev_sloppy's
!>  value; the two criteria are nevertheless not the same pair, since sloppy
!>  pairs it with the tighter gthr 6e-3.)
!> Energy equivalent of the ENERGY criterion: 1e-4 Eh = 0.0628 kcal/mol.
!> The GRADIENT criterion is what actually terminates the stage (ancopt.f90:585
!> `converged = econverged.and.gconverged.and.lowered`, and |dE| drops below
!> ethr long before the gradient drops below gthr), so loosening ethr alone
!> buys very little - for soft-mode systems, whose |dE| lingers above the old
!> 5e-5, it can shave a few cycles; otherwise the cycle count is unchanged.
!> gthr is the knob that controls how much this stage costs.
!> 4e-3 Eh/a0 = 4.74 kcal/mol/A as a vector NORM (0.36 kcal/mol/A RMS per
!> coordinate at N=60), still TIGHTER than this table's olev_crude (1e-2), i.e.
!> the stage stays on the tight side of the crude PRE-optimisation (stage 1).
!> Residual geometry quality, worst case (whole gradient norm on one soft mode):
!>   lambda  0.5 Eh/a0^2 -> 0.010 kcal/mol, dx 0.004 A
!>   lambda  0.1          -> 0.050 kcal/mol, dx 0.021 A
!>   lambda  0.05         -> 0.100 kcal/mol, dx 0.042 A
!> i.e. comfortably inside CREGEN's dedup thresholds (rthr 0.125 A) and far
!> below the error of the MTD-level method itself.  Adding the energy-criterion
!> slop (0.063 kcal/mol) changes none of that.
!> maxcycle is now the global opt_maxcycle_global (120) - the old note that it
!> was pinned to the print_opt_data() clamp of 200 is obsolete: that clamp is
!> gone, and the flat value is now deliberate rather than accidental.  The
!> `acc` quirk below is unaffected.
!> maxcycle/acc no longer follow the olev_normal values here: the budget is
!> global, and acc is inert in any case (see the note above).
    case (olev_tightcs)
      ethr = 1.d-4  !> olev_tightcs
      gthr = 4.d-3  !> olev_tightcs
      maxcycle = n*3
      acc = 1.00d0
!> conformer-search "very tight" (olev_vtightcs): DEFINED as half of
!> olev_tightcs, so that stage 6 of the multilevel chain is always exactly one
!> factor of two below stage 5.  See the long note on the parameter above.
!>   ethr 1e-4 -> 5e-5 Eh      (0.0314 kcal/mol)
!>   gthr 4e-3 -> 2e-3 Eh/a0
!> Kept as its own level rather than reusing olev_vtight (1e-7 / 2e-4), which
!> is the user-facing "vtight" target and must not move.
    case (olev_vtightcs)
      ethr = 5.d-5  !> olev_vtightcs = olev_tightcs / 2
      gthr = 2.d-3  !> olev_vtightcs = olev_tightcs / 2
      maxcycle = n*5
      acc = 0.50d0
    end select
!>--- LOCAL EXTENSION (2026-09-27): one explicit GLOBAL budget replaces the old
!> per-level values and their clamps (`min(maxcycle,10000)` / `max(maxcycle,200)`,
!> the latter of which is what really pinned every multilevel stage to 200 - see
!> the note on opt_maxcycle_global above).  The `maxcycle = n*K` lines in the
!> cases above are now DOCUMENTATION ONLY: they record what CREST intended per
!> level, but they no longer decide anything.  Keep them in sync if the intent
!> ever becomes relevant again, and do not "fix" the resulting warning-free
!> dead store - it is deliberate.
    maxcycle = opt_maxcycle_global
!> user-defined overwrites
    if (calc%maxcycle <= 0) then
      calc%maxcycle = maxcycle
    end if
    if (calc%ethr_opt > 0.0_wp) then
      ethr = calc%ethr_opt
      olev = 99
    end if
    if (calc%gthr_opt > 0.0_wp) then
      gthr = calc%gthr_opt
      olev = 99
    end if
    if (calc%tsopt) then
      calc%hlow_opt = max(calc%hlow_opt,0.250d0)
      calc%iupdat = 1
    end if

    return
  end subroutine get_optthr

!========================================================================================!

  subroutine prdechng(nat3,grad,displ,hess,depred)
!********************************************************************
!* Purpose:
!* Calculates predicted energy change according to the second order
!* model.
!*
!* Input:
!* nat3  - 3*natoms
!* hess  - Hessian matrix stored as lower triangle
!* grad  - Gradient vector
!* displ - Displacement vector
!*
!* Output:
!* depred - Predicted energy change
!******************************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat3
    real(wp),intent(in) :: grad(nat3)
    real(wp),intent(in) :: displ(nat3)
    real(wp),intent(in) :: hess(nat3*(nat3+1)/2)
    !> OUTPUT
    real(wp),intent(out) :: depred
    !> LOCAL
    real(wp),allocatable :: hdx(:)
    real(wp) :: gtmp,htmp
    !> BLAS functions:
    real(wp),external :: ddot
    external :: dspmv
    allocate (hdx(nat3),source=0.0_wp)
    call dspmv('u',nat3,0.5d0,hess,displ,1,0.0d0,hdx,1)
    gtmp = ddot(nat3,displ,1,grad,1)
    htmp = ddot(nat3,displ,1,hdx,1)
    depred = htmp+gtmp
    return
  end subroutine prdechng

  subroutine prdispl(nvar,displ)
    implicit none
    integer,intent(in) :: nvar
    real(wp),intent(in) :: displ(nvar)
    real(wp),allocatable :: er(:)
    integer,allocatable :: merk(:)
    integer :: i,j,ii,k
    integer :: ihilf
    real(wp) :: pp
    allocate (er(nvar),source=0.0_wp)
    allocate (merk(nvar),source=0)

    er = abs(displ)

    do i = 1,nvar
      merk(i) = i
    end do
    do ii = 2,nvar
      i = ii-1
      k = i
      pp = er(i)
      do j = ii,nvar
        if (er(j) .le. pp) cycle
        k = j
        pp = er(j)
      end do
      if (k .eq. i) cycle
      er(k) = er(i)
      er(i) = pp
      ihilf = merk(i)
      merk(i) = merk(k)
      merk(k) = ihilf
    end do

    write (*,'(''Largest |displ|/coords:'',5(f8.4,'' ('',i4,'')''))') &
      (er(i),merk(i),i=1,min(3,nvar))

  end subroutine prdispl

!========================================================================================!
  function bool2string(bool)
    implicit none
    character(len=:),allocatable :: bool2string
    logical :: bool
    if (bool) then
      bool2string = "True"
    else
      bool2string = "False"
    end if
  end function bool2string

!========================================================================================!
  subroutine geoconvav(nc,e,g,val,deriv)
    implicit none
    integer :: nc     !> total number of E/G points
    real(wp) :: e(*)  !> total energy in Eh
    real(wp) :: g(*)  !> norm of Cartesian gradient (in TM: |dE/dxyz|)
    real(wp) :: val   !> av. energy in Eh to be used further
    real(wp) :: deriv !> av. gradient

    integer :: low
    integer :: i,j
    integer,parameter:: nav = 5 !> average over last nav
    real(wp) :: eav,gav

    !> only apply it if sufficient number of points i.e. a "tail" can exist
    !> with the censo blockl = 8 default, this can first be effective in the second
    if (nc .lt. 3*nav) then
      val = e(nc)
      deriv = g(nc)
      return
    end if

    low = max(1,nc-nav+1)
    j = 0
    eav = 0
    do i = nc,low,-1
      j = j+1
      eav = eav+e(i)
      gav = gav+g(i)
    end do
    val = eav/float(j)

    low = max(1,nc-nav+1)
    j = 0
    gav = 0
    do i = nc,low,-1
      j = j+1
      gav = gav+g(i)
    end do
    ! adjust the gradient norm to xtb "conventions" because e.g. a noisy
    ! DCOSMO-RS gradient for large cases can never (even on average)
    ! become lower than the "-opt normal" thresholds
    deriv = gav/float(j)/2.d0
  end subroutine geoconvav

!========================================================================================!
  subroutine rdhess(nat3,h,fname)
    integer,intent(in)  :: nat3
    real(wp),intent(out) :: h(nat3,nat3)
    character(len=*),intent(in) :: fname
    integer  :: iunit,i,j,mincol,maxcol
    character(len=5)  :: adum
    character(len=80) :: a80

    open (newunit=iunit,file=fname)
50  read (iunit,'(a)') a80
    if (index(a80,'$hessian') .ne. 0) then
      do i = 1,nat3
        maxcol = 0
200     mincol = maxcol+1
        maxcol = min(maxcol+5,nat3)
        read (iunit,*) (h(j,i),j=mincol,maxcol)
        if (maxcol .lt. nat3) goto 200
      end do
      close (iunit)
      goto 300
    end if
    goto 50

300 return
  end subroutine rdhess

!========================================================================================!

  subroutine print_optiter(iter)
    implicit none
    integer :: iter
    write (*,'(/,"┌",76("─"),"┐")')
    write (*,'(  "│",32(" ")," CYCLE",i5,1x,32(" "),"│")') iter
    write (*,'(  "└",76("─"),"┘")')
  end subroutine print_optiter

!========================================================================================!

  subroutine print_optsummary(calc,tight,nvar,maxcycle,maxmicro, &
      &                       ethr,gthr,linear,wr)
    implicit none
    !> INPUT
    type(calcdata) :: calc
    integer,intent(in) :: tight,maxcycle,maxmicro,nvar
    real(wp),intent(in) :: ethr,gthr
    logical,intent(in) :: linear,wr
    !> LOCAL
    integer :: iupdat,engine
    real(wp) :: hlow,hmax,maxdispl,s6
    type(mhparam) :: mhset
    !> FORMATS
    character(len=*),parameter :: scifmt = &
                                  '(10x,"│",3x,a,e21.7,1x,a,1x,"│")'
    character(len=*),parameter :: dblfmt = &
                                  '(10x,"│",3x,a,f18.7,5x,a,1x,"│")'
    character(len=*),parameter :: intfmt = &
                                  '(10x,"│",3x,a,i18,      10x,"│")'
    character(len=*),parameter :: chrfmt = &
                                  '(10x,"│",3x,a,a18,      10x,"│")'

!>--- set params
    engine = calc%opt_engine
    iupdat = calc%iupdat
    hlow = calc%hlow_opt !> 0.01 in ancopt, 0.002 too small
    hmax = calc%hmax_opt
    maxdispl = calc%maxdispl_opt
    s6 = mhset%s6 !> slightly better than 30 for various proteins

!>--- print a summary of settings, if desired
    write (*,'(/,10x,"┍",49("━"),"┑")')
    write (*,'(10x,"│",11x,a,11x,"│")') "GEOMETRY OPTIMIZATION SETUP"
    write (*,'(10x,"┝",49("━"),"┥")')
    select case (engine)
    case (0)
      write (*,chrfmt) "algorithm         ","          ANCOPT"
    case (1)
      write (*,chrfmt) "algorithm         ","          L-BFGS"
    case (2)
      write (*,chrfmt) "algorithm         ","Rational Function"
    case (-1)
      write (*,chrfmt) "algorithm         ","Gradient Descent"
    end select

    select case (tight)
    case (-3)
      write (*,chrfmt) "optimization level",'crude (-3)'
    case (-2)
      write (*,chrfmt) "optimization level",'vloose (-2)'
    case (-1)
      write (*,chrfmt) "optimization level",'loose (-1)'
    case (0)
      write (*,chrfmt) "optimization level",'normal (0)'
    case (1)
      write (*,chrfmt) "optimization level",'tight (1)'
    case (2)
      write (*,chrfmt) "optimization level",'vtight (2)'
    case (3)
      write (*,chrfmt) "optimization level",'extreme (3)'
!>--- LOCAL EXTENSION (2026-09-29): the two conformer-search-only levels.
!> Without these cases the log would print "custom" for stage 5/6, which is
!> exactly the kind of nameless level that makes a run hard to read back.
    case (olev_tightcs)
      write (*,chrfmt) "optimization level",'tightcs (4)'
    case (olev_vtightcs)
      write (*,chrfmt) "optimization level",'vtightcs (5)'
    case default
      write (*,chrfmt) "optimization level",'custom'
    end select
    write (*,intfmt) "max. optcycles    ",maxcycle

    select case (engine)
    case (0)
      write (*,intfmt) "ANC micro-cycles  ",maxmicro
    end select

    write (*,intfmt) "degrees of freedom",nvar
    !if (modef > 0) then
    !  write (*,intfmt) "# mode follow     ",modef
    !end if
    write (*,'(10x,"├",49("─"),"┤")')
    select case (engine)
    case (0,2)
      if (calc%exact_rf) then
        write (*,chrfmt) "RF solver         ","spevx"
      else
        write (*,chrfmt) "RF solver         ","davidson"
      end if
    end select
    if (engine .ne. -1) then
      select case (iupdat)
      case (0)
        write (*,chrfmt) "Hessian update    ","bfgs"
      case (1)
        write (*,chrfmt) "Hessian update    ","powell"
      case (2)
        write (*,chrfmt) "Hessian update    ","sr1"
      case (3)
        write (*,chrfmt) "Hessian update    ","bofill"
      case (4)
        write (*,chrfmt) "Hessian update    ","schlegel"
      end select
    end if
    write (*,chrfmt) "write crestopt.log",bool2string(wr)
    if (linear) then
      write (*,chrfmt) "linear (good luck)",bool2string(linear)
    else
      write (*,chrfmt) "linear?           ",bool2string(linear)
    end if
    write (*,scifmt) "energy convergence",ethr,"Eh   "

    select case (engine)
    case (0)
      write (*,scifmt) "grad. convergence ",gthr,"Eh/α "
    case default
      write (*,scifmt) "grad. convergence ",gthr,"Eh/a0"
    end select

    select case (engine)
    case (0)
      write (*,dblfmt) "maximium RF displ.",maxdispl,"    "
      write (*,scifmt) "Hlow (freq-cutoff)",hlow,"     "
      write (*,dblfmt) "Hmax (freq-cutoff)",hmax,"    "
      write (*,dblfmt) "S6 in model hess. ",s6,"    "
    case (2)
      write (*,dblfmt) "maximium RF displ.",maxdispl,"    "
      write (*,'(10x,"│",3x,a,e21.7,1x,a,"│")') &
      & "Hessian guess     ",calc%hguess,"Eh/a0²"
    case (-1)
      write (*,dblfmt) "max step size     ",maxdispl,"    "
    end select
    write (*,'(10x,"└",49("─"),"┘")')

  end subroutine print_optsummary

!========================================================================================!

  subroutine print_convd(econverged,gconverged)
    logical,intent(in) :: econverged,gconverged
    write (*,'(3x,"converged δE/grad :",1x,a5," /",1x,a)') &
    & bool2string(econverged),bool2string(gconverged)
  end subroutine print_convd

!========================================================================================!
!========================================================================================!
end module optimize_utils
