!================================================================================!
! This file is part of CREST.
!
! Copyright (C) 2023 Philipp Pracht
!
! CREST is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CREST is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!> This module implements a simple gradient descend algorithm

module gradientdescent_module
  use crest_parameters
  use crest_calculator
  use axis_module
  use strucrd
  use ls_rmsd

  use optimize_type
  use optimize_maths
  use modelhessian_module
  use hessupdate_module
  use optimize_utils
  implicit none
  private

  public :: gradientdescent

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine gradientdescent(mol,calc,etot,grd,pr,wr,iostatus)
!*************************************************************************
!> subroutine gradientdescent
!> Implementation of a simple gradient descent algorithm
!>
!> Input/Output:
!>      mol  - object containing the molecule,
!>             Cartesian coordinates in Bohrs.
!>             will be overwritten on output
!>     calc  - object containing calculation settings
!>             and optimization thresholds (look for calc% )
!>     etot  - on input initial energy (do a singlepoint before ancopt)
!>             on output final energy
!>      grd  - Cartesian gradient
!>       pr  - printout bool
!>       wr  - logfile (crestopt.log) bool
!>  iostatus - return status of the routine
!>             (success=0, error<0, not converged>0)
!!***********************************************************************
    implicit none
    !> INPUT/OUTPUT
    type(coord),intent(inout) :: mol
    type(calcdata),intent(in) :: calc
    real(wp),intent(inout) :: etot
    real(wp),intent(inout) :: grd(3,mol%nat)
    logical,intent(in) :: pr
    logical,intent(in) :: wr
    integer,intent(out) :: iostatus
    !> LOCAL
    integer  :: tight
    real(wp) :: eel
    real(wp) :: et
    real(wp) :: egap
    logical :: fail
    !> Local objects
    type(coord)   :: molopt
    type(optimizer)  :: anc
    type(mhparam) :: mhset

    real(wp) :: step,amu2au,au2cm,dumi,dumj,damp,hlow,edum,s6,thr
    real(wp) :: maxdispl,gthr,ethr,hmax,energy,rij(3),t1,t0,w1,w0
    real(wp) :: rot(3),maxerise
    integer :: n3,i,j,k,l,jjj,ic,jc,ia,ja,ii,jj,info,nat3
    integer :: nvar,iter,nread,maxcycle,maxmicro,itry,maxopt,iupdat,iii
    integer :: id,ihess,error
    integer :: ilog,nopt
    real(wp) :: lambda,gnorm,dnorm,ddot,eold,xdum,estart,acc,e_in
    real(wp) :: depred,echng,dummy,maxd,alp,gchng,gnold
!>--- LOCAL FIX (2026-10-02), see the notes at the ERISE loop below:
!>  g0norm  = |g| at the point the current iteration STARTS from, i.e. the norm
!>            of the frozen search direction (not of the current trial point).
!>  steptol = smallest displacement that is still worth an energy evaluation
!>            (bohr).  Below it the trial geometry is indistinguishable from the
!>            current one within the precision of the level of theory.
    real(wp) :: g0norm
    real(wp),parameter :: steptol = 1.0e-10_wp

    real(wp),allocatable :: displ(:,:),grdold(:,:),ddispl(:,:)
!>--- LOCAL FIX (2026-10-02): grd0 = frozen gradient of the current point (the
!> search direction of the whole iteration), grdt = scratch buffer for the
!> gradient of a TRIAL point, so that a rejected trial can no longer overwrite
!> the search direction.
    real(wp),allocatable :: grd0(:,:),grdt(:,:)
    real(wp),allocatable :: grmsd(:,:),gdiff(:,:)
    type(convergence_log),allocatable :: avconv
    real(wp) :: U(3,3),x_center(3),y_center(3),rmsdval
    integer :: modef
    logical :: ex,converged,linear,econverged,gconverged,lowered
    real(wp) :: esave

    iostatus = 0
    if (mol%nat .eq. 1) return
!>  defaults
    nopt = mol%nat*3
    tight = calc%optlev
    call get_optthr(mol%nat,tight,calc,ethr,gthr)
!>--- LOCAL FIX (2026-10-02): calc%maxdispl_opt defaults to 1.0_wp
!> (calc_type.f90:254), and that default was translated here into
!> real(mol%nat) - i.e. the trust radius became the ATOM COUNT (87 bohr for a
!> typical CREST system, against the 1.0 bohr that ancopt.f90:366 gives to the
!> very same field).  Keep the value: it is a TRUST RADIUS IN BOHR and is now
!> applied to the resulting displacement below, with exactly the meaning it has
!> in ANCOPT.  Side effect: the "max step size" line of print_optsummary
!> (optutils.f90:586) prints calc%maxdispl_opt and used to be a lie for GD; it
!> is the truth now.
    if(calc%maxdispl_opt .ne. 1.0_wp)then
      maxdispl = calc%maxdispl_opt
    else
      maxdispl = 1.0_wp
    endif

!> max number of iterations
    maxmicro = 100
    maxcycle = calc%maxcycle
    if (maxcycle .lt. maxmicro) maxmicro = maxcycle

    !> check if the molecule is linear
    call axis(mol%nat,mol%at,mol%xyz,rot,dumi)
    linear = (rot(3) .lt. 1.d-10).or.(mol%nat == 2)

    !> set degrees of freedom
    nat3 = 3*mol%nat
    nvar = nat3-6
    if (linear) then
      nvar = nat3-5
    end if
    if (calc%nfreeze .gt. 0) then ! exact fixing
      nvar = nat3-3*calc%nfreeze-3
      if (nvar .le. 0) nvar = 1
    end if

    !$omp critical
    allocate (grmsd(3,mol%nat))
    !$omp end critical

!>--- print a summary of settings, if desired
    if (pr) then
      call print_optsummary(calc,tight,nvar,maxcycle,maxmicro, &
      &                       ethr,gthr,linear,wr)
    end if


!>--- initialize storage
    !$omp critical
    allocate (molopt%at(mol%nat),molopt%xyz(3,mol%nat))
    allocate (displ(3,mol%nat),ddispl(3,mol%nat))
    allocate (grdold(3,mol%nat),gdiff(3,mol%nat))
!>--- LOCAL FIX (2026-10-02)
    allocate (grd0(3,mol%nat),grdt(3,mol%nat))
    !$omp end critical

!>--- backup coordinates, and starting energy
    molopt%nat = mol%nat
    molopt%at = mol%at
    molopt%xyz = mol%xyz
    estart = etot
    gnorm = 0.0_wp
    depred = 0.0_wp
    echng = 0.0_wp
    alp = 1.0_wp
    ddispl = 0.0_wp
    energy = etot

!>--- initialize .log file, if desired
    ilog = 942
    if (wr) then
      open (newunit=ilog,file='crestopt.xyz')
    end if

!>--- The gradient descent iteration loop. "iter"
    iter = 0
    fail = .false.
    converged = .false.

!>--- step 0 (just printout from preceeding SP)
    if (pr) call print_optiter(iter)
    gnorm = norm2(grd)
    if (iostatus .ne. 0) then
      fail = .true.
    end if
    grdold = grd
    gnold = gnorm
    eold = energy
    if (pr) then
      write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') energy
      write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') 0.0_wp
      write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
      write (*,'(2x,"step size",e18.7)',advance='yes') alp
    end if

!>======================================================================
    if (.not.fail) then
      SDiter: do while (iter < maxcycle.and..not.converged)
!>======================================================================
        iter = iter+1
        if (pr) call print_optiter(iter)
        grdold = grd

!>--- LOCAL FIX (2026-10-02): freeze the search direction for THIS iteration.
!> The original code stepped along -alp*grd while grd was overwritten by every
!> single trial (see below).  After the first rejected trial the loop therefore
!> walked along the gradient of a point that had just been REJECTED instead of
!> the gradient of the point it is trying to improve - that is not a line
!> search, the direction changes with every halving.  For a large first step
!> (see the maxdispl note above) the trial point can even sit in a region whose
!> gradient is anti-parallel to the one here, and then NO step size along it
!> can lower the energy, which is exactly how the "FAILED to perform gradient
!> descent" branch was reached.
        grd0 = grd
        g0norm = norm2(grd0)

!>--- check if the energy decreased. If not, make alp smaller and try again
!===================================!
        ERISE: do jjj = 1,maxmicro
!===================================!

!>--- determine a good step size
          if (iter == 1) then
            alp = 1.0_wp !> first step just from gradient
          else
            alp = bbstep(nopt,ddispl,gdiff)
          end if
!>--- LOCAL FIX (2026-10-02): cap the DISPLACEMENT, not alp.  alp multiplies
!> the gradient, so `alp = min(maxdispl,alp)` scales every component of the
!> gradient by the same factor - with the old default maxdispl = real(mol%nat)
!> one scalar of ~87 was enough to blow a stiff component up to several bohr.
!> Dividing by |g0| turns the cap into "no step longer than maxdispl bohr",
!> which keeps the Barzilai-Borwein step intact where the gradient is small
!> (soft modes, the regime BB is actually made for) and is the same quantity
!> that ancopt.f90 applies its trust radius to.
          alp = min(alp,maxdispl/max(g0norm,1.0e-30_wp))
          alp = alp*(0.5_wp**(jjj-1))

!>--- LOCAL FIX (2026-10-02): step-size floor.  Once alp*|g0| falls below the
!> resolution of the geometry, the trial point is - within the precision of the
!> level of theory - the SAME structure, so the energy comes back bitwise
!> identical and echng is EXACTLY 0.0, which `lowered` (a strict "< 0") can
!> never accept.  Every remaining halving is then a guaranteed miss AND a
!> wasted energy+gradient evaluation (~60 of the 100 per failed structure, each
!> one a full DFT gradient when ORCA is the calculator).  Stop the line search
!> here: a step that cannot be represented means this engine has reached a
!> stationary point, so it is converged if the gradient criterion is met, and a
!> plain non-convergence otherwise (handled by the fall-through below).
          if (alp*g0norm .lt. steptol) then
            if (pr) write (*,'(" * step size below geometry resolution. Stopping the line search ...")')
            energy = eold
            grd = grd0
            gnorm = g0norm
            echng = 0.0_wp
            econverged = .true.
            gconverged = gnorm .lt. gthr
            converged = econverged.and.gconverged
            if (converged) etot = energy
            exit SDiter
          end if

!>--- take the gradient descent step
          displ(:,:) = -alp*grd0(:,:)
          molopt%xyz = molopt%xyz+displ

!>--- calculate new energy and gradient
!>--- LOCAL FIX (2026-10-02): the gradient of the TRIAL point goes into a
!> scratch buffer.  grd keeps describing the current point until the step is
!> accepted, so a rejected trial can no longer corrupt the search direction.
          grdt = 0.0_wp
          call engrad(molopt,calc,energy,grdt,iostatus)
          if (iostatus .ne. 0) then
            fail = .true.
            exit SDiter
          end if

!>--- optimization step analysis and printout
          gnorm = norm2(grdt)
          echng = energy-eold
          econverged = abs(echng) .lt. ethr
          gconverged = gnorm .lt. gthr
          lowered = echng .lt. 0.0_wp
          converged = econverged.and.gconverged
          if (pr) then
            write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') energy
            write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') echng
            write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
            write (*,'(2x,"step size",e18.7)',advance='yes') alp
          end if

!>--- for too large energy rise, rewind step. otherwise go to next iteration
!>--- LOCAL FIX (2026-10-02): ".and. .not.converged".  The original test looked
!> at `lowered` only, so a structure that already satisfied BOTH convergence
!> criteria at the trial point (|ΔE| within ethr AND |g| below gthr - the
!> normal situation for a structure handed over from a previous optimisation
!> stage, and unavoidable at a flat minimum with a noisy energy) was rewound
!> anyway and ended up in the jjj == maxmicro branch below, which forces
!> converged = .false. and throws a perfectly good structure away.  Compare
!> ancopt.f90:585, where `lowered` is part of the convergence test instead of
!> being a veto over it.
          if (.not.lowered.and..not.converged) then
            if(pr) write (*,'(" * ENERGY RISE detected. Rewind and try again with smaller step ...")')
            molopt%xyz = molopt%xyz-displ
            !write (*,*) jjj,maxmicro
            if (jjj == maxmicro) then
              write (*,'(" * FAILED to perform gradient descent! Something is very wrong ...")')
              converged = .false.
              exit SDiter
            end if
          else
!>--- accept the step: from here on grd describes the NEW geometry again
            grd = grdt
            exit ERISE
          end if
!===========================!
        end do ERISE !> jjj
!===========================!

!>--- save some things for next iteration
        ddispl = -displ
        gdiff = grd-grdold
        eold = energy
        gnold = gnorm

!>--- dump coordinates after step
        if (wr) then
          call molopt%appendlog(ilog,energy)
        end if

!>--- check convergence
        if (pr) then
          call print_convd(econverged,gconverged)
        end if
        if (converged) then
          converged = .true.
          etot = energy
          exit SDiter
        end if

!>======================================================================
      end do SDiter
    end if
!>======================================================================

!>--- close .log file
    if (wr) then
      close (ilog)
    end if

    if (converged) then
!>--- if the relaxation converged properly do this
      iostatus = 0
      if (pr) then
        call rmsd(mol%nat,mol%xyz,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)
        write (*,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
          "GEOMETRY OPTIMIZATION CONVERGED AFTER",iter,"ITERATIONS"
        write (*,'(72("-"))')
        write (*,'(1x,"total energy gain   :",F18.7,1x,"Eh",F14.4,1x,"kcal/mol")') &
          etot-estart, (etot-estart)*autokcal
        write (*,'(1x,"total RMSD          :",F18.7,1x,"a0",F14.4,1x,"Å")') &
          rmsdval,rmsdval*autoaa
        write (*,'(72("-"))')
      end if
    else if (iostatus .ne. 0) then
!>--- if iostatus =/= 0, something went wrong in the relaxation
      if (pr) then
        write (*,'(/,3x,"***",1x,a,1x,"***",/)') &
          "GEOMETRY RELAXATION FAILED"
      end if
    else
!>--- not converging in the given cycles is considered a FAILURE
      !> some iostatus>0 is selected to signal this
      iostatus = iter
      if (pr) then
        write (*,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
          "FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN",iter,"ITERATIONS"
      end if
    end if

!>--- overwrite input structure with optimized one
    mol%nat = molopt%nat
    mol%at = molopt%at
    mol%xyz = molopt%xyz

!> deallocate data
    !$omp critical
    if (allocated(gdiff)) deallocate (gdiff)
    if (allocated(grdold)) deallocate (grdold)
    if (allocated(ddispl)) deallocate (ddispl)
    if (allocated(displ)) deallocate (displ)
!>--- LOCAL FIX (2026-10-02)
    if (allocated(grd0)) deallocate (grd0)
    if (allocated(grdt)) deallocate (grdt)
    if (allocated(grmsd)) deallocate (grmsd)
    if (allocated(molopt%at)) deallocate (molopt%at)
    if (allocated(molopt%xyz)) deallocate (molopt%xyz)
    !$omp end critical

    return
  end subroutine gradientdescent

!========================================================================================!

  function bbstep(nopt,ddispl,gdiff) result(alp)
!*************************************************************
!* Calculate a good step size via the Barzilai-Borwein method
!*************************************************************
    implicit none
    real(wp) :: alp
    integer,intent(in) :: nopt
    real(wp),intent(in) :: ddispl(nopt)
    real(wp),intent(in) :: gdiff(nopt)
    real(wp) :: dum1,dum2
    dum1 = dot_product(ddispl,gdiff)
    dum2 = dot_product(gdiff,gdiff)
!>--- LOCAL FIX (2026-10-02): guard the division.  dum2 is the squared norm of
!> the gradient CHANGE between two accepted steps.  When that change is below
!> the noise of the level of theory (a loose DFT grid, an SCF that is not
!> converged to the last digit, a numerically flat PES) dum2 is ~0 and alp
!> becomes ±Inf or NaN, and every min() with a NaN argument is undefined, so the
!> poisoned coordinates went straight into the next singlepoint.  dum1 ~ 0 (a step that
!> is effectively orthogonal to its own gradient change, i.e. pure noise) is
!> just as useless: it yields alp ~ 0, i.e. a displacement below the geometry
!> resolution, which the strict `lowered` test can never accept.  Fall back to
!> a plain steepest-descent step (alp = 1) in both cases.
    if (dum2 .lt. 1.0e-14_wp.or.abs(dum1) .lt. 1.0e-14_wp) then
      alp = 1.0_wp
      return
    end if
!>--- abs(): a negative BB step means negative local curvature.  The sign must
!> not reach the displacement (that would turn the descent direction uphill),
!> so only the magnitude is used; the displacement trust radius in
!> gradientdescent() keeps the resulting (possibly large) step in bounds.
    alp = dum1/dum2
    alp = abs(alp)
  end function bbstep

!========================================================================================!
!========================================================================================!
end module gradientdescent_module
