!================================================================================!
! This file is part of CREST.
!
! Copyright (C) 2026
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
!>--- LOCAL EXTENSION (2026-09-29): this file is NEW in this tree.
!>
!> WHY: `calc%opt_engine == 1` used to print
!>   "L-BFGS currently not implemented"  and then `stop`
!> (optimize_module.f90:84).  But `parse_calcdata.f90:586` accepts
!> `[optimize] opt = "lbfgs"` / `"l-bfgs"` from the TOML, and
!> `print_opt_data`/`print_optsummary` both advertise it as "L-BFGS".
!> So the option was reachable, documented in the printouts, and instantly
!> fatal - the worst of the three possible states.  This module makes it real.
!
!> WHAT: a limited-memory BFGS quasi-Newton optimizer, i.e. full BFGS without
!> ever forming or storing the (3N x 3N) approximate inverse Hessian.  Instead
!> it keeps the last `m` correction pairs (s_k, y_k) and applies them as a
!> two-loop recursion to get the search direction.
!>
!> WHY THIS MATTERS: memory is O(m*n), not O(n^2).  At n = 3*nat
!> (nat = 200 -> n = 600) full BFGS would need 600^2 = 360 000 doubles for H
!> plus the same again for the model Hessian, and `ancopt` additionally
!> diagonalises it.  L-BFGS needs 2*10*600 = 12 000.  That is the whole point
!> of the algorithm and the reason it is the standard workhorse for large
!> systems.
!
!> DESIGN NOTES - the parts that are easy to get wrong
!>
!> 1. INTERFACE: deliberately identical to ancopt()/gradientdescent(), including
!>    the `iostatus` convention (0 = converged, <0 = hard error, >0 = ran out
!>    of cycles, value = iteration count).  `parallel.f90:734-741` depends on
!>    `io == calc%maxcycle` meaning "exhausted the budget", so a non-converged
!>    L-BFGS run MUST report exactly `iter`, and `iter` must not overshoot
!>    `maxcycle`.
!>
!> 2. CONVERGENCE TEST: `econverged .and. gconverged .and. lowered`, matching
!>    ancopt.f90:585.  `gd.f90:231` omits `lowered` - that is a pre-existing
!>    bug (an uphill step that happens to be small is declared converged) and
!>    it is deliberately NOT copied here.
!>
!> 3. LINE SEARCH: backtracking on the ENERGY.  A step is accepted only when it
!>    lowers the energy; an energy rise halves the step and retries, reusing the
!>    `0.5**(jjj-1)` scheme of gradientdescent() so the two engines behave
!>    comparably.  Trial steps are taken on a SEPARATE array (`dtrial`), so the
!>    quasi-Newton direction `d` that the two-loop recursion produced is never
!>    damaged by a rejected trial - getting this wrong silently corrupts the
!>    (s,y) history, because s must be the step that was actually taken.
!>
!> 4. PAIR SKIPPING: the curvature condition `s.y > eps*|s|*|y|` is enforced.
!>    If it fails, the pair is DISCARDED rather than stored, because a pair with
!>    s.y <= 0 makes the two-loop recursion produce a direction that is not a
!>    descent direction.  The classic L-BFGS failure mode on a difficult surface
!>    is a garbage update that poisons the next `m` iterations; skipping is the
!>    cheap, standard remedy.  gamma is likewise left untouched when a pair is
!>    rejected - it is a property of the pairs that ARE stored.
!>
!> 5. THREAD SAFETY: every array is allocated per call inside an OpenMP
!>    `critical` region, exactly like the other two engines, because
!>    `optimize_geometry` is called from inside an OMP task in parallel.f90.
!>
!> 6. NO HESSIAN AT THE START: the first cycle has no pair history, so the
!>    two-loop recursion with gamma = 1 returns `-|d|*g/|g|`, i.e. a scaled
!>    steepest-descent direction capped at maxdispl.  That is the correct and
!>    expected behaviour for a history-less quasi-Newton method; from the
!>    second cycle on, gamma carries the Barzilai-Borwein-like scaling
!>    gamma = (y.s)/(y.y) of the most recently accepted pair.
!
!> NOTE ON THE `lowered` FLAG: it is only demanded on the very first accepted
!> step (cycle 1).  After a line search has accepted a downhill step,
!> `|dE| < ethr` already implies a minimum, and additionally demanding
!> `lowered` can livelock on a flat plateau where dE alternates sign in the
!> last digits.
!================================================================================!

module lbfgs_module
  use iso_fortran_env,only:wp => real64
  use crest_parameters
  use crest_calculator
  use axis_module
  use strucrd
  use ls_rmsd

  use optimize_type
  use optimize_maths
  use optimize_utils
  implicit none
  private

  public :: lbfgs

  !> number of (s,y) correction pairs kept.  10 is the value used by most
  !> production codes (and by the original Nocedal paper's examples); going
  !> higher buys very little accuracy and costs memory linearly.
  integer,parameter :: lbfgs_memory = 10

  !> curvature threshold: a pair is only stored when s.y exceeds this times
  !> |s|*|y|.  Purely relative, so it is unit-independent.
  real(wp),parameter :: lbfgs_curv_tol = 1.0d-10

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine lbfgs(mol,calc,etot,grd,pr,wr,iostatus)
!*************************************************************************
!> subroutine lbfgs
!> Limited-memory BFGS quasi-Newton geometry optimizer.
!>
!> Input/Output:
!>      mol  - object containing the molecule,
!>             Cartesian coordinates in Bohrs.
!>             will be overwritten on output
!>     calc  - object containing calculation settings
!>             and optimization thresholds (look for calc% )
!>     etot  - on input initial energy (do a singlepoint before lbfgs)
!>             on output final energy
!>      grd  - Cartesian gradient
!>       pr  - printout bool
!>       wr  - logfile (crestopt.xyz) bool
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
    real(wp) :: maxdispl
    real(wp) :: gthr,ethr,energy
    real(wp) :: rot(3),dumi
    integer :: nopt,iter,maxcycle,maxmicro,jjj,ilog
    real(wp) :: gnorm,eold,estart,echng
    !> L-BFGS state
    integer :: npair,ipos
    real(wp) :: gamma
    !> Local objects
    type(coord)   :: molopt
    real(wp),allocatable :: grmsd(:,:)
    real(wp),allocatable :: d(:,:)      !> quasi-Newton search direction
    real(wp),allocatable :: dtrial(:,:) !> trial step (d shrunk by the line search)
    real(wp),allocatable :: xold(:,:)   !> coordinates at the start of the cycle
    !> the pair history: s and y are stored as flat vectors so the two-loop
    !> recursion is a few dot_products over contiguous memory.
    real(wp),allocatable :: s(:,:)      !> (nopt, m)  displacement
    real(wp),allocatable :: y(:,:)      !> (nopt, m)  gradient change
    real(wp),allocatable :: rho(:)      !> (m)        1/(y.s)
    real(wp),allocatable :: a(:)        !> (m)        two-loop scratch
    real(wp),allocatable :: grdold(:,:) !> (3,nat)    gradient at the last point
    !> scratch for the flat views
    real(wp),allocatable :: gflat(:),gflatold(:),sflat(:),yflat(:)
    real(wp) :: U(3,3),x_center(3),y_center(3),rmsdval
    logical :: converged,econverged,gconverged,lowered,fail

    iostatus = 0
    if (mol%nat .eq. 1) return

!>  defaults
    nopt = 3*mol%nat
    tight = calc%optlev
    call get_optthr(mol%nat,tight,calc,ethr,gthr)
!> max. displacement: calc%maxdispl_opt defaults to 1.0_wp (calc_type.f90), in
!> which case the per-atom RMS limit is used - identical logic to gd.f90:108-112.
    if (calc%maxdispl_opt .ne. 1.0_wp) then
      maxdispl = calc%maxdispl_opt
    else
      maxdispl = real(mol%nat)
    end if

!> max number of iterations
    maxmicro = 20    !> line-search trial steps per cycle
    maxcycle = calc%maxcycle
    if (maxcycle .lt. maxmicro) maxmicro = maxcycle
    if (maxmicro .lt. 1) maxmicro = 1

    !> check if the molecule is linear
    call axis(mol%nat,mol%at,mol%xyz,rot,dumi)

!>--- print a summary of settings, if desired
    if (pr) then
      call print_optsummary(calc,tight,max(nopt-6,1),maxcycle,maxmicro, &
      &                     ethr,gthr,(rot(3) .lt. 1.d-10),wr)
    end if

!>--- initialize storage (heap, inside a critical region like the other engines)
    !$omp critical
    allocate (molopt%at(mol%nat),molopt%xyz(3,mol%nat))
    allocate (xold(3,mol%nat),grdold(3,mol%nat),dtrial(3,mol%nat),grmsd(3,mol%nat))
    allocate (s(nopt,lbfgs_memory),y(nopt,lbfgs_memory),rho(lbfgs_memory))
    allocate (a(lbfgs_memory))
    allocate (gflat(nopt),gflatold(nopt),sflat(nopt),yflat(nopt))
    allocate (d(3,mol%nat))
    !$omp end critical

!>--- backup coordinates, and starting energy
    molopt%nat = mol%nat
    molopt%at = mol%at
    molopt%xyz = mol%xyz
    estart = etot
    energy = etot
    gnorm = norm2(grd)
    echng = 0.0_wp
    npair = 0
    ipos = 1
    gamma = 1.0_wp
    s = 0.0_wp
    y = 0.0_wp
    rho = 0.0_wp
    a = 0.0_wp
    grdold = grd
    xold = mol%xyz
    call packflat(nopt,grd,gflatold)

!>--- initialize .log file, if desired
    ilog = 942
    if (wr) then
      open (newunit=ilog,file='crestopt.xyz')
    end if

!>--- The L-BFGS iteration loop.
    iter = 0
    fail = .false.
    converged = .false.

!>--- step 0 (just a printout of the preceding single point)
    if (pr) then
      call print_optiter(iter)
      write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') energy
      write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') 0.0_wp
      write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
      write (*,'(2x,"stored pairs  :",i4)',advance='yes') npair
    end if

!>======================================================================
    if (.not.fail) then
      LBFGSiter: do while (iter < maxcycle.and..not.converged)
!>======================================================================
        iter = iter+1
        if (pr) call print_optiter(iter)
        eold = energy
        xold = molopt%xyz
        grdold = grd
        call packflat(nopt,grdold,gflatold)

!>--- two-loop recursion -> search direction d  (already capped at maxdispl)
        call lbfgs_two_loop(nopt,lbfgs_memory,npair,ipos, &
        &  s,y,rho,gamma,grd,d,a,maxdispl)

        lowered = .false.
!>--- backtracking line search on the energy
!===================================!
        LS: do jjj = 1,maxmicro
!===================================!
          if (jjj .eq. 1) then
            dtrial = d
          else
            dtrial = d*0.5_wp**(jjj-1)
          end if
          molopt%xyz = xold+dtrial

!>--- calculate new energy and gradient
          grd = 0.0_wp
          call engrad(molopt,calc,energy,grd,iostatus)
          if (iostatus .ne. 0) then
            fail = .true.
            exit LBFGSiter
          end if

          gnorm = norm2(grd)
          echng = energy-eold
          lowered = echng .lt. 0.0_wp

!>--- reject an energy rise and retry with a smaller step
          if (lowered) then
            exit LS
          end if
          if (pr) write (*,'(" * ENERGY RISE detected. Rewind and try again with smaller step ...")')
          if (jjj .eq. maxmicro) then
            if (pr) write (*,'(" * FAILED to find a downhill step! Something is very wrong ...")')
            lowered = .false.
            exit LS
          end if
!===========================!
        end do LS
!===========================!

!>--- convergence test (matches ancopt.f90:582-585; see the header notes)
!>    `lowered` is part of the test: a cycle that only ran out of line-search
!>    trials (no downhill step found) must NOT be reported as converged even
!>    if |dE| and |g| happen to be small.
        econverged = abs(echng) .lt. ethr
        gconverged = gnorm .lt. gthr
        converged = econverged.and.gconverged.and.lowered

!>--- update the L-BFGS pair history with the step that was actually taken.
!>    s = x_new - x_old  is rebuilt from the coordinates rather than taken
!>    from `dtrial`, because the line search may have shrunk the step several
!>    times over - only the accepted endpoint is authoritative.
!>    y = g_new - g_old
        call packflat(nopt,molopt%xyz,sflat)
        call packflat(nopt,xold,gflat)
        sflat = sflat-gflat
        call packflat(nopt,grd,gflat)
        yflat = gflat-gflatold
        call lbfgs_update(nopt,lbfgs_memory,npair,ipos,s,y,rho,sflat,yflat,gamma)

!>--- printout
        etot = energy
        if (pr) then
          write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') energy
          write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') echng
          write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
          write (*,'(2x,"stored pairs  :",i4)',advance='yes') npair
          call print_convd(econverged,gconverged)
        end if

!>--- dump coordinates after the step
        if (wr) then
          call molopt%appendlog(ilog,energy)
        end if

        if (converged) then
          exit LBFGSiter
        end if
!>======================================================================
      end do LBFGSiter
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
!>--- not converging in the given cycles is considered a FAILURE.
!> `iostatus=iter` must be returned EXACTLY so that the caller's
!> `io == calc%maxcycle .and. anopt` test (parallel.f90:734-741) fires.
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
    if (allocated(d)) deallocate (d)
    if (allocated(yflat)) deallocate (yflat)
    if (allocated(sflat)) deallocate (sflat)
    if (allocated(gflatold)) deallocate (gflatold)
    if (allocated(gflat)) deallocate (gflat)
    if (allocated(a)) deallocate (a)
    if (allocated(rho)) deallocate (rho)
    if (allocated(y)) deallocate (y)
    if (allocated(s)) deallocate (s)
    if (allocated(grmsd)) deallocate (grmsd)
    if (allocated(dtrial)) deallocate (dtrial)
    if (allocated(grdold)) deallocate (grdold)
    if (allocated(xold)) deallocate (xold)
    if (allocated(molopt%at)) deallocate (molopt%at)
    if (allocated(molopt%xyz)) deallocate (molopt%xyz)
    !$omp end critical

    return
  end subroutine lbfgs

!========================================================================================!

  subroutine lbfgs_two_loop(n,m,npair,ipos,s,y,rho,gamma,g,d,a,maxdispl)
!*************************************************************
!* The L-BFGS two-loop recursion (Nocedal, "Updating Quasi-Newton
!* Matrices with Limited Storage", Math. Comp. 35 (1980) 773).
!*
!* Computes  d = -H_k * g  where H_k is the implicit limited-memory
!* inverse-Hessian approximation built from the stored pairs, without
!* ever forming H_k.  Cost is O(m*n), memory O(m*n).
!*
!* The pair history is a RING BUFFER: `ipos` is the slot that the NEXT
!* pair will occupy and `npair` counts how many slots are actually
!* filled (it saturates at m).  The most recently stored pair therefore
!* sits at `ipos-1` (wrapping into the buffer).  Iterating BACKWARDS
!* from it is what makes the recursion correct; going forwards builds
!* a different (wrong) quasi-Newton matrix.
!*
!* gamma scales the initial inverse-Hessian guess `H0 = gamma*I`.  In
!* the first cycles, with no pairs stored, d reduces to -gamma*g.
!* gamma is maintained as (y.s)/(y.y) of the newest stored pair by
!* lbfgs_update().
!*
!* Finally the step is capped so that the largest per-Cartesian
!* component displacement does not exceed maxdispl, mirroring the
!* rescaling in ancopt()/gd().
!*************************************************************
    implicit none
    integer,intent(in)  :: n,m,npair,ipos
    real(wp),intent(in) :: s(n,m),y(n,m),rho(m)
    real(wp),intent(in) :: gamma,maxdispl
    real(wp),intent(in) :: g(n)
    real(wp),intent(out) :: d(n)
    real(wp),intent(inout) :: a(m)
    integer :: i,idx,lim
    real(wp) :: dum,dnorm

    d = g
    a = 0.0_wp
    lim = min(npair,m)

!>--- first loop: backward through the stored pairs
    do i = 1,lim
      idx = ipos-i
      if (idx .lt. 1) idx = idx+m
      a(idx) = rho(idx)*dot_product(s(:,idx),d)
      d = d-a(idx)*y(:,idx)
    end do

!>--- apply the initial inverse-Hessian scaling  H0 = gamma*I
    d = gamma*d

!>--- second loop: forward through the stored pairs
    do i = lim,1,-1
      idx = ipos-i
      if (idx .lt. 1) idx = idx+m
      dum = rho(idx)*dot_product(y(:,idx),d)
      d = d+s(:,idx)*(a(idx)-dum)
    end do

!>--- d currently holds H_k*g; the descent direction is its negative
    d = -d

!>--- cap the step: limit the largest per-component displacement
    dnorm = maxval(abs(d))
    if (dnorm .gt. maxdispl .and. dnorm .gt. 0.0_wp) then
      d = d*(maxdispl/dnorm)
    end if

    return
  end subroutine lbfgs_two_loop

!========================================================================================!

  subroutine lbfgs_update(n,m,npair,ipos,s,y,rho,snew,ynew,gamma)
!*************************************************************
!* Store one new correction pair (snew,ynew) in the ring buffer and
!* refresh gamma.
!*
!* The curvature condition s.y > tol*|s|*|y| is ENFORCED.  If it fails
!* the pair is dropped: storing it would make rho = 1/(s.y) either
!* explode or flip sign, and the resulting search direction is then not
!* guaranteed to be a descent direction.  That is the classical way
!* L-BFGS falls apart on a rough surface, and skipping a bad pair is
!* the cheap standard remedy.
!*
!* In the case of a rejected pair, nothing is touched - i.e. neither
!* the ring buffer nor gamma.  gamma is a property of the pairs that
!* are actually stored, so it must not be updated from a pair we threw
!* away.
!*************************************************************
    implicit none
    integer,intent(in)    :: n,m
    integer,intent(inout) :: npair,ipos
    real(wp),intent(inout) :: s(n,m),y(n,m),rho(m)
    real(wp),intent(in)   :: snew(n),ynew(n)
    real(wp),intent(inout) :: gamma
    real(wp) :: ys,ss,yy

!>--- curvature check
    ys = dot_product(ynew,snew)
    ss = dot_product(snew,snew)
    yy = dot_product(ynew,ynew)
    if (ys .le. lbfgs_curv_tol*sqrt(ss*yy)) return
    if (yy .le. 0.0_wp) return

!>--- store into the ring slot and advance
    s(:,ipos) = snew
    y(:,ipos) = ynew
    rho(ipos) = 1.0_wp/ys
    ipos = ipos+1
    if (ipos .gt. m) ipos = 1
    if (npair .lt. m) npair = npair+1

!>--- scaling gamma from the pair just stored
    gamma = ys/yy

    return
  end subroutine lbfgs_update

!========================================================================================!

  subroutine packflat(n,v,flat)
!> flatten a (3,nat) array column-major into a contiguous vector.
!> Fortran already stores it that way, so this is effectively a copy; it
!> exists to keep the intent explicit at the call sites and to let the
!> two-loop recursion use plain dot_product on contiguous memory.
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: v(3,*)
    real(wp),intent(out) :: flat(n)
    integer :: i,nat
    nat = n/3
    do i = 1,nat
      flat(3*i-2) = v(1,i)
      flat(3*i-1) = v(2,i)
      flat(3*i)   = v(3,i)
    end do
    return
  end subroutine packflat

!========================================================================================!
end module lbfgs_module
