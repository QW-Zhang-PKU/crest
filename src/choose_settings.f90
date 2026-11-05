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

!> The routines contained here take care of setting MTD runtimes
!> and how to set corresponding bias parameters.
!> Naturally, this is quite imporant for the overall performance of CREST

!=========================================================================================!
!=========================================================================================!
subroutine md_length_setup(env)
!***********************************************************
!* set the total run time according to flexibility measures
!***********************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use zdata, only:readwbo
  implicit none
  !> IN/OUTPUT
  type(systemdata) :: env    !> MAIN STORAGE OS SYSTEM DATA
  !> LOCAL
  real(wp) :: total,minimum,lenthr
  real(wp) :: flex,av1,rfac,nciflex
  type(coord) :: mol 
  logical :: ex
  integer :: flexdbg,io
  character(len=32) :: envstr
  external :: crest_hard_bonds_setup
!>--- LOCAL EXTENSION (2026-09-29): EXPLICIT INTERFACE for flexi().
!>
!> flexi lives in src/flexi.F90 as an EXTERNAL subroutine (it is not part of a
!> module and there is no interface block anywhere in the tree).  Since
!> 2026-09-27 it takes two OPTIONAL dummies (nfix, fixat) for the bonds declared
!> via --cinp, and Fortran REQUIRES an explicit interface to call a procedure
!> that has an optional dummy.  Without one gfortran rejects the call with
!>
!>   Error: Explicit interface required for 'flexi' at (1): optional argument
!>
!> NOTE ON THE ERROR MESSAGE: gfortran reports that against BOTH call sites
!> below, and for the 4-argument "else" branch it prints
!>
!>   choose_settings.f90:112:57: Error: Missing actual argument for argument
!>   '_formal_7' at (1)
!>
!> which points at the wrong line and names a dummy that does not exist in the
!> source ('_formal_7' is gfortran's internal name for the 7th dummy).  The
!> 4-argument call is perfectly legal - the real problem is the missing
!> interface.  Do not "fix" it by padding the call with dummy arguments.
!>
!> KEEP IN SYNC with the declaration in src/flexi.F90:24-44.  Only the argument
!> TYPES are needed here, but reproducing the intent(in) attributes keeps the
!> two in step and makes a mismatch a compile error rather than silent UB.
  interface
    subroutine flexi(mol,rednat,includeAtom,flex,nfix,fixat)
      use iso_fortran_env,only:wp => real64
      use strucrd
      implicit none
      type(coord),intent(in) :: mol                  !> input molecule
      integer,intent(in)     :: rednat               !> should equal sum(includeAtom(:))
      integer,intent(in)     :: includeAtom(mol%nat) !> 0 if not considered, 1 if considered
      real(wp),intent(out)   :: flex                 !> OUTPUT: flexibility 0..1
      integer,intent(in),optional :: nfix            !> number of declared bonds
      integer,intent(in),optional :: fixat(:,:)      !> (2,nfix) 1-based FULL indices
    end subroutine flexi
  end interface
!> get reference geometry  
  call env%ref%to( mol ) 
!>--- LOCAL EXTENSION: report the rigid (SHAKE) bonds coming from --cinp BEFORE
!> the flexibility banner below.  They are part of the input the flexibility
!> measure is computed from, so the reader needs them first; env_to_mddat()
!> picks the (already filled) list up again later.
  call crest_hard_bonds_setup(env)

!> at least 5ps per MTD
  minimum = 5.0d0
!> Maximum of 200 ps, longer runs can only be conducted by user input
  lenthr = 200.0d0

!>--- optional: explain the flexibility measure BEFORE it is computed
!>    (CREST_FLEXI_DEBUG=1).  Read BY VALUE: a present-but-unparsable value
!>    must fall back to 0, and a failed internal read leaves flexdbg undefined.
  flexdbg = 0
  call get_environment_variable('CREST_FLEXI_DEBUG',value=envstr,status=io)
  if (io == 0) then
    if (len_trim(envstr) > 0) read (envstr,*,iostat=io) flexdbg
    if (io /= 0) flexdbg = 0
  end if
  if (flexdbg > 0) then
    call smallhead('How the MTD length is derived')
    write (stdout,'(1x,a)') 'covalent flexibility   (src/flexi.F90 :: flexi)'
    write (stdout,'(3x,a)') 'for every bond i-j with both atoms CN > 1.2 and neither terminal:'
    write (stdout,'(5x,a)') 'branch  = 2/sqrt(deg_i*deg_j)      branching penalty'
    write (stdout,'(5x,a)') 'ringf   = 1-(1-0.5*(1-exp(-0.06*L)))*rho   L = ring size'
    write (stdout,'(5x,a)') 'rho     = continuous ring membership 0..1, from the'
    write (stdout,'(5x,a)') '          strongest cycle closing the bond (CREST_FLEXI_RINGMODE)'
    write (stdout,'(5x,a)') 'doublef = 1-exp(-4*(WBO-2)^6)      double-bond term, cliff at WBO 1.2-1.5'
    write (stdout,'(5x,a)') 'hybf    = x0.5 for each sp2 carbon (CN < 3.3)'
    write (stdout,'(3x,a)') 'bonds declared with --cinp are ALWAYS in the graph,'
    write (stdout,'(3x,a)') 'even if longer than rcov_i+rcov_j (forming/breaking bonds)'
    write (stdout,'(5x,a)') 'gate: r <= 1.3*(rcov_i+rcov_j), CREST_FLEXI_FORCEFAC,'
    write (stdout,'(5x,a)') '      needed because "atoms: a,b,c" registers ALL pairs'
    write (stdout,'(3x,a)') 'val  = min(branch*ringf*doublef*hybf, cap)   cap = 0.9'
    write (stdout,'(3x,a)') 'flex = ( sum(val^p) / m )^(1/p)    m = bonds, p = 1.5'
    write (stdout,'(3x,a)') '       (CREST_FLEXI_VALCAP / CREST_FLEXI_PMEAN;'// &
    &                       ' cap=1.0 and p=2.0 restore the original)'
    write (stdout,'(1x,a)') 'combination and MD length:'
    write (stdout,'(3x,a)') 'flex   = 0.5*flex + 0.5*nci*sqrt(flex)'
    write (stdout,'(3x,a)') 't(MTD) = 3.0*exp(0.10*flex*max(1,rednat-8))   [ps]'
    write (stdout,'(3x,a)') 'the 5 ps floor is applied BEFORE the --mdlen scaling'
  end if

  call smallhead('Generating MTD length from a flexibility measure')

  if ((env%crestver .ne. crest_solv).and..not.env%NCI) then
    write(stdout,'(1x,a)',advance='no') 'Calculating GFN0-xTB WBOs   ...'
!>-- xtb singlepoint to get WBOs (always GFN0)
    call xtbsp(env,0)
    write (stdout,'(1x,a)') 'done.'
!>-- save those WBOs to the reference
    inquire(file='wbo',exist = ex)
    if(ex)then
    if(.not.allocated(env%ref%wbo)) allocate(env%ref%wbo( mol%nat, mol%nat), source=0.0_wp)   
    call readwbo('wbo',mol%nat, env%ref%wbo)
    endif

!>-- covalent flexibility measure based on WBO and structure only
!>-- pass the bonds the user declared via --cinp: they are inserted into the
!>   bond graph no matter how long they are, so a ring that is still forming
!>   is still recognised as a ring
    if (allocated(env%hcatms)) then
      call flexi( mol, env%rednat, env%includeRMSD, flex, env%nhardcons, env%hcatms)
    else
      call flexi( mol, env%rednat, env%includeRMSD, flex)
    end if
!>-- NCI flexi based on E(HB)/Nat and E(disp)/Nat
    call nciflexi(env,nciflex)
    write (stdout,'(1x,''    covalent flexibility measure :'',f8.3)') flex
    write (stdout,'(1x,''non-covalent flexibility measure :'',f8.3)') nciflex
    if (flexdbg > 0) then
      write (stdout,'(3x,a,f10.6)') 'covalent part                : ',flex
      write (stdout,'(3x,a,f10.6)') 'non-covalent part            : ',nciflex
    end if
!>-- the NCI flex is only relevant if the covalent framework is flexible
    flex = 0.5*flex+0.5*nciflex*sqrt(flex)
    if (flexdbg > 0) then
      write (stdout,'(3x,a,f10.6)') 'combined 0.5*c+0.5*n*sqrt(c) : ',flex
    end if

    if (env%entropic) then
!>--- special case for entropy mode that depends more strongly on the flexibility than on size
!>--- -8 accounts for small systems having no conf.
      av1 = (flex**1.333333)*max(1,env%rednat-8)
!>--- Maximum of 3000 ps, longer runs can only be conducted by user input
      lenthr = 3000.0d0
!>--- minimum is 5 ps set above
      env%tmtd = 4.50d0*exp(0.165d0*av1)
    else
!>--- normal case
      av1 = (flex**1.000000)*max(1,env%rednat-8)
!>--- Maximum of 500 ps
      lenthr = 500.0d0
      env%tmtd = 3.0d0*exp(0.10d0*av1)
    end if
  else
    write (stdout,'(1x,"System flexiblity is set to 1.0 for NCI mode")')
    flex = 1.0d0
    env%tmtd = 0.10*(env%rednat+0.1*env%rednat*env%rednat)
  end if
  env%flexi = flex
  write (stdout,'(1x,''flexibility measure :'',f8.3)') env%flexi
  if (flexdbg > 0) then
    write (stdout,'(3x,a,i0)') 'rednat                       : ',env%rednat
    write (stdout,'(3x,a,f10.3)') 't(MTD) = 3*exp(0.1*f*(r-8))  : ',env%tmtd
  end if

!>-- rfac is used to scale the total MD length according to special runtypes
  total = max(minimum,env%tmtd)
  select case (env%runver)
  case (2,5,6,33)   !> "-quick","-squick","-mquick"
    rfac = 0.5d0
  case (3)          !> "-qcg"
    rfac = 0.25d0
  case (77)
    rfac = 1.50d0
  case (8)
    rfac = 2.0d0
  case (787878)
    rfac = 14.0d0/float(env%nmetadyn)
  case default  !> everything else 1=default, 4=NCI
    rfac = 1.0d0
  end select
!>-- additional user set scaling
  if (env%scallen) then
    write (stdout,'(1x,''t(MTD) based on flexibility :'',  f8.1)') env%tmtd*rfac
    rfac = env%mdlenfac*rfac
    write (stdout,'(1x,''MTD length is scaled by     :'',  f6.3)') env%mdlenfac
  end if

!>-- ONLY use generated MD length if not already set by the user
  if (env%mdtime .le. 0.0d0) then
    if(env%mddat%length_ps > 0.0_wp)then
      total = env%mddat%length_ps
      write(stdout,'(1x,"t(MTD) / ps set via calculator :",  f8.1)') total
    else if (total .gt. lenthr) then
      total = lenthr
      call mtdwarning(lenthr*rfac)
    end if
    env%mdtime = anint(total)*rfac
  else
    write (stdout,'(1x,''t(MTD) / ps set by command line  :'',  f8.1)') env%mdtime
  end if

  write (stdout,'(1x,''t(MTD) / ps    :'',  f8.1)') env%mdtime
  write (stdout,'(1x,''Σ(t(MTD)) / ps :'',  f8.1,'' ('',i0,'' MTDs)'')') &
  & env%mdtime*float(env%nmetadyn),env%nmetadyn

!> A MTD Vbias snapshot is taken every 1 ps
  if(allocated(env%metadlist))then
    env%metadlist(:) = ceiling(env%mdtime)
  endif

  return
end subroutine md_length_setup

!========================================================================================!
!========================================================================================!
subroutine defaultGF(env)
!************************************************************
!* Setmetadynamics default Guiding Force Parameter
!* There are different combinations depending on the runtype
!************************************************************
  use crest_parameters 
  use crest_data
  use filemod
  implicit none
  !> IN/OUTPUT
  type(systemdata) :: env
  !> LOCAL
  integer  :: ia,ik,na,nk,m,nmtdyn,nmtdynmax,nrem
  real(wp) :: alp,k
  real(wp) :: kstart,kinc
  real(wp) :: alpinc
!>--- LOCAL EXTENSION: the SECOND alpha level of the --onebatch 2x2 grid
  real(wp) :: alp2
  character(len=64) :: envstr
  type(filetype) :: biasfile
  logical :: ex
  integer :: i,io
  character(len=:),allocatable :: atmp

  nrem = 0

  if (.not.env%metadynset) then
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    if (env%readbias) then
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      inquire (file='mtdbias',exist=ex)
      if (ex) then
        atmp = ''
        call biasfile%open('mtdbias')
        call biasfile%clearblanks()
        nmtdyn = biasfile%nlines
        call env%allocate(nmtdyn)
        do i = 1,nmtdyn
          atmp = getlarg(biasfile%line(i),1)
          read (atmp,*,iostat=io) k
          if (io == 0) then
            env%metadfac(i) = k*env%rednat
          end if
          atmp = getlarg(biasfile%line(i),2)
          read (atmp,*,iostat=io) alp
          if (io == 0) then
            env%metadexp(i) = alp
          end if
        end do
        !write(*,*) env%metadfac
        !write(*,*) env%metadexp
      else
        error stop "no file 'mtdbias'"
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    else
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      select case (env%runver)
!---------- "-quick","-squick"
      case (2,5) 
        na = 3
        nk = 2
        nmtdyn = na*nk
        alp = 1.2d0 ! start value alpha
        kstart = 0.002d0 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-mquick"
      case (6) 
        na = 3
        nk = 2
        nmtdyn = na*nk
        alp = 1.0d0 ! start value alpha
        kstart = 0.002d0 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-qcg"
      case (3) 
        na = 4
        nk = 3
        nmtdyn = na*nk
        alp = 1.0d0 ! start value alpha
        kstart = 0.00125d0 ! start value k
        alpinc = (3./2.) ! increment
        kinc = (3./2.)   ! increment
!---------- "-nci"
      case (4) 
        na = 3
        nk = 2
        nmtdyn = na*nk
        alp = 1.0d0 ! start value alpha
        kstart = 0.001d0 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-singlerun"
      case (45) 
        na = 1
        nk = 1
        nmtdyn = na*nk
        alp = 1.0d0 ! start value alpha
        kstart = 0.001d0 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-relax"
      case (33)
        na = 1
        nk = 3
        nmtdyn = na*nk
        alp = 0.8d0 ! start value alpha
        kstart = 0.0030 ! start value k
        alpinc = 2.0 ! increment
        kinc = 2.0     ! increment
!---------- "-compress"
      case (77) 
        na = 3
        nk = 3
        nmtdyn = na*nk
        alp = 1.61803 ! start value alpha
        kstart = 0.005d0 ! start value k
        alpinc = 1.61803 ! increment
        kinc = 2.0     ! increment
!--------- "search_1"
      case (crest_s1,crest_mecp) 
        na = 3
        nk = 3
        nmtdyn = (na*nk)
        alp = 1.61803 ! start value alpha
        kstart = 0.0075d0 ! start value k
        alpinc = 1.61803 ! increment
        kinc = 2.0     ! increment
!---------- "-nmtd <x>"
      case (787878)
        call gfdistribute(env%nmetadyn,nk,na,nrem)
        nmtdyn = na*nk !> nrem will be substracted at the end
        alp = 1.3 ! start value alpha
        kstart = 0.0050 ! start value k
        alpinc = (5./3.) ! increment
        kinc = 1.5d0     ! increment
!---------- "-entropy"
      case (111) 
        na = 6
        nk = 4
        nmtdyn = (na*nk)
        alp = 1.61803 ! start value alpha
        kstart = 0.0075d0 ! start value k
        alpinc = 1.61803 ! increment
        kinc = 2.0     ! increment
!---------- default
      case default
        if (env%iterativeV2) then  !for the default iterative mode
          !=======================================================!
          na = 4
          nk = 3
          nrem = 0
          nmtdyn = (na*nk)+2
          call env%allocate(nmtdyn) !> allocate k(Vbias) and α(Vbias)
          alp = 1.3                 !> start value alpha
          kstart = 0.003d0          !> start value k
          alpinc = (5.0d0/3.0d0)    !> increment alpha
          kinc = 2.0d0              !> increment k

          !-- two additional MTDs with extreme values
          env%metadfac(nmtdyn-1) = 0.001*env%rednat
          env%metadexp(nmtdyn-1) = 0.1

          env%metadfac(nmtdyn) = 0.005*env%rednat
          env%metadexp(nmtdyn) = 0.8
          !=======================================================!
        else  !for the non-iterative mode
          !=======================================================!
          na = 6
          nk = 4
          nmtdyn = na*nk
          alp = 1.3          !> start value alpha
          kstart = 0.003d0   !> start value k
          alpinc = (4.0d0/3.0d0)   !> increment alpha
          kinc = (3.0d0/2.0d0)     !> increment k
          !======================================================!
        end if
      end select
!>---- LOCAL EXTENSION ("-onebatch"): override whatever grid the runtype picked
!> with a FIXED 4-point grid.  A single batch has no second MTD round that could
!> compensate for a badly placed Vbias, so the four trajectories are laid out as
!> a 2x2 rectangle in (log k, log alpha), centred on the region that both
!> upstream grids agree on (k ~ 1.5e-3, alpha ~ 0.6), with a 2x spread in each
!> direction:
!>   k     = 0.0020 / 0.0010  the -quick pair.  2.5x gentler than the -nmtd
!>                            maximum (0.005), which matters here: a stronger
!>                            push fights the 10 rigid SHAKE bonds and inflates
!>                            the per-step constraint violation.
!>   alpha = 1.00   / 0.67    bias half-width sqrt(ln2/alpha) = 0.83 / 1.02 Ang.
!>                            0.50 was the original choice and pushed structures
!>                            apart too fast - see the note on alp2 below.
!>                            alpha <= 0.3 stops telling conformers of a ~90
!>                            atom system apart; alpha >= 1.3 is more aggressive
!>                            than one unrepeated batch needs.
!> A 4x1 grid (one k, four alphas) would put all four trajectories at the same
!> push strength - if that k is wrong, all four are wasted.  The 2x2 hedges both
!> dimensions.  An explicit -nmtd <N> is an explicit user request and still wins
!> over this override (runver 787878), so the two never fight.
      if (env%onebatch.and.env%runver.ne.787878) then
        na = 2
        nk = 2
        nrem = 0
        nmtdyn = na*nk
        alp = 1.0d0          !> start value alpha
        kstart = 0.002d0     !> start value k
!>--- LOCAL EXTENSION (2026-09-25): the second alpha level is 0.67, not 0.50.
!> Reason: a hill is k*exp(-alpha*dRMSD^2), so alpha sets how WIDE one hill is -
!> half width at half height = sqrt(ln2/alpha), i.e. 1.18 A at alpha = 0.50 and
!> 1.02 A at 0.67.  A wide hill floods a broad RMSD region after only a few
!> deposits, so the structure is pushed out of the basin it started in before
!> that basin has actually been sampled, and it lands in a high-energy geometry
!> that no cheap optimizer will bring back.  Narrowing the hill by ~14 % buys the
!> trajectory more steps inside each basin at the same hill height k.
!> alpinc is a DIVISOR (the loop below does alp = alp/alpinc), so the number that
!> ends up in the grid is alp/alpinc - hence 1.0d0/alp2 rather than spelling out
!> the reciprocal 1.4925..., which says nothing to a reader.
!> Tunable without recompiling: CREST_ONEBATCH_ALPHA2=<real>, read BY VALUE.
        alp2 = 0.67d0
        call get_environment_variable('CREST_ONEBATCH_ALPHA2',value=envstr,status=io)
        if (io == 0) then
          if (len_trim(envstr) > 0) read (envstr,*,iostat=io) alp2
        end if
!> two separate ifs on purpose: a failed read leaves alp2 UNDEFINED, and .or.
!> does not guarantee short-circuit evaluation
        if (io /= 0) then
          alp2 = 0.67d0
        else if (alp2 <= 0.0d0) then
          alp2 = 0.67d0
        end if
        alpinc = 1.0d0/alp2
        kinc = 2.0d0         !> increment k
      end if
!>---- settings are generated here
      m = 0
      nmtdynmax = nmtdyn-nrem
      call env%allocate(nmtdynmax)   !allocate k(Vbias) and α(Vbias)
      do ia = 1,na
        k = kstart  ! start value k
        do ik = 1,nk
          m = m+1
          if (m > nmtdynmax) cycle !> skip the last nrem setups
          env%metadfac(m) = k*env%rednat*env%mtd_kscal
          env%metadexp(m) = alp
          k = k/kinc     !> increment k
        end do
        alp = alp/alpinc !> increment alpha
      end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++!
    end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++!
  end if
  return
!====================================================================!
contains
!====================================================================!

  subroutine gfdistribute(nsoll,k,a,rem)
!**********************************************************
!* subroutine gfdistribute
!* for a given number of requested MTDs (nsoll),
!* specify the number of different kpush (k) and
!* alpha (a) values, and from the product (k*a),
!* determine how many MTDs have to be neglected (nrem).
!* I.e.,
!*        (k*a)-nrem = nsoll
!*
!**********************************************************
    implicit none
    integer,intent(in) :: nsoll
    integer,intent(out) :: k,a,rem
    real(wp) :: nn,nsq
    k = 1
    a = 1
    rem = 0
    if (nsoll < 1) then
      error stop 'MTD setup failure.'
    end if
    nn = float(nsoll)
    nsq = sqrt(nn)
    a = ceiling(nsq)
    k = nint(nsq)
    rem = abs(nsoll-(a*k))
    return
  end subroutine gfdistribute
end subroutine defaultGF

!=========================================================================================!
!=========================================================================================!
subroutine adjustnormmd(env)
!*************************************************************************
!* Dynamically determine the number of normMDs and settings of staticMTDs
!* Set their number and the different temperatures.
!* Defaults for the static MTDs are more lengthy...
!*************************************************************************
  use crest_parameters
  use crest_data
  implicit none
  !> IN/OUTPUT
  type(systemdata) :: env
  !> LOCAL
  integer :: ndum

  if (env%rotamermds) then
!>--- first the number of normMDs on low conformers
    if (env%nrotammds .le. 0) then !> if no user input was set
      !> multiple short MDs, which has a better parallel efficiency
      !> default is 4 
      env%nrotammds = max(1,nint(float(env%nmetadyn)/4.0d0)) 
    end if

!>--- then the temperature range
    if (env%temps .le. 0) then
      !> at how many different temperatures? 
      !> starting at 400k and increasing 100K for each (200 K for -entropy mode)
      env%temps = 2
      if (env%entropic) then
        env%temps = 1
      end if
    end if
!>--- total number of NORMMDs is temps*nrotammds
  end if

!==============================================!
!>--- settings for static MTDS in entropy mode
!==============================================!
  if (env%entropymd) then 
    env%emtd%iter = 20    !> max number of iterations
    env%emtd%nbias = min(150,nint(env%tmtd/4)) !> max number of bias structures
    env%emtd%nbiasgrow = min(1.4d0,1.2d0+env%tmtd*1.d-3) !> increase of nBias in each cycle
    env%emtd%nMDs = 36          !> number of static MTDs
    env%emtd%lenfac = 0.5d0     !> length (relativ to regular MTDs)
    env%emtd%temperature = env%nmdtemp !> sMTD temperature (default 600 K)
    env%emtd%kpush = 1.d-4+env%tmtd*1.0d-6   !> kpush constant PER ATOM, a bit more for flexible systems 1.d-4+env%tmtd*1.d-6 1.5 zu viel, 0.5 zu wenig
    env%emtd%alpha = 1.0d0        !> some alpha
    env%emtd%mtdramp = 0.015d0    !> parameter to control how "fast" bias is applied in MTD
    if (env%crestver == crest_imtd) then
      if (env%emtd%confthr < 0.0d0) then
        env%emtd%confthr = 0.02d0    !> if we gain less than x% NEW conformers, exit
      end if
      if (env%emtd%sconvthr < 0.0d0) then
        env%emtd%sconvthr = 0.005d0   !> if we gain less than x% NEW entropy, exit
      end if
    end if

    if (env%nmdtemp < 0.d0) then   !> if temperature is not set by the user
      env%nmdtemp = 600.0d0
    end if

!>--- for the new alternative iMTD-sMTD runtype, re-adjust settings
    if (env%crestver == crest_imtd2) then
      ndum = 2
      env%emtd%nklist = ndum
      allocate (env%emtd%klist(ndum))
      env%emtd%klist(1) = env%emtd%kpush
      env%emtd%klist(2) = env%emtd%kpush*2.5d0
      env%emtd%nMDs = 12           !> number of static MTDs
      env%emtd%lenfac = 0.5d0      !> half the length because we have 2 kpush
      if (env%emtd%confthr < 0.0d0) then
        env%emtd%confthr = 0.05d0  !> if we gain less than x% NEW conformers, exit
      end if
      if (env%emtd%sconvthr < 0.0d0) then
        env%emtd%sconvthr = 0.01d0  !> if we gain less than x% NEW entropy, exit
      end if
    end if

!>--- Exclude atoms from static MTD bias
    env%emtd%rmax = 0    !> ignore small rings up to this size in bias
    call mtdatoms(env)
  end if
!==============================================!

  return
end subroutine adjustnormmd

!========================================================================================!
!========================================================================================!
subroutine env_to_mddat(env)
!**********************************************
!* Convert CREST's global MD settings to the
!* MD calculator object that will be used as
!* A basis for the newer calculator routines
!**********************************************
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata) :: env
  real(wp) :: dum
  external :: crest_hard_bonds_setup
!!>--- dont override user-defined settings
!  if(env%mddat%requested) return
!> we will check if any default settings were already set individually, instead
!> the if-statements in the following take care of that

!>--- necessary transfer global settings into mddat object
   if(env%mddat%length_ps <= 0.0_wp)then
   !> total runtime in ps
     env%mddat%length_ps    = env%mdtime
   else
     env%mdtime = env%mddat%length_ps
   endif
   if(env%mddat%tstep <= 0.0_wp)then
   !> time step in fs 
     env%mddat%tstep        = env%mdstep
   endif
   !> simulation steps (would be recovered automatically later, but just to make sure)
   env%mddat%length_steps = nint(env%mddat%length_ps*1000.0_wp / env%mddat%tstep)
   if(env%mddat%tsoll <= 0.0_wp)then
   !> target temperature
     env%mddat%tsoll = env%mdtemp
   endif

   if( env%mddat%dumpstep <= 0.0_wp ) then 
   !> dump frequency in fs
     env%mddat%dumpstep = float(env%mddumpxyz)
   endif
   if(env%mddat%sdump <= 0)then
   !> trajectory structure dump every x steps 
     dum = max(1.0_wp, (env%mddat%dumpstep / env%mddat%tstep))
     env%mddat%sdump = nint(dum)
   endif

!>--- LOCAL EXTENSION: bond constraints read from the constraint input file
!> become RIGID SHAKE bonds.  This must happen before the SHAKE setup below
!> because env%mddat%shk is the template that is copied into every MD
!> (setuptest.f90: MDSTART = env%mddat, parallel.f90: mddat = env%mddat,
!>  dynamics.f90: mddat = env%mddat), so filling it here reaches all of them.
!> The setup+report itself lives in crest_hard_bonds_setup() so that it can
!> also be called from md_length_setup() - there it is printed ABOVE the
!> "Generating MTD length ..." banner.  It is idempotent (nusr guard), so
!> calling it here as well is what keeps the standalone --md path working.
   call crest_hard_bonds_setup(env)

   !> The SHAKE setup (special condition referring to the default)
   env%mddat%shake = env%mddat%shake .and.(env%shake > 0) !> SHAKE algorithm?
!>--- LOCAL EXTENSION: rigid constraints REQUIRE SHAKE to be switched on, even
!> if -shake 0 was given.  Leaving shake_mode at 0 in that case is intended:
!> init_shake() then only auto-detects no bonds at all, so the result is
!> "exactly the user bonds, nothing else".
   if (env%mddat%shk%nusr > 0) env%mddat%shake = .true.
   if( env%mddat%shake .and. env%mddat%shk%shake_mode == 0)then
   env%mddat%shk%shake_mode = env%shake     !> H-only shake =1, all atom =2
   endif 

   if(env%mddat%md_hmass <= 0.0_wp)then
   !> hydrogen mass (to enable longer timesteps)
     env%mddat%md_hmass = env%hmass 
   endif

   ! TODO: WBO reader if shake is applied and wbo file is present

!>--- set flag to signal present settings
  env%mddat%requested = .true.

end subroutine env_to_mddat
!=========================================================================================!
!=========================================================================================!
subroutine crest_hard_bonds_setup(env)
!***********************************************************
!* Turn the distances declared in the --cinp / TOML
!* constraint file into RIGID (SHAKE) bonds and report them.
!*
!* WHY THIS IS A SEPARATE ROUTINE: the report belongs ABOVE the
!* "Generating MTD length from a flexibility measure" banner, i.e. before the
!* flexibility measure the reader is meant to read it next to - seeing the
!* constrained bonds first is what makes the flexibility number (and hence the
!* MTD length) interpretable.  md_length_setup() therefore calls this first.
!* env_to_mddat() calls it too, because the standalone --md / -md paths never
!* go through md_length_setup(); the nusr guard below keeps it from being done
!* (and printed) twice in the paths that do.
!***********************************************************
  use crest_parameters
  use crest_data
  use dynamics_module,only:crest_hard_bonds_to_shake
  implicit none
  type(systemdata) :: env
  if (env%hardcons.and.env%nhardcons > 0) then
    if (allocated(env%hcatms).and.allocated(env%hcdist)) then
!>--- nusr == 0 means "not filled yet".  crest_hard_bonds_to_shake() leaves it
!> at 0 when the list is empty, in which case there is nothing to report anyway.
      if (env%mddat%shk%nusr == 0) then
        call crest_hard_bonds_to_shake(env%mddat%shk,env%nhardcons, &
        &                              env%hcatms,env%hcdist,.true.)
      end if
    end if
  end if
end subroutine crest_hard_bonds_setup


