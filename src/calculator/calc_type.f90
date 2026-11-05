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

module calc_type
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use constraints
  use strucrd,only:coord
!>--- api types
  use tblite_api
  use gfn0_api
  use gfnff_api,only:gfnff_data
  use libpvol_api,only:libpvol_calculator
!>--- other types
  use orca_type
  use lwoniom_module
  implicit none

  character(len=*),public,parameter :: sep = '/'
  character(len=*),public,parameter :: dev0 = ' 2>/dev/null'

!&<
  !> job type enumerator
  type ,private:: enum_jobtype
    integer :: unknown   = 0
    integer :: xtbsys    = 1
    integer :: generic   = 2
    integer :: turbomole = 3
    integer :: orca      = 4
    integer :: terachem  = 5
    integer :: tblite    = 6
    integer :: gfn0      = 7
    integer :: gfn0occ   = 8
    integer :: gfnff     = 9
    integer :: libpvol   = 10
    integer :: lj        = 11
  end type enum_jobtype
  type(enum_jobtype), parameter,public :: jobtype = enum_jobtype()

  character(len=45),parameter,private :: jobdescription(12) = [ &
     & 'Unknown calculation type                    ', &
     & 'xTB calculation via external binary         ', &
     & 'Generic script execution                    ', &
     & 'Systemcall with Turbomole-style in/output   ', &
     & 'Systemcall to the ORCA program package      ', &
     & 'Systemcall to the TeraChem program package  ', &
     & 'xTB calculation via tblite lib              ', &
     & 'GFN0-xTB calculation via GFN0 lib           ', &
     & 'GFN0*-xTB calculation via GFN0 lib          ', &
     & 'GFN-FF calculation via GFNFF lib            ', &
     & 'external pressure calculation via libpvol   ', &
     & 'Lennard-Jones potential calculation         ' ]
!&>

!=========================================================================================!

  public :: calculation_settings
  type :: calculation_settings
!**********************************************************************
!* data object that contains the data for a *SINGLE* calculation level
!**********************************************************************
    integer :: id = 0         !> calculation type (see "jobtype" parameter above)
    integer :: prch = stdout  !> printout channel
    logical :: pr = .false.   !> allow the calculation to produce printout? Results in a lot I/O
    logical :: prappend = .false. !> append printout
    logical :: prstdout = .false. !> special case, fwd some printout to stdout
    integer :: refine_lvl = 0 !> to allow defining different refinement levels

    integer :: chrg = 0          !> molecular charge
    integer :: uhf = 0           !> uhf parameter (xtb) or multiplicity (other)
    logical :: active = .true.   !> active setting to disable the calculation (this is different from weight=0)
    real(wp) :: weight = 1.0_wp  !> calculation weight (when adding them up)

    character(len=:),allocatable :: calcspace  !> subdirectory to perform the calculation in
    character(len=:),allocatable :: calcfile
    character(len=:),allocatable :: gradfile
    character(len=:),allocatable :: path
    character(len=:),allocatable :: other
    character(len=:),allocatable :: binary      !> binary or generic script
    character(len=:),allocatable :: systemcall  !> systemcall for running generic scripts
    character(len=:),allocatable :: description !> see above jobdescription parameter
    character(len=:),allocatable :: shortflag   !> shorter job description

!>--- gradient format specifications
    logical :: numgrad = .false.      !> run numerical gradient (expensive!)
    real(wp) :: gradstep = 0.0005_wp  !> displacement for numerical gradient
    logical :: rdgrad = .true.
    integer :: gradtype = 0
    integer :: gradfmt = 0
    character(len=:),allocatable :: gradkey
    character(len=:),allocatable :: efile

!>--- results/property requests
    real(wp) :: epot = 0.0_wp
    real(wp) :: efix = 0.0_wp
    real(wp) :: etot = 0.0_wp

    !> bond orders
    logical :: rdwbo = .false.
    real(wp),allocatable :: wbo(:,:)

    !> atomic charges
    logical :: rdqat = .false.
    logical :: dumpq = .false.
    real(wp),allocatable :: qat(:)

    !> dipole and dipole gradient
    logical :: rddip = .false.
    real(wp) :: dipole(3) = 0.0_wp
    logical :: rddipgrad = .false.
    real(wp),allocatable :: dipgrad(:,:)

    !> other properties
    logical,allocatable :: getsasa(:) 
    logical :: getlmocent = .false.
    integer :: nprot = 0
    real(wp),allocatable :: protxyz(:,:) 

!>--- API constructs
    integer  :: tblitelvl = 2
    real(wp) :: etemp = 300.0_wp
    real(wp) :: accuracy = 1.0_wp
    logical  :: apiclean = .true.
    integer  :: maxscc = 500
    logical  :: saveint = .false.
    character(len=:),allocatable :: solvmodel
    character(len=:),allocatable :: solvent
!>--- Some optional file name storages
    character(len=:),allocatable :: parametrisation
    logical  :: restart = .false.  !> restart option (some potentials can do this)
    character(len=:),allocatable :: restartfile
    character(len=:),allocatable :: refgeo
    character(len=:),allocatable :: refcharges

!>--- tblite data
    type(tblite_data),allocatable :: tblite
    character(len=:),allocatable :: tbliteparam

!>--- GFN0-xTB data
    type(gfn0_data),allocatable :: g0calc
    integer :: nconfig = 0
    integer,allocatable :: config(:)
    real(wp),allocatable :: occ(:)

!>--- GFN-FF data
    type(gfnff_data),allocatable :: ff_dat

!>--- libpvol data
    integer  :: pvmodel = 1            !> libpvol model type (0=XHCFF, 1=PV)
    integer  :: ngrid = 1202           !> lebedev grid points per atom
    real(wp) :: extpressure = 0.0_wp   !> hydorstatic pressure in Gpa
    real(wp) :: proberad = 1.5_wp      !> proberadius in Angstroem
    integer  :: vdwset = 0             !> Type of VDW radii -> 0 (default) D3, 1 -> Bondi
    real(wp) :: pvradscal = 1.0_wp     !> Scaling factor for SAS radii
    type(libpvol_calculator),allocatable :: libpvol

    !> ONIOM fragment IDs
    integer :: ONIOM_highlowroot = 0
    integer :: ONIOM_id = 0

    !> ORCA job template
    type(orca_input) :: ORCA

!>--- Type procedures
  contains
    procedure :: deallocate => calculation_settings_deallocate
    procedure :: addconfig => calculation_settings_addconfig
    procedure :: autocomplete => calculation_settings_autocomplete
    procedure :: printid => calculation_settings_printid
    procedure :: info => calculation_settings_info
    procedure :: create => create_calclevel_shortcut
    procedure :: norestarts => calculation_settings_norestarts
    procedure :: dumpdipgrad => calculation_dump_dipgrad
  end type calculation_settings

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!

  public :: calcdata
  type :: calcdata
!*****************************************************************************
!* data object that collects settings for *ALL* calculations and constraints.
!*****************************************************************************
    integer :: id = 0  !> this parameter will decide how to return or add up energies and gradients
    integer :: refine_stage = 0 !> to allow iterating different refinement stages

!>--- calculations
    integer :: ncalculations = 0
    type(calculation_settings),allocatable :: calcs(:)
    real(wp),allocatable :: etmp(:)
    real(wp),allocatable :: grdtmp(:,:,:)
    real(wp),allocatable :: eweight(:)
    real(wp),allocatable :: weightbackup(:)
    logical,allocatable  :: activebackup(:)
    real(wp),allocatable :: etmp2(:)
    real(wp),allocatable :: grdtmp2(:,:,:)
    real(wp),allocatable :: eweight2(:)
    real(wp),allocatable :: grdfix(:,:)

!>--- constraints
    integer :: nconstraints = 0
    type(constraint),allocatable :: cons(:)

!>--- scans
    integer :: nscans = 0
    logical :: relaxscan = .true.
    real(wp) :: scansforce = 0.5_wp
    type(scantype),allocatable :: scans(:)

!>--- frozen atoms
    integer :: nfreeze = 0
    logical,allocatable :: freezelist(:)

!>--- results/property requests
    real(wp) :: epot
    real(wp) :: efix
    real(wp) :: etot

!>--- optimization settings
    logical  :: optnewinit = .false.  !> ensure fresh calc/param setup at beginning of opt
!>--- LOCAL EXTENSION (2026-09-27): default raised to .true..
!> An optimisation that runs out of `maxcycle` WITHOUT converging used to be
!> treated as a failed structure: ancopt returns iostatus = iter
!> (ancopt.f90:271-275), crest_oloop checks `io == maxcycle .and. anopt`
!> (parallel.f90:739) and - with anopt = .false. - fell into the else branch
!> `eread(zcopy) = 1.0_wp`, i.e. the structure was given a 1 Eh energy and
!> dropped by the energy window.  Since the global cycle budget is now a
!> deliberately small 120, that would silently throw away every conformer that
!> happens to need 121..120+ cycles.  With anopt = .true. such a structure is
!> KEPT with its partially optimised geometry and its real (unconverged)
!> energy, and CREGEN's RMSD/energy filters decide about it as usual.
!> NOTE: this is the ONLY place anopt is initialised - calc%reset() does not
!> touch it (see the reset block below), and protonate.f90 sets it explicitly.
    logical  :: anopt = .true.  !> allow collecting structures that ran into maxcycle?
    integer  :: optlev = 0
    integer  :: micro_opt = 20
    integer  :: maxcycle = 0
    real(wp) :: maxdispl_opt = 1.0_wp
    real(wp) :: ethr_opt = -1.0_wp  !> ΔE convergence
    real(wp) :: gthr_opt = -1.0_wp  !> RMS force convergence
    real(wp) :: hlow_opt = 0.010_wp
    real(wp) :: hmax_opt = 5.0_wp
    real(wp) :: acc_opt = 1.0_wp
    real(wp) :: maxerise = 1.0e-5_wp
    real(wp) :: hguess = 0.02_wp !> Hessian guess (for some optimizers)
    logical  :: exact_rf = .false.
    logical  :: average_conv = .false.
    logical  :: tsopt = .false.
    integer  :: iupdat = 0  !> 0=BFGS, 1=Powell, 2=SR1, 3=Bofill, 4=Schlegel
    integer  :: opt_engine = 0 !> default: ANCOPT

!>--- GFN0* data, needed for special MECP application
    type(gfn0_data),allocatable  :: g0calc

!>--- printouts and IO
    logical :: pr_energies = .false.
    integer :: eout_unit = stdout
    character(len=:),allocatable :: elog

!>--- LOCAL EXTENSION: "this call does not need a gradient".
!> Set by singlepoint() (i.e. by the refine singlepoint re-ranking,
!> crest_sploop).  A calculator that can skip the gradient then does so - for
!> ORCA that means the template is written WITHOUT the EnGrad keyword and the
!> energy is read back from ORCA.out, saving the CP-SCF/Z-vector solve (and the
!> CPCM gradient term) on every structure of the ensemble.
!> The energy returned is identical; only the gradient stays zero.
!> NOTE: this is a flag on the WHOLE calcdata, not on a single level, because it
!> belongs to "what the caller asked for" and not to "how this level is set up".
    logical :: nogr = .false.

!>--- ONIOM calculator data
    type(lwoniom_data),allocatable :: ONIOM
    type(coord),allocatable :: ONIOMmols(:)
    integer,allocatable :: ONIOMmap(:) !> map ONIOM fragments to calculation_settings
    integer,allocatable :: ONIOMrevmap(:) !> map calculation settings to ONIOM frags (or zero)

!>--- Type procedures
  contains
    procedure :: reset => calculation_reset
    procedure :: init => calculation_init
    generic,public :: add => calculation_add_constraint,calculation_add_settings, &
    & calculation_add_scan,calculation_add_constraintlist
    procedure,private :: calculation_add_constraint,calculation_add_settings, &
    & calculation_add_scan,calculation_add_constraintlist
    procedure :: copy => calculation_copy
    procedure :: printconstraints => calculation_print_constraints
    procedure :: removeconstraint => calculation_remove_constraint
    procedure :: info => calculation_info
    procedure :: ONIOMexpand => calculation_ONIOMexpand
    procedure :: active => calc_set_active
    procedure :: active_restore => calc_set_active_restore
    procedure :: freezegrad => calculation_freezegrad
    procedure :: increase_charge => calculation_increase_charge
    procedure :: decrease_charge => calculation_decrease_charge
    procedure :: dealloc_params => calculation_deallocate_params
  end type calcdata

  public :: get_dipoles

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine calculation_reset(self)
    implicit none
    class(calcdata) :: self

    self%id = 0

    if (allocated(self%calcs)) deallocate (self%calcs)
    self%ncalculations = 0

    if (allocated(self%cons)) deallocate (self%cons)
    self%nconstraints = 0

    self%optlev = 0
    self%micro_opt = 20
    self%maxcycle = 0
    self%maxdispl_opt = 1.0_wp
    self%ethr_opt = -1.0_wp  !> ΔE convergence
    self%gthr_opt = -1.0_wp  !> RMS force convergence
    self%hlow_opt = 0.010_wp
    self%hmax_opt = 5.0_wp
    self%acc_opt = 1.0_wp
    self%maxerise = 1.0e-5_wp
    self%exact_rf = .false.
    self%average_conv = .false.
    self%tsopt = .false.
    self%iupdat = 0  !> 0=BFGS, 1=Powell, 2=SR1, 3=Bofill, 4=Schlegel
    self%opt_engine = 0 !> default: ANCOPT

    self%pr_energies = .false.
    self%eout_unit = stdout
    if (allocated(self%elog)) deallocate (self%elog)

    if (allocated(self%etmp)) deallocate (self%etmp)
    if (allocated(self%grdtmp)) deallocate (self%grdtmp)
    if (allocated(self%etmp2)) deallocate (self%etmp2)
    if (allocated(self%grdtmp2)) deallocate (self%grdtmp2)

    if (allocated(self%g0calc)) deallocate (self%g0calc)

    if (allocated(self%ONIOM)) deallocate (self%ONIOM)
    if (allocated(self%ONIOMmols)) deallocate (self%ONIOMmols)
    if (allocated(self%ONIOMmap)) deallocate (self%ONIOMmap)

    return
  end subroutine calculation_reset

!=========================================================================================!

  subroutine calculation_deallocate_params(self)
!**********************************************
!* Deallocate parametrization/model setup
!* This only applies (or is useful) for a few
!* selected model like GFN-FF, gfn0 and tblite
!**********************************************
    class(calcdata) :: self
    integer :: i,j,k
    if (self%ncalculations > 0) then
      do i = 1,self%ncalculations
        if(allocated(self%calcs(i)%tblite)) deallocate(self%calcs(i)%tblite)
        if(allocated(self%calcs(i)%g0calc)) deallocate(self%calcs(i)%g0calc)
        if(allocated(self%calcs(i)%ff_dat)) deallocate(self%calcs(i)%ff_dat)
        if(allocated(self%calcs(i)%libpvol)) deallocate(self%calcs(i)%libpvol)
      end do
    end if
  end subroutine calculation_deallocate_params

!=========================================================================================!

  subroutine calculation_add_settings(self,cal)
    implicit none
    class(calcdata) :: self
    type(calculation_settings) :: cal
    type(calculation_settings),allocatable :: callist(:)
    integer :: i,j

    if (self%ncalculations < 1) then
      allocate (self%calcs(1))
      self%ncalculations = 1
      self%calcs(1) = cal
    else
      i = self%ncalculations+1
      j = self%ncalculations
      allocate (callist(i))
      callist(1:j) = self%calcs(1:j)
      callist(i) = cal
      call move_alloc(callist,self%calcs)
      self%ncalculations = i
    end if

    return
  end subroutine calculation_add_settings

!=========================================================================================!

  subroutine calculation_add_constraint(self,constr)
    implicit none
    class(calcdata) :: self
    type(constraint) :: constr
    type(constraint),allocatable :: conslist(:)
    integer :: i,j

    if (self%nconstraints < 1) then
      allocate (self%cons(1))
      self%nconstraints = 1
      self%cons(1) = constr
    else
      i = self%nconstraints+1
      j = self%nconstraints
      allocate (conslist(i))
      conslist(1:j) = self%cons(1:j)
      conslist(i) = constr
      call move_alloc(conslist,self%cons)
      self%nconstraints = i
    end if

    return
  end subroutine calculation_add_constraint

  subroutine calculation_add_constraintlist(self,k,constr)
    implicit none
    class(calcdata) :: self
    integer :: k
    type(constraint) :: constr(k)
    type(constraint),allocatable :: conslist(:)
    integer :: i,j

    if (self%nconstraints < 1) then
      allocate (self%cons(k))
      self%nconstraints = k
      self%cons(1:k) = constr(1:k)
    else
      i = self%nconstraints+k
      j = self%nconstraints
      allocate (conslist(i))
      conslist(1:j) = self%cons(1:j)
      conslist(j+1:i) = constr(1:k)
      call move_alloc(conslist,self%cons)
      self%nconstraints = i
    end if

    return
  end subroutine calculation_add_constraintlist

!=========================================================================================!

  subroutine calculation_add_scan(self,scn)
    implicit none
    class(calcdata) :: self
    type(scantype) :: scn
    type(scantype),allocatable :: scnlist(:)
    integer :: i,j

    if (self%nscans < 1) then
      allocate (self%scans(1))
      self%nscans = 1
      self%scans(1) = scn
    else
      i = self%nscans+1
      j = self%nscans
      allocate (scnlist(i))
      scnlist(1:j) = self%scans(1:j)
      scnlist(i) = scn
      call move_alloc(scnlist,self%scans)
      self%nscans = i
    end if

    return
  end subroutine calculation_add_scan

!=========================================================================================!

  subroutine calculation_remove_constraint(self,d)
    implicit none
    class(calcdata) :: self
    type(constraint) :: constr
    type(constraint),allocatable :: conslist(:)
    integer :: i,j,d,d1,d2

    if (self%nconstraints < d) return

    i = self%nconstraints-1
    j = self%nconstraints
    allocate (conslist(i))
    if (d == 1) then
      conslist(1:i) = self%cons(2:j)
    else if (d == j) then
      conslist(1:i) = self%cons(1:i)
    else
      d1 = d-1
      d2 = d+1
      conslist(1:d1) = self%cons(1:d1)
      conslist(d:i) = self%cons(d2:j)
    end if
    call move_alloc(conslist,self%cons)
    self%nconstraints = i

    return
  end subroutine calculation_remove_constraint

!=========================================================================================!
  subroutine calculation_print_constraints(self,chnl)
    implicit none
    class(calcdata) :: self
    integer :: i,j
    integer,optional :: chnl

    if (self%nconstraints < 1) then
      return
    else
      if (present(chnl)) then
        i = chnl
      else
        i = stdout
      end if
      !write(i,*)
      do j = 1,self%nconstraints
        call self%cons(j)%print(i)
      end do
    end if

    return
  end subroutine calculation_print_constraints

!=========================================================================================!

  subroutine calculation_init(self)
    class(calcdata) :: self

    if (allocated(self%elog)) then
      self%pr_energies = .true.
      open (newunit=self%eout_unit,file=self%elog)
    end if

  end subroutine calculation_init

!=========================================================================================!
!> copy a calcdata object from src to self
  subroutine calculation_copy(self,src)
    class(calcdata) :: self
    type(calcdata) :: src
    integer :: i

    self%id = src%id

    self%ncalculations = src%ncalculations
    if (allocated(self%calcs)) deallocate (self%calcs)
    !self%calcs = src%calcs
    do i = 1,self%ncalculations
      call self%add(src%calcs(i))
    end do

    self%nconstraints = src%nconstraints
    if (allocated(self%cons)) deallocate (self%cons)
    !self%cons = src%cons
    do i = 1,self%nconstraints
      call self%add(src%cons(i))
    end do

    self%optlev = src%optlev
    self%micro_opt = src%micro_opt
    self%maxcycle = src%maxcycle
    self%maxdispl_opt = src%maxdispl_opt
    self%hlow_opt = src%hlow_opt
    self%hmax_opt = src%hmax_opt
    self%acc_opt = src%acc_opt
    self%exact_rf = src%exact_rf
    self%average_conv = src%average_conv
    self%tsopt = src%tsopt
    self%iupdat = src%iupdat

    self%pr_energies = src%pr_energies
    self%eout_unit = src%eout_unit
    self%elog = src%elog

    return
  end subroutine calculation_copy

!=========================================================================================!

  subroutine calculation_freezegrad(self,grad)
    class(calcdata) :: self
    real(wp),intent(inout) :: grad(:,:)
    integer :: nat,i
    if (self%nfreeze > 0) then
      nat = size(grad,2)
      if (nat == self%nfreeze) then
        error stop '*** ERROR *** Must not freeze all atoms!'
      end if
      do i = 1,nat
        if (self%freezelist(i)) grad(:,i) = 0.0_wp
      end do
    end if
  end subroutine calculation_freezegrad

!=========================================================================================!

  subroutine calculation_increase_charge(self,dchrg)
!******************************************************************
!* increase the charge of all calculation_settings objects by one
!* or the specified dchrg
!****************************************************************** 
    implicit none
    class(calcdata) :: self
    integer,intent(in),optional :: dchrg
    integer :: i,j
    if(self%ncalculations > 0)then
     if(present(dchrg))then
       j = dchrg
     else
       j = 1
     endif
     do i=1,self%ncalculations
        self%calcs(i)%chrg = self%calcs(i)%chrg + j
     enddo
    endif
    return
  end subroutine calculation_increase_charge

!=========================================================================================!

  subroutine calculation_decrease_charge(self,dchrg)
!******************************************************************
!* decrease the charge of all calculation_settings objects by one
!* or the specified dchrg
!****************************************************************** 
    implicit none
    class(calcdata) :: self
    integer,intent(in),optional :: dchrg
    integer :: i,j
      if(self%ncalculations > 0)then
     if(present(dchrg))then
       j = dchrg
     else
       j = 1
     endif
     do i=1,self%ncalculations
        self%calcs(i)%chrg = self%calcs(i)%chrg - j
     enddo
    endif
    return
  end subroutine calculation_decrease_charge

!=========================================================================================!

  subroutine calc_set_active(self,ids)
    implicit none
    class(calcdata) :: self
    integer,intent(in) :: ids(:)
    integer :: i,j,k,l
    if (allocated(self%weightbackup)) deallocate (self%weightbackup)
    if (allocated(self%activebackup)) deallocate (self%activebackup)
!>--- on-the-fly multiscale definition
    allocate (self%weightbackup(self%ncalculations),source=1.0_wp)
    allocate (self%activebackup(self%ncalculations),source=.true.)
    do i = 1,self%ncalculations
!>--- save backup weights
      self%weightbackup(i) = self%calcs(i)%weight
      self%activebackup(i) = self%calcs(i)%active
!>--- set the weight of all unwanted calculations to 0
      if (.not.any(ids(:) .eq. i)) then
        self%calcs(i)%weight = 0.0_wp
        self%calcs(i)%active = .false.
      else
!>--- and all other to active
        if(self%calcs(i)%weight == 0.0_wp)then
           self%calcs(i)%weight = 1.0_wp
        endif 
        self%calcs(i)%active = .true.
      end if
    end do
  end subroutine calc_set_active

  subroutine calc_set_active_restore(self)
    implicit none
    class(calcdata) :: self
    integer :: i,j,k,l
    if (.not.allocated(self%weightbackup)) return
    if (.not.allocated(self%activebackup)) return
    do i = 1,self%ncalculations
!>--- set all to active and restore saved weights
      self%calcs(i)%weight = self%weightbackup(i)
      self%calcs(i)%active = self%activebackup(i)
      self%eweight(i) = self%weightbackup(i)
    end do
    deallocate (self%weightbackup)
  end subroutine calc_set_active_restore

!==========================================================================================!

  subroutine get_dipoles(calc,diptmp)
!*********************************************
!* Collect the x y and z dipole moments for
!* each calculation level in one array diptmp
!*********************************************
    implicit none
    type(calcdata),intent(inout) :: calc
    real(wp),allocatable,intent(out) :: diptmp(:,:)
    integer :: i,j,k,l,m

    m = calc%ncalculations
    allocate (diptmp(3,m),source=0.0_wp)
    do i = 1,m
      if (calc%calcs(i)%rddip) then
        diptmp(1:3,i) = calc%calcs(i)%dipole(1:3)
      end if
    end do
  end subroutine get_dipoles

!=========================================================================================!

  subroutine calculation_ONIOMexpand(self)
!*******************************************************
!* for an ONIOM calculations some of the calculators
!* have to be duplikated, which is done by this routine
!*******************************************************
    class(calcdata) :: self
    integer :: ncalcs
    integer :: maxid
    integer :: i,j,k,l,newid,j2
    type(calculation_settings) :: calculator
    integer,allocatable :: newids(:,:)
    character(len=40) :: atmp
    if (.not.allocated(self%ONIOM)) return
    ncalcs = self%ONIOM%ncalcs
    maxid = maxval(self%ONIOM%calcids(1,:),1)
    if (maxid > self%ncalculations) then
      write (stdout,'(a)') '**ERROR** in ONIOM setup: not enough calculators defined!'
      error stop
    end if

    write (stdout,'(a)',advance='no') 'Assigning and duplicating calculators for ONIOM setup ...'
    flush (stdout)

    allocate (self%ONIOMmap(ncalcs),source=0)
    allocate (newids(2,self%ONIOM%nfrag),source=0)
    k = 0
    do i = 1,self%ONIOM%nfrag

      do l = 1,2
        !> j is now the ID of the reference calculation_settings object
        j = self%ONIOM%calcids(l,i)
        if (l == 2) then
          !> to exlcude the highest ONIOM layer, we need to cycle
          j2 = self%ONIOM%calcids(1,i)
          if (j == j2) then
            newid = newids(1,i)
            newids(2,i) = newid
            self%calcs(newid)%ONIOM_highlowroot = 3
            self%calcs(newid)%ONIOM_id = i
            cycle
          end if
        end if

        if (any(self%ONIOMmap(:) .eq. j)) then
          !> If one of this type is already in the mapping, duplicate the calculator and add it
          call calculator%deallocate()
          calculator = self%calcs(j)
          !> However, we MUST not use restart I/O options for duplicates!
          call calculator%norestarts()
          call self%add(calculator)
          newid = self%ncalculations
          k = k+1
          self%ONIOMmap(k) = newid
        else
          !> otherwise (i.e. it's not yet present), we can simply add it
          k = k+1
          newid = j
          self%ONIOMmap(k) = newid
        end if
        newids(l,i) = newid

        self%calcs(newid)%ONIOM_highlowroot = l
        self%calcs(newid)%ONIOM_id = i
        select case (l)
        case (1)
          write (atmp,'(a,i0,a)') 'ONIOM.',i,'.high'
        case (2)
          write (atmp,'(a,i0,a)') 'ONIOM.',i,'.low'
        case (3)
          write (atmp,'(a,i0,a)') 'ONIOM.',i,'.root'
        end select
        self%calcs(newid)%calcspace = trim(atmp)

        !> ALWAYS set the weight of ONIOM calcs to 0!
        self%calcs(newid)%weight = 0.0_wp

        !> Check if the ONIOM fragment has a charge attached!
        if (allocated(self%ONIOM%fragment(i)%chrg)) then
          self%calcs(newid)%chrg = self%ONIOM%fragment(i)%chrg
        end if

      end do
    end do
    self%ONIOM%calcids = newids
    deallocate (newids)

    if (allocated(self%ONIOMrevmap)) deallocate (self%ONIOMrevmap)
    allocate (self%ONIOMrevmap(self%ncalculations),source=0)
    do i = 1,self%ncalculations
      do j = 1,self%ONIOM%ncalcs
        if (self%ONIOMmap(j) == i) then
          self%ONIOMrevmap(i) = j
        end if
      end do
    end do

    write (stdout,*) 'done.'
  end subroutine calculation_ONIOMexpand

!========================================================================================!

  subroutine calculation_info(self,iunit)
    implicit none
    class(calcdata) :: self
    integer,intent(in) :: iunit
    integer :: i,j
    character(len=*),parameter :: fmt1 = '(1x,a20," : ",i5)'
    character(len=*),parameter :: fmt2 = '(1x,a20," : ",f12.5)'
    character(len=20) :: atmp
    integer :: constraintype(8)

    write (iunit,'(1x,a)') '----------------'
    write (iunit,'(1x,a)') 'Calculation info'
    write (iunit,'(1x,a)') '----------------'
    if (self%ncalculations <= 0) then
      write (iunit,'("> ",a)') 'No calculation levels set up!'
    else if (self%ncalculations > 1) then
      do i = 1,self%ncalculations
        if (self%calcs(i)%ONIOM_id > 0) then
          write (iunit,'("> ",a,i0,a)') 'Automatic ONIOM setup calculation level ',i,':'
        else
          write (iunit,'("> ",a,i0,a)') 'User-defined calculation level ',i,':'
        end if
        call self%calcs(i)%info(iunit)
      end do
    else
      write (iunit,'("> ",a)') 'User-defined calculation level:'
      call self%calcs(1)%info(iunit)
    end if
    write (iunit,*)

    if (self%nconstraints > 0) then
      write (iunit,'("> ",a)') 'User-defined constraints:'
      if (self%nconstraints <= 20) then
        do i = 1,self%nconstraints
          call self%cons(i)%print(iunit)
        end do
      else
        constraintype(:) = 0
        do i = 1,self%nconstraints
          j = self%cons(i)%type
          if (j > 0.and.j < 9) then
            constraintype(j) = constraintype(j)+1
          end if
        end do
        if (constraintype(1) > 0) write (iunit,'(2x,a,i0)') '# bond constraints    : ',constraintype(1)
        if (constraintype(2) > 0) write (iunit,'(2x,a,i0)') '# angle constraints   : ',constraintype(2)
        if (constraintype(3) > 0) write (iunit,'(2x,a,i0)') '# dihedral constraints: ',constraintype(3)
        if (constraintype(4) > 0) write (iunit,'(2x,a,i0)') '# wall potential      : ',constraintype(4)
        if (constraintype(5) > 0) write (iunit,'(2x,a,i0)') '# wall(logfermi) potential  : ',constraintype(5)
        if (constraintype(8) > 0) write (iunit,'(2x,a,i0)') '# bondrange constraints    : ',constraintype(8)
      end if
      write (iunit,*)
    end if

    if (self%ncalculations > 1) then
      write (iunit,'("> ",a)') 'Energy and gradient processing:'
      select case (self%id)
      case (1:)
        write (iunit,'(1x,a,i0)') 'Energies and gradients of calculation level ',self%id, &
        &  ' will be used'
      case (-1)
        write (iunit,'(1x,a)') 'Special MECP energy and gradients will be constructed'
        write (iunit,'(1x,a)') 'See https://doi.org/10.1021/acs.jpclett.3c00494'
      case default
        write (iunit,'(1x,a)') 'Energies and gradients of all calculation levels will be'// &
        & ' added according to their weights'
      end select
      if (allocated(self%ONIOM)) then
        write (iunit,'(1x,a)') 'ONIOM energy and gradient will be constructed from calculations:'
        do i = 1,self%ncalculations
          if (any(self%ONIOMmap(:) .eq. i)) then
            write (stdout,'(1x,i0)',advance='no') i
          end if
        end do
        write (iunit,*)
      end if
      write (iunit,*)
    end if

    return
  end subroutine calculation_info



!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!> CALCULATION_SETTINGS associated routines
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!

  subroutine calculation_settings_deallocate(self)
    implicit none
    class(calculation_settings) :: self

    if (allocated(self%calcspace)) deallocate (self%calcspace)
    if (allocated(self%calcfile)) deallocate (self%calcfile)
    if (allocated(self%gradfile)) deallocate (self%gradfile)
    if (allocated(self%path)) deallocate (self%path)
    if (allocated(self%other)) deallocate (self%other)
    if (allocated(self%binary)) deallocate (self%binary)
    if (allocated(self%systemcall)) deallocate (self%systemcall)
    if (allocated(self%description)) deallocate (self%description)
    if (allocated(self%wbo)) deallocate (self%wbo)
    if (allocated(self%dipgrad)) deallocate (self%dipgrad)
    if (allocated(self%gradkey)) deallocate (self%gradkey)
    if (allocated(self%efile)) deallocate (self%efile)
    if (allocated(self%solvmodel)) deallocate (self%solvmodel)
    if (allocated(self%solvent)) deallocate (self%solvent)
    if (allocated(self%tblite)) deallocate (self%tblite)
    if (allocated(self%g0calc)) deallocate (self%g0calc)
    if (allocated(self%ff_dat)) deallocate (self%ff_dat)
    if (allocated(self%libpvol)) deallocate (self%libpvol)

    self%id = 0
    self%prch = stdout
    self%chrg = 0
    self%uhf = 0
    self%epot = 0.0_wp
    self%efix = 0.0_wp
    self%etot = 0.0_wp

    self%rdwbo = .false.
    self%rddip = .false.
    self%dipole = 0.0_wp
    self%rddipgrad = .false.
    self%gradtype = 0
    self%gradfmt = 0

    self%tblitelvl = 2
    self%etemp = 300.0_wp
    self%accuracy = 1.0_wp
    self%apiclean = .false.
    self%maxscc = 500
    self%saveint = .false.

    self%ngrid = 230
    self%extpressure = 0.0_wp
    self%proberad = 1.5_wp

    self%ONIOM_highlowroot = 0
    self%ONIOM_id = 0
    return
  end subroutine calculation_settings_deallocate

!=========================================================================================!

  subroutine calculation_settings_addconfig(self,config)
    implicit none
    class(calculation_settings) :: self
    integer,intent(in)  :: config(:)
    integer :: i,j
    integer :: l,lold,lnew,n,nnew
    integer,allocatable :: configtmp(:,:)

    l = size(config,1)
    if (allocated(self%config)) deallocate (self%config)
    allocate (self%config(l))
    self%config = config
  end subroutine calculation_settings_addconfig

!=========================================================================================!

  subroutine calculation_settings_norestarts(self)
!*************************************************
!* remove restart options from this calculation
!*************************************************
    implicit none
    class(calculation_settings) :: self
    self%restart = .false.
    if (allocated(self%refgeo)) deallocate (self%refgeo)
    if (allocated(self%restartfile)) deallocate (self%restartfile)
  end subroutine calculation_settings_norestarts

!=========================================================================================!

!>-- check for missing settings in a calculation_settings object
  subroutine calculation_settings_autocomplete(self,id)
    use iomod,only:crest_scratch_dir
    implicit none
    class(calculation_settings) :: self
    integer,intent(in)  :: id
    integer :: i,j
    character(len=50) :: nmbr
    character(len=:),allocatable :: scratch
    logical,save :: announced = .false.

    !> add a short description
    self%description = trim(jobdescription(self%id+1))
    call calculation_settings_shortflag(self)

    if (.not.allocated(self%calcspace)) then
      !> I've decided to perform all calculations in a separate directory to
      !> avoid accumulation of files in the main workspace
      !
      !> LOCAL EXTENSION: that directory can be moved to node-local storage with
      !> CREST_SCRATCH_DIR (see iomod%crest_scratch_dir).  An explicit
      !> "calcspace = ..." in the TOML always wins.
      !> Motivation: on an NFS / spinning-disk workspace every scratch write of
      !> the external program is a network round trip, which for a cheap level
      !> of theory easily costs more than the calculation itself.
      call crest_scratch_dir(scratch)
      write (nmbr,'(i0)') id
      self%calcspace = trim(scratch)//'calculation.level.'//trim(nmbr)
      if (len(scratch) > 0 .and. .not.announced) then
        announced = .true.
        write (stdout,'(/,1x,a,a)') &
        &  'Calculation scratch directories relocated to ',trim(scratch)
        write (stdout,'(1x,a)') &
        &  '(CREST_SCRATCH_DIR - the ORCA.out / ORCA.err of the last call are there)'
      end if
    end if

    if (self%pr .and. self%prch.ne.stdout) then
      self%prch = self%prch+id
    end if
  end subroutine calculation_settings_autocomplete

!>--- create a short calculation info flag
    subroutine calculation_settings_shortflag(self)
    implicit none
    class(calculation_settings) :: self
    integer :: i,j

    select case( self%id )
    case( jobtype%xtbsys )   
      self%shortflag = 'xtb subprocess'
    case( jobtype%generic ) 
      self%shortflag = 'generic subprocess' 
    case( jobtype%turbomole ) 
      self%shortflag = 'TURBOMOLE subprocess'
    case( jobtype%orca ) 
      self%shortflag = 'ORCA subprocess'
    case( jobtype%terachem ) 
      self%shortflag = 'TeraChem subprocess'
    case( jobtype%tblite )
      select case (self%tblitelvl)
      case (xtblvl%gfn2)
        self%shortflag = 'GFN2-xTB'
      case (xtblvl%gfn1)
        self%shortflag = 'GFN1-xTB'
      case (xtblvl%ipea1)
        self%shortflag = 'IPEA1-xTB'
      case (xtblvl%ceh)
        self%shortflag = 'CEH'
      case (xtblvl%eeq)
        self%shortflag = 'EEQ(D4)'
      case (xtblvl%param)
        self%shortflag = 'parameter file: '//trim(self%tbliteparam)
      end select
    case( jobtype%gfn0 )
      self%shortflag =  'GFN0-xTB' 
    case( jobtype%gfn0occ )
      self%shortflag =  'GFN0-xTB*'
    case( jobtype%gfnff )
      self%shortflag =  'GFN-FF'
    case( jobtype%libpvol )
      self%shortflag =  'LIVPVOL'
    case( jobtype%lj )
      self%shortflag =  'LJ'
    case default
      self%shortflag = 'undefined'
    end select
    if(allocated(self%solvmodel).and.allocated(self%solvent))then
      self%shortflag = self%shortflag//'/'//trim(self%solvmodel)
       self%shortflag = self%shortflag//'('//trim(self%solvent)//')'
    endif
  end subroutine calculation_settings_shortflag

!>-- generate a unique print id for the calculation
  subroutine calculation_settings_printid(self,thread,id)
    implicit none
    class(calculation_settings) :: self
    integer,intent(in)  :: thread,id
    integer :: i,j,dum
    character(len=50) :: nmbr
    dum = 100*(thread+1)
    dum = dum+id
    self%prch = dum
  end subroutine calculation_settings_printid

  subroutine calculation_dump_dipgrad(self,filename)
    implicit none
    class(calculation_settings) :: self
    character(len=*),intent(in) :: filename
    integer :: i,j,ich
    if (.not.allocated(self%dipgrad)) return
    open (newunit=ich,file=filename)
    do i = 1,size(self%dipgrad,2)
      write (ich,'(3f20.10)') self%dipgrad(1:3,i)
    end do
    close (ich)
  end subroutine calculation_dump_dipgrad

!=========================================================================================!

  subroutine calculation_settings_info(self,iunit)
    implicit none
    class(calculation_settings) :: self
    integer,intent(in) :: iunit
    integer :: i,j
    character(len=*),parameter :: fmt1 = '(" :",2x,a20," : ",i0)'
    character(len=*),parameter :: fmt2 = '(" :",2x,a20," : ",f0.5)'
    character(len=*),parameter :: fmt3 = '(" :",2x,a20," : ",a)'
    character(len=*),parameter :: fmt4 = '(" :",1x,a)'
    character(len=20) :: atmp
    logical :: gxtbwarn

    gxtbwarn=.false.

    if (allocated(self%description)) then
      write (iunit,'(" :",1x,a)') trim(self%description)
    else
      write (atmp,*) 'Job type'
      write (iunit,fmt1) atmp,self%id
    end if
    !> more info
    if (self%id == jobtype%tblite) then
      select case (self%tblitelvl)
      case (xtblvl%gfn2)
        write (iunit,fmt4) 'GFN2-xTB level'
      case (xtblvl%gfn1)
        write (iunit,fmt4) 'GFN1-xTB level'
      case (xtblvl%ceh)
        write (iunit,fmt4) 'Charge Extended Hückel (CEH) model'
      end select
    end if
    if (any((/jobtype%orca,jobtype%xtbsys,jobtype%turbomole, &
    &  jobtype%generic,jobtype%terachem/) == self%id)) then
!>--- LOCAL CHANGE (2026-09-30): g-xTB is now detected by its --gxtb flag, not
!> by a binary called 'gxtb'.  The method lives inside the xtb binary as a
!> tblite parameter set, so the old name test no longer fires and the level
!> would be reported as plain "selected binary : xtb", which is misleading in
!> the header of a production run.  Keep accepting the old spelling too
!> (binary containing 'gxtb') so a legacy setup still gets labelled.
      if(index(self%binary,'gxtb').ne.0)then
        write(iunit,fmt4) 'g-xTB method (via xtb subprocess)'
        gxtbwarn = .true.
      else if(allocated(self%other))then
        if(index(self%other,'--gxtb').ne.0)then
          write(iunit,fmt4) 'g-xTB method (via xtb subprocess)'
          gxtbwarn = .true.
        else
          write (iunit,'(" :",3x,a,a)') 'selected binary : ',trim(self%binary)
        end if
      else  
        write (iunit,'(" :",3x,a,a)') 'selected binary : ',trim(self%binary)
      endif
    end if
    if (self%refine_lvl > 0) then
      write (atmp,*) 'refinement stage'
      write (iunit,fmt1) atmp,self%refine_lvl
    end if

    !> system data
    write (atmp,*) 'Molecular charge'
    write (iunit,fmt1) atmp,self%chrg
    if (self%uhf /= 0) then
      write (atmp,*) 'UHF parameter'
      write (iunit,fmt1) atmp,self%uhf
    end if

    if (allocated(self%solvmodel)) then
      write (atmp,*) 'Solvation model'
      write (iunit,fmt3) atmp,trim(self%solvmodel)
    end if
    if (allocated(self%solvent)) then
    write (atmp,*) 'Solvent'
      write (iunit,fmt3) atmp,trim(self%solvent)
    end if

    !> xTB specific parameters
    if (any((/jobtype%tblite,jobtype%xtbsys,jobtype%gfn0,jobtype%gfn0occ/) == self%id)) then
      write (atmp,*) 'Fermi temperature'
      write (iunit,fmt2) atmp,self%etemp
      write (atmp,*) 'Accuracy'
      write (iunit,fmt2) atmp,self%accuracy
      write (atmp,*) 'max SCC cycles'
      write (iunit,fmt1) atmp,self%maxscc
    end if

    write (atmp,*) 'Reset data?'
    if (self%apiclean) write (iunit,fmt3) atmp,'yes'
    write (atmp,*) 'Read WBOs?'
    if (self%rdwbo) write (iunit,fmt3) atmp,'yes'
    write (atmp,*) 'Read dipoles?'
    if (self%rddip) write (iunit,fmt3) atmp,'yes'

    if (self%ONIOM_highlowroot /= 0) then
      select case (self%ONIOM_highlowroot)
      case (1)
        write (atmp,*) 'ONIOM frag ("high")'
      case (2)
        write (atmp,*) 'ONIOM frag ("low")'
      case (3)
        write (atmp,*) 'ONIOM frag ("root")'
      end select
      write (iunit,fmt1) trim(atmp),self%ONIOM_id
    else
      if(self%weight .ne. 1.0_wp)then
        write (atmp,*) 'Weight'
        write (iunit,fmt2) atmp,self%weight
      endif
    end if

    if(gxtbwarn)then
!>--- LOCAL CHANGE (2026-09-30): the old text claimed "Gradients are NUMERICAL".
!> That was true for the FIRST g-xTB release (standalone binary, Turbomole
!> protocol, -grad = finite differences).  In the current xtb 6.7.1 based
!> package g-xTB is a tblite parameter set and the gradient is ANALYTICAL,
!> arriving in the normal .engrad file.  The remaining caveats are the
!> genuinely current ones.
       write(iunit,fmt4) 'NOTE  : g-xTB is a preliminary parametrisation (Z = 1-103).'
       write(iunit,fmt4) 'NOTE  : Solvation (--gbe/--cosmo) is preliminary and may be'
       write(iunit,fmt4) 'NOTE  :   unstable during optimisation.'
    endif

  end subroutine calculation_settings_info

!=========================================================================================!

  subroutine create_calclevel_shortcut(self,levelstring)
!*********************************************************************
!* subroutine create_calclevel_shortcut called with %create(...)
!* Set up a calculation_settings object for a given level of theory
!* More shortcuts can be added as required.
!* Be careful about the intent(out) setting!
!*********************************************************************
    implicit none
    class(calculation_settings),intent(out) :: self
    character(len=*) :: levelstring
    call self%deallocate()
    select case (trim(levelstring))
    case ('gfnff','--gff','--gfnff')
      self%id = jobtype%gfnff
    case ('gfn0','--gfn0')
      self%id = jobtype%gfn0
!>--- LOCAL CHANGE (2026-10-12): GFN2 now runs through the EXTERNAL xtb binary,
!> exactly like g-xTB further below, instead of the in-process tblite library.
!>
!> WHY: tblite opens its own OpenMP team INSIDE the CREST process rather than
!> living within the thread budget CREST sizes for an external child.  With many
!> workers running concurrently this oversubscribes the node, and the MD/MTD
!> phase was observed to end up more than TEN TIMES slower than the very same
!> setup driven through the external binary.  On the subprocess channel the child
!> just inherits OMP_NUM_THREADS like every other external calculator, so the
!> core allocation stays under CREST's control (which is the whole point of
!> CREST_CORES_PER_JOB / CREST_MAX_CORES_PER_JOB).
!>
!> This build is linked WITHOUT tblite (WITH_TBLITE=OFF), so this branch must not
!> point at jobtype%tblite any more - it would abort with "Compiled without
!> tblite support!".  NOTE: '-gfn1' is deliberately NOT rerouted and will abort
!> in this build; only GFN2 (the level actually used) is supported.
!> Solvation: '-alpb <solv>' sets cal%solvmodel/cal%solvent, which ONLY the
!> in-process path ever read.  legacy_wrappers.f90 now translates them into real
!> xtb CLI flags for this channel - see the note there.
    case ('gfn2','--gfn2')
      self%id = jobtype%xtbsys
      self%binary = 'xtb'
      self%other = '--gfn2'
    case ('gfn1','--gfn1')
      self%id = jobtype%tblite
      self%tblitelvl = 1
    case ('gp3')
      self%id = jobtype%turbomole
      self%rdgrad = .false.
      self%binary = 'gp3'
!>--- LOCAL CHANGE (2026-09-30): g-xTB runs through the xtb SUBPROCESS, not
!> through the TURBOMOLE channel.
!> Upstream wires 'gxtb'/'gxtb_dev' to jobtype%turbomole because the FIRST
!> g-xTB release was a standalone binary that only spoke Turbomole's
!> coord/energy/gradient file protocol and had NO analytical gradient (hence
!> the -grad hack and the "numeric gradients" warning everywhere).
!> The current g-xTB is shipped as a full xtb 6.7.1 in which the method is a
!> tblite parameter set (`--gxtb` == "use the g-xTB method via the tblite
!> library", cf. the binary's own help text and the GXTB_H0SPEC /
!> GXTB_REPULSION / GXTB_MULTIPOLE record types in libxtb).  tblite hands
!> back an ANALYTICAL gradient in the ordinary xtb .engrad file, so the whole
!> Turbomole detour is unnecessary - and actively harmful, because the
!> Turbomole path would drive a binary called 'gxtb' with Turbomole-style
!> files instead of the xtb CLI it actually implements.
!> Consequences of moving to xtbsys:
!>   * binary is plain 'xtb' (name resolved via PATH), not 'gxtb'
!>   * '--gxtb' is injected through calc%other
!>   * the gradient is read by rd_grad_engrad (xtbinp.engrad), not rd_grad_tm
!>   * WBOs become available again (rdwbo stays as-is; the tblite path can
!>     write the wbo file, unlike turbom_sc which had it commented out)
!>   * no .CHRG/.UHF side files; charge/uhf go in as --chrg/--uhf
!>
!> LOCAL EXTENSION (2026-09-30) -- g-xTB + implicit solvation
!> ---------------------------------------------------------
!> g-xTB does not use the GFN2 solvation keywords.  Its own package documents
!> two models, selected on the xtb command line:
!>     --gbe   <solvent>   generalized Born, finite epsilon
!>     --cosmo <solvent>   domain-decomposition conductor-like screening
!> Neither is reachable from CREST's CLI: '-g'/'-alpb' write env%gbsa /
!> env%solvent, which legacy_wrappers.f90 turns into cal%solvmodel/cal%solvent
!> -- fields that ONLY the in-process tblite path ever reads
!> (api_engrad.f90 -> tblite_add_solv).  The subprocess path builds its command
!> line from calc%other alone (xtb_sc.f90), so on the xtb channel a '-alpb thf'
!> used to be dropped SILENTLY.
!> So the flag is appended to calc%other here, from a plain environment
!> variable - no recompilation per solvent, and nothing to change at the call
!> sites.  Examples:
!>     export CREST_GXTB_SOLV='--gbe thf'
!>     export CREST_GXTB_SOLV='--cosmo water'
!> Unset/empty => the old, solvent-free behaviour (and the same default for
!> every other level of theory).
!> Whitespace is collapsed so a stray double space cannot reach the shell, and
!> the value is passed through verbatim otherwise -- this is a raw passthrough
!> for an external program, so no attempt is made to validate the solvent name
!> (xtb itself will report an unknown solvent).
!>
!> CAVEAT, read before using it: upstream's g-xTB README calls solvation
!> preliminary.  'gbe' is reported unstable during geometry optimisation and
!> 'cosmo' omits the solute response to the solvent polarisation, i.e. its
!> GRADIENT IS NOT CONSISTENT.  Since CREST optimises and analyses Hessians,
!> treat any number obtained this way as exploratory only.
    case ('gxtb','gxtb_dev')
      self%id = jobtype%xtbsys
      self%binary = 'xtb'
      self%other = '--gxtb'
!>--- LOCAL EXTENSION: --gbe/--cosmo passthrough (see the note above)
      block
        character(len=128) :: sv
        character(len=:),allocatable :: clean
        integer :: k
        sv = ' '
        call getenv('CREST_GXTB_SOLV',sv)
        if (len_trim(sv) > 0) then
!>--- collapse runs of blanks to one, so '--gbe   thf' cannot smuggle an
!>    empty argument onto the shell command line
          clean = ''
          do k = 1,len_trim(sv)
            if (sv(k:k) == ' ') then
              if (len(clean) > 0.and.clean(len(clean):len(clean)) /= ' ') clean = clean//' '
            else
              clean = clean//sv(k:k)
            end if
          end do
          if (len(clean) > 0) self%other = trim(self%other)//' '//trim(clean)
        end if
      end block
      self%rdgrad = .true.
    case ('orca')
      self%id = jobtype%orca

    case ('generic')
      self%id = jobtype%generic

    end select
    call self%autocomplete(self%id)
  end subroutine create_calclevel_shortcut

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!
end module calc_type
