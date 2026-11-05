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

subroutine crest_refine(env,input,output)
!*******************************************************************
!* subroutine crest_refine
!* This subroutine will process the ensemble specified by "input"
!* and either overwrite it or write a new file "output".
!* The routine is intended to be called after geometry optimization
!* to re-rank structures by singlepoint energies, or add other
!* contributions to the energy.
!*******************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use cregen_interface
  use parallel_interface
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: input
  character(len=*),intent(in),optional :: output
!===========================================================!
  integer :: i,j,k,l,io,ich,m,t1,t2
  logical :: pr,wr,ex
  logical :: topopresort  !> LOCAL EXTENSION: topology-only pre-filter before the refinement
!===========================================================!
  character(len=:),allocatable :: outname
  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:),etmp(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  integer :: nrefine,refine_stage
!===========================================================!
!>--- setup
  if (present(output)) then
    outname = output !> new file
  else
    outname = input  !> overwrite
  end if
 
!>--- presorting step, if necessary
  if(env%refine_presort)then
    call newcregen(env,0,input)
    call rename('crest_ensemble.xyz',input)
  endif

!>--- read in
  call rdensemble(input,nat,nall,at,xyz,eread)
  allocate (etmp(nall),source=0.0_wp)
!>--- LOCAL EXTENSION: topology-only pre-filter.  Has to run HERE, i.e. before
!> the conversion below: cregen_topocheck() divides by bohr itself, so it wants
!> exactly the Angstrom coordinates rdensemble() has just produced.
  topopresort = want_topo_presort()
  if (topopresort) call crest_refine_topopresort
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_sploop requires coordinates in Bohrs
    xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!===========================================================!
  DO_REFINE: if (allocated(env%refine_queue)) then
!===========================================================!

    call smallhead('ensemble refinement')

    nrefine = size(env%refine_queue,1)

    do i = 1,nrefine
      refine_stage = env%refine_queue(i)
      !> set the calculator to the correct stage
      env%calc%refine_stage = refine_stage

      select case (refine_stage)
      case (refine%singlepoint)
        write (stdout,'("> Singlepoint re-ranking for ",i0," structures")') nall
        call crest_sploop(env,nat,nall,at,xyz,eread)

      case (refine%correction)
        write (stdout,'("> Additive correction for ",i0," structures")') nall
        call crest_sploop(env,nat,nall,at,xyz,etmp)
        eread(:) = eread(:)+etmp(:)

      case (refine%geoopt)
        write (stdout,'("> Geometry optimization of ",i0," structures")') nall
        call crest_oloop(env,nat,nall,at,xyz,eread,.false.)

      case(refine%confsolv)
        call new_ompautoset(env,'subprocess',1,t1,t2)
        write (stdout,'("> ConfSolv: ΔΔGsoln estimation from 3D directed message passing neural networks (D-MPNN)")')
        call confsolv_request( input, nall, t2, etmp, io)
        if(io == 0)then
        eread(:) = etmp(:)*kcaltoau  !> since CREGEN deals with Eh energies
        endif   
      end select
      write(stdout,*) 
    end do

    !> reset the refinement stage of the calculator
    env%calc%refine_stage = 0

!===========================================================!
  end if DO_REFINE
!===========================================================!

!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: ensemble file must be written in AA
  xyz = xyz / angstrom
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- write output ensemble
  call wrensemble(outname,nat,nall,at,xyz,eread)

!===========================================================!
  deallocate (etmp,eread,xyz,at)
  return
!========================================================================================!
!========================================================================================!
contains
!========================================================================================!
!> LOCAL EXTENSION: should the topology-only pre-filter run?
!>
!> Default is ON (classes.f90: refine_presort_topo = .true.), but it is only
!> ever *effective* when there actually is a refinement stage - otherwise
!> crest_refine() is just rewriting the ensemble and there is nothing to save.
!> env%checktopo is respected, so "-notopo" / TOML notopo = true switches this
!> off as well.
!========================================================================================!
  logical function want_topo_presort() result(yes)
    implicit none
    character(len=64) :: envstr
    integer :: ival,stat,rd
    logical :: tp
    tp = env%refine_presort_topo
!>--- CREST_REFINE_PRESORT_TOPO=0 / =1 overrides the TOML value (and the default)
    call get_environment_variable('CREST_REFINE_PRESORT_TOPO',value=envstr,status=stat)
    if (stat == 0) then
      if (len_trim(envstr) > 0) then
        read (envstr,*,iostat=rd) ival
        if (rd == 0) tp = (ival > 0)
      end if
    end if
    yes = allocated(env%refine_queue).and.tp.and.env%checktopo
    return
  end function want_topo_presort
!========================================================================================!
!> LOCAL EXTENSION: the pre-filter itself.
!>
!> The very same check runs seconds later in sort_and_check -> newcregen ->
!> cregen_topocheck, and everything removed here is removed there too.  Doing it
!> here only means those structures never get an expensive singlepoint.
!> Deliberately NOT done here: energy sorting, the EWIN cut and the RMSD
!> de-duplication - all three depend on the *refinement* energies, which do not
!> exist yet at this point.  That is exactly why the ready-made
!> refine_presort (which runs the FULL cregen) is unsafe to enable by default
!> and this one is safe.
!========================================================================================!
  subroutine crest_refine_topopresort
!>--- LOCAL FIX (2026-10-05): REQUIRED - this routine is the one caller of
!> cregen_topocheck() that runs BEFORE newcregen(), and newcregen() is the only
!> place that arms the reject audit (exclfile_check + exclfile_init).  Without
!> these two calls every TOPO rejection produced here was silently dropped by
!> exclfile_write() (`active` still .false. / `uch` still -1), which is how the
!> log could say "82 of 87 structures survive the topology pre-filter" while
!> crest_excluded_structures.xyz stayed empty: the five dropped structures never
!> reach CREGEN, so its own topocheck never gets a chance to audit them either.
!> Both calls are idempotent - the one in newcregen() stays a no-op, and the
!> file is still truncated exactly once per run.
    use exclfile_module,only:exclfile_check,exclfile_init
    implicit none
    integer :: i,nallnew
    character(len=128),allocatable :: comments(:)
    real(wp),allocatable :: xyzref(:,:,:)
    real(wp),allocatable :: ereadref(:)

    write (stdout,'(1x,a)') '> Topology pre-filter before refinement'
!>--- cregen_topocheck() reorders xyz AND comments with the same permutation,
!> but it knows nothing about the real-valued energy array.  Carry the energies
!> inside the comment lines - that is how CREGEN itself transports them - and
!> read them back with grepenergy() afterwards.
    allocate (comments(nall))
    do i = 1,nall
      write (comments(i),'(f18.8)') eread(i)
    end do
    call exclfile_check()
    call exclfile_init()
    call cregen_topocheck(stdout,env,env%checkiso,nat,nall,at,xyz,comments,nallnew)

    if (nallnew < 1) then
      write (stdout,'(1x,a)') '  **WARNING** no structure passed the topology check - keeping all of them'
      deallocate (comments)
      return
    end if

    if (nallnew < nall) then
      allocate (xyzref(3,nat,nallnew),ereadref(nallnew))
      xyzref(1:3,1:nat,1:nallnew) = xyz(1:3,1:nat,1:nallnew)
      do i = 1,nallnew
        ereadref(i) = grepenergy(comments(i))
      end do
      call move_alloc(xyzref,xyz)
      call move_alloc(ereadref,eread)
      write (stdout,'(1x,a,i0,a,i0,a)') '> ',nallnew,' of ',nall, &
      &   ' structures survive the topology pre-filter'
      nall = nallnew
    else
      write (stdout,'(1x,a)') '> all structures survive the topology pre-filter'
    end if
    deallocate (comments)
    return
  end subroutine crest_refine_topopresort
!========================================================================================!
end subroutine crest_refine
!========================================================================================!
!========================================================================================!
