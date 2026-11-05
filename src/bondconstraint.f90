!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!=========================================================================================!
!> Get the bonds from the topology and write a constraints
!> file with a bond length constrain on each bond.
!>
!> A file "bondlengths" is written with a constraint on each bond.
!> The input argument is the required force constant.
!==========================================================================!
subroutine autoBondConstraint(filename,forceconstant,wbofile)
  use crest_parameters,only:wp
  use zdata
  implicit none
  character(len=*) :: filename
  character(len=*) :: wbofile
  real(wp) :: forceconstant
  type(zmolecule) :: zmol
  !--- get topology
  call simpletopo_file(filename,zmol,.false.,.false.,wbofile)
  !--- get bond matrix and write "bondlengths" file
  call getbmat(zmol,1,forceconstant)
  call zmol%deallocate()
  return
end subroutine autoBondConstraint

subroutine autoBondConstraint_withEZ(filename,forceconstant,wbofile)
  use crest_parameters,only:wp
  use zdata
  implicit none
  character(len=*) :: filename
  character(len=*) :: wbofile
  real(wp) :: forceconstant
  type(zmolecule) :: zmol
  !--- get topology
  call simpletopo_file(filename,zmol,.false.,.false.,wbofile)
  !--- get bond matrix and write "bondlengths" file
  call getbmat(zmol,5,forceconstant)
  call zmol%deallocate()
  return
end subroutine autoBondConstraint_withEZ

!=========================================================================================!
!> Get the bonds from the topology and write a constraints
!> file with a bond length constrain on each bond, but only if it is
!> an M-X bond (M = transition metal atom)
!>
!> A file "bondlengths" is written with a constraint on each bond.
!> The input argument is the required force constant.
!==========================================================================!
subroutine autoMetalConstraint(filename,forceconstant,wbofile)
  use crest_parameters,only:wp
  use zdata
  implicit none
  character(len=*) :: filename
  character(len=*) :: wbofile
  real(wp) :: forceconstant
  type(zmolecule) :: zmol
  !--- get topology
  call simpletopo_file(filename,zmol,.false.,.false.,wbofile)
  !--- get bond matrix and write "bondlengths" file
  call getbmat(zmol,2,forceconstant)
  call zmol%deallocate()
  return
end subroutine autoMetalConstraint

!=========================================================================================!
!> Get the bonds from the topology and write a constraints
!> file with a bond length constrain on each bond, except if the bond
!> is an X-H bond
!>
!> A file "bondlengths" is written with a constraint on each bond.
!> The input argument is the required force constant.
!==========================================================================!
subroutine autoHeavyConstraint(filename,forceconstant)
  use crest_parameters,only:wp
  use zdata
  implicit none
  character(len=*) :: filename
  real(wp) :: forceconstant
  type(zmolecule) :: zmol
  !--- get topology
  call simpletopo_file(filename,zmol,.false.,.false.,'')
  !--- get bond matrix and write "bondlengths" file
  call getbmat(zmol,3,forceconstant)
  call zmol%deallocate()
  return
end subroutine autoHeavyConstraint

!=========================================================================================!
!> Get the bonds from the topology and write a constraints
!> file with a bond length constrain on each bond, but only if the bond
!> is an X-H bond
!>
!> A file "bondlengths" is written with a constraint on each bond.
!> The input argument is the required force constant.
!==========================================================================!
subroutine autoHydrogenConstraint(filename,forceconstant)
  use crest_parameters,only:wp
  use zdata
  implicit none
  character(len=*) :: filename
  real(wp) :: forceconstant
  type(zmolecule) :: zmol
  !--- get topology
  call simpletopo_file(filename,zmol,.false.,.false.,'')
  !--- get bond matrix and write "bondlengths" file
  call getbmat(zmol,4,forceconstant)
  call zmol%deallocate()
  return
end subroutine autoHydrogenConstraint

!=========================================================================================!
!> Some routines related to the Bond matrix (BMAT)
!===========================================================================!
subroutine getbmat(zmol,r,force)
  use crest_parameters,only:wp
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer :: r
  integer :: nat
  real(wp) :: force
  real(wp),allocatable :: xyz(:,:)
  real(wp),allocatable :: bmat(:,:)
  integer,allocatable :: at(:)
  integer,allocatable :: bonds(:,:)
  integer :: i,j,nb

  integer,parameter :: cbonds = 1
  integer,parameter :: cmetal = 2
  integer,parameter :: cheavy = 3
  integer,parameter :: chydro = 4
  integer,parameter :: cistrans = 5
  integer,parameter :: countbonds = 6

  real(wp),parameter :: bohr = 0.52917726_wp

  nat = zmol%nat
  allocate (at(nat),source=0)
  at = zmol%at
  allocate (xyz(3,nat),bmat(nat,nat),source=0.0_wp)

  do i = 1,nat
    xyz(:,i) = zmol%zat(i)%cart(:)
  end do

  nb = 0
  do i = 1,nat
    do j = 1,nat
      if (any(zmol%zat(i)%ngh(:) .eq. j)) then  !only include bonds from the neighbour lists
        bmat(i,j) = zmol%distmat(i,j)
        bmat(i,j) = bmat(i,j)*bohr   !BMAT is in Angstroem!
        nb = nb+1
      else
        cycle
      end if
    end do
  end do

  !--- write the constrain file
  if (r == cbonds) then
    call writeBmatconstr(nat,bmat,force)
  end if
  if (r == cmetal) then
    call writeMetalconstr(nat,at,bmat,force)
  end if
  if (r == cheavy) then
    call writeHeavyconstr(nat,at,bmat,force)
  end if
  if (r == chydro) then
    call writeHydrogenconstr(nat,at,bmat,force)
  end if
  if (r == cistrans) then
    call writeBmatconstr_withEZ(zmol,nat,bmat,force)
  end if

  !--- utility
  if (r == countbonds) then
    allocate (bonds(2,nb),source=0)
    nb = 0
    do i = 1,nat
      do j = 1,i
        if (bmat(i,j) > 1d-6) then
          nb = nb+1
          bonds(1,nb) = j
          bonds(2,nb) = i
        end if
      end do
    end do
    deallocate (bonds)
  end if

  deallocate (bmat,xyz,at)
  return
end subroutine getbmat

!====================================================!
!> Write a constraint file with all bonds
!====================================================!
subroutine writeBmatconstr(nat,bmat,force)
  use crest_parameters,only:wp
  implicit none
  integer :: nat
  real(wp) :: bmat(nat,nat)
  integer :: i,j,ich
  real(wp) :: force
  character(len=20) :: dumm

  open (newunit=ich,file='bondlengths')
  write (ich,'(a)') '$constrain'
  write (dumm,'(f16.4)') force
  write (ich,'(3x,a,a)') 'force constant=',adjustl(trim(dumm))
  do i = 1,nat
    do j = i+1,nat
      if (bmat(i,j) .gt. 0.1_wp) then
        write (ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',i,',',j,',',bmat(i,j)
      end if
    end do
  end do
  write (ich,'(a)') '$end'
  close (ich)
  return
end subroutine writeBmatconstr

!====================================================!
!> Write a constraint file with all TM bonds
!====================================================!
subroutine writeMetalconstr(nat,at,bmat,force)
  use crest_parameters,only:wp
  implicit none
  integer :: nat
  integer :: at(nat)
  real(wp) :: bmat(nat,nat)
  integer :: i,j,ich
  real(wp) :: force
  character(len=20) :: dumm
  logical :: isTMetal  !this is a function

  open (newunit=ich,file='bondlengths')
  write (ich,'(a)') '$constrain'
  write (dumm,'(f16.4)') force
  write (ich,'(3x,a,a)') 'force constant=',adjustl(trim(dumm))
  do i = 1,nat
    do j = i+1,nat
      if (bmat(i,j) .gt. 0.1_wp) then
        if (isTMetal(at(i)).or.isTMetal(at(j))) then
          write (ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',i,',',j,',',bmat(i,j)
        end if
      end if
    end do
  end do
  write (ich,'(a)') '$end'
  close (ich)
  return
end subroutine writeMetalconstr

logical function isTMetal(i)
  implicit none
  integer :: i
  isTMetal = .false.
  if (i .ge. 21.and.i .le. 30) isTMetal = .true.
  if (i .ge. 39.and.i .le. 48) isTMetal = .true.
  if (i .ge. 57.and.i .le. 80) isTMetal = .true.
  return
end function isTMetal

!====================================================!
!> Write a constraint file with all bonds except X-H
!====================================================!
subroutine writeHeavyconstr(nat,at,bmat,force)
  use crest_parameters,only:wp
  implicit none
  integer :: nat
  integer :: at(nat)
  real(wp) :: bmat(nat,nat)
  integer :: i,j,ich
  real(wp) :: force
  character(len=20) :: dumm

  open (newunit=ich,file='bondlengths')
  write (ich,'(a)') '$constrain'
  write (dumm,'(f16.4)') force
  write (ich,'(3x,a,a)') 'force constant=',adjustl(trim(dumm))
  do i = 1,nat
    do j = i+1,nat
      if (bmat(i,j) .gt. 0.1_wp) then
        if ((at(i) .ne. 1).and.(at(j) .ne. 1)) then
          write (ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',i,',',j,',',bmat(i,j)
        end if
      end if
    end do
  end do
  write (ich,'(a)') '$end'
  close (ich)
  return
end subroutine writeHeavyconstr
!====================================================!
!> Write a constraint file with all X-H bonds
!====================================================!
subroutine writeHydrogenconstr(nat,at,bmat,force)
  use crest_parameters,only:wp
  implicit none
  integer :: nat
  integer :: at(nat)
  real(wp) :: bmat(nat,nat)
  integer :: i,j,ich
  real(wp) :: force
  character(len=20) :: dumm

  open (newunit=ich,file='bondlengths')
  write (ich,'(a)') '$constrain'
  write (dumm,'(f16.4)') force
  write (ich,'(3x,a,a)') 'force constant=',adjustl(trim(dumm))
  do i = 1,nat
    do j = i+1,nat
      if (bmat(i,j) .gt. 0.1_wp) then
        if ((at(i) .eq. 1).or.(at(j) .eq. 1)) then
          write (ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',i,',',j,',',bmat(i,j)
        end if
      end if
    end do
  end do
  write (ich,'(a)') '$end'
  close (ich)
  return
end subroutine writeHydrogenconstr

!====================================================!
!> Write a constraint file with all bonds
!====================================================!
subroutine writeBmatconstr_withEZ(zmol,nat,bmat,force)
  use crest_parameters,only:wp,bohr
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer :: nat
  real(wp) :: bmat(nat,nat)
  integer :: i,j,k,l,ich
  real(wp) :: force
  character(len=20) :: dumm

  integer :: ni,nj
  real(wp) :: ndist
  integer :: r
  logical :: smallring

  real(wp),parameter :: springexp = 4.0_wp
  real(wp),parameter :: cclen = 1.384_wp
  integer,parameter  :: ringmax = 5
  integer,parameter  :: ringmax2 = 8

  open (newunit=ich,file='bondlengths')
  write (ich,'(a)') '$constrain'
  write (dumm,'(f16.4)') force
  write (ich,'(3x,a,a)') 'force constant=',adjustl(trim(dumm))
  do i = 1,nat
    do j = i+1,nat
      if (bmat(i,j) .gt. 0.1_wp) then

        smallring = .false.

        if (.not.smallring) then
          write (ich,'(3x,a,i0,a,i0,a,f8.5)') 'distance: ',i,', ',j,', ',bmat(i,j)
        else
          write (ich,'(3x,a,i0,a,i0,a,f8.5,a,f8.5)') 'DISTANCE: ',i,', ',j,', ',bmat(i,j), &
         &       ', ',springexp
        end if
      end if
    end do
  end do
  !write(ich,*) '#E/Z constraints'
  ILOOP: do i = 1,nat
    JLOOP: do j = i+1,nat
      if (bmat(i,j) .gt. 0.1_wp) then
        if ((zmol%at(i) == 6.and.zmol%at(j) == 6).and.  &                     ! ij are Carbon atoms?
        &  (nint(zmol%zat(i)%cn) == 3.and.nint(zmol%zat(j)%cn) == 3).and. &  ! ij are sp2 C?
        &  (bmat(i,j) < cclen)) then                                       ! r_ij could be C=C bond?

          ! check if ij define a C=C bond within a small ring? (cycle, if so)
          do r = 1,zmol%nri
            if (any(zmol%zri(r)%rlist == i).and. &
            &  any(zmol%zri(r)%rlist == j).and. &
            &  (zmol%zri(r)%rs < ringmax2)) then
              cycle JLOOP
            end if
          end do

          do k = 1,zmol%zat(i)%nei
            ni = zmol%zat(i)%ngh(k)
            if (ni == j) cycle
            do l = 1,zmol%zat(j)%nei
              nj = zmol%zat(j)%ngh(l)
              if (nj == i) cycle
              ndist = zmol%dist(ni,nj)*bohr
              !write(ich,'(3x,a,i0,a,i0,a,f8.5)') 'distance: ',ni,', ',nj,', ',ndist
              !--- WARNING: The constraint has to be case-sensitive in order to have an adjustable springexp.
              write (ich,'(3x,a,i0,a,i0,a,f8.5,a,f8.5)') 'DISTANCE: ',ni,', ',nj,', ',ndist, &
              &       ', ',springexp
            end do
          end do
        end if
      end if
    end do JLOOP
  end do ILOOP
  write (ich,'(a)') '$end'
  close (ich)

  return
end subroutine writeBmatconstr_withEZ

!========================================================================================!
!> Routines to create constraints automatically for
!> the internal calculation routines of CREST
!>--------------------------------------------------!
subroutine autoconstraint_internal(env)
  use crest_parameters,only:wp,bohr
  use crest_data
  implicit none
  type(systemdata),intent(inout) :: env
  integer :: i,j,k,l
  integer :: allbonds,allangl,alldihed

  l = env%calc%nconstraints
  allbonds = 0
  allangl = 0
  alldihed = 0

  if (l > 0) then
    do i = 1,l
      k = env%calc%cons(i)%type
      if (k == 11) allbonds = i
      if (k == 21) allangl = i
      if (k == 31) alldihed = i
    end do
  end if

  if (allbonds > 0) then
    call autobond_internal(env,allbonds)
    call env%calc%removeconstraint(allbonds)
  end if

  if (allangl > 0) then
  end if

  if (alldihed > 0) then
  end if

  return
end subroutine autoconstraint_internal
subroutine autobond_internal(env,c)
  use crest_parameters,only:wp,bohr
  use crest_data
  use crest_calculator,only:constraint
  use zdata
  implicit none
  type(systemdata) :: env
  integer :: c
  integer :: i,j,k,l
  real(wp) :: kforce
  type(zmolecule) :: zmol
  real(wp),allocatable :: bmat(:,:)
  integer :: nat,nb
  type(constraint) :: constr

  if (c > 0) then
    kforce = env%calc%cons(c)%fc(1)
  else
    kforce = env%forceconst
  end if
  !>--- get topology
  call simpletopo_file(env%inputcoords,zmol,.false.,.false.,env%wbofile)
  !>--- get bond matrix
  nat = zmol%nat
  allocate (bmat(nat,nat),source=0.0_wp)
  nb = 0
  do i = 1,nat
    do j = 1,nat
      if (any(zmol%zat(i)%ngh(:) .eq. j)) then  !only include bonds from the neighbour lists
        bmat(i,j) = zmol%distmat(i,j)
        bmat(i,j) = bmat(i,j)   !BMAT is not in Angstroem!
        nb = nb+1
      else
        cycle
      end if
    end do
  end do
  !>--- creat constraint objects
  do i = 1,nat
    do j = i+1,nat
      if (bmat(i,j) .gt. 0.1_wp) then
        call constr%deallocate()
        call constr%bondconstraint(i,j,bmat(i,j),kforce)
        call env%calc%add(constr)
      end if
    end do
  end do

  deallocate (bmat)
  call zmol%deallocate()
  return
end subroutine autobond_internal

!=========================================================================================!
!> LOCAL EXTENSION (2026-10-12): generate "bondlengths" for a COMBINATION of the
!> automatic bond-constraint families -cbonds / -cmetal / -cheavy / -clight.
!>
!> WHY THIS EXISTS
!> ---------------
!> The families used to be mutually exclusive: confparse kept a single `ctype`,
!> so passing two of them silently kept only the LAST one.  "-cbonds" also meant
!> ALL bonds, so combining it with -cmetal (or with --cinp) would put two springs
!> on the same distance - and for --cinp even a harmonic term on top of a rigid
!> SHAKE bond.  Here "-cbonds" becomes THE REMAINDER: every bond that no other
!> REQUESTED family and no --cinp constraint already owns.  Used alone it
!> reproduces the old "constrain all bonds" file.
!>
!> PER-FAMILY FORCE CONSTANTS
!> --------------------------
!> One "$constrain" block is written per family, each with its own
!> "force constant=".  That works because parse_xtb_inputfile() loops over ALL
!> blocks of a file and get_xtb_constraint_block() reads the force constant once
!> per block (parse_xtbinput.f90:94-108 and 142-153); one block can only carry
!> one k.
!>
!> --cinp BONDS ARE NEVER TOUCHED HERE.  They have their own file and become
!> harmonic restraints in calc%cons AND rigid (SHAKE) bonds via add_hardcons().
!> A bond counts as owned by --cinp if
!>   * it appears as a "distance:"/"bond:" entry, or
!>   * BOTH of its atoms appear in an "atoms:" entry (a positional restraint
!>     already fixes the distance between two restrained atoms).
!>
!> REMINDER: the reference distance is the bond length of the INPUT structure
!> (see getbmat), so these restraints hold whatever geometry you feed in.  They
!> do not by themselves correct a systematic method bias - feed a geometry whose
!> bond lengths you trust.
!=========================================================================================!
subroutine autoConstraintUnion(filename,wbofile,cinpfile, &
    &          do_bonds,do_metal,do_heavy,do_light, &
    &          fc_bonds,fc_metal,fc_heavy,fc_light,fc_global)
  use crest_parameters,only:wp,bohr
  use zdata
  implicit none
  character(len=*),intent(in) :: filename
  character(len=*),intent(in) :: wbofile
  character(len=*),intent(in) :: cinpfile
  logical,intent(in)  :: do_bonds,do_metal,do_heavy,do_light
  real(wp),intent(in) :: fc_bonds,fc_metal,fc_heavy,fc_light,fc_global

  type(zmolecule) :: zmol
  integer :: nat,i,j,ich
  integer,allocatable :: at(:)
  real(wp),allocatable :: bmat(:,:)
  logical,allocatable :: iscinp(:,:)
  real(wp) :: kb,km,kh,kl

  !>--- topology and neighbour-bond list; distances converted to Angstrom
  !>    (the "bondlengths" file is written in Angstrom - see getbmat())
  call simpletopo_file(filename,zmol,.false.,.false.,wbofile)
  nat = zmol%nat
  allocate (at(nat),source=0)
  at = zmol%at
  allocate (bmat(nat,nat),source=0.0_wp)
  do i = 1,nat
    do j = 1,nat
      if (any(zmol%zat(i)%ngh(:) .eq. j)) bmat(i,j) = zmol%distmat(i,j)*bohr
    end do
  end do

  !>--- bonds that --cinp already owns
  allocate (iscinp(nat,nat),source=.false.)
  if (len_trim(cinpfile) > 0) call read_cinp_bonds()

  !>--- effective force constants; a family without an explicit value falls back
  !>    to the global one, which preserves the old single-family behaviour
  kb = fc_bonds
  km = fc_metal
  kh = fc_heavy
  kl = fc_light
  if (kb < 0.0_wp) kb = fc_global
  if (km < 0.0_wp) km = fc_global
  if (kh < 0.0_wp) kh = fc_global
  if (kl < 0.0_wp) kl = fc_global

  open (newunit=ich,file='bondlengths')
  if (do_metal) call writeblock(ich,1,km)
  if (do_light) call writeblock(ich,2,kl)
  if (do_heavy) call writeblock(ich,3,kh)
  if (do_bonds) call writeblock(ich,0,kb)   !> family 0 = the remainder
  write (ich,'(a)') '$end'
  close (ich)

  deallocate (bmat,at,iscinp)
  call zmol%deallocate()
  return

contains

!>--- fam 1 = metal, 2 = X-H, 3 = heavy-heavy, 0 = the "-cbonds" remainder
  subroutine writeblock(unit,fam,k)
    integer,intent(in) :: unit,fam
    real(wp),intent(in) :: k
    character(len=20) :: dumm
    integer :: a,b,n
    n = 0
    do a = 1,nat
      do b = a+1,nat
        if (bmat(a,b) <= 0.1_wp) cycle
        if (iscinp(a,b)) cycle
        if (infam(a,b,fam)) n = n+1
      end do
    end do
    if (n < 1) return
    write (dumm,'(f16.4)') k
    write (unit,'(a)') '$constrain'
    write (unit,'(3x,a,a)') 'force constant=',adjustl(trim(dumm))
    do a = 1,nat
      do b = a+1,nat
        if (bmat(a,b) <= 0.1_wp) cycle
        if (iscinp(a,b)) cycle
        if (.not.infam(a,b,fam)) cycle
        write (unit,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',a,',',b,',',bmat(a,b)
      end do
    end do
    return
  end subroutine writeblock

!>--- does bond (a,b) belong to family fam?
  logical function infam(a,b,fam) result(res)
    integer,intent(in) :: a,b,fam
    logical :: ism,isxh,ishv
    ism = tmetal(at(a)) .or. tmetal(at(b))
    isxh = (at(a) == 1) .or. (at(b) == 1)
    ishv = (at(a) /= 1) .and. (at(b) /= 1)
    select case (fam)
    case (1)
      res = ism
    case (2)
      res = isxh
    case (3)
      res = ishv
    case default   !> 0: everything not claimed by a REQUESTED other family
      res = .not. ((do_metal.and.ism).or.(do_light.and.isxh).or.(do_heavy.and.ishv))
    end select
    return
  end function infam

!>--- Z ranges of the d- and f-block metals (same criterion as isTMetal())
  logical function tmetal(z) result(res)
    integer,intent(in) :: z
    res = ((z >= 21.and.z <= 30).or.(z >= 39.and.z <= 48).or.(z >= 57.and.z <= 80))
    return
  end function tmetal

!>--- collect the bonds that the --cinp file already constrains
  subroutine read_cinp_bonds()
    character(len=512) :: line,s
    character(len=64)  :: key
    logical,allocatable :: locked(:)
    logical :: ex
    integer :: unit,io,a,b,k,l,ng,iv(2)
    inquire (file=trim(cinpfile),exist=ex)
    if (.not.ex) return
    allocate (locked(nat),source=.false.)
    open (newunit=unit,file=trim(cinpfile),status='old',iostat=io)
    if (io /= 0) then
      deallocate (locked)
      return
    end if
    do
      read (unit,'(a)',iostat=io) line
      if (io /= 0) exit
      s = adjustl(line)
      if (len_trim(s) == 0) cycle
      if (s(1:1) == '$'.or.s(1:1) == '#'.or.s(1:1) == '!') cycle
      k = index(s,':')
      l = index(s,'=')
      if (k > 0.and.(l == 0.or.k < l)) then
        key = s(1:k-1)
        s = s(k+1:)
      else if (l > 0) then
        key = s(1:l-1)
        s = s(l+1:)
      else
        cycle
      end if
      call tolower(key)
      select case (trim(key))
      case ('distance','bond')
        call grab_ints(s,iv,2,ng)
        if (ng >= 2) then
          a = min(iv(1),iv(2))
          b = max(iv(1),iv(2))
          if (a >= 1.and.b <= nat.and.a /= b) iscinp(a,b) = .true.
        end if
      case ('atoms','atom')
        call mark_atoms(s,locked)
      end select
    end do
    close (unit)
    !>--- two atoms under a positional restraint have a fixed distance
    do a = 1,nat
      do b = a+1,nat
        if (locked(a).and.locked(b)) iscinp(a,b) = .true.
      end do
    end do
    deallocate (locked)
    return
  end subroutine read_cinp_bonds

!>--- first nwant integers of a string, commas treated as spaces
  subroutine grab_ints(str,out,nwant,ngot)
    character(len=*),intent(in) :: str
    integer,intent(in) :: nwant
    integer,intent(out) :: out(nwant)
    integer,intent(out) :: ngot
    character(len=512) :: s
    character(len=32) :: tok
    integer :: p,q,n,iv
    out = 0
    ngot = 0
    s = str
    call commatospace(s)
    p = 1
    do while (p <= len_trim(s).and.ngot < nwant)
      do while (p <= len_trim(s).and.s(p:p) == ' ')
        p = p+1
      end do
      if (p > len_trim(s)) exit
      q = p
      do while (q <= len_trim(s).and.s(q:q) /= ' ')
        q = q+1
      end do
      tok = s(p:min(q-1,p+31))
      p = q
      read (tok,*,iostat=n) iv
      if (n == 0) then
        ngot = ngot+1
        out(ngot) = iv
      end if
    end do
    return
  end subroutine grab_ints

!>--- mark every atom index and every "N-M" range found in a string
  subroutine mark_atoms(str,locked)
    character(len=*),intent(in) :: str
    logical,intent(inout) :: locked(nat)
    character(len=512) :: s
    character(len=32) :: tok
    integer :: p,q,n,lo,hi,mid
    s = str
    call commatospace(s)
    p = 1
    do while (p <= len_trim(s))
      do while (p <= len_trim(s).and.s(p:p) == ' ')
        p = p+1
      end do
      if (p > len_trim(s)) exit
      q = p
      do while (q <= len_trim(s).and.s(q:q) /= ' ')
        q = q+1
      end do
      tok = s(p:min(q-1,p+31))
      p = q
      mid = index(tok,'-')
      if (mid > 1.and.mid < len_trim(tok)) then
        read (tok(1:mid-1),*,iostat=n) lo
        if (n /= 0) cycle
        read (tok(mid+1:),*,iostat=n) hi
        if (n /= 0) cycle
      else
        read (tok,*,iostat=n) lo
        if (n /= 0) cycle
        hi = lo
      end if
      if (hi < lo) then
        mid = lo
        lo = hi
        hi = mid
      end if
      lo = max(1,lo)
      hi = min(nat,hi)
      do mid = lo,hi
        locked(mid) = .true.
      end do
    end do
    return
  end subroutine mark_atoms

  subroutine commatospace(s)
    character(len=*),intent(inout) :: s
    integer :: i
    do i = 1,len(s)
      if (s(i:i) == ',') s(i:i) = ' '
    end do
    return
  end subroutine commatospace

  subroutine tolower(s)
    character(len=*),intent(inout) :: s
    integer :: i,c
    do i = 1,len_trim(s)
      c = iachar(s(i:i))
      if (c >= 65.and.c <= 90) s(i:i) = achar(c+32)
    end do
    return
  end subroutine tolower

end subroutine autoConstraintUnion
