!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2025 Philipp Pracht
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

!> This is a collection of various IO subroutines

module iomod
  use iso_fortran_env,wp => real64,stdout => output_unit
  use iso_c_binding
  implicit none
  public

  interface
    function mkdir(path,mode) bind(c,name="mkdir")
      use iso_c_binding
      integer(c_int) :: mkdir
      character(kind=c_char,len=1) :: path(*)
      integer(c_int16_t),value :: mode
    end function mkdir
  end interface

  interface
    !> NOTE: POSIX setenv() takes THREE arguments:
    !>   int setenv(const char *name, const char *value, int overwrite);
    !> The binding used to declare only two, so the "overwrite" flag was
    !> read from an undefined register.  If that happened to be 0, setenv()
    !> silently REFUSED to update an already existing variable.
    !> Consequences: ompenvset() could not change OMP_NUM_THREADS /
    !> MKL_NUM_THREADS / OPENBLAS_NUM_THREADS, so externally started
    !> programs (xtb, orca, ... subprocesses) inherited whatever the
    !> calling shell had exported - e.g. OPENBLAS_NUM_THREADS=1, which
    !> serialises precisely the BLAS work that xtb parallelises with
    !> OMP_NUM_THREADS.
    integer(kind=c_int) function c_setenv(c_name,c_VALUE,c_overwrite) bind(c,name="setenv")
      use iso_c_binding
      character(kind=c_char)   :: c_name(*)
      character(kind=c_char)   :: c_VALUE(*)
      integer(kind=c_int),value :: c_overwrite
    end function
  end interface
  interface setenv
    module procedure setenv_string
    module procedure setenv_int
    module procedure setenv_float
  end interface setenv

  interface
    integer(kind=c_int) function c_symlink(c_from,c_to) bind(c,name="symlink")
      use iso_c_binding
      character(kind=c_char)   :: c_from(*)
      character(kind=c_char)   :: c_to(*)
    end function
  end interface

  interface wrshort
    module procedure wrshort_real
    module procedure wrshort_int
    module procedure wrshort_string
  end interface wrshort

  interface rdshort
    module procedure rdshort_real
    module procedure rdshort_int
    module procedure rdshort_string
  end interface rdshort

  interface uppercase
!    module procedure to_upper
    module procedure to_upper_func
  end interface uppercase

  interface lowercase
!    module procedure to_lower
    module procedure to_lower_func
  end interface lowercase

  interface convert_to_string
    module procedure i8_to_string
    module procedure i16_to_string
    module procedure i32_to_string
    module procedure i64_to_string
    module procedure r32_to_string
    module procedure r64_to_string
    module procedure bool_to_string
  end interface convert_to_string
  interface to_str
    module procedure to_str_i8
    module procedure to_str_i16
    module procedure to_str_i32
    module procedure to_str_i64
    module procedure to_str_r32
    module procedure to_str_r64
    module procedure to_str_bool
  end interface to_str

  interface
    function c_realpath(source_ptr,resolved_ptr) bind(c,name="realpath")
      use iso_c_binding
      implicit none
      type(c_ptr) :: c_realpath  !> function returns a pointer to char

      !> Both arguments must be type(c_ptr), value
      type(c_ptr),value :: source_ptr
      type(c_ptr),value :: resolved_ptr
    end function c_realpath
  end interface

  public :: checkprog
  public :: checkprog_silent
  public :: absolute_path,absolute_filepath
  private :: getpath

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!
!-------------------------------------------------------------
! remove file 'fname'
! if 'fname' contains the wildcard '*', the remove is ignored
!--------------------------------------------------------------
  subroutine remove(fname)
    implicit none
    character(len=*) :: fname
    integer :: ich,och
    logical :: ex
    if (index(fname,'*') .eq. 0) then
      open (newunit=ich,file=fname)
      close (ich,status='delete')
    end if
  end subroutine remove

!-------------------------------------------------------------
! creates an empty file "fname", similar to the shell command "touch"
!--------------------------------------------------------------
  subroutine touch(fname)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    open (newunit=ich,file=fname)
    close (ich)
  end subroutine touch

!-------------------------------------------------------------
! prints the file "fname", similar to the shell command "cat"
!--------------------------------------------------------------
  subroutine cat(fname)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    character(len=512) :: str
    integer :: io
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) str
      if (io < 0) exit
      write (stdout,'(a)') trim(str)
    end do
    close (ich)
  end subroutine cat

!--------------------------------------------------------------------------------------
! prints the file "fname", similar to the shell command "cat" but deletes it at the end
!--------------------------------------------------------------------------------------
  subroutine catdel(fname)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    character(len=512) :: str
    integer :: io
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) str
      if (io < 0) exit
      write (stdout,'(a)') trim(str)
    end do
    close (unit=ich,status='delete')
  end subroutine catdel

!-------------------------------------------------------------------------
! copy a file from path "from" to path "to"
!-------------------------------------------------------------------------
  subroutine copy(from,to)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: from
    character(len=*) :: to
    character(len=512) :: line
    integer :: io
    inquire (file=from,exist=ex)
    if (ex) then
      open (newunit=ich,file=from)
      open (newunit=och,file=to)
      do
        read (ich,'(a)',iostat=io) line
        if (io < 0) exit
        write (och,'(a)') trim(line)
      end do
      close (och)
      close (ich)
    end if
  end subroutine copy

!-------------------------------------------------------------------------
! copy a file from path "from" to a specified sub-directory "to"
!-------------------------------------------------------------------------
  subroutine copysub(from,to)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: from
    character(len=*) :: to
    character(len=512) :: line,dir1
    integer :: io
    inquire (file=from,exist=ex)
    if (ex) then
      call getcwd(dir1)
      inquire (file=from,exist=ex)
      if (ex) then
        open (newunit=ich,file=from)
        call chdir(to)
        open (newunit=och,file=from)
        do
          read (ich,'(a)',iostat=io) line
          if (io < 0) exit
          write (och,'(a)') trim(line)
        end do
        close (och)
        close (ich)
      end if
      call chdir(dir1)
    end if
  end subroutine copysub

!-------------------------------------------------------------------------
! move a file from path "from" to path "to". "from" file is deleted.
!-------------------------------------------------------------------------
  subroutine move(from,to)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: from
    character(len=*) :: to
    inquire (file=from,exist=ex)
    if (ex) then
      call remove(to)
      call rename(from,to)
    else
      write (stdout,'(''file '',a,'' does not exist!'')') trim(from)
    end if
  end subroutine move

!-------------------------------------------------------------------------
! converts a path relative to the current directory (i.e. containing "../"
! and such) into the absolute Path
!-------------------------------------------------------------------------
  subroutine getAbsPath(Path)
    implicit none
    character(len=*),intent(inout) :: Path
    character(len=512) :: thispath
    call getcwd(thispath)  !get current working directory
    call chdir(Path)
    call getcwd(Path)  !get the absolute path which was given as a relative path
    call chdir(thispath) !switch back to original directory
  end subroutine getAbsPath

!-------------------------------------------------------------------------
! the "inquire" command, but in a sub-directory
!-------------------------------------------------------------------------
  subroutine inquiresub(f,sub,bool)
    implicit none
    character(len=*),intent(in) :: f
    character(len=*),intent(in) :: sub
    character(len=512) :: thispath
    logical,intent(out) :: bool
    call getcwd(thispath)  !get current working directory
    call chdir(sub) !go to the subdirectory
    inquire (file=f,exist=bool) !do the "inquire"
    call chdir(thispath) !switch back to original directory
  end subroutine inquiresub

!-------------------------------------------------------------------------
! append content of text file "from" into text file "to", similar to "cat A >> B"
!-------------------------------------------------------------------------
!>--- LOCAL CHANGE (2026-09-27): open the target for APPENDING instead of
!> emulating it with a read-to-EOF scan.  The old code did
!>     open(newunit=ich,file=to)      ! position 0
!>     do; read(ich,*,iostat=io); ...  ! walk to the end
!> which re-read the WHOLE accumulated target on every call.  collectcre()
!> (legacy_algos/confscript2_misc.f90) appends one .cre_* file per restart
!> round through this routine, so merging R files cost O(R^2) reading - and on
!> NFS every one of those passes is a network round trip.  POSITION='APPEND'
!> (F2003, supported by gfortran and ifx) makes each append O(1).
!-------------------------------------------------------------------------
  subroutine appendto(from,to)
    implicit none
    integer :: ich,och
    logical :: ex
    integer :: io
    character(len=*) :: from
    character(len=*) :: to
    character(len=1024) :: str
    open (newunit=ich,file=to,position='append')
    open (newunit=och,file=from)
    do
      read (och,'(a)',iostat=io) str
      if (io < 0) then
        exit
      else
        write (ich,'(a)') trim(str)
      end if
    end do
    close (och)
    close (ich)
  end subroutine appendto

!-------------------------------------------------------------------------
! make a directory via iso_c_binding
!-------------------------------------------------------------------------
  function makedir(str)
    implicit none
    integer :: makedir
    character(len=*) :: str
    makedir = mkdir(str//char(0),int(o'770',c_int16_t)) !create new directory
    return
  end function makedir

!-------------------------------------------------------------------------
! put a environment variable via iso_c_binding
!-------------------------------------------------------------------------
  function setenv_string(env,str)
    implicit none
    integer :: setenv_string
    character(len=*) :: env
    character(len=*) :: str
    setenv_string = c_setenv(env//c_null_char,str//c_null_char,1_c_int) !> 1 = overwrite existing value
    return
  end function setenv_string
  function setenv_int(env,intval)
    implicit none
    integer :: setenv_int
    character(len=*) :: env
    integer :: intval
    character(len=20) :: str
    write (str,'(i0)') intval
    setenv_int = c_setenv(env//c_null_char,trim(str)//c_null_char,1_c_int) !> 1 = overwrite existing value
    return
  end function setenv_int
  function setenv_float(env,floatval)
    implicit none
    integer :: setenv_float
    character(len=*) :: env
    real(wp) :: floatval
    character(len=30) :: str
    write (str,'(f14.6)') floatval
    setenv_float = c_setenv(env//c_null_char,trim(str)//c_null_char,1_c_int) !> 1 = overwrite existing value
    return
  end function setenv_float

!-------------------------------------------------------------------------
! set a symlink from path1 to path2  via iso_c_binding
!-------------------------------------------------------------------------
  function sylnk(path1,path2)
    implicit none
    integer :: sylnk
    character(len=*) :: path1
    character(len=*) :: path2
    sylnk = c_symlink(trim(path1)//c_null_char,trim(path2)//c_null_char) !create new directory
    return
  end function sylnk

!-------------------------------------------------------------------------
! minigrep: a grep subroutine that returns true or false, depending on
! if the substring "str" is present in the file "fil"
!-------------------------------------------------------------------------
  subroutine minigrep(fil,str,bool)
    implicit none
    integer :: ich,och
    logical :: ex
    logical :: bool
    character(len=*) :: fil
    character(len=*) :: str
    character(len=512) :: tmp
    integer :: io
    bool = .false.
    open (newunit=ich,file=fil)
    do
      read (ich,'(a)',iostat=io) tmp
      if (io < 0) exit
      if (index(tmp,str) .ne. 0) then
        bool = .true.
        close (ich)
        exit
      end if
    end do
    close (ich)
    !minigrep = bool
    return
  end subroutine minigrep

!-------------------------------------------------------------------------
! grepval: a grep subroutine that returns the float that follows immeadeatly
! after the substring "str". Works only for the first occurence of "str".
! If str is not present, "bool" returns .false.
!-------------------------------------------------------------------------
  subroutine grepval(fil,str,bool,val)
    implicit none
    integer :: ich,och
    logical :: ex
    logical :: bool
    character(len=*) :: fil
    character(len=*) :: str
    character(len=512) :: tmp
    real(wp) :: val
    integer :: io
    bool = .false.
    val = 0.0d0
    open (newunit=ich,file=fil)
    do
      read (ich,'(a)',iostat=io) tmp
      if (io < 0) exit
      io = index(tmp,str,.true.)
      if (io .ne. 0) then
        io = io+len(str)
        tmp = adjustl(tmp(io+1:))
        read (tmp,*) val
        bool = .true.
        close (ich)
        exit
      end if
    end do
    close (ich)
    !minigrep = bool
    return
  end subroutine grepval

!---------------------------------------------------------------------------
! grepcntx: a grep subroutine that returns the entire line containing "str"
! or the specified line after the occurence of "str", given
! If str is not present, "bool" returns .false.
!---------------------------------------------------------------------------
  subroutine grepcntxt(fil,str,bool,line,context)
    implicit none
    integer :: ich,och
    logical :: ex
    logical,intent(out) :: bool
    logical :: track
    character(len=*),intent(in) :: fil
    character(len=*),intent(in) :: str
    character(len=*),intent(out) :: line
    integer,intent(in),optional :: context
    integer :: cntxt
    character(len=512) :: tmp
    real(wp) :: val
    integer :: io,itrack,reftrack
    bool = .false.
    track = .false.
    val = 0.0d0
    line = ''
    itrack = 0
    if (present(context)) then
      reftrack = max(0,context)
    else
      reftrack = 0
    end if
    open (newunit=ich,file=fil)
    do
      read (ich,'(a)',iostat=io) tmp
      if (io < 0) exit
      io = index(tmp,str,.true.)
      if (io .ne. 0) then
        bool = .true.
        track = .true.
      end if
      if (track) then
        if (itrack == reftrack) then
          line = trim(tmp)
          exit
        else
          itrack = itrack+1
        end if
      end if
    end do
    close (ich)
    return
  end subroutine grepcntxt

!=====================================================================================!
  function getlines(fname,maxwidth) result(n)
    character(len=*) :: fname
    logical :: ex
    integer,intent(out),optional :: maxwidth
    integer :: n,io,w,ich
    character(len=1056) :: atmp
    n = 0
    w = 0
    inquire (file=fname,exist=ex)
    if (.not.ex) return
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit
      n = n+1
      if (len_trim(atmp) > w) w = len_trim(atmp)
    end do
    if (present(maxwidth)) maxwidth = w
    close (ich)
  end function getlines

!--------------------------------------------------------------------------------------
! Write a file with a single line (INT)
!--------------------------------------------------------------------------------------
  subroutine wrshort_int(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    integer :: var
    open (newunit=ich,file=fname)
    write (ich,*) var
    close (ich)
  end subroutine wrshort_int

!--------------------------------------------------------------------------------------
! Write a file with a single line (REAL)
!--------------------------------------------------------------------------------------
  subroutine wrshort_real(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    real(wp) :: var
    open (newunit=ich,file=fname)
    write (ich,*) var
    close (ich)
  end subroutine wrshort_real

!--------------------------------------------------------------------------------------
! Write a file with a single line (STRING)
!--------------------------------------------------------------------------------------
  subroutine wrshort_string(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    character(len=*) :: var
    open (newunit=ich,file=fname)
    write (ich,'(a)') var
    close (ich)
  end subroutine wrshort_string

!--------------------------------------------------------------------------------------
! reads the first line of a file (INT)
!--------------------------------------------------------------------------------------
  subroutine rdshort_int(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    integer :: var
    integer :: io
    open (newunit=ich,file=fname)
    read (ich,*,iostat=io) var
    if (io .ne. 0) var = 0
    close (ich)
  end subroutine rdshort_int

!--------------------------------------------------------------------------------------
! reads the first line of a file (REAL)
!--------------------------------------------------------------------------------------
  subroutine rdshort_real(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    real(wp) :: var
    integer :: io
    open (newunit=ich,file=fname)
    read (ich,*,iostat=io) var
    if (io .ne. 0) var = 0.0d0
    close (ich)
  end subroutine rdshort_real

!--------------------------------------------------------------------------------------
! reads the first line of a file (STRING)
!--------------------------------------------------------------------------------------
  subroutine rdshort_string(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    character(len=*) :: var
    integer :: io
    open (newunit=ich,file=fname)
    read (ich,'(a)',iostat=io) var
    if (io .ne. 0) var = ''
    close (ich)
  end subroutine rdshort_string

!----------------------------------------------------------------------------
! returns entry x of line as string, or as empty if end-of-line is reached

  subroutine glinex(line,x,string)
    implicit none
    integer :: ich,och
    logical :: ex
    integer,intent(in) :: x
    character(len=*),intent(in)  :: line
    character(len=*),intent(out) :: string
    character(len=80) :: stmp
    character(len=1)  :: digit
    integer :: i,k,cs
    logical :: track

    string = ''
    stmp = ''
    cs = 0
    track = .false.
    k = len(trim(line))
    do i = 1,k
      digit = line(i:i)
      if (digit .eq. ' '.or.digit .eq. char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
        if (track) then
          if (cs == x) then
            string = trim(stmp)
            exit
          end if
        end if
        track = .false.
        stmp = ''
        cycle
      else
        if (.not.track) cs = cs+1
        track = .true.
        stmp = trim(stmp)//trim(digit)
        if (i .eq. k.and.x .lt. 0) then  !if x is lt. 0, the last element of the line is returned
          string = trim(stmp)
          exit
        end if
      end if
    end do
    return
  end subroutine glinex

!----------------------------------------------------------------------------
! cut the leading x words (seperated by blanks) of a string, and overwrite it

  subroutine clinex(line,x)
    implicit none
    integer,intent(in) :: x
    character(len=*),intent(inout)  :: line
    character(len=1)  :: digit
    integer :: i,j,k,cs
    logical :: track

    cs = 0
    track = .false.
    k = len(trim(line))
    do i = 1,k
      digit = line(i:i)
      if (digit .eq. ' '.or.digit .eq. char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
        track = .false.
        cycle
      else
        if (.not.track) cs = cs+1
        if (cs .gt. x) then
          j = i-1
          line = line(j:k)
          exit
        end if
        track = .true.
      end if
    end do
    return
  end subroutine clinex

!-------------------------------------------------------------
! convert a string to upper or lower case
!--------------------------------------------------------------
  subroutine to_upper(str)
    character(*),intent(in out) :: str
    integer :: i
    do i = 1,len(str)
      select case (str(i:i))
      case ("a":"z")
        str(i:i) = achar(iachar(str(i:i))-32)
      end select
    end do
  end subroutine to_upper
  function to_upper_func(str) result(out)
    character(*),intent(in) :: str
    character(len=:),allocatable :: out
    out = str
    call to_upper(out)
  end function to_upper_func
  subroutine to_lower(str)
    character(*),intent(in out) :: str
    integer :: i
    do i = 1,len(str)
      select case (str(i:i))
      case ("A":"Z")
        str(i:i) = achar(iachar(str(i:i))+32)
      end select
    end do
  end subroutine to_lower
  function to_lower_func(str) result(out)
    character(*),intent(in) :: str
    character(len=:),allocatable :: out
    out = str
    call to_lower(out)
  end function to_lower_func

!-----------------------------------------------------------------
! a simple routine to lock a feature of the code behind a password
!-----------------------------------------------------------------
  subroutine pwdlock(pwdref)
    implicit none
    character(len=*) :: pwdref    ! check for this password
    character(len=:),allocatable :: pwd
    pwd = pwdref !allocate max length of pwd
    write (stdout,'(/,1x,a)',advance='no') 'Locked feature. Enter password: '
    read (*,*) pwd
    if (pwd == pwdref) then
      write (stdout,'(1x,a,/)') 'Valid. Continue.'
    else
      error stop 'Invalid. Stop.'
    end if
    return
  end subroutine pwdlock

!====================================================================!
  function filechecker(fin,fout) result(have)
    implicit none
    logical :: have
    character(len=*),intent(in) :: fin
    character(len=*),intent(out) :: fout
    have = .false.
    inquire (file=fin,exist=have)
    if (have) then
      fout = fin
    else
      fout = ''
    end if
    return
  end function filechecker

!=========================================================================================!
!=========================================================================================!
!=========================================================================================!

!> Checking directories with Fortran's inquire is not handled by the standard.
!> The interpretation whether or not to report a directory as file is compiler
!> specific and therefore always an extension to the Fortran standard.
  function directory_exist(file) result(exist)
    character(len=*),intent(in) :: file
    logical :: exist
#ifdef __INTEL_COMPILER
    !> Intel provides the directory extension to inquire to handle this case
    inquire (directory=file,exist=exist)
#else
    !> GCC handles directories as files, to make sure we get a directory and
    !> not a file append a path separator and the current dir
    inquire (file=trim(file)//"/.",exist=exist)
#endif
  end function directory_exist

!=========================================================================================!
!=========================================================================================!
!=========================================================================================!

!> a wrapper for the intrinsic isatty function.
!> ifort only seems to work if isatty is declard as external
!> while gfortran does not want that...
  function myisatty(channel) result(term)
    implicit none
    integer,intent(in) :: channel
    logical :: term
#ifdef __INTEL_COMPILER
    logical,external :: isatty
#endif
    term = isatty(channel)
  end function myisatty

!=========================================================================================!
!=========================================================================================!
!=========================================================================================!
!> type conversion routines
!> 8 bit integers
  subroutine i8_to_string(i_in,str_out)
    implicit none
    integer(kind=int8) :: i_in
    character(len=:),allocatable :: str_out
    character(len=4) :: stmp
    write (stmp,*) i_in
    str_out = trim(adjustl(stmp))
  end subroutine i8_to_string

!> 16 bit integers
  subroutine i16_to_string(i_in,str_out)
    implicit none
    integer(kind=int16) :: i_in
    character(len=:),allocatable :: str_out
    character(len=6) :: stmp
    write (stmp,*) i_in
    str_out = trim(adjustl(stmp))
  end subroutine i16_to_string

!> 32 bit integer
  subroutine i32_to_string(i_in,str_out)
    implicit none
    integer(kind=int32) :: i_in
    character(len=:),allocatable :: str_out
    character(len=12) :: stmp
    write (stmp,*) i_in
    str_out = trim(adjustl(stmp))
  end subroutine i32_to_string

!> 64 bit integer
  subroutine i64_to_string(i_in,str_out)
    implicit none
    integer(kind=int64) :: i_in
    character(len=:),allocatable :: str_out
    character(len=60) :: stmp
    write (stmp,*) i_in
    str_out = trim(adjustl(stmp))
  end subroutine i64_to_string

!> 32 bit float
  subroutine r32_to_string(r_in,str_out)
    implicit none
    real(kind=real32) :: r_in
    character(len=:),allocatable :: str_out
    character(len=80) :: stmp
    write (stmp,*) r_in
    str_out = trim(adjustl(stmp))
  end subroutine r32_to_string

!> 64 bit float
  subroutine r64_to_string(r_in,str_out)
    implicit none
    real(kind=wp) :: r_in
    character(len=:),allocatable :: str_out
    character(len=80) :: stmp
    write (stmp,*) r_in
    str_out = trim(adjustl(stmp))
  end subroutine r64_to_string

!> boolean
  subroutine bool_to_string(bool,str_out)
    implicit none
    logical :: bool
    character(len=:),allocatable :: str_out
    if (bool) then
      str_out = 'true'
    else
      str_out = 'false'
    end if
  end subroutine bool_to_string

!> function variants
  function to_str_i8(i8) result(str_out)
    implicit none
    integer(kind=int8) :: i8
    character(len=:),allocatable :: str_out
    call convert_to_string(i8,str_out)
  end function to_str_i8
  function to_str_i16(i16) result(str_out)
    implicit none
    integer(kind=int16) :: i16
    character(len=:),allocatable :: str_out
    call convert_to_string(i16,str_out)
  end function to_str_i16
  function to_str_i32(i32) result(str_out)
    implicit none
    integer(kind=int32) :: i32
    character(len=:),allocatable :: str_out
    call convert_to_string(i32,str_out)
  end function to_str_i32
  function to_str_i64(i64) result(str_out)
    implicit none
    integer(kind=int64) :: i64
    character(len=:),allocatable :: str_out
    call convert_to_string(i64,str_out)
  end function to_str_i64
  function to_str_r32(r32) result(str_out)
    implicit none
    real(kind=real32) :: r32
    character(len=:),allocatable :: str_out
    call convert_to_string(r32,str_out)
  end function to_str_r32
  function to_str_r64(r64) result(str_out)
    implicit none
    real(kind=wp) :: r64
    character(len=:),allocatable :: str_out
    call convert_to_string(r64,str_out)
    str_out = truncate_zeros(str_out)
  end function to_str_r64
  function to_str_bool(bool) result(str_out)
    implicit none
    logical :: bool
    character(len=:),allocatable :: str_out
    call convert_to_string(bool,str_out)
  end function to_str_bool

!> truncate zeros of a float string
  function truncate_zeros(str_in) result(str_out)
    character(len=:),allocatable :: str_in
    character(len=:),allocatable :: str_out
    integer :: i,j,k,l
    str_out = ''
    if (allocated(str_in)) then
      if (index(lowercase(str_in),'e') .ne. 0) then
        str_out = trim(adjustl(str_in))
        return
      end if
      l = len_trim(str_in)
      k = l
      do i = l,1,-1
        k = i
        if (str_in(i:i) .ne. '0') then
          exit
        end if
      end do
      if (str_in(k:k) .eq. '.') k = k+1
      str_out = str_in(1:k)
      str_out = trim(adjustl(str_out))
    end if
  end function truncate_zeros

!========================================================================================!
!========================================================================================!
!========================================================================================!

  subroutine getpath(fname,path)
    implicit none
    character(len=*) :: fname
    character(len=*) :: path
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: rcode,ich,io

    pipe = ' >/dev/null 2>/dev/null'

    checkcall = 'command -v '//trim(fname)//pipe
    call command(checkcall,rcode)

    if (rcode .eq. 0) then
      checkcall = 'command -v '//trim(fname)//' > pathout.tmp 2>/dev/null'
      call command(checkcall,io)
      open (newunit=ich,file='pathout.tmp')
      read (ich,'(a)') path
      close (ich,status='delete')
    end if
    return
  end subroutine getpath
!=========================================================================================!
  subroutine checkprog(fname,r)
    implicit none
    character(len=*) :: fname
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: rcode,r
    character(len=512) :: path

    pipe = ' >/dev/null 2>/dev/null'

    checkcall = 'command -v '//trim(fname)//pipe
    call command(checkcall,rcode)

    write (stdout,'(4x,a,a,a)') 'binary: "',trim(fname),'"'
    if (rcode .ne. 0) then
      write (stdout,'(4x,a)') 'status: not found'
      r = r+1
    else
      write (stdout,'(4x,a)') 'status: present'
      call getpath(fname,path)
      write (stdout,'(4x,a,a)') 'path  : ',trim(path)
    end if

    return
  end subroutine checkprog
!=========================================================================================!
! the same as above but only provide a printout (at all) if the program is not present!
  subroutine checkprog_silent(fname,verbose,iostat)
    implicit none
    character(len=*) :: fname
    logical,intent(in),optional :: verbose
    integer,intent(out),optional :: iostat
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: io

    pipe = ' >/dev/null 2>/dev/null'
    checkcall = 'command -v '//trim(fname)//pipe
    call command(trim(checkcall),exitstat=io)

    if (present(verbose)) then
      if (io .ne. 0.and.verbose) then
        write (stdout,'(4x,a,a,a)') 'binary: "',trim(fname),'"'
        write (stdout,'(4x,a)') 'status: not found!'
      end if
    end if

    if (present(iostat)) then
      iostat = io
    end if

    return
  end subroutine checkprog_silent

!========================================================================================!

  function absolute_path(path_in) result(abs_path)
    use iso_c_binding
    implicit none
    character(len=*),intent(in)  :: path_in
    character(len=:),allocatable :: abs_path
    character(kind=c_char),allocatable,target :: c_path(:)
    character(kind=c_char),target             :: c_buffer(4096)
    type(c_ptr) :: ret_ptr
    integer :: n,i

    !> Trim empty
    if (len_trim(path_in) == 0) then
      allocate (character(len=0) :: abs_path)
      return
    end if

    !> Convert Fortran string to C string (null-terminated)
    n = len_trim(path_in)
    allocate (c_path(n+1))
    c_path(1:n) = transfer(path_in(1:n),c_path(1:n))
    c_path(n+1) = c_null_char

    !> realpath call: use c_loc on the entire arrays, not on slices
    ret_ptr = c_realpath(c_loc(c_path),c_loc(c_buffer))

    !> Check success/failure
    if (.not.c_associated(ret_ptr)) then
      !> realpath failed => return original
      allocate (character(len=n) :: abs_path)
      abs_path = path_in(1:n)
      return
    end if

    !> Find length of C string in c_buffer
    i = 1
    do while (i <= size(c_buffer).and.c_buffer(i) /= c_null_char)
      i = i+1
    end do
    i = i-1
    allocate (character(len=i) :: abs_path)
    abs_path = transfer(c_buffer(1:i),abs_path)
  end function absolute_path

  function absolute_filepath(path_in) result(abs_file)
!****************************************************************************
!* Returns an absolute path even if the final filename does not exist.
!* Splits path_in into directory + file, calls absolute_path(directory),
!* then re-appends the filename.
!* If the input is *just* a file name without any separators, it is returned
!****************************************************************************
    implicit none
    character(len=*),intent(in) :: path_in
    character(len=:),allocatable :: abs_file
    character(len=:),allocatable :: dir_part,file_part,abs_dir
    logical :: has_slash

    !> If empty => return empty
    if (len_trim(path_in) == 0) then
      allocate (character(len=0)::abs_file)
      return
    end if

    !> Split the path into directory part + file part
    call split_path(path_in,dir_part,file_part,has_slash)

    !> If directory part is empty, we don't generate an absolute path!
    if (len_trim(dir_part) == 0) then
      abs_dir = '' !absolute_path(".")
    else
      abs_dir = absolute_path(dir_part)
    end if

    !> If absolute_path failed => just return path_in as is
    if (len_trim(abs_dir) == len_trim(dir_part).and. &
        trim(abs_dir) == trim(dir_part)) then
      !> realpath didn't improve or directory doesn't exist
      allocate (character(len=len_trim(path_in)) :: abs_file)
      abs_file = path_in
      return
    end if

    !> Re-append the file part if any
    if (len_trim(file_part) == 0) then
      !> It's just a directory
      allocate (character(len=len_trim(abs_dir)) :: abs_file)
      abs_file = abs_dir
    else
      !> Ensure abs_dir ends with a slash
      if (abs_dir(len_trim(abs_dir):len_trim(abs_dir)) /= "/") then
        abs_dir = trim(abs_dir)//"/"
      end if

      allocate (character(len=len_trim(abs_dir)+len_trim(file_part)) :: abs_file)
      abs_file = trim(abs_dir)//trim(file_part)
    end if

  end function absolute_filepath

  subroutine split_path(fullpath,dir_part,base_part,has_slash)
!*************************************************************************
!* split_path(fullpath, dir_part, base_part, has_slash):
!* Splits "fullpath" into "dir_part" and "base_part".
!* e.g.  "some/dir/file.txt" => dir_part="some/dir", base_part="file.txt"
!*        "/tmp/" => dir_part="/tmp/", base_part=""
!*        "filename" => dir_part="", base_part="filename"
!*************************************************************************
    implicit none
    character(len=*),intent(in) :: fullpath
    character(len=:),allocatable,intent(out) :: dir_part,base_part
    logical,intent(out) :: has_slash

    integer :: i,n,last_slash
    character(len=:),allocatable :: p

    n = len_trim(fullpath)
    allocate (character(len=n)::p)
    p = fullpath(1:n)

    last_slash = 0
    do i = n,1,-1
      if (p(i:i) == "/") then
        last_slash = i
        exit
      end if
    end do

    if (last_slash == 0) then
      !> No slash at all
      has_slash = .false.
      allocate (character(len=0)::dir_part)
      allocate (character(len=n)::base_part)
      base_part = p
    else
      has_slash = .true.
      allocate (character(len=last_slash)::dir_part)
      dir_part = p(:last_slash)
      if (last_slash < n) then
        allocate (character(len=n-last_slash)::base_part)
        base_part = p(last_slash+1:)
      else
        !> slash is at the end => no file part
        allocate (character(len=0)::base_part)
      end if
    end if

  end subroutine split_path

!========================================================================================!
!========================================================================================!
!========================================================================================!

!> For some reason the behaviour of "call system"  and "call execute_command_line"
!> differs slightly beteen Intel and GNU versions of the program that I have build.
!> I don't understand why this is the case.
!> To circumvent issues for now, I will simply include two different wrappers with
!> precompiler statements.
!> It doesn't matter too much as these are mostly relevant for legacy routines.
!>-----------------------------------------------------------------------------!
!> Scratch root for the per-call calculation directories.
!>
!> CREST runs every single point inside its own directory "calculation.level.<n>"
!> **below the current working directory**.  On a networked / spinning-disk
!> workspace that turns every scratch write of the external program into a
!> network round trip, which easily dominates the runtime of a cheap level of
!> theory.  CREST_SCRATCH_DIR relocates those directories to node-local storage:
!>
!>     export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
!>
!> It is created with "mkdir -p": iomod%makedir() is a bare POSIX mkdir() and
!> would NOT create the parents.  The result is returned WITHOUT a separator
!> when unset and WITH exactly one trailing '/' when set, so callers can simply
!> prepend it to the directory name.
!>
!> NOTE: the root must be unique per job - two concurrent CREST runs sharing
!> one root would both use ".../calculation.level.1".
!>-----------------------------------------------------------------------------!
  subroutine crest_scratch_dir(dir)
    implicit none
    character(len=:),allocatable,intent(out) :: dir
    character(len=512) :: val
    integer :: l,io,estat
    dir = ''
    val = ''
    !> F2003 standard intrinsic (getenv is only a vendor extension)
    call get_environment_variable('CREST_SCRATCH_DIR',value=val,status=estat)
    if (estat /= 0) return
    if (len_trim(val) == 0) return
    dir = trim(val)
    l = len_trim(dir)
    do while (l > 1 .and. dir(l:l) == '/')
      dir = dir(1:l-1)
      l = l-1
    end do
    !> "/" alone (or a string of nothing but slashes) would give the nonsense
    !> path "//calculation.level.N" -> fall back to the default location.
    if (len_trim(dir) == 0 .or. trim(dir) == '/') then
      dir = ''
      return
    end if
    dir = trim(dir)//'/'
    call command('mkdir -p "'//trim(dir)//'"',io)
    if (io /= 0) then
      !> without this the later mkdir(calculation.level.N) fails as well and
      !> the run dies with the unhelpful "Initial geometry optimization failed!"
      write (stdout,'(1x,a,a)') '**WARNING** could not create CREST_SCRATCH_DIR: ',trim(dir)
    end if
  end subroutine crest_scratch_dir

!>-----------------------------------------------------------------------------!
!> Directory for the "cold" intermediate files: written only a handful of times
!> per run, so they gain nothing from node-local storage - but they are of no
!> interest to the user either, so they should not clutter the workspace.
!>
!>   export CREST_WORK_DIR=crest_work    # default
!>   export CREST_WORK_DIR=none          # -> upstream layout, no relocation
!>
!> Same contract as crest_scratch_dir(): '' when disabled, otherwise "<name>/"
!> with exactly one trailing '/'.
!>-----------------------------------------------------------------------------!
  subroutine crest_work_dir(dir)
    implicit none
    character(len=:),allocatable,intent(out) :: dir
    character(len=512) :: val
    integer :: l,io,estat
    dir = ''
    val = ''
    call get_environment_variable('CREST_WORK_DIR',value=val,status=estat)
    if (estat /= 0 .or. len_trim(val) == 0) val = 'crest_work'
    if (trim(val) == 'none' .or. trim(val) == '0') return
    dir = trim(val)
    l = len_trim(dir)
    do while (l > 1 .and. dir(l:l) == '/')
      dir = dir(1:l-1)
      l = l-1
    end do
    if (len_trim(dir) == 0 .or. trim(dir) == '/') then
      dir = ''
      return
    end if
    dir = trim(dir)//'/'
    call command('mkdir -p "'//trim(dir)//'"',io)
    if (io /= 0) then
      write (stdout,'(1x,a,a)') '**WARNING** could not create CREST_WORK_DIR: ',trim(dir)
      dir = ''
    end if
  end subroutine crest_work_dir

!>-----------------------------------------------------------------------------!
!> Relocate CREST's intermediate files out of the workspace.
!>
!> Two classes, matching how often they are touched:
!>   HOT  -> CREST_SCRATCH_DIR (node-local): trajectories, restart dumps, ...
!>   COLD -> CREST_WORK_DIR (subdirectory of the workspace): coord, wbo, ...
!> The user-facing results (crest_conformers.xyz, crest_rotamers.xyz, the
!> ensemble files, ...) are NOT touched and stay where they are.
!>
!> IMPLEMENTATION: symbolic links, not code changes at the open() sites.
!> CREST hard-codes these names in dozens of places ("coord" alone appears ~50
!> times), so rewriting every call site would be a large and unsafe change.
!> A link "coord -> <target>/coord" makes the existing open(file='coord') write
!> to the target while every reader stays untouched.  The decisive advantage is
!> that the failure mode is benign: if a name is missing from the tables below,
!> or if CREST removes a link, the file simply ends up in the workspace again -
!> the run stays correct, it just loses the speed-up / tidiness.
!>
!>   export CREST_STAGE_FILES=0   # escape hatch: disable all of this
!>-----------------------------------------------------------------------------!
  subroutine crest_stage_files()
    implicit none
    character(len=:),allocatable :: hot,work
    character(len=512) :: val
    integer :: i,estat,nhot,ncold
    !>--- written every few MD steps -> node-local scratch
!>--- LOCAL CHANGE (2026-09-27): the list used to hold all five temp directory
!> names CREST has ever used (MDFILES, METADYN, NORMMD, MRMSD, OPTIM - see the
!> array literal removed here), created eagerly at start-up regardless of the
!> run type.  Three of them could never do anything:
!>   MRMSD   - has NO writer left anywhere in this tree.  The only other
!>             occurrence is cleanup.f90's rmrf('MRMSD'), i.e. a leftover from
!>             an older version that removed the directory it no longer fills.
!>   METADYN - the legacy V2 code names its directories METADYN<i>
!>   NORMMD    (confscript2_misc.f90:361 and normmd.f90:233 both do
!>             `write(tmppath,'(a,i0)') 'METADYN',vz`, and set_trj_origins()
!>             appends the same index), and it creates them itself via
!>             `rmrf(tmppath)` + `makedir(trim(tmppath))`
!>             (confscript2_misc.f90:266-267, normmd.f90:48-49).
!>             A symlink of the BARE name can therefore never be matched.
!> What is kept, and why:
!>   MDFILES - the modern MD/MTD worker trajectories (parallel.f90:1269,
!>             `mddats(i)%trajectoryfile = mdir//sep//...`).  It is additionally
!>             relocated by PATH via crest_hot_dir(), which does not depend on
!>             this symlink at all.
!>   OPTIM   - the legacy optimisation directory, genuinely written by
!>             confopt.f90:316 / confscript2_misc.f90:453 / reactor.f90:695 and
!>             the -protonate / -tautomerize / -deprotonate run types
!>             (rename('OPTIM/opt.xyz',...)).
!> Removing the three dead names is behaviour-neutral: nothing ever wrote into
!> them, so the run is identical, there are just three fewer empty directories
!> plus symlinks in the workspace and the scratch root.
    character(len=32),parameter :: hotl(2) = [character(len=32) :: &
    &   'MDFILES','OPTIM']
    !>--- written a few times per run -> subdirectory of the workspace
    character(len=32),parameter :: coldl(7) = [character(len=32) :: &
    &   'coord','wbo','crestopt.xyz','gfnff_adjacency','gfnff_topo', &
    &   '.CHRG','.UHF']

    call get_environment_variable('CREST_STAGE_FILES',value=val,status=estat)
    if (estat == 0) then
      if (trim(val) == '0' .or. trim(val) == 'none') return
    end if

    call crest_scratch_dir(hot)
    call crest_work_dir(work)
    if (len(hot) == 0 .and. len(work) == 0) return

    nhot = size(hotl)
    ncold = size(coldl)
    do i = 1,nhot
      if (len(hot) > 0) call crest_link(hot,trim(hotl(i)),.true.)
    end do
    do i = 1,ncold
      if (len(work) > 0) call crest_link(work,trim(coldl(i)),.false.)
    end do
  end subroutine crest_stage_files

!>-----------------------------------------------------------------------------!
!> Put the symlink "name -> <root><name>" in the working directory.
!> "isdir" also creates the target directory itself.
!> Anything that already exists is left strictly alone.
!>-----------------------------------------------------------------------------!
  subroutine crest_link(root,name,isdir)
    implicit none
    character(len=*),intent(in) :: root
    character(len=*),intent(in) :: name
    logical,intent(in) :: isdir
    character(len=1024) :: cmd
    integer :: io
    logical :: ex

    inquire (file=trim(name),exist=ex)
    if (ex) return
    if (isdir) then
      call command('mkdir -p "'//trim(root)//trim(name)//'"',io)
      if (io /= 0) return
    else
      call command('mkdir -p "'//trim(root)//'"',io)
      if (io /= 0) return
    end if
    !> a dangling link from an earlier run is not "existing" for inquire(), so
    !> remove whatever is there first - never a real file, that was caught above
    cmd = 'rm -f "'//trim(name)//'" ; ln -s "'// &
    &     trim(root)//trim(name)//'" "'//trim(name)//'"'
    call command(trim(cmd),io)
  end subroutine crest_link

!>-----------------------------------------------------------------------------!
!> Hot variant of a scratch subdirectory: "<CREST_SCRATCH_DIR><name>"
!> (or just <name> when no scratch directory is configured).
!>-----------------------------------------------------------------------------!
  subroutine crest_hot_dir(dir,name)
    implicit none
    character(len=:),allocatable,intent(out) :: dir
    character(len=*),intent(in) :: name
    character(len=:),allocatable :: root
    call crest_scratch_dir(root)
    dir = trim(root)//trim(name)
  end subroutine crest_hot_dir

!>-----------------------------------------------------------------------------!
  subroutine command(cmd,exitstat)
    implicit none
    character(len=*),intent(in) :: cmd
    integer,intent(out),optional :: exitstat
    integer :: io
    interface
      !> SIGINT-safe system(): see src/signal.c.  The Fortran SYSTEM intrinsic
      !> and EXECUTE_COMMAND_LINE both end up in POSIX system(3), which is
      !> REQUIRED to ignore SIGINT and SIGQUIT while it waits for the child.
      !> CREST is inside this routine for ~100 % of its run time (waiting for
      !> orca/xtb), so with system(3) Ctrl+C could never terminate CREST -
      !> it only killed the child and CREST sat in waitpid() forever.
      !> crest_runcmd() is fork()+execl("/bin/sh","sh","-c",cmd)+waitpid()
      !> and leaves the signal dispositions alone.
      function crest_runcmd(ccmd) bind(C,name="crest_runcmd") result(stat)
        use iso_c_binding
        character(kind=c_char),dimension(*),intent(in) :: ccmd
        integer(c_int) :: stat
      end function crest_runcmd
    end interface
    !> the external program reads files that were just written: make sure the
    !> console is not still holding a partially filled buffer
    flush (stdout)
    io = crest_runcmd(trim(cmd)//c_null_char)
    if (present(exitstat)) exitstat = io
  end subroutine command

!========================================================================================!
!========================================================================================!
!========================================================================================!
end module iomod

