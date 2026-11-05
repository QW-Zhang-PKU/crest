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

!> module orca_type
!> a minimal implementation for storing an ORCA input file

module orca_type
  use crest_parameters,only:wp,stdout,stderr,autoaa
  use iomod
  use strucrd
  implicit none
  public

  !> Maximum length of one ORCA template line.  This is the same size as the
  !> scratch buffer the reader below reads a line into, so a template line can
  !> never exceed it.
  integer,parameter :: orca_linelen = 1056

  !>--- NOTE on the element type of "input":
  !> "input" is an ALLOCATABLE ARRAY OF FIXED-LENGTH STRINGS on purpose - do NOT
  !> turn it into "character(len=:),allocatable :: input(:)".
  !> A deferred-length character ARRAY component is silently mis-copied when the
  !> surrounding derived type is assigned (gcc PR fortran/89200, a GCC 9
  !> regression; ifort/ifx have the same trap).  After the copy the elements no
  !> longer reference the template but stale heap memory.
  !> CREST hits this because calculation_settings objects get copied around
  !> (tmpcalc = env%calc in setuptest.f90, calc%add() while parsing the TOML),
  !> and this is the ONLY deferred-length character array component in the whole
  !> code base - which is why the damage is confined to the ORCA template.
  !> Symptom: ORCA.in is written full of NUL bytes plus fragments of CREST's own
  !> memory (the banner, the jobdescription() table, the systemcall string).
  !> ORCA then aborts with
  !>   ERROR: non-ASCII character(s) found in words in line ...
  !>   [qcscan1.cpp, line 264]: Please remove all non-ASCII characters
  !> and CREST only reports the opaque "Initial geometry optimization failed!".
  !> Fixed-length elements are copied by a plain memcpy and cannot be mangled.
  type :: orca_input
    character(len=:),allocatable :: cmd
    integer :: nlines = 0
    character(len=orca_linelen),allocatable :: input(:)
    logical :: mpi = .false.
  contains
    procedure :: read => read_orca_input
    procedure :: write => write_orca_input
  end type orca_input

!>=======================================================================================
!> "engrad" handling.  Three things need it:
!>   * read_orca_input()  - must RECOGNISE an existing EnGrad so it does not
!>     append a second one (case-insensitively!)
!>   * write_orca_input() - must be able to REMOVE it again when the caller only
!>     wants an energy (singlepoint(), i.e. calc%nogr)
!> Both work on the raw keyword line, so the matching has to be case-insensitive:
!> ORCA spells it "EnGrad" in its own manual and accepts any capitalisation.
!>=======================================================================================
  character(len=6),parameter :: ENG = 'engrad'

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!---------------------------------------------------------------------------------------!
!> is the (sub)string "engrad" present in line, regardless of capitalisation?
!---------------------------------------------------------------------------------------!
  pure function has_engrad(line) result(found)
    implicit none
    character(len=*),intent(in) :: line
    logical :: found
    integer :: i,n
    found = .false.
    n = len_trim(line)
    do i = 1,n-len(ENG)+1
      if (locase(line(i:i+len(ENG)-1)) == ENG) then
        found = .true.
        return
      end if
    end do
  end function has_engrad

!---------------------------------------------------------------------------------------!
!> lowercase a short string (only A-Z are mapped, everything else is kept)
!---------------------------------------------------------------------------------------!
  pure function locase(str) result(out)
    implicit none
    character(len=*),intent(in) :: str
    character(len=len(str)) :: out
    integer :: i,ic
    out = str
    do i = 1,len_trim(str)
      ic = iachar(out(i:i))
      if (ic >= iachar('A').and.ic <= iachar('Z')) out(i:i) = char(ic+32)
    end do
  end function locase

!---------------------------------------------------------------------------------------!
!> remove the FIRST "engrad" (any capitalisation) from line, plus one adjacent
!> blank.  Used to turn an EnGrad template into a plain single point template.
!---------------------------------------------------------------------------------------!
  function strip_engrad(line) result(out)
    implicit none
    character(len=*),intent(in) :: line
    character(len=len(line)) :: out
    integer :: i,j,n
    out = line
    n = len(ENG)
    do i = 1,len_trim(out)-n+1
      if (locase(out(i:i+n-1)) == ENG) then
!>--- drop one of the two surrounding blanks as well, so "! PBE EnGrad def2"
!> becomes "! PBE def2" and not "! PBE  def2"
        j = i+n
        if (j <= len_trim(out)) then
          if (out(j:j) == ' ') then
            out = out(:i-1)//out(j+1:)
          else
            out = out(:i-1)//out(j:)
          end if
        else
          out = out(:i-1)
        end if
        return
      end if
    end do
  end function strip_engrad

!========================================================================================!
!========================================================================================!

  subroutine read_orca_input(self,fname)
    implicit none
    class(orca_input) :: self
    character(len=*),intent(in) :: fname
    logical :: ex,trackcoord
    integer :: nlines,width
    integer :: io,ich,i,j,k,l
    character(len=1056) :: atmp
    logical :: gotengrad

    inquire (file=fname,exist=ex)
    if (.not.ex) then
      write (stderr,'(3a)') '**ERROR** ORCA input template ',fname,' could not be found!'
      error stop
    else
      write (stdout,'(3a)') 'Reading ORCA input template ',fname,' ...'
!>--- LOCAL EXTENSION: frame the template with the same '*' banner that the
!> TOML echo in parse_datastruct.f90 (root_print2) uses, instead of the '|'
!> bars we used before.  CREST now echoes THREE kinds of input file (TOML,
!> every ORCA template, and the constraining file) and they all have to look
!> alike, otherwise the log is unreadable.
      write (stdout,'(a)') repeat('*',80)
      write (stdout,'("*",1x,a,a,a)') 'INPUT FILE ',trim(fname),' content (without comments):'
      write (stdout,'(a)') repeat('*',80)
    end if
    nlines = getlines(fname,width)
    open (newunit=ich,file=fname)
    nlines = 0
    trackcoord = .false.
    gotengrad = .false.
    !> count lines and check keywords
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit  !> EOF
!>--- LOCAL EXTENSION: echo the template verbatim into the log.
!> Upstream only said "Reading ORCA input template orca.inp ...", so a log file
!> never contained the method/basis/convergence settings that were actually
!> used - and for a multi-level setup (orca.inp vs. orca-high.inp) there was no
!> way to tell afterwards which file held which keywords.  Dumped BEFORE the
!> adjustl()+lowercase() below, i.e. exactly as it is on disk.
      write (stdout,'("*",1x,a)') trim(atmp)
      atmp = adjustl(lowercase(atmp))

      !> ignore coord lines, CREST will write those
      if (trackcoord.and.atmp(1:1) == '*') then
        trackcoord = .false.
        cycle
      end if
      if (atmp(1:1).eq.'*' .and. index(atmp,'xyz ') .ne. 0) then
        trackcoord = .true.
        cycle
      end if
      if (atmp(1:1).eq.'*' .and. index(atmp,'xyzfile') .ne. 0) cycle

      if (.not.trackcoord) nlines = nlines+1

      if (index(atmp,'$new_job') .eq. 1) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': please define only single jobs (no $new_job)!'
        error stop
      end if

      !> check runtypes from the simple input line
!>--- LOCAL FIX: index() is CASE SENSITIVE, so a template written as "EnGrad"
!> (the spelling in the ORCA manual) was not recognised here and CREST appended
!> a SECOND EnGrad below -> "! ... EnGrad EnGrad".  Harmless for ORCA but it
!> makes the generated ORCA.in contradict the template the user wrote.
      if (atmp(1:1) .eq. '!'.and.has_engrad(atmp)) then
        gotengrad = .true.
      end if
      if (atmp(1:1) .eq. '!'.and.index(atmp,'md') .ne. 0) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': only run EnGrad jobs!'
        error stop
      end if
      if (atmp(1:1) .eq. '!'.and.index(atmp,'opt') .ne. 0) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': only run EnGrad jobs!'
        error stop
      end if
      if (atmp(1:1) .eq. '!'.and.index(atmp,'freq') .ne. 0) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': only run EnGrad jobs!'
        error stop
      end if

      !> coordinate input check, e.g.  remove Bohr keyword, if necessary
      if (atmp(1:1) .eq. '%'.and.index(atmp,'coords') .ne. 0) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': please remove %coords block!'
        error stop
      end if

      !> check if there is parallelization
      if (atmp(1:1) .eq. '!'.and.index(atmp,'pal') .ne. 0) then
        self%mpi = .true.
      end if
      if (atmp(1:1) .eq. '%'.and.index(atmp,'pal') .ne. 0) then
        self%mpi = .true.
      end if
    end do
    close (ich)
    write (stdout,'(a)') repeat('*',80)

    !> allocate memory space.  The element length is FIXED (orca_linelen) - see
    !> the note on type(orca_input): a deferred-length character array component
    !> is corrupted by the derived-type copy, which is what produced a binary
    !> ORCA.in.
    !> The SOURCE= length MUST equal orca_linelen: with a shorter SOURCE,
    !> gfortran <= 13 sizes the elements by the SOURCE length instead of the
    !> declared element length (gcc PR fortran/113793), so each later
    !> "self%input(k) = trim(atmp)" writes past the element and the next free()
    !> aborts with "double free or corruption (out)".  GCC 14 fixed this.
    allocate (self%input(nlines),source=repeat(" ",orca_linelen))
    self%nlines = nlines

    !> Open file from the beginning and read into memory
    open (newunit=ich,file=fname)
    k = 0
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit  !> EOF
      atmp = adjustl(atmp)

      !> ignore coord lines, CREST will write those
      if (trackcoord.and.atmp(1:1) == '*') then
        trackcoord = .false.
        cycle
      end if
      if (atmp(1:1).eq.'*' .and. index(atmp,'xyz ') .ne. 0) then
        trackcoord = .true.
        cycle
      end if
      if (atmp(1:1).eq.'*' .and. index(atmp,'xyzfile') .ne. 0) cycle

      if (atmp(1:1) .eq. '!'.and.index(atmp,'bohr') .ne. 0) then
        j = index(atmp,'bohr')
        atmp(j:j+3) = '    '
      end if

      k = k+1
      if (atmp(1:1) .eq. '!'.and..not.gotengrad) then
        atmp = trim(atmp)//' EnGrad'
        gotengrad = .true.
      end if
      self%input(k) = trim(atmp)

    end do
    close (ich)
  end subroutine read_orca_input

!========================================================================================!
!> How many MPI ranks the ORCA child should use
!========================================================================================!

  !>---------------------------------------------------------------------------
  !> The core budget CREST has reserved for THIS call, translated into ORCA's
  !> own parallelisation keyword.
  !>
  !> new_ompautoset() (ompmklset.F90) hands "cores_per_job" to every EXTERNAL
  !> program through OMP_NUM_THREADS, and it already IS the fair share:
  !>   'max'         -> T       : trialOPT (initial optimisation), singlepoints
  !>   'auto',nall   -> T/nall  : ensemble optimisations, MTD/MD workers
  !> ORCA however parallelises over MPI and needs an explicit
  !>     %pal
  !>       nprocs N
  !>     end
  !> block - it never looks at OMP_NUM_THREADS, which is exactly why an
  !> unpatched CREST+ORCA run uses ONE core regardless of -T.
  !> Override with CREST_ORCA_NPROCS (e.g. to cap the ranks per job).
  !>---------------------------------------------------------------------------
  integer function orca_nprocs() result(np)
    implicit none
    character(len=64) :: val
    integer :: r
    call getenv('CREST_ORCA_NPROCS',val)
    read (val,*,iostat=r) np
    if (r /= 0) then
      call getenv('OMP_NUM_THREADS',val)
      read (val,*,iostat=r) np
      if (r /= 0) np = 1
    end if
    np = max(np,1)
  end function orca_nprocs

!========================================================================================!
!>---------------------------------------------------------------------------
!> How many OpenMP/MKL threads the ORCA *child* gets.
!>
!> NOTE: this is deliberately NOT tied to the number of MPI ranks.  Many ranks
!> are only worth it if the per-iteration work is large; for a minimal basis
!> (MINIX) and a few hundred basis functions, 8 MPI processes usually LOSE
!> against 1 process with 8 threads - every SCF iteration pays an Allreduce of
!> the Fock matrix, which is pure overhead at that size.
!>
!> CREST_ORCA_OMP_THREADS=N  -> export OMP/MKL/OPENBLAS_NUM_THREADS=N
!>                              for the ORCA child (default: 1)
!>---------------------------------------------------------------------------
  subroutine orca_omp_threads_string(str)
    implicit none
    character(len=*),intent(out) :: str
    character(len=64) :: val
    integer :: n,r
    n = 1
    call getenv('CREST_ORCA_OMP_THREADS',val)
    read (val,*,iostat=r) n
    if (r /= 0) n = 1
    n = max(n,1)
    write (str,'(a,i0,a,i0,a,i0)') &
    &  ' OMP_NUM_THREADS=',n,' MKL_NUM_THREADS=',n,' OPENBLAS_NUM_THREADS=',n
  end subroutine orca_omp_threads_string

!========================================================================================!

  subroutine write_orca_input(self,fname,mol,chrg,mult,nogr)
    implicit none
    class(orca_input),intent(in) :: self
    character(len=*),intent(in) :: fname
    type(coord),intent(in) :: mol
    integer,intent(in) :: chrg
    integer,intent(in) :: mult
!>--- LOCAL EXTENSION: nogr = .true. writes the template WITHOUT the EnGrad
!> keyword, i.e. turns it into a plain single point.  Used by ORCA_sp() when the
!> caller (singlepoint(), i.e. the refine singlepoint re-ranking) only needs the
!> energy: the gradient costs a CP-SCF/Z-vector solve (plus the CPCM gradient
!> term) that nobody ever reads in that stage.
!> OPTIONAL on purpose: every existing call site keeps the old behaviour.
    logical,intent(in),optional :: nogr
    integer :: ich,i,j,k,l
    integer :: np
    logical :: sponly
    character(len=orca_linelen) :: line

    if(.not.allocated(self%input))then
      write (stderr,'(3a)') '**ERROR** Please provide an ORCA input template!'
      error stop
    endif

    sponly = .false.
    if (present(nogr)) sponly = nogr

    np = orca_nprocs()

    open (newunit=ich,file=fname)
    do i=1,self%nlines
      line = self%input(i)
      if (sponly.and.has_engrad(line)) line = strip_engrad(line)
      write(ich,'(a)') trim(line)
    enddo
    !>--- MPI parallelisation.  Only injected when the template does not already
    !> carry a %pal block (read_orca_input sets self%mpi when it sees one), so an
    !> explicit %pal in orca.inp always wins.
    if (np > 1 .and. .not.self%mpi) then
      write(ich,'(a)') '%pal'
      write(ich,'(a,1x,i0)') '  nprocs',np
      write(ich,'(a)') 'end'
    end if
    write(ich,*)
    write(ich,'(a,1x,i0,1x,i0,a)') '*xyz',chrg,mult,'  # charge and multiplicity (2S+1)'
    do i=1,mol%nat
      write(ich,'(a2,3F25.15)') asym(mol%at(i)),mol%xyz(1:3,i)*autoaa
    enddo
    write(ich,'("*")')
    close (ich)

  end subroutine write_orca_input

!========================================================================================!
!========================================================================================!
end module orca_type
