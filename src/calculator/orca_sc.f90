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

!> module orca_sc
!> A module containing routines for ORCA subrpocess calls

!=========================================================================================!
module orca_sc
  use iso_fortran_env,only:wp => real64,stderr => error_unit
  use strucrd
  use calc_type
  use orca_type,only:orca_nprocs,orca_omp_threads_string
  use iomod,only:makedir,directory_exist,remove,command
  implicit none
  !>--- private module variables and parameters
  private
  integer,parameter :: nf = 2
  character(len=*),parameter :: orcafiles(nf) = [&
          & 'ORCA.in    ','ORCA.engrad']
  character(len=3),parameter :: orca = 'orca'
  character(len=8),parameter :: xyzn = 'ORCA.in'
  character(len=11),parameter :: gf = 'ORCA.engrad'

  public :: ORCA_engrad
  public :: ORCA_sp

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine ORCA_engrad(mol,calc,energy,grad,iostatus)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove

    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    integer :: i,j,k,l,ich,och,io
    logical :: ex

    iostatus = 0

    !>--- setup system call information
    !$omp critical
    call ORCA_setup(mol,calc)
    !$omp end critical

    !>--- do the systemcall
    call initsignal()
    call command(calc%systemcall,iostatus)
    if (iostatus /= 0) return

    !>--- read energy and gradient
    !$omp critical
    call rd_ORCA_engrad(mol,calc,energy,grad,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    !>--- read WBOs?
    !$omp critical
    call rd_ORCA_wbo(mol,calc,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    return
  end subroutine ORCA_engrad

!========================================================================================!
!========================================================================================!
!> ORCA_sp: the SAME run, but WITHOUT the gradient.
!>
!> WHY: crest_sploop() - the "Singlepoint re-ranking for N structures" stage of
!> the refine - only ever uses the ENERGY.  It calls engrad(), ORCA therefore
!> writes an EnGrad job, and the gradient it computes is thrown away: the
!> gradient array in crest_sploop is filled and then deallocated without ever
!> being read.  For a DFT level with a solvation model that is a CP-SCF /
!> Z-vector solve plus the CPCM gradient contribution, on every single structure
!> of the ensemble, for nothing.
!>
!> HOW: ORCA_setup() writes the template with the EnGrad keyword stripped, and
!> the energy is read back from the "FINAL SINGLE POINT ENERGY" line of ORCA.out
!> instead of from ORCA.engrad.  Everything else (calculation space, system
!> call, MPI/%pal block, WBO reading) is identical to ORCA_engrad().
!========================================================================================!
  subroutine ORCA_sp(mol,calc,energy,iostatus)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove

    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    integer,intent(out) :: iostatus

    iostatus = 0

    !>--- setup system call information (input file WITHOUT EnGrad)
    !$omp critical
    call ORCA_setup(mol,calc,nogr=.true.)
    !$omp end critical

    !>--- do the systemcall
    call initsignal()
    call command(calc%systemcall,iostatus)
    if (iostatus /= 0) return

    !>--- read the energy from ORCA.out
    !$omp critical
    call rd_ORCA_energy(mol,calc,energy,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    !>--- read WBOs? (only if requested; see rd_ORCA_wbo)
    !$omp critical
    call rd_ORCA_wbo(mol,calc,iostatus)
    !$omp end critical
    if (iostatus /= 0) return

    return
  end subroutine ORCA_sp

!========================================================================================!

  subroutine ORCA_setup(mol,calc,nogr)
!***********************************
!* Write ORCA inputs, and construct
!* subprocess call
!***********************************
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
!>--- LOCAL EXTENSION: nogr = .true. writes the template WITHOUT the EnGrad
!> keyword (see ORCA_sp).  Optional: the existing call passes nothing and keeps
!> writing EnGrad exactly as before.
    logical,intent(in),optional :: nogr

    integer :: l
    integer :: np
    logical :: sponly
    character(len=:),allocatable :: fname
    character(len=:),allocatable :: cpath
    character(len=10) :: num
    character(len=256) :: ompstr   !> "<N>" thread budget handed to the ORCA child
    character(len=256) :: envstr   !> scratch for reading environment switches
    character(len=256) :: argstr   !> extra args passed straight through to ORCA/MPI
    character(len=512) :: mpiargs  !> ALL extra MPI switches, to be quoted as ONE argv entry
    integer :: i,j,k,ich,och,io,ival
    logical :: ex
    logical,save :: callprinted = .false.  !> CREST_ORCA_PRINT_CALL: print the call once
    logical,save :: optstartprinted = .false.  !> "optimization starting" mark, printed once

    call initsignal()

    !>--- set default binary if not present
    if (.not.allocated(calc%binary)) then
      if(.not.allocated(calc%ORCA%cmd))then
       write(stdout,'(a)') "**WARNING** ORCA cmd line not specified, falling back to 'orca'"
       calc%binary = orca
      else
       calc%binary = calc%ORCA%cmd
      endif
    end if

    !>--- check for the calculation space
    if (allocated(calc%calcspace)) then
      ex = directory_exist(calc%calcspace)
      if (.not.ex) then
        io = makedir(trim(calc%calcspace))
      end if
      cpath = calc%calcspace
    else
      cpath = ''
    end if
    !>--- cleanup old files
    do i = 1,nf
      call remove(trim(cpath)//sep//trim(orcafiles(i)))
    end do
    deallocate (cpath)

    !>--- construct path information and write ORCA.in file
    if (.not.allocated(calc%calcfile)) then
      if (allocated(calc%calcspace)) then
        l = len_trim(calc%calcspace)
        fname = trim(calc%calcspace)
        if (calc%calcspace(l:l) == sep) then
          fname = trim(fname)//xyzn
        else
          fname = trim(fname)//sep//xyzn
        end if
      else
        fname = xyzn
      end if
      calc%calcfile = fname
    else
      fname = calc%calcfile
    end if
    if(calc%uhf < 1) calc%uhf = 1 !> ORCA uses multiplicity, not n_alpha - n_beta!
    sponly = .false.
    if (present(nogr)) sponly = nogr
    call calc%ORCA%write(fname,mol,calc%chrg,calc%uhf,sponly)
    deallocate (fname)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    !>--- if the systemcall was already set up, return
    if (allocated(calc%systemcall)) return
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

    !>--- construct path information for sys-call
    !> ORCA parallelises over MPI (the %pal block that write_orca_input injects).
    !> By default it must NOT additionally pick up the OpenMP/MKL budget CREST
    !> exported for the child: N ranks x N threads would be N*N runnable threads
    !> per job, so the thread layer of the ORCA child is pinned to 1.
    !>
    !> EXCEPTION / tuning knob: for small levels of theory (minimal basis such as
    !> MINIX) MPI hardly pays off - every SCF iteration then spends more time in
    !> Allreduce of the Fock matrix than in actual work.  In that regime
    !>     CREST_ORCA_NPROCS=1 CREST_ORCA_OMP_THREADS=8
    !> (no %pal at all + 8 threads inside one process) can be markedly faster than
    !> 8 MPI ranks.  Try both; the default below stays at 1 thread per rank.
    if (allocated(calc%calcspace)) then
      calc%systemcall = 'cd '//calc%calcspace//' &&'
    else
      calc%systemcall = ''
    end if
    call orca_omp_threads_string(ompstr)
    calc%systemcall = trim(calc%systemcall)//trim(ompstr)
    calc%systemcall = trim(calc%systemcall)//' '//trim(calc%binary)

    calc%systemcall = trim(calc%systemcall)//' ORCA.in'

!>--- LOCAL FIX: ORCA accepts at most TWO command line arguments
!>        orca <input>  "<additional mpi parameters>"
!>    Everything after the input file has to reach the ORCA driver as ONE single
!>    quoted string - that is how ORCA forwards it to the mpirun it spawns
!>    internally (see the ORCA manual: `$(which orca) inp.inp "optional openmpi
!>    arguments"`).  Appending the switches as separate words makes the driver
!>    abort before it ever touches the input file:
!>        "ORCA: too many arguments -
!>         Additional mpi parameters have to be given as a single string!"
!>    which is exactly what CREST_ORCA_MPI_ARGS='--bind-to core' produced
!>    (argv = ORCA.in | --oversubscribe | --bind-to | core  = 3 extra args).
!>    => collect every switch first, then quote the lot exactly once.
    mpiargs = ''
    np = orca_nprocs()
    if((calc%ORCA%mpi .or. np > 1).and.index(calc%binary,'oversubscribe').eq.0)then
      !> OpenMPI disables process binding once it is told to oversubscribe, which
      !> can cost a surprising amount of performance.  CREST_ORCA_NO_OVERSUBSCRIBE
      !> drops the flag if the node has enough free slots anyway.
      call getenv('CREST_ORCA_NO_OVERSUBSCRIBE',envstr)
      ival = 0   !> a failed internal READ leaves the variable UNCHANGED, and
      !> Fortran's .and./.or. do not guarantee short-circuit evaluation, so the
      !> value has to be defined before the test below can look at it
      read (envstr,*,iostat=io) ival
      !> DEFAULT ON.  Only an explicit non-zero value drops the flag; an unset
      !> or unparsable value keeps it.  (The old test was
      !> 'len_trim(envstr) == 0', i.e. ANY non-empty value disabled
      !> --oversubscribe - so CREST_ORCA_NO_OVERSUBSCRIBE=0 did the exact
      !> opposite of what its name says.)
      if (io /= 0.or.ival == 0) then
        mpiargs = '--oversubscribe'
      end if
    endif

    !>--- anything else the site needs passed straight through to ORCA/MPI,
    !> e.g. "--bind-to core" to restore pinning after --oversubscribe.
    !> NOTE: give the switches WITHOUT surrounding quotes now - CREST adds them.
    !> (A value that is already quoted is passed through untouched, so older
    !>  configs like CREST_ORCA_MPI_ARGS='"--bind-to core"' keep working.)
    call getenv('CREST_ORCA_MPI_ARGS',argstr)
    if (len_trim(argstr) > 0) then
      if (len_trim(mpiargs) > 0) then
        mpiargs = trim(mpiargs)//' '//trim(argstr)
      else
        mpiargs = trim(argstr)
      end if
    end if
    if (len_trim(mpiargs) > 0) then
      l = len_trim(mpiargs)
      if (mpiargs(1:1) == '"'.and.mpiargs(l:l) == '"') then
        calc%systemcall = trim(calc%systemcall)//' '//trim(mpiargs)
      else
        calc%systemcall = trim(calc%systemcall)//' "'//trim(mpiargs)//'"'
      end if
    end if

    !>--- add printout information
    !> NOTE: the other subprocess calculators append dev0 (' 2>/dev/null') and
    !> therefore THROW AWAY the child's stderr.  For ORCA that is fatal for
    !> debugging: a broken start-up (LD_LIBRARY_PATH not set -> "error while
    !> loading shared libraries: libmpi.so.40", wrong binary path, license
    !> problem) reports ONLY on stderr and never reaches ORCA.out.  CREST then
    !> just sees a missing ORCA.engrad and aborts with the opaque
    !> "Initial geometry optimization failed!".  Keep the message instead: it is
    !> cheap and it is the first file to look at when ORCA fails.  The file is
    !> created inside calc%calcspace, i.e. next to ORCA.in/ORCA.out.
    calc%systemcall = trim(calc%systemcall)//' '//'> ORCA.out'
    calc%systemcall = trim(calc%systemcall)//' 2> ORCA.err'

!>--- LOCAL EXTENSION: opt-in echo of the command line CREST is about to run.
!> Set CREST_ORCA_PRINT_CALL=1 to see it; 0 or unset suppresses it.  This is the
!> fastest way to settle "what is CREST actually executing?" - a mis-quoted or
!> mis-assembled call is otherwise invisible, because ORCA.out only contains
!> ORCA's own complaints.
!> Printed once per process (the same string is built for every worker anyway).
    call getenv('CREST_ORCA_PRINT_CALL',envstr)
    ival = 0   !> see the note above: a failed READ must not leave ival undefined
    read (envstr,*,iostat=io) ival
    !> DEFAULT OFF, and the switch is the VALUE, not the mere presence of the
    !> variable - the same convention as CREST_DYNAMIC_SHARE in parallel.f90
    !> ('read(val,*,iostat=r); if (r == 0.and.ival == 0) return').  The old test
    !> was 'len_trim(envstr) > 0', which reads any non-empty value as "on", so
    !> CREST_ORCA_PRINT_CALL=0 printed the very line it was meant to suppress.
    !> Only integers are honoured; an unset or unparsable value leaves it off.
    if (io == 0.and.ival /= 0.and..not.callprinted) then
      callprinted = .true.
      write (stdout,'(a)') 'ORCA call: '//trim(calc%systemcall)
      flush (stdout)
    end if
!>--- LOCAL EXTENSION: mark the start of the initial geometry optimization.
!> The FIRST ORCA call of every CREST run is always made from inside that
!> optimization, and it is the one stage that stays completely silent until it
!> is finished (17 min and more with a real DFT level - nothing else in the
!> log tells where the time went).  So print the "starting" line directly
!> behind the echo above; with CREST_ORCA_PRINT_CALL unset it still lands
!> between the "Initial Geometry Optimization" header and the run hashes of
!> ORCA itself, i.e. exactly where the silence used to be.
!> Printed once per process (same reasoning as the echo above).
    if (.not.optstartprinted) then
      optstartprinted = .true.
      write (stdout,'(1x,a)') 'optimization starting ...'
      flush (stdout)
    end if
    return
  end subroutine ORCA_setup

!========================================================================================!

  subroutine rd_ORCA_engrad(mol,calc,energy,grad,iostatus)
!**************************************************
!* Read energy and gradient from ORCA output.
!* ORCA will always write to a .engrad file,
!* which with our naming conventions should always
!* be called ORCA.engrad
!**************************************************
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove
    use gradreader_module,only:rd_grad_engrad

    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus
    integer :: n,c
    real(wp) :: dum
    character(len=128) :: atmp

    integer :: i,j,k,ich,och,io
    logical :: ex

    call initsignal()

    iostatus = 0

    !>--- NOTE: calc%gradfile only holds the file NAME (or a user supplied
    !> name).  The path prefix has to be re-derived from the CURRENT
    !> calc%calcspace on EVERY call - cp. generic_sc.f90.  Caching the full
    !> path breaks parallel runs, where calcspace is switched per worker
    !> thread (see parallel.f90): all workers would read the same stale
    !> gradient file.  See the comment in xtb_sc.f90 for the full story.
    if (.not.allocated(calc%gradfile)) calc%gradfile = gf

    if (allocated(calc%calcspace)) then
      atmp = trim(calc%calcspace)//sep//trim(calc%gradfile)
    else
      atmp = trim(calc%gradfile)
    end if

    inquire (file=atmp,exist=ex)
    if (.not.ex) then
      iostatus = 1
      return
    end if

    c = 0
    open (newunit=ich,file=atmp)
    call rd_grad_engrad(ich,mol%nat,energy,grad,iostatus)
    close (ich)

    return
  end subroutine rd_ORCA_engrad

!========================================================================================!
!> rd_ORCA_energy: read ONLY the energy, from ORCA.out, for a gradient-free run.
!>
!> WHY: with the EnGrad keyword stripped from the template (see write_orca_input
!> and ORCA_sp) ORCA does not write an ORCA.engrad file at all, so there is
!> nothing for rd_ORCA_engrad() to read.  ORCA nevertheless prints the total
!> energy of the converged SCF in the "FINAL SINGLE POINT ENERGY" line of its
!> standard output, and that value is exactly the one it would have put into the
!> .engrad file, so the two are interchangeable.
!>
!> NOTE ON THE PATH: calc%calcspace is switched per worker thread in parallel
!> runs (see parallel.f90), so - exactly as in rd_ORCA_engrad() above and in
!> rd_ORCA_wbo() below - the path prefix is re-derived from the CURRENT
!> calcspace on EVERY call.  Never cache it.
!========================================================================================!
  subroutine rd_ORCA_energy(mol,calc,energy,iostatus)
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type
    use iomod,only:makedir,directory_exist,remove

    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    integer,intent(out) :: iostatus

    character(len=*),parameter :: KEY = 'FINAL SINGLE POINT ENERGY'
    character(len=:),allocatable :: outfile
    character(len=1024) :: atmp
    integer :: ich,io,p,nfound
    logical :: ex
    real(wp) :: dum

    call initsignal()

    iostatus = 0

    !>--- locate the output file of the corresponding ORCA run
    if (allocated(calc%calcspace)) then
      outfile = trim(calc%calcspace)//sep//'ORCA.out'
    else
      outfile = 'ORCA.out'
    end if

    inquire (file=outfile,exist=ex)
    if (.not.ex) then
      write (stderr,'(3a)') '**WARNING** ORCA output ',trim(outfile),&
              & ' not found - cannot read energy'
      iostatus = 1
      return
    end if

    !>--- scan the output.  The LAST occurrence wins: an ORCA run that restarts
    !> or prints several SCFs (e.g. a compound job) repeats the line, and the
    !> final one is the converged result that would have gone into ORCA.engrad.
    nfound = 0
    open (newunit=ich,file=outfile)
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit
      p = index(atmp,KEY)
      if (p == 0) cycle
      read (atmp(p+len(KEY):),*,iostat=io) dum
      if (io /= 0) cycle
      energy = dum
      nfound = nfound+1
    end do
    close (ich)

    if (nfound == 0) then
      write (stderr,'(3a)') '**WARNING** no "',KEY,'" line found in '//trim(outfile)
      write (stderr,'(a)') '            The ORCA run most likely did not converge.'
      iostatus = 1
      return
    end if

    return
  end subroutine rd_ORCA_energy

!========================================================================================!

  subroutine rd_ORCA_wbo(mol,calc,iostatus)
!**************************************************
!* Read Wiberg/Mayer bond orders from the ORCA
!* output file written by ORCA_setup (i.e. ORCA.out).
!*
!* ORCA prints the Mayer bond orders as
!*     Mayer bond orders larger than 0.1
!*     B( 0-C , 1-S ) : 1.0550
!* i.e. with ZERO-based atom indices; those are
!* converted to the one-based CREST convention here.
!* The block is emitted whenever the Mayer population
!* analysis is active, which is the ORCA default and
!* can be forced with "! MAYER" in the input template.
!* If no bond orders are found a warning is issued but
!* the run is NOT aborted: WBOs only influence SHAKE
!* and the topology handling, never the energies.
!**************************************************
    use iso_fortran_env,only:wp => real64
    use strucrd
    use calc_type

    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    integer,intent(out) :: iostatus

    character(len=:),allocatable :: wbofile
    character(len=1024) :: atmp,inside,left,right
    integer :: i,j,p1,p2,pc,p3,ich,io,nwbo
    logical :: ex
    real(wp) :: dum

    call initsignal()

    iostatus = 0
    if (.not.calc%rdwbo) return

    !>--- locate the output file of the corresponding ORCA run
    if (allocated(calc%calcspace)) then
      wbofile = trim(calc%calcspace)//sep//'ORCA.out'
    else
      wbofile = 'ORCA.out'
    end if

    inquire (file=wbofile,exist=ex)
    if (.not.ex) then
      write (stderr,'(3a)') '**WARNING** ORCA output ',trim(wbofile),&
              & ' not found - cannot read bond orders'
      iostatus = 1
      return
    end if

    if (allocated(calc%wbo)) deallocate (calc%wbo)
    allocate (calc%wbo(mol%nat,mol%nat),source=0.0_wp)

    nwbo = 0
    open (newunit=ich,file=wbofile)
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit
      !>--- only consider lines of the form  B( i-X , j-Y ) : value
      if (index(atmp,'B(') == 0) cycle
      if (index(atmp,':') == 0) cycle
      p1 = index(atmp,'(')
      p2 = index(atmp,')')
      if (p1 == 0 .or. p2 <= p1) cycle
      inside = atmp(p1+1:p2-1)
      pc = index(inside,',')
      if (pc <= 1) cycle
      left = adjustl(inside(1:pc-1))
      right = adjustl(inside(pc+1:))
      i = extract_leading_int(left)
      j = extract_leading_int(right)
      if (i < 0 .or. j < 0) cycle
      if (i >= mol%nat .or. j >= mol%nat) cycle
      !>--- numerical value is printed after the colon
      p3 = index(atmp(p2:),':')
      if (p3 == 0) cycle
      read (atmp(p2+p3:),*,iostat=io) dum
      if (io /= 0) cycle
      !>--- ORCA counts atoms from zero, CREST from one
      calc%wbo(i+1,j+1) = dum
      calc%wbo(j+1,i+1) = dum
      nwbo = nwbo+1
    end do
    close (ich)

    if (nwbo == 0) then
      write (stderr,'(3a)') '**WARNING** no Mayer bond orders found in ',&
              & trim(wbofile),' - WBOs set to zero.'
      write (stderr,'(a)') '            Add "! MAYER" to the ORCA input template if WBOs are required.'
    end if

    return
  end subroutine rd_ORCA_wbo

!========================================================================================!

  integer function extract_leading_int(str)
!**************************************************
!* Return the first unsigned integer that appears in
!* a string, or -1 if there is none.  Used to parse
!* the "i-Element" tokens of the ORCA Mayer bond
!* order printout (e.g. "0-C" -> 0, "12-H" -> 12).
!**************************************************
    implicit none
    character(len=*),intent(in) :: str
    integer :: k,p,last

    extract_leading_int = -1

    !>--- position of the first digit
    p = 0
    do k = 1,len_trim(str)
      if (str(k:k) >= '0' .and. str(k:k) <= '9') then
        p = k
        exit
      end if
    end do
    if (p == 0) return

    !>--- end of the digit run
    last = p
    do k = p+1,len_trim(str)
      if (str(k:k) < '0' .or. str(k:k) > '9') exit
      last = k
    end do

    !>--- accumulate the value (no internal I/O / no read -> compiler-proof)
    extract_leading_int = 0
    do k = p,last
      extract_leading_int = extract_leading_int*10+(iachar(str(k:k))-iachar('0'))
    end do

  end function extract_leading_int

!========================================================================================!
!========================================================================================!
end module orca_sc

