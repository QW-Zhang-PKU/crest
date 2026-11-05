!=========================================================================================!
!=========================================================================================!
!>--- LOCAL EXTENSION (2026-09-29): collector for CREGEN-excluded structures.
!>
!> WHY: CREGEN drops structures at four independent gates (broken/clashing,
!> topology mismatch, EWIN cut-off, and the RMSD/rotamer de-duplication).  All
!> four used to be reported as a bare count ("number of removed clashes : 17"),
!> so there was no way to answer the two questions that actually matter when an
!> ensemble comes out smaller than expected:
!>   (a) WHICH structures went away, and
!>   (b) WHY - which gate rejected them.
!> With stage 5/6 running tight thresholds this matters more, because a
!> half-converged structure can be dropped by the RMSD gate and silently
!> disappear from the final ensemble.
!>
!> WHAT: every rejected structure is appended to crest_excluded_structures.xyz
!> together with a machine-readable tag and a human-readable reason encoded in
!> the xyz comment line:
!>
!>   <energy> <TAG> <reason text>
!>
!>   TAG = CLASH   rejected by discardbroken()   (broken geometry / dissociation)
!>         TOPO    rejected by cregen_topocheck() (neighbour list mismatch, E/Z)
!>         EWIN    rejected by cregen_esort()     (outside the energy window)
!>         DUP     rejected by cregen_CRE()       (RMSD/rotamer duplicate)
!>
!> The file is a normal, valid multi-structure xyz file, so it can be fed
!> straight back into crest or opened in any viewer.  The comment line is the
!> only addition; the energy is kept FIRST and in the same format as the input
!> ensemble so that grepenergy() still parses it.
!>
!> DESIGN NOTES:
!>  * Truncated once per run (exclfile_init) and appended to afterwards, so the
!>    file is the complete audit trail of the run rather than just the last
!>    sorting step.
!>  * Never fails the run: every file operation is best-effort.  A missing or
!>    unwritable file must not turn a working conformer search into an error.
!>  * The switch is the environment variable CREST_EXCL_XYZ (default OFF).
!>    Writing the file costs one extra copy of every rejected structure, which
!>    can be a large fraction of the ensemble, so it is opt-in.
!=========================================================================================!

module exclfile_module
  use iso_fortran_env,only:wp => real64,output_unit
  implicit none
  private

  public :: exclfile_check,exclfile_init,exclfile_write,exclfile_free

  character(len=*),parameter :: exclfile_name = 'crest_excluded_structures.xyz'

  logical :: active = .false.
  logical :: initialised = .false.
  integer :: uch = -1     !> unit of the open file, -1 = closed

contains

!>---
!>--- read the environment switch.  Called once per newcregen().
!>--- CREST_EXCL_XYZ=0/none/false forces it off explicitly.
!>---
  subroutine exclfile_check()
    implicit none
    character(len=64) :: val
    integer :: estat

    active = .false.
    call get_environment_variable('CREST_EXCL_XYZ',value=val,status=estat)
    if (estat /= 0) return
    if (trim(val) == '0' .or. trim(val) == 'none' .or. trim(val) == 'false') return
    active = .true.
    return
  end subroutine exclfile_check

!>---
!>--- truncate / open the output file.  Safe to call repeatedly: only the first
!>--- call of a run does anything, so successive newcregen() calls append.
!>---
  subroutine exclfile_init()
    implicit none
    integer :: ios

    if (.not. active) return
    if (initialised) return

    !> start from a clean file so that a re-run does not mix two runs' rejects
    open (newunit=uch,file=exclfile_name,status='replace',action='write',iostat=ios)
    if (ios /= 0) then
      uch = -1
      active = .false.
      return
    end if
    initialised = .true.
    return
  end subroutine exclfile_init

!>---
!>--- append ONE rejected structure.
!>---   nat,at,xyz  - the structure (xyz in whatever unit the caller holds)
!>---   comment     - the original xyz comment line (carries the energy)
!>---   tag,reason  - the gate and why it fired
!>---   factor      - unit conversion to Angstrom (bohr for the internal arrays,
!>---                 1.0 when the caller already holds Angstrom)
!>---
  subroutine exclfile_write(nat,at,xyz,comment,tag,reason,factor)
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    character(len=*),intent(in) :: comment
    character(len=*),intent(in) :: tag,reason
    real(wp),intent(in) :: factor
    integer :: i,ios
    character(len=512) :: line

    if (.not. active) return
    !>--- LOCAL FIX (2026-10-26): open(newunit=uch,...) returns a NEGATIVE
    !> unit by Fortran 2008 design, so `if (uch < 0) return` treated every
    !> successfully opened file as "not open" - this audit never wrote anything.
    if (.not. initialised) return

    write (line,'(a,1x,a,1x,a)') trim(comment),trim(tag),trim(reason)
    write (uch,'(i0)',iostat=ios) nat
    if (ios /= 0) return
    write (uch,'(a)') trim(line)
    do i = 1,nat
      write (uch,'(a2,3(1x,f16.8))',iostat=ios) sym_of(at(i)), &
      &  xyz(1,i)*factor,xyz(2,i)*factor,xyz(3,i)*factor
      if (ios /= 0) return
    end do
    flush (uch)
    return
  end subroutine exclfile_write

!>---
!>--- close the file at the very end of the run.
!>---
  subroutine exclfile_free()
    implicit none
    if (initialised) close (uch)
    uch = -1
    initialised = .false.
    return
  end subroutine exclfile_free

!>---
!>--- minimal element symbol lookup, so this module does not have to depend on
!>--- miscdata / PSE and its (large) transitive import list.  Covers H..Kr,
!>--- which is ample for the light-element conformer searches CREST targets;
!>--- anything else is written as "Xx" rather than crashing.
!>---
  function sym_of(z) result(sym)
    implicit none
    integer,intent(in) :: z
    character(len=2) :: sym
    character(len=2),parameter :: pse(36) = &
    &  ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
    &   'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
    &   'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
    &   'Ga','Ge','As','Se','Br','Kr']

    if (z >= 1 .and. z <= 36) then
      sym = pse(z)
    else
      sym = 'Xx'
    end if
    return
  end function sym_of

end module exclfile_module
