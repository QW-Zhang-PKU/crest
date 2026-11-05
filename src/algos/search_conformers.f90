!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

subroutine crest_search_imtdgc(env,tim)
!*******************************************************************
!* This is the re-implementation of CREST's iMTD-GC default workflow
!* 
!* Compared to the legacy implementation, this version
!* is separated from the entropy algo to keep things clean
!* The entropy algo (sMTD-iMTD) can be found in search_entropy.f90
!*******************************************************************
  use crest_parameters, only: wp,stdout
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod
  use utilities
  use cregen_interface
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,m
  logical :: pr,wr
!===========================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  type(mddata),allocatable :: mddats(:)
  integer :: nsim,nallout

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,ex
  character(len=80) :: atmp,btmp,str
  logical :: multilevel(6)
  logical :: start,lower
!>--- LOCAL EXTENSION (2026-10-26): cumulative MTD-boost factor for the
!>    progressive kpush increase in the MAINLOOP (see apply below).
  real(wp),save :: mtdboost_fac = 1.0_wp
!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,13x,"│")') "CREST iMTD-GC SAMPLING"
  write (stdout,'(10x,"┕",49("━"),"┙")')

!>--- LOCAL CHANGE (2026-10-04): no blank line between the box and the input
!> structure, and none between the coordinates and whatever follows.  The
!> "rigid (SHAKE) bonds" report used to add one of its own as well
!> (dynamics_module.f90, crest_hard_bonds_to_shake), i.e. there were two
!> consecutive empty lines here.
!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)

!>--- saftey termination
  if(mol%nat .le. 2)then
     call catchdiatomic(env)
    return
  endif

!>--- sets the MD length according to a flexibility measure
  call md_length_setup(env) 
!>--- create the MD calculator saved to env
  call env_to_mddat(env)

  if (env%performMTD) then
!>--- (optional) calculate a short 1ps test MTD to check settings
   call tim%start(1,'Trial metadynamics (MTD)')
   call trialmd(env)  
   call tim%stop(1)
   if(env%iostatus_meta .ne. 0 ) return
  end if

!===========================================================!
!>--- Start mainloop 
  env%nreset = 0
  start = .true.
  MAINLOOP : do
    call printiter
    if (.not. start) then
!>--- clean Dir for new iterations, but leave iteration backup files
      call clean_V2i 
      env%nreset = env%nreset + 1
    else 
!>--- at the beginning, wipe directory clean
      call V2cleanup(.false.)
    end if
!===========================================================!
!>--- Meta-dynamics loop
  mtdloop: do i = 1,env%Maxrestart

    write(stdout,*)
    write(stdout,'(1x,a)') '------------------------------'
    write(stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
    write(stdout,'(1x,a)') '------------------------------'

    nsim = -1 !>--- enambles automatic MTD setup in init routines
    call crest_search_multimd_init(env,mol,mddat,nsim)
    allocate (mddats(nsim), source=mddat)
    call crest_search_multimd_init2(env,mddats,nsim)

!>--- LOCAL EXTENSION (2026-10-26): progressive boost of the MTD force constant.
!> A rigid system keeps falling back into the same basin, so a constant bias
!> never reaches a new one.  From the 2nd MTD iteration on, multiply the whole
!> metadfac array by CREST_MTD_BOOST (default 1.25), capped at a cumulative
!> factor of CREST_MTD_BOOST_MAX (default 4.0).  CREST_MTD_BOOST=1 disables it.
!> Note this also scales the -iru static bias, since parallel.f90 uses
!> maxval(env%metadfac) as its k - which is consistent by design.
    block
      character(len=32) :: bstr
      integer :: bstat
      real(wp) :: bfac, bcap
      bfac = 1.25_wp
      bstr = ''
      call get_environment_variable('CREST_MTD_BOOST',value=bstr,status=bstat)
      if (bstat == 0 .and. len_trim(bstr) > 0) then
        read (bstr,*,iostat=bstat) bfac
        if (bstat /= 0) bfac = 1.25_wp
      end if
      bcap = 4.0_wp
      bstr = ''
      call get_environment_variable('CREST_MTD_BOOST_MAX',value=bstr,status=bstat)
      if (bstat == 0 .and. len_trim(bstr) > 0) then
        read (bstr,*,iostat=bstat) bcap
        if (bstat /= 0) bcap = 4.0_wp
      end if
      if (i > 1 .and. bfac > 1.0_wp .and. allocated(env%metadfac)) then
        if (mtdboost_fac * bfac <= bcap + 1.0e-9_wp) then
          env%metadfac(:) = env%metadfac(:) * bfac
          mtdboost_fac = mtdboost_fac * bfac
          write (stdout,'(3x,a,i0,a,f6.3,a,f6.3,a)') &
          &  'MTD kpush boost (round ',i,'): x',bfac,' cumulative x',mtdboost_fac
        end if
      end if
    end block

    call tim%start(2,'Metadynamics (MTD)')
    call crest_search_multimd(env,mol,mddats,nsim)
    call tim%stop(2)
!>--- a file called crest_dynamics.xyz should have been written
    ensnam = 'crest_dynamics.xyz'
!>--- deallocate for next iteration
    if(allocated(mddats))deallocate(mddats)

!==========================================================!
!>--- Reoptimization of trajectories
    call tim%start(3,'Geometry optimization')
    call optlev_to_multilev(env%optlev,multilevel)
    call crest_multilevel_oloop(env,ensnam,multilevel)
    call tim%stop(3)
    if(env%iostatus_meta .ne. 0 ) return

!>--- LOCAL EXTENSION (2026-10-26): keep a STABLE copy of this iteration's full
!> rotamer ensemble before it is renamed away / deleted by clean_V2i.  The
!> -iru static bias in the NEXT iteration reads this file instead of
!> crest_conformers.xyz, which only contains the few structures that survived
!> the ewin/rthr filtering - on rigid systems that can be as little as one.
!>--- save the CRE under a backup name
    call checkname_xyz(crefile,atmp,str)
!>--- LOCAL EXTENSION (2026-10-26): keep a STABLE copy of this iteration's full
!> rotamer ensemble before it is renamed away / deleted by clean_V2i.  The -iru
!> static bias in the NEXT iteration reads this file instead of
!> crest_conformers.xyz, which only holds the few structures that survived the
!> ewin/rthr filtering.  NOTE: checkname_xyz() advances an internal counter, so
!> it must be called exactly ONCE here - that is why the copy sits right after
!> the single existing call.
    if (env%iru) then
      call copy(atmp,'crest_irubias.xyz')
    end if
    call checkname_xyz('.cre',str,btmp)
    call rename(atmp,btmp)
!>--- save cregen output
    call checkname_tmp('cregen',atmp,btmp)
    call rename('cregen.out.tmp',btmp)

!=========================================================!
!>--- cleanup after first iteration and prepare next
    if (i .eq. 1 .and. start) then
      start = .false.
!>-- obtain a first lowest energy as reference
      env%eprivious = env%elowest
!>-- remove the two extreme-value MTDs
      if (.not. env%readbias .and.  env%runver .ne. 33 .and. &
      &   env%runver .ne. 787878 ) then
        env%nmetadyn = env%nmetadyn - 2
      end if
!>-- the cleanup 
      call clean_V2i   
!>-- and always do two cycles of MTDs
      cycle mtdloop 
    endif
!=========================================================!
!>--- Check for lowest energy
    call elowcheck(lower,env)
    if (.not. lower) then
      exit mtdloop
    end if
  enddo mtdloop
!=========================================================!
!>--- collect all ensembles from mtdloop and merge
  write(stdout,*)
  write (stdout,'(''========================================'')')
  write (stdout,'(''           MTD Simulations done         '')')
  write (stdout,'(''========================================'')')
  write (stdout,'(1x,''Collecting ensmbles.'')')
!>-- collecting all ensembles saved as ".cre_*.xyz"
  call collectcre(env)                      
  call newcregen(env,0)
  call checkname_xyz(crefile,atmp,btmp)
!>--- remaining number of structures
  call remaining_in(atmp,env%ewin,nallout) 

!=========================================================!
!>--- (optional) Perform additional MDs on the lowest conformers
  if (env%rotamermds) then
    call tim%start(4,'Molecular dynamics (MD)')
    call crest_rotamermds(env,conformerfile)
    call tim%stop(4)
    if(env%iostatus_meta .ne. 0 ) return

!>--- Reoptimization of trajectories
    call checkname_xyz(crefile,atmp,btmp)
    write(stdout,'('' Appending file '',a,'' with new structures'')')trim(atmp)
    ensnam = 'crest_dynamics.xyz'
    call appendto(ensnam,trim(atmp))
    call tim%start(3,'Geometry optimization')
    call crest_multilevel_wrap(env,trim(atmp),-1)
    call tim%stop(3)
    if(env%iostatus_meta .ne. 0 ) return

    call elowcheck(lower,env)
    if (lower) then
      call checkname_xyz(crefile,atmp,str)
      call checkname_xyz('.cre',str,btmp)
      call rename(atmp,btmp)
      cycle MAINLOOP
    end if
  end if

!=========================================================!
!>--- (optional) Perform GC step
    if (env%performCross) then
      call tim%start(5,'Genetic crossing (GC)')
      call crest_newcross3(env)
      call tim%stop(5)
      if(env%iostatus_meta .ne. 0 ) return

      call confg_chk3(env)
      call elowcheck(lower,env)
      if (lower) then
        call checkname_xyz(crefile,atmp,str)
        call checkname_xyz('.cre',str,btmp)
        call rename(atmp,btmp)
        if (env%iterativeV2) cycle MAINLOOP
      end if
    end if

!==========================================================!
!>--- exit mainloop
   exit MAINLOOP
  enddo MAINLOOP

!==========================================================!
!>--- final ensemble optimization
    write (stdout,'(/)')
    write (stdout,'(3x,''================================================'')')
    write (stdout,'(3x,''|           Final Geometry Optimization        |'')')
    write (stdout,'(3x,''================================================'')')
    call tim%start(3,'Geometry optimization')
    call checkname_xyz(crefile,atmp,str)
    call crest_multilevel_wrap(env,trim(atmp),0) 
    call tim%stop(3)                 
    if(env%iostatus_meta .ne. 0 ) return

!==========================================================!
!>--- final ensemble sorting
!  call newcregen(env,0) 
!> this is actually done within the last crest_multilevel_
!> call, so I comment it out here

!==========================================================!
!>--- print CREGEN results and clean up Directory a bit
    write (stdout,'(/)')
    call smallhead('Final Ensemble Information')
    call V2terminating()

!==========================================================!
  return
end subroutine crest_search_imtdgc

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!

subroutine crest_multilevel_wrap(env,ensnam,level)
!*************************************************
!* wrapper for the multilevel_oloop to select
!* only a single optimization level
!*************************************************
  use crest_parameters, only: wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: ensnam
  integer,intent(in) :: level
  logical :: multilevel(6)
  integer :: k
  multilevel = .false.
  select case(level)
  case( 1: ) !> explicit selection (level is a positie integer)
    k = min(level,6) 
    multilevel(k) =.true.
  case default 
  !>-- map global variable to multilevel selection (level is 0 or negative)
    k = optlevmap_alt(env%optlev) + level
    k = max(1,k)
    multilevel(k) =.true.
  end select
  call crest_multilevel_oloop(env,ensnam,multilevel)
end subroutine crest_multilevel_wrap

!========================================================================================!
subroutine crest_multilevel_oloop(env,ensnam,multilevel_in)
!*******************************************************
!* multilevel optimization loop.
!* construct consecutive optimizations starting with
!* crude thresholds to very tight ones
!*******************************************************
  use crest_parameters, only: wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
!>--- LOCAL EXTENSION (2026-09-27): the "tight" stage of this chain uses its own
!> threshold level, loosened 10x - see the note on olev_tightcs in
!> optimize/optutils.f90.  optimize_module has a bare `private`, so the olev_*
!> constants are NOT re-exported through it and have to be imported here.
!> (2026-09-29: olev_vtightcs, the "half of tight" level used by stage 6, comes
!> from the same module and is imported the same way.)
  use optimize_utils,only:olev_tight,olev_vtight,olev_tightcs,olev_vtightcs
  use utilities
  use crest_restartlog
  use parallel_interface
  implicit none
  type(systemdata) :: env 
  character(len=*),intent(in) :: ensnam
  logical,intent(in) :: multilevel_in(6)
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,pr
  character(len=128) :: inpnam,outnam
  integer :: i,l,k,T,Tn
  real(wp) :: ewinbackup,rthrbackup
  real(wp) :: hlowbackup
  integer :: microbackup
  integer :: optlevelbackup
!>--- LOCAL EXTENSION (2026-09-29): stage 6 (Final Geometry Optimization) now
!> forces its own optimization engine, so it has to be backed up/restored
!> exactly like optlev/hlow_opt/micro_opt - otherwise the stage-6 setting
!> leaks into every later stage of this loop and into a standalone --opt run.
  integer :: enginebackup
  logical :: multilevel(6)

  interface 
    subroutine crest_refine(env,input,output)
      use crest_data
      implicit none
      type(systemdata),intent(inout) :: env
      character(len=*),intent(in) :: input
      character(len=*),intent(in),optional :: output
    end subroutine crest_refine
  end interface

!>--- save backup thresholds
  ewinbackup     = env%ewin
  rthrbackup     = env%rthr
  optlevelbackup = env%calc%optlev
  hlowbackup     = env%calc%hlow_opt
  microbackup    = env%calc%micro_opt
  enginebackup   = env%calc%opt_engine

!>--- set multilevels, or enforce just one
  multilevel(:) = .false.
  if(env%multilevelopt)then
     multilevel(:) = multilevel_in(:)
!>--- LOCAL CHANGE (dml): -onebatch should also get the LAST (vtightcs)
!> stage.  With the default optlev=2 optlev_to_multilev() clears stage 6,
!> so without this the user's vtightcs level would never be reached.
     if (env%onebatch) multilevel(6) = .true.    
  else
    k = optlevmap_alt(env%optlev)
    multilevel(k) = .true.
  endif

  pr = .false.
  l = count(multilevel)
  if( l > 1 )then
  pr = .true.
  write(stdout,*)
  write(stdout,'(1x,a)') '======================================'
  write(stdout,'(1x,a)') '|  Multilevel Ensemble Optimization  |'
  write(stdout,'(1x,a)') '======================================'
  endif
  
!>--- read ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write(stdout,*) '**ERROR** empty ensemble file ',trim(ensnam)
    env%iostatus_meta = status_failed
    return
  endif
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>--- track ensemble for restart
  call trackensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_oloop requires coordinates in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

  write(stdout,'(1x,a,i0,a,a,a)')'Optimizing all ',nall, &
  & ' structures from file "',trim(ensnam),'" ...'

!>--- sequential optimizations of ensembles
  dump = .true. !> optimized structures will be written to crest_ensemble.xyz
  do i=1,6
    if(multilevel(i))then
     !>--- set threads
       call new_ompautoset(env,'auto',nall,T,Tn)
     !>--- set optimization parameters
       call set_multilevel_options(env,i,.true.)
     !>--- run parallel optimizations
       call crest_oloop(env,nat,nall,at,xyz,eread,dump)
       deallocate(eread,at,xyz)
     !>--- rename ensemble and sort
       call checkname_xyz(crefile,inpnam,outnam)
       call rename(ensemblefile,trim(inpnam))
     !>--- check for empty ensemble content
     !> LOCAL CHANGE (2026-09-27): only the SIZE is wanted here, so use the
     !> cheap counter (no coordinate parsing) instead of the full
     !> rdensembleparam() parse.  Identical nat (= largest structure) / nall.
       call count_xyz_structures(trim(inpnam),nat,nall)
       if (nall .lt. 1) then
         write(stdout,*) '**ERROR** empty ensemble file',trim(inpnam)
         env%iostatus_meta = status_failed
         return
       endif

       write(stdout,*)
     !==========================================================!
     !>-- dedicated ensemble refinement step (overwrites inpnam)
      call  crest_refine(env,trim(inpnam))
     !==========================================================!

     !>--- CREGEN sorting
       call sort_and_check(env,trim(inpnam))
       call checkname_xyz(crefile,inpnam,outnam)
     !>--- check for empty ensemble content (again)
     !> LOCAL CHANGE (2026-09-27): size-only -> cheap counter (see above).
       call count_xyz_structures(trim(inpnam),nat,nall)
       if (nall .lt. 1) then
         write(stdout,*) '**ERROR** empty ensemble file',trim(inpnam)
         env%iostatus_meta = status_failed
         return
       endif
     !>--- read new ensemble for next iteration
       allocate (xyz(3,nat,nall),at(nat),eread(nall))
       call rdensemble(trim(inpnam),nat,nall,at,xyz,eread)
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
     !>--- Important: crest_oloop requires coordinates in Bohrs
       xyz = xyz / bohr
     !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
     !>--- restore default sorting thresholds
       env%ewin        = ewinbackup
       env%rthr        = rthrbackup
       env%calc%optlev = optlevelbackup
       env%calc%hlow_opt  = hlowbackup
       env%calc%micro_opt = microbackup
       env%calc%opt_engine = enginebackup
    endif
  enddo

  if(allocated(eread)) deallocate(eread)
  if(allocated(at))  deallocate(at)
  if(allocated(xyz)) deallocate(xyz)
  return
contains
  subroutine set_multilevel_options(env,i,pr)
    implicit none
    type(systemdata) :: env
    integer,intent(in) :: i
    logical,intent(in) :: pr 

    env%calc%hlow_opt  = env%hlowopt
    env%calc%micro_opt = nint(env%microopt)

    select case( i )  
    case( 1 )
     if(pr) call smallhead('crude pre-optimization')
     env%calc%optlev =  -3
     !> larger thresholds
     env%rthr = env%rthr * 2.0d0
     env%ewin = aint(env%ewin * 2.0d0)
    case( 2 )
     if(pr) call smallhead('optimization with very loose thresholds')
     env%calc%optlev =  -2
     env%rthr = env%rthr *1.5d0
     env%ewin = aint(env%ewin * 2.0d0)
    case( 3 )
     if(pr) call smallhead('optimization with loose thresholds')
     env%calc%optlev =  -1
      env%ewin = aint(env%ewin*(10.0d0/6.0d0))
    case( 4 )
     if(pr) call smallhead('optimization with regular thresholds')
     env%calc%optlev =  0
    case( 5 )
     if(pr) call smallhead('optimization with tight thresholds')
!>--- LOCAL EXTENSION (2026-09-27): was olev_tight (1e-6 Eh / 8e-4 Eh/a0).
!> This stage only polishes the ensemble for deduplication; the real refinement
!> happens at a better level afterwards, so 10x looser (1e-5 / 8e-3) is plenty
!> and saves optimisation cycles on every structure.  The criteria actually used
!> are printed by print_opt_data() below, so the log stays truthful.
!>--- LOCAL CHANGE (dml): the loosened 'tightcs' level is now used ONLY by
!> -onebatch.  Every other run type (including --quick) gets the NATIVE
!> olev_tight, so the conformer-search ladder keeps CREST's original meaning
!> unless the user explicitly asks for the onebatch workflow.
     if (env%onebatch) then
       env%calc%optlev = olev_tightcs
     else
       env%calc%optlev = olev_tight
     end if
    case( 6 )
     if(pr) call smallhead('optimization with very tight thresholds')
!>--- LOCAL EXTENSION (2026-09-29): stage 6 is now DEFINED as "half of stage 5",
!> i.e. half of olev_tightcs, instead of the absolute olev_vtight level.
!>
!> WHY: stage 5 was retuned to a separate level (olev_tightcs, 1e-4 Eh / 4e-3
!> Eh/a0) that is deliberately much looser than olev_tight.  Stage 6 used to
!> point at the ABSOLUTE level 2 (olev_vtight = 1e-7 Eh / 2e-4 Eh/a0), which
!> after that retune is 1000x tighter in energy and 20x tighter in gradient
!> than the stage just below it.  A 3-order-of-magnitude jump between two
!> consecutive stages is not a meaningful ladder: stage 6 would then dominate
!> the entire runtime and, as the log in this run shows, fail on most
!> structures ("1 of 5 structures successfully optimized").
!> Deriving stage 6 from olev_tightcs keeps the ladder monotonic and
!> proportional (x2 tighter) no matter how stage 5 is tuned later.
!>
!> The two criteria are halved INDEPENDENTLY from the stage-5 pair:
!>   ethr: 1e-4 -> 5e-5 Eh        (= 0.0314 kcal/mol)
!>   gthr: 4e-3 -> 2e-3 Eh/a0
!> Both the stage-5 pair and this pair are defined in ONE place
!> (optimize/optutils.f90, get_optthr), so a future retune of olev_tightcs
!> cannot desynchronise them.
!> olev_vtight (level 2) itself is left untouched - it remains the level of the
!> user-facing "vtight" optlevel keyword.
!>--- LOCAL CHANGE (dml): same rule as stage 5 - 'vtightcs' only for -onebatch.
     if (env%onebatch) then
       env%calc%optlev = olev_vtightcs
     else
       env%calc%optlev = olev_vtight
     end if
!>--- LOCAL EXTENSION (2026-09-29) end
!>
!>--- LOCAL EXTENSION (2026-09-29): the FINAL Geometry Optimization runs with
!> the gradient-descent engine by default, independent of `calc%opt_engine`.
!>
!> WHY: stage 6 is the last thing that touches the ensemble before CREGEN
!> deduplicates it and the run terminates.  Its only jobs are (a) to shave the
!> residual gradient off structures that stage 5 left slightly above threshold
!> and (b) to never make a geometry worse.  ANCOPT is a quasi-Newton method
!> built on a MODEL Hessian - it is fast when the model is good, but it can and
!> does throw a structure into a nearby saddle or a different basin when the
!> model Hessian is badly conditioned, and at a very tight gradient threshold
!> that shows up as a sporadic hard failure.  Gradient descent has no model
!> Hessian at all: every accepted step is a strictly downhill, capped,
!> Barzilai-Borwein-scaled displacement, so it cannot overshoot on a bad
!> curvature estimate.  It needs more cycles for the same gradient reduction,
!> but it fails far less often, and at this point in the pipeline a slightly
!> larger cycle count is cheaper than losing a structure.
!>
!> The engine is set on the local `env%calc` copy and restored right after
!> every stage (see the save/restore block at the top of this subroutine),
!> so stages 1-5 keep whatever engine the user selected, and a `[optimize]
!> opt = "..."` from the TOML still governs the MTD-level and GC-level
!> optimizations.  `-nogd` does not exist: this is a deliberate default, not a
!> user flag.  Comment out the line below to go back to the global engine.
     env%calc%opt_engine = -1 !> gradient descent
    case default
     if(pr) call smallhead('optimization with default thresholds')
!>--- LOCAL EXTENSION (2026-09-27): these two were hardcoded 6.0_wp / 0.125_wp
!> and would silently undo the tuned defaults from confparse.f90 (8.0 / 0.15)
!> on any multilevel run that reaches this branch.  Keep the single source of
!> truth in confparse.f90 by leaving env%ewin / env%rthr untouched here.
     env%calc%optlev = 0
    end select

    call print_opt_data(env%calc, stdout)

  end subroutine set_multilevel_options
end subroutine crest_multilevel_oloop

!========================================================================================!
subroutine crest_rotamermds(env,ensnam)
!***********************************************************
!* set up and perform several MDs at different temperatures
!* on the lowest few conformers
!***********************************************************
  use crest_parameters, only: wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: ensnam

  integer :: nsim
  type(mddata) :: mddat
  type(mddata),allocatable :: mddats(:)
  type(coord) :: mol 
  type(coord),allocatable :: mols(:)
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  integer :: nstrucs,i,j,k,io
  real(wp) :: temp,newtemp
  character(len=80) :: atmp
  
!>--- coord setup
  call env%ref%to(mol)
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) then
    write(stdout,*) '**ERROR** empty ensemble file',trim(ensnam)
    env%iostatus_meta = status_failed
    return
  endif

!>--- determine how many MDs need to be run and setup
  call adjustnormmd(env)
  nstrucs = min(nall, env%nrotammds)
  nsim = nstrucs * env%temps 
  call crest_search_multimd_init(env,mol,mddat,nsim)
  allocate (mddats(nsim), source=mddat)
  call crest_search_multimd_init2(env,mddats,nsim)
!>--- adjust T's and runtimes
  k = 0
  do i=1,env%temps
    !> each T block 100K higher
    temp = env%nmdtemp + (i-1)*100.0_wp
    do j=1,nstrucs
      k= k + 1
      mddats(k)%tsoll = temp
      !> reduce runtime by 50% compared to MTDs
      mddats(k)%length_ps = mddats(k)%length_ps * 0.5_wp
      call mdautoset(mddats(k),io)
    enddo
  enddo 

!>--- read ensemble and prepare mols
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: mols must be in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  allocate(mols(nsim), source=mol) 
  k = 0
  do i=1,env%temps
    do j=1,nstrucs
      k = k + 1
      mols(k)%at = at
      mols(k)%xyz(:,:) = xyz(:,:,j)
    enddo
  enddo
  deallocate(eread,at,xyz)
  
!>--- print what we are doing
  write(stdout,*)
  write(atmp,'(''Additional regular MDs on lowest '',i0,'' conformer(s)'')')nstrucs
  call smallheadline(trim(atmp)) 

!>--- and finally, run the MDs
  call crest_search_multimd2(env,mols,mddats,nsim)

  if(allocated(mols))deallocate(mols)
  if(allocated(mddats))deallocate(mddats) 
  return
end subroutine crest_rotamermds

!========================================================================================!
subroutine crest_newcross3(env)
!*********************************************************
!* wrapper for the conformational crossing
!* takes the latest crest_rotamers_*, crosses structures
!* and writes the optimized ones back to the file
!*********************************************************
  use crest_parameters
  use crest_data
  use iomod
  use utilities
  implicit none
  type(systemdata) :: env  
  real(wp) :: ewinbackup
  integer  :: i,imax,tmpconf,nremain
  character(len=128) :: inpnam,outnam,refnam
  character(len=512) :: thispath,tmppath
  logical :: multilevel(6)
  real(wp),allocatable :: backupthr(:)

  multilevel = .false.
  call getcwd(thispath)

  do i = 1,1  !>-- technically it would be possible to repeat the crossing
!>-- determine max number of new structures
    imax = min(nint(env%mdtime*50.0d0),5000)
    if (env%setgcmax) then
      imax = nint(env%gcmax)
    else if(imax<0)then
      imax=5000
    end if
    if (env%quick) then
      imax = nint(float(imax)*0.5d0)
    end if

!>-- call the crossing routine
    call checkname_xyz(crefile,refnam,tmppath)
    call touch(trim(tmppath)) 
    call crest_crossing(env,imax,trim(refnam),env%gcmaxparent)
    if (imax .lt. 1) then
      call remove(trim(tmppath))
      return
      exit
    end if

!>-- optimize ensemble
    if (env%gcmultiopt) then !>-- optionally split into two steps
      call optlev_to_multilev(env%optlev,multilevel)
    else
      multilevel(4) = .true.
    end if
!>--- LOCAL EXTENSION (2026-09-29): always send the crossed structures through
!> stage 5 (tight, = olev_tightcs) as well, regardless of how the earlier stages
!> were selected.
!>
!> WHY: GC cuts and splices two conformers at a randomly chosen bond.  The
!> product is a NEW chemical geometry that has never been relaxed at any level -
!> a bond just changed length and two half-molecules were re-joined, so the
!> gradient can be one or two orders of magnitude larger than that of an MTD
!> snapshot.  With the default optlev = 2 the chain's only GC stage used to be
!> multilevel(4) = "regular" (or, with -gcmopt, stages 1 + 5).  A "regular"
!> (optlev 0) relaxation is not enough to bring such a geometry down to where
!> CREGEN can compare it reliably against the MTD-derived conformers, and the
!> subsequent single-sided comparison is exactly where spurious new "conformers"
!> come from.  Running stage 5 unconditionally costs one extra ANCOPT pass over
!> the GC output only (typically a few hundred structures, not the whole
!> ensemble), which is negligible next to the MTD step that produced them.
!>
!> NOTE: stage 5 is idempotent in the sense that it is a normal, monotonically
!> tighter stage - adding it can only improve the geometry handed to the final
!> wrap(...,0), it can never loosen anything.  It also keeps the -nogcmopt path
!> (multilevel(4) alone) meaningful: that flag now selects whether the CRUDE
!> pre-stage runs, not whether any GC optimisation happens at all.
    multilevel(5) = .true.
    call crest_multilevel_oloop(env,'confcross.xyz',multilevel)
    if(env%iostatus_meta .ne. 0 ) return

!>-- append optimized crossed structures and original to a single file
    call checkname_xyz(crefile,inpnam,outnam)
    write(stdout,'(a,a)')'appending new structures to ',trim(refnam)
    call appendto(trim(inpnam),trim(refnam))
    do while(trim(inpnam).ne.trim(refnam))
      call remove(trim(inpnam))
      call checkname_xyz(crefile,inpnam,outnam)
    enddo
  end do
end subroutine crest_newcross3
