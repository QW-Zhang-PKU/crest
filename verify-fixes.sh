#!/usr/bin/env bash
# =============================================================================
#  verify-fixes.sh — confirm that this CREST tree carries every local fix
#
#  Two independent checks:
#    A) reverse-dry-run of each patch in ./patches/  -> "the patch is applied"
#    B) grep for the one distinctive line each fix introduced
#
#  Exit status 0 = all good.  Run it before building if you are unsure whether
#  a source tree already contains the fixes.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

fail=0

# ---------------------------------------------------------------------------
# A) patch bookkeeping
# ---------------------------------------------------------------------------
echo "=== A. patch reverse-check (dry run) ==================================="
if command -v patch >/dev/null 2>&1; then
  if [ -d patches ] && compgen -G "patches/*.patch" >/dev/null; then
    for p in patches/*.patch; do
      # -R reverses the patch; if the tree already contains the fixed lines the
      # reverse dry-run succeeds.  --fuzz=0 keeps it strict.
      if patch -p1 -R --dry-run -f --fuzz=0 <"$p" >/dev/null 2>&1; then
        printf '  APPLIED      %s\n' "$(basename "$p")"
      else
        printf '  NOT-APPLIED  %s\n' "$(basename "$p")"
        fail=$((fail + 1))
      fi
    done
  else
    echo "  (no patches/ directory - skipping)"
  fi
else
  echo "  (patch not installed - skipping)"
fi

# ---------------------------------------------------------------------------
# B) feature spot-checks
# ---------------------------------------------------------------------------
echo
echo "=== B. feature spot-checks ============================================"
# file|pattern|description
CHECKS=(
"src/iomod.F90|c_setenv(c_name,c_VALUE,c_overwrite)|setenv() declared with the 3rd 'overwrite' argument"
"src/iomod.F90|c_setenv(env//c_null_char,trim(str)//c_null_char,1_c_int)|setenv_int() passes overwrite=1"
"src/iomod.F90|1_c_int) !> 1 = overwrite existing value|all three setenv_* call sites updated"
"src/ompmklset.F90|OMP_MAX_ACTIVE_LEVELS|children get OMP_MAX_ACTIVE_LEVELS=1"
"src/calculator/xtb_sc.f90|if (.not.allocated(calc%gradfile)) calc%gradfile = gf|xtb gradient path re-derived per worker"
"src/calculator/orca_sc.f90|calc%gradfile = gf|orca gradient path re-derived per worker"
"src/calculator/orca_sc.f90|extract_leading_int|ORCA Mayer bond-order reader implemented (was an error stub)"
"src/calculator/orca_sc.f90|2> ORCA.err|ORCA stderr captured to ORCA.err (upstream discards it with /dev/null)"
"src/calculator/subprocess_types.f90|character(len=orca_linelen),allocatable :: input(:)|ORCA template is a fixed-length array (gcc PR89200 deferred-length trap)"
"src/calculator/subprocess_types.f90|source=repeat(\" \",orca_linelen)|ORCA allocate SOURCE length == element length (gcc PR113793 heap-corruption trap)"
"src/calculator/subprocess_types.f90|orca_nprocs() result(np)|ORCA core budget read back from OMP_NUM_THREADS (MPI rank count)"
"src/calculator/subprocess_types.f90|np > 1 .and. .not.self%mpi|ORCA.in gets a %pal nprocs N block (ORCA ignores OMP_NUM_THREADS)"
"src/calculator/orca_sc.f90|orca_omp_threads_string|ORCA child thread budget tunable via CREST_ORCA_OMP_THREADS (default 1)"
"src/calculator/orca_sc.f90|CREST_ORCA_NO_OVERSUBSCRIBE|--oversubscribe can be dropped (it disables MPI process binding)"
"src/calculator/orca_sc.f90|CREST_ORCA_MPI_ARGS|extra args (e.g. --bind-to core) can be passed through to ORCA/MPI"
"src/ompmklset.F90|case ('fixed')|new_ompautoset 'fixed': CREST_CORES_PER_JOB cores per job, T/n jobs at once"
"src/algos/parallel.f90|new_ompautoset(env,'fixed',nall,T,Tn)|crest_oloop + crest_sploop reserve a fixed core count per batch job"
"src/ompmklset.F90|case ('share')|new_ompautoset 'share': one job gets exactly the budget a parallel job has"
"src/algos/setuptest.f90|new_ompautoset(env,'share',env%nmetadyn,T,Tn)|trial MTD runs with T/nmetadyn cores (was 'min' = 1 core)"
"src/algos/setuptest.f90|timeestimate(rtime,env%mdtime,env%nmetadyn,env%threads,Tn,trial_md_ps)|MTD runtime estimate uses the measured per-job core budget"
"src/algos/setuptest.f90|trial_md_ps = 0.25_wp|trial MTD shortened to 1/4 (0.25 ps instead of 1 ps)"
"src/algos/setuptest.f90|min(trial_vbdump_fs,trial_md_ps*500.0_wp)|trial Vbias dump is a flat 50 fs, capped at half the trial so >=2 updates always fit (2026-09-23: the old 550*trial_md_ps scaling gave exactly 1 update)"
"src/algos/setuptest.f90|rtime*(mdtime/max(trialps,1.0e-6_wp))|estimate divides out the trial length instead of assuming 1 ps"
"src/confparse.f90|env%mddump = 500|PRODUCTION MTD Vbias dump default halved to 500 fs / 0.5 ps (upstream 1000 fs)"
"src/calculator/turbom_sc.f90|calc%gradfile = gf|turbom gradient path re-derived per worker"
"src/calculator/turbom_sc.f90|if (.not.allocated(calc%efile)) calc%efile = ef|turbom energy path re-derived per worker"
"src/parsing/parse_calcdata.f90|SAFETY GUARD|TOML-before-coordinates SIGSEGV guard"
"src/iomod.F90|subroutine crest_scratch_dir|CREST_SCRATCH_DIR helper: returns + mkdir -p the node-local scratch root"
"src/calculator/calc_type.f90|call crest_scratch_dir(scratch)|calculation.level.N can be relocated off an NFS workspace"
"src/cleanup.f90|calculation.level.*|relocated scratch dirs are cleaned up (never the scratch root itself)"
"src/algos/parallel.f90|crest_calc_is_subprocess|subprocess calculators detected -> no nested OMP team"
"src/dynamics/metadynamics_module.f90|mtd_omp_min_ncur|MTD RMSD OMP region only for long reference lists"
"src/sigterm.f90|signalset|signal handlers installed only once"
"src/sigterm.f90|crest_install_signals|handlers are the async-signal-safe C ones (no more write(*,*)/cs_shutdown in a handler)"
"src/signal.c|crest_runcmd|SIGINT-safe system(): fork+execl+waitpid instead of POSIX system(3)"
"src/signal.c|crest_sighandler|signal handler only uses write(2)/kill()/_exit()"
"src/signal.c|setpgid (0, 0)|external programs run in their own process group (killable as one job)"
"src/signal.c|crest_install_signals|SIGHUP/SIGINT/SIGQUIT/SIGTERM/SIGXCPU installed via sigaction"
"src/iomod.F90|crest_runcmd(trim(cmd)//c_null_char)|command() no longer uses system()/EXECUTE_COMMAND_LINE"
"src/algos/setuptest.f90|subroutine crest_time2string|shared seconds->string helper (renamed: time2string clashes with trialmd_legacy)"
"src/algos/setuptest.f90|call optprof%start(1)|initial geometry optimization is timed"
"src/algos/setuptest.f90|noptcycles = count_xyz_frames('crestopt.xyz',molopt%nat)|initial optimization cycle count read back from crestopt.xyz"
"src/algos/setuptest.f90|Initial geometry optimization: |initial optimization reports runtime + number of cycles"
"src/algos/setuptest.f90|MDSTART%progress = .true.|trial MTD progress marks enabled"
"src/algos/setuptest.f90|MDSTART%progress_tag = 'Trial MTD'|trial MTD progress lines are labelled"
"src/algos/setuptest.f90|'|   MD simulation time   :'|trial MTD prints the settings it actually uses, in the same box layout as the production banner"
"src/algos/setuptest.f90|write (stdout,'(2x,a)') repeat(':',45)|the trial MTD box is closed by a 45-colon rule"
"src/algos/parallel.f90|integer,intent(in),optional :: ntot|parallel_md_block_printout is told the batch size"
"src/algos/parallel.f90|if (closebox) write (stdout,'(2x,a)') repeat(':',45)|only the LAST MTD/MD box of a batch is closed by a 45-colon rule"
"src/choose_settings.f90|subroutine crest_hard_bonds_setup|the --cinp -> rigid SHAKE setup+report is a shared routine"
"src/choose_settings.f90|if (env%mddat%shk%nusr == 0)|it is guarded by nusr, so it is set up and printed exactly once"
"src/choose_settings.f90|call crest_hard_bonds_setup(env)|md_length_setup() reports the SHAKE bonds ABOVE the flexibility banner (env_to_mddat covers the --md path)"
"src/dynamics/dynamics_module.f90|dat%progress_steps = max(1,dat%length_steps/4)|progress interval derived as length_steps/4 -> marks at 25/50/75 %"
"src/dynamics/dynamics_module.f90|call md_progress(stdout,dat,t,pcount0,prate)|MD loop emits coarse progress marks (independent of pr)"
"src/dynamics/dynamics_module.f90|integer :: progress_steps = 0|mddata carries the progress-mark settings"
"src/algos/parallel.f90|mddats(i)%progress = .true.|PRODUCTION MD/MTD also emits progress marks (was only the trial MTD)"
"src/algos/parallel.f90|mddats(i)%progress_tag = atmp|production progress lines are tagged with the job index (*MTD n)"
"src/algos/parallel.f90|call system_clock(count=tb0,count_rate=brate)|phase clock starts before the (sequential) worker setup"
"src/algos/parallel.f90|call system_clock(count=tbmid)|marker splits worker setup from parallel MD sampling"
"src/algos/parallel.f90|TOTAL of this MTD round|true wall time of the whole MTD round is printed (setup + sampling)"
# --- trajectories written as .xyz instead of .trj / optimization log as .xyz
"src/algos/dynamics.f90|character(len=*),parameter :: trjf='crest_dynamics.xyz'|MTD ensemble is crest_dynamics.xyz (was .trj)"
"src/algos/parallel.f90|'crest_',i,'.xyz'|per-worker trajectory inside MDFILES is .xyz (was .trj)"
"src/dynamics/dynamics_module.f90|'crest_',dat%md_index,'.xyz'|dynamics trajectory file is .xyz (was .trj)"
"src/algos/search_entropy.f90|'crest_',i,'.xyz'|entropy MD trajectories are .xyz (was .trj)"
"src/optimize/ancopt.f90|open (newunit=ilog,file='crestopt.xyz')|ANC optimizer writes crestopt.xyz (was crestopt.log)"
"src/optimize/rfo.f90|open (newunit=ilog,file='crestopt.xyz')|RFO optimizer writes crestopt.xyz (was crestopt.log)"
"src/optimize/gd.f90|open (newunit=ilog,file='crestopt.xyz')|GD optimizer writes crestopt.xyz (was crestopt.log)"
# --- intermediate-file staging: hot -> CREST_SCRATCH_DIR, cold -> CREST_WORK_DIR
"src/iomod.F90|subroutine crest_work_dir|CREST_WORK_DIR helper: cold intermediate files get their own subdirectory"
"src/iomod.F90|subroutine crest_stage_files|hot files -> CREST_SCRATCH_DIR, cold files -> CREST_WORK_DIR"
"src/iomod.F90|subroutine crest_hot_dir|helper: hot subdirectory below CREST_SCRATCH_DIR"
"src/iomod.F90|subroutine crest_link|staging is done with symlinks (benign failure: file stays in the workspace)"
"src/iomod.F90|character(len=32),parameter :: hotl(2) = [character(len=32) ::|hot list trimmed to the two names that can actually be written (MDFILES, OPTIM)"
"src/iomod.F90|Removing the three dead names is behaviour-neutral|the removal of METADYN/NORMMD/MRMSD from the hot list is documented where it lives"
"src/iomod.F90|'coord','wbo','crestopt.xyz','gfnff_adjacency','gfnff_topo',|cold list: written a few times per run -> crest_work/"
"src/iomod.F90|export CREST_STAGE_FILES=0|escape hatch: CREST_STAGE_FILES=0 restores the upstream layout"
"src/crest_main.f90|call crest_stage_files()|file staging runs once at start-up"
"src/algos/parallel.f90|call crest_hot_dir(mdirtmp,'MDFILES')|MDFILES relocated by PATH (not symlink: must survive rmrf+makedir)"
"src/algos/setuptest.f90|call crest_hot_dir(dirnam,'TRIALMD')|TRIALMD relocated by PATH into node-local scratch"
# --- log / output format improvements (2026-09-15) --------------------------
"src/printouts.f90|zqw  (locally modified source)|banner names the local maintainer instead of 'usr@localhost.localdomain'"
"src/calculator/subprocess_types.f90|LOCAL EXTENSION: echo the template verbatim|ORCA input template is echoed into the log after 'Reading ORCA input template'"
"src/algos/setuptest.f90|(2x,a,i8,a)|Trial MTD integer settings use a fixed i8 field inside the box (was i0 -> 'steps              :50')"
"src/dynamics/dynamics_module.f90|\", remaining \"|progress line keeps an explicit 'remaining' label (2026-09-24 compact form: was 'expected remaining time', before that the cryptic 'ETA')"
# --- the progress line was 108 columns wide and wrapped on a narrow terminal.
# --- Four jobs interleave their marks inside !$omp critical, so a wrapped line
# --- is not just ugly, it is unreadable.  108 -> 73 columns.
"src/dynamics/dynamics_module.f90|\"% (\",f6.2,\" /\",f6.2,\" ps, elapsed \"|progress line: single spaces and f6.2 ps fields (was f7.2 + double spaces)"
"src/dynamics/dynamics_module.f90|'(i0,\"h\",i2.2,\"m\",i2.2,\"s\")'|compact h/m/s time string: no spaces, zero-padded -> fixed width"
"src/dynamics/dynamics_module.f90|if (nint(seconds) >= 60) then|anint() rounding up to 60 s is normalised before formatting"
"src/dynamics/dynamics_module.f90|t == 3*dat%progress_steps|exactly THREE progress marks (25/50/75 %); no 'mod()'-derived 99/100 % line"
"src/timer.f90|lshowcpu|write_timing can suppress the cpu-time / ratio c/w lines"
# --- the compact timing line used to print "wallmins" AFTER the hours had been
#     subtracted: 2 h 14 min 8.978 s came out as "14 min, 8.978 sec"
"src/timer.f90|walldays > 0.or.wallhours > 0|write_timing no longer reports runtimes modulo one hour"
"src/algos/parallel.f90|.true.,.not.subproc|optimization/singlepoint summary drops cpu-time for external calculators"
"src/algos/singlepoint.f90|.true.,.not.subproc|standalone singlepoint summary drops cpu-time for external calculators"
# --- ORCA command line must pass extra MPI switches as ONE quoted argv entry
"src/calculator/orca_sc.f90|mpiargs = '--oversubscribe'|ORCA MPI switches are collected, not appended one by one"
"src/calculator/orca_sc.f90|mpiargs(1:1)|extra MPI switches are emitted as ONE quoted argv entry (ORCA: 'too many arguments')"
"src/calculator/orca_sc.f90|CREST_ORCA_PRINT_CALL|opt-in echo of the exact ORCA command line"
# --- env switches are read by VALUE, not by mere presence.  The old tests were
# --- 'len_trim(envstr) > 0' / 'len_trim(envstr) == 0', which meant
# --- CREST_ORCA_PRINT_CALL=0 printed the line it was meant to suppress, and
# --- CREST_ORCA_NO_OVERSUBSCRIBE=0 dropped --oversubscribe.  Same convention as
# --- CREST_DYNAMIC_SHARE in parallel.f90.
"src/calculator/orca_sc.f90|read (envstr,*,iostat=io) ival|ORCA env switches are parsed as integers (was: any non-empty value counted as 'set')"
"src/calculator/orca_sc.f90|if (io == 0.and.ival /= 0.and..not.callprinted) then|CREST_ORCA_PRINT_CALL=0 really suppresses the echo"
"src/calculator/orca_sc.f90|if (io /= 0.or.ival == 0) then|CREST_ORCA_NO_OVERSUBSCRIBE=0 keeps --oversubscribe"
# --- "-iru" is actually wired up (upstream declares+parses it but never reads it)
"src/confparse.f90|-iru|'-iru' flag is parsed into env%iru"
"src/algos/parallel.f90|(env%iru.or.env%loadconformer).and.allocated(env%metadfac)|-iru/--loadconformer are READ (upstream: dead switch), guarded by the GF list"
"src/algos/parallel.f90|iruxyz = iruxyz/bohr|static bias coordinates are converted to Bohr"
"src/algos/parallel.f90|mddats(i)%npot = 2|-iru adds a 2nd potential instead of replacing the dynamic MTD"
"src/algos/parallel.f90|mddats(i)%mtd(2)%mtdtype = cv_rmsd_static|the 2nd potential is a STATIC RMSD bias"
"src/algos/parallel.f90|mddats(i)%mtd(2)%kpush = mtds(i)%kpush*irukscal|static bias uses the normal k (scalable via CREST_IRU_KSCALE)"
"src/algos/parallel.f90|CREST_IRU_KSCALE|optional k-scaling knob for the static bias"
"src/algos/parallel.f90|static bias (-iru)|per-job banner reports the static bias"
# --- ONE quoting style for every echoed input file: '*' banner + '* ' line prefix,
#     exactly like parse_datastruct.f90's root_print2() does for the TOML.
"src/calculator/subprocess_types.f90|'INPUT FILE ',trim(fname),' content (without comments):'|ORCA template echo uses the same '*' banner as the TOML echo"
"src/calculator/subprocess_types.f90|write (stdout,'(\"*\",1x,a)') trim(atmp)|ORCA template lines are quoted with '* ' (was '| ')"
"src/confparse.f90|'INPUT FILE ',trim(env%constraints),' content (sorted):'|constrain_input echo uses the same '*' banner as the TOML echo"
"src/confparse.f90|write (stdout,'(\"*\",1x,a)') trim(env%cts%sett(i))|constrain_input lines are quoted with '* ' (was '> ')"
# --- metal bonds are kept out of the CREGEN topology check (same rule as SHAKE)
"src/dynamics/shake_module.f90|crest_is_metal|SHAKE's private metal() table is exposed as a shared 'is this a metal' test"
"src/dynamics/shake_module.f90|if (z == 3.or.z == 4) return|Li/Be exception: treated exactly like init_shake() does"
"src/classes.f90|logical :: excludeTOPOmetals = .true.|metal bonds are excluded from the topology check BY DEFAULT"
"src/cregen.f90|use shake_module,only:crest_is_metal|cregen_topocheck uses the same metal definition as SHAKE"
"src/cregen.f90|call bondtotopo_excl(nat,at,bond,cn,ntopo,toporef,neighmat,topoexcl)|topology is built with the metal exclusion mask (all-false mask == plain bondtotopo)"
"src/cregen.f90|of which metals|the log reports how many atoms were excluded as metals"
"src/confparse.f90|-nometaltopo|CLI switch to go back to judging metal topology"
"src/parsing/parse_maindata.f90|case ('metaltopo')|TOML key metaltopo"
# --- topology-only pre-filter BEFORE the expensive refinement singlepoints
"src/classes.f90|logical :: refine_presort_topo = .true.|topology pre-filter before the refinement is ON by default"
"src/parsing/parse_maindata.f90|case ('refine_presort_topo')|TOML key refine_presort_topo"
"src/cregen.f90|subroutine cregen_topocheck(ch,env,checkez,nat,nall,at,xyz,comments,newnall)|cregen_topocheck is exposed via cregen_interface so crest_refine can call it"
"src/algos/refine.f90|call cregen_topocheck(stdout,env,env%checkiso,nat,nall,at,xyz,comments,nallnew)|crest_refine runs the topology check before the singlepoints"
"src/algos/refine.f90|CREST_REFINE_PRESORT_TOPO|CREST_REFINE_PRESORT_TOPO=0/=1 overrides the TOML value/default"
"src/algos/refine.f90|allocated(env%refine_queue).and.tp.and.env%checktopo|the pre-filter only runs when there IS a refinement (and -notopo disables it)"
# --- '--onebatch': one single MTD batch + MDs, no iterative MTD, EWIN=8
"src/classes.f90|logical :: onebatch = .false.|state flag for the --onebatch mode"
"src/confparse.f90|case ('-onebatch')|CLI switch --onebatch / -onebatch"
"src/confparse.f90|env%Maxrestart   = 1|--onebatch: exactly ONE MTD round"
"src/confparse.f90|env%allowrestart = .false.|--onebatch: elowcheck never reports 'lower', so no MTD re-iteration"
"src/confparse.f90|env%ewin         = 8.0d0|--onebatch: energy window widened to 8 kcal/mol"
"src/parsing/parse_maindata.f90|case ('onebatch')|TOML key onebatch"
"src/printouts.f90|one-batch mode|--onebatch is reported in the settings printout"
"src/printouts.f90|-onebatch          : perform ONE single batch of MTDs|--onebatch is documented in the help text"
# --- --onebatch runs a fixed 2x2 = 4-point k/alpha grid (one level below -quick's 6).
#     A single unrepeated batch has no second round that could compensate for a
#     badly placed Vbias, so the 4 points are a rectangle in (log k, log alpha)
#     rather than 4 alphas at one k: with one k, a wrong k wastes all four.
"src/choose_settings.f90|LOCAL EXTENSION (\"-onebatch\")|the 4-MTD grid override is documented where it lives"
"src/choose_settings.f90|if (env%onebatch.and.env%runver.ne.787878)|the override yields to an explicit -nmtd <N> (runver 787878)"
"src/choose_settings.f90|na = 2|two alpha levels: 1.00 and the second one (alp2)"
"src/choose_settings.f90|0.0020 / 0.0010|k = 0.0020 and 0.0010 (the -quick pair; 2.5x gentler than -nmtd's 0.005)"
"src/choose_settings.f90|0.83 / 1.02 Ang|alpha = 1.00 / 0.67 documented with its bias half-width sqrt(ln2/alpha)"
# --- (2026-09-25) the second alpha level of the onebatch grid: 0.50 -> 0.67.
# --- alpha sets the WIDTH of one Gaussian hill (half width sqrt(ln2/alpha)).
# --- 0.50 -> 1.18 A: a hill that wide floods a broad RMSD region after only a
# --- few deposits and pushes the structure out of a basin before it has been
# --- sampled.  0.67 -> 1.02 A buys more steps inside each basin at the same k.
"src/choose_settings.f90|alp2 = 0.67d0|default of the second alpha level (was 0.50)"
"src/choose_settings.f90|CREST_ONEBATCH_ALPHA2|the second alpha level is tunable without recompiling"
"src/choose_settings.f90|alpinc = 1.0d0/alp2|alpinc is a DIVISOR, so the grid gets alp/alpinc = alp2"
"src/printouts.f90|alpha = 1.00/0.67|-h states the actual grid instead of leaving it implicit"
"src/choose_settings.f90|over this override (runver 787878)|nrem is cleared: a leftover nrem from -nmtd would push nmtdynmax below 1"
"src/classes.f90|a fixed 2x2 = 4-point grid|the onebatch flag documents the 4-MTD grid"
"src/printouts.f90|1 round of 4 MTDs|the settings printout states how many MTDs one batch runs"
# --- '--loadconformer': user-supplied ensemble as a static RMSD bias (the -iru mechanism,
#     but read from readconformers.xyz in the CWD, i.e. active from the FIRST MTD on).
"src/classes.f90|loadconformerfile = 'readconformers.xyz'|the file name is a public parameter next to conformerfile"
"src/classes.f90|logical :: loadconformer = .false.|state flag for the --loadconformer mode"
"src/classes.f90|from the very first MTD|-iru can only bias from iteration 2 on, --loadconformer from the 1st - documented"
"src/confparse.f90|case ('-loadconformer','-loadconf')|CLI switch --loadconformer"
"src/confparse.f90|inquire (file=loadconformerfile,exist=ex)|a missing readconformers.xyz is caught at parse time"
"src/confparse.f90|env%loadconformer = .true.|the flag is only set when the file is actually there"
"src/parsing/parse_maindata.f90|case ('loadconformer','loadconf')|TOML key loadconformer"
"src/algos/parallel.f90|if ((env%iru.or.env%loadconformer).and.allocated(env%metadfac))|the static-bias block is entered for either flag"
"src/algos/parallel.f90|irufile = loadconformerfile|--loadconformer wins over -iru for the bias source"
"src/algos/parallel.f90|call get_environment_variable('CREST_LOADCONF_KSCALE'|a user ensemble can be scaled separately from the auto-generated one"
"src/algos/parallel.f90|call rdensembleparam(irufile,irunat,irunall)|the ensemble size is read from the selected file"
"src/algos/parallel.f90|call rdensemble(irufile,irunat,irunall|the ensemble itself is read from the selected file"
"src/algos/parallel.f90|iruxyz = iruxyz/bohr|cvxyz must still be in Bohr for the static bias"
"src/algos/parallel.f90|cvtype(2) = cv_rmsd_static|the 2nd potential is the static RMSD bias"
"src/printouts.f90|-loadconformer     : read|the option is documented in the help text"
"src/printouts.f90|static RMSD bias|the settings printout says which ensemble biases the run"
# --- the silent initial geometry optimization announces itself
"src/calculator/orca_sc.f90|optstartprinted|'optimization starting' printed once, right after the first external call"
"src/calculator/orca_sc.f90|'optimization starting ...'|the long, previously silent initial optimization is marked in the log"
# --- live output: flush instead of leaving half-written records in the buffer
"src/dynamics/dynamics_module.f90|flush (trj)|MD trajectory is flushed after EVERY dumped frame (no half-written last frame)"
"src/dynamics/dynamics_module.f90|flush (stdout)|stdout flushed on the same cadence, so a redirected crest.out stays live during a silent (pr=.false.) MTD"
"src/dynamics/dynamics_module.f90|Without this the MD table only shows up in blocks of one buffer size|MD table is flushed right after it is printed"
"src/sigterm.f90|flush (stdout)|creststop() flushes stdout before exit(3), otherwise the terminating message can be missing"
# --- --cinp bond constraints are RIGID (SHAKE) by default
"src/classes.f90|logical :: hardcons = .true.|state flag: -cinp bonds are rigid by default"
"src/classes.f90|integer :: nhardcons = 0|bonds read from the constraint FILE are kept separately"
"src/parsing/parse_xtbinput.f90|call add_hardcons(env,i1,i2,dist)|every bond from the constraint file is recorded"
# NOTE: the pairwise branch (atoms:/elements:) builds its constraints much later
# (:300-329) and used to skip the registration entirely, so such files stayed
# harmonic.  Keep this one - a file written with "atoms:" is the common case.
"src/parsing/parse_xtbinput.f90|call add_hardcons(env,j,i,dist)|atoms:/elements: pairwise distance constraints are registered as RIGID too"
"src/parsing/parse_xtbinput.f90|subroutine add_hardcons(env,i,j,dist)|the recorder itself (chunked growth of the list)"
"src/dynamics/dynamics_module.f90|subroutine crest_hard_bonds_to_shake(shk,n,atms,dist,pr)|bond list -> SHAKE user bonds"
# NOTE: these two MUST stay separate.  One "allocate(a,b,source=0)" makes the
# source-expr conformability check apply to BOTH objects, and conslistu is
# integer while distconsu is real(wp).  gfortran then dies with
#   Error: Type of entity at (1) is type incompatible with source-expr at (2)
"src/dynamics/dynamics_module.f90|allocate (shk%conslistu(2,n),source=0)|integer pair list allocated on its own (NOT sharing a source= with the real array)"
"src/dynamics/dynamics_module.f90|allocate (shk%distconsu(n),source=0.0_wp)|real target-distance list allocated on its own"
"src/dynamics/dynamics_module.f90|note: these bonds tie|redundant rigid-bond set is flagged up front (atoms: a,b,c,... is always redundant)"
"src/dynamics/dynamics_module.f90|If SHAKE stalls, narrow the atom list|the redundancy note keeps its actionable hint"
"src/dynamics/dynamics_module.f90|'rigid (SHAKE) bonds taken from the constraint input: '|the SHAKE report no longer prefixes a blank line"
"src/dynamics/shake_module.f90|real(wp),allocatable :: distconsu(:)|user bonds carry an explicit target distance"
"src/dynamics/shake_module.f90|shk%distcons(i) = shk%distconsu(i)*shk%distconsu(i)|explicit target distance wins over the input geometry"
"src/dynamics/shake_module.f90|shk%ncons = shk%nusr+nadd|auto-detected bonds already in the user list are not constrained twice"
"src/choose_settings.f90|call crest_hard_bonds_to_shake(env%mddat%shk,env%nhardcons,|the file bonds are converted into SHAKE bonds (now inside crest_hard_bonds_setup)"
"src/choose_settings.f90|if (env%mddat%shk%nusr > 0) env%mddat%shake = .true.|rigid bonds force SHAKE on, even with -shake 0"
"src/algos/parallel.f90|shk%distconsu = env%mddat%shk%distconsu|the rigid bonds survive 'mddat%shk = shk' in the parallel MTD setup"
"src/algos/setuptest.f90|MD%shake = (shakemode > 0).or.hardcons|the trial MTD keeps SHAKE on for rigid bonds"
"src/algos/setuptest.f90|if (tstep <= 1.0d0.and.hardcons) then|the trial MTD never trades the rigid bonds for a shorter time step"
"src/confparse.f90|case ('-softcons','-soft')|CLI switch to keep -cinp bonds harmonic"
"src/parsing/parse_maindata.f90|case ('hardconstraints','hardcons')|TOML key hardconstraints"
"src/printouts.f90|-softcons          : keep the -cinp bonds as harmonic restraints instead|-softcons is documented in the help text"
# --- (2026-09-22) P1+P2: MD failures must SAY WHAT failed --------------------
# Background: dynamics() has exactly two early exits - engrad() failing
# (dynamics_module.f90:340) and do_shake() failing (:504) - and both were
# reported through the same integer with no text at all.  The trial MTD runs
# with pr=.false. (setuptest.f90:164) while the SHAKE message was behind
# "else if (pr)", so neither diagnosis ever reached the log.
"src/dynamics/shake_module.f90|call shake_failreport(shk,maxdev,maxdev0,min(icyc,shk%maxcyc))|SHAKE non-convergence is reported even with pr=.false. (trial MTD runs silent)"
"src/dynamics/shake_module.f90|subroutine shake_failreport(shk,maxdev,maxdev0,icyc)|per-bond SHAKE report: iterations used, tolerance, worst offenders"
"src/dynamics/shake_module.f90|constraints violated        :|the report says how many constraints are still above the tolerance"
"src/dynamics/dynamics_module.f90|call md_failreport('SINGLEPOINT (energy/gradient)',t,dat%length_steps,term)|a failed singlepoint is named as such, not just a bare 'not converged'"
"src/dynamics/dynamics_module.f90|call md_failreport('SHAKE (rigid constraints)',t,dat%length_steps,term)|a failed SHAKE step is named as such, not just a bare 'not converged'"
"src/dynamics/dynamics_module.f90|subroutine md_failreport(origin,t,nsteps,termstatus)|shared failure report prints the step number and the termination status"
"src/dynamics/dynamics_module.f90|Failing in the FIRST step means the starting geometry|a failure in step 1 points at coord vs coord.ref, not at the time step"
# --- (2026-09-23) the first real failure showed that the FINAL deviation alone
# --- cannot tell an unreachable target from a DIVERGING solver, and that f9.5 /
# --- e12.4 overflow into '*********' and a dropped exponent letter.
"src/dynamics/shake_module.f90|real(wp) :: maxdev0|the deviation BEFORE the first correction is captured"
"src/dynamics/shake_module.f90|if (icyc == 0) maxdev0 = maxdev|captured on the first pass only (icyc is incremented at the end of the body)"
"src/dynamics/shake_module.f90|deviation BEFORE correcting :|the report prints the pre-correction deviation, not only the final one"
"src/dynamics/shake_module.f90|DIAGNOSIS: the step arrived essentially ON the constraints|auto-diagnosis: arrived on-target then blew up == solver instability"
"src/dynamics/shake_module.f90|es12.4,2x,a,es12.4|distances are printed with 'es', never 'f' - f9.5 overflowed to ********* on the first real failure"
# --- (2026-09-23) P4: the SHAKE sweep must be Gauss-Seidel, not Jacobi
"src/dynamics/shake_module.f90|LOCAL EXTENSION (P4, 2026-09-23): Gauss-Seidel instead of Jacobi|P4 rationale is documented where the sweep lives"
"src/dynamics/shake_module.f90|P4: refresh the linearisation at the CURRENT positions|P4: dr is refreshed in the non-frozen branch"
"src/dynamics/shake_module.f90|P4: same refresh here|P4: dr is refreshed in the frozen branch too"
# --- (2026-09-23) P3-a: the budget was the real blocker, not feasibility ------
# Background: measured rho = 0.9437 -> dev(n) = 0.2521*0.9437**n reaches the
# 1e-7 tolerance at n ~ 254 while maxcyc was 250.  The run missed by 4 sweeps.
"src/dynamics/shake_module.f90|LOCAL EXTENSION (P3-a, 2026-09-23): maxcyc 250 -> 1000|P3-a rationale is documented where the default lives"
"src/dynamics/shake_module.f90|integer :: maxcyc = 1000|the SHAKE sweep budget is 1000, not 250 (measured need: ~254)"
"src/dynamics/shake_module.f90|4 SWEEPS SHORT OF THE BUDGET|the P3-a note records the measured 254-vs-250 shortfall"
# --- (2026-09-23) P3-c: the DIAGNOSIS used to blame the time step.  On the real
# --- failure the residual had come down SIX orders of magnitude and sat at
# --- 1.29x the tolerance, yet the maxdev0-gt-1e-3 branch printed: the step is
# --- too large - here a shorter time step DOES help.  Wrong: it ran out of
# --- cycles.  The branches are now keyed on maxdev vs tolshake, and the report
# --- prints the fitted contraction factor plus the sweeps that were missing.
# --- NOTE: every entry in this array MUST end with a closing double quote.
# --- A block where the trailing quote was silently dropped on every line still
# --- LOOKS right in an editor but makes bash fail with "syntax error near
# --- unexpected token |" a few lines further down.  Always `bash -n` the file
# --- after editing the array, and when the error points at a line that looks
# --- perfectly balanced, count the quotes on the lines ABOVE it.
"src/dynamics/shake_module.f90|nmore = log|P3-c: the missing sweep count is extrapolated from that fit"
"src/dynamics/shake_module.f90|contraction per sweep|P3-c: rho is printed, so a geometric decay is visible at a glance"
"src/dynamics/shake_module.f90|more sweeps needed to reach|P3-c: the report says how many sweeps it was short"
"src/dynamics/shake_module.f90|DIAGNOSIS: NOT a feasibility problem|P3-c: came within a factor of tol and ran out of cycles is its own branch"
"src/dynamics/shake_module.f90|out of CYCLES.|P3-c: that branch names the actual cause"
"src/dynamics/shake_module.f90|the residual is NOT shrinking any more|P3-c: a genuine stall is still distinguished from a budget overrun"
# --- (2026-09-23) ORCA trial-MTD default: 0.125 ps instead of 0.25 ps --------
# Background: the trial MTD only calibrates the time step and samples the wall
# time, so for an external SCF it is pure overhead.  The default now depends on
# the calculator, which needs a predicate - hence the new logical function.
"src/algos/parallel.f90|logical function crest_calc_is_orca(calc)|a predicate for driving ORCA, next to its subprocess sibling"
"src/algos/parallel.f90|crest_calc_is_orca = .true.|the predicate actually recognises ORCA"
"src/algos/parallel.f90|if (.not.allocated(calc%calcs)) return|the predicate must not crash on a calculator that is not set up"
"src/algos/setuptest.f90|logical,external :: crest_calc_is_orca|setuptest links the predicate the same way as its sibling"
"src/algos/setuptest.f90|if (crest_calc_is_orca(env%calc)) then|the trial length is chosen on the actual calculator"
"src/algos/setuptest.f90|trial_md_ps = 0.125_wp|the ORCA trial MTD default is 0.125 ps"
"src/algos/setuptest.f90|trial_md_ps = 0.25_wp|every other calculator keeps 0.25 ps"
"src/algos/setuptest.f90|length set by|the settings block prints WHERE the trial length came from"
"src/algos/setuptest.f90|trial_src = 'ORCA default'|the ORCA reason string is recorded (<= 17 chars: it is printed into a fixed 17-column field)"
# --- (2026-09-23) the trial's Vbias dump is a plain 50 fs, no longer scaled ---
# Scaling with the trial length was the wrong lever: the interval answers "how
# often does the bias adapt", which must not depend on the trial's length.
"src/algos/setuptest.f90|trial_vbdump_fs = 50.0_wp|the TRIAL Vbias dump default is 50 fs"
"src/algos/setuptest.f90|CREST_TRIAL_VBDUMP_FS|the trial Vbias dump stays overridable without a rebuild"
"src/algos/setuptest.f90|it is now a plain 50 fs|the rationale for dropping the scaling is recorded"
# --- (2026-09-24) dynamic core sharing: a job that finishes hands its slice to
# --- the jobs that are still running, instead of leaving 1/n of the node idle.
"src/algos/parallel.f90|module crest_dynshare|the dynamic core-sharing state lives next to the MD drivers"
"src/algos/parallel.f90|crest_dynshare_init(nsim,T,env%Threads,Tn)|the round registers its job count, slots, budget and starting share"
"src/algos/parallel.f90|CREST_DYNAMIC_SHARE|opt-out switch for the dynamic sharing"
"src/algos/parallel.f90|CREST_MAX_CORES_PER_JOB|caps what a single job may get (48 MPI ranks can lose against 12)"
"src/algos/parallel.f90|call ompenvset(ds_total)|the env entry is pre-sized so every later setenv() writes in place"
"src/algos/parallel.f90|nrunning = max(1,min(ds_nslots,ds_njobs-ds_nfin))|with waves a finished job is replaced, so only the drained queue frees cores"
"src/algos/parallel.f90|!\$omp critical (crest_dynshare_lock)|the finish counter and the setenv are serialised"
  "src/algos/parallel.f90|call crest_dynshare_release()|every finished job hands its cores over"
  "src/algos/parallel.f90|call crest_dynshare_reset()|the share is restored after the round"
# --- (2026-09-24) --cinp stays RIGID (SHAKE), -cbonds/-cmetal/... stay harmonic
# --- Background: get_xtb_constraint_block() is shared by BOTH readers.  The
# --- --cinp file holds the handful of bonds the user wrote by hand, while the
# --- auto-generated "bondlengths" holds ONE ENTRY PER BOND.  Registering the
# --- latter as rigid froze the whole molecule (or, for -cmetal, the entire
# --- first coordination sphere), so rigidity is now opt-in per file.
"src/parsing/parse_xtbinput.f90|logical,optional,intent(in) :: hard|the rigid-registration flag is opt-in on the public entry point"
"src/parsing/parse_xtbinput.f90|hardcons = .false.|the default is HARMONIC, i.e. plain upstream behaviour"
"src/parsing/parse_xtbinput.f90|if (present(hard)) hardcons = hard|an explicit hard= still wins over the default"
"src/parsing/parse_xtbinput.f90|call get_xtb_constraint_block(env,blk,hardcons)|the flag is threaded down into the \$constrain reader"
"src/parsing/parse_xtbinput.f90|logical,intent(in) :: hard|get_xtb_constraint_block receives it"
"src/parsing/parse_xtbinput.f90|if (hard) call add_hardcons(env,i1,i2,dist)|explicit 'distance:' entries only register as rigid when the file opted in"
"src/parsing/parse_xtbinput.f90|if (hard) call add_hardcons(env,j,i,dist)|the pairwise 'atoms:' path is guarded the same way"
"src/parsing/confparse2.f90|call parse_xtbinputfile(env,'bondlengths',hard=.false.)|-cbonds/-cmetal/-cheavy/-clight stay harmonic restraints"
"src/confparse.f90|call parse_xtbinputfile(env,env%constraints,hard=.true.)|--cinp / TOML 'constraints' opts into RIGID (SHAKE)"
# --- (2026-09-24) the "stop" file: end a long dynamics early and keep the results
# --- Background: a 6 ps MTD on an external calculator can run for days.  Killing
# --- the job loses everything, because nothing has been optimized or written yet.
# --- Creating a file named "stop" in the working directory now ends the sampling
# --- and hands what has been sampled so far to the next stage.
"src/dynamics/dynamics_module.f90|logical :: stopcheck = .true.|mddata gains a per-caller opt-out for the stop-file poll"
"src/dynamics/dynamics_module.f90|logical :: stopfile_initialized = .false.|process-wide 'was the stop file there at startup' reference"
"src/dynamics/dynamics_module.f90|call md_stopcheck_setup(sinterval,stopfile,stoppre)|the poll is armed once per dynamics() call"
"src/dynamics/dynamics_module.f90|call get_environment_variable('CREST_STOP_FILE',value=envstr,status=stat)|the file name is configurable"
"src/dynamics/dynamics_module.f90|call get_environment_variable('CREST_STOP_CHECK_SEC',value=envstr,status=stat)|the interval is configurable (seconds)"
"src/dynamics/dynamics_module.f90|if (rd == 0) interval = val|the interval is read BY VALUE, not by mere presence"
"src/dynamics/dynamics_module.f90|if (interval < 1.0_wp) interval = 1.0_wp|a sub-second interval degenerates into one stat() per step"
"src/dynamics/dynamics_module.f90|stopfile_seen0 = ex|a stop file left over from an EARLIER run is recorded as such"
"src/dynamics/dynamics_module.f90|call md_stopcheck(stopfile,sinterval,scount0,prate,stoppre,stopgo)|the poll sits at the top of the step, so the next engrad() is skipped"
"src/dynamics/dynamics_module.f90|hit = (ex.and.(.not.pre))|only a file that APPEARED while crest runs may trigger"
"src/dynamics/dynamics_module.f90|          term = 3|a deliberate stop needs its own termination status"
"src/dynamics/dynamics_module.f90|call md_stopreport(stopfile,t,dat%length_steps,dat%dumped)|and is reported even for a silent (pr=.false.) parallel MTD"
"src/algos/setuptest.f90|MD%stopcheck = .false.|the trial MTD opts out (it reads io/=0 as 'did not converge')"
"src/algos/dynamics.f90|case (3)|standalone --md: term=3 is a stop on request, not a failure"
"src/algos/parallel.f90|vz,' stopped (\"stop\" file)'|the per-job completion line names the reason"
"src/printouts.f90|CREST_STOP_CHECK_SEC|the two knobs are documented in -h"
# --- (2026-10-05) the honoured stop file is DELETED automatically, so it cannot
# --- cut the next MTD round short as well.
"src/dynamics/dynamics_module.f90|stopfile_hit = .true.|a stop is latched for the whole round, not only for the job that saw it"
"src/dynamics/dynamics_module.f90|subroutine md_stop_cleanup()|the removal happens in the caller, after the parallel region"
"src/dynamics/dynamics_module.f90|close (u,status='delete',iostat=ios)|deleted from Fortran, without pulling iomod's remove() into this module"
"src/dynamics/dynamics_module.f90|stopfile_name = trim(fname)|the name to delete is the one that was actually honoured"
"src/algos/parallel.f90|if (md_stop_seen()) then|both parallel MTD rounds do the housekeeping"
"src/algos/parallel.f90|Structures so far are handed on to the next stage. Stop file is automatically removed.|the round says what happened to the file"
# --- ...and it has to work with -scratch, where crest chdir()s away from the
# --- directory the user is watching.  The launch directory is captured in
# --- parseflags() (getcwd into env%homedir) before any chdir can happen.
"src/dynamics/dynamics_module.f90|character(len=512) :: stopfile_launchdir = ''|the directory crest was STARTED from is remembered"
"src/dynamics/dynamics_module.f90|public :: md_stopcheck_setdir|parseflags() can tell the poll about it"
"src/confparse.f90|call md_stopcheck_setdir(env%homedir)|and does so right after getcwd(), i.e. before scratch.f90 chdir()s"
"src/dynamics/dynamics_module.f90|logical function md_stopfile_exists(fname) result(ex)|one probe covering both directories"
"src/dynamics/dynamics_module.f90|ex = md_stopfile_exists(fname)|both call sites use it (trigger and the 'was it there at startup' reference)"
"src/dynamics/dynamics_module.f90|if (fname(1:1) /= '/') then|an ABSOLUTE name is never prefixed with the launch directory"
# --- (2026-09-26) energy-only single points: a stage that ranks by energy must
# --- not pay for a gradient.  crest_sploop ("> Singlepoint re-ranking for N
# --- structures") called engrad(), so ORCA computed a gradient on every
# --- structure and crest_sploop filled grads() and deallocated it without a
# --- single read.  For a DFT level with a solvent model that is a CP-SCF /
# --- Z-vector solve plus the CPCM gradient term, per structure, for nothing.
"src/calculator/calc_type.f90|logical :: nogr = .false.|calcdata gains an 'this call does not need a gradient' flag"
"src/calculator/calculator.F90|module procedure :: singlepoint_mol|singlepoint() is a real public API (generic interface)"
"src/calculator/calculator.F90|subroutine singlepoint_mol(mol,calc,energy,iostatus)|the energy-only entry point (raises calc%nogr, reuses engrad_mol)"
"src/calculator/calculator.F90|if (calc%nogr) then|potential_core() dispatches on the flag"
"src/calculator/calculator.F90|call ORCA_sp(molptr,calc%calcs(id),calc%etmp(id),iostatus)|ORCA is routed to the gradient-free variant"
"src/calculator/orca_sc.f90|public :: ORCA_sp|ORCA_sp is exported from orca_sc"
"src/calculator/orca_sc.f90|subroutine ORCA_sp(mol,calc,energy,iostatus)|ORCA single point without the gradient"
"src/calculator/orca_sc.f90|call ORCA_setup(mol,calc,nogr=.true.)|ORCA_sp writes the template WITHOUT EnGrad"
"src/calculator/orca_sc.f90|logical,intent(in),optional :: nogr|ORCA_setup() takes the flag as an OPTIONAL arg (every old call site unchanged)"
"src/calculator/orca_sc.f90|subroutine rd_ORCA_energy(mol,calc,energy,iostatus)|the energy-only reader"
"src/calculator/orca_sc.f90|KEY = 'FINAL SINGLE POINT ENERGY'|energy is taken from the FINAL SINGLE POINT ENERGY line of ORCA.out"
"src/calculator/orca_sc.f90|call rd_ORCA_energy(mol,calc,energy,iostatus)|ORCA_sp calls it instead of rd_ORCA_engrad"
"src/calculator/subprocess_engrad.f90|public :: ORCA_sp|the re-export module forwards ORCA_sp"
"src/calculator/subprocess_types.f90|character(len=6),parameter :: ENG = 'engrad'|the EnGrad keyword is handled in one place"
"src/calculator/subprocess_types.f90|pure function has_engrad(line) result(found)|CASE-INSENSITIVE EnGrad recognition"
"src/calculator/subprocess_types.f90|function strip_engrad(line) result(out)|EnGrad can be removed from the template on write"
"src/calculator/subprocess_types.f90|if (sponly.and.has_engrad(line)) line = strip_engrad(line)|write_orca_input() strips EnGrad iff nogr was passed"
"src/calculator/subprocess_types.f90|logical,intent(in),optional :: nogr|write_orca_input() takes it as an OPTIONAL arg"
"src/algos/parallel.f90|calculations(i)%nogr = .true.|crest_sploop marks its per-thread calculator copies as energy-only"
"src/algos/parallel.f90|call singlepoint(mols(job),calculations(job),energy,io)|crest_sploop calls singlepoint(), not engrad()"
"src/algos/parallel.f90|see the calc%nogr note in the setup loop above|the comment says why the gradient is gone"
# --- (2026-09-26) upstream argument-count bug: singlepoint.f90 called
# ---   call crest_sploop(env,nat,nall,at,xyz,eread,.true.)   <- SEVEN actual args
# --- against a routine with SIX dummies (refine.f90 passes six, and the pristine
# --- 3.0.2 tree has the same six-dummy declaration).  It only "worked" because
# --- the x86-64 SysV ABI lets the callee ignore the extra pointer.
"src/algos/singlepoint.f90|LOCAL FIX: this call used to pass a SEVENTH argument|the 7th actual argument of the crest_sploop call is documented and gone"
"src/algos/singlepoint.f90|call crest_sploop(env,nat,nall,at,xyz,eread)|singlepoint.f90 calls crest_sploop with exactly its 6 dummy arguments"
# --- (2026-09-27) performance batch: kill duplicated external calls (T1)
"src/optimize/ancopt.f90|LOCAL CHANGE (2026-09-26): do NOT refetch on the very first iteration|ancopt/relax reuses the caller's initial gradient on its very first pass"
"src/optimize/rfo.f90|LOCAL CHANGE (2026-09-26): same guard as in ancopt/relax|rfo has the same first-pass guard"
"src/algos/parallel.f90|if (.not.subproc) then|the per-slot warm-up engrad is skipped for subprocess calculators"
"src/algos/parallel.f90|call singlepoint(mol,calc,energy,io)|the WBO-only call asks for an energy only (no gradient)"
"src/algos/dynamics.f90|call singlepoint(mol,calc,energy,io)|standalone MD: WBO-only call asks for an energy only"
"src/algos/setuptest.f90|call singlepoint(mol,tmpcalc,energy,io)|trial MD: WBO-only call asks for an energy only"
"src/algos/numhess.f90|call singlepoint(mol,calc,energy,io)|numhess initial single point asks for an energy only"
"src/algos/numhess.f90|letting singlepoint() leave grdtmp zeroed|numhess documents why grdtmp zeroing is safe"
# --- (2026-09-27) performance batch: file/memory round trips (T2)
"src/restartlog.f90|CREST_RESTART_SAVE|the continuous restart snapshot can be switched off"
"src/restartlog.f90|logical function restart_save_enabled|the switch is read once and cached"
"src/iomod.F90|position='append'|appendto() appends instead of scanning to EOF (was O(R^2))"
"src/strucreader.f90|subroutine count_xyz_structures|cheap nat/nall counter that parses no coordinate"
"src/strucreader.f90|public :: count_xyz_structures|the counter is exported"
"src/legacy_algos/confscript2_misc.f90|only:rdensembleparam,rdensemble,count_xyz_structures|sort_and_check imports the counter"
"src/legacy_algos/confscript2_misc.f90|call count_xyz_structures(trim(filename),env%nat,nallin)|sort_and_check counts without a full parse"
"src/legacy_algos/confscript2_misc.f90|call count_xyz_structures(trim(filename),nat,nall)|remaining_in() needs ONE pass instead of two"
"src/algos/search_conformers.f90|call count_xyz_structures(trim(inpnam),nat,nall)|multilevel loop counts without a full parse"
"src/algos/parallel.f90|call remove(atmp)|a merged per-worker trajectory is deleted right away (no duplicate copy on disk)"
# --- (2026-09-27) performance batch: pure-CPU waste (T3)
"src/ls_rmsd.f90|LOCAL CHANGE (2026-09-27): matmul(transpose(U),y(:,i)) depends only on i|the RMSD gradient kernel hoists the matmul out of the j loop"
"src/dynamics/metadynamics_module.f90|LOCAL CHANGE (2026-09-27): the compact copy of the CURRENT structure|the xyzcp subset copy is built once per call, not per reference"
"src/dynamics/dynamics_module.f90|real(wp),allocatable :: grdmtd(:,:)|the MTD bias gradient scratch belongs to dynamics()"
"src/dynamics/dynamics_module.f90|subroutine md_calc_mtd(mol,dat,epot,grd,grdmtd,pr)|md_calc_mtd takes the scratch instead of allocating it"
"src/dynamics/dynamics_module.f90|private to this thread's invocation|the per-step molo copy is no longer inside a global critical"
"src/cregen.f90|logical :: esorted|the energy-window scans verify the sorted-er invariant"
"src/cregen.f90|do while (jmin < i)|the rmap scan uses a monotonically advancing lower bound"
"src/cregen.f90|do while (jmin2 < i)|the duplicate scan uses the same window as the rmap scan"
"src/cregen.f90|real(wp),allocatable :: xyzmask(:,:,:)|heavy-atom sub-structures are cached instead of rebuilt per pair"
"src/cregen.f90|real(wp),allocatable :: xyzmask0(:,:,:)|sub-structure coordinates are cached instead of rebuilt per pair"
# --- (2026-09-27) conformer-search "tight" stage: own threshold level, 10x looser
"src/optimize/optutils.f90|integer,parameter :: olev_tightcs = 4|a separate level for the conformer-search tight stage (keeps the ladder monotonic)"
"src/optimize/optutils.f90|case (olev_tightcs)|get_optthr knows the new level"
"src/optimize/optutils.f90|ethr = 1.d-4  !> olev_tightcs|the new level's energy criterion is 1e-4 Eh (olev_tight is 1e-6)"
"src/optimize/optutils.f90|gthr = 4.d-3  !> olev_tightcs|the new level's gradient criterion is 4e-3 Eh/a0 (olev_tight is 8e-4)"
"src/optimize/optutils.f90|1e-6 / 8e-4 -> 1e-5 / 8e-3 -> 5e-5 / 1.6e-2 -> 5e-5 / 4e-3 -> 1e-4 / 4e-3|the tuning history of the level is recorded (and matches the values above)"
"src/algos/search_conformers.f90|use optimize_utils,only:olev_tightcs|the multilevel chain imports the new level (optimize_module is private)"
"src/algos/search_conformers.f90|env%calc%optlev = olev_tightcs|'optimization with tight thresholds' now uses it"
# --- (2026-09-27) global optimisation budget 120 + keep over-limit structures
"src/optimize/optutils.f90|integer,parameter :: opt_maxcycle_global = 120|the optimisation cycle budget is one explicit global constant"
"src/optimize/optutils.f90|maxcycle = opt_maxcycle_global|get_optthr uses the global budget instead of the per-level n*K values"
"src/calculator/calc_type.f90|logical  :: anopt = .true.|structures that exhaust maxcycle are KEPT, not discarded as failures"
# --- (2026-09-27) tuned CREGEN dedup defaults: EWIN 8, RTHR 0.15, ETHR 0.15
"src/confparse.f90|env%ewin = 8.0d0    !> EWIN - energy window in kcal/mol|default EWIN is 8.0 kcal/mol (was 6.0)"
"src/confparse.f90|env%rthr = 0.15d0   !> RTHR - RMSD thr in Angstroem|default RTHR is 0.15 Ang (was 0.125)"
"src/confparse.f90|env%ethr = 0.15d0   !> ETHR - E threshold in kcal|default ETHR is 0.15 kcal/mol (was 0.05)"
"src/confparse.f90|stage 1 now uses 0.300 Ang|the multiplicative effect on the multilevel rthr factors is recorded in the comment"
"src/algos/search_conformers.f90|Keep the single source of|the multilevel default branch no longer hardcodes ewin/rthr"
# --- (2026-09-28) flexi.F90: NTOP (parameter) vs ntop (counter) were the SAME
# --- symbol -> Intel build error "already has basic type of INTEGER".  The size
# --- parameter is now NTOPLIST so the running counter can keep its name.
"src/flexi.F90|integer,parameter :: NTOPLIST = 20|the bond-list size parameter got a name that cannot collide with the counter"
"src/flexi.F90|integer :: ti(NTOPLIST),tj(NTOPLIST)|the arrays use the renamed parameter"
"src/flexi.F90|if (ntop .lt. NTOPLIST) then|the insert-branch bound uses the renamed parameter"
# --- (2026-09-29) The flexibility diagnostic table must quote FULL (input-file)
# --- atom numbers.  It used to print the REDUCED index, so with nat=97/rednat=92
# --- phosphorus at 50/51 was reported as 47/48 and the offset varied with
# --- position -> impossible to relate the table back to one's own xyz file.
"src/flexi.F90|integer,allocatable::irev(:)|the reverse (reduced -> full) index map exists"
"src/flexi.F90|irev(k) = i|the reverse map is filled from includeAtom()"
"src/flexi.F90|asym(at(i)),irev(i),asym(at(j)),irev(j)|the table prints FULL atom numbers (both endpoints)"
# --- (2026-09-30) g-xTB must reach CREST through the xtb SUBPROCESS with the
# --- --gxtb flag (analytical .engrad gradient), NOT through the TURBOMOLE
# --- channel (which upstream chose for the first, gradient-less g-xTB release).
"src/calculator/calc_type.f90|self%id = jobtype%xtbsys|g-xTB is routed to the xtb subprocess, not to turbomole"
"src/calculator/calc_type.f90|self%other = '--gxtb'|the g-xTB method flag is injected into the xtb call"
"src/calculator/calc_type.f90|self%binary = 'xtb'|g-xTB uses the ordinary xtb binary (name via PATH)"
"src/parsing/parse_calcdata.f90|job%other = '--gxtb'|the TOML method=gxtb_dev path agrees with the CLI path"
"src/confparse.f90|        call gxtb_dev_warning()|-gxtb still prints the notice"
"src/confparse.f90|      case ('-gxtb_dev')|-gxtb_dev still exists as its own case"
"src/printouts.f90|ANALYTICAL gradient) - NOT through the TURBOMOLE interface|the obsolete 'numerical gradient' notice is corrected"
# --- (2026-09-30) g-xTB solvation passthrough ($CREST_GXTB_SOLV -> calc%other).
# --- '-g'/'-alpb' only reach the in-process tblite path; on the xtb subprocess
# --- channel they were dropped SILENTLY, so the flag has to be appended to
# --- calc%other from the environment.  Both doors (CLI shortcut + TOML) must
# --- carry it, hence one anchor per file.
"src/calculator/calc_type.f90|call getenv('CREST_GXTB_SOLV',sv)|the CLI path forwards \$CREST_GXTB_SOLV to the xtb command line"
"src/calculator/calc_type.f90|clean = clean//sv(k:k)|runs of blanks are collapsed before the flags reach the shell"
"src/parsing/parse_calcdata.f90|subroutine append_gxtb_solv(other)|the TOML path has a shared --gbe/--cosmo helper"
"src/parsing/parse_calcdata.f90|call append_gxtb_solv(job%other)|the TOML method=gxtb_dev case actually calls it"
"src/parsing/parse_calcdata.f90|if (index(job%other,'--gxtb') /= 0) call append_gxtb_solv(job%other)|the TOML 'flags' key overwrites job%other, so the passthrough is re-added there"
# --- (2026-09-29) GC output must always go through stage 5, stage 6 is defined
# --- as HALF of stage 5, and CREGEN rejects are collected into an audit file.
"src/algos/search_conformers.f90|multilevel(5) = .true.|GC output always passes stage 5 (tight) before the final wrap"
"src/algos/search_conformers.f90|env%calc%optlev = olev_vtightcs|stage 6 is defined relative to stage 5, not as the absolute olev_vtight"
"src/algos/search_conformers.f90|use optimize_utils,only:olev_tightcs,olev_vtightcs|both conformer-search levels are imported"
"src/optimize/optutils.f90|integer,parameter :: olev_vtightcs = 5|stage 6's level exists as its own constant"
"src/optimize/optutils.f90|case (olev_vtightcs)|get_optthr knows stage 6's level"
"src/optimize/optutils.f90|ethr = 5.d-5  !> olev_vtightcs = olev_tightcs / 2|stage 6 energy criterion is exactly half of stage 5's"
"src/optimize/optutils.f90|gthr = 2.d-3  !> olev_vtightcs = olev_tightcs / 2|stage 6 gradient criterion is exactly half of stage 5's"
"src/optimize/optutils.f90|write (*,chrfmt) \"optimization level\",'tightcs (4)'|stage 5 is named, not printed as 'custom'"
"src/optimize/optutils.f90|write (*,chrfmt) \"optimization level\",'vtightcs (5)'|stage 6 is named, not printed as 'custom'"
# --- CREGEN reject audit (new module + its four call sites + the plumbing)
"src/exclfile.f90|module exclfile_module|the reject-audit module exists"
"src/exclfile.f90|character(len=*),parameter :: exclfile_name = 'crest_excluded_structures.xyz'|the audit file name"
"src/exclfile.f90|call get_environment_variable('CREST_EXCL_XYZ',value=val,status=estat)|the audit is opt-in via CREST_EXCL_XYZ"
"src/exclfile.f90|subroutine exclfile_write(nat,at,xyz,comment,tag,reason,factor)|the single append entry point"
"src/cregen.f90|call exclfile_check()|newcregen honours CREST_EXCL_XYZ"
"src/cregen.f90|call exclfile_init()|the audit file is truncated once per run"
"src/cregen.f90|'CLASH', &|discardbroken() auditors rejections"
"src/cregen.f90|'TOPO', &|cregen_topocheck() auditors rejections"
"src/cregen.f90|'EWIN', &|cregen_esort() auditors rejections"
"src/cregen.f90|'DUP',trim(atmp_e),1.0_wp|cregen_CRE() auditors duplicates"
"src/cregen.f90|logical :: ezfail   !> LOCAL EXTENSION (2026-09-29)|the topology gate can tell E/Z apart in its reason text"
"src/cregen.f90|subroutine cregen_esort(ch,nat,nall,at,xyz,comments,nallout,ewin)|cregen_esort takes the atom list (the audit writer needs it)"
# --- (2026-10-05) the topology PRE-FILTER (algos/refine.f90) calls
# --- cregen_topocheck() BEFORE newcregen() ever arms the audit, so its rejects
# --- were dropped silently: "82 of 87 survive the topology pre-filter" with an
# --- EMPTY crest_excluded_structures.xyz.  It now arms the audit itself.
"src/algos/refine.f90|use exclfile_module,only:exclfile_check,exclfile_init|the topology pre-filter arms the reject audit itself"
"src/algos/refine.f90|call exclfile_init()|and does so BEFORE cregen_topocheck() discards anything"
# --- (2026-10-02) gradient descent: the "FAILED to perform gradient descent!
# --- Something is very wrong ..." bail-out.  The ERISE halving loop can only
# --- succeed while the step is still representable; once alp*|g| underflows the
# --- geometry, ORCA returns a bitwise identical energy, echng is EXACTLY 0.0 and
# --- the strict `lowered` test can never accept it -> ~60 wasted DFT gradients
# --- per failed structure and a good structure thrown away.  Four anchors:
# --- frozen search direction, displacement trust radius, step floor, converged
# --- veto, plus the bbstep division guard.
"src/optimize/gd.f90|real(wp),allocatable :: grd0(:,:),grdt(:,:)|the trial gradient no longer overwrites the search direction"
"src/optimize/gd.f90|grd0 = grd|the search direction is frozen once per iteration"
"src/optimize/gd.f90|g0norm = norm2(grd0)|the trust radius is scaled by the frozen gradient norm"
"src/optimize/gd.f90|alp = min(alp,maxdispl/max(g0norm,1.0e-30_wp))|maxdispl caps the DISPLACEMENT in bohr, not the alp scalar"
"src/optimize/gd.f90|real(wp),parameter :: steptol = 1.0e-10_wp|the halving loop has an explicit step-size floor"
"src/optimize/gd.f90|if (alp*g0norm .lt. steptol) then|the floor is actually tested before every trial"
"src/optimize/gd.f90|if (.not.lowered.and..not.converged) then|a structure that already meets both criteria is not rewound"
"src/optimize/gd.f90|grd = grdt|an accepted step re-points grd at the new geometry"
"src/optimize/gd.f90|if (dum2 .lt. 1.0e-14_wp.or.abs(dum1) .lt. 1.0e-14_wp) then|bbstep cannot divide by a vanishing gradient change"
# --- (2026-10-04) log-format cleanups.  Each anchor pins one of the removals,
# --- because "the line is gone" is exactly what a later careless edit would
# --- put back.
"src/printouts.f90|do i = 1,12|the banner box has no blank line above its bottom border any more"
"src/printouts.f90|3x,''• P.Pracht, F.Bohle|the citation block starts directly under its heading"
"src/algos/search_conformers.f90|write (stdout,*) 'Input structure:'|no blank line between the SAMPLING box and the input structure"
"src/dynamics/dynamics_module.f90|LOCAL CHANGE (2026-10-04): no leading blank line|the SHAKE report does not prefix a blank line"
"src/algos/parallel.f90|subroutine crest_dynshare_barmark()|an open progress bar can be closed before the [cores] line"
"src/algos/parallel.f90|if (.not.ds_announced) then|the core reallocation is reported once per round, not once per job"
"src/algos/parallel.f90|call parallel_md_block_printout(mddats(i),i,nsim)|MTD startup boxes are printed in ascending order, outside the OMP region"
"src/optimize/optimize_module.f90|' (Hessian update type: '|engine and Hessian update share one line"
"src/cleanup.f90|call exclfile_free()|the audit file is closed at run end"
"src/CMakeLists.txt|exclfile.f90|the new file is in the CMake source list"
"src/meson.build|'exclfile.f90'|the new file is in the meson source list"
# --- (2026-09-29) L-BFGS is a real optimizer now, and the Final Geometry
# --- Optimization defaults to the gradient-descent engine.
"src/optimize/lbfgs.f90|module lbfgs_module|the L-BFGS engine exists as its own module"
"src/optimize/lbfgs.f90|subroutine lbfgs(mol,calc,etot,grd,pr,wr,iostatus)|its entry point matches the ancopt/gd/rfo interface exactly"
"src/optimize/lbfgs.f90|subroutine lbfgs_two_loop(n,m,npair,ipos,s,y,rho,gamma,g,d,a,maxdispl)|the two-loop recursion is factored out and testable"
"src/optimize/lbfgs.f90|subroutine lbfgs_update(n,m,npair,ipos,s,y,rho,snew,ynew,gamma)|the pair/ring-buffer update is factored out and testable"
"src/optimize/lbfgs.f90|integer,parameter :: lbfgs_memory = 10|lbfgs keeps m=10 correction pairs"
"src/optimize/lbfgs.f90|real(wp),parameter :: lbfgs_curv_tol = 1.0d-10|the curvature gate is relative, hence unit-independent"
"src/optimize/lbfgs.f90|if (ys .le. lbfgs_curv_tol*sqrt(ss*yy)) return|a non-positive-curvature pair is DROPPED, not stored"
"src/optimize/lbfgs.f90|dtrial = d*0.5_wp**(jjj-1)|the line search shrinks a separate trial step, never d itself"
"src/optimize/lbfgs.f90|converged = econverged.and.gconverged.and.lowered|convergence test matches ancopt (it includes the 'lowered' flag)"
"src/optimize/optimize_module.f90|use lbfgs_module|optimize_module links the L-BFGS engine"
"src/optimize/optimize_module.f90|call lbfgs(molnew,calc,etot,grd,pr,wr,iostatus)|engine 1 dispatches to the real L-BFGS"
"src/algos/search_conformers.f90|env%calc%opt_engine = -1|the Final Geometry Optimization (stage 6) defaults to gradient descent"
"src/algos/search_conformers.f90|integer :: enginebackup|the engine is backed up before the multilevel loop"
"src/algos/search_conformers.f90|enginebackup   = env%calc%opt_engine|... captured from the user's setting"
"src/algos/search_conformers.f90|env%calc%opt_engine = enginebackup|... and restored after every stage, so stages 1-5 are unaffected"
"src/optimize/CMakeLists.txt|lbfgs.f90|the new file is in the optimize CMake source list"
"src/optimize/meson.build|'lbfgs.f90'|the new file is in the optimize meson source list"
# --- (2026-09-29) flexi() gained OPTIONAL dummies (nfix/fixat) but is an
# --- EXTERNAL subroutine, so every call site needs an explicit interface.
# --- Without one gfortran refuses to compile:
# ---   "Explicit interface required for 'flexi' at (1): optional argument"
# --- and mis-reports it as "Missing actual argument for argument '_formal_7'"
# --- at the 4-argument call sites.  Both callers must carry the block.
"src/flexi.F90|integer,intent(in),optional :: nfix|flexi() takes the declared bonds as optional dummies"
"src/choose_settings.f90|subroutine flexi(mol,rednat,includeAtom,flex,nfix,fixat)|md_length_setup has an explicit interface for flexi()"
"src/choose_settings.f90|call flexi( mol, env%rednat, env%includeRMSD, flex, env%nhardcons, env%hcatms)|the --cinp path passes the declared bonds"
"src/minitools.f90|subroutine flexi(mol,rednat,includeAtom,flex,nfix,fixat)|testtopo has an explicit interface for flexi() too"
"src/minitools.f90|call flexi(mol,zmol%nat,inc,flex)|... and the 4-argument call is still legal"
)

for c in "${CHECKS[@]}"; do
  f="${c%%|*}"; rest="${c#*|}"; pat="${rest%%|*}"; desc="${rest#*|}"
  if [ ! -f "$f" ]; then
    printf '  MISSING-FILE %-42s %s\n' "$f" "$desc"; fail=$((fail + 1))
  # -F is REQUIRED: the patterns are literal source snippets and some of them
  # contain '*' (e.g. "550.0_wp*trial_md_ps").  Plain grep would treat that as
  # a regex quantifier and report a false FAIL on a perfectly patched tree.
  elif grep -qF -- "$pat" "$f"; then
    printf '  OK           %-42s %s\n' "$f" "$desc"
  else
    printf '  FAIL         %-42s %s\n' "$f" "$desc"; fail=$((fail + 1))
  fi
done

# --- NEGATIVE checks: things that must be GONE.  The loop above can only assert
# --- presence, and silently reintroducing the old behaviour is exactly the kind
# --- of regression it would never catch.
for c in \
  "src/dynamics/shake_module.f90|'SHAKE did not converge! maxdev='|the silent pr-gated SHAKE message (must be gone: the trial MTD runs with pr=.false.)" \
  "src/dynamics/shake_module.f90|maxcyc = 250|the old 250-sweep budget (must be gone: it was 4 sweeps short on a real system)" \
  "src/dynamics/shake_module.f90|time step DOES help.|the WRONG P3-c advice (must be gone: the residual had come down 6 decades, this was a budget overrun)" \
  "src/algos/setuptest.f90|550.0_wp*trial_md_ps|the Vbias dump must NOT scale with the trial length any more (it gave exactly 1 bias update)" \
  "src/printouts.f90|(1 MTD round, no iteration)|the old onebatch line that did not say how many MTDs (must be gone)" \
  "src/algos/parallel.f90|rdensembleparam(conformerfile|the -iru read must go through irufile, not the hardcoded crest_conformers.xyz (must be gone: it would ignore --loadconformer)" \
  "src/algos/parallel.f90|write (stdout,'(3x,a,a)')    'bias ensemble   : ',conformerfile|the banner must name the file that was ACTUALLY read (must be gone)" \
  "src/dynamics/dynamics_module.f90|\" ps,  elapsed \"|the DOUBLE space before 'elapsed' in the progress line (must be gone; quoted so the old sample in the comment does not match)" \
  "src/dynamics/dynamics_module.f90|i0,\" h \",i0,\" min \",i0,\" sec\"|the spaced time string (must be gone: 17 columns for '3 h 32 min 21 sec')" \
  "src/dynamics/dynamics_module.f90|\", expected remaining time \"|the 26-column label (must be gone: it is what made the progress line wrap)" \
  "src/calculator/orca_sc.f90|if (len_trim(envstr) > 0.and..not.callprinted) then|the presence-only PRINT_CALL test (must be gone: it made =0 print the line)" \
  "src/calculator/orca_sc.f90|if (len_trim(envstr) == 0) then|the presence-only NO_OVERSUBSCRIBE test (must be gone: =0 dropped --oversubscribe)" \
  "src/parsing/parse_xtbinput.f90|              call add_hardcons(env,i1,i2,dist)|the UNGUARDED registration (must be gone: it made EVERY -cbonds/-cmetal bond rigid)" \
  "src/parsing/parse_xtbinput.f90|            call add_hardcons(env,j,i,dist)|the UNGUARDED pairwise registration (must be gone: same reason)" \
  "src/parsing/confparse2.f90|call parse_xtbinputfile(env,'bondlengths')|the unflagged bondlengths call (must be gone: without hard=.false. those bonds became rigid)" \
  "src/algos/setuptest.f90|MD%stopcheck = .true.|the trial MTD must NOT opt into the stop-file poll (io/=0 there means 'shorten dt and retry' and aborts the run)" \
  "src/algos/dynamics.f90|if (io == 0) then|the old two-branch MD status (must be gone: term=3 landed in the else and was reported as ABNORMALLY + status_failed)" \
  "src/choose_settings.f90|alpinc = 2.0d0|the old onebatch alpha divisor (must be gone: that put 0.50 in the grid, which pushed structures apart too fast)" \
  "src/dynamics/dynamics_module.f90|inquire (file=fname,exist=ex)|the CWD-only probe (must be gone: under -scratch that is the scratch dir, not where the user's stop file is)" \
  "src/algos/parallel.f90|call engrad(mols(job),calculations(job),energy,grads(:,:,job),io)|crest_sploop must NOT call engrad() any more (must be gone: it asked ORCA for a gradient nobody read)" \
  "src/calculator/subprocess_types.f90|if (atmp(1:1) .eq. '!'.and.index(atmp,'engrad') .ne. 0) then|the CASE-SENSITIVE EnGrad test (must be gone: 'EnGrad' in a template was not recognised)" \
  "src/algos/singlepoint.f90|call crest_sploop(env,nat,nall,at,xyz,eread,.true.)|the 7-argument crest_sploop call (must be gone: 6 dummies only, it survived on an ABI accident)" \
  "src/algos/dynamics.f90|allocate (grad(3,mol%nat),source=0.0_wp)|the discarded gradient array of the standalone-MD WBO call (must be gone)" \
  "src/algos/setuptest.f90|allocate(grd(3,mol%nat))|the discarded gradient array of the trial-MD WBO call (must be gone)" \
  "src/algos/numhess.f90|ofreq(:),grad0(:,:),energies0(:)|the never-read grad0 array declaration (must be gone)" \
  "src/iomod.F90|backspace (ich)|the read-to-EOF append emulation (must be gone: it re-read the target on every append)" \
  "src/dynamics/dynamics_module.f90|molo%at = mol%at|the never-read molo%at copy (must be gone)" \
  "src/dynamics/dynamics_module.f90|molo%nat = mol%nat|the never-read molo%nat copy (must be gone)" \
  "src/dynamics/dynamics_module.f90|subroutine md_calc_mtd(mol,dat,epot,grd,pr)|the old md_calc_mtd signature (must be gone: the scratch is passed in now)" \
  "src/algos/search_conformers.f90|env%calc%optlev =  1|the old tight-stage level (must be gone: it is olev_tightcs now)" \
  "src/optimize/optutils.f90|maxcycle = max(maxcycle,200)|the old per-level cycle floor (must be gone: it is what pinned every multilevel stage to 200)" \
  "src/optimize/optutils.f90|maxcycle = min(maxcycle,10000)|the old per-level cycle ceiling (must be gone: replaced by the global budget)" \
  "src/calculator/calc_type.f90|logical  :: anopt = .false.|the old 'discard structures that hit maxcycle' default (must be gone)" \
  "src/iomod.F90|'MDFILES','METADYN','NORMMD','MRMSD','OPTIM'|the old 5-entry hot list, including the three names nothing can write to (must be gone)" \
  "src/confparse.f90|env%ewin = 6.0d0|the old 6 kcal EWIN default (must be gone)" \
  "src/confparse.f90|env%rthr = 0.125d0|the old 0.125 Ang RTHR default (must be gone)" \
  "src/confparse.f90|env%ethr = 0.05d0|the old 0.05 kcal ETHR default (must be gone)" \
  "src/algos/search_conformers.f90|env%ewin        = 6.0_wp|the multilevel hardcoded EWIN reset (must be gone: it undid the tuned default)" \
  "src/algos/search_conformers.f90|env%rthr        = 0.125_wp|the multilevel hardcoded RTHR reset (must be gone: it undid the tuned default)" \
  "src/flexi.F90|integer,parameter :: NTOP = 20|the colliding size parameter name (must be gone: NTOP and ntop are one symbol in Fortran)" \
  "src/flexi.F90|ti(NTOP),tj(NTOP)|the arrays on the colliding parameter (must be gone)" \
  "src/flexi.F90|if (ntop .lt. NTOP) then|the bound on the colliding parameter (must be gone)" \
  "src/flexi.F90|asym(at(i)),i,asym(at(j)),j,r*autoaa|the table printing REDUCED indices (must be gone: it reported P at 50/51 as 47/48)" \
  "src/algos/search_conformers.f90|env%calc%optlev =  2|the absolute olev_vtight pin in stage 6 (must be gone: it is olev_vtightcs now)" \
  "src/optimize/optutils.f90|      ethr = 1.d-7  !> olev_vtightcs|stage 6 must NOT use olev_vtight's 1e-7 (must be gone: it is half of stage 5 now)" \
  "src/optimize/optutils.f90|      gthr = 2.d-4  !> olev_vtightcs|stage 6 must NOT use olev_vtight's 2e-4 (must be gone: it is half of stage 5 now)" \
  "src/optimize/optimize_module.f90|L-BFGS currently not implemented|the fatal L-BFGS stub message (must be gone: engine 1 is implemented now)" \
  "src/optimize/lbfgs.f90|subroutine packflat_shift|the placeholder helper from the first draft (must be gone: it was a no-op stub)" \
  "src/optimize/lbfgs.f90|function compute_old|the placeholder helper from the first draft (must be gone: it was a no-op stub)" \
  "src/optimize/lbfgs.f90|1+(i-1)/3*0|the dead '*' + '0' index arithmetic of the first-draft packflat (must be gone)" \
  "src/calculator/calc_type.f90|self%binary = 'gxtb'|a binary literally called 'gxtb' (must be gone: g-xTB is a flag on xtb now)" \
  "src/parsing/parse_calcdata.f90|job%binary = 'gxtb'|the TOML path pointing at a 'gxtb' binary (must be gone)" \
  "src/parsing/parse_calcdata.f90|job%other ='-grad'|the TOML path asking for the turbomole -grad hack (must be gone)" \
  "src/printouts.f90|NUMERICAL gradients which are SLOW and NOISY|the obsolete numerical-gradient warning (must be gone: the gradient is analytical now)"
do
  f="${c%%|*}"; rest="${c#*|}"; pat="${rest%%|*}"; desc="${rest#*|}"
  if grep -qF -- "$pat" "$f"; then
    printf '  FAIL         %-42s %s\n' "$f" "$desc"; fail=$((fail + 1))
  else
    printf '  OK (absent)  %-42s %s\n' "$f" "$desc"
  fi
done

# --- COUNT-based checks.  Presence alone is not enough when a snippet is
# --- legitimately repeated: after P4 the bond vector has to be refreshed in the
# --- convergence-test loop AND in both correction branches, i.e. at least 3
# --- times.  A single occurrence means somebody reverted one branch to the old
# --- Jacobi sweep - which is exactly the kind of half-applied patch that a plain
# --- "does the string exist" check waves through.
want=3
got=$(grep -cF -- "shk%dr(1:3,i) = shk%xyzt(1:3,iat)-shk%xyzt(1:3,jat)" src/dynamics/shake_module.f90)
if [ "$got" -ge "$want" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got" "src/dynamics/shake_module.f90" "P4: bond vector refreshed before every correction, not once per cycle"
else
  printf '  FAIL         %-42s %s\n' "src/dynamics/shake_module.f90" "P4: expected >= $want dr refreshes (test loop + both branches), found $got"
  fail=$((fail + 1))
fi

# --- CREGEN reject audit: all FOUR gates must write, and each exactly the way
# --- it was wired.  Counting the entry point is the only reliable assertion:
# --- a presence check on the subroutine name would still pass if three of the
# --- four calls had been deleted.
want_clash=3   # dissoc / cnorm<1e-6 / distance, three distinct reasons
got_clash=$(grep -cF -- "call exclfile_write(nat,at,xyz(:,:,j),comments(j),'CLASH'" src/cregen.f90)
if [ "$got_clash" -eq "$want_clash" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got_clash" "src/cregen.f90" "CREGEN audit: all three CLASH sub-reasons are recorded"
else
  printf '  FAIL         %-42s %s\n' "src/cregen.f90" "CREGEN audit: expected $want_clash CLASH writers, found $got_clash"
  fail=$((fail + 1))
fi
want_topo=2    # topology mismatch + E/Z isomer
got_topo=$(grep -cF -- "call exclfile_write(nat,at,c1" src/cregen.f90)
if [ "$got_topo" -eq "$want_topo" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got_topo" "src/cregen.f90" "CREGEN audit: both TOPO sub-reasons are recorded (mismatch + E/Z)"
else
  printf '  FAIL         %-42s %s\n' "src/cregen.f90" "CREGEN audit: expected $want_topo TOPO writers, found $got_topo"
  fail=$((fail + 1))
fi
for gate in EWIN DUP; do
  got=$(grep -cF -- "'$gate'" src/cregen.f90)
  if [ "$got" -eq 1 ]; then
    printf '  OK (x%s)      %-42s %s\n' "$got" "src/cregen.f90" "CREGEN audit: the $gate gate writes exactly once"
  else
    printf '  FAIL         %-42s %s\n' "src/cregen.f90" "CREGEN audit: expected 1 $gate writer, found $got"
    fail=$((fail + 1))
  fi
done
# --- the audit must be opt-in: gc enforcement of the default-off contract
if grep -qF -- "if (trim(val) == '0' .or. trim(val) == 'none' .or. trim(val) == 'false') return" src/exclfile.f90; then
  printf '  OK           %-42s %s\n' "src/exclfile.f90" "CREGEN audit: CREST_EXCL_XYZ=0/none/false really disable it"
else
  printf '  FAIL         %-42s %s\n' "src/exclfile.f90" "CREGEN audit: the explicit-off values are missing from exclfile_check()"
  fail=$((fail + 1))
fi

# --- The L-BFGS pair update must be gated TWICE: once on the relative
# --- curvature condition and once on yy>0.  The second guard exists because a
# --- numerically saturated y can slip past the first (sqrt(ss*yy) can become
# --- Inf/NaN on a pathological gradient) and 1/(y.s) would then be stored as a
# --- NaN rho, poisoning every following direction silently.  A presence check
# --- on the subroutine alone would not see the second guard go missing.
# --- NOTE: this is grep -E, so the parentheses must NOT be backslash-escaped -
# --- '\(' matches a literal paren and would find neither guard.
want_guard=2
got_guard=$(grep -cE -- "^    if \(ys \.le\. lbfgs_curv_tol|^    if \(yy \.le\. 0\.0_wp\) return" src/optimize/lbfgs.f90)
if [ "$got_guard" -eq "$want_guard" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got_guard" "src/optimize/lbfgs.f90" "L-BFGS: both curvature gates are present (relative s.y test + yy>0)"
else
  printf '  FAIL         %-42s %s\n' "src/optimize/lbfgs.f90" "L-BFGS: expected $want_guard curvature guards, found $got_guard"
  fail=$((fail + 1))
fi

# --- (2026-09-30) -gxtb and -gxtb_dev must BOTH select the method.  A presence
# --- check on the warn call alone would still pass with -gxtb left as a dead
# --- end (warn + abort), which is the upstream behaviour we deliberately
# --- changed.  Count the assignments instead: exactly two.
want_gfnver_gxtb=2
got_gfnver_gxtb=$(awk '{ l=$0; sub(/^[ \t]+/,"",l); if (l ~ /^!/) next;
                         if (l ~ /^env%gfnver = .gxtb_dev./) n++ }
                       END{ print n+0 }' src/confparse.f90)
if [ "$got_gfnver_gxtb" -eq "$want_gfnver_gxtb" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got_gfnver_gxtb" "src/confparse.f90" "-gxtb AND -gxtb_dev both SELECT gxtb_dev (neither is a dead end)"
else
  printf '  FAIL         %-42s %s\n' "src/confparse.f90" "-gxtb: expected $want_gfnver_gxtb assignments to env%gfnver, found $got_gfnver_gxtb"
  fail=$((fail + 1))
fi

# --- interface blocks for flexi(): the count must EQUAL the number of files that
# --- call it.  A presence check would still pass if the minitools.f90 copy were
# --- dropped, which is exactly the half-application that breaks the build on the
# --- remote (two files, three call sites, all of them need it in scope).
want_flexif=2
got_flexif=$(grep -lF -- "subroutine flexi(mol,rednat,includeAtom,flex,nfix,fixat)" \
              src/choose_settings.f90 src/minitools.f90 2>/dev/null | wc -l)
if [ "$got_flexif" -eq "$want_flexif" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got_flexif" "src/{choose_settings,minitools}.f90" "flexi(): an explicit interface is in scope in BOTH callers (optional dummies require one)"
else
  printf '  FAIL         %-42s %s\n' "src/{choose_settings,minitools}.f90" "flexi(): expected an interface block in $want_flexif callers, found $got_flexif"
  fail=$((fail + 1))
fi

# --- parallel_md_block_printout() got an OPTIONAL third dummy (ntot) on
# --- 2026-09-29, so BOTH callers in src/algos/parallel.f90 now need an explicit
# --- interface as well.  Same failure mode as flexi, hence the same
# --- file-count shape: one interface per caller, the definition carries none.
# --- The two call sites live in DIFFERENT subroutines (crest_search_multimd and
# --- crest_search_multimd2), so two separate blocks are required - a single
# --- host-scoped one would not be visible in the other routine.
want_mdprintf=2
got_mdprintf=$(grep -cF -- "subroutine parallel_md_block_printout(MD,vz,ntot)" \
               src/algos/parallel.f90)
# 3 occurrences expected: 2 interface blocks + the definition itself.
want_mdprintf_total=3
if [ "$got_mdprintf" -eq "$want_mdprintf_total" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$want_mdprintf" "src/algos/parallel.f90" "parallel_md_block_printout(): an explicit interface is in scope in BOTH callers (optional ntot)"
else
  printf '  FAIL         %-42s %s\n' "src/algos/parallel.f90" "parallel_md_block_printout(): expected $want_mdprintf_total declarations (2 interfaces + definition), found $got_mdprintf"
  fail=$((fail + 1))
fi

# --- The optional dummy itself must still be optional, otherwise the 2-argument
# --- and 3-argument call shapes cannot both be legal and the brace of interface
# --- blocks above would be pointless work.
if grep -qF -- "integer,intent(in),optional :: ntot" src/algos/parallel.f90; then
  printf '  OK           %-42s %s\n' "src/algos/parallel.f90" "parallel_md_block_printout(): ntot is still OPTIONAL (both call shapes stay legal)"
else
  printf '  FAIL         %-42s %s\n' "src/algos/parallel.f90" "parallel_md_block_printout(): ntot lost its optional attribute"
  fail=$((fail + 1))
fi

# --- EVERY routine that calls exclfile_*() must `use exclfile_module` IN THAT
# --- ROUTINE.  This is the link-time trap: without the use, the compiler
# --- silently invents an implicit external and the object asks for the
# --- un-mangled `exclfile_write_`, while exclfile.f90 exports
# --- `__exclfile_module_MOD_exclfile_write`.  Both files compile clean; only
# --- the LINK fails, with an error that names cregen.f90 but nothing useful
# --- about WHY.  discardbroken() was missed on 2026-09-29 for exactly this
# --- reason (its calls sit outside the routines the audit patch touched).
# ---
# --- Counting `use exclfile_module` occurrences against the number of routines
# --- that call exclfile_*() is the cheap invariant: 5 callers, 5 use lines.
# --- The count needs awk, not grep: the !>--- comment blocks above the repaired
# --- routine contain prose lines that also start with the word "use", and a
# --- grep pattern broad enough to catch a real use statement (which carries a
# --- trailing "!> LOCAL EXTENSION" comment) also catches that prose.  Real use
# --- statements are statements: strip leading blanks, then require the line to
# --- begin with "use exclfile_module" AND reject anything starting with '!'.
want_exuses=5
got_exuses=$(awk '{ l=$0; sub(/^[ \t]+/,"",l); if (l ~ /^!/) next;
                    if (l ~ /^use[ \t]+exclfile_module([ \t]|$|,)/) n++ }
                  END{ print n+0 }' src/cregen.f90)
if [ "$got_exuses" -eq "$want_exuses" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got_exuses" "src/cregen.f90" "every exclfile-calling routine has 'use exclfile_module' in scope (link-time trap)"
else
  printf '  FAIL         %-42s %s\n' "src/cregen.f90" "expected $want_exuses 'use exclfile_module' (one per calling routine), found $got_exuses"
  fail=$((fail + 1))
fi

# --- And the specific one that broke: discardbroken()'s `use` must sit inside
# --- discardbroken itself, not merely somewhere earlier in the file.  Anchoring
# --- on the line number would be brittle, so instead require the use line to
# --- appear between the routine header and its `implicit none`.
if awk '/^ *subroutine +discardbroken\(/{f=1} f&&/^ *use +exclfile_module/{hit=1} f&&/^ *implicit +none/{print (hit?"YES":"NO"); exit}' src/cregen.f90 | grep -q YES; then
  printf '  OK           %-42s %s\n' "src/cregen.f90" "discardbroken() itself uses exclfile_module (the routine that failed to link)"
else
  printf '  FAIL         %-42s %s\n' "src/cregen.f90" "discardbroken() is missing its OWN 'use exclfile_module' (undefined reference to exclfile_write_)"
  fail=$((fail + 1))
fi

# --- The dynamic-sharing calls must be present in BOTH MD drivers
# --- (crest_search_multimd and crest_search_multimd2).  Wiring only one of them
# --- is the obvious half-application: -nci / -entropy would keep idling.
for snip in "call crest_dynshare_init(nsim,T,env%Threads,Tn)"; do
  want=2
  got=$(grep -cF -- "$snip" src/algos/parallel.f90)
  if [ "$got" -eq "$want" ]; then
    printf '  OK (x%s)      %-42s %s\n' "$got" "src/algos/parallel.f90" "dynamic sharing (MD): '$snip' wired into both MD drivers"
  else
    printf '  FAIL         %-42s %s\n' "src/algos/parallel.f90" "dynamic sharing (MD): expected $want x '$snip' (one per MD driver), found $got"
    fail=$((fail + 1))
  fi
done
# --- ... and in BOTH batch loops (crest_oloop = the ensemble/multilevel
# --- optimization, crest_sploop = the refine singlepoint re-ranking).  Those
# --- are the stages that idle hardest at the end: the queue is drained in waves
# --- of Threads/CREST_CORES_PER_JOB jobs, so the last 1-2 structures would run
# --- on 4 cores while the rest of the machine waits.
for snip in "call crest_dynshare_init(nall,T,env%Threads,Tn)"; do
  want=2
  got=$(grep -cF -- "$snip" src/algos/parallel.f90)
  if [ "$got" -eq "$want" ]; then
    printf '  OK (x%s)      %-42s %s\n' "$got" "src/algos/parallel.f90" "dynamic sharing (batch opt): '$snip' wired into crest_oloop AND crest_sploop"
  else
    printf '  FAIL         %-42s %s\n' "src/algos/parallel.f90" "dynamic sharing (batch opt): expected $want x '$snip' (oloop + sploop), found $got"
    fail=$((fail + 1))
  fi
done
# --- release/reset: 2 MD drivers + crest_oloop + crest_sploop = 4
for snip in "call crest_dynshare_release()" "call crest_dynshare_reset()"; do
  want=4
  got=$(grep -cF -- "$snip" src/algos/parallel.f90)
  if [ "$got" -eq "$want" ]; then
    printf '  OK (x%s)      %-42s %s\n' "$got" "src/algos/parallel.f90" "dynamic sharing: '$snip' in both MD drivers + both batch loops"
  else
    printf '  FAIL         %-42s %s\n' "src/algos/parallel.f90" "dynamic sharing: expected $want x '$snip' (2 MD + oloop + sploop), found $got"
    fail=$((fail + 1))
  fi
done

# --- (2026-09-27) energy-only calls: every converted site must still call
# --- singlepoint() (not engrad()).  A revert of any one of them shows up as a
# --- drop in the count even though the string is still present elsewhere.
for spec in \
  "src/algos/parallel.f90|call singlepoint(|2|crest_sploop + the WBO warm-up" \
  "src/algos/dynamics.f90|call singlepoint(|1|standalone MD WBO" \
  "src/algos/setuptest.f90|call singlepoint(|1|trial MD WBO" \
  "src/algos/numhess.f90|call singlepoint(|1|numhess initial single point"
do
  f="${spec%%|*}"; rest="${spec#*|}"; pat="${rest%%|*}"; rest="${rest#*|}"
  want="${rest%%|*}"; desc="${rest#*|}"
  got=$(grep -cF -- "$pat" "$f")
  if [ "$got" -eq "$want" ]; then
    printf '  OK (x%s)      %-42s %s\n' "$got" "$f" "'$pat' present exactly $want time(s) - $desc"
  else
    printf '  FAIL         %-42s %s\n' "$f" "expected $want x '$pat' ($desc), found $got"
    fail=$((fail + 1))
  fi
done
# --- the restart-snapshot switch must gate all three writers + its own definition
want=4
got=$(grep -cF -- "restart_save_enabled()" src/restartlog.f90)
if [ "$got" -eq "$want" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got" "src/restartlog.f90" "CREST_RESTART_SAVE gates trackrestart + both trackensemble variants"
else
  printf '  FAIL         %-42s %s\n' "src/restartlog.f90" "expected $want uses of restart_save_enabled() (def + 3 gates), found $got"
  fail=$((fail + 1))
fi
# --- the cheap counter must be used at both multilevel sites
want=2
got=$(grep -cF -- "call count_xyz_structures(trim(inpnam),nat,nall)" src/algos/search_conformers.f90)
if [ "$got" -eq "$want" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got" "src/algos/search_conformers.f90" "both multilevel size-counts go through the cheap counter"
else
  printf '  FAIL         %-42s %s\n' "src/algos/search_conformers.f90" "expected $want cheap counter calls, found $got"
  fail=$((fail + 1))
fi
# --- the merged per-worker trajectory must be deleted in BOTH collect routines
want=2
got=$(grep -cF -- "call remove(atmp)" src/algos/parallel.f90)
if [ "$got" -eq "$want" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got" "src/algos/parallel.f90" "both collect() routines delete the merged per-worker trajectory"
else
  printf '  FAIL         %-42s %s\n' "src/algos/parallel.f90" "expected $want 'call remove(atmp)' (one per collect), found $got"
  fail=$((fail + 1))
fi
# --- the MTD scratch is allocated exactly once, in dynamics()
want=1
got=$(grep -cF -- "allocate (grdmtd(3,mol%nat),source=0.0_wp)" src/dynamics/dynamics_module.f90)
if [ "$got" -eq "$want" ]; then
  printf '  OK (x%s)      %-42s %s\n' "$got" "src/dynamics/dynamics_module.f90" "the MTD bias scratch is allocated once per MD run, not per step"
else
  printf '  FAIL         %-42s %s\n' "src/dynamics/dynamics_module.f90" "expected $want allocation of grdmtd, found $got"
  fail=$((fail + 1))
fi

# ---------------------------------------------------------------------------
#  Errors that only gfortran sees.  fparser checks syntax, not types, and this
#  tree is cross-compiled - there is no Fortran compiler here - so without this
#  step a single type error costs a full remote rebuild cycle.
echo
echo "=== C. Fortran semantics: allocate(...) source=/mold= type check ======"
SEMCHECK="../_orig/f90_semantic_check.py"
PYCMD=""
for cand in py python3 python; do
  if command -v "$cand" >/dev/null 2>&1; then PYCMD="$cand"; break; fi
done
if [ -n "$PYCMD" ] && [ -f "$SEMCHECK" ]; then
  # NOTE: pass the DIRECTORY, never a file list.  The derived-type table has to
  # know every type in the tree before any single file is audited, so the
  # checker needs one invocation covering all sources.
  semout=$("$PYCMD" "$SEMCHECK" src 2>&1); semrc=$?
  printf '  %s\n' "$(printf '%s\n' "$semout" | tail -1)"
  if [ "$semrc" -ne 0 ]; then
    printf '%s\n' "$semout" | sed 's/^/  /'
    fail=$((fail + 1))
  fi
else
  echo "  SKIP (no python found, or $SEMCHECK missing)"
fi

# ---------------------------------------------------------------------------
# --- Case-insensitive NAME COLLISIONS in one scoping unit.
# --- Fortran treats NTOP and ntop as the SAME symbol, so
# ---     integer,parameter :: NTOP = 20
# ---     integer :: ...,ntop,...
# --- is a hard error on the Intel compiler ("Symbol 'ntop' already has basic
# --- type of INTEGER") while gfortran 16.2 SILENTLY ACCEPTS it.  That is how it
# --- reached a remote build: neither the syntax nor the semantic checker can
# --- see it.  Cheap standalone check, no .mod files needed.
echo
echo "=== C2. Fortran name collisions (case-insensitive, per scoping unit) ==="
COLLCHECK="../_orig/f90_name_collision_check.py"
if [ -n "$PYCMD" ] && [ -f "$COLLCHECK" ]; then
  colout=$("$PYCMD" "$COLLCHECK" src 2>&1); colrc=$?
  printf '  %s\n' "$(printf '%s\n' "$colout" | tail -1)"
  if [ "$colrc" -ne 0 ]; then
    printf '%s\n' "$colout" | sed 's/^/  /'
    fail=$((fail + 1))
  fi
else
  echo "  SKIP (no python found, or $COLLCHECK missing)"
fi

# ---------------------------------------------------------------------------
# --- Fortran SYNTAX, whole tree.
# --- WHY THIS SECTION EXISTS (2026-10-04): f90_syntax_check.py takes its input
# --- as argv[1:] - a FILE LIST.  Invoked WITHOUT arguments it checks nothing
# --- and still prints "RESULT: 0 failure(s)", so "ran it, it was green" meant
# --- absolutely nothing and the check was in practice never performed.  It is
# --- the only check here that sees a plain typo like "env.legacy" instead of
# --- "env%legacy" (parallel.f90:875), which fparser rejects and which would
# --- otherwise have cost a full remote rebuild cycle.
# --- NOTE: this costs ~2 min - fparser parses all 166 sources.  Drop the
# --- section if that is too slow, but do NOT drop the habit of passing files.
echo
echo "=== C3. Fortran syntax (fparser, all sources) ========================="
SYNCHECK="../_orig/f90_syntax_check.py"
# Known fparser limitations on VALID Fortran that is in this tree (all three
# pre-date 2026-10-04, all three compile with gfortran 16 and ifx):
#   calc_type.f90 / constraints.f90 - "procedure,private :: a, b, &" continued
#                                     over two lines (F2003 type-bound procs)
#   lj.f90                          - a continued expression that starts with
#                                     an operator
SYNFP="calc_type.f90|constraints.f90|lj.f90"
if [ -n "$PYCMD" ] && [ -f "$SYNCHECK" ]; then
  synout=$("$PYCMD" "$SYNCHECK" $(find src -name '*.f90' -o -name '*.F90' | sort | tr '\n' ' ') 2>&1)
  synall=$(printf '%s\n' "$synout" | grep -cE "^  FAIL" || true)
  synreal=$(printf '%s\n' "$synout" | grep -E "^  FAIL" | grep -Ev "$SYNFP" | wc -l)
  echo "  $synall failure(s), $((synall - synreal)) of them in the known-FP list"
  if [ "$synreal" -ne 0 ]; then
    # print only the FAIL blocks whose file is NOT whitelisted (awk keeps the
    # ">>>" detail line together with the file it belongs to)
    printf '%s\n' "$synout" | awk -v fp="$SYNFP" '
      /^  FAIL/ { show = ($0 !~ fp) ? 1 : 0 }
      /^  OK/   { show = 0 }
      show      { print "  " $0 }' | head -40
    fail=$((fail + 1))
  fi
else
  echo "  SKIP (no python found, or $SYNCHECK missing)"
fi

# ---------------------------------------------------------------------------
echo
echo "=== D. dependencies pre-placed (offline build possible) =============="
for d in tblite toml-f gfn0 gfnff lwoniom pvol test-drive; do
  if [ -f "subprojects/$d/CMakeLists.txt" ]; then
    printf '  OK    subprojects/%-11s (%s files)\n' "$d" "$(find "subprojects/$d" -type f 2>/dev/null | wc -l)"
  else
    printf '  MISS  subprojects/%-11s\n' "$d"; fail=$((fail + 1))
  fi
done

echo
if [ "$fail" -eq 0 ]; then
  echo "RESULT: PASS - this tree carries all local fixes and builds offline."
  exit 0
else
  echo "RESULT: FAIL - $fail problem(s). Do NOT build this tree."
  exit 1
fi
