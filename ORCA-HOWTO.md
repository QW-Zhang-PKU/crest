# Using ORCA as the CREST calculator

CREST 3.0.2 **already ships an ORCA subprocess interface** — no new calculator code
is required.  What you configure is:

| piece | file | key |
|---|---|---|
| the driver             | `src/calculator/orca_sc.f90`   | (built in) |
| dispatch               | `src/calculator/calculator.F90` | `case (jobtype%orca)` |
| input template reader  | `src/calculator/subprocess_types.f90` | `type(orca_input)` |
| TOML keys              | `src/parsing/parse_calcdata.f90` | `method` / `orca_cmd` / `orca_template` |

The only functional gap that existed was `rd_ORCA_wbo()`, which used to abort with
`**ERROR** ORCA bond order reader not yet implemented`.  It is implemented now
(see *Bond orders* below).

---

## 1. The three files you need

    your-molecule.xyz      the coordinates
    orca.inp               the calculation level / basis set   (template)
    crest-orca.toml        what ties them together

Both `orca.inp` and `your-system-orca.toml` in this directory are ready-to-use
examples — copy them next to your coordinates and edit three lines.

---

## 2. Recommended setup: ORCA for the refinement, cheap method for the MTD

A conformational search evaluates the energy **and gradient at every MD step** of
every MTD.  That is thousands of single points — a full ORCA run for each one is
usually not affordable.

CREST therefore supports *refinement levels*: the MTD runs on a cheap level and
ORCA is only applied to the final ensemble.  Add **two** calculation levels:

```toml
input = "your-molecule.xyz"

# ---- level 1: cheap, drives the (thousands of) MTD / pre-optimisation steps
[[calculation.level]]
method = "gff"                       # or "gfn2" / "gfn" for external xtb

# ---- level 2: ORCA, used ONLY in the refinement stage
[[calculation.level]]
method        = "orca"
orca_cmd      = "/PATH/TO/orca"
orca_template = "orca.inp"
chrg          = 1
multiplicity  = 1                    # 2S+1  <-- see the warning in section 4
refine        = "opt"                # "opt" = re-optimise the whole ensemble
                                     # "sp"  = single-point re-ranking only
```

The `refine` key is collected into CREST's refinement queue and executed *after*
the ensemble optimisation.  Valid values:

| `refine =` | CREST value | what ORCA then does |
|---|---|---|
| `"sp"`, `"singlepoint"` | 1 | one single point per structure (cheapest) |
| `"add"`, `"correction"` | 2 | one single point per structure, added to the energies |
| `"opt"`, `"optimization"` | 3 | full geometry optimisation of every structure |
| (omitted) | 0 | this level is not used in the refinement |

In the log you will see `> Geometry optimization of N structures` /
`> Singlepoint re-ranking for N structures` when the ORCA level actually fires.

## 3. Simple setup: ORCA drives everything

If you really want ORCA at every MD step (small molecule, cheap method), just make
ORCA the *only* level:

```toml
input = "your-molecule.xyz"

[calculation.level]
method        = "orca"
orca_cmd      = "/PATH/TO/orca"
orca_template = "orca.inp"
chrg          = 1
multiplicity  = 1
```

Then run e.g.

```bash
crest your-molecule-orca.toml -T 40 --quick --noreftopo
```

**Cost warning:** `-T 40` with 6 MTDs gives 6 MTDs x ~6 threads, but each of them
then performs a *serial* chain of ORCA single points.  Expect this to be orders of
magnitude slower than xtb — measure with the trial MTD (0.25 ps) first.

---

## 4. Charge and multiplicity  (read this!)

CREST passes the charge and the multiplicity to ORCA by writing

```
*xyz <charge> <multiplicity>
```

internally, for every single point.  The values come from the **CREST input**, not
from `orca.inp` (any coordinate block you put in the template is stripped).

* `chrg` behaves as expected.
* For the spin, use the key **`multiplicity = 2S+1`**
  (1 = singlet, 2 = doublet, 3 = triplet, ...).

> **Do not use `uhf =` for an open-shell ORCA calculation.**  CREST stores the
> `uhf`/`multiplicity` value in one field and only clamps it to `>= 1` before
> writing it as the ORCA multiplicity.  `multiplicity = 2` therefore gives a
> doublet (correct), whereas `uhf = 1` also becomes multiplicity 1 — i.e. a
> *singlet*, which is wrong.  `uhf` follows the xtb convention (n(alpha)-n(beta))
> and is only safe at its default of 0.

---

## 5. Rules for `orca.inp`

CREST reads the template once, then appends the live geometry before each call.
The reader **rejects** (or silently strips) the following, so don't use them:

| not allowed in `orca.inp` | why |
|---|---|
| a coordinate block (`*xyz ...`, `*xyzfile ...`) | CREST supplies the coordinates |
| a `%coords ... end` block | same |
| `$new_job` | only single jobs are supported |
| `! Opt`, `! Freq` (`! MD` too) | CREST drives the runtype |

If no `! ... EnGrad` line is present, ` EnGrad` is appended to the first `!` line
automatically.  Keeping it explicit (`! r2SCAN-3c EnGrad`) is clearer.

---

## 6. Parallelisation

ORCA parallelises with **MPI**, configured inside the input via

```
%pal
  nprocs N
end
```

It does **not** read `OMP_NUM_THREADS` — and that is the variable CREST uses to
hand a core budget to every external program.  A plain ORCA therefore runs on
**one core no matter what `-T` says**.  This is the most common reason an
otherwise healthy CREST+ORCA setup looks inexplicably slow.

### Fixed here: the block is written for you

`write_orca_input()` reads the budget CREST prepared for the *current* call and
writes the matching `%pal` block.  **Leave `%pal` commented out in `orca.inp`**
— CREST fills it in per call, per worker.

`new_ompautoset()` (ompmklset.F90) already computes a fair share per job and
exports it as `OMP_NUM_THREADS`; that is the number that becomes `nprocs`:

| CREST stage | modus | written as `nprocs` |
|---|---|---|
| initial optimisation (`trialOPT_calculator`), singlepoints | `max` | **T** — every core |
| trial MTD (`trialMD_calculator`) | `share` | **T / nmetadyn** — same as a production MTD |
| MTD / MD (`parallel.f90`) | `auto`, `nsim` | **T / nsim** |
| ensemble / batch optimisation + batch singlepoints (`crest_oloop`, `crest_sploop`) | `fixed` | **CREST_CORES_PER_JOB**, default **4** |
| refinement (`refine.f90`) | `subprocess` | **T** |
| `setuptest` sanity call | `min` | 1 |

So `-T 40` with 6 MTDs gives `%pal nprocs 40` for the initial optimisation and
`%pal nprocs 6` inside each MTD worker — the whole budget is used, and it is
divided evenly once parallelism actually starts.  Nothing to configure.

### Batch stages: fixed 4 cores per job, the rest queue

The batch loops (`crest_oloop`, `crest_sploop`) do **not** use the `T/nall`
scheme.  That scheme degenerates badly for a real conformer ensemble: as soon as
`nall > T` every job is handed exactly **one** core, i.e. hundreds of serial
ORCA gradients.

They use a `fixed` distribution instead:

```
cores per job  = CREST_CORES_PER_JOB          (default 4)
jobs at once   = T / cores_per_job            (integer division)
```

With `-T 40` that is **10 ORCA jobs × 4 cores = 40 cores**, and the remaining
structures wait in the OpenMP task queue — which is exactly the "N in parallel,
the rest queue" behaviour you asked for.  No oversubscription, no thrashing.

When the ensemble is *smaller* than `T / cores_per_job`, the leftover cores go to
the jobs that do exist, so nothing is wasted:

| `-T` | structures | jobs × cores | `nprocs` |
|---|---|---|---|
| 40 | 500 | 10 × 4 = 40 | 4 |
| 40 | 50 | 10 × 4 = 40 | 4 |
| 40 | 5 | 5 × 8 = 40 | 8 |
| 40 | 2 | 2 × 20 = 40 | 20 |
| 16 | 300 | 4 × 4 = 16 | 4 |

Change the per-job core count **without recompiling**:

```bash
CREST_CORES_PER_JOB=8 crest your-system-orca.toml -T 40 --quick --noreftopo
```

(`CREST_ORCA_NPROCS` is a different knob: it overrides only the number written
into `%pal`, without changing how many jobs run at once.  Prefer
`CREST_CORES_PER_JOB` so both stay consistent.)

`--oversubscribe` is appended (once) whenever a `%pal` block is written or the
template already carries one, so OpenMPI does not refuse to start more ranks
than it detected slots.  If you put `oversubscribe` into `orca_cmd` yourself,
CREST sees it and does not add a second one.

### Trial MTD now uses real cores, and the estimate is honest

The trial MTD exists only to calibrate the MD time step and to project the
runtime — it produces no science.  Upstream nevertheless ran it with
`new_ompautoset(env,'min',...)` (comment: *"should be set to 1 to simulate max
parallelization"*), i.e. a fully serial ORCA gradient at every MD step.

Two things were wrong with that:

1. the trial itself was needlessly slow, and
2. `timeestimate()` then **assumed** a perfect `floor(T/nmetadyn)`× linear
   speed-up to convert that serial measurement into a batch prediction — the
   origin of wildly optimistic runtime forecasts.

Now the trial gets `'share'`: exactly the budget one production MTD gets,
`T/nmetadyn` cores (`%pal nprocs 6` for `-T 40` with 6 MTDs).  The measurement is
therefore already a per-job runtime and needs no fictional divisor:

```
Estimated runtime for one MTD (     25.0 ps) on 6 cores: 1 h 12 min 40 sec
Estimated runtime for a batch of 6 MTDs (all parallel) on 40 threads: 1 h 12 min 40 sec
```

When more MTDs are requested than fit at once the estimate accounts for the
waves (`4 at a time, 3 waves`), instead of pretending they are free.

### The trial is also 4× shorter now (0.25 ps, was 1.0 ps)

Same reasoning: it calibrates the time step / SHAKE settings and feeds the
estimator, nothing else.  Two things had to be changed together with the length:

1. **`MTD%cvdump_fs` scales with it** (`550 fs * trial_md_ps`).  The old value was
   hard-coded for a 1 ps trial; keeping it at 550 fs would have given a 0.25 ps
   trial *zero* Vbias (RMSD-CV) updates — the trial would silently degrade into
   a plain MD and stop testing the metadynamics setup.
2. **`timeestimate()` divides out the trial length** (`mdtime/trialps`).  Upstream
   implicitly assumed the measurement covered exactly 1 ps, so a shortened trial
   would have made the forecast 4× too optimistic.

The printed line now shows the real length, and the estimator output is
unaffected by the change:

```
Trial MTD 1 runtime (0.25 ps) :   2 min 30 sec
Estimated runtime for one MTD (     25.0 ps) on 6 cores: 4 h 10 min 0 sec
```

`MDSTART%dumpstep` stays at 20 fs — the trial trajectory is deleted by
`cleanup.f90` and is never read, so the number of dumped frames does not matter.

### What the start of the log looks like now

Both the initial optimisation and the trial MTD used to run with **no printout
at all** (`pr=.false.`), so with a slow external calculator the job looked hung
and there was no number to compare against afterwards.  This tree prints:

```
==================================================
           Initial Geometry Optimization
==================================================
Initial geometry optimization: 37 cycles, 12 sec   (crestopt.log)
Geometry successfully optimized.

==================================================
        Starting trial MTD to test settings
==================================================
Trial MTD settings:
   length             :   0.250 ps
   time step          :    4.00 fs
   steps              : 62
   temperature        :   298.2 K
   SHAKE              : on (mode 2)
   Vbias k / Eh       :   0.108
   Vbias alpha        :   0.5000
   Vbias dump         :   137.5 fs
   cores per job      : 8

 Trial MTD  25% (  0.06 /  0.25 ps, elapsed 12s, remaining 36s)
 Trial MTD  50% (  0.12 /  0.25 ps, elapsed 24s, remaining 24s)
 Trial MTD  75% (  0.19 /  0.25 ps, elapsed 36s, remaining 12s)
 Trial MTD 100% (  0.25 /  0.25 ps, elapsed 48s, remaining 0s)
Trial MTD 1 runtime (0.25 ps) :   48 sec
Estimated runtime for one MTD (     25.0 ps) on 8 cores: 1 h 20 min 0 sec
```

Notes:

- **37 cycles** is the frame count of `crestopt.log`.  Every optimisation engine
  (ANCOPT, RFO, gradient descent) appends exactly one frame per gradient
  evaluation, so this is the number of optimisation cycles — the same number
  ANCOPT prints as `GEOMETRY OPTIMIZATION CONVERGED AFTER n ITERATIONS`, which
  is invisible here because `pr=.false.`.
- **steps = length / time step** = 0.25 ps / 4 fs = 62.  If the trial has to
  restart (it reduces the time step up to 6 times), the later attempts print
  `Reducing the time step to ... fs and trying again ...` and the step count
  shown above belongs to the *first* attempt only.
- **Vbias dump** is `550 fs * trial_md_ps` = 137.5 fs.  It scales with the trial
  length on purpose — see above.
- The progress lines carry a `flush`, otherwise they would sit in the buffer and
  only appear once the job ends, which would defeat the point.  They are emitted
  **independently of `pr`** and only for the trial MTD; the production MTDs are
  untouched.
- `ETA` is `elapsed * (100 - %) / %`, i.e. it assumes the remaining steps cost
  the same as the ones already done.  It is a rough guide, not a forecast.

Do not confuse this with the **production** Vbias interval: that one is
`env%mddump` (`confparse.f90`), read into `mtd%cvdump_fs` in `parallel.f90`, and it
is a separate knob with its own default.  The two were never equal upstream either
(trial 550 fs vs production 1000 fs); production was halved to **500 fs (0.5 ps)**
in this tree.  `-vbdump <ps>` overrides it, and wants **ps**, not fs.

### Overriding

| variable | effect |
|---|---|
| `CREST_ORCA_NPROCS=8` | cap the ranks per ORCA call at 8, regardless of budget. `1` = no `%pal` at all |
| `CREST_ORCA_OMP_THREADS=8` | threads for the ORCA child (default **1**).  Combine with `CREST_ORCA_NPROCS=1` for the threads-instead-of-ranks mode |
| `CREST_ORCA_NO_OVERSUBSCRIBE=1` | do not append `--oversubscribe` (that flag disables MPI process binding). `0`/unset = keep it |
| `CREST_ORCA_PRINT_CALL=1` | echo the exact ORCA command line once per process. `0`/unset = silent (default) |
| `CREST_ORCA_MPI_ARGS="--bind-to core"` | append anything else to the ORCA command line |
| `CREST_TRIAL_MD_PS=0.5` | change the trial MTD length without recompiling (default `0.25`) |
| `CREST_WORK_DIR=crest_work` | subdirectory for the *cold* intermediate files (`coord`, `wbo`, `crestopt.xyz`, `gfnff_*`, `.CHRG`, `.UHF`). `none`/`0` = leave them in the workspace (upstream behaviour) |
| `CREST_STAGE_FILES=0` | disable all intermediate-file staging (`CREST_SCRATCH_DIR` for `calculation.level.N` is *not* affected — that one is independent) |
| `CREST_STOP_CHECK_SEC=300` | how often (seconds) a running dynamics looks for the stop file. `0`/negative = off |
| `CREST_STOP_FILE=stop` | name of that file. Looked up in the working directory **and** in the directory crest was started from |
| unset | use the budget / length above |

#### The `stop` file — cutting a multi-day run short without losing it

With ORCA a 6 ps MTD can run for days, and killing the job loses everything,
because no structure has been optimized or written yet.  Create the file while
crest is running:

```bash
touch stop          # in the directory crest is running in
```

At the next check (default every 300 s) the sampling ends, the trajectory that
exists so far is kept, and the run goes straight on to the optimization /
refinement stage — so you get a smaller ensemble in minutes instead of nothing.
The completion line says `*MTD   2 stopped ("stop" file)`.

Three things to know:

* **A file that was already there when crest started is ignored.**  `stop` has
  to *appear* during the run.  Without that rule a leftover file from last
  week's run would kill today's run after five minutes with no explanation.
  (The report tells you to remove it.)
* **`stop` next to your toml/xyz works in every mode.**  The file is looked up
  in the working directory **and** in the directory crest was started from, so
  it does not matter which of the two it is:
  * plain run (and with `CREST_SCRATCH_DIR`) — the two are the same directory;
  * `-scratch` — crest `chdir`s into the scratch directory
    (`scratch.f90:54`), so the launch directory is probed *in addition*.  This
    is why a relative `CREST_STOP_FILE` no longer needs an absolute path.
  An **absolute** name (`CREST_STOP_FILE=/scratch/$USER/liu/stop`) is used as
  given and is never prefixed with the launch directory.
* **Only the dynamics is affected.**  The multilevel optimization and cregen
  stages do not check it; they also do not need to, because by then structures
  are being written continuously.

All of these are read **by value**, not by mere presence: `FOO=0` means *off*,
it does not mean "set".  (An earlier version tested only whether the variable
was non-empty, so `CREST_ORCA_PRINT_CALL=0` printed the line it was meant to
suppress and `CREST_ORCA_NO_OVERSUBSCRIBE=0` dropped `--oversubscribe`.  Only
integers are honoured; a non-numeric value falls back to the default.)

```bash
CREST_ORCA_NPROCS=8 crest your-system-orca.toml -T 40 --quick --noreftopo
CREST_TRIAL_MD_PS=0.5 crest your-system-orca.toml -T 40 --quick   # a longer trial
```

Worth trying when ORCA's MPI scaling flattens before your core count does — for
a single EnGrad on a small/medium molecule 8–16 ranks is often faster than 40,
because the extra ranks spend their time communicating.

### Two things to watch

1. **Memory.** `%maxcore` is *per MPI process*.  With `nprocs 40` and
   `%maxcore 2000` ORCA can claim 40 × 2 GB = 80 GB.  Size `%maxcore` for the
   **initial optimisation**, i.e. roughly `RAM / T`.
2. **The child's thread layer is pinned to 1.**  The systemcall is built as
   `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 <orca> ORCA.in`,
   because otherwise N ranks × N threads means N² runnable threads per job.

### Verifying

```bash
grep -i "Number of parallel processes" calculation.level.1/ORCA.out
head -20 calculation.level.1/ORCA.in          # should show the %pal block
```

During the run, `htop`/`top` should show ~T busy cores in the initial
optimisation and ~T busy cores spread over `nsim` ORCA process trees during MTD.

### If you insist on `%pal` in the template

An explicit `%pal` always wins: `read_orca_input()` sets `self%mpi` when it sees
one, and the automatic block is then skipped.  In that case drive the launcher
from `orca_cmd`, e.g.

```toml
orca_cmd = "mpirun --oversubscribe -np 8 orca"
```

---

## 7. Constraining part of the structure

Short answer: **exactly the same as with xtb** — no ORCA-specific syntax, no
change to `orca.inp`.

### Why it is identical

Constraints are applied by **CREST**, not by the QM program.  In
`calculator.F90` the generic `engrad()` first calls the calculator
(`potential_core`) and then, for *any* calculator, adds the constraint bias on
top:

```fortran
if (calc%nconstraints > 0) then
  do i = 1,calc%nconstraints
    call calc_constraint(mol%nat,mol%xyz,calc%cons(i),efix,calc%grdfix)
    energy   = energy   + efix
    gradient = gradient + calc%grdfix
  end do
end if
```

`xtb_sc.f90` contains no constraint code at all, and neither does `orca_sc.f90`.
Consequences:

- ORCA always does a plain EnGrad — **constraints cost nothing in ORCA**, which
  matters when a single gradient takes minutes;
- they are enforced in the initial optimisation, the ensemble optimisation *and*
  the MTD/MD (`dynamics_module.f90` calls the same `engrad()`);
- `type(constraint)` holds only integer/real allocatables, so it survives the
  per-worker copies of `calculation_settings` (unlike the ORCA template — that
  one is why the PR89200 fix was needed).

### How to pass it

Command line `-cinp <file>` (aliases `-cfile`, `-c`, `-C`); in the TOML
`constraints = "<file>"` (also accepted: `xtbconstraints`, `xtbinput`).
The file uses xtb-xcontrol syntax, parsed by CREST itself
(`parse_xtbinput.f90`); `#` starts a comment.  A ready-to-adapt example for this
project is `input/your-system.constrains`.

```text
$constrain
   force constant=0.05
   distance: 1, 2, auto
   distance: 1, 36, auto
   angle:    2, 1, 36, auto
$end
```

| key | meaning |
|---|---|
| `force constant=` | Hartree; **default 0.05** when omitted |
| `reference=<file>` | take `auto` values from this geometry, not the input one |
| `distance:` / `bond:` | `i, j, <Å or auto>` |
| `angle:` | `i, j, k, <degrees or auto>` |
| `dihedral:` | `i, j, k, l, <degrees or auto>` |
| `bondrange:` | `i, j, rmin, rmax [, beta [, T]]` |
| `atoms:` | rigid-body restraint — see the warning below |
| `elements:` | select by element symbol |

### `atoms:` is a rigid-body restraint, not a pin

This is the part people get wrong.  `atoms: 1-36` does **not** freeze those
atoms in space.  CREST expands it into *all pairwise distance constraints among
the selected atoms* (`parse_xtbinput.f90` only constrains pairs where **both**
atoms are in the list).  So the fragment

- keeps its internal geometry (rigid),
- but can still **translate and rotate as a whole**,
- and its position relative to the *unselected* atoms is not restrained at all.

Cost grows as n(n-1)/2: 36 atoms → 630 constraints, 50 atoms → 1225.  Prefer a
handful of `distance:` / `angle:` lines unless you really want a rigid block.

There is no user-facing "freeze these atoms" option in the current code —
`calc%freezelist` exists but is only set internally (protonation routines).

### MTD: exclude the constrained part from the RMSD

MTD measures RMSD over **all** atoms by default, so a fragment you just made
rigid still contributes and the MTD bias will push against your own constraint.
Put this outside any `$`-block:

```text
atomlist+ 37-104
```

(`atomlist-` removes the listed atoms instead.)

### Two more gotchas

- A constraining file turns **automatic zsort off** (`zsortwarning2` in
  `confparse.f90`) — intentional: atom numbering has to stay stable.
- A very stiff force constant fights both the optimiser and the MTD bias.  Start
  at the default 0.05 and only raise it if the constraint is visibly violated.

---

## 8. Bond orders (SHAKE, NCI mode, wbo/flexi tools)

`rd_ORCA_wbo()` now parses ORCA's Mayer bond-order block from `ORCA.out`:

```
 Mayer bond orders larger than 0.1
B( 0-C , 1-S ) : 1.0550
B( 0-C , 2-H ) : 0.9152
```

* indices are **0-based in ORCA**, converted to CREST's 1-based convention;
* the matrix is symmetrised;
* if the block is absent, CREST prints a warning and continues with zeroed WBOs
  instead of aborting the run.

ORCA enables the Mayer population analysis by default; add `! MAYER` to
`orca.inp` to be explicit.

> Note: the flexibility measure at the start of a run always uses the **internal
> GFN0** engine, so the MD length is chosen identically for ORCA and xtb runs.
> ORCA WBOs are only needed for `-shake`, NCI mode and the `wbo`/`flexi` tools.

---

## 9. If ORCA does not run: read `ORCA.err`

This build does **not** discard ORCA's stderr.  Every ORCA call is

```
cd calculation.level.<n> && <orca_cmd> ORCA.in > ORCA.out 2> ORCA.err
```

so after a failure there are three files to inspect in `calculation.level.<n>/`:

| file | what it tells you |
|---|---|
| `ORCA.err` | start-up failures: missing shared libraries, wrong binary path, MPI problems |
| `ORCA.out` | ORCA's normal output — input errors, SCF failures, `ORCA TERMINATED NORMALLY` |
| `ORCA.engrad` | the file CREST actually reads; if it is absent, ORCA never finished |

> Upstream CREST (and every other subprocess calculator in it) appends
> `2>/dev/null` to that command line, which silently *throws stderr away*.  A
> broken ORCA start-up then produces no output at all and CREST only reports the
> unhelpful `Initial geometry optimization failed!`.  That redirect is replaced by
> `2> ORCA.err` in this tree (see `src/calculator/orca_sc.f90`, `ORCA_setup`).

### The classic case: `LD_LIBRARY_PATH`

An `orca_*_shared_openmpi418*` package is **dynamically linked** against OpenMPI
4.1.8 and ships its own shared objects.  It only starts if the loader can find
them.  Without that you get, in `ORCA.err`:

```
.../orca_startup_mpi: error while loading shared libraries:
libmpi.so.40: cannot open shared object file: No such file or directory
```

Fix it in the run script, **before** `crest` is called:

```bash
export ORCA_DIR=/PATH/TO/orca_dir
export PATH=$ORCA_DIR:$PATH
export LD_LIBRARY_PATH=$ORCA_DIR:$ORCA_DIR/lib:$LD_LIBRARY_PATH
# only if ORCA uses a separately installed OpenMPI 4.1.8:
# export LD_LIBRARY_PATH=/path/to/openmpi-4.1.8/lib:$LD_LIBRARY_PATH
```

> Why it can work "by hand" but not from a script: a non-interactive
> `bash run.sh` does **not** read `~/.bashrc`, so anything exported only there is
> missing unless the calling shell already had it.

### 60-second standalone test

Do not debug through CREST — test ORCA on its own first:

```bash
cd /tmp
printf '! HF def2-SVP\n*xyz 0 1\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.24\n*\n' > w.inp
/PATH/TO/orca w.inp > w.out 2> w.err
echo "exit=$?"; cat w.err
grep -c 'ORCA TERMINATED NORMALLY' w.out
```

* `exit=0` and `ORCA TERMINATED NORMALLY` -> binary and environment are fine, the
  problem is CREST-side.
* non-empty `w.err`, or no such line -> that is your real error; fix it before
  touching anything else.

---

## 10. Troubleshooting

| symptom | cause / fix |
|---|---|
| `ORCA input template <file> could not be found!` | `orca_template` path is relative to the working directory |
| `please define only single jobs (no $new_job)!` | remove `$new_job` |
| `only run EnGrad jobs!` | remove `Opt` / `Freq` / `MD` |
| `please remove %coords block!` | delete the `%coords` block |
| `ORCA cmd line not specified, falling back to 'orca'` | set `orca_cmd` |
| `Initial geometry optimization failed!` | ORCA returned non-zero / no `ORCA.engrad`.  Read `calculation.level.1/ORCA.err` **first**, then `ORCA.out` |
| `gradient file not found` / very first call fails | ORCA did not run — check `calculation.level.*/ORCA.err` and `ORCA.out` |
| `error while loading shared libraries: libmpi.so.40` | `LD_LIBRARY_PATH` does not cover ORCA + OpenMPI — see section 8 |
| `ORCA bond order reader not yet implemented` | you are on an unpatched tree — this build fixes it |
| every structure identical / MD "terminated EARLY" | stale gradient reuse — needs the `crest-xtbsys-mtd-fix` patch |

The per-call working directory is `calculation.level.<n>` (and
`calculation.level.<n>_<i>` for the `i`-th parallel worker), so `ORCA.in`,
`ORCA.out`, `ORCA.err` and `ORCA.engrad` of the last call are always available
there.  `crestopt.log` (in the main working directory) contains only CREST's own
optimisation record — never ORCA's error messages.

### "ORCA is much slower inside CREST than when I run it by hand"

A factor of 2–3 is common and is almost never CREST's fault.  CREST does not call
ORCA the way you do by hand.  The exact command is

```
cd calculation.level.N && OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
   /path/to/orca ORCA.in --oversubscribe > ORCA.out 2> ORCA.err
```

so there are **four** differences to your manual run, in order of how often they
turn out to be the culprit:

1. **MPI ranks instead of threads.**  CREST injects `%pal nprocs N` (N = the
   per-job core budget).  For a *small* level of theory — MINIX and a few hundred
   basis functions — N MPI ranks are frequently **slower** than 1 process with N
   threads: every SCF iteration pays an Allreduce of the Fock matrix, which is
   pure overhead at that size.  If your manual run used `OMP_NUM_THREADS=8`
   without `%pal`, you compared 1×8 threads against 8×1 ranks.
   → Try `CREST_ORCA_NPROCS=1 CREST_ORCA_OMP_THREADS=8` (rebuild needed, see below)
   → quick, no-rebuild probe: `CREST_ORCA_NPROCS=1` alone (then ORCA is fully
   serial — if that is *still* faster per iteration, MPI is the problem).
2. **`--oversubscribe` disables OpenMPI process binding.**  Ranks can then land
   across sockets / migrate.  Look for the corresponding OpenMPI warning in
   `ORCA.err`.
   → `CREST_ORCA_NO_OVERSUBSCRIBE=1`, and/or
   `CREST_ORCA_MPI_ARGS="--bind-to core"`.
3. **Not the same `%maxcore`.**  The template in this repository ships
   `%maxcore 1000`.  If your manual test used ORCA's default (larger), you were
   not running the same job.  Copy `orca.inp` verbatim, or set
   `export OMP_STACKSIZE`/`%maxcore` identically in both.
4. **Filesystem.**  CREST works in `calculation.level.N` **under the CREST job's
   directory**; if that is an NFS home/share while your manual run used node-local
   scratch, every ORCA scratch write is a network round trip.  Check with
   `df -h .` in the job directory and, if it is a network FS, `cd` to a
   node-local scratch before `sbatch` (or add `cd /scratch/$SLURM_JOB_ID` to the
   job script).

### Move the scratch directories off the network filesystem

Every single point runs inside `calculation.level.N`, created **below CREST's
current working directory**.  On an NFS / spinning-disk workspace every scratch
write of ORCA is a network round trip — for a cheap level of theory that is very
often the dominant cost, and it is invisible in the ORCA timings.

Set

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
```

and the directory becomes `<root>/calculation.level.N` (absolute, so nothing else
changes).  The root is created with `mkdir -p`.  Two rules:

- **The root must be unique per job.**  Two concurrent CREST runs sharing one root
  would both use `.../calculation.level.1`.  `$SLURM_JOB_ID` is the easy way.
- The relocated directory is removed by CREST's clean-up, but only the
  `calculation.level.*` entries — never the root itself, so
  `CREST_SCRATCH_DIR=/tmp` cannot wipe `/tmp`.

The job scripts in `slurm-scripts/` already export it, print `df -h /tmp`, and
delete the root only when CREST finished normally (on failure it is kept, which is
exactly when you need `ORCA.out` / `ORCA.err`).

Two things to check **before** trusting `/tmp`:

1. **Is `/tmp` really disk and not RAM?**  A lot of nodes mount `/tmp` as tmpfs,
   i.e. it comes out of your memory request.  A `%maxcore 1000` ORCA scratch per
   rank plus the trajectory files can then push the job over its memory limit and
   get it killed with an unhelpful OOM message.  Check with
   `df -hT /tmp | tail -1` — if the *type* column says `tmpfs`, use a real local
   disk path instead (`/scratch`, `/local`, `/var/tmp`, ...).  The job script
   prints `df -h /tmp` at start-up precisely so you can see the size.
2. **The path must not contain spaces.**  It is passed to `mkdir -p` quoted, but
   CREST builds further shell commands from `calcspace` that are not.

**No-rebuild alternative:** `calcspace` is a normal TOML key
(`parse_calcdata.f90`, aliases `dir`), so

```toml
[calculation.level]
calcspace = "/tmp/crest-12345/calculation.level.1"
```

works on an unmodified CREST too — but you must create the parent directory
yourself, and you have to name the level index explicitly.

### Where the other intermediate files go

`calculation.level.N` is not the only thing CREST drops in your workspace.  There
are also `coord`, `wbo`, `crestopt.xyz`, `gfnff_adjacency`, `gfnff_topo`, `.CHRG`,
`.UHF`, and the directories `MDFILES/`, `METADYN/`, `NORMMD/`, `MRMSD/`, `OPTIM/`.
They are now split by **write frequency**:

| class | what | where |
|---|---|---|
| **hot** (written every few MD steps) | `MDFILES/`, `METADYN/`, `NORMMD/`, `MRMSD/`, `OPTIM/` | `CREST_SCRATCH_DIR` (node-local) |
| **cold** (written a handful of times per run) | `coord`, `wbo`, `crestopt.xyz`, `gfnff_adjacency`, `gfnff_topo`, `.CHRG`, `.UHF` | `CREST_WORK_DIR`, default `crest_work/` |
| **results** (what you actually want) | `crest_conformers.xyz`, `crest_rotamers.xyz`, `crest_best.xyz`, `crest_ensemble.xyz`, `crest_dynamics.xyz`, ... | unchanged, in the workspace |

The hot/cold names are put in place as **symlinks** at start-up
(`crest_stage_files()` in `iomod.F90`, called once from `crest_main.f90`), so the
layout the rest of CREST sees is unchanged.  Consequences worth knowing:

- **Failure is benign.**  If a name is not in the tables, or CREST removes a link,
  the file simply lands in the workspace again.  The run stays correct, you just
  lose the speed-up / tidiness.
- **Existing files are never touched.**  `crest_link()` returns immediately if
  `inquire()` says the name already exists, so a `coord` you placed by hand is
  safe.
- **`MDFILES/` and `TRIALMD/` are relocated by path, not by symlink.**  CREST
  calls `rmrf()` + `makedir()` on them, which would delete the link; they use
  `crest_hot_dir()` instead.  This also means the scratch root must already exist
  (`iomod%makedir()` is a bare `mkdir()` and cannot create parents) — it does,
  because `crest_scratch_dir()` made it with `mkdir -p`.
- **Everything is off when `CREST_STAGE_FILES=0`**, and only the cold half is off
  with `CREST_WORK_DIR=none`.  `CREST_SCRATCH_DIR` is independent of both.

If you want a completely clean workspace *and* a completely local run, set all
three:

```bash
export CREST_SCRATCH_DIR=/tmp/crest-$SLURM_JOB_ID
export CREST_WORK_DIR=crest_work      # or: none, to keep upstream layout
export CREST_STAGE_FILES=1            # 0 disables the hot/cold relocation
```

**Make the comparison apples-to-apples first.**  Take one frame, copy the *exact*
`ORCA.in` CREST wrote (`calculation.level.N/ORCA.in` — it has the coordinates, the
auto-generated `%pal`, charge and multiplicity in it), and run it by hand **from
that same directory** with the same env prefix and the same `--oversubscribe`.  If
it is then also ~2 s, the difference is the invocation/filesystem, not CREST.

Then diff the two `ORCA.out` files: the header states the number of MPI processes
and the per-process memory limit, and the timing block at the end breaks the wall
time into sections — that tells you which section grew.

---

## 11. `ERROR: non-ASCII character(s) found in words` — CREST wrote a binary `ORCA.in`

If `ORCA.err` is empty and `ORCA.out` ends with

```
ERROR: non-ASCII character(s) found in words on line 24: 64_shared_open?
ERROR: non-ASCII character(s) found in words on line 25: ??═══╝██╔════╝╚══██╔══╝
...
[file orca_tools/Tool-Scanner/qcscan1.cpp, line 264]:
      Please remove all non-ASCII characters from your input file
```

this is **not** an ORCA or environment problem — CREST fed ORCA a corrupted
input file.

### Why

`src/calculator/subprocess_types.f90` declared the template buffer as

```fortran
character(len=:),allocatable :: input(:)     ! deferred-length ARRAY
```

A deferred-length character **array** component is silently mis-copied whenever
the enclosing derived type is assigned — **gcc PR fortran/89200**, a GCC 9
regression (8.2.0 is fine; ifort/ifx have the same class of bug).  After the copy
the elements point at stale heap memory instead of the template.

CREST copies calculation settings in several places
(`setuptest.f90: tmpcalc = env%calc`, `calc%add()` while parsing the TOML), and
this is the **only** deferred-length character array component in the whole code
base — which is exactly why only ORCA is affected while `coord`,
`crest_input_copy.xyz` and the systemcall string all stay correct.

The result is an `ORCA.in` with the correct structure (template block,
`*xyz <chrg> <mult>`, coordinates, `*`) but with the template lines replaced by
NUL bytes and fragments of CREST's own memory — its startup banner, the
`jobdescription()` table, the systemcall string.

### Check

```bash
cd calculation.level.1
file ORCA.in                # reports "data" instead of "ASCII text"
grep -c $'\x00' ORCA.in     # non-zero
```

### Proper fix (needs a rebuild)

The element length is fixed now:

```fortran
integer,parameter :: orca_linelen = 1056
character(len=orca_linelen),allocatable :: input(:)
```

Fixed-length elements are copied by a plain memcpy and cannot be mangled.  Use
the patched tree (`crest-xtbsys-mtd-fix.patch`, or the ready-made `crest-now`
source tree) and recompile.

### Follow-up trap: `double free or corruption (out)` right after "fixing" this

Changing the element to a fixed length but leaving the allocation alone trades
PR89200 for a second bug:

```fortran
width = min(max(width+10,1),orca_linelen)
allocate (self%input(nlines),source=repeat(" ",width))     ! WRONG
```

`width` comes from `getlines()` — the longest template line, typically ~85 — so
the SOURCE length no longer matches the declared element length (1056).  That is
**gcc PR fortran/113793**: gfortran <= 13 does *not* diagnose the mismatch and
sizes the elements by the SOURCE length instead; every later
`self%input(k) = trim(atmp)` then writes 1056 bytes into an ~85-byte element.
The overwrite is silent and the process dies at the next `free()`:

```
double free or corruption (out)
Program received signal SIGABRT
#10  read_orca_input  at src/calculator/subprocess_types.f90:153
```

Fixed in GCC 14 (r14-9996, 2024-04-13).  The patched tree makes the SOURCE
length equal the element length, which is correct on every compiler:

```fortran
allocate (self%input(nlines),source=repeat(" ",orca_linelen))   ! correct
```

**PR89200 and PR113793 must be fixed together.**  If you see `double free`
pointing at `subprocess_types.f90`, you are running a half-applied PR89200 fix.

### Workaround (no rebuild)

Make the template contribute **zero** lines, so there is nothing left to corrupt:
CREST then writes only the coordinate block, which is generated fresh and *is*
correct.  Put the method/basis into a wrapper used as `orca_cmd`.

`orca.zero.inp` — the template; CREST's reader strips a `*xyz` block, so both
lines are skipped and `nlines = 0`:

```
*xyz 0 1
*
```

`orca-with-header.sh` — the wrapper (`chmod +x`):

```bash
#!/bin/bash
# CREST's template lines are corruptible (gcc PR89200), so ORCA.in arrives
# containing only the correct coordinate block.  Prepend the header here.
cat > ORCA.hdr <<'EOF'
! PBE D3 def2-SV(P) def2/J RIJCOSX EnGrad
%maxcore 2000
EOF
cat ORCA.hdr ORCA.in > ORCA.full && mv ORCA.full ORCA.in
exec /PATH/TO/orca "$@"
```

and in the TOML:

```toml
orca_cmd      = "/abs/path/to/orca-with-header.sh"
orca_template = "orca.zero.inp"
```

> Do **not** combine the wrapper with the normal template: every template line
> CREST keeps would still be corrupted.
