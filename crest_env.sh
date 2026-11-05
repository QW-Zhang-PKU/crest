#!/usr/bin/env bash
# =============================================================================
#  CREST runtime environment.  SOURCE this file (do not execute it):
#
#      source ./crest_env.sh
#      source ./crest_env.sh 40        # override the thread count
# =============================================================================

# --- SITE-SPECIFIC SETTINGS --------------------------------------------------
# Edit the two defaults below, or export them before sourcing this file:
#     export CREST_BIN=/path/to/crest/bin
#     export XTB_PREFIX=/path/to/xtb
#     source ./crest_env.sh 40
# build-now.sh installs to $HOME/crest-dist/bin by default.
CREST_BIN="${CREST_BIN:-/PATH/TO/crest/bin}"
XTB_PREFIX="${XTB_PREFIX:-/PATH/TO/xtb}"
export PATH="$CREST_BIN:$PATH"

# --- thread budget -----------------------------------------------------------
# CREST's own -T flag sets this internally as well, but setting it here makes
# sure that every child process (and the BLAS backend) starts with the same
# budget instead of "all cores".
export OMP_NUM_THREADS="${1:-40}"

# --- BLAS/LAPACK backends ----------------------------------------------------
# These are only a SAFE DEFAULT for processes CREST does not manage.
#   * in-process (tblite) runs: CREST pins its own BLAS to one thread through
#     openblas_set_num_threads()/mkl_set_num_threads() in new_ompautoset(),
#     so the value exported here does not matter for CREST itself.
#   * external (subprocess) runs: CREST OVERWRITES OMP_NUM_THREADS,
#     MKL_NUM_THREADS and OPENBLAS_NUM_THREADS for every child with
#     cores_per_job - see ompenvset() in src/ompmklset.F90.  With -T 40 and
#     6 MTDs that is 6, not the 1 exported here.  (Until the 3-argument
#     setenv() binding was fixed in src/iomod.F90 that overwrite silently
#     failed whenever the variable already existed, which is exactly why the
#     child inherited these values instead.)
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# --- OpenMP thread stacks ----------------------------------------------------
# GFN-xTB gradients allocate sizeable private arrays; the default 8 MB stack
# will segfault on larger molecules when running threaded.
export OMP_STACKSIZE=1G
ulimit -s unlimited 2>/dev/null || true

# --- external xtb binary (subprocess mode via TOML method = "gfn") ----------
# Parameters for solvation etc. live in $XTBPATH.  The binary is added to
# PATH below; TOML can also reference it by absolute path.
export XTBPATH="$XTB_PREFIX/share/xtb"
export PATH="$XTB_PREFIX/bin:$PATH"

# --- sanity ------------------------------------------------------------------
if command -v crest >/dev/null 2>&1; then
  echo "crest      : $(command -v crest)"
  echo "threads    : OMP_NUM_THREADS=$OMP_NUM_THREADS  (BLAS pinned to 1)"
  echo "stack      : OMP_STACKSIZE=$OMP_STACKSIZE, ulimit -s $(ulimit -s)"
  if command -v xtb >/dev/null 2>&1; then
    echo "xtb        : $(command -v xtb)"
  else
    echo "xtb        : not in PATH (only needed for subprocess mode)"
  fi
  echo
  echo "in-process :  crest molecule.xyz -T $OMP_NUM_THREADS -gfn2 -chrg 0"
  echo "subprocess :  crest coord your-system-xtb.toml -T $OMP_NUM_THREADS --quick"
else
  echo "WARNING: 'crest' not found in PATH ($CREST_BIN)" >&2
  echo "         fix CREST_BIN in this file, or run: make -C _build install" >&2
fi
