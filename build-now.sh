#!/usr/bin/env bash
# =============================================================================
#  build-now.sh — one-shot build for the fully patched CREST tree in this folder
#
#  Usage (on the Linux machine):
#
#      bash build-now.sh                    # install to $HOME/crest-dist
#      PREFIX=~/opt/crest bash build-now.sh # custom install prefix
#      JOBS=16 bash build-now.sh            # limit compile parallelism
#
#  What it does:
#    1. sanity-checks that all local patches are actually present (verify-fixes.sh)
#    2. if this folder sits on a WSL 9p mount (/mnt/...), stages it to
#       $HOME/crest-now first — compiling 7k files over 9p is painfully slow
#    3. delegates to build_crest.sh (CMake + Ninja/Make, offline, no network)
#
#  Requirements: gfortran >= 10, gcc, cmake (ninja recommended)
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${PREFIX:-$HOME/crest-dist}"
JOBS="${JOBS:-$(nproc)}"

say() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 0) provenance: make sure this tree really carries the local fixes
# ---------------------------------------------------------------------------
say "Checking that the local patches are applied"
if [ -x "$HERE/verify-fixes.sh" ] || [ -f "$HERE/verify-fixes.sh" ]; then
  bash "$HERE/verify-fixes.sh" || die "this tree is NOT fully patched - refusing to build"
else
  echo "  (verify-fixes.sh not found - skipping; not fatal)"
fi

# ---------------------------------------------------------------------------
# 1) avoid building across a WSL 9p / network mount
# ---------------------------------------------------------------------------
SRC="$HERE"
case "$HERE" in
  /mnt/* | /media/* | /run/media/* | /net/* | /Volumes/*)
    STAGE="$HOME/crest-now"
    say "Source is on a slow filesystem ($HERE)"
    echo "  staging to $STAGE and building there instead"
    mkdir -p "$STAGE"
    # cp -a over an existing tree is fine and avoids a bulk rm
    cp -a "$HERE/." "$STAGE/"
    SRC="$STAGE"
    ;;
  *)
    say "Source: $SRC (native filesystem - building in place)"
    ;;
esac

# ---------------------------------------------------------------------------
# 2) toolchain
# ---------------------------------------------------------------------------
say "Toolchain"
for t in gfortran gcc cmake; do
  command -v "$t" >/dev/null 2>&1 || die "'$t' not found in PATH"
done
GFV=$(gfortran -dumpversion | cut -d. -f1)
[ "$GFV" -ge 10 ] || die "gfortran >= 10 required (found $GFV)"
echo "  gfortran   : $(gfortran --version | head -n1)"
echo "  cmake      : $(cmake --version | head -n1)"
if command -v ninja >/dev/null 2>&1; then
  echo "  generator  : Ninja"
else
  echo "  generator  : Unix Makefiles  (install ninja for a faster build)"
fi
echo "  jobs       : $JOBS"
echo "  prefix     : $PREFIX"

# ---------------------------------------------------------------------------
# 3) build (delegates to the shared helper; subprojects are pre-placed)
# ---------------------------------------------------------------------------
say "Building"
NO_SUBMODULE=1 JOBS="$JOBS" bash "$SRC/build_crest.sh" "$SRC" "$PREFIX"

# ---------------------------------------------------------------------------
# 4) what next
# ---------------------------------------------------------------------------
CREST="$PREFIX/bin/crest"
[ -x "$CREST" ] || CREST="$SRC/_build/crest"

cat <<EOF

-----------------------------------------------------------------------------
 Patched CREST is built.

   binary : $CREST

 Run it:

   source $SRC/crest_env.sh 40
   cd <your working directory>          # preferably NOT under /mnt/...
   crest your-system-xtb.toml -T 40 --quick --noreftopo

 Recommended config for the external-xtb MTD path:
   * CREST now exports OMP/MKL/OPENBLAS_NUM_THREADS = cores_per_job (6 with
     -T 40 and 6 MTDs) plus OMP_MAX_ACTIVE_LEVELS=1 for every xtb child.
   * While it runs, verify with:  watch -n 2 $SRC/check-xtb-threads.sh
-----------------------------------------------------------------------------
EOF
