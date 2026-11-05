#!/usr/bin/env bash
# =============================================================================
#  build-single-file.sh — build CREST as ONE self-contained executable
#
#  Usage (on the Linux build machine):
#
#      bash build-single-file.sh                     # Intel + static MKL -> ~/crest-single/crest
#      DIST=~/dist bash build-single-file.sh         # where to drop the final file
#      JOBS=16 bash build-single-file.sh
#      TOOLCHAIN=gnu bash build-single-file.sh       # gfortran + static OpenBLAS instead
#      STRIP=0 bash build-single-file.sh             # keep symbols (default: strip)
#      BIN=/path/to/already/built/crest bash build-single-file.sh   # skip build,
#                                                    # just verify + package
#
#  ---------------------------------------------------------------------------
#  WHAT "SINGLE FILE" MEANS HERE (and what it does NOT mean)
#  ---------------------------------------------------------------------------
#  CREST is already a single ELF: CMake builds
#      crest-exe  ->  links "lib-crest-static" (CMakeLists.txt:281-286)
#  and installs exactly one program, $PREFIX/bin/crest.  Nothing else in
#  $PREFIX is needed to *run* it (the installed .mod/.a files and the license
#  files are for downstream compilation only).
#
#  What the default build does NOT give you is a binary you can scp to a
#  compute node: it is dynamically linked against whatever the build host had,
#  typically
#      libmkl_*.so  libiomp5.so  libifcore.so  libimf.so  libsvml.so  libirc.so
#  (Intel) or libopenblas.so / libgomp.so / libgfortran.so (GNU).  A node
#  without oneAPI installed then fails with "error while loading shared
#  libraries".
#
#  This script therefore pins the *link* mode, not the layout:
#
#      STATIC=intel  ->  -static-intel -static-libgcc -static-libstdc++
#      LA=mkl-static ->  MKL archives instead of .so
#
#  After that the ONLY dynamic dependencies left are glibc's own
#  (libc/libm/libdl/libpthread/ld-linux).  Those exist on every Linux node, so
#  the file is effectively self-contained.  That is what the final gate checks.
#
#  IT IS STILL NOT A "DISTRIBUTE ANYWHERE" BINARY:
#    * glibc is dynamic, so the target node needs glibc >= the build host's.
#      Build on the OLDEST node you want to run on.
#    * xtb / orca are separate programs.  CREST calls them via system()
#      (env%ProgName, src/confparse.f90:137, default 'xtb'); they are NOT
#      linked into the binary and never can be.
#    * Nothing else is read from disk at start-up - the GFN parameters are
#      compiled in via tblite / gfn0 / gfnff.
#
#  ---------------------------------------------------------------------------
#  THE THREE STATIC LEVELS (do not skip this)
#  ---------------------------------------------------------------------------
#    1. STATIC=intel  (DEFAULT, recommended)
#       Intel/GNU runtime + MKL static, glibc dynamic.  Intel supports this for
#       OpenMP.  Result: no libifcore/libimf/libsvml/libirc/libiomp5 needed.
#
#    2. STATIC=all  ->  plain "-static" (glibc included too)
#       Intel explicitly does NOT recommend this for OpenMP codes: libiomp5
#       plus a static glibc can deadlock, and dlopen()/NSS/getaddrinfo() warn
#       or misbehave.  CREST forks external programs, which is exactly the
#       dlopen/exec-heavy case.  This script refuses it unless FORCE_STATIC=1.
#
#    3. upstream -DSTATICBUILD=ON (config/CMakeLists.txt:41)
#       Same "-static" as (2), plus CMAKE_FIND_LIBRARY_SUFFIXES=".a".  Caveat:
#       it also disables the object library and the shared library, and CREST's
#       own target stops linking OpenMP::OpenMP_Fortran directly
#       (CMakeLists.txt:190-200).  The -qopenmp flag then only survives because
#       tblite::tblite propagates it (subprojects/tblite/CMakeLists.txt:63-95,
#       'target_link_libraries(tblite-lib PUBLIC ${lib-deps})').
#       So: STATICBUILD=ON is fine here, but never combine it with a
#       no-tblite build.  This script does not use it.
#
#  ---------------------------------------------------------------------------
#  Environment switches
#  ---------------------------------------------------------------------------
#   TOOLCHAIN   intel (default: ifx/ifort + static MKL) | gnu (gfortran +
#               static OpenBLAS; needs bash build-openblas-static.sh first)
#   DIST        output dir for the finished file      (default ~/crest-single)
#   JOBS        compile parallelism                   (default nproc)
#   PREFIX      internal cmake install prefix         (default $DIST/_install)
#   LA          linear algebra backend                (default mkl-static)
#               mkl-static | openblas-static are the only ones that make sense
#               here; a shared BLAS defeats the whole point and is refused
#               unless FORCE_SHARED_BLAS=1.
#   STATIC      intel (default) | all (refused unless FORCE_STATIC=1) | 0
#   OPENBLAS_DIR  prefix holding lib/libopenblas.a    (default ~/opt/openblas-static)
#   MARCH       extra arch flags, e.g. -xHost        (default: none = portable)
#   OMP_STATIC  1 to also add -qopenmp-link=static (only needed if libiomp5.so
#               still shows up in the gate below; -static-intel normally covers
#               it already, and Intel prefers the dynamic libiomp5)
#   STRIP       1 (default) strip the finished file, 0 keep it
#   SKIP_VERIFY 1 to skip the verify-fixes.sh provenance gate
#   NO_SUBMODULE 1 to skip git submodule init (auto if subprojects/ present)
#   SETVARS     path to Intel oneAPI setvars.sh      (auto-probed if unset)
#   BIN         skip the build, verify + package this binary
#   FORCE_STATIC / FORCE_SHARED_BLAS / FORCE_DEPS
#               1 to turn the corresponding hard error into a warning
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLCHAIN="${TOOLCHAIN:-intel}"
DIST="${DIST:-$HOME/crest-single}"
PREFIX="${PREFIX:-$DIST/_install}"
JOBS="${JOBS:-$( (command -v nproc >/dev/null && nproc) || getconf _NPROCESSORS_ONLN || echo 4 )}"
STRIP="${STRIP:-1}"
SKIP_VERIFY="${SKIP_VERIFY:-0}"
FORCE_STATIC="${FORCE_STATIC:-0}"
FORCE_SHARED_BLAS="${FORCE_SHARED_BLAS:-0}"
FORCE_DEPS="${FORCE_DEPS:-0}"

if [ "$TOOLCHAIN" = "intel" ]; then
  LA="${LA:-mkl-static}"
  STATIC="${STATIC:-intel}"
  OPENBLAS_DIR="${OPENBLAS_DIR:-$HOME/opt/openblas-static}"
else
  LA="${LA:-openblas-static}"
  STATIC="${STATIC:-gnu-runtime}"
  OPENBLAS_DIR="${OPENBLAS_DIR:-$HOME/opt/openblas-static}"
fi

say()  { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33m  !! %s\033[0m\n' "$*"; }
die()  { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 1) preflight - refuse combinations that silently produce a NON-portable file
# ---------------------------------------------------------------------------
say "Preflight"

case "$STATIC" in
  all|full|2)
    if [ "$FORCE_STATIC" = "1" ]; then
      warn "STATIC=all on purpose: full '-static' with OpenMP can deadlock"
      warn "(libiomp5 + static glibc).  You asked for it - not my problem."
    else
      die "STATIC=all means plain '-static' (glibc included).
     Intel does not recommend that for OpenMP: libiomp5 plus a static glibc
     can deadlock, and CREST forks xtb/orca constantly (exec + dlopen heavy).
     Use STATIC=intel (the default) - the remaining glibc dependency is
     harmless because every Linux node has it.
     If you really want it: FORCE_STATIC=1 bash $0"
    fi
    ;;
esac

case "$LA" in
  mkl-static|openblas-static) : ;;
  *)
    if [ "$FORCE_SHARED_BLAS" = "1" ]; then
      warn "LA='$LA' is a SHARED BLAS - the binary will need that .so on the"
      warn "compute node, which defeats the purpose of this script."
    else
      die "LA='$LA' is a shared BLAS/LAPACK.
     The finished file would still need those .so files at run time, i.e. it is
     NOT self-contained.  Use LA=mkl-static (Intel, default) or LA=openblas-static
     (GNU; build it first with 'bash build-openblas-static.sh', then point
     OPENBLAS_DIR at it).
     Proceed anyway:  FORCE_SHARED_BLAS=1 bash $0"
    fi
    ;;
esac

[ "$TOOLCHAIN" = "intel" ] || [ "$TOOLCHAIN" = "gnu" ] \
  || die "unknown TOOLCHAIN='$TOOLCHAIN' (use intel | gnu)"

for t in ldd; do command -v "$t" >/dev/null 2>&1 || warn "'$t' not found - the dependency gate will fall back to readelf"; done

echo "  toolchain  : $TOOLCHAIN"
echo "  BLAS       : $LA"
echo "  static     : $STATIC"
echo "  jobs       : $JOBS"
echo "  dist       : $DIST"

# ---------------------------------------------------------------------------
# 2) build (delegated - do NOT re-implement cmake here)
# ---------------------------------------------------------------------------
mkdir -p "$DIST"
BIN="${BIN:-}"

if [ -n "$BIN" ]; then
  [ -f "$BIN" ] || die "BIN='$BIN' does not exist"
  say "Skipping the build (BIN given) - verifying $BIN"
else
  case "$TOOLCHAIN" in
    intel)
      HELPER="${HELPER:-$HERE/build-now-intel.sh}"
      [ -f "$HELPER" ] || die "'$HELPER' not found.
     This script is a thin front-end: it pins the link mode and then calls
     build-now-intel.sh, so both files must sit in the CREST source root.
     Override with HELPER=/path/to/build-now-intel.sh"
      say "Building via build-now-intel.sh (STATIC=$STATIC LA=$LA)"
      # NOTE: do not add -DCMAKE_EXE_LINKER_FLAGS of your own - it may only be
      # passed once, and build-now-intel.sh already emits it (line ~434).
      PREFIX="$PREFIX" JOBS="$JOBS" LA="$LA" STATIC="$STATIC" \
      WITH_MKL="${WITH_MKL:-0}" SKIP_VERIFY="$SKIP_VERIFY" \
      OMP_STATIC_LINK="${OMP_STATIC:-0}" \
      NO_SUBMODULE="${NO_SUBMODULE:-}" SETVARS="${SETVARS:-}" \
      MARCH="${MARCH:-}" BUILD_TYPE="${BUILD_TYPE:-Release}" \
        bash "$HELPER"
      ;;
    gnu)
      HELPER="${HELPER:-$HERE/build_crest.sh}"
      [ -f "$HELPER" ] || die "'$HELPER' not found (override with HELPER=...)"
      say "Building via build_crest.sh (static OpenBLAS + GNU runtime)"
      [ -f "$OPENBLAS_DIR/lib/libopenblas.a" ] || die \
"libopenblas.a not found in $OPENBLAS_DIR/lib - run 'bash build-openblas-static.sh' first
     (or OPENBLAS_DIR=<prefix> if you already have a static OpenBLAS)."
      # '-static-libgfortran'/'libgcc'/'libstdc++' cover the GNU runtime.
      # libgomp has no -static-libgomp: GCC gives no supported static switch
      # for it, so it may stay dynamic.  The gate below reports the truth.
      STATIC_BLAS=1 OPENBLAS_DIR="$OPENBLAS_DIR" \
      STATIC_BLAS_EXTRA_LIBS="${STATIC_BLAS_EXTRA_LIBS:--static-libgcc -static-libstdc++ -static-libgfortran}" \
      JOBS="$JOBS" \
        bash "$HELPER" "$HERE" "$PREFIX"
      ;;
  esac

  for cand in "$PREFIX/bin/crest" "$HERE/_build-intel/crest" "$HERE/_build/crest"; do
    if [ -x "$cand" ]; then BIN="$cand"; break; fi
  done
fi

[ -n "$BIN" ] && [ -f "$BIN" ] || die "no crest binary found - looked in
     $PREFIX/bin/crest
     $HERE/_build-intel/crest
     $HERE/_build/crest
   Re-run with BIN=/path/to/crest to package an existing build."

say "Binary: $BIN"

# ---------------------------------------------------------------------------
# 3) GATE 1 - is the file actually self-contained?
#    Allowed: glibc's own objects.  Anything compiler/BLAS/OpenMP specific is
#    a hard failure, because that is exactly what is missing on a compute node.
# ---------------------------------------------------------------------------
DEPS_RAW=""
if command -v ldd >/dev/null 2>&1; then
  DEPS_RAW="$(ldd "$BIN" 2>/dev/null || true)"
elif command -v readelf >/dev/null 2>&1; then
  DEPS_RAW="$(readelf -d "$BIN" 2>/dev/null | sed -n 's/.*NEEDED.*\[\(.*\)\].*/\1/p' || true)"
fi

say "Dependency gate"
case "$DEPS_RAW" in
  *"not a dynamic executable"*|*"statically linked"*)
    echo "  fully static - no shared libraries at all"
    DEPS=""
    ;;
  *)
    DEPS="$(printf '%s\n' "$DEPS_RAW" | awk 'NF{print $1}' | grep -E '\.so' || true)"
    ;;
esac

MISSING="$(printf '%s\n' "$DEPS_RAW" | grep 'not found' || true)"
[ -n "$MISSING" ] && { warn "unresolved libraries:"; printf '%s\n' "$MISSING"; }

BAD=""
while IFS= read -r d; do
  [ -n "$d" ] || continue
  b="${d##*/}"
  case "$b" in
    # glibc / the loader: present on every Linux node
    libc.so.*|libm.so.*|libdl.so.*|libdl.so|libpthread.so.*|librt.so.*|ld-linux*.so.*|linux-vdso*.so.*|linux-gate*.so.*) ;;
    # tolerable but note it
    libgcc_s.so.*) echo "  note      : $b (system libgcc, harmless)" ;;
    *)
      BAD="${BAD:+$BAD
}$d"
      ;;
  esac
done <<EOF
$DEPS
EOF

if [ -n "$BAD" ]; then
  warn "the binary still needs these shared libraries:"
  printf '%s\n' "$BAD" | sed 's/^/       /'
  if [ "$FORCE_DEPS" = "1" ]; then
    warn "FORCE_DEPS=1 - continuing anyway.  Copy these .so files along with"
    warn "the binary, or set LD_LIBRARY_PATH on the compute node."
  else
    die "NOT self-contained.
     Typical causes and fixes:
       libiomp5.so / libifcore.so / libimf.so / libsvml.so / libirc.so
           -> STATIC=intel did not take effect.  Check that build-now-intel.sh
              really got STATIC=intel (this script exports it).  If it is only
              libiomp5.so that survives, retry with OMP_STATIC=1 (adds
              -qopenmp-link=static).
       libmkl_*.so
           -> LA=mkl-static did not take effect, or MKLROOT was wrong so CMake
              silently fell back to the shared libs.  Look for
              'BLA_VENDOR=Intel10_64lp' + 'BLA_STATIC=ON' in the configure log.
       libopenblas.so
           -> use OPENBLAS_DIR=<prefix with lib/libopenblas.a> and
              LA=openblas-static.
       libgomp.so / libgfortran.so (GNU)
           -> GCC has no '-static-libgomp'; either accept it (it is in the
              system gcc runtime and usually present) or switch to Intel.
     Override with FORCE_DEPS=1 if you know the node has them."
  fi
else
  echo "  OK - only glibc (libc/libm/libdl/libpthread) remains"
fi

# ---------------------------------------------------------------------------
# 4) GATE 2 - OpenMP really enabled?
#    A statically linked CREST that is silently single-threaded is the worst
#    possible outcome, and nothing in 'ldd' would reveal it.
# ---------------------------------------------------------------------------
say "OpenMP check"
OMP_FLAG_SEEN=0
for bd in "$HERE/_build-intel" "$HERE/_build" "$PREFIX/.."; do
  [ -d "$bd" ] || continue
  FF=$(find "$bd" -name flags.make -path '*crest*' 2>/dev/null | head -n1 || true)
  [ -n "$FF" ] || continue
  if grep -q -E '(^| )-(-qopenmp|fopenmp)( |$)' "$FF" 2>/dev/null; then
    OMP_FLAG_SEEN=1
    echo "  -qopenmp/-fopenmp present in $FF"
    break
  fi
done
if [ "$OMP_FLAG_SEEN" = "0" ]; then
  warn "could not confirm -qopenmp/-fopenmp on any CREST compile line."
  warn "If you built with -DSTATICBUILD=ON, CREST's own target no longer links"
  warn "OpenMP::OpenMP_Fortran (CMakeLists.txt:190-200); the flag then only"
  warn "arrives transitively through tblite::tblite.  Verify with a threaded"
  warn "run before trusting this binary."
fi

# ---------------------------------------------------------------------------
# 5) GATE 3 - does it start with an EMPTY environment?
#    This is the real proof: no LD_LIBRARY_PATH, no Intel vars, nothing.
# ---------------------------------------------------------------------------
say "Clean-environment start-up test"
VTMP="/tmp/_crest_single_version.$$"
env -i "$BIN" --version >"$VTMP" 2>&1 && VRC=0 || VRC=$?
sed 's/^/  /' "$VTMP" || true
# 'crest --version' ends in a bare 'stop' (src/confparse.f90:76-79), so the
# exit code is 0; a non-zero code or a loader complaint is the real signal.
if [ "$VRC" -ne 0 ] || [ ! -s "$VTMP" ] || grep -q 'error while loading shared libraries' "$VTMP"; then
  rm -f "$VTMP"
  die "the binary does not start with an empty environment (exit $VRC).
     That means it still needs something from the environment - see the output
     above.  A truly self-contained crest must run under 'env -i'."
fi
rm -f "$VTMP"
echo "  OK - runs with 'env -i' (no LD_LIBRARY_PATH needed)"

# ---------------------------------------------------------------------------
# 6) package: one file, one place
# ---------------------------------------------------------------------------
say "Packaging"
mkdir -p "$DIST"
OUT="$DIST/crest"
cp -f "$BIN" "$OUT"
chmod 755 "$OUT"
if [ "$STRIP" = "1" ] && command -v strip >/dev/null 2>&1; then
  strip -s "$OUT" 2>/dev/null && echo "  stripped  : yes" || echo "  stripped  : no (strip failed)"
else
  echo "  stripped  : no (STRIP=$STRIP)"
fi

SIZE=$(du -h "$OUT" | cut -f1)
SHA=$( (command -v sha256sum >/dev/null && sha256sum "$OUT" | cut -d' ' -f1) || echo '<no sha256sum>' )
GLIBC=$( (ldd --version 2>/dev/null | head -n1) || echo '<unknown>' )
FTYPE=$( (command -v file >/dev/null && file -b "$OUT") || echo '<no file cmd>' )

echo "  file      : $OUT"
echo "  size      : $SIZE"
echo "  type      : $FTYPE"
echo "  sha256    : $SHA"

# ---------------------------------------------------------------------------
# 7) report
# ---------------------------------------------------------------------------
cat <<EOF

-----------------------------------------------------------------------------
 CREST is now a single self-contained executable.

   $OUT   ($SIZE)

 What travels with it : nothing.  No .so, no module setvars, no LD_LIBRARY_PATH.
 What it still needs  : glibc >= the build host's ($GLIBC)
                        -> build on the OLDEST node you intend to run on.
 What it still calls  : xtb / orca, as separate programs.  They are spawned via
                        system() and can never be linked in
                        (env%ProgName, src/confparse.f90:137).

 Copy it out:

   scp $OUT user@login.cluster:/path/to/bin/crest
   ssh user@login.cluster 'chmod +x /path/to/bin/crest && /path/to/bin/crest --version'

 Run it (ORCA external calculator - remember there is NO CLI flag for ORCA,
 the calculator is chosen by the .toml file only):

   cd <working dir>
   crest struc.xyz crest.toml -T 48 --loadconformer --onebatch

   # crest.toml
   #   method        = "orca"
   #   orca_cmd      = "/path/to/orca"
   #   orca_template = "/path/to/orca.inp"

 Note: crest_env.sh is only needed when the binary is NOT self-contained.  With
 this build you can skip it - but keep MKL_NUM_THREADS=1 in your job script if
 you use crest_env.sh anyway, since CREST re-exports the per-job budget to
 every child (ompenvset).
-----------------------------------------------------------------------------
EOF
