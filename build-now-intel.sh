#!/usr/bin/env bash
# =============================================================================
#  build-now-intel.sh — one-shot Intel (ifx/ifort) build of the patched CREST tree
#
#  Usage (on the Linux machine):
#
#      bash build-now-intel.sh                          # MKL static -> $HOME/crest-dist-intel
#      PREFIX=~/opt/crest-intel bash build-now-intel.sh
#      JOBS=16 bash build-now-intel.sh
#      FC=ifort bash build-now-intel.sh                 # force the classic compiler
#      LA=openblas-static bash build-now-intel.sh       # reuse an existing static OpenBLAS
#      SETVARS=/opt/intel/oneapi/setvars.sh bash build-now-intel.sh
#      WITH_MKL=1 bash build-now-intel.sh               # opt in to WITH_MKL=1 (see below)
#      MARCH="-xHost" bash build-now-intel.sh           # node-specific binary
#
#  Why Intel at all: MKL + the Intel runtime are typically 10-30% faster than the
#  gfortran/OpenBLAS combination, and on a cluster that already ships ORCA you
#  usually have oneAPI installed anyway.
#
#  ---------------------------------------------------------------------------
#  Environment switches (all optional)
#  ---------------------------------------------------------------------------
#   PREFIX         install prefix              (default: $HOME/crest-dist-intel)
#   JOBS           compile parallelism         (default: nproc)
#   FC / CC        Fortran / C compiler        (default: ifx + icx, else ifort + icc)
#   SETVARS        path to Intel oneAPI setvars.sh   (auto-probed if unset)
#   LA             linear algebra backend:
#                    mkl-static    (default) MKL, statically linked
#                                  -> nothing extra needed on the compute node
#                    mkl                     MKL, shared (needs LD_LIBRARY_PATH)
#                    mkl-rt                  MKL via the single mkl_rt shared lib
#                    openblas-static         reuse $OPENBLAS_DIR
#                    openblas                system shared OpenBLAS
#                    auto                    let CMake pick (no forced vendor)
#   OPENBLAS_DIR   prefix holding lib/libopenblas.a   (default $HOME/opt/openblas-static)
#   BUILD_TYPE     cmake build type            (default: Release)
#   OPT            optimisation level          (default: -O3 for Release)
#   MARCH          extra arch flags, e.g. -xHost or -axAVX2  (default: none -> portable
#                  across cluster nodes; -xHost is faster but only runs on the build host)
#   HEAP_ARRAYS    passed through to -heap-arrays, e.g. HEAP_ARRAYS=1 for plain
#                  '-heap-arrays' (all temporaries on the heap) or HEAP_ARRAYS=64
#                  for a size threshold. Use it if threaded runs segfault in the
#                  master thread (stack exhaustion); see 'ifx -help heap-arrays'.
#   FP_PRECISE     1 adds -fp-model=precise    (default 0 = Intel's default fast=1,
#                  which is what upstream's Intel dialect uses)
#   EXTRA_FFLAGS   extra Fortran flags
#   OMP_STATIC_LINK
#                  1 adds -qopenmp-link=static  (only if libiomp5.so is absent on
#                  the compute nodes; NOT Intel's recommended default)
#   STATIC         how much of the runtime to bind into the executable:
#                    0 / none    (default) everything dynamic, as today
#                    intel / 1   -static-intel -static-libgcc -static-libstdc++
#                                Intel + GNU runtime static, libc/libm dynamic.
#                                THIS IS THE ONE YOU WANT for a cluster: the
#                                binary no longer needs libifcore/libimf/libiomp5
#                                at run time, and Intel supports it for OpenMP.
#                    all / 2     -static  (libc included).  Intel does NOT
#                                recommend this for OpenMP codes: libiomp5 plus a
#                                static glibc can deadlock, and dlopen()/
#                                getaddrinfo()/NSS warn or misbehave.  Only for a
#                                compute node with no runtime at all.
#                  Combine with LA=mkl-static (already the default) for a
#                  self-contained binary; a shared BLAS defeats the purpose.
#   WITH_TESTS     ON/OFF                      (default: OFF)
#   WITH_MKL       1 defines WITH_MKL=1        (default 0 = leave it undefined,
#                  exactly like the upstream CMake build).  1 makes CREST call
#                  MKL_Set_Num_Threads/mkl_set_dynamic and tblite call
#                  mkl_free_buffers(); opt-in because that tblite branch is not
#                  exercised by CMake builds.  crest_env.sh already pins MKL to
#                  one thread via MKL_NUM_THREADS=1.
#   SKIP_VERIFY    1 to skip the verify-fixes.sh provenance gate
#   NO_SUBMODULE   1 to skip git submodule init (auto-enabled if subprojects/ exists)
#
#  What it does:
#    1. locates the source root and the Intel toolchain (sources setvars.sh if needed)
#    2. sanity-checks that all local patches are actually present (verify-fixes.sh)
#    3. stages the tree off a WSL 9p / network mount if necessary
#    4. configures with CMake, passing the Intel-specific Fortran flags
#       (config/CMakeLists.txt already supplies -r8 -align array64byte via the
#        MATCHES "Intel" branch - do NOT pass the gfortran flags here)
#    5. builds, installs, smoke-tests, and reports the runtime dependencies
#
#  NOTE on free-form line length: CREST has source lines up to 172 characters.
#  That is fine for Intel (free-form limit is 10000 characters) - the
#  gfortran-only '-ffree-line-length-none' must NOT be passed here, and the
#  fixed-form '-extend-source' is irrelevant for .f90 files.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${PREFIX:-$HOME/crest-dist-intel}"
JOBS="${JOBS:-$( (command -v nproc >/dev/null && nproc) || getconf _NPROCESSORS_ONLN || echo 4 )}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
LA="${LA:-mkl-static}"
OPENBLAS_DIR="${OPENBLAS_DIR:-$HOME/opt/openblas-static}"
FP_PRECISE="${FP_PRECISE:-0}"
WITH_TESTS="${WITH_TESTS:-OFF}"
SKIP_VERIFY="${SKIP_VERIFY:-0}"
OMP_STATIC_LINK="${OMP_STATIC_LINK:-0}"
STATIC="${STATIC:-0}"
WITH_MKL="${WITH_MKL:-0}"
CMAKE="${CMAKE:-cmake}"

say() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33m  !! %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# -1) locate the real source root (the dir that holds CMakeLists.txt + src/)
#     Normally the script sits next to CMakeLists.txt, but accept a subdir
#     invocation (or a wrapped tarball) by searching downward like build_crest.sh.
# ---------------------------------------------------------------------------
SRC_ROOT="${SRC:-$HERE}"
if [ ! -f "$SRC_ROOT/CMakeLists.txt" ] || [ ! -f "$SRC_ROOT/src/CMakeLists.txt" ]; then
  CANDR=$(find "$SRC_ROOT" -maxdepth 3 -name CMakeLists.txt -printf '%h\n' 2>/dev/null | head -n1 || true)
  if [ -n "$CANDR" ] && [ -f "$CANDR/src/CMakeLists.txt" ]; then
    SRC_ROOT="$CANDR"
  else
    # walk upward from wherever the script lives (e.g. src/crest/tools/)
    d="$HERE"; UP=""
    while [ -n "$d" ] && [ "$d" != "/" ]; do
      if [ -f "$d/CMakeLists.txt" ] && [ -f "$d/src/CMakeLists.txt" ]; then UP="$d"; break; fi
      d="$(dirname "$d")"
    done
    [ -n "$UP" ] && SRC_ROOT="$UP"
  fi
fi
{ [ -f "$SRC_ROOT/CMakeLists.txt" ] && [ -f "$SRC_ROOT/src/CMakeLists.txt" ]; } || die \
"no CMakeLists.txt with a src/ subdirectory found near '$HERE'.
     Run this script from the CREST source root, or pass it explicitly:
         SRC=/path/to/crest bash $0"
say "Source root: $SRC_ROOT"

# ---------------------------------------------------------------------------
# 0) provenance: make sure this tree really carries the local fixes
# ---------------------------------------------------------------------------
say "Checking that the local patches are applied"
if [ "$SKIP_VERIFY" = "1" ]; then
  warn "SKIP_VERIFY=1 - provenance gate skipped"
elif [ -f "$SRC_ROOT/verify-fixes.sh" ]; then
  bash "$SRC_ROOT/verify-fixes.sh" || die "this tree is NOT fully patched - refusing to build"
else
  warn "verify-fixes.sh not found - continuing without the provenance gate"
fi

# ---------------------------------------------------------------------------
# 1) avoid building across a WSL 9p / network mount
# ---------------------------------------------------------------------------
HERE="$SRC_ROOT"
SRC="$SRC_ROOT"
case "$HERE" in
  /mnt/* | /media/* | /run/media/* | /net/* | /Volumes/*)
    STAGE="$HOME/crest-now-intel"
    say "Source is on a slow filesystem ($HERE)"
    echo "  staging to $STAGE and building there instead"
    mkdir -p "$STAGE"
    cp -a "$HERE/." "$STAGE/"
    SRC="$STAGE"
    ;;
  *)
    say "Source: $SRC (native filesystem - building in place)"
    ;;
esac

# ---------------------------------------------------------------------------
# 2) toolchain: get the Intel compilers on PATH
# ---------------------------------------------------------------------------
say "Intel toolchain"

if [ -n "${FC:-}" ] && command -v "$FC" >/dev/null 2>&1; then
  echo "  using the already-active toolchain (FC=$FC)"
else
  CAND="${SETVARS:-}"
  if [ -z "$CAND" ]; then
    for c in /opt/intel/oneapi/setvars.sh \
             /mnt/opt/intel/oneapi/setvars.sh \
             /opt/software/intel/oneapi/setvars.sh \
             /mnt/opt/software/intel/oneapi/setvars.sh \
             "$HOME/intel/oneapi/setvars.sh" \
             /opt/intel/bin/compilervars.sh; do
      if [ -f "$c" ]; then CAND="$c"; break; fi
    done
  fi
  if [ -n "$CAND" ]; then
    command -v ifx >/dev/null 2>&1 || command -v ifort >/dev/null 2>&1 || {
      echo "  sourcing $CAND"
      # setvars.sh is noisy and not -u/-e clean; disable those around it
      set +eu
      # shellcheck disable=SC1090
      . "$CAND" >/dev/null 2>&1 || true
      set -eu
    }
  fi
fi

# Fortran compiler: prefer the LLVM-based ifx, fall back to the classic ifort
if [ -n "${FC:-}" ]; then
  command -v "$FC" >/dev/null 2>&1 || die "FC='$FC' not found in PATH"
elif command -v ifx >/dev/null 2>&1; then
  FC=ifx
elif command -v ifort >/dev/null 2>&1; then
  FC=ifort
else
  die "neither 'ifx' nor 'ifort' is in PATH.
     Install Intel oneAPI, or point the script at its environment:
         SETVARS=/path/to/oneapi/setvars.sh bash $0"
fi

# C compiler must match the Fortran compiler family
if [ -n "${CC:-}" ]; then
  command -v "$CC" >/dev/null 2>&1 || die "CC='$CC' not found in PATH"
elif [ "$FC" = "ifx" ] && command -v icx >/dev/null 2>&1; then
  CC=icx
elif command -v icc >/dev/null 2>&1; then
  CC=icc
elif command -v icx >/dev/null 2>&1; then
  CC=icx
else
  die "no Intel C compiler (icx/icc) found in PATH - needed for src/signal.c and src/symmetry_i.c"
fi

command -v "$CMAKE" >/dev/null 2>&1 || die "'$CMAKE' not found in PATH (override with CMAKE=/path/to/cmake)"

FVER=$("$FC" --version 2>/dev/null | head -n1 || echo unknown)
CVER=$("$CC" --version 2>/dev/null | head -n1 || echo unknown)
CMV=$("$CMAKE" --version | head -n1)
echo "  FC         : $FVER"
echo "  CC         : $CVER"
echo "  cmake      : $CMV"

# --- cmake must be new enough to know ifx ----------------------------------
# CMake recognises ifx through the 'IntelLLVM' compiler id, which only exists
# since CMake 3.20.  With an older cmake the identification ends up EMPTY;
# CMake only *warns* ("The Fortran compiler identification is unknown") and
# keeps going, config/CMakeLists.txt then falls into its
#   message(FATAL_ERROR "Please use an Intel or GNU compiler!")
# branch.  Catch it here, where the message can actually be useful.
# Compare major/minor as integers - a float compare gets 3.9 vs 3.20 wrong.
CMNUM=$("$CMAKE" --version | head -n1 | sed -E 's/^[^0-9]*([0-9]+\.[0-9]+).*/\1/')
CMMAJ="${CMNUM%%.*}"
CMMIN="${CMNUM#*.}"
cm_lt() { # true when cmake < $1.$2
  [ -n "$CMNUM" ] || return 1
  [ "$CMMAJ" -lt "$1" ] 2>/dev/null && return 0
  if [ "$CMMAJ" -eq "$1" ] 2>/dev/null && [ "$CMMIN" -lt "$2" ] 2>/dev/null; then
    return 0
  fi
  return 1
}
case "$FC" in
  ifx)
    if cm_lt 3 20; then
      die "cmake $CMNUM is too old for '$FC'.
     ifx identifies itself with the 'IntelLLVM' compiler id, which CMake only
     knows since 3.20.  With cmake $CMNUM the Fortran compiler id comes back
     EMPTY, CMake merely warns ('Fortran compiler identification is unknown')
     and config/CMakeLists.txt aborts with
         'Please use an Intel or GNU compiler!'
     Get a newer cmake, e.g.
         module load cmake                 # or: module load cmake/3.28
         pip install --user --upgrade cmake
         conda install -c conda-forge cmake
     then re-run, or point this script straight at one:
         CMAKE=/path/to/bin/cmake bash $0
     (Or build the GNU way instead:  bash build-now.sh)"
    fi
    cm_lt 3 25 && warn "cmake $CMNUM works with ifx, but >= 3.25 is better tested."
    ;;
esac

if command -v ninja >/dev/null 2>&1; then
  GEN="Ninja"
else
  GEN="Unix Makefiles"
  warn "ninja not found - using Make (slower). Install ninja for a much faster build."
fi

# ---------------------------------------------------------------------------
# 3) subprojects (tblite, toml-f, gfn0, gfnff, lwoniom, pvol, test-drive)
# ---------------------------------------------------------------------------
NO_SUBMODULE="${NO_SUBMODULE:-0}"
if [ -f "$SRC/subprojects/tblite/CMakeLists.txt" ] || [ -f "$SRC/subprojects/toml-f/CMakeLists.txt" ]; then
  NO_SUBMODULE=1
fi
if [ "$NO_SUBMODULE" = "1" ]; then
  say "Using the pre-placed subprojects/ (offline mode)"
  for d in tblite toml-f gfn0 gfnff lwoniom pvol test-drive; do
    if [ -d "$SRC/subprojects/$d" ] && [ ! -f "$SRC/subprojects/$d/CMakeLists.txt" ]; then
      die "subprojects/$d exists but has no CMakeLists.txt - empty submodule placeholder!"
    fi
  done
else
  say "Initialising git subprojects"
  command -v git >/dev/null 2>&1 || die "'git' not found (use NO_SUBMODULE=1 for offline builds)"
  if [ -d "$SRC/.git" ]; then
    git -C "$SRC" submodule update --init --recursive
  else
    die "no .git directory and NO_SUBMODULE=0 - cannot fetch submodules"
  fi
fi

# ---------------------------------------------------------------------------
# 4) linear algebra backend
# ---------------------------------------------------------------------------
CMAKE_EXTRA=()
USE_MKL=0

mkl_dir() { printf '%s' "${MKLROOT:-/opt/intel/oneapi/mkl/latest}"; }

case "$LA" in
  mkl-static|mkl|mkl-rt)
    if [ -z "${MKLROOT:-}" ]; then
      # not fatal: CMake can still find MKL through the compiler driver, but be loud
      warn "MKLROOT is not set - did setvars.sh run? Continuing anyway."
    fi
    MD="$(mkl_dir)"
    case "$LA" in
      mkl-static)
        [ -f "$MD/lib/intel64/libmkl_intel_lp64.a" ] || die \
"static MKL interface library not found under '$MD/lib/intel64'.
     Either source the Intel environment (SETVARS=...) so MKLROOT is right,
     or pick a different backend:
         LA=mkl                 (shared MKL; needs LD_LIBRARY_PATH at run time)
         LA=openblas-static     (reuse \$OPENBLAS_DIR)
         LA=auto                (let CMake decide)"
        CMAKE_EXTRA+=(-DBLA_VENDOR=Intel10_64lp -DBLA_STATIC=ON)
        ;;
      mkl)
        CMAKE_EXTRA+=(-DBLA_VENDOR=Intel10_64lp)
        ;;
      mkl-rt)
        CMAKE_EXTRA+=(-DBLA_VENDOR=Intel10_64_dyn)
        ;;
    esac
    USE_MKL=1
    say "Linear algebra: $LA"
    echo "  MKLROOT    : ${MKLROOT:-<unset>}"
    ;;
  openblas-static)
    [ -f "$OPENBLAS_DIR/lib/libopenblas.a" ] || die \
"libopenblas.a not found in $OPENBLAS_DIR/lib
     Build it first with 'bash build-openblas-static.sh',
     or point OPENBLAS_DIR=<your prefix> at an existing static OpenBLAS."
    # '.a' before '.so' so the static archive wins; '.so' stays in the list so
    # other libraries are still found
    CMAKE_EXTRA+=(
      -DBLA_VENDOR=OpenBLAS
      -DBLA_STATIC=ON
      "-DCMAKE_FIND_LIBRARY_SUFFIXES=.a;.so"
      "-DCMAKE_LIBRARY_PATH=$OPENBLAS_DIR/lib"
    )
    say "Linear algebra: static OpenBLAS ($OPENBLAS_DIR/lib/libopenblas.a)"
    ;;
  openblas)
    CMAKE_EXTRA+=(-DBLA_VENDOR=OpenBLAS)
    say "Linear algebra: shared OpenBLAS"
    ;;
  auto)
    say "Linear algebra: auto (CMake picks whatever it finds)"
    ;;
  *)
    die "unknown LA='$LA' (use mkl-static | mkl | mkl-rt | openblas-static | openblas | auto)"
    ;;
esac

# ---------------------------------------------------------------------------
# 5) Fortran flags
#
#  IMPORTANT: config/CMakeLists.txt sets CMAKE_Fortran_FLAGS itself for Intel
#    dialect = "-g -O2 -r8 -align array64byte -traceback"
#    bounds  = "-check all -fpe0"          (Debug only)
#  so the flags below are only the *extra* layer.  The per-config flags
#  (CMAKE_Fortran_FLAGS_RELEASE) are appended last and therefore win the -O race.
# ---------------------------------------------------------------------------
DEFS=""
# --- WITH_MKL is OPT-IN and deliberately NOT tied to the MKL backend --------
#  Upstream CMake never defines WITH_MKL (only the Meson path does), so every
#  CMake-built CREST has it undefined.  That is a *supported* state:
#    * src/ompmklset.F90 guards with '#ifdef WITH_MKL'  -> code simply omitted;
#      'openblasset' guards with '#ifdef WITH_OPENBLAS' -> no stray OpenBLAS
#      symbol is referenced, so an MKL link stays clean.
#    * tblite/src/tblite/context/type.F90  '#ifndef WITH_MKL / #define
#      WITH_MKL 0' -> self-defaults, so '#if WITH_MKL' is a well-formed 0.
#  Turning it ON defines WITH_MKL=1 (a *value*, so tblite's '#if WITH_MKL'
#  is not left with an empty expression) and activates:
#      MKL_Set_Num_Threads / mkl_set_dynamic  in ompmklset
#      mkl_free_buffers()                     in tblite
#  Use it only on purpose: crest_env.sh already pins MKL to one thread via
#  MKL_NUM_THREADS=1, which is the behaviour the patched thread logic expects.
if [ "$WITH_MKL" = "1" ]; then
  DEFS="$DEFS -DWITH_MKL=1"
  [ "$USE_MKL" = "1" ] || warn "WITH_MKL=1 but the linear-algebra backend is '$LA' - that define only makes sense with MKL"
fi
case "${HEAP_ARRAYS:-}" in
  ""|0) : ;;
  1|all) DEFS="$DEFS -heap-arrays" ;;
  *)   DEFS="$DEFS -heap-arrays $HEAP_ARRAYS" ;;
esac
[ "$FP_PRECISE" = "1" ] && DEFS="$DEFS -fp-model=precise"

BASE_FLAGS="${DEFS# }"
[ -n "${MARCH:-}" ] && BASE_FLAGS="$BASE_FLAGS $MARCH"
[ -n "${EXTRA_FFLAGS:-}" ] && BASE_FLAGS="$BASE_FLAGS $EXTRA_FFLAGS"

if [ "$BUILD_TYPE" = "Debug" ]; then
  CFG_FLAGS="-O0 -g -traceback"
else
  CFG_FLAGS="${OPT:--O3} -DNDEBUG -traceback"
fi

#  CMAKE_EXE_LINKER_FLAGS may only be given ONCE - a second -D for the same
#  variable silently overrides the first one.  Everything therefore has to be
#  accumulated into a single string first.
LINK_FLAGS_ACC=""
if [ "$OMP_STATIC_LINK" = "1" ]; then
  LINK_FLAGS_ACC="$LINK_FLAGS_ACC -qopenmp-link=static"
fi
case "$STATIC" in
  intel|1|runtime)
    LINK_FLAGS_ACC="$LINK_FLAGS_ACC -static-intel -static-libgcc -static-libstdc++"
    ;;
  all|full|2)
    LINK_FLAGS_ACC="$LINK_FLAGS_ACC -static"
    ;;
  0|no|none|"")
    ;;
  *)
    die "unknown STATIC='$STATIC' (use 0 | intel | all)"
    ;;
esac
LINK_FLAGS_ACC="${LINK_FLAGS_ACC# }"
LINK_ARGS=()
if [ -n "$LINK_FLAGS_ACC" ]; then
  LINK_ARGS+=("-DCMAKE_EXE_LINKER_FLAGS=$LINK_FLAGS_ACC")
  case "$LA" in
    mkl-static|openblas-static) : ;;
    auto) : ;;
    *)
      warn "STATIC=$STATIC but the BLAS backend '$LA' is shared - the binary"
      warn "will still need that .so at run time.  Use LA=mkl-static (default)"
      warn "or LA=openblas-static for a self-contained executable."
      ;;
  esac
fi
LINK_FLAGS="${LINK_ARGS[*]:-}"

if [ "$WITH_MKL" = "1" ]; then
  echo "  define     : -DWITH_MKL=1  (CREST calls MKL_Set_Num_Threads/mkl_set_dynamic,"
  echo "               tblite also calls mkl_free_buffers())"
else
  echo "  define     : (none - WITH_MKL left undefined, like the upstream CMake build)"
fi
echo "  base flags : ${BASE_FLAGS:-(none)}"
echo "  ${BUILD_TYPE} flags : $CFG_FLAGS"
[ -n "$LINK_FLAGS" ] && echo "  link flags : $LINK_FLAGS"
case "$STATIC" in
  0|no|none|"") : ;;
  *) echo "  static     : $STATIC   (runtime bound into the binary)" ;;
esac
echo "  jobs       : $JOBS"
echo "  generator  : $GEN"
echo "  prefix     : $PREFIX"

# ---------------------------------------------------------------------------
# 6) configure
# ---------------------------------------------------------------------------
BUILD_DIR="${BUILD_DIR:-$SRC/_build-intel}"
say "Configuring ($BUILD_TYPE, build dir $BUILD_DIR)"
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
CFG_LOG="$BUILD_DIR/_configure.log"

export FC CC
if ! "$CMAKE" -S "$SRC" -B "$BUILD_DIR" -G "$GEN" \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
  -DCMAKE_Fortran_FLAGS="$BASE_FLAGS" \
  -DCMAKE_Fortran_FLAGS_RELEASE="$CFG_FLAGS" \
  -DCMAKE_Fortran_FLAGS_DEBUG="-O0 -g -traceback" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DWITH_OpenMP=TRUE \
  -DWITH_TESTS="$WITH_TESTS" \
  ${LINK_ARGS[@]+"${LINK_ARGS[@]}"} \
  ${CMAKE_EXTRA[@]+"${CMAKE_EXTRA[@]}"} >"$CFG_LOG" 2>&1; then
  echo
  echo "--- the actual cmake errors -------------------------------------------"
  grep -n -i -E 'CMake Error|FATAL_ERROR|Please use an Intel or GNU|identification is unknown|Could NOT find|not found' "$CFG_LOG" | head -n 25 || true
  echo
  echo "--- last 40 lines of $CFG_LOG ---"
  tail -n 40 "$CFG_LOG"
  if grep -q 'Fortran compiler identification is unknown' "$CFG_LOG" 2>/dev/null; then
    echo
    echo ">>> CMake could not identify '$FC' at all ('identification is unknown')."
    echo ">>> With ifx that almost always means cmake is too old: the IntelLLVM"
    echo ">>> compiler id needs cmake >= 3.20, and this run used: $CMV"
    echo ">>> Fix: module load cmake / pip install --user --upgrade cmake, or"
    echo ">>>      CMAKE=/path/to/bin/cmake bash $0"
  fi
  die "configure FAILED (full log: $CFG_LOG)"
fi
tail -n 45 "$CFG_LOG"

# ---------------------------------------------------------------------------
# 6b) did CMake really see an Intel compiler, and did the dialekt flags land?
# ---------------------------------------------------------------------------
CACHE="$BUILD_DIR/CMakeCache.txt"
[ -f "$CACHE" ] || die "configure produced no CMakeCache.txt in $BUILD_DIR - see $CFG_LOG"

# Where does the compiler id actually live?  The authoritative source is the
# generated CMakeFiles/<cmake-version>/CMakeFortranCompiler.cmake, which holds
# 'set(CMAKE_Fortran_COMPILER_ID "IntelLLVM")'.  CMakeCache.txt is only a
# fallback: depending on the CMake version that entry can be absent, or typed
# :INTERNAL/:UNINITIALIZED instead of :STRING.  An earlier version of this
# script parsed only '^CMAKE_Fortran_COMPILER_ID:STRING=' in the cache, found
# nothing, and reported an empty id as if CMake had failed to pick Intel.
FID=""
FVERSION=""
IDFILE=$(find "$BUILD_DIR/CMakeFiles" -maxdepth 2 -name 'CMakeFortranCompiler.cmake' 2>/dev/null | head -n1 || true)
if [ -n "$IDFILE" ]; then
  FID=$(sed -n 's/^set(CMAKE_Fortran_COMPILER_ID[[:space:]]*"\(.*\)")/\1/p' "$IDFILE" | head -n1 || true)
  FVERSION=$(sed -n 's/^set(CMAKE_Fortran_COMPILER_VERSION[[:space:]]*"\(.*\)")/\1/p' "$IDFILE" | head -n1 || true)
fi
[ -n "$FID" ] || FID=$(grep -m1 '^CMAKE_Fortran_COMPILER_ID:' "$CACHE" 2>/dev/null | sed 's/^[^=]*=//' || true)
[ -n "$FVERSION" ] || FVERSION=$(grep -m1 '^CMAKE_Fortran_COMPILER_VERSION:' "$CACHE" 2>/dev/null | sed 's/^[^=]*=//' || true)
FCOMP=$(grep -m1 '^CMAKE_Fortran_COMPILER:' "$CACHE" 2>/dev/null | sed 's/^[^=]*=//' || true)

case "$FID" in
  Intel*)
    echo "  CMake compiler id : $FID $FVERSION  (config/CMakeLists.txt takes the Intel branch)"
    ;;
  '')
    if grep -q 'Fortran compiler identification is unknown' "$CFG_LOG" 2>/dev/null; then
      die "CMake could NOT identify the Fortran compiler ${FCOMP:+('$FCOMP')}.
     It reported 'Fortran compiler identification is unknown' - which is only a
     *warning*, so CMake carried on, but the compiler id stays EMPTY and
     config/CMakeLists.txt cannot take the Intel branch.
     Cause is almost always a cmake too old to know ifx (IntelLLVM needs
     >= 3.20).  This run used: $CMV
     Fix:
         module load cmake          # or: module load cmake/3.28
         pip install --user --upgrade cmake
     or point this script at a newer one:
         CMAKE=/path/to/bin/cmake bash $0"
    fi
    die "could not determine CMAKE_Fortran_COMPILER_ID.
     Looked in: ${IDFILE:-<no CMakeFortranCompiler.cmake found>}
                $CACHE
     Configure claimed success, so inspect $CFG_LOG to see what CMake decided."
    ;;
  *)
    die "CMake identified the Fortran compiler as '$FID', NOT Intel.
     That means config/CMakeLists.txt used the wrong dialect flags (or none).
       compiler CMake picked : ${FCOMP:-unknown}
       FC requested here     : $FC"
    ;;
esac

FLAGFILE=$(find "$BUILD_DIR" -name flags.make -path '*crest-object.dir*' 2>/dev/null | head -n1)
if [ -n "$FLAGFILE" ]; then
  grep -q -- '-r8' "$FLAGFILE" && grep -q -- '-align array64byte' "$FLAGFILE" \
    || die "the Intel dialect flags (-r8 -align array64byte) are missing from $FLAGFILE"
  echo "  Intel dialect     : -r8 -align array64byte present"
  if [ "$WITH_MKL" = "1" ]; then
    grep -q -- '-DWITH_MKL' "$FLAGFILE" \
      && echo "  WITH_MKL          : present" \
      || warn "WITH_MKL=1 was requested but did not reach the compile line"
  fi
else
  warn "could not locate flags.make - skipping the flag self-check"
fi

# ---------------------------------------------------------------------------
# 7) build + install
# ---------------------------------------------------------------------------
say "Building with $JOBS jobs"
if [ "$GEN" = "Ninja" ]; then
  ninja -C "$BUILD_DIR" -j "$JOBS"
else
  make -C "$BUILD_DIR" -j "$JOBS"
fi

say "Installing"
if [ "$GEN" = "Ninja" ]; then
  ninja -C "$BUILD_DIR" install
else
  make -C "$BUILD_DIR" install
fi

CREST="$PREFIX/bin/crest"
[ -x "$CREST" ] || CREST="$BUILD_DIR/crest"

# ---------------------------------------------------------------------------
# 8) smoke test + runtime dependency report
# ---------------------------------------------------------------------------
say "Smoke test"
if [ -x "$CREST" ]; then
  "$CREST" --version || true
else
  die "binary not found - look into $BUILD_DIR"
fi

say "Runtime dependencies of $CREST"
if command -v ldd >/dev/null 2>&1; then
  echo "--- missing libraries ---"
  nf=$(ldd "$CREST" 2>/dev/null | grep 'not found' || true)
  [ -n "$nf" ] && echo "$nf" || echo "  (none)"

  echo "--- BLAS/LAPACK ---"
  bl=$(ldd "$CREST" 2>/dev/null | grep -i -E 'blas|lapack|mkl' || true)
  [ -n "$bl" ] && echo "$bl" || echo "  (none - statically linked)"

  echo "--- OpenMP runtime ---"
  om=$(ldd "$CREST" 2>/dev/null | grep -i -E 'iomp5|gomp' || true)
  [ -n "$om" ] && echo "$om" || echo "  (none - statically linked)"
  case "$om" in
    *iomp5*) echo "  note: libiomp5.so must exist on the compute node. If it does not,"
             echo "        either add the Intel compiler lib dir to LD_LIBRARY_PATH,"
             echo "        or rebuild with OMP_STATIC_LINK=1 / STATIC=intel" ;;
  esac

  echo "--- Intel runtime (libifcore/libimf/libsvml/libirc) ---"
  ir=$(ldd "$CREST" 2>/dev/null \
       | grep -i -E 'libifcore|libifport|libimf|libsvml|libirc|libintlc' || true)
  [ -n "$ir" ] && echo "$ir" || echo "  (none - statically linked)"
  if [ -z "$ir" ] && [ -z "$om" ] && [ -z "$bl" ]; then
    echo "  => no compiler/BLAS/OpenMP runtime left: this binary runs on a node"
    echo "     without the Intel compiler installed."
  fi
fi

cat <<EOF

-----------------------------------------------------------------------------
 Patched CREST is built with the Intel toolchain.

   compiler : $FVER
   BLAS     : $LA
   binary   : $CREST
   build dir: $BUILD_DIR

 Run it:

   source $SRC/crest_env.sh 48
   cd <your working directory>          # preferably NOT under /mnt/...
   crest your-system-orca.toml -T 48 --quick --noreftopo

 Intel-specific notes
   * 'crest --version' prints the compiler that was used, so the metadata
     block doubles as a build receipt.
   * crest_env.sh already exports MKL_NUM_THREADS=1 / OPENBLAS_NUM_THREADS=1,
     and CREST re-exports the per-job budget to every child (ompenvset).
   * If you see "undefined reference to mkl_*" or "iomp5 not found" in a
     threaded run, rebuild with OMP_STATIC_LINK=1.
   * For a node-specific binary add MARCH="-xHost"; the default build is
     portable across the cluster.
-----------------------------------------------------------------------------
EOF
