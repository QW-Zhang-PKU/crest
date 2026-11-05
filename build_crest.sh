#!/usr/bin/env bash
# =============================================================================
#  CREST 3.0.2 build helper  (GNU / CMake)
#  Usage:
#     bash build_crest.sh [SOURCE_DIR] [INSTALL_PREFIX]
#  e.g.
#     bash build_crest.sh $HOME/crest-program $HOME/crest-dist
#
#  Useful environment switches:
#     WITH_TESTS=ON       build the unit tests as well (needs test-drive)
#     NO_SUBMODULE=1      skip `git submodule update --init` (use pre-placed sources)
#     JOBS=16             parallel compile jobs (default: nproc)
#
#  Sizing the build (see section "1b" below):
#     FULL=1              build EVERYTHING (the upstream default feature set)
#     WITH_TBLITE=ON      also re-enable a single optional dependency, e.g.
#     WITH_TBLITE=ON WITH_GFNFF=ON bash build_crest.sh
#
#     Defaults in this tree are SLIM: the in-process engines that a pure
#     external-calculator (xtb / g-xTB / ORCA) conformer search never touches
#     are switched OFF, because tblite alone is ~6200 source files and
#     dominates the build time.
# =============================================================================
set -euo pipefail

SRC="${1:-$HOME/crest-program}"
PREFIX="${2:-$HOME/crest-dist}"
JOBS="${JOBS:-$(nproc)}"
WITH_TESTS="${WITH_TESTS:-OFF}"
NO_SUBMODULE="${NO_SUBMODULE:-0}"

say() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 0) locate the real source root (the dir that holds CMakeLists.txt + src/)
# ---------------------------------------------------------------------------
say "Locating source root under: $SRC"
[ -d "$SRC" ] || die "source directory '$SRC' does not exist"
if [ ! -f "$SRC/CMakeLists.txt" ]; then
  CAND=$(find "$SRC" -maxdepth 3 -name CMakeLists.txt -printf '%h\n' 2>/dev/null | head -n 1 || true)
  [ -n "$CAND" ] && [ -f "$CAND/src/CMakeLists.txt" ] || die "no CMakeLists.txt with src/ found under $SRC"
  SRC="$CAND"
fi
say "Source root: $SRC"

# ---------------------------------------------------------------------------
# 1) toolchain check
# ---------------------------------------------------------------------------
say "Checking toolchain"
for t in gfortran gcc cmake; do
  command -v "$t" >/dev/null 2>&1 || die "'$t' not found in PATH"
done
if [ "$NO_SUBMODULE" != "1" ]; then
  command -v git >/dev/null 2>&1 || die "'git' not found in PATH (use NO_SUBMODULE=1 for offline builds)"
fi
GFV=$(gfortran -dumpversion | cut -d. -f1)
[ "$GFV" -ge 10 ] || die "gfortran >= 10 required (found $GFV)"
command -v ninja >/dev/null 2>&1 && GEN="Ninja" || GEN="Unix Makefiles"
echo "  gfortran   : $(gfortran --version | head -n1)"
echo "  cmake      : $(cmake --version | head -n1)"
echo "  generator  : $GEN"
echo "  jobs       : $JOBS"

# ---------------------------------------------------------------------------
# 1b) optional dependencies (feature switches)
#
#  A conformer search that uses EXTERNAL calculators only - which is the case
#  for the recommended `xtb` (gfn2), `g-xTB` and `ORCA` setups - never calls
#  the in-process engines.  They are therefore OFF by default here:
#
#    WITH_TBLITE=OFF    in-process GFN2-xTB through the tblite library.
#                       ~6200 source files / ~45 MB: this single subproject is
#                       what makes a full CREST build slow.  The external
#                       `xtb` binary covers gfn2 and g-xTB instead.
#    WITH_GFNFF=OFF     in-process GFN-FF.
#    WITH_LIBPVOL=OFF   external-pressure (libpvol) calculations.
#    WITH_LWONIOM=OFF   ONIOM.
#
#  KEPT ON, because they are genuinely required:
#
#    WITH_GFN0=ON       the flexibility/topology measure at the START of every
#                       conformer search uses the internal GFN0 engine, as does
#                       the SHAKE/topology setup and the protonation tools.
#    WITH_TOMLF=ON      the TOML input parser.
#
#  Each switch is a normal CMake cache option (config/CMakeLists.txt:33-39), so
#  it can be re-enabled individually without touching this script:
#
#      WITH_TBLITE=ON bash build_crest.sh              # one dependency back
#      FULL=1 bash build_crest.sh                      # the upstream default
#
#  NOTE: with a dependency OFF its source directory is not needed at all (the
#  code is #ifdef-guarded and compiles to a stub), so a SLIM source tree may
#  simply omit subprojects/tblite, gfnff, lwoniom and pvol.
# ---------------------------------------------------------------------------
if [ "${FULL:-0}" = "1" ]; then
  WITH_TBLITE=ON; WITH_GFNFF=ON; WITH_LIBPVOL=ON; WITH_LWONIOM=ON
  WITH_GFN0=ON;   WITH_TOMLF=ON
else
  WITH_TBLITE="${WITH_TBLITE:-OFF}"
  WITH_GFNFF="${WITH_GFNFF:-OFF}"
  WITH_LIBPVOL="${WITH_LIBPVOL:-OFF}"
  WITH_LWONIOM="${WITH_LWONIOM:-OFF}"
  WITH_GFN0="${WITH_GFN0:-ON}"
  WITH_TOMLF="${WITH_TOMLF:-ON}"
fi

say "Optional dependencies"
for v in WITH_TOMLF WITH_GFN0 WITH_TBLITE WITH_GFNFF WITH_LIBPVOL WITH_LWONIOM; do
  eval "printf '  %-14s : %s\n' \"$v\" \"\${$v}\""
done
[ "$WITH_TOMLF" = "ON" ] || die "WITH_TOMLF=OFF is not supported: CREST needs the TOML parser"

# ---------------------------------------------------------------------------
# 2) subprojects (tblite, toml-f, gfn0, gfnff, lwoniom, libpvol, test-drive)
# ---------------------------------------------------------------------------
# offline: subprojects already present as plain source directories?
# (do not probe tblite specifically - a slim tree may omit it on purpose)
if [ -f "$SRC/subprojects/toml-f/CMakeLists.txt" ] || \
   [ -f "$SRC/subprojects/tblite/CMakeLists.txt" ] || \
   [ -f "$SRC/subprojects/gfn0/CMakeLists.txt" ]; then
  say "Found pre-placed subprojects/ -> offline mode"
  NO_SUBMODULE=1
fi

if [ "$NO_SUBMODULE" = "1" ]; then
  say "Skipping submodule initialisation (NO_SUBMODULE=1)"
  #>--- Only the dependencies that are actually ENABLED have to be present.  A
  #> disabled one is never handed to find_package(), so its source directory may
  #> be missing entirely - that is what makes a slim source tree possible.
  needed_deps=()
  if [ "$WITH_TOMLF"   = "ON" ]; then needed_deps+=(toml-f);   fi
  if [ "$WITH_TBLITE"  = "ON" ]; then needed_deps+=(tblite);   fi
  if [ "$WITH_GFN0"    = "ON" ]; then needed_deps+=(gfn0);     fi
  if [ "$WITH_GFNFF"   = "ON" ]; then needed_deps+=(gfnff);    fi
  if [ "$WITH_LWONIOM" = "ON" ]; then needed_deps+=(lwoniom);  fi
  if [ "$WITH_LIBPVOL" = "ON" ]; then needed_deps+=(pvol);     fi
  if [ "$WITH_TESTS"   = "ON" ]; then needed_deps+=(test-drive); fi
  echo "  required subprojects: ${needed_deps[*]:-<none>}"
  for d in "${needed_deps[@]+"${needed_deps[@]}"}"; do
    if [ ! -d "$SRC/subprojects/$d" ]; then
      die "subprojects/$d is ENABLED but its source directory is missing - provide it or set its WITH_*=OFF switch"
    fi
    if [ ! -f "$SRC/subprojects/$d/CMakeLists.txt" ]; then
      die "subprojects/$d exists but has no CMakeLists.txt - empty submodule placeholder!"
    fi
  done
else
  say "Initialising git subprojects"
  if [ -d "$SRC/.git" ]; then
    git -C "$SRC" submodule update --init --recursive
  elif [ -f "$SRC/.gitmodules" ]; then
    echo "  no .git directory -> cloning submodules by hand"
    while read -r path url branch; do
      [ -n "${path:-}" ] || continue
      if [ ! -d "$SRC/$path" ] || [ -z "$(ls -A "$SRC/$path" 2>/dev/null)" ]; then
        echo "  cloning $url -> $path"
        rm -rf "$SRC/$path"
        if [ -n "${branch:-}" ]; then
          git clone --depth 1 -b "$branch" "$url" "$SRC/$path"
        else
          git clone --depth 1 "$url" "$SRC/$path"
        fi
      fi
    done < <(awk '/^\[submodule/{p="";u="";b=""} \
                  /path *=/{p=$3} /url *=/{u=$3} /branch *=/{b=$3} \
                  /^$/{if(p!="")print p, u, b}' "$SRC/.gitmodules")
  else
    echo "  no .gitmodules -> relying on CMake 'fetch' (needs internet)"
  fi
fi

# ---------------------------------------------------------------------------
# 3a) static-linking options
#
#   STATIC_BLAS=1       只把 BLAS/LAPACK(OpenBLAS) 静态链进去，
#                       libc / libgomp / libgfortran 仍用动态。
#                       -> 推荐：解决 "libopenblas.so.0 not found"，
#                          又不碰 glibc，几乎不会编挂。
#   STATICBUILD=ON      CREST 自带的全静态开关（CMakeLists.txt:48），
#                       连 glibc 一起静态，需要 glibc-static 包，
#                       且静态二进制 + OpenMP 会有警告，能不用就不用。
#   OPENBLAS_DIR=...    静态 OpenBLAS 的安装前缀（内含 lib/libopenblas.a）
# ---------------------------------------------------------------------------
STATIC_BLAS="${STATIC_BLAS:-0}"
STATICBUILD="${STATICBUILD:-OFF}"
OPENBLAS_DIR="${OPENBLAS_DIR:-$HOME/opt/openblas-static}"

CMAKE_EXTRA=()
if [ "$STATIC_BLAS" = "1" ]; then
  say "Static BLAS/LAPACK requested"
  if [ -f "$OPENBLAS_DIR/lib/libopenblas.a" ]; then
    echo "  static OpenBLAS: $OPENBLAS_DIR/lib/libopenblas.a"
    echo "  size           : $(du -h "$OPENBLAS_DIR/lib/libopenblas.a" | cut -f1)"
  else
    die "libopenblas.a not found in $OPENBLAS_DIR/lib — 先跑 bash build-openblas-static.sh，
     或者用 OPENBLAS_DIR=<你的前缀> 指到已有的静态 OpenBLAS 安装上。"
  fi
  # ".a" 排在前面 => 同目录下优先选静态库；".so" 保留，避免其他库找不到
  CMAKE_EXTRA+=(
    -DBLA_VENDOR=OpenBLAS
    -DBLA_STATIC=ON
    "-DCMAKE_FIND_LIBRARY_SUFFIXES=.a;.so"
    "-DCMAKE_LIBRARY_PATH=$OPENBLAS_DIR/lib"
  )
  # A pthread-flavoured libopenblas.a occasionally needs its thread library
  # spelled out (symptom: "undefined reference to `pthread_create'").  Not set
  # by default because on glibc >= 2.34 libpthread is an empty stub and CMake
  # usually appends Threads::Threads itself.
  if [ -n "${STATIC_BLAS_EXTRA_LIBS:-}" ]; then
    echo "  extra link libs: $STATIC_BLAS_EXTRA_LIBS"
    CMAKE_EXTRA+=("-DCMAKE_EXE_LINKER_FLAGS=$STATIC_BLAS_EXTRA_LIBS")
  fi
fi
if [ "$STATICBUILD" = "ON" ]; then
  say "Fully static build (STATICBUILD=ON) — needs glibc-static etc."
  CMAKE_EXTRA+=(-DSTATICBUILD=ON)
fi

# ---------------------------------------------------------------------------
# 3) configure
# ---------------------------------------------------------------------------
say "Configuring (build type Release, install to $PREFIX)"
export FC=gfortran CC=gcc
rm -rf "$SRC/_build"
cmake -S "$SRC" -B "$SRC/_build" -G "$GEN" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS_RELEASE="-O3 -DNDEBUG -fbacktrace -ffree-line-length-none" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DWITH_OpenMP=TRUE \
  -DWITH_TESTS="$WITH_TESTS" \
  -DWITH_TOMLF="$WITH_TOMLF" \
  -DWITH_GFN0="$WITH_GFN0" \
  -DWITH_TBLITE="$WITH_TBLITE" \
  -DWITH_GFNFF="$WITH_GFNFF" \
  -DWITH_LIBPVOL="$WITH_LIBPVOL" \
  -DWITH_LWONIOM="$WITH_LWONIOM" \
  ${CMAKE_EXTRA[@]+"${CMAKE_EXTRA[@]}"} 2>&1 | tail -n 40

# ---------------------------------------------------------------------------
# 4) build
# ---------------------------------------------------------------------------
say "Building with $JOBS jobs"
if [ "$GEN" = "Ninja" ]; then
  ninja -C "$SRC/_build" -j "$JOBS"
else
  make -C "$SRC/_build" -j "$JOBS"
fi

# ---------------------------------------------------------------------------
# 5) install
# ---------------------------------------------------------------------------
say "Installing"
if [ "$GEN" = "Ninja" ]; then
  ninja -C "$SRC/_build" install
else
  make -C "$SRC/_build" install
fi

# ---------------------------------------------------------------------------
# 6) smoke test
# ---------------------------------------------------------------------------
CREST="$PREFIX/bin/crest"
[ -x "$CREST" ] || CREST="$SRC/_build/crest"
say "Smoke test"
if [ -x "$CREST" ]; then
  "$CREST" --version || true
else
  echo "  binary not found, look into $SRC/_build"
fi

# --- 运行期依赖自检：静态编译就是为了这条输出里没有 openblas ----------
if [ -x "$CREST" ]; then
  say "Runtime dependencies"
  if command -v ldd >/dev/null 2>&1; then
    echo "--- 缺失的库 ---"
    nf=$(ldd "$CREST" 2>/dev/null | grep 'not found')
    [ -n "$nf" ] && echo "$nf" || echo "  (无)"
    echo "--- BLAS/LAPACK 相关 ---"
    bl=$(ldd "$CREST" 2>/dev/null | grep -i -E 'blas|lapack')
    if [ -n "$bl" ]; then
      echo "$bl"
      echo "  ^ 仍是动态 BLAS：若目标节点没有该库，请把 .so 拷到 crest 的 ../lib"
    else
      echo "  (无 —— BLAS/LAPACK 已静态链入)"
    fi
  fi
fi

cat <<EOF

-----------------------------------------------------------------------------
 Done.  Binary: $CREST

 Add it to your PATH:
   echo 'export PATH=$PREFIX/bin:\$PATH' >> ~/.bashrc
   source ~/.bashrc

 IMPORTANT before running (OpenBLAS backend spawns extra thread pools):
   export OMP_NUM_THREADS=40
   export OPENBLAS_NUM_THREADS=1
   export MKL_NUM_THREADS=1

 Check which calculator is actually used (tblite API vs. external xtb binary)
 in the "calculation level" block that CREST prints at startup.
-----------------------------------------------------------------------------
EOF
