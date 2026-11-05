#!/usr/bin/env bash
# =============================================================================
#  build-openblas-static.sh — 编译一份【静态 + 通用 CPU】的 OpenBLAS
#
#  为什么需要它：
#    CREST 在编译机上动态链了 libopenblas.so.0，而计算节点(node11)上没有，
#    于是 crest 起不来：
#        crest: error while loading shared libraries: libopenblas.so.0
#    把 OpenBLAS 静态链进 crest 是最彻底的解法。
#
#  两个关键开关：
#    NO_SHARED=1     只产出 libopenblas.a，不产出 .so（避免又链到动态库）
#    DYNAMIC_ARCH=1  【重要】运行时按 CPU 派发内核。
#                    编译机若是 AVX512 而 node11 只有 AVX2，
#                    没有这个选项会直接 illegal instruction。
#
#  用法:
#     bash build-openblas-static.sh                     # 装到 ~/opt/openblas-static
#     PREFIX=/opt/software/openblas-static bash build-openblas-static.sh
#     JOBS=32 bash build-openblas-static.sh
#     THREADS=openmp bash build-openblas-static.sh      # 默认 pthreads
#
#  之后:
#     OPENBLAS_DIR=$PREFIX STATIC_BLAS=1 bash build-now.sh
# =============================================================================
set -euo pipefail

PREFIX="${PREFIX:-$HOME/opt/openblas-static}"
JOBS="${JOBS:-$(nproc)}"
THREADS="${THREADS:-pthreads}"          # pthreads | openmp
VER="${VER:-v0.3.28}"

say() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

command -v gfortran >/dev/null 2>&1 || die "gfortran not found"
command -v make     >/dev/null 2>&1 || die "make not found"
command -v git      >/dev/null 2>&1 || die "git not found (OpenBLAS 源码需要 git 拉取)"

echo "prefix    : $PREFIX"
echo "jobs      : $JOBS"
echo "threads   : $THREADS"
echo "version   : $VER"

# ---------------------------------------------------------------------------
# 1) 源码
# ---------------------------------------------------------------------------
BUILD="$HOME/openblas-build"
mkdir -p "$BUILD"
if [ ! -d "$BUILD/OpenBLAS/.git" ]; then
  say "Cloning OpenBLAS $VER"
  rm -rf "$BUILD/OpenBLAS"
  git clone --depth 1 -b "$VER" https://github.com/OpenMathLib/OpenBLAS "$BUILD/OpenBLAS"
else
  say "OpenBLAS source already present"
fi
cd "$BUILD/OpenBLAS"

# ---------------------------------------------------------------------------
# 2) 线程模型
# ---------------------------------------------------------------------------
case "$THREADS" in
  pthreads|pthread) MAKE_THREAD="USE_OPENMP=0" ; NOTE="pthreads（OpenBLAS 自带线程池）" ;;
  openmp|omp)       MAKE_THREAD="USE_OPENMP=1" ; NOTE="OpenMP（与 crest 共用 libgomp 线程池）" ;;
  *) die "THREADS must be 'pthreads' or 'openmp'" ;;
esac
echo "thread model: $NOTE"

# ---------------------------------------------------------------------------
# 3) 编译（只静态 + 通用 CPU 内核）
# ---------------------------------------------------------------------------
say "Building static OpenBLAS (this takes a few minutes)"
make -j "$JOBS" \
     $MAKE_THREAD \
     NO_SHARED=1 \
     DYNAMIC_ARCH=1 \
     NO_AFFINITY=1 \
     USE_THREAD=1 \
     NUM_THREADS=64

# ---------------------------------------------------------------------------
# 4) 安装
# ---------------------------------------------------------------------------
say "Installing to $PREFIX"
make install PREFIX="$PREFIX" NO_SHARED=1

# ---------------------------------------------------------------------------
# 5) 检查
# ---------------------------------------------------------------------------
say "Result"
A="$PREFIX/lib/libopenblas.a"
if [ -f "$A" ]; then
  echo "  OK: $A"
  echo "  size: $(du -h "$A" | cut -f1)"
else
  die "静态库没有生成，检查上面的编译日志"
fi
if ls "$PREFIX"/lib/libopenblas.so* >/dev/null 2>&1; then
  echo "  警告: 目录下仍有 .so —— 编译时请把 $PREFIX/lib 放在搜索路径最前面，"
  echo "        或直接把 .so 删掉，否则 CMake 可能又链到动态库："
  ls -1 "$PREFIX"/lib/libopenblas.so*
else
  echo "  (目录里没有 .so —— 正确)"
fi

cat <<EOF

-----------------------------------------------------------------------------
 静态 OpenBLAS 就绪: $PREFIX

 接着重编 CREST（只静态化 BLAS/LAPACK，其余保持动态）:

     OPENBLAS_DIR=$PREFIX STATIC_BLAS=1 bash build-now.sh

 编完确认（下面应该没有 blas/lapack，也没有 not found）:

     ldd <prefix>/bin/crest | grep -i -E 'blas|not found'

 若要连 glibc 一起静态（不推荐，需要 glibc-static 且 OpenMP 有警告）:

     STATICBUILD=ON bash build-now.sh
-----------------------------------------------------------------------------
EOF
