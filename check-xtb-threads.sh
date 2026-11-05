#!/usr/bin/env bash
# =============================================================================
#  check-xtb-threads.sh
#
#  Show, for every running xtb process, how many OS threads it actually has and
#  which OMP/BLAS thread budget it inherited.
#
#  WHY: in CREST's "external calculator" (subprocess) mode every MD step starts a
#  fresh xtb process.  CREST hands the per-job core count to that child through
#  OMP_NUM_THREADS / MKL_NUM_THREADS / OPENBLAS_NUM_THREADS (see ompenvset() in
#  src/ompmklset.F90).  If that propagation is broken the child silently keeps
#  whatever the calling shell exported - typically OMP_NUM_THREADS=40 with
#  OPENBLAS_NUM_THREADS=1 - which over-subscribes the OpenMP layer and
#  serialises the BLAS layer underneath.
#
#  USAGE
#     # terminal 1
#     source ~/crest_env.sh 40
#     crest your-system-xtb.toml -T 40 --quick --noreftopo
#
#     # terminal 2
#     watch -n 2 ./check-xtb-threads.sh
#
#  WHAT GOOD LOOKS LIKE (6 MTDs, -T 40 -> cores_per_job = 6)
#     OMP_NUM_THREADS=6  MKL_NUM_THREADS=6  OPENBLAS_NUM_THREADS=6
#     OMP_MAX_ACTIVE_LEVELS=1
#     task count ~= 6-8 per xtb process
#
#  WHAT IS BROKEN
#     OMP_NUM_THREADS=40 (inherited from the shell) or OPENBLAS_NUM_THREADS=1
#     alongside a task count of 1-2  ->  the child thread budget is not being
#     propagated; rebuild CREST with crest-omp-env-fix.patch applied.
# =============================================================================
set -uo pipefail

pids=$(pgrep -x xtb 2>/dev/null || true)
if [ -z "$pids" ]; then
  echo "no xtb process running right now (did the trial/batch finish?)"
  exit 0
fi

printf '%s\n' "----------------------------------------------------------------"
printf '%-8s %-6s  %s\n' "PID" "THREADS" "ENVIRONMENT"
printf '%s\n' "----------------------------------------------------------------"

for p in $pids; do
  nthreads=$(ls "/proc/$p/task" 2>/dev/null | wc -l)
  envs=$(tr '\0' '\n' < "/proc/$p/environ" 2>/dev/null \
         | grep -E '^(OMP|OPENBLAS|MKL|VECLIB|NUMEXPR)_(NUM_THREADS|MAX_ACTIVE_LEVELS)=' \
         | sort | paste -sd' ' -)
  printf '%-8s %-6s  %s\n' "$p" "$nthreads" "${envs:-(unreadable)}"
  # the interesting bits: which binary, and in which calcspace it is running
  cwd=$(readlink "/proc/$p/cwd" 2>/dev/null || echo '?')
  printf '%-8s %-6s  cwd=%s\n' "" "" "$cwd"
done

echo
echo "totals:"
for p in $pids; do ls "/proc/$p/task" 2>/dev/null | wc -l; done \
  | awk '{s+=$1; n++} END {printf "  xtb processes: %d   OS threads total: %d\n", n, s}'
echo "  cores available: $(nproc)"
