#!/usr/bin/env bash
# =============================================================================
#  check-orca-threads.sh  --  CREST + ORCA: does the CPU load match -T ?
#
#  WHY THIS SCRIPT EXISTS
#  ----------------------
#  On a 48-core / 96-thread box, `crest ... -T 48` looks like "every logical
#  core is busy" in htop.  That is expected and harmless IF there are really
#  only ~48 runnable threads (the scheduler spreads them over 96 logical CPUs,
#  each physical core then shows ~50 % on both of its hyperthreads).  It IS a
#  problem if there are ~96 runnable threads: that is 2x oversubscription.
#
#  So the question is never "how many cores look busy" but
#      HOW MANY EXECUTION UNITS ARE ACTUALLY RUNNING?
#  This script answers that, plus the ORCA-specific causes:
#
#    * CREST injects          %pal nprocs <cores_per_job>   into ORCA.in and
#      appends                --oversubscribe               to the orca command
#      (src/calculator/orca_sc.f90 ORCA_setup).
#    * --oversubscribe makes OpenMPI DROP ITS CPU BINDING.  48 unbound MPI ranks
#      then wander over all 96 logical CPUs: cache locality is lost, two ranks
#      can land on the two hyperthreads of the SAME physical core, and memory
#      traffic crosses NUMA nodes.  Nothing is oversubscribed on paper, yet the
#      run is slower than it should be.
#
#  USAGE
#     # terminal 1
#     crest coord your-system-orca.toml -T 48 ...
#     # terminal 2
#     watch -n 5 ./check-orca-threads.sh
#
#  WHAT GOOD LOOKS LIKE (48 physical cores, -T 48, CREST_CORES_PER_JOB=4)
#     busy (equiv. cores) ~= 48
#     orca ranks ~= 48, threads per rank ~= 1-2
#     hyperthread pairs with BOTH siblings >70 %: 0
# =============================================================================
set -uo pipefail

tmp=$(mktemp -d 2>/dev/null || echo /tmp/cot.$$)
mkdir -p "$tmp" 2>/dev/null
trap 'rm -rf "$tmp" 2>/dev/null' EXIT

NCPU=$(nproc 2>/dev/null || echo 0)

# ---------------------------------------------------------------------------
echo "=== 1) CPU topology ==================================================="
TPC=$(lscpu 2>/dev/null | awk -F: '/^Thread\(s\) per core/{gsub(/ /,"",$2); print $2; exit}')
CPS=$(lscpu 2>/dev/null | awk -F: '/^Core\(s\) per socket/{gsub(/ /,"",$2); print $2; exit}')
SOCK=$(lscpu 2>/dev/null | awk -F: '/^Socket\(s\)/{gsub(/ /,"",$2); print $2; exit}')
if [ -n "${TPC:-}" ] && [ -n "${CPS:-}" ] && [ -n "${SOCK:-}" ]; then
  NPHYS=$(( CPS * SOCK ))
else
  NPHYS=$(( NCPU / ${TPC:-1} ))
fi
echo "  logical CPUs    : $NCPU"
echo "  physical cores  : $NPHYS   ($TPC threads/core${SOCK:+ , $SOCK socket(s)})"

# ---------------------------------------------------------------------------
echo
echo "=== 2) how much CPU is REALLY being used (1 s sample) ================="
read_stat() {
  # user nice system idle iowait irq softirq steal  -> busy excludes idle+iowait
  awk '/^cpu[0-9]/{ busy=$2+$3+$4+$7+$8+$9; tot=$2+$3+$4+$5+$6+$7+$8+$9
                    print $1, busy, tot }' /proc/stat 2>/dev/null
}
read_stat > "$tmp/a"; sleep 1; read_stat > "$tmp/b"
awk 'NR==FNR{b[$1]=$2; t[$1]=$3; next}
     ($1 in b){ db=$2-b[$1]; dt=$3-t[$1]; printf "%s %.1f\n", $1, (dt>0)?100*db/dt:0 }
    ' "$tmp/a" "$tmp/b" > "$tmp/pct"
if [ -s "$tmp/pct" ]; then
  BUSYC=$(awk '{s+=$2} END{printf "%.1f", s/100}' "$tmp/pct")
  PCT=$(awk -v n="$NCPU" '{s+=$2} END{printf "%.1f", (n>0)?s/n:0}' "$tmp/pct")
  echo "  busy, in units of a full core : $BUSYC  (of $NCPU logical / $NPHYS physical)"
  echo "  total CPU utilisation         : ${PCT} %"
else
  BUSYC=""
  echo "  (could not sample /proc/stat)"
fi

# ---------------------------------------------------------------------------
echo
echo "=== 3) ORCA processes =================================================="
mapfile -t PIDS < <(pgrep -x orca 2>/dev/null)
if [ "${#PIDS[@]}" -eq 0 ]; then
  mapfile -t PIDS < <(pgrep -f 'orca.*ORCA\.in' 2>/dev/null)
fi
if [ "${#PIDS[@]}" -eq 0 ]; then
  echo "  no orca process running right now"
  echo "  (is a calculation in flight?  try again during the MTD/optimisation phase)"
  ORCA_N=0; ORCA_T=0
else
  printf '  %-8s %-7s %-9s %s\n' "PID" "THREADS" "%CPU" "CWD / nprocs"
  ORCA_T=0
  for p in "${PIDS[@]}"; do
    [ -r "/proc/$p" ] || continue
    nt=$(ls "/proc/$p/task" 2>/dev/null | wc -l); nt=${nt:-0}
    ORCA_T=$(( ORCA_T + nt ))
    cpu=$(ps -o pcpu= -p "$p" 2>/dev/null | tr -d ' ')
    cwd=$(readlink "/proc/$p/cwd" 2>/dev/null | xargs -r basename 2>/dev/null)
    np=$(grep -iA2 '^[[:space:]]*%pal' "/proc/$p/cwd/ORCA.in" 2>/dev/null \
         | grep -io 'nprocs[[:space:]]*[0-9]*' | head -1 | grep -o '[0-9]*')
    printf '  %-8s %-7s %-9s %s%s\n' "$p" "$nt" "${cpu:-?}" "${cwd:-?}" \
           "${np:+  (nprocs $np)}"
  done
  ORCA_N=${#PIDS[@]}
  echo
  echo "  orca ranks total : $ORCA_N"
  echo "  their OS threads : $ORCA_T"
fi

# ---------------------------------------------------------------------------
echo
echo "=== 4) is the orca command line oversubscribing? ======================"
shpid=$(pgrep -f 'sh -c.*orca' 2>/dev/null | head -1)
if [ -n "$shpid" ] && [ -r "/proc/$shpid/cmdline" ]; then
  cmd=$(tr '\0' ' ' < "/proc/$shpid/cmdline" 2>/dev/null)
  echo "  ${cmd:0:300}"
  case "$cmd" in
    *--oversubscribe*) echo "  -> carries --oversubscribe : OpenMPI CPU BINDING IS OFF" ;;
    *)                 echo "  -> no --oversubscribe" ;;
  esac
else
  echo "  (no 'sh -c ... orca' wrapper alive at this instant - re-run during a step)"
fi

# ---------------------------------------------------------------------------
echo
echo "=== 5) hyperthread siblings: are both halves of a core busy? =========="
if [ ! -s "$tmp/pct" ]; then
  echo "  (no per-cpu data)"
else
  pairs_both=0; pairs_any=0; pairs_total=0
  declare -A seen=()
  while read -r name pct; do
    c=${name#cpu}
    sib=$(cat "/sys/devices/system/cpu/cpu$c/topology/thread_siblings_list" 2>/dev/null)
    [ -n "$sib" ] || continue
    key=$(echo "$sib" | tr ',' '\n' | sort -n | paste -sd, -)
    [ -n "${seen[$key]:-}" ] && continue
    seen[$key]=1
    pairs_total=$(( pairs_total + 1 ))
    # sum of the busy% of every sibling in this physical core
    s=$(echo "$sib" | tr ',' '\n' | while read -r s2; do
          awk -v n="cpu$s2" '$1==n{print $2}' "$tmp/pct"
        done | awk '{t+=$1} END{printf "%.1f", t}')
    s_int=${s%.*}
    [ "$s_int" -ge 120 ] 2>/dev/null && pairs_both=$(( pairs_both + 1 ))
    [ "$s_int" -ge 40 ]  2>/dev/null && pairs_any=$(( pairs_any + 1 ))
  done < "$tmp/pct"
  if [ "$pairs_total" -eq 0 ]; then
    echo "  (this kernel does not expose thread_siblings_list - skip)"
    echo "  use instead:  mpstat -P ALL 1 1   and compare cpu N with cpu N+$NPHYS"
  else
    echo "  physical cores examined                    : $pairs_total"
    echo "  cores whose sibling sum is >= 120 %        : $pairs_both   <-- HT contention"
    echo "  cores doing any work (>= 40 %)             : $pairs_any"
  fi
fi

# ---------------------------------------------------------------------------
echo
echo "=== verdict ==========================================================="
TH=$(( NPHYS * 12 / 10 ))
if [ "$ORCA_N" -eq 0 ]; then
  echo "  no ORCA running - nothing to judge yet."
elif [ "$ORCA_T" -gt "$TH" ]; then
  echo "  OVERSUBSCRIBED: $ORCA_T OS threads for $NPHYS physical cores."
  echo "  Every core is running ~2 compute threads.  Expect 20-50 % slowdown."
  echo "  Fix: lower -T, or raise CREST_CORES_PER_JOB so fewer, fatter jobs run."
elif [ "$ORCA_N" -gt "$TH" ] 2>/dev/null; then
  echo "  Many MPI ranks ($ORCA_N) for $NPHYS cores - check the %pal nprocs value"
  echo "  CREST injected and the number of parallel CREST jobs."
else
  echo "  Thread count looks sane (~$ORCA_T units for $NPHYS cores)."
  echo "  If htop shows all $NCPU logical CPUs busy that is the normal"
  echo "  hyperthread spread, NOT oversubscription."
  echo "  The remaining risk is the missing CPU BINDING (see section 4):"
  echo "  try  export CREST_ORCA_MPI_ARGS='--bind-to core'"
fi
echo
