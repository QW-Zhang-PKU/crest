#!/usr/bin/env bash
# =============================================================================
#  bench-nmtd.sh — find the MTD-count that actually uses your cores
#
#  WHY: in external-calculator (subprocess) mode CREST splits the thread budget
#  as
#        parallel_jobs = nmetadyn                    (one CREST thread per MTD)
#        cores_per_job = floor(-T / nmetadyn)        (-> OMP_NUM_THREADS)
#  For a batch of 6 MTDs on -T 40 that is 6 MTDs x 6 threads.  But xtb scales
#  badly with threads while independent xtb processes scale perfectly, so
#  6 x 6 wastes most of the machine.  Shrinking cores_per_job to 1 and running
#  more MTDs in parallel is typically 2-4x faster for the SAME amount of MD.
#
#  USAGE
#      bash bench-nmtd.sh                     # sweep 6 12 20 40 on -T 40
#      bash bench-nmtd.sh 40 6 20 40          # explicit thread count + list
#      TOML=other.toml bash bench-nmtd.sh
#
#  Each value re-runs the whole CREST job (including the ~1.5 min trial MTD),
#  so budget ~5 min per entry — run it in the background.
# =============================================================================
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

T="${1:-40}"
shift 2>/dev/null || true
N=("$@")
[ ${#N[@]} -gt 0 ] || N=(6 12 20 40)

TOML="${TOML:-your-system-xtb.toml}"
[ -f "$TOML" ] || { echo "ERROR: $TOML not found (set TOML=...)"; exit 1; }

# the coordinate file is announced INSIDE the toml (root key: input/coord/...)
COORD="${COORD:-$(grep -iE '^[[:space:]]*(input|coord|coords|structure)[[:space:]]*=' "$TOML" \
                  | head -1 | sed -E 's/.*=[[:space:]]*"([^"]+)".*/\1/')}"
[ -n "$COORD" ] && [ -f "$COORD" ] || { echo "ERROR: coordinate file '$COORD' not found"; exit 1; }

command -v crest >/dev/null 2>&1 || { echo "ERROR: crest not in PATH (source ./crest_env.sh $T)"; exit 1; }

echo "threads      : $T"
echo "MTD counts   : ${N[*]}"
echo "toml / coord : $TOML / $COORD"
echo

RESULTS=()
for n in "${N[@]}"; do
  d="bench-nmtd-$n"
  rm -rf "$d" 2>/dev/null
  mkdir -p "$d"
  cp "$TOML" "$d/" && cp "$COORD" "$d/"

  printf '==============================================================\n'
  printf '  -nmtd %s   (cores_per_job = %d)\n' "$n" "$(( T / n > 0 ? T / n : 1 ))"
  printf '==============================================================\n'

  t0=$SECONDS
  ( cd "$d" && crest "$TOML" -T "$T" -nmtd "$n" --quick --noreftopo ) \
      > "$d/crest.log" 2>&1
  rc=$?
  wall=$(( SECONDS - t0 ))

  # per-MTD completion times from the log
  mtd=$(grep -oE '\*MTD +[0-9]+ completed successfully \.+ +[0-9]+ min, *[0-9.]+ sec' \
        "$d/crest.log" | tail -1 || true)
  est=$(grep -oE 'Estimated runtime for a batch of [0-9]+ MTDs on [0-9]+ threads:.*' \
        "$d/crest.log" | tail -1 || true)
  nstruct=$(grep -oE '[0-9]+ of [0-9]+ structures successfully optimized' \
        "$d/crest.log" | tail -1 || true)

  printf '  -> exit=%s  total wall=%s s\n' "$rc" "$wall"
  [ -n "$est"    ] && printf '  -> %s\n' "$est"
  [ -n "$mtd"    ] && printf '  -> slowest %s\n' "$mtd"
  [ -n "$nstruct" ] && printf '  -> %s\n' "$nstruct"
  [ "$rc" -ne 0 ] && printf '  -> FAILED, see %s/crest.log\n' "$d"
  echo

  RESULTS+=("$n|$(( T / n > 0 ? T / n : 1 ))|$wall|$rc|$mtd")
done

printf '==============================================================\n'
printf '  SUMMARY  (-T %s)\n' "$T"
printf '==============================================================\n'
printf '%-8s %-16s %-12s %s\n' "-nmtd" "threads/MTD" "wall / s" "slowest MTD"
for r in "${RESULTS[@]}"; do
  IFS='|' read -r n c w rc m <<<"$r"
  [ "$rc" -eq 0 ] || w="$w(ERR)"
  printf '%-8s %-16s %-12s %s\n' "$n" "$c" "$w" "$m"
done
echo
echo "Pick the smallest wall time; then re-run it WITHOUT --noreftopo / with your"
echo "real sampling budget. Remember the MTD length per job shrinks as -nmtd grows"
echo "(t(MTD) = mdtime_total / nmetadyn), so use -mdtime to keep it long enough."
