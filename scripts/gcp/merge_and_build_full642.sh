#!/usr/bin/env bash
# Run on janne-pc, after both full642 GCP arms (A, B) have completed
# run_and_pull_pilot.sh and run_veto_and_s0.sh's safety-net pass. Merges both
# arms' tile trees into one combined root and runs a SINGLE dedup/S0 build
# over the full 642-plate set -- dedup must see the whole set at once, same
# as the 9-plate pilot's single pass over all 9 plates, not per-arm.
#
# Usage: ./scripts/gcp/merge_and_build_full642.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
PYTHON_BIN="${PYTHON_BIN:-/home/janne/.micromamba/envs/vasco-py311/bin/python3.11}"

ROOT_A="work/runs/full642_gcp_A"
ROOT_B="work/runs/full642_gcp_B"
COMBINED="work/runs/full642_gcp"
COMBINED_TILES="$COMBINED/tiles"
COMBINED_PLATE_MAP="$COMBINED/plate_map.csv"
COMBINED_EPOCHS="$COMBINED/plate_epochs.csv"

for r in "$ROOT_A" "$ROOT_B"; do
  test -f "$r/plate_map.csv" || { echo "[FATAL] $r/plate_map.csv missing -- arm not finished" >&2; exit 1; }
  test -d "$r/tiles" || { echo "[FATAL] $r/tiles missing -- arm not finished" >&2; exit 1; }
done

mkdir -p "$COMBINED_TILES"

echo "[merge] checking arm A / arm B tile-ID sets are disjoint..."
TILES_A="$(mktemp)"; TILES_B="$(mktemp)"
awk -F, 'NR>1 {print $2}' "$ROOT_A/plate_map.csv" | sort -u > "$TILES_A"
awk -F, 'NR>1 {print $2}' "$ROOT_B/plate_map.csv" | sort -u > "$TILES_B"
OVERLAP="$(comm -12 "$TILES_A" "$TILES_B")"
if [ -n "$OVERLAP" ]; then
  echo "[FATAL] arm A and arm B pulled overlapping tile IDs -- refusing to merge:" >&2
  echo "$OVERLAP" | head -20 >&2
  rm -f "$TILES_A" "$TILES_B"
  exit 1
fi
rm -f "$TILES_A" "$TILES_B"
echo "[merge] OK: tile-ID sets disjoint."

echo "[merge] rsyncing arm A tiles into $COMBINED_TILES ..."
rsync -a "$ROOT_A/tiles/" "$COMBINED_TILES/"
echo "[merge] rsyncing arm B tiles into $COMBINED_TILES ..."
rsync -a "$ROOT_B/tiles/" "$COMBINED_TILES/"

echo "[merge] combining plate_map.csv ..."
{ head -1 "$ROOT_A/plate_map.csv"; tail -n +2 "$ROOT_A/plate_map.csv"; tail -n +2 "$ROOT_B/plate_map.csv"; } > "$COMBINED_PLATE_MAP"
N_A=$(( $(wc -l < "$ROOT_A/plate_map.csv") - 1 ))
N_B=$(( $(wc -l < "$ROOT_B/plate_map.csv") - 1 ))
N_COMBINED=$(( $(wc -l < "$COMBINED_PLATE_MAP") - 1 ))
[ "$N_COMBINED" -eq $((N_A + N_B)) ] || { echo "[FATAL] combined plate_map row count mismatch: $N_COMBINED != $N_A + $N_B" >&2; exit 1; }
echo "[merge] plate_map.csv: $N_A (A) + $N_B (B) = $N_COMBINED rows, confirmed."

echo "[merge] combining plate_epochs.csv ..."
{ head -1 "$ROOT_A/plate_epochs.csv"; tail -n +2 "$ROOT_A/plate_epochs.csv"; tail -n +2 "$ROOT_B/plate_epochs.csv"; } > "$COMBINED_EPOCHS"

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
[ "$BRANCH" = "paper-parity" ] || { echo "[FATAL] repo is on '$BRANCH', not paper-parity" >&2; exit 1; }

RUN_TAG="full642_paper_parity_$(date +%Y%m%d)"
echo "[merge] building S0 (dedup 0.25\", tag $RUN_TAG)..."
"$PYTHON_BIN" scripts/build_run_stage_csvs.py \
  --tiles-root "$COMBINED_TILES" \
  --plate-map-csv "$COMBINED_PLATE_MAP" \
  --dedup-tol-arcsec 0.25 \
  --run-root work/runs \
  --run-tag "$RUN_TAG"

echo "[merge] Done. See work/runs/$RUN_TAG/stage_S0.csv"
