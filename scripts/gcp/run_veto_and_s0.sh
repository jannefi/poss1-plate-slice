#!/usr/bin/env bash
# Run on janne-pc ONLY, after run_and_pull_pilot.sh has pulled every pilot
# plate's tile trees. Runs steps 4-5 (veto + spike mask) locally, using
# janne-pc's existing Gaia/PS1/USNO-B mirrors -- the only reason this split
# is possible at all, since those mirrors exist nowhere else -- then builds
# the pilot's S0.
#
# run_and_pull_pilot.sh already launches step4/5 per plate in the background
# as soon as that plate is pulled (overlapping janne-pc's CPU work with GCP
# slicing the next plate), so by the time this script runs, most or all of
# that work is already done. The per-plate loop below is now a fast,
# idempotent safety-net pass -- run_steps_4_5_parallel.py skips any tile
# whose tile_status.json already shows step4/step5 == ok, so re-running it
# here just mops up anything that did not finish in the background (or lets
# this script be run standalone if the background approach was skipped).
#
# Levers, corrected after this plan was first drafted wrong (see
# BRANCH_PAPER_PARITY.md and 03_NEXT_ACTIONS.md's lever table -- veto and
# spike-mask directions were backwards in the first pass):
#   1 WCSFIX ON       -- pipeline default; must match what the VM used while
#                        slicing, so left alone here too (do not set
#                        VASCO_WCSFIX_DISABLE).
#   2 dedup 0.25"     -- passed explicitly to build_run_stage_csvs.py below;
#                        its own default is 3.0", correct only for raw WCS.
#   3 veto: drop USNO-B (Gaia+PS1 only) -- VASCO_DISABLE_USNOB=1. Pipeline
#                        default is all three ON.
#   4 spike mask: USNO-B -- VASCO_SPIKE_CATALOG=usnob. Only real on the
#                        paper-parity branch; needs VASCO_USNOB_CACHE.
#
# Per-plate epoch: this branch's step4 can resume from a harvested catalogue
# with no tile FITS present, but proper-motion propagation then needs an
# explicit VASCO_PLATE_EPOCH_YEAR (see fetch_plates.sh) or it silently skips
# PM propagation -- a real result change, not a convenience loss.
#
# Usage: PLATES="XE074,XE105,XE002,XE366,XE516,XE585" ./scripts/gcp/run_veto_and_s0.sh
set -euo pipefail

: "${PLATES:?Set PLATES to the same comma-separated plate list used for the GCP arm}"
WORKERS="${WORKERS:-8}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
LOCAL_ROOT="${LOCAL_ROOT:-$REPO_ROOT/work/runs/paper_parity_pilot}"
LOCAL_TILES="$LOCAL_ROOT/tiles"
PLATE_MAP="$LOCAL_ROOT/plate_map.csv"
PLATE_EPOCHS="$LOCAL_ROOT/plate_epochs.csv"
PYTHON_BIN="${PYTHON_BIN:-/home/janne/.micromamba/envs/vasco-py311/bin/python3.11}"

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
[ "$BRANCH" = "paper-parity" ] || { echo "[FATAL] repo is on '$BRANCH', not paper-parity -- VASCO_SPIKE_CATALOG is dead on main" >&2; exit 1; }

test -f "$PLATE_MAP" || { echo "[FATAL] $PLATE_MAP missing -- run run_and_pull_pilot.sh first" >&2; exit 1; }
test -f "$PLATE_EPOCHS" || { echo "[FATAL] $PLATE_EPOCHS missing -- run run_and_pull_pilot.sh first" >&2; exit 1; }

export VASCO_GAIA_CACHE="${VASCO_GAIA_CACHE:-/home/janne/local_cache/gaia}"
export VASCO_PS1_CACHE="${VASCO_PS1_CACHE:-/home/janne/local_cache/ps1}"
export VASCO_USNOB_CACHE="${VASCO_USNOB_CACHE:-/home/janne/local_cache/usnob}"
export VASCO_LDAC_DROP_VIGNET=1
export VASCO_DISABLE_USNOB=1        # lever 3: drop USNO-B from the veto
export VASCO_SPIKE_CATALOG=usnob    # lever 4: spike mask -> USNO-B
unset VASCO_WCSFIX_DISABLE          # lever 1: WCSFIX stays ON, matches the slice step
unset VASCO_CIRCLE_ARCMIN           # square tiles, no circular cut

echo "[CONFIG] veto=gaia+ps1 (usnob dropped)  spike_cat=usnob  wcsfix=ON  dedup=0.25\""

IFS=',' read -ra PLATE_ARR <<< "$PLATES"
for PLATE in "${PLATE_ARR[@]}"; do
  TILES_FILE="$LOCAL_ROOT/tiles_${PLATE}.txt"
  awk -F, -v p="$PLATE" -v root="$LOCAL_TILES" \
    'NR>1 && $1==p {print root "/" $2}' "$PLATE_MAP" > "$TILES_FILE"
  N_TILES=$(wc -l < "$TILES_FILE")
  if [ "$N_TILES" -eq 0 ]; then
    echo "[$PLATE] [WARN] no tiles found in $PLATE_MAP -- skipping"
    continue
  fi

  EPOCH=$(awk -F, -v p="$PLATE" 'NR>1 && $1==p {print $2}' "$PLATE_EPOCHS")
  if [ -z "$EPOCH" ]; then
    echo "[$PLATE] [FATAL] no epoch found in $PLATE_EPOCHS" >&2
    exit 1
  fi
  echo "[$PLATE] $N_TILES tiles, epoch=$EPOCH -- running step4+5..."
  VASCO_PLATE_EPOCH_YEAR="$EPOCH" "$PYTHON_BIN" tools/run_steps_4_5_parallel.py \
    --tiles-file "$TILES_FILE" --workers "$WORKERS"
done

RUN_TAG="paper_parity_pilot_$(date +%Y%m%d)"
echo "[run_veto_and_s0] building S0 (dedup 0.25\", tag $RUN_TAG)..."
"$PYTHON_BIN" scripts/build_run_stage_csvs.py \
  --tiles-root "$LOCAL_TILES" \
  --plate-map-csv "$PLATE_MAP" \
  --dedup-tol-arcsec 0.25 \
  --run-root work/runs \
  --run-tag "$RUN_TAG"

echo "[run_veto_and_s0] Done. See work/runs/$RUN_TAG/stage_S0.csv"
