#!/usr/bin/env bash
# Reusable "shrinking set" post-process funnel: EDGE -> MORPH -> SHAPE -> MAPS.
# Each stage's kept-survivors CSV feeds the next stage as input, so each
# later (more expensive) stage runs against a much smaller population than
# it would standalone against the full stage_S0.csv.
#
# All four stages are already-existing, documented EXPERIMENTAL scripts,
# run here via their normal CLI -- no stage's internal logic is touched by
# this driver, and it does not wire anything into process_one_plate.sh or
# any per-tile pipeline step. See docs/STAGE_MORPH.md, docs/STAGE_SHAPE.md.
#
# Idempotent at the EDGE step: if <run-dir>/stages/stage_S0_EDGE.csv
# already exists, it is reused as-is rather than recomputed -- lets this
# same script be run end-to-end from a bare stage_S0.csv on a future
# campaign, or resumed here on top of already-computed EDGE output.
#
# Usage:
#   VASCO_MAPS_CACHE=<maps_cache> \
#   scripts/run_post_stage_chain.sh \
#       --run-dir work/runs/post_stage_chain_20260802 \
#       --tiles-root data/tiles_archive \
#       --workers 12
set -euo pipefail

RUN_DIR=""
TILES_ROOT="data/tiles_archive"
WORKERS=12
# Default matches stage_edge_post.py's own default (the naive-tessellation
# archive report). A tiles root that isn't in that report -- e.g. a
# supplement root -- MUST pass its own report here, or every one of its
# rows is "missing" from the report and EDGE degenerates into a silent
# pass-through (see the EDGE coverage guard below).
EDGE_REPORT="data/metadata/tile_plate_edge_report.csv"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="$2"; shift 2 ;;
    --tiles-root) TILES_ROOT="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --edge-report-csv) EDGE_REPORT="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$RUN_DIR" ]]; then
  echo "ERROR: --run-dir is required" >&2
  exit 1
fi
if [[ -z "${VASCO_MAPS_CACHE:-}" ]]; then
  echo "ERROR: VASCO_MAPS_CACHE must be set (MAPS stage requires it) -- refusing to run a partial chain" >&2
  exit 1
fi

STAGES_DIR="$RUN_DIR/stages"
mkdir -p "$STAGES_DIR"

echo "[CHAIN] run-dir=$RUN_DIR tiles-root=$TILES_ROOT workers=$WORKERS"
echo "[CHAIN] edge-report=$EDGE_REPORT"

# --- 1. EDGE (S0) -- idempotent ---
EDGE_OUT="$STAGES_DIR/stage_S0_EDGE.csv"
if [[ -f "$EDGE_OUT" ]]; then
  echo "[CHAIN] EDGE: reusing existing $EDGE_OUT"
else
  echo "[CHAIN] EDGE: computing from stage_S0.csv"
  python3 scripts/stage_edge_post.py \
    --run-dir "$RUN_DIR" \
    --input-glob 'stage_S0.csv' \
    --stage S0 \
    --edge-report-csv "$EDGE_REPORT"
fi
N_EDGE=$(( $(wc -l < "$EDGE_OUT") - 1 ))
echo "[CHAIN] EDGE kept: $N_EDGE"

# Coverage guard: rows whose tile_id has no edge-report entry are KEPT by
# stage_edge_post.py, so a report that doesn't cover this tiles root turns
# EDGE into a no-op that looks like a legitimate 0% cut. That silently
# happened to the first supplement pilot (3,276/3,276 rows missing). Warn
# loudly rather than let the number be quietly meaningless.
N_MISSING=$(python3 -c "
import csv,sys
p='$STAGES_DIR/stage_S0_EDGE_flags.csv'
try:
    rows=list(csv.DictReader(open(p)))
except FileNotFoundError:
    sys.exit(0)
print(sum(1 for r in rows if r.get('edge_report_missing')=='1'))
" 2>/dev/null || echo 0)
if [[ "${N_MISSING:-0}" -gt 0 ]]; then
  echo "[CHAIN] WARNING: $N_MISSING EDGE rows had no edge-report entry and were kept unfiltered."
  echo "[CHAIN] WARNING: check --edge-report-csv actually covers tiles-root=$TILES_ROOT"
fi

# --- 2. MORPH (S0M) ---
echo "[CHAIN] MORPH: running on EDGE output"
python3 scripts/stage_morph_post.py \
  --run-dir "$RUN_DIR" \
  --input-glob 'stages/stage_S0_EDGE.csv' \
  --stage S0M \
  --tiles-root "$TILES_ROOT"
MORPH_OUT="$STAGES_DIR/stage_S0M_MORPH.csv"
N_MORPH=$(( $(wc -l < "$MORPH_OUT") - 1 ))
echo "[CHAIN] MORPH kept: $N_MORPH"

# --- 3. SHAPE (S0S) ---
echo "[CHAIN] SHAPE: running on MORPH output"
python3 scripts/stage_shape_post.py \
  --run-dir "$RUN_DIR" \
  --input-glob 'stages/stage_S0M_MORPH.csv' \
  --stage S0S \
  --tiles-root "$TILES_ROOT" \
  --workers "$WORKERS"
SHAPE_OUT="$STAGES_DIR/stage_S0S_SHAPE.csv"
N_SHAPE=$(( $(wc -l < "$SHAPE_OUT") - 1 ))
echo "[CHAIN] SHAPE kept: $N_SHAPE"

# --- 4. MAPS (S1) ---
# PYTHONPATH=. required: stage_maps_post.py imports vasco.maps_cache_query
# at module scope, and Python sets sys.path[0] to scripts/, not the cwd,
# when invoked as `python3 scripts/stage_maps_post.py` (see memory
# maps_gsc_veto_trial_2026_07_27).
echo "[CHAIN] MAPS: running on SHAPE output"
PYTHONPATH=. python3 scripts/stage_maps_post.py \
  --run-dir "$RUN_DIR" \
  --input-glob 'stages/stage_S0S_SHAPE.csv' \
  --stage S1
MAPS_OUT="$STAGES_DIR/stage_S1_MAPS.csv"
N_MAPS=$(( $(wc -l < "$MAPS_OUT") - 1 ))
echo "[CHAIN] MAPS kept: $N_MAPS"

N_S0=$(( $(wc -l < "$RUN_DIR/stage_S0.csv") - 1 ))
echo
echo "[CHAIN] === Funnel summary ==="
echo "[CHAIN] S0 (input):    $N_S0"
echo "[CHAIN] EDGE kept:     $N_EDGE"
echo "[CHAIN] MORPH kept:    $N_MORPH"
echo "[CHAIN] SHAPE kept:    $N_SHAPE"
echo "[CHAIN] MAPS kept:     $N_MAPS"
python3 -c "print(f'[CHAIN] Overall reduction: {$N_S0} -> {$N_MAPS} ({100*(1-$N_MAPS/$N_S0):.1f}%)')"
