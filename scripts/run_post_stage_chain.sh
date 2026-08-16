#!/usr/bin/env bash
# Reusable "shrinking set" post-process funnel: MORPH -> SHAPE -> MAPS.
# Each stage's kept-survivors CSV feeds the next stage as input, so each
# later (more expensive) stage runs against a much smaller population than
# it would standalone against the full stage_S0.csv.
#
# All four stages are already-existing, documented EXPERIMENTAL scripts,
# run here via their normal CLI -- no stage's internal logic is touched by
# this driver, and it does not wire anything into process_one_plate.sh or
# any per-tile pipeline step. See docs/STAGE_MORPH.md, docs/STAGE_SHAPE.md.
#
# The chain starts at MORPH, reading <run-dir>/stage_S0.csv directly. The
# EDGE stage that once ran first was retired on 2026-08-16 (see the note at
# the MORPH step); scripts/stage_edge_post_v2.py supersedes it and is not
# wired in here.
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
while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="$2"; shift 2 ;;
    --tiles-root) TILES_ROOT="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --edge-report-csv) echo "[CHAIN] --edge-report-csv is obsolete: the EDGE stage was retired 2026-08-16; ignoring" >&2; shift 2 ;;
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

# --- 1. MORPH (S0M) ---
# The EDGE stage that used to run here has been RETIRED. scripts/stage_edge_post.py
# was tile-granular -- it classified a whole ~1 deg tile by its worst boundary
# sample point -- and on any tessellation its precomputed report does not cover,
# every row hit the "missing -> keep" default, so the stage exited ok having done
# nothing. On the released tile set the report intersected in 2 tiles out of
# 25,643, i.e. a 99.99% silent no-op. The warning below it was never enough.
# Replaced by scripts/stage_edge_post_v2.py, which is per-row, needs no
# precomputed report, and flags rather than cuts by default. It is deliberately
# NOT wired in here: see docs/PLATE_EDGE_MASK.md.
echo "[CHAIN] MORPH: running on stage_S0.csv"
python3 scripts/stage_morph_post.py \
  --run-dir "$RUN_DIR" \
  --input-glob 'stage_S0.csv' \
  --stage S0M \
  --tiles-root "$TILES_ROOT"
MORPH_OUT="$STAGES_DIR/stage_S0M_MORPH.csv"
N_MORPH=$(( $(wc -l < "$MORPH_OUT") - 1 ))
echo "[CHAIN] MORPH kept: $N_MORPH"

# --- 2. SHAPE (S0S) ---
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

# --- 3. MAPS (S1) ---
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
echo "[CHAIN] MORPH kept:    $N_MORPH"
echo "[CHAIN] SHAPE kept:    $N_SHAPE"
echo "[CHAIN] MAPS kept:     $N_MAPS"
python3 -c "print(f'[CHAIN] Overall reduction: {$N_S0} -> {$N_MAPS} ({100*(1-$N_MAPS/$N_S0):.1f}%)')"
