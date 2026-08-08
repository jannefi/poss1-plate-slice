#!/usr/bin/env bash
# Full pipeline for a single plate: download -> step2/3 -> step4/5 -> archive.
#
# Usage:
#   tools/process_one_plate.sh <PLATE_ID> [--workers N]
#
# --workers: passed through to step2/3 and step4/5 (default 12).
#
# NOTE (2026-07-27): raised from 6 to 12 (= janne-pc's full vCPU count).
# The original 6 was a network-throttling-motivated cap for step4/5,
# back when it hit live VizieR/MAST -- A/B-tested safe for that reason
# specifically (feedback_parallelize_multitile_ops memory). That
# justification is now gone: step4/5 is fully local-cache-backed
# (Gaia/PS1/USNO-B veto AND the spike-mask bright-star fetch, the
# latter fixed 2026-07-27 after being found to cost ~31s/tile via an
# uncached live MAST call). step2/3 (SExtractor+PSFEx, single-threaded
# per process via NTHREADS 1/OPENBLAS_NUM_THREADS=1 below) is genuinely
# CPU-bound and was confirmed under-subscribed at 6 workers (~70%
# system-wide CPU use, 6 of 12 cores idle, verified live via `ps`/
# `mpstat`). Verified with a properly isolated, sequential A/B
# (identical 49 tiles from an already-archived plate, one worker count
# fully completed -- confirmed via exit code, not inferred from file
# timestamps -- before the other started): 12 workers took 616s vs 759s
# at 6 workers -- a real, clean ~19% improvement. (An earlier attempt at
# this same A/B test produced a bogus ~22% figure from two runs that
# had accidentally overlapped in time and contaminated each other --
# retracted; this is the redone, trustworthy number.) Not the naive-
# expected 2x -- step2/3 has real non-parallelizable overhead (per-tile
# I/O, setup) -- but a genuine improvement, not the regression early
# live-campaign numbers seemed to show (that was plate-to-plate field-
# density variance, not the worker count -- see the campaign memory for
# the full story). janne-pc is a dedicated, single-purpose machine
# (Janne's own words: "I don't use it for anything else, not now or in
# general, ever") so there's no interactive-responsiveness reason to
# hold cores back either.
#
# Uses the naive (--include-corners) tessellation plan (repro/mnras-parity
# deviation #2, see context/REPRO_DEVIATIONS.md) -- every plate is a fixed
# 49-tile grid. Any tile already in data/tiles_archive/ (from an earlier,
# smaller-footprint run) is restored to data/tiles/ before step1; step1
# itself then skips any tile whose tile_status.json already shows step1 ok
# (scripts/run_plan.py's built-in resume logic) and downloads only what's
# genuinely new. Stale per-tile SPREAD_MODEL postscore output (keyed to the
# OLD survivor numbering, invalidated by a step2/3 rerun) is deleted as each
# tile is restored. --force on step2/3 and step4/5 applies uniformly to
# every tile, old or new, so a plate is always left in a fully consistent
# state.
#
# Working mode (per Janne): one plate per invocation. This script does
# not loop over plates itself.
set -euo pipefail


REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# All env vars exported once, up front, before step1 -- previously
# VASCO_REPRO_SINGLE_PASS was exported right before step4/5, which meant
# step2/3 (the actual SExtractor stage) always ran in full two-pass mode
# regardless of intent. Keep this block first so every stage sees a
# consistent environment.
#
# VASCO_REPRO_SINGLE_PASS deliberately NOT set (2026-07-26): real
# campaign-scale data showed the single-pass + post-hoc crop-based
# SPREAD_MODEL restore (tools/spread_model_postscore.py) costs 1.65x MORE
# compute than genuine two-pass SExtractor+PSFEx on every tile, since the
# real average survivor density (20.26/tile) is ~2x the crop approach's
# break-even point (~10.5/tile) -- see memory
# spread_model_crop_approach_cost_reversal. Unset (the default) means
# vasco/cli_pipeline.py runs the already-existing run_pass1/run_psfex/
# run_pass2 path; vasco/mnras/filters_mnras.py's SPREAD_MODEL gate
# activates automatically once that column exists, no other code changes
# needed. This only affects plates processed from here on (XE087+) -- the
# existing single-pass archive (XE002-086) is untouched by this switch.
export VASCO_WCSFIX_DISABLE=1
export VASCO_DISABLE_USNOB=1
# Moved from <cache_root> (spinning HDD, sda3) to here
# (SSD, sdb2) 2026-07-27: step4/5 xmatch/veto against these caches was
# measured at ~416s/plate (51% of total per-plate time) on the HDD vs
# 263s/plate for the same stage on EC2's SSD-backed EBS volume against
# equal-sized caches -- random-access lookups against multi-billion-row
# cache files are exactly the workload a spinning disk struggles with.
# Verified byte-identical (file count + rsync dry-run diff) before the
# HDD copy was deleted.
export VASCO_GAIA_CACHE=<gaia_cache>
export VASCO_USNOB_CACHE=<usnob_cache>
export VASCO_PS1_CACHE=<ps1_cache>

# Two-pass mode invokes PSFEx on every tile (not just survivor-bearing
# ones). Cap BLAS threading per subprocess -- without this, each
# psfex/sex call self-multithreads via OpenBLAS and oversubscribes the
# machine (observed load avg 28-44 on a 12-core box while tuning the S0G
# bulk precompute today, for comparatively little real throughput gain).
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Pre-create data/logs before any parallel step1 downloads start. On a
# fresh machine (no prior runs), N parallel step1-download subprocesses
# all race to create this directory simultaneously via
# vasco/downloader.py's configure_logger() -- observed in practice on a
# brand-new EC2 instance as a FileExistsError('data') crash on every tile,
# something janne-pc never surfaced simply because data/logs already
# existed there from months of prior runs.
mkdir -p data/logs

PLAN_CSV="plans/tiles_mnras_plates_naive.csv"
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

TILE_IDS_FILE="$WORKDIR/tile_ids.txt"      # bare tile ids
TILE_PATHS_FILE="$WORKDIR/tile_paths.txt"  # data/tiles/<id> paths

# --- resolve this plate's tile list from the plan CSV ---
python3 - "$PLAN_CSV" "$PLATE_ID" "$TILE_IDS_FILE" "$TILE_PATHS_FILE" <<'PYEOF'
import csv, sys
plan_csv, plate_id, ids_out, paths_out = sys.argv[1:5]
tile_ids = []
with open(plan_csv) as f:
    for row in csv.DictReader(f):
        if row["plate_id"] == plate_id:
            tile_ids.append(row["tile_id"])
if not tile_ids:
    print(f"Plate {plate_id} not found in {plan_csv}", file=sys.stderr)
    sys.exit(1)
with open(ids_out, "w") as f:
    f.write("\n".join(tile_ids) + "\n")
with open(paths_out, "w") as f:
    f.write("\n".join(f"data/tiles/{t}" for t in tile_ids) + "\n")
print(f"{plate_id}: {len(tile_ids)} tiles")
PYEOF

N_TILES=$(wc -l < "$TILE_IDS_FILE")
echo "=== ${PLATE_ID}: ${N_TILES} tiles ==="

# --- step 1: restore any already-archived tiles, then download what's new ---
T_START=$(date +%s)
echo "--- restoring already-archived tiles from data/tiles_archive/ (if any) ---"
RESTORED=0
while read -r tid; do
  [[ -z "$tid" ]] && continue
  if [[ -d "data/tiles/$tid" ]]; then
    continue
  fi
  if [[ -d "data/tiles_archive/$tid" ]]; then
    mv "data/tiles_archive/$tid" "data/tiles/$tid"
    # Stale-data cleanup: a step2/3 rerun regenerates pass2.ldac and the
    # survivor set, which would silently invalidate any leftover
    # SPREAD_MODEL postscore output (keyed to the OLD survivor numbering).
    rm -rf "data/tiles/$tid/postscore" "data/tiles/$tid/catalogs/spread_model_postscore.csv"
    RESTORED=$((RESTORED + 1))
  fi
done < "$TILE_IDS_FILE"
echo "  restored ${RESTORED} tile(s) from archive"
echo "--- step1 download (skips tiles already marked step1 ok, incl. just-restored ones) ---"
python3 scripts/run_plan.py "$PLAN_CSV" --plate "$PLATE_ID"
T_STEP1=$(date +%s)

# --- step 2/3 ---
echo "--- step2/3 (workers=${WORKERS}) ---"
# --force unconditionally: every plate now mixes already-processed tiles
# (whose tile_status.json shows step2+step3 "ok", which
# run_steps_2_3_parallel.py would otherwise skip via its _steps_done()
# resume logic) with brand-new ones. Redundant reprocessing of
# already-correct tiles is cheap (single-pass, ~1-2s/tile) and keeps every
# tile in a plate on the same footing without tracking partial state.
python3 tools/run_steps_2_3_parallel.py --tiles-file "$TILE_PATHS_FILE" --workers "$WORKERS" --force
T_STEP23=$(date +%s)


# --- step 4/5 ---
echo "--- step4/5 (workers=${WORKERS}) ---"
# --force unconditionally, same reasoning as step2/3 above -- ensures the
# (possibly regenerated) pass2.ldac actually gets re-veto'd and re-filtered
# for every tile rather than reusing stale prior outputs.
# Unlike step1/step2/3, this tool's exit code IS meaningful (1 on partial
# failure) -- but we want our own tile_status.json-based verification
# below to run regardless and report exactly which tile/step failed, so
# don't let `set -e` kill the script here.
python3 tools/run_steps_4_5_parallel.py --tiles-file "$TILE_PATHS_FILE" --workers "$WORKERS" --force || true
T_STEP45=$(date +%s)

# --- verify every tile completed all 5 steps ok (exit codes upstream are
#     not fully reliable, see plan risk notes -- re-check tile_status.json
#     directly rather than trusting subprocess exit codes) ---
# A tile whose step1 status is "skip" (no POSS-I coverage at that
# coordinate, e.g. near the celestial pole -- see cmd_step1_download's
# REJECT_NON_FITS handling) is a valid terminal state, not a failure: it
# never has a FITS to process, so step2-5 are correctly never attempted.
echo "--- verifying tile_status.json for all tiles ---"
FAILED=$(python3 - "$TILE_IDS_FILE" <<'PYEOF'
import json, sys
from pathlib import Path
ids_file = sys.argv[1]
tile_ids = [l.strip() for l in open(ids_file) if l.strip()]
failed = []
skipped = []
for tid in tile_ids:
    p = Path("data/tiles") / tid / "tile_status.json"
    if not p.exists():
        failed.append((tid, "no tile_status.json"))
        continue
    steps = json.loads(p.read_text()).get("steps", {})
    if steps.get("step1", {}).get("status") == "skip":
        skipped.append(tid)
        continue
    for s in ("step1", "step2", "step3", "step4", "step5"):
        if steps.get(s, {}).get("status") != "ok":
            failed.append((tid, f"{s}={steps.get(s, {}).get('status')}"))
            break
if skipped:
    print(f"SKIPPED {len(skipped)} tile(s), no POSS-I coverage: {', '.join(skipped)}")
if failed:
    for tid, reason in failed:
        print(f"FAILED {tid}: {reason}")
    sys.exit(1)
print("ALL_OK")
PYEOF
) || { echo "$FAILED"; echo "One or more tiles failed -- NOT archiving. Fix and re-run before continuing."; exit 1; }
echo "$FAILED"

# --- final full_strict candidate count (for the summary) ---
FINAL_COUNT=$(python3 - "$TILE_IDS_FILE" <<'PYEOF'
import csv, sys
from pathlib import Path
tile_ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]
total = 0
for tid in tile_ids:
    p = Path("data/tiles") / tid / "catalogs" / "sextractor_pass2.filtered.csv"
    if p.exists():
        with open(p) as f:
            total += sum(1 for _ in csv.DictReader(f))
print(total)
PYEOF
)

# --- archive to HDD (only reached if verification passed) ---
mkdir -p data/tiles_archive
echo "--- archiving ${N_TILES} tiles to data/tiles_archive/ (HDD) ---"
SSD_BEFORE=$(df --output=avail -B1 . | tail -1)
while read -r tid; do
  [[ -z "$tid" ]] && continue
  mv "data/tiles/$tid" "data/tiles_archive/$tid"
done < "$TILE_IDS_FILE"
SSD_AFTER=$(df --output=avail -B1 . | tail -1)
FREED_MB=$(( (SSD_AFTER - SSD_BEFORE) / 1024 / 1024 ))
T_END=$(date +%s)

echo ""
echo "===================== SUMMARY: ${PLATE_ID} ====================="
echo "Tiles processed:            ${N_TILES}"
echo "Final full_strict candidates: ${FINAL_COUNT}"
echo "SSD space freed by archive:  ${FREED_MB} MB"
echo "Archived to: data/tiles_archive/ -> <archive_dir>/"
echo "--- timing ---"
echo "step1 download:      $(( T_STEP1 - T_START ))s"
echo "step2/3 SExtractor:   $(( T_STEP23 - T_STEP1 ))s"
echo "step4/5 xmatch+veto:  $(( T_STEP45 - T_STEP23 ))s"
echo "verify+archive:       $(( T_END - T_STEP45 ))s"
echo "TOTAL:                $(( T_END - T_START ))s"
echo "=================================================================="
