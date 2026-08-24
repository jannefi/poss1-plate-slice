#!/usr/bin/env bash
# Run ON the GCP VM. Downloads exactly the pilot's plate scans from IRSA's
# bulk, plate-addressed archive -- confirmed live and byte-identical to
# janne-pc's copies (dss1red_XE074.fits: 391985280 bytes both places) -- and
# writes plate_epochs.csv (plate_id,epoch) from each FITS's DATE-OBS, reusing
# vasco.cli_pipeline._plate_epoch_year_from_fits so the epoch math matches
# what the pipeline itself would compute. janne-pc needs this because its
# step4 run never sees the raw plate FITS (only the harvested tile catalogues
# come back), and without an epoch, proper-motion propagation silently skips
# -- a real result change, not a convenience loss.
#
# Usage: PLATES="XE074,XE105,XE002" ./scripts/gcp/fetch_plates.sh
set -euo pipefail

: "${PLATES:?Set PLATES to a comma-separated plate list, e.g. XE074,XE105}"
REPO_DIR="${REPO_DIR:-/home/janne/poss1-plate-slice}"
PLATE_DIR="${PLATE_DIR:-/home/janne/data/plates}"
PYTHON_BIN="${PYTHON_BIN:-/home/janne/.micromamba/envs/vasco-py311/bin/python3.11}"
IRSA_BASE="https://irsa.ipac.caltech.edu/data/DSS/images/dss1red"

mkdir -p "$PLATE_DIR"
IFS=',' read -ra PLATE_ARR <<< "$PLATES"

# Writes plate_epochs.csv INCREMENTALLY, one line per plate right after that
# plate's own download completes (checked via a .done marker, not just the
# .fits file's existence -- curl writes its destination file from the start,
# so mid-download it already "exists" but is truncated). This lets a
# concurrently-running orchestrator (run_and_pull_pilot.sh) start slicing
# plate 1 as soon as it's done, instead of waiting for all N plates to
# download first -- IRSA throttles each connection to ~550-650KB/s, so at
# full 321-plate scale a blocking fetch-everything-then-slice-everything
# design would roughly double total wall clock.
cd "$REPO_DIR"
EPOCHS_CSV="$REPO_DIR/plate_epochs.csv"
[ -f "$EPOCHS_CSV" ] || echo "plate_id,epoch" > "$EPOCHS_CSV"

for P in "${PLATE_ARR[@]}"; do
  DEST="$PLATE_DIR/dss1red_${P}.fits"
  DONE_MARK="${DEST}.done"
  if [ -f "$DONE_MARK" ]; then
    echo "[fetch_plates] $P already present (done marker exists), skipping"
    continue
  fi
  echo "[fetch_plates] fetching $P from IRSA..."
  curl -f -C - -o "$DEST" "$IRSA_BASE/dss1red_${P}.fits"

  EPOCH="$("$PYTHON_BIN" -c "
import sys
sys.path.insert(0, '.')
from vasco.cli_pipeline import _plate_epoch_year_from_fits
e = _plate_epoch_year_from_fits('$DEST')
if e is None:
    raise SystemExit(1)
print(f'{e:.6f}')
")"
  if [ -z "$EPOCH" ]; then
    echo "[fetch_plates] [FATAL] could not read DATE-OBS for $P ($DEST)" >&2
    exit 1
  fi
  echo "${P},${EPOCH}" >> "$EPOCHS_CSV"
  echo "  $P: epoch=$EPOCH"
  touch "$DONE_MARK"
done

echo "[fetch_plates] Done. $(ls "$PLATE_DIR"/*.fits.done 2>/dev/null | wc -l) / ${#PLATE_ARR[@]} plates confirmed complete in $PLATE_DIR"
