#!/usr/bin/env bash
# Run on janne-pc ONLY -- this is the orchestrator for the whole GCP arm of
# the paper-parity pilot. It drives everything over SSH in the janne-pc ->
# VM direction (the VM never needs inbound access to janne-pc's residential
# network, same as the existing EC2 push/pull scripts).
#
# Per plate: slice + 2-pass SExtractor on the VM (steps 1-3 only, no vetoes,
# no catalog mirrors needed there), pull that plate's full tile trees back,
# then kick off janne-pc's own step4/5 (veto + spike mask) for that plate IN
# THE BACKGROUND before moving on to slice the next plate on the VM. This
# overlaps janne-pc's CPU work with GCP's -- previously janne-pc sat idle
# during every slice step and only started its own work after all 9 plates
# were in, which wastes the one machine that's actually free while waiting.
#
# Why one plate at a time on the VM, not a batch: run_fullscale_slice.py's
# scratch_tiles dir is a single fixed path that gets wiped at the START of
# every plate's iteration regardless of --keep-tiles (confirmed by reading
# the loop body) -- so a plate's tiles must be pulled off before the VM is
# told to slice the next one, or they're gone. That constraint is about the
# VM's local disk only; it does not stop janne-pc from working on what it
# already has while the VM moves on.
#
# Levers applied here for step4/5 (see run_veto_and_s0.sh's header for the
# full correction history): veto drops USNO-B (Gaia+PS1 only), spike mask
# switches to USNO-B, WCSFIX stays ON to match the slice step, dedup is NOT
# done here (that needs every plate in, so it stays in run_veto_and_s0.sh).
#
# Usage:
#   GCP_HOST=<external-ip> PLATES="XE074,XE105,XE002,XE366,XE516,XE585" \
#     ./scripts/gcp/run_and_pull_pilot.sh
set -euo pipefail

: "${GCP_HOST:?Set GCP_HOST to the external IP of the VM}"
: "${PLATES:?Set PLATES to a comma-separated plate list}"
GCP_USER="${GCP_USER:-janne}"
KEY_PATH="${KEY_PATH:-$HOME/.ssh/vasco60-gcp-pilot}"
REMOTE_DIR="${REMOTE_DIR:-/home/${GCP_USER}/poss1-plate-slice}"
REMOTE_OUT="${REMOTE_OUT:-/home/${GCP_USER}/work/pilot}"
REMOTE_PLATE_DIR="${REMOTE_PLATE_DIR:-/home/${GCP_USER}/data/plates}"
REMOTE_PYTHON="${REMOTE_PYTHON:-/home/${GCP_USER}/.micromamba/envs/vasco-py311/bin/python3.11}"
SLICE_WORKERS="${SLICE_WORKERS:-14}"
# Lower than a full 12-core budget: two plates' step4/5 can occasionally
# overlap in the background if the VM outpaces janne-pc, and 6+6 still fits
# nproc without starving whichever one is further along.
VETO_WORKERS="${VETO_WORKERS:-6}"
# Minimum free space (GB) on the filesystem backing $LOCAL_TILES before
# starting a plate's pull. One plate's tile tree is ~3-4GB; 50GB is a wide
# margin. Added after the 2026-08-25/26 incident: a full disk crashed both
# arms mid-rsync (silent, 12h unnoticed) instead of stopping cleanly.
MIN_FREE_GB="${MIN_FREE_GB:-50}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
LOCAL_ROOT="${LOCAL_ROOT:-$REPO_ROOT/work/runs/paper_parity_pilot}"
LOCAL_TILES="$LOCAL_ROOT/tiles"
PLATE_MAP="$LOCAL_ROOT/plate_map.csv"
PLATE_EPOCHS="$LOCAL_ROOT/plate_epochs.csv"
LOG_DIR="$LOCAL_ROOT/veto_logs"
mkdir -p "$LOCAL_TILES" "$LOG_DIR"
PYTHON_BIN="${PYTHON_BIN:-/home/janne/.micromamba/envs/vasco-py311/bin/python3.11}"

# Record identity for the watchdog (scripts/gcp/watchdog.sh) so it can find
# and safely restart this process without guessing or pattern-matching.
echo "$$" > "$LOCAL_ROOT/orchestrator.pid"
echo "$GCP_HOST" > "$LOCAL_ROOT/gcp_host.txt"
[ -n "${INSTANCE:-}" ] && echo "$INSTANCE" > "$LOCAL_ROOT/gcp_instance.txt"

check_disk_space() {
  local AVAIL_KB; AVAIL_KB=$(df --output=avail -k "$LOCAL_TILES" | tail -1 | tr -d ' ')
  local AVAIL_GB=$((AVAIL_KB / 1024 / 1024))
  if [ "$AVAIL_GB" -lt "$MIN_FREE_GB" ]; then
    echo "[FATAL] low disk space on $(df --output=target "$LOCAL_TILES" | tail -1 | tr -d ' '): ${AVAIL_GB}GB free, need >= ${MIN_FREE_GB}GB -- stopping before the next pull rather than crashing mid-rsync" >&2
    exit 1
  fi
}

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
[ "$BRANCH" = "paper-parity" ] || { echo "[FATAL] repo is on '$BRANCH', not paper-parity -- VASCO_SPIKE_CATALOG is dead on main" >&2; exit 1; }

# Step4/5 env, set once for every background job this script launches.
export VASCO_GAIA_CACHE="${VASCO_GAIA_CACHE:-/home/janne/local_cache/gaia}"
export VASCO_PS1_CACHE="${VASCO_PS1_CACHE:-/home/janne/local_cache/ps1}"
export VASCO_USNOB_CACHE="${VASCO_USNOB_CACHE:-/home/janne/local_cache/usnob}"
export VASCO_LDAC_DROP_VIGNET=1
export VASCO_DISABLE_USNOB=1        # lever 3: drop USNO-B from the veto
export VASCO_SPIKE_CATALOG=usnob    # lever 4: spike mask -> USNO-B
unset VASCO_WCSFIX_DISABLE          # lever 1: WCSFIX stays ON, matches the slice step
unset VASCO_CIRCLE_ARCMIN
echo "[CONFIG] step4/5 will run with veto=gaia+ps1 (usnob dropped)  spike_cat=usnob  wcsfix=ON"

SSH_OPTS=(-i "$KEY_PATH" -o StrictHostKeyChecking=accept-new
          -o ControlMaster=auto -o "ControlPath=/tmp/vasco-gcp-ssh-%r@%h:%p" -o ControlPersist=10m)
ssh_run() { ssh "${SSH_OPTS[@]}" "${GCP_USER}@${GCP_HOST}" "$@"; }

if [ ! -f "$PLATE_MAP" ]; then
  echo "plate_id,tile_id" > "$PLATE_MAP"
fi

# fetch_plates.sh (run separately, in the background, overlapping with this
# script) downloads plates and writes plate_epochs.csv INCREMENTALLY, one
# line per plate, marking each done with a dss1red_<PLATE>.fits.done file
# only once that plate's curl has actually succeeded (not just started --
# curl writes its destination from byte 0, so a bare .fits file's existence
# does not mean the download is complete). wait_for_plate_fetched polls for
# that marker before this loop tells the VM to slice a plate, so slicing
# stays plate-synchronized with fetching without ever blocking on the WHOLE
# plate list being fetched first -- IRSA's ~550-650KB/s per-connection cap
# means waiting for all 321 plates up front would roughly double wall clock.
wait_for_plate_fetched() {
  local PLATE="$1"
  local MARK="${REMOTE_PLATE_DIR}/dss1red_${PLATE}.fits.done"
  local WAITED=0
  local TIMEOUT=2400   # 40 min -- generous vs. the observed ~11 min/plate fetch rate
  while ! ssh_run "test -f '$MARK'" 2>/dev/null; do
    if [ "$WAITED" -ge "$TIMEOUT" ]; then
      echo "[$PLATE] [FATAL] fetch not done after ${TIMEOUT}s -- check fetch_plates.log on the VM" >&2
      return 1
    fi
    sleep 15
    WAITED=$((WAITED+15))
  done
}

run_step45_for_plate() {
  local PLATE="$1"
  local TILES_FILE="$LOCAL_ROOT/tiles_${PLATE}.txt"
  awk -F, -v p="$PLATE" -v root="$LOCAL_TILES" \
    'NR>1 && $1==p {print root "/" $2}' "$PLATE_MAP" > "$TILES_FILE"
  local N_TILES; N_TILES=$(wc -l < "$TILES_FILE")
  if [ "$N_TILES" -eq 0 ]; then
    echo "[$PLATE] [WARN] no tiles found in plate_map for step4/5" >&2
    return 1
  fi
  local EPOCH; EPOCH=$(awk -F, -v p="$PLATE" 'NR>1 && $1==p {print $2}' "$PLATE_EPOCHS")
  if [ -z "$EPOCH" ]; then
    echo "[$PLATE] [FATAL] no epoch found in $PLATE_EPOCHS" >&2
    return 1
  fi
  VASCO_PLATE_EPOCH_YEAR="$EPOCH" "$PYTHON_BIN" tools/run_steps_4_5_parallel.py \
    --tiles-file "$TILES_FILE" --workers "$VETO_WORKERS"
}

BG_PIDS=()
IFS=',' read -ra PLATE_ARR <<< "$PLATES"
for PLATE in "${PLATE_ARR[@]}"; do
  DONE_MARK="$LOCAL_ROOT/.pulled_${PLATE}"
  if [ -f "$DONE_MARK" ]; then
    echo "[$PLATE] already pulled, skipping slice+pull (still queuing step4/5 if not done)"
  else
    check_disk_space
    echo "[$PLATE] waiting for fetch_plates.sh to finish downloading this plate on the VM..."
    wait_for_plate_fetched "$PLATE" || exit 1
    rsync -az -e "ssh ${SSH_OPTS[*]}" \
      "${GCP_USER}@${GCP_HOST}:${REMOTE_DIR}/plate_epochs.csv" "$PLATE_EPOCHS"

    echo "[$PLATE] slicing + 2-pass SExtractor on GCP (no vetoes, --keep-tiles)..."
    # WCSFIX left at pipeline default (ON, lever 1 = release). Two-pass
    # required for SPREAD_MODEL -- do not set VASCO_REPRO_SINGLE_PASS.
    ssh_run "cd '$REMOTE_DIR' && \
      '$REMOTE_PYTHON' tools/run_fullscale_slice.py \
        --out-dir '$REMOTE_OUT' --plate-dir '$REMOTE_PLATE_DIR' \
        --plates '$PLATE' --crpix-table data/plate_crpix_table.csv \
        --keep-tiles --workers $SLICE_WORKERS"

    echo "[$PLATE] pulling tile trees to $LOCAL_TILES ..."
    rsync -az --progress \
      -e "ssh ${SSH_OPTS[*]}" \
      "${GCP_USER}@${GCP_HOST}:${REMOTE_OUT}/scratch_tiles/" "$LOCAL_TILES/"

    echo "[$PLATE] recording plate_id map for this plate's tiles..."
    # scratch_tiles held exactly this plate's tiles when we rsynced it (it was
    # wiped and re-sliced fresh at the start of this plate's run on the VM),
    # so every tile_RA*_DEC* dir just pulled belongs to $PLATE.
    ssh_run "ls '$REMOTE_OUT/scratch_tiles/'" | while read -r TID; do
      [ -n "$TID" ] && echo "${PLATE},${TID}" >> "$PLATE_MAP"
    done

    echo "[$PLATE] freeing scratch_tiles on the VM before the next plate..."
    ssh_run "rm -rf '$REMOTE_OUT/scratch_tiles'"

    touch "$DONE_MARK"
  fi

  echo "[$PLATE] launching step4/5 (veto+spike) in the background -- GCP moves on to the next plate now"
  ( run_step45_for_plate "$PLATE" > "$LOG_DIR/${PLATE}.log" 2>&1 && \
    echo "[$PLATE] step4/5 OK" >> "$LOG_DIR/${PLATE}.log" || \
    echo "[$PLATE] step4/5 FAILED, see $LOG_DIR/${PLATE}.log" ) &
  BG_PIDS+=($!)
done

echo "[run_and_pull_pilot] all plates sliced+pulled. Waiting for background step4/5 jobs to finish..."
FAILED=0
for PID in "${BG_PIDS[@]}"; do
  wait "$PID" || FAILED=$((FAILED+1))
done
if [ "$FAILED" -gt 0 ]; then
  echo "[run_and_pull_pilot] [WARN] $FAILED step4/5 job(s) failed -- check $LOG_DIR/*.log before running build_run_stage_csvs.py" >&2
else
  echo "[run_and_pull_pilot] all step4/5 jobs completed OK"
fi

echo "[run_and_pull_pilot] Tiles: $LOCAL_TILES  Plate map: $PLATE_MAP"
echo "Next: scripts/gcp/shutdown_vm.sh, then scripts/gcp/run_veto_and_s0.sh (fast -- step4/5 is already done, it just builds S0)"
