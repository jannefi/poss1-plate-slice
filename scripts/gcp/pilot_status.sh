#!/usr/bin/env bash
# Run on janne-pc. Read-only status check across the whole paper-parity pilot
# split: is the VM up, what has it downloaded/sliced, what has janne-pc
# actually pulled and verified, and has veto+S0 run yet. Meant to be safe to
# run at any time, including while run_and_pull_pilot.sh is mid-loop.
set -uo pipefail

PROJECT="${PROJECT:-project-a54d84d6-a5e7-4acc-b49}"
ZONE="${ZONE:-europe-north1-a}"
INSTANCE="${INSTANCE:-vasco-paper-parity-pilot}"
GCP_HOST="${GCP_HOST:-}"
GCP_USER="${GCP_USER:-janne}"
KEY_PATH="${KEY_PATH:-$HOME/.ssh/vasco60-gcp-pilot}"
REMOTE_DIR="${REMOTE_DIR:-/home/${GCP_USER}/poss1-plate-slice}"
REMOTE_OUT="${REMOTE_OUT:-/home/${GCP_USER}/work/pilot}"
REMOTE_PLATE_DIR="${REMOTE_PLATE_DIR:-/home/${GCP_USER}/data/plates}"
PLATES="${PLATES:-XE074,XE105,XE002,XE366,XE516,XE585,XE484,XE540,XE352}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOCAL_ROOT="${LOCAL_ROOT:-$REPO_ROOT/work/runs/paper_parity_pilot}"
LOCAL_TILES="$LOCAL_ROOT/tiles"
PLATE_MAP="$LOCAL_ROOT/plate_map.csv"
PLATE_EPOCHS="$LOCAL_ROOT/plate_epochs.csv"

IFS=',' read -ra PLATE_ARR <<< "$PLATES"
PROBLEMS=0
note_problem() { PROBLEMS=$((PROBLEMS+1)); echo "  [PROBLEM] $*"; }

echo "======================================================================"
echo " PAPER-PARITY PILOT STATUS -- $(date -Is)"
echo "======================================================================"

echo
echo "--- 1. GCP VM ---"
VM_LINE="$(gcloud compute instances list --filter="name=$INSTANCE" \
  --format="value(name,zone,status,networkInterfaces[0].accessConfigs[0].natIP)" \
  --project="$PROJECT" 2>/dev/null)"
if [ -z "$VM_LINE" ]; then
  echo "  VM '$INSTANCE' does not exist (deleted, or never created, or wrong project/zone)."
else
  read -r NAME ZONE_ACTUAL STATUS IP <<< "$VM_LINE"
  echo "  $NAME  zone=$ZONE_ACTUAL  status=$STATUS  ip=${IP:-none}"
  [ "$STATUS" = "RUNNING" ] || note_problem "VM status is '$STATUS', not RUNNING"
  [ -n "$GCP_HOST" ] || GCP_HOST="$IP"
fi

if [ -n "$GCP_HOST" ]; then
  SSH_OPTS=(-i "$KEY_PATH" -o StrictHostKeyChecking=accept-new -o ConnectTimeout=8 -o BatchMode=yes)
  if ssh "${SSH_OPTS[@]}" "${GCP_USER}@${GCP_HOST}" true 2>/dev/null; then
    echo "  SSH reachable at $GCP_HOST"
    echo
    echo "--- 2. GCP-side plate downloads & disk ---"
    ssh "${SSH_OPTS[@]}" "${GCP_USER}@${GCP_HOST}" \
      "ls -la '$REMOTE_PLATE_DIR' 2>/dev/null | grep -c '\.fits\$' || echo 0" \
      | xargs -I{} echo "  plates downloaded: {} / ${#PLATE_ARR[@]}"
    ssh "${SSH_OPTS[@]}" "${GCP_USER}@${GCP_HOST}" "df -h / | tail -1" | \
      awk '{print "  disk: " $3 " used / " $2 " total (" $5 " full)"}'
    ssh "${SSH_OPTS[@]}" "${GCP_USER}@${GCP_HOST}" \
      "ls '$REMOTE_OUT/scratch_tiles' 2>/dev/null | wc -l" | \
      xargs -I{} echo "  scratch_tiles currently resident on VM: {} tile dirs (should be 0 between plates, 49 mid-plate)"
    RUNNING="$(ssh "${SSH_OPTS[@]}" "${GCP_USER}@${GCP_HOST}" \
      "pgrep -af 'run_fullscale_slice.py|sex |psfex ' 2>/dev/null | grep -v grep" 2>/dev/null)"
    if [ -n "$RUNNING" ]; then
      echo "  currently running on VM:"
      echo "$RUNNING" | sed 's/^/    /'
    else
      echo "  nothing running on VM right now"
    fi
  else
    echo "  SSH NOT reachable at $GCP_HOST (VM may still be booting, or already deleted)"
  fi
else
  echo "  (no IP known -- VM absent, skipping remote checks)"
fi

echo
echo "--- 3. janne-pc: per-plate pull status ---"
printf "  %-8s %-8s %-10s %-10s\n" "PLATE" "PULLED" "TILES" "DETECTIONS"
for PLATE in "${PLATE_ARR[@]}"; do
  MARK="$LOCAL_ROOT/.pulled_${PLATE}"
  if [ -f "$MARK" ]; then
    N_TILES=0; N_DET=0
    if [ -f "$PLATE_MAP" ]; then
      N_TILES=$(awk -F, -v p="$PLATE" 'NR>1 && $1==p' "$PLATE_MAP" | wc -l)
    fi
    if [ -d "$LOCAL_TILES" ]; then
      N_DET=$(awk -F, -v p="$PLATE" 'NR>1 && $1==p {print $2}' "$PLATE_MAP" 2>/dev/null | \
        while read -r TID; do
          f="$LOCAL_TILES/$TID/catalogs/sextractor_pass2.csv"
          [ -f "$f" ] && wc -l < "$f"
        done | awk '{s+=$1-1} END {print s+0}')
    fi
    printf "  %-8s %-8s %-10s %-10s\n" "$PLATE" "yes" "$N_TILES" "$N_DET"
    [ "$N_TILES" -eq 49 ] || note_problem "$PLATE: $N_TILES tiles pulled, expected 49"
  else
    printf "  %-8s %-8s %-10s %-10s\n" "$PLATE" "no" "-" "-"
  fi
done

echo
echo "--- 4. plate_epochs.csv ---"
if [ -f "$PLATE_EPOCHS" ]; then
  N_EPOCHS=$(($(wc -l < "$PLATE_EPOCHS") - 1))
  echo "  $N_EPOCHS plates have an epoch recorded"
  # Only a problem if the plate was already pulled but has no epoch --
  # missing epochs for plates not processed yet are just pending, not broken.
  for PLATE in "${PLATE_ARR[@]}"; do
    if [ -f "$LOCAL_ROOT/.pulled_${PLATE}" ] && ! grep -q "^${PLATE}," "$PLATE_EPOCHS"; then
      note_problem "$PLATE was pulled but has no epoch recorded"
    fi
  done
else
  echo "  not pulled yet"
fi

echo
echo "--- 5. veto + S0 build ---"
# Only real dated run-tag dirs (paper_parity_pilot_YYYYMMDD/), not the
# one-off paper_parity_pilot_native_check/ used for the plate-1 parity gate.
LATEST_RUN=$(ls -td "$REPO_ROOT"/work/runs/paper_parity_pilot_[0-9]*/ 2>/dev/null | head -1)
if [ -n "$LATEST_RUN" ]; then
  echo "  latest run dir: $LATEST_RUN"
  if [ -f "${LATEST_RUN}stage_S0.csv" ]; then
    N_S0=$(($(wc -l < "${LATEST_RUN}stage_S0.csv") - 1))
    echo "  stage_S0.csv: $N_S0 rows"
  else
    echo "  stage_S0.csv not built yet"
  fi
else
  echo "  no run_veto_and_s0.sh output yet"
fi

echo
echo "======================================================================"
if [ "$PROBLEMS" -eq 0 ]; then
  echo " OK -- no problems detected."
else
  echo " $PROBLEMS problem(s) flagged above."
fi
echo "======================================================================"
