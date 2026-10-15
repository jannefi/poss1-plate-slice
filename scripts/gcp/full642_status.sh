#!/usr/bin/env bash
# Run on janne-pc, any time, from a fresh shell -- does not depend on any
# particular SSH session or background job being alive. Read-only status
# check across both full642 GCP arms: are the VMs up, how far has each arm's
# fetch/slice/pull/veto gotten, and has the final merged S0 been built yet.
set -uo pipefail

PROJECT="${PROJECT:-project-a54d84d6-a5e7-4acc-b49}"
ZONE="${ZONE:-europe-north1-a}"
GCP_USER="${GCP_USER:-janne}"
KEY_PATH="${KEY_PATH:-$HOME/.ssh/vasco60-gcp-pilot}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

echo "======================================================================"
echo " FULL642 2-ARM STATUS -- $(date -Is)"
echo "======================================================================"

check_arm() {
  local TAG="$1" INSTANCE="$2"
  local LOCAL_ROOT="$REPO_ROOT/work/runs/full642_gcp_${TAG}"
  local PLATES_FILE="$LOCAL_ROOT/plates.txt"
  local PLATE_MAP="$LOCAL_ROOT/plate_map.csv"
  local PLATE_EPOCHS="$LOCAL_ROOT/plate_epochs.csv"
  local N_PLATES_TOTAL=0
  [ -f "$PLATES_FILE" ] && N_PLATES_TOTAL=$(($(tr ',' '\n' < "$PLATES_FILE" | wc -l)))

  echo
  echo "--- ARM $TAG ($INSTANCE) -- $N_PLATES_TOTAL plates assigned ---"

  local VM_LINE IP STATUS
  VM_LINE="$(gcloud compute instances list --filter="name=$INSTANCE" \
    --format="value(name,status,networkInterfaces[0].accessConfigs[0].natIP)" \
    --project="$PROJECT" --zones="$ZONE" 2>/dev/null)"
  if [ -z "$VM_LINE" ]; then
    echo "  VM does not exist (deleted, or not yet created)."
    return
  fi
  read -r _ STATUS IP <<< "$VM_LINE"
  echo "  VM status=$STATUS  ip=${IP:-none}"

  local SSH_OPTS=(-i "$KEY_PATH" -o StrictHostKeyChecking=accept-new -o ConnectTimeout=8 -o BatchMode=yes)
  if [ -n "$IP" ] && ssh "${SSH_OPTS[@]}" "${GCP_USER}@${IP}" true 2>/dev/null; then
    local N_FETCHED
    N_FETCHED="$(ssh "${SSH_OPTS[@]}" "${GCP_USER}@${IP}" \
      "ls /home/${GCP_USER}/data/plates/*.fits 2>/dev/null | wc -l")"
    echo "  plates fetched (VM disk): $N_FETCHED / $N_PLATES_TOTAL"
    local DISK
    DISK="$(ssh "${SSH_OPTS[@]}" "${GCP_USER}@${IP}" "df -h / | tail -1")"
    echo "  VM disk: $(echo "$DISK" | awk '{print $3" used / "$2" total ("$5" full)"}')"
    local RUNNING
    RUNNING="$(ssh "${SSH_OPTS[@]}" "${GCP_USER}@${IP}" \
      "pgrep -af 'fetch_plates.sh|run_fullscale_slice.py|sex |psfex ' 2>/dev/null | grep -v grep")"
    if [ -n "$RUNNING" ]; then
      echo "  currently running on VM:"
      echo "$RUNNING" | sed 's/^/    /' | cut -c1-140
    else
      echo "  nothing running on VM right now"
    fi
  else
    echo "  SSH not reachable (VM may be booting, or fetch/slice already done and this is stale)"
  fi

  local N_PULLED=0
  if [ -f "$PLATES_FILE" ]; then
    while IFS= read -r P; do
      [ -f "$LOCAL_ROOT/.pulled_${P}" ] && N_PULLED=$((N_PULLED+1))
    done < <(tr ',' '\n' < "$PLATES_FILE")
  fi
  echo "  plates pulled to janne-pc: $N_PULLED / $N_PLATES_TOTAL"

  local N_MAP=0 N_EPOCHS=0
  [ -f "$PLATE_MAP" ] && N_MAP=$(($(wc -l < "$PLATE_MAP") - 1))
  [ -f "$PLATE_EPOCHS" ] && N_EPOCHS=$(($(wc -l < "$PLATE_EPOCHS") - 1))
  echo "  tiles recorded in plate_map.csv: $N_MAP   plates with epoch: $N_EPOCHS"

  local N_VETO_LOGS=0 N_VETO_OK=0
  if [ -d "$LOCAL_ROOT/veto_logs" ]; then
    N_VETO_LOGS=$(ls "$LOCAL_ROOT/veto_logs"/*.log 2>/dev/null | wc -l)
    N_VETO_OK=$(grep -l "step4/5 OK" "$LOCAL_ROOT/veto_logs"/*.log 2>/dev/null | wc -l)
  fi
  echo "  step4/5 (veto+spike) logs: $N_VETO_OK OK / $N_VETO_LOGS launched / $N_PLATES_TOTAL total"

  local ORCH_LOG="$LOCAL_ROOT/orchestrator.log"
  if [ -f "$ORCH_LOG" ]; then
    echo "  orchestrator.log tail:"
    tail -3 "$ORCH_LOG" | sed 's/^/    /'
  fi
}

check_arm A vasco-full642-gcp-a
check_arm B vasco-full642-gcp-b

echo
echo "--- FINAL MERGED S0 ---"
LATEST_RUN=$(ls -td "$REPO_ROOT"/work/runs/full642_paper_parity_[0-9]*/ 2>/dev/null | head -1)
if [ -n "$LATEST_RUN" ]; then
  echo "  latest run dir: $LATEST_RUN"
  if [ -f "${LATEST_RUN}stage_S0.csv" ]; then
    N_S0=$(($(wc -l < "${LATEST_RUN}stage_S0.csv") - 1))
    echo "  stage_S0.csv: $N_S0 rows  <-- DONE"
  else
    echo "  stage_S0.csv not built yet"
  fi
else
  echo "  merge_and_build_full642.sh has not been run yet"
fi

echo
echo "======================================================================"
