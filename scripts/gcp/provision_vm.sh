#!/usr/bin/env bash
# Run on janne-pc. Creates one GCP VM for the full642 run, idempotent
# (check-then-create). Same machine type/zone/image as the validated pilot
# (n2-standard-16, europe-north1-a, Debian 13), disk bumped 150GB -> 250GB:
# fetch_plates.sh downloads a whole arm's plate half up front and never
# deletes the FITS, so ~321 plates x ~380MB =~ 122GB must fit alongside the
# env/OS (~15GB) and transient scratch_tiles -- 250GB leaves a comfortable
# margin without touching the fetch/slice scripts.
#
# Usage: INSTANCE=vasco-full642-gcp-a ./scripts/gcp/provision_vm.sh
set -euo pipefail

: "${INSTANCE:?Set INSTANCE to a unique VM name, e.g. vasco-full642-gcp-a}"
PROJECT="${PROJECT:-project-a54d84d6-a5e7-4acc-b49}"
ZONE="${ZONE:-europe-north1-a}"
MACHINE_TYPE="${MACHINE_TYPE:-n2-standard-16}"
DISK_SIZE_GB="${DISK_SIZE_GB:-250}"
IMAGE_FAMILY="${IMAGE_FAMILY:-debian-13}"
IMAGE_PROJECT="${IMAGE_PROJECT:-debian-cloud}"
SSH_KEY_PATH="${SSH_KEY_PATH:-$HOME/.ssh/vasco60-gcp-pilot.pub}"
SSH_USER="${SSH_USER:-janne}"
SSH_KEYS_FILE="$(mktemp)"
trap 'rm -f "$SSH_KEYS_FILE"' EXIT
printf '%s:%s\n' "$SSH_USER" "$(cat "$SSH_KEY_PATH")" > "$SSH_KEYS_FILE"

EXISTING="$(gcloud compute instances list --filter="name=$INSTANCE" \
  --format="value(name)" --project="$PROJECT" --zones="$ZONE" 2>/dev/null)"
if [ -n "$EXISTING" ]; then
  echo "[provision_vm] $INSTANCE already exists in $ZONE, skipping create."
else
  echo "[provision_vm] Creating $INSTANCE ($MACHINE_TYPE, ${DISK_SIZE_GB}GB, $ZONE)..."
  gcloud compute instances create "$INSTANCE" \
    --project="$PROJECT" \
    --zone="$ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --image-family="$IMAGE_FAMILY" \
    --image-project="$IMAGE_PROJECT" \
    --boot-disk-size="${DISK_SIZE_GB}GB" \
    --boot-disk-type=pd-balanced \
    --metadata-from-file="ssh-keys=${SSH_KEYS_FILE}"
fi

echo "[provision_vm] Waiting for external IP..."
IP=""
for _ in $(seq 1 30); do
  IP="$(gcloud compute instances describe "$INSTANCE" --zone="$ZONE" --project="$PROJECT" \
    --format="value(networkInterfaces[0].accessConfigs[0].natIP)" 2>/dev/null)"
  [ -n "$IP" ] && break
  sleep 2
done
[ -n "$IP" ] || { echo "[provision_vm] [FATAL] no external IP after waiting" >&2; exit 1; }

echo "[provision_vm] $INSTANCE ready. External IP: $IP"
echo "  export GCP_HOST=$IP  # for push_to_gcp.sh / run_and_pull_pilot.sh / etc."
