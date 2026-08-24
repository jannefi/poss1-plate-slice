#!/usr/bin/env bash
# Run on janne-pc. Deletes the pilot VM so it never sits idle burning the
# 90-day free-trial credit window -- mirrors the always-power-off discipline
# in scripts/ec2_overnight.sh, except here there's no separate persistent
# data volume to preserve: everything worth keeping was already pulled by
# run_and_pull_pilot.sh, so a full delete (not just a stop) is correct.
set -euo pipefail

PROJECT="${PROJECT:-project-a54d84d6-a5e7-4acc-b49}"
ZONE="${ZONE:-europe-north1-a}"
INSTANCE="${INSTANCE:-vasco-paper-parity-pilot}"

echo "[shutdown_vm] Deleting $INSTANCE in $ZONE (project $PROJECT)..."
gcloud compute instances delete "$INSTANCE" --zone="$ZONE" --project="$PROJECT" --quiet
echo "[shutdown_vm] Done. No further GCP compute cost for this pilot."
