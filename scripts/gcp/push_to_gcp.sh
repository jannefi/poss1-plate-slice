#!/usr/bin/env bash
# Run on janne-pc ONLY. Pushes the repo (code only -- work/ is 23 GB of
# unrelated historical run output and must not travel) and the already-built
# Python 3.11 micromamba env to the GCP pilot VM, over plain SSH/rsync. No
# plate FITS here: the VM fetches its own pilot plates directly from IRSA
# (scripts/gcp/fetch_plates.sh) rather than using janne-pc's upload bandwidth.
# No private V-csv either -- this VM never runs vetoes.
#
# Usage: GCP_HOST=<external-ip> ./scripts/gcp/push_to_gcp.sh
set -euo pipefail

: "${GCP_HOST:?Set GCP_HOST to the external IP of the VM}"
GCP_USER="${GCP_USER:-janne}"
KEY_PATH="${KEY_PATH:-$HOME/.ssh/vasco60-gcp-pilot}"
REMOTE_DIR="${REMOTE_DIR:-/home/${GCP_USER}/poss1-plate-slice}"
MICROMAMBA_ENV_LOCAL="${MICROMAMBA_ENV_LOCAL:-/home/janne/.micromamba/envs/vasco-py311}"
# Same reasoning as the EC2 push script: conda-forge's OpenSSL build bakes
# its cert-file path in at build time, so the remote path must match exactly.
REMOTE_ENV_DIR="${REMOTE_ENV_DIR:-/home/janne/.micromamba/envs/vasco-py311}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REMOTE_ENV_PARENT="$(dirname "$REMOTE_ENV_DIR")"
SSH_OPTS=(-i "$KEY_PATH" -o StrictHostKeyChecking=accept-new)

# Lever 4's code (fetch_bright_usnob, VASCO_SPIKE_CATALOG) only exists on this
# branch -- main just prints the env var in a banner without acting on it.
# Steps 1-3 (what the VM runs) are identical between branches, but pinning
# the SAME commit on both machines removes any doubt.
echo "[push_to_gcp] Ensuring local repo is on paper-parity..."
CUR_BRANCH="$(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD)"
if [ "$CUR_BRANCH" != "paper-parity" ]; then
  echo "[push_to_gcp] Switching local checkout from '$CUR_BRANCH' to 'paper-parity'"
  git -C "$REPO_ROOT" checkout paper-parity
fi

echo "[push_to_gcp] Ensuring remote directories exist..."
ssh "${SSH_OPTS[@]}" "${GCP_USER}@${GCP_HOST}" \
  "mkdir -p '$REMOTE_DIR' '$REMOTE_ENV_PARENT'"

echo "[push_to_gcp] Syncing repo (code only -- excludes work/, .git, data symlink target)..."
rsync -avz --progress \
  -e "ssh ${SSH_OPTS[*]}" \
  --exclude '.git/' \
  --exclude 'work/' \
  --exclude 'results/' \
  --exclude '__pycache__/' \
  --exclude '*.pyc' \
  "$REPO_ROOT/" "${GCP_USER}@${GCP_HOST}:${REMOTE_DIR}/"

echo "[push_to_gcp] Syncing working Python 3.11 env (~1.4GB, one-time)..."
rsync -avz --progress \
  -e "ssh ${SSH_OPTS[*]}" \
  "$MICROMAMBA_ENV_LOCAL/" "${GCP_USER}@${GCP_HOST}:${REMOTE_ENV_DIR}/"

# janne-pc's stilts is NOT the Debian apt package (dpkg -S finds nothing for
# it) -- it's the upstream jar from starlink.ac.uk at /opt/stilts/stilts.jar,
# version "3.5-4" (a completely different numbering scheme than Debian's
# packaged "stilts", which apt on a fresh VM resolved to 3.5.2-1: newer,
# pulling in ant/topcat/skyview-java bloat neither machine needs). Shipping
# the identical jar sidesteps any version-matching question for STILTS
# entirely, the same reasoning as shipping the Python env instead of
# reinstalling packages.
STILTS_JAR_LOCAL="${STILTS_JAR_LOCAL:-/opt/stilts/stilts.jar}"
echo "[push_to_gcp] Syncing stilts.jar (~17MB, ships the exact jar, not an apt package)..."
ssh "${SSH_OPTS[@]}" "${GCP_USER}@${GCP_HOST}" "mkdir -p /home/${GCP_USER}/staging"
rsync -avz --progress \
  -e "ssh ${SSH_OPTS[*]}" \
  "$STILTS_JAR_LOCAL" "${GCP_USER}@${GCP_HOST}:/home/${GCP_USER}/staging/stilts.jar"

echo "[push_to_gcp] Done."
echo "  Repo:   ${REMOTE_DIR}  (branch: paper-parity)"
echo "  Python: ${REMOTE_ENV_DIR}/bin/python3.11"
echo "  Next:   ssh -i $KEY_PATH ${GCP_USER}@${GCP_HOST} 'cd ${REMOTE_DIR} && REPO_DIR=${REMOTE_DIR} bash scripts/gcp/bootstrap.sh'"
