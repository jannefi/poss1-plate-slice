#!/usr/bin/env bash
# Run ON the GCP VM itself, after push_to_gcp.sh has copied the repo and the
# working Python 3.11 env over. Installs SExtractor/PSFEx/STILTS pinned to
# the EXACT versions running on janne-pc, and refuses to continue on a
# mismatch rather than silently drifting -- last time this pipeline was split
# across machines (the EC2 arm), PSFEx came back newer than janne-pc's and
# had to be validated for output parity after the fact. Pinning turns that
# into a loud failure at setup time instead of a surprise discovered later.
set -euo pipefail

REPO_DIR="${REPO_DIR:-/home/janne/poss1-plate-slice}"
PYTHON_BIN="${PYTHON_BIN:-/home/janne/.micromamba/envs/vasco-py311/bin/python3.11}"

# Known-good, read off janne-pc on 2026-08-24 (both machines run Debian 13
# trixie, so these should be exact apt matches from the same archive).
WANT_SEXTRACTOR="2.28.2+ds-1"
WANT_PSFEX="3.24.2-2"
WANT_STILTS_VERSION="3.5-4"   # `stilts -version` output, checked separately below

echo "[bootstrap] apt install pinned source-extractor/psfex + a bare JRE + rsync ..."
sudo apt-get update
sudo apt-get install -y \
    "source-extractor=${WANT_SEXTRACTOR}" \
    "psfex=${WANT_PSFEX}" \
    openjdk-21-jre-headless rsync
# NOT apt-get install stilts: janne-pc's stilts (`dpkg -S` finds nothing for
# it) is not the Debian package at all -- it's the upstream jar from
# starlink.ac.uk at /opt/stilts/stilts.jar. The Debian package resolves to a
# newer, differently-numbered build (3.5.2-1 vs janne-pc's 3.5-4) and drags
# in ant/topcat/skyview-java for no reason. push_to_gcp.sh already shipped
# the identical jar to ~/staging/stilts.jar; install it exactly where
# janne-pc keeps it so the wrapper script matches too.
echo "[bootstrap] Installing the shipped stilts.jar (not the apt package)..."
sudo mkdir -p /opt/stilts
sudo cp "/home/$(whoami)/staging/stilts.jar" /opt/stilts/stilts.jar
printf '#!/usr/bin/env sh\nexec java -jar /opt/stilts/stilts.jar "$@"\n' | sudo tee /usr/local/bin/stilts >/dev/null
sudo chmod +x /usr/local/bin/stilts

echo "[bootstrap] Linking sex -> source-extractor (the apt package does not"
echo "provide this alternative itself; vasco/pipeline_split.py's"
echo "_find_binary(['sex', 'sextractor']) looks for 'sex' specifically)..."
sudo ln -sf /usr/bin/source-extractor /usr/bin/sex

echo "[bootstrap] Verifying installed versions against janne-pc..."
GOT_SEXTRACTOR="$(dpkg-query -W -f='${Version}' source-extractor)"
GOT_PSFEX="$(dpkg-query -W -f='${Version}' psfex)"
GOT_STILTS="$(stilts -version 2>&1 | grep -oP 'STILTS version \K[0-9.-]+' || echo MISSING)"

FAIL=0
[ "$GOT_SEXTRACTOR" = "$WANT_SEXTRACTOR" ] || { echo "[bootstrap][FATAL] source-extractor $GOT_SEXTRACTOR != $WANT_SEXTRACTOR"; FAIL=1; }
[ "$GOT_PSFEX" = "$WANT_PSFEX" ] || { echo "[bootstrap][FATAL] psfex $GOT_PSFEX != $WANT_PSFEX"; FAIL=1; }
[ "$GOT_STILTS" = "$WANT_STILTS_VERSION" ] || { echo "[bootstrap][FATAL] stilts $GOT_STILTS != $WANT_STILTS_VERSION"; FAIL=1; }
if [ "$FAIL" -ne 0 ]; then
  echo "[bootstrap][FATAL] version mismatch against janne-pc -- refusing to" >&2
  echo "continue. Either the apt archive has moved on, or this VM's Debian" >&2
  echo "point release differs. Do not proceed without an explicit decision:" >&2
  echo "pin harder (apt-get install pkg=version=<archived build>), or accept" >&2
  echo "the drift and run the same output-parity check used for the EC2 arm" >&2
  echo "(diff a shared tile's SExtractor/PSFEx output against janne-pc's)." >&2
  exit 1
fi
echo "[bootstrap] Version check OK: source-extractor=$GOT_SEXTRACTOR psfex=$GOT_PSFEX stilts=$GOT_STILTS"

echo "[bootstrap] Verifying shipped Python 3.11 env..."
if [ ! -x "$PYTHON_BIN" ]; then
  echo "Expected python at $PYTHON_BIN not found -- did push_to_gcp.sh finish?" >&2
  exit 1
fi
"$PYTHON_BIN" -c "import galsim, astropy, numpy, pandas, pyarrow, astropy_healpix; print('python env ok:', galsim.__version__, astropy.__version__)"

echo "[bootstrap] Verifying repo layout..."
test -d "$REPO_DIR" || { echo "Expected repo at $REPO_DIR not found" >&2; exit 1; }
test -f "$REPO_DIR/data/plate_crpix_table.csv" || { echo "CRPIX table missing" >&2; exit 1; }
# .git/ is deliberately excluded from push_to_gcp.sh's rsync (no VCS history
# needed here), so there is no branch to check on this machine -- push_to_gcp.sh
# already verified janne-pc itself was on paper-parity before shipping, and
# steps 1-3 (all this VM runs) are identical between branches regardless.

echo "[bootstrap] All checks passed. Ready to fetch plates and slice."
