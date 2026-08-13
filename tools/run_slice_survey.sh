#!/usr/bin/env bash
# Launch the full slice survey with every science-affecting variable set
# explicitly, rather than inherited or left to a default.
#
# Usage:
#   tools/run_slice_survey.sh --plate-dir <dir> --out-dir <dir> [extra args...]
#
# Any extra arguments are passed through to tools/run_fullscale_slice.py, so
# --plates, --workers, --grid and friends work as normal.
#
# Why this script exists. Two variables have silently changed the science in
# real runs of this pipeline, neither raising an error:
#
#   VASCO_CIRCLE_ARCMIN   inherited from an interactive shell, it applied a 30'
#                         circular cut to 106 plates of a full-scale run,
#                         discarding ~21% of detections and reintroducing the
#                         corner gaps that square tiles exist to avoid.
#
#   VASCO_WCSFIX_DISABLE  never set, so WCSFIX -- which defaults to ON -- refit
#                         every tile's astrometry against Gaia for an entire
#                         642-plate campaign. This pipeline is specified on the
#                         raw plate WCS, and the 3.0" dedup applied downstream
#                         is only correct for raw coordinates; WCS-fixed
#                         coordinates require 0.25".
#
# Setting a variable in a launcher is cheap. Discovering months later that a
# survey ran in an undocumented configuration is not.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- explicitly NOT set ---------------------------------------------------
unset VASCO_CIRCLE_ARCMIN        # square tiles, no 30' catalogue cut
unset VASCO_REPRO_SINGLE_PASS    # two-pass SExtractor: single-pass has no
                                 # SPREAD_MODEL, so the morphology gate would
                                 # reject every candidate

# --- explicitly set -------------------------------------------------------
export VASCO_WCSFIX_DISABLE=1    # raw plate WCS; see above
export VASCO_LDAC_DROP_VIGNET=1  # drop the VIGNET column during LDAC->CSV;
                                 # ~6.6x smaller intermediates, no science change

# Local mirrors are optional -- unset, the pipeline queries VizieR/MAST live.
# Set them here if you built mirrors, so every stage sees them: a variable set
# for step 2 but not step 4 falls back to live queries without complaining.
: "${VASCO_GAIA_CACHE:=}"
: "${VASCO_PS1_CACHE:=}"
: "${VASCO_USNOB_CACHE:=}"
[ -n "$VASCO_GAIA_CACHE"  ] && export VASCO_GAIA_CACHE
[ -n "$VASCO_PS1_CACHE"   ] && export VASCO_PS1_CACHE
[ -n "$VASCO_USNOB_CACHE" ] && export VASCO_USNOB_CACHE

cd "$REPO" || exit 1
echo "[LAUNCH] $(date -u +%Y-%m-%dT%H:%M:%SZ)  wcsfix=off circle_cut=off"
echo "[LAUNCH] verify the runner's own [CONFIG] line agrees before walking away"
exec python3 tools/run_fullscale_slice.py "$@"
