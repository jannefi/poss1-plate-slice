# Reproducing this work

End to end, from empty disk to a validated catalogue. Read
[`PARAMETERS.md`](PARAMETERS.md) alongside this — it says which numbers you may
change and which you may not.

## 0. What you are committing to

Budget **~700 GB of input data** before you start. Every figure below is measured
from a real full-survey run over the 642-plate footprint, not estimated.

### Input data, downloaded once

| item | size | notes |
|---|---:|---|
| DSS1-red full plate scans, 642 plates | **252 GB** | uniformly 0.39 GB each |
| STScI cutouts, one per plate | **~6 GB** | 9 MB per 60′ cutout; needed only for the CRPIX table in step 4 |
| Gaia mirror | **58 GB** | |
| Pan-STARRS1 mirror | **308 GB** | by far the largest single item |
| USNO-B1.0 mirror | **61 GB** | |
| **total** | **~685 GB** | |

The three mirrors are 427 GB of that. They are optional in the sense that the
code falls back to live queries, and impractical to omit in practice — see §2.
To reproduce only the **detection** numbers and not the veto chain, the plate
scans and cutouts alone suffice: **~258 GB**.

The plate-scan figure is for this footprint specifically. IRSA's `dss1red`
library holds 932 scans totalling 443 GB, averaging 0.48 GB and reaching 1.06 GB
— but every plate in the 642 is 0.39 GB, so sizing from the library average
over-provisions by about a fifth.

### Working space and outputs

| item | size | notes |
|---|---:|---|
| lean RA/Dec CSVs, whole survey | **~14 GB** | 642 files; the detection-level product |
| survivors tree after the veto chain | **<1 GB** | |
| peak scratch, one plate in flight | **1–8 GB** | 49 tiles, ~30 MB each mid-run and ~150 MB fully populated |
| every tile tree retained (`--keep-tiles`) | **~5 TB** | 31,066 tiles × ~150 MB |

**Disk discipline is load-bearing, not tidiness.** `run_fullscale_slice.py` keeps
exactly one plate's tile tree on disk at a time: slice, extract, write the lean
RA/Dec CSV, delete. That is what holds the peak at single-digit GB instead of
~5 TB. If you pass `--keep-tiles`, size your storage first.

Setting `VASCO_LDAC_DROP_VIGNET=1` drops the `VIGNET` postage stamps during the
LDAC→CSV conversion and is worth it at survey scale; it changes no column the
pipeline consumes.

### Time

| run | cost at 12 workers |
|---|---|
| detection only, no vetoes | **~3 h** |
| full chain, `--with-vetoes` | **~4.7 days**, ~630 s/plate |

Runtime is set by **sky density, not I/O**: `corr(|galactic b|, detections)` is
−0.89, and plates run in `XE` order, which sweeps the galactic plane — so the
printed ETA drifts upward mid-plane and back down afterwards. That is expected,
not a stall. Both runs are **resumable** per plate, so an interruption costs only
the plate in flight.

## 1. Software

Python 3.11+, plus **SExtractor**, **PSFEx** and **STILTS** on `PATH`.

```bash
pip install -r requirements.txt
sex --version && psfex --version && stilts -version
```

## 2. Data

### Plate scans (required)

```
https://irsa.ipac.caltech.edu/data/DSS/images/dss1red/
```

Plate-addressed, so you fetch exactly the plates you intend to use. Files are
named `dss1red_XE*.fits`. The footprint used here is the 642-plate list in
`data/plate_manifest.csv`: every POSS-I red plate whose centre declination is
≥ −3.0°. That list is not hand-curated — regenerate it from the headers with

```bash
python3 tools/build_plate_manifest.py \
    --plate-dir <plate_dir> --out data/plate_manifest.csv
```

The threshold is not a tuned parameter: no plate centre lies between −5° and
−1°, so any cut inside that gap selects the same 642 plates, and the generator
asserts that gap is empty rather than trusting it.

### Archive cutouts (required)

From the STScI/MAST cutout service. Needed for the union arm — and needed
regardless of that, because **cutout headers are the only place the explicit FITS
WCS solution appears**, which step 4 depends on. One cutout per plate is enough
for the correction table.

### Catalogue mirrors (strongly recommended)

```bash
ls scripts/local_cache/         # gaia/  ps1/  usnob/
```

Without them the veto stages hit live VizieR/MAST and you will be rate-limited
long before a survey-scale run finishes.

> The mirrors are **HEALPix `nside=32`, NESTED** ordering. Reading them as RING
> silently returns the wrong sky with no error and no warning — the query
> succeeds and the answer is nonsense. If your cross-match results look
> inexplicable, check this first.

## 3. Configure

```bash
cp config.yaml config.local.yaml     # edit the paths; config.local.yaml is gitignored
python3 -c "from vasco.paths import dump; dump()"
```

Every key must print `OK`. A missing path fails loudly rather than resolving to
someone else's filesystem.

## 4. Build the astrometric correction

**Do this before slicing.** Skipping it leaves a ~2.3″ systematic on a third of
the sky — see [`DSS_WCS_TWO_SOLUTIONS.md`](DSS_WCS_TWO_SOLUTIONS.md).

A prebuilt table for the 642-plate footprint ships as
`data/plate_crpix_table.csv`, so you can skip ahead and slice with it as is.
Rebuilding is how you *verify* it rather than trust it, and that needs one
archive cutout per plate — worth doing before relying on the numbers, since the
correction is derived from headers you can check yourself.

```bash
python3 tools/build_plate_crpix_table.py \
    --plate-dir <plate_dir> --archive-tiles <tiles_dir> \
    --out data/plate_crpix_table.csv
```

Expect roughly **a third of plates** to need correction, at a median ~2.3″, and
the rest to sit near 0.02″. Plates reported with scatter above 0.2 px are
corrected but flagged; a handful is normal, many is a sign your cutout set is
mismatched to your plate set.

## 5. Slice the survey

```bash
# Local mirrors. Optional -- without them the pipeline queries VizieR/MAST live.
export VASCO_GAIA_CACHE=... VASCO_PS1_CACHE=... VASCO_USNOB_CACHE=...

# REQUIRED. WCSFIX refits each tile's astrometry against Gaia and defaults to
# ON, but this pipeline is specified on the raw plate WCS, and the 3.0" dedup in
# step 7 is only correct for raw coordinates -- WCS-fixed coordinates need 0.25".
# Leaving this unset silently produces a different catalogue.
export VASCO_WCSFIX_DISABLE=1

# Must NOT be set: tiles are square here, with no circular catalogue cut.
unset VASCO_CIRCLE_ARCMIN

python3 tools/run_fullscale_slice.py \
    --out-dir work/slice \
    --plate-manifest data/plate_manifest.csv \
    --crpix-table data/plate_crpix_table.csv \
    --workers 12
```

**Read the first log lines before letting it run.** They must say
`circle_cut=off  wcsfix=off`. Both have gone wrong in real runs here: an
inherited `VASCO_CIRCLE_ARCMIN` cut 106 plates to a 30′ circle, and a missing
`VASCO_WCSFIX_DISABLE` ran an entire 642-plate campaign WCS-fixed. Neither
raised an error — the numbers were simply different. That is what the banner is
for. `tools/run_slice_survey.sh` sets every one of these explicitly and is the
safer way to launch.

**Resumable** — a plate whose output CSV exists is skipped, so an interrupted run
continues where it stopped.

Check `work/slice/progress.csv` when it finishes: every plate must report **49
tiles sliced, 49 with catalogues, and 0 skips**. A plate reporting fewer has
failed; the run does not stop for it, because one bad plate should not cost the
other 641. Per-plate slicer output is kept under `work/slice/slice_logs/`.

> An earlier version of this runner discarded the slicer's stdout, and 195 tiles
> lost to a grid walking off the array went unnoticed for weeks. That is why the
> skip count is in the ledger and the logs are kept.

## 6. Veto and filter

Steps 4-5 apply the cross-match vetoes and the candidate-selection filters.

```bash
python3 tools/run_steps_4_5_parallel.py --tiles-file <tiles> --workers 12
```

Set the environment on **every** step, not just the first — `VASCO_WCSFIX_DISABLE`
included, since steps 4-5 are where WCSFIX actually runs. A missing cache
variable silently falls back to live queries rather than failing.

**Verify the vetoes actually covered the sky.** They can be incomplete without
raising anything, so check the result rather than the logs:

```bash
VASCO_GAIA_CACHE=... python3 tools/check_s0_gaia_invariant.py \
    --s0-csv work/slice/runs/<run>/stage_S0.csv
```

Every surviving row passed a Gaia veto, so almost none of them should have a
Gaia source within 5″. If the figure is not ~0, the veto did not see that sky —
whatever the logs say. This check exists because a partial-cone bug in the local
mirror query once left 56% of a catalogue as un-vetoed Gaia stars, with every
stage reporting ok. Run it before trusting any candidate count.

### If that check fails, re-run — do not try to repair the catalogue

It is tempting to fix a bad veto by re-applying the corrected one to the
survivors you already have. **That does not work, and the reason is worth
understanding before you change anything upstream of the filters.**

The morphology gate in `vasco/mnras/filters_mnras.py` applies a 2σ clip to
`FWHM_IMAGE` and `ELONGATION` whose window is `median ± k·1.4826·MAD` **of the
population being filtered** — not a fixed threshold. So the rows that reach the
gate determine which rows the gate keeps. A remainder still full of un-vetoed
stars is tightly clustered in FWHM, which shrinks the MAD, narrows the window,
and cuts genuinely broader objects that a correct run would keep.

The consequence: a veto failure does not merely *add* spurious survivors, it
also *removes* real ones, so the correct catalogue is **not** a subset of the
broken one and cannot be recovered from it. Measured on one plate here, a
post-hoc re-veto recovered 80 of the 84 true survivors; the 4 it could not
recover had passed the veto in both runs and differed only through the clip
window.

Two practical rules follow:

* **Retain `<tile>/catalogs/sextractor_pass2.csv`.** It is ~5 MB per tile
  (~158 GB for a full survey) and it is the last point at which the expensive
  imaging work — slicing, SExtractor, PSFEx — is cheaply reversible. Everything
  downstream of it (veto, filters, spike mask, dedup) can be redone in hours
  *given the detections*. Delete the FITS slices, the LDACs and the PSFs if you
  need the space; keep this file.
* Treat the filter stage as population-dependent whenever you reason about what
  a change upstream of it will do.

## 7. Aggregate the survivors into S0

The veto chain leaves one survivors CSV per tile. Downstream stages need those
collected into a single set, with each row carrying the plate it came from.

Build the tile → plate map first. On this route the runner sliced a plate you
named and stamped `plate_id` onto every row it emitted, so its own output is the
record of what happened:

```bash
python3 tools/build_tile_plate_map_from_radec.py \
    --radec-dir work/slice/radec \
    --out-csv work/slice/tile_plate_map.csv \
    --filtered-dir work/slice/filtered
```

It aborts if any tile maps to two plates, which would mean the run mixed
sources. Do **not** substitute `tools/build_tile_plate_map_from_headers.py`
here — that one resolves a genuine ambiguity in the *archive* route, where a
cutout service picks the plate for you, and it needs tile FITS that this route
deletes as it goes.

Then aggregate and deduplicate:

```bash
python3 scripts/build_run_stage_csvs.py \
    --tiles-root work/slice/filtered \
    --plate-map-csv work/slice/tile_plate_map.csv \
    --catalog-name catalogs/sextractor_pass2.filtered.csv \
    --run-root work/runs --run-tag <tag> \
    --dedup-tol-arcsec 3.0 --accept-empty-file-as-valid-empty --full
```

`--accept-empty-file-as-valid-empty` matters more than it sounds: a tile with no
survivors is written as a 0-byte file, and without the flag those count as
failures rather than as legitimately empty. In the 642-plate run that is 5,063 of
31,458 tiles — 16% of the survey.

**Take the footprint from the run's `tile_manifest.csv`, never from the tile_ids
present in a stage CSV.** Tiles that produced no survivor do not appear in S0 at
all, so deriving coverage from S0 silently shrinks the denominator and inflates
any per-area rate computed from it.

## 8. Post-processing

All 16 stages are in `scripts/stage_*_post.py`; chain them with
`scripts/run_post_stage_chain.sh`. `RESULTS.md` records the chain used for the
published numbers.

**Three stages ship but were not run for those numbers** — SkyBoT, SuperCOSMOS
and VSX. See [`PARAMETERS.md`](PARAMETERS.md); this materially affects the final
row count and comparability.

Deduplication tolerance is **3.0″** under raw plate WCS and is *coupled* to
whether WCS-fix is applied. Read that section before changing either.

## 9. Validate

Two checks, both of which should be run before believing any output.

```bash
python3 tools/audit_independence.py                       # no private data, tree + history
python3 tools/union_parity_fullscale.py \
    --ref-csv <reference.csv> \
    --arm fullplate=work/slice/radec --arm archive=<archive_radec> \
    --combine fullplate+archive --out-dir work/union
```

Sanity checks on the parity output:

- every arm must report a **plausible detection count** in its `[ARM]` line. An
  arm that silently reads zero files reports 0% and looks like a result.
- watch the **1-3″ columns**, not only 5″. The astrometric correction lives
  there; a regression shows up at 1-2″ long before it is visible at 5″. That is
  precisely how the original defect hid for so long.

## Common failure modes

| symptom | cause |
|---|---|
| cross-match returns nonsense, no error | mirror read as RING instead of NESTED |
| veto stage inexplicably slow | cache env var unset on that step; it is querying live |
| everything rejected by the morphology gate | single-pass mode — no PSF model, so no `SPREAD_MODEL` |
| tile count below 49 for a plate | grid walking off the array; check `slice_logs/` |
| ~2.3″ offset on a third of plates | correction table not built, or not passed to the slicer |
| duplicate-inflated candidate count | 0.25″ dedup applied to non-WCS-fixed coordinates |
