# Reproducing this work

End to end, from empty disk to a validated catalogue. Read
[`PARAMETERS.md`](PARAMETERS.md) alongside this — it says which numbers you may
change and which you may not.

## 0. What you are committing to

| resource | requirement |
|---|---|
| plate scans | ~1 GB per plate; 634 plates for the full footprint |
| slice run | ~3 h at 12 workers, whole survey |
| peak disk during slicing | one plate's tiles, a few GB |
| lean RA/Dec output | ~14 GB for the full survey |
| local catalogue mirrors | large, and see the README's warning about bulk access |

**Disk discipline is load-bearing, not tidiness.** `run_fullscale_slice.py` keeps
exactly one plate's tile tree on disk at a time: slice, extract, write the lean
RA/Dec CSV, delete. Retaining every tile tree through the veto stages would need
roughly **7.4 TB**. If you pass `--keep-tiles`, size your storage first.

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
named `dss1red_XE*.fits`. The footprint used here is the 634-plate list in
`data/plate_manifest.csv`.

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
export VASCO_GAIA_CACHE=... VASCO_PS1_CACHE=... VASCO_USNOB_CACHE=...

python3 tools/run_fullscale_slice.py \
    --out-dir work/slice \
    --plate-manifest data/plate_manifest.csv \
    --crpix-table data/plate_crpix_table.csv \
    --workers 12
```

**Resumable** — a plate whose output CSV exists is skipped, so an interrupted run
continues where it stopped.

Check `work/slice/progress.csv` when it finishes: every plate must report **49
tiles sliced, 49 with catalogues, and 0 skips**. A plate reporting fewer has
failed; the run does not stop for it, because one bad plate should not cost the
other 633. Per-plate slicer output is kept under `work/slice/slice_logs/`.

> An earlier version of this runner discarded the slicer's stdout, and 195 tiles
> lost to a grid walking off the array went unnoticed for weeks. That is why the
> skip count is in the ledger and the logs are kept.

## 6. Veto and filter

Steps 4-5 apply the cross-match vetoes and the candidate-selection filters.

```bash
python3 tools/run_steps_4_5_parallel.py --tiles-file <tiles> --workers 12
```

Set the cache environment variables on **every** step, not just the first — a
missing one silently falls back to live queries rather than failing.

## 7. Post-processing

All 16 stages are in `scripts/stage_*_post.py`; chain them with
`scripts/run_post_stage_chain.sh`. `RESULTS.md` records the chain used for the
published numbers.

**Three stages ship but were not run for those numbers** — SkyBoT, SuperCOSMOS
and VSX. See [`PARAMETERS.md`](PARAMETERS.md); this materially affects the final
row count and comparability.

Deduplication tolerance is **3.0″** under raw plate WCS and is *coupled* to
whether WCS-fix is applied. Read that section before changing either.

## 8. Validate

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
