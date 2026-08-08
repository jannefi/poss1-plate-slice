# STAGE_MORPH — Morphology-Based Filtering

**Status:** EXPERIMENTAL — validate rejection rates before promoting to a hard gate.

## Purpose

`stage_morph_post.py` is a run-scoped post-pipeline stage that rejects candidates whose morphology is inconsistent with the local PSF. It is designed to run **first** in the post-pipeline stage chain (before SKYBOT, SCOS, PTF, VSX, GSC) to shrink the candidate set early and reduce load on network-bound stages.

Motivation: POSS-I plates contain photographic artifacts that survive the initial SExtractor gates (FLAGS, SNR, ELONGATION, SPREAD_MODEL thresholds). These artifacts are detectable by comparing each candidate's profile to nearby Gaia-matched reference stars, which characterise the local PSF.

## Method

### PSF reference sample (per tile)

Built from `catalogs/sextractor_pass2.csv` (the full detection catalog, not just surviving candidates):

- FLAGS = 0
- ELONGATION < 1.3
- Positionally matched to a Gaia DR3 star within 3″
- Gaia Gmag in (12, 18) — avoids saturated bright stars and faint low-SNR sources

Typical yield: ~1500 reference stars per tile.

### Rejection criteria

Two independent metrics; a candidate is rejected if **either** fires:

| Metric | Formula | Default threshold | What it catches |
|--------|---------|-------------------|-----------------|
| `fwhm_ratio` | `FWHM_IMAGE / psf_fwhm_median` | > 1.5 | Extended blobs, halos, plate scratches |
| `spread_snr` | `(SPREAD_MODEL − psf_spread_median) / SPREADERR_MODEL` | > 5.0 | Profile deviations vs. local PSFEx model |

**Note:** `CLASS_STAR` is intentionally not used. On photographic (POSS-I) plates, the SExtractor neural-net classifier is unreliable — PSF reference stars themselves score ~0.015 on a 0–1 scale (trained on CCD data, not photographic plates).

### Pass-through cases (no rejection applied)

- Tile has fewer than `--min-psf-stars` (default 5) reference stars → flagged `psf_insufficient`
- Tile catalogs not found on disk → flagged `catalog_missing`
- Candidate NUMBER not found in tile catalog → flagged `candidate_not_found`

All pass-through cases appear in the flags CSV with `reject_flag=0`.

## Calibration

Measured on a 181-tile, 684-candidate run (this dataset):

| Criterion | Rejection |
|-----------|-----------|
| `fwhm_ratio > 1.5` | 17.8% |
| `spread_snr > 5.0` | 49.1% |
| Either (default config) | **50.3%** |

The `spread_snr` threshold is the dominant discriminator. At `spread_snr > 3`, rejection rises to ~54%.

## Usage

```bash
python scripts/stage_morph_post.py \
    --run-dir ./work/runs/run-S1-... \
    --input-glob 'stage_S0.csv' \
    --stage S0M \
    --tiles-root ./data/tiles
```

Feed the output into the rest of the stage chain:

```bash
# After morph:
python scripts/stage_skybot_post.py \
    --run-dir ./work/runs/run-S1-... \
    --input-glob 'stages/stage_S0M_MORPH.csv' \
    --stage S1
```

## Outputs

All outputs land in `<run-dir>/stages/` (override with `--out-dir`).

| File | Content |
|------|---------|
| `stage_S0M_MORPH.csv` | Kept candidates: `src_id, ra, dec` |
| `stage_S0M_MORPH_flags.csv` | All input rows + metrics + reject reason |
| `stage_S0M_MORPH_ledger.json` | Counts, per-tile PSF stats, parameters |

### Flags CSV columns

```
src_id, ra, dec, tile_id, object_id,
fwhm_image, fwhm_ratio,
spread_model, spreaderr_model, spread_snr,
psf_fwhm_median, psf_spread_median, psf_star_count,
reject_flag, reject_reason, source_chunk
```

## CLI reference

```
--run-dir            REQUIRED  Run folder
--input-glob         stage_S0.csv
--stage              S0M
--out-dir            <run-dir>/stages
--tiles-root         ./data/tiles
--verbose            Print per-tile progress

PSF star selection:
--gaia-mag-min       12.0
--gaia-mag-max       18.0
--gaia-match-arcsec  3.0
--elongation-max     1.3
--min-psf-stars      5

Rejection thresholds:
--fwhm-ratio-max     1.5
--spread-snr-max     5.0
```

## Validation

After running, inspect:

1. `ledger.json → totals.rejection_pct` — expect ~50% on a typical polar tile run
2. `ledger.json → tile_psf_summary` — all tiles should show `ok`; investigate any `psf_insufficient`
3. Manually review a sample of `reject_reason` values in the flags CSV
4. Check that `kept + rejected == input` (i.e., no rows dropped silently)

```bash
# Quick check: count rows
wc -l stages/stage_S0M_MORPH.csv stages/stage_S0M_MORPH_flags.csv
```

## Performance

The Gaia proximity lookup uses binary search on declination (O(log N + K) per source, where K ≈ 3–4 stars/query). On a typical tile (3000 SExtractor sources, 4000 Gaia stars), PSF construction takes ~1s. Full 181-tile run completes in ~3 minutes.
