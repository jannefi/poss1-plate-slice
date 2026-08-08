# STAGE_SHAPE — Pixel-level shape/PSF analysis (EXPERIMENTAL)

> **Status: EXPERIMENTAL.** Do not treat rejection by this stage as definitive.
> Validate results manually before promoting to a hard gate.

## Overview

`scripts/stage_shape_post.py` implements pixel-level candidate filtering based
on the method described in:

> Busko (2026) – *Searching for Fast Astronomical Transients in Archival Photographic Plates*
> https://arxiv.org/abs/2603.20407
>
> Reference implementation:
> https://github.com/cuernodegazpacho/plateanalysis

For each candidate in the shrinking set the stage computes two classes of
metrics and applies optional rejection gates.

---

## Metrics

### 1. profile_diff (radial profile comparison)

A neighborhood cutout (default 8 arcmin) centred on the candidate is extracted
from the tile FITS. The pixels are preconditioned:

1. **Invert**: `pixel = 65535 − pixel`
2. **Background subtract**: `photutils.Background2D(box_size=40, filter_size=3,
   SigmaClip(σ=3), MedianBackground)`

A `photutils.RadialProfile` is built for the candidate and for each
flux-matched neighbor star (see *Neighbourhood selection* below) using edge
radii `np.arange(30)/2` (0 … 14.5 px in 0.5 px steps).

Each profile is normalized to [0, 1]. The **profile_diff** RMS is then:

```
averaged_profile = mean(normalized star profiles)
diff = target_profile − averaged_profile
diff = where(averaged_profile ≤ 0.1,  0, diff)   # mask low-signal bins
diff *= averaged_profile                           # weight by profile shape
diff[0:2] *= 0                                     # zero first two bins
profile_diff = sqrt(sum(diff²) / len(diff))
```

**Reject if** `profile_diff > --profile-diff-threshold` (default **0.05**).

### 2. Contour metrics

A tiny cutout (default 21 px) centred on the best pixel centroid is extracted.
It is normalized to uint8 [0–255] and passed through OpenCV thresholds
(default **21** and **45**). For each valid contour (area > 7, perimeter > 0):

| Metric | Formula |
|---|---|
| `circularity` | 4π·A / P² |
| `area` | Contour area A (px²) |
| `shape_defect` | Σ convexity-defect depths / max(w, h) of bounding rect |
| `circle_deviation` | std(distance_to_centre / radius) over contour points |

*Reference implementation parity*: the last valid contour across all threshold iterations
overwrites the stored values (matching the reference implementation).

`shape_confidence = "low"` if `area < 100`, else `"high"`. Low confidence
rows are **never auto-dropped** — the explicit thresholds decide.

**Reject if** `circularity < --circularity-low-limit` (default **0.80**).

### 3. Elongation gate

**Reject if** `ELONGATION > --elongation-limit` (default **1.10**).
Elongation is read from the per-tile `sextractor_pass2.csv`.

---

## Neighbourhood selection

For each candidate, flux-matched neighbour stars are drawn from the tile's
`catalogs/sextractor_pass2.csv`:

1. Star pixel coords (0-indexed) fall inside the neighborhood cutout footprint.
2. `|FLUX_MAX_star − FLUX_MAX_target| / FLUX_MAX_target ≤ flux_range`
   (default `flux_range = 0.1`, i.e. ±10 %).

The number of stars used per candidate is recorded in the ledger.

---

## Pixel coordinate source

For the candidate centroid, the script reads from `sextractor_pass2.csv` in
order of preference:

1. `x_fit` / `y_fit`
2. `XWIN_IMAGE` / `YWIN_IMAGE`

SExtractor pixel coordinates are 1-indexed; the script converts to 0-indexed
for use with numpy/astropy.

---

## Graceful degradation

If a tile asset is missing (FITS, pass2 catalog, pixel coords, invalid ra/dec),
the candidate is marked `shape_failed=1` with a `failure_reason` code. It is
kept in the flags output and counted in the ledger. The stage never crashes on
a single bad tile.

---

## Outputs

All outputs are written to `<run-dir>/stages/`.

| File | Description |
|---|---|
| `stage_{STAGE}_SHAPE.csv` | Survivors: `src_id, tile_id, object_id, ra, dec, profile_diff` |
| `stage_{STAGE}_SHAPE_flags.csv` | All rows with metrics + rejection + failure flags |
| `stage_{STAGE}_SHAPE_ledger.json` | Counts, parameters, neighbourhood QA stats, per-tile summary |

The `ra` and `dec` columns are passed through unchanged from the input stage
CSV (WCS-corrected coordinates). No `NUMBER`/`number` column is emitted.

---

## Usage

```bash
python scripts/stage_shape_post.py \
    --run-dir ./work/runs/run-R3-... \
    --input-glob 'stages/stage_S3PTF.csv' \
    --stage S4 \
    --tiles-root ./data/tiles
```

Parallel execution (one worker per tile, recommended for large runs):

```bash
python scripts/stage_shape_post.py \
    --run-dir ./work/runs/run-R3-... \
    --input-glob 'stages/stage_S3PTF.csv' \
    --stage S4 \
    --tiles-root ./data/tiles \
    --workers 8
```

### Full CLI reference

```
--run-dir                     Run folder (required)
--input-glob                  Glob for input CSV relative to run-dir
--stage                       Stage label (e.g. S4)
--out-dir                     Output directory (default: <run-dir>/stages)
--tiles-root                  Root of tile directories (default: ./data/tiles)
--workers                     Parallel workers (default: 1)
--verbose                     Print per-tile progress

Neighbourhood / radial profile:
--neighborhood-cutout-arcmin  Cutout size in arcmin (default: 8.0)
--edge-radii                  Radial bin edges: 'arange30/2' or CSV floats (default: arange30/2)
--flux-range                  Fractional FLUX_MAX tolerance for neighbour selection (default: 0.1)
--invert-max                  Pixel inversion value (default: 65535)

Contour:
--tiny-cutout-px              Tiny cutout size in px (default: 21)
--opencv-thresholds           Comma-separated OpenCV thresholds (default: 21,45)

Gating:
--profile-diff-threshold      Reject if profile_diff > this (default: 0.05)
--elongation-limit            Reject if elongation > this (default: 1.10)
--circularity-low-limit       Reject if circularity < this (default: 0.80)
```

---

## Acceptance criteria / test plan

1. Outputs (`_SHAPE.csv`, `_SHAPE_flags.csv`, `_SHAPE_ledger.json`) appear under
   `<run-dir>/stages/`.
2. `ledger.totals.input_rows` == number of rows in the input glob.
3. `kept_rows + rejected_rows + failed_rows == input_rows`.
4. No `NUMBER`/`number` column in any output CSV.
5. Spot-check: pick one candidate, confirm `profile_diff` is a non-NaN float
   in a reasonable range [0 … ~1].
6. Graceful degradation: point `--tiles-root` at a partial tile set; all
   missing-asset rows appear in flags with `shape_failed=1`.

---

## Dependencies

| Library | Purpose |
|---|---|
| `astropy` | FITS I/O, WCS, Cutout2D, SigmaClip |
| `photutils` | Background2D, RadialProfile |
| `opencv-python` (`cv2`) | Contour metrics |
| `numpy` | Array math |

Install if missing:

```bash
pip install photutils opencv-python
```
