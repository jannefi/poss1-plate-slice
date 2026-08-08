#!/usr/bin/env python3
"""
EXPERIMENTAL — Post-pipeline pixel-level shape/PSF analysis stage (reference implementation parity).

Status
------
Not an official veto stage. Validate results before promoting to a hard gate.

Goal
----
Given a run directory containing a stage CSV, evaluate each candidate against
pixel-level radial-profile and contour metrics following the method described in:

    Busko (2026) – https://arxiv.org/abs/2603.20407
    Reference implementation: https://github.com/cuernodegazpacho/plateanalysis

Two classes of metrics are computed (reference implementation parity):

1. profile_diff (radial profile comparison)
   - Extract a neighborhood cutout (default 8') centred on the candidate.
   - Invert (65535 − pixel) then subtract photutils Background2D.
   - Build a RadialProfile for the candidate and each flux-matched neighbor star.
   - profile_diff = RMS of weighted difference between target and averaged star
     profiles (see referece implementation ProfileWorker for the exact weighting logic).
   Reject if profile_diff > --profile-diff-threshold (default 0.05).

2. Contour metrics (circularity, area, shape_defect, circle_deviation)
   - Extract a tiny cutout (default 21 px) centred on the candidate centroid.
   - Normalize to uint8, apply OpenCV thresholds [21, 45].
   - For each valid contour (area > 7, perimeter > 0): compute circularity,
     area, convexity-defect shape_defect, and circle_deviation.
   Reject if circularity < --circularity-low-limit (default 0.80).

   shape_confidence = "low" if area < 100, else "high". Never auto-drops.

Elongation gate
---------------
   Reject if elongation > --elongation-limit (default 1.10).
   Elongation is read from the per-tile SExtractor pass2 catalog.

Pixel source access
-------------------
For each candidate:
  - Locate tile FITS under <tiles-root>/<tile_id>/raw/*.fits
  - Locate pixel centroid in catalogs/sextractor_pass2.csv
    (prefer x_fit/y_fit; fallback to XWIN_IMAGE/YWIN_IMAGE)
  - Missing assets → shape_failed flag; row kept in flags output; no crash.

Neighbourhood star selection
-----------------------------
From the tile's pass2 catalog, for each candidate:
  a) Stars whose pixel coords fall inside the neighborhood cutout footprint.
  b) Stars with FLUX_MAX within ±flux_range (default 10 %) of the target.
  Self is excluded. Number of stars used is recorded in ledger.

Outputs (under <run-dir>/stages/)
-----------------------------------
1. stage_<STAGE>_SHAPE.csv
   Surviving candidates: src_id, tile_id, object_id, ra, dec, profile_diff

2. stage_<STAGE>_SHAPE_flags.csv
   All input rows with metrics, reason codes, and failure flags.

3. stage_<STAGE>_SHAPE_ledger.json
   Counts, parameters, neighborhood QA stats, and per-tile summary.

Usage
-----
python scripts/stage_shape_post.py \\
    --run-dir ./work/runs/run-R3-... \\
    --input-glob 'stages/stage_S3PTF.csv' \\
    --stage S4 \\
    --tiles-root ./data/tiles

Parallel execution (time-consuming for large N):
    --workers 8
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import multiprocessing
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cv2
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits as astropy_fits
from astropy.nddata import Cutout2D
from astropy.stats import SigmaClip
from astropy.wcs import WCS

from photutils.background import Background2D, MedianBackground
from photutils.profiles import RadialProfile


# ---------------------------------------------------------------------------
# Format helpers
# ---------------------------------------------------------------------------

_NAN = float("nan")


def _fmt(v) -> str:
    """Format a float for CSV output; empty string for NaN/None."""
    if v is None:
        return ""
    try:
        if math.isnan(float(v)):
            return ""
    except (TypeError, ValueError):
        return ""
    return f"{v:.6g}"


def _now_iso() -> str:
    return datetime.now(tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# CSV / stage helpers (matching stage_morph_post.py conventions)
# ---------------------------------------------------------------------------

def _read_header(path: Path) -> List[str]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.reader(f)
        return [c.strip().lstrip("\ufeff") for c in next(r, [])]


def _detect_src_cols(cols: List[str]) -> Tuple[str, str, str]:
    cset = set(cols)
    if "src_id" in cset:
        src = "src_id"
    elif "row_id" in cset:
        src = "row_id"
    else:
        raise RuntimeError("Input CSV missing required id column 'src_id' (or 'row_id')")
    ra = next((c for c in ["ra", "RA", "RA_corr", "ALPHAWIN_J2000"] if c in cset), None)
    dec = next((c for c in ["dec", "DEC", "Dec_corr", "DELTAWIN_J2000"] if c in cset), None)
    if not ra or not dec:
        raise RuntimeError("Input CSV missing RA/Dec columns")
    return src, ra, dec


def _parse_src_id(src_id: str) -> Tuple[str, str]:
    """Parse 'tile_id:object_id' → (tile_id, object_id)."""
    parts = src_id.split(":", 1)
    return (parts[0], parts[1]) if len(parts) == 2 else (src_id, "")


def _parse_edge_radii(spec: str) -> np.ndarray:
    """Parse --edge-radii argument.

    'arange30/2' → np.arange(30)/2 (default, reference implementation parity)
    '0,0.5,1,...' → np.array([0., 0.5, 1., ...])
    """
    if spec.strip().lower() == "arange30/2":
        return np.arange(30) / 2.0
    try:
        return np.array([float(x) for x in spec.split(",")])
    except ValueError:
        raise ValueError(f"Cannot parse --edge-radii: {spec!r}")


# ---------------------------------------------------------------------------
# Pass2 catalog loader
# ---------------------------------------------------------------------------

def _load_pass2(sex_csv: Path) -> Dict[str, dict]:
    """Load sextractor_pass2.csv keyed by NUMBER string.

    Returns {object_id: {x, y, flux_max, elongation, ra, dec}}.

    Pixel coords (x, y) are SExtractor 1-indexed XWIN/YWIN or x_fit/y_fit,
    converted to 0-indexed (subtract 1) for use with numpy/astropy.
    """
    result: Dict[str, dict] = {}
    try:
        with sex_csv.open(newline="", encoding="utf-8", errors="ignore") as fh:
            reader = csv.DictReader(fh)
            cols = set(reader.fieldnames or [])
            x_col  = next((c for c in ["x_fit", "X_FIT", "XWIN_IMAGE", "X_IMAGE"] if c in cols), None)
            y_col  = next((c for c in ["y_fit", "Y_FIT", "YWIN_IMAGE", "Y_IMAGE"] if c in cols), None)
            ra_col = next((c for c in ["ALPHAWIN_J2000", "ALPHA_J2000", "ra", "RA"] if c in cols), None)
            dec_col = next((c for c in ["DELTAWIN_J2000", "DELTA_J2000", "dec", "DEC"] if c in cols), None)
            flux_col  = next((c for c in ["FLUX_MAX", "flux_max", "FLUX_AUTO", "FLUX_APER"] if c in cols), None)
            elong_col = next((c for c in ["ELONGATION", "elongation"] if c in cols), None)

            for row in reader:
                num = (row.get("NUMBER") or row.get("number") or "").strip()
                if not num:
                    continue
                entry: dict = {}
                try:
                    if x_col and y_col:
                        # SExtractor pixel coords are 1-indexed; convert to 0-indexed
                        entry["x"] = float(row[x_col]) - 1.0
                        entry["y"] = float(row[y_col]) - 1.0
                    if ra_col and dec_col:
                        entry["ra"]  = float(row[ra_col])
                        entry["dec"] = float(row[dec_col])
                    if flux_col:
                        entry["flux_max"] = float(row[flux_col])
                    if elong_col:
                        entry["elongation"] = float(row[elong_col])
                    result[num] = entry
                except (ValueError, KeyError):
                    pass
    except Exception:
        pass
    return result


# ---------------------------------------------------------------------------
# Profile preconditioning (reference implementation parity)
# ---------------------------------------------------------------------------

def _precondition(data: np.ndarray, invert_max: float = 65535.0) -> np.ndarray:
    """Invert and background-subtract a cutout (reference implementation parity).

    Steps:
      1) invert: pixels = invert_max − pixels
      2) subtract Background2D(box_size=40, filter_size=3,
                               SigmaClip(sigma=3), MedianBackground)
    Falls back to inversion only if Background2D raises (e.g. cutout too small).
    """
    inverted = invert_max - data.astype(float)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bkg = Background2D(
                inverted,
                box_size=40,
                filter_size=3,
                sigma_clip=SigmaClip(sigma=3),
                bkg_estimator=MedianBackground(),
            )
        return inverted - bkg.background
    except Exception:
        return inverted


# ---------------------------------------------------------------------------
# Radial profile computation (reference implementation parity)
# ---------------------------------------------------------------------------

def _normalize_profile(profile: np.ndarray) -> np.ndarray:
    """Normalize profile to [0, 1]; returns zeros if flat."""
    pr_max = np.nanmax(profile)
    pr_min = np.nanmin(profile)
    denom = pr_max - pr_min
    if denom == 0.0:
        return np.zeros_like(profile, dtype=float)
    return (profile - pr_min) / denom


def _compute_profile_diff(
    cutout_data: np.ndarray,
    target_xy: Tuple[float, float],
    star_xys: List[Tuple[float, float]],
    edge_radii: np.ndarray,
) -> float:
    """Compute reference implementation profile_diff RMS.

    Returns NaN if target profile fails or no neighborhood star profiles succeed.

    Math (reference implementation ProfileWorker parity):
      averaged_profile = mean(normalized star profiles)
      diff = target_profile − averaged_profile
      diff = where(averaged_profile <= 0.1, 0, diff)   # mask low-signal bins
      diff *= averaged_profile                           # weight by profile
      diff[0:2] *= 0                                     # zero first two bins
      profile_diff = sqrt(sum(diff²) / len(diff))
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rp_target = RadialProfile(cutout_data, target_xy, edge_radii)
        target_prof = _normalize_profile(np.array(rp_target.profile, dtype=float))
    except Exception:
        return _NAN

    if not np.isfinite(target_prof).any():
        return _NAN

    star_profs = []
    for xy in star_xys:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rp = RadialProfile(cutout_data, xy, edge_radii)
            prof = _normalize_profile(np.array(rp.profile, dtype=float))
            if np.isfinite(prof).any():
                star_profs.append(prof)
        except Exception:
            continue

    if not star_profs:
        return _NAN

    averaged_profile = np.nanmean(np.array(star_profs), axis=0)
    diff = target_prof - averaged_profile
    diff = np.where(averaged_profile <= 0.1, 0.0, diff)
    diff *= averaged_profile
    diff[0:2] *= 0.0

    return float(np.sqrt(np.sum(np.square(diff)) / len(diff)))


# ---------------------------------------------------------------------------
# Contour / shape metrics (reference implementation parity)
# ---------------------------------------------------------------------------

def _compute_contour_metrics(
    data: np.ndarray,
    target_xy_orig: Tuple[float, float],
    tiny_cutout_size_px: int,
    opencv_thresholds: List[int],
) -> dict:
    """Compute circularity, area, shape_defect, circle_deviation on a tiny cutout.

    Follows reference implementation ContourWorker: loops over all thresholds and all valid
    contours; the last valid contour processed across all threshold iterations
    determines the stored values (parity with reference implementation).

    Returns dict with keys: circularity, area, shape_defect, circle_deviation,
    shape_confidence. Values are NaN where no valid contour was found.
    """
    out: dict = {
        "circularity": _NAN,
        "area": _NAN,
        "shape_defect": _NAN,
        "circle_deviation": _NAN,
        "shape_confidence": "low",
    }

    try:
        cutout_tiny = Cutout2D(data, position=target_xy_orig, size=tiny_cutout_size_px)
    except Exception:
        return out

    image_float = cutout_tiny.data.astype(float)
    image_uint8 = np.empty(image_float.shape, dtype=np.uint8)
    cv2.normalize(image_float, image_uint8, 0, 255, cv2.NORM_MINMAX, cv2.CV_8U)

    circularity = _NAN
    area = _NAN
    shape_defect = _NAN
    circle_deviation = _NAN

    for t in opencv_thresholds:
        _, thresh = cv2.threshold(image_uint8, t, 255, cv2.THRESH_BINARY)
        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        for contour in contours:
            ar = cv2.contourArea(contour)
            perimeter = cv2.arcLength(contour, True)

            if ar <= 7.0 or perimeter <= 0.0:
                continue

            # Circularity and area (reference implementation parity: last valid overwrites)
            circularity = (4.0 * math.pi * ar) / (perimeter ** 2)
            area = ar

            # Shape defect (cumulative convexity defect depth)
            epsilon = 0.01 * perimeter
            approx = cv2.approxPolyDP(contour, epsilon, True)
            if len(approx) >= 4:
                try:
                    hull = cv2.convexHull(approx, returnPoints=False)
                    defects = cv2.convexityDefects(approx, hull)
                    if defects is not None:
                        sum_depth = sum(defects[i, 0][3] / 256.0 for i in range(defects.shape[0]))
                        x_r, y_r, w_r, h_r = cv2.boundingRect(approx)
                        denom = max(w_r, h_r)
                        if denom > 0:
                            shape_defect = sum_depth / denom
                except cv2.error:
                    pass

            # Circle deviation (RMS of normalised contour-point distances)
            try:
                (x_c, y_c), radius = cv2.minEnclosingCircle(contour)
                if radius > 0:
                    distances = [
                        math.sqrt((pt[0][0] - x_c) ** 2 + (pt[0][1] - y_c) ** 2) / radius
                        for pt in contour
                    ]
                    if distances:
                        circle_deviation = float(np.std(distances))
            except Exception:
                pass

    out["circularity"] = circularity
    out["area"] = area
    out["shape_defect"] = shape_defect
    out["circle_deviation"] = circle_deviation
    out["shape_confidence"] = (
        "low" if (math.isnan(area) or area < 100.0) else "high"
    )
    return out


# ---------------------------------------------------------------------------
# Flags CSV field list
# ---------------------------------------------------------------------------

_FLAGS_FIELDS = [
    "src_id", "tile_id", "object_id", "ra", "dec",
    "profile_diff", "circularity", "area", "shape_defect", "circle_deviation",
    "shape_confidence", "elongation", "stars_used",
    "rej_profile_diff", "rej_elongation", "rej_circularity",
    "reject_flag", "reject_reason",
    "shape_failed", "failure_reason",
    "source_chunk",
]

_KEPT_FIELDS = ["src_id", "tile_id", "object_id", "ra", "dec", "profile_diff"]


def _make_failed_row(cand: dict, reason: str) -> dict:
    """Return a flags-row for a candidate that could not be evaluated."""
    return {
        "src_id": cand["src_id"],
        "tile_id": cand["tile_id"],
        "object_id": cand["object_id"],
        "ra": cand["ra"],
        "dec": cand["dec"],
        "profile_diff": "",
        "circularity": "",
        "area": "",
        "shape_defect": "",
        "circle_deviation": "",
        "shape_confidence": "",
        "elongation": "",
        "stars_used": "",
        "rej_profile_diff": "",
        "rej_elongation": "",
        "rej_circularity": "",
        "reject_flag": "0",
        "reject_reason": "",
        "shape_failed": "1",
        "failure_reason": reason,
        "source_chunk": cand.get("source_chunk", ""),
    }


# ---------------------------------------------------------------------------
# Per-tile worker (called from Pool.map)
# ---------------------------------------------------------------------------

def _process_tile(args: dict) -> List[dict]:
    """Process all candidates from a single tile. Returns list of row dicts."""
    tile_id: str = args["tile_id"]
    candidates: List[dict] = args["candidates"]
    tiles_root = Path(args["tiles_root"])
    params: dict = args["params"]

    tile_dir = tiles_root / tile_id

    # --- Locate FITS ---
    raw_dir = tile_dir / "raw"
    if not raw_dir.exists():
        return [_make_failed_row(c, "raw_dir_missing") for c in candidates]

    fits_files = sorted(raw_dir.glob("*.fits"))
    if not fits_files:
        return [_make_failed_row(c, "fits_not_found") for c in candidates]
    fits_path = fits_files[0]  # use first if multiple; logged in ledger

    # --- Load FITS ---
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with astropy_fits.open(fits_path) as hdul:
                hdr = hdul[0].header
                data = hdul[0].data.astype(float)
            wcs = WCS(hdr)
    except Exception as exc:
        reason = f"fits_load_error:{str(exc)[:80]}"
        return [_make_failed_row(c, reason) for c in candidates]

    # --- Load pass2 ---
    sex_csv = tile_dir / "catalogs" / "sextractor_pass2.csv"
    if not sex_csv.exists():
        return [_make_failed_row(c, "pass2_missing") for c in candidates]

    pass2 = _load_pass2(sex_csv)

    # Build flat list of candidate neighbours from pass2 (stars with full data)
    neigh_pool = [
        {"object_id": oid, **s}
        for oid, s in pass2.items()
        if all(k in s for k in ("x", "y", "flux_max"))
    ]

    edge_radii = params["edge_radii"]
    results: List[dict] = []

    for cand in candidates:
        obj_id = cand["object_id"]
        sex_entry = pass2.get(obj_id)

        if sex_entry is None:
            results.append(_make_failed_row(cand, "sex_row_missing"))
            continue

        if "x" not in sex_entry or "y" not in sex_entry:
            results.append(_make_failed_row(cand, "pixel_coords_missing"))
            continue

        x_px = sex_entry["x"]   # 0-indexed
        y_px = sex_entry["y"]   # 0-indexed
        target_flux_max = sex_entry.get("flux_max", _NAN)
        elongation = sex_entry.get("elongation", _NAN)

        try:
            cand_ra  = float(cand["ra"])
            cand_dec = float(cand["dec"])
        except (ValueError, TypeError):
            results.append(_make_failed_row(cand, "invalid_ra_dec"))
            continue

        # ---- Profile diff -----------------------------------------------
        profile_diff = _NAN
        stars_used = 0

        try:
            nbhd_arcmin = params["neighborhood_cutout_size_arcmin"]
            target_coord = SkyCoord(ra=cand_ra * u.deg, dec=cand_dec * u.deg)
            cutout_size = (nbhd_arcmin * u.arcmin, nbhd_arcmin * u.arcmin)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                nbhd_cutout = Cutout2D(
                    data, target_coord, cutout_size, wcs=wcs, mode="trim"
                )

            nbhd_data = _precondition(nbhd_cutout.data, invert_max=params["invert_max"])

            # Target position in cutout pixel coords
            ox, oy = nbhd_cutout.origin_original   # lower-left (x, y) in original
            tx = x_px - ox
            ty = y_px - oy

            # Select neighbourhood stars: inside cutout footprint + flux filter
            ch, cw = nbhd_cutout.data.shape
            flux_range = params["flux_range"]
            selected_star_xys: List[Tuple[float, float]] = []

            if not math.isnan(target_flux_max):
                flux_lo = target_flux_max * (1.0 - flux_range)
                flux_hi = target_flux_max * (1.0 + flux_range)
                for star in neigh_pool:
                    if star["object_id"] == obj_id:
                        continue  # skip self
                    sf = star["flux_max"]
                    if not (flux_lo <= sf <= flux_hi):
                        continue
                    sx_cut = star["x"] - ox
                    sy_cut = star["y"] - oy
                    if 0.0 <= sx_cut < cw and 0.0 <= sy_cut < ch:
                        selected_star_xys.append((sx_cut, sy_cut))

            stars_used = len(selected_star_xys)
            profile_diff = _compute_profile_diff(
                nbhd_data, (tx, ty), selected_star_xys, edge_radii
            )

        except Exception:
            profile_diff = _NAN

        # ---- Contour metrics --------------------------------------------
        contour = _compute_contour_metrics(
            data,
            (x_px, y_px),
            params["tiny_cutout_size_px"],
            params["opencv_thresholds"],
        )

        # ---- Gating ------------------------------------------------------
        rej_pd = int(
            not math.isnan(profile_diff)
            and profile_diff > params["profile_diff_threshold"]
        )
        rej_el = int(
            not math.isnan(elongation)
            and elongation > params["elongation_limit"]
        )
        circ = contour["circularity"]
        rej_ci = int(
            not math.isnan(circ)
            and circ < params["circularity_low_limit"]
        )

        reject_flag = int(rej_pd or rej_el or rej_ci)

        reasons = []
        if rej_pd:
            reasons.append(f"profile_diff={profile_diff:.4f}>{params['profile_diff_threshold']}")
        if rej_el:
            reasons.append(f"elongation={elongation:.3f}>{params['elongation_limit']}")
        if rej_ci:
            reasons.append(f"circularity={circ:.3f}<{params['circularity_low_limit']}")

        results.append({
            "src_id":            cand["src_id"],
            "tile_id":           tile_id,
            "object_id":         obj_id,
            "ra":                cand["ra"],
            "dec":               cand["dec"],
            "profile_diff":      _fmt(profile_diff),
            "circularity":       _fmt(contour["circularity"]),
            "area":              _fmt(contour["area"]),
            "shape_defect":      _fmt(contour["shape_defect"]),
            "circle_deviation":  _fmt(contour["circle_deviation"]),
            "shape_confidence":  contour["shape_confidence"],
            "elongation":        _fmt(elongation),
            "stars_used":        str(stars_used),
            "rej_profile_diff":  str(rej_pd),
            "rej_elongation":    str(rej_el),
            "rej_circularity":   str(rej_ci),
            "reject_flag":       str(reject_flag),
            "reject_reason":     ";".join(reasons),
            "shape_failed":      "0",
            "failure_reason":    "",
            "source_chunk":      cand.get("source_chunk", ""),
        })

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="[EXPERIMENTAL] Pixel-level shape/PSF analysis stage (reference implementation parity)."
    )
    ap.add_argument("--run-dir", required=True,
                    help="Run folder, e.g. ./work/runs/run-R3-...")
    ap.add_argument("--input-glob", default="stages/stage_S0M_MORPH.csv",
                    help="Glob (relative to run-dir) for input stage CSV.")
    ap.add_argument("--stage", default="S1",
                    help="Stage label for output filenames. Default: S1")
    ap.add_argument("--out-dir", default=None,
                    help="Output directory. Default: <run-dir>/stages")
    ap.add_argument("--tiles-root", default="./data/tiles",
                    help="Root of tile directories. Default: ./data/tiles")
    ap.add_argument("--workers", type=int, default=1,
                    help="Parallel worker processes (one per tile). Default: 1")
    ap.add_argument("--verbose", action="store_true",
                    help="Print per-tile progress.")

    # Neighbourhood / radial profile params
    ap.add_argument("--neighborhood-cutout-arcmin", type=float, default=8.0,
                    help="Neighbourhood cutout size in arcmin. Default: 8.0")
    ap.add_argument("--edge-radii", default="arange30/2",
                    help="Radial bin edges in px: 'arange30/2' or comma list. Default: arange30/2")
    ap.add_argument("--flux-range", type=float, default=0.1,
                    help="Fractional FLUX_MAX range for neighbour selection (±). Default: 0.1")
    ap.add_argument("--invert-max", type=float, default=65535.0,
                    help="Value used to invert pixel intensities. Default: 65535")

    # Contour params
    ap.add_argument("--tiny-cutout-px", type=int, default=21,
                    help="Tiny cutout size in pixels for contour metrics. Default: 21")
    ap.add_argument("--opencv-thresholds", default="21,45",
                    help="Comma-separated OpenCV threshold values. Default: 21,45")

    # Gating thresholds
    ap.add_argument("--profile-diff-threshold", type=float, default=0.05,
                    help="Reject if profile_diff > this. Default: 0.05")
    ap.add_argument("--elongation-limit", type=float, default=1.10,
                    help="Reject if elongation > this. Default: 1.10")
    ap.add_argument("--circularity-low-limit", type=float, default=0.80,
                    help="Reject if circularity < this. Default: 0.80")

    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise SystemExit(f"run-dir not found: {run_dir}")

    tiles_root = Path(args.tiles_root)
    out_dir = Path(args.out_dir) if args.out_dir else (run_dir / "stages")
    out_dir.mkdir(parents=True, exist_ok=True)

    chunks = sorted(run_dir.glob(args.input_glob))
    if not chunks:
        raise SystemExit(f"No inputs matched: {run_dir}/{args.input_glob}")

    src_col, ra_col, dec_col = _detect_src_cols(_read_header(chunks[0]))

    try:
        edge_radii = _parse_edge_radii(args.edge_radii)
    except ValueError as exc:
        raise SystemExit(str(exc))

    try:
        opencv_thresholds = [int(x.strip()) for x in args.opencv_thresholds.split(",")]
    except ValueError:
        raise SystemExit(f"Cannot parse --opencv-thresholds: {args.opencv_thresholds!r}")

    params = {
        "neighborhood_cutout_size_arcmin": args.neighborhood_cutout_arcmin,
        "edge_radii":                      edge_radii,
        "flux_range":                      args.flux_range,
        "invert_max":                      args.invert_max,
        "tiny_cutout_size_px":             args.tiny_cutout_px,
        "opencv_thresholds":               opencv_thresholds,
        "profile_diff_threshold":          args.profile_diff_threshold,
        "elongation_limit":                args.elongation_limit,
        "circularity_low_limit":           args.circularity_low_limit,
    }

    stage = args.stage
    out_kept   = out_dir / f"stage_{stage}_SHAPE.csv"
    out_flags  = out_dir / f"stage_{stage}_SHAPE_flags.csv"
    out_ledger = out_dir / f"stage_{stage}_SHAPE_ledger.json"

    # ------------------------------------------------------------------
    # Pass 1: read all input rows; group by tile
    # ------------------------------------------------------------------
    all_rows: List[dict] = []
    tiles_needed: Dict[str, List[dict]] = {}

    for ch in chunks:
        with ch.open(newline="", encoding="utf-8", errors="ignore") as fh:
            for row in csv.DictReader(fh):
                sid  = (row.get(src_col)  or "").strip()
                ra   = (row.get(ra_col)   or "").strip()
                dec  = (row.get(dec_col)  or "").strip()
                if not sid or not ra or not dec:
                    continue
                tile_id, object_id = _parse_src_id(sid)
                cand = {
                    "src_id": sid, "ra": ra, "dec": dec,
                    "tile_id": tile_id, "object_id": object_id,
                    "source_chunk": ch.name,
                }
                all_rows.append(cand)
                tiles_needed.setdefault(tile_id, []).append(cand)

    total_in = len(all_rows)
    print(f"[SHAPE] input_rows={total_in} tiles={len(tiles_needed)}")

    # ------------------------------------------------------------------
    # Pass 2: process tiles (parallel or serial)
    # ------------------------------------------------------------------
    tile_args = [
        {"tile_id": tid, "candidates": cands,
         "tiles_root": str(tiles_root), "params": params}
        for tid, cands in tiles_needed.items()
    ]

    if args.workers > 1:
        with multiprocessing.Pool(processes=args.workers) as pool:
            tile_results = pool.map(_process_tile, tile_args)
    else:
        tile_results = []
        for i, ta in enumerate(tile_args):
            res = _process_tile(ta)
            tile_results.append(res)
            if args.verbose:
                n_kept = sum(1 for r in res if r["reject_flag"] == "0" and r["shape_failed"] == "0")
                n_fail = sum(1 for r in res if r["shape_failed"] == "1")
                print(f"[SHAPE]   {ta['tile_id']}: in={len(res)} kept={n_kept} failed={n_fail}")

    evaluated = [r for tile_res in tile_results for r in tile_res]

    # ------------------------------------------------------------------
    # Compute totals
    # ------------------------------------------------------------------
    total_rejected = sum(1 for e in evaluated if e["reject_flag"] == "1")
    total_failed   = sum(1 for e in evaluated if e["shape_failed"] == "1")
    total_kept     = sum(1 for e in evaluated
                        if e["reject_flag"] == "0" and e["shape_failed"] == "0")

    rej_pd = sum(1 for e in evaluated if e.get("rej_profile_diff") == "1")
    rej_el = sum(1 for e in evaluated if e.get("rej_elongation")   == "1")
    rej_ci = sum(1 for e in evaluated if e.get("rej_circularity")  == "1")

    # Neighbourhood stars_used stats (excluding failed rows)
    stars_counts = [
        int(e["stars_used"])
        for e in evaluated
        if e.get("stars_used", "") not in ("", None)
    ]
    def _pct(arr, q):
        return float(np.percentile(arr, q)) if arr else _NAN

    nbhd_stats: dict = {}
    if stars_counts:
        zero_stars_rows = sum(1 for c in stars_counts if c == 0)
        nbhd_stats = {
            "min":              min(stars_counts),
            "p50":              float(np.median(stars_counts)),
            "max":              max(stars_counts),
            "zero_stars_rows":  zero_stars_rows,
        }

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    with out_kept.open("w", newline="", encoding="utf-8") as fk, \
         out_flags.open("w", newline="", encoding="utf-8") as ff:

        kept_w  = csv.DictWriter(fk, fieldnames=_KEPT_FIELDS)
        flags_w = csv.DictWriter(ff, fieldnames=_FLAGS_FIELDS)
        kept_w.writeheader()
        flags_w.writeheader()

        for e in evaluated:
            flags_w.writerow(e)
            if e["reject_flag"] == "0" and e["shape_failed"] == "0":
                kept_w.writerow({k: e[k] for k in _KEPT_FIELDS})

    # ------------------------------------------------------------------
    # Per-tile summary for ledger
    # ------------------------------------------------------------------
    per_tile = []
    for tile_id in sorted(tiles_needed.keys()):
        tm = [e for e in evaluated if e["tile_id"] == tile_id]
        t_kept = sum(1 for e in tm if e["reject_flag"] == "0" and e["shape_failed"] == "0")
        t_rej  = sum(1 for e in tm if e["reject_flag"] == "1")
        t_fail = sum(1 for e in tm if e["shape_failed"] == "1")
        failure_reasons = list({e["failure_reason"] for e in tm if e.get("failure_reason")})
        per_tile.append({
            "tile_id":         tile_id,
            "input_rows":      len(tm),
            "kept_rows":       t_kept,
            "rejected_rows":   t_rej,
            "failed_rows":     t_fail,
            "failure_reasons": failure_reasons,
        })

    ledger = {
        "experimental":   True,
        "created":        _now_iso(),
        "stage":          stage,
        "run_dir":        str(run_dir),
        "input_glob":     args.input_glob,
        "tiles_root":     str(tiles_root),
        "workers":        args.workers,
        "parameters": {
            "neighborhood_cutout_size_arcmin": args.neighborhood_cutout_arcmin,
            "edge_radii":                      args.edge_radii,
            "flux_range":                      args.flux_range,
            "invert_max":                      args.invert_max,
            "tiny_cutout_size_px":             args.tiny_cutout_px,
            "opencv_thresholds":               opencv_thresholds,
            "profile_diff_threshold":          args.profile_diff_threshold,
            "elongation_limit":                args.elongation_limit,
            "circularity_low_limit":           args.circularity_low_limit,
        },
        "input_chunks": [p.name for p in chunks],
        "totals": {
            "input_rows":        total_in,
            "kept_rows":         total_kept,
            "rejected_rows":     total_rejected,
            "failed_rows":       total_failed,
            "rejection_pct":     round(100.0 * total_rejected / total_in, 2) if total_in > 0 else 0.0,
        },
        "rejected_by_reason": {
            "rej_profile_diff": rej_pd,
            "rej_elongation":   rej_el,
            "rej_circularity":  rej_ci,
        },
        "neighborhood_stats": nbhd_stats,
        "columns_detected": {
            "src_id_col": src_col,
            "ra_col":     ra_col,
            "dec_col":    dec_col,
        },
        "outputs": {
            "kept_csv":    str(out_kept),
            "flags_csv":   str(out_flags),
            "ledger_json": str(out_ledger),
        },
        "per_tile": per_tile,
    }

    out_ledger.write_text(json.dumps(ledger, indent=2), encoding="utf-8")

    pct = 100.0 * total_rejected / total_in if total_in > 0 else 0.0
    print(
        f"[SHAPE] [EXPERIMENTAL] input={total_in} kept={total_kept} "
        f"rejected={total_rejected} failed={total_failed} ({pct:.1f}% rejection)"
    )
    print(f"[SHAPE] wrote: {out_kept}")
    print(f"[SHAPE] wrote: {out_flags}")
    print(f"[SHAPE] wrote: {out_ledger}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
