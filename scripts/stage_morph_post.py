#!/usr/bin/env python3
"""
EXPERIMENTAL — Post-pipeline morphology-based filtering stage.

Status
------
Not an official veto stage. Validate results before promoting to a hard gate.

Goal
----
Given a run directory containing a stage CSV (e.g. stage_S0.csv), evaluate
each candidate against a per-tile PSF model derived from Gaia-matched
SExtractor sources. Candidates significantly more extended than the local PSF
are rejected. Designed to run first in the post-pipeline stage chain to shrink
the candidate set before network-bound stages (SKYBOT, SCOS, etc.).

Metrics
-------
Two independent rejection criteria (OR logic):

1. fwhm_ratio = FWHM_IMAGE / psf_fwhm_median
   Reject if fwhm_ratio > --fwhm-ratio-max (default 1.5)
   Catches extended blobs, halos, plate scratches.

2. spread_snr = (SPREAD_MODEL - psf_spread_median) / SPREADERR_MODEL
   Reject if spread_snr > --spread-snr-max (default 5.0)
   Catches profile deviations relative to the local PSFEx model.

Note: CLASS_STAR is intentionally NOT used. It is unreliable on photographic
(POSS-I) plates — PSF reference stars themselves score ~0.015 on a 0–1 scale
because the neural net is trained on CCD data.

PSF Reference Sample (per tile)
--------------------------------
Drawn from catalogs/sextractor_pass2.csv (full detection catalog, not just
surviving candidates), so the PSF sample is large (~1500 stars/tile) and
independent of what survived the gates:
  - FLAGS = 0
  - ELONGATION < --elongation-max (default 1.3)
  - Positionally matched to a Gaia DR3 star within --gaia-match-arcsec (3")
  - Gaia Gmag in (--gaia-mag-min, --gaia-mag-max) = (12, 18)

If a tile yields fewer than --min-psf-stars (default 5) reference stars, all
candidates from that tile pass through unchanged (flagged psf_insufficient).

Calibration
-----------
On a 181-tile, 684-candidate test run:
  - fwhm_ratio > 1.5:          17.8% rejection
  - spread_snr > 5.0:          49.1% rejection
  - Either (combined):         50.3% rejection
These are conservative starting points; tune per dataset.

Usage
-----
python scripts/stage_morph_post.py \\
    --run-dir ./work/runs/run-S1-... \\
    --input-glob 'stage_S0.csv' \\
    --stage S0M \\
    --tiles-root ./data/tiles

Outputs (under <run-dir>/stages/ unless --out-dir given)
---------------------------------------------------------
1) stage_<STAGE>_MORPH.csv
   Kept remainder AFTER morphology rejection.
   Columns: src_id, ra, dec

2) stage_<STAGE>_MORPH_flags.csv
   Full audit table for ALL input rows.
   Columns: src_id, ra, dec, tile_id, object_id,
            fwhm_image, fwhm_ratio,
            spread_model, spreaderr_model, spread_snr,
            psf_fwhm_median, psf_spread_median, psf_star_count,
            reject_flag, reject_reason, source_chunk

3) stage_<STAGE>_MORPH_ledger.json
   Counts, per-tile PSF stats, and all parameters used.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Spatial helpers — fast Gaia proximity lookup via sorted dec index
# ---------------------------------------------------------------------------

def _haversine_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Great-circle separation in arcsec (haversine; numerically stable)."""
    ra1r = math.radians(ra1 % 360.0)
    ra2r = math.radians(ra2 % 360.0)
    dec1r = math.radians(dec1)
    dec2r = math.radians(dec2)
    dra = ra2r - ra1r
    ddec = dec2r - dec1r
    a = math.sin(ddec / 2) ** 2 + math.cos(dec1r) * math.cos(dec2r) * math.sin(dra / 2) ** 2
    return 2.0 * math.degrees(math.asin(math.sqrt(min(1.0, max(0.0, a))))) * 3600.0


def _build_sorted_gaia(
    gaia_rows: List[Tuple[float, float]],
) -> Tuple[List[Tuple[float, float]], List[float]]:
    """Sort Gaia (ra, dec) pairs by declination; return (sorted_rows, dec_list)."""
    s = sorted(gaia_rows, key=lambda x: x[1])
    return s, [r[1] for r in s]


def _gaia_has_match(
    ra: float,
    dec: float,
    gaia_sorted: List[Tuple[float, float]],
    gaia_decs: List[float],
    radius_arcsec: float,
) -> bool:
    """Return True if (ra, dec) has any Gaia star within radius_arcsec.

    Uses binary search on dec to restrict candidates, then a quick RA pre-filter
    before the full haversine.  O(log N + K) where K << N.
    """
    radius_deg = radius_arcsec / 3600.0
    lo = bisect.bisect_left(gaia_decs, dec - radius_deg)
    hi = bisect.bisect_right(gaia_decs, dec + radius_deg)
    if lo >= hi:
        return False
    cos_dec = math.cos(math.radians(dec))
    ra_tol = radius_deg / max(cos_dec, 1e-6)
    for i in range(lo, hi):
        gra, gdec = gaia_sorted[i]
        dra = abs(gra - ra)
        if dra > 180.0:
            dra = 360.0 - dra
        if dra > ra_tol:
            continue
        if _haversine_arcsec(ra, dec, gra, gdec) <= radius_arcsec:
            return True
    return False


# ---------------------------------------------------------------------------
# CSV helpers
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


# ---------------------------------------------------------------------------
# Tile PSF model
# ---------------------------------------------------------------------------

@dataclass
class TilePSF:
    tile_id: str
    psf_fwhm_median: float      # NaN when status != "ok"
    psf_spread_median: float    # NaN when status != "ok"
    psf_star_count: int
    status: str                 # "ok" | "psf_insufficient" | "catalog_missing"


def _load_tile_psf(
    tile_id: str,
    tiles_root: Path,
    *,
    gaia_mag_min: float,
    gaia_mag_max: float,
    gaia_match_arcsec: float,
    elongation_max: float,
    min_psf_stars: int,
) -> TilePSF:
    tile_dir = tiles_root / tile_id
    sex_csv = tile_dir / "catalogs" / "sextractor_pass2.csv"
    # Prefer plate-epoch Gaia (Bug #6 fix in vasco/cli_pipeline.py:1306) so PSF
    # reference matching at 3" doesn't silently exclude HPM stars whose J2016
    # positions have drifted from their POSS-epoch detection. Fall through to
    # the raw catalog when the propagated variant is absent or empty.
    gaia_at_plate = tile_dir / "catalogs" / "gaia_neighbourhood_at_plate.csv"
    gaia_raw = tile_dir / "catalogs" / "gaia_neighbourhood.csv"
    if gaia_at_plate.exists() and gaia_at_plate.stat().st_size > 0:
        gaia_csv = gaia_at_plate
    else:
        gaia_csv = gaia_raw

    _nan = float("nan")

    if not sex_csv.exists() or not gaia_csv.exists():
        return TilePSF(tile_id, _nan, _nan, 0, "catalog_missing")

    # Load Gaia reference stars
    gaia_rows: List[Tuple[float, float]] = []
    try:
        with gaia_csv.open(newline="", encoding="utf-8", errors="ignore") as f:
            for row in csv.DictReader(f):
                try:
                    gmag = float(row["Gmag"])
                    if gaia_mag_min < gmag < gaia_mag_max:
                        gaia_rows.append((float(row["ra"]), float(row["dec"])))
                except Exception:
                    pass
    except Exception:
        return TilePSF(tile_id, _nan, _nan, 0, "catalog_missing")

    if not gaia_rows:
        return TilePSF(tile_id, _nan, _nan, 0, "psf_insufficient")

    gaia_sorted, gaia_decs = _build_sorted_gaia(gaia_rows)

    # Scan SExtractor catalog for PSF reference stars
    psf_fwhms: List[float] = []
    psf_spreads: List[float] = []
    try:
        with sex_csv.open(newline="", encoding="utf-8", errors="ignore") as f:
            for row in csv.DictReader(f):
                if row.get("FLAGS", "").strip() != "0":
                    continue
                try:
                    elong = float(row["ELONGATION"])
                    if elong >= elongation_max:
                        continue
                    fwhm = float(row["FWHM_IMAGE"])
                    spread = float(row["SPREAD_MODEL"])
                    ra = float(row["ALPHAWIN_J2000"])
                    dec = float(row["DELTAWIN_J2000"])
                except Exception:
                    continue
                if _gaia_has_match(ra, dec, gaia_sorted, gaia_decs, gaia_match_arcsec):
                    psf_fwhms.append(fwhm)
                    psf_spreads.append(spread)
    except Exception:
        return TilePSF(tile_id, _nan, _nan, 0, "catalog_missing")

    n = len(psf_fwhms)
    if n < min_psf_stars:
        return TilePSF(tile_id, _nan, _nan, n, "psf_insufficient")

    return TilePSF(
        tile_id=tile_id,
        psf_fwhm_median=statistics.median(psf_fwhms),
        psf_spread_median=statistics.median(psf_spreads),
        psf_star_count=n,
        status="ok",
    )


# ---------------------------------------------------------------------------
# Candidate SExtractor data loader
# ---------------------------------------------------------------------------

def _load_sex_candidates(
    sex_csv: Path,
    object_ids: List[str],
) -> Dict[str, dict]:
    """Return {object_id: {fwhm, spread, spreaderr}} for the requested NUMBERs."""
    wanted = set(object_ids)
    out: Dict[str, dict] = {}
    try:
        with sex_csv.open(newline="", encoding="utf-8", errors="ignore") as f:
            for row in csv.DictReader(f):
                num = row.get("NUMBER", "").strip()
                if num not in wanted:
                    continue
                try:
                    out[num] = {
                        "fwhm": float(row["FWHM_IMAGE"]),
                        "spread": float(row["SPREAD_MODEL"]),
                        "spreaderr": float(row["SPREADERR_MODEL"]),
                    }
                except Exception:
                    pass
    except Exception:
        pass
    return out


# ---------------------------------------------------------------------------
# Per-candidate evaluation
# ---------------------------------------------------------------------------

_FLAGS_FIELDS = [
    "src_id", "ra", "dec", "tile_id", "object_id",
    "fwhm_image", "fwhm_ratio",
    "spread_model", "spreaderr_model", "spread_snr",
    "psf_fwhm_median", "psf_spread_median", "psf_star_count",
    "reject_flag", "reject_reason", "source_chunk",
]

_NAN = float("nan")


def _fmt(v: Optional[float]) -> str:
    if v is None or (v != v):  # None or NaN
        return ""
    return f"{v:.6g}"


def _evaluate(
    src_id: str,
    ra: str,
    dec: str,
    tile_id: str,
    object_id: str,
    sex_data: Optional[dict],
    psf: TilePSF,
    fwhm_ratio_max: float,
    spread_snr_max: float,
    source_chunk: str,
) -> dict:
    base = {
        "src_id": src_id, "ra": ra, "dec": dec,
        "tile_id": tile_id, "object_id": object_id,
        "psf_fwhm_median": _fmt(psf.psf_fwhm_median),
        "psf_spread_median": _fmt(psf.psf_spread_median),
        "psf_star_count": str(psf.psf_star_count),
        "source_chunk": source_chunk,
    }

    if psf.status != "ok":
        return {**base,
                "fwhm_image": "", "fwhm_ratio": "",
                "spread_model": "", "spreaderr_model": "", "spread_snr": "",
                "reject_flag": 0, "reject_reason": psf.status}

    if sex_data is None:
        return {**base,
                "fwhm_image": "", "fwhm_ratio": "",
                "spread_model": "", "spreaderr_model": "", "spread_snr": "",
                "reject_flag": 0, "reject_reason": "candidate_not_found"}

    fwhm = sex_data["fwhm"]
    spread = sex_data["spread"]
    spreaderr = sex_data["spreaderr"]

    fwhm_ratio = fwhm / psf.psf_fwhm_median if psf.psf_fwhm_median > 0 else _NAN
    spread_snr = (spread - psf.psf_spread_median) / spreaderr if spreaderr > 0 else _NAN

    reasons = []
    reject = 0

    if fwhm_ratio == fwhm_ratio and fwhm_ratio > fwhm_ratio_max:  # NaN-safe
        reasons.append(f"fwhm_ratio={fwhm_ratio:.3f}>{fwhm_ratio_max}")
        reject = 1
    if spread_snr == spread_snr and spread_snr > spread_snr_max:  # NaN-safe
        reasons.append(f"spread_snr={spread_snr:.1f}>{spread_snr_max}")
        reject = 1

    return {**base,
            "fwhm_image": _fmt(fwhm),
            "fwhm_ratio": _fmt(fwhm_ratio),
            "spread_model": _fmt(spread),
            "spreaderr_model": _fmt(spreaderr),
            "spread_snr": _fmt(spread_snr),
            "reject_flag": reject,
            "reject_reason": ";".join(reasons)}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="[EXPERIMENTAL] Post-pipeline morphology-based filtering stage."
    )
    ap.add_argument("--run-dir", required=True, help="Run folder, e.g. ./work/runs/run-S1-...")
    ap.add_argument(
        "--input-glob", default="stage_S0.csv",
        help="Glob (relative to run-dir) for input stage CSV. Default: stage_S0.csv",
    )
    ap.add_argument("--stage", default="S0M", help="Stage label for output filenames. Default: S0M")
    ap.add_argument("--out-dir", default=None, help="Output directory. Default: <run-dir>/stages")
    ap.add_argument("--tiles-root", default="./data/tiles", help="Root of tile directories. Default: ./data/tiles")
    ap.add_argument("--verbose", action="store_true", help="Print per-tile progress.")

    # PSF reference star selection
    ap.add_argument("--gaia-mag-min", type=float, default=12.0, help="Min Gaia Gmag for PSF stars. Default: 12.0")
    ap.add_argument("--gaia-mag-max", type=float, default=18.0, help="Max Gaia Gmag for PSF stars. Default: 18.0")
    ap.add_argument("--gaia-match-arcsec", type=float, default=3.0, help="Gaia match radius for PSF star ID (arcsec). Default: 3.0")
    ap.add_argument("--elongation-max", type=float, default=1.3, help="Max ELONGATION for PSF reference stars. Default: 1.3")
    ap.add_argument("--min-psf-stars", type=int, default=5, help="Min PSF stars per tile to apply filtering. Default: 5")

    # Rejection thresholds
    ap.add_argument("--fwhm-ratio-max", type=float, default=1.5,
                    help="Reject if FWHM_IMAGE / psf_fwhm_median > this. Default: 1.5")
    ap.add_argument("--spread-snr-max", type=float, default=5.0,
                    help="Reject if (SPREAD_MODEL - psf_spread_median) / SPREADERR_MODEL > this. Default: 5.0")
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

    stage = args.stage
    out_kept   = out_dir / f"stage_{stage}_MORPH.csv"
    out_flags  = out_dir / f"stage_{stage}_MORPH_flags.csv"
    out_ledger = out_dir / f"stage_{stage}_MORPH_ledger.json"

    # Pass 1: read all input rows; group object_ids by tile
    all_rows: List[dict] = []
    tiles_needed: Dict[str, List[str]] = {}

    for ch in chunks:
        with ch.open(newline="", encoding="utf-8", errors="ignore") as f:
            for row in csv.DictReader(f):
                sid = (row.get(src_col) or "").strip()
                ra  = (row.get(ra_col)  or "").strip()
                dec = (row.get(dec_col) or "").strip()
                if not sid or not ra or not dec:
                    continue
                tile_id, object_id = _parse_src_id(sid)
                all_rows.append({
                    "src_id": sid, "ra": ra, "dec": dec,
                    "tile_id": tile_id, "object_id": object_id,
                    "source_chunk": ch.name,
                })
                tiles_needed.setdefault(tile_id, []).append(object_id)

    total_in = len(all_rows)
    print(f"[MORPH] input_rows={total_in} tiles={len(tiles_needed)}")

    # Pass 2: load per-tile PSF models and candidate SExtractor data
    tile_psf: Dict[str, TilePSF] = {}
    tile_sex: Dict[str, Dict[str, dict]] = {}

    n_ok = n_insufficient = n_missing = 0
    for tile_id, obj_ids in tiles_needed.items():
        psf = _load_tile_psf(
            tile_id, tiles_root,
            gaia_mag_min=args.gaia_mag_min,
            gaia_mag_max=args.gaia_mag_max,
            gaia_match_arcsec=args.gaia_match_arcsec,
            elongation_max=args.elongation_max,
            min_psf_stars=args.min_psf_stars,
        )
        tile_psf[tile_id] = psf

        if psf.status == "ok":
            n_ok += 1
            sex_csv = tiles_root / tile_id / "catalogs" / "sextractor_pass2.csv"
            tile_sex[tile_id] = _load_sex_candidates(sex_csv, obj_ids)
        else:
            tile_sex[tile_id] = {}
            if psf.status == "psf_insufficient":
                n_insufficient += 1
            else:
                n_missing += 1

        if args.verbose:
            print(f"[MORPH]   {tile_id}: {psf.status} psf_stars={psf.psf_star_count}"
                  + (f" fwhm_med={psf.psf_fwhm_median:.2f}" if psf.status == "ok" else ""))

    print(f"[MORPH] tiles: ok={n_ok} psf_insufficient={n_insufficient} catalog_missing={n_missing}")

    # Pass 3: evaluate all candidates
    evaluated: List[dict] = []
    for r in all_rows:
        tid = r["tile_id"]
        psf = tile_psf[tid]
        sex_data = tile_sex[tid].get(r["object_id"]) if psf.status == "ok" else None
        evaluated.append(_evaluate(
            src_id=r["src_id"], ra=r["ra"], dec=r["dec"],
            tile_id=tid, object_id=r["object_id"],
            sex_data=sex_data, psf=psf,
            fwhm_ratio_max=args.fwhm_ratio_max,
            spread_snr_max=args.spread_snr_max,
            source_chunk=r["source_chunk"],
        ))

    total_rejected = sum(1 for e in evaluated if e["reject_flag"] == 1)
    total_kept = total_in - total_rejected

    # Write outputs
    with out_kept.open("w", newline="", encoding="utf-8") as fk, \
         out_flags.open("w", newline="", encoding="utf-8") as ff:
        kept_w  = csv.DictWriter(fk, fieldnames=["src_id", "ra", "dec"])
        flags_w = csv.DictWriter(ff, fieldnames=_FLAGS_FIELDS)
        kept_w.writeheader()
        flags_w.writeheader()
        for e in evaluated:
            flags_w.writerow(e)
            if e["reject_flag"] == 0:
                kept_w.writerow({"src_id": e["src_id"], "ra": e["ra"], "dec": e["dec"]})

    # Per-tile summary for ledger
    tile_summaries = []
    for tile_id in sorted(tiles_needed.keys()):
        psf = tile_psf[tile_id]
        tm = [e for e in evaluated if e["tile_id"] == tile_id]
        t_rej = sum(1 for e in tm if e["reject_flag"] == 1)
        tile_summaries.append({
            "tile_id": tile_id,
            "psf_status": psf.status,
            "psf_star_count": psf.psf_star_count,
            "psf_fwhm_median": None if psf.psf_fwhm_median != psf.psf_fwhm_median else psf.psf_fwhm_median,
            "psf_spread_median": None if psf.psf_spread_median != psf.psf_spread_median else psf.psf_spread_median,
            "input_rows": len(tm),
            "rejected_rows": t_rej,
            "kept_rows": len(tm) - t_rej,
        })

    ledger = {
        "experimental": True,
        "run_dir": str(run_dir),
        "input_glob": args.input_glob,
        "stage": stage,
        "tiles_root": str(tiles_root),
        "parameters": {
            "gaia_mag_min": args.gaia_mag_min,
            "gaia_mag_max": args.gaia_mag_max,
            "gaia_match_arcsec": args.gaia_match_arcsec,
            "elongation_max": args.elongation_max,
            "min_psf_stars": args.min_psf_stars,
            "fwhm_ratio_max": args.fwhm_ratio_max,
            "spread_snr_max": args.spread_snr_max,
        },
        "input_chunks": [p.name for p in chunks],
        "totals": {
            "input_rows": total_in,
            "rejected_rows": total_rejected,
            "kept_rows": total_kept,
            "rejection_pct": round(100.0 * total_rejected / total_in, 2) if total_in > 0 else 0.0,
        },
        "tile_psf_summary": {
            "ok": n_ok,
            "psf_insufficient": n_insufficient,
            "catalog_missing": n_missing,
        },
        "per_tile": tile_summaries,
        "outputs": {
            "kept_csv": str(out_kept),
            "flags_csv": str(out_flags),
            "ledger_json": str(out_ledger),
        },
        "columns_detected": {"src_id_col": src_col, "ra_col": ra_col, "dec_col": dec_col},
    }
    out_ledger.write_text(json.dumps(ledger, indent=2), encoding="utf-8")

    pct = 100.0 * total_rejected / total_in if total_in > 0 else 0.0
    print(f"[MORPH] [EXPERIMENTAL] input={total_in} rejected={total_rejected} kept={total_kept} ({pct:.1f}% reduction)")
    print(f"[MORPH] wrote: {out_kept}")
    print(f"[MORPH] wrote: {out_flags}")
    print(f"[MORPH] wrote: {out_ledger}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
