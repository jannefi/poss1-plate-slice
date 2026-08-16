#!/usr/bin/env python3
from __future__ import annotations

# -*- coding: utf-8 -*-
"""
Compute per-tile plate-edge classification for the current naive-
tessellation archive.

Adapted from the original VASCO60 tool (jannefi/vasco commit b1cefc05,
scripts/check_tile_plate_edge.py) -- the core geometry
(radec_to_plate_pixels_gnomonic / enforce_east_left_orientation /
min_edge_distance_px / classify_px / classify_arcsec) is kept as-is;
only tile/plate *resolution* is simplified since this repo's naive-
tessellation plan already gives a clean 1:1 tile_id -> plate_id mapping
(verified 2026-07-27: 31,113 tiles, zero mapped to more than one plate)
-- no REGION-sidecar lookup needed.

For each (tile_id, plate_id) pair in the tile plan:
  1. Read the tile's own WCS from data/tiles_archive/<tile_id>/raw/*.header.json
     (skipped/flagged if the tile isn't archived yet).
  2. Read the plate's own FITS header (PLATERA/PLATEDEC/PLTSCALE/
     XPIXELSZ/YPIXELSZ/NAXIS1/NAXIS2) from the plate header registry.
  3. Project the tile's 8 boundary sample points (4 corners + 4 edge
     midpoints) into the plate's own pixel frame via a gnomonic
     projection centered on the plate's PLATERA/PLATEDEC.
  4. Classify by minimum distance to the plate's [1..NAXIS1]x[1..NAXIS2]
     boundary, in both pixels and arcsec, against configurable
     thresholds (defaults match the original VASCO60 run).

Output: data/metadata/tile_plate_edge_report.csv, same column schema as
the original tool (the consumer that used it, scripts/stage_edge_post.py,
was retired 2026-08-16 in favour of scripts/stage_edge_post_v2.py), so it
needs no changes:
  tile_id, plate_id, plate_nx, plate_ny, tile_center_ra_deg,
  tile_center_dec_deg, plate_center_ra_deg, plate_center_dec_deg,
  sep_arcmin, min_edge_dist_px, min_edge_dist_arcsec, class_px,
  class_arcsec, notes

This is a read-only, non-wired-in demonstration tool -- it does not
touch the primary pipeline, tessellation, or any orchestrator.

Usage:
  python scripts/compute_tile_plate_edge_report.py \\
      --tile-plan-csv plans/tiles_mnras_plates_naive.csv \\
      --tiles-root data/tiles_archive \\
      --headers-dir <plate_headers> \\
      --out-csv data/metadata/tile_plate_edge_report.csv
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from vasco.paths import get as _p


import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from astropy.io.fits import Header
from astropy.wcs import WCS

RAD2AS = 206264.80624709636  # arcsec per radian


# --------------------------- robust JSON loader ---------------------------
def robust_json_load(path: Path):
    try:
        b = path.read_bytes()
    except Exception as e:
        return None, f"read_bytes failed: {e}"
    for enc in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            s = b.decode(enc)
            try:
                return json.loads(s), None
            except Exception:
                try:
                    return json.loads(s, strict=False), None
                except Exception:
                    pass
        except Exception:
            pass
    return None, "json decode failed after fallbacks"


def pick_header_dict(any_json: dict):
    if isinstance(any_json, dict):
        if "header" in any_json and isinstance(any_json["header"], dict):
            return any_json["header"]
        if "selected" in any_json and isinstance(any_json["selected"], dict):
            return any_json["selected"]
    return any_json


# -------------------------- tile WCS (from tile raw) ----------------------
def load_tile_wcs_from_json(tile_dir: Path):
    raw = tile_dir / "raw"
    cands = sorted(raw.glob("*.header.json"))
    if not cands:
        return None, None, None, "", "no .header.json under raw/"
    tj = cands[0]
    data, err = robust_json_load(tj)
    if data is None:
        return None, None, None, tj.name, f"cannot parse JSON: {err}"
    hdr = pick_header_dict(data)
    try:
        H = Header()
        H["NAXIS"] = 2
        H["NAXIS1"] = int(hdr.get("NAXIS1", 0) or 0)
        H["NAXIS2"] = int(hdr.get("NAXIS2", 0) or 0)
        for k in ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2"):
            v = hdr.get(k, None)
            if v is not None:
                H[k] = float(v)
        cd_keys = ("CD1_1", "CD1_2", "CD2_1", "CD2_2")
        if all(hdr.get(k, None) is not None for k in cd_keys):
            for k in cd_keys:
                H[k] = float(hdr[k])
        else:
            for k in ("CDELT1", "CDELT2", "CROTA2"):
                v = hdr.get(k, None)
                if v is not None:
                    H[k] = float(v)
        H["CTYPE1"] = hdr.get("CTYPE1", "RA---TAN")
        H["CTYPE2"] = hdr.get("CTYPE2", "DEC--TAN")
        w = WCS(H)
        return w, H["NAXIS1"], H["NAXIS2"], tj.name, ""
    except Exception as e:
        return None, int(hdr.get("NAXIS1", 0) or 0), int(hdr.get("NAXIS2", 0) or 0), tj.name, f"WCS build failed: {e}"


# --------------------------- plate core (no WCS) --------------------------
def load_plate_core_from_json(path: Path):
    """Plate geometry for the gnomonic projection.

    The tangent point must be the sky position of pixel (cx, cy) -- the centre
    of the *scan*. PLATERA/PLATEDEC is the centre of the *plate*, and on some
    scans those are not the same position.

    Measured over all 932 DSS1-red headers (2026-08-16): the two agree to a
    median 0.062 deg, but **7 plates disagree by ~4.4 deg** -- XE761, XE758,
    XE733, XE574, XE284, XE543, XE541 -- plus XE304 (1.33 deg) and XE293
    (1.16 deg). PLATERA/PLATEDEC there agrees exactly with the DSS PLT*
    plate-solution keywords, so neither value is corrupt; the scan simply is
    not centred on the plate centre while CNPIX still reads (0, 0).

    Using PLATERA as the tangent point while pairing it with pixel (cx, cy)
    silently asserts that the scan IS centred there, which put rows on those
    plates up to 8.8 deg from their nominal centre -- geometrically impossible
    and previously misread as a plate-attribution defect. It is not one: the
    tessellation, the tile->plate map and the detections are all correct
    (those plates' veto survival sits in the 4th-33rd percentile, so their
    astrometry demonstrably works), and only this projection was wrong.

    So: take the centre from the scan's own WCS, and fall back to
    PLATERA/PLATEDEC only if no usable WCS can be built. The offset between
    the two is returned so callers can report it.
    """
    data, err = robust_json_load(path)
    if data is None:
        return None, "cannot parse plate JSON: " + str(err)
    hdr = pick_header_dict(data)
    try:
        nax1 = int(hdr.get("NAXIS1", 0) or 0)
        nax2 = int(hdr.get("NAXIS2", 0) or 0)
        pl_ra = float(hdr["PLATERA"])
        pl_de = float(hdr["PLATEDEC"])
        pltscale = float(hdr["PLTSCALE"])  # arcsec/mm
        xp_um = float(hdr["XPIXELSZ"])     # microns
        yp_um = float(hdr["YPIXELSZ"])
    except Exception:
        return None, "missing required plate keys (NAXIS*, PLATERA/DEC, PLTSCALE, X/YPIXELSZ)"
    cx = nax1 / 2.0
    cy = nax2 / 2.0
    as_per_px_x = pltscale * (xp_um / 1000.0)
    as_per_px_y = pltscale * (yp_um / 1000.0)

    centre_ra, centre_dec = pl_ra, pl_de
    centre_source = "PLATERA/PLATEDEC (no usable WCS)"
    offset_deg = float("nan")
    try:
        w = WCS(Header(hdr.items()), relax=True)
        sky = w.all_pix2world([[cx, cy]], 0)[0]
        wra, wdec = float(sky[0]), float(sky[1])
        if math.isfinite(wra) and math.isfinite(wdec):
            centre_ra, centre_dec = wra, wdec
            centre_source = "scan WCS at (cx, cy)"
            d1 = math.radians(pl_de)
            d2 = math.radians(wdec)
            dr = math.radians(wra - pl_ra)
            s = (math.sin((d2 - d1) / 2.0) ** 2
                 + math.cos(d1) * math.cos(d2) * math.sin(dr / 2.0) ** 2)
            offset_deg = math.degrees(2.0 * math.asin(math.sqrt(min(max(s, 0.0), 1.0))))
    except Exception:
        pass

    return {
        "nax1": nax1, "nax2": nax2,
        "center_ra": centre_ra, "center_dec": centre_dec,
        "plate_ra": pl_ra, "plate_dec": pl_de,
        "center_source": centre_source,
        "center_offset_deg": offset_deg,
        "cx": cx, "cy": cy,
        "as_per_px_x": as_per_px_x, "as_per_px_y": as_per_px_y,
    }, ""


# ------------------------------- geometry ---------------------------------
def radec_to_plate_pixels_gnomonic(ra_deg: np.ndarray, dec_deg: np.ndarray, plate: dict) -> np.ndarray:
    ra0 = math.radians(plate["center_ra"])
    de0 = math.radians(plate["center_dec"])
    as_px_x = plate["as_per_px_x"]
    as_px_y = plate["as_per_px_y"]
    cx, cy = plate["cx"], plate["cy"]
    ra = np.radians(ra_deg)
    de = np.radians(dec_deg)
    dra = (ra - ra0 + np.pi) % (2 * np.pi) - np.pi
    sin_de = np.sin(de)
    cos_de = np.cos(de)
    sin_de0 = math.sin(de0)
    cos_de0 = math.cos(de0)
    cos_dra = np.cos(dra)
    sin_dra = np.sin(dra)
    denom = sin_de0 * sin_de + cos_de0 * cos_de * cos_dra
    denom = np.where(np.abs(denom) < 1e-12, np.nan, denom)
    xi = (cos_de * sin_dra) / denom
    eta = (cos_de0 * sin_de - sin_de0 * cos_de * cos_dra) / denom
    xi_as, eta_as = xi * RAD2AS, eta * RAD2AS
    x = cx + (xi_as / as_px_x)
    y = cy + (eta_as / as_px_y)
    return np.stack([x, y], axis=1)


def enforce_east_left_orientation(plate: dict):
    eps_deg = 5.0 / 3600.0
    ra0 = plate["center_ra"]
    dec0 = plate["center_dec"]
    xy_east = radec_to_plate_pixels_gnomonic(np.array([ra0 + eps_deg]), np.array([dec0]), plate)[0]
    xy_north = radec_to_plate_pixels_gnomonic(np.array([ra0]), np.array([dec0 + eps_deg]), plate)[0]
    cx, cy = plate["cx"], plate["cy"]
    flip_x = xy_east[0] > cx
    flip_y = xy_north[1] < cy
    return {"flip_x": bool(flip_x), "flip_y": bool(flip_y)}


def min_edge_distance_px(points_xy: np.ndarray, plate_nx: int, plate_ny: int) -> float:
    if points_xy.size == 0:
        return float("nan")
    dvals = []
    for x, y in points_xy:
        if np.isnan(x) or np.isnan(y):
            continue
        if x < 1 or x > plate_nx or y < 1 or y > plate_ny:
            d = -min(abs(x - 1), abs(plate_nx - x), abs(y - 1), abs(plate_ny - y))
        else:
            d = min(x - 1, plate_nx - x, y - 1, plate_ny - y)
        dvals.append(d)
    if not dvals:
        return float("nan")
    return float(min(dvals))


def classify_px(mind_px: float, nx: int, ny: int, th_px: float, th_frac: float) -> str:
    if math.isnan(mind_px):
        return "off_plate"
    if mind_px < 0:
        return "edge_touch"
    thresh = max(float(th_px), float(th_frac) * float(min(nx, ny)))
    return "core" if mind_px >= thresh else "near_edge"


def classify_arcsec(mind_as, th_as: float) -> str:
    if mind_as is None or math.isnan(mind_as):
        return "off_plate"
    if mind_as < 0:
        return "edge_touch"
    return "core" if mind_as >= th_as else "near_edge"


def load_tile_plan(csv_path: Path):
    """Yield (tile_id, plate_id) pairs from the naive-tessellation plan."""
    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        for row in csv.DictReader(f):
            tid = (row.get("tile_id") or "").strip()
            pid = (row.get("plate_id") or "").strip()
            if tid and pid:
                yield tid, pid


def main() -> int:
    ap = argparse.ArgumentParser(description="Compute per-tile plate-edge classification report.")
    ap.add_argument("--tile-plan-csv", default="plans/tiles_mnras_plates_naive.csv",
                    help="tile_id -> plate_id mapping (naive-tessellation plan).")
    ap.add_argument("--tiles-root", default="./data/tiles_archive",
                    help="Root containing per-tile dirs with raw/*.header.json.")
    ap.add_argument("--headers-dir", default=str(_p("plate_headers")),
                    help="Plate header registry root (dss1red_<plate_id>.fits.header.json).")
    ap.add_argument("--out-csv", default="./data/metadata/tile_plate_edge_report.csv")
    ap.add_argument("--threshold-px", type=float, default=200.0)
    ap.add_argument("--threshold-frac", type=float, default=0.02)
    ap.add_argument("--threshold-arcsec", type=float, default=60.0)
    args = ap.parse_args()

    tiles_root = Path(args.tiles_root)
    headers_dir = Path(args.headers_dir)
    out_csv = Path(args.out_csv)

    plate_info_cache: dict[str, tuple] = {}  # plate_id -> (plate_info, err)
    rows_out = []
    n_ok = n_no_tile = n_no_plate = 0

    for tile_id, plate_id in load_tile_plan(Path(args.tile_plan_csv)):
        tile_dir = tiles_root / tile_id
        if not tile_dir.exists():
            n_no_tile += 1
            continue  # not yet archived -- silently skip, not an error

        tile_wcs, tnx, tny, tile_json_name, tile_err = load_tile_wcs_from_json(tile_dir)
        if tile_wcs is None or not tnx or not tny:
            rows_out.append({
                "tile_id": tile_id, "plate_id": plate_id,
                "plate_nx": "", "plate_ny": "",
                "tile_center_ra_deg": "", "tile_center_dec_deg": "",
                "plate_center_ra_deg": "", "plate_center_dec_deg": "",
                "sep_arcmin": "", "min_edge_dist_px": "", "min_edge_dist_arcsec": "",
                "class_px": "off_plate", "class_arcsec": "off_plate",
                "notes": f"tile WCS unavailable: {tile_err}",
            })
            continue

        if plate_id not in plate_info_cache:
            plate_json = headers_dir / f"dss1red_{plate_id}.fits.header.json"
            if not plate_json.exists():
                plate_info_cache[plate_id] = (None, f"plate header not found: {plate_json}")
            else:
                plate_info_cache[plate_id] = load_plate_core_from_json(plate_json)
        plate_info, p_err = plate_info_cache[plate_id]

        if plate_info is None:
            n_no_plate += 1
            rows_out.append({
                "tile_id": tile_id, "plate_id": plate_id,
                "plate_nx": "", "plate_ny": "",
                "tile_center_ra_deg": "", "tile_center_dec_deg": "",
                "plate_center_ra_deg": "", "plate_center_dec_deg": "",
                "sep_arcmin": "", "min_edge_dist_px": "", "min_edge_dist_arcsec": "",
                "class_px": "off_plate", "class_arcsec": "off_plate",
                "notes": p_err,
            })
            continue

        orient = enforce_east_left_orientation(plate_info)

        cx, cy = tnx / 2.0, tny / 2.0
        tile_ra_dec = tile_wcs.all_pix2world(np.array([[cx, cy]]), 1)[0]
        tile_cent = (float(tile_ra_dec[0]), float(tile_ra_dec[1]))

        cosd = math.cos(math.radians(plate_info["center_dec"]))
        dra = ((tile_cent[0] - plate_info["center_ra"] + 180.0) % 360.0) - 180.0
        sep_arcmin = math.hypot(dra * cosd, tile_cent[1] - plate_info["center_dec"]) * 60.0

        corners = np.array([[1, 1], [tnx, 1], [tnx, tny], [1, tny], [1, 1]], dtype=float)
        mids = np.array([[tnx / 2, 1], [tnx, tny / 2], [tnx / 2, tny], [1, tny / 2]], dtype=float)
        samples = np.vstack([corners, mids])
        world = tile_wcs.all_pix2world(samples, 1)
        plate_xy = radec_to_plate_pixels_gnomonic(world[:, 0], world[:, 1], plate_info)

        if orient["flip_x"]:
            plate_xy[:, 0] = 2.0 * plate_info["cx"] - plate_xy[:, 0]
        if orient["flip_y"]:
            plate_xy[:, 1] = 2.0 * plate_info["cy"] - plate_xy[:, 1]

        px_margin = min_edge_distance_px(plate_xy, plate_info["nax1"], plate_info["nax2"])
        as_per_px = 0.5 * (plate_info["as_per_px_x"] + plate_info["as_per_px_y"])
        as_margin = px_margin * as_per_px if not math.isnan(px_margin) else float("nan")

        c_px = classify_px(px_margin, plate_info["nax1"], plate_info["nax2"], args.threshold_px, args.threshold_frac)
        c_as = classify_arcsec(as_margin, args.threshold_arcsec)

        rows_out.append({
            "tile_id": tile_id, "plate_id": plate_id,
            "plate_nx": plate_info["nax1"], "plate_ny": plate_info["nax2"],
            "tile_center_ra_deg": f"{tile_cent[0]:.6f}", "tile_center_dec_deg": f"{tile_cent[1]:.6f}",
            "plate_center_ra_deg": f"{plate_info['center_ra']:.6f}",
            "plate_center_dec_deg": f"{plate_info['center_dec']:.6f}",
            "sep_arcmin": f"{sep_arcmin:.3f}",
            "min_edge_dist_px": f"{px_margin:.2f}" if not math.isnan(px_margin) else "",
            "min_edge_dist_arcsec": f"{as_margin:.1f}" if not math.isnan(as_margin) else "",
            "class_px": c_px, "class_arcsec": c_as,
            "notes": f"orient flip_x={orient['flip_x']} flip_y={orient['flip_y']}",
        })
        n_ok += 1

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "tile_id", "plate_id", "plate_nx", "plate_ny",
        "tile_center_ra_deg", "tile_center_dec_deg",
        "plate_center_ra_deg", "plate_center_dec_deg",
        "sep_arcmin", "min_edge_dist_px", "min_edge_dist_arcsec",
        "class_px", "class_arcsec", "notes",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows_out)

    print(f"[OK] wrote {out_csv} rows={len(rows_out)} "
          f"(classified={n_ok}, skipped_not_archived={n_no_tile}, no_plate_header={n_no_plate})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
