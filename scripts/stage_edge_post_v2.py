#!/usr/bin/env python3
"""
EXPERIMENTAL — Post-pipeline plate-core (per-row edge) stage, v2.

Status
------
Not an official veto stage and not wired into the pipeline. Same posture
as scripts/stage_gsc_post.py, scripts/stage_maps_post.py and the v1
scripts/stage_edge_post.py: standalone, read-only, comparable output
contract.

Why v2 exists
-------------
Two problems with v1 (scripts/stage_edge_post.py):

1. **It is tile-granular.** v1 joins each row's tile_id against
   data/metadata/tile_plate_edge_report.csv, which classifies a *tile*
   by the worst of its 8 boundary sample points. A ~1 deg tile whose
   corner pokes over the plate boundary condemns every candidate in it,
   including candidates a full degree inside the plate. Measured on the
   released catalogue, the same 60" threshold applied per-row touches
   0.66% of rows against the ~53% the tile rule removed on the tile set
   it was built for.

2. **On the current catalogue it is a silent no-op.** The edge report
   was built for plans/tiles_mnras_plates_naive.csv (31,064 tiles on a
   regular 1 deg dec grid). The released run uses a different
   tessellation (25,643 distinct tiles, irregular centres). The two tile
   sets intersect in *2* tiles, so 99.99% of rows hit v1's
   "edge_report_missing -> keep" default and the stage drops nothing
   while exiting ok.

What v2 does instead
--------------------
Per-row, no tile granularity, no precomputed report:

    sep = haversine(row_ra, row_dec, PLATERA, PLATEDEC)   # degrees
    in_core = sep <= core_radius_deg

The MAPS/APS documentation restricts reliable POSS-I astrometry to the
central **5.4 deg diameter** of a plate, i.e. a **2.7 deg radius** about
the plate centre, warning of geometric distortion and vignetting outside
it. That 2.7 deg is this stage's default.

Note on MAPS_CORE_RADIUS_DEG = 2.2
----------------------------------
vasco/plan/tessellate_plates.py uses 2.2 deg, defined there as

    2.2 = (5.4 / 2) - 0.5

The 0.5 deg is a **half-tile margin**: it exists so that a whole 1x1 deg
*tile* fits inside the core when only the tile *centre* is tested. A
candidate is a point and has no extent, so re-applying that margin
per-row would discard a 0.5 deg annulus of sky the APS documentation
calls reliable. 2.2 deg is therefore the correct constant for planning
tiles and the wrong one for cutting rows. Use 2.7.

A third value is meaningful if the goal is *footprint parity with old
VASCO60* rather than the APS physics claim. Under VASCO60's rule a tile
centre sits within 2.2 deg and a row within that tile can be up to a
half-diagonal further out, so rows reach

    2.2 + 0.5*sqrt(2) = 2.907 deg

from the plate centre. --vasco60-parity selects that radius.

Shape matters
-------------
The MAPS core is a **circle** about the plate centre. It is not the
square plate boundary (half-width ~3.308 deg, corners at ~4.678 deg), so
this stage cuts on radial separation, not on distance-to-boundary. The
square-frame Chebyshev distance is still recorded per row as a
diagnostic (cheb_own_deg), because it is the right metric for a
different question (how close to the physical plate edge a row sits) --
but it is not what the cut uses.

Default is FLAG, not CUT
------------------------
Because this radius is a very large lever -- on the released catalogue
2.7 deg removes 72.5% of rows and 2.2 deg removes 83.8%, the median row
sitting 3.22 deg from its plate centre -- the stage **records** the
geometry for every row and keeps everything unless --cut is passed. The
ledger always carries a yield curve across candidate radii so the cut
point is an explicit, published decision rather than a buried constant.

Two policies are recorded for every row (--policy selects which one
--cut acts on):

  own : separation from the plate the row was actually detected on.
        The honest quality statement -- our measurement came from that
        plate, at that radius.
  any : minimum separation over every plate in the run that covers the
        position. Keeps sky, but a row kept only because some *other*
        plate could have measured it well is still a row we measured
        badly, unless it is re-extracted from that other plate. Recorded
        for analysis; not recommended as the cut without re-extraction.

Usage
-----
python scripts/stage_edge_post_v2.py \\
    --run-dir ./work/runs/run-S1-... \\
    --input-glob 'stages/stage_S0.csv' \\
    --plate-map-csv /srv/vasco/vasco60/fullscale_veto/tile_plate_map.csv \\
    --stage S1

Outputs (under <run-dir>/stages/)
----------------------------------
1) stage_<STAGE>_EDGE2.csv
   Kept remainder. Identical to the input row set unless --cut.
   Columns: src_id, ra, dec

2) stage_<STAGE>_EDGE2_flags.csv
   Full audit table for ALL input rows. Columns: src_id, ra, dec,
   tile_id, det_plate, sep_own_deg, in_core_own, best_plate,
   sep_best_deg, in_core_best, cheb_own_deg, beyond_corner,
   plate_unresolved, source_chunk

3) stage_<STAGE>_EDGE2_ledger.json
   Counts + parameters, the yield curve, and three integrity counters:
   plate_unresolved_rows, beyond_corner_rows (rows further from their
   own plate centre than that plate's corner distance -- geometrically
   impossible, so a plate-attribution defect) and rows_missing_coords.

How to validate
---------------
Run the self-tests (no pipeline state needed):

    python scripts/stage_edge_post_v2.py --self-test

Against the released catalogue (S0, 122,820 rows), these must hold:
  * --core-radius-deg 2.7 with --policy own --cut keeps 33,649 (27.40%);
    2.2 keeps 19,818 (16.14%); --vasco60-parity (2.907) keeps 42,280
    (34.42%).
  * beyond_corner_rows == 0, and likewise 0 on the PTF terminal set.
    A nonzero count means the plate centre is being taken from the wrong
    keyword again -- see load_plate_core_from_json in
    compute_tile_plate_edge_report.py. Before that fix this read 347 on
    S0 and 159 on PTF, which looked like a plate-attribution defect and
    was not one.
  * plate_sources.centre_from_wcs == 932 (every plate), and
    centre_offset_gt_0p5_deg lists 10 plates, led by XE761/XE758/XE733/
    XE284/XE574/XE543/XE541 at ~4.4 deg.
  * without --cut, kept row count == input row count exactly.
  * sep_best_deg <= sep_own_deg for every row, by construction.

Threshold provenance
--------------------
The 2.7 deg default comes from MAPS/APS documentation, not from tuning.
Do NOT select this radius by what it does to the recall of a comparison
catalogue. A radius chosen that way makes every downstream number circular.
Freeze it on geometric/physical grounds -- or take somebody else's published
criterion, which is what --pasp2025 does -- and report recall as a check
afterwards, never as the selection rule.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

# Geometry is shared with the v1 report generator rather than duplicated --
# two copies of a projection drift apart silently.
_GEN = Path(__file__).resolve().parent / "compute_tile_plate_edge_report.py"


def _load_geometry_module():
    spec = importlib.util.spec_from_file_location("_edge_geom", _GEN)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load geometry helpers from {_GEN}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PLATE_HALF_WIDTH_DEG = 3.308  # square plate half-width in its own frame
CORE_RADIUS_DEFAULT_DEG = 2.7  # 5.4 deg APS reliable-core diameter / 2
VASCO60_PARITY_RADIUS_DEG = 2.2 + 0.5 * math.sqrt(2.0)  # 2.9071...
# Villarroel et al. 2025 (PASP) mask transients >2 deg from the plate centre in
# both the solar-reflection test and the GSO background-density count, naming
# "edge fingerprints or other plate defects" as the reason, and report 22,314
# transients surviving. Confirmed against Doherty's replication repo, whose
# pre-filtered SUPERVIKTIG_HELAVASCO_within2deg_CENTER.csv yields 22,309 on 614
# plates -- and against this project's own geometry, 22,317. Three independent
# routes within 8 rows. Use --pasp2025 to apply the same criterion here.
PASP2025_RADIUS_DEG = 2.0
YIELD_CURVE_RADII = [2.2, 2.4, 2.5, 2.7, VASCO60_PARITY_RADIUS_DEG, 3.0, 3.308]
# Distance to the ARRAY BOUNDARY, in arcmin. A far more targeted lever than the
# radial cut: measured on the released catalogue, removing everything within
# 10' drops 21% of rows at zero cost in matches to the published vanish-possi
# catalogue, because overlapping full-plate coverage retains the sky on a
# neighbour's interior.
EDGE_CURVE_ARCMIN = [1, 2, 3, 5, 10, 15, 20, 30]


def haversine_deg(ra1, dec1, ra2, dec2):
    """Vectorised angular separation in degrees."""
    r1 = np.radians(ra1)
    d1 = np.radians(dec1)
    r2 = np.radians(ra2)
    d2 = np.radians(dec2)
    s = np.sin((d2 - d1) / 2.0) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2.0) ** 2
    return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(s, 0.0, 1.0))))


def _read_header(path: Path) -> List[str]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.reader(f)
        return [c.strip().lstrip("﻿") for c in next(r, [])]


def _detect_cols(cols: List[str], plate_col_opt: Optional[str]) -> Tuple[str, str, str, Optional[str]]:
    cset = {c.strip() for c in cols}

    if "src_id" in cset:
        src = "src_id"
    elif "row_id" in cset:
        src = "row_id"
    else:
        raise RuntimeError("Input CSV missing required id column 'src_id' (or 'row_id')")

    ra = next((c for c in ["ra", "RA", "RA_corr", "ALPHAWIN_J2000", "ALPHA_J2000"] if c in cset), None)
    dec = next((c for c in ["dec", "DEC", "Dec", "Dec_corr", "DELTAWIN_J2000", "DELTA_J2000"] if c in cset), None)
    if not ra or not dec:
        raise RuntimeError("Input CSV missing RA/Dec columns")

    if plate_col_opt:
        if plate_col_opt not in cset:
            raise RuntimeError(f"--plate-col {plate_col_opt!r} not present in input CSV")
        plate = plate_col_opt
    else:
        plate = next((c for c in ["det_plate", "plate_id", "plate"] if c in cset), None)

    return src, ra, dec, plate


def load_plate_centres(headers_dir: Path, geom) -> Dict[str, dict]:
    """plate_id -> plate dict (centre, pixel scale, NAXIS) from the header registry."""
    out: Dict[str, dict] = {}
    for jf in sorted(headers_dir.glob("*.header.json")):
        stem = jf.name
        # dss1red_XE002.fits.header.json -> XE002
        pid = stem.split(".")[0]
        if "_" in pid:
            pid = pid.split("_", 1)[1]
        core, err = geom.load_plate_core_from_json(jf)
        if core is None:
            continue
        # Real WCS, for distance to the ARRAY BOUNDARY. cheb below is a linear
        # gnomonic approximation and drifts ~1' from the GSSS polynomial at the
        # rim -- enough to place an on-plate row outside the array -- so the
        # boundary distance must come from the plate's own solution.
        try:
            import json as _json
            from astropy.io.fits import Header as _Header
            from astropy.wcs import WCS as _WCS
            _h = _Header(geom.pick_header_dict(_json.load(open(jf))).items())
            core["wcs"] = _WCS(_h, relax=True)
        except Exception:
            core["wcs"] = None
        out[pid] = core
    return out


def load_tile_plate_map(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        for row in csv.DictReader(f):
            tid = (row.get("tile_id") or "").strip()
            pid = (row.get("plate_id") or row.get("det_plate") or "").strip()
            if tid and pid:
                out[tid] = pid
    return out


def load_src_plate_map(path: Path) -> Dict[str, str]:
    """src_id -> det_plate, e.g. from the released primary_plate_flags.csv."""
    out: Dict[str, str] = {}
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        for row in csv.DictReader(f):
            sid = (row.get("src_id") or "").strip()
            pid = (row.get("det_plate") or "").strip()
            if sid and pid:
                out[sid] = pid
    return out


@dataclass
class ChunkStats:
    chunk: str
    input_rows: int
    kept_rows: int
    in_core_own: int
    in_core_best: int
    plate_unresolved: int
    beyond_corner: int
    missing_coords: int
    yield_curve_own: Dict[str, int] = field(default_factory=dict)
    yield_curve_best: Dict[str, int] = field(default_factory=dict)
    yield_curve_edge: Dict[str, int] = field(default_factory=dict)


def chebyshev_deg(ra: np.ndarray, dec: np.ndarray, plate: dict, geom) -> np.ndarray:
    """Chebyshev distance from plate centre in the plate's own pixel frame, in degrees.

    Diagnostic only -- the cut uses radial separation. Recorded because it is
    the correct metric for 'how close to the physical plate edge is this row'.
    """
    xy = geom.radec_to_plate_pixels_gnomonic(ra, dec, plate)
    dx = np.abs(xy[:, 0] - plate["cx"]) * (plate["as_per_px_x"] / 3600.0)
    dy = np.abs(xy[:, 1] - plate["cy"]) * (plate["as_per_px_y"] / 3600.0)
    return np.maximum(dx, dy)


def process_chunk(
    path: Path,
    src_col: str,
    ra_col: str,
    dec_col: str,
    plate_col: Optional[str],
    tile_map: Dict[str, str],
    src_map: Dict[str, str],
    plates: Dict[str, dict],
    core_radius: float,
    policy: str,
    cut: bool,
    flags_w: csv.DictWriter,
    kept_w: csv.DictWriter,
) -> ChunkStats:
    src_ids: List[str] = []
    tile_ids: List[str] = []
    ras: List[float] = []
    decs: List[float] = []
    det_plates: List[str] = []
    missing_coords = 0

    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        for row in csv.DictReader(f):
            sid = (row.get(src_col) or "").strip()
            if not sid:
                continue
            try:
                ra = float((row.get(ra_col) or "").strip())
                dec = float((row.get(dec_col) or "").strip())
            except (TypeError, ValueError):
                missing_coords += 1
                continue
            tid = (row.get("tile_id") or "").strip()
            pid = ""
            if plate_col:
                pid = (row.get(plate_col) or "").strip()
            if not pid:
                pid = src_map.get(sid, "")
            if not pid and tid:
                pid = tile_map.get(tid, "")
            src_ids.append(sid)
            tile_ids.append(tid)
            ras.append(ra)
            decs.append(dec)
            det_plates.append(pid)

    n = len(src_ids)
    ra_a = np.asarray(ras, dtype=float)
    dec_a = np.asarray(decs, dtype=float)

    sep_own = np.full(n, np.nan)
    cheb_own = np.full(n, np.nan)
    edge_arcmin = np.full(n, np.nan)
    beyond_corner = np.zeros(n, dtype=bool)
    unresolved = np.zeros(n, dtype=bool)

    by_plate: Dict[str, List[int]] = {}
    for i, pid in enumerate(det_plates):
        if pid and pid in plates:
            by_plate.setdefault(pid, []).append(i)
        else:
            unresolved[i] = True

    for pid, idx in by_plate.items():
        pl = plates[pid]
        ii = np.asarray(idx, dtype=int)
        sep_own[ii] = haversine_deg(ra_a[ii], dec_a[ii], pl["center_ra"], pl["center_dec"])
        cheb_own[ii] = chebyshev_deg(ra_a[ii], dec_a[ii], pl, geom=_GEOM)
        if pl.get("wcs") is not None:
            px, py = pl["wcs"].world_to_pixel_values(ra_a[ii], dec_a[ii])
            d_px = np.minimum(np.minimum(px, pl["nax1"] - px),
                              np.minimum(py, pl["nax2"] - py))
            edge_arcmin[ii] = d_px * pl["as_per_px_x"] / 60.0
        # corner distance of this plate, in its own frame
        half_x = pl["cx"] * (pl["as_per_px_x"] / 3600.0)
        half_y = pl["cy"] * (pl["as_per_px_y"] / 3600.0)
        corner = math.hypot(half_x, half_y)
        beyond_corner[ii] = sep_own[ii] > corner

    # --- best (minimum-separation) plate over the run's plate set ---
    plate_ids = sorted({p for p in det_plates if p and p in plates})
    sep_best = np.full(n, np.nan)
    best_plate = [""] * n
    if plate_ids:
        c_ra = np.array([plates[p]["center_ra"] for p in plate_ids])
        c_dec = np.array([plates[p]["center_dec"] for p in plate_ids])
        CH = 20000
        for s in range(0, n, CH):
            e = min(s + CH, n)
            d = haversine_deg(
                ra_a[s:e, None], dec_a[s:e, None], c_ra[None, :], c_dec[None, :]
            )
            j = np.argmin(d, axis=1)
            sep_best[s:e] = d[np.arange(e - s), j]
            for k, jj in enumerate(j):
                best_plate[s + k] = plate_ids[jj]

    in_core_own = sep_own <= core_radius
    in_core_best = sep_best <= core_radius
    decide = in_core_own if policy == "own" else in_core_best
    # Unresolved rows are never silently cut: they are kept and counted.
    keep = np.ones(n, dtype=bool) if not cut else (decide | unresolved)

    for i in range(n):
        flags_w.writerow(
            {
                "src_id": src_ids[i],
                "ra": f"{ra_a[i]:.10f}",
                "dec": f"{dec_a[i]:.10f}",
                "tile_id": tile_ids[i],
                "det_plate": det_plates[i],
                "sep_own_deg": "" if np.isnan(sep_own[i]) else f"{sep_own[i]:.6f}",
                "in_core_own": int(bool(in_core_own[i])) if not np.isnan(sep_own[i]) else "",
                "best_plate": best_plate[i],
                "sep_best_deg": "" if np.isnan(sep_best[i]) else f"{sep_best[i]:.6f}",
                "in_core_best": int(bool(in_core_best[i])) if not np.isnan(sep_best[i]) else "",
                "cheb_own_deg": "" if np.isnan(cheb_own[i]) else f"{cheb_own[i]:.6f}",
                "edge_dist_arcmin": "" if np.isnan(edge_arcmin[i]) else f"{edge_arcmin[i]:.3f}",
                "beyond_corner": int(bool(beyond_corner[i])),
                "plate_unresolved": int(bool(unresolved[i])),
                "source_chunk": path.name,
            }
        )
        if keep[i]:
            kept_w.writerow({"src_id": src_ids[i], "ra": f"{ra_a[i]:.10f}", "dec": f"{dec_a[i]:.10f}"})

    yc_edge = {f"{t}": int(np.count_nonzero(edge_arcmin < t)) for t in EDGE_CURVE_ARCMIN}
    yc_own = {f"{r:.3f}": int(np.count_nonzero(sep_own > r)) for r in YIELD_CURVE_RADII}
    yc_best = {f"{r:.3f}": int(np.count_nonzero(sep_best > r)) for r in YIELD_CURVE_RADII}

    return ChunkStats(
        chunk=path.name,
        input_rows=n,
        kept_rows=int(np.count_nonzero(keep)),
        in_core_own=int(np.count_nonzero(in_core_own)),
        in_core_best=int(np.count_nonzero(in_core_best)),
        plate_unresolved=int(np.count_nonzero(unresolved)),
        beyond_corner=int(np.count_nonzero(beyond_corner)),
        missing_coords=missing_coords,
        yield_curve_own=yc_own,
        yield_curve_best=yc_best,
        yield_curve_edge=yc_edge,
    )


_GEOM = None  # set in main(); module-level so process_chunk can reach it


def self_test() -> int:
    """Geometry self-tests. No pipeline state, no disk reads beyond this file."""
    ok = True

    def check(name, cond):
        nonlocal ok
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")
        ok = ok and bool(cond)

    print("stage_edge_post_v2 self-tests")

    # haversine basics
    check("haversine identity == 0", abs(float(haversine_deg(10.0, 20.0, 10.0, 20.0))) < 1e-9)
    check("haversine 1 deg in dec", abs(float(haversine_deg(0.0, 0.0, 0.0, 1.0)) - 1.0) < 1e-9)
    check(
        "haversine shrinks with cos(dec) in RA",
        abs(float(haversine_deg(0.0, 60.0, 1.0, 60.0)) - 0.5) < 1e-3,
    )
    check(
        "haversine handles RA wrap at 0/360",
        abs(float(haversine_deg(359.5, 0.0, 0.5, 0.0)) - 1.0) < 1e-9,
    )

    # vectorised broadcasting used for the best-plate search
    d = haversine_deg(np.array([0.0, 0.0])[:, None], np.array([0.0, 0.0])[:, None],
                      np.array([0.0, 1.0])[None, :], np.array([1.0, 0.0])[None, :])
    check("broadcast shape (2,2)", d.shape == (2, 2))

    # the constant relationships this stage exists to get right
    check("2.7 == 5.4/2", abs(CORE_RADIUS_DEFAULT_DEG - 5.4 / 2.0) < 1e-12)
    check(
        "2.2 == 2.7 - half-tile margin",
        abs(2.2 - (5.4 / 2.0 - 0.5)) < 1e-12,
    )
    check(
        "vasco60 parity radius == 2.2 + half-diagonal",
        abs(VASCO60_PARITY_RADIUS_DEG - (2.2 + 0.5 * math.sqrt(2.0))) < 1e-12,
    )
    check("parity radius sits between 2.7 and 3.0", 2.7 < VASCO60_PARITY_RADIUS_DEG < 3.0)
    check(
        "core radius is inside the square plate half-width",
        CORE_RADIUS_DEFAULT_DEG < PLATE_HALF_WIDTH_DEG,
    )

    print("ALL PASS" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


def main() -> int:
    global _GEOM

    ap = argparse.ArgumentParser(
        description="[EXPERIMENTAL] Per-row plate-core (edge) stage, v2."
    )
    ap.add_argument("--self-test", action="store_true", help="Run geometry self-tests and exit")
    ap.add_argument("--run-dir", help="Run folder, e.g. ./work/runs/run-S1-...")
    ap.add_argument("--input-glob", default="stages/stage_S0.csv",
                    help="Glob (relative to run-dir) for input stage CSV")
    ap.add_argument("--stage", default="S1", help="Stage label used in output filenames")
    ap.add_argument("--out-dir", default=None, help="Output directory. Default: <run-dir>/stages")
    ap.add_argument("--headers-dir", default="/srv/vasco/vasco60/metadata/plates/headers",
                    help="Plate header registry (*.header.json with PLATERA/PLATEDEC)")
    ap.add_argument("--plate-map-csv", default=None,
                    help="CSV with tile_id,plate_id (e.g. fullscale_veto/tile_plate_map.csv)")
    ap.add_argument("--src-plate-csv", default=None,
                    help="CSV with src_id,det_plate (e.g. released primary_plate_flags.csv)")
    ap.add_argument("--plate-col", default=None,
                    help="Explicit plate column in the input CSV (default: autodetect)")
    ap.add_argument("--core-radius-deg", type=float, default=CORE_RADIUS_DEFAULT_DEG,
                    help=f"Radial core cut in degrees. Default {CORE_RADIUS_DEFAULT_DEG} "
                         "(APS 5.4 deg reliable-core diameter / 2)")
    ap.add_argument("--pasp2025", action="store_true",
                    help=f"Use {PASP2025_RADIUS_DEG} deg, the plate-centre mask of Villarroel "
                         "et al. 2025 (PASP). Their own criterion, not ours.")
    ap.add_argument("--vasco60-parity", action="store_true",
                    help=f"Use {VASCO60_PARITY_RADIUS_DEG:.4f} deg (2.2 + half-tile diagonal), "
                         "reproducing old VASCO60's effective per-row footprint")
    ap.add_argument("--policy", choices=["own", "any"], default="own",
                    help="Which separation --cut acts on: the detection plate (own) or the "
                         "nearest covering plate (any). Both are always recorded.")
    ap.add_argument("--cut", action="store_true",
                    help="Actually drop out-of-core rows. Default is flag-only (keeps everything).")
    ap.add_argument("--allow-unresolved", action="store_true",
                    help="Exit 0 even if some rows had no resolvable plate (default: exit 3)")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    if not args.run_dir:
        raise SystemExit("--run-dir is required (or use --self-test)")

    _GEOM = _load_geometry_module()

    if args.pasp2025 and args.vasco60_parity:
        raise SystemExit("--pasp2025 and --vasco60-parity are mutually exclusive")
    if args.pasp2025:
        core_radius = PASP2025_RADIUS_DEG
    elif args.vasco60_parity:
        core_radius = VASCO60_PARITY_RADIUS_DEG
    else:
        core_radius = float(args.core_radius_deg)

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise SystemExit(f"run-dir not found: {run_dir}")

    out_dir = Path(args.out_dir) if args.out_dir else (run_dir / "stages")
    out_dir.mkdir(parents=True, exist_ok=True)

    chunks = sorted(run_dir.glob(args.input_glob))
    if not chunks:
        raise SystemExit(f"No inputs matched: {run_dir}/{args.input_glob}")

    headers_dir = Path(args.headers_dir)
    if not headers_dir.exists():
        raise SystemExit(f"headers-dir not found: {headers_dir}")
    plates = load_plate_centres(headers_dir, _GEOM)
    if not plates:
        raise SystemExit(f"no plate headers loaded from {headers_dir}")

    tile_map = load_tile_plate_map(Path(args.plate_map_csv)) if args.plate_map_csv else {}
    src_map = load_src_plate_map(Path(args.src_plate_csv)) if args.src_plate_csv else {}

    hdr = _read_header(chunks[0])
    src_col, ra_col, dec_col, plate_col = _detect_cols(hdr, args.plate_col)
    if not plate_col and not tile_map and not src_map:
        raise SystemExit(
            "No way to resolve a plate per row: input has no plate column, and neither "
            "--plate-map-csv nor --src-plate-csv was given."
        )

    stage = args.stage
    out_kept = out_dir / f"stage_{stage}_EDGE2.csv"
    out_flags = out_dir / f"stage_{stage}_EDGE2_flags.csv"
    out_ledger = out_dir / f"stage_{stage}_EDGE2_ledger.json"

    flags_fields = ["src_id", "ra", "dec", "tile_id", "det_plate", "sep_own_deg", "in_core_own",
                    "best_plate", "sep_best_deg", "in_core_best", "cheb_own_deg",
                    "edge_dist_arcmin", "beyond_corner", "plate_unresolved",
                    "source_chunk"]

    per_chunk: List[ChunkStats] = []
    with out_kept.open("w", newline="", encoding="utf-8") as f_kept, \
         out_flags.open("w", newline="", encoding="utf-8") as f_flags:
        kept_w = csv.DictWriter(f_kept, fieldnames=["src_id", "ra", "dec"])
        flags_w = csv.DictWriter(f_flags, fieldnames=flags_fields)
        kept_w.writeheader()
        flags_w.writeheader()
        for ch in chunks:
            st = process_chunk(ch, src_col, ra_col, dec_col, plate_col, tile_map, src_map,
                               plates, core_radius, args.policy, args.cut, flags_w, kept_w)
            per_chunk.append(st)
            print(f"[EDGE2] {ch.name}: in={st.input_rows} kept={st.kept_rows} "
                  f"in_core_own={st.in_core_own} in_core_best={st.in_core_best} "
                  f"unresolved={st.plate_unresolved} beyond_corner={st.beyond_corner}", flush=True)

    tot_in = sum(s.input_rows for s in per_chunk)
    tot_kept = sum(s.kept_rows for s in per_chunk)
    tot_unres = sum(s.plate_unresolved for s in per_chunk)
    tot_corner = sum(s.beyond_corner for s in per_chunk)

    curve_own: Dict[str, int] = {}
    curve_best: Dict[str, int] = {}
    curve_edge: Dict[str, int] = {}
    for s in per_chunk:
        for k, v in s.yield_curve_own.items():
            curve_own[k] = curve_own.get(k, 0) + v
        for k, v in s.yield_curve_best.items():
            curve_best[k] = curve_best.get(k, 0) + v
        for k, v in s.yield_curve_edge.items():
            curve_edge[k] = curve_edge.get(k, 0) + v

    ledger = {
        "experimental": True,
        "stage": stage,
        "run_dir": str(run_dir),
        "input_glob": args.input_glob,
        "core_radius_deg": core_radius,
        "core_radius_source": (
            "PASP 2025 (Villarroel et al.) plate-centre mask, >2 deg discarded"
            if args.pasp2025 else
            "vasco60-parity: 2.2 + 0.5*sqrt(2), old VASCO60 effective per-row footprint"
            if args.vasco60_parity else
            "APS/MAPS reliable core: 5.4 deg diameter / 2" if abs(core_radius - 2.7) < 1e-9 else
            "operator-specified"
        ),
        "policy": args.policy,
        "cut_applied": bool(args.cut),
        "plate_half_width_deg": PLATE_HALF_WIDTH_DEG,
        "columns_detected": {"src_id": src_col, "ra": ra_col, "dec": dec_col,
                             "plate": plate_col or "(resolved via map)"},
        "plate_sources": {
            "headers_dir": str(headers_dir),
            "plates_loaded": len(plates),
            "plate_map_csv": args.plate_map_csv,
            "src_plate_csv": args.src_plate_csv,
            "centre_from_wcs": sum(
                1 for p in plates.values() if p.get("center_source", "").startswith("scan WCS")
            ),
            "centre_offset_gt_0p5_deg": sorted(
                (
                    {"plate": pid,
                     "offset_deg": round(pl["center_offset_deg"], 4),
                     "plate_ra": round(pl["plate_ra"], 4),
                     "plate_dec": round(pl["plate_dec"], 4),
                     "wcs_ra": round(pl["center_ra"], 4),
                     "wcs_dec": round(pl["center_dec"], 4)}
                    for pid, pl in plates.items()
                    if pl.get("center_offset_deg") == pl.get("center_offset_deg")
                    and pl.get("center_offset_deg", 0.0) > 0.5
                ),
                key=lambda d: -d["offset_deg"],
            ),
        },
        "totals": {
            "input_rows": tot_in,
            "kept_rows": tot_kept,
            "dropped_rows": tot_in - tot_kept,
            "in_core_own": sum(s.in_core_own for s in per_chunk),
            "in_core_best": sum(s.in_core_best for s in per_chunk),
            "plate_unresolved_rows": tot_unres,
            "beyond_corner_rows": tot_corner,
            "rows_missing_coords": sum(s.missing_coords for s in per_chunk),
        },
        "yield_curve_dropped_by_radius": {"own": curve_own, "best": curve_best},
        "yield_curve_within_arcmin_of_array_edge": curve_edge,
        "per_chunk": [s.__dict__ for s in per_chunk],
        "outputs": {
            "kept_csv": str(out_kept),
            "flags_csv": str(out_flags),
            "ledger_json": str(out_ledger),
        },
        "notes": [
            "Radial separation from PLATERA/PLATEDEC; the APS core is a circle, not the "
            "square plate boundary. cheb_own_deg is a diagnostic, not the cut.",
            "MAPS_CORE_RADIUS_DEG=2.2 in vasco/plan/tessellate_plates.py is a tile-CENTRE "
            "constant carrying a 0.5 deg half-tile margin; per-row the equivalent is 2.7.",
            "beyond_corner_rows are geometrically impossible for a correct plate label and "
            "indicate plate attribution defects, not edge effects.",
            "Rows with an unresolvable plate are KEPT and counted, never silently cut.",
        ],
    }
    out_ledger.write_text(json.dumps(ledger, indent=2), encoding="utf-8")

    print(f"[EDGE2] radius={core_radius:.4f} deg policy={args.policy} cut={args.cut}")
    print(f"[EDGE2] in={tot_in} kept={tot_kept} dropped={tot_in - tot_kept} "
          f"({(tot_in - tot_kept) / max(tot_in, 1):.2%})")
    print(f"[EDGE2] wrote: {out_kept}")
    print(f"[EDGE2] wrote: {out_flags}")
    print(f"[EDGE2] wrote: {out_ledger}")

    if tot_unres and not args.allow_unresolved:
        print(f"[EDGE2] ERROR: {tot_unres} rows had no resolvable plate. They were KEPT and "
              f"flagged. Re-run with --allow-unresolved to accept this.", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
