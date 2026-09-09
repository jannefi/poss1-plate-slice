#!/usr/bin/env python3
"""Re-derive per-row (src_id-level) red/blue z-scores from ALREADY-FETCHED
cutouts written by blue_plate_hot_tile_check.py, for the follow-up question:
does the SCOS+PTF+EDGE post-process chain remove the blue-negative rows?

blue_plate_hot_tile_check.py only wrote per-tile summaries. This reads the
same cached FITS (no re-fetch) and the same S0 rows, and writes one row per
S0 survivor with its own red_z/blue_z -- the thing needed to check a
specific row's fate through an unrelated stage.

Usage:
    python3 tools/blue_plate_per_row_export.py \\
      --s0-csv results/s0-642-20260814/stage_S0.csv.gz \\
      --tiles tile_RA...,tile_RA... \\
      --fits-dir work/blue_plate_hot_tile_check_20260828/fits \\
      --out work/blue_plate_hot_tile_check_20260828/per_row.csv
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(REPO))

from tools.blue_plate_hot_tile_check import local_z, resolve_polarity  # noqa: E402


def main():
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--s0-csv", required=True)
    ap.add_argument("--tiles", required=True)
    ap.add_argument("--fits-dir", required=True)
    ap.add_argument("--gaia-cache", default="/home/janne/local_cache/gaia")
    ap.add_argument("--box", type=int, default=15)
    ap.add_argument("--ap", type=int, default=3)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    from astropy.io import fits
    from astropy.wcs import WCS
    from vasco.utils.tile_id import parse_tile_id_center

    tiles = [t.strip() for t in args.tiles.split(",") if t.strip()]
    s0 = pd.read_csv(args.s0_csv)
    fits_dir = Path(args.fits_dir)
    out_rows = []

    for tid in tiles:
        rows = s0[s0.tile_id == tid].reset_index(drop=True)
        if rows.empty:
            print(f"[SKIP] {tid}: 0 rows in S0")
            continue
        tra, tdec = parse_tile_id_center(tid)
        qra = f"{tra:.3f}"
        qdec = f"{tdec:.3f}"
        band_z = {}
        for band, tag in [("red", "dss1-red"), ("blue", "dss1-blue")]:
            cands = list(fits_dir.glob(f"{tag}_{qra}_{qdec}_*arcmin.fits"))
            if not cands:
                print(f"[MISSING] {tid} {band}: no cached FITS matching "
                      f"{tag}_{qra}_{qdec}_*")
                band_z[band] = np.full(len(rows), np.nan)
                continue
            hdr = fits.getheader(cands[0])
            data = fits.getdata(cands[0]).astype(np.float64)
            w = WCS(hdr)
            sign, _, _ = resolve_polarity(data, w, args.gaia_cache, tra, tdec,
                                          args.box, args.ap)
            if sign < 0:
                data = -data
            px, py = w.world_to_pixel_values(rows.ra.to_numpy(float),
                                             rows.dec.to_numpy(float))
            pxr, pyr = np.round(px).astype(int), np.round(py).astype(int)
            nx, ny = int(hdr["NAXIS1"]), int(hdr["NAXIS2"])
            on_array = (px >= 0) & (px < nx) & (py >= 0) & (py < ny)
            z, _ = local_z(data, pxr, pyr, args.box, args.ap)
            z[~on_array] = np.nan
            band_z[band] = z
        for i, r in rows.iterrows():
            out_rows.append({
                "src_id": r.src_id, "tile_id": tid, "object_id": r.object_id,
                "ra": r.ra, "dec": r.dec,
                "red_z": band_z["red"][i], "blue_z": band_z["blue"][i],
            })
        print(f"[OK] {tid}: {len(rows)} rows")

    df = pd.DataFrame(out_rows)
    df.to_csv(args.out, index=False)
    print(f"\n[out] {args.out}: {len(df)} rows")


if __name__ == "__main__":
    main()
