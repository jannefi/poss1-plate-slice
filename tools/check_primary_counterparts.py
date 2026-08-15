#!/usr/bin/env python3
"""Complete the primary-plate partition: does the primary plate's own raw
detection list corroborate each non-primary row?

Adds one column to the flags written by tools/build_primary_plate_flags.py:

    primary_has_det  True if the primary plate holds a raw detection within
                     --radius-arcsec of the row's position

Rows with is_primary=False and primary_has_det=False are single-plate content
in multiply-searched sky -- present on one plate's pixels only. A null control
(every query displaced 6' in RA against the same plate) calibrates the chance
rate; the real rate landing BELOW the null is expected, not anomalous, because
catalogue rows are veto survivors and their sky is therefore star-depleted.

Input --radec-dir is the per-plate detection output of the pipeline's steps
2-3 (one <PLATE>.csv with ra,dec per detection, as written by
tools/run_fullscale_slice.py). It is regenerable from the public plate scans,
so a third party can reproduce this column end to end.

Measured on results/s0-642-20260814 (54,749 non-primary rows):
primary_has_det = 0.22%, null 2.68%.

Usage:
    python3 tools/check_primary_counterparts.py \\
        --flags primary_plate_flags.csv \\
        --catalog results/s0-642-20260814/stage_S0.csv.gz \\
        --radec-dir <steps-2-3 output>/radec \\
        --out primary_plate_flags.csv

How to validate: rerun with a different --null-shift-arcmin; the real rate must
not move, the null rate stays in the same few-percent range.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def _unit(ra, dec):
    ra = np.deg2rad(np.asarray(ra, float))
    dec = np.deg2rad(np.asarray(dec, float))
    return np.column_stack(
        [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--catalog", required=True,
                    help="the catalogue the flags were built from (for ra/dec)")
    ap.add_argument("--radec-dir", required=True,
                    help="per-plate raw detections: <PLATE>.csv with ra,dec")
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--null-shift-arcmin", type=float, default=6.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    fl = pd.read_csv(args.flags)
    cat = pd.read_csv(args.catalog, usecols=["src_id", "ra", "dec"])
    d = fl.merge(cat, on="src_id", how="left")
    if len(d) != len(fl) or d.ra.isna().any():
        sys.exit("[FAIL] flags do not join 1:1 onto the catalogue")

    chord = 2 * np.sin(np.deg2rad(args.radius_arcsec / 3600.0) / 2)
    todo = d[~d.is_primary]
    hit = np.zeros(len(d), bool)
    null = np.zeros(len(d), bool)
    t0, missing = time.time(), []
    for plate, g in todo.groupby("primary_plate"):
        f = Path(args.radec_dir) / f"{plate}.csv"
        if not f.is_file():
            missing.append(plate)
            continue
        det = pd.read_csv(f, usecols=lambda c: c.lower() in ("ra", "dec"))
        cols = {c.lower(): c for c in det.columns}
        tree = cKDTree(_unit(det[cols["ra"]].values, det[cols["dec"]].values))
        loc = d.index.get_indexer(g.index)
        dd, _ = tree.query(_unit(g.ra, g.dec), distance_upper_bound=chord)
        hit[loc] = np.isfinite(dd)
        shift = args.null_shift_arcmin / 60.0 / np.cos(np.deg2rad(g.dec))
        dn, _ = tree.query(_unit(g.ra + shift, g.dec), distance_upper_bound=chord)
        null[loc] = np.isfinite(dn)

    d["primary_has_det"] = hit
    n_np = int((~d.is_primary).sum())
    out = d[["src_id", "tile_id", "det_plate", "primary_plate", "is_primary",
             "primary_has_det", "sep_primary_deg", "sep_margin"]]
    out.to_csv(args.out, index=False)

    single = int(((~d.is_primary) & ~d.primary_has_det).sum())
    print(f"[OK] {len(d):,} rows -> {args.out}  ({time.time()-t0:.0f}s)")
    print(f"     non-primary rows              : {n_np:,}")
    print(f"     primary_has_det               : {hit.sum():,} "
          f"({100*hit.sum()/max(n_np,1):.2f}%)   null {100*null.sum()/max(n_np,1):.2f}%")
    print(f"     single-plate in overlap sky   : {single:,}")
    if missing:
        print(f"[WARN] {len(missing)} primary plates had no radec CSV: "
              f"{missing[:5]} -- their rows carry primary_has_det=False by absence")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
