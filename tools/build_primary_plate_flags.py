#!/usr/bin/env python3
"""Flag catalogue rows by the primary-plate rule: was this row detected on the
plate a per-position query would serve for its coordinates?

WHY THIS EXISTS

This pipeline slices full plate scans, so sky covered by more than one plate is
searched on every plate that covers it. A cutout-based pipeline (the design of
Solano et al. 2022, which fetched 60'x60' cutouts from a DSS service) searches
each queried position once, on the plate the DSS selection function serves.
Full-plate slicing is therefore a *coverage deviation* from that design, and
this tool marks it row by row instead of hiding it or deleting it:

    primary_plate  the plate a per-position query would serve for (ra, dec)
    is_primary     detection plate == primary plate

THE RULE, AND WHY IT HAS NO FREE PARAMETER

primary_plate = the plate with the nearest centre (great-circle; centres are
transcribed from the GSSS headers into the public plate manifest). Nothing is
fitted and nothing is tuned. Validated against 11,727 archive tiles whose
headers record which plate the STScI cutout service actually served at that
position: **nearest-centre matches for 99.04%**, and the ~1% disagreements are
near-equidistant boundary ties (the served plate is the second-nearest centre
in 94% of them).

HOW TO READ THE FLAG -- and how not to

is_primary=False does NOT mean "a cutout pipeline could not have found this".
Cutout pipelines query at *tile centres*, and for sky near the boundary between
two plates, plate(tile centre) routinely differs from plate(source position) --
a measured ~15% effect. So boundary sky is served from effectively either
plate, depending on the grid; only sky deep in a plate's rim, far from any
other plate's centre, is out of a cutout design's reach. The companion tool
`check_primary_counterparts.py` adds the column that completes the picture:

    primary_has_det  the primary plate's raw detections contain a source
                     within 5" of this row

Rows with is_primary=False and primary_has_det=False are **single-plate content
in multiply-searched sky**: present on one plate's pixels only, found by a
full-plate search with certainty and by a cutout design only when its tiling
happens to serve that plate.

THIS IS A PARTITION, NOT A QUALITY CUT. Filtering a catalogue to
is_primary=True discards real content: measured against the public
vanish-possi catalogue (R, 5,399 rows), a hard is_primary filter loses 9.1% of
R matches (98 of 1,072). Quote the partition counts side by side; filter only
with that cost in view.

Usage (self-contained from a release folder):
    python3 tools/build_primary_plate_flags.py \\
        --catalog results/s0-642-20260814/stage_S0.csv.gz \\
        --tile-plate-map <(zcat results/s0-642-20260814/tile_manifest.csv.gz) \\
        --manifest data/plate_manifest.csv \\
        --out primary_plate_flags.csv

How to validate: the output row count must equal the catalogue exactly (the
tool asserts it), and any row can be checked by hand -- its primary_plate is
the manifest plate with the smallest angular distance to (ra, dec).
"""
from __future__ import annotations

import argparse
import sys

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
    ap.add_argument("--catalog", required=True,
                    help="CSV(.gz) with src_id, tile_id, ra, dec")
    ap.add_argument("--tile-plate-map", required=True,
                    help="CSV with tile_id, plate_id -- the detection plate "
                         "(a release's tile_manifest works)")
    ap.add_argument("--manifest", required=True,
                    help="plate manifest with plate_id, ra_deg, dec_deg")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    cat = pd.read_csv(args.catalog)
    for c in ("src_id", "tile_id", "ra", "dec"):
        if c not in cat.columns:
            sys.exit(f"[FAIL] catalogue lacks required column '{c}'")
    n_in = len(cat)

    man = pd.read_csv(args.manifest)
    if not {"plate_id", "ra_deg", "dec_deg"} <= set(man.columns):
        sys.exit("[FAIL] manifest lacks ra_deg/dec_deg -- regenerate it with "
                 "tools/build_plate_manifest.py (centres were added 2026-08-15)")

    tpm = pd.read_csv(args.tile_plate_map)[["tile_id", "plate_id"]] \
        .rename(columns={"plate_id": "det_plate"})
    cat = cat.merge(tpm, on="tile_id", how="left")
    if cat.det_plate.isna().any():
        miss = cat[cat.det_plate.isna()].tile_id.unique()
        sys.exit(f"[FAIL] {len(miss)} tile_ids missing from the tile-plate map, "
                 f"e.g. {miss[:3]} -- refusing to emit partial flags")

    tree = cKDTree(_unit(man.ra_deg.values, man.dec_deg.values))
    dist, idx = tree.query(_unit(cat.ra.values, cat.dec.values), k=2)
    cat["primary_plate"] = man.plate_id.values[idx[:, 0]]
    cat["is_primary"] = cat.det_plate == cat.primary_plate
    # chord -> angle, degrees. sep_margin near 0 marks near-tie boundary sky,
    # where the nearest-centre proxy is least certain (the validated ~1%).
    a0 = np.rad2deg(2 * np.arcsin(np.clip(dist[:, 0] / 2, 0, 1)))
    a1 = np.rad2deg(2 * np.arcsin(np.clip(dist[:, 1] / 2, 0, 1)))
    cat["sep_primary_deg"] = np.round(a0, 4)
    cat["sep_margin"] = np.round(1.0 - a0 / np.maximum(a1, 1e-9), 4)

    out = cat[["src_id", "tile_id", "det_plate", "primary_plate", "is_primary",
               "sep_primary_deg", "sep_margin"]]
    assert len(out) == n_in, "row count changed -- must never happen"
    out.to_csv(args.out, index=False)

    n_p = int(out.is_primary.sum())
    print(f"[OK] {n_in:,} rows -> {args.out}")
    print(f"     is_primary=True : {n_p:,} ({100*n_p/n_in:.2f}%)")
    print(f"     is_primary=False: {n_in-n_p:,} ({100*(n_in-n_p)/n_in:.2f}%)")
    print("[NEXT] complete the partition with tools/check_primary_counterparts.py "
          "(adds primary_has_det from the per-plate raw detections)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
