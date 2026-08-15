#!/usr/bin/env python3
"""Neighbour-plate persistence check: is a candidate's source visible on a
DIFFERENT plate, i.e. at a different epoch?

SHIPPED BUT NOT EXECUTED. This stage is implemented and its impact measured,
but it is deliberately not part of the released catalogue's chain — see
docs/PARAMETERS.md ("Stages that ship but are NOT run") for the measured
numbers and the reason. In short: a row this stage flags is a source seen on
two plates exposed on different nights, which is strong evidence of a
persistent object rather than a transient. Removing such rows would change
what the catalogue claims — from "on the plate and absent from the modern
catalogues we checked" to an assertion about transience — and that is a
different catalogue, to be built deliberately or not at all.

POSS-I red plates overlap their neighbours by roughly half a degree, and
neighbouring fields were photographed on different nights, often different
years. For a candidate inside overlap sky, the neighbouring plate is therefore
a free second epoch. This tool checks every catalogue row against the raw
detections of every OTHER plate whose footprint covers its position, with a
displaced-null control per plate pair.

Only overlap rows are decidable: a row covered by a single plate has no second
epoch to consult, and is reported as 'single_coverage', not as unconfirmed.

Usage:
    python3 tools/stage_neighbor_persistence.py \\
        --catalog results/s0-642-20260814/stage_S0.csv.gz \\
        --flags results/s0-642-20260814/primary_plate_flags.csv.gz \\
        --manifest data/plate_manifest.csv \\
        --radec-dir <steps-2-3 output>/radec \\
        --out neighbor_persistence.csv

How to validate: rerun with a different --null-shift-arcmin; the confirmed
rate must be stable while the null stays in the few-percent range. The
half-width guard can be checked against any plate header (14000 px at
~1.7"/px is ~6.6 deg, half-width ~3.3).
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

PLATE_HALF_WIDTH_DEG = 3.3   # 14000 px * ~1.7"/px / 2; chebyshev, per axis


def _unit(ra, dec):
    ra = np.deg2rad(np.asarray(ra, float))
    dec = np.deg2rad(np.asarray(dec, float))
    return np.column_stack(
        [np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--flags", required=True,
                    help="primary_plate_flags CSV (for det_plate)")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--radec-dir", required=True)
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--null-shift-arcmin", type=float, default=6.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    cat = pd.read_csv(args.catalog, usecols=["src_id", "ra", "dec"])
    fl = pd.read_csv(args.flags, usecols=["src_id", "det_plate"])
    d = cat.merge(fl, on="src_id")
    man = pd.read_csv(args.manifest)

    chord = 2 * np.sin(np.deg2rad(args.radius_arcsec / 3600.0) / 2)
    n_other = np.zeros(len(d), np.int16)      # other plates covering the row
    n_conf = np.zeros(len(d), np.int16)       # of those, plates that detect it
    n_null = np.zeros(len(d), np.int16)
    t0, missing = time.time(), 0

    # One pass per plate: find rows inside THIS plate's footprint that were
    # detected elsewhere, and ask this plate's raw detections about them.
    for _, p in man.iterrows():
        dra = ((d.ra.values - p.ra_deg + 180.0) % 360.0) - 180.0
        cheb = np.maximum(np.abs(dra * np.cos(np.deg2rad(d.dec.values))),
                          np.abs(d.dec.values - p.dec_deg))
        sel = (cheb <= PLATE_HALF_WIDTH_DEG) & (d.det_plate.values != p.plate_id)
        if not sel.any():
            continue
        f = Path(args.radec_dir) / f"{p.plate_id}.csv"
        if not f.is_file():
            missing += 1
            continue
        det = pd.read_csv(f, usecols=lambda c: c.lower() in ("ra", "dec"))
        cols = {c.lower(): c for c in det.columns}
        tree = cKDTree(_unit(det[cols["ra"]].values, det[cols["dec"]].values))
        idx = np.flatnonzero(sel)
        g_ra, g_dec = d.ra.values[idx], d.dec.values[idx]
        n_other[idx] += 1
        dd, _ = tree.query(_unit(g_ra, g_dec), distance_upper_bound=chord)
        n_conf[idx] += np.isfinite(dd).astype(np.int16)
        shift = args.null_shift_arcmin / 60.0 / np.cos(np.deg2rad(g_dec))
        dn, _ = tree.query(_unit(g_ra + shift, g_dec), distance_upper_bound=chord)
        n_null[idx] += np.isfinite(dn).astype(np.int16)

    d["n_other_plates"] = n_other
    d["n_confirming"] = n_conf
    d["persistent"] = n_conf > 0
    d[["src_id", "det_plate", "n_other_plates", "n_confirming",
       "persistent"]].to_csv(args.out, index=False)

    n = len(d)
    dec_able = n_other > 0
    conf = d.persistent.values
    print(f"[OK] {n:,} rows -> {args.out}  ({time.time()-t0:.0f}s)"
          + (f"  [{missing} plates lacked a radec CSV]" if missing else ""))
    print(f"     single_coverage (no other plate)  : {(~dec_able).sum():,} "
          f"({100*(~dec_able).mean():.1f}%)")
    print(f"     overlap rows                      : {dec_able.sum():,}")
    print(f"     confirmed on another epoch        : {conf.sum():,} "
          f"= {100*conf[dec_able].mean():.2f}% of overlap rows")
    exp_null = (n_null[dec_able] > 0).mean()
    print(f"     null (displaced {args.null_shift_arcmin:.0f}')            "
          f"  : {100*exp_null:.2f}%")
    print("[NOTE] NOT part of the release chain. A 'persistent' row is a "
          "candidate seen at two epochs -- removing it changes the catalogue's "
          "claim. See docs/PARAMETERS.md before wiring this in anywhere.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
