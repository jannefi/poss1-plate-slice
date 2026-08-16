#!/usr/bin/env python3
"""Do plate-rim catalogue rows have a counterpart on the neighbouring plate?

The coverage partition rests on a structural claim: rim content is found on one
plate's vignetted edge and is invisible to a pipeline that requests cutouts by
position, because such a pipeline is served the *neighbouring* plate's
well-exposed interior for that sky. The claim is testable — if it holds, rim
rows should almost never have a raw detection on the neighbour.

It is also a robustness check on tools/check_primary_counterparts.py, which
asks the same question under a different geometric rule (is the *primary* plate
-- the one whose centre the position is nearest -- carrying a detection?). Two
independent definitions of "content a per-position design cannot reach" should
give the same rate; if they do not, the partition is an artifact of how it was
drawn. Pass --ref-csv to split the sample by presence in some other catalogue
under the identical rule, which is the only way to check that a
catalogue-conditioned selection did not bias the rate -- comparing across two
different partition rules confounds the rule with the selection.

Definitions, stated because the original script was not retained and its `cheb`
column could not be reproduced (it reaches 3.805°, beyond a plate half-width):

  cheb       distance from the plate centre in **the plate's own pixel frame**,
             max(|x-xc|, |y-yc|) scaled to degrees. A POSS-I plate is square in
             pixel space and rotated with respect to RA/Dec, so a chebyshev
             taken in RA/Dec is not the plate's own geometry. Capped at the
             half-width, ~3.30°, by construction.
  rim        cheb > --rim-deg (default 3.0°)
  neighbour  the nearest plate centre *other than* the detection plate, from
             the public plate manifest
  null       the identical query displaced --null-shift-arcmin in RA against
             the same neighbour plate, which preserves that plate's local
             source density and destroys any real association

The real rate landing BELOW the null is expected, not anomalous: catalogue rows
are veto survivors, so their sky is star-depleted relative to a shifted control.
The discriminating quantity is the real rate, not the difference.

How to validate
---------------
    python3 tools/rim_neighbour_counterparts.py \\
        --catalog results/s0-642-20260814/stage_S0.csv.gz \\
        --flags   results/s0-642-20260814/primary_plate_flags.csv.gz \\
        --plate-manifest data/plate_manifest.csv \\
        --plate-dir <plate-scan-dir> \\
        --radec-dir <run>/radec \\
        --out-dir work/rim_neighbour

Sanity: cheb must not exceed the plate half-width for any row; the tool asserts
it. Re-running with a different --null-shift-arcmin must move the null rate
only within its interval and leave the real rate untouched.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree


def unit(ra, dec) -> np.ndarray:
    ra = np.radians(np.asarray(ra, dtype=np.float64))
    dec = np.radians(np.asarray(dec, dtype=np.float64))
    c = np.cos(dec)
    return np.column_stack((c * np.cos(ra), c * np.sin(ra), np.sin(dec)))


def chord(deg: float) -> float:
    return 2.0 * np.sin(np.radians(deg) / 2.0)


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * (c - h), 100 * (c + h))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--flags", required=True, help="primary_plate_flags.csv[.gz], for det_plate")
    ap.add_argument("--plate-manifest", required=True)
    ap.add_argument("--plate-dir", required=True)
    ap.add_argument("--radec-dir", required=True)
    ap.add_argument("--ref-csv", default=None,
                    help="OPTIONAL comparison catalogue; splits the sample by presence "
                         "in it under the identical rule, to check that conditioning on a "
                         "catalogue does not shift the rate. Never hardcode its path.")
    ap.add_argument("--ref-radius-arcsec", type=float, default=5.0)
    ap.add_argument("--rim-deg", type=float, default=3.0)
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--null-shift-arcmin", type=float, default=6.0)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    cat = pd.read_csv(args.catalog, usecols=["src_id", "ra", "dec"])
    flags = pd.read_csv(args.flags, usecols=["src_id", "det_plate"])
    cat = cat.merge(flags, on="src_id", how="left")
    if cat["det_plate"].isna().any():
        print(f"[FATAL] {int(cat.det_plate.isna().sum())} rows have no det_plate", file=sys.stderr)
        return 2
    pm = pd.read_csv(args.plate_manifest)
    print(f"{len(cat):,} catalogue rows, {cat.det_plate.nunique()} detection plates, "
          f"{len(pm)} plates in the manifest")

    # ---- cheb in each plate's own pixel frame ----
    cheb = np.full(len(cat), np.nan)
    halfwidth = {}
    for plate, g in cat.groupby("det_plate"):
        p = Path(args.plate_dir) / f"dss1red_{plate}.fits"
        if not p.exists():
            continue
        h = fits.getheader(p)
        w = WCS(h)
        nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
        sc = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))  # deg/px
        x, y = w.all_world2pix(g["ra"].to_numpy(), g["dec"].to_numpy(), 0)
        cheb[g.index] = np.maximum(np.abs(x - (nx - 1) / 2), np.abs(y - (ny - 1) / 2)) * sc
        halfwidth[plate] = max(nx, ny) / 2 * sc
    cat = cat.assign(cheb=cheb)

    hw = max(halfwidth.values()) if halfwidth else 3.4
    over = cat["cheb"] > hw + 1e-6
    if over.any():
        print(f"[FATAL] {int(over.sum())} rows have cheb beyond the plate half-width "
              f"({hw:.3f}°); the frame is wrong", file=sys.stderr)
        return 2
    print(f"cheb computed for {int(cat.cheb.notna().sum()):,} rows; "
          f"max {np.nanmax(cat.cheb):.3f}° against a half-width of {hw:.3f}°")

    rim = cat[cat["cheb"] > args.rim_deg].copy()
    print(f"rim rows (cheb > {args.rim_deg:g}°): {len(rim):,} "
          f"({100*len(rim)/len(cat):.1f}% of the catalogue)")

    # ---- nearest plate other than the detection plate ----
    ptree = cKDTree(unit(pm["ra_deg"], pm["dec_deg"]))
    plates = pm["plate_id"].to_numpy()
    d2, j2 = ptree.query(unit(rim["ra"], rim["dec"]), k=2)
    nearest = plates[j2[:, 0]]
    neigh = np.where(nearest == rim["det_plate"].to_numpy(), plates[j2[:, 1]], nearest)
    rim = rim.assign(neigh=neigh)

    # ---- does the neighbour hold a raw detection? real and null ----
    shift = args.null_shift_arcmin / 60.0
    hit = np.zeros(len(rim), dtype=bool)
    nul = np.zeros(len(rim), dtype=bool)
    seen = np.zeros(len(rim), dtype=bool)
    r = chord(args.radius_arcsec / 3600.0)
    groups = list(rim.groupby("neigh"))
    for n, (plate, g) in enumerate(groups, start=1):
        f = Path(args.radec_dir) / f"{plate}.csv"
        if not f.exists():
            continue
        det = pd.read_csv(f, usecols=["ra", "dec"])
        if det.empty:
            continue
        t = cKDTree(unit(det["ra"].to_numpy(), det["dec"].to_numpy()))
        pos = rim.index.get_indexer(g.index)
        seen[pos] = True
        dd, _ = t.query(unit(g["ra"], g["dec"]), k=1)
        hit[pos] = dd <= r
        dn, _ = t.query(unit(g["ra"] + shift / np.cos(np.radians(g["dec"])), g["dec"]), k=1)
        nul[pos] = dn <= r
        if n % 100 == 0 or n == len(groups):
            print(f"  {n}/{len(groups)} neighbour plates", flush=True)
    rim = rim.assign(neigh_searched=seen, neigh_has_det=hit, null_has_det=nul)
    rim.to_csv(out / "rim_rows.csv", index=False)

    def report(label, sub):
        s = sub[sub["neigh_searched"]]
        k, n = int(s["neigh_has_det"].sum()), len(s)
        kn = int(s["null_has_det"].sum())
        lo, hi = wilson(k, n)
        nlo, nhi = wilson(kn, n)
        print(f"\n  {label}")
        print(f"    rim rows searched            {n:>7,}")
        print(f"    neighbour detection <= {args.radius_arcsec:g}\"   {k:>7,}   "
              f"{100*k/n:.3f}%  95% CI [{lo:.3f}, {hi:.3f}]")
        print(f"    null control ({args.null_shift_arcmin:g}' shift)   {kn:>7,}   "
              f"{100*kn/n:.3f}%  95% CI [{nlo:.3f}, {nhi:.3f}]")
        return dict(label=label, n=n, hits=k, pct=100 * k / n, ci=[lo, hi],
                    null=kn, null_pct=100 * kn / n)

    print(f"\n=== NEIGHBOUR-PLATE COUNTERPART RATE (rim > {args.rim_deg:g}°) ===")
    res = [report("ALL rim rows (unconditioned)", rim)]

    if args.ref_csv:
        ref = pd.read_csv(args.ref_csv)
        rc = [c for c in ref.columns if c.lower() in ("ra", "_ra", "raj2000", "ra_deg")][0]
        dc = [c for c in ref.columns if c.lower() in ("dec", "_de", "dej2000", "dec_deg")][0]
        rt = cKDTree(unit(ref[rc], ref[dc]))
        dd, _ = rt.query(unit(rim["ra"], rim["dec"]), k=1)
        absent = dd > chord(args.ref_radius_arcsec / 3600.0)
        rim = rim.assign(ref_absent=absent)
        res.append(report("rim rows ABSENT from the reference (private)", rim[absent]))
        res.append(report("rim rows PRESENT in the reference (private)", rim[~absent]))
        rim.to_csv(out / "rim_rows.csv", index=False)

    (out / "summary.json").write_text(json.dumps(
        {"rim_deg": args.rim_deg, "radius_arcsec": args.radius_arcsec,
         "null_shift_arcmin": args.null_shift_arcmin,
         "n_catalogue": int(len(cat)), "n_rim": int(len(rim)), "results": res}, indent=2))
    print(f"\nledgers -> {out}/rim_rows.csv, summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
