#!/usr/bin/env python3
"""Does catalogue content stay real as you move out to the plate rim?

The coverage partition says rim rows are single-plate content a per-position
pipeline cannot see. A selected sample of outer-rim rows suggested the opposite reading — that much
of the outer rim is scan artifacts — but a sample selected on a correlate of
artifact density cannot carry that extrapolation. This tool settles it on a
random sample stratified by rim depth.

For each sampled row it measures **peak significance at the catalogued position
in our own pixels**, (peak - annulus median) / annulus MAD-sigma on the plate
the row was detected on, and pairs it with the SuperCOSMOS R1 verdict.

Every row also gets a **displaced control** at the same position shifted in RA.
This is what makes "has signal" a measurement rather than an assumption: the
control gives the peak-sigma distribution of nearby blank sky on the same plate
at the same depth, so the excess over control is the real detection rate. Rim
pixels are noisier and more structured than interior pixels, so a fixed sigma
threshold alone would manufacture a depth trend out of nothing.

Interior strata are included deliberately. They anchor the scale: whatever the
method reports for rows everyone agrees are real is the ceiling against which
rim strata should be read.

How to validate
---------------
    python3 tools/rim_depth_profile.py \\
        --catalog results/s0-642-20260814/stage_S0.csv.gz \\
        --flags   results/s0-642-20260814/primary_plate_flags.csv.gz \\
        --plate-dir <plate-scan-dir> \\
        --scos-flags <chain>/stages/scos/scos_flags.csv \\
        --per-stratum 150 --seed 20260816 \\
        --out-dir work/rim_depth

Sanity: the innermost stratum must show a high real rate and a low control
rate. If it does not, the estimator is broken and no stratum means anything.
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

BINS = [0.0, 1.0, 2.0, 2.5, 3.0, 3.05, 3.10, 3.15, 3.20, 3.25, 3.35]
SIGMA = 3.0


def peaks_on_plate(path: Path, ra, dec, r_src=2.0, r_in=8.0, r_out=20.0):
    """Peak significance for many positions on one plate. Opens the file once."""
    out = np.full(len(ra), np.nan)
    if not path.exists():
        return out
    with fits.open(path, memmap=True) as h:
        hdu = h[0]
        w = WCS(hdu.header)
        ny, nx = hdu.shape
        sc = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600
        R = int(r_out / sc) + 3
        try:
            xs, ys = w.all_world2pix(np.asarray(ra), np.asarray(dec), 0)
        except Exception:
            return out
        for i, (x, y) in enumerate(zip(np.atleast_1d(xs), np.atleast_1d(ys))):
            x, y = float(x), float(y)
            x0, x1 = int(x) - R, int(x) + R
            y0, y1 = int(y) - R, int(y) + R
            if x0 < 0 or y0 < 0 or x1 > nx or y1 > ny:
                continue          # clipped: annulus would be biased
            d = np.asarray(hdu.section[y0:y1, x0:x1], float)
            yy, xx = np.mgrid[y0:y1, x0:x1]
            rr = np.hypot(xx - x, yy - y) * sc
            ann, core = d[(rr >= r_in) & (rr <= r_out)], d[rr <= r_src]
            if ann.size < 20 or core.size < 1:
                continue
            bg = np.median(ann)
            sig = 1.4826 * np.median(np.abs(ann - bg))
            if sig > 0:
                out[i] = (core.max() - bg) / sig
    return out


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    hh = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * (c - hh), 100 * (c + hh))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--plate-dir", required=True)
    ap.add_argument("--scos-flags", default=None)
    ap.add_argument("--per-stratum", type=int, default=150)
    ap.add_argument("--null-shift-arcmin", type=float, default=6.0)
    ap.add_argument("--seed", type=int, default=20260816)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    cat = pd.read_csv(args.catalog, usecols=["src_id", "ra", "dec"]).merge(
        pd.read_csv(args.flags, usecols=["src_id", "det_plate"]), on="src_id", how="left")
    print(f"{len(cat):,} catalogue rows")

    # ---- cheb in each plate's own pixel frame, all rows ----
    cheb = np.full(len(cat), np.nan)
    for plate, g in cat.groupby("det_plate"):
        p = Path(args.plate_dir) / f"dss1red_{plate}.fits"
        if not p.exists():
            continue
        h = fits.getheader(p)
        w = WCS(h)
        nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
        sc = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))
        x, y = w.all_world2pix(g["ra"].to_numpy(), g["dec"].to_numpy(), 0)
        cheb[g.index] = np.maximum(np.abs(x - (nx - 1) / 2), np.abs(y - (ny - 1) / 2)) * sc
    cat = cat.assign(cheb=cheb)
    cat["stratum"] = pd.cut(cat["cheb"], BINS)
    print("population per stratum:")
    print(cat["stratum"].value_counts().sort_index().to_string())

    rng = np.random.default_rng(args.seed)
    samp = (cat.dropna(subset=["cheb"]).groupby("stratum", observed=True, group_keys=False)
            .apply(lambda g: g.sample(min(len(g), args.per_stratum), random_state=args.seed)))
    samp = samp.reset_index(drop=True)
    print(f"\nsampled {len(samp):,} rows, seed {args.seed}")

    # ---- peak sigma, real and displaced, one plate open per plate ----
    shift = args.null_shift_arcmin / 60.0
    real = np.full(len(samp), np.nan)
    ctrl = np.full(len(samp), np.nan)
    plates = sorted(samp["det_plate"].dropna().unique())
    for n, plate in enumerate(plates, start=1):
        idx = np.where(samp["det_plate"].to_numpy() == plate)[0]
        p = Path(args.plate_dir) / f"dss1red_{plate}.fits"
        ra = samp["ra"].to_numpy()[idx]
        dec = samp["dec"].to_numpy()[idx]
        real[idx] = peaks_on_plate(p, ra, dec)
        ctrl[idx] = peaks_on_plate(p, ra + shift / np.cos(np.radians(dec)), dec)
        if n % 100 == 0 or n == len(plates):
            print(f"  {n}/{len(plates)} plates", flush=True)
    samp = samp.assign(peak=real, peak_ctrl=ctrl)

    if args.scos_flags:
        sc = pd.read_csv(args.scos_flags, usecols=["src_id", "nmatch_r1"])
        samp = samp.merge(sc, on="src_id", how="left")
        samp["scos_unconf"] = samp["nmatch_r1"] == 0
    samp.to_csv(out / "sample.csv", index=False)

    ok = samp[samp["peak"].notna() & samp["peak_ctrl"].notna()].copy()
    ok["has_sig"] = ok["peak"] > SIGMA
    ok["ctrl_sig"] = ok["peak_ctrl"] > SIGMA

    print(f"\n=== SIGNAL AND SuperCOSMOS BY RIM DEPTH  (n={len(ok):,} measured) ===")
    hdr = (f"  {'cheb':<14}{'n':>5}{'med peak':>10}{'>3s real':>10}{'>3s ctrl':>10}"
           f"{'excess':>9}{'SCOS unconf':>13}")
    print(hdr)
    rows = []
    for st, g in ok.groupby("stratum", observed=True):
        n = len(g)
        r = g["has_sig"].mean() * 100
        c = g["ctrl_sig"].mean() * 100
        u = g["scos_unconf"].mean() * 100 if "scos_unconf" in g else np.nan
        print(f"  {str(st):<14}{n:>5}{g['peak'].median():>10.2f}{r:>9.1f}%{c:>9.1f}%"
              f"{r-c:>8.1f}%{u:>12.1f}%")
        rows.append(dict(stratum=str(st), n=n, median_peak=float(g["peak"].median()),
                         real_pct=r, ctrl_pct=c, excess_pct=r - c, scos_unconf_pct=u))

    if "scos_unconf" in ok:
        print("\n=== IS SCOS AN ACCURATE ARBITER AT EVERY DEPTH? ===")
        print(f"  {'cheb':<14}{'SCOS conf: med peak':>22}{'SCOS unconf: med peak':>24}{'n conf':>9}")
        for st, g in ok.groupby("stratum", observed=True):
            a = g[~g["scos_unconf"]]["peak"]
            b = g[g["scos_unconf"]]["peak"]
            print(f"  {str(st):<14}{(a.median() if len(a) else float('nan')):>22.2f}"
                  f"{(b.median() if len(b) else float('nan')):>24.2f}{len(a):>9}")

    (out / "summary.json").write_text(json.dumps(
        {"seed": args.seed, "per_stratum": args.per_stratum, "sigma": SIGMA,
         "null_shift_arcmin": args.null_shift_arcmin, "strata": rows}, indent=2))
    print(f"\nledgers -> {out}/sample.csv, summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
