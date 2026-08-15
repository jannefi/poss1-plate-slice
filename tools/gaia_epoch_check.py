#!/usr/bin/env python3
"""Are the reference-catalogue rows our Gaia veto removes high-proper-motion
stars -- objects that moved, rather than objects that vanished?

This pipeline vetoes against Gaia positions PROPAGATED TO THE PLATE EPOCH
(~1950s), not against Gaia's own catalogue positions (epoch 2016). Over ~65
years most stars barely move, but a tail crosses the 5" veto radius. A
catalogue built by matching at Gaia's epoch would keep exactly those stars,
while this pipeline removes them.

That yields a specific, falsifiable prediction for the reference rows our Gaia
veto removes:

  IF proper motion explains them, they should be strongly enriched in stars
  within 5" at the PLATE epoch but beyond 5" at the CATALOGUE epoch -- far
  above the chance rate measured by the null control.

  IF NOT -- a Gaia source within 5" at BOTH epochs -- then the reference
  catalogue's own Gaia veto did not remove rows it should have, which is a
  statement about how that catalogue was built, not about epochs.

A row confirmed by this test is a star sitting at the reference position at the
plate epoch: something that moved away, not something that disappeared. That is
exactly the contaminant class a vanishing-source search must remove, and it is
checkable by anyone -- the positions are published and Gaia DR3 is public.

NEEDS NO CATALOGUE QUERY. Both Gaia position sets are already inside every tile
directory: `gaia_neighbourhood.csv` (Gaia's own epoch-2016 positions) and
`gaia_neighbourhood_at_plate.csv` (proper-motion-propagated, the one the veto
reads).

CONTROL ARM: --null-shift-arcmin displaces every reference position in RA
before matching, giving the chance-coincidence rate at the same local Gaia
density. Always read the real rate against the null; Gaia is dense enough that
a bare percentage means little on its own.

Reads the per-row output of tools/funnel_attribution.py.

Usage:
    python3 tools/gaia_epoch_check.py \\
      --funnel-csv funnel_rows.csv \\
      --tiles-root <run>/tiles \\
      --out gaia_epoch.csv

How to validate: rerun with a different --null-shift-arcmin. The real rates
must not move; the null rates should stay in the same range.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

CAT_EPOCH = "gaia_neighbourhood.csv"             # Gaia's own positions, ep 2016
PLATE_EPOCH = "gaia_neighbourhood_at_plate.csv"  # PM-propagated, what we veto on

REPORT_STAGES = ("VETOED_GAIA", "VETOED_PS1", "VETOED_USNOB",
                 "MNRAS_MORPHOLOGY", "SURVIVED_S0")


def _proj(ra, dec, ra0, dec0):
    """Tangent-plane arcsec offsets, RA wrapped."""
    cosd = np.cos(np.deg2rad(dec0))
    dra = ((np.asarray(ra, float) - ra0 + 180.0) % 360.0 - 180.0) * cosd
    return np.column_stack([dra * 3600.0, (np.asarray(dec, float) - dec0) * 3600.0])


def _load(path: Path):
    try:
        if path.stat().st_size == 0:
            return None
        d = pd.read_csv(path, usecols=["ra", "dec", "pmRA", "pmDE", "Gmag"])
    except (FileNotFoundError, ValueError, pd.errors.EmptyDataError):
        return None
    d = d.dropna(subset=["ra", "dec"])
    return d if len(d) else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--funnel-csv", required=True,
                    help="per-row output of tools/funnel_attribution.py")
    ap.add_argument("--tiles-root", required=True)
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--null-shift-arcmin", type=float, default=6.0,
                    help="RA displacement for the chance-rate control arm")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    f = pd.read_csv(args.funnel_csv)
    idx_col = "ref_index" if "ref_index" in f.columns else f.columns[0]
    f = f[f.stage != "NEVER_DETECTED"].copy()
    root = Path(args.tiles_root)

    rows = []
    for tid, g in f.groupby("tile_id"):
        ra0 = float(tid.split("_RA")[1].split("_")[0])
        ds = tid.split("_DEC")[1]
        dec0 = float(ds[1:]) * (1.0 if ds[0] == "p" else -1.0)
        cat = _load(root / tid / "catalogs" / CAT_EPOCH)
        plate = _load(root / tid / "catalogs" / PLATE_EPOCH)
        if cat is None or plate is None:
            continue
        t_cat = cKDTree(_proj(cat.ra.values, cat.dec.values, ra0, dec0))
        t_pl = cKDTree(_proj(plate.ra.values, plate.dec.values, ra0, dec0))

        shift = args.null_shift_arcmin / 60.0 / np.cos(np.deg2rad(dec0))
        for arm, ra_q in (("real", g.ra.values), ("null", g.ra.values + shift)):
            q = _proj(ra_q, g.dec.values, ra0, dec0)
            d_cat, _ = t_cat.query(q)
            d_pl, i_pl = t_pl.query(q)
            for k, (_, r) in enumerate(g.iterrows()):
                j = i_pl[k]
                pm = (float(np.hypot(plate.pmRA.values[j], plate.pmDE.values[j]))
                      if np.isfinite(plate.pmRA.values[j]) else np.nan)
                rows.append({"arm": arm, "ref_index": int(r[idx_col]),
                             "stage": r.stage, "tile_id": tid,
                             "d_cat_epoch": float(d_cat[k]),
                             "d_plate_epoch": float(d_pl[k]),
                             "pm_total_mas_yr": pm,
                             "gmag": float(plate.Gmag.values[j])})

    d = pd.DataFrame(rows)
    if d.empty:
        print("[FATAL] no rows -- do these tiles carry gaia_neighbourhood CSVs?")
        return 1
    R = args.radius_arcsec

    def block(sub, label):
        n = len(sub)
        if n == 0:
            return
        w_pl = sub.d_plate_epoch <= R
        w_cat = sub.d_cat_epoch <= R
        only_pl = w_pl & ~w_cat
        print(f"  {label:34s} n={n:5d}  within{R:.0f}\": "
              f"plate-epoch {100*w_pl.mean():6.2f}%  "
              f"catalogue-epoch {100*w_cat.mean():6.2f}%  "
              f"plate-only {100*only_pl.mean():6.2f}%")
        if only_pl.sum():
            s = sub[only_pl]
            print(f"  {'':34s}        plate-only PM: "
                  f"median {s.pm_total_mas_yr.median():7.1f}  "
                  f"p90 {s.pm_total_mas_yr.quantile(.9):7.1f} mas/yr"
                  f"   (all stars median {sub.pm_total_mas_yr.median():.1f})")

    print(f"\n{'='*104}\nGaia within {R:.0f}\" of the reference position, "
          f"at two epochs\n{'='*104}")
    for arm in ("real", "null"):
        a = d[d.arm == arm]
        print(f"\n[{arm.upper()} arm]"
              + ("" if arm == "real"
                 else f"  (positions displaced {args.null_shift_arcmin:.0f}' in RA "
                      f"-- chance rate at the same local density)"))
        block(a, "ALL DETECTED reference rows")
        for st in REPORT_STAGES:
            block(a[a.stage == st], st)

    real = d[(d.arm == "real") & (d.stage == "VETOED_GAIA")]
    null = d[(d.arm == "null") & (d.stage == "VETOED_GAIA")]
    if len(real) and len(null):
        rc = 100 * (real.d_cat_epoch <= R).mean()
        rp = 100 * (real.d_plate_epoch <= R).mean()
        nc = 100 * (null.d_cat_epoch <= R).mean()
        print(f"\n{'-'*104}\nVERDICT for the 'VETOED_GAIA' rows  (n={len(real)})")
        print(f"  plate epoch (what our veto used) : {rp:6.2f}% have Gaia within {R:.0f}\"")
        print(f"  catalogue epoch (Gaia's own)     : {rc:6.2f}%   null control {nc:6.2f}%")
        if rc < 50 and rp > 80:
            print("  -> PROPER MOTION explains it. These are stars that had moved")
            print("     away from their catalogue positions by the plate epoch, so a")
            print("     catalogue-epoch veto would keep them and this one does not.")
        elif rc > 80:
            print("  -> PROPER MOTION does NOT explain it: a Gaia source sits within")
            print(f"     {R:.0f}\" at BOTH epochs, so the reference catalogue's own Gaia")
            print("     veto did not remove rows it should have.")
        else:
            print("  -> MIXED. Neither explanation is clean; report both numbers.")
        if len(real) < 30:
            print(f"  [SMALL n] {len(real)} rows. Row-exact and individually "
                  f"checkable, but do not quote a tight interval on the rate.")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        d.to_csv(args.out, index=False)
        print(f"\n[OUT] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
