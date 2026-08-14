#!/usr/bin/env python3
"""How far does the WCSFIX refit actually move the released rows?

The release README quotes this displacement to tell readers whether the raw-WCS
and refit variants are interchangeable for positional work. It must therefore be
a measurement, not an estimate: two ~200-row samples of this quantity disagreed by
4x on the p90 and by a factor of 4 on the fraction beyond 2", which is why this
tool reads the whole survey instead of sampling it.

Each per-tile survivor catalogue retains both coordinate systems side by side --
ALPHAWIN_J2000/DELTAWIN_J2000 (the raw plate solution WCSFIX was fitted from) and
RA_corr/Dec_corr (the refit result that reached the catalogue). The displacement is
the angular separation between them, per row.

Reports a bootstrap CI over TILES, not rows: rows within a tile share one refit,
so they are not independent and a row-level CI would be far too narrow.

Also reports the known astrometry-defect ("mode-2") tiles separately, via
--exclude-tiles. They carry ~950 rows each against a survey average near 4.5, so
any sample that happens to include one is dominated by it -- one earlier attempt
here drew 17,851 rows from 319 tiles for that reason and reported a 79" maximum
that belonged entirely to the defect. Those tiles already carry a
"do not use for positional work" warning in the release, so the headline number
must exclude them and state that it does.

Usage:
    python3 tools/measure_wcsfix_shift.py \
      --filtered-root work/slice/filtered --out work/wcsfix_shift.csv
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

COLS = ["ALPHAWIN_J2000", "DELTAWIN_J2000", "RA_corr", "Dec_corr"]


def do_tile(job):
    tile_id, root = job
    p = Path(root) / tile_id / "catalogs" / "sextractor_pass2.filtered.csv"
    try:
        if p.stat().st_size == 0:
            return tile_id, None
        d = pd.read_csv(p, usecols=COLS).dropna()
    except Exception:
        return tile_id, None
    if d.empty:
        return tile_id, None
    # Wrap the RA difference. A plain subtraction reads ~360 deg for any row
    # within a few arcsec of the meridian, and the defect tiles happen to cluster
    # near RA 0, which produced a spurious 270 deg maximum on the first run.
    dra = ((d.RA_corr.values - d.ALPHAWIN_J2000.values + 180.0) % 360.0) - 180.0
    sep = np.hypot(dra * np.cos(np.deg2rad(d.DELTAWIN_J2000.values)),
                   d.Dec_corr.values - d.DELTAWIN_J2000.values) * 3600.0
    return tile_id, sep.astype(np.float32)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--filtered-root", required=True,
                    help="Directory of <tile_id>/catalogs/sextractor_pass2.filtered.csv")
    ap.add_argument("--tile-list", default=None,
                    help="Optional CSV with a tile_id column. Avoids listing a "
                         "31k-entry directory on a spinning disk, which is the "
                         "slowest part of this measurement.")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--exclude-tiles", default=None,
                    help="CSV with a tile_id column -- the known astrometry-defect "
                         "tiles. Reported separately rather than dropped silently.")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--out", default=None, help="Per-tile summary CSV.")
    args = ap.parse_args()

    root = Path(args.filtered_root)
    if args.tile_list:
        tiles = pd.read_csv(args.tile_list).tile_id.astype(str).tolist()
    else:
        tiles = sorted(p.name for p in root.iterdir() if p.is_dir())
    print(f"[CONFIG] {len(tiles)} tiles under {root}", flush=True)

    per_tile, all_sep, done = [], [], 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for f in as_completed([ex.submit(do_tile, (t, str(root))) for t in tiles]):
            tid, sep = f.result()
            done += 1
            if sep is not None:
                all_sep.append(sep)
                per_tile.append({"tile_id": tid, "rows": len(sep),
                                 "median_as": float(np.median(sep)),
                                 "max_as": float(sep.max())})
            if done % 2000 == 0:
                print(f"  {done}/{len(tiles)}  tiles with rows={len(per_tile)}",
                      flush=True)

    if not all_sep:
        print("[FATAL] no rows measured")
        return 1
    df = pd.DataFrame(per_tile)

    excl = set()
    if args.exclude_tiles:
        e = pd.read_csv(args.exclude_tiles)
        col = "tile_id" if "tile_id" in e.columns else e.columns[0]
        excl = set(e[col].astype(str))
    keep_i = [i for i, r in enumerate(per_tile) if r["tile_id"] not in excl]
    drop_i = [i for i, r in enumerate(per_tile) if r["tile_id"] in excl]
    if drop_i:
        d = np.concatenate([all_sep[i] for i in drop_i])
        print(f"\n[DEFECT TILES] {len(drop_i)} tiles, {len(d)} rows "
              f"({100.0*len(d)/sum(len(x) for x in all_sep):.1f}% of all rows) "
              f"-- median {np.median(d):.2f}\", max {d.max():.2f}\"")
        print("  Excluded from the headline below. These already carry a "
              "do-not-use-positionally warning.")
        all_sep = [all_sep[i] for i in keep_i]
        df = df[~df.tile_id.isin(excl)]
    a = np.concatenate(all_sep)

    print(f"\n{'='*72}\nWCSFIX displacement of released rows\n{'='*72}")
    print(f"  {len(a)} rows over {len(df)} tiles with survivors "
          f"({len(tiles)} tiles examined)")
    print(f"\n  {'quantile':<10} {'arcsec':>9}")
    for q in (25, 50, 75, 90, 95, 99):
        print(f"  p{q:<9} {np.percentile(a, q):8.3f}")
    print(f"  {'max':<10} {a.max():8.3f}")
    print(f"\n  {'threshold':<12} {'frac beyond':>12}")
    for t in (0.25, 1.0, 2.0, 5.0):
        print(f"  {t:>6.2f}\"      {100*(a > t).mean():10.2f}%")

    # Bootstrap over tiles -- rows in a tile share one refit.
    rng = np.random.RandomState(0)
    idx = np.arange(len(all_sep))
    meds, f2 = [], []
    for _ in range(args.bootstrap):
        s = rng.choice(idx, len(idx), replace=True)
        c = np.concatenate([all_sep[i] for i in s])
        meds.append(np.median(c))
        f2.append(100 * (c > 2.0).mean())
    print(f"\n  bootstrap over {len(all_sep)} tiles, {args.bootstrap} resamples:")
    print(f"    median      {np.mean(meds):6.3f}\"  95% CI "
          f"{np.percentile(meds,2.5):.3f}-{np.percentile(meds,97.5):.3f}\"")
    print(f"    frac > 2\"   {np.mean(f2):6.2f}%  95% CI "
          f"{np.percentile(f2,2.5):.2f}-{np.percentile(f2,97.5):.2f}%")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, index=False)
        print(f"\n[OUT] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
