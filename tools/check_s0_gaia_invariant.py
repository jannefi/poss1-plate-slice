#!/usr/bin/env python3
"""Standing invariant: S0 must contain almost no Gaia stars.

The veto chain removes any detection with a Gaia source within 5". So the
fraction of S0 rows that HAVE a Gaia counterpart within 5" must be ~0. If it
is not, the veto did not actually run over that sky -- whatever the logs say.

This is the check that would have caught the 2026-08 partial-cone bug the day
it was introduced. In run full642-20260812 the healthy tiles sat at 0.0% while
277 pathological tiles sat at a median 98.3%, holding 56% of the entire
catalogue. Nothing else in the pipeline noticed: the veto stage ran, wrote its
xmatch file, and reported ok.

Reports per tile and overall, and exits non-zero if the overall rate or any
single tile exceeds the threshold, so it can gate a release.

Run:
    VASCO_GAIA_CACHE=/path/to/gaia python3 tools/check_s0_gaia_invariant.py \\
        --s0-csv <run>/stage_S0.csv --out work/s0_gaia_invariant.csv

Options:
    --max-tile-frac   fail if any tile exceeds this (default 0.20)
    --max-total-frac  fail if the catalogue overall exceeds this (default 0.02)
    --sample-per-tile cap rows checked per tile, for speed (default 300)
"""
from __future__ import annotations

import argparse
import csv
import functools
import glob
import os
import sys
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy_healpix import HEALPix

MATCH_ARCSEC = 5.0
_HP = HEALPix(nside=32, order="nested")


# Pixel cache bounded by BYTES, not by entry count.
#
# This was an lru_cache(maxsize=512), which bounds how many pixels are held but
# not how big they are -- and nside=32 Gaia pixels vary enormously: the median
# is ~1 MB on disk while the largest holds 4.29 M sources, 68.6 MB resident as
# ra+dec float64. 512 of those is 34.3 GB, i.e. more RAM than a typical box has,
# so a run over crowded sky could take the machine down. Counting bytes makes
# the ceiling something you choose rather than something the sky chooses.
_PIX_CACHE: "OrderedDict[tuple[str, int], tuple]" = OrderedDict()
_PIX_BYTES = 0
_PIX_BUDGET = 2 * 1024 ** 3  # 2 GiB; override with --cache-budget-gb


def _pixel(cache: str, p: int):
    global _PIX_BYTES
    key = (cache, p)
    hit = _PIX_CACHE.get(key)
    if hit is not None:
        _PIX_CACHE.move_to_end(key)
        return hit[0]

    files = glob.glob(f"{cache}/parquet/healpix_5={p}/*.parquet")
    if not files:
        val, nbytes = None, 0
    else:
        df = pd.concat([pd.read_parquet(f, columns=["ra", "dec"]) for f in files])
        ra, dec = df.ra.values, df.dec.values
        val, nbytes = (ra, dec), int(ra.nbytes + dec.nbytes)

    _PIX_CACHE[key] = (val, nbytes)
    _PIX_BYTES += nbytes
    while _PIX_BYTES > _PIX_BUDGET and len(_PIX_CACHE) > 1:
        _, (_, evicted) = _PIX_CACHE.popitem(last=False)
        _PIX_BYTES -= evicted
    return val


def _gaia_near(cache: str, ra: np.ndarray, dec: np.ndarray, pad_deg: float = 0.6):
    """Gaia rows around a set of positions, via the pixels they touch."""
    from vasco.local_cache_query import _cone_pixels
    pix = set()
    step = max(1, len(ra) // 40)
    for r, d in zip(ra[::step], dec[::step]):
        pix.update(_cone_pixels(_HP, float(r), float(d), pad_deg * 60.0))
    R, D = [], []
    for p in pix:
        g = _pixel(cache, int(p))
        if g is not None:
            R.append(g[0])
            D.append(g[1])
    if not R:
        return None
    return SkyCoord(np.concatenate(R) * u.deg, np.concatenate(D) * u.deg)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--s0-csv", required=True)
    ap.add_argument("--out", default=None, help="per-tile CSV ledger")
    ap.add_argument("--max-tile-frac", type=float, default=0.20)
    ap.add_argument("--max-total-frac", type=float, default=0.02)
    ap.add_argument("--sample-per-tile", type=int, default=300)
    ap.add_argument("--cache-budget-gb", type=float, default=2.0,
                    help="RAM ceiling for the Gaia pixel cache. Pixel sizes vary "
                         "~70x, so this is bounded by bytes rather than by count.")
    args = ap.parse_args()

    global _PIX_BUDGET
    _PIX_BUDGET = int(args.cache_budget_gb * 1024 ** 3)

    cache = os.getenv("VASCO_GAIA_CACHE")
    if not cache:
        print("[FATAL] VASCO_GAIA_CACHE is unset -- this check needs the Gaia mirror")
        return 2

    tiles = defaultdict(list)
    for row in csv.DictReader(open(args.s0_csv)):
        tiles[row["tile_id"]].append((float(row["ra"]), float(row["dec"])))
    print(f"[IN] {args.s0_csv}: {sum(len(v) for v in tiles.values()):,} rows "
          f"in {len(tiles):,} tiles")

    rng = np.random.RandomState(0)
    ledger, n_checked, n_matched = [], 0, 0
    for i, (tid, pts) in enumerate(sorted(tiles.items())):
        arr = np.array(pts)
        if len(arr) > args.sample_per_tile:
            arr = arr[rng.choice(len(arr), args.sample_per_tile, replace=False)]
        gc = _gaia_near(cache, arr[:, 0], arr[:, 1])
        if gc is None:
            frac = float("nan")
        else:
            sc = SkyCoord(arr[:, 0] * u.deg, arr[:, 1] * u.deg)
            _, d2d, _ = sc.match_to_catalog_sky(gc)
            hit = int((d2d.arcsec < MATCH_ARCSEC).sum())
            frac = hit / len(arr)
            n_checked += len(arr)
            n_matched += hit
        ledger.append({"tile_id": tid, "s0_rows": len(pts),
                       "checked": len(arr), "gaia_frac": frac})
        if i % 2000 == 0:
            print(f"  ...{i}/{len(tiles)}", flush=True)

    df = pd.DataFrame(ledger)
    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        df.to_csv(args.out, index=False)
        print(f"[OUT] {args.out}")

    total = n_matched / max(n_checked, 1)
    bad = df[df.gaia_frac > args.max_tile_frac].sort_values("s0_rows", ascending=False)
    print(f"\n[RESULT] overall S0 rows with a Gaia source within {MATCH_ARCSEC:.0f}\": "
          f"{100*total:.2f}%  (threshold {100*args.max_total_frac:.2f}%)")
    print(f"[RESULT] tiles over {100*args.max_tile_frac:.0f}%: {len(bad)} "
          f"holding {int(bad.s0_rows.sum()):,} S0 rows")
    if len(bad):
        print("\n  worst offenders:")
        for _, r in bad.head(15).iterrows():
            print(f"    {r.tile_id:34s} s0={int(r.s0_rows):6d}  gaia={100*r.gaia_frac:5.1f}%")

    failed = (total > args.max_total_frac) or len(bad) > 0
    print("\nRESULT:", "FAIL -- the veto did not cover this sky" if failed else "PASS")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    sys.exit(main())
