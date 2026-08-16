#!/usr/bin/env python3
"""Prove that dropping VIGNET during LDAC->CSV conversion changes nothing else.

`sex_default.param` asks SExtractor for VIGNET(45,45) -- a 45x45 postage stamp,
2025 floats per detection. `stilts tcopy` serialises every one of them into the
CSV: 72.8 MB for 3,187 rows, of which the 35 columns anything actually reads are
1.2 MB. Nothing in this repo reads VIGNET from the CSV (PSFEx reads it from the
LDAC, and single-pass runs no PSFEx), so the conversion spends ~94% of its time
on a column with no consumer.

This script runs both conversions over real tiles and requires them to agree
exactly. The bar is deliberately absolute rather than "close enough":

    same row count, same schema minus VIGNET, and max relative difference
    EXACTLY 0 on every retained numeric column.

Anything less means the change is not a pure speed change and should not ship.
A weaker tolerance would also mask the failure mode found in the alternative
implementation this replaced: astropy's CSV writer honours each column's FITS
display format and silently truncated THETA_IMAGE from 40.433807 to 40.43 on
3,185 of 3,187 rows, which "close enough" would have waved through.

Tiles are only read, never written -- output goes to a scratch directory.

How to validate
---------------
    python3 tools/verify_ldac_vignet_drop.py \
        --tiles-root <tiles_dir> --n 200 --workers 4 \
        --out work/ldac_conversion_equivalence.csv

Exit status is 0 only if every tile passed. The summary prints the worst
relative difference seen across all tiles and columns; it must be 0.0.
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from vasco.paths import get as _p

import argparse
import random
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HDU_TRIES = ["#LDAC_OBJECTS", "#2", "#1", "#0", "#3", "#4", "#5", "#6", "#7", "#8", ""]


def convert(ldac: Path, out: Path, drop_vignet: bool) -> float | None:
    """Run one conversion the way the pipeline does. Returns seconds, or None."""
    for ext in HDU_TRIES:
        in_arg = f"in={ldac}{ext}" if ext else f"in={ldac}"
        argv = (["stilts", "tpipe", in_arg, "cmd=delcols VIGNET", f"out={out}", "ofmt=csv"]
                if drop_vignet else
                ["stilts", "tcopy", in_arg, f"out={out}", "ofmt=csv"])
        t0 = time.time()
        try:
            subprocess.run(argv, check=True, capture_output=True)
        except Exception:
            continue
        if out.exists() and out.stat().st_size > 0:
            return time.time() - t0
    return None


def compare(tile: Path) -> dict:
    """Convert one tile both ways and measure disagreement."""
    r = {"tile": tile.name, "status": "", "rows_full": 0, "rows_drop": 0,
         "cols_full": 0, "cols_drop": 0, "max_rel_diff": np.nan,
         "worst_col": "", "sec_full": np.nan, "sec_drop": np.nan,
         "bytes_full": 0, "bytes_drop": 0}
    ldac = tile / "pass2.ldac"
    if not ldac.exists():
        r["status"] = "no_ldac"
        return r

    tmp = Path(tempfile.mkdtemp(prefix="ldaccmp_"))
    try:
        f_csv, d_csv = tmp / "full.csv", tmp / "drop.csv"
        r["sec_full"] = convert(ldac, f_csv, drop_vignet=False)
        if r["sec_full"] is None:
            r["status"] = "convert_failed"
            return r
        r["bytes_full"] = f_csv.stat().st_size
        a = pd.read_csv(f_csv, low_memory=False)
        r["rows_full"] = len(a)

        # Tiles cut by the two-pass config use default.param, which never asks
        # for VIGNET. `delcols VIGNET` is a hard error there, the pipeline falls
        # back to plain tcopy, and the output is identical by construction --
        # there is nothing to compare, and it is not a failure.
        if not any(c.startswith("VIGNET") for c in a.columns):
            r["status"] = "no_vignet"
            r["cols_full"] = len(a.columns)
            return r

        r["sec_drop"] = convert(ldac, d_csv, drop_vignet=True)
        if r["sec_drop"] is None:
            r["status"] = "convert_failed"
            return r
        r["bytes_drop"] = d_csv.stat().st_size
        b = pd.read_csv(d_csv, low_memory=False)
        r["rows_drop"] = len(b)
        kept = [c for c in a.columns if not c.startswith("VIGNET")]
        r["cols_full"], r["cols_drop"] = len(kept), len(b.columns)

        if len(a) != len(b):
            r["status"] = "row_count_mismatch"
            return r
        if kept != list(b.columns):
            r["status"] = "schema_mismatch"
            return r

        worst, wcol = 0.0, ""
        for col in kept:
            x, y = a[col].to_numpy(), b[col].to_numpy()
            if not (np.issubdtype(x.dtype, np.number) and np.issubdtype(y.dtype, np.number)):
                if not (a[col].astype(str) == b[col].astype(str)).all():
                    worst, wcol = np.inf, col
                continue
            # NaNs must sit in the same places, then compare the finite values.
            if not np.array_equal(np.isnan(x), np.isnan(y)):
                worst, wcol = np.inf, col
                break
            m = ~np.isnan(x)
            if not m.any():
                continue
            rel = np.abs(x[m] - y[m]) / np.maximum(np.abs(x[m]), 1e-30)
            if rel.size and rel.max() > worst:
                worst, wcol = float(rel.max()), col
        r["max_rel_diff"], r["worst_col"] = worst, wcol
        r["status"] = "ok" if worst == 0.0 else "DIFFERS"
        return r
    except Exception as e:                                   # noqa: BLE001
        r["status"] = f"error:{type(e).__name__}"
        return r
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tiles-root", default=str(_p("tiles_dir")))
    ap.add_argument("--n", type=int, default=200, help="tiles to sample [200]")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20260808)
    ap.add_argument("--out", default="work/ldac_conversion_equivalence.csv")
    args = ap.parse_args()

    root = Path(args.tiles_root)
    print(f"[SCAN] {root}", flush=True)
    tiles = [p for p in root.iterdir() if p.is_dir() and (p / "pass2.ldac").exists()]
    print(f"[SCAN] {len(tiles):,} tiles carry a pass2.ldac")
    random.Random(args.seed).shuffle(tiles)
    tiles = tiles[:args.n]
    print(f"[PLAN] comparing {len(tiles)} tiles, {args.workers} workers\n", flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(compare, t): t for t in tiles}
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            rows.append(r)
            if r["status"] not in ("ok", "no_vignet") or i % 25 == 0:
                print(f"  [{i}/{len(tiles)}] {r['tile']:<34} {r['status']:<18} "
                      f"rows={r['rows_full']:>6} rel={r['max_rel_diff']}", flush=True)

    d = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(args.out, index=False)

    ok = d[d.status == "ok"]
    print(f"\n=== {len(d)} tiles: {d.status.value_counts().to_dict()} ===")
    print(f"  ({int((d.status == 'no_vignet').sum())} carry no VIGNET -- two-pass config, "
          f"pipeline falls back to tcopy, nothing to compare)")
    if len(ok):
        print(f"  rows compared      {int(ok.rows_full.sum()):,}")
        print(f"  worst rel diff     {ok.max_rel_diff.max():.3e}   (must be 0.000e+00)")
        print(f"  seconds  full      {ok.sec_full.sum():8.1f}   median {ok.sec_full.median():.2f}")
        print(f"  seconds  no VIGNET {ok.sec_drop.sum():8.1f}   median {ok.sec_drop.median():.2f}")
        print(f"  speedup            {ok.sec_full.sum()/max(ok.sec_drop.sum(), 1e-9):.1f}x")
        print(f"  size     full      {ok.bytes_full.sum()/1e9:8.2f} GB")
        print(f"  size     no VIGNET {ok.bytes_drop.sum()/1e9:8.2f} GB "
              f"({100*(1-ok.bytes_drop.sum()/max(ok.bytes_full.sum(),1)):.1f}% smaller)")
    bad = d[~d.status.isin(["ok", "no_vignet"])]
    print(f"\nwrote {args.out}")
    if len(bad):
        print(f"\n[FAIL] {len(bad)} tiles did not match exactly:")
        print(bad[["tile", "status", "max_rel_diff", "worst_col"]].to_string(index=False))
        return 1
    print("\n[PASS] every tile identical on every retained column.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
