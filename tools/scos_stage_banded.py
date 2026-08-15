#!/usr/bin/env python3
"""SuperCOSMOS cross-scan stage, resolved by band rather than by whole table.

WHY THIS EXISTS INSTEAD OF scripts/stage_supercosmos_post.py

The shipped stage matches against `supercosmos.sources` as a whole. That table
is a MERGE over four photographic bands -- B, R1, R2 and I -- and I is the
SERC/POSS-II photographic I band (IV-N emulsion, roughly 715-900 nm), which
reaches into the near-infrared. Measured on a 1 deg^2 patch at dec +29: of
21,078 sources, 4,957 have no R1 detection and 573 (2.7%) are I-only. So a
whole-table match can be satisfied by an I-band-only source -- a near-infrared
detection reaching the result by way of a merge, which is not what a
cross-scan check of red POSS-I plates should rest on.

There is a resolution that is better on both axes at once, so this tool takes
it. It asks for three match counts in ONE query per chunk:

    nmatch_any   any band            -- what the shipped stage would have done
    nmatch_noI   B, R1 or R2         -- optical only, I excluded
    nmatch_r1    R1 only             -- POSS-I E: the same plate material

R1 is the scientifically correct arm for the stated purpose. MNRAS 2022 uses
SuperCOSMOS as a CROSS-SCAN consistency check: an independent digitization of
the same photographic material should see a real source, while a scan artifact
(a dust speck, a plate flaw on one scan) will not be reproduced. R1 IS POSS-I
E -- our plates, digitized by someone else. R2 is POSS-II, a different epoch,
so matching it would test persistence rather than scan consistency, which is a
different question and one this project measures separately.

Costing all three in one query means the I-band contribution is measured rather
than argued about, at no extra service load.

RESUMABLE. Each chunk writes its own TAP output; re-running skips completed
chunks. Kind to a service known to cancel large jobs.

Usage:
    python3 tools/scos_stage_banded.py \\
      --in-csv <run-dir>/stages/stage_S0_upload.csv \\
      --out-dir <run-dir>/stages/scos \\
      --chunk-size 5000

How to validate: --limit-chunks 1 first and read the printed rates; the three
counts must satisfy nmatch_r1 <= nmatch_noI <= nmatch_any for every row (the
tool asserts it).
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd

TAPURL = "https://dc.g-vo.org/__system__/tap/run"
TABLE = "supercosmos.sources"


def adql(radius_deg: float, table: str) -> str:
    return (
        "SELECT u.row_id AS row_id, "
        "COUNT(*) AS nmatch_any, "
        "SUM(CASE WHEN s.objidb IS NOT NULL OR s.objidr1 IS NOT NULL "
        "OR s.objidr2 IS NOT NULL THEN 1 ELSE 0 END) AS nmatch_noi, "
        "SUM(CASE WHEN s.objidr1 IS NOT NULL THEN 1 ELSE 0 END) AS nmatch_r1 "
        "FROM TAP_UPLOAD.t1 AS u "
        f"JOIN {table} AS s "
        "ON 1 = CONTAINS( POINT('ICRS', s.raj2000, s.dej2000), "
        f"CIRCLE('ICRS', u.ra, u.dec, {radius_deg:.10f}) ) "
        "GROUP BY u.row_id"
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-csv", required=True, help="src_id, ra, dec upload view")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--chunk-size", type=int, default=5000)
    ap.add_argument("--limit-chunks", type=int, default=0,
                    help="stop after N chunks (0 = all); use 1 for a pilot")
    ap.add_argument("--tapurl", default=TAPURL)
    ap.add_argument("--table", default=TABLE)
    ap.add_argument("--sleep", type=float, default=2.0,
                    help="seconds between chunks; be a good citizen")
    args = ap.parse_args(argv)

    stilts = shutil.which("stilts")
    if not stilts:
        raise SystemExit("[FAIL] stilts not on PATH")
    out = Path(args.out_dir)
    (out / "parts").mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(args.in_csv)
    if "src_id" not in d.columns:
        raise SystemExit("[FAIL] input needs src_id, ra, dec")
    n = len(d)
    radius_deg = args.radius_arcsec / 3600.0
    q = adql(radius_deg, args.table)
    nchunks = (n + args.chunk_size - 1) // args.chunk_size
    todo = nchunks if not args.limit_chunks else min(nchunks, args.limit_chunks)
    print(f"[SCOS] {n:,} rows, {nchunks} chunks of {args.chunk_size}, "
          f"running {todo}; radius {args.radius_arcsec}\"")
    print(f"[SCOS] ADQL: {q}\n")

    t0 = time.time()
    for i in range(todo):
        part = out / "parts" / f"part_{i+1:04d}.csv"
        if part.exists() and part.stat().st_size:
            print(f"[SCOS] part {i+1}/{todo}: cached")
            continue
        sub = d.iloc[i * args.chunk_size:(i + 1) * args.chunk_size]
        up_csv = out / "parts" / f"_up_{i+1:04d}.csv"
        up_vot = out / "parts" / f"_up_{i+1:04d}.vot"
        # The upload column is row_id: 'NUMBER'/'number' are forbidden by the
        # CSV contract and 'src_id' is kept as the join key on our side.
        sub.rename(columns={"src_id": "row_id"})[["row_id", "ra", "dec"]] \
           .to_csv(up_csv, index=False)
        subprocess.run([stilts, "tcopy", f"in={up_csv}", "ifmt=csv",
                        f"out={up_vot}", "ofmt=votable"], check=True)
        try:
            subprocess.run([stilts, "tapquery", f"tapurl={args.tapurl}",
                            "nupload=1", f"upload1={up_vot}", "upname1=t1",
                            "ufmt1=votable", f"adql={q}",
                            f"out={part}", "ofmt=csv"], check=True)
        except subprocess.CalledProcessError as e:
            print(f"[SCOS] part {i+1} FAILED (exit {e.returncode}). "
                  f"Re-run to resume from here.")
            return 1
        up_csv.unlink(missing_ok=True)
        up_vot.unlink(missing_ok=True)
        el = time.time() - t0
        print(f"[SCOS] part {i+1}/{todo}  {len(sub)} rows  "
              f"{el:.0f}s elapsed  ~{el/(i+1)*(todo-i-1):.0f}s left", flush=True)
        time.sleep(args.sleep)

    parts = sorted((out / "parts").glob("part_*.csv"))
    got = [pd.read_csv(p) for p in parts if p.stat().st_size]
    m = (pd.concat(got, ignore_index=True) if got
         else pd.DataFrame(columns=["row_id", "nmatch_any", "nmatch_noi", "nmatch_r1"]))
    done_rows = min(todo * args.chunk_size, n)
    base = d.iloc[:done_rows][["src_id", "ra", "dec"]].copy()
    f = base.merge(m.rename(columns={"row_id": "src_id"}), on="src_id", how="left")
    for c in ("nmatch_any", "nmatch_noi", "nmatch_r1"):
        f[c] = f[c].fillna(0).astype(int)
    assert (f.nmatch_r1 <= f.nmatch_noi).all() and (f.nmatch_noi <= f.nmatch_any).all(), \
        "band counts are not nested -- the query is wrong"

    flags = out / "scos_flags.csv"
    f.to_csv(flags, index=False)
    res = {"rows_queried": int(len(f)), "radius_arcsec": args.radius_arcsec,
           "chunks_done": todo, "chunks_total": nchunks, "table": args.table}
    for name, col in (("any_band", "nmatch_any"), ("optical_noI", "nmatch_noi"),
                      ("r1_possI_E", "nmatch_r1")):
        hit = int((f[col] > 0).sum())
        res[name] = {"matched": hit, "matched_pct": round(100 * hit / len(f), 3),
                     "unmatched": int(len(f) - hit)}
    (out / "scos_ledger.json").write_text(json.dumps(res, indent=2) + "\n")

    print(f"\n[SCOS] {len(f):,} rows queried")
    print(f"{'arm':14s} {'matched':>9s} {'%':>7s}   {'UNMATCHED (would be dropped)':>30s}")
    for name in ("any_band", "optical_noI", "r1_possI_E"):
        r = res[name]
        print(f"{name:14s} {r['matched']:9,} {r['matched_pct']:6.2f}%   "
              f"{r['unmatched']:>12,}")
    print(f"\n[OUT] {flags}\n[OUT] {out/'scos_ledger.json'}")
    print("[NOTE] Nothing is dropped here. This writes flags only; the shrink "
          "is a separate, deliberate anti-join.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
