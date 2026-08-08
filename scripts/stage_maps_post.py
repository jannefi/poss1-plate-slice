#!/usr/bin/env python3
"""
EXPERIMENTAL — Post-pipeline MAPS (Minnesota Automated Plate Scanner)
reduction stage.

Status
------
Not an official veto stage. Use for exploration and candidate reduction
only. Results are not yet validated against the full pipeline. Do not
use as a hard gate without further testing.

Goal
----
Given a run directory containing a stage CSV (e.g. stage_S0.csv), cross-
match against the local MAPS POSS-I catalog mirror (vasco.maps_cache_query,
see that module and README.data_format.md in the mirror root for format
details) within a configurable radius (default 5 arcsec) and write a
shrinking-set stage output. Mirrors scripts/stage_gsc_post.py's exact
CLI/output contract so the two are directly comparable, but is fully
local (no VizieR/network dependency) since the MAPS mirror lives on disk.

Why MAPS is a good veto candidate for POSS-I red detections
-------------------------------------------------------------
MAPS is pure POSS-I material -- the same imaging generation as both our
own SExtractor detections and V's source plates -- unlike GSC (which
mixes POSS-I/POSS-II/other-epoch plates). Only objects independently
matched on *both* the O and E plate exposures are in the catalog at all
(unmatched images are scratches/defects/noise, excluded before this
stage ever sees them), so a MAPS match is a comparatively strong "this
position was seen by a real, independently-confirmed exposure" signal.

Requires VASCO_MAPS_CACHE set to the mirror root (containing
maps_plate_index.parquet and parquet_icrs_by_plate/P###.parquet).

Usage
-----
VASCO_MAPS_CACHE=<maps_cache> \\
python scripts/stage_maps_post.py \\
    --run-dir ./work/runs/run-S1-... \\
    --input-glob 'stages/stage_S0.csv' \\
    --stage S1

Outputs (under <run-dir>/stages/)
----------------------------------
1) stage_<STAGE>_MAPS.csv
   Kept remainder AFTER MAPS elimination (rows WITHOUT a MAPS match).
   Columns: src_id, ra, dec

2) stage_<STAGE>_MAPS_flags.csv
   Full audit table for ALL input rows.
   Columns: src_id, ra, dec, has_maps_match, best_sep_arcsec,
   maps_field, maps_starnum, maps_galnode, maps_flag, maps_zero_coverage

3) stage_<STAGE>_MAPS_ledger.json
   Counts + parameters used, including zero_coverage_rows (rows where no
   MAPS plate's bounding box overlapped the query at all -- a data
   availability gap, e.g. near P003/P926 or between plate footprints --
   distinct from "plate(s) overlapped, zero objects within radius").

Notes
-----
- Input must have src_id (or row_id), ra, dec columns.
- flag_ok_only (default on) restricts MAPS matches to flag % 100 == 0
  (Oflag == 0 and Eflag == 0 -- no moments/classifier error, not a
  scratch, not stripe-clipped, has background coverage, non-negative
  sky). NOT flag == 0 -- see vasco/maps_cache_query.py's docstring for
  why (Eduplicates == 1 is the normal case, not a defect).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def _read_header(path: Path) -> List[str]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.reader(f)
        return [c.strip().lstrip("﻿") for c in next(r, [])]


def _detect_cols(cols: List[str]) -> Tuple[str, str, str]:
    cset = {c.strip() for c in cols}

    if "src_id" in cset:
        src = "src_id"
    elif "row_id" in cset:
        src = "row_id"
    else:
        raise RuntimeError("Input CSV missing required id column 'src_id' (or 'row_id')")

    ra_candidates = ["ra", "RA", "RA_corr", "ALPHAWIN_J2000", "ALPHA_J2000"]
    dec_candidates = ["dec", "DEC", "Dec", "Dec_corr", "DELTAWIN_J2000", "DELTA_J2000"]

    ra = next((c for c in ra_candidates if c in cset), None)
    dec = next((c for c in dec_candidates if c in cset), None)
    if not ra or not dec:
        raise RuntimeError(
            "Input CSV missing RA/Dec columns (expected one of ra/dec, RA/DEC, RA_corr/Dec_corr, etc.)"
        )
    return src, ra, dec


def _count_data_rows(path: Path) -> int:
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            n = sum(1 for _ in f)
        return max(0, n - 1)
    except Exception:
        return 0


def _iter_input_rows(path: Path, src_col: str, ra_col: str, dec_col: str) -> Iterable[Tuple[str, str, str]]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.DictReader(f)
        for row in r:
            sid = (row.get(src_col) or "").strip()
            ra = (row.get(ra_col) or "").strip()
            dec = (row.get(dec_col) or "").strip()
            if not sid or not ra or not dec:
                continue
            yield sid, ra, dec


@dataclass
class ChunkStats:
    chunk: str
    input_rows: int
    matched_rows: int
    kept_rows: int
    zero_coverage_rows: int


def main() -> int:
    ap = argparse.ArgumentParser(
        description="[EXPERIMENTAL] Post-pipeline MAPS reduction stage."
    )
    ap.add_argument("--run-dir", required=True, help="Run folder, e.g. ./work/runs/run-S1-...")
    ap.add_argument(
        "--input-glob",
        default="stages/stage_S0.csv",
        help="Glob (relative to run-dir) for input stage CSV. Default: stages/stage_S0.csv",
    )
    ap.add_argument("--stage", default="S1", help="Stage label used in output filenames. Default: S1")
    ap.add_argument("--out-dir", default=None, help="Output directory. Default: <run-dir>/stages")
    ap.add_argument("--radius-arcsec", type=float, default=5.0, help="Match radius in arcsec. Default: 5")
    ap.add_argument("--flag-ok-only", dest="flag_ok_only", action="store_true", default=True,
                     help="Restrict MAPS matches to flag %% 100 == 0 (default: on)")
    ap.add_argument("--no-flag-ok-only", dest="flag_ok_only", action="store_false",
                     help="Match against MAPS regardless of flag value")
    args = ap.parse_args()

    if not os.getenv("VASCO_MAPS_CACHE"):
        raise SystemExit("VASCO_MAPS_CACHE must be set to the MAPS mirror root")

    from vasco.maps_cache_query import query_maps, plates_overlapping_count

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise SystemExit(f"run-dir not found: {run_dir}")

    out_dir = Path(args.out_dir) if args.out_dir else (run_dir / "stages")
    out_dir.mkdir(parents=True, exist_ok=True)

    chunks = sorted(run_dir.glob(args.input_glob))
    if not chunks:
        raise SystemExit(f"No inputs matched: {run_dir}/{args.input_glob}")

    hdr = _read_header(chunks[0])
    src_col, ra_col, dec_col = _detect_cols(hdr)

    stage = args.stage
    out_kept = out_dir / f"stage_{stage}_MAPS.csv"
    out_flags = out_dir / f"stage_{stage}_MAPS_flags.csv"
    out_ledger = out_dir / f"stage_{stage}_MAPS_ledger.json"

    flags_fields = ["src_id", "ra", "dec", "has_maps_match", "best_sep_arcsec",
                    "maps_field", "maps_starnum", "maps_galnode", "maps_flag",
                    "maps_zero_coverage", "source_chunk"]

    total_in = total_match = total_kept = total_zero_cov = 0
    per_chunk: List[ChunkStats] = []

    with out_kept.open("w", newline="", encoding="utf-8") as f_kept, \
         out_flags.open("w", newline="", encoding="utf-8") as f_flags:

        kept_w = csv.DictWriter(f_kept, fieldnames=["src_id", "ra", "dec"])
        flags_w = csv.DictWriter(f_flags, fieldnames=flags_fields)
        kept_w.writeheader()
        flags_w.writeheader()

        for ch in chunks:
            in_rows = _count_data_rows(ch)
            if in_rows == 0:
                per_chunk.append(ChunkStats(ch.name, 0, 0, 0, 0))
                continue

            mcount = kcount = zcount = 0

            for sid, ra_s, dec_s in _iter_input_rows(ch, src_col, ra_col, dec_col):
                ra, dec = float(ra_s), float(dec_s)

                n_overlap = plates_overlapping_count(ra, dec, radius_arcsec=args.radius_arcsec)
                zero_coverage = (n_overlap == 0)
                if zero_coverage:
                    zcount += 1

                df = query_maps(ra, dec, radius_arcsec=args.radius_arcsec,
                                flag_ok_only=args.flag_ok_only)
                has_match = df is not None and len(df) > 0

                if has_match:
                    best = df.iloc[0]
                    sep = f"{best['sep_arcsec']:.6f}"
                    maps_field = str(int(best["POSS_field"]))
                    maps_starnum = str(int(best["starnumO"]))
                    maps_galnode = str(int(best["galnodO_x1000"]))
                    maps_flag = str(int(best["flag"]))
                    mcount += 1
                else:
                    sep = maps_field = maps_starnum = maps_galnode = maps_flag = ""

                flags_w.writerow({
                    "src_id": sid,
                    "ra": ra_s,
                    "dec": dec_s,
                    "has_maps_match": 1 if has_match else 0,
                    "best_sep_arcsec": sep,
                    "maps_field": maps_field,
                    "maps_starnum": maps_starnum,
                    "maps_galnode": maps_galnode,
                    "maps_flag": maps_flag,
                    "maps_zero_coverage": 1 if zero_coverage else 0,
                    "source_chunk": ch.name,
                })

                if not has_match:
                    kept_w.writerow({"src_id": sid, "ra": ra_s, "dec": dec_s})
                    kcount += 1

            per_chunk.append(ChunkStats(ch.name, in_rows, mcount, kcount, zcount))
            total_in += in_rows
            total_match += mcount
            total_kept += kcount
            total_zero_cov += zcount

    ledger = {
        "experimental": True,
        "run_dir": str(run_dir),
        "input_glob": args.input_glob,
        "stage": stage,
        "radius_arcsec": float(args.radius_arcsec),
        "flag_ok_only": bool(args.flag_ok_only),
        "maps_cache": os.getenv("VASCO_MAPS_CACHE"),
        "input_chunks": [p.name for p in chunks],
        "totals": {
            "input_rows": total_in,
            "matched_rows": total_match,
            "kept_rows": total_kept,
            "zero_coverage_rows": total_zero_cov,
        },
        "per_chunk": [cs.__dict__ for cs in per_chunk],
        "outputs": {"kept_csv": str(out_kept), "flags_csv": str(out_flags), "ledger_json": str(out_ledger)},
        "columns_detected": {"src_id_col": src_col, "ra_col": ra_col, "dec_col": dec_col},
    }
    out_ledger.write_text(json.dumps(ledger, indent=2), encoding="utf-8")

    print(f"[MAPS] [EXPERIMENTAL] chunks={len(chunks)} input_rows={total_in} "
          f"matched={total_match} kept={total_kept} zero_coverage={total_zero_cov}")
    print(f"[MAPS] wrote: {out_kept}")
    print(f"[MAPS] wrote: {out_flags}")
    print(f"[MAPS] wrote: {out_ledger}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
