#!/usr/bin/env python3
"""
EXPERIMENTAL — Post-pipeline plate-edge exclusion stage (demonstration).

Status
------
Not an official veto stage, and not a proposal to revert naive
tessellation (a deliberate, locked choice for this project -- see
context/02_DECISIONS.md). This reproduces old VASCO60's plate-edge
exclusion policy against our *existing* naive-tessellation output,
purely as a standalone comparison/demonstration, same posture as
scripts/stage_gsc_post.py and scripts/stage_maps_post.py.

Goal
----
Given a run directory containing a stage CSV with tile_id per row (e.g.
stage_S0.csv), join each row's tile_id against
data/metadata/tile_plate_edge_report.csv (see
scripts/compute_tile_plate_edge_report.py) and keep only "core" rows --
i.e. reproduce old VASCO60's is_edge_core() rule:
    is_core = (class_px == 'core') OR (class_arcsec == 'core')

Unlike GSC/MAPS, this cut is entirely local and tile-granular (every
candidate from the same tile gets the same verdict) -- no per-row
positional cross-match needed, since each candidate already carries its
own originating tile_id.

Usage
-----
python scripts/stage_edge_post.py \\
    --run-dir ./work/runs/run-S1-... \\
    --input-glob 'stages/stage_S0.csv' \\
    --edge-report-csv data/metadata/tile_plate_edge_report.csv \\
    --stage S1

Outputs (under <run-dir>/stages/)
----------------------------------
1) stage_<STAGE>_EDGE.csv
   Kept remainder (rows whose tile is edge-core).
   Columns: src_id, ra, dec

2) stage_<STAGE>_EDGE_flags.csv
   Full audit table for ALL input rows.
   Columns: src_id, ra, dec, tile_id, is_core, class_px, class_arcsec,
   min_edge_dist_px, min_edge_dist_arcsec, plate_id, edge_report_missing

3) stage_<STAGE>_EDGE_ledger.json
   Counts + parameters used, including rows whose tile_id has no entry
   at all in the edge report (edge_report_missing) -- distinct from a
   real "not core" verdict.

Notes
-----
- Input must have src_id (or row_id), tile_id, ra, dec columns.
- Rows whose tile_id isn't in the edge report are KEPT by default
  (--drop-missing to instead drop them) since a missing entry means "we
  don't know" (e.g. tile not covered by the edge-report run), not "this
  is an edge tile."
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def _read_header(path: Path) -> List[str]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.reader(f)
        return [c.strip().lstrip("﻿") for c in next(r, [])]


def _detect_cols(cols: List[str]) -> Tuple[str, str, str, str]:
    cset = {c.strip() for c in cols}

    if "src_id" in cset:
        src = "src_id"
    elif "row_id" in cset:
        src = "row_id"
    else:
        raise RuntimeError("Input CSV missing required id column 'src_id' (or 'row_id')")

    if "tile_id" not in cset:
        raise RuntimeError("Input CSV missing required 'tile_id' column")

    ra_candidates = ["ra", "RA", "RA_corr", "ALPHAWIN_J2000", "ALPHA_J2000"]
    dec_candidates = ["dec", "DEC", "Dec", "Dec_corr", "DELTAWIN_J2000", "DELTA_J2000"]
    ra = next((c for c in ra_candidates if c in cset), None)
    dec = next((c for c in dec_candidates if c in cset), None)
    if not ra or not dec:
        raise RuntimeError("Input CSV missing RA/Dec columns")
    return src, "tile_id", ra, dec


def _count_data_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        n = sum(1 for _ in f)
    return max(0, n - 1)


def load_edge_report(csv_path: Path) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        for row in csv.DictReader(f):
            tid = (row.get("tile_id") or "").strip()
            if not tid:
                continue
            out[tid] = row
    return out


def is_edge_core(edge_rec: dict) -> bool:
    return (edge_rec.get("class_px") == "core") or (edge_rec.get("class_arcsec") == "core")


def _iter_input_rows(path: Path, src_col: str, tile_col: str, ra_col: str, dec_col: str) -> Iterable[Tuple[str, str, str, str]]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        for row in csv.DictReader(f):
            sid = (row.get(src_col) or "").strip()
            tid = (row.get(tile_col) or "").strip()
            ra = (row.get(ra_col) or "").strip()
            dec = (row.get(dec_col) or "").strip()
            if not sid or not tid or not ra or not dec:
                continue
            yield sid, tid, ra, dec


@dataclass
class ChunkStats:
    chunk: str
    input_rows: int
    core_rows: int
    kept_rows: int
    missing_rows: int


def main() -> int:
    ap = argparse.ArgumentParser(description="[EXPERIMENTAL] Plate-edge exclusion demonstration stage.")
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--input-glob", default="stages/stage_S0.csv")
    ap.add_argument("--edge-report-csv", default="data/metadata/tile_plate_edge_report.csv")
    ap.add_argument("--stage", default="S1")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--drop-missing", action="store_true",
                    help="Drop rows whose tile_id has no edge-report entry (default: keep them).")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise SystemExit(f"run-dir not found: {run_dir}")

    edge_map = load_edge_report(Path(args.edge_report_csv))
    if not edge_map:
        raise SystemExit(f"edge report empty or not found: {args.edge_report_csv}")

    out_dir = Path(args.out_dir) if args.out_dir else (run_dir / "stages")
    out_dir.mkdir(parents=True, exist_ok=True)

    chunks = sorted(run_dir.glob(args.input_glob))
    if not chunks:
        raise SystemExit(f"No inputs matched: {run_dir}/{args.input_glob}")

    hdr = _read_header(chunks[0])
    src_col, tile_col, ra_col, dec_col = _detect_cols(hdr)

    stage = args.stage
    out_kept = out_dir / f"stage_{stage}_EDGE.csv"
    out_flags = out_dir / f"stage_{stage}_EDGE_flags.csv"
    out_ledger = out_dir / f"stage_{stage}_EDGE_ledger.json"

    flags_fields = ["src_id", "ra", "dec", "tile_id", "is_core", "class_px", "class_arcsec",
                    "min_edge_dist_px", "min_edge_dist_arcsec", "plate_id",
                    "edge_report_missing", "source_chunk"]

    total_in = total_core = total_kept = total_missing = 0
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

            ccount = kcount = mcount = 0

            for sid, tid, ra, dec in _iter_input_rows(ch, src_col, tile_col, ra_col, dec_col):
                edge_rec = edge_map.get(tid)
                missing = edge_rec is None
                if missing:
                    mcount += 1
                    core = not args.drop_missing
                    class_px = class_arcsec = min_px = min_as = plate_id = ""
                else:
                    core = is_edge_core(edge_rec)
                    class_px = edge_rec.get("class_px", "")
                    class_arcsec = edge_rec.get("class_arcsec", "")
                    min_px = edge_rec.get("min_edge_dist_px", "")
                    min_as = edge_rec.get("min_edge_dist_arcsec", "")
                    plate_id = edge_rec.get("plate_id", "")

                if core:
                    ccount += 1

                flags_w.writerow({
                    "src_id": sid, "ra": ra, "dec": dec, "tile_id": tid,
                    "is_core": 1 if core else 0,
                    "class_px": class_px, "class_arcsec": class_arcsec,
                    "min_edge_dist_px": min_px, "min_edge_dist_arcsec": min_as,
                    "plate_id": plate_id,
                    "edge_report_missing": 1 if missing else 0,
                    "source_chunk": ch.name,
                })

                if core:
                    kept_w.writerow({"src_id": sid, "ra": ra, "dec": dec})
                    kcount += 1

            per_chunk.append(ChunkStats(ch.name, in_rows, ccount, kcount, mcount))
            total_in += in_rows
            total_core += ccount
            total_kept += kcount
            total_missing += mcount

    ledger = {
        "experimental": True,
        "run_dir": str(run_dir),
        "input_glob": args.input_glob,
        "stage": stage,
        "edge_report_csv": args.edge_report_csv,
        "drop_missing": bool(args.drop_missing),
        "input_chunks": [p.name for p in chunks],
        "totals": {
            "input_rows": total_in,
            "core_rows": total_core,
            "kept_rows": total_kept,
            "edge_report_missing_rows": total_missing,
        },
        "per_chunk": [cs.__dict__ for cs in per_chunk],
        "outputs": {"kept_csv": str(out_kept), "flags_csv": str(out_flags), "ledger_json": str(out_ledger)},
        "columns_detected": {"src_id_col": src_col, "tile_id_col": tile_col, "ra_col": ra_col, "dec_col": dec_col},
    }
    out_ledger.write_text(json.dumps(ledger, indent=2), encoding="utf-8")

    print(f"[EDGE] [EXPERIMENTAL] chunks={len(chunks)} input_rows={total_in} "
          f"core={total_core} kept={total_kept} edge_report_missing={total_missing}")
    print(f"[EDGE] wrote: {out_kept}")
    print(f"[EDGE] wrote: {out_flags}")
    print(f"[EDGE] wrote: {out_ledger}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
