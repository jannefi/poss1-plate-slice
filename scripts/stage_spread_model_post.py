#!/usr/bin/env python3
"""stage_spread_model_post.py

Run-scoped post-pipeline stage S0G: SPREAD_MODEL gate.

The primary pipeline runs single-pass SExtractor (no PSFEx, no SPREAD_MODEL
column) to keep candidate generation cheap. `tools/spread_model_postscore.py`
restores a real, PSF-based SPREAD_MODEL measurement for each tile's survivors
via a validated crop + synthetic-PSF remeasurement, and writes it to
`catalogs/spread_model_postscore.csv` per tile. This stage is a thin CSV-join
that reads that already-computed file and applies the locked gate
(`SPREAD_MODEL > --spread-model-min`, default -0.002, context/02_DECISIONS.md)
-- it does NOT run PSFEx/SExtractor itself, so the multi-hour measurement
cost is paid once (by `spread_model_postscore.py`), not on every
post-process run.

Prerequisite
------------
Run `tools/spread_model_postscore.py --tile-ids-file ... --tiles-root ...
[--skip-existing]` against every tile referenced by the input stage CSV
BEFORE running this stage, so `catalogs/spread_model_postscore.csv` exists
per tile. This stage does not fail if that file is missing for some tiles --
see "Conservative-keep policy" below -- but the resulting kept count is then
only an upper bound, not a final number (see `postscore_coverage` in the
ledger).

Conservative-keep policy
-------------------------
Following the same conservative-keep-and-flag convention as
`stage_scope_dec_post.py` (dec_parse_error) and `stage_morph_post.py`
(psf_insufficient) -- never silent-drop, never hard-error:

  - Tile has no spread_model_postscore.csv yet  -> KEEP, reason=postscore_missing
  - object_id not found in that tile's postscore rows -> KEEP, reason=postscore_row_not_found
  - Row found but matched=False (no two-pass detection found) -> KEEP, reason=postscore_unmatched
  - Row found, matched=True, SPREAD_MODEL <= --spread-model-min -> REJECT, reason=spread_model_below_gate
  - Row found, matched=True, SPREAD_MODEL >  --spread-model-min -> KEEP, reason=(none)

The gate is recomputed here from the raw SPREAD_MODEL value against
--spread-model-min, rather than trusting the persisted gate_pass column
verbatim, so a future threshold change is a cheap stage re-run rather than
a multi-hour recompute of spread_model_postscore.py.

Inputs
------
- One or more CSVs (relative to --run-dir, via --input-glob) containing:
  src_id, ra, dec  (src_id = "tile_id:object_id", the repo-wide convention)
  Typically 'stage_S0.csv' -- this stage is designed to run first, directly
  on S0's output.

Outputs (written under <run-dir>/stages by default)
-----------------------------------------------------
1) stage_<STAGE>_SPREAD.csv
   Kept survivors after the SPREAD_MODEL gate.
   Columns: src_id, ra, dec

2) stage_<STAGE>_SPREAD_flags.csv
   Audit table for ALL input rows.
   Columns: src_id, ra, dec, tile_id, object_id, spread_model,
            spreaderr_model, sep_px, postscore_status, gate_pass,
            reject_reason, is_rejected

3) stage_<STAGE>_SPREAD_ledger.json
   Parameters, totals, rejected_by_reason breakdown, and a
   postscore_coverage block (tiles_referenced, tiles_with_postscore,
   tiles_missing_postscore, plus per-reason keep counts) so a run's
   SPREAD_MODEL coverage -- and therefore whether its kept count is final
   or still an upper bound -- is auditable from the ledger alone.

Usage
-----
    python scripts/stage_spread_model_post.py \\
        --run-dir ./work/runs/run-S0G-... \\
        --input-glob 'stage_S0.csv' \\
        --stage S0G \\
        --tiles-root ./data/tiles_archive
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def _read_header(path: Path) -> List[str]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.reader(f)
        return [c.strip().lstrip("\ufeff") for c in next(r, [])]


def _detect_src_cols(cols: List[str]) -> Tuple[str, str, str]:
    cset = set(cols)
    if "src_id" in cset:
        src = "src_id"
    elif "row_id" in cset:
        src = "row_id"
    else:
        raise RuntimeError("Input CSV missing required id column 'src_id' (or 'row_id')")
    ra = next((c for c in ["ra", "RA", "RA_corr", "ALPHAWIN_J2000"] if c in cset), None)
    dec = next((c for c in ["dec", "DEC", "Dec_corr", "DELTAWIN_J2000"] if c in cset), None)
    if not ra or not dec:
        raise RuntimeError("Input CSV missing RA/Dec columns")
    return src, ra, dec


def _parse_src_id(src_id: str) -> Tuple[str, str]:
    """Parse 'tile_id:object_id' -> (tile_id, object_id)."""
    parts = src_id.split(":", 1)
    return (parts[0], parts[1]) if len(parts) == 2 else (src_id, "")


def _norm_object_id(object_id: str) -> str:
    """Normalize an object id to a plain int string ('9436', not '9436.0')."""
    try:
        return str(int(float(object_id)))
    except Exception:
        return object_id.strip()


def _load_postscore(tile_id: str, tiles_root: Path, postscore_name: str) -> Optional[Dict[str, dict]]:
    """Return {object_id_str: row} for a tile's postscore CSV, or None if absent/empty."""
    p = tiles_root / tile_id / postscore_name
    if not p.exists() or p.stat().st_size == 0:
        return None
    out: Dict[str, dict] = {}
    try:
        with p.open(newline="", encoding="utf-8", errors="ignore") as f:
            for row in csv.DictReader(f):
                num = _norm_object_id(row.get("survivor_number", ""))
                if num:
                    out[num] = row
    except Exception:
        return None
    return out


_FLAGS_FIELDS = [
    "src_id", "ra", "dec", "tile_id", "object_id",
    "spread_model", "spreaderr_model", "sep_px",
    "postscore_status", "gate_pass", "reject_reason", "is_rejected",
]


def _fmt(v) -> str:
    if v is None or v == "":
        return ""
    try:
        fv = float(v)
    except Exception:
        return str(v)
    if fv != fv:  # NaN
        return ""
    return f"{fv:.6g}"


def _evaluate(src_id: str, ra: str, dec: str, tile_id: str, object_id: str,
              postscore: Optional[Dict[str, dict]], spread_model_min: float) -> dict:
    base = {"src_id": src_id, "ra": ra, "dec": dec, "tile_id": tile_id, "object_id": object_id}

    if postscore is None:
        return {**base, "spread_model": "", "spreaderr_model": "", "sep_px": "",
                "postscore_status": "postscore_missing", "gate_pass": "",
                "reject_reason": "postscore_missing", "is_rejected": 0}

    row = postscore.get(_norm_object_id(object_id))
    if row is None:
        return {**base, "spread_model": "", "spreaderr_model": "", "sep_px": "",
                "postscore_status": "row_not_found", "gate_pass": "",
                "reject_reason": "postscore_row_not_found", "is_rejected": 0}

    matched = str(row.get("matched", "")).strip().lower() in ("true", "1")
    if not matched:
        return {**base, "spread_model": "", "spreaderr_model": "", "sep_px": _fmt(row.get("sep_px")),
                "postscore_status": "unmatched", "gate_pass": "",
                "reject_reason": "postscore_unmatched", "is_rejected": 0}

    try:
        sm = float(row["SPREAD_MODEL"])
    except Exception:
        return {**base, "spread_model": "", "spreaderr_model": "", "sep_px": _fmt(row.get("sep_px")),
                "postscore_status": "unmatched", "gate_pass": "",
                "reject_reason": "postscore_unmatched", "is_rejected": 0}

    gate_pass = sm > spread_model_min
    return {**base,
            "spread_model": _fmt(sm),
            "spreaderr_model": _fmt(row.get("SPREADERR_MODEL")),
            "sep_px": _fmt(row.get("sep_px")),
            "postscore_status": "ok",
            "gate_pass": gate_pass,
            "reject_reason": "" if gate_pass else "spread_model_below_gate",
            "is_rejected": 0 if gate_pass else 1}


def main() -> int:
    ap = argparse.ArgumentParser(description="S0G: SPREAD_MODEL gate (thin join over precomputed postscore CSVs).")
    ap.add_argument("--run-dir", required=True, help="Run folder, e.g. ./work/runs/run-S0G-...")
    ap.add_argument("--input-glob", default="stage_S0.csv",
                     help="Input CSV glob, relative to run-dir. Default: stage_S0.csv")
    ap.add_argument("--stage", default="S0G", help="Stage label used in output filenames. Default: S0G")
    ap.add_argument("--out-dir", default=None, help="Output directory. Default: <run-dir>/stages")
    ap.add_argument("--tiles-root", default="./data/tiles",
                     help="Root of tile directories. Default: ./data/tiles (use ./data/tiles_archive "
                          "for tiles already archived off the SSD).")
    ap.add_argument("--postscore-name", default="catalogs/spread_model_postscore.csv",
                     help="Relative path under each tile dir to the precomputed postscore CSV "
                          "written by tools/spread_model_postscore.py.")
    ap.add_argument("--spread-model-min", type=float, default=-0.002,
                     help="Locked SPREAD_MODEL gate threshold (context/02_DECISIONS.md). Default: -0.002")
    ap.add_argument("--verbose", action="store_true", help="Print per-tile coverage as it's loaded.")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise SystemExit(f"run-dir not found: {run_dir}")

    tiles_root = Path(args.tiles_root)
    out_dir = Path(args.out_dir) if args.out_dir else (run_dir / "stages")
    out_dir.mkdir(parents=True, exist_ok=True)

    chunks = sorted(run_dir.glob(args.input_glob))
    if not chunks:
        raise SystemExit(f"No inputs matched: {run_dir}/{args.input_glob}")

    src_col, ra_col, dec_col = _detect_src_cols(_read_header(chunks[0]))

    stage = args.stage
    out_kept = out_dir / f"stage_{stage}_SPREAD.csv"
    out_flags = out_dir / f"stage_{stage}_SPREAD_flags.csv"
    out_ledger = out_dir / f"stage_{stage}_SPREAD_ledger.json"

    # Pass 1: read all input rows, group object_ids by tile.
    all_rows: List[dict] = []
    tiles_needed: Dict[str, List[str]] = {}
    for ch in chunks:
        with ch.open(newline="", encoding="utf-8", errors="ignore") as f:
            for row in csv.DictReader(f):
                sid = (row.get(src_col) or "").strip()
                ra = (row.get(ra_col) or "").strip()
                dec = (row.get(dec_col) or "").strip()
                if not sid or not ra or not dec:
                    continue
                tile_id, object_id = _parse_src_id(sid)
                all_rows.append({"src_id": sid, "ra": ra, "dec": dec,
                                  "tile_id": tile_id, "object_id": object_id})
                tiles_needed.setdefault(tile_id, []).append(object_id)

    total_in = len(all_rows)
    print(f"[S0G-SPREAD] input_rows={total_in} tiles={len(tiles_needed)}")

    # Pass 2: lazy-load each referenced tile's postscore CSV once.
    tile_postscore: Dict[str, Optional[Dict[str, dict]]] = {}
    n_with_postscore = 0
    for tile_id in tiles_needed:
        ps = _load_postscore(tile_id, tiles_root, args.postscore_name)
        tile_postscore[tile_id] = ps
        if ps is not None:
            n_with_postscore += 1
        if args.verbose:
            print(f"[S0G-SPREAD]   {tile_id}: "
                  f"{'postscore ok (' + str(len(ps)) + ' rows)' if ps is not None else 'MISSING'}")

    n_missing_tiles = len(tiles_needed) - n_with_postscore
    print(f"[S0G-SPREAD] tiles: with_postscore={n_with_postscore} missing={n_missing_tiles}")

    # Pass 3: evaluate every row.
    evaluated: List[dict] = [
        _evaluate(r["src_id"], r["ra"], r["dec"], r["tile_id"], r["object_id"],
                  tile_postscore[r["tile_id"]], args.spread_model_min)
        for r in all_rows
    ]

    total_rejected = sum(1 for e in evaluated if e["is_rejected"] == 1)
    total_kept = total_in - total_rejected

    rejected_by_reason: Dict[str, int] = {}
    kept_by_reason: Dict[str, int] = {}
    for e in evaluated:
        reason = e["reject_reason"] or "ok"
        bucket = rejected_by_reason if e["is_rejected"] == 1 else kept_by_reason
        bucket[reason] = bucket.get(reason, 0) + 1

    # Write outputs.
    with out_kept.open("w", newline="", encoding="utf-8") as fk, \
         out_flags.open("w", newline="", encoding="utf-8") as ff:
        kept_w = csv.DictWriter(fk, fieldnames=["src_id", "ra", "dec"])
        flags_w = csv.DictWriter(ff, fieldnames=_FLAGS_FIELDS)
        kept_w.writeheader()
        flags_w.writeheader()
        for e in evaluated:
            flags_w.writerow(e)
            if e["is_rejected"] == 0:
                kept_w.writerow({"src_id": e["src_id"], "ra": e["ra"], "dec": e["dec"]})

    ledger = {
        "run_dir": str(run_dir),
        "input_glob": args.input_glob,
        "stage": stage,
        "tiles_root": str(tiles_root),
        "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "parameters": {
            "spread_model_min": args.spread_model_min,
            "postscore_name": args.postscore_name,
        },
        "input_chunks": [p.name for p in chunks],
        "totals": {
            "input_rows": total_in,
            "kept_rows": total_kept,
            "rejected_rows": total_rejected,
            "rejection_pct": round(100.0 * total_rejected / total_in, 2) if total_in > 0 else 0.0,
        },
        "rejected_by_reason": rejected_by_reason,
        "kept_by_reason": kept_by_reason,
        "postscore_coverage": {
            "tiles_referenced": len(tiles_needed),
            "tiles_with_postscore": n_with_postscore,
            "tiles_missing_postscore": n_missing_tiles,
            "rows_kept_via_postscore_missing": kept_by_reason.get("postscore_missing", 0),
            "rows_kept_via_row_not_found": kept_by_reason.get("postscore_row_not_found", 0),
            "rows_kept_via_unmatched": kept_by_reason.get("postscore_unmatched", 0),
            "note": "kept_rows is an UPPER BOUND on the SPREAD_MODEL-adjusted count "
                    "unless tiles_missing_postscore == 0 and the three rows_kept_via_* "
                    "counts above are all 0 -- those rows were never actually evaluated "
                    "against the gate.",
        },
        "outputs": {
            "kept_csv": str(out_kept),
            "flags_csv": str(out_flags),
            "ledger_json": str(out_ledger),
        },
        "columns_detected": {"src_id_col": src_col, "ra_col": ra_col, "dec_col": dec_col},
    }
    out_ledger.write_text(json.dumps(ledger, indent=2), encoding="utf-8")

    pct = 100.0 * total_rejected / total_in if total_in > 0 else 0.0
    print(f"[S0G-SPREAD] input={total_in} rejected={total_rejected} kept={total_kept} ({pct:.1f}% reduction)")
    print(f"[S0G-SPREAD] postscore coverage: {n_with_postscore}/{len(tiles_needed)} tiles "
          f"({n_missing_tiles} missing -- kept count is an upper bound if > 0)")
    print(f"[S0G-SPREAD] wrote: {out_kept}")
    print(f"[S0G-SPREAD] wrote: {out_flags}")
    print(f"[S0G-SPREAD] wrote: {out_ledger}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
