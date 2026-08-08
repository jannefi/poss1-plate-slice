#!/usr/bin/env python3
"""stage_scope_dec_post.py

Run-scoped post-pipeline stage S6: science-grade subset by declination.

Sources with Dec < --dec-min are outside the science scope (default: northern
hemisphere, Dec ≥ 0°) and removed from the surviving set.

Motiviation: the VASCO60 primary comparison target (MNRAS 107k list) covers
Dec ≥ 0°.  An optional --dec-min -3.0 view includes sources that may appear
in the MNRAS list despite slightly negative declination.

If a source's Dec cannot be parsed it is conservatively KEPT and flagged with
reject_reason="dec_parse_error".

Inputs
------
- One CSV (relative to --run-dir) containing: src_id, ra, dec
  (typically stages/stage_S5_VSX.csv)

Outputs (written under <run-dir>/stages by default)
---------------------------------------------------
1) stage_S6_SCOPE_DEC.csv
   Kept survivors after declination scope gate.
   Columns: src_id, ra, dec

2) stage_S6_SCOPE_DEC_flags.csv
   Audit table for ALL input rows.
   Columns: src_id, ra, dec, dec_value, reject_reason, is_rejected

3) stage_S6_SCOPE_DEC_ledger.json
   Parameters + totals.

Usage
-----
    python scripts/stage_scope_dec_post.py --run-dir ./work/runs/run-R1-...
    python scripts/stage_scope_dec_post.py --run-dir ./work/runs/run-R1-... --dec-min -3.0
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Tuple


def main() -> int:
    ap = argparse.ArgumentParser(
        description="S6: Declination scope gate — restrict to science-grade Dec range."
    )
    ap.add_argument("--run-dir", required=True,
                    help="Run folder, e.g. ./work/runs/run-R1-...")
    ap.add_argument("--input-glob", default="stages/stage_S5_VSX.csv",
                    help="Input CSV (relative to run-dir). Default: stages/stage_S5_VSX.csv")
    ap.add_argument("--stage", default="S6",
                    help="Stage label used in output filenames. Default: S6")
    ap.add_argument("--out-dir", default=None,
                    help="Output directory. Default: <run-dir>/stages")
    ap.add_argument("--dec-min", type=float, default=0.0,
                    help="Minimum declination (degrees) to keep a source. "
                         "Sources below this threshold are rejected as out-of-scope. "
                         "Default: 0.0 (northern hemisphere)")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise SystemExit(f"run-dir not found: {run_dir}")

    out_dir = Path(args.out_dir) if args.out_dir else (run_dir / "stages")
    out_dir.mkdir(parents=True, exist_ok=True)

    inputs = sorted(run_dir.glob(args.input_glob))
    if not inputs:
        raise SystemExit(f"No inputs matched: {run_dir}/{args.input_glob}")

    stage = args.stage
    dec_min = args.dec_min

    out_kept = out_dir / f"stage_{stage}_SCOPE_DEC.csv"
    out_flags = out_dir / f"stage_{stage}_SCOPE_DEC_flags.csv"
    out_ledger = out_dir / f"stage_{stage}_SCOPE_DEC_ledger.json"

    total_in = 0
    total_kept = 0
    total_rejected = 0
    total_parse_errors = 0

    with out_kept.open("w", newline="", encoding="utf-8") as f_kept, \
         out_flags.open("w", newline="", encoding="utf-8") as f_flags:

        kept_w = csv.DictWriter(f_kept, fieldnames=["src_id", "ra", "dec"])
        kept_w.writeheader()
        flags_w = csv.DictWriter(
            f_flags,
            fieldnames=["src_id", "ra", "dec", "dec_value", "reject_reason", "is_rejected"],
        )
        flags_w.writeheader()

        for p in inputs:
            with p.open(newline="", encoding="utf-8", errors="ignore") as f:
                for row in csv.DictReader(f):
                    sid = (row.get("src_id") or "").strip()
                    if not sid:
                        continue
                    ra_str = (row.get("ra") or "").strip()
                    dec_str = (row.get("dec") or "").strip()
                    total_in += 1

                    try:
                        ra = float(ra_str)
                        dec = float(dec_str)
                    except Exception:
                        # Keep conservatively on parse error
                        kept_w.writerow({"src_id": sid, "ra": ra_str, "dec": dec_str})
                        flags_w.writerow({
                            "src_id": sid, "ra": ra_str, "dec": dec_str,
                            "dec_value": "", "reject_reason": "dec_parse_error",
                            "is_rejected": 0,
                        })
                        total_kept += 1
                        total_parse_errors += 1
                        continue

                    if dec < dec_min:
                        flags_w.writerow({
                            "src_id": sid, "ra": f"{ra:.12f}", "dec": f"{dec:.12f}",
                            "dec_value": f"{dec:.6f}",
                            "reject_reason": "dec_below_scope",
                            "is_rejected": 1,
                        })
                        total_rejected += 1
                    else:
                        kept_w.writerow({"src_id": sid, "ra": f"{ra:.12f}", "dec": f"{dec:.12f}"})
                        flags_w.writerow({
                            "src_id": sid, "ra": f"{ra:.12f}", "dec": f"{dec:.12f}",
                            "dec_value": f"{dec:.6f}", "reject_reason": "", "is_rejected": 0,
                        })
                        total_kept += 1

    ledger = {
        "run_dir": str(run_dir),
        "input_glob": args.input_glob,
        "stage": stage,
        "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "parameters": {
            "dec_min": dec_min,
        },
        "totals": {
            "input_rows": total_in,
            "kept_rows": total_kept,
            "rejected_rows": total_rejected,
            "parse_errors_kept": total_parse_errors,
        },
        "rejected_by_reason": {
            "dec_below_scope": total_rejected,
        },
        "outputs": {
            "kept_csv": str(out_kept),
            "flags_csv": str(out_flags),
            "ledger_json": str(out_ledger),
        },
    }

    out_ledger.write_text(json.dumps(ledger, indent=2), encoding="utf-8")

    print(f"[S6-DEC] input={total_in} kept={total_kept} rejected={total_rejected} "
          f"parse_errors_kept={total_parse_errors}")
    print(f"[S6-DEC] dec_min={dec_min}")
    print(f"[S6-DEC] wrote: {out_kept}")
    print(f"[S6-DEC] wrote: {out_flags}")
    print(f"[S6-DEC] wrote: {out_ledger}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
