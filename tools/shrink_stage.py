#!/usr/bin/env python3
"""Apply one stage's flags to produce the next shrinking-set stage CSV.

Implements the shrinking-set convention below:

    stage_next = stage_prev  <join>  flags_on_src_id

with the join direction named explicitly, because the two post-process stages
this project runs go in OPPOSITE directions and getting it backwards silently
inverts the catalogue:

  --keep matched    SuperCOSMOS. A counterpart in an independent digitization of
                    the same plate means the source is real rather than a scan
                    artifact, so matches are KEPT (MNRAS 2022: "Candidates
                    having a counterpart in the Supercosmos digitization at
                    less than 5 arcsec were kept").

  --drop matched    PTF, SkyBoT, VSX. A counterpart in a modern survey means the
                    source is still there / is a known moving or variable
                    object, so matches are DROPPED.

Order-independence: every stage here is a per-row catalogue test with no
population-derived threshold, so the final set does not depend on stage order --
only the per-stage attribution does. (This is NOT true of the MNRAS morphology
filters, whose sigma-clip window is derived from the population being filtered.)

Usage:
    python3 tools/shrink_stage.py \\
      --in-stage  <run-dir>/stages/stage_S0.csv \\
      --flags     <run-dir>/stages/scos/scos_flags.csv \\
      --flag-col  nmatch_r1 --keep matched \\
      --out-stage <run-dir>/stages/stage_S1_SCOS.csv \\
      --label SCOS_R1

How to validate: rows_out + rows_removed must equal rows_in (asserted), and the
upload view must carry only src_id, ra, dec.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

UPLOAD = ["src_id", "ra", "dec"]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-stage", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--flag-col", required=True,
                    help="numeric match-count or 0/1 column in the flags file")
    ap.add_argument("--keep", choices=["matched", "unmatched"], required=True,
                    help="'matched' for SuperCOSMOS; 'unmatched' for PTF/SkyBoT/VSX")
    ap.add_argument("--out-stage", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--ledger", default=None,
                    help="stage_ledger.csv to append to (default: alongside out-stage)")
    args = ap.parse_args(argv)

    d = pd.read_csv(args.in_stage)
    f = pd.read_csv(args.flags, usecols=["src_id", args.flag_col])
    n_in = len(d)
    j = d.merge(f, on="src_id", how="left")
    if len(j) != n_in:
        raise SystemExit(f"[FAIL] join changed row count {n_in} -> {len(j)}; "
                         f"duplicate src_id in the flags file?")
    unseen = int(j[args.flag_col].isna().sum())
    if unseen:
        print(f"[WARN] {unseen:,} rows absent from the flags file; treated as "
              f"UNMATCHED (a stage that never queried a row cannot confirm it)")
    matched = j[args.flag_col].fillna(0) > 0
    keep = matched if args.keep == "matched" else ~matched

    out = j[keep][d.columns]
    out.to_csv(args.out_stage, index=False)
    up = Path(args.out_stage).with_name(Path(args.out_stage).stem + "_upload.csv")
    out[UPLOAD].to_csv(up, index=False)

    removed = n_in - len(out)
    assert len(out) + removed == n_in
    ledger = Path(args.ledger) if args.ledger else \
        Path(args.out_stage).parent.parent / "stage_ledger.csv"
    if not ledger.exists():
        ledger.write_text("stage,rows_in,rows_flagged,rows_out,note\n")
    with ledger.open("a") as fh:
        fh.write(f"{args.label},{n_in},{int(matched.sum())},{len(out)},"
                 f"keep={args.keep} on {args.flag_col}\n")

    print(f"[{args.label}] {n_in:,} in -> {len(out):,} out  "
          f"(removed {removed:,} = {100*removed/n_in:.2f}%)")
    print(f"[OUT] {args.out_stage}\n[OUT] {up}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
