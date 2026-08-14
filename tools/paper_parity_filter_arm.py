#!/usr/bin/env python3
"""How much does the third veto (USNO-B) change the catalogue?

Solano et al. (2022) veto against two catalogues, Gaia and PS1. This pipeline adds
USNO-B, which is a deliberate deviation. Anyone comparing against a paper-parity
implementation needs to know its size.

It cannot be read off the released catalogue. The vetoes run BEFORE the MNRAS
filters, and `_robust_sigma_clip` takes its median and MAD from the population
being filtered -- so dropping a veto changes the population, changes the clip
window, and changes which rows the FILTERS cut. Subtracting the USNO-B removals
from the final catalogue would give the wrong answer.

What this tool does instead is exact: it re-runs the filter stage from
`sextractor_pass2.after_ps1_veto.csv` -- the population as it stood before the
USNO-B veto -- so the clip window is re-derived from the correct population, as
the pipeline would have derived it.

VALIDATION IS THE POINT, not an afterthought. The same code path is first run on
`after_usnob_veto.csv`, which must reproduce the tile's actual
`sextractor_pass2.filtered.csv` exactly. A tile that fails that control is
reported and excluded rather than quietly averaged in: if the reimplementation
does not reproduce the pipeline, its paper-parity arm means nothing.

Reads only. Never writes into a tile.

Requires tiles that retained the full per-stage chain (`--keep-tiles`).

Usage:
    python3 tools/paper_parity_filter_arm.py \
      --tiles-root work/slice/tiles --out work/paper_parity_arm.csv
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

# Exactly the configuration cli_pipeline.py's _apply_mnras_filters_and_spikes
# uses. Divergence here silently invalidates the comparison, which is what the
# control arm below exists to catch.
EXTRACT_CFG = {"flags_equal": 0, "snr_win_min": 30.0}
MORPH_CFG = {
    "fwhm_lower": 2.0, "fwhm_upper": 7.0, "elongation_lt": 1.3,
    "spread_model_min": -0.002, "sigma_clip": True, "sigma_k": 2.0,
    "extent_delta_lt": 2.0, "extent_min": 1.0,
}


def _filter_chain(csv_path: Path, bright_cache: Path):
    """Run extract -> morphology -> spikes and return the surviving NUMBERs."""
    from astropy.table import Table
    from vasco.mnras.filters_mnras import (apply_extract_filters,
                                           apply_morphology_filters)
    from vasco.mnras.spikes import (SpikeConfig, SpikeRuleConst, SpikeRuleLine,
                                    apply_spike_cuts)
    from vasco.cli_pipeline import _read_bright_cache

    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return set()
    try:
        tab = Table.read(str(csv_path), format="ascii.csv")
    except Exception:
        return set()
    if len(tab) == 0:
        return set()

    tab = apply_extract_filters(tab, cfg=dict(EXTRACT_CFG))
    tab = apply_morphology_filters(tab, cfg=dict(MORPH_CFG))
    if len(tab) == 0:
        return set()

    bright = []
    if bright_cache.exists() and bright_cache.stat().st_size > 0:
        try:
            bright = _read_bright_cache(bright_cache)
        except Exception:
            bright = []
    rows = [dict(zip(tab.colnames, r)) for r in tab.as_array()]
    kept, _ = apply_spike_cuts(
        rows, bright,
        SpikeConfig(search_radius_arcmin=1.5,
                    rules=[SpikeRuleConst(const_max_mag=12.4),
                           SpikeRuleLine(a=-0.09, b=15.3)]))
    return {int(r["NUMBER"]) for r in kept}


def _numbers(path: Path) -> set[int]:
    try:
        if path.stat().st_size == 0:
            return set()
        return set(pd.read_csv(path, usecols=["NUMBER"]).NUMBER.astype(int))
    except Exception:
        return set()


def do_tile(job):
    tile_id, root = job
    d = Path(root) / tile_id / "catalogs"
    bright = d / "ps1_bright_stars_r16_rad45.csv"
    actual = _numbers(d / "sextractor_pass2.filtered.csv")

    # Control: reproduce the pipeline from the same input it used.
    control = _filter_chain(d / "sextractor_pass2.after_usnob_veto.csv", bright)
    faithful = (control == actual)

    # Paper-parity: same chain, but from before the USNO-B veto.
    paper = _filter_chain(d / "sextractor_pass2.after_ps1_veto.csv", bright)

    return {
        "tile_id": tile_id,
        "released_rows": len(actual),
        "control_rows": len(control),
        "control_faithful": bool(faithful),
        "paper_parity_rows": len(paper),
        "added_by_dropping_usnob": len(paper - actual),
        "lost_by_dropping_usnob": len(actual - paper),
        "in_after_ps1": len(_numbers(d / "sextractor_pass2.after_ps1_veto.csv")),
        "in_after_usnob": len(_numbers(d / "sextractor_pass2.after_usnob_veto.csv")),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles-root", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = Path(args.tiles_root)
    tiles = sorted(p.name for p in root.iterdir() if p.is_dir())
    print(f"[CONFIG] {len(tiles)} tiles under {root}", flush=True)
    for v in ("VASCO_DISABLE_MNRAS_FILTERS", "VASCO_DISABLE_USNOB",
              "VASCO_CIRCLE_ARCMIN"):
        if os.getenv(v):
            print(f"[CONFIG][WARN] {v}={os.getenv(v)} set -- the filter config "
                  f"here is hardcoded to the released settings, so a leaked env "
                  f"var will make the control arm disagree.")

    rows, done = [], 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(do_tile, (t, str(root))) for t in tiles]
        for f in as_completed(futs):
            rows.append(f.result())
            done += 1
            if done % 100 == 0:
                print(f"  {done}/{len(tiles)}", flush=True)

    df = pd.DataFrame(rows)
    ok = df[df.control_faithful]
    bad = df[~df.control_faithful]

    print(f"\n{'='*78}\nCONTROL: does this reproduce the pipeline?\n{'='*78}")
    print(f"  tiles reproducing filtered.csv exactly : {len(ok)}/{len(df)}"
          f"  ({100.0*len(ok)/len(df):.2f}%)")
    if len(bad):
        print(f"  tiles FAILING the control              : {len(bad)}  "
              f"-- excluded from the numbers below")
        print(f"     released vs control row totals: {bad.released_rows.sum()}"
              f" vs {bad.control_rows.sum()}")
    if len(ok) == 0:
        print("\n[FATAL] no tile reproduced the pipeline. The paper-parity arm "
              "cannot be trusted and no number is reported.")
        return 1

    rel, pap = int(ok.released_rows.sum()), int(ok.paper_parity_rows.sum())
    print(f"\n{'='*78}\nTHE DEVIATION: USNO-B as the third veto\n{'='*78}")
    print(f"  {len(ok)} validated tiles")
    print(f"  released      (Gaia + PS1 + USNO-B) : {rel:8d} rows")
    print(f"  paper parity  (Gaia + PS1)          : {pap:8d} rows")
    print(f"  difference                          : {pap-rel:+8d} rows"
          f"  ({100.0*(pap-rel)/rel:+.2f}%)")
    print(f"  ratio paper/released                : {pap/rel:8.3f}x")
    print(f"\n  rows added by dropping USNO-B  : {int(ok.added_by_dropping_usnob.sum()):7d}")
    print(f"  rows LOST by dropping USNO-B   : {int(ok.lost_by_dropping_usnob.sum()):7d}"
          f"   <- non-zero only because the clip window moves")
    lost = int(ok.lost_by_dropping_usnob.sum())
    if lost:
        print(f"\n  That {lost} is the whole reason this needed a re-run rather than a")
        print(f"  subtraction: dropping a veto does not only ADD rows. The wider")
        print(f"  population shifts the sigma-clip window and removes rows the")
        print(f"  released catalogue kept.")
    print(f"\n  veto-stage removals for reference: after_ps1 "
          f"{int(ok.in_after_ps1.sum())} -> after_usnob "
          f"{int(ok.in_after_usnob.sum())} "
          f"({100.0*(1-ok.in_after_usnob.sum()/max(1,ok.in_after_ps1.sum())):.2f}% cut)")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, index=False)
        print(f"\n[OUT] {args.out}")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    raise SystemExit(main())
