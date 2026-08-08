#!/usr/bin/env python3
"""Dry-run vasco.mnras.filters_mnras against already-vetoed candidates.

Answers: "how much would the real (unmodified) extract+morphology filter
set remove, on top of the Gaia+PS1 veto, given single-pass has no
SPREAD_MODEL?" Reads catalogs/sextractor_pass2.after_usnob_veto.csv per
tile (the same remainder _apply_mnras_filters_and_spikes consumes in the
real pipeline) and applies apply_extract_filters + apply_morphology_filters
with the exact cfg used there. Read-only — writes nothing back to the
tile directories.

Usage:
  micromamba run -n vasco-py311 python tools/apply_mnras_filters_dryrun.py \
    --tiles-root ./data/tiles --tile-list plans/tiles_mnras_top_plate_XE524_naive.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from astropy.table import Table

from vasco.mnras.filters_mnras import apply_extract_filters, apply_morphology_filters

EXTRACT_CFG = {'flags_equal': 0, 'snr_win_min': 30.0}
MORPH_CFG = {
    'fwhm_lower': 2.0, 'fwhm_upper': 7.0, 'elongation_lt': 1.3,
    'spread_model_min': -0.002, 'sigma_clip': True, 'sigma_k': 2.0,
    'extent_delta_lt': 2.0, 'extent_min': 1.0,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--tiles-root', type=Path, default=Path('./data/tiles'))
    ap.add_argument('--tile-list', type=Path, required=True,
                    help='Plan CSV with a tile_id column')
    args = ap.parse_args()

    with args.tile_list.open(newline='', encoding='utf-8') as f:
        tiles = [row['tile_id'] for row in csv.DictReader(f)]

    has_spread = None
    tot_n0 = tot_extract = tot_morph = 0
    print(f"{'tile':28s} {'remainder':>10s} {'after_extract':>14s} {'after_morph':>12s}")
    for tid in tiles:
        p = args.tiles_root / tid / 'catalogs' / 'sextractor_pass2.after_usnob_veto.csv'
        if not p.exists() or p.stat().st_size == 0:
            continue
        try:
            tab = Table.read(str(p), format='ascii.csv')
        except Exception as e:
            print(f"{tid}: read failed ({e})")
            continue
        n0 = len(tab)
        if has_spread is None:
            has_spread = 'SPREAD_MODEL' in tab.colnames
        t1 = apply_extract_filters(tab, cfg=EXTRACT_CFG)
        n1 = len(t1)
        t2 = apply_morphology_filters(t1, cfg=MORPH_CFG)
        n2 = len(t2)
        tot_n0 += n0; tot_extract += n1; tot_morph += n2
        print(f"{tid:28s} {n0:10d} {n1:14d} {n2:12d}")

    print(f"\nSPREAD_MODEL column present: {has_spread}  (single-pass -> should be False; "
          f"that check is silently skipped by both functions if absent)")
    print(f"\nTOTAL  remainder={tot_n0}  after_extract_filters={tot_extract}  "
          f"after_morphology_filters={tot_morph}")
    if tot_n0:
        print(f"  extract_filters kept {100.0*tot_extract/tot_n0:.1f}%")
    if tot_extract:
        print(f"  morphology_filters (of extract survivors) kept {100.0*tot_morph/tot_extract:.1f}%")
    if tot_n0:
        print(f"  combined kept {100.0*tot_morph/tot_n0:.1f}% of remainder")


if __name__ == '__main__':
    main()
