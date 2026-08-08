#!/usr/bin/env python3
"""Choose a dedup tolerance valid for the coordinates the pipeline actually has.

The problem this answers
------------------------
`context/02_DECISIONS.md` locks dedup at 0.25 arcsec *under WCS-fixed
coordinates*. WCS-fix is a deliberately ACTIVE deviation
(`context/REPRO_DEVIATIONS.md` item 6) -- raw plate WCS is used as-is, because
MNRAS 2022 describes no coordinate-correction step and enabling the fix measured
WORSE parity project-wide. Only ~1.5% of archive tiles carry any wcsfix
artifact.

So the locked tolerance is being applied outside its stated domain. Nothing is
broken; a precondition was dropped when WCS-fix was disabled and the tolerance
was never revisited. The consequence is measurable: the full-sky S0 holds 16,055
pairs within 5 arcsec, and 100% of them are CROSS-TILE -- the same source seen
through two overlapping tiles, each carrying its own unaligned WCS solution.

Why widening is safe here, and how this checks it
-------------------------------------------------
Genuine close pairs essentially do not exist at POSS-I resolution: of 2,385,582
raw SExtractor detections across 400 tiles, only 96 (0.004%) have a neighbour
within 5 arcsec. So a wider tolerance merges duplicates rather than distinct
sources. This sweep verifies that rather than assuming it, by tracking:

  cross_tile_pairs   should RISE then PLATEAU once the tolerance exceeds typical
                     inter-tile astrometric scatter -- the plateau knee is the
                     radius to pick, not a round number
  intra_tile_pairs   must stay ~0. Intra-tile pairs are the signature of merging
                     genuinely distinct sources, since SExtractor does not emit
                     two detections that close within one tile
  cluster sizes      clusters of 3+ become common at wider radii; the
                     representative-selection rule must still behave

Read-only: computes counts, writes no run outputs.

How to validate
---------------
    sudo systemd-run --scope -p MemoryMax=24G -p MemorySwapMax=0 \\
      --uid=janne --gid=janne --working-directory="$PWD" --quiet \\
      python3 \\
      tools/dedup_radius_sweep.py \\
        --stage-csv work/runs/full635-20260802-dedupfix/stage_S0.csv \\
        --out-dir work/dedup_sweep

Pick the radius where cross-tile pairs flatten AND intra-tile pairs are still
~0. If intra-tile climbs before the plateau, widening is unsafe and the tolerance
must stay tight regardless of the duplicate count.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def xyz(ra, dec):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def chord(arcsec: float) -> float:
    return 2.0 * math.sin(math.radians(arcsec / 3600.0) / 2.0)


def components(n, pairs):
    """Union-find over the pair graph -> cluster sizes and kept count."""
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in pairs:
        ra_, rb = find(a), find(b)
        if ra_ != rb:
            parent[rb] = ra_
    sizes = Counter()
    for i in range(n):
        sizes[find(i)] += 1
    return Counter(sizes.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stage-csv", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--radii", default="0.25,0.5,1,1.5,2,2.5,3,4,5,6,8")
    args = ap.parse_args()

    ra, dec, tile = [], [], []
    with Path(args.stage_csv).open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                ra.append(float(r["ra"])); dec.append(float(r["dec"]))
            except (TypeError, ValueError, KeyError):
                continue
            tile.append(r.get("tile_id", ""))
    ra, dec = np.asarray(ra), np.asarray(dec)
    tile = np.asarray(tile)
    n = len(ra)
    print(f"[IN] {n} rows from {Path(args.stage_csv).name}")

    tree = cKDTree(xyz(ra, dec))
    radii = [float(x) for x in args.radii.split(",") if x.strip()]
    rows = []
    print(f"\n{'radius':>8} {'pairs':>9} {'cross':>9} {'intra':>7} "
          f"{'kept':>9} {'removed':>9} {'max_clust':>10}")
    print("-" * 68)
    prev_cross = None
    for R in radii:
        pairs = list(tree.query_pairs(r=chord(R)))
        cross = sum(1 for a, b in pairs if tile[a] != tile[b])
        intra = len(pairs) - cross
        sizes = components(n, pairs)
        kept = sum(sizes.values())          # one representative per component
        maxc = max(sizes) if sizes else 1
        growth = "" if prev_cross is None else f"  (+{cross - prev_cross})"
        prev_cross = cross
        rows.append(dict(radius_arcsec=R, pairs=len(pairs), cross_tile=cross,
                         intra_tile=intra, kept=kept, removed=n - kept,
                         max_cluster=maxc,
                         cluster_hist={str(k): v for k, v in sorted(sizes.items())}))
        print(f"{R:>8.2f} {len(pairs):>9d} {cross:>9d} {intra:>7d} "
              f"{kept:>9d} {n-kept:>9d} {maxc:>10d}{growth}")

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    (Path(args.out_dir) / "dedup_radius_sweep.json").write_text(json.dumps({
        "input": str(args.stage_csv), "n_rows": n, "sweep": rows,
        "guidance": "Pick the radius where cross_tile plateaus AND intra_tile is "
                    "still ~0. intra_tile rising before the plateau means genuine "
                    "close pairs are being merged and widening is unsafe.",
    }, indent=2))
    print(f"\nwrote {args.out_dir}/dedup_radius_sweep.json")


if __name__ == "__main__":
    main()
