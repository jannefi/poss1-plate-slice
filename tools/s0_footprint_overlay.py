#!/usr/bin/env python3
"""RA/Dec overlay of the tile footprint and the S0 candidate positions.

A release figure. Publishing a footprint overlay alongside each release is the
convention this project inherited from VASCO60, and it is what the VASCO team
asks for: it shows at a glance which sky the run actually covered, where tiles
produced candidates and where they produced none.

BUILT ENTIRELY FROM RELEASED ARTIFACTS -- `tile_manifest.csv.gz` and
`stage_S0.csv.gz` from the release directory itself, nothing from /srv and
nothing private. So anyone who downloads the release can regenerate this exact
figure, which is the point of shipping it. There is no reference-catalogue layer
here by design; a comparison against an unpublished list would not be
reproducible by a reader and does not belong in a release figure.

Tile centres come from the tile_id via vasco.utils.tile_id.parse_tile_id_center,
the same parser the pipeline uses, so the plotted grid cannot drift from the one
that was processed.

Usage:
    python3 tools/s0_footprint_overlay.py                      # latest release
    python3 tools/s0_footprint_overlay.py --release-dir results/s0-642-20260814
    python3 tools/s0_footprint_overlay.py --out-png /tmp/x.png --dpi 200

How to validate: the printed counts must match the release's own RUN_SUMMARY.txt
(tiles, plates, S0 rows). If they do not, you are plotting a different run.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from vasco.utils.tile_id import parse_tile_id_center  # noqa: E402

# Colour-blind-safe pair, and deliberately ordered: the tile grid is context and
# sits underneath in a muted tone, the candidates are the subject and sit on top.
C_TILE = "#6699cc"
C_S0 = "#cc5500"


def stat_box(ax, x, y, ra, dec, label, colour, ha):
    n = len(ra)
    npos = int((dec > 0).sum())
    txt = (
        f"{label}\n"
        f"N = {n:,}\n"
        f"Dec > 0: {npos:,}  ({100.0 * npos / max(n, 1):.1f}%)\n"
        f"Dec ≤ 0: {n - npos:,}  ({100.0 * (n - npos) / max(n, 1):.1f}%)\n"
        f"RA: {ra.min():.1f}–{ra.max():.1f} deg\n"
        f"Dec: {dec.min():.1f}–{dec.max():.1f} deg"
    )
    ax.text(
        x, y, txt, transform=ax.transAxes, fontsize=9, color=colour,
        ha=ha, va="bottom", multialignment=ha,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                  edgecolor=colour, alpha=0.9),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-dir", default="results/s0-642-20260814",
                    help="Release directory holding tile_manifest.csv.gz and stage_S0.csv.gz")
    ap.add_argument("--out-png", default=None,
                    help="Default: <release-dir>/s0_footprint_overlay.png")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    rel = Path(args.release_dir)
    man_p, s0_p = rel / "tile_manifest.csv.gz", rel / "stage_S0.csv.gz"
    for p in (man_p, s0_p):
        if not p.exists():
            print(f"[FATAL] missing release artifact: {p}")
            return 2
    out = Path(args.out_png) if args.out_png else rel / "s0_footprint_overlay.png"

    man = pd.read_csv(man_p)
    s0 = pd.read_csv(s0_p)

    centres = man.tile_id.map(parse_tile_id_center)
    bad = int(centres.isna().sum())
    if bad:
        # Never silently drop tiles from a footprint figure: a missing tile is
        # exactly the thing the figure exists to show.
        print(f"[FATAL] {bad} of {len(man)} tile_ids did not parse to a centre; "
              f"the footprint would be understated. Fix the parser or the manifest.")
        return 2
    man["ra"] = [c[0] for c in centres]
    man["dec"] = [c[1] for c in centres]

    empty = man.rows_emitted_to_S0.fillna(0) == 0
    n_plates = man.plate_id.nunique()
    print(f"[INFO] tiles {len(man):,}  plates {n_plates}  "
          f"tiles with survivors {int((~empty).sum()):,}  empty {int(empty.sum()):,}")
    print(f"[INFO] S0 rows {len(s0):,}")
    print("[CHECK] compare the three counts above against the release RUN_SUMMARY.txt")

    # Two panels, for two different readers.
    #
    # (a) the overlay: what the run covered and where candidates fell -- the
    #     information a footprint figure exists to carry.
    # (b) S0 alone, plotted the way the VASCO team plots V: one colour, large
    #     opaque markers, near-square panel.
    #
    # Panel (b) is not decoration. A wide panel with 0.35pt markers spreads the
    # same points over ~3.7x the area with ~57x smaller marks than their
    # published style, so the identical catalogue looks perhaps an order of
    # magnitude sparser. Ours has ~14% MORE rows than V, and a reader comparing
    # our figure against theirs by eye would conclude the opposite. Matching
    # their style in one panel removes that false impression at no cost.
    fig, (ax, axb) = plt.subplots(
        2, 1, figsize=(15, 13.4), height_ratios=[1.0, 0.72])
    # Asymmetric encoding on purpose. There are 4x as many candidates as tiles,
    # so drawing both at one marker size buries the footprint under the
    # candidates -- i.e. buries the thing the figure is named for. The tile grid
    # is therefore drawn as a broad, semi-transparent FIELD (large marker, low
    # alpha) and the candidates as fine points on top, so coverage reads as an
    # area and the candidate structure reads as detail against it.
    ax.scatter(man.ra, man.dec, s=6.0, c=C_TILE, alpha=0.22, linewidths=0,
               rasterized=True, label=f"Tile centres (N={len(man):,})")
    ax.scatter(s0.ra, s0.dec, s=0.35, c=C_S0, alpha=0.45, linewidths=0,
               rasterized=True, label=f"S0 candidate positions (N={len(s0):,})")

    ax.set_xlim(0, 360)
    # Crop to the data rather than drawing the empty southern hemisphere: this
    # release is northern POSS-I by construction, and half a blank panel just
    # shrinks the part a reader needs to see.
    dlo = float(min(man.dec.min(), s0.dec.min()))
    dhi = float(max(man.dec.max(), s0.dec.max()))
    pad = 0.03 * (dhi - dlo)
    ax.set_ylim(dlo - pad, min(90.0, dhi + pad))
    ax.set_xticks(range(0, 361, 30))
    ax.set_xlabel("RA (deg)", fontsize=12)
    ax.set_ylabel("Dec (deg)", fontsize=12)
    ax.set_title(
        f"POSS-I tile footprint and S0 candidate positions — "
        f"{n_plates} plates, {len(man):,} tiles",
        fontsize=13, fontweight="bold",
    )
    ax.grid(alpha=0.25, linewidth=0.5)
    leg = ax.legend(loc="upper right", fontsize=10, markerscale=8, framealpha=0.95)
    for h in leg.legend_handles:
        h.set_alpha(1.0)

    stat_box(ax, 0.012, 0.03, man.ra.values, man.dec.values,
             "Tile centres", C_TILE, "left")
    stat_box(ax, 0.988, 0.03, s0.ra.values, s0.dec.values,
             "S0 candidates", C_S0, "right")
    ax.text(0.5, 0.03,
            f"tiles with ≥1 survivor: {int((~empty).sum()):,}\n"
            f"tiles with none: {int(empty.sum()):,}\n"
            f"mean S0 per plate: {len(s0) / max(n_plates, 1):.0f}",
            transform=ax.transAxes, fontsize=9, color="#333333",
            ha="center", va="bottom", multialignment="center",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                      edgecolor="#999999", alpha=0.9))

    # --- panel (b): S0 alone, comparable-by-eye with published V figures -------
    axb.scatter(s0.ra, s0.dec, s=14, c=C_S0, alpha=1.0, linewidths=0,
                rasterized=True)
    axb.set_xlim(0, 360)
    axb.set_ylim(dlo - pad, min(90.0, dhi + pad))
    axb.set_xticks(range(0, 361, 30))
    axb.set_xlabel("RA (deg)", fontsize=12)
    axb.set_ylabel("Dec (deg)", fontsize=12)
    axb.set_title(
        f"S0 candidate positions alone (N={len(s0):,}) — large opaque markers, "
        f"for like-for-like comparison with published transient-sample figures",
        fontsize=11,
    )
    ax.set_xlabel("")

    # Stamp the release into the image. Figures get screenshotted and pasted
    # away from the directory that explains them -- this project has already had
    # one number travel without its qualifier. A detached copy of this PNG should
    # still say which release it describes and what stage "S0" means.
    fig.text(0.995, 0.005,
             f"{rel.name} · S0 = post-veto, post-filter, deduplicated at 0.25″ · "
             f"regenerate: tools/s0_footprint_overlay.py",
             fontsize=7.5, color="#777777", ha="right", va="bottom")

    fig.tight_layout()
    fig.savefig(out, dpi=args.dpi)
    print(f"[OUT] {out}  ({out.stat().st_size / 1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
