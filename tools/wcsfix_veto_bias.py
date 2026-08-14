#!/usr/bin/env python3
"""Does the Gaia-bootstrapped WCSFIX refit inflate the Gaia veto?

The published catalogue's coordinates come from a per-tile degree-2 refit
bootstrapped off Gaia, and Gaia is then the first 5" veto. The release README
concedes that a coherent tile-wide pull toward Gaia's frame could raise Gaia-veto
rates, and that it had not been quantified. This quantifies it.

Measuring the Gaia rate alone cannot answer the question, because a rise is
expected either way: the raw plate solution has real residual error, so aligning
to Gaia genuinely improves the astrometry and *should* find more true
counterparts. "Circular" would mean the refit snaps toward Gaia BEYOND that
genuine improvement.

THE DISCRIMINATOR: catalogues the refit never saw. WCSFIX bootstraps off Gaia
only. PS1 and USNO-B sit on essentially the same ICRS frame but played no part in
the fit, so:

  * genuine astrometric improvement -> match rates rise by a COMPARABLE amount
    against all three catalogues, and median separations fall for all three.
  * Gaia-specific overfitting      -> the Gaia rise EXCEEDS the PS1/USNO-B rises,
    and the excess is the size of the circularity.

Both arms are also run through a null control (detections displaced in RA) to
show the chance rate is unchanged by the refit, which it must be.

Everything needed is already inside each tile: the wcsfix catalogue carries the
raw (ALPHAWIN_J2000) and refit (RA_corr) positions side by side, and the tile
holds its own Gaia / PS1 / USNO-B neighbourhood catalogues.

Usage:
    python3 tools/wcsfix_veto_bias.py \
      --tiles-root work/slice/tiles \
      --max-tiles 250 --out work/wcsfix_veto_bias.csv
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# Reference catalogues. Gaia is what WCSFIX bootstrapped off; the other two are
# the independent controls. All three are matched at the plate epoch where a
# plate-epoch version exists, matching what the veto itself does.
REFS = {
    "gaia":  ("gaia_neighbourhood_at_plate.csv",  "bootstrap catalogue"),
    "ps1":   ("ps1_neighbourhood.csv",            "control, unseen by the fit"),
    "usnob": ("usnob_neighbourhood_at_plate.csv", "control, unseen by the fit"),
}


def _proj(ra, dec, ra0, dec0):
    cosd = np.cos(np.deg2rad(dec0))
    dra = ((np.asarray(ra, float) - ra0 + 180.0) % 360.0 - 180.0) * cosd
    return np.column_stack([dra * 3600.0, (np.asarray(dec, float) - dec0) * 3600.0])


def do_tile(job):
    tile_id, root, radius, cap, null_arcmin, seed = job
    d = Path(root) / tile_id / "catalogs"
    ra0 = float(tile_id.split("_RA")[1].split("_")[0])
    ds = tile_id.split("_DEC")[1]
    dec0 = float(ds[1:]) * (1.0 if ds[0] == "p" else -1.0)

    try:
        det = pd.read_csv(d / "sextractor_pass2.wcsfix.csv",
                          usecols=["ALPHAWIN_J2000", "DELTAWIN_J2000",
                                   "RA_corr", "Dec_corr"]).dropna()
    except Exception:
        return None
    if len(det) < 50:
        return None
    if len(det) > cap:
        det = det.iloc[np.random.RandomState(seed).choice(len(det), cap,
                                                          replace=False)]

    shift = null_arcmin / 60.0 / np.cos(np.deg2rad(dec0))
    arms = {
        "raw":   (det.ALPHAWIN_J2000.values, det.DELTAWIN_J2000.values),
        "refit": (det.RA_corr.values, det.Dec_corr.values),
    }

    out = {"tile_id": tile_id, "dec": dec0, "n_det": len(det)}
    # How far the refit actually moved this tile -- the covariate to check the
    # bias against.
    out["shift_median_as"] = float(np.median(np.hypot(
        (det.RA_corr.values - det.ALPHAWIN_J2000.values)
        * np.cos(np.deg2rad(det.DELTAWIN_J2000.values)),
        det.Dec_corr.values - det.DELTAWIN_J2000.values) * 3600.0))

    for ref, (fname, _) in REFS.items():
        try:
            r = pd.read_csv(d / fname, usecols=["ra", "dec"]).dropna()
        except Exception:
            continue
        if len(r) < 10:
            continue
        tree = cKDTree(_proj(r.ra.values, r.dec.values, ra0, dec0))
        for arm, (ra_v, dec_v) in arms.items():
            for tag, off in (("", 0.0), ("_null", shift)):
                dist, _ = tree.query(_proj(ra_v + off, dec_v, ra0, dec0))
                out[f"{ref}_{arm}{tag}_frac"] = float(np.mean(dist <= radius))
                if not tag:
                    out[f"{ref}_{arm}_medsep"] = float(np.median(dist))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles-root", required=True)
    ap.add_argument("--radius-arcsec", type=float, default=5.0,
                    help="The veto radius. Default matches the pipeline.")
    ap.add_argument("--max-tiles", type=int, default=250)
    ap.add_argument("--max-detections", type=int, default=1500,
                    help="Per-tile sample cap, for I/O.")
    ap.add_argument("--null-shift-arcmin", type=float, default=6.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = Path(args.tiles_root)
    tiles = sorted(p.name for p in root.iterdir() if p.is_dir())
    if len(tiles) > args.max_tiles:
        step = len(tiles) // args.max_tiles
        tiles = tiles[::step][:args.max_tiles]
    print(f"[CONFIG] {len(tiles)} tiles, radius {args.radius_arcsec}\", "
          f"cap {args.max_detections} det/tile, null shift "
          f"{args.null_shift_arcmin}'", flush=True)

    rows, done = [], 0
    jobs = [(t, str(root), args.radius_arcsec, args.max_detections,
             args.null_shift_arcmin, 0) for t in tiles]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for f in as_completed([ex.submit(do_tile, j) for j in jobs]):
            r = f.result()
            if r:
                rows.append(r)
            done += 1
            if done % 50 == 0:
                print(f"  {done}/{len(tiles)}", flush=True)

    df = pd.DataFrame(rows)
    if df.empty:
        print("[FATAL] no tiles measured")
        return 1

    print(f"\n{'='*90}")
    print(f"WCSFIX veto bias -- {len(df)} tiles, {int(df.n_det.sum())} detections")
    print(f"refit displacement: median of per-tile medians "
          f"{df.shift_median_as.median():.3f}\"")
    print(f"{'='*90}")
    print(f"{'catalogue':<8} {'role':<28} {'raw':>8} {'refit':>8} {'delta':>8} "
          f"{'medsep raw':>11} {'refit':>8}")
    print("-" * 90)
    deltas = {}
    for ref, (_, role) in REFS.items():
        c_raw, c_ref = f"{ref}_raw_frac", f"{ref}_refit_frac"
        if c_raw not in df.columns:
            continue
        raw = 100 * df[c_raw].mean()
        ref_ = 100 * df[c_ref].mean()
        deltas[ref] = ref_ - raw
        print(f"{ref:<8} {role:<28} {raw:7.2f}% {ref_:7.2f}% {ref_-raw:+7.2f} "
              f"{df[f'{ref}_raw_medsep'].median():10.2f}\" "
              f"{df[f'{ref}_refit_medsep'].median():7.2f}\"")
    print("-" * 90)
    print(f"{'':<8} {'NULL control (chance rate)':<28}", end="")
    for ref in REFS:
        c = f"{ref}_raw_null_frac"
        if c in df.columns:
            print(f"  {ref}: {100*df[c].mean():.2f}% raw / "
                  f"{100*df[f'{ref}_refit_null_frac'].mean():.2f}% refit", end="")
    print()

    if "gaia" in deltas and ("ps1" in deltas or "usnob" in deltas):
        ctrl = [deltas[k] for k in ("ps1", "usnob") if k in deltas]
        ctrl_mean = float(np.mean(ctrl))
        excess = deltas["gaia"] - ctrl_mean
        print(f"\n{'-'*90}\nVERDICT")
        print(f"  Gaia match rate rises      {deltas['gaia']:+.2f} points")
        print(f"  unseen controls rise       {ctrl_mean:+.2f} points "
              f"(mean of {', '.join(f'{k} {deltas[k]:+.2f}' for k in ('ps1','usnob') if k in deltas)})")
        print(f"  GAIA-SPECIFIC EXCESS       {excess:+.2f} points")
        if abs(excess) < 0.5:
            print("  -> The refit is a genuine astrometric improvement. It raises")
            print("     matches against catalogues it never saw by the same amount,")
            print("     so the Gaia rise is real counterparts, not circularity.")
        elif excess > 0:
            print("  -> Gaia-specific inflation. The refit pulls toward Gaia beyond")
            print("     the astrometric improvement the controls confirm. The excess")
            print("     above is the size of the circularity and must be quoted.")
        else:
            print("  -> Gaia rises LESS than the controls; no Gaia-specific bias.")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, index=False)
        print(f"\n[OUT] {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
