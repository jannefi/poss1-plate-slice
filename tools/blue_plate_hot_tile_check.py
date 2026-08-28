#!/usr/bin/env python3
"""Do hot-tile S0 survivors show real flux on an INDEPENDENT second exposure?

[[hot_tiles_four_hypotheses_eliminated]] ruled out un-vetoed Gaia, displaced
catalogued stars, scan artifacts (survive SuperCOSMOS at 55.8%), and
USNO-B-catalogued sources for the released catalogue's hot tiles. The one
hypothesis left standing is a plate-emulsion defect -- but SuperCOSMOS
re-scans the SAME physical glass, so it can't separate "real feature of this
one piece of glass" from "real astronomical source". The POSS-I O (blue)
plate is a genuinely different exposure (different glass, same approximate
epoch) -- a defect specific to the E emulsion will not appear on it; a real
star or galaxy generally will, at least partially.

Method, reused from tools/v_flux_vs_null_on_plate.py (same repo family,
different repo): local-background z-score per position (box/aperture,
median/MAD background, off-array -> NaN), polarity resolved per-image
against bright Gaia stars rather than assumed. Run on BOTH bands per tile:
red is the positive control (S0 survivors are red-plate detections by
construction, so they must show real flux there), blue is the actual
question.

Fully S0-only. No V anywhere in this script -- public data (STScI cutouts,
Gaia for polarity resolution only) and our own released S0. Publishable
regardless of the V-publication question.

Usage:
    python3 tools/blue_plate_hot_tile_check.py \\
      --s0-csv results/s0-642-20260814/stage_S0.csv.gz \\
      --tile-manifest results/s0-642-20260814/tile_manifest.csv.gz \\
      --tiles tile_RA24.686_DECp33.529,tile_RA306.161_DECp63.564,... \\
      --out-dir work/blue_plate_hot_tile_check_20260828
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(REPO))


def local_z(data: np.ndarray, xs, ys, box: int, ap: int):
    """Local-background z-score of the pixel peak in an aperture at each (x,y).
    Identical method to v_flux_vs_null_on_plate.py::local_z (vasco60-repro-parity)."""
    from astropy.stats import mad_std

    ny, nx = data.shape
    z = np.full(len(xs), np.nan)
    peak = np.full(len(xs), np.nan)
    for i, (x, y) in enumerate(zip(xs, ys)):
        x, y = int(x), int(y)
        if x - box < 0 or x + box >= nx or y - box < 0 or y + box >= ny:
            continue
        win = data[y - box:y + box + 1, x - box:x + box + 1].copy()
        core = win[box - ap:box + ap + 1, box - ap:box + ap + 1].copy()
        win[box - ap:box + ap + 1, box - ap:box + ap + 1] = np.nan
        bg = win[np.isfinite(win)]
        if bg.size < 20:
            continue
        med, sig = np.median(bg), mad_std(bg)
        peak[i] = core.max()
        z[i] = (peak[i] - med) / max(sig, 1e-6)
    return z, peak


def resolve_polarity(data, w, gaia_cache, plate_ra, plate_dec, box, ap, n_stars=40):
    from vasco.local_cache_query import _cone_query

    g = _cone_query(gaia_cache, plate_ra, plate_dec, 30.0,
                    columns=["ra", "dec", "phot_g_mean_mag"])
    if not len(g):
        return 1, 0.0, 0
    g = g.sort_values("phot_g_mean_mag").head(n_stars)
    px, py = w.world_to_pixel_values(g["ra"].to_numpy(float), g["dec"].to_numpy(float))
    px, py = np.round(px).astype(int), np.round(py).astype(int)
    z, _ = local_z(data, px, py, box, ap)
    med = float(np.nanmedian(z)) if np.isfinite(z).any() else 0.0
    sign = 1 if med >= 0 else -1
    return sign, med, int(np.isfinite(z).sum())


def summarize(z):
    ok = np.isfinite(z)
    if not ok.sum():
        return {"n": 0}
    return {"n": int(ok.sum()), "median_z": round(float(np.nanmedian(z[ok])), 2),
            "frac_z_gt3": round(100 * float(np.mean(z[ok] > 3)), 1),
            "frac_z_gt5": round(100 * float(np.mean(z[ok] > 5)), 1)}


def process_tile(tile_id, tile_ra, tile_dec, rows, args, gaia_cache):
    from astropy.io import fits
    from astropy.wcs import WCS
    from vasco.downloader import fetch_skyview_dss

    out_fits = Path(args.out_dir) / "fits"
    out_fits.mkdir(parents=True, exist_ok=True)

    result = {"tile_id": tile_id, "tile_ra": tile_ra, "tile_dec": tile_dec,
              "n_survivors": len(rows)}
    bands_data = {}
    for band, survey in [("red", "dss1-red"), ("blue", "dss1-blue")]:
        p = fetch_skyview_dss(tile_ra, tile_dec, size_arcmin=args.size_arcmin,
                              survey=survey, out_dir=str(out_fits))
        hdr = fits.getheader(p)
        data = fits.getdata(p).astype(np.float64)
        w = WCS(hdr)
        nx, ny = int(hdr["NAXIS1"]), int(hdr["NAXIS2"])

        sign, star_z, n_star = resolve_polarity(data, w, gaia_cache, tile_ra, tile_dec,
                                                 args.box, args.ap)
        if sign < 0:
            data = -data

        px, py = w.world_to_pixel_values(rows.ra.to_numpy(float), rows.dec.to_numpy(float))
        pxr, pyr = np.round(px).astype(int), np.round(py).astype(int)
        on_array = (px >= 0) & (px < nx) & (py >= 0) & (py < ny)

        rng = np.random.default_rng(20260828)
        margin = args.box + 5
        n_null = max(len(rows), 100)
        null_x = rng.integers(margin, nx - margin, size=n_null)
        null_y = rng.integers(margin, ny - margin, size=n_null)

        sz, speak = local_z(data, pxr, pyr, args.box, args.ap)
        nz, npeak = local_z(data, null_x, null_y, args.box, args.ap)

        band_result = {
            "survey_header": str(hdr.get("SURVEY", "")),
            "on_array": int(on_array.sum()), "on_array_pct": round(100 * on_array.sum() / max(len(rows), 1), 1),
            "polarity_sign": sign, "polarity_star_median_z": round(star_z, 2), "polarity_n_stars": n_star,
            "survivors": summarize(sz), "null": summarize(nz),
        }
        result[band] = band_result
        bands_data[band] = (data, w, pxr, pyr, sz, on_array)
        print(f"  [{tile_id}][{band}] SURVEY={band_result['survey_header']} "
              f"on_array={band_result['on_array']}/{len(rows)} "
              f"survivors={band_result['survivors']} null={band_result['null']}")

    if args.render:
        render_tile(tile_id, bands_data, args)

    return result


def render_tile(tile_id, bands_data, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.visualization import ZScaleInterval

    out_png = Path(args.out_dir) / "png"
    out_png.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.6))
    zs = ZScaleInterval()
    for ax, band in zip(axes, ["red", "blue"]):
        data, w, px, py, z, on_array = bands_data[band]
        vmin, vmax = zs.get_limits(data)
        ax.imshow(data, cmap="gray", origin="lower", vmin=vmin, vmax=vmax)
        ok = on_array
        ax.scatter(px[ok], py[ok], s=18, facecolors="none",
                  edgecolors=np.where(z[ok] > 3, "lime", "red"), linewidths=1.0)
        ax.set_title(band, fontsize=13)
        ax.set_xticks([]); ax.set_yticks([])
    # One shared title + one shared legend line, not per-panel -- long
    # per-panel titles used to overlap where the two subplots meet.
    fig.suptitle(tile_id, fontsize=13, y=0.99)
    fig.text(0.5, 0.015, "green = z>3 (significant local flux)   "
             "red = z<=3 or off-array", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fig.savefig(out_png / f"{tile_id}.png", dpi=110)
    plt.close(fig)
    print(f"  [{tile_id}] wrote {out_png}/{tile_id}.png")


def main():
    warnings.filterwarnings("ignore")
    ap_ = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap_.add_argument("--s0-csv", required=True)
    ap_.add_argument("--tiles", required=True, help="comma-separated tile_id list")
    ap_.add_argument("--gaia-cache", default="/home/janne/local_cache/gaia")
    ap_.add_argument("--size-arcmin", type=float, default=60.0)
    ap_.add_argument("--box", type=int, default=15)
    ap_.add_argument("--ap", type=int, default=3)
    ap_.add_argument("--render", action="store_true", default=True)
    ap_.add_argument("--out-dir", default=str(REPO / "work" / "blue_plate_hot_tile_check_20260828"))
    args = ap_.parse_args()

    from vasco.utils.tile_id import parse_tile_id_center

    tiles = [t.strip() for t in args.tiles.split(",") if t.strip()]
    s0 = pd.read_csv(args.s0_csv)
    print(f"[CONFIG] {len(tiles)} tiles, S0={len(s0)} rows, box={args.box} ap={args.ap} "
          f"size={args.size_arcmin}'")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    results = []
    for tid in tiles:
        rows = s0[s0.tile_id == tid]
        if rows.empty:
            print(f"[SKIP] {tid}: 0 rows in S0")
            continue
        tra, tdec = parse_tile_id_center(tid)
        print(f"[TILE] {tid}  centre {tra:.4f} {tdec:.4f}  {len(rows)} survivors")
        r = process_tile(tid, tra, tdec, rows, args, args.gaia_cache)
        results.append(r)
        (out / "results.json").write_text(json.dumps(results, indent=2))

    # Aggregate
    def agg(band, pop):
        allz_n = sum(r[band][pop]["n"] for r in results if r[band][pop]["n"])
        if not allz_n:
            return {"n": 0}
        w_med = np.average([r[band][pop]["median_z"] for r in results if r[band][pop]["n"]],
                           weights=[r[band][pop]["n"] for r in results if r[band][pop]["n"]])
        w_gt3 = np.average([r[band][pop]["frac_z_gt3"] for r in results if r[band][pop]["n"]],
                           weights=[r[band][pop]["n"] for r in results if r[band][pop]["n"]])
        return {"n": allz_n, "weighted_median_z": round(float(w_med), 2),
                "weighted_frac_z_gt3": round(float(w_gt3), 1)}

    summary = {
        "n_tiles": len(results),
        "red_survivors": agg("red", "survivors"), "red_null": agg("red", "null"),
        "blue_survivors": agg("blue", "survivors"), "blue_null": agg("blue", "null"),
    }
    print("\n" + "=" * 70)
    print("[SUMMARY]")
    print(json.dumps(summary, indent=2))
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[out] {out}/results.json, summary.json"
         + (f", png/*.png" if args.render else ""))


if __name__ == "__main__":
    main()
