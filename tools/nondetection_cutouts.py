#!/usr/bin/env python3
"""Side-by-side cutouts for reference rows we never detected.

`archive_slice_parity.py` finds a handful of reference S0 rows with no detection
of ours within 30". Statistics cannot say whether those are real sources we
missed, sources on pixels we never saw, or spurious rows in the reference
catalogue. Pixels can.

Each row gets two stamps of the same sky at the same scale:

  left   the reference run's own tile, as served by the STScI cutout service
  right  our plate scan, cut from the IRSA full-plate FITS we slice from

Overlaid: the reference position (circle) and every one of our own raw
detections inside the stamp (crosses), so "we detected nothing here" is visible
rather than asserted.

The FITS `REGION` keyword of each side is printed, because plate identity is the
first thing to rule out — a plan-derived `plate_id` is a label, whereas REGION is
what the archive actually served. If the two sides disagree, the miss is plate
selection and not detection.

DSS scans carry no photometric calibration, so the stretch is per-stamp zscale
and brightness is not comparable between panels.

How to validate
---------------
    python3 tools/nondetection_cutouts.py \
        --rows work/archive_slice_parity/never_detected_58.csv \
        --bucket D \
        --ref-tiles <vasco60>/data/tiles \
        --plate-dir <plate-scan-dir> \
        --radec-dir <run>/radec \
        --out work/archive_slice_parity/nondetection_cutouts.png

Sanity: both panels must show the same star field. If they do not, the WCS or
the plate is wrong, and nothing else in the figure means anything.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astropy.coordinates import SkyCoord
import astropy.units as u


def stamp(path: str, ra: float, dec: float, half_px: int):
    """Read a square subarray centred on (ra, dec). Returns (data, wcs, header)."""
    with fits.open(path, memmap=True) as hdul:
        hdu = hdul[0]
        w = WCS(hdu.header)
        try:
            x, y = w.all_world2pix(ra, dec, 0)
        except Exception:
            x, y = w.wcs_world2pix(ra, dec, 0)
        x, y = int(round(float(x))), int(round(float(y)))
        ny, nx = hdu.shape
        x0, x1 = max(0, x - half_px), min(nx, x + half_px)
        y0, y1 = max(0, y - half_px), min(ny, y + half_px)
        if x1 <= x0 or y1 <= y0:
            return None, None, dict(hdu.header)
        data = np.asarray(hdu.section[y0:y1, x0:x1], dtype=float)
        sub = w.deepcopy()
        sub = sub[y0:y1, x0:x1]
        return data, sub, dict(hdu.header)


def show(ax, data, wcs, ra, dec, dets, title, sub):
    if data is None or data.size == 0:
        ax.text(0.5, 0.5, "no pixel coverage", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title, fontsize=8)
        return
    lo, hi = ZScaleInterval().get_limits(data)
    ax.imshow(data, origin="lower", cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
    try:
        cx, cy = wcs.all_world2pix(ra, dec, 0)
    except Exception:
        cx, cy = wcs.wcs_world2pix(ra, dec, 0)
    ax.add_patch(plt.Circle((float(cx), float(cy)), 6, fill=False, color="#ff3b30", lw=1.4))
    if dets is not None and len(dets):
        try:
            dx, dy = wcs.all_world2pix(dets["ra"].to_numpy(), dets["dec"].to_numpy(), 0)
            ok = (dx > 0) & (dx < data.shape[1]) & (dy > 0) & (dy < data.shape[0])
            ax.plot(dx[ok], dy[ok], "+", color="#32d74b", ms=9, mew=1.3, ls="none")
        except Exception:
            pass
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=8)
    ax.text(0.02, 0.02, sub, transform=ax.transAxes, fontsize=6.5, color="#ffd60a",
            va="bottom", ha="left")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--bucket", default="D", help="leading letter of the bucket column to render")
    ap.add_argument("--ref-tiles", required=True)
    ap.add_argument("--plate-dir", required=True)
    ap.add_argument("--radec-dir", required=True)
    ap.add_argument("--sex-csv", default="work/archive_slice_parity/rows_with_published_sex.csv")
    ap.add_argument("--box-arcmin", type=float, default=2.5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.rows)
    if "bucket" in df.columns:
        df = df[df["bucket"].astype(str).str.startswith(args.bucket)]
    df = df.reset_index(drop=True)
    if df.empty:
        print("[FATAL] no rows selected", file=sys.stderr)
        return 2
    print(f"rendering {len(df)} rows (bucket {args.bucket})")

    sex = None
    if Path(args.sex_csv).exists():
        sex = pd.read_csv(args.sex_csv, usecols=lambda c: c in
                          ("src_id", "SNR_WIN", "MAG_AUTO", "SPREAD_MODEL", "ELONGATION", "FLAGS"))

    fig, axes = plt.subplots(len(df), 2, figsize=(7.2, 3.5 * len(df)))
    axes = np.atleast_2d(axes)
    report = []

    for i, r in df.iterrows():
        ra, dec = float(r["ra"]), float(r["dec"])
        plate = str(r["our_nearest_plate"])

        # our own raw detections inside a generous box around the position
        det = pd.read_csv(Path(args.radec_dir) / f"{plate}.csv", usecols=["ra", "dec"])
        m = ((np.abs(det["dec"] - dec) < 0.05) &
             (np.abs(((det["ra"] - ra + 180) % 360) - 180) * np.cos(np.radians(dec)) < 0.05))
        det = det[m]

        # ---- reference side: the STScI-served tile ----
        cand = glob.glob(str(Path(args.ref_tiles) / str(r["ref_tile_id"]) / "raw" / "*.fits"))
        ref_region = ref_date = ref_plt = "?"
        if cand:
            # REGION lives in the FITS, not in header.json's "selected" block.
            rh = fits.getheader(cand[0])
            ref_region = str(rh.get("REGION", "?")).strip()
            ref_plt = str(rh.get("PLTLABEL", "?")).strip()
            ref_date = str(rh.get("DATE-OBS", "?"))[:10]
            scale = abs(rh.get("CD2_2") or rh.get("CDELT2") or 4.7e-4) * 3600
            half = int(args.box_arcmin * 60 / 2 / max(scale, 0.1))
            d1, w1, _ = stamp(cand[0], ra, dec, half)
        else:
            d1 = w1 = None

        # ---- our side: the IRSA plate scan we slice from ----
        pf = Path(args.plate_dir) / f"dss1red_{plate}.fits"
        our_region = "?"
        our_plt = "?"
        if pf.exists():
            ph = fits.getheader(pf)
            our_region = str(ph.get("REGION", "?")).strip()
            our_plt = str(ph.get("PLTLABEL", "?")).strip()
            d2, w2, _ = stamp(str(pf), ra, dec, int(args.box_arcmin * 60 / 2 / 1.7))
        else:
            d2 = w2 = None

        ann = ""
        if sex is not None:
            row = sex[sex["src_id"] == r["src_id"]]
            if len(row):
                q = row.iloc[0]
                ann = (f"SNR {q.get('SNR_WIN', float('nan')):.0f}  "
                       f"mag {q.get('MAG_AUTO', float('nan')):.2f}  "
                       f"elong {q.get('ELONGATION', float('nan')):.2f}")

        show(axes[i, 0], d1, w1, ra, dec, det,
             f"REFERENCE (STScI)  {r['ref_tile_id']}\nREGION {ref_region}  {ref_date}", ann)
        show(axes[i, 1], d2, w2, ra, dec, det,
             f"OURS (IRSA plate slice)  plate {plate}\nREGION {our_region}   "
             f"nearest our detection {r['dist_raw_arcsec']:.1f}\"", "")
        report.append(dict(src_id=r["src_id"], ra=ra, dec=dec, ref_plate=r["ref_plate"],
                           our_plate=plate, ref_REGION=ref_region, our_REGION=our_region,
                           ref_PLTLABEL=ref_plt, our_PLTLABEL=our_plt,
                           region_agrees=(ref_region == our_region),
                           dist_raw_arcsec=r["dist_raw_arcsec"],
                           our_dets_in_box=int(len(det))))

    fig.suptitle(f"Reference S0 rows with no detection of ours within 30\"  "
                 f"({args.box_arcmin:g}' stamps; red = reference position, green + = our raw detections)",
                 fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")

    rep = pd.DataFrame(report)
    csv = str(Path(args.out).with_suffix(".csv"))
    rep.to_csv(csv, index=False)
    print("\n=== plate identity: what each side actually served ===")
    print(rep[["src_id", "ref_REGION", "ref_PLTLABEL", "our_REGION", "our_PLTLABEL",
               "region_agrees", "dist_raw_arcsec", "our_dets_in_box"]].to_string(index=False))
    print(f"\nREGION agrees on {int(rep.region_agrees.sum())}/{len(rep)} rows")
    print(f"wrote {csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
