#!/usr/bin/env python3
"""Is the ~2.25" full-plate slice offset in the PIXELS or in the WCS?

The whole-sky union test showed the full-plate arm stepping from 62.76% at 2" to
94.85% at 3" -- a systematic. Measured per plate it is sharply bimodal: either
0.017" or ~2.25", with 25 of 66 sampled plates in the bad mode, and against Gaia
it is our slice that is wrong, not the archive cutout (XE291: 2.075" vs 0.378").

Two candidate causes, with opposite remedies:

  (A) WCS  -- the images agree pixel-for-pixel, but the plate solution as astropy
      evaluates it disagrees with the one STScI applied to its own copy. A bug on
      our side, fixable in slice_plate_tiles.py, and no calibration table needed.

  (B) PIXELS -- IRSA's full-plate scan is registered differently from what STScI
      serves, so the same sky sits at a different place in the array. Nothing in
      the slicer can fix that, and a per-plate Gaia correction becomes the right
      answer rather than a fudge factor.

The test separates them. Take an archive cutout, find the matching region of the
full plate BY WCS, and cross-correlate the pixels:

  - offset ~0 px  -> images agree, so the disagreement is in the WCS      => (A)
  - offset ~1.3 px in the direction of the sky offset -> scans differ     => (B)

1.3 px is what 2.25" comes to at 1.7"/px. Always run a known-good plate as a
control in the same invocation: a method that reports "no shift" on a bad plate
is only meaningful if it reports "no shift" on a good one too, and the same
comparison on XE524 previously gave correlation 0.999998.

How to validate
---------------
    python3 tools/diagnose_slice_offset.py \
        --bad XE162 --bad XE140 --control XE524

Reads only; writes nothing outside the scratch directory.
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from vasco.paths import get as _p
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PLATE_DIR = str(_p("plate_dir"))
TILES_ROOT = str(_p("tiles_dir"))


def archive_tile_for(plate: str, t2p: pd.DataFrame) -> tuple[str, Path] | None:
    for _, r in t2p[t2p.plate_id == plate].iterrows():
        p = Path(TILES_ROOT) / str(r.tile_id) / "raw" / str(r.tile_fits)
        if p.exists():
            return str(r.tile_id), p
    return None


def phase_shift(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """Sub-pixel (dy, dx) shift taking a onto b, plus peak correlation."""
    a = a.astype(float) - np.median(a)
    b = b.astype(float) - np.median(b)
    # Taper to stop edge discontinuities from dominating the transform.
    wy = np.hanning(a.shape[0])[:, None]
    wx = np.hanning(a.shape[1])[None, :]
    A = np.fft.fft2(a * wy * wx)
    B = np.fft.fft2(b * wy * wx)
    R = A * np.conj(B)
    R /= np.maximum(np.abs(R), 1e-30)
    c = np.fft.ifft2(R).real
    peak = float(c.max() / max(c.std(), 1e-30))
    iy, ix = np.unravel_index(np.argmax(c), c.shape)

    def refine(idx, axis_len, line):
        m = line[idx]
        l = line[(idx - 1) % axis_len]
        r = line[(idx + 1) % axis_len]
        d = 2 * m - l - r
        return idx + (0.5 * (r - l) / d if abs(d) > 1e-30 else 0.0)

    fy = refine(iy, c.shape[0], c[:, ix])
    fx = refine(ix, c.shape[1], c[iy, :])
    if fy > c.shape[0] / 2:
        fy -= c.shape[0]
    if fx > c.shape[1] / 2:
        fx -= c.shape[1]
    return fy, fx, peak


def run(plate: str, label: str, t2p: pd.DataFrame) -> dict | None:
    from astropy.io import fits
    from astropy.nddata import Cutout2D
    from astropy.wcs import WCS

    got = archive_tile_for(plate, t2p)
    if not got:
        print(f"  [{plate}] no archive tile FITS on disk -- skipped")
        return None
    tile_id, tile_path = got

    th = fits.open(tile_path)
    tdata, thdr = th[0].data, th[0].header
    tw = WCS(thdr)

    ph = fits.open(f"{PLATE_DIR}/dss1red_{plate}.fits", memmap=False)[0]
    pw = WCS(ph.header)
    scale = ph.header["XPIXELSZ"] / 1000.0 * ph.header["PLTSCALE"]

    ny, nx = tdata.shape
    # Centre of the archive tile, in sky, then in full-plate pixels.
    ra_c, dec_c = [float(v) for v in tw.pixel_to_world_values(nx / 2.0, ny / 2.0)]
    px, py = [float(v) for v in pw.world_to_pixel_values(ra_c, dec_c)]

    side = min(ny, nx, 1024)
    try:
        cut = Cutout2D(ph.data, position=(px, py), size=(side, side), wcs=pw, copy=True)
    except Exception as e:
        print(f"  [{plate}] cutout failed: {e}")
        return None
    tcut = Cutout2D(tdata, position=(nx / 2.0, ny / 2.0), size=(side, side), copy=True)

    dy, dx, peak = phase_shift(cut.data, tcut.data)
    shift_px = float(np.hypot(dy, dx))
    print(f"  [{plate}] {label:<8} tile {tile_id}")
    print(f"      pixel shift  dx={dx:+.3f}  dy={dy:+.3f}  |shift|={shift_px:.3f} px "
          f"= {shift_px * scale:.3f}\"   (corr peak {peak:.1f} sigma)")
    return {"plate": plate, "label": label, "dx_px": dx, "dy_px": dy,
            "shift_px": shift_px, "shift_arcsec": shift_px * scale, "peak_sigma": peak}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bad", action="append", default=[], help="Plate in the ~2.25\" mode")
    ap.add_argument("--control", action="append", default=[], help="Known-good plate")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if not args.bad or not args.control:
        raise SystemExit("[FATAL] pass at least one --bad and one --control; "
                         "a result on bad plates alone proves nothing")

    t2p = pd.read_csv("data/metadata/tile_to_plate.csv")
    t2p = t2p[t2p.plate_id.notna() & t2p.tile_fits.notna()]

    print("Cross-correlating full-plate slice against archive cutout of the same sky.")
    print("  ~0 px      -> images agree, disagreement is in the WCS   (our bug)")
    print("  ~1.3 px    -> scans registered differently               (needs calibration)\n")
    rows = []
    for p in args.bad:
        r = run(p, "BAD", t2p)
        if r:
            rows.append(r)
    for p in args.control:
        r = run(p, "CONTROL", t2p)
        if r:
            rows.append(r)

    d = pd.DataFrame(rows)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        d.to_csv(args.out, index=False)
        print(f"\nwrote {args.out}")
    if len(d) and d.label.nunique() > 1:
        bad = d[d.label == "BAD"].shift_px.median()
        ctl = d[d.label == "CONTROL"].shift_px.median()
        print(f"\nmedian |shift|:  BAD {bad:.3f} px   CONTROL {ctl:.3f} px")
        print("VERDICT:", "PIXELS differ -- scans are registered differently"
              if bad - ctl > 0.5 else
              "PIXELS agree -- the offset is in the WCS, and it is ours to fix")
    return 0


if __name__ == "__main__":
    sys.exit(main())
