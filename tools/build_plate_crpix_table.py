#!/usr/bin/env python3
"""Per-plate CRPIX correction, derived from public headers alone.

Why this exists
---------------
STScI's DSS headers carry TWO independent astrometric solutions:

  1. the GSSS plate polynomial  (AMDX/AMDY, PPO1-6, CNPIX1/2, XPIXELSZ, PLTSCALE)
  2. an explicit FITS-standard  CRPIX / CRVAL / CD set

**astropy/wcslib prefers (1). SExtractor reads only (2).** On 207 of 633 POSS-I
plates the two disagree by ~2.34"; on the remaining 426 they agree to 0.125".
The explicit set is the correct one -- it is what matches Gaia, at ~0.06".

IRSA's full-plate scans carry ONLY solution (1), so tiles sliced from them
inherit the wrong answer on a third of plates. This table restores solution (2)'s
convention by measuring the offset between them, per plate, and the slicer
applies it as a CRPIX shift.

See docs/DSS_WCS_TWO_SOLUTIONS.md for the full diagnosis.

Why it needs no catalogue and no fitting
----------------------------------------
The offset is a property of the two header solutions, not of the sky. It is pure
header arithmetic: take an archive cutout that carries both solutions, ask each
where the same data pixel points, and difference them. No Gaia, no reference
catalogue, no fitted parameter, no threshold to tune -- so the table is exactly
reproducible by anyone with the same public headers, and auditable line by line.

An earlier version instead fitted a discrete -1/0 flag against Gaia. That worked
(it agrees with this table on 99.0% of plates) but it had to be *measured*
against a catalogue and could in principle be classified wrong. This cannot.

Sign convention -- get this wrong and you double the error
----------------------------------------------------------
`delta` is measured as (pixel where the GSSS solution places a sky position)
minus (pixel where that data actually sits, per CNPIX bookkeeping). A WCS whose
CRPIX is reduced by d evaluates at pixel p exactly as the original would at p+d.
So the correction is:

    CRPIX_corrected = CRPIX - delta

Applying +delta instead moves the tile the wrong way and roughly doubles the
residual rather than removing it.

Accuracy
--------
The offset is near-constant across a plate, but not perfectly: ~0.1 px of real
scatter remains (~0.17"), so this is a large improvement on the uncorrected
2.34" and a modest one on the discrete flag's 0.282". Do not expect 0.000".

How to validate
---------------
    python3 tools/build_plate_crpix_table.py \\
        --plate-dir <dir of dss1red_XE*.fits> \\
        --archive-tiles <dir of archive cutout FITS> \\
        --out data/plate_crpix_table.csv

Flagged plates should come out near |delta| ~ 1.0 px and unflagged near 0.0, with
`scatter_px` well under 0.2 for both. A plate whose scatter is large is not
trustworthy and is reported rather than silently used.
"""
from __future__ import annotations

import argparse
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# Keywords that make astropy prefer the GSSS solution. Removing them all is what
# leaves the explicit CRPIX/CD set standing on its own.
DSS_KEYS = (
    ["CNPIX1", "CNPIX2", "PPO1", "PPO2", "PPO3", "PPO4", "PPO5", "PPO6",
     "XPIXELSZ", "YPIXELSZ", "PLTSCALE", "PLTRAH", "PLTRAM", "PLTRAS",
     "PLTDECSN", "PLTDECD", "PLTDECM", "PLTDECS"]
    + [f"AMD{a}{i}" for a in ("X", "Y") for i in range(1, 21)]
    + [f"AMDRE{a}{i}" for a in ("X", "Y") for i in range(1, 21)]
)


def index_archive_tiles(root: Path) -> dict[str, list[Path]]:
    """Map REGION -> cutout files that carry BOTH solutions."""
    from astropy.io import fits

    by_plate: dict[str, list[Path]] = defaultdict(list)
    for p in sorted(root.rglob("*.fits")):
        try:
            h = fits.getheader(p)
        except Exception:
            continue
        if "CRPIX1" in h and "AMDX1" in h and h.get("REGION"):
            by_plate[str(h["REGION"]).strip()].append(p)
    return by_plate


def measure(plate: str, plate_fits: Path, tiles: list[Path], max_tiles: int) -> dict:
    from astropy.io import fits
    from astropy.wcs import WCS

    r = {"plate": plate, "delta_x_px": np.nan, "delta_y_px": np.nan,
         "scatter_px": np.nan, "offset_arcsec": np.nan, "n_tiles": 0, "status": ""}
    if not plate_fits.exists():
        r["status"] = "no_plate_fits"
        return r
    if not tiles:
        r["status"] = "no_archive_tile"
        return r

    ph = fits.getheader(plate_fits)
    pw = WCS(ph)                                   # GSSS solution (what we must fix)
    scale = ph["XPIXELSZ"] / 1000.0 * ph["PLTSCALE"]

    dxs, dys = [], []
    for tp in tiles[:max_tiles]:
        th = fits.getheader(tp)
        hs = th.copy()
        for k in DSS_KEYS:
            hs.remove(k, ignore_missing=True)
        try:
            ws = WCS(hs)                           # explicit CRPIX/CD -- the truth
        except Exception:
            continue
        nx, ny = int(th["NAXIS1"]), int(th["NAXIS2"])
        g = np.linspace(0, nx - 1, 5)
        x, y = [a.ravel() for a in np.meshgrid(g, np.linspace(0, ny - 1, 5))]

        ra, dec = ws.all_pix2world(x, y, 1)        # FITS 1-based, as SExtractor reads
        px, py = pw.world_to_pixel_values(ra, dec)
        # Where that data really sits on the plate, from CNPIX bookkeeping (0-based).
        tx = x - 1 + th["CNPIX1"] - ph["CNPIX1"]
        ty = y - 1 + th["CNPIX2"] - ph["CNPIX2"]
        dxs.append(np.median(px - tx))
        dys.append(np.median(py - ty))

    if not dxs:
        r["status"] = "no_usable_tile"
        return r
    dx, dy = float(np.median(dxs)), float(np.median(dys))
    r.update(delta_x_px=dx, delta_y_px=dy, n_tiles=len(dxs),
             scatter_px=float(np.hypot(np.std(dxs), np.std(dys))),
             offset_arcsec=float(np.hypot(dx, dy) * scale), status="ok")
    return r


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plate-dir", required=True)
    ap.add_argument("--archive-tiles", required=True,
                    help="Directory tree of archive cutout FITS carrying both solutions.")
    ap.add_argument("--out", default="data/plate_crpix_table.csv")
    ap.add_argument("--max-tiles", type=int, default=4,
                    help="Cutouts averaged per plate. More costs time, buys little.")
    args = ap.parse_args()

    plate_dir = Path(args.plate_dir)
    print(f"[SCAN] indexing archive cutouts under {args.archive_tiles} ...", flush=True)
    by_plate = index_archive_tiles(Path(args.archive_tiles))
    print(f"[SCAN] {len(by_plate)} plates have a cutout carrying both solutions")

    plates = sorted({p.stem.replace("dss1red_", "")
                     for p in plate_dir.glob("dss1red_XE*.fits")})
    rows = [measure(pl, plate_dir / f"dss1red_{pl}.fits", by_plate.get(pl, []),
                    args.max_tiles) for pl in plates]

    d = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(args.out, index=False)

    ok = d[d.status == "ok"]
    big = ok[ok.offset_arcsec > 1.5]
    print(f"\n=== {len(d)} plates: {d.status.value_counts().to_dict()} ===")
    print(f"  needing correction (>1.5\"): {len(big)} / {len(ok)} "
          f"({100 * len(big) / max(len(ok), 1):.1f}%)")
    if len(big):
        print(f"  their median offset        : {big.offset_arcsec.median():.3f}\"")
        print(f"  their median |delta|       : "
              f"{np.hypot(big.delta_x_px, big.delta_y_px).median():.3f} px")
    rest = ok[ok.offset_arcsec <= 1.5]
    if len(rest):
        print(f"  the rest, median offset    : {rest.offset_arcsec.median():.3f}\"")
    shaky = ok[ok.scatter_px > 0.2]
    print(f"  high scatter (>0.2 px, not trustworthy): {len(shaky)}")
    if len(shaky):
        print(shaky[["plate", "delta_x_px", "delta_y_px", "scatter_px", "n_tiles"]]
              .to_string(index=False))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
