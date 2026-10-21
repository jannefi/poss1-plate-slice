#!/usr/bin/env python3
"""Which plate does a catalogue's rows sit on, when more than one plate covers them?

POSS-I plates overlap. On overlapping sky a source can be measured on any of the
plates that cover it, and which one a catalogue used is a property of how that
catalogue was made rather than of the sky. This tool measures it: for every row,
it finds every plate whose array covers the position, computes the distance to
the nearest array edge on each, and asks whether the row's own plate is the one
that puts it furthest from any edge.

The answer separates two designs cleanly. A tessellation assigns a position to
whichever tile it falls in, so the plate it lands on is whatever the grid
happened to choose -- the furthest-from-edge plate only by luck. A per-position
request to a cutout service asks for a named position and gets back the plate
the service picks, and STScI's `getimage` picks FURTHEST_FROM_EDGE by default
(see docs/SCAN_SOURCE_SENSITIVITY.md).

WHAT THIS DOES NOT TELL YOU
---------------------------
It measures the geometry of a finished catalogue. It does not show how any
catalogue was produced, who ran what, or in which order. A high rate is
consistent with per-position plate selection and inconsistent with a
position-blind grid; it is not a reconstruction of anyone's pipeline, and this
tool does not attempt one.

INPUT
-----
`--ref-csv` is a LOCAL path and nothing is downloaded. Reference catalogues are
other people's data under their own terms; this repository ships none and
fetches none. SVO-style exports are read as-is, leading `#` included.

`--footprints` defaults to the table shipped here: all 932 POSS-I red plates in
the IRSA archive, carrying the plate-solution keywords from the public FITS
headers. A reader's astropy rebuilds the same WCS this tool used, so the result
does not depend on our arithmetic -- and it needs 168 KB rather than a few
hundred gigabytes of plate scans.

Note the 932 against the 642 plates this project processed: coverage is
computed over every archive plate, including ones we never ran, because a
position can be covered by a plate nobody chose to process.

How to validate
---------------
Run it on this repository's own released catalogue first. A tessellation should
score low; if it scores high, the footprint table or the matching is wrong:

    python3 tools/plate_choice_test.py \\
        --ref-csv results/s0-642-paper-parity-20260828/stage_S0.csv.gz \\
        --manifest results/s0-642-paper-parity-20260828/tile_manifest.csv.gz

That reports 36.3% -- near what a grid should give, and the control that tells
you the tool is measuring what it claims.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
from astropy.io.fits import Header
from astropy.wcs import WCS

RA_NAMES = ("ra", "ra_deg", "raj2000", "_raj2000", "alpha", "alpha_j2000")
DEC_NAMES = ("dec", "de", "dec_deg", "dej2000", "_dej2000", "delta", "delta_j2000")
PLATE_NAMES = ("name", "plate", "plate_id", "det_plate", "field")
META = ("plate_id", "center_ra", "center_dec", "as_per_px")


def read_table(path: Path) -> pd.DataFrame:
    """Read a CSV, tolerating gzip and the leading '#' SVO puts on its header."""
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt") as f:
        text = f.read()
    lines = text.splitlines()
    if lines and lines[0].lstrip().startswith("#"):
        lines[0] = lines[0].lstrip().lstrip("#").strip()
    return pd.read_csv(io.StringIO("\n".join(lines)))


def find_col(df, candidates):
    lower = {c.lower().strip(): c for c in df.columns}
    for c in candidates:
        if c in lower:
            return lower[c]
    return None


def load_plates(path: Path) -> dict:
    fp = read_table(path)
    keys = [c for c in fp.columns if c not in META]
    out = {}
    for row in fp.to_dict("records"):
        h = Header()
        for k in keys:
            h[k.upper()] = row[k]
        out[row["plate_id"]] = {
            "wcs": WCS(h, relax=True),
            "nax1": int(row["naxis1"]), "nax2": int(row["naxis2"]),
            "as_per_px": float(row["as_per_px"]),
            "cra": float(row["center_ra"]), "cdec": float(row["center_dec"]),
        }
    return out


def _edge_px(pl, ra, dec):
    px, py = pl["wcs"].world_to_pixel_values(ra, dec)
    d = np.minimum(np.minimum(px, pl["nax1"] - px),
                   np.minimum(py, pl["nax2"] - py))
    return d * pl["as_per_px"] / 60.0


def edge_on_own_plate(ra, dec, pid, plates):
    out = np.full(len(ra), np.nan)
    for p in pd.unique(pid):
        if p not in plates:
            continue
        ii = np.flatnonzero(pid == p)
        out[ii] = _edge_px(plates[p], ra[ii], dec[ii])
    return out


def coverage(ra, dec, plates, chunk_deg=5.0):
    """Number of covering plates, and the best edge distance any of them offers."""
    n = len(ra)
    ncov = np.zeros(n, np.int32)
    best = np.full(n, -np.inf)
    cra, cdec = np.radians(ra), np.radians(dec)
    for pl in plates.values():
        sep = np.degrees(np.arccos(np.clip(
            np.sin(cdec) * np.sin(np.radians(pl["cdec"])) +
            np.cos(cdec) * np.cos(np.radians(pl["cdec"])) *
            np.cos(cra - np.radians(pl["cra"])), -1, 1)))
        ii = np.flatnonzero(sep < chunk_deg)
        if not len(ii):
            continue
        e = _edge_px(pl, ra[ii], dec[ii])
        k = e > 0
        jj = ii[k]
        ncov[jj] += 1
        best[jj] = np.maximum(best[jj], e[k])
    best[~np.isfinite(best)] = np.nan
    return ncov, best


def main() -> int:
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ref-csv", required=True, type=Path,
                    help="LOCAL path to the catalogue under test. Never downloaded.")
    ap.add_argument("--manifest", type=Path, default=None,
                    help="tile_id,plate_id map, when the catalogue names tiles "
                         "rather than plates")
    ap.add_argument("--footprints", type=Path,
                    default=here / "data" / "plate_footprints.csv.gz")
    ap.add_argument("--tolerance-arcsec", type=float, default=0.6,
                    help="how close to the best edge distance still counts as "
                         "'the furthest-from-edge plate' (default 0.6\")")
    ap.add_argument("--sweep", action="store_true",
                    help="report the rate across a range of tolerances")
    ap.add_argument("--out-summary", type=Path, default=None)
    args = ap.parse_args()

    plates = load_plates(args.footprints)
    print(f"[PLATES] {len(plates)} footprints from {args.footprints}")

    df = read_table(args.ref_csv)
    rac, decc = find_col(df, RA_NAMES), find_col(df, DEC_NAMES)
    if not rac or not decc:
        print(f"[FATAL] no RA/Dec columns in {args.ref_csv} (saw {list(df.columns)})")
        return 2

    pid_col = find_col(df, PLATE_NAMES)
    if args.manifest is not None:
        man = read_table(args.manifest)
        tcol, pcol = find_col(man, ("tile_id",)), find_col(man, ("plate_id",))
        tref = find_col(df, ("tile_id",))
        if not (tcol and pcol and tref):
            print("[FATAL] --manifest needs tile_id/plate_id, and the catalogue a tile_id")
            return 2
        df = df.assign(_plate=df[tref].map(man.set_index(tcol)[pcol]))
        pid_col = "_plate"
    if not pid_col:
        print(f"[FATAL] no plate column in {args.ref_csv}; pass --manifest if it "
              f"names tiles (saw {list(df.columns)})")
        return 2

    df = df[df[pid_col].notna()]
    ra = pd.to_numeric(df[rac], errors="coerce").values
    dec = pd.to_numeric(df[decc], errors="coerce").values
    pid = df[pid_col].astype(str).values
    good = np.isfinite(ra) & np.isfinite(dec)
    ra, dec, pid = ra[good], dec[good], pid[good]
    print(f"[REF]    {args.ref_csv}  rows={len(ra)}  coords={rac}/{decc}  plate={pid_col}")

    ncov, best = coverage(ra, dec, plates)
    edge = edge_on_own_plate(ra, dec, pid, plates)
    m = (ncov > 1) & np.isfinite(best) & np.isfinite(edge)
    short = (best - edge)[m]
    tol = args.tolerance_arcsec / 60.0
    rate = float(np.mean(short < tol))

    print(f"\n  covered by >1 plate : {int(m.sum())} of {len(ra)} "
          f"({100*m.mean():.1f}%)")
    print(f"  own plate IS the furthest-from-edge plate: {100*rate:.1f}%"
          f"   (tolerance {args.tolerance_arcsec:g}\")")
    rest = short[short >= tol]
    if len(rest):
        print(f"  median shortfall when it is not: {np.median(rest):.2f}'")
    print(f"  shortfall percentiles: p50 {np.percentile(short,50):.3f}'  "
          f"p90 {np.percentile(short,90):.3f}'  p99 {np.percentile(short,99):.2f}'")

    summary = {"reference": str(args.ref_csv), "rows": int(len(ra)),
               "multi_coverage": int(m.sum()),
               "tolerance_arcsec": args.tolerance_arcsec,
               "furthest_from_edge_fraction": round(rate, 5)}

    if args.sweep:
        print("\n  tolerance sweep:")
        summary["sweep"] = {}
        for t_as in (0.3, 0.6, 1.2, 2.4, 5.0, 10.0, 30.0):
            r = float(np.mean(short < t_as / 60.0))
            summary["sweep"][f"{t_as:g}"] = round(r, 5)
            print(f"    {t_as:>6.1f}\"   {100*r:6.1f}%")
        print("\n  The two DSS archives' plate solutions disagree by up to ~2.4\"")
        print("  (docs/DSS_WCS_TWO_SOLUTIONS.md). If the rate is flat across that")
        print("  band, the astrometric choice is not what produced it.")

    if args.out_summary:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        args.out_summary.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"\n[WROTE] {args.out_summary}")

    print("\n  A high rate is consistent with per-position plate selection and")
    print("  inconsistent with a position-blind grid. It is not a reconstruction")
    print("  of anyone's pipeline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
