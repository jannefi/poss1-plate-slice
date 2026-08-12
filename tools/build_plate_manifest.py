#!/usr/bin/env python3
"""Generate the plate manifest from plate headers alone.

The footprint this project processes is defined by a rule, not by a list:

    every POSS-I red plate whose centre declination is >= -3.0 deg

Applied to IRSA's dss1red library that selects 642 plates, XE002..XE643.

Why a rule and not a list: a list has to come from somewhere, and a list whose
origin cannot be stated cannot support a reproducibility claim. Anyone with the
IRSA scans can re-run this script and get the same 642 plates without access to
any private catalogue.

Note the limit applies to plate CENTRES while the survey limit it derives from
concerns data coverage. That is an interpretive step, and it is harmless here:
the northern and southern plate sets are separated by a 4.4 deg gap (lowest
northern centre -0.81, highest southern -5.19), so every threshold in that gap
selects exactly the same 642 plates. The script asserts this, so a future edit
to DEC_MIN cannot silently change the footprint.

Usage:
    python tools/build_plate_manifest.py --plate-dir /path/to/dss1red \\
        --out data/plate_manifest.csv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
from astropy.io import fits  # noqa: E402

DEC_MIN = -3.0
# The empty band between the northern and southern plate sets. Any threshold
# inside it yields the same footprint; if a future library breaks that, the
# footprint is no longer threshold-insensitive and the claim above is void.
GAP = (-5.0, -1.0)
PLATE_RX = re.compile(r"dss1red_(XE\d+)\.fits$")


def plate_centre_dec(header) -> float:
    """Plate centre declination from the GSSS keywords, in degrees.

    PLTDECSN carries the sign as a separate '+'/'-' field, so the sign must be
    applied to the assembled magnitude -- reading PLTDECD alone silently loses
    it for southern plates.
    """
    sign = -1.0 if str(header.get("PLTDECSN", "+")).strip() == "-" else 1.0
    return sign * (
        abs(header.get("PLTDECD", 0))
        + header.get("PLTDECM", 0) / 60.0
        + header.get("PLTDECS", 0) / 3600.0
    )


def scan(plate_dir: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for f in sorted(plate_dir.glob("dss1red_XE*.fits")):
        m = PLATE_RX.search(f.name)
        if not m:
            continue
        out[m.group(1)] = plate_centre_dec(fits.getheader(f))
    if not out:
        sys.exit(f"[FAIL] no dss1red_XE*.fits found under {plate_dir}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plate-dir", required=True, type=Path,
                    help="directory of IRSA dss1red_XE*.fits scans")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--dec-min", type=float, default=DEC_MIN)
    args = ap.parse_args()

    decs = scan(args.plate_dir)
    selected = sorted(p for p, d in decs.items() if d >= args.dec_min)

    # Threshold insensitivity rests on one fact: no plate centre lies inside the
    # gap between the northern and southern sets. Check that fact directly --
    # comparing the caller's selection against probes would also reject a
    # deliberately different --dec-min, which is a legitimate thing to ask for.
    inside = sorted(p for p, d in decs.items() if GAP[0] < d < GAP[1])
    if inside:
        sys.exit(
            f"[FAIL] {len(inside)} plate(s) have centre declinations inside the "
            f"{GAP[0]}..{GAP[1]} gap ({', '.join(inside[:5])}), so the northern and "
            f"southern sets are no longer cleanly separated and no single "
            f"threshold defines this footprint unambiguously. Re-derive the rule "
            f"before trusting the output."
        )
    if not (GAP[0] < args.dec_min < GAP[1]):
        print(f"[WARN] --dec-min {args.dec_min} lies outside the {GAP[0]}..{GAP[1]} "
              f"gap, so this is NOT the published footprint rule.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["plate_id"])
        for p in selected:
            w.writerow([p])

    below = len(decs) - len(selected)
    print(f"[OK] scanned {len(decs)} plates in {args.plate_dir}")
    print(f"     dec >= {args.dec_min}: {len(selected)} selected, {below} below the limit")
    print(f"     range {selected[0]}..{selected[-1]}  ->  {args.out}")


if __name__ == "__main__":
    main()
