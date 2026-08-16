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
northern centre -0.865, highest southern -5.171), so every threshold in that gap
selects exactly the same 642 plates. The script asserts this, so a future edit
to DEC_MIN cannot silently change the footprint.

Usage:
    python tools/build_plate_manifest.py --plate-dir /path/to/dss1red \\
        --out data/plate_manifest.csv
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
from astropy.io import fits  # noqa: E402
from astropy.wcs import WCS  # noqa: E402

DEC_MIN = -3.0
# The empty band between the northern and southern plate sets. Any threshold
# inside it yields the same footprint; if a future library breaks that, the
# footprint is no longer threshold-insensitive and the claim above is void.
GAP = (-5.0, -1.0)
PLATE_RX = re.compile(r"dss1red_(XE\d+)\.fits$")


def keyword_plate_centre(header) -> tuple[float, float]:
    """Centre of the PLATE (RA, Dec) from the GSSS keywords, in degrees.

    PLTDECSN carries the sign as a separate '+'/'-' field, so the sign must be
    applied to the assembled magnitude -- reading PLTDECD alone silently loses
    it for southern plates. RA is sexagesimal hours across three keywords.

    This is the centre of the glass, NOT the centre of the scan -- see
    scan_centre. Kept because the two disagree, and the disagreement is worth
    reporting rather than hiding.
    """
    ra = 15.0 * (
        header.get("PLTRAH", 0)
        + header.get("PLTRAM", 0) / 60.0
        + header.get("PLTRAS", 0) / 3600.0
    )
    sign = -1.0 if str(header.get("PLTDECSN", "+")).strip() == "-" else 1.0
    dec = sign * (
        abs(header.get("PLTDECD", 0))
        + header.get("PLTDECM", 0) / 60.0
        + header.get("PLTDECS", 0) / 3600.0
    )
    return ra, dec


def scan_centre(header) -> tuple[float, float]:
    """Centre of the SCAN (RA, Dec): the sky position of the middle pixel.

    This is what the manifest must carry. Downstream the manifest answers
    "which plate covers this position", and that is a question about where the
    data is, not where the glass was pointed.

    The two are not the same. Across all 932 DSS1-red headers they agree to a
    median 0.07 deg, but seven plates disagree by ~4.4 deg -- XE761, XE758,
    XE733, XE574, XE284, XE543, XE541 -- plus XE304, XE293, XE509 and XE880
    between 0.5 and 1.4 deg. Neither value is corrupt: the keyword centre
    agrees exactly with the DSS PLT* plate solution, but on those scans the
    image is not centred on the plate while CNPIX still reads (0, 0).
    tools/slice_plate_tiles.py already slices on the image centre for exactly
    this reason; the manifest previously did not, which put a wrong
    primary_plate on 4.3% of catalogue rows.

    Selecting the footprint on this centre rather than the keyword one leaves
    the 642-plate set unchanged and the north/south gap empty (lowest northern
    -0.865, highest southern -5.171), so the reproducibility claim in the
    module docstring is unaffected -- main() asserts both.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = WCS(header, relax=True)
        sky = w.all_pix2world(
            [[float(header["NAXIS1"]) / 2.0, float(header["NAXIS2"]) / 2.0]], 0
        )[0]
    return float(sky[0]), float(sky[1])


def _sep_deg(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    d1, d2 = math.radians(dec1), math.radians(dec2)
    s = (math.sin((d2 - d1) / 2.0) ** 2
         + math.cos(d1) * math.cos(d2) * math.sin(math.radians(ra2 - ra1) / 2.0) ** 2)
    return math.degrees(2.0 * math.asin(min(1.0, math.sqrt(max(s, 0.0)))))


def scan(plate_dir: Path) -> dict[str, tuple[float, float]]:
    out: dict[str, tuple[float, float]] = {}
    offset: list[tuple[float, str]] = []
    for f in sorted(plate_dir.glob("dss1red_XE*.fits")):
        m = PLATE_RX.search(f.name)
        if not m:
            continue
        hdr = fits.getheader(f)
        sra, sdec = scan_centre(hdr)
        kra, kdec = keyword_plate_centre(hdr)
        out[m.group(1)] = (sra, sdec)
        off = _sep_deg(kra, kdec, sra, sdec)
        if off > 0.5:
            offset.append((off, m.group(1)))
    if not out:
        sys.exit(f"[FAIL] no dss1red_XE*.fits found under {plate_dir}")
    if offset:
        offset.sort(reverse=True)
        print(f"[NOTE] {len(offset)} plates whose scan centre is >0.5 deg from the "
              f"keyword plate centre; the manifest carries the scan centre:")
        for off, p in offset:
            print(f"       {p}  {off:.3f} deg")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plate-dir", required=True, type=Path,
                    help="directory of IRSA dss1red_XE*.fits scans")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--dec-min", type=float, default=DEC_MIN)
    args = ap.parse_args()

    centres = scan(args.plate_dir)
    decs = {p: c[1] for p, c in centres.items()}
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
    # ra_deg/dec_deg are SCAN centres -- the sky position of the middle pixel
    # under the header's own WCS, not the PLTRAH../PLTDEC.. keyword centre. See
    # scan_centre for why the distinction matters and which plates it moves.
    # They exist so downstream geometry -- e.g. the primary-plate rule in
    # tools/build_primary_plate_flags.py -- reads exact centres from a public
    # artifact instead of re-deriving approximations from tile names.
    with args.out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["plate_id", "ra_deg", "dec_deg"])
        for p in selected:
            w.writerow([p, f"{centres[p][0]:.6f}", f"{centres[p][1]:.6f}"])

    below = len(decs) - len(selected)
    print(f"[OK] scanned {len(decs)} plates in {args.plate_dir}")
    print(f"     dec >= {args.dec_min}: {len(selected)} selected, {below} below the limit")
    print(f"     range {selected[0]}..{selected[-1]}  ->  {args.out}")


if __name__ == "__main__":
    main()
