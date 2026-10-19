#!/usr/bin/env python3
"""Match a published reference catalogue against this repository's released ones.

The cheap half of the comparison. Given a reference catalogue and one or more
of this repository's released candidate catalogues, this answers one question:
**what fraction of the reference rows have one of ours within N arcseconds?**

It needs nothing but two CSV files, runs in seconds on a laptop, and is the
tool to reach for when checking a published recall figure. For *which stage*
consumed the rows that do not match, use `tools/funnel_attribution.py` instead
-- that one needs the per-tile stage catalogues, which are tens of gigabytes
and are not shipped with this repository.

THE ONE THING THIS TOOL CANNOT TELL YOU
---------------------------------------
A non-match here has two completely different causes and this tool cannot
separate them:

  * we never detected anything at that position, or
  * we detected it and a later stage removed it.

Those mean opposite things about the pipeline, and the headline figures differ
by a factor of five depending on which you mean -- raw detection runs in the
high nineties, candidate-level survival around a fifth. **Quote a number from
this tool as "has a row in <named catalogue> within N arcsec", never as
"detected" or "recovered".** Separating the two causes is what
`funnel_attribution.py` is for.

Positional match is also not identity: a reference row and one of ours at the
same position may have been selected by different pipelines for different
reasons. Agreement here is agreement of position, nothing more.

INPUT
-----
`--ref-csv` is a local path. This tool never downloads anything. Reference
catalogues are other people's data, published under their own terms; fetching
and caching them inside a repository would redistribute them, so this tool
reads a file you obtained yourself and this repository ships none.

It reads SVO-style exports directly, including the leading `#` on the header
line that `cs.php?...&format=ascii` emits, and finds the coordinate columns
under any of the usual spellings. If the reference carries a plate column it is
used to break the results down per plate.

How to validate
---------------
Against a published number, before trusting it on anything else:

    python3 tools/reference_match.py \\
        --ref-csv <your copy of the 5,399-row vanish-possi export> \\
        --catalog results/s0-642-20260814/stage_S0.csv.gz \\
        --match-radius-arcsec 5

`docs/FUNNEL_ATTRIBUTION.md` reports 1,072 of 5,399 (19.86%) for exactly this
comparison. If this tool does not reproduce that, do not use its output.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

RA_NAMES = ("ra", "ra_deg", "raj2000", "_raj2000", "alpha", "alpha_j2000")
DEC_NAMES = ("dec", "de", "dec_deg", "dej2000", "_dej2000", "delta", "delta_j2000")
PLATE_NAMES = ("name", "plate", "plate_id", "det_plate", "field")


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "r")


def read_table(path: Path) -> pd.DataFrame:
    """Read a CSV, tolerating the leading '#' SVO puts on its header line."""
    with _open(path) as f:
        first = f.readline()
    header = 0
    names = None
    if first.lstrip().startswith("#"):
        names = [c.strip() for c in first.lstrip().lstrip("#").strip().split(",")]
        header = 0
        with _open(path) as f:
            return pd.read_csv(f, skiprows=1, names=names)
    with _open(path) as f:
        return pd.read_csv(f, header=header)


def find_col(df: pd.DataFrame, candidates) -> str | None:
    lower = {c.lower().strip(): c for c in df.columns}
    for c in candidates:
        if c in lower:
            return lower[c]
    return None


def unit_vectors(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.radians(ra_deg.astype(float))
    dec = np.radians(dec_deg.astype(float))
    cd = np.cos(dec)
    return np.column_stack((cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)))


def chord(radius_arcsec: float) -> float:
    """Chord length on the unit sphere subtending the given angle."""
    return 2.0 * np.sin(np.radians(radius_arcsec / 3600.0) / 2.0)


def match(ref_xyz: np.ndarray, cat_xyz: np.ndarray, radius_arcsec: float):
    """Nearest catalogue row per reference row, and whether it is within radius."""
    tree = cKDTree(cat_xyz)
    dist, idx = tree.query(ref_xyz, k=1)
    sep_arcsec = np.degrees(2.0 * np.arcsin(np.clip(dist / 2.0, 0, 1))) * 3600.0
    return sep_arcsec, idx, sep_arcsec <= radius_arcsec


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Positional match of a reference catalogue against released catalogues.")
    ap.add_argument("--ref-csv", required=True, type=Path,
                    help="LOCAL path to the reference catalogue. Never downloaded.")
    ap.add_argument("--catalog", required=True, type=Path, action="append",
                    help="a released catalogue (repeatable); e.g. "
                         "results/s0-642-20260814/stage_S0.csv.gz")
    ap.add_argument("--match-radius-arcsec", type=float, default=5.0)
    ap.add_argument("--sweep", action="store_true",
                    help="also report 1/2/3/5/10 arcsec, to show radius sensitivity")
    ap.add_argument("--per-plate", action="store_true",
                    help="break results down by the reference's plate column")
    ap.add_argument("--out", type=Path, default=None,
                    help="per-reference-row CSV. Contains reference rows -- write it "
                         "OUTSIDE this repository and do not commit it.")
    ap.add_argument("--out-summary", type=Path, default=None, help="summary JSON")
    args = ap.parse_args()

    ref = read_table(args.ref_csv)
    rra, rdec = find_col(ref, RA_NAMES), find_col(ref, DEC_NAMES)
    if not rra or not rdec:
        print(f"[FATAL] no RA/Dec columns in {args.ref_csv} (saw {list(ref.columns)})")
        return 2
    ref = ref[np.isfinite(pd.to_numeric(ref[rra], errors="coerce"))
              & np.isfinite(pd.to_numeric(ref[rdec], errors="coerce"))].copy()
    ref_xyz = unit_vectors(pd.to_numeric(ref[rra]).values, pd.to_numeric(ref[rdec]).values)
    plate_col = find_col(ref, PLATE_NAMES) if args.per_plate else None

    print(f"[REF] {args.ref_csv}  rows={len(ref)}  coords={rra}/{rdec}"
          + (f"  plate={plate_col}" if plate_col else ""))

    summary = {"reference": str(args.ref_csv), "reference_rows": int(len(ref)),
               "match_radius_arcsec": args.match_radius_arcsec, "catalogs": {}}
    out = ref[[c for c in (rra, rdec, plate_col) if c]].copy()

    for cat_path in args.catalog:
        cat = read_table(cat_path)
        cra, cdec = find_col(cat, RA_NAMES), find_col(cat, DEC_NAMES)
        if not cra or not cdec:
            print(f"[FATAL] no RA/Dec columns in {cat_path} (saw {list(cat.columns)})")
            return 2
        cat_xyz = unit_vectors(pd.to_numeric(cat[cra]).values,
                               pd.to_numeric(cat[cdec]).values)
        sep, _, hit = match(ref_xyz, cat_xyz, args.match_radius_arcsec)
        n, tot = int(hit.sum()), len(ref)
        label = cat_path.parent.name or cat_path.stem
        print(f"\n[CATALOG] {cat_path}  rows={len(cat)}")
        print(f"  reference rows with a row in this catalogue within "
              f"{args.match_radius_arcsec:g}\": {n} of {tot} ({100*n/tot:.2f}%)")
        entry = {"path": str(cat_path), "rows": int(len(cat)),
                 "matched": n, "fraction": round(n / tot, 6)}

        if args.sweep:
            entry["sweep"] = {}
            print("    radius sweep:")
            for r in (1.0, 2.0, 3.0, 5.0, 10.0):
                k = int((sep <= r).sum())
                entry["sweep"][f"{r:g}"] = {"matched": k, "fraction": round(k / tot, 6)}
                print(f"      {r:>5.1f}\"  {k:>7}  {100*k/tot:6.2f}%")

        if plate_col:
            g = pd.DataFrame({"plate": ref[plate_col].values, "hit": hit})
            per = g.groupby("plate")["hit"].agg(["size", "sum"])
            per["frac"] = per["sum"] / per["size"]
            entry["per_plate"] = {"plates": int(len(per)),
                                  "median_fraction": round(float(per["frac"].median()), 6),
                                  "min_fraction": round(float(per["frac"].min()), 6),
                                  "max_fraction": round(float(per["frac"].max()), 6)}
            print(f"    per plate ({len(per)} plates): median {100*per['frac'].median():.2f}%"
                  f"  range {100*per['frac'].min():.2f}-{100*per['frac'].max():.2f}%")

        out[f"sep_arcsec__{label}"] = np.round(sep, 4)
        out[f"matched__{label}"] = hit
        summary["catalogs"][label] = entry

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, index=False)
        print(f"\n[WROTE] {args.out}  ({len(out)} rows) -- contains reference rows; "
              f"do not commit it")
    if args.out_summary:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        args.out_summary.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"[WROTE] {args.out_summary}")

    print("\nA match means a row of ours lies within the radius. It does NOT mean "
          "'detected':\n  a non-match may be a row we never detected or one a later "
          "stage removed.\n  Use tools/funnel_attribution.py to separate those.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
