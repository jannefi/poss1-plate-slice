#!/usr/bin/env python3
"""Which PLATE EPOCHS back a USNO-B match? Decides whether the match is circular.

Context
-------
A blunt "57% of our S0 has a USNO-B source within 5 arcsec" was read as making
USNO-B unusable, on the grounds that USNO-B1.0 is built from digitised POSS-I
plates -- including the POSS-I E (red) plates this pipeline runs SExtractor on --
so a match could be the same photons on the same glass.

That is true only for matches backed *solely* by POSS-I red. USNO-B records
which plates contributed to each source:

    B1mag  POSS-I  O  (blue)   ~1950s, SAME survey, DIFFERENT exposure
    R1mag  POSS-I  E  (red)    ~1950s, the plates we detect on -> CIRCULAR
    B2mag  POSS-II J  (blue)   ~1990s, independent later epoch
    R2mag  POSS-II F  (red)    ~1990s, independent later epoch
    Imag   POSS-II N  (near-IR)~1990s, independent later epoch

So the match population splits three ways:
  R1-only            purely circular, carries no independent information
  R1+B1, no POSS-II  POSS-I both colours -- the MAPS-like case. The blue plate is
                     a separate exposure, so this DOES add information (an
                     emulsion defect on one plate does not repeat on the other),
                     but it is still 1950s-epoch only.
  any POSS-II        genuinely independent modern epoch, ~40 years later. This is
                     the subset that could serve as a real "still there later"
                     veto.

Note the near-IR Imag: per the standing no-IR-veto rule, POSS-II N must not be
counted on its own as optical evidence, so it is reported separately.

Fully local. Run under a cgroup cap (see tools/catalog_xmatch_local.py).
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from vasco.paths import get as _p

import argparse
import csv
import gc
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds
from astropy import units as u
from astropy_healpix import HEALPix
from scipy.spatial import cKDTree

MAGS = ["B1mag", "R1mag", "B2mag", "R2mag", "Imag"]


def xyz(ra, dec):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def chord(arcsec):
    return 2.0 * math.sin(math.radians(arcsec / 3600.0) / 2.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--catalog-csv", required=True)
    ap.add_argument("--cache", default=str(_p("usnob_cache") / "parquet"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--nside", type=int, default=32)
    ap.add_argument("--dataset-refresh-pixels", type=int, default=100)
    args = ap.parse_args()

    ra, dec = [], []
    with Path(args.catalog_csv).open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                ra.append(float(r["ra"])); dec.append(float(r["dec"]))
            except (TypeError, ValueError, KeyError):
                continue
    ra, dec = np.asarray(ra), np.asarray(dec)
    n = len(ra)
    print(f"[CAT] {n} rows")

    hp = HEALPix(nside=args.nside, order="nested")   # mirrors are NESTED
    tree = cKDTree(xyz(ra, dec))
    best = np.full(n, np.inf)
    flags = np.zeros((n, len(MAGS)), dtype=bool)

    px = np.asarray(hp.lonlat_to_healpix(ra * u.deg, dec * u.deg))
    by_pix = {}
    for i, p in enumerate(px):
        by_pix.setdefault(int(p), []).append(i)
    by_pix = {k: np.asarray(v) for k, v in by_pix.items()}

    def nbrs(p):
        v = np.atleast_1d(hp.neighbours(int(p))).ravel()
        return [int(x) for x in v if np.isfinite(x) and x >= 0]

    to_read = sorted({q for p in by_pix for q in [int(p)] + nbrs(p)})
    print(f"[SCAN] {len(to_read)} pixels", flush=True)

    dset = ds.dataset(args.cache, format="parquet", partitioning="hive")
    lim = chord(args.radius_arcsec)
    for i, Q in enumerate(to_read, 1):
        if args.dataset_refresh_pixels and i > 1 and i % args.dataset_refresh_pixels == 1:
            del dset; gc.collect()
            dset = ds.dataset(args.cache, format="parquet", partitioning="hive")
        cand = [Q] + nbrs(Q)
        sel = [by_pix[c] for c in cand if c in by_pix]
        if not sel:
            continue
        idx = np.concatenate(sel)
        for sub in dset.scanner(columns=["ra", "dec"] + MAGS,
                                filter=pc.field("healpix_5") == Q).to_batches():
            if sub.num_rows == 0:
                continue
            bt = cKDTree(xyz(sub.column("ra").to_numpy(zero_copy_only=False),
                             sub.column("dec").to_numpy(zero_copy_only=False)))
            d, j = bt.query(tree.data[idx], k=1, distance_upper_bound=lim)
            ok = np.isfinite(d)
            if ok.any():
                s = np.degrees(2.0 * np.arcsin(np.clip(d[ok] / 2.0, 0, 1))) * 3600.0
                hit, jj = idx[ok], j[ok]
                better = s < best[hit]
                if better.any():
                    hh, jm = hit[better], jj[better]
                    best[hh] = s[better]
                    for k, mc in enumerate(MAGS):
                        col = sub.column(mc).to_numpy(zero_copy_only=False)
                        v = np.asarray(col, dtype=float)[jm]
                        flags[hh, k] = np.isfinite(v)
            del bt
        if i % 500 == 0:
            print(f"  [{i}/{len(to_read)}]", flush=True)

    matched = np.isfinite(best)
    m = flags[matched]
    r1 = m[:, MAGS.index("R1mag")]
    b1 = m[:, MAGS.index("B1mag")]
    p2 = m[:, MAGS.index("B2mag")] | m[:, MAGS.index("R2mag")]
    p2_opt = p2                                  # POSS-II optical only
    p2_ir = m[:, MAGS.index("Imag")]             # near-IR, reported separately

    nm = int(matched.sum())
    cls = Counter()
    for k in range(nm):
        if p2_opt[k]:
            cls["POSS-II optical (independent modern epoch)"] += 1
        elif b1[k]:
            cls["POSS-I O+E only (MAPS-like, 1950s both colours)"] += 1
        elif r1[k]:
            cls["POSS-I E only (CIRCULAR with our own detections)"] += 1
        else:
            cls["other/none of the above"] += 1

    print(f"\n=== matched within {args.radius_arcsec}\": {nm}/{n} ({100*nm/n:.2f}%) ===")
    for k, v in cls.most_common():
        print(f"  {v:8d}  {100*v/nm:6.2f}%   {k}")
    print(f"\n  (near-IR POSS-II N present on {int(p2_ir.sum())} = "
          f"{100*p2_ir.sum()/nm:.2f}% -- never counted as optical evidence)")
    print(f"  of ALL {n} catalogue rows: "
          f"{100*(np.array([p2_opt[k] for k in range(nm)]).sum())/n:.2f}% have a "
          f"POSS-II-optical-backed USNO-B match")

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    (Path(args.out_dir) / "usnob_epoch_breakdown.json").write_text(json.dumps({
        "n_catalog": n, "n_matched": nm, "radius_arcsec": args.radius_arcsec,
        "classes": dict(cls), "poss2_nir_present": int(p2_ir.sum()),
        "note": "R1-only is circular with POSS-I E detections. POSS-II optical is "
                "the only subset usable as an independent later-epoch veto.",
    }, indent=2))
    print(f"\nwrote {args.out_dir}/usnob_epoch_breakdown.json")


if __name__ == "__main__":
    main()
