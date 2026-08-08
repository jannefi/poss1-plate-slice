#!/usr/bin/env python3
"""Whole-sky union raw parity against a public reference catalogue.

Measures what fraction of a reference catalogue's rows have any detection within
a given radius, for each arm and for unions of arms. The reference used here is
the published POSS-I vanishing-source catalogue (SVO `vanish-possi`, 5,399 rows),
so the entire measurement is reproducible from public data.

Arms are directories of lean RA/Dec CSVs at any granularity (per plate, per
tile). Every arm is scored on the SAME reference rows.

Memory design -- read before changing the loop
----------------------------------------------
The full-plate arm alone is ~186 million detections. Treeing that at once costs
well over 10 GB and there are several arms. Instead each arm is streamed in
bounded chunks: tree the chunk, query all reference rows against it with a
distance cap, and keep a running minimum. Taking a per-chunk k=1 nearest
neighbour and then the minimum across chunks is exact for "distance to the
nearest detection", and peak memory is O(one chunk) + O(reference), independent
of arm size.

Querying from the detection side instead would NOT be exact -- two reference rows
can share a nearest detection and only one would learn about it.

Interpreting the result
-----------------------
`fullplate` is the locally-sliced arm (IRSA full-plate scans + the published
plate list); `archive` is the cutout-service arm. Their union is the method this
repository documents. Watch the 1-3" columns as well as 5": the per-plate CRPIX
correction (docs/DSS_WCS_TWO_SOLUTIONS.md) is what makes the tight radii
meaningful, and a regression there shows up at 1-2" long before it shows at 5".

How to validate
---------------
    python3 tools/union_parity_fullscale.py --ref-csv <reference.csv> \
        --arm fullplate=<dir> --arm archive=<dir> \
        --combine fullplate+archive \
        --out-dir work/union

Each arm must report a plausible detection count in its [ARM] line; an arm that
silently reads zero files would otherwise show as 0% and look like a result.
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

RADII = (1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 30.0, 60.0)


def xyz(ra, dec):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def chord(arcsec):
    return 2.0 * math.sin(math.radians(arcsec / 3600.0) / 2.0)


def load_ref(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        head = csv.DictReader(f).fieldnames or []
    rk = next(c for c in head if c.strip().lower() in ("ra", "_ra", "raj2000"))
    dk = next(c for c in head if c.strip().lower() in ("dec", "_dec", "dej2000"))
    d = pd.read_csv(path, usecols=[rk, dk])
    return d[rk].to_numpy(float), d[dk].to_numpy(float)


def stream_arm(directory: Path, v_xyz, max_arcsec, chunk_rows, label):
    """Running minimum separation from every reference row to this arm."""
    best = np.full(len(v_xyz), np.inf)
    lim = chord(max_arcsec)
    files = sorted(directory.glob("*.csv"))
    if not files:
        raise SystemExit(f"[FATAL] no CSVs under {directory}")
    buf_ra, buf_dec, n_buf, n_tot, n_chunks = [], [], 0, 0, 0

    def flush():
        nonlocal buf_ra, buf_dec, n_buf, n_chunks
        if not n_buf:
            return
        ra = np.concatenate(buf_ra); dec = np.concatenate(buf_dec)
        tree = cKDTree(xyz(ra, dec))
        d, _ = tree.query(v_xyz, k=1, distance_upper_bound=lim)
        np.minimum(best, d, out=best)
        n_chunks += 1
        del tree, ra, dec
        buf_ra, buf_dec, n_buf = [], [], 0
        gc.collect()

    for i, f in enumerate(files, 1):
        try:
            d = pd.read_csv(f, usecols=lambda c: c.strip().lower() in ("ra", "dec"))
        except Exception:
            continue
        cols = {c.strip().lower(): c for c in d.columns}
        if "ra" not in cols or "dec" not in cols or not len(d):
            continue
        buf_ra.append(d[cols["ra"]].to_numpy(float))
        buf_dec.append(d[cols["dec"]].to_numpy(float))
        n_buf += len(d); n_tot += len(d)
        if n_buf >= chunk_rows:
            flush()
        if i % 5000 == 0:
            print(f"    [{label}] {i}/{len(files)} files, {n_tot:,} detections", flush=True)
    flush()
    print(f"  [ARM] {label:<12} {len(files):>6} files, {n_tot:>12,} detections, "
          f"{n_chunks} chunks", flush=True)
    # chord -> arcsec
    out = np.full(len(v_xyz), np.inf)
    ok = np.isfinite(best)
    out[ok] = np.degrees(2.0 * np.arcsin(np.clip(best[ok] / 2.0, 0, 1))) * 3600.0
    return out, n_tot


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ref-csv", required=True,
                    help="Public reference catalogue with RA/Dec columns.")
    ap.add_argument("--arm", action="append", required=True, help="LABEL=DIR_OF_CSVS")
    ap.add_argument("--combine", action="append", default=[], help="LABEL_A+LABEL_B")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-arcsec", type=float, default=60.0)
    ap.add_argument("--chunk-rows", type=int, default=4_000_000)
    args = ap.parse_args()

    ref_ra, ref_dec = load_ref(Path(args.ref_csv))
    refx = xyz(ref_ra, ref_dec)
    print(f"[IN] reference rows {len(ref_ra):,}")

    dists, counts = {}, {}
    for spec in args.arm:
        label, d = spec.split("=", 1)
        print(f"[ARM] streaming {label} from {d}", flush=True)
        dists[label], counts[label] = stream_arm(Path(d), refx, args.max_arcsec,
                                                 args.chunk_rows, label)
    for spec in args.combine:
        labels = spec.split("+")
        miss = [l for l in labels if l not in dists]
        if miss:
            print(f"[WARN] skipping {spec}: unknown arm(s) {miss}")
            continue
        dists[spec] = np.minimum.reduce([dists[l] for l in labels])
        counts[spec] = sum(counts[l] for l in labels)

    print(f"\n=== WHOLE-SKY UNION RAW PARITY on {len(ref_ra):,} reference rows ===")
    hdr = "arm".ljust(22) + "".join(f"{r:>8.0f}\"" for r in RADII)
    print(hdr); print("-" * len(hdr))
    out = {"n_rows": int(len(ref_ra)), "by_arm": {}}
    for label, dv in dists.items():
        row = [100.0 * np.mean(dv <= r) for r in RADII]
        out["by_arm"][label] = {"detections": int(counts.get(label, 0)),
                                **{str(r): round(x, 3) for r, x in zip(RADII, row)}}
        print(label.ljust(22) + "".join(f"{x:>8.2f}%" for x in row))

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    (Path(args.out_dir) / "union_fullscale.json").write_text(json.dumps(out, indent=2))
    np.save(Path(args.out_dir) / "dists.npy",
            np.vstack([dists[k] for k in dists]))
    (Path(args.out_dir) / "arm_order.json").write_text(json.dumps(list(dists)))
    print(f"\nwrote {args.out_dir}/union_fullscale.json")


if __name__ == "__main__":
    main()
