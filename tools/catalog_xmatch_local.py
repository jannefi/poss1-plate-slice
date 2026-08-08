#!/usr/bin/env python3
"""Crossmatch a positional catalogue against a local catalogue mirror (Gaia DR3,
USNO-B1.0, ...), with a null control, and report matched fraction vs radius.

Why the null control is mandatory
---------------------------------
These reference catalogues are dense. Against Gaia DR3, roughly 6% of
*arbitrary* sky positions have a source within 30 arcsec purely by chance, and
at wider radii it approaches certainty. A bare "N% of rows have a counterpart"
number is therefore uninterpretable on its own. This tool runs the identical
match against the same positions with RA shifted by --null-shift degrees, which
preserves the local reference surface density while destroying any real
association, and reports both.

This exists because a headline project claim -- that most rows of a reference
catalogue had tight Gaia counterparts, supposedly contradicting that catalogue's
stated "no counterpart within 5 arcsec" construction -- turned out to conflate
two different measurements. Re-running this end to end refuted it: the
catalogue's own published coordinates had zero matches inside 5 arcsec, with a
step function at exactly the stated veto radius, which is the signature of a
correctly applied veto. Any similar claim should be checked with this tool
before it is published.

Memory design -- read this before changing the loop
---------------------------------------------------
The first version built the KD-tree on the *reference* side: pull a batch of
HEALPix pixels into memory, tree them, query the catalogue against it. On
USNO-B1.0 (1.05B rows, ~85k rows per healpix-5 pixel) a batch of 120 core pixels
plus neighbours reached 24.9 GB RSS on a 30 GB machine and drove it into swap.
It survived Gaia only by luck.

This version streams instead. The reference mirror is read in bounded row
batches; each batch is treed and both the real and the null catalogue positions
are queried against it in the same pass, keeping a running minimum separation.
Peak memory is O(one batch) + O(catalogue), independent of mirror size, and one
pass serves both real and null so the mirror is read once rather than twice.

Exactness note: taking a per-batch k=1 nearest neighbour and then the running
minimum across batches is exact for "distance to the nearest reference source".
Querying from the reference side instead would NOT be -- two catalogue rows can
share a nearest reference row, and only one of them would learn about it.

Fully local: nothing is sent over the network, so a private catalogue's
coordinates never leave the machine.

Run it under a cgroup cap
-------------------------
On janne-pc (30GB, dedicated to this project) always launch via:

    sudo systemd-run --scope -p MemoryMax=24G -p MemorySwapMax=0 \
        --uid=janne --gid=janne --working-directory="$PWD" --quiet \
        python3 tools/...

`MemorySwapMax=0` is the setting that matters. With 31GB of swap the kernel
prefers swapping to killing, so a runaway grows past RAM and thrashes the box
unresponsive rather than dying -- that is what froze this machine three times on
2026-08-05. The cap is not resource rationing on a dedicated box; it converts
"machine needs a hard reset" into "one process exits 137".

How to validate
---------------
    python3 tools/catalog_xmatch_local.py \
        --catalog-csv /path/to/catalogue.csv \
        --cache <gaia_cache>/parquet --cache-label gaia \
        --out-dir work/gaia_xmatch

Read the 1 arcsec binned table, not just the cumulative one. A veto applied at
radius R shows as exactly zero real matches below R and a sharp step at R; a
genuine physical association shows as a broad excess of real over null; and
"real approximately equals null" means no signal at all, however large the raw
percentage looks.
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import math
from pathlib import Path

import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds
from astropy import units as u
from astropy_healpix import HEALPix
from scipy.spatial import cKDTree

RADII = [1, 2, 3, 5, 10, 15, 20, 30]


def xyz(ra, dec):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def chord(arcsec):
    return 2.0 * math.sin(math.radians(arcsec / 3600.0) / 2.0)


def mem_available_gb() -> float:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024 / 1024
    except Exception:
        pass
    return float("inf")


def load_positions(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        head = csv.DictReader(f).fieldnames or []
    ra_col = next((c for c in head if c.strip().lower() in ("ra", "_ra", "raj2000", "ra_deg")), None)
    dec_col = next((c for c in head if c.strip().lower() in ("dec", "_dec", "dej2000", "dec_deg")), None)
    if not ra_col or not dec_col:
        raise SystemExit(f"[FATAL] no RA/Dec columns in {path}; saw {head[:12]}")
    ra, dec = [], []
    with path.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                ra.append(float(r[ra_col])); dec.append(float(r[dec_col]))
            except (TypeError, ValueError):
                continue
    return np.asarray(ra), np.asarray(dec), ra_col, dec_col


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--catalog-csv", required=True,
                    help="Positional catalogue. Path is never hardcoded.")
    ap.add_argument("--cache", required=True,
                    help="Parquet mirror dir, hive-partitioned on healpix_5, with ra/dec.")
    ap.add_argument("--cache-label", default="ref", help="Name used in output labels.")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--nside", type=int, default=32, help="HEALPix nside (level 5).")
    ap.add_argument("--hp-order", default="nested", choices=("nested", "ring"),
                    help="MUST match how the mirror's healpix_5 column was built. "
                         "The local mirrors are NESTED (see "
                         "vasco/local_cache_query.py::_get_hp). Using 'ring' here "
                         "silently reads the WRONG sky region -- the ids are valid "
                         "either way, so nothing errors, you just get a different "
                         "patch of sky and a near-empty match list.")
    ap.add_argument("--rows-per-batch", type=int, default=8_000_000,
                    help="Reference rows per streamed batch. ~8M keeps peak near "
                         "1.5GB; lower it on a small machine.")
    ap.add_argument("--null-shift", type=float, default=5.0,
                    help="RA shift in degrees for the null control.")
    ap.add_argument("--max-arcsec", type=float, default=max(RADII))
    ap.add_argument("--filter-col", default=None,
                    help="Optional mirror column to threshold, e.g. nDetections. "
                         "Use it to reproduce a production veto's own selectivity: "
                         "a raw PS1 match at 5 arcsec is near-saturated (78%% of "
                         "random positions hit something), so without the same cut "
                         "the pipeline applies, a real veto signature is invisible.")
    ap.add_argument("--filter-min", type=float, default=None,
                    help="Keep mirror rows with --filter-col >= this value.")
    ap.add_argument("--min-free-gb", type=float, default=4.0,
                    help="Refuse to start below this much available RAM.")
    ap.add_argument("--dataset-refresh-pixels", type=int, default=100,
                    help="Rebuild the pyarrow Dataset every N pixels. A "
                         "long-lived Dataset accumulates per-fragment state as "
                         "more of the sky is touched, so RSS grows with pixels "
                         "PROCESSED, not rows read -- measured on USNO-B: "
                         "monotonic climb to an 8GB OOM kill, versus a flat "
                         "~2.2GB plateau when refreshed every 100 pixels. "
                         "0 disables refreshing (do not, on a large mirror).")
    args = ap.parse_args()

    free = mem_available_gb()
    if free < args.min_free_gb:
        raise SystemExit(f"[FATAL] only {free:.1f}GB RAM available, need "
                         f"{args.min_free_gb}GB. Lower --rows-per-batch or free memory.")
    print(f"[MEM] {free:.1f}GB available, batch={args.rows_per_batch:,} rows")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ra, dec, rc, dcol = load_positions(Path(args.catalog_csv))
    n = len(ra)
    print(f"[CAT] {n} rows from {Path(args.catalog_csv).name} (cols {rc}/{dcol})")

    null_ra = (ra + args.null_shift) % 360.0
    real_tree = cKDTree(xyz(ra, dec))
    null_tree = cKDTree(xyz(null_ra, dec))
    sep_real = np.full(n, np.inf)
    sep_null = np.full(n, np.inf)
    lim = chord(args.max_arcsec)

    # One pass, one read per pixel.
    #
    # A first attempt streamed the whole filtered mirror in fixed row batches and
    # queried the ENTIRE catalogue against every batch. pyarrow ignores
    # batch_size once a filter is applied and emits one batch per row group, so
    # that became ~42,000 batches of ~800 rows, each triggering a full
    # 107k-point query -- hours of work for no reason.
    #
    # Instead: index the catalogue by pixel, then for each pixel Q read Q alone
    # and query only the catalogue points sitting in Q or in a neighbour of Q.
    # That is exactly the set of points a source in Q could match, because the
    # neighbour relation is symmetric and max_arcsec is far smaller than a
    # pixel. Every mirror row is read once, and every tree is one pixel wide.
    hp = HEALPix(nside=args.nside, order=args.hp_order)
    print(f"[HP] nside={args.nside} order={args.hp_order}")

    def index_by_pixel(a, d):
        px = np.asarray(hp.lonlat_to_healpix(a * u.deg, d * u.deg))
        out = {}
        for i, p in enumerate(px):
            out.setdefault(int(p), []).append(i)
        return {k: np.asarray(v) for k, v in out.items()}

    real_by_pix = index_by_pixel(ra, dec)
    null_by_pix = index_by_pixel(null_ra, dec)

    def nbrs(p):
        vals = np.atleast_1d(hp.neighbours(int(p))).ravel()
        return [int(x) for x in vals if np.isfinite(x) and x >= 0]

    to_read = set()
    for src in (real_by_pix, null_by_pix):
        for p in src:
            to_read.add(int(p))
            to_read.update(nbrs(p))
    to_read = sorted(to_read)
    print(f"[SCAN] {len(to_read)} of {hp.npix} healpix-5 pixels to read", flush=True)

    dset = ds.dataset(args.cache, format="parquet", partitioning="hive")
    seen = 0
    for i, Q in enumerate(to_read, 1):
        # Periodically drop and rebuild the Dataset. See --dataset-refresh-pixels:
        # pyarrow accumulates per-fragment state across the sky, so without this
        # RSS climbs with pixel count until the process is OOM-killed.
        if args.dataset_refresh_pixels and i > 1 and \
                i % args.dataset_refresh_pixels == 1:
            del dset
            gc.collect()
            dset = ds.dataset(args.cache, format="parquet", partitioning="hive")
        cand = [Q] + nbrs(Q)
        r_idx = np.concatenate([real_by_pix[c] for c in cand if c in real_by_pix]) \
            if any(c in real_by_pix for c in cand) else None
        n_idx = np.concatenate([null_by_pix[c] for c in cand if c in null_by_pix]) \
            if any(c in null_by_pix for c in cand) else None
        if r_idx is None and n_idx is None:
            continue
        # Read the pixel INCREMENTALLY, never as one table. A single healpix-5
        # pixel is small on average (~85k rows for USNO-B) but Galactic-plane
        # pixels are far larger, and materialising one whole plus its KD-tree is
        # what put this machine into swap. Row groups are tens of MB, and the
        # candidate set per pixel is tiny (tens of points), so treeing each
        # sub-batch and re-querying costs almost nothing.
        _flt = pc.field("healpix_5") == Q
        if args.filter_col and args.filter_min is not None:
            _flt = _flt & (pc.field(args.filter_col) >= args.filter_min)
        for sub in dset.scanner(columns=["ra", "dec"], filter=_flt).to_batches():
            if sub.num_rows == 0:
                continue
            bt = cKDTree(xyz(sub.column("ra").to_numpy(zero_copy_only=False),
                             sub.column("dec").to_numpy(zero_copy_only=False)))
            seen += sub.num_rows
            for idx, tree, sep in ((r_idx, real_tree, sep_real),
                                   (n_idx, null_tree, sep_null)):
                if idx is None or idx.size == 0:
                    continue
                d, _ = bt.query(tree.data[idx], k=1, distance_upper_bound=lim)
                ok = np.isfinite(d)
                if ok.any():
                    s = np.degrees(2.0 * np.arcsin(np.clip(d[ok] / 2.0, 0, 1))) * 3600.0
                    hit = idx[ok]
                    sep[hit] = np.minimum(sep[hit], s)
            del bt
        if i % 250 == 0:
            print(f"  [{i}/{len(to_read)} pixels] {seen:,} reference rows, "
                  f"{mem_available_gb():.1f}GB free", flush=True)
    print(f"[SCAN] done, {seen:,} reference rows read")

    print(f"\n{'radius':>8} | {'REAL matched':>17} | {'NULL (chance)':>17} | {'excess':>9}")
    print("-" * 64)
    cum = []
    for R_ in RADII:
        if R_ > args.max_arcsec:
            continue
        r_n, n_n = int((sep_real <= R_).sum()), int((sep_null <= R_).sum())
        cum.append(dict(radius_arcsec=R_, real=r_n, null=n_n, excess=r_n - n_n,
                        real_pct=100 * r_n / n, null_pct=100 * n_n / n))
        print(f"{R_:>6}\" | {r_n:9d} {100*r_n/n:6.2f}% | {n_n:9d} {100*n_n/n:6.2f}% | {r_n-n_n:+9d}")

    print("\n[1-arcsec bins]  real / null  (excess)")
    binned = []
    for lo in range(0, int(min(20, args.max_arcsec))):
        r_n = int(((sep_real > lo) & (sep_real <= lo + 1)).sum())
        n_n = int(((sep_null > lo) & (sep_null <= lo + 1)).sum())
        binned.append(dict(lo=lo, hi=lo + 1, real=r_n, null=n_n, excess=r_n - n_n))
        print(f"  {lo:2d}-{lo+1:2d}\": {r_n:6d} / {n_n:6d}   ({r_n-n_n:+6d})")

    np.save(out_dir / "sep_real.npy", sep_real)
    np.save(out_dir / "sep_null.npy", sep_null)
    (out_dir / "xmatch_summary.json").write_text(json.dumps(dict(
        n_rows=n, catalog=str(args.catalog_csv), cache=str(args.cache),
        cache_label=args.cache_label, null_shift_deg=args.null_shift,
        hp_nside=args.nside, hp_order=args.hp_order,
        filter_col=args.filter_col, filter_min=args.filter_min,
        reference_rows_streamed=seen, cumulative=cum, binned_1arcsec=binned,
        note="A step function at radius R with zero real matches below it is the "
             "signature of a veto applied at R, not of absent counterparts. "
             "Compare real against null, never real alone.",
    ), indent=2))
    print(f"\nwrote {out_dir}/xmatch_summary.json")


if __name__ == "__main__":
    main()
