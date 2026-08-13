#!/usr/bin/env python3
"""
Post-pipeline USNO-B1.0 reduction stage, backed by the LOCAL mirror.

Relationship to scripts/stage_usno_post.py
------------------------------------------
That script predates the local catalogue mirrors: it drives STILTS
`cdsskymatch` against VizieR I/284/out over the network, and consumes the old
chunked `upload_positional_*.csv` convention. Both are obsolete for this
purpose -- a live CDS crossmatch of a 200k+ row stage is slow, rate-limited and
non-reproducible. This is a local-mirror replacement with the same
kept/flags/ledger output contract as `stage_gsc_post.py` / `stage_maps_post.py`,
reading `stage_S0.csv` (or any `src_id,ra,dec`-bearing stage CSV) directly.

The old script is left in place untouched; nothing references either.

Why USNO-B is worth testing at all
----------------------------------
The pipeline's live vetoes are Gaia + PS1 only. USNO-B is wired into
`vasco/cli_pipeline.py` but switched OFF for the whole archive via
`VASCO_DISABLE_USNOB=1`, a deliberate documented deviation
(context/REPRO_DEVIATIONS.md item 4) -- verified in the archive itself:
`usnob_veto_enabled: false` and zero rows removed in 300/300 sampled tiles.
USNO-B is also the catalogue MNRAS 2022 used for its spike check, and it is
purely optical, so unlike GSC 2.4.2 it does not fall foul of the standing
"no IR-derived veto stages" rule.

The open question this answers: does USNO-B remove anything the existing
Gaia+PS1 chain has not already removed, or is it redundant the way MAPS turned
out to be?

EXPERIMENTAL. Not wired into any orchestrator, same posture as GSC/MAPS/EDGE.

MEMORY: why this script froze a 30 GB machine twice, and what changed
--------------------------------------------------------------------
Both freezes came from a *scattered* input (rows sampled across the whole
642-plate footprint), and the row count was never the problem -- a 6,490-row
control killed the box while a 6,490-row single plate would have been fine.
The input rows only choose WHICH pixels are read, never how much is read.

Three compounding defects, all fixed here:

1. **The neighbour halo multiplied scattered input by ~9.** The old loop took
   every core pixel and unconditionally added all 8 `hp.neighbours(p)`. On
   compact sky those halos overlap and `need` is barely larger than the batch;
   on scattered sky nothing overlaps and `need` reaches 9x the batch. At the
   old default `--batch 120` that is ~1,080 of 12,288 pixels -- **8.8% of the
   entire mirror in one call.**

   The halo was also unnecessary. The match radius is 5" against a ~1.8 deg
   nside=32 pixel, so a neighbour can only matter when a query point lies
   within 5" of a pixel boundary -- a ~0.1% minority. Pixel selection is now
   per query point and exact (see `_pixels_for_points`).

2. **`to_table` materialised the whole batch at once**, and `col()` did
   `.to_numpy(zero_copy_only=False).astype(float)` -- two further full copies
   per column, for up to six columns. Reading is now streamed through
   `scanner().to_batches()`, so peak memory is set by one record batch, not by
   how much sky the batch happens to cover.

3. **The batch was bounded by pixel COUNT while pixels vary ~20x in size**
   (measured on the mirror: median 46,845 rows, p90 182,129, p99 605,981, max
   923,532). Same failure mode as an `lru_cache(maxsize=N)` over
   wildly-varying values. Batches are now bounded by estimated ROWS, from a
   one-off parquet-footer scan cached beside the mirror.

Old worst case, scattered + p99 sky: >31 GB for the Arrow table alone, before
copies. That is the freeze.

`--legacy-halo` restores the old pixel selection for A/B comparison. It is
there to prove the two agree on veto decisions, not for production use.

Outputs
-------
1) stage_<STAGE>_USNOB.csv        kept remainder (rows with NO USNO-B match)
2) stage_<STAGE>_USNOB_flags.csv  audit row per input: match flag + separation
3) stage_<STAGE>_USNOB_ledger.json counts + parameters

How to validate
---------------
    VASCO_USNOB_CACHE=<usnob_cache> \\
    python3 scripts/stage_usnob_post.py \\
        --stage-csv work/runs/<run>/stage_S0.csv \\
        --out-dir work/runs/<run>/stages --stage S0

Check the ledger's `matched_rows` against a null run (`--null-shift 5.0`, which
offsets every position in RA to destroy real association while preserving local
USNO-B density). USNO-B is dense; a raw match rate is not interpretable without
that comparison.

Prefer driving this **per plate**. Every batch is then compact by construction,
which bounds the working set independently of anything above.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from astropy import units as u
from astropy_healpix import HEALPix
from scipy.spatial import cKDTree

_AZIMUTHS = np.radians(np.arange(0.0, 360.0, 45.0))   # 8 samples on the disc rim


def xyz(ra, dec):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def chord(arcsec: float) -> float:
    return 2.0 * math.sin(math.radians(arcsec / 3600.0) / 2.0)


def _offset(ra_deg, dec_deg, dist_rad, az_rad):
    """Positions `dist_rad` from (ra,dec) along azimuth `az_rad`, on the sphere."""
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    sd, cd = np.sin(dec), np.cos(dec)
    sD, cD = math.sin(dist_rad), math.cos(dist_rad)
    dec2 = np.arcsin(np.clip(sd * cD + cd * sD * np.cos(az_rad), -1.0, 1.0))
    ra2 = ra + np.arctan2(np.sin(az_rad) * sD * cd, cD - sd * np.sin(dec2))
    return np.degrees(ra2) % 360.0, np.degrees(dec2)


def _pixels_for_points(hp, ra, dec, radius_arcsec, own_pix):
    """Pixels that can hold a source within `radius_arcsec` of each point.

    Exact, and cheap because the radius is ~1/1300 of a pixel. Each point's own
    pixel is always required. A neighbour is required only when the point sits
    within the radius of a boundary, which is detected by sampling 8 positions
    on the rim of the disc; any point whose rim leaves its own pixel is then
    resolved exactly with the pipeline's own overlap test, so a sliver thinner
    than the rim sampling cannot be missed.

    Returns (all_pixels, extra_by_index) -- the union, and per-point extras for
    points that straddle a boundary.
    """
    d = math.radians(radius_arcsec / 3600.0)
    allpix = set(int(p) for p in own_pix)
    extra = defaultdict(set)
    straddle = np.zeros(len(ra), dtype=bool)

    for az in _AZIMUTHS:
        r2, d2 = _offset(ra, dec, d, az)
        p2 = np.asarray(hp.lonlat_to_healpix(r2 * u.deg, d2 * u.deg), dtype=np.int64)
        diff = p2 != own_pix
        if diff.any():
            straddle |= diff
            for i in np.flatnonzero(diff):
                extra[int(i)].add(int(p2[i]))
                allpix.add(int(p2[i]))

    # Exact resolution for the rare boundary cases. `_cone_pixels` is the
    # pipeline's own provably-complete overlap test -- the one whose absence
    # caused the 2026-08 partial-cone bug -- so this cannot under-select.
    n_str = int(straddle.sum())
    if n_str:
        try:
            from vasco.local_cache_query import _cone_pixels
            for i in np.flatnonzero(straddle):
                for p in _cone_pixels(hp, float(ra[i]), float(dec[i]),
                                      radius_arcsec / 60.0):
                    extra[int(i)].add(int(p))
                    allpix.add(int(p))
        except ImportError:
            # Fall back to the blanket halo for the straddlers only. Still far
            # smaller than haloing everything, and never under-selects.
            for i in np.flatnonzero(straddle):
                for p in np.atleast_1d(hp.neighbours(int(own_pix[i]))):
                    if p >= 0:
                        extra[int(i)].add(int(p))
                        allpix.add(int(p))
    return allpix, extra, n_str


def _pixel_row_counts(cache: str, rebuild: bool = False) -> dict[int, int]:
    """Rows per healpix_5 partition, from parquet footers. Cached beside the mirror.

    Needed to bound a batch by rows rather than by pixel count. The scan reads
    only footers (no data) but walks ~12k directories, so it is cached; the
    mirror is immutable in practice.
    """
    import glob
    import re
    cache_file = Path(cache).parent / "pixel_row_counts.json"
    if cache_file.exists() and not rebuild:
        try:
            return {int(k): int(v) for k, v in json.loads(cache_file.read_text()).items()}
        except Exception:
            pass
    counts: dict[int, int] = {}
    for f in glob.glob(f"{cache}/healpix_5=*/*.parquet"):
        m = re.search(r"healpix_5=(\d+)", f)
        if not m:
            continue
        p = int(m.group(1))
        counts[p] = counts.get(p, 0) + pq.ParquetFile(f).metadata.num_rows
    try:
        cache_file.write_text(json.dumps({str(k): v for k, v in counts.items()}))
    except OSError:
        pass
    return counts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--stage-csv", required=True,
                    help="Input stage CSV with src_id,ra,dec (e.g. stage_S0.csv).")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--stage", default="S0", help="Label used in output filenames.")
    ap.add_argument("--cache", default=None,
                    help="USNO-B parquet dir. Defaults to $VASCO_USNOB_CACHE/parquet.")
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--null-shift", type=float, default=0.0,
                    help="If non-zero, shift every RA by this many degrees before "
                         "matching -- a chance-rate control, not a real run.")
    ap.add_argument("--allow-poss1-only", action="store_true",
                    help="Match against ALL USNO-B entries, including those whose only "
                         "evidence is POSS-I. CIRCULAR with our own detections -- for "
                         "comparison runs only, never for a published veto.")
    ap.add_argument("--count-nir", action="store_true",
                    help="Also accept POSS-II N (Imag) as second-epoch evidence. Off by "
                         "default: it is a near-IR band, and no veto stage here may rest "
                         "on infrared evidence. Measured cost of excluding it is ~1%% of "
                         "the cut.")
    ap.add_argument("--plate-epoch", type=float, default=1953.0,
                    help="Epoch to propagate USNO-B positions back to, in years. POSS-I "
                         "spans roughly 1949-1958, so the default is its midpoint; the "
                         "residual from per-plate variation is ~0.5\" for a fast star, "
                         "against the ~4.7\" the propagation itself corrects.")
    ap.add_argument("--no-pm", action="store_true",
                    help="Skip proper-motion propagation. USNO-B positions are epoch "
                         "2000.0 and our detections are ~1950, so this leaves a "
                         "magnitude-correlated error; for comparison runs only.")
    ap.add_argument("--nside", type=int, default=32)
    ap.add_argument("--hp-order", default="nested", choices=("nested", "ring"),
                    help="Ordering of the mirror's healpix_5 column. The mirrors are "
                         "NESTED; 'ring' reads the wrong sky and fails silently.")
    ap.add_argument("--batch", type=int, default=120,
                    help="Maximum CORE pixels per batch. Also bounded by --max-batch-rows.")
    ap.add_argument("--max-batch-rows", type=int, default=40_000_000,
                    help="Cap on reference rows pulled per batch, from cached per-pixel "
                         "counts. ~40M rows x 6 float64 cols is ~1.9 GB. This, not "
                         "--batch, is what actually bounds memory: pixels vary ~20x in "
                         "size, so a count-based limit bounds nothing.")
    ap.add_argument("--scan-batch-rows", type=int, default=2_000_000,
                    help="Record-batch size for streaming. Peak memory scales with this.")
    ap.add_argument("--rebuild-pixel-counts", action="store_true",
                    help="Re-scan parquet footers instead of using the cached counts.")
    ap.add_argument("--legacy-halo", action="store_true",
                    help="Old pixel selection: blanket 8-neighbour halo on every core "
                         "pixel. For A/B verification only -- this is what made a "
                         "scattered input read ~9x more sky than it needed.")
    args = ap.parse_args()

    cache = args.cache
    if not cache:
        base = os.getenv("VASCO_USNOB_CACHE")
        if not base:
            print("[USNOB][FATAL] set VASCO_USNOB_CACHE or pass --cache", file=sys.stderr)
            return 2
        cache = str(Path(base) / "parquet")
    if not Path(cache).exists():
        print(f"[USNOB][FATAL] cache not found: {cache}", file=sys.stderr)
        return 2

    rows = []
    with Path(args.stage_csv).open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                rows.append((r.get("src_id", ""), float(r["ra"]), float(r["dec"])))
            except (TypeError, ValueError, KeyError):
                continue
    if not rows:
        print(f"[USNOB][FATAL] no usable rows in {args.stage_csv}", file=sys.stderr)
        return 2
    ra = np.array([r[1] for r in rows])
    dec = np.array([r[2] for r in rows])
    n = len(rows)
    print(f"[USNOB] input rows={n} radius={args.radius_arcsec}\" cache={cache}"
          + (f" NULL-SHIFT={args.null_shift} deg" if args.null_shift else ""), flush=True)

    q_ra = (ra + args.null_shift) % 360.0 if args.null_shift else ra
    # The mirror's `healpix_5` column is nside=32 NESTED. Reading it as RING
    # silently returns a DIFFERENT REGION OF SKY -- no error, no warning, just a
    # match rate that looks like a weak real signal. This was hardcoded to "ring"
    # and every run of this script before 2026-08-08 matched against the wrong
    # sky: 9.7% removal against a 6.2% null, where the correct value is 75%
    # against a 16% null. Print the order so the mistake is visible in any log.
    hp = HEALPix(nside=args.nside, order=args.hp_order)
    print(f"[HP] nside={args.nside} order={args.hp_order}", flush=True)
    dset = ds.dataset(cache, format="parquet", partitioning="hive")

    own = np.asarray(hp.lonlat_to_healpix(q_ra * u.deg, dec * u.deg), dtype=np.int64)
    by_pix = defaultdict(list)
    for i, p in enumerate(own):
        by_pix[int(p)].append(i)
    pixels = sorted(by_pix)

    # --- pixel selection ----------------------------------------------------
    if args.legacy_halo:
        extra_by_i = {}
        n_straddle = -1
        print("[USNOB][WARN] --legacy-halo: blanket 8-neighbour halo. A/B use only.",
              flush=True)
    else:
        _, extra_by_i, n_straddle = _pixels_for_points(
            hp, q_ra, dec, args.radius_arcsec, own)
        print(f"[PIX] {len(pixels):,} core pixels; {n_straddle:,}/{n:,} points "
              f"({100.0*n_straddle/max(n,1):.3f}%) lie within {args.radius_arcsec}\" "
              f"of a pixel boundary and pull in a neighbour", flush=True)

    counts = _pixel_row_counts(cache, rebuild=args.rebuild_pixel_counts)

    # --- batching, bounded by ROWS as well as by pixel count -----------------
    batches: list[list[int]] = []
    cur: list[int] = []
    cur_rows = 0
    for p in pixels:
        est = counts.get(p, 0)
        if cur and (len(cur) >= args.batch or cur_rows + est > args.max_batch_rows):
            batches.append(cur)
            cur, cur_rows = [], 0
        cur.append(p)
        cur_rows += est
    if cur:
        batches.append(cur)
    print(f"[BATCH] {len(batches)} batches (<= {args.batch} pixels, "
          f"<= {args.max_batch_rows:,} rows each)", flush=True)

    sep = np.full(n, np.inf)
    n_poss2_kept = n_poss2_dropped = 0
    P = xyz(q_ra, dec)
    n_zero_cov = 0
    peak_need = 0

    cols = ["ra", "dec"]
    if not args.allow_poss1_only:
        cols += ["B2mag", "R2mag"] + (["Imag"] if args.count_nir else [])
    if not args.no_pm:
        cols += ["pmRA", "pmDE"]

    for bi, core in enumerate(batches, 1):
        idx = np.concatenate([by_pix[p] for p in core])
        need = set(core)
        if args.legacy_halo:
            for p in core:
                need.update(int(x) for x in np.atleast_1d(hp.neighbours(p)) if x >= 0)
        else:
            for i in idx:
                e = extra_by_i.get(int(i))
                if e:
                    need.update(e)
        peak_need = max(peak_need, len(need))

        scanner = dset.scanner(columns=cols,
                               filter=pc.field("healpix_5").isin(sorted(need)),
                               batch_size=args.scan_batch_rows)

        n_ref_seen = 0
        for rb in scanner.to_batches():
            if rb.num_rows == 0:
                continue

            def col(name):
                return rb.column(name).to_numpy(zero_copy_only=False).astype(float)

            c_ra, c_dec = col("ra"), col("dec")

            # --- second-epoch restriction -----------------------------------
            # USNO-B is built partly FROM POSS-I, so matching against all of it
            # would veto our own detections against themselves. Its 2+-survey
            # merge rule means a POSS-I-only source never becomes an entry at
            # all, but that rule can also be satisfied by POSS-I O *and* POSS-I
            # E -- two 1950s Palomar emulsions, not an independent epoch.
            # Requiring POSS-II (B2/R2) is therefore the whole argument, not a
            # refinement. Missing magnitudes are stored as NaN, verified, so
            # isfinite() is a valid test here; a catalogue encoding them as
            # 99.99 would silently mark every row second-epoch and quietly
            # restore the circular veto.
            keep2 = None
            if not args.allow_poss1_only:
                has2 = np.isfinite(col("B2mag")) | np.isfinite(col("R2mag"))
                if args.count_nir:
                    has2 |= np.isfinite(col("Imag"))
                keep2 = np.flatnonzero(has2)
                n_poss2_kept += int(keep2.size)
                n_poss2_dropped += int(has2.size - keep2.size)
                if keep2.size == 0:
                    continue
                c_ra, c_dec = c_ra[keep2], c_dec[keep2]

            # --- proper motion to the plate epoch ---------------------------
            # USNO-B positions are epoch 2000.0; the plates are ~1950. Over ~47
            # yr a 100 mas/yr star moves 4.7", comparable to the whole veto
            # radius. The errors are not random either: Monet+2003 documents
            # magnitude-dependent fixed-pattern astrometric residuals on POSS-I.
            if not args.no_pm:
                pm_ra, pm_de = col("pmRA"), col("pmDE")
                if keep2 is not None:
                    pm_ra, pm_de = pm_ra[keep2], pm_de[keep2]
                dt = args.plate_epoch - 2000.0
                pm_ra = np.nan_to_num(pm_ra)
                pm_de = np.nan_to_num(pm_de)
                # pmRA is mas/yr and already includes cos(dec); undo it for a
                # delta in RA. Applying it without dividing shrinks the
                # correction by cos(dec) -- invisible near the equator, 10x
                # wrong near the pole.
                c_dec = c_dec + pm_de * dt / 3.6e6
                c_ra = c_ra + (pm_ra * dt / 3.6e6) / np.cos(
                    np.radians(np.clip(c_dec, -89.9, 89.9)))

            n_ref_seen += len(c_ra)
            tree = cKDTree(xyz(c_ra, c_dec))
            d, _ = tree.query(P[idx], k=1,
                              distance_upper_bound=chord(args.radius_arcsec))
            ok = np.isfinite(d)
            if ok.any():
                s = np.degrees(2.0 * np.arcsin(np.clip(d[ok] / 2.0, 0, 1))) * 3600.0
                cur_sep = sep[idx]
                cur_sep[ok] = np.minimum(cur_sep[ok], s)
                sep[idx] = cur_sep

        if n_ref_seen == 0:
            n_zero_cov += len(idx)
        if bi % 10 == 0 or bi == len(batches):
            print(f"  [batch {bi}/{len(batches)}  peak_need={peak_need} pixels]",
                  flush=True)

    matched = sep <= args.radius_arcsec
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"stage_{args.stage}_USNOB"
    kept_p, flags_p, led_p = (out_dir / f"{tag}.csv", out_dir / f"{tag}_flags.csv",
                              out_dir / f"{tag}_ledger.json")

    with kept_p.open("w", newline="", encoding="utf-8") as fk, \
            flags_p.open("w", newline="", encoding="utf-8") as ff:
        kw = csv.DictWriter(fk, fieldnames=["src_id", "ra", "dec"])
        fw = csv.DictWriter(ff, fieldnames=["src_id", "ra", "dec",
                                            "has_usnob_match", "best_sep_arcsec"])
        kw.writeheader(); fw.writeheader()
        for i, (sid, a, dcl) in enumerate(rows):
            fw.writerow({"src_id": sid, "ra": f"{a:.8f}", "dec": f"{dcl:.8f}",
                         "has_usnob_match": bool(matched[i]),
                         "best_sep_arcsec": "" if not np.isfinite(sep[i]) else f"{sep[i]:.3f}"})
            if not matched[i]:
                kw.writerow({"src_id": sid, "ra": f"{a:.8f}", "dec": f"{dcl:.8f}"})

    n_match = int(matched.sum())
    led = {
        "stage": "USNOB", "experimental": True, "wired_into_orchestrator": False,
        "input_csv": str(args.stage_csv), "cache": cache,
        "radius_arcsec": args.radius_arcsec, "null_shift_deg": args.null_shift,
        "hp_nside": args.nside, "hp_order": args.hp_order,
        # What this run actually vetoed against. A number quoted without these
        # is not interpretable: the POSS-II restriction is the difference
        # between an independent second-epoch veto and a circular one.
        "second_epoch": {
            "poss2_restricted": not args.allow_poss1_only,
            "nir_counted": bool(args.count_nir),
            "reference_rows_kept": n_poss2_kept,
            "reference_rows_dropped_poss1_only": n_poss2_dropped,
        },
        "proper_motion": {
            "propagated": not args.no_pm,
            "plate_epoch": None if args.no_pm else args.plate_epoch,
            "source_epoch": None if args.no_pm else 2000.0,
        },
        # Memory shape of the run, so a freeze is diagnosable after the fact.
        "reads": {
            "pixel_selection": "legacy_halo" if args.legacy_halo else "exact_per_point",
            "core_pixels": len(pixels),
            "batches": len(batches),
            "peak_pixels_per_batch": peak_need,
            "points_straddling_a_boundary": None if args.legacy_halo else n_straddle,
            "max_batch_rows": args.max_batch_rows,
            "scan_batch_rows": args.scan_batch_rows,
        },
        "counts": {"in_rows": n, "matched_rows": n_match,
                   "kept_rows": n - n_match, "zero_coverage_rows": n_zero_cov},
        "removal_pct": round(100.0 * n_match / n, 4),
        "outputs": {"kept_csv": str(kept_p), "flags_csv": str(flags_p),
                    "ledger_json": str(led_p)},
        "note": ("USNO-B is dense; interpret matched_rows against a --null-shift run. "
                 "poss2_restricted=false is CIRCULAR with POSS-I detections and must "
                 "not be used for a published veto."),
    }
    led_p.write_text(json.dumps(led, indent=2), encoding="utf-8")
    print(f"[USNOB] in={n} matched={n_match} ({100*n_match/n:.2f}%) "
          f"kept={n-n_match} zero_coverage={n_zero_cov}")
    print(f"[USNOB] peak pixels in one batch: {peak_need}")
    print(f"[USNOB] wrote {kept_p}\n[USNOB] wrote {flags_p}\n[USNOB] wrote {led_p}")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    raise SystemExit(main())
