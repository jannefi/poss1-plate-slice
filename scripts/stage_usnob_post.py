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
from astropy import units as u
from astropy_healpix import HEALPix
from scipy.spatial import cKDTree


def xyz(ra, dec):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def chord(arcsec: float) -> float:
    return 2.0 * math.sin(math.radians(arcsec / 3600.0) / 2.0)


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
    ap.add_argument("--batch", type=int, default=120)
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

    pix = np.asarray(hp.lonlat_to_healpix(q_ra * u.deg, dec * u.deg))
    by_pix = defaultdict(list)
    for i, p in enumerate(pix):
        by_pix[int(p)].append(i)
    pixels = sorted(by_pix)
    sep = np.full(n, np.inf)
    n_poss2_kept = n_poss2_dropped = 0
    P = xyz(q_ra, dec)
    n_zero_cov = 0

    for b0 in range(0, len(pixels), args.batch):
        core = pixels[b0:b0 + args.batch]
        need = set(core)
        for p in core:
            need.update(int(x) for x in np.atleast_1d(hp.neighbours(p)) if x >= 0)
        cols = ["ra", "dec"]
        if not args.allow_poss1_only:
            cols += ["B2mag", "R2mag"] + (["Imag"] if args.count_nir else [])
        if not args.no_pm:
            cols += ["pmRA", "pmDE"]
        tbl = dset.to_table(columns=cols,
                            filter=pc.field("healpix_5").isin(sorted(need)))
        idx = np.concatenate([by_pix[p] for p in core])
        if tbl.num_rows == 0:
            n_zero_cov += len(idx)
            continue

        def col(name):
            return tbl.column(name).to_numpy(zero_copy_only=False).astype(float)

        c_ra, c_dec = col("ra"), col("dec")

        # --- second-epoch restriction -------------------------------------
        # USNO-B is built partly FROM POSS-I, so matching against all of it
        # would veto our own detections against themselves. Its 2+-survey merge
        # rule means a POSS-I-only source never becomes an entry at all, but
        # that rule can also be satisfied by POSS-I O *and* POSS-I E -- two
        # 1950s Palomar emulsions, not an independent epoch. Requiring POSS-II
        # (B2/R2) is therefore the whole argument, not a refinement.
        # Missing magnitudes are stored as NaN, verified, so notna() is a valid
        # test here; a catalogue encoding them as 99.99 would silently mark
        # every row second-epoch and quietly restore the circular veto.
        if not args.allow_poss1_only:
            has2 = np.isfinite(col("B2mag")) | np.isfinite(col("R2mag"))
            if args.count_nir:
                has2 |= np.isfinite(col("Imag"))
            keep2 = np.flatnonzero(has2)
            if keep2.size == 0:
                n_zero_cov += len(idx)
                continue
            c_ra, c_dec = c_ra[keep2], c_dec[keep2]
            n_poss2_kept += int(keep2.size)
            n_poss2_dropped += int(has2.size - keep2.size)

        # --- proper motion to the plate epoch ------------------------------
        # USNO-B positions are epoch 2000.0; the plates are ~1950. Over ~47 yr a
        # 100 mas/yr star moves 4.7", comparable to the whole veto radius. The
        # errors are not random either: Monet+2003 documents magnitude-dependent
        # fixed-pattern astrometric residuals on POSS-I specifically.
        if not args.no_pm:
            pm_ra, pm_de = col("pmRA"), col("pmDE")
            if not args.allow_poss1_only:
                pm_ra, pm_de = pm_ra[keep2], pm_de[keep2]
            dt = args.plate_epoch - 2000.0
            pm_ra = np.nan_to_num(pm_ra)
            pm_de = np.nan_to_num(pm_de)
            # pmRA is mas/yr and already includes cos(dec); undo it for a
            # delta in RA. Applying it without dividing shrinks the correction
            # by cos(dec) -- invisible near the equator, 10x wrong near the pole.
            c_dec = c_dec + pm_de * dt / 3.6e6
            c_ra = c_ra + (pm_ra * dt / 3.6e6) / np.cos(np.radians(np.clip(c_dec, -89.9, 89.9)))

        tree = cKDTree(xyz(c_ra, c_dec))
        d, _ = tree.query(P[idx], k=1, distance_upper_bound=chord(args.radius_arcsec))
        ok = np.isfinite(d)
        s = np.degrees(2.0 * np.arcsin(np.clip(d[ok] / 2.0, 0, 1))) * 3600.0
        cur = sep[idx]
        cur[ok] = np.minimum(cur[ok], s)
        sep[idx] = cur
        done = b0 + len(core)
        if (b0 // args.batch) % 10 == 0 or done >= len(pixels):
            print(f"  [{done}/{len(pixels)} pixels]", flush=True)

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
    print(f"[USNOB] wrote {kept_p}\n[USNOB] wrote {flags_p}\n[USNOB] wrote {led_p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
