#!/usr/bin/env python3
"""Which pipeline stage removed each reference-catalogue row that we detected?

This pipeline finds nearly every source in a published reference catalogue at
the raw-detection stage, and keeps a much smaller fraction in the final
candidate list. This tool attributes that difference to the exact stage that
consumed each reference row -- veto by veto, gate by gate -- instead of
reasoning about it from aggregate counts.

EXACT, NOT RE-MATCHED PER STAGE. Every per-tile stage catalogue carries the
SExtractor `NUMBER` column, which is unique within a tile and strictly NESTED
down the chain: each later stage is a subset of the previous one. So a
reference row is matched to its nearest raw detection exactly ONCE, and every
later question -- "did it survive the PS1 veto?" -- is answered by set
membership on that NUMBER. Buckets therefore sum to the number of reference
rows in footprint by construction (asserted), with no radius ambiguity
accumulating stage by stage.

Pipeline order, and the file that proves each step:

    sextractor_pass2.wcsfix.csv                     raw (WCS-fixed)
      -> .after_gaia_veto.csv                       Gaia 5"
      -> .after_ps1_veto.csv                        PS1 5"
      -> .after_usnob_veto.csv                      USNO-B 5"
      -> .after_usnob_veto.rejected_extract.csv     FLAGS==0 & SNR_WIN>=30
      -> .after_usnob_veto.rejected_morphology.csv  FWHM/ELONG/SPREAD_MODEL/extent
      -> sextractor_spike_rejected.csv              bright-star spike mask
      -> .filtered.csv                              survivors
      -> (global dedup)                             -> stage_S0.csv

THREE THINGS TO READ BEFORE QUOTING THE OUTPUT

1. **The ordering caveat.** The vetoes run BEFORE the MNRAS filters, and the
   morphology filter's sigma-clip window is derived from the population being
   filtered. The attribution is therefore what the pipeline ACTUALLY did, in
   its actual order. It is NOT a counterfactual: "would the filters have cut
   this row anyway?" has no clean answer, because removing a veto changes the
   population and hence the clip window.
2. **NEVER_DETECTED is truncated unless you restrict to plate cores.** A
   reference row's detection often lives on an overlapping tile outside the
   sampled set, and reads here as undetected. Restrict to rows well inside a
   sampled plate before quoting any detection rate. The conditional
   attribution among detected rows is unaffected.
3. **Coordinates are RA_corr/Dec_corr** -- the WCS-fixed positions the veto
   matches on and that S0 carries. The tool refuses a tile whose raw catalogue
   lacks those columns rather than silently falling back to the plate WCS,
   which differs by up to tens of arcsec on some tiles.

A reference row inside more than one tile footprint (tessellation overlap) is
credited with its BEST outcome across tiles -- the furthest stage it reached
anywhere -- so overlap cannot manufacture a loss.

Usage:
    python3 tools/funnel_attribution.py \\
      --tiles-root <run>/tiles \\
      --ref-csv vanish_possi.csv \\
      --s0-csv results/s0-642-20260814/stage_S0.csv.gz \\
      --label R --out funnel_rows.csv --out-summary funnel.json

How to validate: the printed bucket table must sum to the in-footprint row
count (the tool asserts this and exits non-zero otherwise), and re-running
with a different --match-radius-arcsec should move NEVER_DETECTED without
reshuffling the relative shares of the later stages.
"""
from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# Stage ranks: higher == survived further. Used to reduce a reference row seen
# by several overlapping tiles to its single best outcome.
STAGES = [
    "NEVER_DETECTED",
    "VETOED_GAIA",
    "VETOED_PS1",
    "VETOED_USNOB",
    "MNRAS_EXTRACT",
    "MNRAS_MORPHOLOGY",
    "SPIKE_MASKED",
    "MNRAS_UNATTRIBUTED",
    "DEDUP_DROPPED",
    "SURVIVED_S0",
]
RANK = {s: i for i, s in enumerate(STAGES)}

CAT = "catalogs/"
FILES = {
    "raw": CAT + "sextractor_pass2.wcsfix.csv",
    "gaia": CAT + "sextractor_pass2.after_gaia_veto.csv",
    "ps1": CAT + "sextractor_pass2.after_ps1_veto.csv",
    "usnob": CAT + "sextractor_pass2.after_usnob_veto.csv",
    "rej_extract": CAT + "sextractor_pass2.after_usnob_veto.rejected_extract.csv",
    "rej_morph": CAT + "sextractor_pass2.after_usnob_veto.rejected_morphology.csv",
    "rej_spike": CAT + "sextractor_spike_rejected.csv",
    "filtered": CAT + "sextractor_pass2.filtered.csv",
}

# Reference catalogues in this field use several spellings for the same two
# columns; accept them rather than making the caller pre-rename a public file.
RA_NAMES = ("Ra", "ra", "RA", "RAJ2000", "ra_deg", "_RAJ2000")
DEC_NAMES = ("Dec", "dec", "DEC", "DEJ2000", "dec_deg", "_DEJ2000")


def read_reference(path: str) -> tuple[np.ndarray, np.ndarray]:
    """(ra, dec) arrays from a reference catalogue, column naming tolerated."""
    d = pd.read_csv(path)
    ra_col = next((c for c in RA_NAMES if c in d.columns), None)
    dec_col = next((c for c in DEC_NAMES if c in d.columns), None)
    if ra_col is None or dec_col is None:
        raise SystemExit(
            f"[FAIL] {path}: no recognised RA/Dec columns. Have "
            f"{list(d.columns)[:8]}; accept any of {RA_NAMES} / {DEC_NAMES}")
    d = d[[ra_col, dec_col]].dropna()
    return d[ra_col].to_numpy(float), d[dec_col].to_numpy(float)


def tile_center(tile_id: str) -> tuple[float, float]:
    ra = float(tile_id.split("_RA")[1].split("_")[0])
    d = tile_id.split("_DEC")[1]
    return ra, float(d[1:]) * (1.0 if d[0] == "p" else -1.0)


def _numbers(path: Path) -> set[int]:
    """NUMBER column as a set. Empty or absent file -> empty set, which is a
    real pipeline state: a stage that keeps nothing writes a 0-byte file."""
    try:
        if path.stat().st_size == 0:
            return set()
        return set(pd.read_csv(path, usecols=["NUMBER"]).NUMBER.astype(int))
    except (FileNotFoundError, ValueError, pd.errors.EmptyDataError):
        return set()


def _reasons(path: Path) -> dict[int, str]:
    """NUMBER -> reject_reason. The MNRAS reject files record which gate fired,
    so a filter loss is attributable to FWHM vs elongation vs SPREAD_MODEL vs
    the sigma clip without re-deriving any of them."""
    try:
        if path.stat().st_size == 0:
            return {}
        d = pd.read_csv(path, usecols=["NUMBER", "reject_reason"])
        return dict(zip(d.NUMBER.astype(int), d.reject_reason.astype(str)))
    except (FileNotFoundError, ValueError, pd.errors.EmptyDataError):
        return {}


def _load_raw(path: Path):
    """(NUMBER, RA_corr, Dec_corr) for the WCS-fixed raw catalogue."""
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        d = pd.read_csv(path, usecols=["NUMBER", "RA_corr", "Dec_corr"])
    except ValueError:
        return "NO_WCSFIX"      # refuse rather than use the plate WCS silently
    d = d.dropna(subset=["RA_corr", "Dec_corr"])
    if d.empty:
        return None
    return (d.NUMBER.to_numpy(int), d.RA_corr.to_numpy(float),
            d.Dec_corr.to_numpy(float))


def _proj(ra, dec, ra0, dec0):
    """Local tangent-plane arcsec offsets; exact enough over one ~1 deg tile.
    The RA difference is wrapped, so a tile at RA 0.3 works."""
    cosd = np.cos(np.deg2rad(dec0))
    dra = ((np.asarray(ra, float) - ra0 + 180.0) % 360.0 - 180.0) * cosd
    return np.column_stack([dra * 3600.0, (np.asarray(dec, float) - dec0) * 3600.0])


def do_tile(job):
    (tile_id, tiles_root, r_ra, r_dec, radius, foot_deg, s0_by_tile) = job
    root = Path(tiles_root) / tile_id
    ra0, dec0 = tile_center(tile_id)

    # Reference rows geometrically inside this tile. foot_deg=0.5 is the
    # INSCRIBED circle of a 60'x60' tile -- deliberately conservative, since a
    # larger radius pulls in rows that fall outside the pixels and scores them
    # as detection failures.
    cosd = np.cos(np.deg2rad(dec0))
    d = np.hypot(((r_ra - ra0 + 180.0) % 360.0 - 180.0) * cosd, r_dec - dec0)
    inside = np.flatnonzero(d <= foot_deg)
    if inside.size == 0:
        return tile_id, [], None

    raw = _load_raw(root / FILES["raw"])
    if raw == "NO_WCSFIX":
        return tile_id, [], "no wcsfix columns"
    if raw is None:
        return tile_id, [(int(i), "NEVER_DETECTED", float("nan"), -1, "")
                         for i in inside], "no raw detections"

    num, rra, rdec = raw
    tree = cKDTree(_proj(rra, rdec, ra0, dec0))
    dist, idx = tree.query(_proj(r_ra[inside], r_dec[inside], ra0, dec0),
                           distance_upper_bound=radius)

    sets = {k: _numbers(root / FILES[k])
            for k in ("gaia", "ps1", "usnob", "rej_extract", "rej_morph",
                      "rej_spike", "filtered")}
    why: dict[int, str] = {}
    for k in ("rej_extract", "rej_morph"):
        why.update(_reasons(root / FILES[k]))
    s0 = s0_by_tile.get(tile_id)
    s0_tree = cKDTree(_proj(s0[0], s0[1], ra0, dec0)) if s0 is not None else None

    out = []
    for k, ri in enumerate(inside):
        if not np.isfinite(dist[k]):
            out.append((int(ri), "NEVER_DETECTED", float("nan"), -1, ""))
            continue
        n = int(num[idx[k]])
        sep = float(dist[k])
        if n not in sets["gaia"]:
            st = "VETOED_GAIA"
        elif n not in sets["ps1"]:
            st = "VETOED_PS1"
        elif n not in sets["usnob"]:
            st = "VETOED_USNOB"
        elif n in sets["rej_extract"]:
            st = "MNRAS_EXTRACT"
        elif n in sets["rej_morph"]:
            st = "MNRAS_MORPHOLOGY"
        elif n in sets["rej_spike"]:
            st = "SPIKE_MASKED"
        elif n not in sets["filtered"]:
            st = "MNRAS_UNATTRIBUTED"
        elif s0_tree is None:
            st = "SURVIVED_S0"          # no S0 supplied; cannot see dedup
        else:
            dd, _ = s0_tree.query(_proj([r_ra[ri]], [r_dec[ri]], ra0, dec0),
                                  distance_upper_bound=radius)
            st = "SURVIVED_S0" if np.isfinite(dd[0]) else "DEDUP_DROPPED"
        out.append((int(ri), st, sep, n, why.get(n, "")))
    return tile_id, out, None


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiles-root", required=True,
                    help="run directory holding <tile_id>/catalogs/...")
    ap.add_argument("--ref-csv", required=True,
                    help="reference catalogue with RA/Dec columns")
    ap.add_argument("--s0-csv", default=None,
                    help="final dedup'd catalogue (src_id,tile_id,ra,dec). "
                         "Without it DEDUP_DROPPED cannot be separated from "
                         "SURVIVED_S0 and dedup losses read as survivors.")
    ap.add_argument("--tile-plate-map", default=None,
                    help="tile_id,plate_id -- adds plate_id to the row output")
    ap.add_argument("--match-radius-arcsec", type=float, default=5.0)
    ap.add_argument("--footprint-radius-deg", type=float, default=0.5)
    ap.add_argument("--workers", type=int, default=4,
                    help="local file reads only; modest values are kind to "
                         "spinning disks")
    ap.add_argument("--label", default="run")
    ap.add_argument("--out", default=None, help="per-reference-row CSV")
    ap.add_argument("--out-summary", default=None, help="summary JSON")
    args = ap.parse_args()

    # A stray pipeline env var does not affect this read-only tool, but it does
    # mean the tiles being read may not be the configuration you assume.
    for var in ("VASCO_CIRCLE_ARCMIN", "VASCO_REPRO_SINGLE_PASS",
                "VASCO_DISABLE_USNOB", "VASCO_WCSFIX_DISABLE"):
        if os.getenv(var):
            print(f"[CONFIG][WARN] {var}={os.getenv(var)} is set in this shell; "
                  f"the tiles you are reading may not be what you expect.")

    root = Path(args.tiles_root)
    tiles = sorted(p.name for p in root.iterdir() if p.is_dir())
    r_ra, r_dec = read_reference(args.ref_csv)
    print(f"[CONFIG] tiles={len(tiles)}  reference rows={len(r_ra)}  "
          f"match={args.match_radius_arcsec}\"  "
          f"footprint={args.footprint_radius_deg} deg")
    print("[CONFIG] coords=RA_corr/Dec_corr (WCS-fixed)")

    s0_by_tile: dict[str, tuple] = {}
    if args.s0_csv:
        s0 = pd.read_csv(args.s0_csv, usecols=["tile_id", "ra", "dec"])
        s0 = s0[s0.tile_id.isin(set(tiles))]
        for tid, g in s0.groupby("tile_id"):
            s0_by_tile[tid] = (g.ra.to_numpy(float), g.dec.to_numpy(float))
        print(f"[CONFIG] S0 rows on these tiles: {len(s0)} across "
              f"{len(s0_by_tile)} tiles")
    else:
        print("[CONFIG][WARN] no --s0-csv: dedup losses will read as SURVIVED_S0")

    best: dict[int, tuple] = {}
    notes: dict[str, str] = {}
    jobs = [(t, str(root), r_ra, r_dec, args.match_radius_arcsec,
             args.footprint_radius_deg, s0_by_tile) for t in tiles]
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(do_tile, j): j[0] for j in jobs}
        for f in as_completed(futs):
            tid, rows, note = f.result()
            if note:
                notes[tid] = note
            for ri, st, sep, num, reason in rows:
                rk = RANK[st]
                if ri not in best or rk > best[ri][0]:
                    best[ri] = (rk, st, sep, tid, num, reason)
            done += 1
            if done % 100 == 0:
                print(f"  {done}/{len(tiles)} tiles  rows seen={len(best)}",
                      flush=True)

    if not best:
        print("[FATAL] no reference rows inside any tile footprint")
        return 1

    df = pd.DataFrame(
        [{"ref_index": ri, "stage": st, "sep_arcsec": sep, "tile_id": tid,
          "det_number": num, "reject_reason": reason,
          "ra": r_ra[ri], "dec": r_dec[ri]}
         for ri, (_, st, sep, tid, num, reason) in best.items()]
    ).sort_values("ref_index")
    if args.tile_plate_map:
        m = pd.read_csv(args.tile_plate_map)
        col = "plate_id" if "plate_id" in m.columns else m.columns[1]
        df = df.merge(m[["tile_id", col]].rename(columns={col: "plate_id"}),
                      on="tile_id", how="left")

    n = len(df)
    counts = df.stage.value_counts()
    detected = n - int(counts.get("NEVER_DETECTED", 0))

    print(f"\n{'='*74}\n[{args.label}] reference rows in footprint: {n}")
    print(f"{'stage':24s} {'rows':>7s} {'% of all':>9s} {'% of detected':>14s}")
    print("-" * 74)
    for s in STAGES:
        c = int(counts.get(s, 0))
        if c == 0:
            continue
        print(f"{s:24s} {c:7d} {100.0*c/n:8.2f}% "
              + ("           n/a" if s == "NEVER_DETECTED"
                 else f"{100.0*c/detected:13.2f}%"))
    print("-" * 74)
    assert int(counts.sum()) == n, "buckets do not sum -- attribution is broken"
    print(f"{'TOTAL':24s} {n:7d} {100.0:8.2f}%")

    raw_recall = 100.0 * detected / n
    s0_recall = 100.0 * int(counts.get("SURVIVED_S0", 0)) / n
    veto = sum(int(counts.get(s, 0))
               for s in ("VETOED_GAIA", "VETOED_PS1", "VETOED_USNOB"))
    mnras = sum(int(counts.get(s, 0))
                for s in ("MNRAS_EXTRACT", "MNRAS_MORPHOLOGY", "SPIKE_MASKED",
                          "MNRAS_UNATTRIBUTED"))
    dedup = int(counts.get("DEDUP_DROPPED", 0))
    lost = detected - int(counts.get("SURVIVED_S0", 0))

    print(f"\n[{args.label}] raw detection recall : {raw_recall:6.2f}%  "
          f"(truncated unless restricted to plate cores -- see the docstring)")
    print(f"[{args.label}] S0 recall            : {s0_recall:6.2f}%")
    print(f"[{args.label}] gap to attribute     : {raw_recall - s0_recall:6.2f} "
          f"points ({lost} rows)\n")
    if lost:
        print("  share of the gap:")
        print(f"    veto chain (Gaia+PS1+USNO-B) {veto:7d}  {100.0*veto/lost:6.2f}%")
        print(f"    MNRAS filters + spikes       {mnras:7d}  {100.0*mnras/lost:6.2f}%")
        print(f"    global dedup                 {dedup:7d}  {100.0*dedup/lost:6.2f}%")

    reason_counts = (df[df.reject_reason.astype(str) != ""]
                     .reject_reason.value_counts().to_dict())
    if reason_counts:
        nr = sum(reason_counts.values())
        print("\n  which MNRAS gate fired (exact, from reject_reason):")
        for k, c in sorted(reason_counts.items(), key=lambda x: -x[1]):
            print(f"    {k:28s} {c:6d}  {100.0*c/nr:6.2f}%")

    summary = {
        "label": args.label, "tiles_root": str(root), "tiles": len(tiles),
        "match_radius_arcsec": args.match_radius_arcsec,
        "footprint_radius_deg": args.footprint_radius_deg,
        "ref_rows_total": int(len(r_ra)), "ref_rows_in_footprint": n,
        "raw_recall_pct": round(raw_recall, 4),
        "s0_recall_pct": round(s0_recall, 4),
        "gap_points": round(raw_recall - s0_recall, 4), "gap_rows": int(lost),
        "counts": {s: int(counts.get(s, 0)) for s in STAGES},
        "gap_share_pct": ({"veto_chain": round(100.0 * veto / lost, 4),
                           "mnras_filters_and_spikes": round(100.0 * mnras / lost, 4),
                           "dedup": round(100.0 * dedup / lost, 4)} if lost else {}),
        "mnras_reject_reasons": {k: int(v) for k, v in reason_counts.items()},
        "tile_notes": notes,
    }
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, index=False)
        print(f"\n[OUT] {args.out}")
    if args.out_summary:
        Path(args.out_summary).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out_summary).write_text(json.dumps(summary, indent=2))
        print(f"[OUT] {args.out_summary}")
    if notes:
        print(f"[NOTE] {len(notes)} tiles flagged: {sorted(set(notes.values()))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
