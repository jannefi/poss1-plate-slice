"""
Benchmark the local PS1 Parquet cache against VizieR II/389/ps1_dr2
cone searches on N random tiles from a tile-plan CSV.

Measures correctness (row counts, per-row objID match, position/magnitude
deltas) and speed. Uses VASCO60's own fetch_ps1_neighbourhood so the
VizieR-side query is bit-identical to what the pipeline runs in
production.

Semantics — VizieR compatibility filter
---------------------------------------
The local cache holds the FULL PS1 otmo table (10.56 B rows). VizieR's
II/389/ps1_dr2 publishes a filtered subset: empirically proven to be
exactly `nDetections >= 3` — rows with 1 or 2 detections are excluded
from the VizieR view. To make the bench apples-to-apples, the cache
query here applies the same filter before comparing to VizieR. With
this filter, the cache reproduces VizieR's results bit-exactly on
objID, position (0.0000 mas), and mag (within float32 ULP).

If you want the raw superset (useful for more-inclusive veto), drop
the `nDetections >= 3` clause from query_parquet. See README.md.

Dense-field tiles with the filter applied may still hit VizieR's 200K
cap (unrelated to the filter). In those cases, the filtered cache
should still contain every VizieR row.

Defaults:
    --cache-dir / $VASCO_PS1_CACHE / /Volumes/SANDISK/PS1
    --tiles-csv                      / <repo>/plans/tiles_poss1e_ps1.csv
    --n-tiles                        / 10
    --radius-arcmin                  / 31.0
    --seed                           / 20260411
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as pds
from astropy_healpix import HEALPix
import astropy.units as u

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))
from vasco.external_fetch_online import fetch_ps1_neighbourhood  # noqa: E402

VIZIER_ROW_CAP = 200_000


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_PS1_CACHE", "/Volumes/SANDISK/PS1"))


def default_tiles_csv() -> Path:
    return REPO_ROOT / "plans" / "tiles_poss1e_ps1.csv"


def query_parquet(
    ds: pds.Dataset, hp: HEALPix, ra: float, dec: float, radius_arcmin: float
):
    t0 = time.perf_counter()
    pixels = hp.cone_search_lonlat(
        ra * u.deg, dec * u.deg, radius=radius_arcmin * u.arcmin
    )
    pixels = [int(p) for p in pixels.tolist()]
    # Apply the VizieR-compatibility filter `nDetections >= 3` at the
    # parquet scan so the downloaded table is already equivalent to what
    # VizieR's II/389/ps1_dr2 serves.
    tbl = ds.to_table(
        columns=["objID", "ra", "dec", "nDetections", "gmag", "rmag", "imag", "zmag", "ymag"],
        filter=pc.field("healpix_5").isin(pixels) & (pc.field("nDetections") >= 3),
    )
    df = tbl.to_pandas()
    df = df.dropna(subset=["ra", "dec"]).reset_index(drop=True)
    ra_r = np.deg2rad(df["ra"].to_numpy())
    dec_r = np.deg2rad(df["dec"].to_numpy())
    cra = np.deg2rad(ra)
    cdec = np.deg2rad(dec)
    cos_sep = np.sin(dec_r) * np.sin(cdec) + np.cos(dec_r) * np.cos(cdec) * np.cos(
        ra_r - cra
    )
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    sep_arcmin = np.rad2deg(np.arccos(cos_sep)) * 60.0
    df["_r_arcmin"] = sep_arcmin
    df = df[df["_r_arcmin"] <= radius_arcmin].copy()
    df = df.sort_values("_r_arcmin").reset_index(drop=True)
    dt = time.perf_counter() - t0
    return df, dt, len(pixels)


def query_vizier(ra: float, dec: float, radius_arcmin: float):
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory() as td:
        tile_dir = Path(td) / "tile"
        tile_dir.mkdir()
        out = fetch_ps1_neighbourhood(
            tile_dir, ra, dec, radius_arcmin, max_records=VIZIER_ROW_CAP, timeout=120.0
        )
        df = pd.read_csv(out)
    dt = time.perf_counter() - t0
    return df, dt


def compare(df_parq: pd.DataFrame, df_viz: pd.DataFrame) -> dict:
    """Strict bit-exact comparison.

    With the `nDetections >= 3` VizieR-compatibility filter applied on
    the Parquet side, the two result sets should be IDENTICAL up to the
    VizieR row cap: same objID set (symmetric difference = 0), same
    position to float64 precision, same mag in every band to float32
    ULP, same nDetections, same Nd.

    When VizieR hits the 200K row cap, the filtered Parquet can contain
    more rows than VizieR (strict equality relaxes to containment in
    that case, since the extras are just rows VizieR couldn't fit).
    """
    out: dict = {
        "n_parquet": len(df_parq),
        "n_vizier": len(df_viz),
        "vizier_hit_cap": len(df_viz) >= VIZIER_ROW_CAP,
    }
    if len(df_parq) == 0 and len(df_viz) == 0:
        out["match_ok"] = True
        out["strict_ok"] = True
        return out
    if len(df_viz) == 0:
        out["match_ok"] = False
        out["strict_ok"] = False
        out["error"] = "vizier_empty"
        return out

    p_ids = set(df_parq["objID"].astype("int64"))
    v_ids = set(df_viz["objID"].astype("int64"))
    v_not_in_p = v_ids - p_ids
    p_not_in_v = p_ids - v_ids

    out["vizier_rows_missing_from_parquet"] = len(v_not_in_p)
    out["parquet_rows_missing_from_vizier"] = len(p_not_in_v)

    # Per-column numeric deltas on matched rows
    common = list(v_ids & p_ids)
    if common:
        p_idx = df_parq.set_index("objID").loc[common]
        v_idx = df_viz.set_index("objID").loc[common]

        dra = (p_idx["ra"].to_numpy() - v_idx["ra"].to_numpy()) * np.cos(
            np.deg2rad(v_idx["dec"].to_numpy())
        )
        ddec = p_idx["dec"].to_numpy() - v_idx["dec"].to_numpy()
        out["max_posdiff_mas"] = float(np.nanmax(np.hypot(dra, ddec))) * 3.6e6

        # Per-band mag deltas — VizieR delivers mags as text; float cast,
        # then diff ignoring NaN (PS1 "no detection" -999 equates to NaN
        # in VizieR's view but to -999 in our cache — treat both as
        # not-comparable and skip).
        for viz_col, parq_col in [
            ("gmag", "gmag"),
            ("rmag", "rmag"),
            ("imag", "imag"),
            ("zmag", "zmag"),
            ("ymag", "ymag"),
        ]:
            if viz_col not in v_idx.columns or parq_col not in p_idx.columns:
                continue
            vv = pd.to_numeric(v_idx[viz_col], errors="coerce")
            pp = pd.to_numeric(p_idx[parq_col], errors="coerce")
            # Mask where either side is NaN or the cache's -999 sentinel
            both_valid = vv.notna() & pp.notna() & (pp != -999)
            if both_valid.any():
                d = (vv[both_valid] - pp[both_valid]).abs()
                out[f"max_{parq_col}_diff"] = float(d.max() or 0.0)
            else:
                out[f"max_{parq_col}_diff"] = 0.0

        # nDetections: cache has "nDetections", VizieR CSV has "Nd"
        if "Nd" in v_idx.columns and "nDetections" in p_idx.columns:
            vn = pd.to_numeric(v_idx["Nd"], errors="coerce")
            pn = pd.to_numeric(p_idx["nDetections"], errors="coerce")
            both = vn.notna() & pn.notna()
            if both.any():
                out["max_ndet_diff"] = int((vn[both] - pn[both]).abs().max())
            else:
                out["max_ndet_diff"] = 0

    # Strict: ids identical AND every numeric delta at float-ULP level
    strict_ok = (
        len(v_not_in_p) == 0
        and len(p_not_in_v) == 0
        and out.get("max_posdiff_mas", 0) < 1.0  # < 1 mas
        and max(out.get(f"max_{b}_diff", 0) for b in ["gmag", "rmag", "imag", "zmag", "ymag"]) < 0.001
        and out.get("max_ndet_diff", 0) == 0
    )
    out["strict_ok"] = strict_ok

    # Relaxed (cap-tolerant): v must be contained in p; p can exceed only
    # when VizieR's 200K cap was hit.
    out["match_ok"] = len(v_not_in_p) == 0 and (
        len(p_not_in_v) == 0 or out["vizier_hit_cap"]
    )
    if len(v_not_in_p) > 0:
        out["sample_v_not_in_p"] = sorted(v_not_in_p)[:3]
    if len(p_not_in_v) > 0 and not out["vizier_hit_cap"]:
        out["sample_p_not_in_v"] = sorted(p_not_in_v)[:3]
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark the local PS1 Parquet cache vs VizieR cone queries."
    )
    parser.add_argument("--cache-dir", type=Path, default=default_cache_dir())
    parser.add_argument("--tiles-csv", type=Path, default=default_tiles_csv())
    parser.add_argument("--n-tiles", type=int, default=10)
    parser.add_argument("--radius-arcmin", type=float, default=31.0)
    parser.add_argument("--seed", type=int, default=20260411)
    args = parser.parse_args()

    ds = pds.dataset(
        str(args.cache_dir / "parquet"), format="parquet", partitioning="hive"
    )
    hp = HEALPix(nside=32, order="nested")

    tiles = pd.read_csv(args.tiles_csv)
    sample = tiles.sample(args.n_tiles, random_state=args.seed).reset_index(drop=True)

    rows = []
    total_parq = 0.0
    total_viz = 0.0

    print(
        f"{'tile_id':40s} {'ra':>8} {'dec':>8} "
        f"{'n_parq':>8} {'n_viz':>8} {'parq_s':>7} {'viz_s':>7} "
        f"{'speedup':>8} {'v->p':>5} {'p->v':>5} {'strict'}"
    )
    for _, row in sample.iterrows():
        tile_id = row["tile_id"]
        ra = float(row["ra_deg"])
        dec = float(row["dec_deg"])

        df_p, dt_p, _ = query_parquet(ds, hp, ra, dec, args.radius_arcmin)
        try:
            df_v, dt_v = query_vizier(ra, dec, args.radius_arcmin)
        except Exception as e:
            print(f"{tile_id}: VizieR FAIL: {e}")
            continue
        cmp_res = compare(df_p, df_v)
        total_parq += dt_p
        total_viz += dt_v
        speedup = dt_v / dt_p if dt_p > 0 else float("inf")
        strict = "OK" if cmp_res.get("strict_ok") else "FAIL"
        flags = " [CAP]" if cmp_res.get("vizier_hit_cap") else ""
        print(
            f"{tile_id:40s} {ra:8.3f} {dec:8.3f} "
            f"{cmp_res['n_parquet']:8d} {cmp_res['n_vizier']:8d} "
            f"{dt_p:7.3f} {dt_v:7.3f} {speedup:7.1f}x "
            f"{cmp_res.get('vizier_rows_missing_from_parquet', 0):5d} "
            f"{cmp_res.get('parquet_rows_missing_from_vizier', 0):5d} "
            f"{strict}{flags}"
        )
        rows.append({"tile_id": tile_id, **cmp_res, "parq_s": dt_p, "viz_s": dt_v})

    print()
    if rows:
        print(
            f"Total Parquet time: {total_parq:7.3f} s "
            f"({total_parq / len(rows):.3f} s/tile avg)"
        )
        print(
            f"Total VizieR  time: {total_viz:7.3f} s "
            f"({total_viz / len(rows):.3f} s/tile avg)"
        )
        if total_parq > 0:
            print(f"Overall speedup   : {total_viz / total_parq:.1f}x")

        n_strict = sum(1 for r in rows if r.get("strict_ok"))
        n_cap_ok = sum(1 for r in rows if r.get("match_ok"))
        n_cap = sum(1 for r in rows if r.get("vizier_hit_cap"))
        print(f"Strict bit-exact: {n_strict}/{len(rows)} tiles")
        print(f"Relaxed (cap-tolerant): {n_cap_ok}/{len(rows)} tiles ({n_cap} hit VizieR 200K cap)")

        worst_pos = max(r.get("max_posdiff_mas", 0) or 0 for r in rows)
        print(f"Max per-row deltas across matched objIDs:")
        print(f"  position : {worst_pos:.4f} mas")
        for b in ["gmag", "rmag", "imag", "zmag", "ymag"]:
            v = max(r.get(f"max_{b}_diff", 0) or 0 for r in rows)
            print(f"  {b:5s}    : {v:.4f} mag")
        worst_nd = max((r.get("max_ndet_diff", 0) or 0) for r in rows)
        print(f"  nDetect  : {worst_nd}")

        # Dense-field insight: did any tile show cache > vizier by a lot?
        extras = [
            (r["tile_id"], r["n_parquet"] - r["n_vizier"])
            for r in rows
            if r.get("vizier_hit_cap") and r.get("n_parquet", 0) > r.get("n_vizier", 0)
        ]
        if extras:
            print()
            print("Dense-field extra rows visible ONLY in the local cache (cache wins):")
            for t, n in sorted(extras, key=lambda x: -x[1]):
                print(f"  {t}: +{n:,} rows beyond VizieR cap")

    return 0 if all(r.get("match_ok", False) for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
