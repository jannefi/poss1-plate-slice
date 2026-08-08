"""
Benchmark the local Gaia Parquet cache against VizieR I/355/gaiadr3 cone
searches, on N random tiles from a tile-plan CSV. Measures correctness
(row counts, per-row position / magnitude / proper-motion deltas) and
speed (Parquet seconds vs VizieR seconds vs speedup factor).

Uses VASCO60's own fetch_gaia_neighbourhood so the VizieR-side query is
bit-identical to what the pipeline runs in production.

Defaults:
    --cache-dir / $VASCO_GAIA_CACHE / /Volumes/SANDISK/Gaia
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
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as pds
from astropy_healpix import HEALPix
import astropy.units as u

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))
from vasco.external_fetch_online import fetch_gaia_neighbourhood  # noqa: E402


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_GAIA_CACHE", "/Volumes/SANDISK/Gaia"))


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
    tbl = ds.to_table(
        columns=["source_id", "ra", "dec", "pmra", "pmdec", "phot_g_mean_mag"],
        filter=pc.field("healpix_5").isin(pixels),
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


def query_vizier(
    ra: float, dec: float, radius_arcmin: float, max_rows: int = 200_000
):
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory() as td:
        tile_dir = Path(td) / "tile"
        tile_dir.mkdir()
        out = fetch_gaia_neighbourhood(
            tile_dir, ra, dec, radius_arcmin, max_rows=max_rows, timeout=120.0
        )
        df = pd.read_csv(out)
    dt = time.perf_counter() - t0
    return df, dt


def compare(
    df_parq: pd.DataFrame, df_viz: pd.DataFrame, max_rows: int
) -> dict:
    out: dict = {
        "n_parquet": len(df_parq),
        "n_vizier": len(df_viz),
        "row_count_match": len(df_parq) == len(df_viz),
        "vizier_hit_cap": len(df_viz) >= max_rows,
    }
    if len(df_parq) == 0 or len(df_viz) == 0:
        out["match_ok"] = len(df_parq) == len(df_viz)
        return out

    tol_deg = 1e-6
    p_ra = df_parq["ra"].to_numpy()
    p_dec = df_parq["dec"].to_numpy()
    p_idx_sort = np.argsort(p_ra)
    p_ra_sorted = p_ra[p_idx_sort]

    matched_parquet = np.zeros(len(df_parq), dtype=bool)
    unmatched_vizier = 0
    max_posdiff_mas = 0.0
    max_gdiff = 0.0
    max_pmradiff = 0.0
    max_pmdediff = 0.0

    v_ra = df_viz["ra"].to_numpy()
    v_dec = df_viz["dec"].to_numpy()
    v_g = df_viz.get("Gmag")
    v_pmra = df_viz.get("pmRA")
    v_pmde = df_viz.get("pmDE")

    for i in range(len(df_viz)):
        ra_i = v_ra[i]
        dec_i = v_dec[i]
        tol_ra = tol_deg / max(np.cos(np.deg2rad(dec_i)), 1e-6)
        lo = np.searchsorted(p_ra_sorted, ra_i - tol_ra)
        hi = np.searchsorted(p_ra_sorted, ra_i + tol_ra)
        candidates = p_idx_sort[lo:hi]
        best_j = -1
        best_sep = float("inf")
        for j in candidates:
            dra = (p_ra[j] - ra_i) * np.cos(np.deg2rad(dec_i))
            ddec = p_dec[j] - dec_i
            s = np.hypot(dra, ddec)
            if s < best_sep:
                best_sep = s
                best_j = j
        if best_j < 0 or best_sep > tol_deg:
            unmatched_vizier += 1
            continue
        matched_parquet[best_j] = True
        max_posdiff_mas = max(max_posdiff_mas, best_sep * 3.6e6)
        if v_g is not None and not pd.isna(v_g.iloc[i]):
            g_parq = df_parq["phot_g_mean_mag"].iloc[best_j]
            if not pd.isna(g_parq):
                max_gdiff = max(max_gdiff, abs(float(v_g.iloc[i]) - float(g_parq)))
        if v_pmra is not None and not pd.isna(v_pmra.iloc[i]):
            pmra_p = df_parq["pmra"].iloc[best_j]
            if not pd.isna(pmra_p):
                max_pmradiff = max(
                    max_pmradiff, abs(float(v_pmra.iloc[i]) - float(pmra_p))
                )
        if v_pmde is not None and not pd.isna(v_pmde.iloc[i]):
            pmde_p = df_parq["pmdec"].iloc[best_j]
            if not pd.isna(pmde_p):
                max_pmdediff = max(
                    max_pmdediff, abs(float(v_pmde.iloc[i]) - float(pmde_p))
                )

    out["matched"] = int(matched_parquet.sum())
    out["unmatched_vizier"] = unmatched_vizier
    out["unmatched_parquet"] = int((~matched_parquet).sum())
    out["max_posdiff_mas"] = max_posdiff_mas
    out["max_Gmag_diff"] = max_gdiff
    out["max_pmRA_diff_masyr"] = max_pmradiff
    out["max_pmDE_diff_masyr"] = max_pmdediff
    out["match_ok"] = unmatched_vizier == 0 and (
        out["unmatched_parquet"] == 0 or out["vizier_hit_cap"]
    )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark the local Gaia Parquet cache vs VizieR cone queries."
    )
    parser.add_argument("--cache-dir", type=Path, default=default_cache_dir())
    parser.add_argument("--tiles-csv", type=Path, default=default_tiles_csv())
    parser.add_argument("--n-tiles", type=int, default=10)
    parser.add_argument("--radius-arcmin", type=float, default=31.0)
    parser.add_argument("--seed", type=int, default=20260411)
    parser.add_argument("--max-rows", type=int, default=200_000)
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
        f"{'speedup':>8} {'pos_mas':>9} {'Gdiff':>7} {'ok'}"
    )
    for _, row in sample.iterrows():
        tile_id = row["tile_id"]
        ra = float(row["ra_deg"])
        dec = float(row["dec_deg"])

        df_p, dt_p, n_pixels = query_parquet(ds, hp, ra, dec, args.radius_arcmin)
        df_v, dt_v = query_vizier(ra, dec, args.radius_arcmin, args.max_rows)
        cmp_res = compare(df_p, df_v, args.max_rows)
        total_parq += dt_p
        total_viz += dt_v
        speedup = dt_v / dt_p if dt_p > 0 else float("inf")
        ok = "OK" if cmp_res["match_ok"] else "FAIL"
        print(
            f"{tile_id:40s} {ra:8.3f} {dec:8.3f} "
            f"{cmp_res['n_parquet']:8d} {cmp_res['n_vizier']:8d} "
            f"{dt_p:7.3f} {dt_v:7.3f} {speedup:7.1f}x "
            f"{cmp_res.get('max_posdiff_mas', 0):9.4f} "
            f"{cmp_res.get('max_Gmag_diff', 0):7.4f} "
            f"{ok}"
            + (" [CAPPED]" if cmp_res["vizier_hit_cap"] else "")
        )
        rows.append(
            {
                "tile_id": tile_id,
                **cmp_res,
                "parq_s": dt_p,
                "viz_s": dt_v,
                "n_hp_pixels": n_pixels,
            }
        )

    print()
    print(
        f"Total Parquet time: {total_parq:7.3f} s  "
        f"({total_parq / args.n_tiles:.3f} s/tile avg)"
    )
    print(
        f"Total VizieR  time: {total_viz:7.3f} s  "
        f"({total_viz / args.n_tiles:.3f} s/tile avg)"
    )
    if total_parq > 0:
        print(f"Overall speedup   : {total_viz / total_parq:.1f}x")

    n_ok = sum(1 for r in rows if r["match_ok"])
    n_capped = sum(1 for r in rows if r["vizier_hit_cap"])
    print(f"Match: {n_ok}/{len(rows)} tiles ok, {n_capped} hit VizieR cap")

    worst_pos = max(r.get("max_posdiff_mas", 0) for r in rows)
    worst_g = max(r.get("max_Gmag_diff", 0) for r in rows)
    worst_pmra = max(r.get("max_pmRA_diff_masyr", 0) for r in rows)
    worst_pmde = max(r.get("max_pmDE_diff_masyr", 0) for r in rows)
    print("Max per-row deltas across all matches:")
    print(f"  position : {worst_pos:.4f} mas")
    print(f"  Gmag     : {worst_g:.4f} mag")
    print(f"  pmRA     : {worst_pmra:.4f} mas/yr")
    print(f"  pmDE     : {worst_pmde:.4f} mas/yr")
    return 0 if n_ok == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
