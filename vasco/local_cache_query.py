"""
Query local HEALPix-5 Parquet caches for Gaia DR3, PS1 DR2, and USNO-B1.0.

Used by fetch_gaia_neighbourhood, fetch_ps1_neighbourhood, and
fetch_usnob_neighbourhood as a drop-in replacement for VizieR when a
local cache is available.

Activation: set VASCO_GAIA_CACHE, VASCO_PS1_CACHE, or VASCO_USNOB_CACHE
to the root directory of the corresponding cache (e.g. /Volumes/SANDISK/Gaia).
When unset, the pipeline falls through to VizieR as before.

Each function writes a CSV with the same column names and format as the
VizieR-backed fetcher it replaces, so downstream code (xmatch, wcsfix,
analysis) is unaffected.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np

_HP = None  # lazy-loaded
_DS_CACHE: dict[str, object] = {}  # cache_dir → pyarrow.dataset.Dataset


def _get_hp():
    global _HP
    if _HP is None:
        from astropy_healpix import HEALPix
        _HP = HEALPix(nside=32, order="nested")
    return _HP


def _get_dataset(cache_dir: str):
    """Return a cached pyarrow Dataset for the given cache directory.

    Partition discovery (scanning 12k+ directories) is expensive on USB
    drives. Caching the Dataset object makes subsequent queries instant.
    """
    import pyarrow.dataset as pds
    if cache_dir not in _DS_CACHE:
        parquet_dir = Path(cache_dir) / "parquet"
        _DS_CACHE[cache_dir] = pds.dataset(
            str(parquet_dir), format="parquet", partitioning="hive"
        )
    return _DS_CACHE[cache_dir]


def _cone_query(cache_dir: str, ra: float, dec: float, radius_arcmin: float,
                columns: list[str], *, parquet_filter=None):
    """Load rows from a local HP5 cache within a cone.

    Returns a pandas DataFrame with the requested columns plus '_r' (sep in arcmin).
    """
    import astropy.units as u
    import pyarrow.compute as pc

    hp = _get_hp()
    ds = _get_dataset(cache_dir)

    pixels = hp.cone_search_lonlat(ra * u.deg, dec * u.deg,
                                   radius=radius_arcmin * u.arcmin)
    pixels = [int(p) for p in pixels.tolist()]

    filt = pc.field("healpix_5").isin(pixels)
    if parquet_filter is not None:
        filt = filt & parquet_filter

    tbl = ds.to_table(columns=columns, filter=filt)
    df = tbl.to_pandas()

    if len(df) == 0:
        df["_r"] = []
        return df

    ra_r = np.deg2rad(df["ra"].to_numpy())
    dec_r = np.deg2rad(df["dec"].to_numpy())
    cra, cdec = np.deg2rad(ra), np.deg2rad(dec)
    cos_sep = (np.sin(dec_r) * np.sin(cdec) +
               np.cos(dec_r) * np.cos(cdec) * np.cos(ra_r - cra))
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    sep_arcmin = np.rad2deg(np.arccos(cos_sep)) * 60.0

    df["_r"] = sep_arcmin
    df = df[df["_r"] <= radius_arcmin].copy()
    df = df.sort_values("_r").reset_index(drop=True)
    return df


# ── Gaia DR3 ─────────────────────────────────────────────────────────

def query_gaia(tile_dir: Path, ra: float, dec: float,
               radius_arcmin: float) -> Path | None:
    """Query local Gaia cache → write gaia_neighbourhood.csv.

    Returns the output path on success, None if VASCO_GAIA_CACHE is unset.
    """
    cache = os.getenv("VASCO_GAIA_CACHE")
    if not cache:
        return None

    out = Path(tile_dir) / "catalogs" / "gaia_neighbourhood.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    df = _cone_query(cache, ra, dec, radius_arcmin,
                     columns=["ra", "dec", "phot_g_mean_mag", "pmra", "pmdec"])

    with out.open("w", newline="") as f:
        w = csv.writer(f)
        # Match the VizieR fetcher's CSV header exactly:
        # ra, dec, Gmag, BPmag, RPmag, pmRA, pmDE, Plx, _r
        w.writerow(["ra", "dec", "Gmag", "BPmag", "RPmag", "pmRA", "pmDE", "Plx", "_r"])
        for _, r in df.iterrows():
            w.writerow([
                f"{r['ra']:.8f}",
                f"{r['dec']:.8f}",
                f"{r['phot_g_mean_mag']:.4f}" if not np.isnan(r["phot_g_mean_mag"]) else "",
                "",  # BPmag — not in cache
                "",  # RPmag — not in cache
                f"{r['pmra']:.3f}" if not np.isnan(r["pmra"]) else "",
                f"{r['pmdec']:.3f}" if not np.isnan(r["pmdec"]) else "",
                "",  # Plx — not in cache
                f"{r['_r']:.6f}",
            ])
    return out


# ── PS1 DR2 ──────────────────────────────────────────────────────────

def query_ps1(tile_dir: Path, ra: float, dec: float,
              radius_arcmin: float) -> Path | None:
    """Query local PS1 cache → write ps1_neighbourhood.csv.

    Applies nDetections >= 3 filter for VizieR II/389 compatibility.
    Returns the output path on success, None if VASCO_PS1_CACHE is unset.
    """
    import pyarrow.compute as pc

    cache = os.getenv("VASCO_PS1_CACHE")
    if not cache:
        return None

    out = Path(tile_dir) / "catalogs" / "ps1_neighbourhood.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    df = _cone_query(
        cache, ra, dec, radius_arcmin,
        columns=["objID", "ra", "dec", "nDetections", "gmag", "rmag", "imag", "zmag", "ymag"],
        parquet_filter=pc.field("nDetections") >= 3,
    )

    # PS1 cache stores -999 sentinel for non-detections; VizieR sends empty.
    # Convert -999 → empty for CSV output.
    mag_cols = ["gmag", "rmag", "imag", "zmag", "ymag"]

    with out.open("w", newline="") as f:
        w = csv.writer(f)
        # Match VizieR fetcher header: objID,ra,dec,Nd,gmag,rmag,imag,zmag,ymag,_r
        w.writerow(["objID", "ra", "dec", "Nd", "gmag", "rmag", "imag", "zmag", "ymag", "_r"])
        for _, r in df.iterrows():
            row = [
                str(int(r["objID"])),
                f"{r['ra']:.8f}",
                f"{r['dec']:.8f}",
                str(int(r["nDetections"])),
            ]
            for mc in mag_cols:
                v = r[mc]
                if np.isnan(v) or v <= -999:
                    row.append("")
                else:
                    row.append(f"{v:.4f}")
            row.append(f"{r['_r']:.6f}")
            w.writerow(row)
    return out


def query_ps1_bright_stars(ra: float, dec: float, radius_arcmin: float = 45.0,
                           rmag_max: float = 16.0, mindetections: int = 2):
    """Query local PS1 cache for bright stars → list[BrightStar].

    Local-cache equivalent of vasco.mnras.spikes.fetch_bright_ps1 (live
    MAST call): same PS1 DR2 catalog (build_cache.py maps rMeanPSFMag ->
    rmag, nDetections as-is 1:1), same nDetections/rmag filters, same
    valid-magnitude sanity check. The live call measured at ~31s/tile;
    this exists to avoid paying that on every tile (2026-07-27).

    Returns None if VASCO_PS1_CACHE is unset (caller falls through to
    the live fetch, matching every other function in this file).
    """
    import pyarrow.compute as pc
    from vasco.mnras.spikes import BrightStar, _is_valid_mag

    cache = os.getenv("VASCO_PS1_CACHE")
    if not cache:
        return None

    df = _cone_query(
        cache, ra, dec, radius_arcmin,
        columns=["ra", "dec", "nDetections", "rmag"],
        parquet_filter=pc.field("nDetections") >= mindetections,
    )

    out = []
    for _, r in df.iterrows():
        rmag = r["rmag"]
        if np.isnan(rmag) or rmag <= -999:
            continue
        rmag = float(rmag)
        if (not _is_valid_mag(rmag)) or (not (rmag <= rmag_max)):
            continue
        out.append(BrightStar(ra=float(r["ra"]), dec=float(r["dec"]), rmag=rmag))
    return out


# ── USNO-B1.0 ────────────────────────────────────────────────────────

def query_usnob(tile_dir: Path, ra: float, dec: float,
                radius_arcmin: float) -> Path | None:
    """Query local USNO-B cache → write usnob_neighbourhood.csv.

    Returns the output path on success, None if VASCO_USNOB_CACHE is unset.
    """
    cache = os.getenv("VASCO_USNOB_CACHE")
    if not cache:
        return None

    out = Path(tile_dir) / "catalogs" / "usnob_neighbourhood.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    df = _cone_query(
        cache, ra, dec, radius_arcmin,
        columns=["ra", "dec", "B1mag", "R1mag", "B2mag", "R2mag", "Imag", "pmRA", "pmDE"],
    )

    mag_cols = ["R1mag", "R2mag", "B1mag", "B2mag", "Imag"]
    pm_cols = ["pmRA", "pmDE"]

    with out.open("w", newline="") as f:
        w = csv.writer(f)
        # Match actual astroquery/VizieR USNO-B output exactly (confirmed from
        # prod CSV header in 5b48711): ra,dec,R1mag,R2mag,B1mag,B2mag,Imag,pmRA,pmDE,_r
        # No 'USNO-B1.0' designation column; _r is the separation from query centre.
        w.writerow(["ra", "dec", "R1mag", "R2mag", "B1mag", "B2mag", "Imag",
                    "pmRA", "pmDE", "_r"])
        for _, r in df.iterrows():
            row = [
                f"{r['ra']:.8f}",
                f"{r['dec']:.8f}",
            ]
            for mc in mag_cols:
                v = r[mc]
                if v is None or (isinstance(v, float) and np.isnan(v)):
                    row.append("")
                else:
                    row.append(f"{v:.2f}")
            for pc_ in pm_cols:
                v = r[pc_]
                if v is None or (isinstance(v, float) and np.isnan(v)):
                    row.append("")
                else:
                    row.append(str(int(v)))
            row.append(f"{r['_r']:.6f}")
            w.writerow(row)
    return out
