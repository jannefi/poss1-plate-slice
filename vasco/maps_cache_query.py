"""
Query the local MAPS (Minnesota Automated Plate Scanner) POSS-I catalog
mirror for a positional cross-match/veto.

Unlike vasco/local_cache_query.py's Gaia/PS1/USNO-B caches (built by this
project as a global HEALPix-partitioned Parquet dataset), the MAPS mirror
is *already* naturally chunked by POSS-I plate -- 632 per-plate Parquet
files plus a small bounding-box index. There's no benefit to repartitioning
~90M rows into HEALPix just to match that convention; instead this module
prunes to the 1-3 overlapping plates via the bbox index, then reads only
those plate files. This is a deliberate choice: an earlier ad hoc attempt
to query this mirror with DuckDB SQL directly (no bbox pre-filtering)
repeatedly ran out of memory -- almost certainly from scanning/joining
across all 632 plate files at once. Never do a full-catalog scan here.

Activation: set VASCO_MAPS_CACHE to the root directory of the mirror
(e.g. <maps_cache>), which must
contain `maps_plate_index.parquet` and `parquet_icrs_by_plate/P###.parquet`.
When unset, query_maps() returns None (matching every function in
local_cache_query.py).

Format reference: README.data_format.md in the cache root (fetched from
https://aps.umn.edu/MAPS/README.data_format). Key facts this module
relies on:
  - flag = 100*Eduplicates + 10*Oflag + Eflag; Oflag/Eflag 0 = OK, 1-6 =
    various defects. A clean reference set is flag % 100 == 0 (Oflag ==
    Eflag == 0) -- NOT flag == 0 outright: empirically, 99.3% of objects
    have Eduplicates == 1 (checked against P004, 2026-07-27), so requiring
    flag == 0 would reject nearly the whole catalog. The README's "set to
    zero" phrasing for Eduplicates does not match the real data.
  - galnodO_x1000 > 500 -> galaxy, < 500 -> star (ANN classifier,
    O-plate preferred).
  - magdO_x1000 (diameter-magnitude) is the star convention, magiO_x1000
    (integrated) is the galaxy convention.
  - P003 and P926 are absent from the catalog entirely -- a query whose
    only overlapping plate would have been one of those instead finds
    zero overlapping plates. This is a data-availability gap, not the
    same thing as "plate(s) found, zero objects within radius" -- callers
    that need to distinguish the two should check `overlapping_plates` in
    the returned dict-like result (see query_maps docstring).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

_INDEX_CACHE: dict[str, "pd.DataFrame"] = {}
_PLATE_TABLE_CACHE: dict[tuple[str, int], "pd.DataFrame"] = {}

_OBJ_COLUMNS = [
    "POSS_field", "starnumO", "ra_icrs_deg", "dec_icrs_deg",
    "galnodO_x1000", "magdO_x1000", "magiO_x1000", "flag",
]


def _plate_index(cache_dir: str) -> "pd.DataFrame":
    """Load and cache maps_plate_index.parquet (632 rows, trivial size)."""
    if cache_dir not in _INDEX_CACHE:
        import pyarrow.parquet as pq
        idx_path = Path(cache_dir) / "maps_plate_index.parquet"
        _INDEX_CACHE[cache_dir] = pq.read_table(str(idx_path)).to_pandas()
    return _INDEX_CACHE[cache_dir]


def _plates_overlapping(cache_dir: str, ra: float, dec: float,
                        radius_deg: float) -> list[int]:
    """Return POSS_field numbers whose bounding box overlaps the query cone.

    Empirically verified (2026-07-27) against all 632 plates in
    maps_plate_index.parquet: every plate's ra_min_deg < ra_max_deg (no
    plate is itself stored wrapped across 0/360 deg) -- near-pole plates
    genuinely span nearly the full 0-360 deg RA range (real POSS-I
    circumpolar field footprints, not an index bug), so a plain interval
    overlap test against ra_min_deg/ra_max_deg is correct for the plate
    side. The only wraparound handling needed is for the tiny query
    margin itself, when the query point sits within `radius_deg` of
    0/360 deg.
    """
    idx = _plate_index(cache_dir)

    dec_ok = (idx["dec_max_deg"] >= dec - radius_deg) & (idx["dec_min_deg"] <= dec + radius_deg)

    # RA half-width grows toward the pole (great-circle radius -> RA extent).
    cosd = max(np.cos(np.deg2rad(dec)), 0.05)
    ra_radius_deg = min(radius_deg / cosd, 180.0)
    ra_lo = (ra - ra_radius_deg) % 360.0
    ra_hi = (ra + ra_radius_deg) % 360.0

    if ra_lo <= ra_hi:
        ra_ok = (idx["ra_max_deg"] >= ra_lo) & (idx["ra_min_deg"] <= ra_hi)
    else:
        # Query margin itself straddles 0/360 deg: overlap if the plate
        # touches either of the two wrapped pieces [ra_lo, 360) or [0, ra_hi].
        ra_ok = (idx["ra_max_deg"] >= ra_lo) | (idx["ra_min_deg"] <= ra_hi)

    return idx.loc[dec_ok & ra_ok, "POSS_field"].tolist()


def _plate_table(cache_dir: str, field: int) -> "pd.DataFrame | None":
    """Load (and cache) one plate's object table, columns-projected."""
    key = (cache_dir, field)
    if key not in _PLATE_TABLE_CACHE:
        import pyarrow.parquet as pq
        path = Path(cache_dir) / "parquet_icrs_by_plate" / f"P{field:03d}.parquet"
        if not path.exists():
            _PLATE_TABLE_CACHE[key] = None
        else:
            _PLATE_TABLE_CACHE[key] = pq.read_table(
                str(path), columns=_OBJ_COLUMNS
            ).to_pandas()
    return _PLATE_TABLE_CACHE[key]


def query_maps(ra: float, dec: float, radius_arcsec: float = 5.0, *,
              flag_ok_only: bool = True) -> "pd.DataFrame | None":
    """Query the local MAPS mirror within radius_arcsec of (ra, dec).

    Returns None if VASCO_MAPS_CACHE is unset. Otherwise returns a
    DataFrame (possibly empty) with columns:
        ra, dec (ICRS deg), galnodO_x1000, magdO, magiO (real magnitudes,
        already /1000), flag, POSS_field, starnumO, sep_arcsec
    sorted by sep_arcsec ascending.

    If flag_ok_only is True (default), rows are restricted to
    flag % 100 == 0 -- i.e. Oflag == 0 and Eflag == 0 (no moments/
    classifier error, not a scratch, not stripe-clipped, has background
    coverage, non-negative sky) -- before the radius filter is applied.
    The Eduplicates component (flag // 100) is deliberately ignored here:
    empirically (checked against P004, 2026-07-27) 99.3% of objects have
    Eduplicates == 1, not 0 -- README.data_format.md's "otherwise it is
    set to zero" phrasing does not match the real data, so Eduplicates==1
    is treated as the normal/expected case, not a defect.

    Distinguishing "no MAPS coverage here" from "MAPS covered it, no
    match": call _plates_overlapping directly, or check whether the
    returned DataFrame is empty *and* no plate file existed for any
    overlapping POSS_field (P003/P926 gaps) vs. genuinely zero objects
    within radius on a plate that does exist. Callers that need this
    distinction (e.g. stage_maps_post.py's ledger) should call
    plates_overlapping_count() alongside this function.
    """
    cache = os.getenv("VASCO_MAPS_CACHE")
    if not cache:
        return None

    radius_deg = radius_arcsec / 3600.0
    fields = _plates_overlapping(cache, ra, dec, radius_deg)

    frames = []
    for field in fields:
        tbl = _plate_table(cache, field)
        if tbl is None or len(tbl) == 0:
            continue
        frames.append(tbl)

    if not frames:
        return pd.DataFrame(columns=[
            "ra", "dec", "galnodO_x1000", "magdO", "magiO", "flag",
            "POSS_field", "starnumO", "sep_arcsec",
        ])

    df = pd.concat(frames, ignore_index=True)

    if flag_ok_only:
        df = df[df["flag"] % 100 == 0]

    if len(df) == 0:
        return df.rename(columns={"ra_icrs_deg": "ra", "dec_icrs_deg": "dec"})

    ra_r = np.deg2rad(df["ra_icrs_deg"].to_numpy())
    dec_r = np.deg2rad(df["dec_icrs_deg"].to_numpy())
    cra, cdec = np.deg2rad(ra), np.deg2rad(dec)
    cos_sep = (np.sin(dec_r) * np.sin(cdec) +
              np.cos(dec_r) * np.cos(cdec) * np.cos(ra_r - cra))
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    sep_arcsec = np.rad2deg(np.arccos(cos_sep)) * 3600.0

    df = df.assign(sep_arcsec=sep_arcsec)
    df = df[df["sep_arcsec"] <= radius_arcsec].copy()

    df["magdO"] = df["magdO_x1000"] / 1000.0
    df["magiO"] = df["magiO_x1000"] / 1000.0
    df = df.rename(columns={"ra_icrs_deg": "ra", "dec_icrs_deg": "dec"})
    df = df[["ra", "dec", "galnodO_x1000", "magdO", "magiO", "flag",
            "POSS_field", "starnumO", "sep_arcsec"]]
    df = df.sort_values("sep_arcsec").reset_index(drop=True)
    return df


def plates_overlapping_count(ra: float, dec: float,
                             radius_arcsec: float = 5.0) -> int | None:
    """Number of MAPS plates whose bbox overlaps this query (0 = coverage gap).

    Returns None if VASCO_MAPS_CACHE is unset.
    """
    cache = os.getenv("VASCO_MAPS_CACHE")
    if not cache:
        return None
    radius_deg = radius_arcsec / 3600.0
    return len(_plates_overlapping(cache, ra, dec, radius_deg))
