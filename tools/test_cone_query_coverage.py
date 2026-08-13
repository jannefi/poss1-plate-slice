#!/usr/bin/env python3
"""Regression test: the local-cache cone query must load the WHOLE cone.

Guards the 2026-08 partial-cone bug. `HEALPix.cone_search_lonlat` returns only
pixels whose CENTRES lie inside the radius, so at nside=32 (pixels ~1.8 deg
across, veto cone ~0.76 deg) it silently drops pixels that overlap the cone
with their centre outside. Nothing raises; the veto catalogue is just short,
and everything in the missing sky survives the veto. On the worst tiles of the
642-plate run 45-48% of the Gaia cone was never read, which put ~173k
un-vetoed Gaia stars into S0.

Two independent checks:

  A. GEOMETRY (no data needed, always runs) -- `_cone_pixels` must return a
     superset of the pixels a dense sampling of the cone actually lands in.
     Includes the tile centres that were observed to fail in production and a
     declination sweep across the HEALPix polar-cap boundary at +/-41.81 deg.

  B. DATA (skipped unless VASCO_GAIA_CACHE is set) -- the row count the query
     returns must equal a brute-force scan over a generously padded pixel set.
     This is the end-to-end statement: no catalogue row inside the radius is
     missed.

Run:
    python3 tools/test_cone_query_coverage.py
    VASCO_GAIA_CACHE=/path/to/gaia python3 tools/test_cone_query_coverage.py

Exit code 0 = pass, 1 = fail. Against the pre-fix code, check A fails on the
production tile centres below.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import astropy.units as u
from astropy_healpix import HEALPix

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from vasco.local_cache_query import _cone_pixels, _HP_NSIDE  # noqa: E402

# step4 uses radius = size_arcmin * sqrt(2)/2 + 3.0; 60' tiles -> 45.4'
RADIUS_ARCMIN = 45.4

# Tile centres measured to lose 28-51% of their Gaia cone in run
# full642-20260812. All sit just above the polar-cap boundary.
PRODUCTION_FAILURES = [
    (1.662, 59.887),    # XE078, 2913 survivors, 47.0% of cone missing
    (7.775, 62.552),    # XE078, 2764 survivors, 48.2%
    (357.938, 59.840),  # XE078, 1204 survivors, 45.1%
    (1.657, 58.951),    # XE078,  143 survivors, 11.1%
    (85.014, 43.904),   # XE203, 1592 survivors, 47.3%
    (82.457, 43.885),   # XE203,  978 survivors, 36.4%
    (345.332, 77.776),  # XE027, 2073 survivors, 51.0%
    (355.553, 79.502),  # XE027,  836 survivors, 33.0%
]

# Tile centres that were already complete -- these must STAY complete, i.e.
# the fix must not be a blanket widening that hides a real regression.
PRODUCTION_HEALTHY = [
    (3.469, 58.936),    # XE078, 7 survivors, 0.0% missing
    (5.386, 59.828),    # XE078, 0 survivors
    (81.127, 44.775),   # XE203, 8 survivors
    (0.006, 29.192),    # XE347, 2 survivors
]


def truth_pixels(hp, ra, dec, radius_arcmin, n=120_000, seed=11):
    """Pixels a dense sampling of the cone actually falls into."""
    rng = np.random.RandomState(seed)
    rad = np.radians(radius_arcmin / 60.0)
    cth = 1 - rng.rand(n) * (1 - np.cos(rad))
    th = np.arccos(cth)
    ph = rng.rand(n) * 2 * np.pi
    # dense rim: boundary pixels are the ones at risk
    m = 20_000
    th = np.concatenate([th, np.full(m, rad)])
    ph = np.concatenate([ph, np.linspace(0, 2 * np.pi, m)])

    x = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    v0 = np.array([np.cos(np.radians(dec)) * np.cos(np.radians(ra)),
                   np.cos(np.radians(dec)) * np.sin(np.radians(ra)),
                   np.sin(np.radians(dec))])
    z = np.array([0.0, 0.0, 1.0])
    ax = np.cross(z, v0)
    s = np.linalg.norm(ax)
    if s < 1e-12:
        rot = np.eye(3)
    else:
        ax = ax / s
        a = np.arccos(np.clip(np.dot(z, v0), -1, 1))
        k = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        rot = np.eye(3) + np.sin(a) * k + (1 - np.cos(a)) * k @ k
    p = rot @ x
    lat = np.degrees(np.arcsin(np.clip(p[2], -1, 1)))
    lon = np.degrees(np.arctan2(p[1], p[0])) % 360.0
    return set(int(q) for q in hp.lonlat_to_healpix(lon * u.deg, lat * u.deg).tolist())


def check_geometry() -> int:
    hp = HEALPix(nside=_HP_NSIDE, order="nested")
    cases = [("production-failure", ra, dec) for ra, dec in PRODUCTION_FAILURES]
    cases += [("production-healthy", ra, dec) for ra, dec in PRODUCTION_HEALTHY]
    # sweep across the polar-cap boundary and up to the pole
    for dec in [0, 10, 20, 30, 38, 41.0, 41.81, 42.5, 45, 50, 55, 60, 65, 70, 75, 80, 85, 88]:
        for ra in np.arange(0.0, 360.0, 11.0):
            cases.append(("sweep", float(ra), float(dec)))

    bad = []
    npix = []
    for label, ra, dec in cases:
        got = set(_cone_pixels(hp, ra, dec, RADIUS_ARCMIN))
        npix.append(len(got))
        missing = truth_pixels(hp, ra, dec, RADIUS_ARCMIN) - got
        if missing:
            bad.append((label, ra, dec, sorted(missing)))

    print(f"[A] geometry: {len(cases)} cone centres, "
          f"median {np.median(npix):.0f} pixels/cone, max {max(npix)}")
    if bad:
        print(f"[A] FAIL -- {len(bad)} centres miss an overlapping pixel:")
        for label, ra, dec, miss in bad[:12]:
            print(f"      {label:20s} ra={ra:7.3f} dec={dec:7.3f}  missing {miss[:6]}")
        return 1
    print("[A] PASS -- no overlapping pixel is ever dropped")
    return 0


def check_data() -> int:
    cache = os.getenv("VASCO_GAIA_CACHE")
    if not cache:
        print("[B] SKIP -- set VASCO_GAIA_CACHE to run the data check")
        return 0

    import glob
    import pandas as pd
    from vasco.local_cache_query import _cone_query

    hp = HEALPix(nside=_HP_NSIDE, order="nested")

    def brute(ra, dec):
        """Every mirror row within the radius, from a generously padded set."""
        cand = hp.cone_search_lonlat(
            ra * u.deg, dec * u.deg,
            radius=(RADIUS_ARCMIN + 200.0) * u.arcmin)
        frames = []
        for p in np.unique(np.asarray(cand, dtype=np.int64)):
            for f in glob.glob(f"{cache}/parquet/healpix_5={int(p)}/*.parquet"):
                frames.append(pd.read_parquet(f, columns=["ra", "dec"]))
        if not frames:
            return 0
        df = pd.concat(frames)
        a = np.deg2rad(df.ra.values)
        d = np.deg2rad(df.dec.values)
        cr, cd = np.deg2rad(ra), np.deg2rad(dec)
        cs = np.clip(np.sin(d) * np.sin(cd) + np.cos(d) * np.cos(cd) * np.cos(a - cr), -1, 1)
        return int((np.rad2deg(np.arccos(cs)) * 60.0 <= RADIUS_ARCMIN).sum())

    fails = 0
    for ra, dec in PRODUCTION_FAILURES[:4] + PRODUCTION_HEALTHY[:2]:
        got = len(_cone_query(cache, ra, dec, RADIUS_ARCMIN,
                              columns=["ra", "dec", "phot_g_mean_mag", "pmra", "pmdec"]))
        want = brute(ra, dec)
        ok = got == want
        fails += (not ok)
        print(f"[B] ra={ra:7.3f} dec={dec:7.3f}  query={got:7d}  brute={want:7d}  "
              f"{'OK' if ok else 'FAIL (short by %d)' % (want - got)}")
    if fails:
        print(f"[B] FAIL -- {fails} cone(s) returned fewer rows than the sky contains")
        return 1
    print("[B] PASS -- every catalogue row inside the radius is returned")
    return 0


if __name__ == "__main__":
    rc = check_geometry() | check_data()
    print("\nRESULT:", "PASS" if rc == 0 else "FAIL")
    sys.exit(rc)
