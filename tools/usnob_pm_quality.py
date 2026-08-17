#!/usr/bin/env python3
"""Are USNO-B1.0's high proper motions real? Measured against Gaia DR3.

Why this exists
---------------
A veto stage that matches a catalogue against 1950s plate detections has to
decide what to do about proper motion. Matching at the catalogue's own epoch
misses genuine fast stars, which sit 10-25" from their J2000 position on a
POSS-I plate. Propagating the catalogue back to the plate epoch fixes that --
but only if the proper motions are trustworthy.

This measures whether USNO-B1.0's are, using nothing but the two catalogues.
It is independent of any pipeline: give it two local mirrors and a list of
plate centres and it reproduces the numbers in
docs/FUNNEL_ATTRIBUTION.md.

Method
------
Random cones inside the plate footprint. Every USNO-B entry gets a SYMMETRIC
pair of tests against Gaia DR3, because the obvious one-sided test is biased:
a real 500 mas/yr star has moved ~8" between J2000 and Gaia's J2016 epoch and
would fail a "is Gaia at the J2000 position" check purely by moving, leaving
only stationary stars in the associated set and guaranteeing the conclusion.

  FABRICATED     a Gaia source within --assoc-arcsec of the USNO-B *catalogued*
                 J2000 position, itself moving slower than --lo-pm-masyr.
                 The entry is a standing star and its catalogued PM is not real.
  GENUINE        a Gaia source within --assoc-arcsec of the USNO-B position
                 PROPAGATED BY ITS OWN PM to Gaia's epoch, with a Gaia PM within
                 --pm-tol-factor of the claimed value. The star really moves,
                 and roughly as USNO-B says.
  INDETERMINATE  neither.

Read the high-PM stratum against the low-PM control, which establishes the rate
at which the association step finds a counterpart at all.

The companion measurement, --density-check, simply counts how many fast stars
Gaia finds in the same sky as USNO-B claims. That comparison needs no
classification at all and is the harder number to argue with.

Inputs
------
Two local parquet mirrors, hive-partitioned on a `healpix_5` column built with
nside=32 in **NESTED** order (the same convention as
tools/catalog_xmatch_local.py -- reading them as RING silently queries the
wrong sky), each carrying ra/dec plus proper motions:

  USNO-B1.0   ra, dec, pmRA, pmDE, B1mag, R1mag, B2mag, R2mag, Imag
  Gaia DR3    ra, dec, pmra, pmdec

and a CSV of plate centres with plate_id, ra_deg, dec_deg (this repository's
data/plate_manifest.csv).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as ds
from astropy import units as u
from astropy_healpix import HEALPix
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[1]

# Widen the pixel search by more than any nside=32 pixel's circumradius, so no
# overlapping pixel can be missed before the explicit prune below.
PIXEL_MAX_RADIUS_ARCMIN = 120.0

USNOB_COLS = ["ra", "dec", "pmRA", "pmDE", "B1mag", "R1mag", "B2mag", "R2mag", "Imag"]
GAIA_COLS = ["ra", "dec", "pmra", "pmdec"]


def xyz(ra, dec):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def chord(arcsec):
    return 2.0 * math.sin(math.radians(arcsec / 3600.0) / 2.0)


def cone_pixels(hp, ra, dec, radius_arcmin):
    """Every pixel that OVERLAPS the cone, not merely those centred inside it.

    `HEALPix.cone_search_lonlat` returns pixels whose *centres* fall within the
    radius. At nside=32 a pixel spans ~1.8 deg, so a pixel can overlap a
    sub-degree cone while its centre lies well outside -- and it is then
    silently dropped, leaving part of the cone unread with no error and no
    empty file to notice. Searching wide and pruning back is the fix; the
    widened search costs nothing because the prune runs before any I/O.
    """
    cand = hp.cone_search_lonlat(ra * u.deg, dec * u.deg,
                                 radius=(radius_arcmin + PIXEL_MAX_RADIUS_ARCMIN) * u.arcmin)
    cand = np.unique(np.asarray(cand, dtype=np.int64))
    if cand.size == 0:
        return []
    r_deg = radius_arcmin / 60.0
    cra, cdec = math.radians(ra), math.radians(dec)

    def sep_deg(lon, lat):
        a, d = np.radians(np.asarray(lon, float)), np.radians(np.asarray(lat, float))
        cs = np.clip(np.sin(d) * math.sin(cdec) + np.cos(d) * math.cos(cdec) * np.cos(a - cra),
                     -1.0, 1.0)
        return np.degrees(np.arccos(cs))

    # These three are exhaustive: if the cone is not wholly inside the pixel
    # and the pixel is not wholly inside the cone, any intersection must make
    # the pixel boundary cross the cone.
    #
    # 1) pixel centre inside the cone
    clon, clat = hp.healpix_to_lonlat(cand)
    keep = sep_deg(clon.to_value(u.deg), clat.to_value(u.deg)) <= r_deg
    # 2) cone centre inside the pixel (cone smaller than one pixel)
    keep |= cand == int(hp.lonlat_to_healpix(ra * u.deg, dec * u.deg))
    # 3) pixel boundary enters the cone. step=32 samples every ~0.055 deg; the
    #    0.05 deg slack covers grazing contact where the chord inside the cone
    #    is shorter than the sample spacing.
    blon, blat = hp.boundaries_lonlat(cand, 32)
    bsep = sep_deg(blon.to_value(u.deg), blat.to_value(u.deg))
    keep |= (bsep <= r_deg + 0.05).any(axis=1)
    return [int(p) for p in cand[keep]]


def cone(dataset, hp, ra, dec, radius_arcmin, columns):
    """Rows of a healpix-partitioned mirror inside a cone, as numpy columns."""
    pixels = cone_pixels(hp, ra, dec, radius_arcmin)
    if not pixels:
        return None
    tbl = dataset.to_table(columns=columns, filter=pc.field("healpix_5").isin(pixels))
    if tbl.num_rows == 0:
        return None
    out = {c: tbl.column(c).to_numpy(zero_copy_only=False).astype(float) for c in columns}
    r1, d1 = math.radians(ra), math.radians(dec)
    r2, d2 = np.radians(out["ra"]), np.radians(out["dec"])
    cs = np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2), -1, 1)
    keep = np.degrees(np.arccos(cs)) * 60.0 <= radius_arcmin
    if not keep.any():
        return None
    return {c: v[keep] for c, v in out.items()}


def sample_fields(centres, n_fields, offset_deg, seed):
    """Random points on plate sky: uniform inside a disc around each centre."""
    rng = np.random.default_rng(seed)
    pick = rng.choice(len(centres), size=min(n_fields, len(centres)), replace=False)
    for i in pick:
        pid, pra, pdec = centres[int(i)]
        rad = offset_deg * math.sqrt(float(rng.random()))
        th = 2 * math.pi * float(rng.random())
        dec = pdec + rad * math.sin(th)
        ra = (pra + (rad * math.cos(th)) / max(math.cos(math.radians(dec)), 1e-6)) % 360.0
        yield pid, ra, dec


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--usnob-mirror", required=True,
                    help="parquet dir, hive-partitioned on healpix_5 (nside=32 NESTED)")
    ap.add_argument("--gaia-mirror", required=True)
    ap.add_argument("--plate-manifest", default=str(REPO / "data" / "plate_manifest.csv"))
    ap.add_argument("--plate-subset", default=None,
                    help="optional CSV with a plate_id column, to restrict the draw")
    ap.add_argument("--n-fields", type=int, default=60)
    ap.add_argument("--cone-arcmin", type=float, default=12.0)
    ap.add_argument("--offset-deg", type=float, default=2.0)
    ap.add_argument("--hi-pm-masyr", type=float, default=150.0)
    ap.add_argument("--lo-pm-masyr", type=float, default=20.0)
    ap.add_argument("--assoc-arcsec", type=float, default=2.0)
    ap.add_argument("--pm-tol-factor", type=float, default=2.0)
    ap.add_argument("--gaia-epoch", type=float, default=2016.0)
    ap.add_argument("--nside", type=int, default=32)
    ap.add_argument("--seed", type=int, default=20260817)
    ap.add_argument("--density-check", action="store_true",
                    help="also count Gaia's own fast stars in the same fields")
    ap.add_argument("--out-dir", default="usnob_pm_quality")
    args = ap.parse_args()

    keep = None
    if args.plate_subset:
        keep = {r["plate_id"] for r in csv.DictReader(open(args.plate_subset))}
    centres = [
        (r["plate_id"], float(r["ra_deg"]), float(r["dec_deg"]))
        for r in csv.DictReader(open(args.plate_manifest))
        if keep is None or r["plate_id"] in keep
    ]
    print(f"[CAT] {len(centres)} plate centres available")

    hp = HEALPix(nside=args.nside, order="nested")   # MUST match the mirrors
    dsu = ds.dataset(args.usnob_mirror, format="parquet", partitioning="hive")
    dsg = ds.dataset(args.gaia_mirror, format="parquet", partitioning="hive")

    pm_all, r1_all, gpm0_all, gpm1_all = [], [], [], []
    gaia_total = gaia_fast = 0
    n_done = 0

    for k, (pid, ra, dec) in enumerate(sample_fields(centres, args.n_fields,
                                                     args.offset_deg, args.seed), 1):
        uu = cone(dsu, hp, ra, dec, args.cone_arcmin, USNOB_COLS)
        if uu is None:
            continue
        # pad the Gaia cone so entries near the rim can still find a counterpart
        gg = cone(dsg, hp, ra, dec, args.cone_arcmin + 1.0, GAIA_COLS)
        n_done += 1

        pmra, pmde = np.nan_to_num(uu["pmRA"]), np.nan_to_num(uu["pmDE"])
        pm = np.hypot(pmra, pmde)
        gpm0 = np.full(len(pm), np.nan)
        gpm1 = np.full(len(pm), np.nan)

        if gg is not None:
            tree = cKDTree(xyz(gg["ra"], gg["dec"]))
            gpm = np.hypot(np.nan_to_num(gg["pmra"]), np.nan_to_num(gg["pmdec"]))

            def probe(qra, qdec):
                d, j = tree.query(xyz(qra, qdec), k=1,
                                  distance_upper_bound=chord(args.assoc_arcsec))
                out = np.full(len(pm), np.nan)
                ok = np.isfinite(d)
                if ok.any():
                    out[ok] = gpm[j[ok]]
                return out

            gpm0 = probe(uu["ra"], uu["dec"])
            dt = args.gaia_epoch - 2000.0
            gpm1 = probe(uu["ra"] + (pmra * dt / 3.6e6) / np.cos(np.radians(uu["dec"])),
                         uu["dec"] + pmde * dt / 3.6e6)

            if args.density_check:
                # count inside the UNPADDED cone only -- gg is padded so that
                # rim entries can find a counterpart, but the padded area is
                # not what `area` below measures
                inside = cone(dsg, hp, ra, dec, args.cone_arcmin, GAIA_COLS)
                if inside is not None:
                    ipm = np.hypot(np.nan_to_num(inside["pmra"]), np.nan_to_num(inside["pmdec"]))
                    gaia_total += len(ipm)
                    gaia_fast += int((ipm >= args.hi_pm_masyr).sum())

        pm_all.append(pm)
        r1_all.append(np.isfinite(uu["R1mag"]))
        gpm0_all.append(gpm0)
        gpm1_all.append(gpm1)
        if k % 10 == 0:
            print(f"  [{k}] cumulative USNO-B entries: {sum(len(x) for x in pm_all)}",
                  flush=True)

    pm = np.concatenate(pm_all)
    r1 = np.concatenate(r1_all)
    gpm0, gpm1 = np.concatenate(gpm0_all), np.concatenate(gpm1_all)

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(gpm1 > 0, pm / gpm1, np.inf)
    genuine = (np.isfinite(gpm1) & (ratio <= args.pm_tol_factor)
               & (ratio >= 1.0 / args.pm_tol_factor))
    fabricated = np.isfinite(gpm0) & (gpm0 < args.lo_pm_masyr) & ~genuine

    def stratum(mask, name):
        n = int(mask.sum())
        if not n:
            return {"name": name, "n": 0}
        fa, ge = int((mask & fabricated).sum()), int((mask & genuine).sum())
        return {"name": name, "n": n,
                "fabricated": fa, "fabricated_pct": round(100 * fa / n, 2),
                "genuine": ge, "genuine_pct": round(100 * ge / n, 2),
                "indeterminate": n - fa - ge,
                "indeterminate_pct": round(100 * (n - fa - ge) / n, 2)}

    hi, lo = pm >= args.hi_pm_masyr, pm < args.lo_pm_masyr
    area = n_done * math.pi * (args.cone_arcmin / 60.0) ** 2
    summary = {
        "n_fields": n_done, "cone_arcmin": args.cone_arcmin,
        "area_deg2": round(area, 3), "n_usnob_entries": int(len(pm)), "seed": args.seed,
        "thresholds": {"hi_pm_masyr": args.hi_pm_masyr, "lo_pm_masyr": args.lo_pm_masyr,
                       "assoc_arcsec": args.assoc_arcsec, "pm_tol_factor": args.pm_tol_factor},
        "strata": [stratum(hi, f"USNO-B PM >= {args.hi_pm_masyr:.0f} mas/yr"),
                   stratum(hi & r1, "  ... and R1-backed (POSS-I red)"),
                   stratum(lo, f"CONTROL: USNO-B PM < {args.lo_pm_masyr:.0f} mas/yr")],
    }
    if args.density_check:
        n_hi = int(hi.sum())
        summary["density_check"] = {
            "gaia_sources": gaia_total,
            "gaia_pm_ge_threshold": gaia_fast,
            "gaia_per_deg2": round(gaia_fast / area, 2) if area else None,
            "usnob_pm_ge_threshold": n_hi,
            "usnob_per_deg2": round(n_hi / area, 2) if area else None,
            "usnob_overclaim_factor": round(n_hi / gaia_fast, 1) if gaia_fast else None,
        }

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "usnob_pm_quality.json").write_text(json.dumps(summary, indent=2))

    print(f"\n=== {len(pm)} USNO-B entries over {n_done} fields ({area:.2f} deg^2) ===")
    print(f"    {'stratum':46} {'n':>7} {'fabricated':>15} {'genuine':>15} {'indet':>15}")
    for s in summary["strata"]:
        if s["n"]:
            print(f"    {s['name']:46} {s['n']:>7} "
                  f"{s['fabricated']:>7} ({s['fabricated_pct']:5.1f}%) "
                  f"{s['genuine']:>7} ({s['genuine_pct']:5.1f}%) "
                  f"{s['indeterminate']:>7} ({s['indeterminate_pct']:5.1f}%)")
    if args.density_check:
        d = summary["density_check"]
        print(f"\n  density check over {area:.2f} deg^2, PM >= {args.hi_pm_masyr:.0f} mas/yr:")
        print(f"    Gaia DR3 : {d['gaia_pm_ge_threshold']:6d}  ({d['gaia_per_deg2']}/deg^2)")
        print(f"    USNO-B   : {d['usnob_pm_ge_threshold']:6d}  ({d['usnob_per_deg2']}/deg^2)")
        print(f"    USNO-B overclaims by {d['usnob_overclaim_factor']}x")
    print(f"\nwrote {out}/usnob_pm_quality.json")


if __name__ == "__main__":
    main()
