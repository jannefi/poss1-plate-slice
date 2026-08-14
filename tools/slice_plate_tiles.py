#!/usr/bin/env python3
"""Slice a local full-plate POSS-I scan into pipeline-ready tiles.

This is a drop-in replacement for `step1-download`. Instead of asking an archive
for the sky at a position and accepting whichever plate it decides to serve, it
takes a named plate off local disk and cuts the tiles out of that plate. Plate
identity stops being a lottery outcome and becomes an input.

That distinction is the whole point. The archive's cutout services are addressed
by coordinates, never by plate, so `plate(tile centre)` and `plate(exact object
position)` disagree for a large fraction of tiles -- measured at 36% on the
635-plate campaign -- and no client-side tiling scheme can see or control it.
Slicing locally removes the selection step entirely.

The archive's 60'x60' cap is a bandwidth policy, not a scientific constraint: its
cutout is itself a subarray extraction from this same plate. Verified directly --
a local slice matches the served cutout of the same sky at pixel correlation
1.000000, and SExtractor run on both finds 99.64% of sources in common within 1"
with a median MAG_AUTO difference of 0.0004.

The WCS is the whole difficulty -- read this before changing it
---------------------------------------------------------------
Full-plate DSS FITS carry the GSSS polynomial plate solution (AMDX1..20,
AMDY1..20, PPO1..6, XPIXELSZ, PLTSCALE) and NO standard CTYPE/CD/CRVAL.

astropy reads that solution correctly *in memory* -- confirmed against Gaia with
a null control: real positions land 1.0 px from a stellar peak, positions shifted
by 0.05 deg land 20 px away. `Cutout2D` also carries it through a slice.

But `WCS.to_header()` serialises it as PC + CDELT, and SExtractor's WCS reader
ignores PC matrices. It then falls back on CDELT -- which for a DSS solution is
in degrees per *millimetre*, roughly 40x the pixel scale. The resulting file is
astropy-valid and round-trips perfectly through astropy, so nothing looks wrong;
SExtractor silently produces sky coordinates that are wrong by degrees. In the
first gate run this showed up as two catalogues with near-identical source counts
sharing 0.05% of their positions.

So this module does two things that are NOT optional:

  1. Refits a clean TAN from points sampled across the slice, using the working
     in-memory WCS as truth (`fit_wcs_from_points`). Residual is ~0.10" median,
     0.33" max -- well inside the 3.0" dedup tolerance.
  2. Serialises as an explicit CD matrix, stripping PC/CDELT/CUNIT.

Do not "simplify" either step back to `cut.wcs.to_header()`. Validate any change
by running SExtractor on both a served tile and a local slice of the same sky and
comparing sky positions -- comparing headers, or round-tripping through astropy,
will not catch this class of failure.

The grid is laid out in pixel space, and must stay that way
-----------------------------------------------------------
Tile centres are stepped across the pixel array and converted to world
coordinates, never the reverse. Stepping in RA/Dec -- even with the RA step
divided by cos(dec), which keeps the sky angle right -- fails near the pole,
because a rectangle in RA/Dec is not a rectangle in the plate's tangent plane.

On XE002 (dec +84.75) the RA step is 0.9346/cos(87.56) = 21.9 deg per column, so
the top-row corners sat +/-65.8 deg away in RA and fell off the array entirely:
two tiles lost to "Arrays do not overlap", and only 90.1% of the plate covered.
Measured over 634 plates, the RA/Dec grid covered 90-93% of the ten plates above
dec 83 and ~96% of the sixteen between 77 and 79, against ~99% elsewhere.

Cutout2D is likewise given its position and size in pixels. Passing an angular
size makes it re-derive the pixel size from the local WCS, which reintroduces the
same distortion in the tile dimensions.

A third of plates need a one-pixel CRPIX correction
---------------------------------------------------
Pass --crpix-table (from tools/build_plate_crpix_table.py) or ~33% of plates come
out ~2.4" from Gaia. The plate solution's pixel origin is ambiguous by one pixel
and which way it falls is plate-dependent: 207 of 633 plates need CRPIX-1 in both
axes, and the other 426 are already at ~0.1" and are RUINED by the same shift
(XE185 measured 0.229" -> 2.211"). So it must be a per-plate lookup, never a
global constant.

The table is discrete -- 0 or -1, chosen by which hypothesis leaves less residual
against Gaia, no fitted parameters -- and derived only from Gaia and public plate
headers, so anyone can regenerate and audit it. Survey-wide it takes the mean
offset from 0.867" to 0.163".

The root cause is still open: the correction is measured per plate, not predicted
from the header. AMDX8/AMDY8 separate the two groups at only 86% accuracy, so
they are a lead, not a rule. Russell et al. 1990 (AJ 99, 2059) defines the GSSS
astrometric model and is the place to look.

Safety
------
Refuses to write into a tiles root containing 'tiles_archive' unless
--allow-production is passed. The production archive and its shared metadata
directory have been damaged by concurrent writes before; pilots get their own
tree.

How to validate
---------------
    python3 tools/slice_plate_tiles.py \
        --plate-fits <plate_dir>/dss1red_XE524.fits \
        --tiles-dir work/runs/plate_slice_pilot/tiles \
        --tiles-file-out work/runs/plate_slice_pilot/tiles_XE524.txt

Then run steps 2-5 over the emitted tiles file exactly as for downloaded tiles.
Every tile's dss1red_title.txt must show the REGION of the plate it was cut
from; if any disagrees, the grid or the plate library is wrong.
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from vasco.utils.tile_id import format_tile_id  # noqa: E402


# Any tile whose TAN refit lands above this is not written. Healthy tiles sit at
# 0.04-0.17"; the diverged ones were 29" and 144". Nothing legitimate is near 1".
MAX_REFIT_RESID_ARCSEC = 1.0


def plate_centre(h) -> tuple[float, float]:
    ra = (h["PLTRAH"] + h["PLTRAM"] / 60.0 + h["PLTRAS"] / 3600.0) * 15.0
    sgn = -1.0 if str(h.get("PLTDECSN", "+")).strip().startswith("-") else 1.0
    dec = sgn * (abs(h["PLTDECD"]) + h["PLTDECM"] / 60.0 + h["PLTDECS"] / 3600.0)
    return ra, dec


def clean_tan_header(cut, plate_header, n_fit: int = 25):
    """Refit a TAN from sampled points and serialise it as a CD matrix.

    Both halves matter -- see the module docstring. Returns (header, residual).
    """
    from astropy.io import fits
    from astropy.wcs.utils import fit_wcs_from_points

    ny, nx = cut.data.shape
    fx, fy = np.meshgrid(np.linspace(0, nx - 1, n_fit), np.linspace(0, ny - 1, n_fit))
    fx, fy = fx.ravel(), fy.ravel()
    sky = cut.wcs.pixel_to_world(fx, fy)

    # proj_point is REQUIRED, not cosmetic. Left at its default of "center",
    # fit_wcs_from_points derives the fiducial from lon.min()/lon.max(), which for
    # a field crossing RA 0 are both *at* the wrap rather than at the field's
    # edges -- so the fiducial lands near the meridian instead of on the tile. On
    # an easy field the fit absorbs that into CD/CRPIX and nothing is lost, which
    # is why this survived a full survey unnoticed. Near the pole, where the GSSS
    # solution is genuinely hard to represent as a plain TAN, the fit instead
    # DIVERGES: two tiles of XE011 came out at median residuals of 143.9" and
    # 29.2" against ~0.10" for the other 47. Passing the cutout's own centre fixes
    # both to ~0.10" and changes nothing on any other tile.
    ctr = cut.wcs.pixel_to_world((nx - 1) / 2.0, (ny - 1) / 2.0)
    fw = fit_wcs_from_points((fx, fy), sky, proj_point=ctr,
                             projection="TAN", sip_degree=None)
    resid = fw.pixel_to_world(fx, fy).separation(sky).arcsec

    hdr = fw.to_header(relax=True)
    cd = {}
    for i in (1, 2):
        d = hdr.get(f"CDELT{i}", 1.0)
        for j in (1, 2):
            cd[(i, j)] = d * hdr.get(f"PC{i}_{j}", 1.0 if i == j else 0.0)
    for i in (1, 2):
        for j in (1, 2):
            hdr[f"CD{i}_{j}"] = cd[(i, j)]
            hdr.remove(f"PC{i}_{j}", ignore_missing=True)
        hdr.remove(f"CDELT{i}", ignore_missing=True)
        hdr.remove(f"CUNIT{i}", ignore_missing=True)
    for k in ("REGION", "PLTLABEL", "PLATEID", "DATE-OBS", "SURVEY", "EQUINOX"):
        if k in plate_header:
            hdr[k] = plate_header[k]
    hdr["VASCOSRC"] = ("local-plate-slice", "cut from full plate, not archive cutout")
    return fits.Header(hdr), float(np.median(resid)), float(resid.max())


def write_tile(tile_dir: Path, name: str, data, hdr, plate_header, size_arcmin):
    from astropy.io import fits

    raw = tile_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    fpath = raw / name
    fits.PrimaryHDU(data.astype(np.int16), header=hdr).writeto(fpath, overwrite=True)

    sel_keys = ("SURVEY", "PLATEID", "PLATE-ID", "PLATE", "DATE-OBS", "RA", "DEC",
                "EQUINOX", "MJD-OBS", "NAXIS1", "NAXIS2", "CD1_1", "CD1_2", "CD2_1",
                "CD2_2", "CDELT1", "CDELT2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2")
    full = fits.getheader(fpath)
    sidecar = {
        "fits_file": name,
        "selected": {k: (str(full[k]) if k in full else None) for k in sel_keys},
        "header": {k: str(full[k]) for k in full},
    }
    (raw / f"{name}.header.json").write_text(json.dumps(sidecar, indent=2))

    (raw / "dss1red_title.txt").write_text(
        "\n".join([
            f"PLTLABEL: {plate_header.get('PLTLABEL','')}",
            f"PLATEID: {plate_header.get('PLATEID','')}",
            f"REGION: {plate_header.get('REGION','')}",
            f"DATE-OBS: {plate_header.get('DATE-OBS','')}",
            f"FITS: {name}",
            f"SOURCE: {name}.header.json",
            "SEP_DEG:",
        ]) + "\n")

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    (tile_dir / "tile_status.json").write_text(json.dumps(
        {"tile_id": tile_dir.name,
         "steps": {"step1": {"status": "ok", "ts": ts, "via": "slice_plate_tiles"}}},
        indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plate-fits", required=True)
    ap.add_argument("--tiles-dir", required=True, help="Tiles ROOT (tile folders created under it).")
    ap.add_argument("--tiles-file-out", default=None, help="Write emitted tile dirs, one per line.")
    ap.add_argument("--grid", type=int, default=7,
                    help="grid x grid tiles per plate [7]. Step is derived in PIXELS so "
                         "the grid spans the array exactly, so grid size sets the overlap: 7 gives "
                         "~6.5%% (and 49 tiles, matching the campaign's per-plate count, "
                         "which keeps before/after comparisons like-for-like); 8 gives "
                         "~20%%. Even 6.5%% is ~3.9 arcmin, far larger than any object, so "
                         "seam splitting is not a concern at either setting.")
    ap.add_argument("--size-arcmin", type=float, default=60.0)
    ap.add_argument("--survey-prefix", default="dss1-red")
    ap.add_argument("--allow-production", action="store_true",
                    help="Permit writing into a production tiles_archive tree. Do not use.")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--crpix-table", default=None,
                    help="CSV from tools/build_plate_crpix_table.py. Plates it flags get "
                         "CRPIX shifted by -1 in both axes. Without it, ~33%% of plates "
                         "come out ~2.4\" from Gaia; see the module docstring. A plate the "
                         "table does not cover is a fatal error, not a warning.")
    ap.add_argument("--allow-missing-crpix", action="store_true",
                    help="Slice plates the CRPIX table does not cover, uncorrected. Only for "
                         "deliberately working against a partial table -- the resulting tiles "
                         "may sit ~2.3\" off with nothing downstream to flag it.")
    args = ap.parse_args()

    tiles_dir = Path(args.tiles_dir)
    if "tiles_archive" in str(tiles_dir.resolve()) and not args.allow_production:
        raise SystemExit(f"[FATAL] refusing to write into production archive: {tiles_dir}\n"
                         f"        pilots must use their own tiles root.")

    from astropy.io import fits
    from astropy.nddata import Cutout2D
    from astropy.wcs import WCS

    # Read the plate fully into RAM (392 MB) rather than memmapping it. With
    # memmap, each of the 49 Cutout2D calls re-reads from the HDD and the pass
    # costs ~83 s; resident it is a few seconds. Memory is bounded and small.
    hdu = fits.open(args.plate_fits, memmap=False)[0]
    _ = hdu.data.shape
    ph = hdu.header
    pw = WCS(ph)
    region = str(ph.get("REGION", "")).strip()
    scale = ph["XPIXELSZ"] / 1000.0 * ph["PLTSCALE"]          # micron -> mm -> arcsec/px
    extent_deg = ph["NAXIS1"] * scale / 3600.0
    # Centre the grid on the IMAGE centre, NOT on the PLTRAH/PLTDEC keywords.
    # Those name the plate's nominal sky centre, which is not where the scanned
    # array is centred: measured across 634 plates the two differ by a median
    # 0.071 deg, by >0.25 deg on 12 plates, and by ~4.4 deg on five (XE509,
    # XE543, XE574, XE284, XE541). Grid-from-keywords therefore walks off the
    # array -- XE284 lost 33 of its 49 tiles to "Arrays do not overlap" -- and
    # even where every tile succeeds it shifts coverage off one edge by the
    # offset. The image centre is by construction the centre of the data.
    ra_c, dec_c = [float(v) for v in pw.pixel_to_world_values(
        ph["NAXIS1"] / 2.0, ph["NAXIS2"] / 2.0)]
    kra, kdec = plate_centre(ph)
    _off = np.hypot((ra_c - kra) * np.cos(np.radians(dec_c)), dec_c - kdec)

    # Per-plate one-pixel origin correction. Roughly a third of plates carry a
    # ~2.4" offset against Gaia that CRPIX-1 in both axes removes, while the
    # rest already sit at ~0.1" and would be ruined by the same shift -- so it
    # is looked up per plate and never applied globally.
    # A plate that the table does not cover is FATAL, not a warning. The older
    # behaviour printed one line into this plate's own slicer log and carried
    # on. That line reaches nothing anyone watches -- not progress.csv, not the
    # console summary, not run.log -- while the tiles come out looking entirely
    # normal and sitting ~2.3" off, inside the 5" veto radius where no pass/fail
    # test can see it. It fired six times in the 2026-08-12 run, when the
    # manifest grew to 642 plates and this table still held 634, and was
    # harmless only because those six happened to need corrections <0.21".
    crpix_dx = crpix_dy = 0.0
    if args.crpix_table:
        import pandas as pd
        _t = pd.read_csv(args.crpix_table)
        _r = _t[_t.plate.astype(str) == region]
        _why = None
        if not len(_r):
            _why = f"{region} is not in {args.crpix_table}"
        elif str(_r.status.iloc[0]) != "ok":
            _why = (f"{region} has status={str(_r.status.iloc[0])!r} in "
                    f"{args.crpix_table}, so no usable correction was measured")
        if _why and not args.allow_missing_crpix:
            raise SystemExit(
                f"[FATAL] {_why}.\n"
                f"        Slicing it would inherit the GSSS solution, which is ~2.3\" off on\n"
                f"        about a third of plates -- inside the 5\" veto radius, so nothing\n"
                f"        downstream would flag it. Extend the table to cover this plate:\n"
                f"          python3 tools/build_plate_crpix_table.py \\\n"
                f"              --plate-dir <plate_dir> --archive-tiles <tiles_dir> \\\n"
                f"              --out {args.crpix_table}\n"
                f"        (it needs one archive cutout carrying both WCS solutions for this\n"
                f"        plate), or pass --allow-missing-crpix to slice it uncorrected on\n"
                f"        purpose."
            )
        if _why:
            print(f"[CRPIX] {_why}")
            print(f"[CRPIX] {region} sliced UNCORRECTED by explicit --allow-missing-crpix")
        else:
            crpix_dx = float(_r.delta_x_px.iloc[0])
            crpix_dy = float(_r.delta_y_px.iloc[0])
            if float(_r.scatter_px.iloc[0]) > 0.2:
                print(f"[CRPIX] {region} scatter {float(_r.scatter_px.iloc[0]):.3f} px "
                      f"-- correction applied but less certain than usual")
    print(f"[CRPIX] correction dx={crpix_dx:+.3f} dy={crpix_dy:+.3f} px")

    # The grid is laid out in PIXEL space, not in RA/Dec -- see the module
    # docstring. Centres run from half a tile in to half a tile from the far
    # edge, so the grid spans the array exactly and `grid` alone sets the
    # overlap.
    n = args.grid
    tile_px = args.size_arcmin * 60.0 / scale
    nx, ny = int(ph["NAXIS1"]), int(ph["NAXIS2"])
    tw, th = int(round(tile_px)), int(round(tile_px))

    def centres(span: int) -> np.ndarray:
        if n == 1:
            return np.array([span / 2.0])
        lo, hi = tile_px / 2.0, span - tile_px / 2.0
        if hi <= lo:                                   # tile wider than the plate
            return np.full(n, span / 2.0)
        return np.linspace(lo, hi, n)

    cx, cy = centres(nx), centres(ny)
    step_px = float(cx[1] - cx[0]) if n > 1 else 0.0
    overlap = 100.0 * (1.0 - step_px / tile_px) if tile_px else 0.0
    print(f"[PLATE] {region}  {nx}x{ny} px  {scale:.3f}\"/px  "
          f"extent {extent_deg:.3f} deg  image centre {ra_c:.5f} {dec_c:.5f}  "
          f"(PLTRA/DEC keywords differ by {_off:.3f} deg)")
    print(f"[GRID]  {n}x{n} tiles of {args.size_arcmin:.0f}' ({tw} px)  "
          f"step {step_px:.1f} px = {step_px * scale / 3600.0:.4f} deg  "
          f"overlap {overlap:.1f}%")

    emitted, resid_all, bad_wcs = [], [], []
    for iy, py in enumerate(cy):
        for ix, px in enumerate(cx):
            ra, dec = [float(v) for v in pw.pixel_to_world_values(px, py)]
            ra %= 360.0
            tid = format_tile_id(ra, dec)
            tdir = tiles_dir / tid
            name = f"{args.survey_prefix}_{ra:.3f}_{dec:.3f}_{args.size_arcmin:.0f}arcmin.fits"
            if (tdir / "raw" / name).exists() and not args.overwrite:
                emitted.append(str(tdir))
                continue
            try:
                # Position and size in pixels: exact, and free of the angular
                # -> pixel conversion that distorts near the pole.
                cut = Cutout2D(hdu.data, position=(px, py), size=(th, tw),
                               wcs=pw, copy=True, mode="trim")
            except Exception as exc:                       # off-plate corner
                print(f"  [SKIP] {tid}: {exc}")
                continue
            hdr, rmed, rmax = clean_tan_header(cut, ph)
            # A tile whose WCS refit did not converge is worse than a missing
            # tile: its detections are real but land at wrong coordinates, no
            # veto can match them, and they survive to the catalogue in their
            # hundreds. Two such tiles contributed 2,096 rows to an earlier
            # release before anyone looked at a PER-TILE residual. Refuse to
            # write them.
            if rmed > MAX_REFIT_RESID_ARCSEC:
                print(f"  [BAD-WCS] {tid}: TAN refit residual {rmed:.3f}\" "
                      f"(max {rmax:.3f}\") exceeds {MAX_REFIT_RESID_ARCSEC}\" "
                      f"-- tile NOT written")
                bad_wcs.append((tid, rmed, rmax))
                continue
            # CRPIX_corrected = CRPIX - delta. A WCS whose CRPIX is reduced by d
            # evaluates at pixel p exactly as the original would at p+d, and delta
            # is measured as (where GSSS puts a position) - (where it really is).
            # The opposite sign roughly doubles the error instead of removing it.
            if crpix_dx or crpix_dy:
                hdr["CRPIX1"] = hdr["CRPIX1"] - crpix_dx
                hdr["CRPIX2"] = hdr["CRPIX2"] - crpix_dy
                hdr["VASCOCPX"] = (f"{crpix_dx:+.4f},{crpix_dy:+.4f}",
                                   "per-plate CRPIX correction, header-derived")
            resid_all.append(rmed)
            write_tile(tdir, name, cut.data, hdr, ph, args.size_arcmin)
            emitted.append(str(tdir))
        print(f"  row {iy+1}/{n} done ({len(emitted)} tiles)", flush=True)

    if resid_all:
        # Report the worst tile alongside the median. The median alone hides a
        # diverged fit completely -- 47 good tiles average a 143.9" failure down
        # to 0.1041", which is exactly how this defect reached a release.
        print(f"[WCS]   TAN refit residual: median of medians "
              f"{np.median(resid_all):.4f}\"  worst tile {max(resid_all):.4f}\"")
    if bad_wcs:
        print(f"[WCS][FAIL] {len(bad_wcs)} tile(s) exceeded "
              f"{MAX_REFIT_RESID_ARCSEC}\" and were not written:")
        for tid, rmed, rmax in bad_wcs:
            print(f"           {tid}  median {rmed:.3f}\"  max {rmax:.3f}\"")
    print(f"[OUT]   {len(emitted)} tiles under {tiles_dir}")
    if args.tiles_file_out:
        Path(args.tiles_file_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.tiles_file_out).write_text("\n".join(emitted) + "\n")
        print(f"[OUT]   tiles file {args.tiles_file_out}")


if __name__ == "__main__":
    main()
