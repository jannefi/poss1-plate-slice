#!/usr/bin/env python3
"""Validate a cutout-based fast path for SPREAD_MODEL against real ground truth.

spread_model_postscore.py computes SPREAD_MODEL by re-running SExtractor's
pass2 on the FULL field (~88s/tile), then crossmatching survivors against it.
This script checks whether running pass2 on a small crop centered on each
survivor's known (X_IMAGE, Y_IMAGE) -- reusing the tile's already-built
pass1.psf -- gives close-enough SPREAD_MODEL values and, more importantly,
the same SPREAD_MODEL>-0.002 gate classification (context/02_DECISIONS.md).

Ground truth for each survivor comes from data already on disk from the
374-tile postscore run: catalogs/spread_model_postscore.csv (matched
SPREAD_MODEL from the full-field two-pass run) joined to
catalogs/sextractor_pass2.filtered.csv (X_IMAGE/Y_IMAGE, keyed by NUMBER).
No new full-field processing is performed -- only small crops.

Usage:
  python tools/validate_spread_model_cutout.py --tile-ids-file /path/to/ids.txt \
      [--tiles-root ./data/tiles_archive] [--halfwidths 150,300]
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from vasco.pipeline_split import run_pass2  # noqa: E402

SPREAD_MODEL_MIN = -0.002  # locked, context/02_DECISIONS.md
CONFIG_ROOT = str(REPO_ROOT / "configs")

_psfex_cache: dict[Path, object] = {}


def build_synthetic_psf(psf_path: Path, x: float, y: float, out_path: Path) -> None:
    """Evaluate the tile's real (spatially-varying) PSFEx model at (x, y) via
    GalSim's DES_PSFEx reader, and write it out as a minimal degree-0
    (POLNAXIS=0, single PSF_MASK component) PSFEx-format .psf file.

    This sidesteps the crop-vs-full-field coordinate-lookup bug: a degree-0
    PSF has no position dependence, so it doesn't matter that SExtractor
    computes X_IMAGE relative to the crop's own (wrong) local grid when it
    looks the PSF up -- there's only one PSF to find. The absolute flux
    normalization of the written mask does not matter: verified empirically
    (this session) that SPREAD_MODEL is invariant to an overall scale factor
    on the PSF template to ~1e-7, so GalSim's native flux=1 rendering is used
    as-is, no calibration against PSFEx's own internal normalization needed.
    """
    import galsim
    import galsim.des

    if psf_path not in _psfex_cache:
        _psfex_cache[psf_path] = galsim.des.DES_PSFEx(str(psf_path))
    des_psfex = _psfex_cache[psf_path]

    with fits.open(psf_path) as f:
        hdr = f[1].header
        nx, ny = hdr["PSFAXIS1"], hdr["PSFAXIS2"]
        samp = hdr["PSF_SAMP"]
        fwhm = hdr["PSF_FWHM"]
        loaded = hdr.get("LOADED", 0)
        accepted = hdr.get("ACCEPTED", 0)

    psf = des_psfex.getPSF(galsim.PositionD(x, y))
    img = psf.drawImage(nx=nx, ny=ny, scale=samp)
    mask = img.array.astype(">f4")

    col = fits.Column(name="PSF_MASK", format=f"{nx * ny}E",
                       array=mask.reshape(1, -1), dim=f"({nx}, {ny}, 1)")
    hdu = fits.BinTableHDU.from_columns([col], name="PSF_DATA")
    hdu.header["LOADED"] = loaded
    hdu.header["ACCEPTED"] = accepted
    hdu.header["CHI2"] = 0.0
    hdu.header["POLNAXIS"] = 0
    hdu.header["POLNGRP"] = 0
    hdu.header["PSF_FWHM"] = fwhm
    hdu.header["PSF_SAMP"] = samp
    hdu.header["PSFNAXIS"] = 3
    hdu.header["PSFAXIS1"] = nx
    hdu.header["PSFAXIS2"] = ny
    hdu.header["PSFAXIS3"] = 1
    fits.HDUList([fits.PrimaryHDU(), hdu, hdu]).writeto(out_path, overwrite=True)


def _extract_ldac_to_csv(ldac_path: Path, out_csv: Path) -> bool:
    for ext in ("#LDAC_OBJECTS", "#2", "#1", "#0"):
        try:
            subprocess.run(
                ["stilts", "tcopy", f"in={ldac_path}{ext}", f"out={out_csv}", "ofmt=csv"],
                check=True, capture_output=True,
            )
        except Exception:
            continue
        if out_csv.exists() and out_csv.stat().st_size > 0:
            with open(out_csv, newline="") as f:
                hdr = next(csv.reader(f), [])
            if len(hdr) > 2:
                return True
    return False


def load_survivors_with_ground_truth(tile_dir: Path) -> list[dict]:
    gt_csv = tile_dir / "catalogs" / "spread_model_postscore.csv"
    filt_csv = tile_dir / "catalogs" / "sextractor_pass2.filtered.csv"
    if not gt_csv.exists() or not filt_csv.exists():
        return []

    with open(filt_csv, newline="") as f:
        by_number = {r["NUMBER"]: r for r in csv.DictReader(f)}

    out = []
    with open(gt_csv, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("matched") != "True" or not row.get("SPREAD_MODEL"):
                continue
            num = row.get("survivor_number")
            src = by_number.get(num)
            if src is None:
                continue
            out.append({
                "survivor_number": num,
                "x_image": float(src["X_IMAGE"]),
                "y_image": float(src["Y_IMAGE"]),
                "gt_spread_model": float(row["SPREAD_MODEL"]),
            })
    return out


def crop_and_measure(raw_fits: Path, psf_path: Path, x: float, y: float,
                      halfwidth: int, scratch_dir: Path, synth_psf: bool = False) -> float | None:
    scratch_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = scratch_dir / "raw"
    raw_dir.mkdir(exist_ok=True)

    with fits.open(raw_fits) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        ny, nx = data.shape
        # X_IMAGE/Y_IMAGE are 1-indexed FITS pixel coords; array is 0-indexed.
        cx, cy = int(round(x)) - 1, int(round(y)) - 1
        x0, x1 = max(0, cx - halfwidth), min(nx, cx + halfwidth)
        y0, y1 = max(0, cy - halfwidth), min(ny, cy + halfwidth)
        stamp = data[y0:y1, x0:x1]
        stamp_hdr = header.copy()
        stamp_hdr["NAXIS1"] = stamp.shape[1]
        stamp_hdr["NAXIS2"] = stamp.shape[0]
        crop_path = raw_dir / "crop.fits"
        fits.PrimaryHDU(stamp, header=stamp_hdr).writeto(crop_path, overwrite=True)

    local_psf = scratch_dir / psf_path.name
    if synth_psf:
        build_synthetic_psf(psf_path, x, y, local_psf)
    elif not local_psf.exists():
        shutil.copy2(psf_path, local_psf)

    try:
        p2 = run_pass2(crop_path, scratch_dir, local_psf, config_root=CONFIG_ROOT)
    except RuntimeError:
        return None

    out_csv = scratch_dir / "crop_pass2.csv"
    if not _extract_ldac_to_csv(p2, out_csv):
        return None
    with open(out_csv, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None

    # crop center in crop-local (1-indexed) pixel coords
    local_cx = (x - 1 - x0) + 1
    local_cy = (y - 1 - y0) + 1
    best = min(rows, key=lambda r: (float(r["X_IMAGE"]) - local_cx) ** 2 +
               (float(r["Y_IMAGE"]) - local_cy) ** 2)
    d = ((float(best["X_IMAGE"]) - local_cx) ** 2 + (float(best["Y_IMAGE"]) - local_cy) ** 2) ** 0.5
    if d > 5:  # px; nearest detection isn't actually our survivor
        return None
    return float(best["SPREAD_MODEL"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tile-ids-file", type=Path, required=True)
    ap.add_argument("--tiles-root", type=Path, default=Path("./data/tiles_archive"))
    ap.add_argument("--halfwidths", type=str, default="150,300")
    ap.add_argument("--synth-psf", action="store_true",
                     help="Build a per-survivor position-correct constant PSF via GalSim "
                          "instead of reusing the tile's shared spatially-varying PSF.")
    ap.add_argument("--scratch-root", type=Path,
                     default=Path("/tmp/claude-1000/-home-janne-code-vasco60-repro-parity/"
                                  "e0081d61-bb22-4b2a-852e-fce4112158d2/scratchpad/cutout_validation"))
    args = ap.parse_args()

    halfwidths = [int(x) for x in args.halfwidths.split(",")]
    tile_ids = [l.strip() for l in args.tile_ids_file.read_text().splitlines() if l.strip()]

    results = []  # list of dicts: tile, survivor, halfwidth, gt, crop_val, edge_dist
    for i, tid in enumerate(tile_ids, 1):
        tile_dir = args.tiles_root / tid
        survivors = load_survivors_with_ground_truth(tile_dir)
        if not survivors:
            print(f"[{i}/{len(tile_ids)}] {tid}: no survivors with ground truth, skipping")
            continue

        raw_candidates = sorted((tile_dir / "raw").glob("*.fits"))
        psf_path = tile_dir / "postscore" / "pass1.psf"
        if not raw_candidates or not psf_path.exists():
            print(f"[{i}/{len(tile_ids)}] {tid}: missing raw fits or pass1.psf, skipping")
            continue
        raw_fits = raw_candidates[0]

        with fits.open(raw_fits) as hdul:
            ny, nx = hdul[0].data.shape

        print(f"[{i}/{len(tile_ids)}] {tid}: {len(survivors)} survivors, image {nx}x{ny}", flush=True)

        for s in survivors:
            edge_dist = min(s["x_image"], s["y_image"], nx - s["x_image"], ny - s["y_image"])
            for hw in halfwidths:
                scratch = args.scratch_root / tid / f"{s['survivor_number']}_{hw}"
                val = crop_and_measure(raw_fits, psf_path, s["x_image"], s["y_image"], hw, scratch,
                                        synth_psf=args.synth_psf)
                results.append({
                    "tile": tid, "survivor": s["survivor_number"], "halfwidth": hw,
                    "gt": s["gt_spread_model"], "crop": val, "edge_dist": edge_dist,
                })
                shutil.rmtree(scratch, ignore_errors=True)

    print("")
    print("===================== CUTOUT VALIDATION SUMMARY =====================")
    for hw in halfwidths:
        rows = [r for r in results if r["halfwidth"] == hw and r["crop"] is not None]
        n_total = len([r for r in results if r["halfwidth"] == hw])
        n_ok = len(rows)
        if not rows:
            print(f"halfwidth={hw}px: 0/{n_total} produced a match, skipping stats")
            continue
        deltas = np.array([r["crop"] - r["gt"] for r in rows])
        gt_pass = np.array([r["gt"] > SPREAD_MODEL_MIN for r in rows])
        crop_pass = np.array([r["crop"] > SPREAD_MODEL_MIN for r in rows])
        agree = (gt_pass == crop_pass)
        print(f"--- halfwidth={hw}px ({n_ok}/{n_total} matched) ---")
        print(f"  delta (crop-gt): mean={deltas.mean():.5f} median={np.median(deltas):.5f} "
              f"|delta| 95th pct={np.percentile(np.abs(deltas), 95):.5f} max={np.abs(deltas).max():.5f}")
        print(f"  gate agreement: {agree.sum()}/{len(agree)} ({100*agree.mean():.1f}%)")
        if (~agree).any():
            print("  DISAGREEMENTS:")
            for r in [r for r in rows if (r["crop"] > SPREAD_MODEL_MIN) != (r["gt"] > SPREAD_MODEL_MIN)]:
                print(f"    {r['tile']} survivor={r['survivor']} edge_dist={r['edge_dist']:.0f}px "
                      f"gt={r['gt']:.5f} crop={r['crop']:.5f}")
        # bucket by edge distance to check the PSF-variation-lookup theory
        ed = np.array([r["edge_dist"] for r in rows])
        near_edge = ed < np.median(ed)
        for label, mask in [("near-edge half", near_edge), ("near-center half", ~near_edge)]:
            if mask.sum() == 0:
                continue
            print(f"    {label}: mean|delta|={np.abs(deltas[mask]).mean():.5f} "
                  f"gate agreement={100*agree[mask].mean():.1f}% (n={mask.sum()})")
    print("=======================================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
