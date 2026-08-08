#!/usr/bin/env python3
"""Post-hoc SPREAD_MODEL scoring for single-pass survivor candidates.

The primary pipeline runs single-pass SExtractor (no PSFEx, no SPREAD_MODEL
column) so candidate generation stays permissive. This tool restores a real
SPREAD_MODEL value for just the small number of survivors in each tile's
catalogs/sextractor_pass2.filtered.csv, without paying for a full-field
second SExtractor pass (which costs ~88s/tile re-measuring thousands of
field objects to get SPREAD_MODEL for ~10-20 survivors).

Approach (validated via tools/validate_spread_model_cutout.py, 96.7% gate
agreement against real full-field ground truth across 60 survivors/10
plates -- see memory spread_model_cutout_speedup_validated):
  1. run_pass1 + run_psfex (unchanged, full field) build the tile's real,
     spatially-varying PSF model -- still needed as the source for step 2.
  2. For each survivor: evaluate that PSF model at the survivor's own
     (X_IMAGE, Y_IMAGE) via GalSim's DES_PSFEx reader, write it out as a
     small constant (degree-0) PSFEx-format .psf file, and run SExtractor's
     pass2 on a small crop centered on the survivor using that file.
     A naive crop reusing the tile's *shared* spatially-varying PSF file
     gets the wrong answer for ~5% of survivors, because SExtractor
     evaluates PSF-variation at the crop-local pixel position, not the
     true field position -- the per-survivor constant PSF sidesteps that
     bug entirely (confirmed: a flat/degree-0 PSF shrinks the crop-vs-
     full-field delta 7x). SPREAD_MODEL is scale-invariant to the PSF
     template's absolute flux normalization (verified empirically), so no
     calibration against PSFEx's internal normalization is needed.
  3. run_pass2 itself is untouched -- SExtractor still does the actual
     SPREAD_MODEL measurement; GalSim is only used to read/evaluate the
     PSF PSFEx already built.

Residual ~3.3% disagreement in validation was traced to crowded-field
deblending (a bright neighbor within ~22px), a different failure mode
this approach doesn't target -- accepted as a known limitation.

Purely additive: writes catalogs/spread_model_postscore.csv per tile and
never modifies sextractor_pass2.filtered.csv. The -0.002 threshold is the
locked SPREAD_MODEL gate from context/02_DECISIONS.md, applied here as a
reportable overlay, not a filter that changes the candidate list.

Usage:
  python tools/spread_model_postscore.py --tile-ids-file /path/to/ids.txt \
      [--tiles-root ./data/tiles] [--workers 4] [--crop-halfwidth 300]
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from vasco.pipeline_split import run_pass1, run_psfex, run_pass2  # noqa: E402
from validate_spread_model_cutout import build_synthetic_psf  # noqa: E402

SPREAD_MODEL_MIN = -0.002  # locked, context/02_DECISIONS.md
CONFIG_ROOT = str(REPO_ROOT / "configs")
CROP_HALFWIDTH_DEFAULT = 300  # px; validated size, see module docstring


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


def crop_and_measure(data, header, psf_path: Path, x: float, y: float,
                      halfwidth: int, scratch_dir: Path) -> tuple[float, float, float] | None:
    """Crop a stamp around (x, y), build a synthetic constant PSF at that
    exact position, run real SExtractor pass2 on the crop, and return
    (SPREAD_MODEL, SPREADERR_MODEL, sep_px) for the nearest detection to
    the crop center, or None if nothing usable was found."""
    scratch_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = scratch_dir / "raw"
    raw_dir.mkdir(exist_ok=True)

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

    local_psf = scratch_dir / "synth.psf"
    build_synthetic_psf(psf_path, x, y, local_psf)

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
    d_px = ((float(best["X_IMAGE"]) - local_cx) ** 2 + (float(best["Y_IMAGE"]) - local_cy) ** 2) ** 0.5
    if d_px > 5:  # px; nearest detection isn't actually our survivor
        return None
    sme = best.get("SPREADERR_MODEL", "")
    return float(best["SPREAD_MODEL"]), float(sme) if sme not in ("", None) else float("nan"), d_px


def postscore_tile(tile_id: str, tiles_root: Path, crop_halfwidth: int = CROP_HALFWIDTH_DEFAULT,
                    skip_existing: bool = False) -> dict:
    tile_dir = tiles_root / tile_id
    out_csv = tile_dir / "catalogs" / "spread_model_postscore.csv"
    if skip_existing and out_csv.exists() and out_csv.stat().st_size > 0:
        return {"tile_id": tile_id, "status": "skipped_existing", "survivors": 0}

    filtered_csv = tile_dir / "catalogs" / "sextractor_pass2.filtered.csv"
    if not filtered_csv.exists():
        return {"tile_id": tile_id, "status": "no_filtered_csv", "survivors": 0}

    with open(filtered_csv, newline="") as f:
        survivors = list(csv.DictReader(f))
    if not survivors:
        return {"tile_id": tile_id, "status": "no_survivors", "survivors": 0}

    raw_candidates = sorted((tile_dir / "raw").glob("*.fits"))
    if not raw_candidates:
        return {"tile_id": tile_id, "status": "no_raw_fits", "survivors": len(survivors)}
    raw_fits = raw_candidates[0]

    postscore_dir = tile_dir / "postscore"
    postscore_raw = postscore_dir / "raw"
    postscore_raw.mkdir(parents=True, exist_ok=True)
    link = postscore_raw / raw_fits.name
    if not link.exists():
        link.symlink_to(raw_fits.resolve())

    p1_ldac, _ = run_pass1(link, postscore_dir, config_root=CONFIG_ROOT)
    psf = run_psfex(p1_ldac, postscore_dir, config_root=CONFIG_ROOT)

    with fits.open(link) as hdul:
        data = hdul[0].data
        header = hdul[0].header

        out_rows = []
        n_matched = n_gate_pass = n_gate_fail = 0
        crops_dir = postscore_dir / "crops"
        for s in survivors:
            x_image = float(s["X_IMAGE"])
            y_image = float(s["Y_IMAGE"])
            s_ra = float(s.get("ALPHA_J2000") or s.get("ALPHAWIN_J2000"))
            s_dec = float(s.get("DELTA_J2000") or s.get("DELTAWIN_J2000"))
            row = {
                "survivor_number": s.get("NUMBER", ""),
                "ra": s_ra,
                "dec": s_dec,
                "matched": False,
                "sep_px": "",
                "SPREAD_MODEL": "",
                "SPREADERR_MODEL": "",
                "gate_pass": "",
            }
            scratch = crops_dir / str(s.get("NUMBER", ""))
            result = crop_and_measure(data, header, psf, x_image, y_image, crop_halfwidth, scratch)
            shutil.rmtree(scratch, ignore_errors=True)
            if result is not None:
                sm, sme, sep_px = result
                gate_pass = sm > SPREAD_MODEL_MIN
                row.update({
                    "matched": True,
                    "sep_px": round(sep_px, 3),
                    "SPREAD_MODEL": sm,
                    "SPREADERR_MODEL": sme,
                    "gate_pass": gate_pass,
                })
                n_matched += 1
                n_gate_pass += int(gate_pass)
                n_gate_fail += int(not gate_pass)
            out_rows.append(row)

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)

    return {
        "tile_id": tile_id,
        "status": "ok",
        "survivors": len(survivors),
        "matched": n_matched,
        "gate_pass": n_gate_pass,
        "gate_fail": n_gate_fail,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tile-ids-file", type=Path, required=True)
    ap.add_argument("--tiles-root", type=Path, default=Path("./data/tiles"))
    ap.add_argument("--crop-halfwidth", type=int, default=CROP_HALFWIDTH_DEFAULT,
                     help=f"Crop half-width in px around each survivor. Default {CROP_HALFWIDTH_DEFAULT}, "
                          "the size validated in tools/validate_spread_model_cutout.py.")
    ap.add_argument("--workers", type=int, default=1,
                     help="Tiles processed concurrently. Each tile's own PSFEx work doesn't "
                          "benefit from more threads (ATLAS/OpenBLAS single-threaded PSFEx build "
                          "in this environment), but independent tiles do parallelize. Default: 1")
    ap.add_argument("--skip-existing", action="store_true",
                     help="Skip a tile if catalogs/spread_model_postscore.csv already exists and "
                          "is non-empty, without touching PSFEx/SExtractor. Makes a large batch "
                          "resumable after an interruption without redoing completed tiles.")
    args = ap.parse_args()

    tile_ids = [l.strip() for l in args.tile_ids_file.read_text().splitlines() if l.strip()]

    t_start = time.time()
    results = []
    if args.workers <= 1:
        for i, tid in enumerate(tile_ids, 1):
            print(f"[{i}/{len(tile_ids)}] {tid} ...", flush=True)
            r = postscore_tile(tid, args.tiles_root, args.crop_halfwidth, args.skip_existing)
            results.append(r)
            print(f"  -> {r}", flush=True)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(postscore_tile, tid, args.tiles_root, args.crop_halfwidth,
                                    args.skip_existing): tid
                       for tid in tile_ids}
            for i, fut in enumerate(as_completed(futures), 1):
                tid = futures[fut]
                r = fut.result()
                results.append(r)
                print(f"[{i}/{len(tile_ids)}] {tid} -> {r}", flush=True)
    elapsed = time.time() - t_start

    total_survivors = sum(r.get("survivors", 0) for r in results)
    total_matched = sum(r.get("matched", 0) for r in results)
    total_pass = sum(r.get("gate_pass", 0) for r in results)
    total_fail = sum(r.get("gate_fail", 0) for r in results)
    skipped_tiles = [r["tile_id"] for r in results if r["status"] == "skipped_existing"]
    failed_tiles = [r["tile_id"] for r in results if r["status"] not in ("ok", "skipped_existing")]

    print("")
    print("===================== SPREAD_MODEL POSTSCORE SUMMARY =====================")
    print(f"Tiles processed:                  {len(tile_ids)} (workers={args.workers}, "
          f"crop_halfwidth={args.crop_halfwidth}px)")
    print(f"Wall time:                        {elapsed:.1f}s")
    if skipped_tiles:
        print(f"Tiles skipped (--skip-existing):  {len(skipped_tiles)}")
    if failed_tiles:
        print(f"Tiles with issues (see above):    {len(failed_tiles)} -> {failed_tiles}")
    print(f"Total survivors (pre-gate):        {total_survivors}")
    print(f"Matched to a two-pass detection:   {total_matched}")
    print(f"Would PASS SPREAD_MODEL>{SPREAD_MODEL_MIN} gate: {total_pass}")
    print(f"Would FAIL SPREAD_MODEL>{SPREAD_MODEL_MIN} gate: {total_fail}")
    print("============================================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
