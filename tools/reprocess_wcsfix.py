#!/usr/bin/env python3
"""Reprocess already-downloaded ESO tiles (from tools/compare_eso_vs_stsci.py)
through step2-5 with WCSFIX enabled, to test whether re-deriving each tile's
astrometry against the local Gaia cache (instead of trusting ESO's own
unrefined plate-solution WCS) fixes the veto-radius-miss mechanism found in
eso_vs_stsci_sextractor_comparison memory.

Reuses the already-fetched, already-POSS-I-validated raw FITS via symlink --
no new ESO network fetches. WCSFIX is exercised via its existing env-var
toggle (VASCO_WCSFIX_DISABLE unset) -- no pipeline code touched.

Usage:
  python3 tools/reprocess_wcsfix.py \
      --source-root work/eso_vs_stsci/20260802_unbiased100 \
      --out-dir work/eso_vs_stsci/20260802_unbiased100_wcsfix \
      --workers 12
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_eso_vs_stsci import (  # noqa: E402
    ARCHIVE,
    MATCH_ARCSEC,
    PIPELINE_ENV,
    REPO,
    load_radec,
    match_stats,
    tile_center,
)

WCSFIX_ENV = {k: v for k, v in PIPELINE_ENV.items() if k != "VASCO_WCSFIX_DISABLE"}


def run_pipeline_step(step: str, workdir: Path, extra_args: list[str] | None = None) -> None:
    import subprocess

    env = dict(os.environ)
    env.pop("VASCO_WCSFIX_DISABLE", None)  # ensure it can't leak through from the parent shell
    env.update(WCSFIX_ENV)
    cmd = [sys.executable, "-m", "vasco.cli_pipeline", step, "--workdir", str(workdir)]
    if extra_args:
        cmd += extra_args
    r = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        raise RuntimeError(f"{step} failed for {workdir.name}: {r.stderr[-1500:]}")


def process_one_tile(tile_id: str, source_root: Path, out_dir: Path) -> dict:
    ra, dec = tile_center(tile_id)
    result = {"tile_id": tile_id, "ra": ra, "dec": dec}

    src_fits = list((source_root / tile_id / "raw").glob("*.fits"))
    if not src_fits:
        result["error"] = "no_source_fits"
        return result
    src_fits = src_fits[0]

    stsci_dir = ARCHIVE / tile_id
    s_ra, s_dec = load_radec(stsci_dir / "catalogs" / "sextractor_pass2.csv")
    sf_ra, sf_dec = load_radec(stsci_dir / "catalogs" / "sextractor_pass2.filtered.csv")
    result["stsci_raw_n"] = len(s_ra)
    result["stsci_filtered_n"] = len(sf_ra)

    eso_orig_dir = source_root / tile_id / "catalogs"
    eo_ra, eo_dec = load_radec(eso_orig_dir / "sextractor_pass2.csv")
    eof_ra, eof_dec = load_radec(eso_orig_dir / "sextractor_pass2.filtered.csv")
    result["eso_orig_raw_n"] = len(eo_ra)
    result["eso_orig_filtered_n"] = len(eof_ra)

    sandbox = out_dir / tile_id
    raw_dir = sandbox / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    link_path = raw_dir / src_fits.name
    if not link_path.exists():
        link_path.symlink_to(src_fits.resolve())
    (sandbox / "RUN_INDEX.json").write_text(
        json.dumps([{"tile": src_fits.stem}]), encoding="utf-8"
    )

    try:
        run_pipeline_step("step2-pass1", sandbox)
        run_pipeline_step("step3-psf-and-pass2", sandbox)
        run_pipeline_step("step4-xmatch", sandbox, ["--size-arcmin", "60"])
        run_pipeline_step("step5-filter-within5", sandbox)
    except Exception as e:
        result["error"] = f"pipeline_failed: {e}"
        return result

    status_path = sandbox / "catalogs" / "wcsfix_status.json"
    if status_path.exists():
        try:
            status = json.loads(status_path.read_text())
            result["wcsfix_ok"] = bool(status.get("ok"))
            result["wcsfix_reason"] = status.get("reason", "")
            fit = status.get("fit") or {}
            result["wcsfix_sigma_arcsec"] = fit.get("sigma_arcsec")
            result["wcsfix_tie_points"] = status.get("tie_points")
        except Exception:
            result["wcsfix_ok"] = None
    else:
        result["wcsfix_ok"] = False
        result["wcsfix_reason"] = "no_status_file"

    # sextractor_pass2.csv is the PRE-WCSFIX raw file (no RA_corr/Dec_corr) --
    # the corrected coordinates only exist in sextractor_pass2.wcsfix.csv
    # (when WCSFIX succeeded) or fall back to the raw file otherwise.
    wcsfix_raw_path = sandbox / "catalogs" / "sextractor_pass2.wcsfix.csv"
    raw_path = wcsfix_raw_path if wcsfix_raw_path.exists() else sandbox / "catalogs" / "sextractor_pass2.csv"
    ew_ra, ew_dec = load_radec(raw_path)
    ewf_ra, ewf_dec = load_radec(sandbox / "catalogs" / "sextractor_pass2.filtered.csv")
    result["eso_wcsfix_raw_n"] = len(ew_ra)
    result["eso_wcsfix_filtered_n"] = len(ewf_ra)

    n_corr, sep_corr = match_stats(ew_ra, ew_dec, s_ra, s_dec, MATCH_ARCSEC, offset_correct=False)
    result["raw_matched_to_stsci"] = n_corr
    result["raw_median_sep_to_stsci_arcsec"] = sep_corr

    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-root", default="work/eso_vs_stsci/20260802_unbiased100")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--tiles", nargs="*", default=None, help="Subset of tile_ids (default: all in source-root)")
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    source_root = Path(args.source_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tiles = args.tiles if args.tiles else sorted(
        p.name for p in source_root.glob("tile_*") if p.is_dir()
    )
    print(f"[INFO] reprocessing {len(tiles)} tiles with WCSFIX enabled", file=sys.stderr)

    results = []
    if args.workers <= 1:
        for tid in tiles:
            res = process_one_tile(tid, source_root, out_dir)
            results.append(res)
            print(f"[RESULT] {res}", file=sys.stderr)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(process_one_tile, tid, source_root, out_dir): tid for tid in tiles}
            for fut in as_completed(futs):
                res = fut.result()
                results.append(res)
                print(f"[RESULT] {res}", file=sys.stderr)

    out_csv = out_dir / "results.csv"
    fieldnames = sorted({k for r in results for k in r.keys()})
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f"[INFO] wrote {out_csv}", file=sys.stderr)

    ok = [r for r in results if "error" not in r]
    errs = [r for r in results if "error" in r]
    tot_stsci = sum(r.get("stsci_filtered_n", 0) for r in ok)
    tot_orig = sum(r.get("eso_orig_filtered_n", 0) for r in ok)
    tot_wcsfix = sum(r.get("eso_wcsfix_filtered_n", 0) for r in ok)
    n_wcsfix_ok = sum(1 for r in ok if r.get("wcsfix_ok"))

    lines = [
        "# ESO WCSFIX reprocessing summary\n",
        f"Tiles: {len(results)} ({len(ok)} ok, {len(errs)} errors)",
        f"WCSFIX succeeded (ok=true): {n_wcsfix_ok}/{len(ok)}\n",
        "| metric | STScI | ESO original | ESO+WCSFIX |",
        "|---|---|---|---|",
        f"| total filtered survivors | {tot_stsci} | {tot_orig} | {tot_wcsfix} |",
    ]
    excess_orig = tot_orig - tot_stsci
    excess_wcsfix = tot_wcsfix - tot_stsci
    if excess_orig:
        reduction = 100 * (1 - excess_wcsfix / excess_orig)
        lines.append(f"\nExcess vs STScI: {excess_orig} (orig) -> {excess_wcsfix} (WCSFIX), "
                     f"{reduction:.1f}% reduction")
    if errs:
        lines.append(f"\nErrors: {len(errs)}")
        for r in errs:
            lines.append(f"- {r['tile_id']}: {r.get('error')}")

    out_md = out_dir / "SUMMARY.md"
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[INFO] wrote {out_md}", file=sys.stderr)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
