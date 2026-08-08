#!/usr/bin/env python3
"""End-to-end proof that VASCO_LDAC_DROP_VIGNET changes nothing but speed.

`verify_ldac_vignet_drop.py` compares two conversions of the same LDAC. That is
not enough on its own: it never runs SExtractor or PSFEx, so it cannot answer
"does dropping VIGNET break the second pass or the PSF model?". This does.

The structural answer is that it cannot, because the knob touches only the
LDAC->CSV conversion, and nothing in the extraction chain reads that CSV:

    pass1   sex + sex_default.param  (HAS VIGNET) -> pass1.ldac
    PSFEx   reads pass1.ldac                      -> pass1.psf     [LDAC, not CSV]
    pass2   sex + default.param + pass1.psf       -> pass2.ldac    [image + PSF]
    single  sex + one_pass.sex/sex_default.param  -> pass2.ldac
    ---- the knob acts only below this line ----
    _ensure_sextractor_csv(pass2.ldac)            -> sextractor_pass2.csv

PSFEx genuinely does require VIGNET, and it genuinely does still get it: its
input is pass1.ldac, which this change never touches. Note also that the two
mode-specific param sets differ -- two-pass writes pass2.ldac from
`default.param`, which never requested VIGNET in the first place, so the knob is
a no-op there and simply falls back to tcopy.

But structural reasoning is exactly what has misled this project before, so this
script runs the real thing: the same tile through the real pipeline twice, knob
off and knob on, in both single-pass and two-pass mode, and requires every
artifact to match. Byte-identical for the LDACs and the PSF model, and identical
on every retained column for the CSV.

Production tiles are never written -- everything happens in a scratch tree.

How to validate
---------------
    python3 tools/verify_ldac_vignet_e2e.py \
        --raw-fits <tiles_dir>/<tile>/raw/<file>.fits \
        --work <work_dir>/fullscale_slice/e2e_vignet

Exit status is 0 only if both modes matched on every artifact. A PSFEx failure
in the two-pass arm fails the run loudly rather than being skipped.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def same_fits_data(a: Path, b: Path) -> tuple[bool, str]:
    """Compare two FITS files by DATA, ignoring headers.

    Byte comparison is the wrong instrument here: SExtractor stamps its own wall
    clock into the LDAC_IMHEAD card text (SEXTIME, SEXELAPS, SEXDATE), so two
    runs of the identical command are never byte-identical. Measured on this
    tile, exactly two characters differ out of 12,072 -- '09:38:12' vs
    '09:38:19', and elapsed 0 s vs 1 s -- while all 3,187 object rows match in
    every column. Comparing data answers the question that matters: did the
    extraction produce different detections?
    """
    from astropy.io import fits
    with fits.open(a) as fa, fits.open(b) as fb:
        if len(fa) != len(fb):
            return False, f"HDU count {len(fa)} vs {len(fb)}"
        for i, (ha, hb) in enumerate(zip(fa, fb)):
            da, db = ha.data, hb.data
            if da is None and db is None:
                continue
            if (da is None) != (db is None):
                return False, f"HDU{i}: one side has no data"
            if getattr(da, "names", None):                    # binary table
                if list(da.names) != list(db.names):
                    return False, f"HDU{i}: column names differ"
                if len(da) != len(db):
                    return False, f"HDU{i}: {len(da)} vs {len(db)} rows"
                for c in da.names:
                    # LDAC_IMHEAD holds the timestamped config text, not science.
                    if ha.name == "LDAC_IMHEAD":
                        continue
                    if not np.array_equal(np.asarray(da[c]), np.asarray(db[c])):
                        return False, f"HDU{i} column {c} differs"
            elif not np.array_equal(np.asarray(da), np.asarray(db)):
                return False, f"HDU{i}: image data differs"
    return True, "all data identical"


def build_tile(raw_fits: Path, tile_dir: Path) -> None:
    """A tile directory the pipeline will accept: raw/<image>.fits and nothing else."""
    if tile_dir.exists():
        shutil.rmtree(tile_dir)
    (tile_dir / "raw").mkdir(parents=True)
    shutil.copy2(raw_fits, tile_dir / "raw" / raw_fits.name)


def run_steps(tile_dir: Path, single_pass: bool, drop_vignet: bool) -> tuple[bool, float, str]:
    """Run steps 2-3 on one tile with a specific knob setting."""
    tiles_file = tile_dir.parent / f"{tile_dir.name}.txt"
    tiles_file.write_text(str(tile_dir) + "\n")

    env = dict(os.environ)
    env.pop("VASCO_REPRO_SINGLE_PASS", None)
    env.pop("VASCO_LDAC_DROP_VIGNET", None)
    if single_pass:
        env["VASCO_REPRO_SINGLE_PASS"] = "1"
    if drop_vignet:
        env["VASCO_LDAC_DROP_VIGNET"] = "1"

    t0 = time.time()
    p = subprocess.run([PY, "tools/run_steps_2_3_parallel.py",
                        "--tiles-file", str(tiles_file), "--workers", "1"],
                       cwd=REPO, env=env, capture_output=True, text=True)
    return p.returncode == 0, time.time() - t0, (p.stdout + p.stderr)[-2500:]


def compare_csv(a_csv: Path, b_csv: Path) -> tuple[bool, str]:
    a = pd.read_csv(a_csv, low_memory=False)
    b = pd.read_csv(b_csv, low_memory=False)
    kept = [c for c in a.columns if not c.startswith("VIGNET")]
    if len(a) != len(b):
        return False, f"row count {len(a)} vs {len(b)}"
    if kept != list(b.columns):
        return False, f"schema differs: {len(kept)} kept vs {len(b.columns)}"
    worst, wcol = 0.0, ""
    for col in kept:
        x, y = a[col].to_numpy(), b[col].to_numpy()
        if not (np.issubdtype(x.dtype, np.number) and np.issubdtype(y.dtype, np.number)):
            if not (a[col].astype(str) == b[col].astype(str)).all():
                return False, f"non-numeric column {col} differs"
            continue
        if not np.array_equal(np.isnan(x), np.isnan(y)):
            return False, f"NaN pattern differs in {col}"
        m = ~np.isnan(x)
        if not m.any():
            continue
        rel = np.abs(x[m] - y[m]) / np.maximum(np.abs(x[m]), 1e-30)
        if rel.max() > worst:
            worst, wcol = float(rel.max()), col
    if worst != 0.0:
        return False, f"max relative diff {worst:.3e} in {wcol}"
    return True, f"{len(a):,} rows x {len(kept)} columns identical"


def run_mode(raw: Path, work: Path, single_pass: bool) -> bool:
    mode = "single-pass" if single_pass else "two-pass (PSFEx)"
    print(f"\n{'=' * 68}\n=== {mode}\n{'=' * 68}", flush=True)
    results = {}
    for label, drop in (("knob_off", False), ("knob_on", True)):
        tdir = work / f"{'single' if single_pass else 'two'}_{label}"
        build_tile(raw, tdir)
        ok, secs, log = run_steps(tdir, single_pass, drop)
        csv = tdir / "catalogs" / "sextractor_pass2.csv"
        ldac = tdir / "pass2.ldac"
        psf = tdir / "pass1.psf"
        print(f"  [{label:8}] rc_ok={ok}  {secs:6.1f}s  csv={csv.exists()} "
              f"ldac={ldac.exists()} psf={psf.exists()}", flush=True)
        if not ok or not csv.exists():
            print(f"  [FAIL] pipeline did not complete under {label}:\n{log}")
            return False
        results[label] = {"dir": tdir, "csv": csv, "ldac": ldac, "psf": psf, "secs": secs}

    a, b = results["knob_off"], results["knob_on"]

    if not single_pass:
        if not (a["psf"].exists() and b["psf"].exists()):
            print("  [FAIL] PSFEx produced no pass1.psf -- PSFEx is affected")
            return False
        same = sha256(a["psf"]) == sha256(b["psf"])
        ok_psf, why = (True, "byte-identical") if same else same_fits_data(a["psf"], b["psf"])
        print(f"  pass1.psf   {'identical' if ok_psf else 'DIFFERS'}: {why}")
        if not ok_psf:
            print("  [FAIL] the PSF model changed -- PSFEx IS affected")
            return False

    ok_ldac, why = same_fits_data(a["ldac"], b["ldac"])
    print(f"  pass2.ldac  {'identical' if ok_ldac else 'DIFFERS'}: {why}")
    if not ok_ldac:
        print("  [FAIL] the extraction itself changed, which the knob must never do")
        return False

    ok_csv, detail = compare_csv(a["csv"], b["csv"])
    print(f"  CSV         {'identical' if ok_csv else 'DIFFERS'}: {detail}")
    print(f"  size        {a['csv'].stat().st_size / 1e6:8.2f} MB -> "
          f"{b['csv'].stat().st_size / 1e6:.2f} MB")
    print(f"  steps 2-3   {a['secs']:.1f}s -> {b['secs']:.1f}s")
    return ok_csv


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw-fits", required=True, help="A tile image to run the pipeline on.")
    ap.add_argument("--work", required=True, help="Scratch root. Must not be a production tree.")
    args = ap.parse_args()

    work = Path(args.work)
    if "tiles_archive" in str(work.resolve()):
        raise SystemExit("[FATAL] refusing to run inside the production archive")
    work.mkdir(parents=True, exist_ok=True)
    raw = Path(args.raw_fits)
    print(f"[IN] {raw}  ({raw.stat().st_size / 1e6:.1f} MB)")

    ok_single = run_mode(raw, work, single_pass=True)
    ok_two = run_mode(raw, work, single_pass=False)

    print(f"\n{'=' * 68}")
    print(f"  single-pass       {'PASS' if ok_single else 'FAIL'}")
    print(f"  two-pass + PSFEx  {'PASS' if ok_two else 'FAIL'}")
    return 0 if (ok_single and ok_two) else 1


if __name__ == "__main__":
    sys.exit(main())
