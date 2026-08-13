#!/usr/bin/env python3
"""Full-scale plate slicing: all 635 plates, lean output, bounded disk.

Answers the independent-replication question at scale. The 8-plate pilot showed
the union of full-plate tiles and archive tiles beats the archive alone by +7.1
pts @5" and +15.6 @10", and that ~50% of its residual gap was simply the
lower-bound artifact of processing 8 plates out of 635. This processes all of
them so the artifact disappears and the answer is real.

Why it does not run steps 4-5
------------------------------
Raw detection parity needs only `sextractor_pass2.csv`. Steps 4-5 build the
Gaia/PS1/USNO veto chain, which writes four more near-copies of the detection
table plus cached neighbourhood catalogues -- 251 MB per tile, which over 31,115
tiles is **7.44 TB**, more than the free space on the data disk. Steps 2-3 alone
plus immediate extraction keeps the whole run in a few GB. Veto stages can be
re-run later on whatever subset matters; they are not needed to measure recall.

Disk discipline
---------------
One plate is resident at a time: slice -> SExtractor -> extract lean RA/Dec ->
delete the tile tree. Peak usage is one plate's tiles (~4 GB); the surviving
output is ~1.5 MB per plate. Without the delete this run cannot fit on disk at
all, so the cleanup is load-bearing, not tidiness.

Resumable: a plate whose output CSV already exists is skipped, so an interrupted
run continues where it stopped.

How to validate
---------------
    python3 tools/run_fullscale_slice.py --out-dir work/runs/fullscale_slice \
        --plate-dir <plate_dir> --workers 12

Check `progress.csv`: every plate must report 49 tiles sliced, 49 with catalogs,
and **0 skips**. A plate reporting fewer has failed and should be re-run; the run
does not stop for it, because one bad plate should not cost the other 634. Each
plate's slicer output is kept at `slice_logs/<PLATE>.log` -- the first version of
this runner discarded it, and 195 tiles lost to a grid that walked off the array
went unnoticed for weeks as a result.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable
# Caches must be present on EVERY step. A missing one does not fail -- it
# silently falls back to live VizieR/MAST, which at survey scale means a run
# that appears to work and takes days.
REQUIRED_ENV_VETO = ["VASCO_GAIA_CACHE", "VASCO_PS1_CACHE", "VASCO_USNOB_CACHE"]


def plate_list(plate_dir: Path, only: str | None):
    if only:
        return [p.strip() for p in only.split(",") if p.strip()]
    rx = re.compile(r"dss1red_(XE\d+)\.fits$")
    out = []
    for f in sorted(plate_dir.glob("dss1red_XE*.fits")):
        m = rx.search(f.name)
        if m:
            out.append(m.group(1))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--plate-dir", default="<plate_dir>")
    ap.add_argument("--plates", default=None, help="Comma-separated subset; default all on disk.")
    ap.add_argument("--plate-manifest", default=None,
                    help="CSV with a plate_id column restricting to VASCO's own plate set.")
    ap.add_argument("--grid", type=int, default=7)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--crpix-table", default=None,
                    help="Per-plate CRPIX correction table; passed through to the slicer. "
                         "Without it ~33%% of plates land ~2.4\" from Gaia.")
    ap.add_argument("--with-vetoes", action="store_true",
                    help="Run steps 4-5 per plate while its tiles are still on disk, and "
                         "keep each tile's filtered survivors. Needed to produce "
                         "candidates; a recall-only run does not want this.")
    ap.add_argument("--keep-tiles", action="store_true",
                    help="Do NOT delete tile trees. Needs ~7 TB for a full run; do not use.")
    args = ap.parse_args()

    # Pass mode is explicit, never inherited. Single-pass produces no PSF model
    # and therefore no SPREAD_MODEL, and the MNRAS morphology gate then rejects
    # every candidate -- fine for a recall-only run, useless for producing one.
    single = bool(os.environ.get("VASCO_REPRO_SINGLE_PASS"))
    if args.with_vetoes and single:
        raise SystemExit(
            "[FATAL] --with-vetoes needs two-pass SExtractor, but "
            "VASCO_REPRO_SINGLE_PASS is set. Single-pass has no SPREAD_MODEL, so "
            "the morphology gate would reject every candidate. Unset it.")
    if args.with_vetoes:
        missing = [e for e in REQUIRED_ENV_VETO if not os.environ.get(e)]
        if missing:
            raise SystemExit(f"[FATAL] --with-vetoes needs local caches; unset: {missing}")
    print(f"[MODE] {'single-pass' if single else 'two-pass'}, "
          f"vetoes {'ON' if args.with_vetoes else 'OFF'}")

    # Echo the environment-driven settings. These change the science silently:
    # a stray VASCO_CIRCLE_ARCMIN inherited from an interactive shell once cut
    # every tile to a 30' circle for 106 plates of a full-scale run, discarding
    # ~21% of detections and reintroducing the corner gaps the square-tile
    # design exists to avoid. Nothing errored; the numbers were simply wrong.
    circle = os.environ.get("VASCO_CIRCLE_ARCMIN", "") or "off"
    # WCSFIX defaults to ON inside the pipeline, so its absence from this banner
    # is exactly how the 642-plate run of 2026-08 ran WCS-fixed on all 642
    # plates while the docs described raw plate WCS. Print the state; never
    # leave it to be inferred from a default.
    wcsfix = "off" if os.environ.get("VASCO_WCSFIX_DISABLE") else "ON"
    print(f"[CONFIG] circle_cut={circle}  wcsfix={wcsfix}  "
          f"crpix_table={args.crpix_table or 'NONE'}  "
          f"drop_vignet={os.environ.get('VASCO_LDAC_DROP_VIGNET', '0')}")
    if args.with_vetoes:
        for e in REQUIRED_ENV_VETO:
            print(f"[CONFIG] {e}={os.environ.get(e)}")
    if circle != "off":
        print("[CONFIG][WARN] a circular cut is ACTIVE -- README documents square "
              "tiles with no 30' cut. Unset VASCO_CIRCLE_ARCMIN unless you mean it.",
              flush=True)
    if wcsfix == "ON":
        print("[CONFIG][WARN] WCSFIX is ACTIVE -- deviation #6 documents raw plate "
              "WCS, and 02_DECISIONS.md pairs WCS-fixed coordinates with a 0.25\" "
              "dedup, not the 3.0\" in use. Export VASCO_WCSFIX_DISABLE=1 unless "
              "you mean it.", flush=True)

    out = Path(args.out_dir)
    if "tiles_archive" in str(out.resolve()):
        raise SystemExit("[FATAL] refusing to run inside the production archive")
    radec = out / "radec"
    radec.mkdir(parents=True, exist_ok=True)
    filtered = out / "filtered"
    if args.with_vetoes:
        filtered.mkdir(parents=True, exist_ok=True)
    scratch = out / "scratch_tiles"

    plates = plate_list(Path(args.plate_dir), args.plates)
    if args.plate_manifest:
        want = set(pd.read_csv(args.plate_manifest).plate_id.astype(str))
        plates = [p for p in plates if p in want]
    print(f"[PLAN] {len(plates)} plates, grid {args.grid}x{args.grid}, workers {args.workers}")

    # The slicer's stdout used to be captured and thrown away, so 195 tiles lost
    # to "Arrays do not overlap" passed silently and only turned up in a tile
    # count weeks later. Keep it: one log per plate, and the skip count in the
    # ledger so a shortfall is visible without reading logs at all.
    logs = out / "slice_logs"
    logs.mkdir(parents=True, exist_ok=True)

    prog = out / "progress.csv"
    if not prog.exists():
        prog.write_text("plate,tiles_sliced,tiles_with_catalogs,detections,skips,survivors,seconds,status\n")

    t_start = time.time()
    for i, plate in enumerate(plates, 1):
        dest = radec / f"{plate}.csv"
        if dest.exists():
            print(f"[{i}/{len(plates)}] {plate} SKIP (done)", flush=True)
            continue
        t0 = time.time()
        if scratch.exists():
            shutil.rmtree(scratch, ignore_errors=True)
        tiles_file = out / "_tiles.txt"

        proc = subprocess.run(
            [PY, "tools/slice_plate_tiles.py",
             "--plate-fits", f"{args.plate_dir}/dss1red_{plate}.fits",
             "--tiles-dir", str(scratch), "--grid", str(args.grid),
             "--tiles-file-out", str(tiles_file)]
            + (["--crpix-table", args.crpix_table] if args.crpix_table else []),
            cwd=REPO, capture_output=True, text=True)
        (logs / f"{plate}.log").write_text(proc.stdout + proc.stderr)
        n_skip = proc.stdout.count("[SKIP]")
        if proc.returncode != 0:
            print(f"[{i}/{len(plates)}] {plate} SLICE FAILED rc={proc.returncode} "
                  f"(see {logs / (plate + '.log')})", flush=True)
            with prog.open("a") as f:
                f.write(f"{plate},0,0,0,0,0,{time.time()-t0:.0f},slice_failed\n")
            continue
        n_sliced = len([l for l in tiles_file.read_text().split() if l])
        if n_skip or n_sliced != args.grid ** 2:
            print(f"[{i}/{len(plates)}] {plate} [WARN] {n_sliced}/{args.grid**2} tiles, "
                  f"{n_skip} skipped -- grid is walking off the array", flush=True)

        # Keep the step runners' output. It used to be captured and dropped, so
        # a plate that failed every tile reported "0 cat" with no reason -- the
        # traceback that said why (a missing import, in one real case) existed
        # and was thrown away. One log per plate costs kilobytes.
        steps_log = logs / f"{plate}.steps.log"
        with steps_log.open("w") as sf:
            sf.write(f"=== steps 2+3 :: {plate} ===\n")
            sf.flush()
            subprocess.run([PY, "tools/run_steps_2_3_parallel.py",
                            "--tiles-file", str(tiles_file), "--workers", str(args.workers)],
                           cwd=REPO, stdout=sf, stderr=subprocess.STDOUT, text=True)
            if args.with_vetoes:
                sf.write(f"\n=== steps 4+5 :: {plate} ===\n")
                sf.flush()
                subprocess.run([PY, "tools/run_steps_4_5_parallel.py",
                                "--tiles-file", str(tiles_file), "--workers", str(args.workers)],
                               cwd=REPO, stdout=sf, stderr=subprocess.STDOUT, text=True)

        # --- extract lean RA/Dec, then drop the heavy tree -----------------
        rows, n_cat = [], 0
        for tdir in [Path(l) for l in tiles_file.read_text().split() if l]:
            c = tdir / "catalogs" / "sextractor_pass2.csv"
            if not c.exists():
                continue
            try:
                d = pd.read_csv(c, usecols=["ALPHA_J2000", "DELTA_J2000", "MAG_AUTO"])
            except Exception:
                continue
            d = d.rename(columns={"ALPHA_J2000": "ra", "DELTA_J2000": "dec", "MAG_AUTO": "mag"})
            d["tile_id"] = tdir.name
            rows.append(d)
            n_cat += 1
        # Survivors are harvested into a tiles_root-shaped tree so that
        # scripts/build_run_stage_csvs.py can run over it unchanged after the
        # heavy trees are gone. Only the filtered CSV is kept: a few tens of KB
        # per tile against ~250 MB for the tile itself.
        n_surv = 0
        if args.with_vetoes:
            for tdir in [Path(l) for l in tiles_file.read_text().split() if l]:
                src = tdir / "catalogs" / "sextractor_pass2.filtered.csv"
                if not src.exists():
                    continue
                dst = filtered / tdir.name / "catalogs" / "sextractor_pass2.filtered.csv"
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, dst)
                for extra in ("tile_status.json",):
                    e = tdir / extra
                    if e.exists():
                        shutil.copyfile(e, filtered / tdir.name / extra)
                n_surv += 1

        n_det = 0
        if rows:
            allr = pd.concat(rows, ignore_index=True)
            allr["plate_id"] = plate
            allr.to_csv(dest, index=False)
            n_det = len(allr)
        if not args.keep_tiles:
            shutil.rmtree(scratch, ignore_errors=True)

        dt = time.time() - t0
        status = "ok" if (n_cat == n_sliced == args.grid ** 2 and n_det) else "partial"
        if args.with_vetoes and n_surv != n_sliced:
            status = "partial"
        with prog.open("a") as f:
            f.write(f"{plate},{n_sliced},{n_cat},{n_det},{n_skip},{n_surv},{dt:.0f},{status}\n")
        done = i
        eta = (time.time() - t_start) / max(done, 1) * (len(plates) - done) / 3600.0
        print(f"[{i}/{len(plates)}] {plate} {status}: {n_sliced} sliced, {n_cat} cat, "
              f"{n_det:,} det, "
              + (f"{n_surv} filt, " if args.with_vetoes else "")
              + f"{dt:.0f}s   ETA {eta:.1f}h", flush=True)

    print(f"\n[DONE] {time.time()-t_start:.0f}s total; lean output under {radec}")


if __name__ == "__main__":
    main()
