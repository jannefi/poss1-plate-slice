#!/usr/bin/env python3
import argparse
import json
import os
import subprocess
import sys
import glob
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

def discover_tiles(root: Path) -> list[Path]:
    pat = str(root / "tile_RA*_DEC*")
    return [Path(p) for p in glob.glob(pat, recursive=True)]

def read_tiles_file(path: Path) -> list[Path]:
    tiles = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tiles.append(Path(line))
    return tiles

def rm_paths(paths: list[Path], dry_run: bool = False) -> None:
    for p in paths:
        if dry_run:
            print(f"[DRY] rm -f {p}")
        else:
            try:
                p.unlink(missing_ok=True)
            except Exception:
                pass

def rm_globs(tile: Path, patterns: list[str], dry_run: bool = False) -> None:
    for pat in patterns:
        for s in glob.glob(str(tile / pat)):
            rm_paths([Path(s)], dry_run=dry_run)

def run_cmd(cmd: list[str], env: dict[str,str]) -> int:
    p = subprocess.run(cmd, env=env)
    return p.returncode

def _read_tile_status(tile: Path) -> dict:
    try:
        p = tile / "tile_status.json"
        return json.loads(p.read_text(encoding="utf-8")).get("steps", {}) if p.exists() else {}
    except Exception:
        return {}


def _step_done(tile: Path, step: str) -> bool:
    """Return True if a single step is complete per tile_status.json."""
    return _read_tile_status(tile).get(step, {}).get("status", "") == "ok"


def _steps_done(tile: Path) -> bool:
    """Return True if step4 and step5 are both complete.

    Primary signal: tile_status.json step4.status==ok and step5.status==ok.
    Fallback: legacy marker files for tiles processed before tile_status.json.
    """
    steps = _read_tile_status(tile)
    if steps.get("step4", {}).get("status") == "ok" and steps.get("step5", {}).get("status") == "ok":
        return True
    # Fallback: legacy marker files
    xdir = tile / "xmatch"
    return (xdir / ".ok_step4_local").exists() and (xdir / ".ok_step5_within5").exists()


def process_one(tile: Path, args, env: dict[str,str]) -> tuple[str, Path, int, str]:
    if Path(".STOP").exists():
        return ("STOP", tile, 0, "STOP file present")

    xdir = tile / "xmatch"
    cdir = tile / "catalogs"
    xdir.mkdir(parents=True, exist_ok=True)
    cdir.mkdir(parents=True, exist_ok=True)

    # Skip logic (unless --force)
    if not args.force:
        if _steps_done(tile):
            return ("SKIP", tile, 0, "step4+5 complete")

    # Optional cleanup (avoid stale skip-fast behavior)
    if args.clean:
        if args.force or not _steps_done(tile):
            rm_globs(tile, [
                "xmatch/sex_*_xmatch*.csv",
                "catalogs/sextractor_pass2.after_*_veto.csv",
                "catalogs/sextractor_pass2.filtered.csv",
                "catalogs/sextractor_pass2.wcsfix.csv",
                "catalogs/wcsfix_status.json",
                "catalogs/_wcsfix_bootstrap_gaia.csv",
            ], dry_run=args.dry_run)

    # Step4
    step4_done = not args.force and _step_done(tile, "step4")
    if not step4_done:
        rc = run_cmd([
            sys.executable, "-m", "vasco.cli_pipeline", "step4-xmatch",
            "--workdir", str(tile),
            "--xmatch-radius-arcsec", str(args.radius_arcsec),
            "--size-arcmin", str(args.size_arcmin),
        ], env)
        if rc != 0:
            return ("FAIL_STEP4", tile, rc, "step4-xmatch failed")

    # Step5
    step5_done = not args.force and _step_done(tile, "step5")
    if not step5_done:
        if args.force:
            pass  # within5arcsec files no longer generated; nothing to clean
        rc = run_cmd([
            sys.executable, "-m", "vasco.cli_pipeline", "step5-filter-within5",
            "--workdir", str(tile),
        ], env)
        if rc != 0:
            return ("FAIL_STEP5", tile, rc, "step5-filter-within5 failed")

    return ("OK", tile, 0, "done")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiles-file", default="", help="Optional list of tile directories (one per line)")
    ap.add_argument("--tiles-root", default="./data/tiles", help="Tile root for auto-discovery")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--size-arcmin", type=float, default=60.0)
    ap.add_argument("--clean", action="store_true", help="remove stale outputs before processing")
    ap.add_argument("--force", action="store_true", help="re-run step4/5 even if markers exist")
    ap.add_argument("--only-missing", action="store_true", help="process tiles missing either marker")
    ap.add_argument("--dry-run", action="store_true", help="print deletes + commands without running")
    ap.add_argument("--wcsfix-fallback", action="store_true",
                    help="set env vars for WCSFIX fallback attempt (bootstrap 15\", min 10, deg 1)")
    args = ap.parse_args()

    root = Path(args.tiles_root)
    if args.tiles_file:
        tiles = read_tiles_file(Path(args.tiles_file))
    else:
        tiles = discover_tiles(root)

    tiles = [t for t in tiles if t.exists()]
    if not tiles:
        print("No tiles found.", file=sys.stderr)
        return 2

    env = os.environ.copy()

    # IMPORTANT: This only helps if you implement fallback in the pipeline
    # OR you run with a single config. For true "retry on failure" logic,
    # we need a code change in cli_pipeline / ensure_wcsfix_catalog wrapper.
    if args.wcsfix_fallback:
        env["VASCO_WCSFIX_BOOTSTRAP_ARCSEC"] = "15"
        env["VASCO_WCSFIX_MIN_MATCHES"] = "10"
        env["VASCO_WCSFIX_DEGREE"] = "1"

    print(f"[info] tiles={len(tiles)} workers={args.workers} radius={args.radius_arcsec}\" size={args.size_arcmin}' clean={args.clean} force={args.force}")
    if args.tiles_file:
        print(f"[info] tiles-file={args.tiles_file}")
    else:
        print(f"[info] auto-discovery from {root}")
    if args.wcsfix_fallback:
        print("[info] wcsfix-fallback env set: bootstrap=15\" min_matches=10 degree=1")

    ok = fail = skip = stop = 0

    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(process_one, t, args, env) for t in tiles]
        for fut in as_completed(futs):
            status, tile, rc, msg = fut.result()
            if status == "OK":
                ok += 1
            elif status == "SKIP":
                skip += 1
            elif status == "STOP":
                stop += 1
            else:
                fail += 1
            print(f"[{status}] {tile} rc={rc} {msg}")

            if status == "STOP":
                print("[info] STOP requested; exiting dispatcher loop.")
                break

    print(f"[summary] ok={ok} skip={skip} fail={fail} stop={stop} total={len(tiles)}")
    return 0 if fail == 0 else 1

if __name__ == "__main__":
    raise SystemExit(main())
