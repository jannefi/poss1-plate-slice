
#!/usr/bin/env python3
"""
Run step2 and step3 per tile concurrently across tiles.

Usage:
  python scripts/run_steps_2_3_parallel.py --tiles-file /tmp/tiles.txt --workers 6

Tips:
- Increase --workers if CPU and disk allow; decrease if the disk becomes the bottleneck.
- This script prints concise per-tile status and a final summary.
- Drop a ".STOP" file in the cwd to pause between tiles; rerunning later skips
  tiles whose tile_status.json already shows step2+step3 as "ok" (mirrors the
  resume pattern in run_steps_4_5_parallel.py). Pass --force to redo anyway.
"""
import argparse, json, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

def run(cmd: list[str]) -> int:
    # Stream minimal output; rely on tile logs for detailed info
    try:
        return subprocess.run(cmd, check=False).returncode
    except Exception:
        return 1

def _read_tile_status(tile: Path) -> dict:
    try:
        p = tile / "tile_status.json"
        return json.loads(p.read_text(encoding="utf-8")).get("steps", {}) if p.exists() else {}
    except Exception:
        return {}

def _steps_done(tile: Path) -> bool:
    """True if step2 and step3 are both already "ok" per tile_status.json."""
    steps = _read_tile_status(tile)
    return steps.get("step2", {}).get("status") == "ok" and steps.get("step3", {}).get("status") == "ok"

def process_tile(tile: str, force: bool) -> tuple[str, bool, str]:
    if Path(".STOP").exists():
        return (tile, True, f"{Path(tile).name}: STOP file present, skipped")
    tile_path = Path(tile)
    if not force and _steps_done(tile_path):
        return (tile, True, f"{Path(tile).name}: SKIP (step2+step3 already ok)")
    t0 = time.time()
    tile = str(tile_path.resolve())
    # Step 2
    rc2 = run([sys.executable,"-u","-m","vasco.cli_pipeline","step2-pass1","--workdir",tile])
    if rc2 not in (0,):  # 2 can be "missing raw" or similar; treat non-zero as soft fail
        # Still attempt step 3; gating may skip it
        pass
    # Step 3
    rc3 = run([sys.executable,"-u","-m","vasco.cli_pipeline","step3-psf-and-pass2","--workdir",tile])
    ok = (rc3 == 0)
    dt = time.time() - t0
    msg = f"{Path(tile).name}: step2 rc={rc2}, step3 rc={rc3}, {dt:.1f}s"
    return tile, ok, msg

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiles-file", required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--force", action="store_true", help="re-run step2/3 even if tile_status.json shows ok")
    args = ap.parse_args()

    tiles = [line.strip() for line in Path(args.tiles_file).read_text().splitlines() if line.strip()]
    print(f"[2+3] Tiles: {len(tiles)}, workers={args.workers}")
    ok_n = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_tile, t, args.force): t for t in tiles}
        for i, fut in enumerate(as_completed(futs), 1):
            tile, ok, msg = fut.result()
            print(f"[{i:>5}/{len(tiles)}] {msg}")
            if ok: ok_n += 1
    print(f"[2+3] Done. OK tiles ~{ok_n}/{len(tiles)}. See per-tile logs for details.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
