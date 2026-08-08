#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_plan.py — Plan-driven Step 1 orchestrator for VASCO60.

Reads a tile plan CSV (produced by vasco.plan.tessellate_plates) and calls
step1-download sequentially for each tile.

Resume behaviour
----------------
A tile is considered done if its tile_status.json records step1.status of
"ok" or "skip" (fallback: RUN_COUNTS.json exists for tiles downloaded before
tile_status.json was introduced).  Done tiles are logged as SKIP and do NOT
count toward --limit.

Usage
-----
    python scripts/run_plan.py plans/tiles_poss1e_ps1.csv
    python scripts/run_plan.py plans/tiles_smoke.csv --limit 1
    python scripts/run_plan.py plans/tiles_smoke.csv --limit 1 --dry-run
    python scripts/run_plan.py plans/tiles_poss1e_ps1.csv --plate XE309
    python scripts/run_plan.py plans/tiles_poss1e_ps1.csv --plate XE309 --dry-run
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import logging
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Logging setup — append to ./logs/run_plan.log
# ---------------------------------------------------------------------------
LOG_PATH = Path("./logs/run_plan.log")

def _setup_logger() -> logging.Logger:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("run_plan")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s %(message)s", datefmt="%Y-%m-%dT%H:%M:%S")

    fh = logging.FileHandler(LOG_PATH, mode="a", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    return logger


# ---------------------------------------------------------------------------
# Resume check
# ---------------------------------------------------------------------------

def _is_done(tile_id: str, tiles_dir: Path) -> bool:
    """Return True if step1 already completed successfully for this tile.

    Primary signal: tile_status.json has step1.status in {"ok", "skip"}.
    Fallback signal: RUN_COUNTS.json exists (tiles downloaded before
    tile_status.json was introduced).
    """
    tile_dir = tiles_dir / tile_id
    status_path = tile_dir / "tile_status.json"
    if status_path.exists():
        try:
            data = json.loads(status_path.read_text(encoding="utf-8"))
            s = data.get("steps", {}).get("step1", {}).get("status", "")
            if s in {"ok", "skip"}:
                return True
        except Exception:
            pass
    return (tile_dir / "RUN_COUNTS.json").exists()


# ---------------------------------------------------------------------------
# Step1 invocation
# ---------------------------------------------------------------------------

def _run_step1(ra: str, dec: str, size_arcmin: str, survey: str,
               tiles_dir: Path) -> int:
    """Call step1-download as a subprocess. Returns exit code."""
    cmd = [
        sys.executable, "-m", "vasco.cli_pipeline", "step1-download",
        "--ra", ra,
        "--dec", dec,
        "--size-arcmin", size_arcmin,
        "--survey", survey,
        "--workdir", str(tiles_dir),
    ]
    result = subprocess.run(cmd)
    return result.returncode


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def run(plan_path: Path, tiles_dir: Path, limit: int | None, dry_run: bool,
        logger: logging.Logger, plate_filter: str | None = None) -> None:
    if not plan_path.exists():
        logger.error(f"Plan file not found: {plan_path}")
        sys.exit(1)

    with plan_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if plate_filter is not None:
        rows = [r for r in rows if r.get("plate_id") == plate_filter]
        if not rows:
            logger.error(f"Plate '{plate_filter}' not found in {plan_path}")
            sys.exit(1)

    total = len(rows)
    downloaded = 0
    skipped    = 0
    failed     = 0

    start_msg = f"START plan={plan_path} total_rows={total} limit={limit} dry_run={dry_run}"
    if plate_filter is not None:
        start_msg += f" plate={plate_filter}"
    logger.info(start_msg)

    try:
        for i, row in enumerate(rows, start=1):
            tile_id    = row["tile_id"]
            plate_id   = row["plate_id"]
            ra         = row["ra_deg"]
            dec        = row["dec_deg"]
            size       = row.get("size_arcmin", "60")
            survey     = row.get("survey", "dss1-red")
            prefix     = f"[{i}/{total}] {tile_id}  plate={plate_id}"

            # Resume check
            if _is_done(tile_id, tiles_dir):
                logger.info(f"SKIP  {prefix}")
                skipped += 1
                continue

            # Limit applies to downloads only (not skips)
            if limit is not None and downloaded >= limit:
                logger.info(f"LIMIT reached ({limit}); stopping.")
                break

            if dry_run:
                logger.info(f"DRY   {prefix}  ra={ra} dec={dec}")
                downloaded += 1
                continue

            rc = _run_step1(ra, dec, size, survey, tiles_dir)
            if rc == 0:
                logger.info(f"OK    {prefix}")
                downloaded += 1
            else:
                logger.warning(f"FAIL  {prefix}  exit_code={rc}")
                failed += 1

    except KeyboardInterrupt:
        logger.info(
            f"INTERRUPT  downloaded={downloaded}  skipped={skipped}  "
            f"failed={failed}  visited={downloaded + skipped + failed}"
        )
        sys.exit(130)

    logger.info(
        f"DONE  downloaded={downloaded}  skipped={skipped}  "
        f"failed={failed}  total_visited={downloaded + skipped + failed}"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Drive step1-download from a tile plan CSV."
    )
    p.add_argument("plan", help="Path to plan CSV (e.g. plans/tiles_poss1e_ps1.csv)")
    p.add_argument("--tiles-dir", default="./data/tiles",
                   help="Tiles root directory passed as --workdir to step1 [./data/tiles]")
    p.add_argument("--limit", type=int, default=None,
                   help="Stop after N successful downloads (skips don't count)")
    p.add_argument("--dry-run", action="store_true",
                   help="Print actions without running step1-download")
    p.add_argument("--plate", metavar="PLATE_ID", default=None,
                   help="Download only tiles attributed to this plate_id in the plan "
                        "(e.g. --plate XE309); error if not found")
    args = p.parse_args(argv)

    logger = _setup_logger()
    run(
        plan_path=Path(args.plan),
        tiles_dir=Path(args.tiles_dir),
        limit=args.limit,
        dry_run=args.dry_run,
        logger=logger,
        plate_filter=args.plate,
    )


if __name__ == "__main__":
    main()
