#!/usr/bin/env python3
"""
run_detection_totals.py — summarise raw detection counts for a run.

Reads tile_manifest.csv from the run folder, then for each tile reads
MNRAS_SUMMARY.json to get veto_start_rows (= pass2 unfiltered detections).
Also sums the manifest columns rows_in_tile_filtered_csv and rows_emitted_to_S0
so you can see the full funnel from raw detections → post-pipeline input.

Also reports per-plate survivor counts from the final post-pipeline stage CSV
(auto-detected as the last stage_S*.csv in stages/, or set with --final-stage).

Usage:
    python tools/run_detection_totals.py --run-dir work/runs/<run>
    python tools/run_detection_totals.py --run-dir work/runs/<run> --tiles-root data/tiles
    python tools/run_detection_totals.py --run-dir work/runs/<run> --final-stage stage_S5_VSX.csv
"""

import argparse
import collections
import csv
import json
import sys
from pathlib import Path


def _detect_final_stage(stages_dir):
    """Return the last stage_S*.csv (excluding _flags files) sorted by name."""
    candidates = sorted(
        f for f in stages_dir.glob("stage_S*.csv")
        if "_flags" not in f.name and "_ledger" not in f.name
    )
    return candidates[-1] if candidates else None


def main():
    ap = argparse.ArgumentParser(description="Summarise raw detection totals for a run")
    ap.add_argument("--run-dir", default=None, help="Run folder (contains tile_manifest.csv). Not required for --scan-all-tiles.")
    ap.add_argument(
        "--tiles-root",
        default="data/tiles",
        help="Root of tile folders (default: data/tiles)",
    )
    ap.add_argument(
        "--final-stage",
        default=None,
        help="Final stage CSV filename in stages/ (default: auto-detect last stage_S*.csv)",
    )
    ap.add_argument(
        "--list-ps1-truncated",
        action="store_true",
        help=(
            "Print tile dirs (one per line, to stdout) where catalogs/ps1_neighbourhood.csv "
            "has exactly 50000 rows — the old cap before it was raised to 200K. "
            "Pipe to a re-fetch script to update affected tiles."
        ),
    )
    ap.add_argument(
        "--clear-post1-for-truncated",
        action="store_true",
        help=(
            "Clear post1.status from tile_status.json for all PS1-truncated tiles "
            "(data rows == 50000). This lets the next delta run of build_run_stage_csvs.py "
            "pick them up as active rather than skipping them."
        ),
    )
    ap.add_argument(
        "--scan-all-tiles",
        action="store_true",
        help=(
            "Scan tiles-root directly (no manifest) for PS1-truncated tiles. "
            "Prints tile dirs to stdout. Does not require --run-dir. "
            "Use this to find truncated tiles that were delta-skipped in all run manifests."
        ),
    )
    ap.add_argument(
        "--clear-post1-from-stdin",
        action="store_true",
        help=(
            "Read tile dirs from stdin (one per line) and clear post1.status from each "
            "tile_status.json. Use after re-fetching, when ps1_neighbourhood.csv is no "
            "longer at the old 50K cap and --clear-post1-for-truncated would find nothing. "
            "Example: cat /tmp/truncated_tiles.txt /tmp/remaining_truncated.txt | sort -u | "
            "python tools/run_detection_totals.py --clear-post1-from-stdin"
        ),
    )
    args = ap.parse_args()

    run_dir = Path(args.run_dir) if args.run_dir else None
    tiles_root = Path(args.tiles_root)

    # --scan-all-tiles: scan tiles_root directly, no manifest needed
    if args.scan_all_tiles:
        OLD_CAP = 50000
        truncated = []
        all_tile_dirs = sorted(p for p in tiles_root.glob("tile_RA*") if p.is_dir())
        for td in all_tile_dirs:
            ps1_csv = td / "catalogs" / "ps1_neighbourhood.csv"
            if not ps1_csv.exists():
                continue
            with open(ps1_csv) as f:
                n_lines = sum(1 for _ in f)
            if n_lines - 1 == OLD_CAP:
                truncated.append(td)
        print(f"# PS1-truncated tiles (data rows == {OLD_CAP}): {len(truncated)} of {len(all_tile_dirs)} total", file=sys.stderr)
        for td in truncated:
            print(td)
        sys.exit(0)

    # --clear-post1-from-stdin: read tile dirs from stdin, clear post1 status
    if args.clear_post1_from_stdin:
        tile_dirs = [line.strip() for line in sys.stdin if line.strip()]
        cleared = 0
        for td in tile_dirs:
            status_path = Path(td) / "tile_status.json"
            if not status_path.exists():
                print(f"  SKIP (no tile_status.json): {td}", file=sys.stderr)
                continue
            try:
                data = json.loads(status_path.read_text(encoding="utf-8"))
                data.get("steps", {}).pop("post1", None)
                status_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
                print(f"  cleared: {Path(td).name}")
                cleared += 1
            except Exception as e:
                print(f"  ERROR {td}: {e}", file=sys.stderr)
        print(f"Done. {cleared}/{len(tile_dirs)} tiles cleared. Re-run build_run_stage_csvs.py in delta mode.")
        sys.exit(0)

    # Shared helper: find tiles where PS1 cache was truncated at the old 50K cap
    def _find_ps1_truncated(manifest_path, tiles_root):
        if not manifest_path.exists():
            print(f"ERROR: manifest not found: {manifest_path}", file=sys.stderr)
            sys.exit(1)
        with open(manifest_path, newline="") as f:
            all_rows = list(csv.DictReader(f))

        def _int_local(val):
            return int(val) if val and val.strip() else 0

        active = [r for r in all_rows if _int_local(r.get("skipped_delta", "0")) == 0]
        OLD_CAP = 50000
        truncated = []
        for row in active:
            ps1_csv = tiles_root / row["tile_id"] / "catalogs" / "ps1_neighbourhood.csv"
            if not ps1_csv.exists():
                continue
            with open(ps1_csv) as f:
                n_lines = sum(1 for _ in f)
            if n_lines - 1 == OLD_CAP:
                truncated.append(tiles_root / row["tile_id"])
        return truncated, len(active)

    # --list-ps1-truncated: emit tile dirs to stdout (pipeable)
    if args.list_ps1_truncated:
        truncated, n_active = _find_ps1_truncated(run_dir / "tile_manifest.csv", tiles_root)
        print(f"# PS1-truncated tiles (data rows == 50000): {len(truncated)} of {n_active} active", file=sys.stderr)
        for td in truncated:
            print(td)
        sys.exit(0)

    # --clear-post1-for-truncated: remove post1.status so next delta run reprocesses them
    if args.clear_post1_for_truncated:
        truncated, n_active = _find_ps1_truncated(run_dir / "tile_manifest.csv", tiles_root)
        print(f"Clearing post1 status for {len(truncated)} of {n_active} active PS1-truncated tiles...")
        cleared = 0
        for td in truncated:
            status_path = td / "tile_status.json"
            if not status_path.exists():
                print(f"  SKIP (no tile_status.json): {td.name}")
                continue
            try:
                data = json.loads(status_path.read_text(encoding="utf-8"))
                data.get("steps", {}).pop("post1", None)
                status_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
                print(f"  cleared: {td.name}")
                cleared += 1
            except Exception as e:
                print(f"  ERROR {td.name}: {e}", file=sys.stderr)
        print(f"Done. {cleared}/{len(truncated)} tiles cleared. Re-run build_run_stage_csvs.py in delta mode.")
        sys.exit(0)

    if run_dir is None:
        print("ERROR: --run-dir is required for this mode", file=sys.stderr)
        sys.exit(1)

    manifest_path = run_dir / "tile_manifest.csv"
    if not manifest_path.exists():
        print(f"ERROR: manifest not found: {manifest_path}", file=sys.stderr)
        sys.exit(1)

    with open(manifest_path, newline="") as f:
        rows = list(csv.DictReader(f))

    def _int(val):
        """Parse int, treating blank/missing as 0."""
        return int(val) if val and val.strip() else 0

    # tile_id → plate_id lookup (all tiles, including skipped)
    tile_to_plate = {r["tile_id"]: r.get("plate_id", "").strip() for r in rows}

    active_rows = [r for r in rows if _int(r.get("skipped_delta", "0")) == 0]
    skipped_rows = [r for r in rows if _int(r.get("skipped_delta", "0")) != 0]

    n_tiles = len(rows)
    n_active = len(active_rows)
    n_skipped = len(skipped_rows)

    pass2_filtered_total = sum(_int(r["rows_in_tile_filtered_csv"]) for r in active_rows)
    emitted_to_s0_total = sum(_int(r["rows_emitted_to_S0"]) for r in active_rows)

    active_plates = sorted({tile_to_plate.get(r["tile_id"], "") for r in active_rows if tile_to_plate.get(r["tile_id"], "")})
    n_plates = len(active_plates)

    pass2_unfiltered_total = 0
    missing_summary = []
    tiles_with_summary = 0

    for row in active_rows:
        tile_id = row["tile_id"]
        summary_path = tiles_root / tile_id / "MNRAS_SUMMARY.json"

        if not summary_path.exists():
            missing_summary.append(tile_id)
            continue

        with open(summary_path) as f:
            summary = json.load(f)

        pass2_unfiltered_total += summary.get("veto_start_rows", 0)
        tiles_with_summary += 1

    # --- Print header ---
    print(f"\nRun: {run_dir.name}")
    print(f"Tiles in manifest : {n_tiles}  ({n_active} active, {n_skipped} delta-skipped)")
    print(f"Plates covered    : {n_plates}  ({', '.join(active_plates)})")
    print(f"Tiles with summary: {tiles_with_summary}")
    if missing_summary:
        print(f"  Missing MNRAS_SUMMARY.json: {len(missing_summary)} tile(s)")
        for t in missing_summary[:10]:
            print(f"    {t}")
        if len(missing_summary) > 10:
            print(f"    ... and {len(missing_summary) - 10} more")

    # --- Funnel ---
    print()
    print("Detection funnel:")
    print(f"  pass2 unfiltered (veto_start_rows)  : {pass2_unfiltered_total:>10,}")
    print(f"  pass2 filtered (post-MNRAS gates)   : {pass2_filtered_total:>10,}")
    print(f"  emitted to S0 (stage input)         : {emitted_to_s0_total:>10,}")

    if pass2_unfiltered_total > 0:
        gate_pct = 100.0 * (1 - pass2_filtered_total / pass2_unfiltered_total)
        s0_pct = 100.0 * (1 - emitted_to_s0_total / pass2_unfiltered_total)
        print()
        print(f"  MNRAS gates removed               : {gate_pct:.3f}% of pass2 unfiltered")
        print(f"  Surviving to S0                   : {100 - s0_pct:.4f}% of pass2 unfiltered")

    print()
    print("Note: pass1 detection count is not available (stored as binary .ldac, no CSV).")

    # --- Per-plate survivor breakdown from final stage CSV ---
    stages_dir = run_dir / "stages"
    if not stages_dir.exists():
        return

    if args.final_stage:
        final_csv = stages_dir / args.final_stage
    else:
        final_csv = _detect_final_stage(stages_dir)

    if final_csv is None or not final_csv.exists():
        print(f"\n(No final stage CSV found in {stages_dir})", file=sys.stderr)
        return

    # Count survivors per plate; src_id format: tile_id:object_id
    plate_survivors = collections.Counter()
    plate_tiles_seen = collections.defaultdict(set)
    total_survivors = 0

    with open(final_csv, newline="") as f:
        for src_row in csv.DictReader(f):
            src_id = src_row.get("src_id", "")
            tile_id = src_id.split(":")[0] if ":" in src_id else src_row.get("tile_id", "")
            plate_id = tile_to_plate.get(tile_id, "unknown")
            plate_survivors[plate_id] += 1
            plate_tiles_seen[plate_id].add(tile_id)
            total_survivors += 1

    print(f"\nFinal stage       : {final_csv.name}")
    print(f"Total survivors   : {total_survivors}")
    print()
    print(f"  {'Plate':<10}  {'Tiles':>6}  {'Survivors':>10}  {'Surv/tile':>10}")
    print(f"  {'-'*10}  {'-'*6}  {'-'*10}  {'-'*10}")
    for plate in sorted(plate_survivors, key=lambda p: -plate_survivors[p]):
        n_surv = plate_survivors[plate]
        n_t = len(plate_tiles_seen[plate])
        ratio = n_surv / n_t if n_t else 0.0
        print(f"  {plate:<10}  {n_t:>6}  {n_surv:>10}  {ratio:>10.2f}")


if __name__ == "__main__":
    main()
