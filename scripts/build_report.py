#!/usr/bin/env python3
"""
build_report.py — cross-run funnel report for vasco60 post-pipeline stages.

Reads all run folders under --runs-dir, collects counts from RUN_SUMMARY.txt
and stage ledger JSON files, and produces:

  <out-dir>/report-<timestamp>/
    funnel.txt          — ASCII funnel table with rejection %
    funnel.json         — same data, machine-readable
    survivors.csv       — deduplicated survivors: src_id,tile_id,ra,dec,plate_id,obs_date,run_id
    tile_coverage.csv   — tiles seen: tile_id,plate_id,run_id,rows_to_S0
    report_index.txt    — one-page human summary

Usage:
    python scripts/build_report.py --runs-dir ./work/runs --out-dir ./work/reports
"""

import argparse
import csv
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# Stage pipeline order definition
# ---------------------------------------------------------------------------
# Each entry: (ledger_glob, display_label, section, short_label)
# section: "pipeline" | "postprocess"
STAGES = [
    # These are synthesised from RUN_SUMMARY.txt, not ledgers
    ("__pass2_raw__",  "pass2 raw detections",    "pipeline",    "RAW"),
    ("__S0__",         "S0  post-MNRAS+dedup",    "pipeline",    "S0"),
    # Post-process stages in order
    ("stage_S0M_MORPH_ledger.json",  "S0M MORPH",   "postprocess", "S0M"),
    ("stage_S0S_SHAPE_ledger.json",  "S4S SHAPE",   "postprocess", "S0S"),
    ("stage_S1_GSC_ledger.json",     "S1  GSC",     "postprocess", "S1"),
    ("stage_S2_SKYBOT_ledger.json",  "S2  SKYBOT",  "postprocess", "S2"),
    ("stage_S3_SCOS_ledger.json",    "S3  SCOS",    "postprocess", "S3"),
    ("stage_S4_PTF_ledger.json",     "S4  PTF",     "postprocess", "S4"),
    ("stage_S5_VSX_ledger.json",          "S5  VSX",        "postprocess", "S5"),
    ("stage_S6_SCOPE_DEC_ledger.json",    "S6  SCOPE_DEC",  "postprocess", "S6"),
]

# Checked in order; first found wins.  Runs without S6 fall back to S5.
FINAL_STAGE_CANDIDATES = ["stage_S6_SCOPE_DEC.csv", "stage_S5_VSX.csv"]


def _get_final_stage_csv(stages_dir: Path) -> Path | None:
    for name in FINAL_STAGE_CANDIDATES:
        p = stages_dir / name
        if p.exists():
            return p
    return None


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _compute_pass2_raw(manifest_path: Path, tiles_root: Path) -> int | None:
    """Compute total pass2 unfiltered detections by summing veto_start_rows from
    each active tile's MNRAS_SUMMARY.json — same logic as run_detection_totals.py.
    Returns None if manifest is missing; returns 0 if no MNRAS_SUMMARY.json found.
    """
    if not manifest_path.exists():
        return None
    with manifest_path.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    def _int(val):
        return int(val) if val and val.strip() else 0

    active = [r for r in rows if _int(r.get("skipped_delta", "0")) == 0]
    total = 0
    found = 0
    for row in active:
        summary_path = tiles_root / row["tile_id"] / "MNRAS_SUMMARY.json"
        if not summary_path.exists():
            continue
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            total += summary.get("veto_start_rows", 0)
            found += 1
        except Exception:
            continue
    return total if found > 0 else None


def _parse_run_summary(path: Path) -> dict:
    """Extract structured fields from RUN_SUMMARY.txt.
    Only reads lines up to and including the 'upload_skybot:' line.
    """
    result = {"pass2_raw": None, "S0_rows": None, "mode": None,
              "tiles_processed": None, "run_dir": None}
    if not path.exists():
        return result

    # Read the script-generated portion only; stop before the per-plate table
    # (user may have appended free text after that table).
    stop_pattern = re.compile(r"^\s+Plate\s+Tiles\s+Survivors")
    lines = []
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if stop_pattern.match(line):
                break
            lines.append(line)

    text = "".join(lines)

    # pass2 unfiltered
    m = re.search(r"pass2 unfiltered.*?:\s+([\d,]+)", text)
    if m:
        result["pass2_raw"] = int(m.group(1).replace(",", ""))

    # dedup_rows (S0 input after dedup)
    m = re.search(r"dedup_rows:\s+(\d+)", text)
    if m:
        result["S0_rows"] = int(m.group(1))
    else:
        # fallback: final_rows_for_stage_and_uploads
        m = re.search(r"final_rows_for_stage_and_uploads:\s+(\d+)", text)
        if m:
            result["S0_rows"] = int(m.group(1))

    m = re.search(r"^mode:\s+(\S+)", text, re.MULTILINE)
    if m:
        result["mode"] = m.group(1)

    m = re.search(r"^tiles_processed:\s+(\d+)", text, re.MULTILINE)
    if m:
        result["tiles_processed"] = int(m.group(1))

    m = re.search(r"^run_dir:\s+(\S+)", text, re.MULTILINE)
    if m:
        result["run_dir"] = m.group(1)

    return result


def _load_ledger(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _kept(ledger: dict | None) -> int | None:
    if ledger is None:
        return None
    t = ledger.get("totals", {})
    for key in ("kept_rows", "kept"):
        if key in t:
            return int(t[key])
    return None


def _input_rows(ledger: dict | None) -> int | None:
    if ledger is None:
        return None
    t = ledger.get("totals", {})
    for key in ("input_rows", "total_rows"):
        if key in t:
            return int(t[key])
    return None


def _run_id(run_dir: Path) -> str:
    """Extract e.g. 'R1' from 'run-R1-20260327_165043'."""
    m = re.search(r"run-([^-]+)-", run_dir.name)
    return m.group(1) if m else run_dir.name


def _pct(new_val: int | None, old_val: int | None) -> str:
    if new_val is None or old_val is None or old_val == 0:
        return ""
    rej = old_val - new_val
    pct = 100.0 * rej / old_val
    # Avoid rounding artefact: only show -100.0% when truly all rows removed
    if new_val > 0 and round(pct, 1) >= 100.0:
        return "-99.9%"
    return f"-{pct:.1f}%"


# ---------------------------------------------------------------------------
# Main collection logic
# ---------------------------------------------------------------------------

def _load_tile_date_obs(tile_to_plate_csv: Path) -> dict:
    """Return {tile_id: tile_date_obs} from tile_to_plate.csv."""
    result = {}
    if not tile_to_plate_csv.exists():
        return result
    with tile_to_plate_csv.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            tid = row.get("tile_id", "").strip()
            obs = row.get("tile_date_obs", "").strip()
            if tid:
                result[tid] = obs
    return result


def collect_run(run_dir: Path, tile_date_obs_map: dict | None = None,
                tiles_root: Path | None = None) -> dict:
    """Return a dict with all counts for one run folder."""
    rid = _run_id(run_dir)
    summary = _parse_run_summary(run_dir / "RUN_SUMMARY.txt")
    if tiles_root is None:
        tiles_root = Path("data/tiles")

    # Load all ledgers
    ledgers = {}
    stages_dir = run_dir / "stages"
    for ledger_glob, _, _, *_ in STAGES:
        if ledger_glob.startswith("__"):
            continue
        p = stages_dir / ledger_glob
        ledgers[ledger_glob] = _load_ledger(p)

    # Build counts list in STAGES order
    counts = []
    for ledger_glob, label, section, *_ in STAGES:
        if ledger_glob == "__pass2_raw__":
            counts.append(_compute_pass2_raw(run_dir / "tile_manifest.csv", tiles_root))
        elif ledger_glob == "__S0__":
            counts.append(summary["S0_rows"])
        else:
            counts.append(_kept(ledgers.get(ledger_glob)))

    # Load tile_manifest for plate coverage
    tiles = {}
    manifest = run_dir / "tile_manifest.csv"
    if manifest.exists():
        with manifest.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("rows_emitted_to_S0", "0") not in ("0", ""):
                    tiles[row["tile_id"]] = {
                        "plate_id": row.get("plate_id", ""),
                        "rows_to_S0": int(row.get("rows_emitted_to_S0", 0)),
                        "run_id": rid,
                    }

    # Load final survivors
    survivors = []
    final_csv = _get_final_stage_csv(stages_dir)
    if final_csv is not None:
        # Build plate_id lookup from manifest
        plate_map = {t: d["plate_id"] for t, d in tiles.items()}
        with final_csv.open(encoding="utf-8") as f:
            for row in csv.DictReader(f):
                src_id = row.get("src_id", "")
                tile_id = src_id.split(":")[0] if ":" in src_id else row.get("tile_id", "")
                survivors.append({
                    "src_id":   src_id,
                    "tile_id":  tile_id,
                    "ra":       row.get("ra", ""),
                    "dec":      row.get("dec", ""),
                    "plate_id": plate_map.get(tile_id, ""),
                    "obs_date": (tile_date_obs_map or {}).get(tile_id, ""),
                    "run_id":   rid,
                })

    return {
        "run_id":          rid,
        "run_dir":         str(run_dir),
        "mode":            summary.get("mode", ""),
        "tiles_processed": summary.get("tiles_processed"),
        "counts":          counts,   # parallel to STAGES list
        "tiles":           tiles,
        "survivors":       survivors,
    }


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _fmt_count(v: int | None) -> str:
    if v is None:
        return "-"
    if v >= 1_000_000:
        return f"{v/1_000_000:.2f}M"
    if v >= 1_000:
        return f"{v:,}"
    return str(v)


def build_funnel_text(runs: list[dict]) -> str:
    """Transposed funnel: one row per run, one column per stage."""
    n_stages = len(STAGES)
    short_labels = [s[3] for s in STAGES]

    # Pre-compute cells[run_idx][stage_idx] = display string
    # Also compute a TOTAL row (sum across runs per stage)
    all_cells: list[list[str]] = []
    total_counts: list[int | None] = [None] * n_stages
    total_denoms: list[int | None] = [None] * n_stages

    for r in runs:
        run_cells: list[str] = []
        prev_count: int | None = None
        for s_idx, (ledger_glob, _label, _section, _short) in enumerate(STAGES):
            count = r["counts"][s_idx]

            if ledger_glob.startswith("__"):
                denom = prev_count
            else:
                ledger_path = Path(r["run_dir"]) / "stages" / ledger_glob
                led = _load_ledger(ledger_path)
                ledger_input = _input_rows(led)
                denom = ledger_input if ledger_input is not None else prev_count

            pct = _pct(count, denom)
            cell = _fmt_count(count)
            if pct:
                cell = f"{cell} ({pct})"
            run_cells.append(cell)

            # Accumulate totals
            if count is not None:
                total_counts[s_idx] = (total_counts[s_idx] or 0) + count
            if denom is not None:
                total_denoms[s_idx] = (total_denoms[s_idx] or 0) + denom

            prev_count = count
        all_cells.append(run_cells)

    # Build TOTAL row
    total_cells: list[str] = []
    for s_idx in range(n_stages):
        tc = total_counts[s_idx]
        pct = _pct(tc, total_denoms[s_idx]) if s_idx > 0 else ""
        cell = _fmt_count(tc)
        if pct:
            cell = f"{cell} ({pct})"
        total_cells.append(cell)

    # Column widths: max of short label, all run cells, total cell
    run_id_w = max(len(r["run_id"]) for r in runs)
    run_id_w = max(run_id_w, 3)  # at least "Run"
    col_widths = [
        max(len(short_labels[s]),
            max(len(all_cells[r][s]) for r in range(len(runs))),
            len(total_cells[s]))
        for s in range(n_stages)
    ]

    # Section separator positions: insert "|" between pipeline and postprocess
    def _row(label: str, cells: list[str]) -> str:
        row = f"{label:<{run_id_w}}"
        last_section = None
        for s_idx, (_, _lbl, section, _short) in enumerate(STAGES):
            if last_section == "pipeline" and section == "postprocess":
                row += "  |"
            row += f"  {cells[s_idx]:>{col_widths[s_idx]}}"
            last_section = section
        return row

    # Header row
    header_cells = short_labels[:]
    header = f"{'Run':<{run_id_w}}"
    last_section = None
    for s_idx, (_, _lbl, section, short) in enumerate(STAGES):
        if last_section == "pipeline" and section == "postprocess":
            header += "  |"
        header += f"  {short:>{col_widths[s_idx]}}"
        last_section = section

    sep = "-" * len(_row("", ["—" * col_widths[s] for s in range(n_stages)]))

    lines = [
        "VASCO60 — Detection Funnel",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f"Runs: {', '.join(r['run_id'] for r in runs)}",
        f"Columns: RAW S0 = pipeline  |  S0M..S6 = post-pipeline stages",
        "",
        header,
        sep,
    ]
    for r_idx, r in enumerate(runs):
        lines.append(_row(r["run_id"], all_cells[r_idx]))
    lines.append(sep)
    lines.append(_row("TOTAL", total_cells))
    lines.append("")
    return "\n".join(lines)


def build_funnel_json(runs: list[dict]) -> dict:
    rows = []
    for i, (ledger_glob, label, section, *_) in enumerate(STAGES):
        entry = {"stage": label, "section": section, "runs": {}}
        for r in runs:
            entry["runs"][r["run_id"]] = r["counts"][i]
        rows.append(entry)
    return {"generated": datetime.now(timezone.utc).isoformat(), "funnel": rows}


def build_survivors(runs: list[dict]) -> list[dict]:
    seen = {}
    result = []
    for r in runs:
        for s in r["survivors"]:
            sid = s["src_id"]
            if sid not in seen:
                seen[sid] = True
                result.append(s)
    return result


def build_tile_coverage(runs: list[dict]) -> list[dict]:
    rows = []
    for r in runs:
        for tile_id, td in sorted(r["tiles"].items()):
            rows.append({
                "tile_id":    tile_id,
                "plate_id":   td["plate_id"],
                "run_id":     r["run_id"],
                "rows_to_S0": td["rows_to_S0"],
            })
    return rows


def build_index(runs: list[dict], survivors: list[dict],
                tile_coverage: list[dict], report_dir: Path) -> str:
    total_tiles = len(tile_coverage)
    total_surv = len(survivors)

    plates = sorted({s["plate_id"] for s in survivors if s["plate_id"]})
    all_plates_tiles = sorted({r["plate_id"] for r in tile_coverage if r["plate_id"]})

    cross_run_dupes = sum(
        len(r["survivors"]) for r in runs
    ) - total_surv

    lines = [
        "VASCO60 — Report Index",
        f"Generated : {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        f"Report dir: {report_dir}",
        "",
        f"Runs processed : {len(runs)}  ({', '.join(r['run_id'] for r in runs)})",
        "Run modes      : " + ", ".join(r["run_id"] + "=" + r["mode"] for r in runs),
        "",
        f"Tiles in coverage : {total_tiles}",
        f"Plates in coverage: {len(all_plates_tiles)}  ({', '.join(all_plates_tiles)})",
        "",
        f"Final survivors (deduplicated) : {total_surv}",
        f"Cross-run duplicates dropped   : {cross_run_dupes}",
        f"Plates with survivors          : {len(plates)}  ({', '.join(plates)})",
        "",
        "Per-run survivor counts:",
    ]
    for r in runs:
        lines.append(f"  {r['run_id']:6s}  {len(r['survivors'])} survivors")

    lines += [
        "",
        "Output files:",
        "  funnel.txt        — stage-by-stage rejection funnel",
        "  funnel.json       — machine-readable funnel data",
        "  survivors.csv     — deduplicated final candidates with plate_id",
        "  tile_coverage.csv — all tiles processed, by run",
        "  report_index.txt  — this file",
    ]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Build cross-run report for vasco60.")
    parser.add_argument("--runs-dir", default="./work/runs",
                        help="Directory containing run-R* folders")
    parser.add_argument("--out-dir", default="./work/reports",
                        help="Parent directory for report output")
    parser.add_argument("--run-glob", default="run-R*",
                        help="Glob to match run folder names (default: run-R*)")
    parser.add_argument("--tile-to-plate", default="./data/metadata/tile_to_plate.csv",
                        help="tile_to_plate.csv path (default: ./data/metadata/tile_to_plate.csv)")
    parser.add_argument("--tiles-root", default="./data/tiles",
                        help="Root of tile folders containing MNRAS_SUMMARY.json (default: ./data/tiles)")
    args = parser.parse_args()

    runs_dir = Path(args.runs_dir)
    out_parent = Path(args.out_dir)

    def _run_sort_key(p: Path) -> tuple:
        m = re.search(r"run-([A-Za-z]*)(\d+)-", p.name)
        return (m.group(1), int(m.group(2))) if m else (p.name, 0)

    run_dirs = sorted(runs_dir.glob(args.run_glob), key=_run_sort_key)
    if not run_dirs:
        print(f"[REPORT] ERROR: no run folders found matching {runs_dir}/{args.run_glob}",
              file=sys.stderr)
        sys.exit(1)

    print(f"[REPORT] Found {len(run_dirs)} run folder(s): "
          f"{', '.join(d.name for d in run_dirs)}")

    tile_date_obs_map = _load_tile_date_obs(Path(args.tile_to_plate))
    if tile_date_obs_map:
        print(f"[REPORT] Loaded obs dates for {len(tile_date_obs_map)} tiles from {args.tile_to_plate}")
    else:
        print(f"[REPORT] WARNING: tile_to_plate.csv not found or empty at {args.tile_to_plate} — obs_date will be blank")

    runs = []
    for rd in run_dirs:
        print(f"[REPORT]   collecting {rd.name} ...", end=" ", flush=True)
        r = collect_run(rd, tile_date_obs_map, tiles_root=Path(args.tiles_root))
        print(f"{len(r['survivors'])} survivors, {len(r['tiles'])} tiles")
        runs.append(r)

    # --- Build outputs ---
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    report_dir = out_parent / f"report-{ts}"
    report_dir.mkdir(parents=True, exist_ok=True)

    funnel_txt = build_funnel_text(runs)
    funnel_json = build_funnel_json(runs)
    survivors = build_survivors(runs)
    tile_coverage = build_tile_coverage(runs)
    index_txt = build_index(runs, survivors, tile_coverage, report_dir)

    # funnel.txt
    p = report_dir / "funnel.txt"
    p.write_text(funnel_txt, encoding="utf-8")
    print(f"[REPORT] wrote: {p}")

    # funnel.json
    p = report_dir / "funnel.json"
    p.write_text(json.dumps(funnel_json, indent=2), encoding="utf-8")
    print(f"[REPORT] wrote: {p}")

    # survivors.csv
    p = report_dir / "survivors.csv"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["src_id", "tile_id", "ra", "dec", "plate_id", "obs_date", "run_id"])
        w.writeheader()
        w.writerows(survivors)
    print(f"[REPORT] wrote: {p}  ({len(survivors)} rows)")

    # tile_coverage.csv
    p = report_dir / "tile_coverage.csv"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["tile_id", "plate_id", "run_id", "rows_to_S0"])
        w.writeheader()
        w.writerows(tile_coverage)
    print(f"[REPORT] wrote: {p}  ({len(tile_coverage)} rows)")

    # report_index.txt
    p = report_dir / "report_index.txt"
    p.write_text(index_txt, encoding="utf-8")
    print(f"[REPORT] wrote: {p}")

    print(f"\n[REPORT] Done — {len(survivors)} deduplicated survivors across "
          f"{len(runs)} runs.")
    print(f"[REPORT] Report dir: {report_dir}")


if __name__ == "__main__":
    main()
