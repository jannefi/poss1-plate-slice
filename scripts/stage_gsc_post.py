#!/usr/bin/env python3
"""
EXPERIMENTAL — Post-pipeline GSC 2.4.2 reduction stage.

Status
------
Not an official veto stage. Use for exploration and candidate reduction only.
Results are not yet validated against the full pipeline. Do not use as a
hard gate without further testing.

Goal
----
Given a run directory containing a stage CSV (e.g. stage_S4_VSX.csv), cross-match
against the Guide Star Catalog 2.4.2 (VizieR table I/353/gsc242) within a
configurable radius (default 5 arcsec) and write a shrinking-set stage output.

In testing against vasco60 datasets, this stage has produced 10-20% reductions.
The GSC 2.4.2 covers the full sky and includes multi-epoch optical photometry,
making it a good complement to the existing veto chain.

Why GSC is a good veto candidate for POSS-I red detections
----------------------------------------------------------
In GSC-II (GSC 2.3), export rules require objects to be detected on at least one
IIIaF/IIIaJ/IV-N plate. Single-plate defects are discarded and unmatched artefacts
are excluded. The published classification is the mode across plates, ignoring
artifact/defect classifications for matched objects under the assumption that
objects matched on multiple plates are real.

Reference: https://iopscience.iop.org/article/10.1088/0004-6256/136/2/735


Usage
-----
python scripts/stage_gsc_post.py \\
    --run-dir ./work/runs/run-S1-... \\
    --input-glob 'stages/stage_S4_VSX.csv' \\
    --stage S5

Outputs (under <run-dir>/stages/)
----------------------------------
1) stage_<STAGE>_GSC.csv
   Kept remainder AFTER GSC elimination (rows WITHOUT a GSC match).
   Columns: src_id, ra, dec

2) stage_<STAGE>_GSC_flags.csv
   Full audit table for ALL input rows.
   Columns: src_id, ra, dec, has_gsc_match, best_sep_arcsec, gsc2_id

3) stage_<STAGE>_GSC_ledger.json
   Counts + parameters used.

Notes
-----
- Requires STILTS in PATH.
- Uses STILTS cdsskymatch against VizieR (network required).
- Input must have src_id (or row_id), ra, dec columns.
- GSC 2.4.2 Class column: 0=star, 3=non-star, 2=blend; all classes are matched
  (no class filtering) since the goal is positional vetoing.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

DEFAULT_CDSTABLE = "I/353/gsc242"


def _ensure_tool(name: str) -> str:
    path = shutil.which(name)
    if path:
        return path
    raise RuntimeError(f"Required tool '{name}' not found in PATH")


def _read_header(path: Path) -> List[str]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.reader(f)
        return [c.strip().lstrip("\ufeff") for c in next(r, [])]


def _detect_cols(cols: List[str]) -> Tuple[str, str, str]:
    cset = {c.strip() for c in cols}

    if "src_id" in cset:
        src = "src_id"
    elif "row_id" in cset:
        src = "row_id"
    else:
        raise RuntimeError("Input CSV missing required id column 'src_id' (or 'row_id')")

    ra_candidates = ["ra", "RA", "RA_corr", "ALPHAWIN_J2000", "ALPHA_J2000"]
    dec_candidates = ["dec", "DEC", "Dec", "Dec_corr", "DELTAWIN_J2000", "DELTA_J2000"]

    ra = next((c for c in ra_candidates if c in cset), None)
    dec = next((c for c in dec_candidates if c in cset), None)
    if not ra or not dec:
        raise RuntimeError(
            "Input CSV missing RA/Dec columns (expected one of ra/dec, RA/DEC, RA_corr/Dec_corr, etc.)"
        )
    return src, ra, dec


def _count_data_rows(path: Path) -> int:
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            n = sum(1 for _ in f)
        return max(0, n - 1)
    except Exception:
        return 0


def _try_import_vasco_cdsskymatch():
    try:
        from vasco.utils.cdsskymatch import cdsskymatch  # type: ignore
        return cdsskymatch
    except Exception:
        return None


def _run_stilts_cdsskymatch(
    stilts_bin: str,
    in_csv: Path,
    out_csv: Path,
    *,
    ra_col: str,
    dec_col: str,
    cdstable: str,
    radius_arcsec: float,
    find: str = "best",
    blocksize: Optional[int] = None,
) -> None:
    cmd = [
        stilts_bin,
        "cdsskymatch",
        "ifmt=csv",
        f"in={str(in_csv)}",
        f"ra={ra_col}",
        f"dec={dec_col}",
        f"cdstable={cdstable}",
        f"radius={radius_arcsec}",
        f"find={find}",
        "ofmt=csv",
        f"out={str(out_csv)}",
    ]
    if blocksize is not None:
        cmd.append(f"blocksize={int(blocksize)}")

    subprocess.run(cmd, check=True, text=True, capture_output=True)


@dataclass
class ChunkStats:
    chunk: str
    input_rows: int
    matched_rows: int
    kept_rows: int


def _load_matches(match_csv: Path, src_col: str) -> Dict[str, Tuple[float, str]]:
    """Return mapping src_id -> (best_sep_arcsec, gsc2_id) from STILTS output."""
    if not match_csv.exists() or match_csv.stat().st_size == 0:
        return {}
    with match_csv.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.DictReader(f)
        if not r.fieldnames:
            return {}
        fields = set(r.fieldnames)
        dist_field = "angDist" if "angDist" in fields else None
        gsc2_field = "GSC2" if "GSC2" in fields else None

        out: Dict[str, Tuple[float, str]] = {}
        for row in r:
            sid = row.get(src_col, "") or row.get("src_id", "") or row.get("row_id", "")
            if not sid:
                continue
            sep = float("nan")
            if dist_field:
                try:
                    sep = float(row.get(dist_field, "nan"))
                except Exception:
                    pass
            gsc2_id = (row.get(gsc2_field, "") or "") if gsc2_field else ""
            out[sid] = (sep, gsc2_id)
    return out


def _normalize_ra_chunk(path: Path, ra_col: str, out_path: Path) -> None:
    """Rewrite a chunk CSV with ra_col wrapped into [0, 360).

    Some naive-tessellation tiles straddle the RA=0/360 wrap boundary, so
    their candidates' RA can come out slightly negative (e.g. -1.68 instead
    of 358.32) -- the exact same sky position, just outside the [0,360)
    range CDS's xmatch service requires ("field RA < 0 or > 360!", a hard
    server-side rejection of the whole batch, not just the offending rows).
    """
    with path.open(newline="", encoding="utf-8", errors="ignore") as f_in:
        r = csv.DictReader(f_in)
        fieldnames = r.fieldnames or []
        rows = list(r)
    for row in rows:
        try:
            row[ra_col] = f"{float(row[ra_col]) % 360.0:.8f}"
        except Exception:
            continue
    with out_path.open("w", newline="", encoding="utf-8") as f_out:
        w = csv.DictWriter(f_out, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def _iter_input_rows(path: Path, src_col: str, ra_col: str, dec_col: str) -> Iterable[Tuple[str, str, str]]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.DictReader(f)
        for row in r:
            sid = (row.get(src_col) or "").strip()
            ra = (row.get(ra_col) or "").strip()
            dec = (row.get(dec_col) or "").strip()
            if not sid or not ra or not dec:
                continue
            yield sid, ra, dec


def main() -> int:
    ap = argparse.ArgumentParser(
        description="[EXPERIMENTAL] Post-pipeline GSC 2.4.2 reduction stage."
    )
    ap.add_argument("--run-dir", required=True, help="Run folder, e.g. ./work/runs/run-S1-...")
    ap.add_argument(
        "--input-glob",
        default="stages/stage_S4_VSX.csv",
        help="Glob (relative to run-dir) for input stage CSV. Default: stages/stage_S4_VSX.csv",
    )
    ap.add_argument("--stage", default="S5", help="Stage label used in output filenames. Default: S5")
    ap.add_argument("--out-dir", default=None, help="Output directory. Default: <run-dir>/stages")
    ap.add_argument("--cdstable", default=DEFAULT_CDSTABLE, help=f"CDS table id. Default: {DEFAULT_CDSTABLE}")
    ap.add_argument("--radius-arcsec", type=float, default=5.0, help="Match radius in arcsec. Default: 5")
    ap.add_argument("--blocksize", type=int, default=1000, help="STILTS blocksize. Default: 1000")
    ap.add_argument("--keep-temp", action="store_true", help="Keep per-chunk match files (debug).")
    args = ap.parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.exists():
        raise SystemExit(f"run-dir not found: {run_dir}")

    out_dir = Path(args.out_dir) if args.out_dir else (run_dir / "stages")
    out_dir.mkdir(parents=True, exist_ok=True)

    chunks = sorted(run_dir.glob(args.input_glob))
    if not chunks:
        raise SystemExit(f"No inputs matched: {run_dir}/{args.input_glob}")

    hdr = _read_header(chunks[0])
    src_col, ra_col, dec_col = _detect_cols(hdr)

    cdsskymatch_wrapper = _try_import_vasco_cdsskymatch()
    stilts_bin = _ensure_tool("stilts")

    stage = args.stage
    out_kept   = out_dir / f"stage_{stage}_GSC.csv"
    out_flags  = out_dir / f"stage_{stage}_GSC_flags.csv"
    out_ledger = out_dir / f"stage_{stage}_GSC_ledger.json"

    flags_fields = ["src_id", "ra", "dec", "has_gsc_match", "best_sep_arcsec", "gsc2_id", "source_chunk"]

    total_in = total_match = total_kept = 0
    per_chunk: List[ChunkStats] = []

    with out_kept.open("w", newline="", encoding="utf-8") as f_kept, \
         out_flags.open("w", newline="", encoding="utf-8") as f_flags:

        kept_w  = csv.DictWriter(f_kept,  fieldnames=["src_id", "ra", "dec"])
        flags_w = csv.DictWriter(f_flags, fieldnames=flags_fields)
        kept_w.writeheader()
        flags_w.writeheader()

        with tempfile.TemporaryDirectory(prefix="gsc_stage_") as tdir:
            tdir_p = Path(tdir)

            for ch in chunks:
                in_rows = _count_data_rows(ch)
                if in_rows == 0:
                    per_chunk.append(ChunkStats(ch.name, 0, 0, 0))
                    continue

                match_out = tdir_p / f"{ch.stem}__gsc_xmatch.csv"
                ch_norm = tdir_p / f"{ch.stem}__ra_normalized.csv"
                _normalize_ra_chunk(ch, ra_col, ch_norm)

                if cdsskymatch_wrapper is not None:
                    cdsskymatch_wrapper(
                        str(ch_norm),
                        str(match_out),
                        ra=ra_col,
                        dec=dec_col,
                        cdstable=args.cdstable,
                        radius_arcsec=float(args.radius_arcsec),
                        find="best",
                        ofmt="csv",
                        omode="out",
                        blocksize=int(args.blocksize) if args.blocksize else None,
                    )
                else:
                    _run_stilts_cdsskymatch(
                        stilts_bin, ch_norm, match_out,
                        ra_col=ra_col,
                        dec_col=dec_col,
                        cdstable=args.cdstable,
                        radius_arcsec=float(args.radius_arcsec),
                        find="best",
                        blocksize=int(args.blocksize) if args.blocksize else None,
                    )

                matches = _load_matches(match_out, src_col)
                mcount = len(matches)
                kcount = 0

                for sid, ra, dec in _iter_input_rows(ch_norm, src_col, ra_col, dec_col):
                    has = 1 if sid in matches else 0
                    sep, gsc2_id = matches.get(sid, (None, ""))

                    flags_w.writerow({
                        "src_id": sid,
                        "ra": ra,
                        "dec": dec,
                        "has_gsc_match": has,
                        "best_sep_arcsec": (f"{sep:.6f}" if (sep is not None and sep == sep) else ""),
                        "gsc2_id": gsc2_id,
                        "source_chunk": ch.name,
                    })

                    if has == 0:
                        kept_w.writerow({"src_id": sid, "ra": ra, "dec": dec})
                        kcount += 1

                per_chunk.append(ChunkStats(ch.name, in_rows, mcount, kcount))
                total_in    += in_rows
                total_match += mcount
                total_kept  += kcount

                if args.keep_temp:
                    dbg_dir = out_dir / "tmp_gsc_matches"
                    dbg_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(match_out, dbg_dir / match_out.name)

    ledger = {
        "experimental": True,
        "run_dir": str(run_dir),
        "input_glob": args.input_glob,
        "stage": stage,
        "cdstable": args.cdstable,
        "radius_arcsec": float(args.radius_arcsec),
        "blocksize": int(args.blocksize) if args.blocksize else None,
        "input_chunks": [p.name for p in chunks],
        "totals": {"input_rows": total_in, "matched_rows": total_match, "kept_rows": total_kept},
        "per_chunk": [cs.__dict__ for cs in per_chunk],
        "outputs": {"kept_csv": str(out_kept), "flags_csv": str(out_flags), "ledger_json": str(out_ledger)},
        "columns_detected": {"src_id_col": src_col, "ra_col": ra_col, "dec_col": dec_col},
    }
    out_ledger.write_text(json.dumps(ledger, indent=2), encoding="utf-8")

    print(f"[GSC] [EXPERIMENTAL] chunks={len(chunks)} input_rows={total_in} matched={total_match} kept={total_kept}")
    print(f"[GSC] wrote: {out_kept}")
    print(f"[GSC] wrote: {out_flags}")
    print(f"[GSC] wrote: {out_ledger}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
