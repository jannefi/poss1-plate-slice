#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_run_stage_csvs.py

Create run-scoped CSV artifacts for Post-pipeline "shrinking set" fetchers.

Outputs (under ./work/runs/run-<date>/ by default):
 - source_extractor_final_filtered.csv       (master S0, canonical schema + annotations)
 - source_extractor_final_filtered__dedup.csv (derived; astronomical dedup across tiles)
 - stage_S0.csv                              (FINAL stage CSV; driven from dedup set)
 - upload_positional.csv                     (FINAL: src_id,ra,dec) + chunked variants
 - tile_manifest.csv                         (per-tile accounting; includes delta skips)

Optional debug/audit:
 - stage_S0__raw.csv
 - upload_positional__raw.csv (+ chunks)

Contract (canonical schema):
 src_id    = tile_id + ":" + object_id
 tile_id   = tile folder name (tile_RA..._DEC[pm]...)
 object_id = internal NUMBER (never emit raw NUMBER in upload CSVs)
 ra/dec    = prefer WCS-fixed coords when present; else fallbacks
 plate_id  = from tile_plan.csv (via --plate-map-csv)

Delta policy:
 Default: skip tiles that already have post1.status==ok in tile_status.json.
 Use --full to reprocess all tiles regardless.
 post1.status is written per tile on completion (including empty tiles).

Dedup policy (science-grade; WCSFIX-ready):
 Duplicates defined globally (no plate_id partition) by true angular separation:
   sep_arcsec(ra,dec) <= dedup_tol_arcsec
 Robust spatial hashing in unit-sphere XYZ coordinates.
 Deterministic representative selection per cluster: tie-break by src_id.
 (Previously partitioned by plate_id before comparing; that silently missed
 true duplicates whose two detections carried different plan-assigned
 plate_id labels near plate junctions/overlaps -- fixed 2026-08-02.)
"""

import argparse
import csv
import datetime as _dt
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


# ----------------------------
# Tile discovery (flat: ./data/tiles/<tile-id>/)
# ----------------------------
_PATTERN = "tile_RA*_DEC*"


def iter_tile_dirs(tiles_root: Path) -> Iterable[Path]:
    """Yield tile dirs under tiles_root (flat layout only)."""
    tiles_root = Path(tiles_root)
    if tiles_root.exists():
        for p in sorted(tiles_root.glob(_PATTERN)):
            if p.is_dir():
                yield p


# ----------------------------
# Plate map (tile_id -> plate_id / REGION)
# ----------------------------
def load_plate_map(csv_path: Path) -> Dict[str, str]:
    """
    Expect columns: tile_id + one of (plate_id, irsa_region, REGION, region).
    Returns dict: tile_id -> plate_id.
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return {}

    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        r = csv.reader(f)
        hdr = next(r, [])
        cols = [c.strip() for c in hdr]

    region_col = None
    for cand in ("plate_id", "irsa_region", "REGION", "region"):
        if cand in cols:
            region_col = cand
            break
    if region_col is None or "tile_id" not in cols:
        return {}

    idx_tile = cols.index("tile_id")
    idx_reg = cols.index(region_col)

    out: Dict[str, str] = {}
    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        r = csv.reader(f)
        _ = next(r, None)
        for row in r:
            if not row or len(row) <= max(idx_tile, idx_reg):
                continue
            tid = str(row[idx_tile]).strip()
            reg = str(row[idx_reg]).strip()
            if tid:
                out[tid] = reg
    return out


# ----------------------------
# tile_status.json helpers
# ----------------------------
def _read_tile_steps(tile_dir: Path) -> dict:
    """Return the 'steps' dict from tile_status.json, or {} on any error."""
    try:
        p = tile_dir / "tile_status.json"
        if p.exists():
            return json.loads(p.read_text(encoding="utf-8")).get("steps", {})
    except Exception:
        pass
    return {}


def _is_post1_done(tile_dir: Path) -> bool:
    """Delta gate: only 'ok' is considered done (failed tiles are NOT skipped)."""
    return _read_tile_steps(tile_dir).get("post1", {}).get("status") == "ok"


def _set_post1_status(tile_dir: Path, status: str, reason: str = "") -> None:
    """Merge post1 status into tile_status.json (atomic write).

    status: 'ok' or 'failed'
    reason: optional human-readable reason (persisted for auditability)

    Note: delta mode only skips when status == 'ok'.
    """
    status = str(status).strip().lower()
    if status not in {"ok", "failed"}:
        status = "failed"

    p = tile_dir / "tile_status.json"
    try:
        data = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    except Exception:
        data = {}

    post1 = {"status": status}
    if reason:
        post1["reason"] = str(reason)

    data.setdefault("steps", {})["post1"] = post1

    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(data, indent=2), encoding="utf-8")
    tmp.replace(p)

# ----------------------------
# RA/Dec column picking
# ----------------------------
_RA_CANDS = ["RA_corr", "RA_CORR", "ALPHAWIN_J2000", "ALPHA_J2000", "RA", "X_WORLD", "RAJ2000"]
_DEC_CANDS = ["Dec_corr", "DEC_corr", "DEC_CORR", "DELTAWIN_J2000", "DELTA_J2000", "DEC", "Y_WORLD", "DEJ2000"]


def detect_header_cols(csv_path: Path) -> List[str]:
    with csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        r = csv.reader(f)
        return [c.strip() for c in next(r, [])]


def pick_radec_cols(cols: List[str]) -> Optional[Tuple[str, str]]:
    colset = set(cols)
    ra = next((c for c in _RA_CANDS if c in colset), None)
    dec = next((c for c in _DEC_CANDS if c in colset), None)
    if ra and dec:
        return ra, dec
    for a, b in [("ra", "dec"), ("RA_ICRS", "DE_ICRS")]:
        if a in colset and b in colset:
            return a, b
    return None


def pick_object_id_col(cols: List[str]) -> Optional[str]:
    for cand in ("NUMBER", "number", "object_id", "objectnumber", "objID"):
        if cand in cols:
            return cand
    return None


# ----------------------------
# Chunk writer
# ----------------------------
def write_chunks(rows: List[dict], out_path: Path, fieldnames: List[str], chunk_size: int, chunk_prefix: str):
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    chunks = []
    if len(rows) <= chunk_size:
        return chunks

    for i in range(0, len(rows), chunk_size):
        chunk = rows[i:i + chunk_size]
        idx = (i // chunk_size) + 1
        chunk_path = out_path.with_name(f"{chunk_prefix}_{idx:07d}.csv")
        with chunk_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(chunk)
        chunks.append(chunk_path)
    return chunks


# ----------------------------
# Dedup helpers
# ----------------------------
def angsep_arcsec(ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float) -> float:
    """Great-circle separation in arcsec (haversine; stable for tiny angles)."""
    ra1 = math.radians(ra1_deg % 360.0)
    ra2 = math.radians(ra2_deg % 360.0)
    dec1 = math.radians(dec1_deg)
    dec2 = math.radians(dec2_deg)
    dra = ra2 - ra1
    ddec = dec2 - dec1
    s1 = math.sin(ddec / 2.0)
    s2 = math.sin(dra / 2.0)
    a = s1 * s1 + math.cos(dec1) * math.cos(dec2) * s2 * s2
    a = min(1.0, max(0.0, a))
    c = 2.0 * math.asin(math.sqrt(a))
    return math.degrees(c) * 3600.0


def radec_to_unit_xyz(ra_deg: float, dec_deg: float) -> Tuple[float, float, float]:
    ra = math.radians(ra_deg % 360.0)
    dec = math.radians(dec_deg)
    cosd = math.cos(dec)
    return cosd * math.cos(ra), cosd * math.sin(ra), math.sin(dec)


def tol_arcsec_to_chord(tol_arcsec: float) -> float:
    tol_rad = (tol_arcsec / 3600.0) * (math.pi / 180.0)
    return 2.0 * math.sin(tol_rad / 2.0)


class UnionFind:
    def __init__(self, n: int):
        self.p = list(range(n))
        self.rank = [0] * n

    def find(self, a: int) -> int:
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a: int, b: int):
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.p[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.p[rb] = ra
        else:
            self.p[rb] = ra
            self.rank[ra] += 1


def dedup_rows_by_plate_radius_xyz(rows: List[dict], tol_arcsec: float) -> Tuple[List[dict], int]:
    """
    Science-grade global dedup by true angular separation <= tol_arcsec.
    Robust neighbor search uses XYZ unit-sphere binning (3D spatial hash).
    Representative selection: deterministic tie-break by src_id (lexicographic).
    Output ordering follows original input order of chosen representatives.

    Previously partitioned by plate_id before comparing, which silently missed
    true duplicates whose two detections carried different plan-assigned
    plate_id labels (common near naive-tessellation plate junctions/overlaps,
    e.g. two nearby real plates' fixed-size grids landing near-identical tiles).
    Dedup is now global -- spatial hash binning already limits comparisons to
    geometrically nearby points regardless of total row count, so removing the
    partition adds negligible cost (~1s at 255,921 rows) while fixing the miss.
    """
    if not rows:
        return rows, 0
    tol_arcsec = float(tol_arcsec)
    if tol_arcsec <= 0:
        return rows, 0

    cell = tol_arcsec_to_chord(tol_arcsec)
    if cell <= 0:
        return rows, 0

    xyz = [radec_to_unit_xyz(float(r["ra"]), float(r["dec"])) for r in rows]

    def key(x: float, y: float, z: float) -> Tuple[int, int, int]:
        return int(x / cell), int(y / cell), int(z / cell)

    bins: Dict[Tuple[int, int, int], List[int]] = {}
    uf = UnionFind(len(rows))

    for i, r in enumerate(rows):
        x, y, z = xyz[i]
        ix, iy, iz = key(x, y, z)

        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    cand = bins.get((ix + dx, iy + dy, iz + dz))
                    if not cand:
                        continue
                    for j in cand:
                        rj = rows[j]
                        if angsep_arcsec(float(r["ra"]), float(r["dec"]),
                                         float(rj["ra"]), float(rj["dec"])) <= tol_arcsec:
                            uf.union(i, j)

        bins.setdefault((ix, iy, iz), []).append(i)

    comps: Dict[int, List[int]] = {}
    for i in range(len(rows)):
        root = uf.find(i)
        comps.setdefault(root, []).append(i)

    chosen_indices: Set[int] = set()
    for members in comps.values():
        if len(members) == 1:
            chosen_indices.add(members[0])
            continue

        def rep_key(i: int) -> str:
            return str(rows[i].get("src_id") or "")

        rep_i = min(members, key=rep_key)
        chosen_indices.add(rep_i)

    out = [rows[i] for i in range(len(rows)) if i in chosen_indices]
    return out, len(rows) - len(out)


# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(description="Build run-scoped stage CSVs for shrinking-set fetchers.")
    ap.add_argument("--tiles-root", default="./data/tiles",
                    help="Tiles root directory (flat layout: ./data/tiles/<tile-id>/).")
    ap.add_argument("--plate-map-csv", default="./plans/tiles_poss1e_ps1.csv",
                    help="Plan CSV with tile_id + plate_id columns (from tessellate_plates).")
    ap.add_argument("--run-root", default="./work/runs",
                    help="Root for run output folders.")
    ap.add_argument("--run-tag", default="",
                    help="Optional run tag. Default: timestamp (run-YYYYMMDD_HHMMSS).")
    ap.add_argument("--chunk-size", type=int, default=2000,
                    help="Chunk size for upload/stage files.")
    ap.add_argument("--catalog-name", default="catalogs/sextractor_pass2.filtered.csv",
                    help="Relative path under tile dir to read survivors from.")
    ap.add_argument("--accept-empty-file-as-valid-empty", action="store_true", help="Treat a 0-byte survivors CSV as an explicitly valid-empty tile (mark post1 ok). Default is conservative: 0-byte file is treated as failure.")
    ap.add_argument("--full", action="store_true",
                    help="Reprocess all tiles regardless of post1 status in tile_status.json. "
                         "Default (delta mode): skip tiles already marked post1.status=ok.")

    # Dedup controls
    ap.add_argument("--dedup-tol-arcsec", type=float, default=3.0,
                    help="Astronomical dedup tolerance in arcsec (default 3.0). "
                         "The historical 0.25 was correct only UNDER WCS-FIXED "
                         "COORDINATES (context/02_DECISIONS.md), but WCS-fix is a "
                         "deliberately active deviation (REPRO_DEVIATIONS item 6) "
                         "so the archive's coordinates are unaligned and each tile "
                         "carries its own raw plate WCS. Measured on the full-sky "
                         "S0 (tools/dedup_radius_sweep.py): cross-tile duplicate "
                         "pairs plateau by 3\" (16,020 of the 16,058 found by 6\"), "
                         "while intra-tile pairs -- the signature of merging "
                         "genuinely distinct sources -- stay at exactly 0 until 8\". "
                         "3.0 sits at the knee with ~2.7x margin. Restore 0.25 only "
                         "if WCS-fix is re-enabled.")
    ap.add_argument("--no-dedup", dest="dedup_enable", action="store_false", default=True,
                    help="Disable astronomical dedup (not recommended).")
    ap.add_argument("--dedup-round-digits", type=int, default=6,
                    help="(deprecated/ignored) old rounding-based dedup parameter.")

    ap.add_argument("--write-raw-stage-and-uploads", action="store_true",
                    help="Also write stage/uploads for the raw (non-dedup) set.")

    args = ap.parse_args()

    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_tag = args.run_tag.strip() or f"run-{ts}"
    run_dir = Path(args.run_root) / run_tag
    run_dir.mkdir(parents=True, exist_ok=True)

    plate_map = load_plate_map(Path(args.plate_map_csv))

    tiles = list(iter_tile_dirs(Path(args.tiles_root)))
    seen_tile_ids: Set[str] = set()
    uniq_tiles: List[Path] = []
    for td in tiles:
        if td.name.startswith("tile_RA") and td.name not in seen_tile_ids:
            seen_tile_ids.add(td.name)
            uniq_tiles.append(td)

    manifest_rows: List[dict] = []
    out_rows: List[dict] = []
    delta_skipped = 0

    for td in uniq_tiles:
        tile_id = td.name
        plate_id = plate_map.get(tile_id, "")

        # Delta check
        if not args.full and _is_post1_done(td):
            manifest_rows.append({
                "tile_id": tile_id,
                "tile_path": str(td),
                "plate_id": plate_id,
                "rows_in_tile_filtered_csv": "",
                "rows_emitted_to_S0": "",
                "skipped_delta": 1,
                "notes": "delta skip: post1 already done",
            })
            delta_skipped += 1
            continue

        cat_path = td / args.catalog_name
        empty_file_valid = False
        n_in = 0
        n_out = 0
        note = ""

        if not cat_path.exists():
            note = "missing survivors csv"  # set later if needed
        elif cat_path.stat().st_size == 0:
            # A zero-byte file is the pipeline's deliberate "no survivors"
            # artifact: cli_pipeline writes it both at the morphology early-exit
            # and as veto carry-forward. It is a valid empty result, not a lost
            # tile, and must stay distinguishable from an absent file -- that is
            # the only signal a genuinely failed tile leaves behind.
            empty_file_valid = True
        else:
            cols = detect_header_cols(cat_path)
            radec = pick_radec_cols(cols)
            objcol = pick_object_id_col(cols)
            if not radec:
                note = "missing RA/Dec columns"
            elif not objcol:
                note = "missing object id column (NUMBER)"
            else:
                ra_col, dec_col = radec
                with cat_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
                    for row in csv.DictReader(f):
                        n_in += 1
                        try:
                            object_id = int(float(row.get(objcol, "")))
                        except Exception:
                            continue
                        try:
                            ra = float(row.get(ra_col, "nan"))
                            dec = float(row.get(dec_col, "nan"))
                        except Exception:
                            continue
                        out_rows.append({
                            "src_id": f"{tile_id}:{object_id}",
                            "tile_id": tile_id,
                            "object_id": object_id,
                            "ra": ra,
                            "dec": dec,
                            "plate_id": plate_id,
                        })
                        n_out += 1

        # Advance delta state only on success (or explicitly valid-empty).
        # Failure cases must NOT be marked ok, otherwise later delta runs can silently skip tiles.
        if note:
            _set_post1_status(td, "failed", note)
        else:
            if empty_file_valid:
                _set_post1_status(td, "ok", "valid_empty_emptyfile")
            elif n_in == 0:
                _set_post1_status(td, "ok", "valid_empty")
            elif n_out > 0:
                _set_post1_status(td, "ok")
            else:
                _set_post1_status(td, "failed", "no_valid_rows_emitted")

        manifest_rows.append({
            "tile_id": tile_id,
            "tile_path": str(td),
            "plate_id": plate_id,
            "rows_in_tile_filtered_csv": n_in,
            "rows_emitted_to_S0": n_out,
            "skipped_delta": 0,
            "notes": note or ("valid empty (no survivors)" if empty_file_valid else ""),
        })

    # manifest
    manifest_path = run_dir / "tile_manifest.csv"
    mf_fields = ["tile_id", "tile_path", "plate_id",
                 "rows_in_tile_filtered_csv", "rows_emitted_to_S0", "skipped_delta", "notes"]
    with manifest_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=mf_fields)
        w.writeheader()
        w.writerows(manifest_rows)

    # dedup src_id within base set
    seen: Set[str] = set()
    base_rows: List[dict] = []
    dup_srcid_dropped = 0
    for r in out_rows:
        sid = r["src_id"]
        if sid in seen:
            dup_srcid_dropped += 1
            continue
        seen.add(sid)
        base_rows.append(r)

    # astronomical dedup
    dedup_rows = base_rows
    dedup_dropped = 0
    if args.dedup_enable:
        dedup_rows, dedup_dropped = dedup_rows_by_plate_radius_xyz(base_rows, args.dedup_tol_arcsec)

    final_rows = dedup_rows

    master_fields = ["src_id", "tile_id", "object_id", "ra", "dec", "plate_id"]

    # write master S0 (base)
    master_base_path = run_dir / "source_extractor_final_filtered.csv"
    with master_base_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=master_fields)
        w.writeheader()
        w.writerows(base_rows)

    # write dedup master
    master_dedup_path = run_dir / "source_extractor_final_filtered__dedup.csv"
    with master_dedup_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=master_fields)
        w.writeheader()
        w.writerows(dedup_rows)

    # FINAL stage/uploads from dedup set
    stage_fields = ["src_id", "tile_id", "object_id", "ra", "dec"]
    stage_path = run_dir / "stage_S0.csv"
    with stage_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=stage_fields)
        w.writeheader()
        for r in final_rows:
            w.writerow({k: r[k] for k in stage_fields})

    upload_pos_fields = ["src_id", "ra", "dec"]
    upload_pos_rows = [{k: r[k] for k in upload_pos_fields} for r in final_rows]
    upload_pos_path = run_dir / "upload_positional.csv"
    write_chunks(upload_pos_rows, upload_pos_path, upload_pos_fields, args.chunk_size, "upload_positional_chunk")

    if args.write_raw_stage_and_uploads:
        stage_raw_path = run_dir / "stage_S0__raw.csv"
        with stage_raw_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=stage_fields)
            w.writeheader()
            for r in base_rows:
                w.writerow({k: r[k] for k in stage_fields})

        upload_pos_raw_path = run_dir / "upload_positional__raw.csv"
        upload_pos_raw_rows = [{k: r[k] for k in upload_pos_fields} for r in base_rows]
        write_chunks(upload_pos_raw_rows, upload_pos_raw_path, upload_pos_fields,
                     args.chunk_size, "upload_positional__raw_chunk")

    # summary
    mode = "full" if args.full else "delta"
    summary = run_dir / "RUN_SUMMARY.txt"
    summary.write_text(
        "\n".join([
            f"run_dir: {run_dir}",
            f"mode: {mode}",
            f"tiles_scanned: {len(uniq_tiles)}",
            f"tiles_delta_skipped: {delta_skipped}",
            f"tiles_processed: {len(uniq_tiles) - delta_skipped}",
            f"S0_rows_raw: {len(out_rows)}",
            f"S0_rows_unique_src_id: {len(base_rows)}",
            f"S0_src_id_duplicates_dropped: {dup_srcid_dropped}",
            f"dedup_enabled: {bool(args.dedup_enable)}",
            f"dedup_method: global angular_sep <= {args.dedup_tol_arcsec:.3f}\" (no plate_id partition; XYZ spatial hash; deterministic rep: src_id)",
            f"dedup_rows: {len(dedup_rows)} (dropped={dedup_dropped})",
            f"final_rows_for_stage_and_uploads: {len(final_rows)} (dedup set)",
            f"plate_map_csv: {args.plate_map_csv}",
            f"master_csv: {master_base_path.name}",
            f"master_dedup_csv: {master_dedup_path.name}",
            f"stage_csv: {stage_path.name}",
            f"upload_positional: {upload_pos_path.name}",
            "upload_skybot: (disabled) not produced by this script",
        ]) + "\n",
        encoding="utf-8"
    )

    print(f"[OK] run folder: {run_dir}  mode={mode}")
    print(f"[OK] tiles: scanned={len(uniq_tiles)} processed={len(uniq_tiles) - delta_skipped} skipped(delta)={delta_skipped}")
    print(f"[OK] master(base): {master_base_path.name} rows={len(base_rows)} (dup_src_id_dropped={dup_srcid_dropped})")
    print(f"[OK] dedup: {master_dedup_path.name} rows={len(dedup_rows)} (dropped={dedup_dropped})")
    print(f"[OK] FINAL stage/uploads (dedup): rows={len(final_rows)}")
    print(f"[OK] manifest: {manifest_path.name} rows={len(manifest_rows)}")


if __name__ == "__main__":
    main()
