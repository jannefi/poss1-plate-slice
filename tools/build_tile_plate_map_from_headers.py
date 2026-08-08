#!/usr/bin/env python3
"""
Derive an authoritative tile_id -> plate_id map by reading each archived
tile's own FITS header sidecar (raw/*.fits.header.json -> header.REGION),
instead of trusting a metadata CSV or a tessellation plan.

Why this exists: two other sources of the same mapping are unreliable.

  1. data/metadata/tiles_registry.csv and tile_to_plate.csv are written by
     vasco.utils.tile_metadata during step1-download. Before the locking
     fix in that module, concurrent step1 fetches clobbered each other's
     rows, so any tiles root fetched in parallel can be missing most of
     its rows (observed: 3,350 of 16,951 tiles present).
  2. A naive-tessellation plan's plate_id is the plate the *plan assumed*
     covers the tile center, which is not always the plate the archive
     actually served (see the plate-selection notes in
     context/03_NEXT_ACTIONS.md).

The FITS header is the archive's own answer to "which plate is this?", so
it is the source of truth -- consistent with vasco.utils.tile_metadata's
documented rule that plate_id is frozen to header REGION.

Output schema is tile_id,plate_id, which both consumers already accept
without modification:
  - scripts/build_run_stage_csvs.py   --plate-map-csv
  - scripts/compute_tile_plate_edge_report.py --tile-plan-csv

No catalog-specific logic: this script takes a tiles root and knows
nothing about any candidate catalog. Safe to commit.

Usage:
  python3 tools/build_tile_plate_map_from_headers.py \
      --tiles-root data/tiles \
      --out-csv work/tile_plate_map.csv

  # optionally also repair that root's own metadata CSVs (single writer,
  # so no race), filling in rows for tiles present on disk but missing:
  python3 tools/build_tile_plate_map_from_headers.py \
      --tiles-root data/tiles \
      --out-csv work/tile_plate_map.csv \
      --repair-metadata-dir data/metadata
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from vasco.utils.tile_metadata import (  # noqa: E402
    _csv_write_lock,
    _read_csv_rows,
    _read_header_sidecar,
    _utc_now_iso,
    _write_csv_atomic,
)

# Same schemas as vasco.utils.tile_metadata's own updaters.
T2P_FIELDS = [
    'tile_id', 'plate_id', 'tile_region', 'tile_survey', 'tile_date_obs', 'tile_fits',
    'irsa_region', 'irsa_filename', 'irsa_survey', 'irsa_platelabel', 'irsa_plateid',
    'irsa_date_obs', 'irsa_center_sep_deg',
]
REGISTRY_FIELDS = [
    'tile_id', 'ra_deg', 'dec_deg', 'survey', 'size_arcmin', 'pixel_scale_arcsec',
    'status', 'downloaded_utc', 'source', 'plate_id', 'notes',
]


def iter_tile_dirs(tiles_root: Path):
    for td in sorted(tiles_root.iterdir()):
        if td.is_dir() and td.name.startswith("tile_RA"):
            yield td


def tile_fits(tile_dir: Path):
    """Return the tile's raw FITS path, or None if it was never fetched."""
    cands = sorted((tile_dir / "raw").glob("*.fits"))
    return cands[0] if cands else None


def read_tile_header(tile_dir: Path):
    """Return (header_dict, fits_path). Empty dict if unreadable/absent."""
    fits_path = tile_fits(tile_dir)
    if fits_path is None:
        return {}, None
    hdr, _ = _read_header_sidecar(fits_path)
    return hdr, fits_path


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--tiles-root", required=True,
                    help="Root containing per-tile dirs (tile_RA..._DEC[pm]...).")
    ap.add_argument("--out-csv", required=True,
                    help="Output tile_id,plate_id CSV.")
    ap.add_argument("--repair-metadata-dir", default=None,
                    help="If set, also fill in missing rows in that dir's "
                         "tile_to_plate.csv and tiles_registry.csv.")
    ap.add_argument("--unresolved-csv", default=None,
                    help="If set, write the tiles whose REGION could not be "
                         "resolved here (tile_id,reason) for auditing.")
    args = ap.parse_args()

    tiles_root = Path(args.tiles_root)
    if not tiles_root.is_dir():
        raise SystemExit(f"tiles root not found: {tiles_root}")

    resolved: list[tuple[str, str]] = []
    unresolved: list[tuple[str, str]] = []
    plates: set[str] = set()

    for td in iter_tile_dirs(tiles_root):
        hdr, fits_path = read_tile_header(td)
        if fits_path is None:
            unresolved.append((td.name, "no raw FITS"))
            continue
        region = (hdr.get("REGION") or "").strip()
        if not region:
            unresolved.append((td.name, "no REGION in header sidecar"))
            continue
        resolved.append((td.name, region))
        plates.add(region)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["tile_id", "plate_id"])
        w.writerows(resolved)

    print(f"[PLATEMAP] tiles scanned:  {len(resolved) + len(unresolved)}")
    print(f"[PLATEMAP] resolved:       {len(resolved)}")
    print(f"[PLATEMAP] unresolved:     {len(unresolved)}")
    print(f"[PLATEMAP] distinct plates:{len(plates)}")
    print(f"[PLATEMAP] wrote {out_csv}")

    if args.unresolved_csv:
        up = Path(args.unresolved_csv)
        up.parent.mkdir(parents=True, exist_ok=True)
        with up.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["tile_id", "reason"])
            w.writerows(unresolved)
        print(f"[PLATEMAP] wrote {up}")

    if not args.repair_metadata_dir:
        return

    # --- optional metadata repair -------------------------------------
    # One bulk read-modify-write per file, under the module's own lock.
    # Calling update_tiles_registry()/update_tile_to_plate_csv() per tile
    # would be O(n^2): each rewrites the entire CSV.
    meta_dir = Path(args.repair_metadata_dir)
    n_t2p = _bulk_repair(meta_dir / "tile_to_plate.csv", T2P_FIELDS, resolved,
                         tiles_root, _t2p_row)
    n_reg = _bulk_repair(meta_dir / "tiles_registry.csv", REGISTRY_FIELDS, resolved,
                         tiles_root, _registry_row)
    print(f"[PLATEMAP] backfilled tile_to_plate.csv rows:  {n_t2p}")
    print(f"[PLATEMAP] backfilled tiles_registry.csv rows: {n_reg}")


def _bulk_repair(out: Path, fieldnames, resolved, tiles_root: Path, make_row) -> int:
    """Add a row for every resolved tile missing from `out`. Existing rows win."""
    with _csv_write_lock(out):
        rows = _read_csv_rows(out)
        added = 0
        for tile_id, plate_id in resolved:
            if tile_id in rows:
                continue
            hdr, fits_path = read_tile_header(tiles_root / tile_id)
            rows[tile_id] = make_row(tile_id, plate_id, hdr, fits_path)
            added += 1
        if added:
            _write_csv_atomic(out, fieldnames, rows)
    return added


def _t2p_row(tile_id, plate_id, hdr, fits_path) -> dict:
    survey = (hdr.get("SURVEY") or "").strip()
    date_obs = (hdr.get("DATE-OBS") or "").strip()
    fits_name = fits_path.name if fits_path else ""
    return {
        "tile_id": tile_id, "plate_id": plate_id, "tile_region": plate_id,
        "tile_survey": survey, "tile_date_obs": date_obs, "tile_fits": fits_name,
        "irsa_region": plate_id, "irsa_filename": fits_name, "irsa_survey": survey,
        "irsa_platelabel": (hdr.get("PLTLABEL") or "").strip(),
        "irsa_plateid": (hdr.get("PLATEID") or "").strip(),
        "irsa_date_obs": date_obs, "irsa_center_sep_deg": "",
    }


def _registry_row(tile_id, plate_id, hdr, fits_path) -> dict:
    ra, dec = _center_from_header(hdr)
    nx = _as_float(hdr.get("NAXIS1"))
    scale = abs(_as_float(hdr.get("CD1_1")) or 0.0) * 3600.0
    return {
        "tile_id": tile_id,
        "ra_deg": f"{ra:.6f}", "dec_deg": f"{dec:.6f}",
        "survey": (hdr.get("SURVEY") or "").strip(),
        "size_arcmin": f"{(nx * scale / 60.0) if (nx and scale) else 0.0:.3f}",
        "pixel_scale_arcsec": f"{scale:.3f}",
        "status": "ok", "downloaded_utc": _utc_now_iso(),
        "source": "build_tile_plate_map_from_headers",
        "plate_id": plate_id,
        "notes": "backfilled from FITS header sidecar",
    }


def _as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _center_from_header(hdr) -> tuple[float, float]:
    """Tile center from CRVAL (tiles are fetched centered on their position)."""
    return (_as_float(hdr.get("CRVAL1")) or 0.0,
            _as_float(hdr.get("CRVAL2")) or 0.0)


if __name__ == "__main__":
    main()
