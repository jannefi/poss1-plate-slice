# -*- coding: utf-8 -*-
"""vasco.utils.tile_metadata

Utilities to keep dataset metadata up to date *during* Step1-download.

Retires separate post-download scripts that were previously required to:
  A) maintain data/metadata/tile_to_plate.csv (tile -> plate_id/REGION)
  B) write per-tile raw/dss1red_title.txt sidecar
  C) maintain data/metadata/tiles_registry.csv

Notes
-----
- plate_id is frozen to FITS header REGION.
- We read REGION/PLTLABEL/PLATEID/DATE-OBS from the local FITS header JSON sidecar
  written next to the FITS by Step1.
"""

from __future__ import annotations

import csv
import json
import os
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Tuple

try:
    import fcntl  # type: ignore
except Exception:  # pragma: no cover
    fcntl = None


@dataclass
class TilePlateRow:
    tile_id: str
    plate_id: str  # == REGION
    irsa_region: str = ''
    irsa_platelabel: str = ''
    irsa_plateid: str = ''
    irsa_date_obs: str = ''
    tile_survey: str = ''
    tile_date_obs: str = ''
    tile_fits: str = ''
    irsa_filename: str = ''
    irsa_center_sep_deg: str = ''


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _lock_exclusive(fp):
    if fcntl is None:
        return
    try:
        fcntl.flock(fp.fileno(), fcntl.LOCK_EX)
    except Exception:
        return


def _unlock(fp):
    if fcntl is None:
        return
    try:
        fcntl.flock(fp.fileno(), fcntl.LOCK_UN)
    except Exception:
        return


@contextmanager
def _csv_write_lock(out: Path):
    """Serialise the *whole* read-modify-write of a shared metadata CSV.

    Previously these updaters locked only the temp file they were about to
    write, which gives no mutual exclusion: concurrent step1 downloads each
    read the same snapshot, then each replaced `out` with its own
    snapshot-plus-one-row, silently dropping every row added in between.
    They also raced on a shared `.tmp` path, so the loser's `tmp.replace()`
    could raise FileNotFoundError outright. The lock must be held across
    read -> modify -> write -> replace, and must live on a file that is
    never itself replaced.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    lock_path = out.with_suffix(out.suffix + '.lock')
    with lock_path.open('a+', encoding='utf-8') as lf:
        _lock_exclusive(lf)
        try:
            yield
        finally:
            _unlock(lf)


def _read_csv_rows(out: Path) -> Dict[str, dict]:
    """Read an existing metadata CSV into {tile_id: row}. Call under the lock."""
    rows: Dict[str, dict] = {}
    if out.exists() and out.stat().st_size > 0:
        try:
            with out.open('r', newline='', encoding='utf-8') as f:
                for r in csv.DictReader(f):
                    tid = (r.get('tile_id') or '').strip()
                    if tid:
                        rows[tid] = r
        except Exception:
            rows = {}
    return rows


def _write_csv_atomic(out: Path, fieldnames, rows: Dict[str, dict]) -> None:
    """Rewrite `out` atomically. Call under the lock.

    The temp name carries the pid so that a failure to acquire the lock
    (fcntl unavailable) degrades to lost rows rather than a hard crash.
    """
    tmp = out.with_suffix(out.suffix + f'.tmp.{os.getpid()}')
    with tmp.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for tid in sorted(rows.keys()):
            w.writerow({k: rows[tid].get(k, '') for k in fieldnames})
    tmp.replace(out)


def _read_header_sidecar(fits_path: Path) -> Tuple[dict, Optional[Path]]:
    sidecar = fits_path.with_suffix(fits_path.suffix + '.header.json')
    if not sidecar.exists() or sidecar.stat().st_size == 0:
        return {}, None
    try:
        payload = json.loads(sidecar.read_text(encoding='utf-8'))
        if isinstance(payload, dict):
            hdr = payload.get('header', payload)
            if isinstance(hdr, dict):
                return hdr, sidecar
    except Exception:
        pass
    return {}, sidecar


def ensure_metadata_dirs(data_root: Path) -> Path:
    meta = data_root / 'metadata'
    meta.mkdir(parents=True, exist_ok=True)
    return meta


def update_tile_to_plate_csv(meta_dir: Path, row: TilePlateRow, filename: str = 'tile_to_plate.csv') -> Path:
    out = meta_dir / filename
    fieldnames = [
        'tile_id', 'plate_id', 'tile_region', 'tile_survey', 'tile_date_obs', 'tile_fits',
        'irsa_region', 'irsa_filename', 'irsa_survey', 'irsa_platelabel', 'irsa_plateid',
        'irsa_date_obs', 'irsa_center_sep_deg'
    ]

    with _csv_write_lock(out):
        rows = _read_csv_rows(out)
        rows[row.tile_id] = {
            'tile_id': row.tile_id,
            'plate_id': row.plate_id,
            'tile_region': row.plate_id,
            'tile_survey': row.tile_survey,
            'tile_date_obs': row.tile_date_obs,
            'tile_fits': row.tile_fits,
            'irsa_region': row.irsa_region or row.plate_id,
            'irsa_filename': row.irsa_filename or row.tile_fits,
            'irsa_survey': row.tile_survey,
            'irsa_platelabel': row.irsa_platelabel,
            'irsa_plateid': row.irsa_plateid,
            'irsa_date_obs': row.irsa_date_obs or row.tile_date_obs,
            'irsa_center_sep_deg': row.irsa_center_sep_deg,
        }
        _write_csv_atomic(out, fieldnames, rows)
    return out


def update_tiles_registry(meta_dir: Path, *, tile_id: str, ra_deg: float, dec_deg: float,
                          survey: str, size_arcmin: float, pixel_scale_arcsec: float,
                          status: str = 'ok', source: str = 'step1-download', vasco_plateid: str, notes: str = '') -> Path:
    out = meta_dir / 'tiles_registry.csv'
    fieldnames = [
        'tile_id','ra_deg','dec_deg','survey','size_arcmin','pixel_scale_arcsec',
        'status','downloaded_utc','source', 'plate_id', 'notes'
    ]

    with _csv_write_lock(out):
        rows = _read_csv_rows(out)
        rows[tile_id] = {
            'tile_id': tile_id,
            'ra_deg': f'{ra_deg:.6f}',
            'dec_deg': f'{dec_deg:.6f}',
            'survey': survey,
            'size_arcmin': f'{float(size_arcmin):.3f}',
            'pixel_scale_arcsec': f'{float(pixel_scale_arcsec):.3f}',
            'status': status,
            'downloaded_utc': _utc_now_iso(),
            'source': source,
            'plate_id': vasco_plateid,
            'notes': notes,
        }
        _write_csv_atomic(out, fieldnames, rows)
    return out


def write_dss1red_title(tile_dir: Path, row: TilePlateRow, *, prefer_local_header: bool = True) -> Path:
    raw = tile_dir / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    title_path = raw / 'dss1red_title.txt'

    src_rel = ''
    if prefer_local_header and row.tile_fits:
        local_json = raw / f'{row.tile_fits}.header.json'
        if local_json.exists():
            try:
                src_rel = os.path.relpath(local_json, raw)
            except Exception:
                src_rel = local_json.name
    if not src_rel:
        src_rel = row.irsa_filename or row.tile_fits or ''

    content_lines = [
        f'PLTLABEL: {row.irsa_platelabel}',
        f'PLATEID: {row.irsa_plateid}',
        f'REGION: {row.plate_id}',
        f'DATE-OBS: {row.irsa_date_obs or row.tile_date_obs}',
        f'FITS: {row.irsa_filename or row.tile_fits}',
        f'SOURCE: {src_rel}',
        f'SEP_DEG: {row.irsa_center_sep_deg}',
    ]
    title_path.write_text("\n".join(content_lines) + '\n', encoding='utf-8', newline="\n")
    return title_path


def update_all_after_download(*, tile_dir: Path, fits_path: Path, tile_id: str,
                              ra_deg: float, dec_deg: float, survey: str,
                              size_arcmin: float, pixel_scale_arcsec: float,
                              data_root: Path, prefer_local_header: bool = True) -> dict:
    hdr, sidecar = _read_header_sidecar(fits_path)

    region = str(hdr.get('REGION','') or '').strip()
    platelabel = str(hdr.get('PLTLABEL','') or '').strip()
    plateid = str(hdr.get('PLATEID','') or '').strip()
    date_obs = str(hdr.get('DATE-OBS','') or '').strip()

    meta_dir = ensure_metadata_dirs(data_root)

    row = TilePlateRow(
        tile_id=tile_id,
        plate_id=region,
        irsa_region=region,
        irsa_platelabel=platelabel,
        irsa_plateid=plateid,
        irsa_date_obs=date_obs,
        tile_survey=survey,
        tile_date_obs=date_obs,
        tile_fits=fits_path.name,
        irsa_filename=fits_path.name,
        irsa_center_sep_deg='',
    )

    out_map = update_tile_to_plate_csv(meta_dir, row)
    out_reg = update_tiles_registry(meta_dir,
                                    tile_id=tile_id,
                                    ra_deg=ra_deg,
                                    dec_deg=dec_deg,
                                    survey=survey,
                                    size_arcmin=size_arcmin,
                                    pixel_scale_arcsec=pixel_scale_arcsec,
                                    status='ok',
                                    source='step1-download',
                                    vasco_plateid=region)
    out_title = write_dss1red_title(tile_dir, row, prefer_local_header=prefer_local_header)

    return {
        'tile_to_plate_csv': str(out_map),
        'tiles_registry_csv': str(out_reg),
        'dss1red_title_txt': str(out_title),
        'plate_id': region,
        'header_sidecar': str(sidecar) if sidecar else '',
    }
