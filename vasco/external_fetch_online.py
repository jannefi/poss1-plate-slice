
from __future__ import annotations
from pathlib import Path
import csv
import requests
import os
import time


__all__ = [
    'fetch_gaia_neighbourhood',
    'fetch_ps1_neighbourhood',
]

# --- USNO-B (I/284) via VizieR TSV -------------------------------------------
def fetch_usnob_neighbourhood(tile_dir: Path | str,
                              ra_deg: float, dec_deg: float,
                              radius_arcmin: float,
                              *, max_rows: int = 200000,
                              timeout: float = 60.0) -> Path:
    """
    Fetch USNO-B1.0 neighborhood (I/284) via VizieR TSV -> CSV.
    Writes catalogs/usnob_neighbourhood.csv

    If VASCO_USNOB_CACHE is set, queries the local Parquet cache instead.

    Env vars:
        VASCO_USNOB_ATTEMPTS  number of retry attempts (default 3)
    """
    from vasco.local_cache_query import query_usnob
    cached = query_usnob(tile_dir, ra_deg, dec_deg, radius_arcmin)
    if cached is not None:
        return cached
    tile_dir = Path(tile_dir)
    out = tile_dir / 'catalogs' / 'usnob_neighbourhood.csv'
    out.parent.mkdir(parents=True, exist_ok=True)
    base = 'https://vizier.u-strasbg.fr/viz-bin/asu-tsv'
    cols = ['RAJ2000','DEJ2000','R1mag','R2mag','B1mag','B2mag','Imag','pmRA','pmDE','_r']
    params = {
        '-source': 'I/284/out',
        '-c': f'{ra_deg:.8f} {dec_deg:.8f}',
        '-c.r': f'{radius_arcmin:.6f}',
        '-out.max': str(int(max_rows)),
        '-out.add': '_r',
        '-out.form': 'dec',
        '-sort': '_r',
        '-out': ','.join(cols),
    }
    import requests, csv
    attempts = int(os.getenv('VASCO_USNOB_ATTEMPTS', '3'))
    last_exc = None
    for attempt in range(1, attempts + 1):
        try:
            r = requests.get(base, params=params, timeout=timeout)
            r.raise_for_status()
            lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith('#')]
            if not lines:
                out.write_text('ra,dec\n', encoding='utf-8')
                return out
            header = lines[0].split('\t')
            colmap = {name: idx for idx, name in enumerate(header)}
            keep = [c for c in cols if c in colmap]
            if 'RAJ2000' not in colmap or 'DEJ2000' not in colmap:
                raise ValueError(
                    f"VizieR USNO-B response missing expected columns; "
                    f"got header: {header[:5]}"
                )
            data_rows = [ln.split('\t') for ln in lines[1:]]
            with out.open('w', newline='', encoding='utf-8') as f:
                w = csv.writer(f)
                w.writerow(['ra' if c=='RAJ2000' else ('dec' if c=='DEJ2000' else c) for c in keep])
                for row in data_rows:
                    try:
                        float(row[colmap['RAJ2000']]); float(row[colmap['DEJ2000']])
                    except Exception:
                        continue
                    w.writerow([row[colmap[c]] for c in keep])
            return out
        except Exception as e:
            last_exc = e
            if attempt < attempts:
                time.sleep(min(10, 1.5 ** attempt))
    raise last_exc

# --- PS1 cap override (optional) ----------------------------------------------
def _ps1_effective_radius_deg(radius_arcmin: float) -> float:
    # Previous cap was 0.5 deg; allow override via env
    cap = float(os.getenv('VASCO_PS1_MAX_RADIUS_DEG', '0.5'))
    return min(float(radius_arcmin) / 60.0, cap)

# ---------------------------------------------------------------
# Gaia DR3 via CDS/VizieR (normalize to ra/dec; skip unit rows)
# ---------------------------------------------------------------
_VIZIER_TSV = 'https://vizier.u-strasbg.fr/viz-bin/asu-tsv'
_DEF_GAIA_COLS = [
    'RA_ICRS', 'DE_ICRS', 'Gmag', 'BPmag', 'RPmag', 'pmRA', 'pmDE', 'Plx', '_r'
]

def fetch_gaia_neighbourhood(tile_dir: Path | str,
                             ra_deg: float, dec_deg: float,
                             radius_arcmin: float,
                             *, max_rows: int = 200000,
                             timeout: float = 60.0) -> Path:
    """Fetch Gaia DR3 neighborhood via VizieR TSV -> CSV; map RA/Dec and filter unit rows.

    If VASCO_GAIA_CACHE is set, queries the local Parquet cache instead.
    """
    from vasco.local_cache_query import query_gaia
    cached = query_gaia(tile_dir, ra_deg, dec_deg, radius_arcmin)
    if cached is not None:
        return cached
    tile_dir = Path(tile_dir)
    out = tile_dir / 'catalogs' / 'gaia_neighbourhood.csv'
    out.parent.mkdir(parents=True, exist_ok=True)

    params = {
        '-source': 'I/355/gaiadr3',
        '-c': f'{ra_deg} {dec_deg}',
        '-c.r': f'{radius_arcmin:.6f}',
        '-out.max': str(max_rows),
        '-out.add': '_r',
        '-out.form': 'dec',
        '-sort': '_r',
        '-out': ','.join(_DEF_GAIA_COLS),
    }
    r = requests.get(_VIZIER_TSV, params=params, timeout=timeout)
    r.raise_for_status()

    # Keep non-empty, non-comment lines
    lines = [row for row in r.text.splitlines() if row and not row.startswith('#')]
    if not lines:
        out.write_text('ra,dec')
        return out

    header = lines[0].split('	')
    data_rows = [ln.split('	') for ln in lines[1:]]
    colmap = {name: idx for idx, name in enumerate(header)}
    keep = [c for c in _DEF_GAIA_COLS if c in colmap]

    with out.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['ra' if c=='RA_ICRS' else ('dec' if c=='DE_ICRS' else c) for c in keep])
        for row in data_rows:
            try:
                # Ensure RA/DEC are numeric (skip unit rows like 'deg')
                float(row[colmap['RA_ICRS']]); float(row[colmap['DE_ICRS']])
            except Exception:
                continue
            w.writerow([row[colmap[c]] for c in keep])
    return out

# ---------------------------------------------------------------
# PS1 DR2 via Vizier Catalogs API (mean table, explicit columns)
# ---------------------------------------------------------------

def fetch_ps1_neighbourhood(tile_dir: Path | str,
                             ra_deg: float, dec_deg: float,
                             radius_arcmin: float,
                             *, max_records: int = 200000,
                             timeout: float = 60.0) -> Path:
    """Fetch PS1 DR2 neighborhood via CDS/VizieR API with explicit columns and progress logs.

    If VASCO_PS1_CACHE is set, queries the local Parquet cache instead.

    Honors environment variables:
        VASCO_PS1_CACHE (local Parquet cache — bypasses VizieR entirely)
        VASCO_PS1_RADIUS_DEG
        VASCO_PS1_TIMEOUT
        VASCO_PS1_ATTEMPTS
        VASCO_PS1_COLUMNS
    """
    from vasco.local_cache_query import query_ps1
    cached = query_ps1(tile_dir, ra_deg, dec_deg, radius_arcmin)
    if cached is not None:
        return cached
    tile_dir = Path(tile_dir)
    out = tile_dir / 'catalogs' / 'ps1_neighbourhood.csv'
    out.parent.mkdir(parents=True, exist_ok=True)

    # Radius (degrees). Cap to 0.5 deg; allow override for dev/testing.
    radius_deg = None
    _r = os.getenv('VASCO_PS1_RADIUS_DEG')
    if _r:
        try:
            radius_deg = float(_r)
        except Exception:
            radius_deg = None
    if radius_deg is None:
        radius_deg = min(float(radius_arcmin) / 60.0, 0.5)

    base = 'https://vizier.u-strasbg.fr/viz-bin/asu-tsv'
    url = base

    _timeout = float(os.getenv('VASCO_PS1_TIMEOUT', timeout))
    _attempts = int(os.getenv('VASCO_PS1_ATTEMPTS', '3'))

    default_cols = [
        'objID', 'RAJ2000', 'DEJ2000', 'Nd', 'gmag', 'rmag', 'imag', 'zmag', 'ymag', '_r'
    ]
    _cols_override = os.getenv('VASCO_PS1_COLUMNS')
    cols = [c.strip() for c in _cols_override.split(',')] if _cols_override else default_cols

    params = {
        '-source': 'II/389/ps1_dr2',
        '-c': f'{ra_deg:.8f} {dec_deg:.8f}',
        '-c.rd': f'{radius_deg:.8f}',  # degrees; use -c.rd not -c.r (arcmin)
        '-out.max': str(int(max_records)),
        '-out.add': '_r',
        '-out.form': 'dec',
        '-sort': '_r',
        '-out': ','.join(cols)
    }

    def _try_once():
        t0 = time.time()
        print(f'[POST][PS1] GET {url} (timeout={_timeout}s, radius={radius_deg}deg)')
        r = requests.get(url, params=params, timeout=_timeout)
        r.raise_for_status()
        dt = time.time() - t0
        print(f'[POST][PS1] OK in {dt:.2f}s ({len(r.content)} bytes)')
        return r

    last_exc = None
    for k in range(1, _attempts + 1):
        try:
            r = _try_once()
            # Parse VizieR TSV response into clean CSV (strip comment/unit/separator rows)
            lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith('#')]
            if not lines:
                out.write_text('ra,dec\n', encoding='utf-8')
                return out
            header = lines[0].split('\t')
            colmap = {name.strip(): idx for idx, name in enumerate(header)}
            keep = [c for c in cols if c in colmap]
            with out.open('w', newline='', encoding='utf-8') as f:
                w = csv.writer(f)
                w.writerow(['ra' if c == 'RAJ2000' else ('dec' if c == 'DEJ2000' else c) for c in keep])
                for ln in lines[1:]:
                    row = ln.split('\t')
                    try:
                        float(row[colmap['RAJ2000']]); float(row[colmap['DEJ2000']])
                    except Exception:
                        continue  # skip unit row ('deg'), separator row ('---'), malformed
                    w.writerow([row[colmap[c]] if colmap[c] < len(row) else '' for c in keep])
            return out
        except Exception as e:
            last_exc = e
            print(f'[POST][WARN] PS1 attempt {k}/{_attempts} failed: {e}')
            if k < _attempts:
                time.sleep(min(10, 1.5 ** k))
    raise last_exc
