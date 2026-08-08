#!/usr/bin/env python3
"""Quick SExtractor config sweep for repro/mnras-parity tuning.

Runs one or more single-pass SExtractor configs against a few
already-downloaded tiles and reports, per (config, tile):

  raw          detections from SExtractor
  circle30     detections within the 30' inscribed circle (pipeline PRE_S0 cut)
  after_gaia   circle30 minus Gaia matches   (5", find=best1, cached catalog)
  survivors    after_gaia minus PS1 matches  (5", find=best1, cached catalog)
  rec_raw      % of reference-catalog rows in the tile with a raw detection <=5"
  rec_off      same, after subtracting the per-tile median systematic offset

Uses only local data: tile FITS + cached gaia/ps1 neighbourhood CSVs.
The reference catalog path is a runtime argument on purpose — never
hardcode it (see context/REPRO_DEVIATIONS.md privacy note).

Example:
  micromamba run -n vasco-py311 python tools/sweep_sex_configs.py \
    --config configs/one_pass.sex --config configs/one_pass_v2.sex \
    --tile tile_RA359.795_DECp4.791 --tile tile_RA359.795_DECp8.791 \
    --ref-catalog /path/to/reference.csv --ref-plate XE524 \
    --ref-name-col Name --ref-ra-col Ra --ref-dec-col Dec
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree

from vasco.mnras.filters_mnras import apply_extract_filters, apply_morphology_filters

REPO = Path(__file__).resolve().parents[1]
MATCH_ARCSEC = 5.0
CIRCLE_ARCMIN = 30.0

# Same filter cfg as tools/apply_mnras_filters_dryrun.py -- keep in sync.
EXTRACT_CFG = {'flags_equal': 0, 'snr_win_min': 30.0}
MORPH_CFG = {
    'fwhm_lower': 2.0, 'fwhm_upper': 7.0, 'elongation_lt': 1.3,
    'spread_model_min': -0.002, 'sigma_clip': True, 'sigma_k': 2.0,
    'extent_delta_lt': 2.0, 'extent_min': 1.0,
}


def wrap_dra(ra, center):
    return ((np.asarray(ra) - center + 180.0) % 360.0) - 180.0


def run_sex(config: Path, fits: Path, workdir: Path) -> Path:
    """Run SExtractor in an isolated workdir; return path to output CSV."""
    workdir.mkdir(parents=True, exist_ok=True)
    # Stage the chosen config plus every support file any config might name.
    for f in [config, *REPO.glob('configs/*.conv'), REPO / 'configs' / 'default.nnw',
              REPO / 'configs' / 'sex_default.param']:
        shutil.copy2(f, workdir / f.name)
    ldac = workdir / 'sweep.ldac'
    cmd = [
        'sex', str(fits.resolve()),
        '-c', config.name,
        '-CATALOG_NAME', ldac.name,
        '-CATALOG_TYPE', 'FITS_LDAC',
        '-PSF_NAME', '',
        # Force the single-pass param set regardless of what the config says:
        # default.param includes SPREAD_MODEL, which needs a PSF model.
        '-PARAMETERS_NAME', 'sex_default.param',
    ]
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"sex failed ({config.name}): {r.stderr[-500:]}")
    out_csv = workdir / 'sweep.csv'
    for ext in ('#LDAC_OBJECTS', '#2', '#1'):
        p = subprocess.run(
            ['stilts', 'tcopy', f'in={ldac}{ext}', f'out={out_csv}', 'ofmt=csv'],
            capture_output=True, text=True)
        if p.returncode == 0 and out_csv.exists() and out_csv.stat().st_size > 0:
            return out_csv
    raise RuntimeError(f"stilts tcopy failed on {ldac}")


def read_radec(path: Path, ra_col: str, dec_col: str):
    ras, decs = [], []
    with path.open(newline='', encoding='utf-8', errors='ignore') as f:
        for row in csv.DictReader(f):
            try:
                ras.append(float(row[ra_col])); decs.append(float(row[dec_col]))
            except (ValueError, KeyError, TypeError):
                continue
    return np.array(ras), np.array(decs)


def stilts_veto(cand_csv: Path, catalog_csv: Path, out_csv: Path,
                cat_ra: str, cat_dec: str) -> int:
    """1not2 best1 xmatch (mirrors the pipeline veto); returns surviving rows."""
    r = subprocess.run(
        ['stilts', 'tskymatch2', f'in1={cand_csv}', f'in2={catalog_csv}',
         f'out={out_csv}', 'ra1=ALPHA_J2000', 'dec1=DELTA_J2000',
         f'ra2={cat_ra}', f'dec2={cat_dec}', f'error={MATCH_ARCSEC}',
         'join=1not2', 'find=best1', 'ofmt=csv'],
        capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"stilts veto failed: {r.stderr[-300:]}")
    with out_csv.open(errors='ignore') as f:
        return max(0, sum(1 for _ in f) - 1)


def tile_center_from_id(tile_id: str):
    # tile_RA359.795_DECp4.791
    ra = float(tile_id.split('_RA')[1].split('_')[0])
    d = tile_id.split('_DEC')[1]
    dec = float(d[1:]) * (1 if d[0] == 'p' else -1)
    return ra, dec


def evaluate(config: Path, tile_id: str, tiles_root: Path,
             ref_ra, ref_dec, tmp_root: Path) -> dict:
    tile_dir = tiles_root / tile_id
    fits = next((tile_dir / 'raw').glob('*.fits'))
    work = tmp_root / config.stem / tile_id
    sex_csv = run_sex(config, fits, work)

    sra, sdec = read_radec(sex_csv, 'ALPHA_J2000', 'DELTA_J2000')
    tra, tdec = tile_center_from_id(tile_id)
    cosd = np.cos(np.deg2rad(tdec))

    res = {'config': config.name, 'tile': tile_id, 'raw': len(sra)}

    # 30' circle cut (pipeline PRE_S0 semantics)
    d_center = np.hypot(wrap_dra(sra, tra) * cosd, sdec - tdec)
    keep = d_center <= CIRCLE_ARCMIN / 60.0
    res['circle30'] = int(keep.sum())
    circle_csv = work / 'circle.csv'
    with sex_csv.open(errors='ignore') as fi, circle_csv.open('w', newline='') as fo:
        rd = csv.DictReader(fi)
        wr = csv.DictWriter(fo, fieldnames=rd.fieldnames)
        wr.writeheader()
        for i, row in enumerate(rd):
            if i < len(keep) and keep[i]:
                wr.writerow(row)

    # Cached-catalog vetoes, same order as the pipeline
    catdir = tile_dir / 'catalogs'
    gaia = catdir / 'gaia_neighbourhood_at_plate.csv'
    if not (gaia.exists() and gaia.stat().st_size > 0):
        gaia = catdir / 'gaia_neighbourhood.csv'
    after_gaia = work / 'after_gaia.csv'
    res['after_gaia'] = stilts_veto(circle_csv, gaia, after_gaia, 'ra', 'dec')
    after_ps1_csv = work / 'after_ps1.csv'
    res['survivors'] = stilts_veto(after_gaia, catdir / 'ps1_neighbourhood.csv',
                                   after_ps1_csv, 'ra', 'dec')

    # Real MNRAS extract+morphology filters on top of the veto survivors --
    # answers "does this config change reach the final candidate list",
    # not just raw/veto-stage counts.
    if res['survivors'] > 0:
        try:
            tab = Table.read(str(after_ps1_csv), format='ascii.csv')
            t1 = apply_extract_filters(tab, cfg=EXTRACT_CFG)
            t2 = apply_morphology_filters(t1, cfg=MORPH_CFG)
            res['after_mnras'] = len(t2)
        except Exception as e:
            print(f"MNRAS filter step failed for {config.name}/{tile_id}: {e}", file=sys.stderr)
            res['after_mnras'] = -1
    else:
        res['after_mnras'] = 0

    # Recovery vs reference rows inside this tile's circle
    d_ref = np.hypot(wrap_dra(ref_ra, tra) * cosd, ref_dec - tdec)
    in_tile = d_ref <= CIRCLE_ARCMIN / 60.0
    rra, rdec = ref_ra[in_tile], ref_dec[in_tile]
    res['ref_rows'] = len(rra)
    if len(rra) == 0 or len(sra) == 0:
        res['rec_raw'] = res['rec_off'] = float('nan')
        return res

    tree = cKDTree(np.column_stack([wrap_dra(sra, tra) * cosd, sdec - tdec]))
    q = np.column_stack([wrap_dra(rra, tra) * cosd, rdec - tdec])
    dist, idx = tree.query(q, k=1)
    res['rec_raw'] = 100.0 * float((dist * 3600.0 <= MATCH_ARCSEC).sum()) / len(rra)

    close = dist * 3600.0 <= 30.0
    if close.sum() >= 3:
        off = np.median(q[close] - tree.data[idx[close]], axis=0)
        dist2, _ = tree.query(q - off, k=1)
        res['rec_off'] = 100.0 * float((dist2 * 3600.0 <= MATCH_ARCSEC).sum()) / len(rra)
        res['off_arcsec'] = (round(off[0] * 3600, 2), round(off[1] * 3600, 2))
    else:
        res['rec_off'] = res['rec_raw']
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--config', action='append', required=True, type=Path)
    ap.add_argument('--tile', action='append', required=True)
    ap.add_argument('--tiles-root', type=Path, default=Path('./data/tiles'))
    ap.add_argument('--ref-catalog', required=True, type=Path,
                    help='Reference catalog CSV (runtime arg only; never commit the path)')
    ap.add_argument('--ref-plate', required=True)
    ap.add_argument('--ref-name-col', default='Name')
    ap.add_argument('--ref-ra-col', default='Ra')
    ap.add_argument('--ref-dec-col', default='Dec')
    ap.add_argument('--workdir', type=Path, default=None,
                    help='Where to put sweep outputs (default: temp dir)')
    ap.add_argument('--workers', type=int, default=3,
                    help='Parallel (config, tile) evaluations (default: 3)')
    args = ap.parse_args()

    ras, decs = [], []
    with args.ref_catalog.open(newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            if row.get(args.ref_name_col) == args.ref_plate:
                try:
                    ras.append(float(row[args.ref_ra_col]))
                    decs.append(float(row[args.ref_dec_col]))
                except (ValueError, KeyError):
                    continue
    ref_ra, ref_dec = np.array(ras), np.array(decs)
    print(f"Reference rows for {args.ref_plate}: {len(ref_ra)}")

    tmp_root = args.workdir or Path(tempfile.mkdtemp(prefix='sex_sweep_'))
    print(f"Sweep workdir: {tmp_root}\n")

    jobs = [(cfg, tile) for cfg in args.config for tile in args.tile]
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(evaluate, cfg.resolve(), tile, args.tiles_root, ref_ra, ref_dec, tmp_root): (cfg, tile)
            for cfg, tile in jobs
        }
        for fut in as_completed(futures):
            cfg, tile = futures[fut]
            try:
                r = fut.result()
            except Exception as e:
                print(f"FAILED {cfg.name} x {tile}: {e}", file=sys.stderr)
                continue
            results.append(r)
            print(f"done: {cfg.name} x {tile}")

    order = {(cfg.name, tile): i for i, (cfg, tile) in enumerate(jobs)}
    results.sort(key=lambda r: order[(r['config'], r['tile'])])

    hdr = (f"\n{'config':22s} {'tile':26s} {'raw':>7s} {'circ30':>7s} {'aftGaia':>8s} "
           f"{'surviv':>7s} {'mnras':>6s} {'ref':>5s} {'rec_raw':>8s} {'rec_off':>8s}  offset\"(ra,dec)")
    print(hdr)
    print('-' * len(hdr))
    for r in results:
        print(f"{r['config']:22s} {r['tile']:26s} {r['raw']:7d} {r['circle30']:7d} "
              f"{r['after_gaia']:8d} {r['survivors']:7d} {r['after_mnras']:6d} {r['ref_rows']:5d} "
              f"{r['rec_raw']:7.1f}% {r['rec_off']:7.1f}%  {r.get('off_arcsec', '')}")

    print(f"\n{'config':22s} {'sum_raw':>8s} {'sum_surv':>9s} {'sum_mnras':>10s}")
    by_cfg = {}
    for r in results:
        d = by_cfg.setdefault(r['config'], {'raw': 0, 'survivors': 0, 'after_mnras': 0})
        d['raw'] += r['raw']; d['survivors'] += r['survivors']; d['after_mnras'] += r['after_mnras']
    for cfg in dict.fromkeys(r['config'] for r in results):
        d = by_cfg[cfg]
        print(f"{cfg:22s} {d['raw']:8d} {d['survivors']:9d} {d['after_mnras']:10d}")


if __name__ == '__main__':
    main()
