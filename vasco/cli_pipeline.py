# -*- coding: utf-8 -*-
from __future__ import annotations
import argparse, json, time, subprocess, os, shutil, math
from pathlib import Path
try:
    from dotenv import load_dotenv; load_dotenv()  # reads .env if present (won't override explicit exports)
except ImportError:
    pass  # python-dotenv not installed; set env vars explicitly or: pip install python-dotenv
from typing import List, Tuple
import warnings
from vasco.utils.tile_id import parse_tile_id_center, format_tile_id
from vasco.utils.tile_metadata import update_all_after_download

# Silence pyerfa/ERFA warnings that clutter runs (must run before importing astropy/erfa)
warnings.filterwarnings(
    "ignore",
    message=r'ERFA function ".*" yielded .*dubious year.*',
    module=r"erfa\.core",
)
warnings.filterwarnings(
    "ignore",
    message=r'ERFA function ".*" yielded .*distance overridden.*',
    module=r"erfa\.core",
)

from . import downloader as dl

_PIXEL_SCALE_ARCSEC = 1.7  # DSS POSS-I E native scale; not user-configurable
from .exporter3 import export_and_summarize
from .utils.coords import parse_ra as _parse_ra, parse_dec as _parse_dec
from .pipeline_split import run_pass1, run_psfex, run_pass2, run_single_pass
from vasco.external_fetch_online import (fetch_gaia_neighbourhood, fetch_ps1_neighbourhood)
from vasco.external_fetch_usnob_vizier import fetch_usnob_neighbourhood
from vasco.mnras.xmatch_stilts import (xmatch_sextractor_with_gaia, xmatch_sextractor_with_ps1)
from vasco.utils.stilts_wrapper import stilts_xmatch

# --- NEW imports for MNRAS modules (filters, spikes, HPM, buckets/report) ---
from astropy.table import Table
from vasco.mnras.filters_mnras import apply_extract_filters, apply_morphology_filters
from vasco.mnras.spikes import (
    fetch_bright_ps1, fetch_bright_usnob, apply_spike_cuts,
    SpikeConfig, SpikeRuleConst, SpikeRuleLine
)
from vasco.mnras.hpm import backprop_gaia_row
from vasco.mnras.buckets import init_buckets, finalize
from vasco.mnras.report import write_summary
import csv as _csv
from vasco.wcsfix_early import ensure_wcsfix_catalog, WcsFixConfig



# --- helpers ---
def _ensure_tool_cli(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool '{tool}' not found in PATH.")

# --- POSSI-E enforcement & header export ---
def _fits_survey(path: Path) -> str:
    from astropy.io import fits
    try:
        with fits.open(path, memmap=False) as hdul:
            hdr = hdul[0].header if hdul and hdul[0].header else {}
            return str(hdr.get('SURVEY','')).strip()
    except Exception:
        return ''

def _write_fits_header_json(fits_path: Path) -> Path:
    """Write a JSON sidecar with selected keys and full header dump next to FITS."""
    from astropy.io import fits
    import json as _json
    fits_path = Path(fits_path)
    sidecar = fits_path.with_suffix(fits_path.suffix + '.header.json')
    try:
        with fits.open(fits_path, memmap=False) as hdul:
            hdr = hdul[0].header if hdul and hdul[0].header else {}
            sel_keys = ['SURVEY','PLATEID','PLATE-ID','PLATE','DATE-OBS','RA','DEC','EQUINOX','MJD-OBS',
                        'NAXIS1','NAXIS2','CD1_1','CD1_2','CD2_1','CD2_2','CDELT1','CDELT2',
                        'CRPIX1','CRPIX2','CRVAL1','CRVAL2']
            selected = {k: (str(hdr.get(k)) if hdr.get(k) is not None else None) for k in sel_keys}
            full = {str(k): (str(hdr.get(k)) if hdr.get(k) is not None else None) for k in hdr.keys()}
            payload = {'fits_file': fits_path.name, 'selected': selected, 'header': full}
            sidecar.write_text(_json.dumps(payload, indent=2), encoding='utf-8')
    except Exception:
        sidecar.write_text(json.dumps({'fits_file': fits_path.name, 'error': 'header_read_failed'}),
                           encoding='utf-8')
    return sidecar

def _enforce_possi_e_or_skip(fits_path: Path, logger) -> None:
    """Check FITS SURVEY and skip (delete + raise) if not POSSI-E; else write header sidecar."""
    survey = _fits_survey(fits_path)
    if survey != 'POSSI-E':
        try:
            msg = f"[STEP1][FILTER] Non-POSS plate; SURVEY={survey!r} — file will be discarded"
            if logger: logger.info(msg)
            else: print(msg)
        finally:
            try:
                Path(fits_path).unlink(missing_ok=True)
            except Exception:
                pass
        raise RuntimeError(f"Non-POSS plate returned by STScI: SURVEY={survey!r}")
    else:
        sidecar = _write_fits_header_json(Path(fits_path))
        if logger: logger.info(f"[STEP1][HEADER] Wrote FITS header sidecar: {sidecar.name}")
        else: print(f"[STEP1][HEADER] Wrote FITS header sidecar: {sidecar.name}")

# --- small helpers ---
def _expected_stem(ra: float, dec: float, survey: str, size_arcmin: float) -> str:
    sv_name = dl.SURVEY_ALIASES.get(survey.lower(), survey)
    tag = sv_name.lower().replace(' ', '-')
    return f"{tag}_{ra:.6f}_{dec:.6f}_{int(round(size_arcmin))}arcmin"

def _to_float_ra(val: str | float) -> float:
    try:
        return float(val)
    except Exception:
        return float(_parse_ra(str(val)))

def _to_float_dec(val: str | float) -> float:
    try:
        return float(val)
    except Exception:
        return float(_parse_dec(str(val)))

def _read_bright_cache(path: Path):
    import csv as _csv
    from vasco.mnras.spikes import BrightStar
    out = []
    if not path.exists() or path.stat().st_size == 0:
        return out
    with path.open(newline='', encoding='utf-8') as f:
        r = _csv.DictReader(f)
        for row in r:
            try:
                out.append(BrightStar(
                    ra=float(row['ra']),
                    dec=float(row['dec']),
                    rmag=float(row['rmag']),
                ))
            except Exception:
                continue
    return out

def _write_bright_cache(path: Path, bright):
    import csv as _csv
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = _csv.DictWriter(f, fieldnames=['ra','dec','rmag'])
        w.writeheader()
        for b in bright:
            w.writerow({'ra': b.ra, 'dec': b.dec, 'rmag': b.rmag})

# --- NEW: reporting helpers for veto-first Step4a/4b ---

def _csv_data_rows(path: Path) -> int:
    """Return number of data rows (excluding header) for a CSV file; 0 if missing/empty; -1 if unreadable."""
    try:
        p = Path(path)
        if not p.exists() or p.stat().st_size == 0:
            return 0
        with p.open('r', newline='', encoding='utf-8', errors='ignore') as f:
            n = sum(1 for _ in f)
        return max(0, n - 1)
    except Exception:
        return -1

def _csv_has_data_row_fast(path: Path) -> bool:
    """
    True if CSV has at least one data row (not just header).
    Reads only the first two lines -> O(1) and avoids scanning big files.
    """
    try:
        p = Path(path)
        if not p.exists() or p.stat().st_size == 0:
            return False
        with p.open("r", encoding="utf-8", errors="ignore") as f:
            hdr = f.readline()
            if not hdr:
                return False
            second = f.readline()
            return bool(second)
    except Exception:
        return False

def _augment_summary_json(tile_dir: Path, extra: dict) -> None:
    """Merge extra fields into <tile>/MNRAS_SUMMARY.json if it exists (best-effort)."""
    try:
        tile_dir = Path(tile_dir)
        p = tile_dir / 'MNRAS_SUMMARY.json'
        if not p.exists() or p.stat().st_size == 0:
            return
        base = json.loads(p.read_text(encoding='utf-8'))
        if not isinstance(base, dict):
            return
        base.update(extra)
        p.write_text(json.dumps(base, indent=2), encoding='utf-8')
    except Exception:
        return

def _analyze_rejection_reasons(
    src_csv: Path,
    out_rej_extract: Path,
    out_rej_morph: Path,
) -> dict:
    """
    Analyze remainder rows against the SAME extract + morphology rules as filters_mnras.py
    and write rejection CSVs with a reject_reason column.

    Rules mirrored from filters_mnras.py:
      Extract:
        - FLAGS == 0
        - SNR_WIN > 30   (note: strict '>' per filters_mnras.py)
      Morphology:
        - optional robust sigma-clip (default True): FWHM_IMAGE and ELONGATION within k*MAD
        - SPREAD_MODEL > -0.002
        - 2 < FWHM_IMAGE < 7
        - ELONGATION < 1.3
        - extent (if XMIN/XMAX/YMIN/YMAX exist):
            abs((XMAX-XMIN) - (YMAX-YMIN)) < 2
            (XMAX-XMIN) > 1 and (YMAX-YMIN) > 1
    """
    import csv
    import numpy as np

    src_csv = Path(src_csv)
    out_rej_extract = Path(out_rej_extract)
    out_rej_morph = Path(out_rej_morph)

    counts = {
        "late_reject_flags": 0,
        "late_reject_snr": 0,
        "late_reject_sigma_clip_fwhm": 0,
        "late_reject_sigma_clip_elongation": 0,
        "late_reject_spread_model": 0,
        "late_reject_fwhm": 0,
        "late_reject_elongation": 0,
        "late_reject_extent_delta": 0,
        "late_reject_extent_min": 0,
        "late_kept_hard_gates": 0,
    }

    if (not src_csv.exists()) or src_csv.stat().st_size == 0:
        out_rej_extract.write_text("", encoding="utf-8")
        out_rej_morph.write_text("", encoding="utf-8")
        return counts

    with src_csv.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.DictReader(f)
        rows = list(r)
        fieldnames = r.fieldnames or []

    if not fieldnames:
        out_rej_extract.write_text("", encoding="utf-8")
        out_rej_morph.write_text("", encoding="utf-8")
        return counts

    def _f(row, key, default=None):
        try:
            v = row.get(key, "")
            if v is None or v == "":
                return default
            return float(v)
        except Exception:
            return default

    # --- robust sigma clip (mirrors filters_mnras._robust_sigma_clip) ---
    def _robust_sigma_clip(x: np.ndarray, k: float = 2.0) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        med = np.nanmedian(x)
        mad = np.nanmedian(np.abs(x - med))
        sigma = 1.4826 * mad
        if not np.isfinite(sigma) or sigma <= 0:
            return np.isfinite(x)
        return np.isfinite(x) & (np.abs(x - med) <= k * sigma)

    sigma_clip = True  # matches filters_mnras default
    sigma_k = 2.0
    # repro/mnras-parity deviation #1: single-pass has no SPREAD_MODEL column;
    # skip that reason so this audit CSV doesn't mislabel rows as rejected
    # for a gate that was never evaluated.
    single_pass_mode = bool(os.getenv('VASCO_REPRO_SINGLE_PASS'))

    # Precompute sigma-clip masks over the whole remainder, if columns exist
    sc_fwhm = None
    sc_elong = None
    if sigma_clip and ("FWHM_IMAGE" in fieldnames):
        fwhm_all = np.array([_f(row, "FWHM_IMAGE", np.nan) for row in rows], dtype=float)
        sc_fwhm = _robust_sigma_clip(fwhm_all, k=sigma_k)
    if sigma_clip and ("ELONGATION" in fieldnames):
        elong_all = np.array([_f(row, "ELONGATION", np.nan) for row in rows], dtype=float)
        sc_elong = _robust_sigma_clip(elong_all, k=sigma_k)

    have_extent = {"XMAX_IMAGE", "XMIN_IMAGE", "YMAX_IMAGE", "YMIN_IMAGE"}.issubset(set(fieldnames))

    rej_extract = []
    rej_morph = []

    for i, row in enumerate(rows):
        reasons_extract = []
        reasons_morph = []

        # --- Extract rules (mirrors filters_mnras.apply_extract_filters) ---
        flags = _f(row, "FLAGS", None)
        snr = _f(row, "SNR_WIN", None)

        if flags is None or flags != 0:
            reasons_extract.append("flags")
        # NOTE: filters_mnras uses '>' not '>='
        if snr is None or not (snr > 30.0):
            reasons_extract.append("snr")

        if reasons_extract:
            if "flags" in reasons_extract:
                counts["late_reject_flags"] += 1
            if "snr" in reasons_extract:
                counts["late_reject_snr"] += 1
            rr = dict(row)
            rr["reject_reason"] = "|".join(reasons_extract)
            rej_extract.append(rr)
            continue

        # --- Sigma clipping (mirrors filters_mnras behavior) ---
        if sc_fwhm is not None and (not bool(sc_fwhm[i])):
            reasons_morph.append("sigma_clip_fwhm")
        if sc_elong is not None and (not bool(sc_elong[i])):
            reasons_morph.append("sigma_clip_elongation")

        # --- Hard morphology rules ---
        if not single_pass_mode:
            spread = _f(row, "SPREAD_MODEL", None)
            if spread is None or not (spread > -0.002):
                reasons_morph.append("spread_model")

        fwhm = _f(row, "FWHM_IMAGE", None)
        if fwhm is None or not (fwhm > 2.0 and fwhm < 7.0):
            reasons_morph.append("fwhm")

        elong = _f(row, "ELONGATION", None)
        if elong is None or not (elong < 1.3):
            reasons_morph.append("elongation")

        # Extent rules (correct rule: abs(dx - dy) < 2, dx>1, dy>1)
        if have_extent:
            xmax = _f(row, "XMAX_IMAGE", None)
            xmin = _f(row, "XMIN_IMAGE", None)
            ymax = _f(row, "YMAX_IMAGE", None)
            ymin = _f(row, "YMIN_IMAGE", None)
            if None not in (xmax, xmin, ymax, ymin):
                dx = xmax - xmin
                dy = ymax - ymin
                if not (np.isfinite(dx) and np.isfinite(dy)):
                    reasons_morph.append("extent_nan")
                else:
                    if not (dx > 1.0 and dy > 1.0):
                        reasons_morph.append("extent_min")
                    if not (abs(dx - dy) < 2.0):
                        reasons_morph.append("extent_delta")

        # --- Tally ---
        if reasons_morph:
            # counts by reason
            for rname in reasons_morph:
                if rname == "sigma_clip_fwhm":
                    counts["late_reject_sigma_clip_fwhm"] += 1
                elif rname == "sigma_clip_elongation":
                    counts["late_reject_sigma_clip_elongation"] += 1
                elif rname == "spread_model":
                    counts["late_reject_spread_model"] += 1
                elif rname == "fwhm":
                    counts["late_reject_fwhm"] += 1
                elif rname == "elongation":
                    counts["late_reject_elongation"] += 1
                elif rname == "extent_min":
                    counts["late_reject_extent_min"] += 1
                elif rname == "extent_delta":
                    counts["late_reject_extent_delta"] += 1
            rr = dict(row)
            rr["reject_reason"] = "|".join(reasons_morph)
            rej_morph.append(rr)
        else:
            counts["late_kept_hard_gates"] += 1

    out_rej_extract.parent.mkdir(parents=True, exist_ok=True)
    out_rej_morph.parent.mkdir(parents=True, exist_ok=True)

    def _write(path, out_rows):
        if not out_rows:
            path.write_text("", encoding="utf-8")
            return
        fns = list(out_rows[0].keys())
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fns)
            w.writeheader()
            w.writerows(out_rows)

    _write(out_rej_extract, rej_extract)
    _write(out_rej_morph, rej_morph)
    return counts

# --- NEW: MNRAS integration helpers ---

def _apply_mnras_filters_and_spikes(tile_dir: Path, sex_csv: Path, buckets: dict) -> Path:
    """
    Apply SNR/FLAGS + morphology filters and bright-star/diffraction spike cuts
    to a COPY of the SExtractor CSV before any cross-matching.
    Returns: catalogs/sextractor_pass2.filtered.csv
    """
    from astropy.table import Table
    import csv as _csv

    tile_dir = Path(tile_dir)
    sex_csv = Path(sex_csv)
    out_csv = sex_csv.with_name('sextractor_pass2.filtered.csv')

    # Load original (unfiltered)
    try:
        tab = Table.read(str(sex_csv), format='ascii.csv')
    except Exception:
        out_csv.write_text('', encoding='utf-8')
        return out_csv

    n0 = len(tab)

    if os.getenv('VASCO_DISABLE_MNRAS_FILTERS'):
        # repro/mnras-parity deviation #5 (see context/REPRO_DEVIATIONS.md):
        # pass every remainder row through unchanged (spike removal below
        # still runs).
        print('[POST][INFO]', tile_dir.name,
              'MNRAS extract/morphology filters disabled by env (VASCO_DISABLE_MNRAS_FILTERS) -> pass-through')
    else:
        # 1) FLAGS==0 & SNR_WIN>=30
        tab = apply_extract_filters(tab, cfg={'flags_equal': 0, 'snr_win_min': 30.0})

        # 2) Morphology gates (PSF-aware)
        tab = apply_morphology_filters(
            tab,
            cfg={
                'fwhm_lower': 2.0,
                'fwhm_upper': 7.0,
                'elongation_lt': 1.3,
                'spread_model_min': -0.002,
                # keep the paper-style robust clipping enabled (it defaults True in filters_mnras.py)
                'sigma_clip': True,
                'sigma_k': 2.0,
                # enable the paper’s pixel-extent guards when columns exist (your configs do output them)
                'extent_delta_lt': 2.0,
                'extent_min': 1.0,
            }
        )

    n1 = len(tab)
    buckets['morphology_rejected'] += max(0, n0 - n1)
    # Early exit: if nothing survives morphology, don't fetch bright-star catalogs
    if n1 == 0:
        # Preserve a valid CSV artifact (empty) and a valid rejected artifact (empty)
        out_csv.write_text('', encoding='utf-8')
        rej_path = tile_dir / 'catalogs' / 'sextractor_spike_rejected.csv'
        rej_path.write_text('', encoding='utf-8')
        return out_csv

    # Write intermediate
    tab.write(str(out_csv), format='ascii.csv', overwrite=True)

    # 3) Bright-star spike removal via PS1 (r<=16, fetch covers full tile + spike margin)
    # Spike physics (per 02_DECISIONS.md): reject candidates within 90″ of a bright star.
    # The fetch is anchored at tile center and must reach every candidate that could be
    # spike-rejected. For a 60×60′ square tile the far corner is 30·√2 ≈ 42.43′ from
    # center; adding the 90″ (1.5′) spike radius gives ≈43.9′ → rounded up to 45′.
    # search_radius_arcmin=1.5 enforces the 90″ paper rule as the per-candidate gate.

    # Which catalogue supplies the bright stars. PS1 is the historical default
    # (an engineering substitution, deviation #4); 'usnob' restores the
    # catalogue Solano et al. (2022) actually use, and makes the magnitude
    # constants below photometrically consistent -- they are the paper's, and
    # the paper calibrated them on USNO-B R, not PS1 r.
    # The 90" association radius and the rules themselves do NOT change with
    # this switch: exactly one variable moves.
    spike_cat = (os.getenv('VASCO_SPIKE_CATALOG') or 'ps1').strip().lower()
    if spike_cat not in ('ps1', 'usnob'):
        raise SystemExit(f"[FATAL] VASCO_SPIKE_CATALOG={spike_cat!r} is not one of ps1|usnob")
    # Checked HERE, not inside the try below: that try ends in
    # `except Exception: bright = []`, which would swallow the mirror-missing
    # error and silently run with NO spike mask at all -- a worse outcome than
    # either catalogue. Fail before entering it.
    if spike_cat == 'usnob' and not os.getenv('VASCO_USNOB_CACHE'):
        raise SystemExit(
            "[FATAL] VASCO_SPIKE_CATALOG=usnob needs VASCO_USNOB_CACHE set "
            "(the USNO-B1.0 Parquet mirror). Refusing to run: the fallback "
            "path would silently produce an unmasked catalogue.")

    center = _tile_center_from_index_or_name(tile_dir)
    bright = []
    if center:
        # Cache key carries the catalogue, so a PS1 cache can never be read as
        # USNO-B (or vice versa) after a switch.
        cache_path = (tile_dir / 'catalogs' /
                      f'{spike_cat}_bright_stars_r16_rad45.csv')
        try:
            # Use cache if present
            if cache_path.exists() and cache_path.stat().st_size > 0:
                bright = _read_bright_cache(cache_path)
            elif spike_cat == 'usnob':
                # Mirror-only by design; raises rather than degrading silently.
                bright = fetch_bright_usnob(
                    center[0], center[1], radius_arcmin=45.0, rmag_max=16.0
                )
                _write_bright_cache(cache_path, bright)
            else:
                bright = fetch_bright_ps1(
                    center[0], center[1],
                    radius_arcmin=45.0, rmag_max=16.0, mindetections=2
                )
                # Save cache for fast reruns
                _write_bright_cache(cache_path, bright)
        except Exception:
            bright = []

        with open(out_csv, newline='') as f:
            rdr = _csv.DictReader(f)
            rows = list(rdr)

        kept, rejected = apply_spike_cuts(
            rows, bright,
            SpikeConfig(
                search_radius_arcmin=1.5,  # 90″ per paper (02_DECISIONS.md)
                rules=[
                    SpikeRuleConst(const_max_mag=12.4),
                    SpikeRuleLine(a=-0.09, b=15.3),  # slope per arcsec
                ],
            )
        )
        buckets['spikes_rejected'] += len(rejected)
        # The spike stage used to log nothing at all, so after a run there was no
        # record of whether the mask ran or which catalogue supplied its stars --
        # recoverable only by inferring it from the direction of the row count.
        # An empty star list is called out separately: `except Exception:
        # bright = []` above degrades to NO MASK, which removes nothing and so
        # reads downstream as a permissive catalogue rather than as a failure.
        if not bright:
            print('[POST][WARN]', tile_dir.name,
                  f'spike mask ran with ZERO bright stars (catalogue={spike_cat}) '
                  f'-- nothing was masked; treat this tile as unmasked, not as a result',
                  flush=True)
        else:
            print('[POST][INFO]', tile_dir.name,
                  f'spike mask: catalogue={spike_cat} bright={len(bright)} '
                  f'rejected={len(rejected)} kept={len(kept)}', flush=True)

        # Write final filtered rows (may be empty, header preserved)
        fieldnames = (kept[0].keys() if kept else (rows[0].keys() if rows else []))
        with open(out_csv, 'w', newline='') as fo:
            if fieldnames:
                w = _csv.DictWriter(fo, fieldnames=fieldnames)
                w.writeheader(); w.writerows(kept)
            else:
                fo.write('')

        # Diagnostics
        rej_path = tile_dir / 'catalogs' / 'sextractor_spike_rejected.csv'
        with rej_path.open('w', newline='') as fo:
            if rejected:
                w = _csv.DictWriter(fo, fieldnames=rejected[0].keys())
                w.writeheader(); w.writerows(rejected)
            else:
                fo.write('')

    return out_csv

    # HPM filtering helpers unchanged

def _sep_arcsec(ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float) -> float:
    ra1 = math.radians(ra1_deg); dec1 = math.radians(dec1_deg)
    ra2 = math.radians(ra2_deg); dec2 = math.radians(dec2_deg)
    s = 2*math.asin(math.sqrt(
        math.sin((dec2-dec1)/2)**2 +
        math.cos(dec1)*math.cos(dec2)*math.sin((ra2-ra1)/2)**2
    ))
    return math.degrees(s) * 3600.0

def _filter_hpm_gaia(xdir: Path, buckets: dict, poss_sep_arcsec: float = 5.0) -> None:
    """After Gaia xmatch, back-propagate Gaia positions to POSS epoch; flag HPM mismatches."""
    gx = xdir / 'sex_gaia_xmatch.csv'
    if not gx.exists():
        return
    with gx.open(newline='') as f:
        rdr = _csv.DictReader(f)
        rows = list(rdr)
    kept, flagged = [], []
    for row in rows:
        try:
            ra_bp, dec_bp = backprop_gaia_row(row, target_epoch=1950.0)
            poss_ra = float(row.get('ALPHA_J2000', row.get('ra', 'nan')))
            poss_dec = float(row.get('DELTA_J2000', row.get('dec', 'nan')))
            if any(map(lambda x: isinstance(x, float) and math.isnan(x),
                       [ra_bp, dec_bp, poss_ra, poss_dec])):
                kept.append(row); continue
            sep = _sep_arcsec(poss_ra, poss_dec, ra_bp, dec_bp)
            if sep <= poss_sep_arcsec:
                kept.append(row)
            else:
                r2 = dict(row); r2['hpm_sep_arcsec'] = f"{sep:.3f}"
                flagged.append(r2)
        except Exception:
            kept.append(row)
    buckets['hpm_objects'] += len(flagged)
    # Write only the small flagged-HPM audit file; hpmclean is not used downstream.
    out_flag = xdir / 'sex_gaia_hpm_flagged.csv'
    if flagged:
        with out_flag.open('w', newline='') as fo:
            w = _csv.DictWriter(fo, fieldnames=flagged[0].keys())
            w.writeheader(); w.writerows(flagged)
    else:
        out_flag.write_text('', encoding='utf-8')
    # Delete the large source xmatch file now that HPM is done (not used downstream).
    try:
        gx.unlink()
    except Exception:
        pass

# --- commands ---

def cmd_one(args: argparse.Namespace) -> int:
    """
    One-shot pipeline: 1+2+3 + export + xmatch, with deferred tile creation in step 1.
    Avoids materializing data/tiles/<tileid> on download errors (non-FITS/non-WCS/non-POSS).
    """
    # DO NOT pre-create tile directory; keep logger central (optional)
    lg = dl.configure_logger(Path('./data/logs'))
    ra = _to_float_ra(args.ra)
    dec = _to_float_dec(args.dec)
    tile_id = format_tile_id(ra, dec)
    run_dir = Path(args.workdir) / tile_id  # tile folder name is always computed, never caller-supplied

    def _cache_ok(path: Path) -> bool:
        """Return True if cache CSV exists, is non-empty, and has RA/Dec columns."""
        try:
            p = Path(path)
            return p.exists() and p.stat().st_size > 0 and _csv_has_radec(p)
        except Exception:
            return False

    # --- STEP 1: download with deferral (downloader will stage & promote only on success)
    try:
        fits = dl.fetch_skyview_dss(
            ra, dec,
            size_arcmin=args.size_arcmin,
            survey=args.survey,
            pixel_scale_arcsec=_PIXEL_SCALE_ARCSEC,
            out_dir=run_dir / 'raw',  # downloader creates this only on success
            logger=lg
        )
        # Enforce POSSI-E post-promotion; may unlink & raise if not POSSI-E
        _enforce_possi_e_or_skip(Path(fits), lg)
    except RuntimeError as e:
        # Non-POSS enforcement path keeps your original bookkeeping — but only
        # write RUN_* artifacts if the tile folder already exists.
        if 'Non-POSS plate returned by STScI' in str(e):
            print('[SKIP]', f'RA={ra:.6f}', f'Dec={dec:.6f}', '-> non-POSS; tile omitted.')
            counts = {'planned': 1, 'downloaded': 0, 'processed': 0, 'filtered_non_poss': 1}
            missing = [{
                'ra': float(ra),
                'dec': float(dec),
                'expected_stem': _expected_stem(ra, dec, args.survey, args.size_arcmin)
            }]
            if run_dir.exists():
                _write_json(run_dir / 'RUN_COUNTS.json', counts)
                _write_json(run_dir / 'RUN_MISSING.json', missing)
                _write_json(run_dir / 'RUN_INDEX.json', [])
                _write_overview(run_dir, counts, [], missing)
            return 0
        # For non-FITS / non-WCS / other failures: downloader already wrote error artifacts
        print('[STEP1][ERROR]', str(e))
        return 1

    # --- SUCCESS PATH ---
    # STEP 2 + STEP 3
    if os.getenv('VASCO_REPRO_SINGLE_PASS'):
        # repro/mnras-parity deviation #1: single-pass SExtractor, no PSFEx,
        # no SPREAD_MODEL. Requires VASCO_DISABLE_MNRAS_FILTERS (deviation #5)
        # to also be set, or the SPREAD_MODEL morphology gate rejects everything.
        if not os.getenv('VASCO_DISABLE_MNRAS_FILTERS'):
            print('[STEP2/3][WARN] VASCO_REPRO_SINGLE_PASS set without '
                  'VASCO_DISABLE_MNRAS_FILTERS -> SPREAD_MODEL gate will reject all candidates')
        p2 = run_single_pass(fits, run_dir, config_root='configs')
    else:
        p1, _ = run_pass1(fits, run_dir, config_root='configs')
        psf = run_psfex(p1, run_dir, config_root='configs')
        p2 = run_pass2(fits, run_dir, psf, config_root='configs')

    # LDAC isn't directly readable by most tools; extract to CSV right away
    # (stilts tcopy, multi-HDU probing) rather than waiting for step4-xmatch.
    sex_csv = _ensure_sextractor_csv(run_dir, p2)
    print('[STEP2/3] catalog CSV ->', sex_csv)

    # Exports & QA
    export_and_summarize(p2, run_dir, export=args.export, histogram_col=args.hist_col)

    # STEP 4: xmatch (local)
    # +3' margin beyond circumscribed circle to prevent edge leakage (e.g. source 2319 X=2031)
    radius_arcmin = args.size_arcmin * (2 ** 0.5) * 0.5 + 3.0

    # Cache-aware fetch-on-miss: avoid network calls if cache already exists
    gaia_cache = run_dir / 'catalogs' / 'gaia_neighbourhood.csv'
    try:
        if os.getenv('VASCO_FORCE_FETCH_GAIA'):
            fetch_gaia_neighbourhood(run_dir, ra, dec, radius_arcmin)
        elif _cache_ok(gaia_cache):
            print('[POST][INFO]', run_dir.name, 'Gaia cache present — skipping fetch')
        else:
            fetch_gaia_neighbourhood(run_dir, ra, dec, radius_arcmin)
    except Exception as e:
        print('[POST][WARN]', run_dir.name, 'Gaia fetch failed:', e)

    ps1_cache = run_dir / 'catalogs' / 'ps1_neighbourhood.csv'
    try:
        if os.getenv('VASCO_DISABLE_PS1'):
            print('[POST][INFO]', run_dir.name, 'PS1 disabled by env — skipping fetch')
        else:
            if os.getenv('VASCO_FORCE_FETCH_PS1'):
                fetch_ps1_neighbourhood(run_dir, ra, dec, radius_arcmin)
            elif _cache_ok(ps1_cache):
                print('[POST][INFO]', run_dir.name, 'PS1 cache present — skipping fetch')
            else:
                fetch_ps1_neighbourhood(run_dir, ra, dec, radius_arcmin)
    except Exception as e:
        print('[POST][WARN]', run_dir.name, 'PS1 fetch failed:', e)

    usnob_cache = run_dir / 'catalogs' / 'usnob_neighbourhood.csv'
    try:
        if os.getenv('VASCO_DISABLE_USNOB'):
            print('[POST][INFO]', run_dir.name, 'USNO-B disabled by env — skipping fetch')
        else:
            if os.getenv('VASCO_FORCE_FETCH_USNOB'):
                fetch_usnob_neighbourhood(run_dir, ra, dec, radius_arcmin)
                print('[POST]', run_dir.name, 'USNO-B (VizieR) -> catalogs/usnob_neighbourhood.csv')
            elif _cache_ok(usnob_cache):
                print('[POST][INFO]', run_dir.name, 'USNO-B cache present — skipping fetch')
            else:
                fetch_usnob_neighbourhood(run_dir, ra, dec, radius_arcmin)
                print('[POST]', run_dir.name, 'USNO-B (VizieR) -> catalogs/usnob_neighbourhood.csv')
    except Exception as e:
        print('[POST][WARN]', run_dir.name, 'USNO-B fetch failed:', e)

    try:
        _post_xmatch_tile(run_dir, p2, radius_arcsec=float(args.xmatch_radius_arcsec))
    except Exception as e:
        print('[POST][WARN] xmatch failed for', run_dir.name, ':', e)

    # Final run bookkeeping & overview (tile dir exists on success)
    results = [{'tile': Path(fits).stem, 'pass1': str(p1), 'psf': str(psf), 'pass2': str(p2)}]
    counts = {'planned': 1, 'downloaded': 1, 'processed': 1, 'filtered_non_poss': 0}
    _write_json(run_dir / 'RUN_INDEX.json', results)
    _write_json(run_dir / 'RUN_COUNTS.json', counts)
    _write_json(run_dir / 'RUN_MISSING.json', [])
    _write_overview(run_dir, counts, results, [])
    print('Run directory:', run_dir)
    return 0


# --- helpers for run bookkeeping ---
def _build_run_dir(base: str | Path | None = None) -> Path:
    base = Path(base) if base else Path('data') / 'runs'
    base.mkdir(parents=True, exist_ok=True)
    return base

def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding='utf-8')

def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2), encoding='utf-8')

import datetime as _datetime

def _update_tile_status(run_dir: Path, step: str, status: str) -> None:
    """Merge-update tile_status.json with the result of one pipeline step.

    Status values: "ok", "skip", "fail".
    Reads existing file (if any), updates only the given step key, then
    re-writes.  Failures are caught and printed as warnings so they never
    abort a pipeline step.
    """
    path = run_dir / 'tile_status.json'
    try:
        if path.exists():
            data = json.loads(path.read_text(encoding='utf-8'))
        else:
            data = {'tile_id': run_dir.name, 'steps': {}}
        data['steps'][step] = {
            'status': status,
            'ts': _datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
        }
        _write_json(path, data)
    except Exception as e:
        print(f'[WARN] tile_status update failed for {step}:', e)

# --- step2, step3, post-xmatch, cds-xmatch, step4, step5, step6 ---

def cmd_step2_pass1(args: argparse.Namespace) -> int:
    run_dir = _build_run_dir(Path(args.workdir) if args.workdir else None)
    raw = run_dir / 'raw'
    fits = next((p for p in sorted(raw.glob('*.fits'))), None)
    if not fits:
        print('[STEP2][ERROR] No FITS in raw/. Run step1-download first.')
        _update_tile_status(run_dir, 'step2', 'fail')
        return 2
    if os.getenv('VASCO_REPRO_SINGLE_PASS'):
        # repro/mnras-parity deviation #1: skip pass1/PSFEx/pass2 split;
        # step3-psf-and-pass2 becomes a no-op when this is set (see there).
        p2 = run_single_pass(fits, run_dir, config_root='configs')
        _write_json(run_dir / 'RUN_INDEX.json', [{'tile': Path(fits).stem, 'single_pass': str(p2)}])
        print('[STEP2] single-pass ->', p2)
        # LDAC isn't directly readable; extract to CSV right away (stilts tcopy).
        sex_csv = _ensure_sextractor_csv(run_dir, p2)
        print('[STEP2] catalog CSV ->', sex_csv)
        _update_tile_status(run_dir, 'step2', 'ok')
        return 0
    p1, _ = run_pass1(fits, run_dir, config_root='configs')
    _write_json(run_dir / 'RUN_INDEX.json', [{'tile': Path(fits).stem, 'pass1': str(p1)}])
    print('[STEP2] pass1 ->', p1)
    _update_tile_status(run_dir, 'step2', 'ok')
    return 0

def cmd_step3_psf_and_pass2(args: argparse.Namespace) -> int:
    run_dir = _build_run_dir(Path(args.workdir) if args.workdir else None)
    if os.getenv('VASCO_REPRO_SINGLE_PASS'):
        # step2-pass1 already produced pass2.ldac directly; nothing to do here.
        print('[STEP3] single-pass mode -> pass2.ldac already produced by step2-pass1, skipping')
        _update_tile_status(run_dir, 'step3', 'ok')
        return 0
    p1 = run_dir / 'pass1.ldac'
    raw = run_dir / 'raw'
    fits = next((p for p in sorted(raw.glob('*.fits'))), None)
    if not p1.exists() or not fits:
        print('[STEP3][ERROR] pass1.ldac or FITS missing. Run step2-pass1 first.')
        _update_tile_status(run_dir, 'step3', 'fail')
        return 2
    psf = run_psfex(p1, run_dir, config_root='configs')
    p2 = run_pass2(fits, run_dir, psf, config_root='configs')
    print('[STEP3] psf ->', psf, '; pass2 ->', p2)
    # LDAC isn't directly readable; extract to CSV right away (stilts tcopy).
    sex_csv = _ensure_sextractor_csv(run_dir, p2)
    print('[STEP3] catalog CSV ->', sex_csv)
    _update_tile_status(run_dir, 'step3', 'ok')
    return 0

def _csv_has_radec(csv_path: Path) -> bool:
    import csv
    try:
        with open(csv_path, newline='') as f:
            hdr = next(csv.reader(f))
        cols = {h.strip() for h in hdr}
        for a,b in [('RA_corr','Dec_corr'), ('RA_corr','DEC_corr'),
                ('ra','dec'), ('RA_ICRS','DE_ICRS'), ('RAJ2000','DEJ2000'),
                ('RA','DEC'), ('lon','lat'), ('raMean','decMean'), ('RAMean','DecMean'),
                ('ALPHA_J2000','DELTA_J2000'), ('ALPHAWIN_J2000','DELTAWIN_J2000'),
                ('X_WORLD','Y_WORLD')]:
            if a in cols and b in cols:
                return True
        return False
    except Exception:
        return False

def _detect_radec_columns(csv_path: Path) -> Tuple[str, str] | None:
    import csv
    try:
        with open(csv_path, newline='') as f:
            hdr = next(csv.reader(f))
        cols = [h.strip() for h in hdr]
        pairs = [('RA_corr','Dec_corr'), ('RA_corr','DEC_corr'),
             ('ALPHA_J2000','DELTA_J2000'), ('ALPHAWIN_J2000','DELTAWIN_J2000'),
             ('X_WORLD','Y_WORLD'), ('RAJ2000','DEJ2000'), ('RA_ICRS','DE_ICRS'),
             ('ra','dec'), ('RA','DEC')]
        for a,b in pairs:
            if a in cols and b in cols:
                return a,b
        return None
    except Exception:
        return None

# NEW March 2026
# circular filter

def _wrap_deg_pm180(x: float) -> float:
    return (x + 180.0) % 360.0 - 180.0

def _sep_arcsec_haversine(ra1_deg: float, dec1_deg: float, ra2_deg: float, dec2_deg: float) -> float:
    """Great-circle separation in arcsec, wrap-aware in RA."""
    ra1 = math.radians(ra1_deg)
    dec1 = math.radians(dec1_deg)
    ra2 = math.radians(ra2_deg)
    dec2 = math.radians(dec2_deg)

    d_ra = math.radians(_wrap_deg_pm180(ra2_deg - ra1_deg))
    d_dec = dec2 - dec1

    a = math.sin(d_dec/2.0)**2 + math.cos(dec1) * math.cos(dec2) * math.sin(d_ra/2.0)**2
    c = 2.0 * math.asin(min(1.0, math.sqrt(a)))
    return math.degrees(c) * 3600.0

def _circle_filter_csv_inplace(
    csv_path: Path,
    *,
    center_ra_deg: float,
    center_dec_deg: float,
    radius_arcmin: float,
    before_copy: Path,
    status_json: Path,
) -> None:
    """
    Filter csv_path to rows within radius_arcmin of center, overwriting csv_path.
    Writes before_copy (once) and status_json marker.
    Streaming implementation (does not load full CSV into memory).
    """
    radius_arcsec = float(radius_arcmin) * 60.0

    # If already filtered with the same radius, skip
    if status_json.exists():
        try:
            st = json.loads(status_json.read_text(encoding="utf-8"))
            if float(st.get("radius_arcmin", -1)) == float(radius_arcmin) and st.get("state") == "done":
                return
        except Exception:
            pass

    # Preserve original square once
    if (not before_copy.exists()) and csv_path.exists() and csv_path.stat().st_size > 0:
        before_copy.write_bytes(csv_path.read_bytes())

    # Detect RA/Dec columns from header
    with csv_path.open("r", newline="", encoding="utf-8", errors="ignore") as fi:
        rdr = _csv.reader(fi)
        header = next(rdr, [])
    cols = [h.strip().lstrip("﻿") for h in header]
    colset = set(cols)

    candidates = [
        ("RA_corr", "Dec_corr"),
        ("ALPHAWIN_J2000", "DELTAWIN_J2000"),
        ("ALPHA_J2000", "DELTA_J2000"),
        ("X_WORLD", "Y_WORLD"),
    ]
    ra_col = dec_col = None
    for a, b in candidates:
        if a in colset and b in colset:
            ra_col, dec_col = a, b
            break
    if ra_col is None:
        raise RuntimeError("circle filter: could not find RA/Dec columns in sextractor_pass2.csv")

    tmp = csv_path.with_suffix(csv_path.suffix + ".circle.tmp")
    kept = 0
    seen = 0

    with csv_path.open("r", newline="", encoding="utf-8", errors="ignore") as fi, \
         tmp.open("w", newline="", encoding="utf-8") as fo:
        dr = _csv.DictReader(fi)
        if not dr.fieldnames:
            return
        dw = _csv.DictWriter(fo, fieldnames=dr.fieldnames)
        dw.writeheader()

        for row in dr:
            seen += 1
            try:
                ra = float(row.get(ra_col, "nan"))
                dec = float(row.get(dec_col, "nan"))
                if not (math.isfinite(ra) and math.isfinite(dec)):
                    continue
                sep_arcsec = _sep_arcsec_haversine(center_ra_deg, center_dec_deg, ra, dec)
                if sep_arcsec <= radius_arcsec:
                    dw.writerow(row)
                    kept += 1
            except Exception:
                continue

    tmp.replace(csv_path)

    status = {
        "state": "done",
        "radius_arcmin": float(radius_arcmin),
        "center_ra_deg": float(center_ra_deg),
        "center_dec_deg": float(center_dec_deg),
        "rows_in": int(seen),
        "rows_out": int(kept),
    }
    status_json.write_text(json.dumps(status, indent=2), encoding="utf-8")
    
# UPDATED ensure: schema + rows aware

def _ensure_sextractor_csv(tile_dir: Path, pass2_ldac: str | Path) -> Path:
    """
    Ensure catalogs/sextractor_pass2.csv exists and is usable:
      - non-empty
      - contains RA/Dec
      - has at least one data row
      - includes required columns used by filters
    If invalid, re-extract from LDAC using STILTS (multiple HDUs).

    Set VASCO_LDAC_DROP_VIGNET=1 to drop the VIGNET postage-stamp column during
    the conversion: ~16x faster and ~98% smaller, bit-identical on every retained
    column, with automatic fall-through to the legacy full copy. See the comment
    at the fast path below before changing it.
    """
    _ensure_tool_cli('stilts')
    import csv

    tile_dir = Path(tile_dir)
    pass2_ldac = Path(pass2_ldac)

    cat_dir = tile_dir / 'catalogs'
    cat_dir.mkdir(parents=True, exist_ok=True)
    sex_csv = cat_dir / 'sextractor_pass2.csv'
    probe = cat_dir / '_probe.csv'

    REQUIRED_COLS = {
        'ALPHA_J2000', 'DELTA_J2000',   # coordinates
        'FLAGS', 'SNR_WIN',             # extract-time
        'FWHM_IMAGE', 'ELONGATION',     # morphology
    }
    if not os.getenv('VASCO_REPRO_SINGLE_PASS'):
        # SPREAD_MODEL is PSF-aware and only exists after PSFEx (pass2).
        # repro/mnras-parity deviation #1 (single-pass) has no PSF model.
        REQUIRED_COLS.add('SPREAD_MODEL')

    def _header(path: Path) -> List[str]:
        try:
            with path.open(newline='') as f:
                return [c.strip() for c in next(csv.reader(f), [])]
        except Exception:
            return []

    def _valid(path: Path) -> bool:
        try:
            if not path.exists() or path.stat().st_size == 0:
                return False
        except Exception:
            return False

        hdr = _header(path)
        if len(hdr) <= 2:
            return False
        if not _csv_has_radec(path):
            return False
        # at least one data row?
        try:
            with path.open(newline='') as f:
                rdr = csv.reader(f)
                next(rdr, None)
                if next(rdr, None) is None:
                    return False
        except Exception:
            return False
        # required schema
        cols = set(hdr)
        if not REQUIRED_COLS.issubset(cols):
            return False
        return True

    # Fast path — but never reuse a CSV older than the LDAC it derives from:
    # a rerun SExtractor pass regenerates pass2.ldac, and silently keeping the
    # previous extraction would feed stale detections to everything downstream
    # (same staleness class as the _stage_to_run_folder config bug).
    def _stale(csv_path: Path) -> bool:
        try:
            return (pass2_ldac.exists()
                    and pass2_ldac.stat().st_mtime > csv_path.stat().st_mtime)
        except Exception:
            return False

    if _valid(sex_csv) and not _stale(sex_csv):
        return sex_csv

    # Re-extract with multi-HDU probing. Retried up to 3x: LDAC->CSV
    # extraction via stilts has shown intermittent failures under
    # concurrent (parallel-tile) load even when the LDAC itself is valid
    # and non-empty -- a fresh retry of the identical extraction reliably
    # succeeds (see repro-parity XE002 investigation, 2026-07-16).
    hdu_tries = ['#LDAC_OBJECTS', '#2', '#1', '#0', '#3', '#4', '#5', '#6', '#7', '#8', '']
    max_attempts = 3

    def _run_ladder(argv_for, attempts: int) -> bool:
        """Sweep every HDU selector, retrying the whole sweep. True if sex_csv written."""
        for attempt in range(1, attempts + 1):
            for ext in hdu_tries:
                in_arg = f"in={str(pass2_ldac)}{ext}" if ext else f"in={str(pass2_ldac)}"
                try:
                    subprocess.run(argv_for(in_arg), check=True, capture_output=True)
                except Exception:
                    continue
                if _valid(probe):
                    try:
                        probe.replace(sex_csv)
                    except Exception:
                        shutil.copyfile(probe, sex_csv)
                    return True
            if attempt < attempts:
                time.sleep(1.5)
        return False

    # Fast path: drop VIGNET during the conversion instead of after it.
    #
    # sex_default.param asks for VIGNET(45,45), i.e. 2025 floats per detection,
    # and tcopy serialises every one of them to CSV text -- 72.8 MB for 3,187
    # rows, of which the retained 35 columns are 1.2 MB. Nothing in this repo
    # reads VIGNET from the CSV; PSFEx reads it from the LDAC, and single-pass
    # runs no PSFEx at all. Measured on a real pass2.ldac: 3.20 s -> 0.20 s, and
    # the 35 retained columns come out bit-identical (max relative difference
    # 0.000e+00), so this is a speed change only.
    #
    # Deliberately NOT a Python/astropy reimplementation: that measured 0.70 s,
    # 3.5x slower than this, and astropy's CSV writer honours each column's FITS
    # display format, silently truncating THETA_IMAGE from 40.433807 to 40.43 on
    # 3,185 of 3,187 rows.
    #
    # One attempt only, then fall through to the untouched tcopy ladder below:
    # `delcols` is a hard error ('No column matching "VIGNET"') on any HDU that
    # lacks the column, so the fallback is load-bearing, not decorative. The
    # worst case of enabling this knob is exactly today's behaviour.
    if os.getenv('VASCO_LDAC_DROP_VIGNET'):
        if _run_ladder(lambda a: ['stilts', 'tpipe', a, 'cmd=delcols VIGNET',
                                  f'out={str(probe)}', 'ofmt=csv'], attempts=1):
            return sex_csv

    if _run_ladder(lambda a: ['stilts', 'tcopy', a,
                              f'out={str(probe)}', 'ofmt=csv'], max_attempts):
        return sex_csv

    # Every attempt across every HDU selector failed to produce a valid,
    # non-empty catalog. Do NOT silently write an empty placeholder here:
    # that made a transient extraction race indistinguishable from a
    # genuinely successful, genuinely-empty tile (both looked like
    # tile_status.json step=="ok" with a 0-row catalog). Raise instead so
    # the failure is loud and the tile is correctly marked "fail" by the
    # caller, per this project's audit-standard (every stage must record
    # in_rows/out_rows/failure reasons -- a silent empty catalog violates
    # that).
    try:
        ldac_exists = pass2_ldac.exists()
        ldac_size = pass2_ldac.stat().st_size if ldac_exists else 'n/a'
    except Exception:
        ldac_exists, ldac_size = 'unknown', 'unknown'
    raise RuntimeError(
        f"_ensure_sextractor_csv: failed to extract a valid catalog from "
        f"{pass2_ldac} after {max_attempts} attempt(s) x {len(hdu_tries)} "
        f"HDU selectors ({hdu_tries}). LDAC exists={ldac_exists} "
        f"size={ldac_size}."
    )


# ---------------------------------------------------------------------------
# Bug #6 fix: epoch-propagate catalog positions to plate epoch before matching.
# Without this, high-PM stars drift beyond the 5" cone gate over the ~50-60 yr
# plate-to-catalog epoch gap and leak through Gaia/USNO-B vetoes entirely.
# ---------------------------------------------------------------------------

def _plate_epoch_year_from_fits(fits_path: Path) -> float | None:
    """Return decimal year of DATE-OBS for the given FITS, or None on failure.

    DSS plates sometimes have overflowed time components (e.g. T11:77:00 where
    minutes=77).  Same normalization as stage_skybot_post.parse_dateobs.
    """
    import re as _re
    try:
        from astropy.io import fits as _fits
        from astropy.time import Time as _Time
        hdr = _fits.getheader(str(fits_path))
        date_obs = hdr.get('DATE-OBS') or hdr.get('DATE-AVG')
        if not date_obs:
            return None
        t = str(date_obs).strip()
        if ' ' in t and 'T' not in t:
            t = t.replace(' ', 'T')
        # Fast path: standard ISO
        try:
            return float(_Time(t).decimalyear)
        except Exception:
            pass
        # Slow path: normalize overflowed hh/mm/ss (common in DSS FITS headers)
        m = _re.match(r'^(\d{4}-\d{2}-\d{2})T(\d{2}):(\d{2}):(\d{2})', t)
        if not m:
            return None
        ymd = m.group(1)
        hh, mm, ss = int(m.group(2)), int(m.group(3)), int(m.group(4))
        mm_add, ss = divmod(ss, 60)
        mm += mm_add
        hh_add, mm = divmod(mm, 60)
        hh = (hh + hh_add) % 24
        return float(_Time(f'{ymd}T{hh:02d}:{mm:02d}:{ss:02d}').decimalyear)
    except Exception:
        return None


def _propagate_catalog_epoch(
    in_csv: Path,
    out_csv: Path,
    catalog_epoch_year: float,
    plate_epoch_year: float,
    ra_col: str,
    dec_col: str,
    pmra_col: str,
    pmde_col: str,
) -> bool:
    """Epoch-propagate a catalog CSV's RA/Dec to `plate_epoch_year` using PMs.

    Reads `in_csv`, writes `out_csv` with ra/dec columns updated in place for
    rows that have finite proper motions. Rows with missing or invalid PMs are
    written unchanged (conservative fall-through to raw catalog epoch).
    Returns True on success; leaves out_csv untouched and returns False on error.

    `pmra_col` is expected to be already multiplied by cos(dec) — this is the
    VizieR/Gaia convention where pmRA in mas/yr is the great-circle rate in the
    RA direction, not dRA/dt.
    """
    try:
        import pandas as pd
        import numpy as np
        df = pd.read_csv(in_csv)
        if not len(df):
            df.to_csv(out_csv, index=False)
            return True

        dt_yr = plate_epoch_year - catalog_epoch_year  # negative for POSS vs Gaia/USNO-B

        def _sf(x):
            try:
                return float(x)
            except (TypeError, ValueError):
                return float('nan')

        ra_deg = df[ra_col].apply(_sf).values
        dec_deg = df[dec_col].apply(_sf).values
        pmra = df[pmra_col].apply(_sf).values
        pmde = df[pmde_col].apply(_sf).values

        has_pm = (
            np.isfinite(pmra) & np.isfinite(pmde)
            & np.isfinite(ra_deg) & np.isfinite(dec_deg)
        )

        # mas/yr * yr / (1000 mas/arcsec * 3600 arcsec/deg) = deg
        ddec_deg = np.where(has_pm, pmde * dt_yr / 3.6e6, 0.0)
        # pmRA is cos(dec)-corrected, so divide by cos(dec) to get pure RA change
        cos_dec = np.cos(np.radians(np.where(has_pm, dec_deg, 0.0)))
        dra_deg = np.where(has_pm, (pmra * dt_yr / 3.6e6) / cos_dec, 0.0)

        df[ra_col] = ra_deg + dra_deg
        df[dec_col] = dec_deg + ddec_deg
        # Drop rows that ended up with NaN positions (unparseable input) so
        # downstream consumers like WCSFIX don't hit matmul warnings on NaNs.
        valid = np.isfinite(df[ra_col].astype(float).values) & np.isfinite(df[dec_col].astype(float).values)
        df = df.loc[valid].reset_index(drop=True)
        df.to_csv(out_csv, index=False)
        return True
    except Exception as e:
        print(f'[EPOCH-PROP][WARN] failed {in_csv.name} -> {out_csv.name}: {e}')
        return False


def _post_xmatch_tile(tile_dir, pass2_ldac, *, radius_arcsec: float = 5.0) -> None:
    """ Step4 (LOCAL backend) — veto-first ordering (Step4a/Step4b).

    Step4a (optical veto, elimination semantics):
      - Gaia veto: write xmatch/sex_gaia_xmatch.csv (matched pairs) and
        catalogs/sextractor_pass2.after_gaia_veto.csv (carry-forward unmatched)
      - PS1 veto on Gaia-unmatched: write xmatch/sex_ps1_xmatch.csv and
        catalogs/sextractor_pass2.after_ps1_veto.csv
      - USNO-B veto on PS1-unmatched: write xmatch/sex_usnob_xmatch.csv and
        catalogs/sextractor_pass2.after_usnob_veto.csv

    Step4b (late filters):
      - Apply extract/morphology filters + spike cuts on the post-veto remainder.

    Extra outputs added by this implementation:
      - catalogs/sextractor_pass2.after_usnob_veto.rejected_extract.csv
      - catalogs/sextractor_pass2.after_usnob_veto.rejected_morphology.csv
      - Additional fields merged into MNRAS_SUMMARY.json (veto counts + rejection counts)
    """
    tile_dir = Path(tile_dir)
    xdir = tile_dir / 'xmatch'
    xdir.mkdir(parents=True, exist_ok=True)
    catdir = tile_dir / 'catalogs'
    catdir.mkdir(parents=True, exist_ok=True)

    # 1) Ensure base SExtractor CSV exists
    sex_csv = _ensure_sextractor_csv(tile_dir, pass2_ldac)
    # 1b) Optional MNRAS-style circular cut at catalogue level (<= VASCO_CIRCLE_ARCMIN)
    # Rationale: run SExtractor on a square cutout (e.g., 60x60) but keep only detections within a
    # circular region (e.g., 30 arcmin radius) before any veto/xmatch.
    try:
        circle_arcmin = float(os.getenv('VASCO_CIRCLE_ARCMIN', '0') or '0')
    except Exception:
        circle_arcmin = 0.0

    if circle_arcmin > 0.0:
        center0 = _tile_center_from_index_or_name(tile_dir)
        if center0 is None:
            print('[STEP4][WARN]', tile_dir.name, 'VASCO_CIRCLE_ARCMIN set but tile center unavailable; skipping circle cut')
        else:
            before_copy = catdir / 'sextractor_pass2_before_circle_filter.csv'
            status_json = catdir / 'circle_filter_status.json'
            try:
                _circle_filter_csv_inplace(
                    sex_csv,
                    center_ra_deg=float(center0[0]),
                    center_dec_deg=float(center0[1]),
                    radius_arcmin=float(circle_arcmin),
                    before_copy=before_copy,
                    status_json=status_json,
                )
                # If filtering resulted in an empty catalog, keep it (parity behaviour) but warn loudly.
                try:
                    out_rows = _csv_data_rows(sex_csv)
                    print('[STEP4][INFO]', tile_dir.name, f'circle cut <= {circle_arcmin:.1f}′ applied -> rows={out_rows}')
                except Exception:
                    print('[STEP4][INFO]', tile_dir.name, f'circle cut <= {circle_arcmin:.1f}′ applied')
            except Exception as e:
                print('[STEP4][WARN]', tile_dir.name, 'circle cut failed; using square catalog:', e)

    # 2) Neighbourhood caches (fetch stage happens before this in cmd_step4_xmatch).
    # Prefer epoch-propagated variants (Bug #6 fix) with fall-through to the raw
    # catalogs — both WCSFIX bootstrap and the cone-match veto chain benefit from
    # matching against plate-epoch positions instead of raw J2016/J2000.
    def _prefer_plate(name: str) -> Path:
        p = catdir / f'{name}_at_plate.csv'
        r = catdir / f'{name}.csv'
        return p if (p.exists() and p.stat().st_size > 0) else r
    gaia_csv = _prefer_plate('gaia_neighbourhood')
    ps1_csv = catdir / 'ps1_neighbourhood.csv'   # PS1 has no PMs in this fetch
    usnob_csv = _prefer_plate('usnob_neighbourhood')

    buckets = init_buckets()

    # 3) Early canonical coordinates (WCSFIX) using local Gaia cache
    sex_for_veto = sex_csv
    if os.getenv('VASCO_WCSFIX_DISABLE'):
        # repro/mnras-parity deviation #6: use raw plate FITS WCS as-is.
        # sex_for_veto already defaults to sex_csv (raw coords); this is the
        # same self-degrading fallback path WCSFIX already uses on failure,
        # just made explicit/auditable instead of abusing threshold knobs.
        print('[POST][INFO]', tile_dir.name, 'WCSFIX disabled by env (VASCO_WCSFIX_DISABLE) -> using raw coords')
    else:
        try:
            import html as _html

            center = _tile_center_from_index_or_name(tile_dir)

            # Primary config (existing behavior; env-controlled)
            cfg = WcsFixConfig(
                bootstrap_radius_arcsec=float(os.getenv('VASCO_WCSFIX_BOOTSTRAP_ARCSEC', '5.0')),
                degree=int(os.getenv('VASCO_WCSFIX_DEGREE', '2')),
                min_matches=int(os.getenv('VASCO_WCSFIX_MIN_MATCHES', '20')),
            )

            # Optional fallback config (second-chance) for the specific failure mode:
            # "too few tie points (x < min_matches)"
            fallback_enabled = os.getenv('VASCO_WCSFIX_FALLBACK', '1').strip().lower() not in ('0', 'false', 'no')
            cfg_fb = WcsFixConfig(
                bootstrap_radius_arcsec=float(os.getenv('VASCO_WCSFIX_FALLBACK_BOOTSTRAP_ARCSEC', '15.0')),
                degree=int(os.getenv('VASCO_WCSFIX_FALLBACK_DEGREE', '1')),
                min_matches=int(os.getenv('VASCO_WCSFIX_FALLBACK_MIN_MATCHES', '10')),
            )

            if gaia_csv.exists() and gaia_csv.stat().st_size > 0:
                out_wcs, status = ensure_wcsfix_catalog(
                    tile_dir,
                    sex_csv,
                    gaia_csv,
                    center=center,
                    cfg=cfg,
                    force=bool(os.getenv('VASCO_WCSFIX_FORCE', '').strip()),
                )

                if status.get('ok'):
                    sex_for_veto = out_wcs
                    print('[POST]', tile_dir.name, 'WCSFIX OK ->', out_wcs.name)
                else:
                    reason = status.get('reason') or ''
                    reason_norm = _html.unescape(str(reason)).lower()
                    print('[POST][INFO]', tile_dir.name, 'WCSFIX skipped/failed -> using raw coords:', reason)

                    # Second-chance retry: only for the known failure mode
                    if fallback_enabled and ('too few tie points' in reason_norm):
                        print('[POST][INFO]', tile_dir.name,
                            f'WCSFIX fallback retry: bootstrap={cfg_fb.bootstrap_radius_arcsec}" '
                            f'min_matches={cfg_fb.min_matches} degree={cfg_fb.degree}')

                        out_wcs2, status2 = ensure_wcsfix_catalog(
                            tile_dir,
                            sex_csv,
                            gaia_csv,
                            center=center,
                            cfg=cfg_fb,
                            # If user explicitly set FORCE, honor it; otherwise do not force unnecessarily.
                            force=bool(os.getenv('VASCO_WCSFIX_FORCE', '').strip()),
                        )

                        if status2.get('ok'):
                            sex_for_veto = out_wcs2
                            print('[POST]', tile_dir.name, 'WCSFIX OK (fallback) ->', out_wcs2.name)
                        else:
                            print('[POST][INFO]', tile_dir.name,
                                'WCSFIX fallback failed -> using raw coords:', status2.get('reason'))
            else:
                print('[POST][INFO]', tile_dir.name, 'WCSFIX skipped: gaia_neighbourhood.csv missing/empty')
        except Exception as e:
            print('[POST][WARN]', tile_dir.name, 'WCSFIX error -> using raw coords:', e)

    # Candidate RA/Dec columns (prefer RA_corr/Dec_corr when present)
    cand_cols = _detect_radec_columns(sex_for_veto) or ('ALPHA_J2000', 'DELTA_J2000')
    ra1, dec1 = cand_cols

    def _copy_or_empty(src: Path, dst: Path) -> None:
        try:
            if src.exists() and src.stat().st_size > 0:
                shutil.copyfile(src, dst)
            else:
                dst.write_text('', encoding='utf-8')
        except Exception:
            try:
                dst.write_text('', encoding='utf-8')
            except Exception:
                pass

    def _veto(stage: str,
          in_candidates: Path,
          catalog: Path,
          out_match: Path,
          out_unmatched: Path,
          default_cat_cols: tuple[str, str],
          disable_env: str | None = None) -> bool:

        # Env disable: treat as skip, carry-forward unchanged
        if disable_env and os.getenv(disable_env):
            print('[POST][INFO]', tile_dir.name, f'{stage} disabled by env ({disable_env}) -> skipping veto')
            _copy_or_empty(in_candidates, out_unmatched)
            return False

        # If candidate input is empty/header-only, nothing to veto; carry-forward empty
        if not _csv_has_data_row_fast(in_candidates):
            print('[POST][INFO]', tile_dir.name, f'{stage} veto skipped: candidates empty (0 data rows)')
            out_unmatched.write_text('', encoding='utf-8')
            return False

        # Catalog must exist, have RA/Dec *and* have at least one data row
        if not (catalog.exists() and _csv_has_radec(catalog) and _csv_has_data_row_fast(catalog)):
            # This is the key fix for header-only ps1_neighbourhood/usnob_neighbourhood
            print('[POST][INFO]', tile_dir.name, f'{stage} veto skipped: catalogue missing/invalid/empty -> carry-forward unchanged')
            buckets.setdefault('missing_inputs', 0)
            buckets['missing_inputs'] += 1
            _copy_or_empty(in_candidates, out_unmatched)
            return False

        cat_cols = _detect_radec_columns(catalog) or default_cat_cols
        ra2, dec2 = cat_cols

        # Matched pairs (find=best1: each SExtractor source independently finds its best match)
        stilts_xmatch(
            str(in_candidates), str(catalog), str(out_match),
            ra1=ra1, dec1=dec1, ra2=ra2, dec2=dec2,
            radius_arcsec=radius_arcsec,
            join_type='1and2',
            find='best1',
            ofmt='csv',
        )

        # Carry-forward unmatched (find=best1: prevents two sources competing for same catalog star)
        stilts_xmatch(
            str(in_candidates), str(catalog), str(out_unmatched),
            ra1=ra1, dec1=dec1, ra2=ra2, dec2=dec2,
            radius_arcsec=radius_arcsec,
            join_type='1not2',
            find='best1',
            ofmt='csv',
        )

        print('[POST]', tile_dir.name, f'{stage} veto ->', out_match)
        return True

    # Step4a outputs
    after_gaia = catdir / 'sextractor_pass2.after_gaia_veto.csv'
    after_ps1 = catdir / 'sextractor_pass2.after_ps1_veto.csv'
    after_usnob = catdir / 'sextractor_pass2.after_usnob_veto.csv'

    out_gaia = xdir / 'sex_gaia_xmatch.csv'
    out_ps1 = xdir / 'sex_ps1_xmatch.csv'
    out_usnob = xdir / 'sex_usnob_xmatch.csv'

    # Veto accounting
    veto_start_rows = _csv_data_rows(sex_for_veto)

    gaia_ok = _veto('Gaia', sex_for_veto, gaia_csv, out_gaia, after_gaia, default_cat_cols=('ra', 'dec'))
    veto_after_gaia_rows = _csv_data_rows(after_gaia)

    ps1_ok = _veto('PS1', after_gaia, ps1_csv, out_ps1, after_ps1,
                   default_cat_cols=('raMean', 'decMean'), disable_env='VASCO_DISABLE_PS1')
    veto_after_ps1_rows = _csv_data_rows(after_ps1)

    usno_ok = _veto('USNO-B', after_ps1, usnob_csv, out_usnob, after_usnob,
                    default_cat_cols=('RAJ2000', 'DEJ2000'), disable_env='VASCO_DISABLE_USNOB')
    veto_after_usnob_rows = _csv_data_rows(after_usnob)

    # “Eliminated” (based on carry-forward rowcounts)
    veto_gaia_eliminated = max(0, veto_start_rows - veto_after_gaia_rows)
    veto_ps1_eliminated = max(0, veto_after_gaia_rows - veto_after_ps1_rows)
    veto_usnob_eliminated = max(0, veto_after_ps1_rows - veto_after_usnob_rows)

    # Pair counts (matched-pairs outputs; not unique candidates)
    veto_gaia_pairs = _csv_data_rows(out_gaia)
    veto_ps1_pairs = _csv_data_rows(out_ps1)
    veto_usnob_pairs = _csv_data_rows(out_usnob)

    # Step4b: late filters/spikes on the post-veto remainder
    remainder = after_usnob
    try:
        if (not remainder.exists()) or remainder.stat().st_size == 0:
            # Nothing survived the vetoes. Write the empty survivors file anyway:
            # _apply_mnras_filters_and_spikes is the only writer of filtered.csv
            # and is skipped by this early return, so without this a tile emptied
            # by a veto leaves NO file, while one emptied by the morphology gate
            # leaves a 0-byte one. Two representations of the same outcome, and
            # the absent one is indistinguishable from a tile lost to a crash.
            (catdir / 'sextractor_pass2.filtered.csv').write_text('', encoding='utf-8')
            write_summary(tile_dir, finalize(buckets), md_path='MNRAS_SUMMARY.md', json_path='MNRAS_SUMMARY.json')
            _augment_summary_json(tile_dir, {
                'veto_start_rows': veto_start_rows,
                'veto_after_gaia_rows': veto_after_gaia_rows,
                'veto_after_ps1_rows': veto_after_ps1_rows,
                'veto_after_usnob_rows': veto_after_usnob_rows,
                'veto_gaia_eliminated': veto_gaia_eliminated,
                'veto_ps1_eliminated': veto_ps1_eliminated,
                'veto_usnob_eliminated': veto_usnob_eliminated,
                'veto_gaia_pairs': veto_gaia_pairs,
                'veto_ps1_pairs': veto_ps1_pairs,
                'veto_usnob_pairs': veto_usnob_pairs,
                'veto_gaia_ok': bool(gaia_ok),
                'veto_ps1_ok': bool(ps1_ok),
                'veto_usnob_ok': bool(usno_ok),
            })
            return
    except Exception:
        write_summary(tile_dir, finalize(buckets), md_path='MNRAS_SUMMARY.md', json_path='MNRAS_SUMMARY.json')
        _augment_summary_json(tile_dir, {
            'veto_start_rows': veto_start_rows,
            'veto_after_gaia_rows': veto_after_gaia_rows,
            'veto_after_ps1_rows': veto_after_ps1_rows,
            'veto_after_usnob_rows': veto_after_usnob_rows,
            'veto_gaia_eliminated': veto_gaia_eliminated,
            'veto_ps1_eliminated': veto_ps1_eliminated,
            'veto_usnob_eliminated': veto_usnob_eliminated,
            'veto_gaia_pairs': veto_gaia_pairs,
            'veto_ps1_pairs': veto_ps1_pairs,
            'veto_usnob_pairs': veto_usnob_pairs,
            'veto_gaia_ok': bool(gaia_ok),
            'veto_ps1_ok': bool(ps1_ok),
            'veto_usnob_ok': bool(usno_ok),
        })
        return

    # Write rejection reasons for the post-veto remainder (should be small)
    rej_extract = catdir / 'sextractor_pass2.after_usnob_veto.rejected_extract.csv'
    rej_morph = catdir / 'sextractor_pass2.after_usnob_veto.rejected_morphology.csv'
    reason_counts = _analyze_rejection_reasons(remainder, rej_extract, rej_morph)

    # Apply the actual filters/spikes pipeline
    _apply_mnras_filters_and_spikes(tile_dir, remainder, buckets)

    # Keep HPM late and conservative (flagging only)
    try:
        _filter_hpm_gaia(xdir, buckets)
    except Exception:
        pass

    write_summary(tile_dir, finalize(buckets), md_path='MNRAS_SUMMARY.md', json_path='MNRAS_SUMMARY.json')

    # Augment summary with veto stage counts + hard-gate rejection counts + env flags
    extra = {
        'ps1_veto_enabled': not bool(os.getenv('VASCO_DISABLE_PS1')),
        'usnob_veto_enabled': not bool(os.getenv('VASCO_DISABLE_USNOB')),
        'veto_start_rows': veto_start_rows,
        'veto_after_gaia_rows': veto_after_gaia_rows,
        'veto_after_ps1_rows': veto_after_ps1_rows,
        'veto_after_usnob_rows': veto_after_usnob_rows,
        'veto_gaia_eliminated': veto_gaia_eliminated,
        'veto_ps1_eliminated': veto_ps1_eliminated,
        'veto_usnob_eliminated': veto_usnob_eliminated,
        'veto_gaia_pairs': veto_gaia_pairs,
        'veto_ps1_pairs': veto_ps1_pairs,
        'veto_usnob_pairs': veto_usnob_pairs,
        'veto_gaia_ok': bool(gaia_ok),
        'veto_ps1_ok': bool(ps1_ok),
        'veto_usnob_ok': bool(usno_ok),
    }
    extra.update(reason_counts)
    _augment_summary_json(tile_dir, extra)

# --- CDS logging & helpers ---

def _csv_row_count(path: Path) -> int:
    import csv
    try:
        with open(path, newline='') as f:
            r = csv.reader(f)
            next(r, None)
            return sum(1 for _ in r)
    except Exception:
        return -1

def _coords_from_tile_dirname(name: str) -> tuple[float, float] | None:
    return parse_tile_id_center(name)

def _tile_center_from_index_or_name(run_dir: Path) -> tuple[float, float] | None:
    try:
        recs = json.loads((Path(run_dir) / 'RUN_INDEX.json').read_text(encoding='utf-8'))
        if recs:
            stem = Path(recs[0].get('tile','')).name
            parts = stem.split('_')
            return float(parts[1]), float(parts[2])
    except Exception:
        pass
    return _coords_from_tile_dirname(Path(run_dir).name)

# --- CLI ---

def cmd_step4_xmatch(args: argparse.Namespace) -> int:
    run_dir = _build_run_dir(Path(args.workdir) if args.workdir else None)
    p2 = run_dir / 'pass2.ldac'
    if not p2.exists():
        # The LDAC's ONLY role downstream is producing sextractor_pass2.csv via
        # _ensure_sextractor_csv, which already declines to regenerate when the
        # CSV is newer. So a tile that still has a valid CSV needs no LDAC: the
        # veto, filters, spike mask and dedup all read the CSV.
        #
        # This matters because it makes an already-extracted run re-runnable
        # under a different configuration -- WCSFIX on/off, a different veto set,
        # a different spike catalogue -- without repeating SExtractor + PSFEx,
        # which is the only expensive stage (~25 min/plate). The harvested
        # catalogue is pre-veto and pre-spike-mask, so it is independent of all
        # of those choices.
        #
        # Announced loudly rather than silently: with no LDAC there is nothing
        # to check the CSV's freshness against, so the operator must know the
        # detections came from a prior extraction.
        fallback_csv = run_dir / 'catalogs' / 'sextractor_pass2.csv'
        if fallback_csv.exists() and fallback_csv.stat().st_size > 0:
            print(f'[STEP4][INFO] pass2.ldac absent; proceeding from the existing '
                  f'{fallback_csv.name} ({fallback_csv.stat().st_size} bytes). '
                  f'Detections come from a PRIOR extraction and cannot be '
                  f'freshness-checked against an LDAC.')
        else:
            print('[STEP4][ERROR] pass2.ldac missing and no usable '
                  'catalogs/sextractor_pass2.csv to fall back on. '
                  'Run step3-psf-and-pass2 first.')
            return 2

    import os

    def _truthy_env(name: str) -> bool:
        v = os.getenv(name, '').strip().lower()
        return v not in ('', '0', 'false', 'no')

    NO_FETCH = _truthy_env('VASCO_STEP4_NO_FETCH')

    def _cache_ok(path: Path) -> bool:
        """Return True if cache CSV exists, is non-empty, and has RA/Dec columns."""
        try:
            p = Path(path)
            return p.exists() and p.stat().st_size > 0 and _csv_has_radec(p)
        except Exception:
            return False

    # Header-only sentinel for PS1 mean.csv (keeps "fetch ran / 0 rows" semantics)
    PS1_HEADER = (
        "objID,raMean,decMean,nDetections,ng,nr,ni,nz,ny,"
        "gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag\n"
    )

    def _ensure_ps1_sentinel(ps1_path: Path) -> None:
        """Ensure ps1_neighbourhood.csv exists and is header-only sentinel (non-empty)."""
        try:
            ps1_path.parent.mkdir(parents=True, exist_ok=True)
            # Only write if missing or empty; do not overwrite real data.
            if (not ps1_path.exists()) or ps1_path.stat().st_size == 0:
                ps1_path.write_text(PS1_HEADER, encoding='utf-8')
        except Exception:
            pass

    # Best-effort tile center (used by neighbourhood fetchers)
    try:
        stem = Path(json.loads((run_dir / 'RUN_INDEX.json').read_text(encoding='utf-8'))[0]['tile']).name
        parts = stem.split('_')
        ra_t = float(parts[1])
        dec_t = float(parts[2])
    except Exception:
        ra_t, dec_t = 0.0, 0.0

    # +3' margin beyond circumscribed circle to prevent edge leakage (e.g. source 2319 X=2031)
    radius_arcmin = args.size_arcmin * (2 ** 0.5) * 0.5 + 3.0

    gaia_cache = run_dir / 'catalogs' / 'gaia_neighbourhood.csv'
    ps1_cache  = run_dir / 'catalogs' / 'ps1_neighbourhood.csv'
    usnob_cache = run_dir / 'catalogs' / 'usnob_neighbourhood.csv'

    # -----------------------------
    # Cache-aware fetch-on-miss (with optional NO_FETCH)
    # -----------------------------
    if NO_FETCH:
        print('[STEP4][INFO]', run_dir.name, 'VASCO_STEP4_NO_FETCH=1 -> skipping Gaia/PS1/USNO fetch; using existing caches only')
    else:
        # ---- Gaia ----
        try:
            if os.getenv('VASCO_FORCE_FETCH_GAIA'):
                fetch_gaia_neighbourhood(run_dir, ra_t, dec_t, radius_arcmin)
            elif _cache_ok(gaia_cache):
                print('[STEP4][INFO]', run_dir.name, 'Gaia cache present — skipping fetch')
            else:
                fetch_gaia_neighbourhood(run_dir, ra_t, dec_t, radius_arcmin)
        except Exception as e:
            print('[STEP4][WARN]', run_dir.name, 'Gaia fetch failed:', e)

        # ---- PS1 ----
        try:
            if os.getenv('VASCO_DISABLE_PS1'):
                print('[STEP4][INFO]', run_dir.name, 'PS1 disabled by env')
            else:
                if dec_t < -30.0:
                    print('[STEP4][INFO]', run_dir.name, f'PS1 skipped (Dec={dec_t:.3f} < -30°, outside coverage)')
                    _ensure_ps1_sentinel(ps1_cache)
                else:
                    if os.getenv('VASCO_FORCE_FETCH_PS1'):
                        fetch_ps1_neighbourhood(run_dir, ra_t, dec_t, radius_arcmin)
                    elif _cache_ok(ps1_cache):
                        print('[STEP4][INFO]', run_dir.name, 'PS1 cache present — skipping fetch')
                    else:
                        fetch_ps1_neighbourhood(run_dir, ra_t, dec_t, radius_arcmin)

                    # If PS1 returned 0 bytes, normalize to header-only sentinel
                    try:
                        if ps1_cache.exists() and ps1_cache.stat().st_size == 0:
                            _ensure_ps1_sentinel(ps1_cache)
                    except Exception:
                        pass
        except Exception as e:
            print('[STEP4][WARN]', run_dir.name, 'PS1 fetch failed:', e)

        # ---- USNO-B ----
        try:
            if os.getenv('VASCO_DISABLE_USNOB'):
                print('[STEP4][INFO]', run_dir.name, 'USNO-B disabled by env')
            else:
                if os.getenv('VASCO_FORCE_FETCH_USNOB'):
                    fetch_usnob_neighbourhood(run_dir, ra_t, dec_t, radius_arcmin)
                    print('[STEP4]', run_dir.name, 'USNO-B -> catalogs/usnob_neighbourhood.csv')
                elif _cache_ok(usnob_cache):
                    print('[STEP4][INFO]', run_dir.name, 'USNO-B cache present — skipping fetch')
                else:
                    fetch_usnob_neighbourhood(run_dir, ra_t, dec_t, radius_arcmin)
                    print('[STEP4]', run_dir.name, 'USNO-B -> catalogs/usnob_neighbourhood.csv')
        except Exception as e:
            print('[STEP4][WARN]', run_dir.name, 'USNO-B fetch failed:', e)

    # ------------------------------------------------------------------
    # Bug #6 fix: epoch-propagate Gaia (J2016.0) and USNO-B (J2000.0)
    # catalog positions to the plate epoch using catalogued PMs, so the
    # 5" cone match correctly vetoes high-PM stars at their plate-epoch
    # positions instead of letting them leak through.
    # ------------------------------------------------------------------
    try:
        raw_fits = next(iter(sorted((run_dir / 'raw').glob('*.fits'))), None)
        plate_epoch = _plate_epoch_year_from_fits(raw_fits) if raw_fits else None
        if plate_epoch is None:
            # Post-hoc re-runs work from a harvested catalogue with no tile FITS,
            # so the epoch cannot be read from a header. It is a per-PLATE
            # constant (one plate, one DATE-OBS), so an explicit override is
            # exact rather than approximate -- see the plate epoch lookup built
            # for the epoch delivery.
            #
            # Worth overriding rather than tolerating: without an epoch the run
            # SKIPS proper-motion propagation entirely, and over a ~65-year
            # baseline a 15 mas/yr star moves 1.0", a fifth of the veto radius.
            # The veto then silently fails to remove real moving stars, which is
            # a change in results, not a degradation in convenience.
            env_epoch = os.getenv('VASCO_PLATE_EPOCH_YEAR', '').strip()
            if env_epoch:
                try:
                    plate_epoch = float(env_epoch)
                    print(f'[STEP4][INFO] {run_dir.name} plate epoch from '
                          f'VASCO_PLATE_EPOCH_YEAR={plate_epoch:.3f} (no tile FITS)')
                except ValueError:
                    print(f'[STEP4][WARN] {run_dir.name} VASCO_PLATE_EPOCH_YEAR='
                          f'{env_epoch!r} is not a number; ignoring')
        if plate_epoch is not None:
            print(f'[STEP4][INFO] {run_dir.name} plate epoch = {plate_epoch:.3f}')
            catdir = run_dir / 'catalogs'
            if gaia_cache.exists() and gaia_cache.stat().st_size > 0:
                ok = _propagate_catalog_epoch(
                    gaia_cache, catdir / 'gaia_neighbourhood_at_plate.csv',
                    catalog_epoch_year=2016.0, plate_epoch_year=plate_epoch,
                    ra_col='ra', dec_col='dec', pmra_col='pmRA', pmde_col='pmDE',
                )
                if ok:
                    print(f'[STEP4][INFO] {run_dir.name} Gaia epoch-propagated (dt={plate_epoch-2016.0:+.1f} yr)')
            if usnob_cache.exists() and usnob_cache.stat().st_size > 0:
                ok = _propagate_catalog_epoch(
                    usnob_cache, catdir / 'usnob_neighbourhood_at_plate.csv',
                    catalog_epoch_year=2000.0, plate_epoch_year=plate_epoch,
                    ra_col='ra', dec_col='dec', pmra_col='pmRA', pmde_col='pmDE',
                )
                if ok:
                    print(f'[STEP4][INFO] {run_dir.name} USNO-B epoch-propagated (dt={plate_epoch-2000.0:+.1f} yr)')
        else:
            print(f'[STEP4][WARN] {run_dir.name} plate epoch unavailable — skipping epoch propagation')
    except Exception as e:
        print(f'[STEP4][WARN] {run_dir.name} epoch propagation failed: {e}')

    # Proceed to xmatch using whatever caches exist
    _post_xmatch_tile(run_dir, p2, radius_arcsec=float(args.xmatch_radius_arcsec))
    _update_tile_status(run_dir, 'step4', 'ok')
    return 0

def cmd_step5_filter_within5(args: argparse.Namespace) -> int:
    # within5arcsec CSVs are legacy artifacts from the old vasco pipeline (used for
    # coordinate fixing and reporting). In vasco60 the veto chain uses tskymatch2
    # join=1not2 directly; nothing downstream reads within5arcsec files.
    # Step5 is kept as a no-op so the CLI and parallel runner don't break.
    run_dir = _build_run_dir(Path(args.workdir) if args.workdir else None)
    xdir = run_dir / 'xmatch'
    if not xdir.exists():
        print('[STEP5][ERROR] xmatch/ missing. Run step4-xmatch first.')
        _update_tile_status(run_dir, 'step5', 'fail')
        return 2
    # Clean up any stale within5arcsec files left from previous runs.
    removed = 0
    for stale in sorted(xdir.glob('*_within5arcsec.csv')):
        try:
            stale.unlink()
            removed += 1
        except Exception:
            pass
    if removed:
        print(f'[STEP5] Removed {removed} stale within5arcsec file(s).')
    print('[STEP5] Done (within5arcsec generation removed; veto chain unaffected).')
    _update_tile_status(run_dir, 'step5', 'ok')
    return 0

def cmd_step6_summarize(args: argparse.Namespace) -> int:
    # LEGACY / OPTIONAL — not part of the standard vasco60 pipeline run.
    # Exports raw pass2.ldac (pre-veto, pre-circle-cut) to final_catalog.csv/.ecsv
    # and writes QA plots + RUN_SUMMARY.md.  The production path reads
    # catalogs/sextractor_pass2.filtered.csv (step5 output) via build_run_stage_csvs.py;
    # that script does not consume any output produced here.
    run_dir = _build_run_dir(Path(args.workdir) if args.workdir else None)
    p2 = run_dir / 'pass2.ldac'
    if not p2.exists():
        print('[STEP6][ERROR] pass2.ldac missing. Run step3-psf-and-pass2 first.')
        _update_tile_status(run_dir, 'step6', 'fail')
        return 2
    export_and_summarize(p2, run_dir, export=args.export, histogram_col=args.hist_col)
    _write_text(run_dir / 'RUN_SUMMARY.md', '# Summary written')
    print('[STEP6] Summary + exports written.')
    _update_tile_status(run_dir, 'step6', 'ok')
    return 0

# argparse + main

def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog='vasco.cli_pipeline',
        description='VASCO pipeline orchestrator (split workflow + POSSI-E guard; LDAC)'
    )
    sub = p.add_subparsers(dest='cmd')
    one = sub.add_parser('one2pass', help='One RA/Dec -> 1+2+3 + xmatch + summarize')
    one.add_argument('--ra', type=str, required=True)
    one.add_argument('--dec', type=str, required=True)
    one.add_argument('--size-arcmin', type=float, default=60.0)
    one.add_argument('--survey', default='dss1-red')
    one.add_argument('--export', choices=['none','csv','parquet','both'], default='csv')
    one.add_argument('--hist-col', default='FWHM_IMAGE')
    one.add_argument('--workdir', required=True, help='Tiles root directory; tile folder computed from --ra/--dec')
    one.add_argument('--xmatch-radius-arcsec', type=float, default=5.0)
    one.set_defaults(func=cmd_one)

    s1 = sub.add_parser('step1-download', help='Download tile FITS to raw/ (POSSI-E enforced; header sidecar)')
    s1.add_argument('--ra', type=str, required=True)
    s1.add_argument('--dec', type=str, required=True)
    s1.add_argument('--size-arcmin', type=float, default=60.0)
    s1.add_argument('--survey', default='dss1-red')
    s1.add_argument('--workdir', required=True, help='Tiles root directory; tile folder computed from --ra/--dec')
    s1.set_defaults(func=cmd_step1_download)

    s2 = sub.add_parser('step2-pass1', help='Run SExtractor pass 1')
    s2.add_argument('--workdir', required=True)
    s2.set_defaults(func=cmd_step2_pass1)

    s3 = sub.add_parser('step3-psf-and-pass2', help='Run PSFEx and PSF-aware pass 2')
    s3.add_argument('--workdir', required=True)
    s3.set_defaults(func=cmd_step3_psf_and_pass2)

    s4 = sub.add_parser('step4-xmatch', help='Cross-match (local; PS1 coverage guard)')
    s4.add_argument('--workdir', required=True)
    s4.add_argument('--xmatch-radius-arcsec', type=float, default=5.0)
    s4.add_argument('--size-arcmin', type=float, default=60.0)
    s4.set_defaults(func=cmd_step4_xmatch)

    s5 = sub.add_parser('step5-filter-within5', help='Filter xmatch to <= 5 arcsec')
    s5.add_argument('--workdir', required=True)
    s5.set_defaults(func=cmd_step5_filter_within5)

    s6 = sub.add_parser('step6-summarize', help='Export final CSV/ECSV + QA + RUN_*')
    s6.add_argument('--workdir', required=True)
    s6.add_argument('--export', choices=['none','csv','parquet','both'], default='csv')
    s6.add_argument('--hist-col', default='FWHM_IMAGE')
    s6.set_defaults(func=cmd_step6_summarize)

    args = p.parse_args(argv)
    if hasattr(args, 'func'):
        return args.func(args)
    p.print_help()
    return 0

# step1-download implementation (unchanged structure)

def cmd_step1_download(args: argparse.Namespace) -> int:
    lg = dl.configure_logger(Path('./data/logs'))
    ra = _to_float_ra(args.ra)
    dec = _to_float_dec(args.dec)
    tile_id = format_tile_id(ra, dec)
    run_dir = Path(args.workdir) / tile_id

    try:
        # Downloader will mkdir raw/ only on success (late promotion)
        fits = dl.fetch_skyview_dss(
            ra, dec,
            size_arcmin=args.size_arcmin,
            survey=args.survey,
            pixel_scale_arcsec=_PIXEL_SCALE_ARCSEC,
            out_dir=run_dir / 'raw',   # created only on success
            logger=lg,
        )
        # POSSI‑E enforcement (may delete and raise)
        _enforce_possi_e_or_skip(Path(fits), lg)
        # ---- vasco60 metadata maintenance (no separate post-step scripts) ----
        try:
            data_root = Path('./data')  # repo-root symlink to dataset root
            info = update_all_after_download(
                tile_dir=run_dir,
                fits_path=Path(fits),
                tile_id=run_dir.name,
                ra_deg=float(ra),
                dec_deg=float(dec),
                survey=str(args.survey),
                size_arcmin=float(args.size_arcmin),
                pixel_scale_arcsec=_PIXEL_SCALE_ARCSEC,
                data_root=data_root,
                prefer_local_header=True,
            )
            print('[STEP1][META]', run_dir.name, 'plate_id=', info.get('plate_id',''),
                  '->', Path(info.get('tile_to_plate_csv','')).name)
        except Exception as _e:
            print('[STEP1][WARN] metadata update failed:', _e)

        # ---- SUCCESS PATH ----
        print('[STEP1] Downloaded FITS ->', fits)
        counts = {'planned': 1, 'downloaded': 1, 'processed': 0, 'filtered_non_poss': 0}
        _write_json(run_dir / 'RUN_COUNTS.json', counts)
        _write_json(run_dir / 'RUN_INDEX.json', [{'tile': Path(fits).stem}])
        _write_json(run_dir / 'RUN_MISSING.json', [])
        _write_overview(run_dir, counts, [{'tile': Path(fits).stem}], [])
        _update_tile_status(run_dir, 'step1', 'ok')
        return 0

    except RuntimeError as e:
        # Known non‑POSS case: treated as a skip
        if 'Non-POSS plate returned by STScI' in str(e):
            print('[SKIP]', f'RA={ra:.6f}', f'Dec={dec:.6f}', '-> non-POSS; tile omitted.')
            counts = {'planned': 1, 'downloaded': 0, 'processed': 0, 'filtered_non_poss': 1}
            missing = [{
                'ra': float(ra),
                'dec': float(dec),
                'expected_stem': _expected_stem(ra, dec, args.survey, args.size_arcmin),
            }]
            # Create the tile folder if needed so the skip status actually
            # persists -- _write_json/_update_tile_status don't create
            # parent dirs themselves, so a tile that's never been touched
            # before would otherwise silently record nothing at all (caught
            # 2026-07-19 alongside the REJECT_NON_FITS case below: a fresh
            # tile hitting either skip path left no tile_status.json,
            # so a later step2 attempt would run anyway and fail with its
            # own confusing "no FITS" error instead of a clean, resumable
            # "already checked, no coverage" record).
            run_dir.mkdir(parents=True, exist_ok=True)
            _write_json(run_dir / 'RUN_COUNTS.json', counts)
            _write_json(run_dir / 'RUN_MISSING.json', missing)
            _write_json(run_dir / 'RUN_INDEX.json', [])
            _write_overview(run_dir, counts, [], missing)
            _update_tile_status(run_dir, 'step1', 'skip')
            return 0

        # STScI returned a non-FITS response (e.g. an HTML error page) --
        # observed near the celestial pole, where POSS-I simply has no
        # coverage. Same "no data available" outcome as the non-POSS case
        # above, so treated the same way (skip, not a failure needing
        # investigation). Any *other* RuntimeError falls through below and
        # is a real failure -- previously this whole except block had no
        # trailing `return`, so any unhandled RuntimeError fell off the end
        # and implicitly returned None, which sys.exit(None) turns into
        # exit code 0 -- silently logging "OK" for a download that never
        # happened (caught 2026-07-19 via 7 near-pole tiles that got a
        # correct [REJECT_NON_FITS] error but were still marked OK).
        if 'STScI returned non-FITS' in str(e):
            print('[SKIP]', f'RA={ra:.6f}', f'Dec={dec:.6f}', '-> no POSS-I coverage; tile omitted.')
            counts = {'planned': 1, 'downloaded': 0, 'processed': 0, 'filtered_non_poss': 1}
            missing = [{
                'ra': float(ra),
                'dec': float(dec),
                'expected_stem': _expected_stem(ra, dec, args.survey, args.size_arcmin),
            }]
            run_dir.mkdir(parents=True, exist_ok=True)
            _write_json(run_dir / 'RUN_COUNTS.json', counts)
            _write_json(run_dir / 'RUN_MISSING.json', missing)
            _write_json(run_dir / 'RUN_INDEX.json', [])
            _write_overview(run_dir, counts, [], missing)
            _update_tile_status(run_dir, 'step1', 'skip')
            return 0

        print('[STEP1][ERROR]', str(e))
        return 1

# overview writer unchanged

def _write_overview(run_dir: Path, counts: dict, results: list, missing: list[dict] | None = None) -> None:
    lines = ['# Run Overview','',
             f"**Planned**: {counts.get('planned', 0)}",
             f"**Downloaded**: {counts.get('downloaded', 0)}",
             f"**Processed**: {counts.get('processed', 0)}",
             f"**Non-POSS filtered**: {counts.get('filtered_non_poss', 0)}",
             '']
    if results:
        lines.append('## Tiles (first 10)')
        for rec in results[:10]:
            t = rec.get('tile','?')
            p2 = Path(rec.get('pass2','pass2.ldac')).name
            lines.append(f"- `{t}` → `{p2}`")
        if len(results) > 10:
            lines.append(f"… and {len(results)-10} more tiles.")
        lines.append('')
    if missing:
        lines.append('## Missing tiles (planned but not processed) — first 15')
        for rec in missing[:15]:
            ra = rec.get('ra'); dec = rec.get('dec'); stem = rec.get('expected_stem')
            lines.append(f"- RA={ra:.6f} Dec={dec:.6f} → expected `{stem}`")
        if len(missing) > 15:
            lines.append(f"… and {len(missing)-15} more missing tiles.")
        lines.append('')
    _write_text(run_dir / 'RUN_OVERVIEW.md', ''.join(lines))

if __name__ == '__main__':
    raise SystemExit(main())
