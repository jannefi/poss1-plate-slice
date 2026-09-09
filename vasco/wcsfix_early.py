from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, List
from vasco.utils.tile_id import parse_tile_id_center

import numpy as np

from vasco.utils.stilts_wrapper import stilts_xmatch


@dataclass
class WcsFixConfig:
    # Bootstrap match radius (arcsec) to build tie points from local Gaia cache
    bootstrap_radius_arcsec: float = 5.0
    # Polynomial degree (2 => 2D quadratic)
    degree: int = 2
    # Minimum tie points required
    min_matches: int = 20
    # Robust fit iterations
    robust_iters: int = 3
    # Residual clipping (arcsec)
    clip_arcsec_min: float = 1.5
    clip_arcsec_sigma: float = 3.0
    clip_arcsec_max: float = 10.0
    # Output filename
    out_name: str = "sextractor_pass2.wcsfix.csv"
    # Status JSON name
    status_name: str = "wcsfix_status.json"

    # --- guard (docs/WCSFIX_GUARD.md) -------------------------------------
    # A degree-2 fit is well constrained where its tie points are and can be
    # badly wrong where they are not; the XE296 corner tile fitted 312 tie
    # points, 155 of them in one 4x4 cell, reported ok, and displaced real
    # stars ~7" so that they cleared the 5" vetoes. Three layers:
    #  1. fit-support gates -> degree-1 fallback / "unsupported"
    #  2. a post-fit self-check on the corrected catalogue
    #  3. `suspect` in the status, which the S0 build quarantines
    guard_enabled: bool = True
    min_tie_points_deg2: int = 500          # fewer tie points -> fit degree 1
    max_sigma_deg2_arcsec: float = 0.4      # worse at degree 2 -> refit degree 1
    max_sigma_deg1_arcsec: float = 0.6      # worse at degree 1 -> unsupported
    concentration_grid: int = 4             # tie-point coverage grid over the tile
    concentration_max_frac: float = 0.40    # > this share in one cell -> degree 1
    sparse_cell_min: int = 10               # a cell with fewer tie points is sparse
    sparse_cells_max: int = 3               # more sparse cells than this -> degree 1
    # self-check: among corrected detections with NO Gaia within the veto
    # radius, the fraction with a Gaia star at [lo, hi) minus the same on
    # positions shifted north by `shift` (chance). Displaced stars show as an
    # excess; healthy tiles sit at ~0 regardless of field density.
    selfcheck_enabled: bool = True
    selfcheck_veto_arcsec: float = 5.0
    selfcheck_lo_arcsec: float = 5.0
    selfcheck_hi_arcsec: float = 10.0
    selfcheck_shift_arcsec: float = 60.0
    selfcheck_min_snr: float = 10.0
    selfcheck_min_unmatched: int = 50
    selfcheck_excess_thr: float = 0.20


def _wrap_deg_pm180(x: np.ndarray) -> np.ndarray:
    """Wrap degrees to [-180, +180)."""
    return (x + 180.0) % 360.0 - 180.0


def _wrap_deg_0_360(x: np.ndarray) -> np.ndarray:
    """Wrap degrees to [0, 360)."""
    return x % 360.0


def _deg_to_arcsec(d: np.ndarray) -> np.ndarray:
    return d * 3600.0


def _sep_arcsec(ra1_deg: np.ndarray, dec1_deg: np.ndarray, ra2_deg: np.ndarray, dec2_deg: np.ndarray) -> np.ndarray:
    """
    Great-circle separation in arcsec (vectorized).
    """
    ra1 = np.radians(ra1_deg)
    dec1 = np.radians(dec1_deg)
    ra2 = np.radians(ra2_deg)
    dec2 = np.radians(dec2_deg)
    s = 2.0 * np.arcsin(
        np.sqrt(
            np.sin((dec2 - dec1) / 2.0) ** 2
            + np.cos(dec1) * np.cos(dec2) * np.sin((ra2 - ra1) / 2.0) ** 2
        )
    )
    return np.degrees(s) * 3600.0


def _pick_sex_radec_cols(header: List[str]) -> Tuple[str, str]:
    cols = set(header)
    # Prefer windowed world coords for stability; raw as fallback
    if "ALPHAWIN_J2000" in cols and "DELTAWIN_J2000" in cols:
        return "ALPHAWIN_J2000", "DELTAWIN_J2000"
    if "ALPHA_J2000" in cols and "DELTA_J2000" in cols:
        return "ALPHA_J2000", "DELTA_J2000"
    if "X_WORLD" in cols and "Y_WORLD" in cols:
        return "X_WORLD", "Y_WORLD"
    raise ValueError("Could not find SExtractor RA/Dec columns in header.")


def _read_csv_header(path: Path) -> List[str]:
    with path.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.reader(f)
        hdr = next(r, [])
    return [h.strip().lstrip("﻿") for h in hdr]


def _bootstrap_match(tile_dir: Path,
                     sex_csv: Path,
                     gaia_csv: Path,
                     ra_col: str,
                     dec_col: str,
                     cfg: WcsFixConfig) -> Path:
    """
    Create a lightweight bootstrap match file between SExtractor CSV and Gaia neighbourhood cache.
    Uses STILTS to avoid pulling big data into memory.
    """
    out = tile_dir / "catalogs" / "_wcsfix_bootstrap_gaia.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    # Join 1and2 with "best" match within bootstrap radius
    stilts_xmatch(
        str(sex_csv),
        str(gaia_csv),
        str(out),
        ra1=ra_col,
        dec1=dec_col,
        ra2="ra",
        dec2="dec",
        radius_arcsec=float(cfg.bootstrap_radius_arcsec),
        join_type="1and2",
        find="best",
        ofmt="csv",
    )
    return out


def _poly_features(dra_deg: np.ndarray, ddec_deg: np.ndarray, degree: int) -> np.ndarray:
    """
    Build 2D polynomial feature matrix up to degree 2 (or 1).
    For degree=2: [1, x, y, x^2, x*y, y^2]
    For degree=1: [1, x, y]
    """
    x = dra_deg
    y = ddec_deg
    if degree <= 1:
        return np.column_stack([np.ones_like(x), x, y])
    # degree 2
    return np.column_stack([np.ones_like(x), x, y, x * x, x * y, y * y])


def _tie_point_coverage(rows: List[dict], cfg: WcsFixConfig) -> dict:
    """How the tie points spread over the tile, from the bootstrap rows' pixel
    coordinates: the share held by the fullest cell of a grid, and how many
    cells are sparse. Falls back to an empty result if no pixel columns."""
    xs, ys = [], []
    for row in rows:
        try:
            xs.append(float(row.get("X_IMAGE", row.get("XWIN_IMAGE"))))
            ys.append(float(row.get("Y_IMAGE", row.get("YWIN_IMAGE"))))
        except (TypeError, ValueError):
            continue
    if len(xs) < 4:
        return {"grid": cfg.concentration_grid, "max_cell_frac": None, "sparse_cells": None}
    x = np.asarray(xs)
    y = np.asarray(ys)
    g = int(cfg.concentration_grid)
    # the tile's own extent; 2118 px is the standard 60' cutout, never smaller
    hi = max(2118.0, float(x.max()), float(y.max()))
    ix = np.clip(((x - 1.0) / hi * g).astype(int), 0, g - 1)
    iy = np.clip(((y - 1.0) / hi * g).astype(int), 0, g - 1)
    counts = np.zeros((g, g), dtype=int)
    np.add.at(counts, (iy, ix), 1)
    return {
        "grid": g,
        "max_cell_frac": float(counts.max() / max(len(x), 1)),
        "sparse_cells": int(np.count_nonzero(counts < cfg.sparse_cell_min)),
        "cell_counts": counts.tolist(),
    }


def _unit_xyz(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _self_check(out_csv: Path, gaia_csv: Path, ra_col_out: str, dec_col_out: str,
                cfg: WcsFixConfig) -> dict:
    """Post-fit self-check on the corrected catalogue (docs/WCSFIX_GUARD.md).

    Among detections with no Gaia counterpart inside the veto radius -- the
    ones that will survive the Gaia veto -- measure the fraction with a Gaia
    star at [lo, hi). Do the same on positions shifted north by `shift` for
    the chance level. A displaced-star tile shows a large positive excess
    (XE296's corner tile: +47 points; XE084's: +36 in a dense field); healthy
    tiles sit within a few points of zero.
    """
    try:
        from scipy.spatial import cKDTree
    except Exception as e:  # pragma: no cover - environment without scipy
        return {"enabled": False, "reason": f"scipy unavailable: {e}"}

    ra, de = [], []
    with out_csv.open(newline="", encoding="utf-8", errors="ignore") as f:
        r = csv.DictReader(f)
        fields = set(r.fieldnames or [])
        has_snr, has_flags = "SNR_WIN" in fields, "FLAGS" in fields
        for row in r:
            try:
                a = float(row[ra_col_out])
                d = float(row[dec_col_out])
            except (TypeError, ValueError, KeyError):
                continue
            if has_snr:
                try:
                    if not float(row["SNR_WIN"]) > cfg.selfcheck_min_snr:
                        continue
                except (TypeError, ValueError):
                    continue
            if has_flags:
                try:
                    if int(float(row["FLAGS"])) != 0:
                        continue
                except (TypeError, ValueError):
                    pass
            ra.append(a)
            de.append(d)
    gra, gde = [], []
    with gaia_csv.open(newline="", encoding="utf-8", errors="ignore") as f:
        for row in csv.DictReader(f):
            try:
                gra.append(float(row["ra"]))
                gde.append(float(row["dec"]))
            except (TypeError, ValueError, KeyError):
                continue
    n = len(ra)
    if n == 0 or len(gra) == 0:
        return {"enabled": True, "n_detections": n, "n_gaia": len(gra), "suspect": False,
                "reason": "nothing to check"}

    tree = cKDTree(_unit_xyz(np.asarray(gra), np.asarray(gde)))
    ra_a, de_a = np.asarray(ra), np.asarray(de)
    ub = 2.0 * np.sin(np.radians(cfg.selfcheck_hi_arcsec / 3600.0) / 2.0)

    def nearest_arcsec(dec_arr):
        d, _ = tree.query(_unit_xyz(ra_a, dec_arr), k=1, distance_upper_bound=ub)
        sep = np.full(n, np.inf)
        ok = np.isfinite(d)
        sep[ok] = 2.0 * np.degrees(np.arcsin(np.minimum(d[ok] / 2.0, 1.0))) * 3600.0
        return sep

    out = {"enabled": True, "n_detections": n, "n_gaia": len(gra),
           "veto_arcsec": cfg.selfcheck_veto_arcsec,
           "window_arcsec": [cfg.selfcheck_lo_arcsec, cfg.selfcheck_hi_arcsec]}
    for lab, dec_arr in (("", de_a), ("_chance", de_a + cfg.selfcheck_shift_arcsec / 3600.0)):
        sep = nearest_arcsec(dec_arr)
        unmatched = ~(sep < cfg.selfcheck_veto_arcsec)
        k = int(unmatched.sum())
        in_win = (sep[unmatched] >= cfg.selfcheck_lo_arcsec) & (sep[unmatched] < cfg.selfcheck_hi_arcsec)
        out[f"n_unmatched{lab}"] = k
        out[f"frac_window{lab}"] = float(in_win.mean()) if k else None
    if out["frac_window"] is None or out["frac_window_chance"] is None:
        out.update({"excess": None, "suspect": False, "reason": "too few unmatched detections"})
        return out
    out["excess"] = float(out["frac_window"] - out["frac_window_chance"])
    out["suspect"] = bool(out["n_unmatched"] >= cfg.selfcheck_min_unmatched
                          and out["excess"] > cfg.selfcheck_excess_thr)
    return out


def _robust_fit_offsets(dra_det: np.ndarray,
                        ddec_det: np.ndarray,
                        dra_off: np.ndarray,
                        ddec_off: np.ndarray,
                        degree: int,
                        cfg: WcsFixConfig) -> Tuple[np.ndarray, np.ndarray, dict]:
    """
    Robust least squares fit of offsets (Gaia - detection) as a polynomial in (dra_det, ddec_det).
    Returns coefficients for RA offset and Dec offset.
    """
    X = _poly_features(dra_det, ddec_det, degree=degree)

    # Start with all rows valid
    mask = np.isfinite(dra_off) & np.isfinite(ddec_off) & np.all(np.isfinite(X), axis=1)

    info = {"iters": 0, "kept": int(mask.sum()), "dropped": int((~mask).sum())}

    if mask.sum() < cfg.min_matches:
        raise RuntimeError(f"too few usable matches after NaN filter ({mask.sum()} < {cfg.min_matches})")

    for it in range(cfg.robust_iters):
        Xm = X[mask]
        y_ra = dra_off[mask]
        y_de = ddec_off[mask]

        # Fit via least squares
        coef_ra, *_ = np.linalg.lstsq(Xm, y_ra, rcond=None)
        coef_de, *_ = np.linalg.lstsq(Xm, y_de, rcond=None)

        # Compute residuals in arcsec
        pred_ra = X @ coef_ra
        pred_de = X @ coef_de
        res_arcsec = _deg_to_arcsec(np.sqrt((dra_off - pred_ra) ** 2 + (ddec_off - pred_de) ** 2))

        # Robust sigma estimate (MAD)
        r = res_arcsec[np.isfinite(res_arcsec)]
        if r.size == 0:
            break
        med = np.median(r)
        mad = np.median(np.abs(r - med))
        sigma = 1.4826 * mad if mad > 0 else (np.std(r) if r.size > 10 else 0.0)

        # Clip threshold
        thr = max(cfg.clip_arcsec_min, cfg.clip_arcsec_sigma * sigma) if sigma > 0 else cfg.clip_arcsec_max
        thr = min(thr, cfg.clip_arcsec_max)

        new_mask = mask & (res_arcsec <= thr)
        info.update({"iters": it + 1, "sigma_arcsec": float(sigma), "clip_thr_arcsec": float(thr),
                     "kept": int(new_mask.sum()), "dropped": int((~new_mask).sum())})

        # Stop if stable
        if new_mask.sum() == mask.sum():
            mask = new_mask
            break

        mask = new_mask

        if mask.sum() < cfg.min_matches:
            raise RuntimeError(f"too few matches after robust clipping ({mask.sum()} < {cfg.min_matches})")

    # Final fit with final mask
    Xm = X[mask]
    y_ra = dra_off[mask]
    y_de = ddec_off[mask]
    coef_ra, *_ = np.linalg.lstsq(Xm, y_ra, rcond=None)
    coef_de, *_ = np.linalg.lstsq(Xm, y_de, rcond=None)

    return coef_ra, coef_de, info


def ensure_wcsfix_catalog(tile_dir: Path,
                          sex_csv: Path,
                          gaia_csv: Path,
                          *,
                          center: Optional[Tuple[float, float]] = None,
                          cfg: Optional[WcsFixConfig] = None,
                          force: bool = False) -> Tuple[Path, dict]:
    """
    Ensure a WCSFIX-augmented SExtractor catalog exists for this tile.

    Inputs:
      - sex_csv: catalogs/sextractor_pass2.csv (big)
      - gaia_csv: catalogs/gaia_neighbourhood.csv (cache)

    Output:
      - catalogs/sextractor_pass2.wcsfix.csv containing all original columns + RA_corr,Dec_corr
      - catalogs/wcsfix_status.json describing success/failure and fit diagnostics

    Returns: (path_used_for_downstream, status_dict)
    """
    tile_dir = Path(tile_dir)
    sex_csv = Path(sex_csv)
    gaia_csv = Path(gaia_csv)
    cfg = cfg or WcsFixConfig()

    out_csv = tile_dir / "catalogs" / cfg.out_name
    status_path = tile_dir / "catalogs" / cfg.status_name
    status: dict = {
        "ok": False,
        "reason": None,
        "sex_csv": str(sex_csv),
        "gaia_csv": str(gaia_csv),
        "out_csv": str(out_csv),
        "bootstrap_radius_arcsec": cfg.bootstrap_radius_arcsec,
        "degree": cfg.degree,
        "min_matches": cfg.min_matches,
    }

    # If already exists and not forcing, use it
    if out_csv.exists() and out_csv.stat().st_size > 0 and not force:
        status.update({"ok": True, "reason": "cached"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return out_csv, status

    # Preconditions
    if not sex_csv.exists() or sex_csv.stat().st_size == 0:
        status.update({"ok": False, "reason": "missing sextractor_pass2.csv"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    if not gaia_csv.exists() or gaia_csv.stat().st_size == 0:
        status.update({"ok": False, "reason": "missing gaia_neighbourhood.csv"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    # Determine center if not provided
    if center is None:
        # Parse from tile folder name (supports both legacy and vasco60 naming)
        center = parse_tile_id_center(tile_dir.name)
    if center is None:
        status.update({"ok": False, "reason": "missing tile center"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    ra0, dec0 = float(center[0]), float(center[1])
    status.update({"tile_center_ra": ra0, "tile_center_dec": dec0})

    # Determine sextractor coordinate columns
    try:
        hdr = _read_csv_header(sex_csv)
        ra_col, dec_col = _pick_sex_radec_cols(hdr)
    except Exception as e:
        status.update({"ok": False, "reason": f"cannot pick sextractor radec cols: {e}"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    status.update({"sex_ra_col": ra_col, "sex_dec_col": dec_col})

    # Bootstrap match
    try:
        boot = _bootstrap_match(tile_dir, sex_csv, gaia_csv, ra_col, dec_col, cfg)
    except Exception as e:
        status.update({"ok": False, "reason": f"bootstrap match failed: {e}"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    # Read tie points
    try:
        with boot.open(newline="", encoding="utf-8", errors="ignore") as f:
            r = csv.DictReader(f)
            rows = list(r)
    except Exception as e:
        status.update({"ok": False, "reason": f"cannot read bootstrap csv: {e}"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    if not rows:
        status.update({"ok": False, "reason": "bootstrap match produced 0 rows"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    # Extract numeric arrays
    det_ra = []
    det_de = []
    g_ra = []
    g_de = []

    for row in rows:
        try:
            dra = float(row[ra_col])
            dde = float(row[dec_col])
            gra = float(row["ra"])
            gde = float(row["dec"])
        except Exception:
            continue
        det_ra.append(dra)
        det_de.append(dde)
        g_ra.append(gra)
        g_de.append(gde)

    det_ra = np.asarray(det_ra, dtype=float)
    det_de = np.asarray(det_de, dtype=float)
    g_ra = np.asarray(g_ra, dtype=float)
    g_de = np.asarray(g_de, dtype=float)

    n = det_ra.size
    status.update({"bootstrap_rows": int(len(rows)), "tie_points": int(n)})

    if n < cfg.min_matches:
        status.update({"ok": False, "reason": f"too few tie points ({n} < {cfg.min_matches})"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    # Work in local (wrapped) coordinates around tile center
    dra_det = _wrap_deg_pm180(det_ra - ra0)
    ddec_det = (det_de - dec0)

    # Offsets to Gaia (wrapped small RA offset)
    dra_off = _wrap_deg_pm180(g_ra - det_ra)
    ddec_off = (g_de - det_de)

    # --- guard, layer 1: is a degree-2 fit supported by these tie points? ---
    guard: dict = {"enabled": bool(cfg.guard_enabled), "degree_requested": int(cfg.degree),
                   "gate_reasons": [], "unsupported": False, "suspect": False}
    degree_used = int(cfg.degree)
    coverage = _tie_point_coverage(rows, cfg)
    guard["coverage"] = coverage
    if cfg.guard_enabled and degree_used > 1:
        if n < cfg.min_tie_points_deg2:
            guard["gate_reasons"].append(f"tie points {n} < {cfg.min_tie_points_deg2}")
        if coverage.get("max_cell_frac") is not None and coverage["max_cell_frac"] > cfg.concentration_max_frac:
            guard["gate_reasons"].append(
                f"{coverage['max_cell_frac']:.2f} of tie points in one {coverage['grid']}x{coverage['grid']} cell")
        if coverage.get("sparse_cells") is not None and coverage["sparse_cells"] > cfg.sparse_cells_max:
            guard["gate_reasons"].append(
                f"{coverage['sparse_cells']} cells with < {cfg.sparse_cell_min} tie points")
        if guard["gate_reasons"]:
            degree_used = 1

    # Robust fit (degree-1 fallback if the degree-2 fit is unsupported or too loose)
    try:
        coef_ra, coef_de, fit_info = _robust_fit_offsets(
            dra_det, ddec_det, dra_off, ddec_off,
            degree=degree_used, cfg=cfg
        )
        if (cfg.guard_enabled and degree_used > 1
                and fit_info.get("sigma_arcsec") is not None
                and fit_info["sigma_arcsec"] > cfg.max_sigma_deg2_arcsec):
            guard["gate_reasons"].append(
                f"degree-2 sigma {fit_info['sigma_arcsec']:.3f}\" > {cfg.max_sigma_deg2_arcsec}\"")
            degree_used = 1
            coef_ra, coef_de, fit_info = _robust_fit_offsets(
                dra_det, ddec_det, dra_off, ddec_off,
                degree=degree_used, cfg=cfg
            )
    except Exception as e:
        status.update({"ok": False, "reason": f"fit failed: {e}", "guard": guard})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status

    guard["degree_used"] = degree_used
    if (cfg.guard_enabled and fit_info.get("sigma_arcsec") is not None
            and fit_info["sigma_arcsec"] > cfg.max_sigma_deg1_arcsec):
        # Raw plate coordinates are not a safe fallback (on XE296's corner they
        # are what failed), so the best available fit is still applied and the
        # tile is marked for quarantine instead.
        guard["unsupported"] = True
        guard["unsupported_reason"] = (
            f"degree-{degree_used} sigma {fit_info['sigma_arcsec']:.3f}\" > {cfg.max_sigma_deg1_arcsec}\"")

    status.update({
        "fit": fit_info,
        "coef_ra": [float(x) for x in coef_ra.tolist()],
        "coef_de": [float(x) for x in coef_de.tolist()],
        "guard": guard,
    })

    # Apply to full SExtractor catalog and write output
    try:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        tmp = out_csv.with_suffix(out_csv.suffix + ".tmp")

        with sex_csv.open(newline="", encoding="utf-8", errors="ignore") as fin, tmp.open("w", newline="", encoding="utf-8") as fout:
            rdr = csv.DictReader(fin)
            fieldnames = list(rdr.fieldnames or [])
            # Avoid duplicate columns if re-running
            for extra in ("RA_corr", "Dec_corr"):
                if extra not in fieldnames:
                    fieldnames.append(extra)

            w = csv.DictWriter(fout, fieldnames=fieldnames)
            w.writeheader()

            for row in rdr:
                try:
                    ra_det_row = float(row.get(ra_col, "nan"))
                    dec_det_row = float(row.get(dec_col, "nan"))
                except Exception:
                    ra_det_row = float("nan")
                    dec_det_row = float("nan")

                if not (math.isfinite(ra_det_row) and math.isfinite(dec_det_row)):
                    row["RA_corr"] = ""
                    row["Dec_corr"] = ""
                    w.writerow(row)
                    continue

                dra = _wrap_deg_pm180(np.array([ra_det_row - ra0], dtype=float))[0]
                dde = (dec_det_row - dec0)

                X = _poly_features(np.array([dra], dtype=float), np.array([dde], dtype=float), degree=degree_used)
                off_ra = (X @ coef_ra).item()
                off_de = (X @ coef_de).item()


                ra_corr = _wrap_deg_0_360(ra_det_row + off_ra)
                dec_corr = dec_det_row + off_de

                row["RA_corr"] = f"{ra_corr:.10f}"
                row["Dec_corr"] = f"{dec_corr:.10f}"
                w.writerow(row)

        tmp.replace(out_csv)

        # --- guard, layer 2: does the corrected catalogue look displaced? ---
        if cfg.guard_enabled and cfg.selfcheck_enabled:
            try:
                guard["selfcheck"] = _self_check(out_csv, gaia_csv, "RA_corr", "Dec_corr", cfg)
            except Exception as e:  # never let the check itself break the step
                guard["selfcheck"] = {"enabled": True, "error": f"{type(e).__name__}: {e}", "suspect": False}
            guard["suspect"] = bool(guard["unsupported"] or guard["selfcheck"].get("suspect"))
        else:
            guard["suspect"] = bool(guard["unsupported"])
        status["guard"] = guard

        status.update({"ok": True, "reason": "wrote", "out_rows": "unknown"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass

        return out_csv, status

    except Exception as e:
        status.update({"ok": False, "reason": f"write failed: {e}"})
        try:
            status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        except Exception:
            pass
        return sex_csv, status