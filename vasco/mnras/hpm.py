from __future__ import annotations
from typing import Tuple, Dict, Any
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time


def _pm(row: Dict[str, Any], *keys: str) -> float:
    """Return the first present, parseable PM value among `keys`, else 0.0.

    Gaia neighbourhood CSVs (VizieR + local cache) write `pmRA`/`pmDE`;
    empty string means NaN. Older callers may pass `pmra`/`pmdec`.
    """
    for k in keys:
        v = row.get(k)
        if v in (None, ""):
            continue
        try:
            return float(v)
        except (TypeError, ValueError):
            continue
    return 0.0


def backprop_gaia_row(row: Dict[str, Any], target_epoch: float = 1950.0) -> Tuple[float, float]:
    """Back-propagate a Gaia row (ra, dec, pmRA/pmra, pmDE/pmdec) to target_epoch (Julian year)."""
    try:
        pm_ra_cosdec = _pm(row, "pmRA", "pmra")
        pm_dec = _pm(row, "pmDE", "pmdec")
        c = SkyCoord(ra=float(row['ra']) * u.deg, dec=float(row['dec']) * u.deg,
                     pm_ra_cosdec=pm_ra_cosdec * u.mas / u.yr,
                     pm_dec=pm_dec * u.mas / u.yr,
                     obstime=Time(2016.0, format='jyear'))
        c2 = c.apply_space_motion(Time(target_epoch, format='jyear'))
        return c2.ra.deg, c2.dec.deg
    except Exception:
        return float('nan'), float('nan')


# How to validate:
#   With a real Gaia row {ra: 85.2658575289, dec: 28.78976672882,
#                         pmRA: -59.53, pmDE: -47.884} and target 1950.0,
#   the returned (ra_bp, dec_bp) must differ from (ra, dec) by ~5.04″
#   (66 yr × ~75 mas/yr total). Pre-fix, the shift was 0.000″.
