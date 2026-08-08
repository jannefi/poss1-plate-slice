# -*- coding: utf-8 -*-
"""vasco.utils.tile_id

Tile ID helpers for the vasco60 naming convention:

  tile_RA130.013_DECp33.081   (Dec >= 0)
  tile_RA130.013_DECm33.081   (Dec < 0)

Returns (ra, dec) normalized to:
  ra in [0, 360)
  dec in [-90, +90]

Used by multiple modules (pipeline, wcsfix, tiles adapters).
"""

from __future__ import annotations

import re
from typing import Optional, Tuple

_RE_NEW = re.compile(r"^tile_RA(?P<ra>[0-9]+(?:\.[0-9]+)?)_DEC(?P<sign>[pm])(?P<dec>[0-9]+(?:\.[0-9]+)?)$")


def _norm_ra_dec(ra: float, dec: float) -> Tuple[float, float]:
    ra = float(ra) % 360.0
    dec = float(dec)
    if dec > 90.0:
        dec = 90.0
    if dec < -90.0:
        dec = -90.0
    return ra, dec


def parse_tile_id_center(tile_id: str) -> Optional[Tuple[float, float]]:
    """Parse (ra, dec) center from a tile directory name.

    Returns None if the name doesn't match known schemes.
    """
    s = (tile_id or "").strip()

    m = _RE_NEW.match(s)
    if m:
        ra = float(m.group('ra'))
        dec_abs = float(m.group('dec'))
        sign = m.group('sign')
        dec = dec_abs if sign == 'p' else -dec_abs
        return _norm_ra_dec(ra, dec)

    return None


def format_tile_id(ra_deg: float, dec_deg: float, ndp: int = 3) -> str:
    """Format a tile_id in the new vasco60 naming style."""
    ra, dec = _norm_ra_dec(ra_deg, dec_deg)
    sign = 'p' if dec >= 0 else 'm'
    return f"tile_RA{ra:.{ndp}f}_DEC{sign}{abs(dec):.{ndp}f}"
