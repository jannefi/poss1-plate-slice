"""
Build the USNO-B1.0 10-column HEALPix-5 Parquet cache from VizieR TAP.

Source : https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync  (ADQL)
Table  : "I/284/out" (Monet+ 2003; 1,045,913,669 rows)

Strategy: paginate the catalog in DEC zones (0.1 deg wide, 1800 zones),
submit a synchronous ADQL query per zone projecting only the 10 columns
VASCO60 needs, parse the CSV response, compute healpix_5 from (ra,dec),
write one or more Parquet files to the local HP5-partitioned cache.

Columns kept (matching fetch_usnob_neighbourhood in
vasco/external_fetch_usnob_vizier.py):
    USNO-B1_0 (id, string),
    RAJ2000->ra (float64 deg),
    DEJ2000->dec (float64 deg),
    B1mag, R1mag, B2mag, R2mag, Imag (float32 mag),
    pmRA, pmDE (int16 mas/yr)

Resumable via per-zone .done markers.

Defaults:
    --cache-dir / $VASCO_USNOB_CACHE / /Volumes/SANDISK/USNOB
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
import threading
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.dataset as pds
from astropy_healpix import HEALPix
import astropy.units as u

TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
TAP_TABLE = '"I/284/out"'

# ADQL column list (in order). The id column's canonical name is
# "USNO-B1.0" — VizieR's CSV output renders it as "USNO-B1_0" but the
# ADQL parser requires the dotted form in double quotes. Aliasing to
# `id` also breaks the parser (it splits the quoted name on the dash),
# so we don't alias; the CSV header comes back as "USNO-B1_0" and we
# rename it at parse time.
COL_SOURCE = [
    '"USNO-B1.0"',
    "RAJ2000",
    "DEJ2000",
    "B1mag",
    "R1mag",
    "B2mag",
    "R2mag",
    "Imag",
    "pmRA",
    "pmDE",
]
# Column name in the returned CSV header for the id column
CSV_ID_COL = "USNO-B1_0"

# Output schema for the Parquet cache
OUTPUT_SCHEMA = pa.schema(
    [
        ("id", pa.string()),
        ("ra", pa.float64()),
        ("dec", pa.float64()),
        ("B1mag", pa.float32()),
        ("R1mag", pa.float32()),
        ("B2mag", pa.float32()),
        ("R2mag", pa.float32()),
        ("Imag", pa.float32()),
        ("pmRA", pa.int16()),
        ("pmDE", pa.int16()),
        ("healpix_5", pa.int32()),
    ]
)

PUBLISHED_ROWS = 1_045_913_669  # from I/284 ReadMe
ZONE_SIZE_DEG = 0.1  # 1800 zones spanning [-90, +90]
N_ZONES = 1800

_log_lock = threading.Lock()
_write_lock = threading.Lock()
_log_path: Path | None = None

_HP5 = HEALPix(nside=2**5, order="nested")


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_USNOB_CACHE", "/Volumes/SANDISK/USNOB"))


def log(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    with _log_lock:
        print(line, flush=True)
        if _log_path is not None:
            try:
                with open(_log_path, "a") as f:
                    f.write(line + "\n")
            except OSError:
                pass


def zone_bounds(idx: int) -> tuple[float, float]:
    """dec bounds [lo, hi) for zone idx, where idx=0 is south pole zone."""
    lo = -90.0 + idx * ZONE_SIZE_DEG
    hi = lo + ZONE_SIZE_DEG
    return lo, hi


def tap_query_zone(idx: int, attempts: int = 5) -> bytes:
    lo, hi = zone_bounds(idx)
    # Use < hi for all but the north-polar zone, which uses <= hi to
    # include dec == +90 sources (if any).
    op = "<=" if idx == N_ZONES - 1 else "<"
    adql = (
        f"SELECT {','.join(COL_SOURCE)} "
        f"FROM {TAP_TABLE} "
        f"WHERE DEJ2000 >= {lo} AND DEJ2000 {op} {hi}"
    )
    params = {
        "REQUEST": "doQuery",
        "LANG": "ADQL",
        "FORMAT": "csv",
        "MAXREC": "10000000",
        "QUERY": adql,
    }
    url = TAP_URL + "?" + urllib.parse.urlencode(params)
    last_err: Exception | None = None
    for k in range(attempts):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "vasco60-usnob-build/1.0"}
            )
            with urllib.request.urlopen(req, timeout=600) as r:
                data = r.read()
            # TAP errors come back with HTTP 200 but a VOTable error payload.
            # CSV format wraps errors differently — detect by checking for an
            # error string in the header of the response.
            head = data[:4096].decode("utf-8", errors="replace").lower()
            if "error" in head and ("vot" in head or "<votable" in head or head.startswith("<?xml")):
                raise RuntimeError(f"TAP returned error for zone {idx}: {head[:500]}")
            return data
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_err = e
            log(f"TAP retry {k+1}/{attempts} zone {idx}: {e}")
            time.sleep(2 ** k)
        except RuntimeError as e:
            last_err = e
            log(f"TAP error {k+1}/{attempts} zone {idx}: {e}")
            time.sleep(2 ** k)
    raise RuntimeError(f"zone {idx} failed after {attempts} attempts: {last_err}")


def parse_csv(data: bytes) -> pa.Table:
    """Parse VizieR TAP CSV with explicit types + null handling."""
    # Nulls in VizieR CSV are empty fields.
    convert = pcsv.ConvertOptions(
        column_types={
            CSV_ID_COL: pa.string(),
            "RAJ2000": pa.float64(),
            "DEJ2000": pa.float64(),
            "B1mag": pa.float32(),
            "R1mag": pa.float32(),
            "B2mag": pa.float32(),
            "R2mag": pa.float32(),
            "Imag": pa.float32(),
            "pmRA": pa.int16(),
            "pmDE": pa.int16(),
        },
        null_values=["", "NaN", "nan", "NULL", "null", "--"],
        strings_can_be_null=True,
    )
    return pcsv.read_csv(io.BytesIO(data), convert_options=convert)


def transform_and_write(tbl: pa.Table, out_parquet: Path, zone_idx: int) -> tuple[int, int]:
    ra = tbl.column("RAJ2000").to_numpy(zero_copy_only=False)
    dec = tbl.column("DEJ2000").to_numpy(zero_copy_only=False)
    valid = ~(np.isnan(ra) | np.isnan(dec))
    if not valid.all():
        tbl = tbl.filter(pa.array(valid))
        ra = ra[valid]
        dec = dec[valid]
    if len(ra) == 0:
        return 0, 0
    hp5 = _HP5.lonlat_to_healpix(ra * u.deg, dec * u.deg).astype(np.int32)

    out = pa.table(
        {
            "id": tbl.column(CSV_ID_COL),
            "ra": tbl.column("RAJ2000"),
            "dec": tbl.column("DEJ2000"),
            "B1mag": tbl.column("B1mag"),
            "R1mag": tbl.column("R1mag"),
            "B2mag": tbl.column("B2mag"),
            "R2mag": tbl.column("R2mag"),
            "Imag": tbl.column("Imag"),
            "pmRA": tbl.column("pmRA"),
            "pmDE": tbl.column("pmDE"),
            "healpix_5": pa.array(hp5, type=pa.int32()),
        }
    )

    partitioning = pds.partitioning(pa.schema([("healpix_5", pa.int32())]), flavor="hive")
    basename = f"USNOB-zone{zone_idx:04d}"
    with _write_lock:
        pds.write_dataset(
            out,
            base_dir=str(out_parquet),
            partitioning=partitioning,
            format="parquet",
            basename_template=f"{basename}-part-{{i}}.parquet",
            existing_data_behavior="overwrite_or_ignore",
            file_options=pds.ParquetFileFormat().make_write_options(
                compression="zstd", compression_level=3
            ),
        )
    n_parts = int(np.unique(hp5).size)
    return out.num_rows, n_parts


def process_zone(idx: int, out_parquet: Path, done_dir: Path) -> dict:
    t0 = time.time()
    marker = done_dir / f"zone_{idx:04d}.done"
    if marker.exists():
        return {"idx": idx, "skipped": True}
    try:
        data = tap_query_zone(idx)
        tbl = parse_csv(data)
        n_rows, n_parts = transform_and_write(tbl, out_parquet, idx)
        lo, hi = zone_bounds(idx)
        marker.write_text(f"{lo:.4f} {hi:.4f} {n_rows} {n_parts}\n")
        return {"idx": idx, "rows": n_rows, "parts": n_parts, "sec": time.time() - t0}
    except Exception as e:
        log(f"FAIL zone {idx}: {e}")
        traceback.print_exc()
        return {"idx": idx, "error": str(e)}


def main() -> int:
    global _log_path
    parser = argparse.ArgumentParser(
        description="Build the USNO-B1.0 10-column HEALPix-5 Parquet cache from VizieR TAP."
    )
    parser.add_argument("--cache-dir", type=Path, default=default_cache_dir())
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, default=0, help="process only N zones (0 = all)")
    parser.add_argument("--start", type=int, default=0, help="start at zone index")
    parser.add_argument("--end", type=int, default=N_ZONES, help="end before zone index (exclusive)")
    args = parser.parse_args()

    cache: Path = args.cache_dir
    parquet_dir = cache / "parquet"
    done_dir = cache / ".done"
    _log_path = cache / "progress.log"
    for d in (parquet_dir, done_dir):
        d.mkdir(parents=True, exist_ok=True)

    log(f"cache_dir={cache}")
    log(f"TAP endpoint: {TAP_URL}")

    zones = list(range(args.start, args.end))
    if args.limit:
        zones = zones[: args.limit]

    pending = [i for i in zones if not (done_dir / f"zone_{i:04d}.done").exists()]
    log(f"{len(pending)}/{len(zones)} zones to process")

    total_rows = 0
    n_ok = 0
    n_fail = 0
    t_start = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(process_zone, i, parquet_dir, done_dir): i for i in pending}
        for k, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            if "error" in res:
                n_fail += 1
                continue
            if res.get("skipped"):
                continue
            n_ok += 1
            total_rows += res["rows"]
            if k % 20 == 0 or k == len(pending):
                elapsed = time.time() - t_start
                rate_per_hr = k / elapsed * 3600 if elapsed > 0 else 0
                eta_s = (len(pending) - k) / (k / elapsed) if k > 0 and elapsed > 0 else 0
                log(
                    f"{k}/{len(pending)} rows={total_rows:,} "
                    f"rate={rate_per_hr:.0f}/hr eta={eta_s/60:.0f}min"
                )
    dt = time.time() - t_start
    log(f"DONE ok={n_ok} fail={n_fail} rows={total_rows:,} elapsed={dt:.0f}s")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
