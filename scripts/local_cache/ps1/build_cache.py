"""
Build the PS1 DR2 9-column HEALPix-5 Parquet cache from STScI's AWS
HATS dataset.

Source : s3://stpubdata/panstarrs/ps1/public/hats/otmo/otmo/dataset/
         (anonymous read; STScI public dataset mirror)

Strategy: pyarrow.fs.S3FileSystem with Parquet column projection.
For each HATS source parquet file, we read ONLY the 9 columns we
need (~9% of the source bytes) — never download the full 135-column
file. No local staging.

Columns kept (VASCO60 PS1 query set):
    objID, raMean->ra, decMean->dec, nDetections,
    gMeanPSFMag->gmag, rMeanPSFMag->rmag, iMeanPSFMag->imag,
    zMeanPSFMag->zmag, yMeanPSFMag->ymag
    + derived healpix_5 (computed from ra/dec via astropy_healpix
    since PS1 objID does NOT encode HEALPix position)

Output layout: matches scripts/local_cache/gaia/ — Hive-partitioned
Parquet with `healpix_5=<N>/` directories, level 5 nested (12,288 pixels).

Resumable via per-source-file .done markers.

Defaults (in priority order):
    --cache-dir / $VASCO_PS1_CACHE / /Volumes/SANDISK/PS1

Precision policy: float64 ra/dec (matching Gaia cache), float32 for the
5 PSF mags (float32 ULP << PS1 photometric precision). See README.md
for the full rationale.
"""

from __future__ import annotations

import argparse
import os
import sys
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.dataset as pds
import pyarrow.fs as pafs
import pyarrow.parquet as pq
from astropy_healpix import HEALPix
import astropy.units as u

S3_BUCKET = "stpubdata"
S3_PREFIX = "panstarrs/ps1/public/hats/otmo/otmo/dataset/"
S3_REGION = "us-east-1"

PUBLISHED_HATS_ROWS = 10_560_724_292  # from otmo/properties hats_nrows

SOURCE_COLS = [
    "objID",
    "raMean",
    "decMean",
    "nDetections",
    "gMeanPSFMag",
    "rMeanPSFMag",
    "iMeanPSFMag",
    "zMeanPSFMag",
    "yMeanPSFMag",
]

OUTPUT_SCHEMA = pa.schema(
    [
        ("objID", pa.int64()),
        ("ra", pa.float64()),
        ("dec", pa.float64()),
        ("nDetections", pa.int16()),
        ("gmag", pa.float32()),
        ("rmag", pa.float32()),
        ("imag", pa.float32()),
        ("zmag", pa.float32()),
        ("ymag", pa.float32()),
        ("healpix_5", pa.int32()),
    ]
)

_log_lock = threading.Lock()
_write_lock = threading.Lock()
_log_path: Path | None = None

_HP5 = HEALPix(nside=2**5, order="nested")


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_PS1_CACHE", "/Volumes/SANDISK/PS1"))


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


def list_s3_parquet(fs: pafs.S3FileSystem) -> list[tuple[str, int]]:
    """Walk the HATS dataset prefix and return [(key, size), ...] of all .parquet leaves."""
    out: list[tuple[str, int]] = []
    selector = pafs.FileSelector(f"{S3_BUCKET}/{S3_PREFIX}", recursive=True)
    for info in fs.get_file_info(selector):
        if info.is_file and info.path.endswith(".parquet") and not info.path.endswith(
            "/_metadata"
        ) and not info.path.endswith("/_common_metadata"):
            # info.path is "stpubdata/panstarrs/ps1/public/hats/otmo/otmo/dataset/..."
            out.append((info.path, info.size))
    return out


def _marker_name(s3_path: str) -> str:
    """Make a filesystem-safe .done marker name from the S3 path."""
    tail = s3_path.split(S3_PREFIX, 1)[1]  # "Norder=7/Dir=110000/Npix=115100.parquet"
    return tail.replace("/", "_").replace(".parquet", ".done")


def _basename(s3_path: str) -> str:
    """Unique basename for output parquet files from this source file."""
    tail = s3_path.split(S3_PREFIX, 1)[1]
    # "Norder=7/Dir=110000/Npix=115100.parquet" -> "PS1-N7-D110000-P115100"
    parts = tail.replace(".parquet", "").split("/")
    norder = parts[0].split("=")[1]
    d = parts[1].split("=")[1] if len(parts) > 1 else "0"
    npix = parts[2].split("=")[1] if len(parts) > 2 else "0"
    return f"PS1-N{norder}-D{d}-P{npix}"


def transform_one(
    s3: pafs.S3FileSystem, s3_path: str, out_parquet: Path
) -> tuple[int, int]:
    """Stream 9 columns from one HATS parquet file to local HP5 partitioned output."""
    with s3.open_input_file(s3_path) as f:
        tbl = pq.read_table(f, columns=SOURCE_COLS)

    ra = tbl.column("raMean").to_numpy(zero_copy_only=False)
    dec = tbl.column("decMean").to_numpy(zero_copy_only=False)
    # HATS data always has valid ra/dec; if any NaN, exclude
    valid = ~(np.isnan(ra) | np.isnan(dec))
    if not valid.all():
        tbl = tbl.filter(pa.array(valid))
        ra = ra[valid]
        dec = dec[valid]

    hp5 = _HP5.lonlat_to_healpix(ra * u.deg, dec * u.deg).astype(np.int32)

    # Rename + cast in one shot
    out_tbl = pa.table(
        {
            "objID": tbl.column("objID").cast(pa.int64()),
            "ra": tbl.column("raMean").cast(pa.float64()),
            "dec": tbl.column("decMean").cast(pa.float64()),
            "nDetections": tbl.column("nDetections").cast(pa.int16()),
            "gmag": tbl.column("gMeanPSFMag").cast(pa.float32()),
            "rmag": tbl.column("rMeanPSFMag").cast(pa.float32()),
            "imag": tbl.column("iMeanPSFMag").cast(pa.float32()),
            "zmag": tbl.column("zMeanPSFMag").cast(pa.float32()),
            "ymag": tbl.column("yMeanPSFMag").cast(pa.float32()),
            "healpix_5": pa.array(hp5, type=pa.int32()),
        }
    )

    partitioning = pds.partitioning(
        pa.schema([("healpix_5", pa.int32())]), flavor="hive"
    )
    base = _basename(s3_path)
    with _write_lock:
        pds.write_dataset(
            out_tbl,
            base_dir=str(out_parquet),
            partitioning=partitioning,
            format="parquet",
            basename_template=f"{base}-part-{{i}}.parquet",
            existing_data_behavior="overwrite_or_ignore",
            file_options=pds.ParquetFileFormat().make_write_options(
                compression="zstd", compression_level=3
            ),
        )
    n_parts = int(np.unique(hp5).size)
    return out_tbl.num_rows, n_parts


def process_file(
    s3: pafs.S3FileSystem,
    s3_path: str,
    size: int,
    out_parquet: Path,
    done_dir: Path,
) -> dict:
    t0 = time.time()
    marker = done_dir / _marker_name(s3_path)
    if marker.exists():
        return {"path": s3_path, "skipped": True}
    try:
        n_rows, n_parts = transform_one(s3, s3_path, out_parquet)
        marker.write_text(f"{size}  {n_rows}  {n_parts}\n")
        dt = time.time() - t0
        return {"path": s3_path, "rows": n_rows, "parts": n_parts, "sec": dt}
    except Exception as e:
        log(f"FAIL {s3_path}: {e}")
        traceback.print_exc()
        return {"path": s3_path, "error": str(e)}


def main() -> int:
    global _log_path
    parser = argparse.ArgumentParser(
        description="Build the PS1 DR2 9-column HEALPix-5 Parquet cache."
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=default_cache_dir(),
        help=f"Output cache directory (default: {default_cache_dir()})",
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--limit", type=int, default=0, help="process only N files (0 = all)"
    )
    args = parser.parse_args()

    cache_dir: Path = args.cache_dir
    parquet_dir = cache_dir / "parquet"
    done_dir = cache_dir / ".done"
    _log_path = cache_dir / "progress.log"
    for d in (parquet_dir, done_dir):
        d.mkdir(parents=True, exist_ok=True)

    log(f"cache_dir={cache_dir}")
    log(f"opening S3 filesystem (anonymous, region={S3_REGION})")
    s3 = pafs.S3FileSystem(anonymous=True, region=S3_REGION)

    log("listing HATS source files...")
    t_list = time.time()
    files = list_s3_parquet(s3)
    files.sort(key=lambda x: x[0])
    log(f"listed {len(files):,} source files in {time.time() - t_list:.1f}s")

    if args.limit:
        files = files[: args.limit]

    pending = [
        (p, s) for p, s in files if not (done_dir / _marker_name(p)).exists()
    ]
    log(f"{len(pending):,} files to process ({len(files) - len(pending):,} already done)")

    total_rows = 0
    n_ok = 0
    n_fail = 0
    t_start = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(process_file, s3, p, sz, parquet_dir, done_dir): p
            for p, sz in pending
        }
        for i, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            if "error" in res:
                n_fail += 1
                continue
            if res.get("skipped"):
                continue
            n_ok += 1
            total_rows += res["rows"]
            if i % 50 == 0 or i == len(pending):
                elapsed = time.time() - t_start
                rate_per_hr = i / elapsed * 3600 if elapsed > 0 else 0
                eta_s = (len(pending) - i) / (i / elapsed) if i > 0 and elapsed > 0 else 0
                log(
                    f"{i}/{len(pending)} rows={total_rows:,} "
                    f"rate={rate_per_hr:.0f}/hr eta={eta_s/60:.0f}min"
                )

    dt = time.time() - t_start
    log(f"DONE ok={n_ok} fail={n_fail} rows={total_rows:,} elapsed={dt:.0f}s")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
