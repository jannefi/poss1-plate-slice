"""
Build the Gaia DR3 6-column HEALPix-5 Parquet cache from ESA's bulk CDN.

This is the exact script used to build /Volumes/SANDISK/Gaia on 2026-04-11.
Resumable via per-file .done markers, MD5-verified, stream-deletes staged
CSVs after transform so disk usage stays bounded.

Defaults match the original build host; override with env vars or flags
for replication on another machine.

Paths (in priority order):
    --cache-dir      / $VASCO_GAIA_CACHE   / /Volumes/SANDISK/Gaia
    --staging-dir    / $VASCO_GAIA_STAGING / /private/tmp/claude/gaia_staging

Columns kept: source_id, ra, dec, pmra, pmdec, phot_g_mean_mag
Partition key: healpix_5 = source_id >> 49  (Gaia source_id encodes HEALPix-12
              nested; shifting 14 more bits drops to nested level 5.)

Note on proper motion precision: pmra/pmdec are deliberately stored as
float32. See MANIFEST.json "precision" section for the full rationale.
Float32 ULP at realistic pm magnitudes is 4-7 orders of magnitude below
Gaia's own measurement uncertainty.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import os
import sys
import threading
import time
import traceback
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.dataset as pds

BASE_URL = "http://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/"
CHECKSUM_FILE_NAME = "_MD5SUM.txt"
REQUIRED_COLS = ["source_id", "ra", "dec", "pmra", "pmdec", "phot_g_mean_mag"]
COLUMN_TYPES = {
    "source_id": pa.int64(),
    "ra": pa.float64(),
    "dec": pa.float64(),
    "pmra": pa.float32(),
    "pmdec": pa.float32(),
    "phot_g_mean_mag": pa.float32(),
}

_log_lock = threading.Lock()
_write_lock = threading.Lock()
_log_path: Path | None = None


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_GAIA_CACHE", "/Volumes/SANDISK/Gaia"))


def default_staging_dir() -> Path:
    return Path(os.environ.get("VASCO_GAIA_STAGING", "/private/tmp/claude/gaia_staging"))


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


def load_checksums(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        parts = raw.split()
        if len(parts) == 2:
            md5, fname = parts
            checksums[fname] = md5
    # _MD5SUM.txt lists itself with a mathematically unsatisfiable hash
    # (a file cannot embed its own MD5 — any value embedded changes the file).
    # See verify_cache.py for the full cosmetic-proof. Filter it out.
    checksums.pop(CHECKSUM_FILE_NAME, None)
    return checksums


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_one(fname: str, staging: Path, attempts: int = 3) -> Path:
    dst = staging / fname
    tmp = dst.with_suffix(dst.suffix + ".part")
    url = BASE_URL + fname
    last_err: Exception | None = None
    for i in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            tmp.rename(dst)
            return dst
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_err = e
            log(f"download retry {i+1}/{attempts} for {fname}: {e}")
            time.sleep(2**i)
    raise RuntimeError(f"download failed for {fname}: {last_err}")


def strip_ecsv_header(path: Path) -> tuple[bytes, bytes]:
    """Return (column_header_line, data_bytes) with `#`-prefixed ECSV metadata removed."""
    with gzip.open(path, "rb") as f:
        while True:
            line = f.readline()
            if not line:
                raise RuntimeError(f"unexpected EOF in {path}")
            if not line.startswith(b"#"):
                header = line
                break
        data = f.read()
    return header, data


def transform_one(csv_gz: Path, out_parquet: Path) -> tuple[int, int]:
    """Read 6 columns from an ECSV.gz, partition by healpix_5, write Parquet.

    Returns (n_rows, n_partitions_written).
    """
    header, data = strip_ecsv_header(csv_gz)
    buf = io.BytesIO(header + data)
    table = pcsv.read_csv(
        buf,
        read_options=pcsv.ReadOptions(use_threads=True),
        convert_options=pcsv.ConvertOptions(
            include_columns=REQUIRED_COLS,
            column_types=COLUMN_TYPES,
            null_values=["", "null", "NaN", "nan", "NULL"],
            strings_can_be_null=True,
        ),
    )
    sid = table.column("source_id").to_numpy(zero_copy_only=False)
    hp5 = (sid >> 49).astype(np.int32)
    table = table.append_column("healpix_5", pa.array(hp5, type=pa.int32()))

    partitioning = pds.partitioning(
        pa.schema([("healpix_5", pa.int32())]), flavor="hive"
    )
    basename = csv_gz.name.replace(".csv.gz", "")
    with _write_lock:
        pds.write_dataset(
            table,
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
    return table.num_rows, n_parts


def process_file(
    fname: str,
    expected_md5: str,
    staging: Path,
    out_parquet: Path,
    done_dir: Path,
) -> dict:
    t0 = time.time()
    marker = done_dir / (fname + ".done")
    if marker.exists():
        return {"fname": fname, "skipped": True}
    try:
        csv_gz = download_one(fname, staging)
        got = md5_of(csv_gz)
        if got != expected_md5:
            csv_gz.unlink(missing_ok=True)
            raise RuntimeError(
                f"MD5 mismatch for {fname}: got {got}, want {expected_md5}"
            )
        n_rows, n_parts = transform_one(csv_gz, out_parquet)
        csv_gz.unlink(missing_ok=True)
        marker.write_text(f"{expected_md5}  {n_rows}  {n_parts}\n")
        dt = time.time() - t0
        return {"fname": fname, "rows": n_rows, "parts": n_parts, "sec": dt}
    except Exception as e:
        log(f"FAIL {fname}: {e}")
        traceback.print_exc()
        return {"fname": fname, "error": str(e)}


def fetch_checksum_file(meta_dir: Path) -> Path:
    meta_dir.mkdir(parents=True, exist_ok=True)
    dst = meta_dir / CHECKSUM_FILE_NAME
    if not dst.exists():
        url = BASE_URL + CHECKSUM_FILE_NAME
        with urllib.request.urlopen(url, timeout=60) as r:
            dst.write_bytes(r.read())
    return dst


def fetch_metadata_files(meta_dir: Path) -> None:
    meta_dir.mkdir(parents=True, exist_ok=True)
    for fname in ["_citation.txt", "_disclaimer.txt", "_license.txt", "_readme.txt"]:
        dst = meta_dir / fname
        if dst.exists():
            continue
        url = BASE_URL.rsplit("/", 2)[0] + "/" + fname
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                dst.write_bytes(r.read())
        except Exception as e:
            log(f"metadata fetch failed for {fname}: {e}")


def main() -> int:
    global _log_path
    parser = argparse.ArgumentParser(
        description="Build the Gaia DR3 6-column HEALPix-5 Parquet cache."
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=default_cache_dir(),
        help=f"Output cache directory (default: {default_cache_dir()})",
    )
    parser.add_argument(
        "--staging-dir",
        type=Path,
        default=default_staging_dir(),
        help=f"Staging directory for in-flight .csv.gz (default: {default_staging_dir()})",
    )
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument(
        "--limit", type=int, default=0, help="process only N files (0 = all)"
    )
    parser.add_argument(
        "--only",
        nargs="*",
        help="process only these filenames (space separated)",
    )
    args = parser.parse_args()

    cache_dir: Path = args.cache_dir
    staging: Path = args.staging_dir
    parquet_dir = cache_dir / "parquet"
    done_dir = cache_dir / ".done"
    meta_dir = cache_dir / "metadata"
    _log_path = cache_dir / "progress.log"

    for d in (staging, parquet_dir, done_dir, meta_dir):
        d.mkdir(parents=True, exist_ok=True)

    log(f"cache_dir={cache_dir}")
    log(f"staging_dir={staging}")
    log("fetching checksum + metadata files")
    checksum_path = fetch_checksum_file(meta_dir)
    fetch_metadata_files(meta_dir)
    checksums = load_checksums(checksum_path)
    log(f"loaded {len(checksums)} data-file checksums")

    if args.only:
        todo = [f for f in args.only if f in checksums]
    else:
        todo = sorted(checksums.keys())
    if args.limit:
        todo = todo[: args.limit]

    pending = [f for f in todo if not (done_dir / (f + ".done")).exists()]
    log(
        f"{len(pending)} files to process "
        f"({len(todo) - len(pending)} already done)"
    )

    total_rows = 0
    n_ok = 0
    n_fail = 0
    t_start = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {
            ex.submit(
                process_file, f, checksums[f], staging, parquet_dir, done_dir
            ): f
            for f in pending
        }
        for i, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            fname = res["fname"]
            if "error" in res:
                n_fail += 1
                continue
            if res.get("skipped"):
                continue
            n_ok += 1
            total_rows += res["rows"]
            rate = (i / (time.time() - t_start)) * 3600
            log(
                f"{i}/{len(pending)} {fname} rows={res['rows']:,} "
                f"parts={res['parts']} sec={res['sec']:.1f} ~{rate:.0f} files/hr"
            )

    dt = time.time() - t_start
    log(f"DONE ok={n_ok} fail={n_fail} rows={total_rows:,} elapsed={dt:.0f}s")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
