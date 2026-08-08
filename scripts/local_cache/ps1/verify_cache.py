"""
Verify the PS1 DR2 local Parquet cache.

Unlike Gaia (which has an ESA-published _MD5SUM.txt manifest), the PS1
HATS source on AWS has no per-file checksum sidecar. Integrity is
enforced implicitly by Parquet's own file-level structure (read
failures raise) and by row-count cross-checks.

Checks performed:
  A. Published HATS row count vs. sum of per-file .done markers
  B. Published HATS row count vs. Parquet dataset recount
  C. Internal consistency: .done marker sum == Parquet recount
  D. Source file coverage: .done marker count == expected HATS file count
  E. Partition coverage: every HP5 pixel with dec > -30° has at least
     one Parquet file (Pan-STARRS is a northern-hemisphere survey and
     does not cover dec < ~-30°, so ~3,000 southern HP5 pixels are
     expected to be empty)
  F. Schema correctness: all 10 expected columns present with the
     expected types

Defaults:
    --cache-dir / $VASCO_PS1_CACHE / /Volumes/SANDISK/PS1
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.dataset as pds
import pyarrow.fs as pafs
from astropy_healpix import HEALPix
import astropy.units as u

S3_BUCKET = "stpubdata"
S3_PREFIX = "panstarrs/ps1/public/hats/otmo/otmo/dataset/"
PUBLISHED_HATS_ROWS = 10_560_724_292
EXPECTED_SOURCE_FILES = 27_161
EXPECTED_PARTITIONS = 12_288

EXPECTED_SCHEMA = {
    "objID": pa.int64(),
    "ra": pa.float64(),
    "dec": pa.float64(),
    "nDetections": pa.int16(),
    "gmag": pa.float32(),
    "rmag": pa.float32(),
    "imag": pa.float32(),
    "zmag": pa.float32(),
    "ymag": pa.float32(),
    "healpix_5": pa.int32(),
}


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_PS1_CACHE", "/Volumes/SANDISK/PS1"))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Defence-in-depth verification of the PS1 local Parquet cache."
    )
    parser.add_argument("--cache-dir", type=Path, default=default_cache_dir())
    parser.add_argument(
        "--skip-network",
        action="store_true",
        help="Skip check D (no S3 listing call).",
    )
    args = parser.parse_args()

    cache = args.cache_dir
    done_dir = cache / ".done"
    parquet_dir = cache / "parquet"

    results: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str) -> None:
        mark = "PASS" if ok else "FAIL"
        results.append((name, ok, detail))
        print(f"[{mark}] {name}: {detail}")

    # Sum rows from .done markers
    total_rows_from_markers = 0
    bad_markers = 0
    done_files = [p for p in done_dir.iterdir() if p.suffix == ".done"]
    for p in done_files:
        try:
            parts = p.read_text().strip().split()
            # Format: "size  rows  partitions"
            total_rows_from_markers += int(parts[1])
        except Exception:
            bad_markers += 1

    # ---------- A ----------
    check(
        "A. Sum of .done marker rows == published HATS row count",
        total_rows_from_markers == PUBLISHED_HATS_ROWS and bad_markers == 0,
        f"markers={total_rows_from_markers:,} published={PUBLISHED_HATS_ROWS:,} "
        f"bad_markers={bad_markers}",
    )

    # ---------- B ----------
    ds = pds.dataset(str(parquet_dir), format="parquet", partitioning="hive")
    parquet_rows = ds.count_rows()
    check(
        "B. Parquet dataset recount == published HATS row count",
        parquet_rows == PUBLISHED_HATS_ROWS,
        f"parquet={parquet_rows:,} published={PUBLISHED_HATS_ROWS:,}",
    )

    # ---------- C ----------
    check(
        "C. Parquet rows == sum(.done markers) (internal consistency)",
        parquet_rows == total_rows_from_markers,
        f"parquet={parquet_rows:,} markers={total_rows_from_markers:,}",
    )

    # ---------- D ----------
    if args.skip_network:
        check(
            "D. .done marker count == expected HATS file count",
            len(done_files) == EXPECTED_SOURCE_FILES,
            f"done={len(done_files)} expected={EXPECTED_SOURCE_FILES} (S3 list skipped)",
        )
    else:
        try:
            s3 = pafs.S3FileSystem(anonymous=True, region="us-east-1")
            selector = pafs.FileSelector(
                f"{S3_BUCKET}/{S3_PREFIX}", recursive=True
            )
            infos = s3.get_file_info(selector)
            s3_parquet_count = sum(
                1
                for i in infos
                if i.is_file
                and i.path.endswith(".parquet")
                and not i.path.endswith("/_metadata")
                and not i.path.endswith("/_common_metadata")
            )
            check(
                "D. .done marker count == live S3 HATS file count",
                len(done_files) == s3_parquet_count,
                f"done={len(done_files)} s3={s3_parquet_count}",
            )
        except Exception as e:
            check("D. .done marker count == live S3 HATS file count", False, f"S3 list failed: {e}")

    # ---------- E ----------
    # PS1 is a northern-hemisphere survey (dec > ~-30°). HP5 pixels whose
    # centers are south of that are expected to have no data; only check
    # that pixels north of the PS1 limit are populated, and that any
    # missing pixel has its center in the exclusion zone.
    partitions_on_disk = {
        p.name.split("=")[1]
        for p in parquet_dir.iterdir()
        if p.name.startswith("healpix_5=")
    }
    partition_ids = {int(x) for x in partitions_on_disk if x.isdigit()}
    full = set(range(EXPECTED_PARTITIONS))
    missing = sorted(full - partition_ids)
    hp = HEALPix(nside=32, order="nested")
    if missing:
        _, lats = hp.healpix_to_lonlat(np.array(missing))
        missing_decs = lats.to(u.deg).value
        non_polar_missing = [
            m for m, d in zip(missing, missing_decs) if d > -30.0
        ]
    else:
        non_polar_missing = []
    check(
        "E. All HP5 pixels with dec > -30° populated (PS1 has no southern data)",
        len(non_polar_missing) == 0,
        f"present={len(partition_ids)}/{EXPECTED_PARTITIONS} "
        f"missing={len(missing)} (all_below_-30={len(missing)-len(non_polar_missing)}) "
        + (
            f"unexpected_missing_sample={non_polar_missing[:5]}"
            if non_polar_missing
            else "expected_southern_gap"
        ),
    )

    # ---------- F ----------
    schema_ok = True
    schema_detail = []
    for name, expected_type in EXPECTED_SCHEMA.items():
        field = ds.schema.field(name) if name in ds.schema.names else None
        if field is None:
            schema_ok = False
            schema_detail.append(f"{name}: MISSING")
        elif field.type != expected_type:
            schema_ok = False
            schema_detail.append(f"{name}: got {field.type}, want {expected_type}")
    check(
        "F. Schema matches expected (10 columns, correct types)",
        schema_ok,
        "all 10 columns correct" if schema_ok else "; ".join(schema_detail),
    )

    print()
    n_pass = sum(1 for _, ok, _ in results if ok)
    n_fail = sum(1 for _, ok, _ in results if not ok)
    print(f"TOTAL: {n_pass} pass, {n_fail} fail")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
