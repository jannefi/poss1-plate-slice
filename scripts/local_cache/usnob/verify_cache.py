"""
Verify the USNO-B1.0 local Parquet cache.

Checks:
  A. Sum of per-zone .done marker rows == Parquet dataset row count.
  B. Row count within 1% of published USNO-B1.0 total (1,045,913,669).
     VizieR TAP may exclude a small number of rows via internal filters;
     a tight tolerance would be misleading.
  C. All 1800 dec zones have a .done marker (full-sky coverage).
  D. Schema correctness.
  E. Partition coverage: every HP5 partition populated (USNO-B is
     all-sky so no exclusion zones, unlike PS1).

Defaults:
    --cache-dir / $VASCO_USNOB_CACHE / /Volumes/SANDISK/USNOB
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.dataset as pds

PUBLISHED_ROWS = 1_045_913_669
EXPECTED_ZONES = 1800
EXPECTED_PARTITIONS = 12288

EXPECTED_SCHEMA = {
    "id": pa.string(),
    "ra": pa.float64(),
    "dec": pa.float64(),
    "B1mag": pa.float32(),
    "R1mag": pa.float32(),
    "B2mag": pa.float32(),
    "R2mag": pa.float32(),
    "Imag": pa.float32(),
    "pmRA": pa.int16(),
    "pmDE": pa.int16(),
    "healpix_5": pa.int32(),
}


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_USNOB_CACHE", "/Volumes/SANDISK/USNOB"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify the USNO-B1.0 local Parquet cache.")
    parser.add_argument("--cache-dir", type=Path, default=default_cache_dir())
    args = parser.parse_args()

    cache = args.cache_dir
    parquet_dir = cache / "parquet"
    done_dir = cache / ".done"

    results: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str) -> None:
        mark = "PASS" if ok else "FAIL"
        results.append((name, ok, detail))
        print(f"[{mark}] {name}: {detail}")

    # ---------- A + C prep ----------
    marker_rows = 0
    zones_seen = 0
    bad_markers = 0
    for p in done_dir.iterdir():
        if not p.name.startswith("zone_") or p.suffix != ".done":
            continue
        try:
            parts = p.read_text().strip().split()
            marker_rows += int(parts[2])  # lo hi rows parts
            zones_seen += 1
        except Exception:
            bad_markers += 1

    ds = pds.dataset(str(parquet_dir), format="parquet", partitioning="hive")
    parquet_rows = ds.count_rows()

    # ---------- A ----------
    check(
        "A. Sum of .done marker rows == Parquet recount",
        marker_rows == parquet_rows and bad_markers == 0,
        f"markers={marker_rows:,} parquet={parquet_rows:,} bad={bad_markers}",
    )

    # ---------- B ----------
    pct = 100.0 * parquet_rows / PUBLISHED_ROWS
    check(
        "B. Row count within 1% of published USNO-B1.0 total",
        abs(parquet_rows - PUBLISHED_ROWS) < 0.01 * PUBLISHED_ROWS,
        f"parquet={parquet_rows:,} published={PUBLISHED_ROWS:,} ({pct:.2f}%)",
    )

    # ---------- C ----------
    check(
        "C. All 1800 dec zones processed",
        zones_seen == EXPECTED_ZONES,
        f"zones_seen={zones_seen}/{EXPECTED_ZONES}",
    )

    # ---------- D ----------
    schema_ok = True
    schema_detail = []
    for name, expected_type in EXPECTED_SCHEMA.items():
        if name not in ds.schema.names:
            schema_ok = False
            schema_detail.append(f"{name}: MISSING")
            continue
        got = ds.schema.field(name).type
        if got != expected_type:
            schema_ok = False
            schema_detail.append(f"{name}: got {got}, want {expected_type}")
    check(
        "D. Schema matches expected (11 columns incl. healpix_5)",
        schema_ok,
        "all columns correct" if schema_ok else "; ".join(schema_detail),
    )

    # ---------- E ----------
    partition_ids = {
        int(p.name.split("=")[1])
        for p in parquet_dir.iterdir()
        if p.name.startswith("healpix_5=")
    }
    missing = sorted(set(range(EXPECTED_PARTITIONS)) - partition_ids)
    check(
        "E. All 12,288 HP5 partitions populated (USNO-B is all-sky)",
        len(missing) == 0,
        f"present={len(partition_ids)}/{EXPECTED_PARTITIONS}"
        + (f" missing_sample={missing[:5]}" if missing else ""),
    )

    print()
    n_pass = sum(1 for _, ok, _ in results if ok)
    n_fail = sum(1 for _, ok, _ in results if not ok)
    print(f"TOTAL: {n_pass} pass, {n_fail} fail")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
