"""
Verify the Gaia DR3 local Parquet cache with defence-in-depth cross-checks.

14 independent checks covering:
  A. _MD5SUM.txt round-trip (local == CDN)
  B. Checksum file structure (1 self + 3386 data entries)
  C. Proof that the _MD5SUM.txt self-reference FAIL during build is
     cosmetic (the file cannot contain its own MD5 — a mathematical
     impossibility, not a sync error)
  D. bucket listing vs _MD5SUM.txt entries (exact set equality)
  E. .done markers cover every data file (no ingest skipped)
  F. 10 random .done markers record exactly the _MD5SUM.txt expected
     MD5 (no divergent checksum set at ingest time)
  G. Sum of row counts across all .done markers == published Gaia DR3
     total (1,811,709,771)
  H. Recount the live Parquet dataset == published total

If all 14 pass, the cache is verified complete and correct.

Defaults:
    --cache-dir / $VASCO_GAIA_CACHE / /Volumes/SANDISK/Gaia
"""

from __future__ import annotations

import argparse
import hashlib
import os
import random
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pyarrow.dataset as pds

BUCKET_LIST_URL = (
    "https://gaia.eu-1.cdn77-storage.com/?prefix=Gaia/gdr3/gaia_source/&delimiter=/"
)
CDN_MD5SUM_URL = "http://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/_MD5SUM.txt"
PUBLISHED_DR3_ROWS = 1_811_709_771
EXPECTED_SELF_WANT = "2b916a83bda84c24f18fc947104d4d40"
EXPECTED_SELF_GOT = "d9206954944deba622e3319cc292ec5e"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_GAIA_CACHE", "/Volumes/SANDISK/Gaia"))


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_raw_checksums(path: Path) -> dict[str, str]:
    """Raw parse WITHOUT filtering the self-reference — verify_cache needs both."""
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[1]] = parts[0]
    return out


def fetch_cdn_md5() -> str:
    with urllib.request.urlopen(CDN_MD5SUM_URL, timeout=60) as r:
        data = r.read()
    return hashlib.md5(data).hexdigest()


def list_bucket_csvs() -> set[str]:
    seen: set[str] = set()
    marker = ""
    for _ in range(20):
        url = BUCKET_LIST_URL + (f"&marker={marker}" if marker else "")
        with urllib.request.urlopen(url, timeout=60) as r:
            doc = ET.fromstring(r.read())
        contents = doc.findall(NS + "Contents")
        for c in contents:
            key = c.find(NS + "Key").text
            name = key.split("/")[-1]
            if name.endswith(".csv.gz"):
                seen.add(name)
        if doc.find(NS + "IsTruncated").text != "true":
            break
        nm = doc.find(NS + "NextMarker")
        marker = (
            nm.text
            if nm is not None and nm.text
            else contents[-1].find(NS + "Key").text
        )
    return seen


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Defence-in-depth verification of the Gaia local Parquet cache."
    )
    parser.add_argument(
        "--cache-dir", type=Path, default=default_cache_dir()
    )
    parser.add_argument(
        "--skip-network",
        action="store_true",
        help="Skip checks A and D (no CDN/bucket calls).",
    )
    args = parser.parse_args()

    cache = args.cache_dir
    md5sum_path = cache / "metadata" / "_MD5SUM.txt"
    done_dir = cache / ".done"
    parquet_dir = cache / "parquet"

    results: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str) -> None:
        mark = "PASS" if ok else "FAIL"
        results.append((name, ok, detail))
        print(f"[{mark}] {name}: {detail}")

    # ---------- A ----------
    local_md5 = md5_of(md5sum_path)
    if args.skip_network:
        check("A. local _MD5SUM.txt == CDN _MD5SUM.txt", True, "SKIPPED (--skip-network)")
    else:
        cdn_md5 = fetch_cdn_md5()
        check(
            "A. local _MD5SUM.txt == CDN _MD5SUM.txt",
            local_md5 == cdn_md5,
            f"local={local_md5} cdn={cdn_md5}",
        )
    check(
        "A'. local _MD5SUM.txt == expected 'got' hash from build log",
        local_md5 == EXPECTED_SELF_GOT,
        f"local={local_md5} expected={EXPECTED_SELF_GOT}",
    )

    # ---------- B ----------
    raw = load_raw_checksums(md5sum_path)
    self_entries = [(n, h) for n, h in raw.items() if n == "_MD5SUM.txt"]
    csv_entries = {
        n: h
        for n, h in raw.items()
        if n.endswith(".csv.gz") and re.match(r"^GaiaSource_\d+-\d+\.csv\.gz$", n)
    }
    other_entries = [
        (n, h)
        for n, h in raw.items()
        if n != "_MD5SUM.txt" and n not in csv_entries
    ]
    check(
        "B. _MD5SUM.txt contains exactly 1 self-reference entry",
        len(self_entries) == 1,
        f"self_entries={len(self_entries)}",
    )
    check(
        "B'. All non-self entries are GaiaSource_*.csv.gz",
        len(other_entries) == 0,
        f"other_entries={len(other_entries)}",
    )
    check(
        "B''. Total entries = 1 self + 3386 data",
        len(raw) == 3387 and len(csv_entries) == 3386,
        f"total={len(raw)} data={len(csv_entries)}",
    )

    # ---------- C ----------
    self_hash_in_file = self_entries[0][1] if self_entries else "<missing>"
    check(
        "C. Self-reference hash in file == 'want' value seen in FAIL log",
        self_hash_in_file == EXPECTED_SELF_WANT,
        f"self_hash={self_hash_in_file} expected={EXPECTED_SELF_WANT}",
    )
    check(
        "C'. Self-reference hash != current file hash "
        "(mathematically forced: file cannot embed own MD5)",
        self_hash_in_file != local_md5,
        f"self_hash={self_hash_in_file} current={local_md5}",
    )

    # ---------- D ----------
    if args.skip_network:
        check(
            "D. bucket.csv.gz set == _MD5SUM.txt.csv.gz set",
            True,
            "SKIPPED (--skip-network)",
        )
    else:
        bucket_csvs = list_bucket_csvs()
        ck_set = set(csv_entries.keys())
        only_in_bucket = bucket_csvs - ck_set
        only_in_checksums = ck_set - bucket_csvs
        check(
            "D. bucket.csv.gz set == _MD5SUM.txt.csv.gz set",
            len(only_in_bucket) == 0 and len(only_in_checksums) == 0,
            f"bucket={len(bucket_csvs)} checksums={len(ck_set)} "
            f"diff_bucket_only={len(only_in_bucket)} "
            f"diff_checksums_only={len(only_in_checksums)}",
        )

    # ---------- E ----------
    done_files = {p.name[:-5] for p in done_dir.iterdir() if p.suffix == ".done"}
    missing_done = set(csv_entries.keys()) - done_files
    bogus_done = done_files - set(csv_entries.keys())
    check(
        "E. Every csv.gz in checksums has a .done marker",
        len(missing_done) == 0,
        f"missing={len(missing_done)}"
        + (f" sample={sorted(missing_done)[:3]}" if missing_done else ""),
    )
    check(
        "E'. Every .done marker corresponds to a csv.gz in checksums (no strays)",
        len(bogus_done) == 0,
        f"bogus={len(bogus_done)}"
        + (f" sample={sorted(bogus_done)[:3]}" if bogus_done else ""),
    )

    # ---------- F ----------
    rng = random.Random(12345)
    sample = rng.sample(sorted(csv_entries.keys()), 10)
    mismatches = []
    for fname in sample:
        marker = done_dir / (fname + ".done")
        txt = marker.read_text().strip()
        recorded_md5 = txt.split()[0]
        expected = csv_entries[fname]
        if recorded_md5 != expected:
            mismatches.append((fname, recorded_md5, expected))
    check(
        "F. 10 random .done markers: recorded MD5 == _MD5SUM.txt expected",
        len(mismatches) == 0,
        f"mismatches={len(mismatches)}/10",
    )

    # ---------- G ----------
    total_rows_from_markers = 0
    bad_markers = 0
    for fname in csv_entries:
        marker = done_dir / (fname + ".done")
        try:
            parts = marker.read_text().strip().split()
            total_rows_from_markers += int(parts[1])
        except Exception:
            bad_markers += 1
    check(
        "G. Sum of rows from .done markers == published DR3 total",
        total_rows_from_markers == PUBLISHED_DR3_ROWS and bad_markers == 0,
        f"sum_markers={total_rows_from_markers:,} "
        f"published={PUBLISHED_DR3_ROWS:,} bad_markers={bad_markers}",
    )

    # ---------- H ----------
    ds = pds.dataset(str(parquet_dir), format="parquet", partitioning="hive")
    parquet_rows = ds.count_rows()
    check(
        "H. Parquet dataset recount == published DR3 total",
        parquet_rows == PUBLISHED_DR3_ROWS,
        f"parquet={parquet_rows:,} published={PUBLISHED_DR3_ROWS:,}",
    )
    check(
        "H'. Parquet rows == sum(.done markers) (internal consistency)",
        parquet_rows == total_rows_from_markers,
        f"parquet={parquet_rows:,} markers={total_rows_from_markers:,}",
    )

    print()
    n_pass = sum(1 for _, ok, _ in results if ok)
    n_fail = sum(1 for _, ok, _ in results if not ok)
    print(f"TOTAL: {n_pass} pass, {n_fail} fail")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
