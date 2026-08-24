#!/usr/bin/env python3
"""Split data/plate_manifest.csv into two roughly-equal, alternating plate
lists for the 2-VM full642 GCP run -- alternating (even/odd row index) rather
than first-half/second-half so any regional clustering in the manifest's
ordering averages out between the two arms.

Writes work/runs/full642_gcp_A/plates.txt and .../full642_gcp_B/plates.txt,
each a single line of comma-joined plate IDs -- the exact format the existing
scripts/gcp/*.sh already expect via `PLATES=$(cat plates.txt)`.

Usage: python3 scripts/gcp/split_manifest.py
"""
import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MANIFEST = REPO_ROOT / "data" / "plate_manifest.csv"

with MANIFEST.open() as f:
    plates = [row["plate_id"] for row in csv.DictReader(f)]

arm_a = plates[0::2]
arm_b = plates[1::2]

for tag, arm in (("A", arm_a), ("B", arm_b)):
    out_dir = REPO_ROOT / "work" / "runs" / f"full642_gcp_{tag}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "plates.txt"
    out_path.write_text(",".join(arm) + "\n")
    print(f"arm {tag}: {len(arm)} plates -> {out_path}")

assert len(arm_a) + len(arm_b) == len(plates)
assert set(arm_a).isdisjoint(arm_b)
print(f"total: {len(plates)} plates, split confirmed disjoint")
