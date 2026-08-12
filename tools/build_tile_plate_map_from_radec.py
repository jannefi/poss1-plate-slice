#!/usr/bin/env python3
"""Derive tile_id -> plate_id from the slice runner's own radec output.

Downstream stages (the shrinking-set S0/S1/S2 chain) need to know which plate
each tile came from. There are two routes through this project and they have
different sources of truth for that:

  * **Archive route** -- tiles fetched from the STScI cutout service. The
    service chooses which plate covers the requested position, and it is not
    always the plate a plan assumed, so the answer must be read back out of
    each tile's own FITS header. That is
    `tools/build_tile_plate_map_from_headers.py`.

  * **Slice route (this tool)** -- tiles cut locally out of a named IRSA
    full-plate scan. Here there is no ambiguity to resolve: the runner sliced a
    specific plate and stamped `plate_id` onto every row it emitted. The radec
    CSV *is* the record of what was actually done.

Do not substitute one for the other. Reading slice-route provenance out of FITS
headers would work but is slower and needs tiles that no longer exist -- the
runner keeps one plate's tile tree on disk at a time and deletes it.

Reads only the two columns it needs, so cost is set by file count rather than by
the ~190M detection rows in a full-survey run.

Usage:
    python3 tools/build_tile_plate_map_from_radec.py \\
        --radec-dir work/slice/radec --out-csv work/slice/tile_plate_map.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--radec-dir", required=True, type=Path,
                    help="Directory of per-plate radec CSVs from run_fullscale_slice.py.")
    ap.add_argument("--out-csv", required=True, type=Path)
    ap.add_argument("--filtered-dir", type=Path, default=None,
                    help="Optional survivors tree. If given, reports tile dirs present "
                         "on disk but absent from the map, which would otherwise vanish "
                         "silently from any footprint derived downstream.")
    args = ap.parse_args()

    files = sorted(args.radec_dir.glob("*.csv"))
    if not files:
        sys.exit(f"[FAIL] no radec CSVs under {args.radec_dir}")

    mapping: dict[str, str] = {}
    conflicts: list[tuple[str, str, str]] = []
    for f in files:
        with f.open() as fh:
            rdr = csv.reader(fh)
            try:
                hdr = next(rdr)
            except StopIteration:
                continue
            try:
                ti, pi = hdr.index("tile_id"), hdr.index("plate_id")
            except ValueError:
                sys.exit(f"[FAIL] {f} lacks tile_id/plate_id columns: {hdr}")
            n = max(ti, pi)
            for row in rdr:
                if len(row) <= n:
                    continue
                t, p = row[ti], row[pi]
                prev = mapping.get(t)
                if prev is None:
                    mapping[t] = p
                elif prev != p:
                    # A tile attributed to two plates means the run mixed sources;
                    # downstream provenance would be wrong, so surface it loudly.
                    conflicts.append((t, prev, p))

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tile_id", "plate_id"])
        for t in sorted(mapping):
            w.writerow([t, mapping[t]])

    print(f"[OK] {len(files)} radec files -> {len(mapping):,} tiles, "
          f"{len(set(mapping.values()))} plates")
    print(f"     wrote {args.out_csv}")

    if conflicts:
        print(f"[FAIL] {len(conflicts)} tile(s) attributed to more than one plate, "
              f"e.g. {conflicts[:3]}", file=sys.stderr)
        sys.exit(1)

    if args.filtered_dir and args.filtered_dir.is_dir():
        on_disk = {d.name for d in args.filtered_dir.iterdir() if d.is_dir()}
        unmapped = sorted(on_disk - set(mapping))
        print(f"     survivors tree: {len(on_disk):,} tile dirs, "
              f"{len(unmapped)} not in the map")
        if unmapped:
            # A tile with zero detections emits no radec rows, so it cannot appear
            # here -- but it still occupies footprint. Downstream footprint counts
            # must come from the stage builder's tile_manifest.csv, never from the
            # tile_ids present in a stage CSV.
            print(f"     e.g. {unmapped[:5]}")


if __name__ == "__main__":
    main()
