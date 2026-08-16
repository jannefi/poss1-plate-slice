#!/usr/bin/env python3
"""One-way parity of a published S0 catalogue against this project's S0.

The question this answers: of the candidates a *previous* pipeline published,
how many does the current pipeline also produce? It is deliberately one-way.
The reference run here (VASCO60 final_release_v1) applied a 30' circular cut
inside each 60'x60' tile and covered a wider plate set, so the two catalogues
are neither nested nor co-extensive; recall in the reference -> ours direction
is the only well-posed statement. The reverse ratio is a candidate-count
excess, not parity, and this tool does not report it.

Scoping is reported at three levels side by side so the denominator is never
hidden:

  all           every reference row
  plate-scoped  reference rows whose own plate_id is one this run processed
  footprint     reference rows geometrically inside this run's tile footprint,
                using the parity funnel's 0.75 deg membership radius

The footprint comes from the run's ``tile_manifest.csv``, never from the tile
ids appearing in the stage CSV -- the latter silently drops every tile that
produced no survivor and shrinks the footprint (see the funnel's own notes).

Two arms:

  A  reference -> our S0 (the deduplicated survivor set). Headline.
  B  reference -> our *raw* SExtractor detections, per plate. Splits every
     arm-A miss into "never detected" and "detected, then removed by the veto
     and filter chain", which is the part that actually distinguishes two
     pipelines rather than two sky coverages.

Arm B's radec cache carries ALPHA_J2000, i.e. the raw plate WCS, while the S0
catalogue carries the WCS-fixed RA_corr/Dec_corr. On the ~1/3 of plates where
WCSFIX moves coordinates the two differ by up to ~2.34", so arm B is reported
as a radius sweep and never as a single number.

How to validate
---------------
    python3 tools/archive_slice_parity.py \
        --ref-s0 <vasco60>/releases/final_release_v1/run/stage_S0.csv \
        --ref-manifest <vasco60>/releases/final_release_v1/run/tile_manifest.csv \
        --our-s0 results/s0-642-20260814/stage_S0.csv.gz \
        --our-manifest results/s0-642-20260814/tile_manifest.csv.gz \
        --radec-dir <run>/radec \
        --out-dir work/archive_slice_parity

Self-tests run automatically and abort on failure:
  * every reference row below our catalogue's declination floor must fall
    outside the footprint;
  * per-plate row counts must sum to the footprint-scoped total;
  * arm-B distances must never exceed arm-A distances for the same row
    (S0 is a subset of the raw detections, so the nearest raw source can only
    be closer or equal -- modulo the WCS caveat above, which is reported).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

MEMBERSHIP_RADIUS_DEG = 0.75  # > 60' tile half-diagonal (~0.707 deg), same as the parity funnel
RADII_ARCSEC = (1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 30.0)
HEADLINE_RADIUS_ARCSEC = 5.0

TILE_ID_RE = re.compile(r"^tile_RA(?P<ra>[0-9.]+)_DEC(?P<sign>[pm])(?P<dec>[0-9.]+)$")


# --------------------------------------------------------------------------
# geometry helpers -- everything is done on unit vectors so RA wrap and the
# poles need no special cases.
# --------------------------------------------------------------------------

def unit_vectors(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.radians(np.asarray(ra_deg, dtype=np.float64))
    dec = np.radians(np.asarray(dec_deg, dtype=np.float64))
    cd = np.cos(dec)
    return np.column_stack((cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)))


def chord(angle_deg: float) -> float:
    return 2.0 * np.sin(np.radians(angle_deg) / 2.0)


def chord_to_arcsec(d: np.ndarray) -> np.ndarray:
    d = np.clip(np.asarray(d, dtype=np.float64) / 2.0, -1.0, 1.0)
    return np.degrees(2.0 * np.arcsin(d)) * 3600.0


def tile_centres(tile_ids) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse tile centres out of the tile_id naming convention.

    Returns (ra, dec, ok_mask). The convention is locked in 02_DECISIONS, so a
    tile id that does not parse is a real problem and is surfaced, not skipped
    silently.
    """
    ra = np.full(len(tile_ids), np.nan)
    dec = np.full(len(tile_ids), np.nan)
    for i, tid in enumerate(tile_ids):
        m = TILE_ID_RE.match(str(tid).strip())
        if not m:
            continue
        ra[i] = float(m.group("ra"))
        d = float(m.group("dec"))
        dec[i] = d if m.group("sign") == "p" else -d
    return ra, dec, np.isfinite(ra) & np.isfinite(dec)


# --------------------------------------------------------------------------
# arm B worker -- one plate per call, so each process holds one plate's
# detections and nothing else. Never build a global tree over the whole radec
# cache; it is ~230M rows.
# --------------------------------------------------------------------------

_ARM_B_CTX: dict = {}


def _arm_b_init(ref_ra: np.ndarray, ref_dec: np.ndarray, radec_dir: str) -> None:
    _ARM_B_CTX["ref_xyz"] = unit_vectors(ref_ra, ref_dec)
    _ARM_B_CTX["radec_dir"] = radec_dir


def _arm_b_plate(task: tuple[str, np.ndarray]):
    plate, row_idx = task
    path = Path(_ARM_B_CTX["radec_dir"]) / f"{plate}.csv"
    if not path.exists():
        return plate, row_idx, None, 0
    det = pd.read_csv(path, usecols=["ra", "dec"], dtype={"ra": np.float64, "dec": np.float64})
    if det.empty:
        return plate, row_idx, None, 0
    tree = cKDTree(unit_vectors(det["ra"].to_numpy(), det["dec"].to_numpy()))
    d, _ = tree.query(_ARM_B_CTX["ref_xyz"][row_idx], k=1)
    return plate, row_idx, chord_to_arcsec(d), int(len(det))


# --------------------------------------------------------------------------

def recall_table(dist_arcsec: np.ndarray, mask: np.ndarray) -> dict:
    n = int(mask.sum())
    out = {"n": n}
    if n == 0:
        out.update({f"{r:g}": None for r in RADII_ARCSEC})
        return out
    d = dist_arcsec[mask]
    finite = np.isfinite(d)
    for r in RADII_ARCSEC:
        out[f"{r:g}"] = float(100.0 * np.count_nonzero(finite & (d <= r)) / n)
    return out


def fmt_recall(label: str, row: dict) -> str:
    cells = " ".join(
        f"{row[f'{r:g}']:6.2f}" if row[f"{r:g}"] is not None else "    --" for r in RADII_ARCSEC
    )
    return f"  {label:<26} {row['n']:>7,}  {cells}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-s0", required=True, help="published stage_S0.csv to test recall OF")
    ap.add_argument("--ref-manifest", required=True, help="published tile_manifest.csv (tile_id -> plate_id)")
    ap.add_argument("--our-s0", required=True, help="this project's stage_S0.csv[.gz]")
    ap.add_argument("--our-manifest", required=True, help="this project's tile_manifest.csv[.gz]")
    ap.add_argument("--radec-dir", default=None, help="per-plate raw detection CSVs; enables arm B")
    ap.add_argument("--membership-deg", type=float, default=MEMBERSHIP_RADIUS_DEG)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---------------- load ----------------
    ref = pd.read_csv(args.ref_s0)
    ref_man = pd.read_csv(args.ref_manifest, usecols=["tile_id", "plate_id"])
    ref = ref.merge(ref_man.drop_duplicates("tile_id"), on="tile_id", how="left")
    ref = ref.rename(columns={"plate_id": "ref_plate"})

    ours = pd.read_csv(args.our_s0, usecols=["ra", "dec"])
    our_man = pd.read_csv(args.our_manifest, usecols=["tile_id", "plate_id"])

    print(f"reference S0 : {len(ref):,} rows, {ref['ref_plate'].nunique()} plates")
    print(f"our S0       : {len(ours):,} rows")
    print(f"our manifest : {len(our_man):,} tiles, {our_man['plate_id'].nunique()} plates")

    t_ra, t_dec, ok = tile_centres(our_man["tile_id"].to_numpy())
    if not ok.all():
        print(f"[FATAL] {int((~ok).sum())} tile ids in our manifest do not parse", file=sys.stderr)
        return 2
    our_man = our_man.assign(t_ra=t_ra, t_dec=t_dec)

    our_plates = set(our_man["plate_id"].astype(str))
    ref_plates = set(ref["ref_plate"].dropna().astype(str))
    print(f"plates: reference {len(ref_plates)}, ours {len(our_plates)}, "
          f"ours-not-in-reference {len(our_plates - ref_plates)}")

    # ---------------- footprint ----------------
    tile_tree = cKDTree(unit_vectors(our_man["t_ra"].to_numpy(), our_man["t_dec"].to_numpy()))
    ref_xyz = unit_vectors(ref["ra"].to_numpy(), ref["dec"].to_numpy())
    d_tile, i_tile = tile_tree.query(ref_xyz, k=1)
    dist_tile_deg = np.degrees(2.0 * np.arcsin(np.clip(d_tile / 2.0, -1.0, 1.0)))
    in_footprint = dist_tile_deg <= args.membership_deg
    nearest_plate = our_man["plate_id"].to_numpy()[i_tile]

    in_plates = ref["ref_plate"].astype(str).isin(our_plates).to_numpy()

    # ---------------- arm A ----------------
    s0_tree = cKDTree(unit_vectors(ours["ra"].to_numpy(), ours["dec"].to_numpy()))
    d_s0, _ = s0_tree.query(ref_xyz, k=1)
    dist_s0 = chord_to_arcsec(d_s0)

    # ---------------- arm B ----------------
    dist_raw = np.full(len(ref), np.inf)
    plate_of_raw = np.array([""] * len(ref), dtype=object)
    det_counts: dict[str, int] = {}
    if args.radec_dir:
        # A reference row is searched against every plate whose tile footprint
        # contains it, not just the nearest -- plates overlap heavily and the
        # nearest tile centre is not necessarily on the plate that shows the
        # source.
        neighbours = tile_tree.query_ball_point(ref_xyz, r=chord(args.membership_deg))
        plate_arr = our_man["plate_id"].to_numpy()
        by_plate: dict[str, list[int]] = {}
        for row_i, tiles in enumerate(neighbours):
            if not tiles:
                continue
            for p in set(plate_arr[tiles]):
                by_plate.setdefault(str(p), []).append(row_i)
        tasks = [(p, np.asarray(sorted(idx), dtype=np.int64)) for p, idx in sorted(by_plate.items())]
        print(f"arm B: {len(tasks)} plates to scan, "
              f"{sum(len(t[1]) for t in tasks):,} (row, plate) pairs")
        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=_arm_b_init,
            initargs=(ref["ra"].to_numpy(), ref["dec"].to_numpy(), args.radec_dir),
        ) as ex:
            for n_done, (plate, row_idx, dists, n_det) in enumerate(
                ex.map(_arm_b_plate, tasks, chunksize=1), start=1
            ):
                det_counts[plate] = n_det
                if dists is not None:
                    better = dists < dist_raw[row_idx]
                    plate_of_raw[row_idx[better]] = plate
                    dist_raw[row_idx[better]] = dists[better]
                if n_done % 50 == 0 or n_done == len(tasks):
                    print(f"    {n_done}/{len(tasks)} plates", flush=True)

    # ---------------- self-tests ----------------
    failures = []
    dec_floor = float(pd.read_csv(args.our_s0, usecols=["dec"])["dec"].min())
    below = ref["dec"].to_numpy() < (dec_floor - args.membership_deg)
    if below.any() and in_footprint[below].any():
        failures.append(
            f"{int(in_footprint[below].sum())} reference rows below our declination floor "
            f"({dec_floor:.3f}) were scored in-footprint")
    if args.radec_dir:
        # S0 is a subset of the raw detections, so the nearest raw source can
        # only be closer -- except where WCSFIX moved the S0 coordinate.
        both = np.isfinite(dist_s0) & np.isfinite(dist_raw) & in_footprint
        viol = both & (dist_raw > dist_s0 + 2.5)
        if viol.any():
            failures.append(
                f"{int(viol.sum())} in-footprint rows have arm-B distance > arm-A + 2.5\" "
                f"(max excess {float((dist_raw - dist_s0)[viol].max()):.2f}\")")

    # ---------------- per-row ledger ----------------
    rows = pd.DataFrame({
        "src_id": ref["src_id"],
        "ra": ref["ra"],
        "dec": ref["dec"],
        "ref_tile_id": ref["tile_id"],
        "ref_plate": ref["ref_plate"],
        "ref_plate_processed": in_plates,
        "dist_nearest_tile_deg": np.round(dist_tile_deg, 6),
        "in_footprint": in_footprint,
        "our_nearest_plate": nearest_plate,
        "dist_s0_arcsec": np.round(dist_s0, 4),
        "dist_raw_arcsec": np.where(np.isfinite(dist_raw), np.round(dist_raw, 4), np.nan),
        "raw_plate": plate_of_raw,
    })
    rows.to_csv(out_dir / "rows.csv", index=False)

    # ---------------- per-plate ----------------
    fp = rows[rows["in_footprint"]].copy()
    if len(fp) != int(in_footprint.sum()):
        failures.append("footprint row count disagrees between ledger and mask")
    grp = fp.groupby("ref_plate", dropna=False)
    per_plate = pd.DataFrame({
        "ref_rows_in_footprint": grp.size(),
        "s0_match_5as": grp["dist_s0_arcsec"].apply(lambda s: int((s <= HEADLINE_RADIUS_ARCSEC).sum())),
        "raw_match_5as": grp["dist_raw_arcsec"].apply(lambda s: int((s <= HEADLINE_RADIUS_ARCSEC).sum())),
        "median_dist_s0_arcsec": grp["dist_s0_arcsec"].median().round(3),
    }).reset_index()
    per_plate["s0_recall_pct"] = (100.0 * per_plate["s0_match_5as"] / per_plate["ref_rows_in_footprint"]).round(2)
    per_plate["raw_recall_pct"] = (100.0 * per_plate["raw_match_5as"] / per_plate["ref_rows_in_footprint"]).round(2)
    per_plate = per_plate.sort_values("s0_recall_pct")
    per_plate.to_csv(out_dir / "per_plate.csv", index=False)
    if int(per_plate["ref_rows_in_footprint"].sum()) != len(fp):
        failures.append("per-plate rows do not sum to the footprint-scoped total")

    # ---------------- report ----------------
    scopes = {
        "all": np.ones(len(ref), dtype=bool),
        "plate_scoped": in_plates,
        "footprint": in_footprint,
    }
    summary = {
        "ref_s0": os.path.abspath(args.ref_s0),
        "our_s0": os.path.abspath(args.our_s0),
        "radec_dir": os.path.abspath(args.radec_dir) if args.radec_dir else None,
        "membership_radius_deg": args.membership_deg,
        "n_ref_rows": int(len(ref)),
        "n_our_rows": int(len(ours)),
        "n_our_tiles": int(len(our_man)),
        "our_dec_floor": dec_floor,
        "scopes": {k: int(v.sum()) for k, v in scopes.items()},
        "arm_a_recall_pct": {k: recall_table(dist_s0, v) for k, v in scopes.items()},
        "self_test_failures": failures,
    }
    if args.radec_dir:
        summary["arm_b_recall_pct"] = {k: recall_table(dist_raw, v) for k, v in scopes.items()}
        fpm = in_footprint
        matched_s0 = fpm & (dist_s0 <= HEADLINE_RADIUS_ARCSEC)
        matched_raw = fpm & (dist_raw <= HEADLINE_RADIUS_ARCSEC)
        summary["decomposition_at_5as"] = {
            "in_footprint": int(fpm.sum()),
            "in_our_s0": int(matched_s0.sum()),
            "detected_but_not_in_s0": int((matched_raw & ~matched_s0).sum()),
            "not_detected": int((fpm & ~matched_raw).sum()),
        }
        summary["raw_detections_scanned"] = int(sum(det_counts.values()))

    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    hdr = " ".join(f"{r:>6g}" for r in RADII_ARCSEC)
    print("\n=== ARM A: reference S0 -> our S0, recall %% at radius (arcsec) ===")
    print(f"  {'scope':<26} {'n':>7}  {hdr}")
    for k, v in scopes.items():
        print(fmt_recall(k, summary["arm_a_recall_pct"][k]))
    if args.radec_dir:
        print("\n=== ARM B: reference S0 -> our RAW detections (pre-WCSFIX coords) ===")
        print(f"  {'scope':<26} {'n':>7}  {hdr}")
        for k, v in scopes.items():
            print(fmt_recall(k, summary["arm_b_recall_pct"][k]))
        d = summary["decomposition_at_5as"]
        n = d["in_footprint"]
        print(f"\n=== DECOMPOSITION at {HEADLINE_RADIUS_ARCSEC:g}\", footprint-scoped (n={n:,}) ===")
        for key, label in (("in_our_s0", "in our S0"),
                           ("detected_but_not_in_s0", "detected, chain removed it"),
                           ("not_detected", "never detected")):
            print(f"  {label:<28} {d[key]:>7,}  {100.0*d[key]/n:6.2f}%")

    print("\n=== WORST PLATES by S0 recall (>=20 reference rows in footprint) ===")
    worst = per_plate[per_plate["ref_rows_in_footprint"] >= 20].head(15)
    print(worst.to_string(index=False))

    print(f"\nledgers -> {out_dir}/rows.csv, per_plate.csv, summary.json")
    if failures:
        print("\n[SELF-TEST FAILURES]", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1
    print("self-tests: all passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
