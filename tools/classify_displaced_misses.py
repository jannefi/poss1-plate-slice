#!/usr/bin/env python3
"""Is the nearby detection the SAME object, or did we miss this one?

For reference S0 rows with no detection of ours within 5" but one within 30"
(buckets B and C of `archive_slice_parity.py` §4.2), distance alone cannot say
whether our pipeline saw the same source with a displaced centroid or saw a
different source entirely and missed this one. Two independent tests decide it.

**Test 1 — is there signal at the reference position in our own pixels?**
Peak significance (peak - annulus median) / annulus MAD-sigma, measured on the
plate scan we slice from. Signal present means the source is in our data.

**Test 2 — does the reference run itself have a separate detection where our
nearest detection sits?** Read from its retained full `sextractor_pass2.csv`.
If the reference extracted a source at the reference position AND another at
our detection's position, the two are distinct objects, so our nearest
detection is a neighbour and we genuinely missed the reference row. If the
reference has nothing at our detection's position, the likelier reading is one
object whose centroid the two runs place differently.

Classification:

  SEPARATE   signal at ref position + reference has its own detection at our
             detection's position -> two objects, we missed one
  MISSED     signal at ref position, no reference detection at ours
             -> source is in our pixels and we did not extract it
  DISPLACED  no signal at ref position -> nothing there in our data; our
             nearby detection is the plausible counterpart, positions differ
  NOSIGNAL   no signal in either run's pixels at the reference position

Reading the reference tiles costs care: `sextractor_pass2.csv` is ~65 MB/tile
because of VIGNET, so always pass `usecols`. One tile at a time, never a
concatenated frame.

How to validate
---------------
    python3 tools/classify_displaced_misses.py \
        --rows work/archive_slice_parity/never_detected_58.csv \
        --buckets BC \
        --ref-tiles <vasco60>/data/tiles \
        --plate-dir <plate-scan-dir> \
        --radec-dir <run>/radec \
        --out-dir work/archive_slice_parity/displaced

Sanity: the recovered nearest-detection distance must reproduce
`dist_raw_arcsec` from the parity run; the tool asserts agreement to 0.01".
"""

from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from scipy.spatial import cKDTree

TILE_ID_RE = re.compile(r"^tile_RA(?P<ra>[0-9.]+)_DEC(?P<sign>[pm])(?P<dec>[0-9.]+)$")
SIGNAL_SIGMA = 3.0
SAME_DET_ARCSEC = 2.0


def uv(ra, dec) -> np.ndarray:
    ra = np.radians(np.asarray(ra, dtype=np.float64))
    dec = np.radians(np.asarray(dec, dtype=np.float64))
    c = np.cos(dec)
    return np.column_stack((c * np.cos(ra), c * np.sin(ra), np.sin(dec)))


def sep_arcsec(d) -> np.ndarray:
    return np.degrees(2.0 * np.arcsin(np.clip(np.asarray(d) / 2.0, -1.0, 1.0))) * 3600.0


def peak_sigma(path, ra, dec, r_src=2.0, r_in=8.0, r_out=20.0):
    with fits.open(path, memmap=True) as h:
        hdu = h[0]
        w = WCS(hdu.header)
        try:
            x, y = w.all_world2pix(ra, dec, 0)
        except Exception:
            x, y = w.wcs_world2pix(ra, dec, 0)
        x, y = float(x), float(y)
        sc = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600
        R = int(r_out / sc) + 3
        ny, nx = hdu.shape
        x0, x1 = max(0, int(x) - R), min(nx, int(x) + R)
        y0, y1 = max(0, int(y) - R), min(ny, int(y) + R)
        if (x1 - x0) < 2 * R or (y1 - y0) < 2 * R:
            return np.nan  # clipped at the frame edge; annulus would be biased
        d = np.asarray(hdu.section[y0:y1, x0:x1], float)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    rr = np.hypot(xx - x, yy - y) * sc
    ann, core = d[(rr >= r_in) & (rr <= r_out)], d[rr <= r_src]
    if ann.size < 20 or core.size < 1:
        return np.nan
    bg = np.median(ann)
    sig = 1.4826 * np.median(np.abs(ann - bg))
    return (core.max() - bg) / sig if sig > 0 else np.nan


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--buckets", default="BC")
    ap.add_argument("--ref-tiles", required=True)
    ap.add_argument("--plate-dir", required=True)
    ap.add_argument("--radec-dir", required=True)
    ap.add_argument("--our-manifest", required=True,
                    help="released tile_manifest.csv[.gz]; defines which plates cover a position")
    ap.add_argument("--control-n", type=int, default=150)
    ap.add_argument("--box-arcmin", type=float, default=1.5)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.rows)
    df = df[df["bucket"].astype(str).str[0].isin(list(args.buckets))].reset_index(drop=True)
    print(f"{len(df)} rows in buckets {args.buckets}")

    # ---- our nearest detection, with its coordinates this time ----
    # Must search EVERY plate whose tiles cover the position, exactly as arm B
    # did -- 14 of these rows are covered by 2-4 plates and the nearest
    # detection is not necessarily on the nearest plate. Searching only
    # `our_nearest_plate` disagreed with the parity run on 9 rows.
    man = pd.read_csv(args.our_manifest, usecols=["tile_id", "plate_id"])
    t_ra, t_dec = [], []
    for t in man["tile_id"]:
        m = TILE_ID_RE.match(str(t).strip())
        t_ra.append(float(m.group("ra")))
        d = float(m.group("dec"))
        t_dec.append(d if m.group("sign") == "p" else -d)
    tile_tree = cKDTree(uv(t_ra, t_dec))
    near = tile_tree.query_ball_point(uv(df["ra"], df["dec"]),
                                      r=2 * np.sin(np.radians(0.75) / 2))
    plate_arr = man["plate_id"].to_numpy()
    by_plate: dict[str, list[int]] = {}
    for i, ts in enumerate(near):
        for p in set(plate_arr[ts]):
            by_plate.setdefault(str(p), []).append(i)

    best = np.full(len(df), np.inf)
    bra = np.full(len(df), np.nan)
    bdec = np.full(len(df), np.nan)
    bplate = np.array([""] * len(df), dtype=object)
    for plate, idx in sorted(by_plate.items()):
        f = Path(args.radec_dir) / f"{plate}.csv"
        if not f.exists():
            continue
        det = pd.read_csv(f, usecols=["ra", "dec"])
        if det.empty:
            continue
        idx = np.asarray(idx)
        d, j = cKDTree(uv(det["ra"].to_numpy(), det["dec"].to_numpy())).query(
            uv(df["ra"].to_numpy()[idx], df["dec"].to_numpy()[idx]), k=1)
        s = sep_arcsec(d)
        better = s < best[idx]
        sel = idx[better]
        best[sel] = s[better]
        bra[sel] = det["ra"].to_numpy()[j[better]]
        bdec[sel] = det["dec"].to_numpy()[j[better]]
        bplate[sel] = plate
    df = df.assign(our_det_ra=bra, our_det_dec=bdec, our_det_sep=best, our_det_plate=bplate)

    bad = (df["our_det_sep"] - df["dist_raw_arcsec"]).abs() > 0.01
    if bad.any():
        print(f"[FATAL] nearest-detection distance disagrees with the parity run on "
              f"{int(bad.sum())} rows", file=sys.stderr)
        return 2
    print("nearest-detection distances reproduce dist_raw_arcsec exactly")

    # ---- pixels + the reference run's own extraction ----
    COLS = ["NUMBER", "RA_corr", "Dec_corr"]
    out = []
    for i, r in df.iterrows():
        # "Is there signal in OUR data" must consider every plate covering the
        # position, not the plate our nearest detection happened to come from
        # -- those differ on 14 of 45 rows, and measuring on a plate that
        # covers the position only at its vignetted edge invents a non-detection.
        cover = sorted({str(p) for p in plate_arr[near[i]]})
        our_sig = np.nan
        for p in cover:
            pf = Path(args.plate_dir) / f"dss1red_{p}.fits"
            if not pf.exists():
                continue
            s = peak_sigma(str(pf), r["ra"], r["dec"])
            if np.isfinite(s) and (not np.isfinite(our_sig) or s > our_sig):
                our_sig, plate = s, p
        if not np.isfinite(our_sig):
            plate = str(r["our_det_plate"]) or str(r["our_nearest_plate"])

        cand = glob.glob(str(Path(args.ref_tiles) / str(r["ref_tile_id"]) / "raw" / "*.fits"))
        ref_sig = peak_sigma(cand[0], r["ra"], r["dec"]) if cand else np.nan

        # RA_corr/Dec_corr are written by the wcsfix stage, so the WCS-fixed
        # coordinates live in pass2.wcsfix.csv; plain pass2.csv carries only
        # ALPHA_J2000/DELTA_J2000. Prefer the fixed frame, since that is what
        # the reference S0 and our S0 both use.
        cats = Path(args.ref_tiles) / str(r["ref_tile_id"]) / "catalogs"
        ref_at_ours = np.nan
        for fname, racol, deccol in (("sextractor_pass2.wcsfix.csv", "RA_corr", "Dec_corr"),
                                     ("sextractor_pass2.csv", "ALPHA_J2000", "DELTA_J2000")):
            p2 = cats / fname
            if not p2.exists():
                continue
            want = ("NUMBER", racol, deccol)
            rd = pd.read_csv(p2, usecols=lambda c: c in want)   # usecols: VIGNET is 65 MB/tile
            if racol not in rd.columns:
                continue
            rd = rd.dropna(subset=[racol, deccol])
            if len(rd):
                t = cKDTree(uv(rd[racol].to_numpy(), rd[deccol].to_numpy()))
                dd, _ = t.query(uv([r["our_det_ra"]], [r["our_det_dec"]]), k=1)
                ref_at_ours = float(sep_arcsec(dd)[0])
            break

        sig_here = (our_sig >= SIGNAL_SIGMA)
        ref_has_own = np.isfinite(ref_at_ours) and ref_at_ours <= SAME_DET_ARCSEC
        if not sig_here and not (ref_sig >= SIGNAL_SIGMA):
            cls = "NOSIGNAL"
        elif not sig_here:
            cls = "DISPLACED"
        elif ref_has_own:
            cls = "SEPARATE"
        else:
            cls = "MISSED"
        out.append(dict(src_id=r["src_id"], bucket=str(r["bucket"])[0], ra=r["ra"], dec=r["dec"],
                        ref_plate=r["ref_plate"], dist_raw_arcsec=r["dist_raw_arcsec"],
                        our_peak_sigma=our_sig, ref_peak_sigma=ref_sig,
                        ref_det_near_our_det_arcsec=ref_at_ours, classification=cls,
                        our_det_ra=r["our_det_ra"], our_det_dec=r["our_det_dec"],
                        our_det_plate=r["our_det_plate"], our_sig_plate=plate))
        print(f"  {i+1}/{len(df)} {r['src_id'][:38]:<38} {cls}", flush=True)

    res = pd.DataFrame(out)
    res.to_csv(out_dir / "classification.csv", index=False)

    print("\n=== CLASSIFICATION ===")
    for k, v in res["classification"].value_counts().items():
        print(f"  {k:<10} {v:>3}  {100*v/len(res):5.1f}%")
    print("\n=== by bucket ===")
    print(pd.crosstab(res["bucket"], res["classification"]).to_string())
    print("\n=== peak significance at the reference position ===")
    print(res.groupby("classification")[["our_peak_sigma", "ref_peak_sigma", "dist_raw_arcsec"]]
          .median().round(2).to_string())
    print(f"\npixel-source agreement: correlation "
          f"{res['our_peak_sigma'].corr(res['ref_peak_sigma']):.4f}")

    # ---- built-in control: the estimator must light up on rows we DID reproduce ----
    allrows = pd.read_csv(args.rows if "rows.csv" in args.rows else
                          str(Path(args.rows).parent / "rows.csv"))
    ctrl = allrows[(allrows.in_footprint) & (allrows.dist_s0_arcsec <= 5)]
    if len(ctrl) > args.control_n:
        ctrl = ctrl.sample(args.control_n, random_state=16)
    cv = []
    for _, q in ctrl.iterrows():
        c = glob.glob(str(Path(args.ref_tiles) / str(q.ref_tile_id) / "raw" / "*.fits"))
        if c:
            cv.append(peak_sigma(c[0], q.ra, q.dec))
    cv = pd.Series(cv).dropna()
    print(f"\n=== CONTROL: {len(cv)} reference rows we DID reproduce, same estimator ===")
    print(f"  median {cv.median():.1f} sigma, {100*(cv>SIGNAL_SIGMA).mean():.1f}% above "
          f"{SIGNAL_SIGMA:g} sigma   (vs {res.ref_peak_sigma.median():.1f} / "
          f"{100*(res.ref_peak_sigma>SIGNAL_SIGMA).mean():.1f}% for these rows)")
    if cv.median() < 10:
        print("  [WARN] control is weak -- the estimator, not the data, may be the problem")

    # ---- compact grid, our pixels, both positions marked ----
    n = len(res)
    ncol, nrow = 5, int(np.ceil(n / 5))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.3 * ncol, 2.55 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax in axes[n:]:
        ax.axis("off")
    for k, r in res.iterrows():
        ax = axes[k]
        pf = Path(args.plate_dir) / f"dss1red_{r['our_det_plate']}.fits"
        ax.set_xticks([]); ax.set_yticks([])
        if not pf.exists():
            ax.axis("off"); continue
        with fits.open(pf, memmap=True) as h:
            hdu = h[0]; w = WCS(hdu.header)
            x, y = w.all_world2pix(r["ra"], r["dec"], 0)
            sc = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600
            R = int(args.box_arcmin * 60 / 2 / sc)
            ny, nx = hdu.shape
            x0, x1 = max(0, int(x) - R), min(nx, int(x) + R)
            y0, y1 = max(0, int(y) - R), min(ny, int(y) + R)
            d = np.asarray(hdu.section[y0:y1, x0:x1], float)
            sub = w[y0:y1, x0:x1]
        lo, hi = ZScaleInterval().get_limits(d)
        ax.imshow(d, origin="lower", cmap="gray", vmin=lo, vmax=hi, interpolation="nearest")
        cx, cy = sub.all_world2pix(r["ra"], r["dec"], 0)
        ax.add_patch(plt.Circle((float(cx), float(cy)), 4, fill=False, color="#ff3b30", lw=1.3))
        ox, oy = sub.all_world2pix(r["our_det_ra"], r["our_det_dec"], 0)
        ax.plot([float(ox)], [float(oy)], "+", color="#0a84ff", ms=11, mew=1.6)
        ax.set_title(f"{r['classification']}  {r['dist_raw_arcsec']:.1f}\"  "
                     f"{r['our_peak_sigma']:.1f}$\\sigma$", fontsize=7)
    fig.suptitle(f"Reference rows with our nearest detection at 5-30\" "
                 f"({args.box_arcmin:g}' stamps, our plate scan; red = reference position, "
                 f"blue + = our nearest detection)", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    png = out_dir / "displaced_grid.png"
    fig.savefig(png, dpi=140)
    print(f"\nwrote {out_dir}/classification.csv and {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
