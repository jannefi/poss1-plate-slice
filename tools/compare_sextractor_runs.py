#!/usr/bin/env python3
import sys, math
from pathlib import Path
import pandas as pd

print("[compare] script loaded", flush=True)

# R points in your field (edit if needed)
R_POINTS = [
    ("1544", 130.44295, 32.908028),
    ("1546", 130.54050, 33.181760),
    ("1551", 130.23894, 32.718590),
    ("1554", 130.13702, 32.670372),
    ("1558", 130.14957, 32.877388),
    ("1579", 129.58102, 33.028930),
    ("1581", 129.51256, 32.996216),
]

def wrap_pm180(x):
    return (x + 180.0) % 360.0 - 180.0

def sep_arcsec(ra1, dec1, ra2, dec2):
    ra1r, dec1r = math.radians(ra1), math.radians(dec1)
    ra2r, dec2r = math.radians(ra2), math.radians(dec2)
    dra = math.radians(wrap_pm180(ra2 - ra1))
    ddec = dec2r - dec1r
    a = math.sin(ddec/2)**2 + math.cos(dec1r)*math.cos(dec2r)*math.sin(dra/2)**2
    c = 2*math.asin(min(1.0, math.sqrt(a)))
    return math.degrees(c) * 3600.0

def load_sex_csv(tile_dir: Path):
    p = tile_dir / "catalogs" / "sextractor_pass2.csv"
    if not p.exists():
        raise FileNotFoundError(f"Missing {p}")
    df = pd.read_csv(p)

    for ra_col, dec_col in [
        ("RA_corr","Dec_corr"),
        ("ALPHAWIN_J2000","DELTAWIN_J2000"),
        ("ALPHA_J2000","DELTA_J2000"),
        ("X_WORLD","Y_WORLD"),
    ]:
        if ra_col in df.columns and dec_col in df.columns:
            return df, ra_col, dec_col, p

    raise RuntimeError("No usable RA/Dec columns found in sextractor_pass2.csv")

def nearest(df, ra_col, dec_col, ra, dec):
    seps = df.apply(lambda r: sep_arcsec(ra, dec, float(r[ra_col]), float(r[dec_col])), axis=1)
    idx = seps.idxmin()
    return idx, float(seps.loc[idx])

def gate_row(row):
    ok = {}
    ok["flags"]  = (float(row.get("FLAGS", 1)) == 0)
    ok["snr"]    = (float(row.get("SNR_WIN", float("nan"))) > 30.0)
    ok["elong"]  = (float(row.get("ELONGATION", float("nan"))) < 1.3)
    fwhm = float(row.get("FWHM_IMAGE", float("nan")))
    ok["fwhm"]   = (fwhm > 2.0 and fwhm < 7.0)
    ok["spread"] = (float(row.get("SPREAD_MODEL", float("nan"))) > -0.002)
    return ok

def summarize(tile_dir: Path, label: str):
    df, ra_col, dec_col, path = load_sex_csv(tile_dir)

    print(f"\n== {label} ==", flush=True)
    print(f"Tile: {tile_dir}", flush=True)
    print(f"Catalog: {path}", flush=True)
    print(f"Using coords: {ra_col}/{dec_col}", flush=True)

    cols = ["FLAGS","SNR_WIN","ELONGATION","FWHM_IMAGE","SPREAD_MODEL"]
    present = [c for c in cols if c in df.columns]
    print("Columns present:", present, flush=True)

    if "SPREAD_MODEL" in df.columns:
        print("SPREAD_MODEL describe:", flush=True)
        print(df["SPREAD_MODEL"].describe(percentiles=[0.05,0.5,0.95]), flush=True)

    if "FWHM_IMAGE" in df.columns:
        print("FWHM_IMAGE describe:", flush=True)
        print(df["FWHM_IMAGE"].describe(percentiles=[0.05,0.5,0.95]), flush=True)

    rows_out = []
    for mid, ra, dec in R_POINTS:
        idx, sep = nearest(df, ra_col, dec_col, ra, dec)

        out = {"mid": mid, "sep_arcsec": sep}

        if sep > 5.0:
            out.update({
                "match_within5": False,
                "FLAGS": None,
                "SNR_WIN": None,
                "ELONGATION": None,
                "FWHM_IMAGE": None,
                "SPREAD_MODEL": None,
                "gate_all": False,
                "gate_flags": None,
                "gate_snr": None,
                "gate_elong": None,
                "gate_fwhm": None,
                "gate_spread": None,
                "note": 'NO within-5" match (nearest is farther)',
            })
        else:
            r = df.loc[idx]
            out.update({
                "match_within5": True,
                "FLAGS": r.get("FLAGS", None),
                "SNR_WIN": r.get("SNR_WIN", None),
                "ELONGATION": r.get("ELONGATION", None),
                "FWHM_IMAGE": r.get("FWHM_IMAGE", None),
                "SPREAD_MODEL": r.get("SPREAD_MODEL", None),
                "note": "",
            })
            g = gate_row(r)
            out["gate_flags"]  = g["flags"]
            out["gate_snr"]    = g["snr"]
            out["gate_elong"]  = g["elong"]
            out["gate_fwhm"]   = g["fwhm"]
            out["gate_spread"] = g["spread"]
            out["gate_all"]    = all(g.values())

        rows_out.append(out)

    outdf = pd.DataFrame(rows_out)
    print("\nR-nearest summary:", flush=True)
    print(outdf.sort_values("mid").to_string(index=False), flush=True)

def main():
    if len(sys.argv) != 3:
        print("Usage: compare_sextractor_runs.py <OLD_TILE_DIR> <NEW_TILE_DIR>", flush=True)
        return 2
    old_dir = Path(sys.argv[1])
    new_dir = Path(sys.argv[2])
    summarize(old_dir, "OLD")
    summarize(new_dir, "NEW")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
