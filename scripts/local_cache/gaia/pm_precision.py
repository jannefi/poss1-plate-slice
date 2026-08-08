"""
Empirical analysis of proper motion range and float32 precision loss in
the Gaia local Parquet cache.

Answers the question "is float32 storage of pmra/pmdec good enough?" by
comparing float32 ULP at realistic |pm| magnitudes against Gaia DR3's
own published pm uncertainty (~0.02 mas/yr bright, ~1.5 mas/yr faint).

Sample: ~2 million rows drawn from 16 random HEALPix-5 partitions.

Defaults:
    --cache-dir / $VASCO_GAIA_CACHE / /Volumes/SANDISK/Gaia
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pyarrow.compute as pc
import pyarrow.dataset as pds


def default_cache_dir() -> Path:
    return Path(os.environ.get("VASCO_GAIA_CACHE", "/Volumes/SANDISK/Gaia"))


def f32_ulp(x: float) -> float:
    x32 = np.float32(x)
    return float(
        np.abs(np.nextafter(x32, np.float32(np.inf), dtype=np.float32) - x32)
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Empirical float32 precision analysis for Gaia pmra/pmdec."
    )
    parser.add_argument("--cache-dir", type=Path, default=default_cache_dir())
    parser.add_argument("--n-partitions", type=int, default=16)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    ds = pds.dataset(
        str(args.cache_dir / "parquet"), format="parquet", partitioning="hive"
    )

    rng = np.random.default_rng(args.seed)
    pixels = rng.choice(12288, size=args.n_partitions, replace=False).tolist()
    tbl = ds.to_table(
        columns=["pmra", "pmdec"],
        filter=pc.field("healpix_5").isin([int(p) for p in pixels]),
    )

    pmra = tbl.column("pmra").to_numpy(zero_copy_only=False)
    pmde = tbl.column("pmdec").to_numpy(zero_copy_only=False)
    valid = ~(np.isnan(pmra) | np.isnan(pmde))
    pmra = pmra[valid]
    pmde = pmde[valid]
    pm = np.hypot(pmra, pmde)

    print(
        f"Sampled sources         : {len(pm):,} "
        f"(from {args.n_partitions} random HEALPix-5 partitions)"
    )
    print(
        f"pmra range              : {pmra.min():+12.3f}  .. "
        f"{pmra.max():+12.3f}  mas/yr"
    )
    print(
        f"pmdec range             : {pmde.min():+12.3f}  .. "
        f"{pmde.max():+12.3f}  mas/yr"
    )
    print(
        f"|pm| range              : {pm.min():12.3f}  .. "
        f"{pm.max():12.3f}  mas/yr"
    )
    print()
    print("|pm| distribution (mas/yr):")
    for p in [50, 90, 99, 99.9, 99.99, 100]:
        print(f"  {p:6.2f}th percentile : {np.percentile(pm, p):10.3f}")
    print()

    print("float32 ULP at representative |pm| magnitudes:")
    for x in [0.01, 0.1, 1.0, 10.0, 100.0, 1_000.0, 10_000.0]:
        print(
            f"  |pm|={x:10.3f} mas/yr  ->  float32 ULP = {f32_ulp(x):.3e} mas/yr"
        )
    print()

    sample_idx = rng.choice(len(pm), size=min(1_000_000, len(pm)), replace=False)
    sampled_pm = pm[sample_idx]
    for unc_thresh in [0.01, 0.05, 0.1, 1.0]:
        exceeds = sum(1 for x in sampled_pm if f32_ulp(float(x)) > unc_thresh)
        print(
            f"Fraction of rows where float32 ULP > {unc_thresh:6.3f} mas/yr: "
            f"{exceeds:,}/{len(sampled_pm):,} = "
            f"{exceeds / len(sampled_pm) * 100:.4f}%"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
