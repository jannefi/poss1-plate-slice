# scripts/local_cache/gaia/

Build a local, HEALPix-partitioned, 6-column Parquet cache of Gaia DR3
from ESA's public CDN bulk dump. Replaces per-tile VizieR cone-search
queries in the VASCO60 pipeline for the veto and cross-match paths.

## What the cache holds

- **1,811,709,771 rows** — every Gaia DR3 source with a `source_id`.
- **6 columns** — `source_id`, `ra`, `dec`, `pmra`, `pmdec`, `phot_g_mean_mag`
  (plus a derived `healpix_5` partition key).
- **~60 GB on disk** as zstd Parquet, partitioned Hive-style by HEALPix
  nested level 5 (12,288 partitions, 1-8 files each, ~15,600 files total).
- Epoch: **J2016.0** (Gaia DR3 reference epoch). Apply `pmra/pmdec` to
  bring positions to another epoch.

## Precision note

`pmra` and `pmdec` are stored as **float32**, not float64. This is a
deliberate trade:

- float32 ULP at 100 mas/yr = 7.6e-6 mas/yr
- float32 ULP at 10,000 mas/yr (Barnard's star) = 9.8e-4 mas/yr
- Gaia's own published pm uncertainty: ~0.02-1.5 mas/yr

The representation error is 4-7 orders of magnitude below Gaia's own
measurement uncertainty at every realistic pm magnitude. Re-run
`pm_precision.py` to confirm the numbers against this machine's cache.

If you need bit-exact pm for regression diffing against VizieR, change
the two `pa.float32()` entries for pmra/pmdec in `build_cache.py` to
`pa.float64()`, delete `parquet/` and `.done/`, and re-run. Adds ~10 GB
to the cache.

## Replication recipe

Prerequisites: Python 3.11+, `pyarrow`, `numpy`, `pandas`,
`astropy`, `astropy-healpix`. On the canonical build host (macOS,
`/Volumes/SANDISK` mounted) no arguments are needed. On any other host,
set `VASCO_GAIA_CACHE` and `VASCO_GAIA_STAGING` first.

```
# 0. (optional) set cache and staging dirs for this host
export VASCO_GAIA_CACHE=/path/to/gaia_cache
export VASCO_GAIA_STAGING=/path/to/staging   # ~1 TB free is plenty; 700 MB peak

# 1. build
python3 scripts/local_cache/gaia/build_cache.py --workers 3

# 2. verify (14 defence-in-depth checks; no network: add --skip-network)
python3 scripts/local_cache/gaia/verify_cache.py

# 3. confirm float32 pm precision claim empirically
python3 scripts/local_cache/gaia/pm_precision.py

# 4. benchmark against VizieR on 10 random tiles (~90 s of VizieR load)
python3 scripts/local_cache/gaia/bench_vs_vizier.py
```

Expected outcomes:

| Step | Expected result |
|---|---|
| build | `DONE ok=3386 fail=0 rows=1,811,709,771 elapsed≈12000s` (3 workers, ~3.5 hr on a typical broadband link) |
| verify | `TOTAL: 14 pass, 0 fail` |
| pm_precision | `Fraction of rows where float32 ULP > 0.01 mas/yr: 0/1,000,000 = 0.0000%` |
| bench | `Match: 10/10 tiles ok`, `Overall speedup: ~35x` |

## Known cosmetic issue: `_MD5SUM.txt` FAIL during build

ESA's upstream `_MD5SUM.txt` lists itself as one of its own entries,
with a hash that cannot mathematically match the file's content (a file
cannot contain its own MD5 — any value embedded changes the file, so no
hash you put in is the hash of the file you just wrote). The build
script now filters the self-reference from its work queue, so this
never appears in a fresh run. If you see it in an old log, `verify_cache.py`
Check C proves it is cosmetic.

## Build outputs on disk

```
$VASCO_GAIA_CACHE/
  parquet/healpix_5=<0..12287>/*.parquet    60 GB    the cache itself
  metadata/                                 ~1 MB   ESA _MD5SUM.txt + _license + _citation + _disclaimer + _readme
  .done/                                    ~few MB per-file ingest markers (resumable)
  MANIFEST.json                                     build provenance + schema + precision section
  LICENSE.txt                                       CC BY 4.0 + full ESA attribution text
  progress.log                                      human-readable build log
```

`MANIFEST.json` and `LICENSE.txt` are not written by `build_cache.py`
itself (they were drafted once at cache creation and live in the cache
directory). If you build from scratch on a new host you should copy
them from the canonical cache or regenerate them; the build process
does not depend on them existing.

## Upstream

- Bulk dump: http://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/
- Archive (authoritative): http://archives.esac.esa.int/gaia
- License: https://www.cosmos.esa.int/web/gaia-users/license — CC BY 4.0
- Citation: https://www.cosmos.esa.int/web/gaia-users/credits
