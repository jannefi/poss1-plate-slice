# scripts/local_cache/usnob/

Build a local, HEALPix-partitioned, 10-column Parquet cache of USNO-B1.0
(Monet et al. 2003) from VizieR TAP. Replaces per-tile VizieR
`I/284/out` cone-search queries in the VASCO60 veto / cross-match path.

## What the cache holds

- **~1,045,913,669 rows** — every USNO-B1.0 source.
- **10 columns** — `id` (USNO-B designation), `ra`, `dec`, `B1mag`,
  `R1mag`, `B2mag`, `R2mag`, `Imag`, `pmRA`, `pmDE`, plus the derived
  `healpix_5` partition key. These are the same 10 columns VASCO60's
  `fetch_usnob_neighbourhood` pulls today.
- **~25–40 GB on disk** as zstd Parquet, partitioned Hive-style by
  HEALPix nested level 5 (12,288 partitions, identical layout to the
  Gaia and PS1 caches).

## Source

**VizieR TAP** at `https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync`,
table `"I/284/out"` (the VizieR-mirrored copy of USNO-B1.0 described in
Monet et al. 2003, AJ, 125, 984).

Unlike Gaia (which has a static CDN bulk tree) and PS1 (which has a
ready-made HATS Parquet mirror on AWS), USNO-B has no canonical bulk
file distribution — the original USNO hosts at `nofs.navy.mil` are
offline, the CDS `ftp/I/284/` tree holds only the ReadMe and a 1000-row
sample, and STScI's AWS public bucket does not mirror it. The only live
bulk source is VizieR's own query service.

The build script paginates the catalog in **1800 declination zones of
0.1° each**, submits a synchronous ADQL query per zone that projects
only the 10 required columns, parses the CSV response, computes
`healpix_5` from the J2000 position, and writes one Parquet file per
zone into the local HP5-partitioned tree.

## Column name mapping

| VizieR native (I/284) | Local cache |
|---|---|
| `USNO-B1_0` (alias of `USNO-B1.0`) | `id` |
| `RAJ2000` | `ra` |
| `DEJ2000` | `dec` |
| `B1mag` | `B1mag` |
| `R1mag` | `R1mag` |
| `B2mag` | `B2mag` |
| `R2mag` | `R2mag` |
| `Imag` | `Imag` |
| `pmRA` | `pmRA` |
| `pmDE` | `pmDE` |

## Precision policy

- `ra`/`dec`: stored as **float64** (matching Gaia and PS1 caches for a
  uniform position type across catalogs).
- 5 mags: stored as **float32**. USNO-B's nominal photometric precision
  is ~0.3 mag — float32 ULP (<1e-5 mag at these ranges) is ~5 orders of
  magnitude below that.
- `pmRA`/`pmDE`: stored as **int16** matching USNO-B's published type
  (units mas/yr, range ±32,767 covers the full catalog).
- `id`: stored as UTF-8 string (12 characters). USNO-B encodes the dec
  zone in the first 4 characters, so string representation preserves
  the full catalog identifier.

## Replication recipe

Prerequisites: Python 3.11+, `pyarrow`, `numpy`, `pandas`, `astropy`,
`astropy-healpix`.

```
# 0. (optional) set cache dir for this host
export VASCO_USNOB_CACHE=/path/to/usnob_cache

# 1. build (~1800 dec-zone TAP queries, 4 workers, resumable)
python3 scripts/local_cache/usnob/build_cache.py --workers 4

# 2. verify
python3 scripts/local_cache/usnob/verify_cache.py

# 3. benchmark against VizieR on 10 random tiles
python3 scripts/local_cache/usnob/bench_vs_vizier.py
```

Expected:

| Step | Expected result |
|---|---|
| build | `DONE ok=1800 fail=0 rows≈1,045,913,669` |
| verify | `TOTAL: 5 pass, 0 fail` |
| bench | 10/10 strict bit-exact match on matched IDs |

## Runtime

VizieR TAP is a shared public service and the zone-by-zone query pattern
is sensitive to server load. At 4 workers, expected runtime is in the
**1–4 hour** range depending on CDS's current responsiveness. The build
is resumable via `.done` markers per zone, so interruption costs at
most one in-flight zone.

## Disk layout

```
$VASCO_USNOB_CACHE/
  parquet/healpix_5=<0..12287>/*.parquet    ~25-40 GB    the cache itself
  .done/                                    1,800 per-zone markers
  progress.log                              build log
  MANIFEST.json                             provenance + schema (post-build)
  LICENSE.txt                               Monet+ 2003 attribution (post-build)
```

## Upstream

- VizieR catalog page: https://cdsarc.cds.unistra.fr/viz-bin/cat/I/284
- Paper: Monet et al. 2003, AJ, 125, 984 — https://ui.adsabs.harvard.edu/abs/2003AJ....125..984M
- ReadMe: https://cdsarc.cds.unistra.fr/ftp/I/284/ReadMe

USNO-B1.0 is in the public domain. When used in publications, cite
Monet et al. 2003.
