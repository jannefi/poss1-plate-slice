# scripts/local_cache/ps1/

Build a local, HEALPix-partitioned, 9-column Parquet cache of Pan-STARRS1
DR2 from STScI's public AWS HATS dataset.
Replaces per-tile VizieR `II/389/ps1_dr2` cone-search queries in the
VASCO60 veto and cross-match path.

## About the `otmo` table name

`otmo` is the name STScI/LSDB gave to the HATS directory on AWS
(`s3://stpubdata/panstarrs/ps1/public/hats/otmo/`). It is **not** a
standard name in STScI's documented PS1 DR2 schema (which names the
object tables `ObjectThin`, `MeanObject`, `StackObjectThin`, etc.).

From the schema — position + uncertainties + `nDetections` from
ObjectThin, plus the full set of per-band Mean*PSFMag / Mean*KronMag /
Mean*ApMag fields from MeanObject — it is effectively **ObjectThin
joined to MeanObject on objID**. The abbreviation therefore appears to
stand for **O**bject**T**hin + **M**ean**O**bject. The HATS team has
not documented this in a README on the bucket, so this reading is
inferred from the schema, not from an official source.

## How the cache compares to VizieR II/389/ps1_dr2

The **cache is a strict superset of VizieR's II/389/ps1_dr2**. In a
5-tile bit-exact comparison, VizieR's row set = `otmo` filtered to
`nDetections >= 3`. Objects with 1 or 2 detections are present in
`otmo` and in our cache, but absent from VizieR.

Counts from one test tile at (RA 137.355, Dec +43.418), 30′ cone:

| source | rows | filter |
|---|---|---|
| cache (unfiltered) | 237,204 | none |
| cache (`nDetections >= 3`) | **17,768** | `nDet >= 3` |
| VizieR `II/389/ps1_dr2` | 17,768 | (server-side, opaque to caller) |

Symmetric difference on objID sets: **0** on all 5 tiles tested.

**Consequences for downstream use:**

- **For a VizieR-compatible query** — apply `nDetections >= 3` to every
  cache query. `bench_vs_vizier.py` does this by default so its
  comparison against the live VizieR service is apples-to-apples.
- **For a more-inclusive veto path** (probably desirable for VASCO60) —
  omit the filter and use the full cache. This will reject more
  candidates than the current VizieR-backed veto, including ones
  supported only by 1 or 2 PS1 detections. That's arguably more
  conservative and safer, but is a semantic change from the current
  pipeline behaviour, and should be noted if you wire the cache in
  unfiltered.
- **For raw-detection counting / differential-imaging style work** —
  use unfiltered; you want the full detection list.

The cache holds both populations (via the `nDetections` column), so
either semantic is one line of query away.

## What the cache holds

- **10,560,724,292 rows** — every PS1 DR2 otmo source.
- **9 columns** — `objID`, `ra`, `dec`, `nDetections`, `gmag`, `rmag`,
  `imag`, `zmag`, `ymag` (plus a derived `healpix_5` partition key).
- **~180–250 GB on disk** as zstd Parquet, partitioned Hive-style by
  HEALPix nested level 5 (12,288 partitions, identical layout to the
  Gaia cache).

## Source

STScI's Pan-STARRS DR2 HATS mirror, hosted anonymously on AWS:

```
s3://stpubdata/panstarrs/ps1/public/hats/otmo/otmo/dataset/
```

- Native format: HATS (HEALPix Adaptive Tiling Scheme) Parquet, with
  adaptive partition orders 2/5/6/7 (27,161 leaf files, ~3 TB total).
- Built by the PS1/STScI team using `hats-import`, published 2024-09-18.
- All 9 columns VASCO60 needs are present in the source schema under
  STScI's native names (`raMean`, `decMean`, `nDetections`,
  `gMeanPSFMag`, …). The build script projects and renames them.

Key insight: **we never download the full 3 TB**. Each HATS leaf
parquet file has 135 columns, but we use `pyarrow.fs.S3FileSystem` +
Parquet column projection to read only the 9 columns we want. This
transfers ~9% of the source bytes (~274 GB over the network) instead of
the full 3 TB. No local staging — each file is streamed directly from
S3 through the transform into the local partitioned output.

## Column name mapping

| STScI native (otmo) | VizieR (II/389) | Local cache |
|---|---|---|
| `objID` | `objID` | `objID` |
| `raMean` | `RAJ2000` | `ra` |
| `decMean` | `DEJ2000` | `dec` |
| `nDetections` | `Nd` | `nDetections` |
| `gMeanPSFMag` | `gmag` | `gmag` |
| `rMeanPSFMag` | `rmag` | `rmag` |
| `iMeanPSFMag` | `imag` | `imag` |
| `zMeanPSFMag` | `zmag` | `zmag` |
| `yMeanPSFMag` | `ymag` | `ymag` |

## Precision and sentinel values

- `ra`/`dec`: stored as **float64** (matching Gaia cache for a
  uniform position type across catalogs).
- 5 PSF mags: stored as **float32**. float32 ULP at PS1 mag ranges
  (10-25 mag) is ~1e-6 mag, far below PS1's own photometric precision
  (~5-50 mmag for well-measured sources).
- `nDetections`: **int16** (PS1 publishes as int16; real values typically
  1-400).
- `objID`: **int64** (PS1's native 64-bit object ID).
- **PS1's `-999` sentinel is preserved as-is** for sources with no
  detection in a given band. If you need NaN semantics, filter at
  query time.

## Partitioning

Unlike the HATS source (adaptive Norder 2-7), the local cache is
flattened to **fixed HEALPix level 5 nested** — 12,288 partitions,
1.83° across each, matching the Gaia cache exactly. This is what makes
query code uniform across both caches.

The trade: dense-field HP5 partitions can be large (up to ~10 MB
compressed, ~3 M rows for extreme galactic-plane pixels), but each
partition still loads in well under a second for a cone search. The
dominant cost at query time is the final angular-distance filter, not
partition I/O.

**`healpix_5` is computed via `astropy_healpix.lonlat_to_healpix(ra,
dec)`** — PS1's `objID` does NOT encode HEALPix position (unlike
Gaia's `source_id`), so the bit-shift trick used in the Gaia build
does not apply here.

## Replication recipe

Prerequisites: Python 3.11+, `pyarrow` (with S3 support),
`numpy`, `pandas`, `astropy`, `astropy-healpix`. On the canonical build
host (macOS, `/Volumes/SANDISK` mounted), no arguments needed. On any
other host, set `VASCO_PS1_CACHE` first.

```
# 0. (optional) set cache dir for this host
export VASCO_PS1_CACHE=/path/to/ps1_cache

# 1. build
python3 scripts/local_cache/ps1/build_cache.py --workers 8

# 2. verify (6 cross-checks; no network: add --skip-network)
python3 scripts/local_cache/ps1/verify_cache.py

# 3. benchmark against VizieR on 10 random tiles
python3 scripts/local_cache/ps1/bench_vs_vizier.py
```

Expected outcomes:

| Step | Expected result |
|---|---|
| build | `DONE ok=27161 fail=0 rows=10,560,724,292` |
| verify | `TOTAL: 6 pass, 0 fail` |
| bench | `Match: 10/10 tiles ok`; dense-field tiles should show the cache returning *more* rows than VizieR's 200K cap |

## Runtime and bandwidth

- **Expected**: ~2 hours (my initial estimate based on benchmarking one
  S3 read in isolation).
- **Actual on a ~90 Mbps residential link**: ~8 hours. The S3 download
  throughput saturates home broadband before the 8 worker threads
  saturate the CPU. This is a **bandwidth-limited** workload; more
  workers on a slow pipe helps very little.
- If running from a cloud VM in `us-east-1`, expect ~20 minutes
  (same bucket region, multi-gigabit S3 bandwidth).
- **Runtime is not a concern for correctness** — the build is resumable
  via `.done` markers, so interruption just means continuing where you
  left off.

## Disk layout

```
$VASCO_PS1_CACHE/
  parquet/healpix_5=<0..12287>/*.parquet    180-250 GB    the cache itself
  .done/                                    per-source-file ingest markers (27,161 files, ~few MB)
  progress.log                              build log
  MANIFEST.json                             provenance + schema + row count (written post-build)
  LICENSE.txt                               STScI PS1 DR2 attribution (written post-build)
```

## Upstream and licensing

- Public AWS mirror: https://registry.opendata.aws/pan-starrs/ (STScI
  AWS Public Dataset — anonymous access)
- S3 location: `s3://stpubdata/panstarrs/ps1/public/`
- Original PS1 DR2 release: https://outerspace.stsci.edu/display/PANSTARRS/
- Pan-STARRS1 is a facility of the University of Hawaii Institute for
  Astronomy, STScI, and partners. Data products are publicly released.
- Citation: see https://outerspace.stsci.edu/display/PANSTARRS/PS1+News
  (Chambers et al. 2016 for the PS1 survey paper, Flewelling et al.
  2020 for DR2 processing).
