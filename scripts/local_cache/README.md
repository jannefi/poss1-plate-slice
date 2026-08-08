# scripts/local_cache/

Offline replication tooling for building local Parquet caches of the
large external catalogs VASCO60 queries during pipeline runs. Each
subdirectory builds one catalog's cache from its upstream bulk dump.

## Why

The pipeline normally fetches per-tile neighbourhoods from VizieR/TAP
and CDS endpoints. For the common vetoes and cross-matches, this hits
row caps in dense fields and has highly variable tail latency. A local
HEALPix-partitioned Parquet cache fixes both: no caps, no network, no
tail, ~20-50x faster per tile for the vetoes we actually run.

## Catalogs

| Catalog | Status | Source | Cache dir default |
|---|---|---|---|
| [Gaia DR3](gaia/) | done | `http://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/` | `/Volumes/SANDISK/Gaia` |
| [PS1 DR2](ps1/) | done | STScI HATS `s3://stpubdata/panstarrs/ps1/public/hats/otmo/` | `/Volumes/SANDISK/PS1` |
| [USNO-B 1.0](usnob/) | done | VizieR TAP `I/284/out` | `/Volumes/SANDISK/USNOB` |

## Conventions for adding a new catalog

1. **One subdirectory per catalog** (`gaia/`, `ps1/`, `usnob/`). Each is
   self-contained: build script, verify script, optional bench and
   analysis scripts, and a per-catalog README with the replication recipe.
2. **Paths are configurable** via `--cache-dir` flag with an environment
   variable fallback (`VASCO_<NAME>_CACHE`) and a hardcoded default
   matching the canonical build host. Scripts must run with zero
   arguments on the canonical host and run anywhere else by setting the
   env var.
3. **Resumable builds.** Writes a per-upstream-file `.done` marker so
   re-runs after an interrupt pick up where they left off.
4. **MD5 verification** against the upstream checksum manifest before
   any bytes are accepted into the cache.
5. **Partitioning by HEALPix level 5** (12,288 pixels, ~1.8° across) is
   the default for source catalogs used in tile-scale veto / cross-match
   work. Justified in `gaia/README.md`.
6. **Parquet with zstd** is the default on-disk format.
7. **No pipeline code changes live here.** These scripts produce data
   that other code reads. Wiring a cache into `vasco/` is a separate
   change done under the regular Plan Mode protocol.
8. **License and citation:** each cache writes a `LICENSE.txt` and a
   `MANIFEST.json` with build provenance, schema, row count, and the
   upstream attribution required by that catalog's license.

## What does NOT belong here

- Anything that runs as part of a normal pipeline run — that's `vasco/`.
- Per-tile neighbourhood fetch code — that's `vasco/external_fetch_online.py`.
- Pre-warming scripts that walk already-fetched tile caches —
  that's `scripts/prewarm_neighbourhood_cache.py`.
- Replication of specific VASCO60 results (plots, rejection renders) —
  that's `scripts/replication/`.
