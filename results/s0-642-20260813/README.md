# S0 — POSS-I candidate catalogue, 642 plates (2026-08-13)

**135,066 rows.** Detections on POSS-I red plates that survive the MNRAS 2022
filter chain, the Gaia / PS1 / USNO-B 5″ vetoes, the diffraction-spike mask, and
global deduplication at 0.25″.

Produced entirely from public inputs. Nothing here derives from any unpublished
catalogue — see [`tools/audit_independence.py`](../../tools/audit_independence.py).

## Files

| file | rows | what |
|---|---:|---|
| `stage_S0.csv.gz` | 135,066 | the catalogue |
| `tile_manifest.csv.gz` | 31,458 | every tile processed, with its plate |
| `RUN_SUMMARY.txt` | — | parameters and counts as the build recorded them |
| `verification_s0_gaia_invariant.csv.gz` | 25,642 | per-tile Gaia-contamination check |
| `verification_dedup_radius_sweep.json` | — | the measurement behind the 0.25″ tolerance |
| `known_astrometry_defect_tiles.csv` | 14 | tiles with a known ~12″ position error — exclude from crossmatches |
| `SHA256SUMS` | — | integrity of the files as shipped |
| `SHA256SUMS.uncompressed` | — | integrity of the *contents* |

`stage_S0.csv` columns: `src_id, tile_id, object_id, ra, dec`. Positions are
ICRS degrees on the DSS/GSSS plate solution with the per-plate CRPIX correction
applied (`docs/DSS_WCS_TWO_SOLUTIONS.md`).

### Content hashes — quote these, not the `.gz` ones

SHA-256 of the **decompressed** files:

| file | rows | bytes | sha256 |
|---|---:|---:|---|
| `stage_S0.csv` | 135,066 | 12,103,728 | `9c788c30cd7c9c16ef99d3b6184a0aca27c385bf057f0b10229d879438d4bc73` |
| `tile_manifest.csv` | 31,458 | 3,882,323 | `94c342cd2cc0dbd29171a5096d3d083579001952e0d0a96a1587a50bc9663563` |
| `verification_s0_gaia_invariant.csv` | 25,642 | 864,546 | `486ce8c50b3faaaa4b552e3f0d949ec0cef1ce098707eb18b62da118306b4dc4` |

**Why both.** `SHA256SUMS` covers the `.gz` files and verifies the *transfer*.
It is not a stable identifier of the data: gzip output depends on the
implementation, the compression level and an embedded timestamp, so recompressing
these exact bytes on this same machine already produces a different `.gz`. Anyone
who repacks, mirrors, or regenerates the archive will fail that check while
holding identical data. The uncompressed hashes above are what identify the
catalogue, and they are what a citation should pin.

```bash
sha256sum -c SHA256SUMS               # the files as shipped
zcat stage_S0.csv.gz | sha256sum      # must equal the table above
```

## How it was produced

642 POSS-I red plates (centre δ ≥ −3.0°) from IRSA's `dss1red` library, sliced
locally into 7×7 tiles per plate, then SExtractor two-pass with PSFEx, the
Solano et al. (2022) §2 filters, and the veto chain. Full recipe:
[`docs/REPRODUCING.md`](../../docs/REPRODUCING.md). Every threshold and its
provenance: [`docs/PARAMETERS.md`](../../docs/PARAMETERS.md).

Three implemented stages — SkyBoT, SuperCOSMOS, VSX — were **not run**. In the
predecessor pipeline they removed ~0.6% of survivors combined, and SkyBoT and
VSX removed nothing at all.

## This catalogue is a repair. Read this before using it.

An earlier build of the same run gave **310,700 rows and was wrong**.

The mirror cone query selected HEALPix pixels with `cone_search_lonlat`, which
returns only pixels whose *centres* fall inside the radius. At nside=32 a pixel
spans ~1.8° against a ~0.76° veto cone, so pixels overlapping the cone with
their centre outside were silently dropped. The veto catalogue was partial and
real stars survived into the candidate list. Nothing errored — the stage ran,
wrote its crossmatch file, and reported success.

It acts only where the veto cone crosses the HEALPix polar-cap boundary at
δ = arcsin(2/3) = 41.81°. **2.18× inflation; 56% of the old catalogue was
un-vetoed Gaia stars.**

Four things worth knowing about the repair, because you cannot check them from
the code alone:

1. **It was found by chasing an unexplained excess, not by a test.** The
   standing check that would have caught it on day one
   (`tools/check_s0_gaia_invariant.py`) was written afterwards. It reports
   10.02% on the old build and **0.02%** on this one — that ledger is included
   above so you can verify the claim rather than take it.
2. **The cheap repair failed.** Re-vetoing the affected rows in place does not
   work: the morphology filter derives its FWHM window from the population being
   filtered, so the extra stars narrowed that window and the bug *removed* real
   rows too. Those cannot be recovered post hoc. 504 tiles were fully re-run.
3. **Two obvious scoping shortcuts are wrong.** Declination alone is not
   sufficient (1,503 tiles above the boundary were provably unaffected), and
   "re-vetoing this tile changed nothing" is not sufficient either — it would
   have left 90 corrupt tiles in place.
4. **An independent catalogue corroborates the repair.** USNO-B, which played no
   part in it, finds POSS-II second-epoch counterparts for **98.94%** of the
   removed rows (7.9× above its own null-shifted chance rate) versus **0.09%**
   of the rows retained (95× below chance), agreeing with the pipeline's own
   98.76% removal to 0.18 points. That signal was present and recorded as
   unexplained on earlier runs, and was not followed up.

**Residual uncertainty**: 90 tiles holding ~504 rows sit in a class the scoping
check flags but the re-veto could not confirm. Their contribution could move in
either direction.

## Known defect: 9.8% of rows have unreliable coordinates

**14 tiles carry 13,270 rows (9.8%) whose positions are wrong by roughly 12″.**
They are listed in `known_astrometry_defect_tiles.csv`. **Exclude them from any
positional crossmatch**, or you will measure this defect rather than whatever
you set out to measure.

The detections themselves are real. On the affected tiles they are uniformly
distributed, sit on genuine flux peaks, and are indistinguishable in FWHM (2.60
vs 2.61 px) and instrumental magnitude (8.90 vs 8.85) from the rows that behave
normally. What fails is the coordinate solution: on a control tile the brightest
Gaia stars land a median **2.24 px** from a flux peak against **7.07 px** for a
deliberately-shifted null, while on an affected tile they land **7.14 px** away
against a **7.21 px** null — i.e. the astrometry there is statistically
indistinguishable from random.

That is also why these rows survived the pipeline in such numbers. The vetoes
match Gaia, PS1 and USNO-B within 5″; a row displaced by ~12″ matches nothing,
so nothing removes it. A normal tile yields ~3 rows here, these yield hundreds
to thousands.

**Ruled out** as causes: the partial-cone veto bug (their cone queries were
complete), plate-edge geometry (their distance from plate centre matches the
control distribution), the plate solutions themselves (other tiles on the same
plates are normal, median 1–13 rows), the tile WCS refit (residual 0.10″,
identical to healthy plates), and spurious detections. The cause is **not yet
identified**, and the rows are expected to be **repairable** rather than
discardable — most are probably ordinary stars that the vetoes will remove once
the positions are corrected, which would make a corrected catalogue *smaller*
than this one.

Being explicit because it is the first thing this defect will affect: any
crossmatch against another catalogue will be depressed by these rows, and the
size of that effect on your numbers depends on your match radius.

## What this catalogue is not

A row here is a source that is **on the plate and not in the modern catalogues
we checked**. That is a statement about two catalogues, not about the universe.
It is not a vanishing-object list, not a transient list, and no row should be
called a detection of anything until it has been inspected individually.

Deduplication is global at 0.25″, which is the tolerance this project pairs with
WCS-fixed coordinates. The included sweep shows that choice is **conservative**:
it leaves ~258 real cross-tile duplicate pairs unmerged, while merging genuinely
distinct sources does not begin until 4″. A 3.0″ tolerance would give 134,821
rows — a 0.18% difference.

## Verify it

```bash
sha256sum -c SHA256SUMS
VASCO_GAIA_CACHE=<gaia_mirror> python3 ../../tools/check_s0_gaia_invariant.py \
    --s0-csv <(zcat stage_S0.csv.gz) --out /tmp/s0_gaia.csv
```

The invariant must PASS: the fraction of rows with a Gaia source within 5″ must
be ~0, because the veto removed exactly those. If it fails, the veto did not
cover that sky — whatever the logs say.

## Licence and attribution

Catalogue: see [`LICENSE-DATA`](../../LICENSE-DATA). Plate scans are DSS/POSS-I
via IRSA; the acknowledgements those require are in
[`ACKNOWLEDGMENTS.md`](../../ACKNOWLEDGMENTS.md) and must travel with any reuse.
