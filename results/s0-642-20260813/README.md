# S0 — POSS-I candidate catalogue, 642 plates (2026-08-13)

> **SUPERSEDED 2026-08-14 by [`../s0-642-20260814/`](../s0-642-20260814/) — do not
> use this catalogue.** Eight tiles listed below as a known astrometry defect are
> now repaired: the cause was the tile WCS refit diverging where a tile crosses
> RA 0 near the pole. Those tiles held 12,273 rows here and hold 30 with correct
> astrometry, so the current catalogue is 122,820 rows. This folder is kept
> because its README carries the diagnosis and the two corrections that led to it.

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
| `known_astrometry_defect_tiles.csv` | 14 | 8 tiles with a ~11″ position error (`verdict=defect`) plus 6 previously mis-flagged and now cleared |
| `SHA256SUMS` | — | integrity of the files as shipped |
| `SHA256SUMS.uncompressed` | — | integrity of the *contents* |

`stage_S0.csv` columns: `src_id, tile_id, object_id, ra, dec`. Positions are
ICRS degrees on the DSS/GSSS plate solution with the per-plate CRPIX correction
applied (`docs/DSS_WCS_TWO_SOLUTIONS.md`), **and then a per-tile astrometric
refit against Gaia** — see "The coordinates carry a Gaia refit" below before
using them for anything positional. That refit is not what
[`docs/REPRODUCING.md`](../../docs/REPRODUCING.md) specifies, and reading that
document alone will not reproduce this file.

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

The veto chain is **Gaia + PS1 + USNO-B**. Solano et al. (2022) use two, Gaia and
PS1; USNO-B is this project's addition and is a deliberate deviation, not an
implementation of the paper. It is the third-largest single stage by removals. If
you are comparing against a paper-parity pipeline, that difference is yours to
account for.

## The coordinates carry a Gaia refit

**This is a disclosure of a step that should not have run, and did.**

On top of the plate solution, the pipeline applies a per-tile astrometric refit
("WCSFIX"): it matches that tile's own SExtractor detections to Gaia positions
propagated to the plate epoch, within 5″, then fits a **degree-2 polynomial**
correction in each axis, σ-clipped at 1.5″ over two iterations, requiring ≥20 tie
points (fallback: 15″ bootstrap, degree 1, ≥10 points). A typical tile fits ~6,000
tie points to a residual σ of ~0.11″, dropping ~2% as outliers.

**It is not cosmetic.** Measured over all **311,915 survivor rows on 26,115
tiles** (`tools/measure_wcsfix_shift.py`), the refit displaces positions by:

| p25 | median | p75 | p90 | p95 | p99 | max |
|---:|---:|---:|---:|---:|---:|---:|
| 0.234″ | **0.476″** | 0.929″ | 1.708″ | 2.400″ | 3.406″ | 6.835″ |

**73%** of rows move further than the 0.25″ dedup tolerance and **7.6%** further
than 2″. Bootstrapped over tiles — the right unit, since every row in a tile
shares one refit — the median is 0.476″ (95% CI 0.426–0.525) and the fraction
beyond 2″ is 7.58% (95% CI 5.97–9.36).

Two things about that population. It is the per-tile survivor rows *before*
global deduplication, so sources appearing in several overlapping tiles are
counted once per tile; those sit near tile edges, where a degree-2 field is
largest, which may bias the figure slightly high. And the 14 known
astrometry-defect tiles are excluded (13,304 rows, 4.1%; their own displacement
is median 0.41″, max 21.31″).

**It defaults ON and was left on unintentionally.**
[`docs/REPRODUCING.md`](../../docs/REPRODUCING.md) specifies this pipeline on the
*raw* plate WCS and instructs setting `VASCO_WCSFIX_DISABLE=1`. That document even
records "a missing `VASCO_WCSFIX_DISABLE` ran an entire 642-plate campaign
WCS-fixed" as a cautionary note. **That campaign is this catalogue.** The
instruction and the artifact disagreed, and the artifact was released without the
discrepancy being noticed. It is stated here rather than quietly corrected because
a dataset whose argument is reproducibility cannot ship instructions that do not
reproduce it.

**To reproduce this file**, leave `VASCO_WCSFIX_DISABLE` *unset* and keep the
0.25″ dedup tolerance. Following `REPRODUCING.md` as written gives the raw-WCS
variant with a 3.0″ tolerance — a legitimate catalogue, but a different one.

### The circularity, and its actual reach

Gaia informs the astrometry and Gaia is then the first veto. Say that plainly.

What limits it: the correction is a **single smooth degree-2 field per tile fitted
over thousands of tie points**, so it cannot move an individual source onto an
individual Gaia star — every source in the tile moves coherently. The residual σ
of ~0.11″ is at the measurement noise of the detections, i.e. the field is
describing a real plate-scale distortion rather than absorbing per-source offsets.

**This has now been measured, and the circularity does not bite.** The test is
`tools/wcsfix_veto_bias.py`: rerun the 5″ veto criterion on the same detections
under raw and refit coordinates, against Gaia *and* against two catalogues the fit
never saw. Over 250 tiles and 372,957 detections:

| catalogue | role in the fit | raw | refit | Δ | median sep, raw → refit |
|---|---|---:|---:|---:|---|
| Gaia | the bootstrap catalogue | 93.25% | 93.26% | **+0.01** | 0.50″ → 0.16″ |
| PS1 | never seen by the fit | 96.84% | 96.89% | +0.05 | 0.67″ → 0.37″ |
| USNO-B | never seen by the fit | 96.98% | 97.00% | +0.02 | 0.53″ → 0.27″ |

**Gaia-specific excess: −0.03 points** — Gaia rises *less* than the controls, not
more. Null-shifted chance rates are unchanged (10.19 → 10.22%, 20.96 → 20.99%,
10.43 → 10.43%), as they must be if nothing pathological is happening.

The median separation roughly halves against **all three** frames, including the
two the refit never touched. That is what a real astrometric correction looks
like; snapping to Gaia would improve Gaia alone.

There is a structural reason the result is this clean: **the veto threshold is 5″
and the refit displaces sources by a median 0.48″, with 99% inside 3.4″.** A
correction an order of magnitude below the threshold can only flip membership for
sources sitting almost exactly at 5″. The refit buys astrometric precision; it
cannot buy veto decisions.

Caveat: measured on the 504 tiles that retain full per-stage catalogues, δ 41.5 to
86.6. The argument above is geometric and should not be declination-sensitive, but
the measurement is not survey-wide.

### If you want raw plate coordinates

The per-tile survivor catalogues retained **both** coordinate systems
(`ALPHA_J2000`/`DELTA_J2000` alongside `RA_corr`/`Dec_corr`), so a raw-WCS-reported
variant of this catalogue is derivable without re-running anything.

One caveat makes that variant weaker than it looks: **the vetoes matched on the
refit coordinates**, so which rows are here was decided WCS-fixed. Swapping the
reported positions does not undo that. A catalogue genuinely free of the refit
requires a re-run with `VASCO_WCSFIX_DISABLE=1`, because the per-tile detection
catalogues from this campaign were not retained.

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

## Known defect: 9.1% of rows have unreliable coordinates

**8 tiles carry 12,273 rows (9.1%) whose positions are wrong by roughly 11″.**
They are the rows marked `verdict=defect` in
`known_astrometry_defect_tiles.csv`. **Exclude them from any positional
crossmatch**, or you will measure this defect rather than whatever you set out to
measure.

**Corrected 2026-08-14.** An earlier version of that file listed **14** tiles and
13,270 rows, and this section claimed 9.8%. Six of those tiles were mis-flagged:
they were identified by an anomalous survivor count, which is a symptom of bad
astrometry but also of other things. Measured directly, their positions are
**fine** — median separation to Gaia of 0.19–0.33″, against 0.16″ for normal sky.
Telling users to discard 997 good rows was itself an error, so the file now keeps
all 14 rows with an explicit `verdict` column (`defect` / `cleared`) and the
measured separations behind it, rather than silently shrinking.

### What triggers it: the tile crosses RA 0

The eight defective tiles are exactly the eight whose sky footprint **straddles
the RA 0/360 meridian**. The six cleared ones are exactly the six that do not. The
separation is complete — 8 of 8 against 0 of 6 — and the severity scales with
declination, because a 1° tile spans 1/cos(δ) degrees of right ascension and at
δ 87 that is 19°, making a crossing both likelier and wider.

The per-tile astrometric refit's own fit residual records the same gradient:
1.30″ at δ 41.8, 2.56″ at δ 58.8, 6.58″ at δ 76.6, and 28.6″ at δ 86.5, against
0.19–0.48″ on the tiles that do not cross.

**The refit is not the cause — it is a failed rescue.** On these tiles it improves
the median separation to Gaia from 11.31″ to 11.29″, because it bootstraps within
5″ and cannot find true counterparts at 11″. The damage is already present in the
tile's base WCS.

**The mechanism is identified, and fixed in the code (2026-08-14).** The slicer
refits a clean TAN for each tile with astropy's `fit_wcs_from_points`, and did not
pass `proj_point`. At its default of `"center"` that function derives the
projection's fiducial from `lon.min()`/`lon.max()`, which for a field crossing
RA 0 are both *at* the wrap rather than at the field's edges — so the fiducial
lands near the meridian instead of on the tile.

On an ordinary field the fit simply absorbs that into `CD`/`CRPIX` and nothing is
lost, which is why it survived a whole survey unnoticed. Near the pole, where the
GSSS plate solution is genuinely hard to represent as a plain TAN, the fit
**diverges** instead. Measured across all 49 tiles of XE011:

| tile | refit residual, as released | with `proj_point` passed |
|---|---:|---:|
| RA349.417 δ+86.5 | **143.933″** | **0.082″** |
| RA11.593 δ+87.0 | **29.165″** | **0.130″** |
| the other 47 | 0.04–0.16″ | unchanged |

Meridian crossing alone is not sufficient — nine of XE011's tiles cross it and
only the two nearest the pole diverge. It is the combination.

**The pipeline measured this failure and averaged it away.** The residual was
computed for every tile at slice time, but only a plate-level `median of medians`
was ever reported, and 47 healthy tiles reduce a 143.9″ failure to `0.1041″`. The
slicer now refuses to write any tile whose refit exceeds 1″, and prints the worst
tile alongside the median.

**These rows are repairable, and this catalogue has not yet been rebuilt.** The
fix is in the code but the released `stage_S0.csv` is unchanged and still contains
the 12,273 affected rows — that is why they are still listed as `defect` above.
Repair means re-slicing the five affected plates and re-running those eight tiles,
which will move positions by ~11″ and let the vetoes finally act on them. Expect
the corrected catalogue to be *smaller*: most of these are ordinary stars that
only survived because nothing could match them.

The detections themselves are real. On the affected tiles they are uniformly
distributed, sit on genuine flux peaks, and are indistinguishable in FWHM (2.60
vs 2.61 px) and instrumental magnitude (8.90 vs 8.85) from rows that behave
normally. What fails is the coordinate solution: on a control tile the brightest
Gaia stars land a median **2.24 px** from a flux peak against **7.07 px** for a
deliberately-shifted null, while on an affected tile they land **7.14 px** away
against a **7.21 px** null — the astrometry there is statistically
indistinguishable from random.

That is also why these rows survive the pipeline in such numbers. The vetoes match
Gaia, PS1 and USNO-B within 5″; a row displaced by ~11″ matches nothing, so
nothing removes it. A normal tile yields ~3 rows here, these yield hundreds to
thousands.

**Ruled out** as causes: the partial-cone veto bug (their cone queries were
complete), plate-edge geometry (their distance from plate centre matches the
control distribution), the plate solutions themselves (other tiles on the same
plates are normal — the eight come from five plates of 49 tiles each), and
spurious detections. The rows are expected to be **repairable** rather than
discardable — most are probably ordinary stars the vetoes will remove once the
positions are corrected, which would make a corrected catalogue *smaller* than
this one.

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
