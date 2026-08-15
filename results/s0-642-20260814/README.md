# S0 — POSS-I candidate catalogue, 642 plates (2026-08-14)

**122,820 rows.** Detections on POSS-I red plates that survive the MNRAS 2022
filter chain, the Gaia / PS1 / USNO-B 5″ vetoes, the diffraction-spike mask, and
global deduplication at 0.25″.

**Supersedes [`../s0-642-20260813/`](../s0-642-20260813/)** (135,066 rows), which
should no longer be used. That build's ~11″ astrometry defect is repaired here.
The older folder is kept because its README carries the diagnosis.

Produced entirely from public inputs. Nothing here derives from any unpublished
catalogue — see [`tools/audit_independence.py`](../../tools/audit_independence.py).

## Files

| file | rows | what |
|---|---:|---|
| `stage_S0.csv.gz` | 122,820 | the catalogue |
| `tile_manifest.csv.gz` | 31,458 | every tile processed, with its plate |
| `verification_s0_gaia_invariant.csv.gz` | 25,643 | per-tile Gaia-contamination check |
| `verification_dedup_radius_sweep.json` | — | the measurement behind the 0.25″ tolerance |
| `RUN_SUMMARY.txt` | — | parameters and counts as the build recorded them |
| `repaired_astrometry_tiles.csv` | 14 | the 8 repaired tiles and the 6 that were never defective |
| `primary_plate_flags.csv.gz` | 122,820 | per-row coverage partition (added 2026-08-15, see below) |
| `SHA256SUMS` | — | integrity of the files as shipped |
| `SHA256SUMS.uncompressed` | — | integrity of the *contents* |

`stage_S0.csv` columns: `src_id, tile_id, object_id, ra, dec`. Positions are ICRS
degrees on the DSS/GSSS plate solution with the per-plate CRPIX correction applied
(`docs/DSS_WCS_TWO_SOLUTIONS.md`), **and then a per-tile astrometric refit against
Gaia** — see "The coordinates carry a Gaia refit" below.

### Content hashes — quote these, not the `.gz` ones

| file | rows | bytes | sha256 |
|---|---:|---:|---|
| `stage_S0.csv` | 122,820 | 10,995,353 | `2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0` |
| `tile_manifest.csv` | 31,458 | 3,850,873 | `5dcb90dc5d98550e5a60246aced2b097922a267c69e81f27d45d16a288142a99` |
| `verification_s0_gaia_invariant.csv` | 25,643 | 864,551 | `18f1ba101b1b6752f36ce20fdc94e31ab8c3639523ee62fb901fc0ea0b0427e2` |
| `primary_plate_flags.csv` | 122,820 | 11,400,373 | `ae7a599408a8694439b5150bb1320232ad9bdc6d70ce7373cbcd3f58b6b2debf` |

`SHA256SUMS` covers the `.gz` files and verifies the *transfer* only — gzip output
is not reproducible across implementations, so anyone who repacks or mirrors this
will fail that check while holding identical data. The uncompressed hashes above
are what identify the catalogue and what a citation should pin.

```bash
sha256sum -c SHA256SUMS               # the files as shipped
zcat stage_S0.csv.gz | sha256sum      # must equal the table above
```

## What changed from the 2026-08-13 build

**8 tiles had a diverged WCS. Repaired, and the catalogue got smaller.**

| | rows |
|---|---:|
| 2026-08-13 release | 135,066 |
| rows removed on the 8 repaired tiles | −12,273 |
| rows they yield with correct astrometry | +30 |
| net dedup difference | −3 |
| **this release** | **122,820** |

The 8 tiles held **12,273 rows and now hold 30** — 0.2%. That is the expected
result and the strongest confirmation the diagnosis was right: they were ordinary
stars that survived only because an ~11″ position error meant no veto could match
them. Repairing the astrometry let Gaia, PS1 and USNO-B remove them, as they had
already removed their neighbours everywhere else.

Per-tile detail is in `repaired_astrometry_tiles.csv`. The astrometric refit's
own residual on those tiles fell from **1.3–28.6″ to 0.11–1.45″**, and its usable
tie points rose from 28–6,360 to 2,629–14,074.

**Six tiles previously flagged were never defective.** The 2026-08-13 release
listed 14 tiles and told users to exclude all of them. Six — five on XE074 plus
one on XE296 — have median separations to Gaia of 0.19–0.33″ against 0.16″ for
normal sky. They were flagged by an anomalous survivor count, which is a symptom
of bad astrometry but also of other things, and that inference was never checked
against the astrometry. Those 997 rows are ordinary and are unchanged here.

### The root cause

The slicer refits a clean TAN for each tile with astropy's `fit_wcs_from_points`
and did not pass `proj_point`. At its default of `"center"` that function derives
the projection's fiducial from `lon.min()`/`lon.max()`, which for a field crossing
RA 0 are both *at the wrap* rather than at the field's edges — so the fiducial
lands near the meridian instead of on the tile.

On an ordinary field the fit absorbs that into `CD`/`CRPIX` and nothing is lost,
which is why it survived a whole survey. Near the pole, where the GSSS plate
solution is genuinely hard to represent as a plain TAN, the fit **diverges**
instead. Across all 49 tiles of XE011:

| tile | as released 2026-08-13 | with `proj_point` passed |
|---|---:|---:|
| RA349.417 δ+86.5 | **143.933″** | **0.082″** |
| RA11.593 δ+87.0 | **29.165″** | **0.130″** |
| the other 47 | 0.04–0.16″ | unchanged |
| plate median of medians | 0.1041″ | 0.1045″ |

Meridian crossing alone is not sufficient — nine of XE011's tiles cross it and
only the two nearest the pole diverge. It is the combination, which is why the
affected tiles are exactly the high-declination crossers.

**The pipeline measured this failure and averaged it away.** The residual was
computed for every tile at slice time, but only a plate-level *median of medians*
was reported, and 47 healthy tiles reduce a 143.9″ failure to `0.1041″` — a
number that looks perfect and is in the slice log. The slicer now refuses to write
any tile whose refit exceeds 1″ (healthy tiles sit at 0.04–0.17″) and prints the
worst tile beside the median.

One honest note on how this was found: an earlier candidate was tested and
rejected first. A *synthetic* meridian-crossing field with a displaced fiducial
fits to 0.002–0.022″, so the displaced fiducial looked harmless. The synthetic was
too easy — a TAN fitted to a pure TAN. Only the real GSSS-distorted near-pole
cutout reproduces the divergence. Diagnosing this class of failure needs the
actual plate.

## This catalogue is also a repair of a veto bug

An earlier build of the same run gave **310,700 rows and was wrong**, for an
unrelated reason. The mirror cone query selected HEALPix pixels with
`cone_search_lonlat`, which returns only pixels whose *centres* fall inside the
radius. At nside=32 a pixel spans ~1.8° against a ~0.76° veto cone, so pixels
overlapping the cone with their centre outside were silently dropped: the veto
catalogue was partial and real stars survived. Nothing errored.

It acts only where the veto cone crosses the HEALPix polar-cap boundary at
δ = arcsin(2/3) = 41.81°. **2.18× inflation; 56% of that build was un-vetoed Gaia
stars.**

Four things you cannot check from the code alone:

1. **It was found by chasing an unexplained excess, not by a test.** The standing
   check that would have caught it on day one
   (`tools/check_s0_gaia_invariant.py`) was written afterwards. It reports 10.02%
   on the old build and **0.02%** on this one — the ledger is included above.
2. **The cheap repair failed.** Re-vetoing in place does not work: the morphology
   filter derives its FWHM window from the population being filtered, so the extra
   stars narrowed that window and the bug *removed* real rows too. 504 tiles were
   fully re-run.
3. **Two obvious scoping shortcuts are wrong.** Declination alone is insufficient
   (1,503 tiles above the boundary were provably unaffected), and "re-vetoing this
   tile changed nothing" is insufficient too — it would have left 90 corrupt tiles.
4. **An independent catalogue corroborates it.** USNO-B, which played no part in
   the repair, finds POSS-II second-epoch counterparts for **98.94%** of the
   removed rows (7.9× above its own null) versus **0.09%** of rows retained (95×
   below), agreeing with the pipeline's own 98.76% removal to 0.18 points. That
   signal was present and recorded as unexplained on earlier runs, and not
   followed up.

**Residual uncertainty**: 90 tiles holding ~504 rows sit in a class the scoping
check flags but the re-veto could not confirm.

## The coordinates carry a Gaia refit

On top of the plate solution, the pipeline applies the per-tile refit described
above: detections matched to proper-motion-propagated Gaia within 5″, a degree-2
polynomial per axis, σ-clipped at 1.5″ over two iterations, ≥20 tie points
(fallback 15″/degree-1/≥10). A typical tile fits ~6,000 tie points to ~0.11″.

**It defaults ON and was left on unintentionally.**
[`docs/REPRODUCING.md`](../../docs/REPRODUCING.md) specifies this pipeline on the
*raw* plate WCS and instructs `VASCO_WCSFIX_DISABLE=1`. **To reproduce this file**,
leave that variable *unset* and keep the 0.25″ dedup tolerance; following
`REPRODUCING.md` as written gives the raw-WCS variant with a 3.0″ tolerance — a
legitimate catalogue, but a different one.

Measured over 311,915 survivor rows on 26,115 tiles of the previous build, the
refit displaces positions by a median **0.476″** (p90 1.708″, p99 3.406″, max
6.835″); 73% move further than the dedup tolerance and 7.6% further than 2″.

### The circularity, and its measured reach

Gaia informs the astrometry and Gaia is then the first veto. Stated plainly.

What limits it: the correction is a single smooth degree-2 field per tile fitted
over thousands of tie points, so it cannot move an individual source onto an
individual Gaia star — every source in the tile moves coherently.

**Measured, and it does not bite.** `tools/wcsfix_veto_bias.py` reruns the 5″ veto
criterion on the same detections under raw and refit coordinates, against Gaia
*and* two catalogues the fit never saw. Over 250 tiles, 372,957 detections:

| catalogue | role | raw | refit | Δ | median sep raw → refit |
|---|---|---:|---:|---:|---|
| Gaia | bootstrap | 93.25% | 93.26% | **+0.01** | 0.50″ → 0.16″ |
| PS1 | unseen by the fit | 96.84% | 96.89% | +0.05 | 0.67″ → 0.37″ |
| USNO-B | unseen by the fit | 96.98% | 97.00% | +0.02 | 0.53″ → 0.27″ |

**Gaia-specific excess: −0.03 points** — Gaia rises *less* than the controls.
Median separation halves against all three frames, including the two the refit
never touched, which is what a real astrometric correction looks like.

The structural reason: **the veto threshold is 5″ and the refit displaces sources
by a median 0.48″.** A correction an order of magnitude below the threshold can
only flip membership for sources sitting almost exactly at 5″.

## Two deviations from Solano et al. (2022)

**The veto chain is three catalogues, not two.** USNO-B is this project's
addition. Measured on the 504 tiles retaining a full per-stage chain, a
paper-parity catalogue (Gaia + PS1 only) is **+16.2%** larger. That cannot be
obtained by subtracting USNO-B's removals — the vetoes run before the filters and
the σ-clip window is population-derived, so dropping the veto also *removes* rows
the released build kept. See `tools/paper_parity_filter_arm.py`.

**Three implemented stages were not run** — SkyBoT, SuperCOSMOS, VSX. In the
predecessor pipeline they removed ~0.6% of survivors combined, and SkyBoT and VSX
removed nothing.

## The coverage partition — `primary_plate_flags.csv.gz` (added 2026-08-15)

Full-plate slicing searches sky on **every** plate that covers it. A
cutout-based pipeline — the design of Solano et al. (2022) — queries a DSS
service per position and searches each position once, on the plate the service
selects. That is a third deviation, and unlike the other two it concerns
*where* we looked, not how we filtered. This sidecar marks it row by row; the
catalogue itself is unchanged and its hashes above remain valid.

Per row: `det_plate` (where the detection was made), `primary_plate` (the plate
a per-position query would serve — nearest plate centre, centres transcribed
from the GSSS headers into the public manifest; **no fitted or tuned
parameter**), `is_primary`, `primary_has_det` (the primary plate's own raw
detections contain a source within 5″ of this row), and `sep_margin` (how far
the row sits from the plate-selection boundary).

The nearest-centre rule was validated against 11,727 archive tiles whose
headers record the plate the STScI cutout service actually served:
**99.04% agreement**, with the ~1% disagreements being near-equidistant
boundary ties.

| partition | rows |
|---|---:|
| whole catalogue | 122,820 |
| `is_primary` — sky a per-position design searches on the same plate | 68,071 (55.4%) |
| non-primary with a primary-plate counterpart | 122 |
| **single-plate content in multiply-searched sky** | **54,627 (44.5%)** |

Those 54,627 rows exist on one plate's pixels only — the primary plate shows
no raw detection within 5″ (0.22%, against a 2.68% shifted-null; landing below
the null is expected, since catalogue rows are veto survivors and their sky is
star-depleted). A full-plate search finds such content with certainty; a
cutout design finds it only when its tile grid happens to serve that plate —
for boundary sky that is close to a coin flip (plate-at-tile-centre vs
plate-at-source-position differ for a measured ~15% of positions), and for
deep-rim sky effectively never.

**This is a partition, not a quality cut.** Filtering to `is_primary` discards
real content: measured against the public vanish-possi catalogue, it loses
**9.1% of R matches** (98 of 1,072). Quote the partition counts side by side;
anyone filtering should do it with that cost in view.

Two honest limits: 743 rows have a primary plate outside the 634-plate
detection library, so their `primary_has_det=False` is by absence of data, not
measurement; and a genuine single-epoch transient on a plate rim would sit in
the single-plate partition too — the flag says where a row was findable, not
what it is.

**An independent catalogue now says more about that partition.** SuperCOSMOS —
a separate digitization of the same POSS-I E plates, and one that played no part
in this rule, which is pure geometry — fails to confirm **60.5%** of the
single-plate rows against **23.6%** of the `is_primary` rows. PTF, testing
present-day persistence instead of scan quality, separates the two not at all
(3.46% vs 3.24%). So the single-plate partition is disproportionately **scan
artifacts**, not merely sky that one pipeline searched and another did not.
Method, confound tests and caveats:
[`docs/POSTPROCESS_STAGES.md`](../../docs/POSTPROCESS_STAGES.md). This changes
nothing in the catalogue or in this sidecar; it is evidence about how to read
them.

Regenerate with `tools/build_primary_plate_flags.py` followed by
`tools/check_primary_counterparts.py`; both run from this folder's own files
plus the public plate manifest and the (regenerable) per-plate detection CSVs.

## What this catalogue is not

A row here is a source that is **on the plate and not in the modern catalogues we
checked**. That is a statement about two catalogues, not about the universe. It is
not a vanishing-object list, not a transient list, and no row should be called a
detection of anything until it has been inspected individually.

## Verify it

```bash
sha256sum -c SHA256SUMS
VASCO_GAIA_CACHE=<gaia_mirror> python3 ../../tools/check_s0_gaia_invariant.py \
    --s0-csv <(zcat stage_S0.csv.gz) --out /tmp/s0_gaia.csv
```

The invariant must PASS: the fraction of rows with a Gaia source within 5″ must be
~0, because the veto removed exactly those. It reports **0.02%** here. If it fails,
the veto did not cover that sky — whatever the logs say.

## Licence and attribution

Catalogue: see [`LICENSE-DATA`](../../LICENSE-DATA). Plate scans are DSS/POSS-I via
IRSA; the acknowledgements those require are in
[`ACKNOWLEDGMENTS.md`](../../ACKNOWLEDGMENTS.md) and must travel with any reuse.
