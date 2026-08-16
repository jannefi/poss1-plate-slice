# Results

All figures here are measured against the **published** POSS-I vanishing-source
catalogue (SVO `vanish-possi`, 5,399 rows). Nothing on this page depends on data
that is not publicly available, so every number can be reproduced independently —
see [`docs/REPRODUCING.md`](docs/REPRODUCING.md).

## Footprint

Every POSS-I red plate whose centre declination is ≥ −3.0°: **642** of the 932
scans in IRSA's `dss1red` library. Three independent public checks support that
boundary, and none of them requires a catalogue this project cannot redistribute.

**The threshold is not a tuned parameter.** No plate centre lies between −5° and
−1°, so every cut inside that 4.4° gap selects the same 642 plates.
`tools/build_plate_manifest.py` asserts the gap is empty rather than assuming it,
so a future library that broke the separation would fail loudly instead of
quietly returning a different footprint.

**The southern edge is where the published catalogue ends.** The Solano et al.
(2022) vanishing-object catalogue (SVO `vanish-possi`, 5,399 rows) spans
−3.32° ≤ δ ≤ +87.97°, with 76 rows below δ = 0 and exactly one below −3.

**The plates are populated across essentially the whole set.** The publicly
released NeoWISE-proximate subset (SVO `vanish-neowise`, 171,753 rows — set *W*
of Watters et al. 2026, characterised there as a spatially uniform-random sample
of the parent catalogue) has sources on **633 of the 642**.

For reference, VASCO's published analyses used 635 plates (Villarroel et al.
2026, Fig. 1). This footprint is 642. Which plates the two lists share is not
something this repository can establish, because the 635-plate list has not been
published.

## Raw detection recall

Fraction of published catalogue rows with a detection within a given radius,
over the full 642-plate footprint. 7×7 tiles per plate, per-plate CRPIX
correction applied.

| arm | 1″ | 2″ | 3″ | 5″ | 10″ | 30″ |
|---|---:|---:|---:|---:|---:|---:|
| `fullplate` — local slices of IRSA full-plate scans | 83.39% | 94.26% | 96.91% | **97.24%** | 97.85% | 99.04% |
| `archive` — STScI cutout service | 67.61% | 81.66% | 84.78% | 85.42% | 86.59% | 90.68% |
| **`fullplate+archive`** | 86.22% | 96.06% | 98.30% | **98.59%** | 99.00% | 99.56% |

Arm sizes: `fullplate` 189,241,189 detections over 642 files; `archive`
179,887,397 over 31,004 tiles.

Extending the footprint from the 634 plates of the first campaign to all 642
moved these figures by at most 0.02 points — `fullplate` at 1″ rose from 83.37%
to 83.39%, the union at 2″ fell from 96.07% to 96.06%, and every other cell is
unchanged. That is the expected result rather than a disappointing one: the
reference catalogue was built over its authors' own footprint, so the plates
added here contribute coverage where it has almost no rows to recover. The
recall numbers and the footprint claim are therefore close to independent —
widening the footprint does not flatter the recall.

**The locally-sliced arm alone exceeds the archive arm by ~12 points at 5″.**
That is the central result: naming the plate yourself outperforms letting a
cutout service choose one for you. The two arms miss different sources, so their
union beats either.

Reproduce with:

```bash
python3 tools/union_parity_fullscale.py \
    --ref-csv <vanish_possi.csv> \
    --arm fullplate=<slice_radec> --arm archive=<archive_radec> \
    --combine fullplate+archive --out-dir work/union_R
```

### Provenance of the archive arm

The `archive` arm's detections were produced by this project's predecessor
pipeline at identical settings, and were **not** regenerated for this release —
doing so requires re-downloading 31,004 cutouts from the STScI service. The
`fullplate` arm was generated end to end by the code in this repository.

## Astrometric correction

| | plates | median offset |
|---|---:|---:|
| needing correction | 209 of 642 (32.6%) | 2.29″ → ~0.09″ |
| already correct | 433 | 0.02″ |

Derived from public headers alone, with no catalogue and no fitted parameter —
see [`docs/DSS_WCS_TWO_SOLUTIONS.md`](docs/DSS_WCS_TWO_SOLUTIONS.md). Its effect
on recall is confined to the tight radii: it is worth roughly 29 points at 1–2″
and nothing at 3″ or beyond, because a ~2.3″ systematic sits *inside* a 5″ match
radius and is invisible to a 5″ pass/fail test.

Two plates (XE285, XE311) exceed the 0.2 px scatter threshold and are corrected
with a warning rather than silently.

## Candidate catalogue

**S0 — 122,820 rows** over the 642-plate footprint (31,458 tiles, 0 skips).
Detections surviving the MNRAS filters, the Gaia / PS1 / USNO-B 5″ vetoes, the
spike mask, and global deduplication at 0.25″.

Two deviations from Solano et al. (2022) that anyone comparing pipelines needs up
front. **The veto chain is three catalogues, not two** — USNO-B is this project's
addition. And **the coordinates carry a per-tile degree-2 astrometric refit
against Gaia** ("WCSFIX"), which defaults on and was left on unintentionally:
[`docs/REPRODUCING.md`](docs/REPRODUCING.md) specifies the raw plate WCS, so
following it as written builds a *different* catalogue whose hash will not match.
The refit moves released rows by a median 0.476″ (p90 1.71″, measured over
311,915 rows). Both are set out in
the release [README](results/s0-642-20260814/README.md), including the
Gaia-in-the-astrometry-then-Gaia-in-the-veto circularity — measured against PS1
and USNO-B controls the refit never saw, and found not to bite: a Gaia-specific
excess of **−0.03 points**, because a sub-arcsecond correction cannot move a
source across a 5″ threshold.

**Released**: [`results/s0-642-20260814/`](results/s0-642-20260814/) — the
catalogue gzipped, the tile manifest, the per-tile Gaia-contamination ledger and
the dedup sweep that justify the numbers below, and two hash manifests. Its
[README](results/s0-642-20260814/README.md) carries the method, the caveats and
the bug disclosure, so the folder stands alone if it is ever separated from this
page.

Cite the **uncompressed** content hash, not the `.gz` one — gzip output is not
reproducible across implementations. `stage_S0.csv` is
`2ff92f2210acb387ef9ef4b88d561595d3883e9aab27065042627272b96590f0`.

The 2026-08-13 build (135,066 rows) is **superseded**: eight tiles carried a
diverged WCS refit, and repairing it removed 12,273 rows that had survived only
because an ~11″ position error meant no veto could match them. Those tiles now
yield 30 rows. `results/s0-642-20260813/` is retained for its diagnosis.

### The coverage partition (added 2026-08-15)

There is a third deviation from the cutout-based design, and it concerns
*where* this pipeline looked rather than how it filtered: full-plate slicing
searches sky on **every** plate that covers it, while a per-position cutout
pipeline (Solano et al. 2022's design) searches each position once, on the
plate a DSS service selects. The release now carries
`primary_plate_flags.csv.gz`, partitioning the catalogue by a
**zero-parameter rule** — nearest plate centre, validated at **99.04%**
against 11,727 tiles whose headers record the plate STScI actually served:

| partition | rows |
|---|---:|
| whole catalogue | 122,820 |
| `is_primary` — sky a per-position design searches on the same plate | 68,071 (55.4%) |
| **single-plate content in multiply-searched sky** | **54,627 (44.5%)** |

The single-plate rows exist on one plate's pixels only — the primary plate's
own raw detections show nothing within 5″ of them (0.22%, against a 2.68%
shifted null). A full-plate search finds such content with certainty; a cutout
design only when its tile grid happens to serve that plate. This is the
measured explanation for why a full-plate catalogue is larger than a
cutout-based one at identical filters.

**The flags partition; they do not judge.** Filtering to `is_primary` costs a
measured **9.1% of matches to the published vanish-possi catalogue** (98 of
1,072), so both counts are quoted side by side and any filtering is the
consumer's decision, made with that cost in view. Rule, validation and limits:
the release [README](results/s0-642-20260814/README.md); tools:
`tools/build_primary_plate_flags.py`, `tools/check_primary_counterparts.py`.

### How much does the third veto cost?

USNO-B is the deviation most likely to matter to anyone comparing against a
paper-parity implementation, so its size is measured rather than asserted
(`tools/paper_parity_filter_arm.py`).

It cannot be read off the released catalogue by subtraction. The vetoes run
*before* the MNRAS filters, and the morphology σ-clip derives its window from the
population being filtered — so removing a veto changes the population, changes the
window, and changes which rows the **filters** cut. The measurement therefore
re-runs the filter stage from the pre-USNO-B population, letting the clip window be
re-derived as the pipeline would have derived it.

On the 504 tiles that retain the full per-stage chain, and on which the same code
path reproduces every tile's actual `filtered.csv` **exactly (504/504)**:

| arm | rows |
|---|---:|
| released — Gaia + PS1 + USNO-B | 2,341 |
| paper parity — Gaia + PS1 | **2,721** |
| difference | **+380 (+16.2%)** |

Dropping USNO-B **adds 390 rows and removes 10**. Those 10 are the point: a veto
removal is not purely additive downstream, which is exactly why a subtraction would
have given the wrong answer.

Note the asymmetry — USNO-B cuts 7.70% at the veto stage (47,043 → 43,420) but the
finished catalogue grows by 16.2%. The rows it removes are disproportionately
*good* ones, which is expected: they are real POSS-II second-epoch stars and so
pass the MNRAS quality gates at well above the average rate.

Caveat: these are the 504 re-run tiles, δ 41.5–86.6, not a survey-wide sample.

Three implemented stages — SkyBoT, SuperCOSMOS and VSX — were **not applied** to
this catalogue. Two have since been measured against it at full scale:
**SuperCOSMOS would remove 40.0%** (49,139 rows) and PTF 2.75%, while SkyBoT is
bounded at under 0.76% yield and was not run. SuperCOSMOS is therefore by far
the largest unapplied stage, and it also corroborates the coverage partition
from the outside — see
[`docs/POSTPROCESS_STAGES.md`](docs/POSTPROCESS_STAGES.md) for the measurements
and controls, and [`docs/PARAMETERS.md`](docs/PARAMETERS.md) for the
parameters.

### The veto had a bug, and this number is the repaired one

The first build of this catalogue was wrong, and the correction is large enough
that it would be misleading to present the result without it.

**The defect.** The local-mirror cone query built its HEALPix pixel list with
`astropy_healpix.cone_search_lonlat`, which returns only pixels whose *centres*
fall inside the radius. At nside=32 a pixel spans ~1.8° against a ~0.76° veto
cone, so pixels overlapping the cone with their centre outside were silently
dropped. The veto catalogue was partial, and real stars in that sky survived
into the candidate list. Nothing raised: the stage ran, wrote its crossmatch
file, and reported success. All three vetoes shared the query, so a Gaia leak
was not caught by PS1 or USNO-B either.

**Where it bites.** Only where the veto cone crosses the HEALPix polar-cap
boundary at δ = arcsin(2/3) = 41.81°. Below it the pixelisation is a regular
grid and nothing was lost. Declination predicts where the bug *can* act, never
where it *did*: 1,503 tiles above that boundary were provably unaffected.

**Scale.** The first build gave 310,700 rows. **2.18× inflation** — 56% of it
was un-vetoed Gaia stars.

**How it was found, and how late.** By chasing an unexplained excess in the
candidate count, not by a test. The standing check that would have caught it on
day one — S0 rows with a Gaia source within 5″ must be ~0, now
`tools/check_s0_gaia_invariant.py` — was written *afterwards*. On the published
build it reports 10.02%; on this one, 0.02%.

**How it was repaired.** A post-hoc re-veto of the affected rows was tried first
and **falsified**: the MNRAS morphology filter derives its FWHM window from the
population being filtered, so the extra stars narrowed that window and the bug
*removed* real rows too. No post-process can put those back. The repair is
therefore a genuine partial re-run of 504 tiles — the tiles measured, one by
one, to have had a catalogue source in a dropped pixel. Two cheaper scoping
heuristics were tested and both are wrong: declination (see above), and "the
re-veto removed nothing from this tile", which would have left 90 corrupt tiles
in place.

**Independent corroboration.** USNO-B — which played no part in the repair —
finds POSS-II second-epoch counterparts for **98.94%** of the removed rows
(7.9× above its own null-shifted chance rate), against **0.09%** of the rows
that remain (95× below chance). It agrees with the pipeline's own 98.76%
removal on those tiles to 0.18 points. This signal was present, and recorded as
unexplained, on earlier runs — an outside instrument flagged the defect long
before it was diagnosed, and the signal was not followed up.

**Residual uncertainty**, stated rather than smoothed: 90 tiles holding ~504
rows sit in a class the scoping check flags but the re-veto could not confirm,
and their contribution could move in either direction.

The fix is [`vasco/local_cache_query.py::_cone_pixels`](vasco/local_cache_query.py)
— search with a margin wide enough that no overlapping pixel can be missed, then
prune with an explicit overlap test. `tools/test_cone_query_coverage.py` passes
against it and **fails** against the old code.

## Which stage removes the published sources?

Recall says how many reference sources the pipeline finds; it does not say what
happens to the rest. [`docs/FUNNEL_ATTRIBUTION.md`](docs/FUNNEL_ATTRIBUTION.md)
attributes each loss to the exact stage that consumed it, by row identity
rather than by re-matching. Headline, on a survey-representative 18-plate
sample: detection inside plate cores is **98.6%**, about a fifth survive to S0
(consistent with a whole-survey crossmatch at 19.86%), and of the losses
**83% are the MNRAS 2022 quality filters** against 17% for the whole veto
chain — two thirds of everything lost is the paper's own `SPREAD_MODEL` and
extraction gates.

The same page carries a result that needs none of this pipeline to check:
**nine published reference rows have a Gaia star at their position at the
plate epoch and none at Gaia's catalogue epoch** — objects that moved rather
than objects that vanished. Sample sizes there are small and stated as such.

## Does the image source change the catalogue?

This pipeline slices plate-addressed IRSA scans; the obvious alternative is
position-addressed archive cutouts.
[`docs/ARCHIVE_VS_SLICE_PARITY.md`](docs/ARCHIVE_VS_SLICE_PARITY.md) measures
the difference against an earlier catalogue of ours built the other way
(15,303 rows, STScI cutouts, 725 plates), one-way and footprint-scoped to
13,642 rows. The footprints are **identical above dec 0** — 100.0% coverage in
every band, no geometric edge cases.

**Raw detection agreement is 99.57% at 5″ and 99.93% at 30″.** Of the 58 rows
unmatched at 5″, cutout inspection accounts for all of them: at most **25 are
real sources not reproduced (≤0.183%)**, while **27 have no source above 3σ at
the catalogued position in either run's pixels (0.198%)**. Roughly half the
residual is reference-side, so image sourcing costs this pipeline essentially
nothing at the detection stage.

**Candidate agreement is lower, 88.61%**, and that gap is not the image source.
Read from per-row `reject_reason`, **85.9%** of the removals are the MNRAS
quality filters against **9.4%** for the veto chain — independently reproducing
the 83% / 17% split measured above against the unrelated SVO catalogue. The
cause is threshold sensitivity: the rows this pipeline drops are the rows that
barely passed the other one, sitting at `SNR_WIN` 31.07 against a gate of 30 at
the 10th percentile. **88.61% therefore measures how reproducible the MNRAS
2022 filter chain is across pixel realisations**, not how much the pixels
differ.

## What these numbers do and do not say

Recall measures whether this pipeline *finds* the published sources. It says
nothing about what those sources are. See *What this pipeline claims* in the
README: a detection surviving the full chain is a source that is on the plate and
not in modern catalogues, which is a statement about two catalogues rather than
about the universe.
