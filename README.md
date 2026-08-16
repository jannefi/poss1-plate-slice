# poss1-plate-slice

An independent, reproducible search for vanishing sources in POSS-I, built by
slicing tiles locally from full-plate DSS scans instead of requesting cutouts
from an archive service.

Everything runs from **public inputs only**: the published POSS-I plate list,
IRSA's bulk full-plate scans, and public reference catalogues.
`tools/audit_independence.py` exists so you can verify that rather than take my
word for it — it scans the working tree *and* the entire git history.

```bash
python3 tools/audit_independence.py
```

## Background — why this exists

This project is derived from **VASCO60**, my earlier pipeline for reproducing the
POSS-I vanishing-source search of Solano et al. (2022). VASCO60 was built to follow
the published methodology closely. It works, but it ran into a ceiling that no
amount of parameter tuning could lift, and the ceiling turned out to be
architectural rather than scientific.

Two goals drove the split into this repository:

1. **Reach parity with the reference catalogue.** VASCO60 could not, and
   diagnosis showed why — see below.
2. **Be reproducible by someone else.** A result that depends on data only I hold
   is not a result anyone can check.

The published headline numbers here are measured against the **public** POSS-I
vanishing-source catalogue, so a third party can reproduce both the method and
its validation end to end.

### The archive-cutout ceiling

Cutout services decide *for you* which plate covers a requested position, and
that decision is **position-dependent**: request the same patch of sky at two
slightly different coordinates and you can be served two different plates.

For a survey-scale search this becomes a systematic you cannot correct after the
fact, because plate identity was never yours to choose. Denser tiling does not
fix it — the probability that a service picks a different plate plateaus at
arcminute scales, so a hundredfold increase in tile count buys only a few points.
It is a property of the interface, not of the sampling.

Slicing from full-plate scans inverts the relationship: **you name the plate, and
the tiling is yours.** This repository does both and unions the two arms, because
they miss different things.

Measured head-to-head against an earlier catalogue of ours built the archive way,
the slice route loses essentially nothing at the detection stage — **99.57%
agreement at 5″**, with at most 25 real sources in 13,642 not reproduced. Where
the two catalogues do diverge, at the candidate level, the cause is the MNRAS
2022 filter chain rather than the pixels. See
[`docs/ARCHIVE_VS_SLICE_PARITY.md`](docs/ARCHIVE_VS_SLICE_PARITY.md).

### The astrometric finding

DSS headers carry **two independent astrometric solutions**, and common tools
disagree about which to use. `astropy`/`wcslib` evaluates the GSSS plate
polynomial; SExtractor reads the explicit FITS `CRPIX`/`CD` keywords. On **209 of
642** POSS-I plates these disagree by **~2.3 arcsec**.

The explicit solution is correct. IRSA's full-plate scans carry *only* the GSSS
solution, so naive slicing inherits the wrong answer on a third of the sky — a
systematic that sits *inside* a 5-arcsec match radius and is therefore invisible
to the obvious test.

`tools/build_plate_crpix_table.py` measures and corrects it from public headers
alone: no catalogue, no fitting, no tuned threshold. Full write-up in
**[`docs/DSS_WCS_TWO_SOLUTIONS.md`](docs/DSS_WCS_TWO_SOLUTIONS.md)** — useful to
anyone slicing DSS plate scans, whatever they are looking for.

## What this pipeline claims

This pipeline reproduces what can be called the **Palomar transient effect**:
point sources that are present on POSS-I plates and have no counterpart in modern
deep surveys, after every veto here has been applied.

**It makes no claim about what those objects are, and offers no hypothesis.**

Two explanations remain fully consistent with the data, and nothing in this
pipeline can distinguish between them:

- genuine astrophysical sources that were there and are not now, and
- defects in the photographic emulsion, its handling, or the scanning process.

Both may even be true of different rows in the same output. A detection that
survives the full chain is a source that is *on the plate* and *not in modern
catalogues*. That is the entire claim, and it is a statement about two
catalogues, not about the universe.

Settling the question almost certainly requires physically examining the original
glass plates — under a microscope, in the archive at Caltech. No amount of
reprocessing the digital scans can substitute for that, because every digital
copy inherits whatever the emulsion and the scanner did. This repository is built
to make the *selection* reproducible and auditable, so that whoever does go to the
plates has a defensible list to take with them.

## Deviations from VASCO60

Deliberate simplifications, listed so they stay auditable rather than silent.

| deviation | rationale |
|---|---|
| **Simplified tessellation** — plain 7×7 grid of 60′×60′ tiles per plate, no corner-avoidance gate | The gate excluded plate corners on astrometric-quality grounds. Dropping it keeps full plate coverage and makes the grid trivially describable, which matters more here than trimming the worst corners. |
| **No 30′ circular cut after download** | VASCO60 treated each 60′ tile as a 30′-radius circle. That discards ~21% of every tile and, combined with the grid, leaves gaps. Full square tiles are kept instead. |
| **Grid laid out in pixel space, not RA/Dec** | Stepping in RA/Dec breaks near the pole — at dec ≈ +85 the RA step exceeds 20° per column and tiles walk off the array. Pixel-space layout is uniform at every declination. |
| **Per-plate CRPIX correction** | New; has no VASCO60 equivalent. See above. |

## Getting the data

**Budget ~700 GB of input data.** Measured over the 642-plate footprint: 252 GB
of plate scans, ~6 GB of cutouts, and 427 GB of catalogue mirrors. Working space
is modest next to that — ~14 GB of lean output and 1–8 GB of scratch, because the
runner holds one plate on disk at a time. A **detection-only** reproduction that
skips the veto chain needs the first two items alone, **~258 GB**. Full breakdown
and time budgets in [`docs/REPRODUCING.md`](docs/REPRODUCING.md).

### DSS1 full-plate scans (required)

The recommended source is IRSA's bulk, **plate-addressed** archive:

```
https://irsa.ipac.caltech.edu/data/DSS/images/dss1red/
```

This is what makes the approach reproducible: files are addressed *by plate*, so
you fetch exactly the plate you intend to use. No cutout service offers that. One
file per plate, `dss1red_XE*.fits`.

### Archive cutouts (required for the union arm and the correction table)

From the STScI/MAST cutout service. These are also the **only** place the
explicit FITS WCS solution appears, so the correction table needs at least one
cutout per plate even if you do not use the archive arm.

### Local catalogue mirrors (strongly recommended)

The veto stages cross-match against Gaia (**58 GB**), Pan-STARRS1 (**308 GB**) and
USNO-B1.0 (**61 GB**). Builders are in `scripts/local_cache/`. Without local
mirrors the stages fall back to live VizieR/MAST queries, which is far slower and
will get you rate-limited well before a survey-scale run completes.

> **Practical warning.** Building these mirrors is a bulk-download exercise, and
> for at least one or two of the surveys it is impractical outside a
> well-connected, US-based host — I built mine on an AWS EC2 instance in a US
> region. Budget for that. Check each survey's current bulk-access terms and
> egress costs before starting; the specifics change.

## Prerequisites

Python 3.11+, plus **SExtractor**, **PSFEx** and **STILTS** on `PATH`.
See `requirements.txt`.

## Quick start

```bash
cp config.yaml config.local.yaml          # then edit the paths
python3 -c "from vasco.paths import dump; dump()"     # check they resolve
```

Build the per-plate astrometric correction, then slice the survey:

```bash
python3 tools/build_plate_crpix_table.py \
    --plate-dir <plate_dir> --archive-tiles <tiles_dir> \
    --out data/plate_crpix_table.csv

tools/run_slice_survey.sh \
    --out-dir work/slice \
    --plate-manifest data/plate_manifest.csv \
    --crpix-table data/plate_crpix_table.csv \
    --workers 12
```

Use the wrapper, not `run_fullscale_slice.py` directly: it sets
`VASCO_WCSFIX_DISABLE=1` and clears `VASCO_CIRCLE_ARCMIN`, both of which change
the science silently if left to a default. Check the runner's `[CONFIG]` line
reads `circle_cut=off  wcsfix=off`.

The slice runner is **resumable** — a plate whose output already exists is
skipped, so an interrupted run continues where it stopped.

Recall against a reference catalogue:

```bash
python3 tools/union_parity_fullscale.py \
    --ref-csv <reference.csv> \
    --arm fullplate=work/slice/radec --arm archive=<archive_radec> \
    --combine fullplate+archive --out-dir work/union
```

## Repository layout

| path | |
|---|---|
| `tools/slice_plate_tiles.py` | cut tiles from one full-plate scan |
| `tools/build_plate_crpix_table.py` | per-plate astrometric correction |
| `tools/run_slice_survey.sh` | **launch the survey** — sets the environment explicitly |
| `tools/run_fullscale_slice.py` | whole-survey runner, resumable |
| `tools/check_s0_gaia_invariant.py` | verify the vetoes covered the sky |
| `tools/test_cone_query_coverage.py` | regression test for the mirror cone query |
| `tools/union_parity_fullscale.py` | recall against a reference catalogue |
| `tools/archive_slice_parity.py` | archive-cutout vs plate-slice agreement, both arms |
| `tools/nondetection_cutouts.py` | side-by-side stamps for unmatched rows |
| `tools/classify_displaced_misses.py` | same object, neighbour, or nothing there? |
| `tools/rim_neighbour_counterparts.py` | second partition rule, drawn on plate radius |
| `tools/rim_depth_profile.py` | signal and SuperCOSMOS agreement vs plate radius |
| `tools/audit_independence.py` | proves no private data, tree and history |
| `scripts/stage_*_post.py` | 16 post-processing / veto stages |
| `scripts/local_cache/` | build local Gaia / PS1 / USNO-B mirrors |
| `docs/PARAMETERS.md` | every threshold, its origin, and whether you may change it |
| `docs/REPRODUCING.md` | end-to-end reproduction guide |
| `docs/` | method, astrometry, stage documentation |
| `results/` | derived catalogues, gzipped |

All 16 post-processing stages ship, not only those used for the headline result,
so you can chain them differently and see what changes. `RESULTS.md` records
which chain produced the published numbers.

**Three of them are not applied** to the current results — SkyBoT, SuperCOSMOS
and VSX. Two have since been measured against this catalogue at full scale, and
one of them is large:

| stage | applied? | measured effect on this catalogue |
|---|---|---|
| SuperCOSMOS | no | **would remove 40.0%** (49,139 rows) |
| PTF | no | would remove 2.75% |
| SkyBoT | no | not run; ~23 h for a full pass, yield bounded under 0.76% |
| VSX | no | not measured here |

The SuperCOSMOS figure deserves the emphasis: an independent digitization of the
same POSS-I E plates does not confirm two fifths of this catalogue, and that
result survives testing against coverage, declination, plate-edge vignetting and
magnitude. It also corroborates the release's coverage partition from the
outside, failing 60.5% of the single-plate rows against 23.6% of the rest —
while PTF, testing present-day persistence, separates the two not at all. That
split marks where an independent digitization stops agreeing; it is **not** a
quality score for either partition, and a direct check on the pixels shows most
rim rows carry real flux that SuperCOSMOS nonetheless does not confirm. Method,
controls and caveats:
[`docs/POSTPROCESS_STAGES.md`](docs/POSTPROCESS_STAGES.md).

An earlier version of this section put the cost of these omissions at **~0.6% of
survivors**, citing the predecessor pipeline. That figure was measured on a pool
of 11,027 rows already filtered by morphology and shape stages and **does not
transfer** to the raw post-veto catalogue released here; it is retained, with its
scope stated, in [`docs/PARAMETERS.md`](docs/PARAMETERS.md).

## References

**The search this reproduces**

- **Solano, E., et al. (2022)**, *Discovering vanishing objects in POSS I red
  images using the Virtual Observatory*, MNRAS **515**(1), 1380–1391.
  [doi:10.1093/mnras/stac1552](https://doi.org/10.1093/mnras/stac1552)

  This is the paper the pipeline reproduces. The candidate-selection filters
  implemented in `vasco/mnras/filters_mnras.py` — and referred to throughout this
  repository as "the MNRAS filters" — are those of its **Section 2, Candidate
  selection**.

- Villarroel, B., et al. (2025), *Aligned, Multiple-transient Events in the First
  Palomar Sky Survey*, PASP **137**(10), 104504.
  [doi:10.1088/1538-3873/ae0afe](https://doi.org/10.1088/1538-3873/ae0afe)
- Villarroel, B., et al., *Is there a background population of high-albedo
  objects in geosynchronous orbits around Earth?*, arXiv:2204.06091.
- Bruehl, S., & Villarroel, B. (2025), *Transients in the Palomar Observatory Sky
  Survey (POSS-I) may be associated with nuclear testing and reports of
  unidentified anomalous phenomena*, *Scientific Reports* **15**, 34125.
  [doi:10.1038/s41598-025-21620-3](https://doi.org/10.1038/s41598-025-21620-3)

**The debate over whether these sources are real**

Listed because a reproduction is only useful if the reader can see what is
contested. This repository takes no side; see *What this pipeline claims*.

- Hambly, N. C., & Blair, C. (2024), *On the nature of apparent transient sources
  on the National Geographic Society–Palomar Observatory Sky Survey glass copy
  plates*, RAS Techniques and Instruments **3**(1), 73 — attributes the profiles
  to emulsion flaws.
- Villarroel, B., Solano, E., & Marcy, G. W. (2025), *On the Image Profiles of
  Transients in the Palomar Sky Survey*, arXiv:2507.15896 — the response.
- Watters, W. A., Dominé, L., Little, S., Pratt, C., Knuth, K. H., & Szenher, M.
  (2026), *Critical Evaluation of Studies Alleging Evidence for Technosignatures
  in the POSS1-E Photographic Plates*, PASA (accepted), arXiv:2601.21946.
- Villarroel, B., Streblyanska, A., Bruehl, S., & Geier, S. (2026), *A Response to
  ... Watters et al. (2026)*, arXiv:2602.15171.

**Other independent reproductions**

- **Hayes, Z. (2026)**, *Independent Recovery of Vanishing Sources on POSS-I
  Photographic Plates Using Automated Source Detection and Cross-Epoch Matching*,
  arXiv:2604.04810 — the closest work to this repository: an independent
  automated pipeline, reporting a 63.9% catalogue-level match rate.
- Doherty, B. (2026), *Independent Replication of Nuclear Test-Transient
  Correlations and Earth Shadow Deficit in POSS-I Photographic Plates*,
  arXiv:2604.00056.
- Doherty, B. (2026), *Statistically Significant Linear Alignments Among
  High-Confidence Transient Candidates on POSS-I Photographic Plates*,
  arXiv:2605.01190.
- Busko, I. (2026), *Searching for Fast Astronomical Transients in Archival
  Photographic Plates*, arXiv:2603.20407.
- Busko, I. (2026), *Fast Astronomical Transients in Archival Photographic
  Plates: Using optical aberrations as a tool for discerning real images, from
  plate artifacts*, arXiv:2606.08319.

  > **Scope.** This paper is sometimes invoked as lending support to the POSS-I
  > transient case. It does not bear on POSS-I one way or the other, and the
  > reason is in the paper itself. This analysis uses plates from the
  > Hamburger Sternwarte 0.6-m Doppel-Reflektor via the APPLAUSE archive — **not
  > POSS-I**, and not the scans used here. The paper describes its findings as
  > "preliminary only" and notes that APPLAUSE, drawing on many different
  > telescopes, "does not lend itself easily to the construction of large and
  > homogeneous data sets".
  >
  > The technique is also optics-specific: it discriminates real images from
  > artifacts using a telescope's characteristic aberration signature, and a
  > parabolic reflector and the Schmidt camera that took POSS-I do not share one.
  > Neither its conclusions nor its method should be read as applying to POSS-I or
  > to this work.

**Astrometry**

- Lasker, B. M., et al. (1990), AJ **99**, 2019 — GSSS plate model and scanning.
- Russell, J. L., et al. (1990), AJ **99**, 2059 — GSSS astrometric model.
- Morrison, J. E., et al. (2001), AJ **121**, 1752 — GSC 1.2; DSS headers remain
  GSC 1.1-based.

## Acknowledgements

Special thanks to [Beatriz Villarroel](https://orcid.org/0000-0002-4101-237X) and
[Alina Streblyanska](https://orcid.org/0000-0001-8876-9102) for guidance and
support.

Special thanks to Ivo Busko for his
[plateanalysis](https://github.com/cuernodegazpacho/plateanalysis) software and
the related [arXiv:2603.20407](https://arxiv.org/abs/2603.20407) publication, and
for his help with that approach.

Many thanks to [Mick West](https://www.metabunk.org/members/mick-west.1/) for
finding, fixing and reporting several bugs and improving the pipeline. The local
Gaia / PS1 / USNO-B mirror builders in `scripts/local_cache/` come from his
[vasco60 fork](https://github.com/MickWest/vasco60), and are a large part of why
this pipeline can run without hammering live catalogue services.

Survey data acknowledgements — which are **conditions of use, not courtesies** —
are in [`ACKNOWLEDGMENTS.md`](ACKNOWLEDGMENTS.md). Cite them if you use these
products.

## Licence

Code: **MIT** (`LICENSE`). Derived catalogues: **CC BY 4.0** (`LICENSE-DATA`).
Neither covers the underlying survey imagery.

## Independence

This repository is not affiliated with, endorsed by, or produced in collaboration
with the VASCO collaboration or any other group, and contains no unpublished data
from any of them. It is an independent reproduction, validated against published
data.
