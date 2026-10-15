# A POSS-II red control — is the candidate rate a property of the plates or of us?

This pipeline turns POSS-I red plates into transient candidates at a rate of
roughly 200 per plate. That number only means something if it is a property of
POSS-I. If the same chain produced the same rate on any photographic survey, it
would be a property of the chain instead, and the right reading of the catalogue
would be "this is what our filters do to a plate".

So: run the identical chain over POSS-II red plates covering **the same sky**,
and compare per-plate yields against the POSS-I plate that overlaps each one.

**The short answer**: POSS-II yields about **an eighth** of POSS-I's candidates
on the same sky — per-plate ratio **0.1243 ± 0.0152** over **30 plate pairs**,
pooled 644 / 5,939 = 0.1084. The suppression does not come from any one part of
the sky or any one class of plate: split by declination, by galactic latitude
and by plate grade, no stratum departs from it (Kruskal–Wallis p = 0.45, 0.38,
0.23 respectively). Whatever drives POSS-I's rate, this pipeline does not
produce it from photographic plates in general.

**What this does not settle** is *why*, and that limit is not a footnote — it is
in the test's construction. POSS-II red is a different emulsion at a different
depth and a different epoch, sampled 1.68× finer, and "fewer candidates on a
cleaner, deeper plate" is a perfectly ordinary explanation that this comparison
cannot exclude. See [Limits](#limits) before quoting anything here. The
independent point of comparison that *does* control for era and instrument class
is Busko (2026), arXiv:2603.20407, discussed at the end.

Everything below uses public inputs: IRSA's full-plate POSS-II red scans, this
repository's own released POSS-I catalogues, and the pipeline's unmodified
entry points. Per-plate counts ship as
[`docs/poss2_red_control_per_plate.csv`](poss2_red_control_per_plate.csv).

## Getting the right plates: the naming trap at IRSA

IRSA's `dss2red` directory holds 1,907 files named `dss2red_<REGION>.fits`, and
they are **four different red surveys** under one directory name. Only one is
Palomar:

| prefix | files | survey | telescope | site latitude |
|---|---:|---|---|---:|
| **XP** | **897** | **POSSII-F** | **Oschin Schmidt, Palomar** | **+33.356** |
| XS | 606 | AAO-SES | UK Schmidt | −31.277 |
| ER | 288 | SERC-ER | UK Schmidt | −31.277 |
| GR | 116 | AAO-GR | UK Schmidt | −31.277 |

1,010 of the 1,907 files are southern UK Schmidt plates. Taking the directory at
face value would have put a different telescope at a different hemisphere into
the control arm. Take `XP` only. (The POSS-I red equivalents in `dss1red` are
`XE`, which is how this repository addresses them already.)

An index built from the first 40 kB of each file — headers only, by HTTP range
request — gives the usable set: 897 XP plates, declination −0.3° to +87.5°,
epochs 1986–1999, median exposure 75 min, grades A 248 / B 401 / C1 248. Two
data-quality exclusions came out of it: **10 plates carry an identical
placeholder `DATE-OBS` of 1952-08-22**, four decades before POSS-II began, and
two more carry 6–7 minute exposures. All twelve are excluded.

## How pairs were chosen

A POSS-II plate is paired with the POSS-I plate whose centre lies within **0.6°**
of its own — close enough that the two plates, each about 6.5° square, cover
substantially the same sky. 36 of the 897 XP plates pair that way against the 642
POSS-I plates in this repository's paper-parity run; **30 of those 36 were run**.

Pairing on sky is what makes the comparison work, because it cancels the things
that would otherwise dominate a survey-wide count: field density, galactic
latitude, and our own gate choices all apply identically to both arms of a pair.
The unit of comparison is the pair, not the survey.

The 30 pairs span declination −0.2° to +64.9°, galactic latitude 1.0° to 86.6°,
all three plate grades, POSS-II exposures 47–110 min and POSS-II epochs
1987–1999. POSS-I yields on the partner plates run from 17 to 462 rows, so the
set is not drawn from the easy end of either survey.

## The pixel scale differs by 1.68×, and that changes every pixel-valued gate

This is the part of the test that had to be got right, and it is the reason a
naïve "run the same config on both" would have produced a suppression for purely
instrumental reasons — in the same direction as the hypothesis under test.

| | POSS-I red (XE103) | POSS-II red (XP140) |
|---|---|---|
| array | 14000 × 13999 | 23040 × 23040 |
| pixel size | 25.2845 µm | 15.0295 / 15.0000 µm |
| plate scale | 67.2 ″/mm | 67.2 ″/mm |
| **sampling** | **1.699 ″/px** | **1.010 ″/px** |
| plate | 6.61° | 6.46° |
| emulsion / filter | 103aE / plexi | IIIaF / RG610 |

POSS-II is sampled 1.68× finer, so every parameter expressed in pixels means
something different on the sky:

| parameter | on POSS-I | on POSS-II if left unchanged |
|---|---|---|
| `2 < FWHM_IMAGE < 7` px | 3.4″ – 11.9″ | **2.0″ – 7.1″** |
| `extent_delta < 2.0` px | 3.4″ | 2.0″ |
| `DETECT_MINAREA 5` px | 14.4 □″ | **5.1 □″** |
| `PHOT_APERTURES 5` px | 8.5″ | 5.0″ |
| `BACK_SIZE 64` px | 109″ | 65″ |

Left alone, the morphology gate is 1.68× tighter on the sky and the detection
minimum area 2.8× smaller. Every pixel-valued parameter was therefore **scaled
by 1.682 for the POSS-II arm** — detection, photometry, background and the
morphology gates — so that both arms gate on the same angular quantities.

Two consequences worth stating, because they bound how much of this is a choice:

- **It matters a lot.** On the first pair, the four plausible regimes (identical
  pixel parameters, gates scaled only, detection scaled only, fully scaled) give
  99, 132, 24 and 33 candidates respectively against POSS-I's 315 — a 5.5× spread
  that is wider than the difference between "a factor of three" and "a factor of
  ten". The *direction* is the same in all four; the magnitude is not.
- **The PSF model barely matters.** A second pass scaled PSFEx's own pixel
  parameters as well (`SAMPLE_FWHMRANGE`, `PSF_SIZE`), which changed the yield by
  one source, 33 → 32, on an identical raw detection count of 268,957. All 30
  pairs use the fully scaled configuration.

Things that derive the scale rather than assume it needed no change and got
none: the tile slicer sets tile size in arcmin and converts with the plate's own
scale (POSS-II tiles come out at 3564 px × 1.0076″ = 59.85′ against POSS-I's
2119 px × 1.7002″ = 60.04′, so the tiles differ in sky area by 0.6%); SExtractor
reads the scale from the FITS WCS because `PIXEL_SCALE` is absent from both
configs; PSFEx samples automatically. `SEEING_FWHM 1.2` is wrong for both
surveys — POSS seeing is 3–5″ — but it feeds `CLASS_STAR`, which this chain does
not gate on, and it is the same value on both sides, so it cancels. Flagged, not
changed.

**Proper motion is handled by construction, and it matters here more than
anywhere else.** The candidate criterion is "no Gaia or Pan-STARRS counterpart
within 5″", and a star moves about 2.5× further between 1950 and Gaia's epoch
than between 1990 and Gaia's. Left uncorrected, that alone would make POSS-II
yield fewer candidates for reasons having nothing to do with transients. Both
arms veto against Gaia **propagated to that plate's own epoch**, which is what
makes the ratio interpretable at all.

## The funnel, both surveys, 30 pairs

The POSS-I side is this repository's released catalogues restricted to the 30
paired plates: `results/s0-642-paper-parity-20260828/` at the gated stage, and
the 15′ edge-cut `results/s1-642-paper-parity-edge15-20260908/` after the edge
stage (both on the `paper-parity` branch). The POSS-II side is the same chain
run from the IRSA scans, with the scaling above, deduplicated at 0.25″ so the
counts sit on the same footing as the released POSS-I catalogue.

| stage | POSS-II | POSS-I | ratio |
|---|---:|---:|---:|
| gated, angular-matched, deduplicated | 644 | 5,939 | **0.108** |
| − densitometry spots | 642 | 5,939 | 0.108 |
| − edge stage at 15′ | 412 | 4,083 | 0.101 |
| − SHAPE | **9** | **475** | **0.019** |

Per plate: **198 candidates on POSS-I against 21.5 on POSS-II**; after the edge
stage, 136 against 13.7.

**The ratio barely moves through three of the four stages.** Neither of the
cleaning stages is doing the work:

- **Densitometry spots** are the opaque calibration wedges DPOSS used to
  linearise its plates (Gal et al., `astro-ph/0210298`). They are detected per
  plate rather than masked at fixed positions, because the "same size, same
  place" assumption is only half right — the size holds at 300–350 px square,
  the position wanders several hundred px between plates, and **one plate in six
  has none at all**. Cost: 0.2% of plate area. Effect: 2 of 644 rows.
- **The edge stage** is `scripts/stage_edge_post_v2.py` with the released S1's
  own parameters (`--core-radius-deg 5.0 --policy own --min-edge-arcmin 15`),
  read from that release's ledger. It removes proportionately from both sides
  (35.8% of POSS-II, 31.3% of POSS-I) and moves the ratio from 0.108 to 0.101.
  An early single-plate reading that POSS-II survivors were an edge population
  did not survive the other five plates of the six-pair chain, where the pooled
  median edge distance came out at ~31′ on both surveys.

**SHAPE is the one stage that separates them**, and that is discussed under
[What survives, and what it looks like](#what-survives-and-what-it-looks-like).

## The headline number, clustered by plate

Candidate rows cluster hard by plate, so treating 6,583 rows as independent
observations would understate the uncertainty by a large factor. The effective
sample is **30 plates**. The headline is therefore the mean of the per-plate
ratios with a plate-clustered standard error, the same unit this project uses
elsewhere.

| | six pairs (earlier run) | **thirty pairs** |
|---|---:|---:|
| per-plate mean ratio | 0.115 ± 0.017 | **0.1243 ± 0.0152** |
| pooled | 0.107 | **0.1084** |
| per-plate median | 0.108 | 0.098 |
| per-plate range | 0.070 – 0.188 | **0.000 – 0.353** |

The six-plate estimate held: 0.124 against 0.115, a third of a standard error
apart, with five times the plates. The interesting question was never the
distance from parity (57σ, and never in doubt) but whether six plates had drawn
a lucky patch of sky. They had not.

### Per-plate table

| POSS-II | POSS-I | grade | exp (min) | POSS-II epoch | dec | \|b\| | POSS-II gated | POSS-I S0 | ratio |
|---|---|:--|---:|---|---:|---:|---:|---:|---:|
| XP460 | XE337 | A | 65 | 1995-08 | +30.2 | 4.9 | 0 | 17 | 0.000 |
| XP864 | XE618 | A | 75 | 1997-05 | −0.2 | 60.0 | 11 | 292 | 0.038 |
| XP792 | XE558 | C1 | 87 | 1993-05 | +4.8 | 64.6 | 16 | 294 | 0.054 |
| XP852 | XE608 | B | 55 | 1997-01 | −0.2 | 37.1 | 10 | 181 | 0.055 |
| XP090 | XE060 | A | 50 | 1999-03 | +64.9 | 37.2 | 12 | 214 | 0.056 |
| XP828 | XE588 | B | 65 | 1991-09 | +0.3 | 59.9 | 26 | 462 | 0.056 |
| XP433 | XE313 | A | 60 | 1998-04 | +29.8 | 42.8 | 21 | 372 | 0.056 |
| XP882 | XE633 | B | 85 | 1991-06 | +0.2 | 11.3 | 2 | 34 | 0.059 |
| XP888 | XE638 | A | 90 | 1990-07 | +0.3 | 37.0 | 18 | 258 | 0.070 |
| XP780 | XE548 | B | 47 | 1998-12 | +4.8 | 39.9 | 11 | 155 | 0.071 |
| XP424 | XE305 | C1 | 80 | 1988-12 | +30.0 | 1.0 | 8 | 111 | 0.072 |
| XP858 | XE613 | B | 60 | 1998-04 | −0.2 | 58.0 | 11 | 128 | 0.086 |
| XP451 | XE329 | B | 110 | 1989-05 | +29.9 | 47.0 | 25 | 287 | 0.087 |
| XP870 | XE623 | C1 | 95 | 1993-06 | −0.1 | 40.8 | 28 | 291 | 0.096 |
| XP468 | XE344 | A | 80 | 1991-09 | +30.3 | 24.5 | 30 | 309 | 0.097 |
| XP140 | XE103 | B | 60 | 1994-07 | +60.1 | 26.8 | 31 | 315 | 0.098 |
| XP416 | XE298 | B | 75 | 1991-10 | +30.3 | 26.6 | 19 | 166 | 0.114 |
| XP097 | XE066 | B | 75 | 1990-01 | +64.8 | 50.9 | 37 | 317 | 0.117 |
| XP452 | XE330 | B | 90 | 1989-05 | +29.9 | 42.1 | 22 | 173 | 0.127 |
| XP840 | XE598 | B | 80 | 1990-11 | +0.1 | 15.3 | 25 | 190 | 0.132 |
| XP846 | XE603 | C1 | 47 | 1997-02 | −0.1 | 11.3 | 18 | 130 | 0.138 |
| XP786 | XE553 | B | 105 | 1991-01 | +4.7 | 62.2 | 37 | 261 | 0.142 |
| XP459 | XE336 | A | 90 | 1987-06 | +30.1 | 9.2 | 6 | 38 | 0.158 |
| XP442 | XE321 | B | 75 | 1990-01 | +29.8 | 86.6 | 36 | 214 | 0.168 |
| XP834 | XE593 | B | 80 | 1990-10 | +0.2 | 40.7 | 24 | 128 | 0.188 |
| XP894 | XE643 | B | 65 | 1987-09 | +0.3 | 58.0 | 50 | 227 | 0.220 |
| XP425 | XE306 | A | 70 | 1989-11 | +30.0 | 5.3 | 20 | 82 | 0.244 |
| XP876 | XE628 | C1 | 85 | 1991-06 | +0.0 | 15.3 | 21 | 77 | 0.273 |
| XP124 | XE089 | B | 95 | 1989-12 | +59.9 | 34.0 | 45 | 148 | 0.304 |
| XP117 | XE083 | B | 65 | 1992-10 | +60.2 | 4.8 | 24 | 68 | 0.353 |

### A methodological trap, found and worth repeating

The batch runner skipped a plate that produced **zero** gated candidates before
the per-plate table was built — a natural-looking `if not len(c): continue`.
XP460 is exactly such a plate: 1.4 M raw detections, none surviving the gates.

**Silently dropping a genuine zero biases the mean of per-plate ratios upward.**
With XP460 missing the analysis reported n = 29 and a mean of 0.1286; restoring
it at its true ratio of 0/17 = 0.000 gives n = 30 and **0.1243**. The pooled
figure is what catches this, because it must reproduce the chain's own totals
(644 / 5,939) and does so only when every run plate is present. The chain totals
were never wrong — POSS-II contributed no rows either way — but the per-plate
statistic was, by 3% of its value and in the flattering direction.

## No stratum behaves differently

This is what six pairs could not test, and it is the main thing 30 buys.

| stratum | plates | pooled | per-plate |
|---|---:|---:|---:|
| dec < 0° | 5 | 0.076 | 0.083 ± 0.017 |
| dec 0–20° | 10 | 0.110 | 0.126 ± 0.025 |
| dec 20–40° | 10 | 0.106 | 0.112 ± 0.021 |
| dec > 40° | 5 | 0.140 | 0.186 ± 0.060 |
| \|b\| < 20° | 9 | 0.166 | 0.159 ± 0.038 |
| \|b\| 20–40° | 8 | 0.101 | 0.108 ± 0.029 |
| \|b\| > 40° | 13 | 0.100 | 0.110 ± 0.016 |
| grade A | 8 | 0.075 | 0.090 ± 0.027 |
| grade B | 17 | 0.126 | 0.140 ± 0.021 |
| grade C1 | 5 | 0.101 | 0.127 ± 0.039 |

Kruskal–Wallis on the per-plate ratios: declination **p = 0.45**, galactic
latitude **p = 0.38**, plate grade **p = 0.23**. Spearman against the continuous
variables is null throughout — declination ρ = +0.17 (p = 0.37), |b| ρ = −0.21
(p = 0.28), exposure ρ = +0.24 (p = 0.21). The largest apparent trend, low
galactic latitude running high, is 1.2σ.

**The suppression is uniform.** It is not carried by a declination band, by the
Galactic plane, by plate quality or by exposure time. That matters more than the
central value does: a ratio driven by a single stratum would point at a selection
or processing artifact somewhere in the comparison, and there is no such stratum.

### Seven high-background plates, and what they are not

Seven POSS-II plates yield 0.6–1.5 M raw detections where a typical plate gives
60–200 k, and take 4–14× longer to process. Their POSS-I partners are also the
poorest yielders in the set (17–111 rows against a typical 150–460), so **both
surveys go quiet on the same sky**. Their ratios scatter in both directions
(0.000 to 0.353) and the Spearman against raw detection count is null, so they
neither drive nor dilute the result. Worth knowing for run-time budgeting, not
for interpretation.

## What survives, and what it looks like

`scripts/stage_shape_post.py` — [STAGE_SHAPE](STAGE_SHAPE.md), which implements
Busko (2026), arXiv:2603.20407 — is the one stage in the chain that asks a
different question from the gates: it compares each candidate's radial profile
against **flux-matched neighbour stars in the same cutout**. Its reference is
empirical and local, the stars on that plate, so it carries no assumption about
the telescope's optics and transfers between surveys by construction. Its own
pixel-defined parameters were scaled 1.682× for POSS-II like everything else.

It is also the only stage that moves the ratio: 0.101 → **0.019**. Rejection
runs at 97.8% on POSS-II against 88.4% on POSS-I. Measured on the earlier
six-pair chain, where the contingency table is the one that was tested, the
difference is significant (Fisher exact p = 0.0042), and the per-metric
comparison separates on all three:

| metric | POSS-II | POSS-I | Mann–Whitney |
|---|---:|---:|---:|
| profile_diff (vs neighbour stars) | 0.043 | 0.028 | **5.8 × 10⁻⁸** |
| circularity | 0.405 | 0.719 | **3.3 × 10⁻⁹** |
| elongation | 1.190 | 1.157 | 0.0084 |

**This does not say nine real transients survive on POSS-II.** SHAPE is
EXPERIMENTAL and removes 88% of this repository's own POSS-I catalogue here; the
survivor count is the output of a chain that guts both surveys, not an estimate
of anything physical. What is meaningful is that the ratio *moved*, and in which
direction: through the gates, dedup, the densitometry mask and the edge cut it
sat flat at 0.10–0.11, and the one stage that asks "does this look like a star on
this plate" takes it to a fiftieth. The POSS-II remainder is not merely fewer, it
is systematically less star-like — which is the opposite of what a population of
real transients should look like beside POSS-I's.

Inspected at the cutout level on the six-pair chain, the survivors on **both**
surveys come out indistinguishable from magnitude-matched stars on their own
plates (elongation 1.06 on both survivor sets against 1.07–1.08 for the
magnitude-matched star samples), which is what the stage selects for. The surviving difference between the surveys is one of number.

A cross-epoch check that sounds decisive and is not, recorded so it is not cited
as though it were: of 120 POSS-I survivors on the six-pair chain, **0 have a
POSS-II raw detection within 2″** and 2 within 5″. The useful part is the bound —
~1.7% are persistent sources that the Gaia and Pan-STARRS vetoes missed, and a
third epoch agrees with them. The useless part is everything else, because a
transient and a POSS-I plate defect make the *same* prediction: absence on
POSS-II.

## Limits

**The emulsion-and-depth objection is not answered by this test, and this is the
limit that matters.** POSS-II F is a IIIaF emulsion behind an RG610 filter,
reaching about 1.5 magnitudes deeper than POSS-I's 103aE, at a finer grain, at an
epoch 35–45 years later, sampled 1.68× finer. Several entirely mundane mechanisms
predict fewer unmatched sources on such a plate regardless of transients: fainter
real sources are likelier to have a Pan-STARRS counterpart, finer grain means
fewer grain clumps survive as point sources, and the plates are graded and
processed differently. **This comparison is same-sky, but it is not
same-instrument.** It cannot separate "far fewer transients" from "a cleaner,
deeper plate".

- **The angular matching is a choice, and its magnitude shows.** The four
  plausible ways of handling the 1.68× scale difference span 5.5× in yield on the
  tested pair. The direction is the same in all four regimes, and the
  fully-matched regime is the defensible one, but the *size* of the suppression
  depends on a decision made here rather than on anything in the data.
- **Two SHAPE thresholds are not exposed and were not scaled** — a contour
  validity cut at `area > 7` px and a confidence split at `area < 100` px. Neither
  auto-drops a row, but neither was angularly matched either.
- **SHAPE is EXPERIMENTAL.** Nothing downstream of it should be read as a rate.
  Nine survivors carries a Poisson error near ±3.
- **The POSS-I side carries no densitometry mask.** The wedges are a DPOSS
  calibration feature and none appeared in the POSS-I corner boxes checked, so the
  asymmetry should be small, but it is an asymmetry. POSS-I's own corner furniture
  — printed and handwritten labels, which are mid-density text rather than opaque
  patches — is unhandled on both sides and needs a different detector.
- **Depth and exposure pull in opposite directions** and the net sign was not
  measured. POSS-II reaches deeper, which helps a faint flash clear detection; a
  75-minute exposure dilutes a fixed-duration flash about 0.4 magnitudes more than
  POSS-I's ~50 minutes.
- **One survey pair, one pipeline.** Thirty plates establish that the ratio is
  uniform across our own strata. That says nothing about whether POSS-I's rate is
  a property of POSS-I or of this pipeline's interaction with POSS-I specifically.
- **Pair overlap is not exact.** Plate centres agree within 0.6° and the tiles
  differ in sky area by 0.6%, so a pair covers substantially but not identically
  the same sky.

## The independent comparison that does control for era: Busko (2026)

Busko, *Searching for Fast Astronomical Transients in Archival Photographic
Plates*, [arXiv:2603.20407](https://arxiv.org/abs/2603.20407) — the paper
[STAGE_SHAPE](STAGE_SHAPE.md) implements — ran a transient search on mid-1950s
plates from the Hamburg Observatory Großer Schmidtspiegel, taken in pairs about
30 minutes apart, with a candidate defined as present on one exposure and absent
on the next. He reports **41 plates, 70 candidates, 35 after visual vetting** —
1.7 per plate raw, 0.85 vetted. He describes the instrument as having "very
similar optical and mechanical characteristics as the Samuel Oschin telescope …
similar plate scales, plate sizes, resolutions, and limiting magnitudes", so
per-plate is a fair unit.

| | per plate |
|---|---:|
| POSS-I, gated, matched, deduplicated (30 pairs) | **198** |
| POSS-I, after the edge stage | 136 |
| POSS-II, gated, matched, deduplicated (30 pairs) | 21.5 |
| POSS-II, after the edge stage | 13.7 |
| Busko (2026), Hamburg, raw / vetted | **1.7 / 0.85** |

**This is why the Hamburg comparison carries weight the POSS-II one cannot.**
Busko's plates are *contemporary* with POSS-I, on a telescope he describes as
comparable, at comparable limiting magnitudes — so they control for exactly the
era-and-instrument axis that POSS-II does not, and they still show a rate one to
two orders of magnitude below POSS-I's. POSS-I is the outlier on two independent
axes at once.

**The caveat that must travel with that comparison**: his selection is not ours.
He requires a source to vanish between two exposures 30 minutes apart; this
pipeline requires absence from modern catalogues plus the MNRAS 2022 filter
chain. Those select different things and the rates are not interchangeable. What
survives the caveat is the order of magnitude, not the number.

One further distinction, because the two are easy to conflate: a separate, later
paper of Busko's argues from coma aberration on a different, coma-prone
non-Schmidt telescope. That argument is flagged elsewhere in this repository as
untransferable to POSS, and it is **not** the paper above. The shape method used
here is the earlier, empirical one.

## How to validate

The per-plate counts are in
[`docs/poss2_red_control_per_plate.csv`](poss2_red_control_per_plate.csv),
one row per pair. Both headline figures come straight out of it, and the pooled
one is the consistency check — it must reproduce the chain's own totals:

```bash
python3 - <<'PY'
import numpy as np, pandas as pd
g = pd.read_csv("docs/poss2_red_control_per_plate.csv")
r = g.ratio
print("pooled              : %d / %d = %.4f" % (g.poss2_gated.sum(), g.poss1_s0.sum(),
                                                g.poss2_gated.sum() / g.poss1_s0.sum()))
print("per-plate mean +/- se: %.4f +/- %.4f (n=%d)"
      % (r.mean(), r.std(ddof=1) / np.sqrt(len(r)), len(r)))
PY
```

Expect `644 / 5939 = 0.1084` and `0.1243 +/- 0.0152 (n=30)`. If n comes out as
29, a zero-candidate plate has been dropped somewhere upstream — see the trap
above.

Rebuilding the POSS-II arm from scratch needs the public scans and the same
entry points as any run on this pipeline, with the pixel-valued parameters
scaled by 1.682 and two ingest accommodations: the slicer's filename pattern is
`dss1red_XE*`, so POSS-II plates are passed to `tools/slice_plate_tiles.py`
directly rather than through `tools/run_fullscale_slice.py`, and
`--allow-missing-crpix` is required because the CRPIX table has no XP rows (with
WCSFIX on, the per-tile Gaia refit absorbs the offset).

```bash
tools/slice_plate_tiles.py --plate-fits <dir>/dss2red_XP140.fits \
    --tiles-dir <out> --allow-missing-crpix
python -m vasco.cli_pipeline step2-pass1 --workdir <tile>
python -m vasco.cli_pipeline step3-psf-and-pass2 --workdir <tile>
VASCO_PLATE_EPOCH_YEAR=<poss2 epoch> python -m vasco.cli_pipeline step4-xmatch --workdir <tile>
```

POSS-II red full-plate scans: `https://irsa.ipac.caltech.edu/data/DSS/images/dss2red/`
— `XP` prefix only.

The POSS-I side needs no rerun: it is the released
`results/s0-642-paper-parity-20260828/` and
`results/s1-642-paper-parity-edge15-20260908/` restricted to the 30 partner
plates in the CSV, and both releases ship their own hashes.
