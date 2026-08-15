# Post-process stages — SuperCOSMOS, PTF, SkyBoT

What happens to the released catalogue when the post-process stages of Solano
et al. (2022) are applied to it, measured rather than assumed. **None of this is
applied to `results/s0-642-20260814/`** — every stage here writes flags, the
released catalogue is unchanged, and its content hashes remain valid.

All inputs are public: the released catalogue, the public plate scans, GAVO's
`supercosmos.sources`, and IRSA's `ptf_objects`.

Tools: [`tools/scos_stage_banded.py`](../tools/scos_stage_banded.py),
[`tools/ptf_stage_coverage.py`](../tools/ptf_stage_coverage.py),
[`tools/shrink_stage.py`](../tools/shrink_stage.py).

## The chain

| stage | rows in | rows out | removed |
|---|---:|---:|---:|
| S0 (released) | — | 122,820 | — |
| SuperCOSMOS (R1) | 122,820 | 73,681 | **−40.0%** |
| PTF (`ngoodobs > 0`) | 73,681 | 71,654 | −2.75% |
| SkyBoT | — | — | *not run — see below* |

Order does not affect the result. Each stage is a per-row catalogue test with no
population-derived threshold, so the final set is order-independent; only the
per-stage attribution shifts. (This is **not** true of the MNRAS morphology
filters, whose σ-clip window is derived from the population being filtered.)

The two stages join in opposite directions, and getting that backwards silently
produces the complement of the catalogue. SuperCOSMOS **keeps** matches — Solano
et al. (2022): *"Candidates having a counterpart in the Supercosmos digitization
at less than 5 arcsec were kept."* PTF **drops** them: a modern-epoch counterpart
means the source is still there.

## SuperCOSMOS: 40% of the catalogue is not confirmed by a second scan

An independent digitization of the same POSS-I E plates fails to confirm
**49,139 of 122,820 rows** at 5″.

### Which band, and why it matters

`supercosmos.sources` is a merge over four photographic bands — B, R1, R2 and I.
The I band is SERC/POSS-II photographic I (IV-N emulsion, ~715–900 nm), which
reaches into the near-infrared, and in a 1 deg² probe at δ+29, 573 of 21,078
sources (2.7%) are I-only. Matching the merged table therefore lets a
near-infrared-only detection satisfy a check that is supposed to be about red
POSS-I plates. The tool returns all three arms from one query so the choice is
measured, not asserted:

| arm | matched | unmatched |
|---|---:|---:|
| any band | 61.92% | 46,764 |
| optical only, I excluded | 61.19% | 47,669 |
| **R1 = POSS-I E only (used here)** | **59.99%** | **49,139** |

Excluding the I band costs **0.74%** (905 rows). R1 is also the right arm on the
merits: R1 *is* POSS-I E, the same photographic material digitized independently
from different glass copies, so this is the cross-scan consistency check the
method intends. R2 is POSS-II — a different epoch, which would test persistence
instead.

**What a non-match can and cannot mean.** For R1, a genuine source on the plate
must appear in both digitizations, modulo detection thresholds. So a non-match
means an artifact on one scan (dust, emulsion flaw, plate defect present on one
copy and not the other) *or* a rejection by SuperCOSMOS's own detection and
merge pipeline. It cannot mean "a real source only one scan saw". These
measurements do not separate those two cases.

### Four confounds, tested and excluded

| confound | test | result |
|---|---|---|
| no SuperCOSMOS coverage there | 800 unmatched positions, count sources within 60″ | **0 holes**; median 20 sources nearby |
| coverage gradient with declination | match rate per declination band | flat 56–70% |
| plate-edge vignetting depresses any scan | ordinary raw detections, rim vs core | 95.7% vs 99.5% — **3.8 pts** |
| catalogue rows are simply fainter | magnitude-matched, same plates | below |

The magnitude control is decisive. Catalogue rows and ordinary raw detections
have an identical median magnitude of −12.05, and in every bin:

| mag bin | catalogue rows | ordinary detections |
|---|---:|---:|
| −13.5..−12.5 | 63.7% (n=193) | 98.5% (n=1,189) |
| −12.5..−11.5 | 70.2% (n=554) | 98.2% (n=1,865) |
| −11.5..−10.5 | 72.4% (n=152) | 98.1% (n=1,273) |

**~30 points of deficit at matched brightness.** Catalogue rows are
qualitatively different from ordinary detections of the same magnitude on the
same plates — not a fainter version of them. (n=924 rows with recovered
magnitudes over 6 plates; modest, but the effect is far larger than its own
sampling noise and consistent across all three bins.)

## PTF: a null, and the null is the useful part

PTF reproduces the paper's quality gate, `COALESCE(ngoodobs,0) > 0`.

| | inside PTF coverage | match, of covered rows |
|---|---:|---:|
| `is_primary` | 82.9% | 3.46% |
| single-plate | 75.9% | 3.24% |
| all (73,681) | 80.9% | 3.40% |

**Coverage is measured, not assumed.** PTF is not all-sky: 19.1% of rows lie
where it never observed, falling to 65.8% coverage at δ 60–90. This creates no
false removals — for a drop-matches veto an untested row is kept, the
conservative direction — but a naive run would quote 2.75% while including
untestable sky in the denominator. 3.40% over covered rows is the interpretable
figure. (IRSA's ADQL rejects conditional aggregates, so coverage costs a second
simple query per chunk.)

## What the two stages together establish

The release ships a coverage partition
([`primary_plate_flags.csv.gz`](../results/s0-642-20260814/README.md)) splitting
the catalogue into sky a per-position pipeline also searches and single-plate
content in multiply-searched sky. SuperCOSMOS played no part in that rule, which
is pure geometry from published plate centres — so it is an independent test of
it:

| partition | SuperCOSMOS unconfirmed | PTF match (of covered) |
|---|---:|---:|
| `is_primary` (68,071) | **23.6%** | 3.46% |
| single-plate (54,749) | **60.5%** | 3.24% |

**SuperCOSMOS separates the partition sharply; PTF does not separate it at all.**
That is what should happen if the two measure different things: SuperCOSMOS asks
whether a row is a real feature on the plate, so an artifact-rich population
fails it; PTF asks whether a real object is still visible today, which has no
reason to track scan quality. The rim deficit (37 points) is ten times the
measured vignetting sensitivity effect (3.8 points).

So the single-plate partition is **disproportionately scan artifacts**, and the
evidence for that comes from a catalogue with no connection to the geometry that
defined it.

## SkyBoT: not run, with a measured rate

SkyBoT is epoch-aware, so cones cannot be batched by position — each plate epoch
needs its own query. Screening the stage set means **17,654 cone searches** on an
80′ grid over 641 distinct plate epochs. Calibrated on 100 randomly sampled
cells covering 393 rows:

| | |
|---|---|
| latency | median **4.67 s**, p90 14.85 s, max 32.8 s |
| reliability | 99/100 HTTP 200; **18 of 100 required a retry** |
| yield | **0 solar-system objects, 0 rows matched at 5″** |
| full pass | **~23 h** at median latency, ~73 h at p90 |

**Not run.** ~23–35 h against a service needing a retry on 18% of requests, for
a yield bounded above by **0.76%** (rule of three on 0/393) — under ~546 rows.
That is an upper bound, not a measured zero, and the sample covers 0.55% of the
stage set.

A plausible reason the yield is low: asteroids trail on a photographic exposure
and are likely removed by the `ELONGATION < 1.3` gate long before a solar-system
check sees them. That is inference, not measurement, though it is consistent
with the predecessor pipeline measuring exactly zero for this stage.

Reversible: per-plate epochs are derivable from the scan headers (all 932 carry
a clean `DATE-OBS`), and the runner caches every response so a pass can be
resumed rather than restarted.

## On comparing with the published funnel

Solano et al. (2022) report ~298K candidates after their pipeline, then
−288,770 for "sources in other catalogues", then SkyBoT 189, PTF 35,
**SuperCOSMOS 3,592**, high-proper-motion 178, visual 2 — leaving 5,234 against
the published 5,399. SuperCOSMOS is 39.9% of the pool it saw, close to the 40.0%
here.

**That agreement is not validation, and should not be read as any.** Two
reasons. Their "other catalogues" cut includes infrared data whose use was
subsequently retracted, so the 9,230 rows their SuperCOSMOS stage saw are
downstream of a step the authors themselves later stepped away from. And the
stages are not comparable: this catalogue corresponds to their **~298K** point —
post-pipeline, post-Gaia, post-PS1, no infrared — not to their 9,230. Two ~40%
figures on pools differing 13× in size and by an entire retracted veto.

What establishes the number here is the confound testing above, not the
coincidence.
