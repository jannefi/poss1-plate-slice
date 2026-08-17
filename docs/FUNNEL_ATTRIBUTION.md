# Funnel attribution — which stage removes the published sources?

This pipeline detects almost every source in the published POSS-I
vanishing-source catalogue (SVO `vanish-possi`, Solano et al. 2022, 5,399 rows)
at the raw-detection stage, and keeps about a fifth of them in the final
candidate list. This page measures *which stage* removes the rest — veto by
veto, gate by gate — rather than inferring it from aggregate counts.

Everything here uses public inputs only: public plate scans, the public
reference catalogue, and this repository's own pipeline output. Tools:
[`tools/funnel_attribution.py`](../tools/funnel_attribution.py),
[`tools/gaia_epoch_check.py`](../tools/gaia_epoch_check.py) and
[`tools/usnob_pm_quality.py`](../tools/usnob_pm_quality.py).

**Sample sizes are small, and that governs how these numbers may be used.** The
reference catalogue has 5,399 rows spread over 642 plates, so any subset of
tiles holds few of them: 218 rows in the main sample here, 112 of which are
losses to attribute. Shares below are reported to the precision the counts
support and no further. Nothing on this page should be quoted with a tight
interval.

## Method

Every per-tile stage catalogue carries the SExtractor `NUMBER` column, unique
within a tile and strictly nested down the chain. So each reference row is
matched to its nearest raw detection **once**, at 5″, and every later question
is set membership on that identifier. Buckets sum to the in-footprint row count
by construction, and the tool asserts it.

Two measurements, both on tile sets this project retained with their full
per-stage chain:

* **18-plate sample** — 882 tiles, plates chosen to span declination and
  galactic latitude evenly (XE025, XE053, XE089, XE116, XE137, XE215, XE235,
  XE241, XE311, XE320, XE340, XE435, XE449, XE452, XE539, XE567, XE582, XE623).
  Median plate centre declination 30.0° against the survey's 29.8°, so this is
  the survey-representative arm.
* **504-tile high-declination set** — δ 41.5–86.6, retained from an unrelated
  repair. Supplementary, and useful precisely because it is *not*
  representative: see the epoch result below.

## Detection: ~99% inside plate cores

| arm | rows in core | raw detection | survives to S0 |
|---|---:|---:|---:|
| 18-plate | 72 | **98.61%** | 19.44% |
| 504-tile | 17 | 94.12% | 17.65% |

"Core" means rows at least 1° inside every edge of their plate (≤2.25° from the
plate centre, whose exact position comes from the public manifest). **The
restriction is necessary, not cosmetic.** Measured over each arm's full
footprint the same detection rates read 66.97% and 63.08% — but those are an
artifact of sampling tiles rather than the whole sky: a reference row near the
edge of a sampled plate often has its detection on an overlapping tile that
was not sampled, and scores here as undetected. The conditional attribution
among *detected* rows is unaffected by this, which is why the sections below
use the full footprint while the detection rate uses cores.

The ~20% survival rate agrees with an independent whole-survey crossmatch of
the released catalogue against the reference catalogue: **1,072 of 5,399
(19.86%)** have an S0 row within 5″.

## Attribution: the quality filters dominate, not the vetoes

18-plate arm, 112 losses:

| bucket | rows | share of losses |
|---|---:|---:|
| Gaia veto | 2 | 1.8% |
| PS1 veto | 2 | 1.8% |
| USNO-B veto | 15 | 13.4% |
| **veto chain total** | **19** | **17.0%** |
| MNRAS extraction gates (`FLAGS`, `SNR_WIN`) | 40 | 35.7% |
| MNRAS morphology gates | 53 | 47.3% |
| **MNRAS filters total** | **93** | **83.0%** |
| spike mask / global dedup | 0 | 0% |

**The USNO-B row does not survive scrutiny.** It is a defect in this project's
own veto stage rather than a property of the reference catalogue — see
[the USNO-B section](#the-usno-b-veto-fires-on-fabricated-proper-motions)
below. Read the veto-chain total as **4 rows (3.6%)**, not 19, and the filter
share as correspondingly larger.

Which gate actually fired, read from the `reject_reason` column rather than
re-derived:

| reason | rows | share of filter losses |
|---|---:|---:|
| `spread_model` | 34 | 36.6% |
| `flags` | 22 | 23.7% |
| `snr` | 18 | 19.4% |
| `spread_model` + `fwhm` | 9 | 9.7% |
| `elongation` | 6 | 6.5% |
| others (combinations) | 4 | 4.3% |

**Two thirds of the losses are gates published in Solano et al. (2022)** —
`SPREAD_MODEL > −0.002`, which rejects sources sharper than the PSF, plus the
extraction requirements `FLAGS == 0` and `SNR_WIN ≥ 30`. That the extraction
gates alone remove 43% of the filter losses is the part worth investigating:
these are rows that passed the reference pipeline's own SNR threshold and fail
ours. Different cutout pixels, a different background estimate, or blending
under this project's tessellation are all plausible; 58 rows is a lead, not a
finding.

**The ordering caveat travels with every share above.** The vetoes run before
the filters, and the morphology sigma-clip window is derived from the
population being filtered — so this is what the pipeline did in its actual
order, not a counterfactual. Removing a veto would change the clip window and
therefore change which rows the *filters* cut.

## The Gaia losses are stars that moved

The pipeline vetoes against Gaia positions propagated to the plate epoch
(~1950s), not Gaia's catalogue positions (2016). On the high-declination arm,
where the Gaia bucket is large enough to test (11 rows):

| | within 5″ |
|---|---:|
| **plate epoch** — what the veto used | **81.82%** |
| **catalogue epoch** — Gaia's own positions | **0.00%** |
| null control (positions displaced 6′ in RA) | 0.00% |

Nine published reference rows have a Gaia star sitting at their position at the
plate epoch and none at the catalogue epoch. Those are **objects that moved,
not objects that vanished** — precisely the contaminant class a
vanishing-source search must remove, and a mechanism that a catalogue matched
at Gaia's own epoch would not catch.

On the survey-representative arm the Gaia bucket is 2 rows and neither is
proper-motion-explained, consistent with high-proper-motion stars being
relatively more prominent at high declination and high galactic latitude.
Whether the effect generalises across the survey is not settled by 2 rows.

**n = 11 is small, and the claim does not rest on the percentage.** Each of the
nine is an individually checkable published coordinate: take the reference
position, query Gaia DR3, propagate the proper motion back ~65 years, and the
star is there. No part of this needs this pipeline.

## The USNO-B veto fires on fabricated proper motions

USNO-B is **this project's own addition** to the veto chain — Solano et al.
(2022) use Gaia and Pan-STARRS — and it removes 15 reference rows, 13.4% of all
losses and 79% of the veto chain's total. An earlier version of this page left
that as an open question. It has now been inspected, and the answer is that
**the removals are spurious and the stage is circular with the plates it runs
on.**

### What the epoch propagation does

The pipeline propagates USNO-B positions to the plate epoch before the 5″ veto
cone. That is correct in principle: a genuine high-proper-motion star sits
10–25″ from its J2000 position on a 1950s plate, and matching at J2000 would
miss it.

On the 15 rows (plus one from a second tile set) the propagation is the entire
story:

| | at epoch 2000.0 | at plate epoch |
|---|---:|---:|
| matched within 5″ | **0 / 16** | **15 / 16** |
| RA+5° null control | 3 / 16 | 3 / 16 |

Nearest USNO-B source at 2000.0 is 5.0–25.5″ away; after propagation the median
separation is **0.62″**. The null control does not move between epochs at all
(3.13→3.13, 3.23→3.23, 4.84→4.84), because chance neighbours have no proper
motion — real firings converge, chance firings do not.

### The proper motions are not real

The veto invoked motions of **186–589 mas/yr**. Across all 4,549 Gaia DR3
sources within 3′ of those sixteen positions, the fastest star moves at
**101 mas/yr**. For 12 of the 16, Gaia has a source within 5″ of the USNO-B
*catalogued* position moving at **0–13 mas/yr** — a discrepancy of 26–238×.

That generalises, and not narrowly. Sixty random 12′ fields inside the plate
footprint, 7.54 deg², 150,333 USNO-B entries
([`tools/usnob_pm_quality.py`](../tools/usnob_pm_quality.py)):

| | count in 7.54 deg² | per deg² |
|---|---:|---:|
| Gaia DR3 sources with PM ≥ 150 mas/yr | **31** | 4.1 |
| USNO-B entries claiming PM ≥ 150 mas/yr | **9,200** | 1,220.2 |

**USNO-B overclaims fast stars by 297×.** Testing each entry symmetrically —
*fabricated* if a slow Gaia star sits at its catalogued J2000 position,
*genuine* if Gaia sits where its own proper motion predicts for 2016.0 at a
comparable rate — only **0.3%** of the high-PM entries are confirmed genuine,
against 30.7% demonstrably stationary. The 29 confirmed genuine correspond
closely to the 31 real fast stars Gaia finds. A low-PM control stratum leaves
only 18.6% indeterminate against 69.0% at high PM, so this is not the
association test failing.

The symmetry matters: a naive "is Gaia at the J2000 position" test is biased,
because a real 500 mas/yr star has moved 8″ by Gaia's epoch and would fail it
purely by moving, leaving only stationary stars in the associated set.

### Why the matches land so precisely — the circularity

USNO-B records which surveys back each entry. **All 16 carry an `R1` leg:
POSS-I red, the same glass this pipeline runs SExtractor on.**

USNO-B has paired a POSS-I detection with an unrelated POSS-II detection and
fitted a proper motion to the pair. Back-propagating by that fitted motion
sweeps the entry 8–25″ across the sky and reproduces the R1 position — the
POSS-I detection — to sub-arcsecond precision, because that detection is an
input to the fit. **The stage discards a POSS-I detection on the evidence that
USNO-B recorded that same POSS-I detection.**

Across 9,052 real veto firings on 1,250 tiles retained with their full
per-stage chain:

| | firings | share |
|---|---:|---:|
| static match at J2000 — propagation irrelevant | 5,121 | 56.6% |
| **fired only because of propagation** | **3,931** | **43.4%** |

Of the propagation-dependent firings, **0.5% are genuine**, 28.0% are
demonstrably fabricated, 71.5% are indeterminate, and **89.6% are `R1`-backed**.
The 0.5% independently reproduces the 0.3% above. Median motion invoked:
309 mas/yr.

### Two limits on these numbers

The 1,250 tiles were selected around published reference positions, so their
firing rate is **not a survey rate**. And 71.5% indeterminate means the
fabricated fraction is a floor, not an estimate — what bounds the other side is
the 0.3–0.5% confirmation rate, which makes it unlikely the indeterminates are
real fast stars.

### What this means for the released catalogue

The released catalogue was produced with this veto active and propagating, so
these removals are inside it. The per-tile intermediate state for that run has
since been deleted, so the affected rows cannot be recovered by inspection and
would have to be recomputed by re-running.

Nothing has been altered in response. The released files and their hashes
stand, and this section is the disclosure rather than a correction to them.

**USNO-B is not part of the method being reproduced.** It was an optional extra
stage, and the planned paper-compliant run drops it — along with the other
deviations from Solano et al. (2022) — which removes the defect at the source
rather than patching it. Until that run exists, treat the USNO-B row of the
attribution table above as spurious and the veto chain as 4 rows, not 19.

## Reproducing

```bash
python3 tools/funnel_attribution.py \
    --tiles-root <run>/tiles \
    --ref-csv vanish_possi.csv \
    --s0-csv results/s0-642-20260814/stage_S0.csv.gz \
    --tile-plate-map <tile_plate_map.csv> \
    --label R --out funnel_rows.csv --out-summary funnel.json

python3 tools/gaia_epoch_check.py \
    --funnel-csv funnel_rows.csv --tiles-root <run>/tiles --out gaia_epoch.csv
```

The USNO-B proper-motion measurement needs no pipeline output at all — two
local catalogue mirrors and this repository's plate manifest are enough, so it
is checkable by anyone holding USNO-B1.0 and Gaia DR3:

```bash
python3 tools/usnob_pm_quality.py \
    --usnob-mirror <usnob parquet dir> \
    --gaia-mirror <gaia parquet dir> \
    --density-check
```

It is seed-stable (`--seed`, default 20260817). Both mirrors must be
hive-partitioned on a `healpix_5` column built at nside=32 in **nested** order;
read as ring, the queries silently return the wrong sky.

The tile sets used here are pipeline output, regenerable from the public plate
scans by the documented run procedure — a run must retain its per-tile
`catalogs/` directories, since the attribution reads the per-stage files
directly. The reference catalogue is public at SVO (`vanish-possi`).
