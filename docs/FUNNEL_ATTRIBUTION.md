# Funnel attribution — which stage removes the published sources?

This pipeline detects almost every source in the published POSS-I
vanishing-source catalogue (SVO `vanish-possi`, Solano et al. 2022, 5,399 rows)
at the raw-detection stage, and keeps about a fifth of them in the final
candidate list. This page measures *which stage* removes the rest — veto by
veto, gate by gate — rather than inferring it from aggregate counts.

Everything here uses public inputs only: public plate scans, the public
reference catalogue, and this repository's own pipeline output. Tools:
[`tools/funnel_attribution.py`](../tools/funnel_attribution.py) and
[`tools/gaia_epoch_check.py`](../tools/gaia_epoch_check.py).

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

## An open question about the USNO-B veto

USNO-B is this project's addition to the veto chain, and it removes 15
reference rows — 13.4% of all losses, and 79% of the veto chain's total. Worth
checking directly against USNO-B rather than assuming: this pipeline propagates
proper motion to the plate epoch before matching, which a match at USNO-B's
own epoch 2000.0 would not do, so the two need not agree. Fifteen rows is small
enough to inspect individually.

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

The tile sets used here are pipeline output, regenerable from the public plate
scans by the documented run procedure — a run must retain its per-tile
`catalogs/` directories, since the attribution reads the per-stage files
directly. The reference catalogue is public at SVO (`vanish-possi`).
