# The plate-edge mask — an instrument, applied to nothing

`scripts/stage_edge_post_v2.py` measures where each catalogue row sits on its
plate and records it. It can also cut on that, but **nothing in this repository
has been cut**. The released catalogue is unchanged and its hashes stand.

That posture is deliberate, and it copies the source. Villarroel et al. 2025
(PASP) apply a 2° plate-centre mask inside two specific analyses while leaving
their catalogue at full size. We do the same: publish the instrument and the
measurement, leave the data alone, and let anyone who wants the cut apply it
themselves.

---

## Why a plate-edge mask exists at all

MNRAS 2022 says nothing about plate edges. **PASP 2025 does**, describing the
transient sample as

> expected to contain a substantial number of false positives, including
> clustered artifacts such as edge fingerprints or other plate defects that
> contaminate our sample

and acting on it twice: masking "edge transients (>2° from plate center)" in
the solar-reflection test, and counting "the transients within 2° from the
center to avoid potential defects on the edge of the plate (22,314
transients)" for the background-density estimate.

So 2° is **their** number, not one we tuned. `--pasp2025` applies exactly it.

### The 2° reading: 22,314 is what they kept

The sentence admits two readings — 22,314 could be the count retained inside 2°
or the count discarded outside it — and the difference matters, so it is worth
pinning down. **It is the retained subset**, and the paper settles this on its
own terms: the solar-reflection test describes masking "edge transients (>2°
from plate center)", which places the *discarded* population outside 2°, and
the background-density passage then counts "the transients within 2° from the
center". The parenthetical count belongs to the within-2° set. A discarded
reading would also leave roughly 85,000 transients for a background-density
estimate explicitly described as avoiding plate edges, which is not coherent.

Corroboration comes from an independent replication of the VASCO shadow
statistics — Doherty, *Independent Replication of Nuclear Test-Transient
Correlations and Earth Shadow Deficit in POSS-I Photographic Plates*,
[arXiv:2604.00056](https://arxiv.org/abs/2604.00056) — which analysed both the
full catalogue and a 2° centre-of-plate subset, stating that the restriction
"eliminates potential edge artifacts from plate scanning". Its supporting code,
cited in the paper's Data Availability section as
`github.com/dca-doherty/VASCO-Replication.git`, includes a run over a dataset
named `SUPERVIKTIG_HELAVASCO_within2deg_CENTER.csv`, described in the script
header as "pre-filtered to transients within 2 degrees of their plate center",
whose results file records `Plate field radius: 2.0 deg`, 614 plates and
`Total transients: 22309` — five rows from the published figure, and
unambiguously a retained subset.

**Two caveats on that corroboration, both worth stating.** First, the
repository is no longer publicly accessible: it returns 404 and no Internet
Archive snapshot exists, so a reader cannot check it today. The paper itself
remains citable and is the primary reference here. Second, the 22,309 file was
pre-filtered by the original team rather than recomputed independently, so it
evidences *what they meant* rather than providing a separate measurement.

**Consequence, on their own published figures:** those 22,314 are what survives
the 2° criterion out of the full transient list, which is roughly five times
larger — so their own edge criterion sets aside about **four fifths** of the
catalogue whenever an analysis depends on it.

### "Within 2° of plate center" is under-specified

Doherty's paper computes its own centre-of-plate subset by "unit-vector
averaging of source positions per plate" and retains **31,525 transients
(29.2%)** — against the ~20.7% implied by the published 22,314. The same
stated criterion therefore yields materially different subsets depending on
whether "plate center" means **the centre of the plate** or **the centroid of
the transients found on it**. Those coincide only if the transients are
distributed symmetrically across the plate, which is precisely what is in
question when the concern is edge artifacts.

`--pasp2025` uses the plate's actual centre, taken from the scan's WCS. Anyone
reproducing a published edge-masked count should state which definition they
used; it is worth roughly nine percentage points of yield here.

---

## What the mask is measured against

Two distances per row, both recorded, neither one privileged:

- **`sep_own_deg`** — angular separation from the centre of the plate the row
  was detected on. This is what a plate-centre mask means, and what
  `--pasp2025` cuts on.
- **`edge_dist_arcmin`** — distance to the physical boundary of the scanned
  array, in the plate's own pixel grid via its real WCS. A different and much
  more targeted quantity; see below.

The plate centre comes from the **scan's WCS at its middle pixel**, not from
the `PLATERA`/`PLTRAH` keywords. On eleven of 932 DSS1-red scans those differ,
seven of them by ~4.4°, because the image is not centred on the plate while
`CNPIX` still reads (0,0). Using the keyword there silently places on-plate
rows outside the array.

Note also that a plate is **square** (half-width ~3.31°, corners at ~4.68°)
while a centre-distance mask is **circular**. The two describe different
regions; `--pasp2025` is circular because that is what the source specifies.

---

## Running it

Flag only — writes the geometry for every row and keeps them all:

```
python scripts/stage_edge_post_v2.py \
    --run-dir <run> --input-glob 'stages/stage_S0.csv' \
    --plate-map-csv <tile_id,plate_id CSV> --stage S1
```

Apply the PASP 2025 criterion, this time actually cutting:

```
python scripts/stage_edge_post_v2.py \
    --run-dir <run> --input-glob 'stages/stage_S0.csv' \
    --plate-map-csv <tile_id,plate_id CSV> \
    --pasp2025 --policy own --cut --stage PASP
```

Other presets: `--vasco60-parity` (2.907° = 2.2° tile-centre limit plus a
half-tile diagonal), or `--core-radius-deg` for anything else. The ledger
always carries a yield curve across candidate radii and across distance to the
array edge, so the cut point stays an explicit decision. `--self-test` runs the
geometry checks with no pipeline state.

---

## Measured impact

### The PASP 2° mask on this catalogue

| stage | rows | after 2° mask | retained |
|---|---:|---:|---:|
| `stage_S0` (released) | 122,820 | 15,798 | **12.9%** |
| after SuperCOSMOS + PTF | 71,654 | 11,322 | **15.8%** |

For scale, the same criterion retains 22.9% of the published vanish-possi
catalogue (R), and about a fifth of the VASCO transient list per the figures
above.

**This catalogue is rim-heavy by construction, and that is not a quality
statement.** Full-plate slicing tiles the entire plate including its corners;
a cutout-based pipeline queries one position at a time and never generates
those tiles. A lower core fraction here is the tessellation, not the
candidates.

### SuperCOSMOS already removes rim rows preferentially

The core fraction rising from 12.9% to 15.8% across the post-process chain is
the interesting part. Splitting what the chain removed:

| | rows removed by SCOS+PTF | removal rate |
|---|---:|---:|
| outside 2° | 46,690 | **43.6%** |
| inside 2° | 4,476 | **28.3%** |

A 1.54× differential. An independent digitisation of the same glass discards
rim rows half again as often as core rows — arrived at with no geometry in the
loop, since SuperCOSMOS knows nothing about where our plate centres are.

### Cost in published-catalogue recall

| set | matches to R inside 2° |
|---|---:|
| S0 | 242 of 1,234 (19.61%) |
| S0 + PASP 2° mask | 241 of 1,234 (**19.53%**) |
| after SCOS+PTF | 235 (19.04%) |
| after SCOS+PTF + PASP 2° | 234 (**18.96%**) |

Inside the region the mask keeps, applying it costs **one** R match. Core rows
are matched by core rows. Recall against *all* of R falls from 19.86% to 4.46%,
but that is the mask doing its job — 77% of R lies outside 2°, and removing
that sky necessarily removes those matches.

### A sharper cut than radius-from-centre

Distance to the array boundary is far more targeted than distance from the
centre:

| cut | removes from S0 | cost in R matches |
|---|---:|---:|
| radius > 2.7° from plate centre | 65% | substantial |
| within 10′ of the array edge | 21% | **zero** |
| within 5′ of the array edge | 7.6% | **zero** |

The edge cut is nearly free because overlapping full-plate coverage keeps the
sky on a neighbouring plate's interior — the detection is dropped, the sky is
not. On the post-chain set only 0.36% of rows sit within 5′ of an edge, against
7.6% in S0, so **SuperCOSMOS has already removed that population**; the
geometric cut is a cross-check on it rather than an additional filter.

---

## Conclusions

1. **2° is the source's own criterion**, confirmed by two independent sources,
   and on their own published figures it sets aside about four fifths of their
   catalogue whenever it is applied.
2. **Nothing here has been cut.** The instrument is published; the release is
   untouched. Matching the source's posture is the point.
3. **Rim-heaviness is tessellation, not quality.** Full-plate slicing generates
   corner sky that a cutout pipeline never visits.
4. **SuperCOSMOS and plate geometry agree without sharing information.** The
   chain removes rim rows 1.54× as often as core rows, and the near-edge
   population it removes is the same one a geometric cut targets.
5. **Distance to the array boundary is the better lever** — 21% of rows for
   zero cost in published-catalogue recall, against 65% for a comparable
   radial cut.

## Caveats

- A centre-distance mask is a **statement about astrometric reliability**, not
  about whether a detection is real. The APS documentation restricts reliable
  POSS-I astrometry to the central 5.4°; it does not say sources outside it are
  false. Treating it as a veto over-reads the source, which is why the default
  here is to flag.
- The rows a mask removes are **not uniformly junk**. Sub-arcminute-from-edge
  rows are unambiguous plate fog, but a random sample of the 2.7–3.0° band
  looks like ordinary sky, and at 1′ zoom most markers sit on a compact round
  source. Dust specks and grain clumps also look stellar at 1.7″/px, so visual
  inspection cannot settle this either way.
- `cheb_own_deg` in the flags output is a linear-gnomonic approximation and
  drifts ~1′ from the GSSS polynomial at the rim. Use `edge_dist_arcmin` for
  boundary questions; it is computed from the plate's own WCS.
