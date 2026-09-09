# S1 — paper-parity S0 with the plate-edge cut applied (2026-09-08)

**95,966 rows.** The paper-parity S0 catalogue
([`../s0-642-paper-parity-20260828/`](../s0-642-paper-parity-20260828/),
134,976 rows) with every row lying within **15 arcminutes of its own plate's
array boundary** removed — 39,010 rows, 28.90%.

This is a **derived, cleaned** product. It is published because the cut is
cheap to describe and expensive to reproduce: regenerating it needs the full
plate scans and a per-plate WCS pass, which not everyone has the disk or the
time for. The numbers behind it were already published as figures; this is
the actual file.

**It does not supersede anything.** S0 remains the catalogue to cite. This
is S0 minus a geometric region, offered for anyone who wants the
edge-cleaned set without rebuilding it.

## Why 15′, and what the cut removes

Distance to the plate's *array boundary* is a far sharper lever than
distance from the plate centre: a comparable radial cut removes ~65% of S0,
while the edge cut removes 21% at 10′ and 28.9% at 15′. Full derivation in
[`../../docs/PLATE_EDGE_MASK.md`](../../docs/PLATE_EDGE_MASK.md).

**What the evidence does and does not show, per annulus** (all 134,976 S0
rows; SuperCOSMOS from the same build's post-process flags, R1 arm):

| distance to array edge | share of S0 | median FWHM | median SPREAD_MODEL | SuperCOSMOS-unconfirmed | match to the public `vanish-possi` list |
|---|---:|---:|---:|---:|---:|
| <5′ | 6.8% | **3.87** | **0.011** | **97%** | 0.00% |
| 5–10′ | 12.5% | 3.01 | 0.005 | 91% | 0.01% |
| 10–15′ | 9.6% | 3.18 | 0.003 | 52% | 0.00% |
| 15–20′ | 7.8% | 3.14 | 0.003 | 22% | 0.15% |
| 20–30′ | 13.0% | 3.07 | 0.003 | 17% | 0.98% |
| >30′ | 50.3% | 3.04–3.11 | 0.0035 | 17–22% | 1.5% |

1. **The innermost 5′ is instrumental.** It is the only band where the
   survivors themselves look different — fuzzier (FWHM 3.87 vs 3.0–3.2) with
   three times the SPREAD_MODEL — and it holds the plate-fog and label-text
   detections documented in `PLATE_EDGE_MASK.md`.
2. **The rim carries a real excess.** Detections that survive both catalogue
   vetoes are 2.16× more frequent within 15′ of an edge than in the interior
   (Fisher p = 3×10⁻⁷), and before the morphology gate that rim population is
   measurably fuzzier and more elongated. A real transient cannot know where
   a plate edge is, so the *excess* is instrumental. It is an excess over the
   interior rate, though: a geometric cut removes the interior-like rows in
   the same band with it, and after the gate morphology no longer separates
   them.
3. **Beyond ~5′ the survivors look like the interior.** At 10–15′ FWHM,
   ELONGATION and SPREAD_MODEL match the plate interior and SNR_WIN is
   higher. Cutting that band is a **conservative margin**, not a finding that
   its rows are instrumental.

## What the recall and SuperCOSMOS figures can and cannot say

Two figures that look like independent support are not, for the rows this
cut targets:

- **Recall against comparison catalogues is uninformative at the rim by
  construction.** A position-addressed cutout service (STScI's `getimage`;
  ESO's server showed the same hand-off in our probes) hands any position near
  one plate's edge to the neighbouring plate that keeps it furthest from an
  array boundary —
  see the `FURTHEST_FROM_EDGE` note in
  [`../../docs/SCAN_SOURCE_SENSITIVITY.md`](../../docs/SCAN_SOURCE_SENSITIVITY.md).
  A catalogue built that way never examined any plate's outer ~18′, so a
  genuine single-plate event on our plate's rim would not be in it either.
  The public list's match rate steps from ~0% to 0.15% to 1% across
  10–15′ / 15–20′ / 20–30′ exactly where such a hand-off would put it; the
  "zero recall loss" inside 15′ measures that geometry, not our rows.
- **SuperCOSMOS is not an independent scan of the same glass at the rim.**
  The merged `supercosmos.sources` table carries no plate id, but its
  single-plate (R1-only) sources reveal the plate epoch: along the array edge
  of XE181 and XE002 the catalogued POSS-I E detections inside ~10′ belong to
  the *neighbouring* plate (its epoch), switching to the plate's own epoch
  from ~15′ inward. "SuperCOSMOS-unconfirmed" inside ~10′ therefore means
  "not a persistent source", which a single-plate event also is. Two plates
  probed; the switch distance is from the mid-latitude one.
- **"Absent from the overlapping neighbour"** describes every single-plate
  candidate, so it does not separate an artifact from an event.

Two alternative explanations of the rim excess were tested and **rejected**:
our own cross-plate dedup tie-break (the pre-dedup catalogue shows the same
pattern) and a quality gradient among survivors (morphology is flat outside
5′ because the MNRAS gate has already flattened it).

## Why 15′ and not 10′

The threshold was read off the comparison-list match-rate curve at the point
where its step completes. Given the above, that point is the cutout service's
hand-off distance rather than a property of this catalogue, and the
morphology and SuperCOSMOS figures would have supported 10′ as readily. The
release keeps 15′ as a conservative margin; `edge_flags.csv.gz` carries the
edge distance of every S0 row, so the cut can be re-drawn at 10′ (21%) or 5′
(6.8%) without the plate scans. Nothing here changes S0.

## Files

| file | rows | what |
|---|---:|---|
| `stage_S1.csv.gz` | 95,966 | the cleaned catalogue |
| `edge_flags.csv.gz` | 134,976 | per-row audit for **all** S0 rows: plate and edge distance, so the cut can be checked or re-cut at another threshold |
| `stage_S1_ledger.json` | — | parameters and counts as the stage recorded them |
| `SHA256SUMS` | — | integrity of the files as shipped |

`stage_S1.csv.gz` columns: `src_id, tile_id, object_id, ra, dec, plate_id,
edge_dist_arcmin`. The first five match S0 exactly; the last two are added so
the cut is self-documenting. Minimum `edge_dist_arcmin` in this file is
15.000 by construction.

## Reproducing it

```
python scripts/stage_edge_post_v2.py \
    --run-dir <run> --input-glob 'stages/stage_S0.csv' \
    --plate-map-csv <tile_id,plate_id> \
    --stage S1 --cut --core-radius-deg 5.0 --min-edge-arcmin 15
```

`--core-radius-deg 5.0` is deliberate: it opens the radial core cut past the
plate corner distance (~4.68°) so it removes nothing, leaving the edge
distance as the only active criterion. The stage is **off by default and not
wired into the pipeline** — applying it is an explicit per-run decision, and
it prints a warning when active, because the papers describe no such
cleaning step.
