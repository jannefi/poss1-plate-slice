# S1 — paper-parity S0 with the plate-edge cut applied (2026-09-08; revision 2, 2026-09-09)

**95,826 rows.** The paper-parity S0 catalogue
([`../s0-642-paper-parity-20260828/`](../s0-642-paper-parity-20260828/),
134,976 rows) with every row lying within **15 arcminutes of its own plate's
array boundary** removed — 39,010 rows, 28.90% — **and, since revision 2,
every row of two tiles whose astrometric refit is known to be broken** — a
further 140 rows.

## Revision 2 (2026-09-09) — two tiles excluded

Revision 1 (2026-09-08, 95,966 rows) carried 137 rows of
`tile_RA24.686_DECp33.529` (plate XE296) and 3 of `tile_RA63.809_DECp57.375`
(plate XE084). On both tiles the per-tile Gaia refit that precedes the vetoes
fitted badly — 312 tie points at σ 0.48″ on the first, σ 0.88″ on the second,
both reported `ok` — and extrapolates ~7″ wrong across the tile, so real
stars there clear the 5″ vetoes and survive as rows: **104 of the 137 kept
XE296 rows sit 3–8″ from a Gaia star** (13 would by chance). The whole tile
is excluded, not only the matched rows, because the polynomial is wrong for
every row on it. Mechanism, the survey-wide check (this is one bad tile, not
a population) and the pipeline guard now being added:
[`../../docs/WCSFIX_GUARD.md`](../../docs/WCSFIX_GUARD.md). The defect is
inherited from S0, whose README carries the same appendix; S0's files are
unchanged.

The release was unannounced, so it is replaced in place rather than re-issued
under a new name. `edge_flags.csv.gz` is byte-identical to revision 1 (it
audits all 134,976 S0 rows and is unaffected by the exclusion);
`stage_S1.csv.gz`, `stage_S1_ledger.json`, `RUN_SUMMARY.txt` and both sums
files changed. Revision 1's content hash was
`8ef7e593dd2b66d5c3427d07d712ab3d17724d914d02b6560f73df02434730fe`; the
command under "Reproducing it" reproduces it exactly when `--exclude-tiles`
is omitted.

| file | rows | sha256 of the uncompressed content |
|---|---:|---|
| `stage_S1.csv` (revision 2) | 95,826 | `57f15791fe945d23feb6e2715ec85e5fd69883dc3ddc0b78dd463db02bda350f` |

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
| `stage_S1.csv.gz` | 95,826 | the cleaned catalogue (revision 2) |
| `edge_flags.csv.gz` | 134,976 | per-row audit for **all** S0 rows: plate and edge distance, so the cut can be checked or re-cut at another threshold; unchanged since revision 1 |
| `stage_S1_ledger.json` | — | parameters and counts as the stage recorded them, including the excluded tiles and their row counts |
| `SHA256SUMS` | — | integrity of the files as shipped |
| `SHA256SUMS.uncompressed` | — | identity of the catalogue content |

`stage_S1.csv.gz` columns: `src_id, tile_id, object_id, ra, dec, plate_id,
edge_dist_arcmin`. The first five match S0 exactly; the last two are added so
the cut is self-documenting. Minimum `edge_dist_arcmin` in this file is
15.000 by construction, and no row carries either excluded `tile_id`.

## Reproducing it

```
python scripts/stage_edge_post_v2.py \
    --run-dir <run> --input-glob 'stages/stage_S0.csv' \
    --plate-map-csv <tile_id,plate_id> \
    --stage S1 --cut --core-radius-deg 5.0 --min-edge-arcmin 15 \
    --exclude-tiles tile_RA24.686_DECp33.529,tile_RA63.809_DECp57.375
```

The stage writes `src_id, ra, dec` for the kept rows plus a flags file; the
shipped `stage_S1.csv` is those `src_id`s joined back to S0 for `tile_id` and
`object_id`, with `plate_id` (the flags' `det_plate`) and `edge_dist_arcmin`
appended. Without `--exclude-tiles` the same command reproduces revision 1.

`--core-radius-deg 5.0` is deliberate: it opens the radial core cut past the
plate corner distance (~4.68°) so it removes nothing, leaving the edge
distance as the only active criterion. The stage is **off by default and not
wired into the pipeline** — applying it is an explicit per-run decision, and
it prints a warning when active, because the papers describe no such
cleaning step.
