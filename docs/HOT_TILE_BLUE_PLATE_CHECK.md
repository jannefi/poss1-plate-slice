# Do hot-tile survivors show up on an independent exposure?

A small number of tiles (244 of 642 plates, at a ≥20-survivors/tile threshold) carry a disproportionate
share of catalogue survivors — the "hot tile" concentration. Four
explanations have been checked and ruled out for these tiles: un-vetoed Gaia
stars, astrometrically displaced catalogued stars, scan artifacts (they
survive an independent SuperCOSMOS digitization of the same plates), and
USNO-B-catalogued sources. One explanation was left standing without direct
evidence: a plate-emulsion or scan defect specific to the red (E) glass.
SuperCOSMOS can't test that hypothesis — it re-scans the *same* physical
glass, so a defect in the emulsion itself would survive that check too.

The POSS-I O (blue) plate is a genuinely independent exposure: different
glass, same night, exposed 8–10 minutes against the red plate's 45–60,
typically within about half an hour of it. A *persistent* astronomical
source should generally show at least some flux there; a defect specific to
the red emulsion should not — and, as the revisit below spells out, so
should nothing that lasts less than the gap between the two exposures. This
applies to both `results/s0-642-20260814/` and
`results/s0-642-paper-parity-20260828/` — the tiles checked here are drawn
from the primary release's own hot-tile ranking, and the mechanism being
tested is a property of the plate material, not of either build's veto
configuration.

## Method

`tools/blue_plate_hot_tile_check.py`. For each tile: fetch one 60′ STScI
`dss1-red` and one 60′ `dss1-blue` cutout at the tile centre
(`vasco.downloader.fetch_skyview_dss`). Local-background z-score per
position — a small aperture for the peak, a surrounding annulus for a
median/MAD background estimate, off-array positions excluded. Per image,
**polarity is resolved against the 40 brightest Gaia stars in the footprint**
rather than assumed. Every catalogue survivor in the tile, and a matched
population of random on-array null positions, gets a z-score on both bands.
Red is the positive control — survivors are red-plate detections by
construction, so they must show real flux there; blue is the actual
question. Gaia is used only to fix which pixel-value sign means "star" in
each image, never as a match or veto criterion.

**Positive-control cross-check, checked before trusting anything else**: the
40 bright Gaia stars used for polarity resolution show robust, strong blue
detection in every tile checked (median z 9.86–25.67, all resolved) —
roughly 40–45% of their own red-band strength, but nowhere near the null
floor. This rules out registration error, aperture mis-sizing, or "the blue
plate is just too shallow to matter" as explanations for what follows: real
stars are easily detected by this method on these exact images.

## Sample and results

Two passes. First, 6 tiles hand-picked for diversity (top tile per plate,
one per distinct plate, from the released catalogue's top hot tiles):
**XE296** (291 survivors), **XE074** (186), **XE407** (97), **XE491** (91),
**XE049** (85), **XE366** (83). Hand-picking meant the ratio it produced
couldn't be trusted as a population estimate on its own, so a second pass
drew 30 tiles genuinely at random (fixed seed, from the same 671-tile /
244-plate pool defined by the `rows_emitted_to_S0 ≥ 20` threshold above),
excluding the 6 already checked.

**Pilot (6 tiles)**: one tile (XE296) is a clear blue-positive exception —
63.5% of survivors show blue flux above a 12.7% null rate, and both bands
independently show a genuinely dense star field there. One tile (XE407) is
inconclusive — even the red positive control is weak. The remaining 4 tiles
combined: survivors show blue flux at **9.0%** (z>3) against a **16.0%**
null rate — less than what random field positions show, the opposite of
what an ordinary star population should do.

**Random 30-tile draw**: 9 of 30 tiles (30%) fail the red positive control
and are excluded — a real finding on its own, since the hand-picked pilot's
lower inconclusive rate (1 of 6) had understated how often even the red
signal is marginal once tiles aren't cherry-picked toward clear cases. Of
the 21 valid tiles: **17 (81%) show survivors below their own tile's null
blue-detection rate**, only 4 show an excess — a two-sided sign test on that
split gives **p = 0.0072**. Combined: blue survivors 9.04% (z>3) against a
14.90% null — reproducing the pilot's ratio almost exactly, now backed by an
unbiased sample and a real significance test rather than 6 hand-picked
tiles.

## Visual examples

![XE491: a moiré scan artifact under the survivor cluster](figures/hot_tile_blue_check/XE491_moire_artifact.png)

**XE491** — the cleanest case. The red image shows an obvious rippled
interference pattern running through part of the tile, clearly not a
stellar structure. Nearly every survivor sits inside that band. The blue
image shows an ordinary field in the same region, and almost none of those
positions show blue flux.

![XE296: real stars in a dense field, the exception](figures/hot_tile_blue_check/XE296_real_stars_and_stamp.png)

**XE296** — the exception, and since 2026-09-09 a known one. Red shows a
printed plate-label/stamp region and what looks like a genuinely denser
field; blue independently confirms a dense, ordinary field in the same
area, and most survivors show real blue flux. The first version of this
page read that as sky density. It is not: this tile
(`tile_RA24.686_DECp33.529`) is the one whose per-tile astrometric refit
failed — 312 tie points, σ 0.48″, reported `ok` — so that ~85% of its
survivors are ordinary Gaia stars displaced by a coherent ~7″ past the 5″
vetoes. They show blue flux because they *are* stars. See the "Known
defect" appendix of either S0 README and
[`WCSFIX_GUARD.md`](WCSFIX_GUARD.md).

![XE366: the known title-text tile](figures/hot_tile_blue_check/XE366_titletext_no_blue.png)

**XE366** — a printed plate-title-text band runs across the tile. Survivors
sit on what look like real point sources in red; almost none show blue flux
at the same positions.

![XE622: a scan/mosaic-seam strip](figures/hot_tile_blue_check/XE622_edge_seam_zero_blue.png)

**XE622** — most survivors sit in the main field on plausible real sources,
but several sit in a distinct bright vertical strip at the tile's edge (a
scan/mosaic-seam region), all showing zero blue flux.

![XE050: a crowded field, zero blue](figures/hot_tile_blue_check/XE050_crowded_field_zero_blue.png)

**XE050** — a very dense field. Survivors sit on plausible stars in red;
zero show blue flux, despite blue showing an equally dense field nearby.
Flagged below as the source of a real, unresolved caveat.

## Revisited 2026-09-09 — the blue result is not specific to hot tiles

The check above compared hot-tile survivors with random *positions* in the
same images. It never compared them with *ordinary* survivors. That
comparison has now been made, with the same method (5″ core, 25″ background
box, four 45″-offset nulls, polarity from each image and confirmed against
Gaia), on catalogue rows that are not hot-tile rows:

| rows | n | red z>5 (control) | blue z>3 vs null | blue z>5 vs null | red-only |
|---|---:|---:|---|---|---:|
| hot-tile survivors (this page, 21 random tiles) | — | — | 9.0% vs 14.9% | — | — |
| S1 rows that also appear in the published 5,399-row list (Solano et al. 2022) | 294 | 100% | 22.4% vs 25.5% | **4.1% vs 10.2%** | 77.6% |
| S1 rows sampled across coverage classes | 889 | 99–100% | 19.5% vs 20.1% | **2.2% vs 6.5%** | 77–81% |

Same field and same night verified per row from the served plate's header
(99% and 96–100%). **Every class of survivor is at or below its own null in
blue** — hot tiles, rows the published papers themselves list, and the rest
of the catalogue alike. "Fewer blue counterparts than chance" is therefore a
property of essentially everything that passes the vetoes and the MNRAS
gates, and it cannot be evidence about hot tiles in particular. (Rows sit
*below* chance because vetoed positions are star-depleted by construction —
the same effect seen against Gaia neighbours.)

The inference also over-reached in a second way. "A real source should show
some blue flux" is true of persistent sources — which the Gaia/PS1 vetoes
have already removed. For what remains, blue absence is exactly what the
source papers' own target class produces: Solano et al. (2022) report that
only 2 of their 5,399 candidates were seen in both colours and attribute the
scarcity to colour and to the 8–10 vs 45–60 min exposures. An emulsion or
scan defect on the red glass and an event shorter than the gap between the
two exposures predict the same thing, present in red and absent in blue, and
this test cannot tell them apart.

## Conclusion (revised 2026-09-09)

The measurement stands: on a random draw, 81% of checkable hot tiles show
fewer independent-exposure counterparts than random positions, p = 0.0072.
What it means has narrowed. The blue plate excludes *persistent* sources —
already excluded by the vetoes — and separates neither hot tiles from the
rest of the catalogue nor defects from sub-exposure events. The direct
evidence that some hot tiles are defects is the **visual** material above:
XE491's moiré band, XE366's printed text, XE622's seam — pixel structure
that is not sky, with the survivors sitting in it. The one "exception" is a
broken astrometric refit, not a dense field. What distinguishes hot tiles
remains their concentration, not their colour behaviour, and the nature of
the red-only population — for hot tiles, for the rest of this catalogue and
for the published lists alike — is not decided by any single-plate
photometric test.

## Caveats

- **Sample scale**: 36 tiles checked against 244 affected plates
  survey-wide. The direction and rough magnitude are now backed by a real
  significance test, but the exact ratio should not be over-generalised
  without a larger draw.
- POSS-I O is intrinsically shallower than POSS-I E — every comparison here
  is survivor-vs-null on the *same* image, never an absolute blue-z
  threshold read on its own.
- The null itself tracks field density (its own z>3 rate varies roughly
  10–44% across tiles), which is why survivors are always compared against
  their own tile's null rather than a fixed threshold.
- **An unresolved confound**: in a sufficiently crowded field, the
  local-background noise estimate itself rises, which can suppress
  z-scores for genuine sources too (XE050 above). This isn't yet separated
  from the emulsion-defect signal — a de-blended or PSF-fit measurement
  would be needed to fully disentangle the two in dense tiles. Both
  readings point the same direction (survivors under-represented on the
  independent exposure); they differ on how much of that is defect versus
  method sensitivity.

## Files

- `tools/blue_plate_hot_tile_check.py` — the check.
- `vasco/downloader.py` — `fetch_skyview_dss`'s `dss1-blue` support (a real
  bug was found and fixed while building this: the STScI request parameter
  that selects plate colour was hardcoded to red regardless of which colour
  was requested; unexercised by any prior caller, so no earlier finding was
  affected).
