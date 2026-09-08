# Two archives of the same plates — does the scan source change the catalogue?

DSS1 plate scans are available from more than one archive. This pipeline
pulls full-plate FITS from IRSA, addressed by plate name (see
[the archive-cutout ceiling](../README.md#the-archive-cutout-ceiling)). A
second public archive — STScI's own `FitsArchive/XEsurvey` folder — hosts
scans of the same physical plates, packaged independently. Both are public,
both trace to the same photographic originals, and a third party reproducing
this project's method might reasonably pull from either.

This page measures whether that choice matters.

**The short answer**: it depends entirely on where in the plate a tile sits.
Near a plate's centre, the two archives are indistinguishable — same stars,
same background, same final candidate list. Near a plate's physical edge,
the STScI copy carries a **real, structural pixel defect** in a band close
to the array boundary — not decompression noise, not a rounding artefact — a
sharp region where its pixel values run strongly negative where IRSA's run
strongly positive. On the 6 tiles checked here (2 plates, both radial
zones, full pipeline including the veto and MNRAS filter chain), that
defect was present on **3 of 3 rim tiles checked** and **absent on 2 of 2
centre tiles**, with no exceptions. Its effect on the final candidate list
per tile ranges from none to a 46% count difference, depending on whether
any given tile's real sources happen to sit near the defective region.

**This project is not affected** — every result in this repository uses one
consistently-sourced set of IRSA scans throughout. This page exists for
anyone comparing catalogues built from different scan sources, or
attempting an independent reproduction that pulls plates from STScI
instead.

Everything here uses public inputs and the pipeline's own unmodified code —
`tools/slice_plate_tiles.py`, `vasco/cli_pipeline.py` `step2`–`step4`, run
directly against both archives' pixels with identical parameters (the
`paper-parity` branch's veto/spike-mask configuration, WCSFIX on). No code
was changed to produce anything below.

## What is compared

| | IRSA | STScI XEsurvey |
|---|---|---|
| source | `irsa.ipac.caltech.edu/data/DSS/images/dss1red` | `archive.stsci.edu/missions/dss/FitsArchive/XEsurvey` |
| addressing | plate name | plate name |
| GSSS astrometric keywords | identical on ~78% of plates; **different on ~22%** — see below | |
| provenance match to IRSA plate IDs | — | 926/932 (99.4%) |
| pipeline / parameters | identical both sides — same `.sex` configs (byte-diffed), same veto set (Gaia+PS1), same spike catalogue (USNO-B), same WCSFIX correction, same per-plate epoch | |

> **Correction (2026-09-08).** An earlier version of this page stated that
> the two archives' GSSS astrometric keywords are byte-identical, and that
> astrometric identity was "established and solid". **That is true only for
> about 78% of plates.** It was verified on two plates, both of which turned
> out to belong to the agreeing majority. A 250-plate survey since then
> found the archives carry **different plate solutions on ~22% of plates**.
> The pixel findings below are unaffected — they were measured on plates
> where the two solutions agree exactly — but the astrometric claim was too
> strong and is corrected here rather than quietly edited away.

## A second, independent difference: the plate solution itself

Sampled 250 plates at random from the 932 present in both archives and
compared the two solutions at nine points spanning each array. Pixel
geometry was controlled first, because on a small number of plates the two
archives hold *physically different scans* under one plate name (different
`PLATEID`, array size and pixel scale — e.g. IRSA serving a finer 23040²
rescan where STScI serves 14000²); comparing those at matched fractional
pixel positions measures the geometry difference, not the astrometry, and
they are excluded.

Among the 248 geometry-matched plates:

| | plates | share |
|---|---:|---:|
| identical plate solution | 193 | 77.8% |
| **different plate solution** | **55** | **22.2%** |

For the divergent 55: median offset **2.39″**, and on **5 plates the
disagreement exceeds 5″** somewhere on the array. The `AMDX`/`AMDY`
polynomial coefficients differ in 55 of 55, so it is the plate solution
proper, not a derived keyword.

**Which archive changed:** STScI. Using STScI's own 1994 CD-ROM headers
(`DSS/cdheaders/dss1/<plate>.hhh`) as an independent third reference, IRSA
sits the same distance from the 1994 solution on divergent and agreeing
plates alike (2.58″ vs 2.54″), while STScI roughly doubles its distance on
exactly the divergent ones (2.53″ → 4.80″). Two cutout services agree with
IRSA: STScI's own `dss_search` returns coefficients byte-identical to
IRSA's on all four plates spot-checked, and ESO's DSS mirror lands within
0.15″ of IRSA on the divergent plates while sitting 2.5″ from the STScI
full-plate files. **The divergence is specific to the full-plate
`FitsArchive` FITS, not to STScI as an archive.**

**This does not establish which solution is more accurate.** The historical
GSSS solution is itself known to sit ~2.3″ from modern reference frames on
many plates (see [the astrometric
finding](../README.md#the-astrometric-finding)), so a re-solve could well be
an improvement. What is established is that the two differ, and where the
difference was introduced. Anyone cross-matching a catalogue built from one
archive against one built from the other should expect ~2.4″ systematic
disagreement on about a fifth of plates — comparable to a typical match
radius.

What follows is about the *pixel values themselves*, and was measured on
plates whose solutions agree exactly, so it is independent of the above.

## The mechanism: a hard-edged clip near each plate's physical boundary

Traced one anomalous candidate — bright, real, `FLAGS=0` on the IRSA side —
that survives on IRSA and is dropped on STScI with a corrupted measurement
(flux 4.1x too high, signal-to-noise flipped negative). The star's own
point-spread core is pixel-identical between archives; the corruption comes
entirely from the *background estimation window* SExtractor uses around it,
which straddles a sharp rectangular region where STScI's pixel values
invert sign. Mapped through the tile's WCS into the full plate's own pixel
frame, that region sits only **~6–9 arcmin from the plate array's own
edge** — a thin clipped border, not a smooth gradient across the plate.

![Rim vs. clean tile, both archives, independently z-scaled](figures/scan_source_sensitivity/rim_vs_clean_scan_source_comparison.png)

*Top row*: a tile ~4° from its plate's centre (near the physical edge).
The two archives are visibly different at a glance — STScI's whole
background tone is shifted, and the bright vignette strip that reads white
in IRSA reads dark in STScI. *Bottom row*: a tile ~1.9° from centre. The
two panels are visually indistinguishable, down to a faint ring-shaped
plate defect present identically in both (a real, shared feature of the
photographic plate, not an artefact of either archive).

![One anomalous source, traced to its background window straddling the clip boundary](figures/scan_source_sensitivity/flux_auto_anomaly_edge_clip.png)

## Results across 6 tiles, 2 plates, both zones

Full pipeline run both ways per tile (Gaia+PS1 veto, USNO-B spike mask,
WCSFIX on, MNRAS filter chain) — not just raw detection counts:

| tile | plate | distance from plate centre | IRSA survivors | STScI survivors | matched (≤2″) | STScI negative-pixel fraction |
|---|---|---:|---:|---:|---:|---:|
| `RA313.291_DECp86.291` | XE002 | 3.97° | 24 | 21 | 75% | 24.2% |
| `RA15.438_DECp81.733` | XE002 | 3.37° | 10 | 10 | **100%** | 16.5% |
| `RA21.998_DECp84.430` | XE002 | 1.87° | 11 | 11 | 90.9% | **0.0%** |
| `RA268.426_DECp50.755` | XE181 | 3.97° | 22 | 14 | 54.5% | 34.0% |
| `RA272.877_DECp48.983` | XE181 | 1.87° | 7 | 4 | 57.1% | **0.0%** |

**The pixel defect itself is deterministic**: every rim tile (>2.8° from
plate centre) shows it, every centre tile (<2.1°) does not — confirmed on
two plates. **Its effect on the final candidate list is not** — one rim
tile with 16.5% corrupted pixels still matched 100%, because none of its
real sources happened to sit near the defective region; another matched
only 54.5%. It depends on where individual sources fall, not just on
whether the defect is present in the tile.

The two centre-tile match rates (90.9%, 57.1%) look large in percentage
terms but come from only 4–11 survivors per tile — ordinary near-threshold
photometric jitter between two independently-processed pixel copies (both
archives round-trip the same photographic original through different
digitisation/compression pipelines), not the structural defect above:
directly confirmed by checking pixel data at both centre tiles — **0.0%
negative pixels, no exceptions**.

## Does the defect manufacture candidates that shouldn't pass?

A natural follow-up: could the corrupted background inflate a source's
measured significance enough to push it *over* the survey's SNR gate when
it would otherwise fail? The mechanism can clearly move flux/SNR
substantially — up to a full sign flip in the case above, a ~40% inflation
on a second, independently-checked real star. Tracing every candidate that
passed on STScI without an IRSA counterpart (across all 6 tiles, matched
against IRSA's raw catalogue at a wide radius, not just the tight tolerance
used for the headline numbers) found the mechanism doing exactly this kind
of inflation on one more real star — but in every case checked, the
inflated source was already above the gate on the clean side too, just by
a smaller margin. One genuine gate-flip (fail on IRSA, pass on STScI) was
found, but its background window was independently confirmed to be
completely uncorrupted — ordinary cross-archive photometric noise, not
this mechanism. **Mechanistically capable of manufacturing a false pass;
not caught doing so in the tiles checked so far.**

## A third, independent archive: same answer

If the defect were a property of the plate itself — something any faithful
digitization would reproduce — a third archive should show it too. Fetched
the same sky position from **ESO's own DSS mirror**
(`archive.eso.org/dss/dss/image`), a fully independent archive.

![Three independent archives, same plate](figures/scan_source_sensitivity/eso_third_archive_check.png)

Confirmed same plate first, then compared: **ESO shows 0.0% negative
pixels — matching IRSA, not STScI.** Two independent archives agree with
each other and disagree with STScI specifically, on the one plate where a
clean three-way comparison was possible.

A second attempt, on the other plate checked in this note, ran directly
into the archive-cutout ceiling described earlier in this README: querying
ESO's position-addressed service at that tile's exact coordinates returned
a **different plate entirely** than the one IRSA and STScI both serve for
that position — an unplanned, live example of exactly that problem, not a
hypothetical one.

## Caveats

- **6 tiles, 2 plates.** Both POSS-I red, both near/mid-northern
  declination (XE002 near-polar, XE181 +48°). No claim about the exact
  distribution of impact across a full survey, and no check yet at low
  declination or on a different plate generation.
- **Centre-tile match percentages are noisy at these small sample sizes**
  (4–11 rows) — the qualitative finding (no pixel defect, smaller absolute
  impact than rim tiles) is solid; the specific percentages are not
  precise estimates of a survey-wide rate.
- **The exact margin from the plate edge (~6–9′) and the transition's
  sharpness were characterised on one tile in detail** and corroborated on
  five more via presence/absence of the pixel signature; not yet mapped
  at sub-tile resolution on every plate edge.
- **What causes the STScI-side pixel inversion at the edge is not
  identified** — only that it is present, structural, and archive-specific.
  It is a pixel-processing difference, not a geometry one: on the plates
  where it was measured, both archives' GSSS keywords are identical, so the
  coordinate system is not involved. (On ~22% of plates *overall* the
  solutions do differ — see the section above — but those are a separate,
  independent finding and not the plates used for the pixel work here.)
- **The pixel defect is present in STScI's 2005 raw scan, not introduced by
  the 2014 FITS packaging.** The `ScanArchive/XEsurvey/*.pim` files are raw
  headerless 16-bit arrays, and two sampled rows (one in the border band,
  one at plate centre) reproduce the 2014 FITS values exactly — the FITS
  conversion is a byte-order repackaging. Two rows of one plate, so a strong
  indication rather than a whole-array proof.
- **A position-addressed cutout service does not appear able to reach the
  defect.** STScI's `getimage` defaults to a `FURTHEST_FROM_EDGE` plate
  choice — among plates covering a position it serves whichever keeps the
  extraction furthest from any array boundary — and each plate's edge is a
  neighbour's interior. Across 16 probes on 4 edges of 3 plates, the closest
  any served extraction came to its own plate's array edge was 560 px,
  against a defect band at ≤330 px. So the defect is reachable mainly by
  slicing full-plate scans directly. Probes on 3 plates only; survey-boundary
  plates, where a plate has no neighbour on one side, are untested.
- **The 22% solution-divergence rate is a proportion from a 250-plate random
  sample** (≈ ±3 points of sampling error) — read it as "about a fifth of
  plates", not as 22.2% exactly.
- **Which plate solution is more accurate is not established.** The
  comparison above is header-level and provenance-level only. Deciding it
  requires cross-matching real detections against a modern reference frame
  under each solution.
- **The candidate-manufacturing question is checked, not settled**: no
  confirmed case of the defect flipping a genuine fail to a pass, only
  that it moves measured SNR by enough (40-100%+) that it plausibly could.
- **The third-archive cross-check is one plate.** The second attempt hit a
  different-plate mismatch from ESO's own service before a same-plate
  comparison was possible.

## Reproducing

```
tools/slice_plate_tiles.py --plate-fits <plate>.fits --tiles-dir <out> \
    --crpix-table data/plate_crpix_table.csv
python -m vasco.cli_pipeline step2-pass1 --workdir <tile>
python -m vasco.cli_pipeline step3-psf-and-pass2 --workdir <tile>
VASCO_DISABLE_USNOB=1 VASCO_SPIKE_CATALOG=usnob VASCO_PLATE_EPOCH_YEAR=<epoch> \
    python -m vasco.cli_pipeline step4-xmatch --workdir <tile>
```

Same entry points as any normal run — nothing here required pipeline
changes. STScI XEsurvey plates: `https://archive.stsci.edu/missions/dss/FitsArchive/XEsurvey/`.
