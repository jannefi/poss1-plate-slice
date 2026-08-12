# Results

All figures here are measured against the **published** POSS-I vanishing-source
catalogue (SVO `vanish-possi`, 5,399 rows). Nothing on this page depends on data
that is not publicly available, so every number can be reproduced independently —
see [`docs/REPRODUCING.md`](docs/REPRODUCING.md).

## Footprint

Every POSS-I red plate whose centre declination is ≥ −3.0°: **642** of the 932
scans in IRSA's `dss1red` library. Three independent public checks support that
boundary, and none of them requires a catalogue this project cannot redistribute.

**The threshold is not a tuned parameter.** No plate centre lies between −5° and
−1°, so every cut inside that 4.4° gap selects the same 642 plates.
`tools/build_plate_manifest.py` asserts the gap is empty rather than assuming it,
so a future library that broke the separation would fail loudly instead of
quietly returning a different footprint.

**The southern edge is where the published catalogue ends.** The Solano et al.
(2022) vanishing-object catalogue (SVO `vanish-possi`, 5,399 rows) spans
−3.32° ≤ δ ≤ +87.97°, with 76 rows below δ = 0 and exactly one below −3.

**The plates are populated across essentially the whole set.** The publicly
released NeoWISE-proximate subset (SVO `vanish-neowise`, 171,753 rows — set *W*
of Watters et al. 2026, characterised there as a spatially uniform-random sample
of the parent catalogue) has sources on **633 of the 642**.

For reference, VASCO's published analyses used 635 plates (Villarroel et al.
2026, Fig. 1). This footprint is 642. Which plates the two lists share is not
something this repository can establish, because the 635-plate list has not been
published.

## Raw detection recall

Fraction of published catalogue rows with a detection within a given radius,
over the full 642-plate footprint. 7×7 tiles per plate, per-plate CRPIX
correction applied.

| arm | 1″ | 2″ | 3″ | 5″ | 10″ | 30″ |
|---|---:|---:|---:|---:|---:|---:|
| `fullplate` — local slices of IRSA full-plate scans | 83.39% | 94.26% | 96.91% | **97.24%** | 97.85% | 99.04% |
| `archive` — STScI cutout service | 67.61% | 81.66% | 84.78% | 85.42% | 86.59% | 90.68% |
| **`fullplate+archive`** | 86.22% | 96.06% | 98.30% | **98.59%** | 99.00% | 99.56% |

Arm sizes: `fullplate` 189,241,189 detections over 642 files; `archive`
179,887,397 over 31,004 tiles.

Extending the footprint from the 634 plates of the first campaign to all 642
moved these figures by at most 0.02 points — `fullplate` at 1″ rose from 83.37%
to 83.39%, the union at 2″ fell from 96.07% to 96.06%, and every other cell is
unchanged. That is the expected result rather than a disappointing one: the
reference catalogue was built over its authors' own footprint, so the plates
added here contribute coverage where it has almost no rows to recover. The
recall numbers and the footprint claim are therefore close to independent —
widening the footprint does not flatter the recall.

**The locally-sliced arm alone exceeds the archive arm by ~12 points at 5″.**
That is the central result: naming the plate yourself outperforms letting a
cutout service choose one for you. The two arms miss different sources, so their
union beats either.

Reproduce with:

```bash
python3 tools/union_parity_fullscale.py \
    --ref-csv <vanish_possi.csv> \
    --arm fullplate=<slice_radec> --arm archive=<archive_radec> \
    --combine fullplate+archive --out-dir work/union_R
```

### Provenance of the archive arm

The `archive` arm's detections were produced by this project's predecessor
pipeline at identical settings, and were **not** regenerated for this release —
doing so requires re-downloading 31,004 cutouts from the STScI service. The
`fullplate` arm was generated end to end by the code in this repository.

## Astrometric correction

| | plates | median offset |
|---|---:|---:|
| needing correction | 209 of 642 (32.6%) | 2.29″ → ~0.09″ |
| already correct | 433 | 0.02″ |

Derived from public headers alone, with no catalogue and no fitted parameter —
see [`docs/DSS_WCS_TWO_SOLUTIONS.md`](docs/DSS_WCS_TWO_SOLUTIONS.md). Its effect
on recall is confined to the tight radii: it is worth roughly 29 points at 1–2″
and nothing at 3″ or beyond, because a ~2.3″ systematic sits *inside* a 5″ match
radius and is invisible to a 5″ pass/fail test.

Two plates (XE285, XE311) exceed the 0.2 px scatter threshold and are corrected
with a warning rather than silently.

## Candidate catalogue

**Not yet released.** The veto chain and post-processing run is in progress; this
section will carry the catalogue, the exact stage chain that produced it, and the
per-stage row counts.

When it lands, note in advance that three implemented stages — SkyBoT,
SuperCOSMOS and VSX — were **not executed**, for the reasons and with the
measured impact given in [`docs/PARAMETERS.md`](docs/PARAMETERS.md).

## What these numbers do and do not say

Recall measures whether this pipeline *finds* the published sources. It says
nothing about what those sources are. See *What this pipeline claims* in the
README: a detection surviving the full chain is a source that is on the plate and
not in modern catalogues, which is a statement about two catalogues rather than
about the universe.
