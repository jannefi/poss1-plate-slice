# Results

All figures here are measured against the **published** POSS-I vanishing-source
catalogue (SVO `vanish-possi`, 5,399 rows). Nothing on this page depends on data
that is not publicly available, so every number can be reproduced independently —
see [`docs/REPRODUCING.md`](docs/REPRODUCING.md).

## Raw detection recall

Fraction of published catalogue rows with a detection within a given radius.
634 plates, 7×7 tiles per plate, per-plate CRPIX correction applied.

| arm | 1″ | 2″ | 3″ | 5″ | 10″ |
|---|---:|---:|---:|---:|---:|
| `fullplate` — local slices of IRSA full-plate scans | 83.37% | 94.26% | 96.91% | **97.24%** | 97.85% |
| `archive` — STScI cutout service | 67.61% | 81.66% | 84.77% | 85.42% | 86.59% |
| **`fullplate+archive`** | 86.22% | 96.07% | 98.30% | **98.59%** | 99.00% |

Arm sizes: `fullplate` 186,480,066 detections over 634 files; `archive`
179,887,397 over 31,004 tiles.

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
| needing correction | 207 of 634 (32.6%) | 2.29″ → ~0.09″ |
| already correct | 427 | 0.02″ |

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
