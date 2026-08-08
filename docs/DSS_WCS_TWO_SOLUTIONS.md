# DSS headers carry two astrometric solutions, and they disagree

**Summary.** STScI's Digitized Sky Survey headers contain two independent
astrometric solutions for the same plate. `astropy`/`wcslib` uses one; SExtractor
uses the other. On 207 of 633 POSS-I plates they disagree by ~2.34 arcsec. The
explicit FITS-standard solution is the correct one.

This is a property of the archive data and of how common tools read it, not of
this pipeline. Anyone slicing DSS full-plate scans will hit it, and it is not
described in the DSS literature we could find.

## The two solutions

A DSS cutout header carries both of:

1. **The GSSS plate polynomial** — `AMDX1-20`, `AMDY1-20`, `PPO1-6`, `CNPIX1/2`,
   `XPIXELSZ`, `YPIXELSZ`, `PLTSCALE`, `PLTRAH/M/S`, `PLTDECSN/D/M/S`. This is the
   plate model described by Russell et al. (1990) and evaluated by the classic
   `dsspos` routine in wcstools and IRAF.

2. **An explicit FITS-standard solution** — `CTYPE1/2`, `CRPIX1/2`, `CRVAL1/2`,
   `CD1_1 … CD2_2`, alongside `WCSNAME = 'DSS'`.

**`astropy.wcs` prefers (1) even when (2) is present** — wcslib recognises the DSS
keyword set and evaluates the polynomial, silently overriding the standard
keywords. **SExtractor reads only (2)**, because it implements the FITS WCS
standard and knows nothing about GSSS.

Both tools are self-consistent. On a header where the two solutions agree,
nothing is wrong. On a header where they disagree, the two tools return
positions ~2.3 arcsec apart for the same pixel, with no warning from either.

## Measurement

Across 633 POSS-I plates, comparing the two solutions on the same header:

| group | plates | median disagreement |
|---|---:|---:|
| divergent | 207 (32.7%) | **2.344"** |
| consistent | 426 (67.3%) | 0.125" |

Discrimination is essentially clean: 99.0% of the divergent group exceeds 1.5",
against 0.2% of the consistent group.

**The explicit solution is correct.** Against Gaia (G<18, proper motions
propagated to the plate epoch), detections positioned by the explicit solution
land 0.06" away; the same detections positioned by the GSSS polynomial land
2.1" away.

## Why it matters for full-plate scans

IRSA's bulk full-plate scans carry **only** solution (1) — they have no
`CRPIX`/`CRVAL`/`CD` at all. Cutouts from the STScI service carry both.

So a tile cut from a full-plate scan inherits the GSSS answer, while the
equivalent archive cutout is read by SExtractor via the explicit answer. The two
disagree by 2.3" on a third of plates even though **the pixel data is identical
and aligned** — we verified this by brute force, recovering an exact match at
zero offset with the `CNPIX` bookkeeping intact.

## What this is not

- **Not** a scanner-epoch artifact. `XPIXELSZ` takes two values (25.2845 and
  24.9956) marking pre/post-1985 scanning, per Lasker et al. (1990) Table IV note
  d, but 628 of 633 plates share one value and are divergent at 32.8% — a clean
  negative.
- **Not** predictable from any full-plate header keyword. A scan of all 138
  keywords across all 633 plates found nothing above `PLTDECD` at 79% accuracy, a
  declination correlate rather than a cause. The discriminator is not *in* the
  full-plate header; it is the *difference between two solutions*, and only the
  cutout headers carry both.
- **Not** the GSC 1.1 vs 1.2 issue. Morrison et al. (2001) note that DSS header
  astrometry is consistent with GSC 1.1 rather than 1.2, with radial and
  magnitude-dependent residuals growing beyond 2.7 degrees from plate centre.
  That is a real and separate defect — it is a smooth radial term, near zero at
  plate centre, whereas this one is a uniform whole-plate offset measurable *at*
  plate centre.

## The correction

`tools/build_plate_crpix_table.py` measures the per-plate offset between the two
solutions using public headers only, and `tools/slice_plate_tiles.py` applies it
via `--crpix-table`. It is header arithmetic: no catalogue, no fitting, no tuned
threshold.

The offset is near-constant across a plate but not perfectly — about 0.1 px of
real scatter remains, so expect ~0.17" residual rather than zero.

## Effect

Correcting it moves the tight-radius recall of the locally-sliced arm sharply
and leaves radii of 3" and beyond untouched -- exactly what a ~2.4" systematic
predicts, since it was always *inside* a 5" match radius. That is also why it
went unnoticed for so long: no 5" pass/fail test could see it. It was found by
looking at the *shape* of a match-distance distribution, where a population
jumping between the 2" and 3" columns is unmistakable.

Quantified figures against the public reference catalogue are reported in
`RESULTS.md`, alongside the run that produced them.

> Numbers are deliberately not repeated here. Any recall figure depends on the
> reference catalogue it was measured against, and a percentage carries no
> marker of its provenance once copied into prose -- so they live in one place,
> next to the command that generates them.

## References

- Lasker, B. M., et al. 1990, AJ, 99, 2019 — GSSS plate model and scanning
- Russell, J. L., et al. 1990, AJ, 99, 2059 — GSSS astrometric model
- Morrison, J. E., et al. 2001, AJ, 121, 1752 — GSC 1.2 astrometry; DSS headers
  remain GSC 1.1-based
