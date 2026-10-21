# Which plate, where plates overlap

**2026-10-21.**

POSS-I plates overlap. On overlapping sky a source can be measured on any of
the plates covering it, and which one a catalogue used is a property of how
that catalogue was assembled, not of the sky. It is also measurable from the
finished catalogue, with nothing but the positions, the plate names and the
public plate geometry.

This page reports that measurement for this repository's own catalogue and for
the published VASCO 2025 sample, and sets out what it does and does not
support. Tool: [`tools/plate_choice_test.py`](../tools/plate_choice_test.py).
Everything here is reproducible from public inputs in a few seconds.

## The measurement

For every row: find every plate whose array covers the position, compute the
distance to the nearest array edge on each, and ask whether the row's own plate
is the one that puts it furthest from any edge. Coverage is computed over all
**932** POSS-I red plates in the IRSA archive, not only the 642 this project
processed — a position can be covered by a plate nobody chose to run.

| | rows | covered by >1 plate | own plate is the furthest-from-edge plate |
|---|---:|---:|---:|
| published catalogue of 107,875 candidates (doi:10.1038/s41598-026-66148-2) | 107,875 | 45,985 (42.6%) | **99.9%** |
| this repository's paper-parity S0 | 134,976 | 88,667 (65.7%) | **36.3%** |

The gap is not a matter of degree. Where the published catalogue is not on the
furthest-from-edge plate, the median shortfall is **0.04′**; ours is **20.59′**.
The shortfall distribution is the clearest way to see it:

| shortfall percentile | published | ours |
|---|---:|---:|
| p50 | 0.000′ | 10.138′ |
| p90 | 0.000′ | 38.984′ |
| p99 | 0.00′ | 92.37′ |

The published catalogue does not merely tend towards the most interior plate.
Out to the 99th percentile it is on exactly that plate, with zero shortfall.

## It is not an astrometric artifact

The two DSS archives' plate solutions disagree by up to about 2.4″ on some
plates ([`DSS_WCS_TWO_SOLUTIONS.md`](DSS_WCS_TWO_SOLUTIONS.md)), so the obvious
objection is that the result turns on which solution was used. It does not.
Sweeping the tolerance from well inside that band to well outside it:

| tolerance | 0.3″ | 0.6″ | 1.2″ | 2.4″ | 5″ | 10″ | 30″ |
|---|---:|---:|---:|---:|---:|---:|---:|
| published | 99.9% | 99.9% | 99.9% | 99.9% | 100.0% | 100.0% | 100.0% |
| ours | 36.3% | 36.3% | 36.3% | 36.3% | 36.4% | 36.5% | 36.9% |

Flat across two orders of magnitude. An arcsecond-scale astrometric choice
cannot produce a sixty-point difference.

## The control is inside the same catalogues

On sky covered by only one plate there is no choice to make, so if the
measurement is picking up plate *selection* rather than something else, the two
catalogues should agree there and diverge only where a choice exists. They do:

| edge distance on the row's own plate | min | 1st pct | median | within 15′ |
|---|---:|---:|---:|---:|
| published, single-coverage sky | 17.63′ | 35.32′ | 76.01′ | 0.00% |
| ours, single-coverage sky | 5.02′ | 35.09′ | 76.52′ | 0.02% |
| published, all rows | 8.58′ | 18.44′ | 52.64′ | 0.06% |
| ours, all rows | 0.14′ | 1.62′ | 30.26′ | 28.90% |

Where no choice exists the two are the same catalogue to the percentile. Where
a choice exists, ours lands wherever the grid put it — 28.9% of our rows sit
within 15′ of an edge — and the published one does not.

**That also explains its edge distribution without any edge cut.** The
published catalogue has often been described as having a floor around 18′ from
the plate edge; its 1st percentile here is 18.44′. No cut is needed to produce
that. It is what taking the most interior covering plate gives you, because the
sky rarely offers a most-interior plate that is itself near an edge.

## What this supports, and what it does not

**Supports:** the rows were placed on plates by a rule evaluated *per position*.
STScI's `getimage` picks `FURTHEST_FROM_EDGE` by default
([`SCAN_SOURCE_SENSITIVITY.md`](SCAN_SOURCE_SENSITIVITY.md)), and that is the
rule the geometry shows. It is inconsistent with a position-blind tessellation,
which is what our own 36.3% is: a grid cannot prefer the interior plate,
because it does not know where the source is when it chooses the tile.

**Does not support — and we are not claiming it:** any account of who ran what,
in what order, or with which software. A per-position rule cannot have been
applied at detection, because the positions are not known until after
detection; that much follows from the geometry. Anything beyond it — how the
first pass was done, what was dropped, what was re-measured — does not follow
from this measurement, and this page does not assert it.

**Why it matters anyway.** It bears directly on what an independent
reproduction can reach. A reproducer working from a position-blind tiling is
not doing the same thing, whatever its pitch, and the recall ceiling documented
in [`REFERENCE_CATALOGUE_REPRODUCIBILITY.md`](REFERENCE_CATALOGUE_REPRODUCIBILITY.md)
§1 is the quantitative form of that. This page is the mechanism behind that
ceiling rather than a separate claim.

## Reproducing it

```bash
# the control first: a tessellation should score low
python3 tools/plate_choice_test.py \
    --ref-csv results/s0-642-paper-parity-20260828/stage_S0.csv.gz \
    --manifest results/s0-642-paper-parity-20260828/tile_manifest.csv.gz

# then any catalogue you have locally, with --sweep for the tolerance scan
python3 tools/plate_choice_test.py --ref-csv <your copy> --sweep
```

[`data/plate_footprints.csv.gz`](../data/plate_footprints.csv.gz) carries the
plate-solution keywords of all 932 archive plates, read from the public IRSA
FITS headers. Your astropy rebuilds the same WCS this tool used, so the result
does not rest on our arithmetic — and it is 168 KB instead of a few hundred
gigabytes of plate scans.

This repository ships no copy of any reference catalogue. `--ref-csv` is a
local path; the tool downloads nothing.

## Limits

- It measures a finished catalogue's geometry. It is evidence about a
  selection rule, not a reconstruction of a pipeline.
- Single-coverage rows carry no information about plate choice, which is why
  they are reported separately and used as the control.
- Coverage uses array boundaries from each plate's own solution. Vignetting,
  emulsion edges and furniture are not modelled; "on the plate" here means
  inside the array.
- The 932-plate footprint set is the IRSA archive as it stands. STScI has
  announced corrections to DSS headers, so anyone repeating this later should
  record which version they used.
