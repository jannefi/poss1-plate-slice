# The reference catalogue was published — what that changes here

**2026-10-15.**

## What happened

On **2026-10-14** the code repository accompanying Bruehl et al. (2026) was
made public at `github.com/dca-doherty/VASCO-ML`. It contains the
full **107,875-row VASCO 2025 sample** (Bruehl et al. 2026) — the catalogue
this pipeline's output has been compared against throughout — as
`data/SUPERVIKTIG_HELAVASCO_validated_v4.csv`, described in that repository as
"107,875 transient candidates from the Solano et al. detection pipeline" (Bruehl et al. 2026).

The paper is published: Bruehl, Doherty, Streblyanska & Villarroel (2026),
*Use of machine learning to enhance detection of transient astronomical
phenomena in historical observatory images*, *Scientific Reports* **16**,
31890, [doi:10.1038/s41598-026-66148-2](https://doi.org/10.1038/s41598-026-66148-2),
published 10 August 2026 and previously circulated as arXiv:2604.18799.

Until that date the catalogue was available only on request. Every statement in
this repository that describes it as unavailable should be read with this date
attached.

## What the released file contains

Twenty columns: `ra`, `dec`, `obs_date`, then catalogue-level measurements
(`psf_ratio`, `elongation`, `compactness`, `sharpness`, `snr`,
`n_comparison_stars`, `candidate_score`), classification flags, and three
machine-learning columns.

Two properties are worth knowing before anyone plans work against it, because
both are easy to assume otherwise:

- **There is no plate identifier and no time of day** *in the GitHub copy*.
  The epoch is given to the day only, and on the published POSS-I observing
  record a single date commonly covers several plates, so the exposure a given
  row came from is not determined by those columns alone. **Corrected
  2026-10-19:** the same rows are also published through the Spanish Virtual
  Observatory at
  [`vanish-bruehl26`](https://svocats.cab.inta-csic.es/vanish-bruehl26/), and
  that version does carry both — `name` (the plate) and `ut_observation` (a
  full timestamp).
- **The `ml_probability` column is not the per-candidate score behind the
  paper's figures** *in the GitHub copy*. (**Corrected 2026-10-19:** the SVO
  version carries the real scored vector as `prob` — see
  [`REFERENCE_CATALOGUE_REPRODUCIBILITY.md`](REFERENCE_CATALOGUE_REPRODUCIBILITY.md) §4.) It takes 287 distinct values across the
  107,875 rows (doi:10.1038/s41598-026-66148-2), and 99.70% of them carry one
  identical value. The repository's own committed
  decile tables span a continuous range in ten bins of roughly equal size, so
  the scored vector those tables were built from is a different, finer-grained
  one; it is not among the committed files. The 250 training labels and the
  image-derived features are likewise listed as not included.

These are observations about what a reader will find in a public file. **They
are not a claim about the paper's results**, which this repository takes no
position on, consistent with *What this pipeline claims* in the README.

## What it changes for this repository

**Nothing about the method, and nothing about any published number here.**

What changes is what a reader can now check. This pipeline's claim has always
been that it needs nothing but public inputs — the published plate list, IRSA's
full-plate scans, and public catalogues. While the reference catalogue was
unavailable, that claim was nearly self-evident: it could not have been used,
because it could not be obtained. Now that anyone can download it, the claim
needs demonstrating rather than assuming, and two things do that:

1. **`tools/audit_independence.py`** scans both the working tree and the entire
   git history for any trace of the reference catalogue, and exits non-zero on
   a hit. It was widened on 2026-10-15 after an audit found five prose
   references that its earlier patterns had not anticipated; all five were
   comments, none was data. Those already in the published history, which
   editing cannot reach, are recorded in
   `tools/audit_independence_baseline.json` with a reason for each, so the gate
   fails on anything new rather than failing permanently. See
   `docs/DATASET_NAMING.md`.
2. **The dates are checkable and they are not close.** Every commit that
   produced a result here — every pipeline change, every stage, every released
   catalogue — falls between **2026-08-08** and **2026-09-11**. The reference
   catalogue became available on **2026-10-14**, at least a month after the
   last of them. Nothing published here could have been tuned against a file
   that did not yet exist publicly, and `git log` is the evidence.

   Commits dated after 2026-10-14 are documentation only, this note among them;
   `git log --stat` separates the two, and any future commit that touches the
   pipeline will fall on the other side of that line and cannot rely on this
   argument.

That second point is the stronger of the two, and it has a shelf life: it holds
for the existing history and stops holding the moment new work is done with the
file in hand.

## What this repository will do with it

- **Cite it.** The paper is cited in the README's reference list, with its DOI.
- **Not redistribute it.** The repository's `LICENSE` grants MIT for code; its
  README states separately that the research datasets "are made available for
  reproducibility of the associated publication". **Corrected 2026-10-19:** a
  third repository published alongside the challenge states instead "Data: CC
  BY 4.0", which would permit redistribution with attribution. The two
  statements conflict and we do not rely on either: this repository ships no
  copy of the file and no derivative containing its rows, and anyone intending
  to redistribute should settle the terms with the authors. The SVO archive
  asks separately to be acknowledged.
- **Declare any future use.** If any result here is ever derived with that file
  in the loop, it will say so on the result, not in a footnote. The point of
  the audit above is to make the absence checkable; the point of this paragraph
  is that the absence is a choice that could change, and would be disclosed if
  it did.

## What is not claimed here

This note records a fact about availability and its consequences for
reproducibility. It makes no claim about whether the sources in that catalogue
are real, which remains contested and is surveyed in the README's reference
list; and it makes no claim about the correctness of the analysis built on it.
