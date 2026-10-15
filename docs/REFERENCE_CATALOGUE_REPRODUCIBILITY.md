# Four things the published release lets anyone check

**2026-10-15.** On 2026-10-14 the 107,875-row VASCO 2025 sample (Bruehl et al. 2026)
became public, in the code repository accompanying *Use of machine
learning to enhance detection of transient astronomical phenomena in historical
observatory images*, Sci Rep **16**, 31890,
[doi:10.1038/s41598-026-66148-2](https://doi.org/10.1038/s41598-026-66148-2).
Background on the release itself is in
[`REFERENCE_CATALOGUE_RELEASE_20261015.md`](REFERENCE_CATALOGUE_RELEASE_20261015.md).

This note records four things that can now be verified by anyone, because the
inputs to all four are public. It is a reproducibility note and nothing more:
**it takes no position on whether the sources in that catalogue are real**,
which is contested and surveyed in the README's reference list, and it makes no
claim about the correctness of any analysis built on them.

**Declaration.** This repository now has access to that file. Nothing in this
repository was derived using it — every commit that produced a result here
predates its release by at least a month, and `tools/audit_independence.py`
checks that claim over the tree and the whole history. The measurements cited
below were made before the release, from public inputs, and are unchanged by
it. Any future result that uses the file will say so on the result.

---

## 1. The catalogue cannot be derived from the method the 2022 paper describes

MNRAS 515(1):1380 describes tessellating the sky in circular regions of 30
arcmin radius, fed by an image cutout service. **No tessellation of that kind
can produce this catalogue**, and the reason is measurable without any access
to it.

Which plate an archive cutout service returns for a position is deterministic
in the position, and it changes on **arcminute** scales rather than degree
scales. We measured the curve directly — 2,025 probes across 15 offset bins —
and it does not decay to zero at small offset: the probability that the plate
chosen for a tile centre differs from the plate chosen for a source inside that
tile rises within a few arcminutes and then plateaus at roughly 17%.

That sets a ceiling on how much of the catalogue any position-blind grid can
recover, whatever its pitch:

| tile pitch | predicted raw recall | cost relative to 60′ |
|---|---:|---:|
| 60′ | 85.0% | 1× |
| 12′ | 90.1% | 25× |
| 6′ | 93.8% | 100× |
| 2′ | 97.0% | 800× |

A hundredfold increase in downloads buys 8.8 points. Adaptive subdivision fails
for the same reason — it is also position-blind. The model was validated three
ways, including against two held-out controls, and it predicts our own observed
raw recall (85.02% predicted against 84.28% measured).

**This is a statement about the published method, not about anyone's
competence.** Two independent facts support it. Hambly & Blair (2024)
independently confirmed that this repository's implementation of the 2022
candidate-selection filters matches the paper, so the gap is not a filter bug.
And Hayes (2026, arXiv:2604.04810), a separate automated pipeline built
independently of this one, reports a 63.9% catalogue-level match — a different
route reaching the same wall.

The honest summary is that **the published text under-determines the result**.
That is a normal and common finding about a method description, and it is worth
saying plainly because several parties have now spent real effort on it.

## 2. What the coordinates do, and do not, show about how it was built

The ceiling above implies a negative that is worth stating precisely: **the
catalogue was not produced by a position-blind tessellation of the sky**, grid
or circular, because no such tessellation reaches it.

**What we could not show.** We looked for a positive signature of tiling in the
released coordinates themselves and did not find one, and the reason is
instructive rather than suggestive: our own tessellation is a 7×7 grid in each
plate's *pixel* frame, so even a catalogue that was built on a grid carries no
periodicity in RA/Dec to detect. A phase test therefore has no working positive
control and cannot support a conclusion either way. We also checked for the
coverage gaps a missing-tile pattern would leave, at 5° and at 1° resolution,
and found none: both catalogues occupy every 5° RA bin, and the cells where
this repository's catalogue has rows and the released one does not are spread
across all 360° of RA rather than forming stripes.

So the claim in this section rests on the recall ceiling, not on a signature in
the coordinates. **Anyone repeating this should not expect the coordinates to
show it directly.**

**What the released file cannot settle at all.** It carries `ra`, `dec` and
`obs_date` and no plate identifier or exposure time. A single observing date
commonly covers several plates, so for most rows the exposure a row came from
is not determined by the released columns. Any account of the construction that
depends on knowing the plate per row is therefore not checkable against the
public file, and we make none here.

## 3. Why this repository's catalogue is not the same list

Two separate reasons, and conflating them is the main way this gets misreported.

**Stage the number before quoting it.** Agreement between two catalogues of
this kind is radically different at different stages of the pipeline, and a
figure quoted without its stage is close to meaningless:

| stage | agreement |
|---|---|
| raw detection — is there a source at that position at all | high, mid-90s percent |
| candidate level — after the published quality gates | roughly 19% |

A previously published figure from this repository, 94–98%, is **raw detection
recall**. It has been quoted externally as a replication rate. It is not one.

**The gates are supposed to remove rows.** Most of the difference is the 2022
paper's own quality criteria doing what the paper says they do — `FLAGS = 0`,
`SNR_WIN > 30`, elongation and `SPREAD_MODEL` thresholds, the spike mask,
deduplication. A catalogue that applies them strictly is *expected* to be
smaller than one that does not. The remaining difference, in both directions, is
implementation: detection thresholds, deblending, and choices the published text
does not specify. Per-stage counts are in
[`FUNNEL_ATTRIBUTION.md`](FUNNEL_ATTRIBUTION.md) and `RESULTS.md`; the
thresholds themselves, with the origin of each, are in
[`PARAMETERS.md`](PARAMETERS.md).

The fair comparison is therefore not total against total. It is like against
like, at a named stage, with the gate behaviour of both sides stated.

## 4. The machine-learning results cannot be independently reproduced

The repository's README states: *"Every paper number can be regenerated from the
raw inputs with a single command."* We tested that on a fresh clone, in the
exact pinned environment its `requirements.txt` specifies, on 2026-10-15.

The documented verification command succeeds:

```
python reproduce.py --verify-only      ->  22 passed, 0 failed, exit 0
```

But it recomputes nothing. All 22 checks read values out of result files
committed to the repository and compare them against expected values also
committed to the repository. It confirms that two files in the same clone agree
with each other.

The actual reproduction does not succeed. Running it as documented:

```
python reproduce.py
  Stage 1/4: Core pipeline (train + score)
  [ERROR] Feature cache not found: results/ml_pipeline_red_allcandidates/feature_cache.parquet
  [ERROR] Run with extraction first, or provide an existing cache.
  [FAIL] Core pipeline (train + score) failed (exit 1, 0.5s)
```

It stops in half a second at the first of four stages, and the orchestrator
invokes the pipeline with `--no-extract`, so the documented command never
attempts to build the cache it requires. Working around that by hand, each
obstacle in turn:

1. **The catalogue is not where the configuration looks for it** — the config
   names `data/raw/`, the file ships at `data/`.
2. **The pinned requirements are incomplete** — the training stage imports
   `seaborn`, which `requirements.txt` does not list.
3. **Feature extraction reports success without the images.** With the path
   fixed it completes in under a second and writes a cache of all
   107,875 candidates (doi:10.1038/s41598-026-66148-2) in 29 columns — the released file's own 20 columns plus nine derived from them.
   **None of the 16 image-derived features** the paper describes are in it:
   no FWHM, aperture flux, ellipticity, connected-pixel count, gradient
   magnitude or radial symmetry. The configuration points at
   `fits_cutouts: D:/plate_data_3deg`, a local path on an author's machine, and
   with no images present the stage produced a catalogue-only cache and
   reported `[OK]`.
4. **Training then stops, and this one cannot be worked around:**
   `ValueError: No label files loaded. Check config.yaml paths.` The 250 expert
   labels are listed in the README as not included, and without them the
   classifier cannot be trained.

The scored catalogue itself is also absent — it is excluded by `.gitignore` and
there are no GitHub releases — and the one scoring column in the released file
is not the one behind the paper's figures: it takes 287 distinct values over
the 107,875 rows (doi:10.1038/s41598-026-66148-2) and 99.70% of rows carry a
single identical value, while the committed decile tables span a continuous
range in ten bins of roughly equal size.

**What this does and does not mean.** It means the published results cannot be
regenerated from the published artifacts, contrary to what the README states,
and that a reader who runs the documented command will get a pass that verifies
nothing. **It is not evidence that the results are wrong.** They may be entirely
correct; the point is that no third party can currently check, which is a
property of the release rather than of the analysis.

## What would resolve each point

- For §1 and §2: the procedure that actually generated the list, in enough
  detail to re-run — in particular how sky positions were chosen, since that is
  the step the published description does not determine.
- For §3: nothing from anyone else. Our side is published, and the comparison
  is now runnable by anyone.
- For §4: the 250 labels, the feature matrix or the images it was built from,
  and the scored catalogue. Any one of them moves the situation; all three
  would close it.

## Reproducing this note

The ceiling measurement in §1 uses only public inputs and is in this
repository. The §4 results come from a clone of the public repository at commit
`406ac08` in the environment its own `requirements.txt` pins. This repository
ships **no copy of the released catalogue and no derivative containing its
rows** — that file is available from its authors, under the terms they set, and
their README asks that the paper be cited by anyone using it.
