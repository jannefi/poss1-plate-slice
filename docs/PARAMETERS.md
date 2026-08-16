# Parameters

Every threshold this pipeline applies, where it comes from, and whether you may
change it.

The distinction that matters throughout: **published values are not free
parameters.** A threshold taken from Solano et al. (2022) reproduces that paper
only at its published value. Changing one is a deviation from the method, to be
recorded as such — not tuning.

| origin | meaning |
|---|---|
| **paper** | From Solano et al. (2022). Changing it means you are no longer reproducing that method. |
| **engineering** | Our choice, made for cost or correctness. Change freely, but read the rationale first. |
| **measured** | Derived from data, with a tool in this repository that regenerates it. |

## Detection

| parameter | value | origin |
|---|---|---|
| `DETECT_THRESH` | 5 | paper — extraction time, `configs/sex_pass*.sex` |
| tile geometry | 60′ × 60′ square | paper |
| tiles per plate | 7 × 7 grid, pixel space | engineering |
| SExtractor passes | two (PSFEx model required for `SPREAD_MODEL`) | paper |

A single-pass mode exists (`VASCO_REPRO_SINGLE_PASS=1`) and is **not** used for
production: without a PSF model there is no `SPREAD_MODEL` column, so the
morphology gate below rejects everything.

## Candidate selection — Solano et al. (2022), Section 2

Implemented in `vasco/mnras/filters_mnras.py`. **All paper values.**

| parameter | value |
|---|---|
| `FLAGS` | `== 0` |
| `SNR_WIN` | `> 30` |
| `ELONGATION` | `< 1.3` |
| `SPREAD_MODEL` | `> -0.002` (absolute threshold) |
| σ-clipping on `FWHM_IMAGE`, `ELONGATION` | optional, 2σ via median + MAD |

## Cross-match veto

| parameter | value | origin |
|---|---|---|
| veto radius, Gaia / PS1 / USNO-B1.0 | 5″ | paper |
| match mode | STILTS `find=best` — one match per candidate row | engineering |
| proper motion | propagated to the plate epoch before matching | engineering |

Proper-motion propagation is not optional in practice. Over a ~65-year baseline a
15 mas/yr star moves 1.0″, a fifth of the veto radius, so matching Gaia at its
own epoch silently fails to veto real moving stars.

## Bright-star / diffraction-spike mask

| parameter | value | origin |
|---|---|---|
| catalogue | Pan-STARRS1, `rMeanPSFMag <= 16` | **deviation** |
| fetch radius | 90′ | engineering |
| association distance | magnitude-dependent rules, `SpikeConfig` | paper |

**Deviation:** Solano et al. use USNO-B1.0 for the spike mask; this pipeline uses
PS1. This is a **deliberate, test-driven engineering choice**, made in the
original `vasco` pipeline and carried forward — see its CHANGELOG entry
[`[0.06.9]`, 2025-11-22](https://github.com/jannefi/vasco/blob/b1cefc05d054d5b235603a80433fa3f16157f186/CHANGELOG.md?plain=1#L23).
It is reversible.

The two have not been compared head to head at survey scale. A local USNO-B
mirror now makes that practical, and it is an open item.

The 90′ **fetch** radius is not an association distance. A 60′ square tile's far
corner is 30·√2 ≈ 42.4′ from centre, so 90′ covers the whole tile plus margin for
stars outside it whose spikes reach in. How close a candidate must lie to a
bright star to be *called* a spike is set by the magnitude-dependent rules.

## Stages that ship but are NOT run

Three stages of the published method are **implemented here and deliberately not
executed** for the current results, plus one stage of our own (last row).

| stage | script | status |
|---|---|---|
| SkyBoT — known solar-system objects | `scripts/stage_skybot_post.py` | not run, **rate + yield bound measured** |
| SuperCOSMOS — plate-artifact discrimination | `scripts/stage_supercosmos_post.py` | not run, **measured: −40.0%** |
| VSX — known variable stars | `scripts/stage_vsx_post.py` | not run |
| Neighbour-plate persistence — second epoch from overlapping plates | `tools/stage_neighbor_persistence.py` | not run, measured |

**Neighbour-plate persistence (ours, not in the paper — measured 2026-08-15,
not applied).** POSS-I plates overlap, and neighbouring fields were exposed on
different nights, so for the 64.8% of catalogue rows inside overlap sky the
neighbouring plate is a free second epoch. Measured on the released catalogue:
**298 of 79,548 overlap rows (0.37%) have a detection on another plate within
5″**, against a 2.87% displaced-null — *below* chance, the expected signature
of veto-survivor sky. Two consequences, stated in both directions: running the
stage would remove almost nothing (~0.2% of the catalogue), and the released
overlap rows are almost never persistent-across-epochs objects. It is not run
because a row it flags is a source seen at two epochs — removing it changes the
catalogue's claim from "on the plate and absent from the modern catalogues we
checked" to an assertion about transience, which is a different catalogue, to
be built deliberately or not at all.

**Why:** all three are online services. At survey scale they are extremely slow,
and their back ends routinely cancel large submitted jobs, so a full-catalogue
pass is a project in itself rather than a pipeline step. A local VSX mirror
exists; SkyBoT and SuperCOSMOS have no practical local equivalent.

> **SuperCOSMOS has since been measured on THIS catalogue and the figure below
> does not transfer.** Run against all 122,820 released rows it removes
> **40.0%**, not 0.5% — see [`POSTPROCESS_STAGES.md`](POSTPROCESS_STAGES.md).
> The two numbers are not in conflict: 0.5% was measured on a pool of 11,027
> that had already passed morphology and shape stages, while S0 is the raw
> post-veto population. Read the 0.5% below as a property of that filtered pool
> only. The same page measures PTF (−2.75%) and reports SkyBoT as not run with
> a measured per-query rate and a yield upper bound.

**What this means for the numbers — measured, not assumed.** Solano et al. (2022)
indicate stages of this kind remove of order 3,000-4,000 objects. When they were
actually run in the predecessor pipeline (VASCO60), against **11,027 survivors**,
the measured removal was far smaller:

| stage | measured removal | rows |
|---|---:|---:|
| SuperCOSMOS | 0.5% | ~55 |
| PTF | 0.1% | ~11 |
| VSX | **0%** | 0 |
| SkyBoT | **0%** | 0 |

**~63 rows in total, about 0.6% of survivors** — and all of it from SuperCOSMOS
and PTF. **SkyBoT and VSX each removed nothing whatsoever.**

That is worth stating as a measured zero rather than an absence of evidence: on
**that** candidate population, neither the solar-system check nor the
variable-star check had any yield.

**The "half a percent" reading that used to follow here has been retired.** It
generalised a 0.5% SuperCOSMOS removal from an 11,027-row morphology-filtered
pool to this catalogue, where the measured figure is 40.0%. The zeros for SkyBoT
and VSX may well carry over — SkyBoT's own calibration here found no
solar-system objects in 393 sampled rows — but the SuperCOSMOS number plainly
does not.

(A plausible reading of the SkyBoT zero is that moving objects trail on a
photographic exposure and are already removed by the `ELONGATION < 1.3` gate long
before the solar-system check sees them. That is inference, not measurement.)

The two figures are not necessarily in conflict: the paper's number applies at a
different point in its chain, against a much larger input set, so the denominators
differ. But on a post-veto survivor pool of this size, the measured effect is
small.

So the honest statement is: our counts omit a cut worth **roughly half a percent**
on comparable input, effectively all of it SuperCOSMOS. Our survivor pool here is bigger, so the
absolute number would rise and the rate could differ — but this is a minor
caveat, not a reason to treat the catalogue as provisional.

This is a resource decision, not a methodological objection, and it is
reversible: the stages are present and wired, and running them later changes only
the tail of the chain.

Cone radius, when SkyBoT is run: **60′** (paper).

## Deduplication

| parameter | value | origin |
|---|---|---|
| tolerance, raw plate WCS | **3.0″, applied globally** | measured |
| tolerance, WCS-fixed coordinates | 0.25″ | measured |

> **These two are coupled and must never be decoupled.** The 0.25″ figure is
> only valid for WCS-fixed coordinates, where every tile has been aligned to
> Gaia. On the raw plate WCS each tile carries its own unaligned solution and
> the same source seen through two overlapping tiles lands arcseconds apart —
> applying 0.25″ there leaves the great majority of cross-tile duplicates
> undetected and inflates the candidate count by several percent.
>
> The WCS-fix stage defaults **on**, and the released catalogue
> (`results/s0-642-20260814/`) is WCS-fixed with the 0.25″ tolerance; the
> raw-WCS variant in `docs/REPRODUCING.md` uses 3.0″. Whichever catalogue you
> build, **change both settings together.**

Regenerate the justification with `tools/dedup_radius_sweep.py`: cross-tile
duplicate pairs should rise and then plateau, while intra-tile pairs — the
signature of merging genuinely distinct sources — stay at zero well past the
chosen value. Confirm the tolerance actually in force in the stage configuration
before a production run; do not assume the default.

## Astrometry

| parameter | value | origin |
|---|---|---|
| per-plate CRPIX correction | per-plate, `data/plate_crpix_table.csv` | measured |
| scatter rejection | warn above 0.2 px | engineering |
| per-tile refit distortion | none — pure TAN (`sip_degree=None`) | engineering, see below |
| refit residual guard | refuse tile above 1.0″ (`MAX_REFIT_RESID_ARCSEC`) | engineering |

Regenerate with `tools/build_plate_crpix_table.py`. Pure header arithmetic — no
catalogue, no fitting, no tuned threshold. Residual after correction is ~0.09″ on
plates that need it; see `docs/DSS_WCS_TWO_SOLUTIONS.md`.

### The ~0.1″ TAN representation term — a documented limitation, not an oversight

Each tile's WCS is a plain TAN refit against the plate's GSSS polynomial
solution (`tools/slice_plate_tiles.py`, `fit_wcs_from_points` with `proj_point`
at the tile centre). A pure TAN cannot represent the GSSS distortion exactly
over a ~1° tile, which leaves a **median ~0.10″ systematic (0.04–0.17″ across
healthy tiles)** in every tile's coordinates relative to the plate solution.
Fitting SIP distortion terms (`sip_degree=3`) removes it almost entirely
(~0.0003″, measured during the XE011 divergence diagnosis). The decision **not**
to do so is deliberate, for three reasons:

1. **SExtractor ignores SIP.** The coordinates every downstream stage consumes
   (`ALPHAWIN_J2000`/`DELTAWIN_J2000`) are computed from CTYPE/CRPIX/CD alone.
   SIP terms in the header would improve no catalogue coordinate — they would
   only make astropy-based readers disagree with the catalogue by ~0.1″, the
   same tool-divergence class documented in `docs/DSS_WCS_TWO_SOLUTIONS.md`.
   Realising the gain would require a coordinate-recompute stage after
   detection, or a refit in the TPV convention that SExtractor does read.
2. **The pipeline cannot feel 0.1″.** A measured 0.48″ coordinate change moves
   the 5″ veto rates by ≤0.05 points (the WCS-fix bias measurement in the
   release README), so ~0.1″ scales to ~0.01 points — tens of rows in a
   ~123k-row catalogue, below run-to-run noise. The optional Gaia refit fits
   ~6,000 tie points per tile to ~0.11″ downstream and absorbs any smooth error
   field anyway. And the GSSS solution itself carries ~1.15″ against Gaia:
   this term is under a tenth of the error it rides on.
3. **The term is bounded, not assumed.** The slicer computes the refit residual
   for every tile and refuses to write any tile above 1.0″, printing the worst
   tile beside the median. A regression of this representation cannot pass
   silently.

Revisit only for a science case that needs sub-0.2″ *absolute* astrometry tied
to the plate solution itself (e.g. plate-epoch proper motions, or
characterising reference-catalogue coordinate errors at the 0.1″ level). Then
do it properly: SIP plus a recompute step (or TPV), a full re-run, and its own
release — not a header-only change.

## Deviations from Solano et al. (2022)

Collected here so they stay auditable rather than silent. Rationale for the first
three is in the README.

1. **Simplified tessellation** — plain 7×7 grid, no corner-avoidance gate.
2. **No 30′ circular cut** — full square tiles retained.
3. **Grid in pixel space, not RA/Dec** — the latter breaks near the pole.
4. **PS1 rather than USNO-B1.0** for the spike mask.
5. **WCS-fix (per-tile Gaia refit) in the released catalogue** — no equivalent
   in the original; defaults on and was left on unintentionally, disclosed in
   the release README with the circularity measurement. Coupled to the dedup
   tolerance above; the raw-WCS build in `docs/REPRODUCING.md` is the
   paper-faithful variant.
6. **Per-plate CRPIX correction** — no equivalent in the original; corrects a
   defect specific to slicing full-plate scans.
7. **SkyBoT, SuperCOSMOS and VSX not applied** — see above. Measured against
   this catalogue, SuperCOSMOS alone would remove **40.0%**, making it by a wide
   margin the deviation with the largest effect on the final row count
   ([`POSTPROCESS_STAGES.md`](POSTPROCESS_STAGES.md)).
8. **Full-plate slicing searches sky on every covering plate** — a cutout
   pipeline searches each position once, on the plate a DSS service selects.
   Flagged per row, not removed: `tools/build_primary_plate_flags.py` +
   `tools/check_primary_counterparts.py` partition the released catalogue by a
   zero-parameter nearest-centre rule (validated at 99.04% against 11,727
   recorded STScI plate selections). 44.4% of the released rows are
   single-plate content in multiply-searched sky. Filtering the partition
   costs a measured 6.5% of matches to the public vanish-possi catalogue, so
   both counts are always quoted together.

## Conventions

- Tile identifiers: `tile_RA<ra>_DECp<dec>` / `tile_RA<ra>_DECm<dec>` — `p`/`m`
  rather than `+`/`-`, which are awkward in paths.
- Every tile carries `tile_status.json`, enabling delta runs and audit.
- Every stage writes a `.json` ledger recording `in_rows`, `out_rows` and failure
  reasons. A stage that cannot say how many rows it consumed and emitted is not
  auditable, and this pipeline's central claim is auditability.
