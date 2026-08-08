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
executed** for the current results.

| stage | script | status |
|---|---|---|
| SkyBoT — known solar-system objects | `scripts/stage_skybot_post.py` | not run |
| SuperCOSMOS — plate-artifact discrimination | `scripts/stage_supercosmos_post.py` | not run |
| VSX — known variable stars | `scripts/stage_vsx_post.py` | not run |

**Why:** all three are online services. At survey scale they are extremely slow,
and their back ends routinely cancel large submitted jobs, so a full-catalogue
pass is a project in itself rather than a pipeline step. A local VSX mirror
exists; SkyBoT and SuperCOSMOS have no practical local equivalent.

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
this candidate population, neither the solar-system check nor the variable-star
check has any yield. Of the three stages omitted here, two are known to cut
nothing, so the real cost of the omission is SuperCOSMOS alone — of order half a
percent.

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

> **These two are coupled and must never be decoupled.** The 0.25″ figure is only
> valid for WCS-fixed coordinates. This pipeline does **not** apply WCS-fix, so
> each tile carries its own unaligned solution and the same source seen through
> two overlapping tiles lands arcseconds apart. Applying 0.25″ to unaligned
> coordinates leaves the great majority of cross-tile duplicates undetected and
> inflates the candidate count by several percent.
>
> **If you enable WCS-fix, revert the tolerance to 0.25″ in the same change.**

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

Regenerate with `tools/build_plate_crpix_table.py`. Pure header arithmetic — no
catalogue, no fitting, no tuned threshold. Residual after correction is ~0.09″ on
plates that need it; see `docs/DSS_WCS_TWO_SOLUTIONS.md`.

## Deviations from Solano et al. (2022)

Collected here so they stay auditable rather than silent. Rationale for the first
three is in the README.

1. **Simplified tessellation** — plain 7×7 grid, no corner-avoidance gate.
2. **No 30′ circular cut** — full square tiles retained.
3. **Grid in pixel space, not RA/Dec** — the latter breaks near the pole.
4. **PS1 rather than USNO-B1.0** for the spike mask.
5. **No WCS-fix** — coupled to the dedup tolerance above.
6. **Per-plate CRPIX correction** — no equivalent in the original; corrects a
   defect specific to slicing full-plate scans.
7. **SkyBoT, SuperCOSMOS and VSX not executed** — see above. This is the
   deviation with the largest effect on the final row count.

## Conventions

- Tile identifiers: `tile_RA<ra>_DECp<dec>` / `tile_RA<ra>_DECm<dec>` — `p`/`m`
  rather than `+`/`-`, which are awkward in paths.
- Every tile carries `tile_status.json`, enabling delta runs and audit.
- Every stage writes a `.json` ledger recording `in_rows`, `out_rows` and failure
  reasons. A stage that cannot say how many rows it consumed and emitted is not
  auditable, and this pipeline's central claim is auditability.
