# This branch is an EXPERIMENT, not a recommended configuration

`paper-parity` exists to measure one thing: how many more candidates the
pipeline yields when configured to follow the *text* of Solano et al. (2022)
more literally than `main` does. It is published so the work is durable and
auditable. **Do not treat it as the configuration to reproduce VASCO with.**

## What it changes relative to `main`

| # | lever | this branch | `main` |
|---|---|---|---|
| 1 | coordinate correction | WCSFIX **off** | WCSFIX **on** |
| 2 | dedup tolerance | 3.0″ (paired with lever 1) | 0.25″ |
| 3 | cross-match veto | Gaia + PS1 | Gaia + PS1 + USNO-B |
| 4 | diffraction-spike catalogue | USNO-B | PS1 |

Levers 1 and 2 are locked together and **must never be decoupled**: a dedup
tolerance is only meaningful relative to the astrometric frame it operates in.

## What the experiment found

Turning WCSFIX off produced a headline **+21.1%** excess over the published
release — and roughly **half of that was veto failure, not signal**. The paper
specifies a 5″ cross-match radius, which presupposes astrometry better than 5″;
raw plate WCS does not always deliver that. Where the raw-WCS error approaches
5″, the Gaia veto stops matching catalogued stars and they survive as
candidates. On the worst tile measured, 272 new "candidates" were catalogued
stars: **0.4%** matched Gaia at raw coordinates, **94.5%** matched once shifted
into the corrected frame. Tiles with ≥3″ error were 2.8% of tiles but 58% of the
excess.

Re-running the same 17 plates with WCSFIX **kept** — so that failure mode is
absent by construction — gives **+11.8%** (4,075 vs 3,644 rows), and that is the
defensible paper-parity figure. 94.5% of release rows are retained, and the
excess is broad (16 of 17 plates, spread 3-42%) rather than concentrated in a
few pathological tiles.

**Conclusion — corrected 2026-08-21.** An earlier revision of this file argued
that the paper's method *implied* an unstated alignment step. That inference went
beyond the data and is withdrawn: it generalized this pipeline's astrometric
error tail to the original pipeline. The defensible statement is narrower.

A 5″ cross-match veto presupposes astrometry comfortably inside 5″. Whether the
original pipeline met that with its native plate solutions alone, or by some
unstated step, cannot be determined from the paper. What *is* measured here
applies to this pipeline: slicing tiles from full-plate scans leaves a
field-dependent astrometric residual (median ~0.5″, p90 ~1.7″, max ~6.8″ against
Gaia), and on the small tail of tiles where that residual approaches 5″ a
raw-coordinate veto silently passes catalogued stars — the failure measured
above. WCSFIX exists to remove that failure mode, and the rows it keeps out are
Gaia-catalogued stars. `main` keeps it on as an engineering requirement of this
pipeline's own astrometry, not as a reconstruction of anything the original
pipeline is claimed to have done.

Per-lever attribution (2026-08-21) reinforces this: regenerating the release
configuration post-hoc from the same detections shows the +11.8% is essentially
all lever 3 — **88.9% of the gained rows sit on a USNO-B source within 5″**,
against 0.0% of the rows both configurations share — with lever 4 net negative.
The excess measures a veto stage, not astrometry, and not missed transients.

Anything quoted from this branch must name both the arm and the stage.
