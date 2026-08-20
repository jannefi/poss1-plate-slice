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

**Conclusion: "no coordinate correction" and "a 5″ veto radius" are internally
inconsistent.** An alignment step is implied by the method even though the paper
never states one, which is why `main` keeps WCSFIX on.

Anything quoted from this branch must name both the arm and the stage.
