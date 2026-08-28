# S0 — paper-parity build, 642 plates (2026-08-28)

**134,976 rows.** A second candidate catalogue, built from the same public
POSS-I plate inputs as [`../s0-642-20260814/`](../s0-642-20260814/) (the
primary release, 122,820 rows), but with two deliberate configuration
changes intended to track the literal method described in Solano et al.
(2022) / MNRAS 2022 more closely than the primary release's own operational
choices do. **This does not supersede the primary release** — the two are
companion builds, kept side by side, and the primary release remains the
one to cite unless you specifically want the paper-parity variant.

Produced entirely from public inputs, same as the primary release — see
[`tools/audit_independence.py`](../../tools/audit_independence.py).

**Hot-tile survivors, in both this build and the primary release, have now
been checked against an independent second exposure** (the POSS-I O/blue
plate) — direct, statistically significant evidence for a plate-emulsion or
scan-defect origin on a majority of the tiles sampled, not just the last
hypothesis left standing by elimination:
[`docs/HOT_TILE_BLUE_PLATE_CHECK.md`](../../docs/HOT_TILE_BLUE_PLATE_CHECK.md).

## Files

| file | rows | what |
|---|---:|---|
| `stage_S0.csv.gz` | 134,976 | the catalogue |
| `tile_manifest.csv.gz` | 31,458 | every tile processed, with its plate |
| `verification_s0_gaia_invariant.csv.gz` | 26,277 | per-tile Gaia-contamination check — **PASS**, 0.02% (see below) |
| `RUN_SUMMARY.txt` | — | parameters and counts as the build recorded them |
| `SHA256SUMS` | — | integrity of the files as shipped |
| `SHA256SUMS.uncompressed` | — | integrity of the *contents* |

`stage_S0.csv` columns: `src_id, tile_id, object_id, ra, dec` — identical
schema to the primary release.

**Not shipped here, unlike the primary release**: `primary_plate_flags.csv.gz`
(the coverage-partition sidecar) and its accompanying analysis, and
`repaired_astrometry_tiles.csv`. The astrometry defect that file documents
was a one-time bug in the shared slicer, fixed upstream well before this
build — nothing here needed a separate repair pass. The coverage-partition
analysis is a known gap, not reproduced for this release; see the primary
release's own README for that methodology if you want to run it against
this catalogue yourself.

## What this is: two configuration differences from the primary release

Both builds share the same plate inputs, the same slicing code (steps 1–3
are byte-identical between the `main` and `paper-parity` branches), the
same MNRAS filter chain, WCS-fix, and 0.25″ dedup tolerance. They differ in
exactly two places — see
[`BRANCH_PAPER_PARITY.md`](../../BRANCH_PAPER_PARITY.md) for the full lever
history:

| | primary release | this build (paper-parity) |
|---|---|---|
| veto catalogues | Gaia + PS1 + USNO-B | **Gaia + PS1 only** |
| spike mask | PS1 | **USNO-B** |

The primary release's choice to include USNO-B in the veto and use PS1 for
the spike mask is VASCO60's own operational decision, not something either
paper specifies. This build drops USNO-B from the veto and uses it for the
spike mask instead, tracking the MNRAS 2022 description more literally.

### The dedup tie-break fix

Both builds dedup at 0.25″, but this build carries a correctness fix to
*which* row survives inside a duplicate cluster. The prior rule picked the
lexicographically smallest `tile_id:object_id` string — arbitrary, no
relationship to detection quality. This build picks the cluster member with
the largest distance to its own plate's array edge, falling back to the old
rule only when a member's plate can't be resolved. Row *count* after dedup
is unaffected (same clusters, same drop count) — only which representative
survives changes.

At this run's scale, the fix has a real, non-trivial effect: **7,976
multi-member clusters, of which 4,096 (51.4%) picked a different
representative than the old rule would have.** The primary release predates
this fix.

## Results

| | rows |
|---|---:|
| primary release (`../s0-642-20260814/`) | 122,820 |
| this build | 134,976 |
| difference | **+9.90%** (+12,156 rows) |

The excess is attributable to dropping USNO-B from the veto (more rows
survive that would otherwise be vetoed as catalogued stars), partially
offset by the spike-mask switch to USNO-B (a net-negative contribution).

## Comparison to the primary release — what actually changed underneath

The +9.90% headline is a net figure and hides substantially more churn than
it suggests. Positional match between the two catalogues (stable across a
0.25″–5″ tolerance sweep, so not a matching-radius artifact):

| | count | % |
|---|---:|---|
| in common (both catalogues) | 115,946 | Jaccard 81.7% |
| primary-only (dropped in this build) | ~6,875 | 5.60% of the primary release |
| this-build-only (genuinely added) | ~19,029 | 14.10% of this build |

So the net +9.90% sits on top of roughly 26,000 rows (19–21%) of underlying
churn — rows added and dropped that don't show up in the row-count delta
alone.

**Attribution**: a controlled experiment on a 9-plate sample (holding raw
detections and tiling fixed, varying only the two levers) found the lever
difference alone accounts for essentially all of the churn — comparing the
primary release's own rows against the same lever configuration rebuilt on
freshly-sliced tiles matched at **100.0%** (varies only in whether it's
old-release-build vs newly-sliced tiles), while switching only the levers
(new tiles, old levers vs new tiles, paper-parity levers) reproduced almost
exactly the same churn rate as comparing the two full releases directly.
In short: **the two configuration changes above explain the difference
completely** — there is no unexplained mechanism, no code-drift effect, and
no tiling/astrometry difference doing meaningful work here.

## Independent post-process measurement (SuperCOSMOS, PTF)

Measured, not applied — same "measure only" posture as the primary
release. SCOS and PTF are pure per-row tests (no population-derived
threshold), so results are order-independent between them.

| | primary release | this build |
|---|---:|---:|
| SCOS unconfirmed (R1 arm) | 40.0% | **36.44%** (49,186 rows) |
| PTF match, of covered rows | 3.40% | **3.64%** |
| PTF coverage | ~80.9% | 80.7% (108,932 rows) |

SCOS confirmation is slightly *better* on this build, consistent with the
excess being ordinary USNO-B-catalogued stars rather than scan artifacts.

### Chained funnel (SCOS → PTF → edge distance)

Three variants, all measurement only — nothing here is applied to
`stage_S0.csv` or shipped as a separate file.

| stage | rows in | rows out | % of this build's S0 |
|---|---:|---:|---:|
| S0 | — | 134,976 | 100% |
| SCOS (keep matched) | 134,976 | 85,790 | 63.6% |
| + PTF (drop matched) | 85,790 | 83,198 | 61.6% |
| + edge distance ≥10′ from array edge | 83,198 | 81,409 | 60.3% |

The primary release's equivalent SCOS→PTF chain: 122,820 → 73,681 → 71,654
(58.3%) — same shape. The edge-distance step removes only 2.15% at this
point in the chain (far below the 19.30% it removes from the whole S0
directly), because SCOS has already removed most of the near-edge
population as a side effect.

A fourth variant — a 2° radial cut from plate centre, matching PASP2025's
own criterion — was also measured (12,903 rows, 9.56% of S0) but is **not
included above**: that cut is dominated by tile-area geometry rather than
being selective for quality, and neither paper applies it as a blanket
catalogue filter, only inside specific analyses. Included here only as a
caveat against being asked for it, not as a recommended figure.

## What this catalogue is not

Not a replacement for the primary release, not independently re-verified
against SuperCOSMOS/PTF beyond the measurement above, and not carrying the
primary release's coverage-partition (`primary_plate_flags.csv.gz`)
analysis. Treat it as a companion build for anyone specifically interested
in how closely the published method (rather than VASCO60's own operational
choices) reproduces at full survey scale.

## Verify it

```bash
sha256sum -c SHA256SUMS               # the files as shipped
zcat stage_S0.csv.gz | sha256sum      # must equal the uncompressed hash below
VASCO_GAIA_CACHE=<gaia_mirror> python3 ../../tools/check_s0_gaia_invariant.py \
    --s0-csv <(zcat stage_S0.csv.gz) --out /tmp/s0_gaia.csv
```

The invariant must PASS: the fraction of rows with a Gaia source within 5″
must be ~0, because the veto removed exactly those. It reports **0.02%**
here — identical to the primary release's own figure. 0 of the 9,781 tiles
with ≥5 rows exceed the per-tile threshold.

### Content hashes — quote these, not the `.gz` ones

`SHA256SUMS` covers the `.gz` files and verifies the *transfer* only — gzip
output is not reproducible across implementations. The hashes below (also
in `SHA256SUMS.uncompressed`) are what identify the catalogue and what a
citation should pin.

| file | rows | bytes | sha256 |
|---|---:|---:|---|
| `stage_S0.csv` | 134,976 | 12,087,327 | `cdd8e5500819eff752fc0d7a630a7909f58e57f53eb3391022ecc3bb976902fa` |
| `tile_manifest.csv` | 31,458 | 3,043,231 | `5c7d8848b6333dae54b6f95cc5985148a786c01a93f727c42b788445504e0e8a` |
| `verification_s0_gaia_invariant.csv` | 26,277 | 886,818 | `0841d5e251ff06d4ce8628b59a03161df7834962ac5ac7916d7aa7ccd867ed24` |

## Licence and attribution

Same terms as the primary release — see
[`../s0-642-20260814/README.md`](../s0-642-20260814/README.md#licence-and-attribution).
