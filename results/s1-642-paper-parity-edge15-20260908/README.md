# S1 — paper-parity S0 with the plate-edge cut applied (2026-09-08)

**95,966 rows.** The paper-parity S0 catalogue
([`../s0-642-paper-parity-20260828/`](../s0-642-paper-parity-20260828/),
134,976 rows) with every row lying within **15 arcminutes of its own plate's
array boundary** removed — 39,010 rows, 28.90%.

This is a **derived, cleaned** product. It is published because the cut is
cheap to describe and expensive to reproduce: regenerating it needs the full
plate scans and a per-plate WCS pass, which not everyone has the disk or the
time for. The numbers behind it were already published as figures; this is
the actual file.

**It does not supersede anything.** S0 remains the catalogue to cite. This
is S0 minus a geometric region, offered for anyone who wants the
edge-cleaned set without rebuilding it.

## Why 15′, and what the cut removes

Distance to the plate's *array boundary* is a far sharper lever than
distance from the plate centre: a comparable radial cut removes ~65% of S0,
while the edge cut removes 21% at 10′ and 28.9% at 15′ — at **essentially
zero cost in recall against published comparison catalogues**. Full
derivation in [`../../docs/PLATE_EDGE_MASK.md`](../../docs/PLATE_EDGE_MASK.md).

Three independent lines of evidence say the removed population is
instrumental, not sky:

1. **It is nearly absent from published comparison catalogues.** Rows within
   15′ of an array edge match a published comparison list at **0.17%**,
   against 55–64% for rows more than 30′ in — a ~300× difference. Dropping
   the whole zone costs **0.00%** recall against the public 5,399-row
   `vanish-possi` list.
2. **SuperCOSMOS has already removed it.** Only 0.36% of the post-SuperCOSMOS
   set sits within 5′ of an edge, against 7.6% in S0 — an independent
   catalogue reached the same conclusion without being asked to.
3. **The detections are plate-specific.** They appear on one plate's rim and
   not on the overlapping neighbour, where the same sky sits ≥18′ inside the
   array — the signature of vignetting, boundary structure and emulsion
   damage at a physical plate edge rather than of anything on the sky.

Two explanations were tested and **rejected**: that the effect is an artifact
of our own cross-plate dedup tie-break (refuted — the pre-dedup catalogue
shows the identical step), and that it is a quality gradient among survivors
(not supported — `ELONGATION` and `FWHM` are flat across every edge bin,
because the MNRAS gate has already flattened them).

## Honest limitation

The recall argument is **partly circular**: the comparison catalogues may
themselves under-sample plate edges, so there is little there to lose by
construction. That is why points 2 and 3 above matter — they are the
independent legs. Recall alone could not have settled this, and is not what
selected the threshold.

The threshold itself is a **round number chosen from a documented yield
curve**, not fitted: 5′ / 10′ / 15′ remove 6.8% / 19.3% / 28.9% with no
recall cost at any of them. 15′ was selected as the point where the
match-rate step has fully completed (0.17% inside, 8.8% at 15–20′, 37.5% at
20–30′), not by optimising an outcome.

## Files

| file | rows | what |
|---|---:|---|
| `stage_S1.csv.gz` | 95,966 | the cleaned catalogue |
| `edge_flags.csv.gz` | 134,976 | per-row audit for **all** S0 rows: plate and edge distance, so the cut can be checked or re-cut at another threshold |
| `stage_S1_ledger.json` | — | parameters and counts as the stage recorded them |
| `SHA256SUMS` | — | integrity of the files as shipped |

`stage_S1.csv.gz` columns: `src_id, tile_id, object_id, ra, dec, plate_id,
edge_dist_arcmin`. The first five match S0 exactly; the last two are added so
the cut is self-documenting. Minimum `edge_dist_arcmin` in this file is
15.000 by construction.

## Reproducing it

```
python scripts/stage_edge_post_v2.py \
    --run-dir <run> --input-glob 'stages/stage_S0.csv' \
    --plate-map-csv <tile_id,plate_id> \
    --stage S1 --cut --core-radius-deg 5.0 --min-edge-arcmin 15
```

`--core-radius-deg 5.0` is deliberate: it opens the radial core cut past the
plate corner distance (~4.68°) so it removes nothing, leaving the edge
distance as the only active criterion. The stage is **off by default and not
wired into the pipeline** — applying it is an explicit per-run decision, and
it prints a warning when active, because the papers describe no such
cleaning step.
