## Post-pipeline stages (run-scoped shrinking set)

**Goal:** Produce a run-scoped folder containing the current shrinking survivor set (S0 → S0G → S1 → S2 …), plus provenance and artifacts for each stage.  
**Canonical run folder:** `./work/runs/run-S1-<date>/`

### Overview (what this produces)
- A **master audit CSV** for the run
- A **science-grade dedup CSV**
- A **current survivor set** (“edge-core”) that shrinks forward at each stage
- Per-stage artifacts:
  - `stage_SX_*.csv` (carry-forward survivors)
  - `stage_SX_*_flags.csv` (per-row flags for auditing)
  - `stage_SX_*_ledger.json` (counts + parameters + stats)

> **Invariant:** Every stage outputs a smaller “current survivors” CSV used as the input for the next stage. 

---

### 1) Build run-scoped stage CSVs (S0)
**Script:** `scripts/build_run_stage_csvs.py`  
**Purpose:** Create the initial run folder with master CSVs, the current survivor set and upload chunks.

**Run (typical):**
```sh
RUN=./work/runs/run-S1-$(date +%Y%m%d_%H%M%S)
python scripts/build_run_stage_csvs.py --run-tag "$(basename "$RUN")"
```

**Key outputs (under $RUN/):**

- source_extractor_final_filtered.csv (master, audit)
- source_extractor_final_filtered__dedup.csv (science-grade dedup)
- stage_S0.csv (minimal stage view)
- upload_positional.csv and upload_positional_chunk_*.csv (S0 upload view for next fetcher)
- tile_manifest.csv, RUN_SUMMARY.txt, allow/exclude list copies


### 2) SPREAD_MODEL gate (S0G)

**Script:** `scripts/stage_spread_model_post.py`  
**Purpose:** Reject candidates whose restored, PSF-based SPREAD_MODEL falls at or below the
locked `-0.002` threshold (`context/02_DECISIONS.md`). The primary pipeline runs single-pass
SExtractor (no PSFEx, no SPREAD_MODEL column) to keep candidate generation cheap, so this gate
restores a real morphology signal at the post-process stage instead.

This is a **thin CSV-join stage**, not a compute stage — it does not run PSFEx/SExtractor
itself. It reads each tile's `catalogs/spread_model_postscore.csv`, already written by a
separate, explicit run of `tools/spread_model_postscore.py` (a validated crop + synthetic-PSF
remeasurement, ~96.7% agreement with full-field ground truth, ~8.4s/survivor — scales with
survivor count, not a fixed per-tile cost, so a full corpus can take hours and is deliberately
kept out of the fast stage-chain path). Run that tool first, over every tile the input stage CSV
references, before running this stage:

```sh
python tools/spread_model_postscore.py \
  --tile-ids-file /path/to/tile_ids.txt \
  --tiles-root ./data/tiles_archive \
  --workers 10 \
  --skip-existing
```

`--skip-existing` makes the (potentially multi-hour) precompute resumable: a tile is skipped
without touching PSFEx/SExtractor if its `catalogs/spread_model_postscore.csv` already exists
and is non-empty.

**Run the gate:**
```sh
python scripts/stage_spread_model_post.py \
  --run-dir "$RUN" \
  --input-glob 'stage_S0.csv' \
  --stage S0G \
  --tiles-root ./data/tiles_archive \
  --spread-model-min -0.002
```

**Expected outputs (under $RUN/stages/):**

- stage_S0G_SPREAD.csv (carry forward)
- stage_S0G_SPREAD_flags.csv
- stage_S0G_SPREAD_ledger.json

> **Note (conservative-keep policy):** a row whose tile has no precomputed postscore file yet,
> whose object wasn't found in that file, or whose crop remeasurement didn't match a two-pass
> detection is **kept**, not dropped or rejected — the same conservative-keep-and-flag
> convention used by the S6 declination gate below. This means `stage_S0G_SPREAD.csv`'s row
> count is an **upper bound** until the ledger's `postscore_coverage.tiles_missing_postscore`
> is `0` and its three `rows_kept_via_*` counts are all `0` — check the ledger before treating a
> run's count as final rather than an estimate.


### 3) SkyBoT stage (run once, keep artifacts, shrink forward)

SkyBoT is typically a small cutter. Run it once, keep artifacts, then shrink forward without requerying.

Foreground run:
```sh
RUN=./work/runs/run-S1... \
STAGE=S1 \
INPUT='stages/stage_S0G_SPREAD.csv' \
bash scripts/run_skybot_stage_bg.sh start
```

**Expected outputs (under $RUN/stages/):**

- stage_S1_SKYBOT.csv (carry forward)
- stage_S1_SKYBOT_flags.csv
- stage_S1_SKYBOT_ledger.json


### 4) SuperCOSMOS stage (shrink to S2)

Script: `scripts/stage_supercosmos_post.py`

Run
```sh
python scripts/stage_supercosmos_post.py \
  --run-dir "$RUN" \
  --input-glob 'stages/stage_S1_SKYBOT.csv' \
  --stage S2 \
  --radius-arcsec 5 \
  --chunk-size 5000 \
  --mode keep_matches
```

**Expected outputs (under $RUN/stages/):**

- stage_S2_SCOS.csv (carry forward)
- stage_S2_SCOS_flags.csv
- stage_S2_SCOS_ledger.json

### 5) PTF stage (shrink to S3)

Script: `scripts/stage_ptf_post.py`
Run:
```sh
python scripts/stage_ptf_post.py \
  --run-dir "$RUN" \
  --input-glob 'stages/stage_S2_SCOS.csv' \
  --stage S3 \
  --radius-arcsec 5 \
  --ptf-table ptf_objects
```

**Expected outputs (under $RUN/stages/):**

- stage_S3_PTF.csv
- stage_S3_PTF_flags.csv
- stage_S3_PTF_ledger.json

### 6) VSX stage (local mirror; shrink forward)

**Prerequisite:** This stage requires a local mirror of the VSX (Variable Star Index) catalogue. Use the bootstrap script from the vasco tooling to fetch it:

```bash
bash tools/fetchers/bootstrap_vsx_catalog.sh
```

Source: [bootstrap_vsx_catalog.sh](https://github.com/jannefi/vasco/blob/main/tools/fetchers/bootstrap_vsx_catalog.sh)

Script: `scripts/stage_vsx_post.py`
Run:
```sh
python scripts/stage_vsx_post.py \
  --run-dir "$RUN" \
  --input-glob 'stages/stage_S3_PTF.csv' \
  --stage S4 \
  --radius-arcsec 5
```

**Expected outputs (under $RUN/stages/):**

- stage_S4_VSX.csv
- stage_S4_VSX_flags.csv
- stage_S4_VSX_ledger.json

### 7) Declination scope gate — S6 (science-grade subset)

**Script:** `scripts/stage_scope_dec_post.py`  
**Purpose:** Restrict the surviving set to the declared science scope by declination. The default threshold (Dec ≥ 0°) aligns with the MNRAS 107k comparison list, which covers the northern sky. Sources below the threshold are flagged as `dec_below_scope` and removed; sources whose Dec cannot be parsed are kept conservatively.

Run (default northern scope):
```sh
python scripts/stage_scope_dec_post.py \
  --run-dir "$RUN" \
  --dec-min 0.0
```

Run (MNRAS-comparison view, allows slightly negative Dec):
```sh
python scripts/stage_scope_dec_post.py \
  --run-dir "$RUN" \
  --dec-min -3.0
```

**Expected outputs (under $RUN/stages/):**

- `stage_S6_SCOPE_DEC.csv` — carry-forward survivors
- `stage_S6_SCOPE_DEC_flags.csv` — per-row flags (`dec_value`, `reject_reason`, `is_rejected`)
- `stage_S6_SCOPE_DEC_ledger.json` — parameters (`dec_min`), totals, rejection breakdown

> **Note:** This stage replaces the earlier CLASS_STAR morphology gate (S6_MORPH_CLASS), which was superseded by the SPREAD_MODEL classifier already applied upstream at the tile level.

---

### Consolidated reporting (all runs)

After delta runs exist, maintain an all-up report that spans all run folders (initial big run + later deltas).The consolidated report should also materialize a single current-survivors view across the union of runs so downstream fetchers always consume one canonical shrinking set.

**Recommended inputs:**

- $RUN/STAGE_LEDGER.csv (rows_in / rows_flagged / rows_out per stage)
- Per-run manifest: run tag, tile selection/range, created timestamp, schema/version

**Recommended outputs:**

- All-up counts per stage across runs
- Current survivors (union view) exported as CSV + upload_positional chunks for the next fetcher




