
# VASCO60 Prewarmers (Unified)

This bundle provides a **unified neighbourhood prewarmer** and thin **background wrappers** for PS1 and USNO-B.

## Files

- `scripts/prewarm_neighbourhood_cache.py` — single runner for **Gaia/PS1/USNO-B**.
- `scripts/prewarm_ps1_neighbourhood_bg.sh` — nohup/PID/STOP wrapper calling the unified runner with `--catalog ps1`.
- `scripts/prewarm_usnob_neighbourhood_bg.sh` — nohup/PID/STOP wrapper calling the unified runner with `--catalog usnob`.

> Gaia already has an existing background wrapper in your repo; keep it or adapt its call to the unified runner.

## Tile naming contract
Matches directory names of the form `tile_RA<ra>_DECp<dec>` or `tile_RA<ra>_DECm<dec>` (e.g., `tile_RA130.013_DECp33.081`). The regex in the runner is:

```python
_TILE_RE = re.compile(r"^tile_RA(?P<ra>\d+(?:\.\d+)?)_DEC(?P<hem>[pm])(?P<dec>\d+(?:\.\d+)?)$")
```

## Tile geometry vs analysis circle

Vasco60 downloads **60×60 arcmin** tiles (step1). From step2 onward, sources are
restricted to a **≤30′ circle** around the tile center (`VASCO_CIRCLE_ARCMIN=30`).
Prewarm radii are sized for this **analysis circle**, not the download tile:

- A source at 30′ from center needs a catalogue match within 5″ → catalogue must
  cover 30′ + 5″ ≈ **31′**.
- Sources in the 30′–42.4′ annulus (corners of the 60×60 square) are cut before
  any xmatch result is used, so they do not drive the prewarm radius requirement.

## Radius policy

| Catalog | Shell default | Effective radius | Rationale |
|---|---|---|---|
| Gaia | 31′ | 31′ | 30′ circle + 5″ xmatch safety |
| USNO-B | 31′ | 31′ | same |
| PS1 | 35′ | 34.8′ | spike halo buffer beyond the 30′ edge; `VASCO_PS1_RADIUS_DEG=0.58` exported to bypass the 0.5° fetcher cap |

**PS1 cap note:** `fetch_ps1_neighbourhood` caps the VizieR cone at 0.5° (30′) by
default. The background wrapper exports `VASCO_PS1_RADIUS_DEG=0.58` so the Python
child process sees it and the full 35′ request reaches VizieR. If you call the
runner directly, either pass a small enough radius (≤30′) or export the variable
yourself:
```bash
export VASCO_PS1_RADIUS_DEG=0.58
python scripts/prewarm_neighbourhood_cache.py --catalog ps1 --radius-arcmin 35 ...
```

## Examples
```bash
# Gaia (31′)
python scripts/prewarm_neighbourhood_cache.py --catalog gaia  --tiles-root ./data/tiles --radius-arcmin 31 --workers 6

# PS1 (35′ — export env var to bypass 0.5° fetcher cap)
export VASCO_PS1_RADIUS_DEG=0.58
python scripts/prewarm_neighbourhood_cache.py --catalog ps1   --tiles-root ./data/tiles --radius-arcmin 35 --workers 4

# USNO‑B (31′)
python scripts/prewarm_neighbourhood_cache.py --catalog usnob --tiles-root ./data/tiles --radius-arcmin 31 --workers 4
```

## Notes
- The unified runner writes per‑tile CSVs under `<tile>/catalogs/` with catalog‑specific names.
- USNO‑B requires `fetch_usnob_neighbourhood()` to be present in `vasco/external_fetch_online.py`. If not, the runner will raise a helpful error.
- All scripts are **resumable** and write a progress JSON into `./logs/`.
