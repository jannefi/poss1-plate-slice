"""Resolve data roots without hardcoding anyone's filesystem.

Every path this project needs comes from here, in a fixed precedence:

    1. environment variable  POSS1_<KEY>    -- per-invocation override
    2. config.local.yaml                    -- your machine, untracked
    3. config.yaml                          -- the tracked default
    4. a built-in fallback                  -- almost certainly wrong for you

The point of (1) and (2) is that a clone should run without editing a tracked
file, and that a wrong path fails loudly rather than silently reading nothing.
Use `require()` when an empty result would be indistinguishable from a genuine
zero -- a cache directory that does not exist otherwise looks exactly like a
cache miss, and the pipeline will quietly fall back to live catalogue queries.

How to validate
---------------
    python3 -c "from vasco.paths import dump; dump()"

prints every key, its resolved value, and whether it exists on disk.
"""
from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]

# No machine-specific defaults: a path that is not configured must fail loudly
# rather than resolve to somebody else's filesystem and read nothing.
KEYS = ("plate_dir", "tiles_dir", "work_dir", "gaia_cache", "ps1_cache",
        "usnob_cache", "maps_cache", "plate_headers")
FALLBACK = {"work_dir": "work"}


@lru_cache(maxsize=1)
def _config() -> dict:
    cfg: dict = {}
    for name in ("config.yaml", "config.local.yaml"):   # local wins
        p = _ROOT / name
        if not p.exists():
            continue
        try:
            import yaml
            cfg.update(yaml.safe_load(p.read_text()) or {})
        except ImportError:
            # PyYAML absent: parse the flat `key: value` subset we actually use
            # rather than failing outright, so a bare clone still runs.
            for line in p.read_text().splitlines():
                line = line.split("#", 1)[0].strip()
                if ":" in line:
                    k, v = line.split(":", 1)
                    if v.strip():
                        cfg[k.strip()] = v.strip()
    return cfg


def get(key: str) -> Path:
    """Resolved path for `key`. Never raises; may point at something absent."""
    env = os.environ.get(f"POSS1_{key.upper()}")
    if key not in KEYS:
        raise KeyError(f"unknown path key {key!r}; known: {sorted(KEYS)}")
    raw = env or _config().get(key) or FALLBACK.get(key)
    if raw is None or str(raw).startswith("<"):
        raise SystemExit(
            f"[FATAL] path {key!r} is not configured.\n"
            f"        Set it in config.local.yaml, or export POSS1_{key.upper()}.\n"
            f"        See config.yaml for what it should point at.")
    p = Path(str(raw)).expanduser()
    return p if p.is_absolute() else (_ROOT / p)


def require(key: str) -> Path:
    """As `get`, but fails loudly if it is not on disk."""
    p = get(key)
    if not p.exists():
        raise SystemExit(
            f"[FATAL] {key} -> {p} does not exist.\n"
            f"        Set it in config.local.yaml or export POSS1_{key.upper()}."
        )
    return p


def dump() -> None:
    for k in sorted(KEYS):
        try:
            p = get(k)
            print(f"  {k:15s} {'OK ' if p.exists() else 'MISSING'}  {p}")
        except SystemExit:
            print(f"  {k:15s} UNSET")


if __name__ == "__main__":
    dump()
