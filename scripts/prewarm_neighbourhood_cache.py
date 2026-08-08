
#!/usr/bin/env python3
from __future__ import annotations
import os, re, json, time, signal, logging, sys
from logging.handlers import RotatingFileHandler
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---------------- repo import helper ----------------
def _add_repo_to_syspath() -> Path:
    here = Path(__file__).resolve()
    add_path = None
    for p in [here] + list(here.parents):
        if (p / "vasco" / "__init__.py").exists():
            add_path = p; break
        if (p / "src" / "vasco" / "__init__.py").exists():
            add_path = p / "src"; break
    if add_path is None:
        cwd = Path.cwd().resolve()
        if str(cwd) not in sys.path:
            sys.path.insert(0, str(cwd))
        raise ModuleNotFoundError("Could not locate 'vasco' package.")
    if str(add_path) not in sys.path:
        sys.path.insert(0, str(add_path))
    return add_path

_DETECTED_REPO = _add_repo_to_syspath()

# ----------- fetchers (VizieR-backed) -----------
from vasco.external_fetch_online import (
    fetch_gaia_neighbourhood,      # I/355 via VizieR
    fetch_ps1_neighbourhood,       # II/389 via VizieR
)
try:
    from vasco.external_fetch_online import fetch_usnob_neighbourhood  # I/284 via VizieR
    _HAS_USNOB = True
except Exception:
    _HAS_USNOB = False

_TILE_RE = re.compile(r"^tile_RA(?P<ra>\d+(?:\.\d+)?)_DEC(?P<hem>[pm])(?P<dec>\d+(?:\.\d+)?)$")

CAT = {
    "gaia":  ("gaia_neighbourhood.csv",  "gaia"),
    "ps1":   ("ps1_neighbourhood.csv",   "ps1"),
    "usnob": ("usnob_neighbourhood.csv", "usnob"),
}

class StopRequested(Exception):
    pass

# ---------------- tile helpers ----------------
def iter_tile_dirs(tiles_root: Path):
    for root, dirs, _ in os.walk(tiles_root):
        for d in dirs:
            if _TILE_RE.match(d):
                yield Path(root) / d

def parse_center(name: str) -> tuple[float, float] | None:
    m = _TILE_RE.match(name)
    if not m: return None
    ra = float(m.group('ra'))
    dec = float(m.group('dec'))
    if m.group('hem') == 'm':
        dec = -dec
    return ra, dec

# ---------------- commons ----------------
def atomic_write_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)

def cache_exists(path: Path) -> bool:
    try:
        return path.exists() and path.stat().st_size > 0
    except Exception:
        return False

# ---------------- main ----------------
def main():
    import argparse
    ap = argparse.ArgumentParser(description="Unified prewarmer for Gaia/PS1/USNO-B neighbourhood caches.")
    ap.add_argument("--catalog", choices=CAT.keys(), required=True)
    ap.add_argument("--tiles-root", default="./data/tiles")
    ap.add_argument("--logs-dir", default="./logs")
    ap.add_argument("--workers", type=int, default=4)

    ap.add_argument("--radius-arcmin", type=float, default=31.0,
                    help="Neighbourhood radius (arcmin). Default 31′ for 30′ circle + 5″ safety.")
    ap.add_argument("--max-rows", type=int, default=200000,
                    help="Gaia/USNO-B row cap; PS1 uses --max-records.")
    ap.add_argument("--max-records", type=int, default=50000,
                    help="PS1 record cap.")
    ap.add_argument("--timeout", type=float, default=90.0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--retry", type=int, default=3)
    ap.add_argument("--progress-every", type=int, default=100)
    args = ap.parse_args()

    tiles_root = Path(args.tiles_root)
    logs_dir = Path(args.logs_dir); logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / f"prewarm_{args.catalog}_neighbourhood_cache.log"
    progress_path = logs_dir / f"prewarm_{args.catalog}_neighbourhood_progress.json"
    stop_file = logs_dir / f"PREWARM_{args.catalog.upper()}_NEIGH_STOP"

    logger = logging.getLogger(f"prewarm_{args.catalog}_neigh"); logger.setLevel(logging.INFO)
    handler = RotatingFileHandler(log_path, maxBytes=10_000_000, backupCount=5, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    sh = logging.StreamHandler(); sh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(sh)

    stop = {"flag": False}
    def _sig_handler(_s, _f):
        stop["flag"] = True
        logger.warning("Stop signal received; exiting after in-flight tasks complete.")
    signal.signal(signal.SIGINT, _sig_handler)
    signal.signal(signal.SIGTERM, _sig_handler)

    out_name, _ = CAT[args.catalog]
    counters = {
        "catalog": args.catalog,
        "tiles_found": 0, "tiles_scheduled": 0,
        "tiles_cached_skip": 0, "tiles_fetched": 0, "tiles_failed": 0,
        "tiles_no_center": 0, "tiles_zero_rows": 0,
        "started_at": time.strftime("%Y-%m-%d %H:%M:%S"), "last_tile": None,
        "radius_arcmin": float(args.radius_arcmin),
    }

    def write_progress():
        atomic_write_text(progress_path, json.dumps(counters, indent=2))

    def do_one(tile_dir: Path):
        if stop["flag"] or stop_file.exists():
            raise StopRequested()
        tile_id = tile_dir.name
        counters["last_tile"] = tile_id
        out = tile_dir / "catalogs" / out_name
        if cache_exists(out):
            return ("cached", tile_id, 0)
        ctr = parse_center(tile_id)
        if not ctr:
            return ("no_center", tile_id, 0)
        ra, dec = ctr
        last_err = None
        for attempt in range(1, args.retry + 1):
            try:
                if args.catalog == "ps1":
                    fetch_ps1_neighbourhood(tile_dir, ra, dec, float(args.radius_arcmin),
                                             max_records=int(args.max_records), timeout=float(args.timeout))
                elif args.catalog == "gaia":
                    fetch_gaia_neighbourhood(tile_dir, ra, dec, float(args.radius_arcmin),
                                              max_rows=int(args.max_rows), timeout=float(args.timeout))
                else:  # usnob
                    if not _HAS_USNOB:
                        raise RuntimeError("fetch_usnob_neighbourhood not available; add it to external_fetch_online.py")
                    fetch_usnob_neighbourhood(tile_dir, ra, dec, float(args.radius_arcmin),
                                              max_rows=int(args.max_rows), timeout=float(args.timeout))
                # Quick row count (skip header)
                rows = 0
                try:
                    with out.open('r', encoding='utf-8', errors='ignore') as f:
                        next(f, None)
                        for _ in f:
                            rows += 1
                            if rows > 0: break
                except Exception:
                    rows = -1
                return ("fetched", tile_id, rows)
            except Exception as e:
                last_err = str(e)
                time.sleep(min(2 * attempt, 10))
        return ("failed", tile_id, last_err or "unknown_error")

    tiles = list(iter_tile_dirs(tiles_root))
    counters["tiles_found"] = len(tiles)
    logger.info(f"prewarm start: catalog={args.catalog} tiles_root={tiles_root} tiles_found={len(tiles)} workers={args.workers}")
    write_progress()

    to_run = []
    for td in tiles:
        if args.limit and len(to_run) >= args.limit:
            break
        to_run.append(td)
    counters["tiles_scheduled"] = len(to_run)
    logger.info(f"prewarm scheduled: {len(to_run)} tiles")
    write_progress()

    done = 0
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = {ex.submit(do_one, td): td for td in to_run}
        for fut in as_completed(futs):
            td = futs[fut]
            try:
                status, tile_id, meta = fut.result()
                if status == "cached":
                    counters["tiles_cached_skip"] += 1
                elif status == "fetched":
                    counters["tiles_fetched"] += 1
                    if int(meta) == 0:
                        counters["tiles_zero_rows"] += 1
                    logger.info(f"[OK] {tile_id} {args.catalog}_neighbourhood ready")
                elif status == "no_center":
                    counters["tiles_no_center"] += 1
                    logger.warning(f"[SKIP] {tile_id} no_center")
                else:
                    counters["tiles_failed"] += 1
                    logger.warning(f"[FAIL] {tile_id} err={meta}")
            except StopRequested:
                logger.warning("StopRequested: exiting loop.")
                stop["flag"] = True
                break
            except Exception as e:
                counters["tiles_failed"] += 1
                logger.warning(f"[FAIL] {td.name} unexpected={e}")
            done += 1
            if done % args.progress_every == 0:
                write_progress()
                logger.info(
                    f"progress: done={done}/{len(to_run)} "
                    f"fetched={counters['tiles_fetched']} cached={counters['tiles_cached_skip']} "
                    f"failed={counters['tiles_failed']}"
                )
    write_progress()
    logger.info("prewarm done: " + json.dumps(counters))

if __name__ == "__main__":
    main()
