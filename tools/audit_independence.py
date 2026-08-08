#!/usr/bin/env python3
"""Prove this repository is free of the private reference catalogue.

The whole claim of this project is that the method needs nothing but public
inputs -- the published plate list, IRSA's full-plate scans, and public
catalogues. That claim is only worth something if it is checkable, so this script
is the check, and anyone can run it.

It scans two places, because they fail differently:

  * the working tree -- what you get when you clone
  * the full git history (`git log -p`) -- what a deleted file leaves behind

A file removed in a later commit is still in the pack, still fetched by every
clone, and still findable. That is why this repository was built by copying into
a fresh `git init` rather than by pruning a fork: there is no history to leak.
This script keeps it that way.

Exit status is 0 only if both scans are clean, so it works as a pre-commit hook
and in CI.

How to validate
---------------
    python3 tools/audit_independence.py            # tree + history
    python3 tools/audit_independence.py --tree-only

Deliberately verify it can fail: add a file containing one of the banned tokens
and confirm a non-zero exit. A gate that has never failed has not been tested.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

# Tokens that would indicate the private catalogue, its filesystem location, its
# row count, or an analysis that consumed it. Matched case-insensitively.
BANNED = [
    r"vasco_main_list",
    r"vasco-cats",
    r"107[,_]?875",
    r"\bv[-_]csv\b",
    r"tiles_supplement_v",
    r"radec_cache_supplement",
    r"\bv[-_]parity\b",
    r"\bbackfill\b",
    # Prose references to the private catalogue. The filename patterns above miss
    # these entirely -- "V-in-footprint" in a comment was found by eye, not by
    # this script, which is why they are here now.
    r"\bV[- ]?(in[- ]footprint|rows|catalog|catalogue|list|coords|coordinates|positions)\b",
    r"\bthe V\b",
]

# Paths that are allowed to mention the tokens, because their job is to name them.
ALLOW = {"tools/audit_independence.py"}

SKIP_DIRS = {".git", "__pycache__", ".venv", "venv", "node_modules"}
SKIP_SUFFIX = {".pyc", ".gz", ".fits", ".parquet", ".png", ".jpg", ".pdf"}


def _candidate_files(root: Path):
    """Files a clone would actually contain.

    Scanning the raw working tree is the wrong question: it sweeps in gitignored
    scratch -- test subsets, run outputs, caches -- which are never distributed
    and which will trip this check constantly. A gate that cries wolf gets
    ignored, and then it protects nothing. Ask git what is tracked or
    tracked-able instead, and fall back to a filesystem walk outside a repo.
    """
    try:
        out = subprocess.run(
            ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
            cwd=root, capture_output=True, text=True, check=True)
        rels = [l for l in out.stdout.splitlines() if l.strip()]
        if rels:
            return [(root / r, r) for r in rels]
    except Exception:
        pass
    return [(p, str(p.relative_to(root))) for p in root.rglob("*")]


def scan_tree(root: Path) -> list[tuple[str, int, str]]:
    pat = re.compile("|".join(BANNED), re.I)
    hits = []
    for p, rel in _candidate_files(root):
        if not p.is_file() or p.suffix in SKIP_SUFFIX:
            continue
        if SKIP_DIRS & set(p.parts):
            continue
        if rel in ALLOW:
            continue
        try:
            text = p.read_text(errors="ignore")
        except Exception:
            continue
        for i, line in enumerate(text.splitlines(), 1):
            if pat.search(line):
                hits.append((rel, i, line.strip()[:120]))
    return hits


def scan_history(root: Path) -> list[str]:
    """Grep every blob ever committed. Empty history is trivially clean."""
    try:
        n = subprocess.run(["git", "rev-list", "--count", "HEAD"], cwd=root,
                           capture_output=True, text=True)
        if n.returncode != 0 or n.stdout.strip() in ("", "0"):
            return []
        out = subprocess.run(
            ["git", "log", "-p", "--all", "--no-color"],
            cwd=root, capture_output=True, text=True, errors="ignore")
    except Exception as e:
        return [f"history scan failed: {e}"]
    pat = re.compile("|".join(BANNED), re.I)
    hits, cur = [], ""
    for ln in out.stdout.splitlines():
        # Track which file each hunk belongs to, so ALLOW applies to history the
        # same way it applies to the tree. Without this the audit flags its own
        # banned-token list the moment it is committed, and the only way to get a
        # green run is to stop reading the output -- which defeats the gate.
        if ln.startswith("+++ b/"):
            cur = ln[6:].strip()
            continue
        if ln.startswith("+") and cur not in ALLOW and pat.search(ln):
            hits.append(f"{cur}: {ln.strip()[:100]}")
    return hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--tree-only", action="store_true")
    args = ap.parse_args()
    root = Path(args.root)

    tree = scan_tree(root)
    hist = [] if args.tree_only else scan_history(root)

    if tree:
        print(f"[FAIL] {len(tree)} working-tree hit(s):")
        for f, i, l in tree:
            print(f"   {f}:{i}: {l}")
    else:
        print("[OK]   working tree clean")

    if not args.tree_only:
        if hist:
            print(f"\n[FAIL] {len(hist)} history hit(s) -- these cannot be fixed by "
                  f"editing files; the history itself carries them:")
            for l in hist[:40]:
                print(f"   {l}")
        else:
            print("[OK]   git history clean")

    ok = not tree and not hist
    print("\nINDEPENDENCE AUDIT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
