#!/usr/bin/env python3
"""Prove this repository is free of the reference catalogue it is compared to.

The whole claim of this project is that the method needs nothing but public
inputs -- the published plate list, IRSA's full-plate scans, and public
catalogues. That claim is only worth something if it is checkable, so this script
is the check, and anyone can run it.

Updated 2026-10-15. The reference catalogue is no longer unpublished: its
107,875 rows were released on 2026-10-14 in the code repository accompanying
Bruehl et al. (2026), Scientific Reports 16, 31890
(doi:10.1038/s41598-026-66148-2). That changes the reason this gate exists but
not whether it is needed, in two ways, and both make it stricter rather than
looser:

  * Independence is now harder to demonstrate, not easier. While the catalogue
    was unavailable, "we could not have used it" was nearly self-evident. Now
    that anyone can download it, only an explicit, checkable audit separates a
    pipeline that reproduces a result from one that was tuned against it. This
    script is that audit.
  * The naming convention is now a separate requirement of its own. In
    anything outward-facing this project calls the dataset the "VASCO 2025
    sample" -- the owning team's own term for it -- and never a single-letter
    label. See docs/DATASET_NAMING.md. The bare-token rule below enforces that.

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
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

# Tokens that would indicate the private catalogue, its filesystem location, its
# row count, or an analysis that consumed it. Matched case-insensitively.
BANNED = [
    r"vasco_main_list",
    r"vasco-cats",
    r"\bv[-_]csv\b",
    r"tiles_supplement_v",
    r"radec_cache_supplement",
    r"\bv[-_]parity\b",
    r"\bbackfill\b",
    # Prose references to the reference catalogue. The filename patterns above
    # miss these entirely -- "V-in-footprint" in a comment was found by eye, not
    # by this script, which is why they are here now.
    #
    # Widened 2026-10-15 after an audit found four live occurrences that the
    # earlier patterns did not catch: "V's source plates", "No V anywhere",
    # "the V-publication question", and "published V figures". Each was a
    # grammatical form the list did not anticipate, which is the argument for
    # the bare-token rule in BARE_V rather than for yet more hand-written
    # alternatives.
]

# Naming forms. These are about what the dataset is CALLED, so unlike the list
# above they are exempt on a line that is plainly about photometry -- Johnson V
# is a real band and "the V band" is not a leak. Keeping the two lists separate
# is the point: a filename or an artifact name is never excused by context,
# a prose label sometimes is.
BANNED_NAME = [
    r"\bV[- ]?(in[- ]footprint|rows|catalog|catalogue|list|coords|coordinates|positions)\b",
    r"\bV[- ]?(parity|recall|match|matched|free|hostile|neutral|supplement|publication|side|only)\b",
    r"\bthe V\b",
    r"\bno V\b",
    r"\bset V\b",
]

# The catch-all. Any bare "V" or "V's" used as a name, in any grammatical form
# the list above failed to anticipate. This is what makes the gate robust
# against rephrasing instead of against a fixed vocabulary.
BARE_V = r"\bV[\u2019']?s?\b"

# ...except where a capital V is a real astronomical term rather than a label.
# Johnson V is a photometric band and appears legitimately in magnitude talk,
# so a line that is clearly about photometry is not a naming violation. Kept
# deliberately narrow: it exempts the line, so anything added to it weakens the
# gate, and nothing belongs here that is not a standard term of art.
EXEMPT = r"\bV[- ]?(band|mag|magnitude|filter)\b|\bJohnson\b|\bB[- ]?V\b"

# Tokens that are PUBLIC, but only safe with their source attached.
#
# The transient-catalogue row count was treated as a private marker until
# 2026-08-16. It is not: it is stated in Doherty, "Independent Replication of
# Nuclear Test-Transient Correlations and Earth Shadow Deficit in POSS-I
# Photographic Plates", arXiv:2604.00056, alongside the 635-plate figure, and
# appears in the wider literature. Blanket-banning it forced circumlocutions
# ("roughly five times larger") that made published arithmetic harder to check,
# which is the opposite of what this audit is for.
#
# It stays conditional rather than simply allowed. Writing the number is fine
# when the citation travels with it on the same line; writing it bare is what
# this audit should still catch, because a bare row count reads as a property of
# a catalogue we hold rather than a figure someone else published.
CONDITIONAL = [
    (
        r"107[,_]?875",
        r"2604\.00056|s41598-026-66148|10\.1038|arXiv|PASP|Doherty|Bruehl|published",
        "public via doi:10.1038/s41598-026-66148-2 (and arXiv:2604.00056) -- "
        "keep the citation on the same line",
    ),
]

# Paths that are allowed to mention the tokens, because their job is to name them.
ALLOW = {"tools/audit_independence.py", "docs/DATASET_NAMING.md"}


_NAME_RE = re.compile("|".join(BANNED_NAME), re.I)
_BARE_V_RE = re.compile(BARE_V)          # case-sensitive: only a capital V is a label
_EXEMPT_RE = re.compile(EXEMPT, re.I)


def _bare_v_hit(line: str) -> bool:
    """True for a naming violation on this line.

    Covers both the spelled-out forms in BANNED_NAME and any bare "V"/"V's"
    used as a label. The bare check is case-sensitive on purpose: lowercase "v"
    turns up constantly in version strings and abbreviations and is never the
    label; only the capital is.

    Either way the line is excused if it is plainly about photometry, because
    Johnson V is a band and not this catalogue.
    """
    if not (_NAME_RE.search(line) or _BARE_V_RE.search(line)):
        return False
    return not _EXEMPT_RE.search(line)


def _conditional_hit(line: str) -> str | None:
    """Return a reason if a conditional token appears without its context."""
    for tok, ctx, why in CONDITIONAL:
        if re.search(tok, line, re.I) and not re.search(ctx, line, re.I):
            return why
    return None

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
                continue
            if _bare_v_hit(line):
                hits.append((rel, i, f"[bare name] {line.strip()[:100]}  <- use "
                                     f"\"VASCO 2025 sample\"; see docs/DATASET_NAMING.md"))
                continue
            why = _conditional_hit(line)
            if why:
                hits.append((rel, i, f"[uncited] {line.strip()[:100]}  <- {why}"))
    return hits


BASELINE_PATH = Path(__file__).with_name("audit_independence_baseline.json")


def _key(path: str, line: str) -> str:
    """Stable identifier for one historical hit.

    Deliberately a hash rather than the text. A baseline that quoted the
    offending lines would reintroduce into this repository exactly the tokens
    the audit exists to keep out -- the file would be its own violation, and it
    would have to be allow-listed to avoid tripping the gate, which is how
    allow-lists quietly grow until they cover everything. The hash pins the hit
    precisely without reproducing it: change the line and the key changes, so a
    baselined entry cannot be stretched to cover something new.
    """
    return hashlib.sha256(f"{path}\n{line.strip()}".encode()).hexdigest()[:16]


def load_baseline() -> dict:
    """Known pre-existing history hits, accepted once and tracked thereafter.

    Git history cannot be edited in place; removing a line from it means
    rewriting every commit after it and force-pushing, which breaks every
    existing clone. For hits that are comments rather than data, that cure is
    worse than the disease -- but leaving the gate permanently red is worse
    still, because a check that always fails stops being read, and then it
    protects nothing. So known hits are recorded here, with a reason, and the
    gate fails only on hits that are not in the record.

    This weakens the gate by exactly the number of entries in the file, which
    is why the file is small, reviewed, and prints its count on every run.
    """
    if not BASELINE_PATH.exists():
        return {}
    try:
        data = json.loads(BASELINE_PATH.read_text())
    except Exception as e:
        print(f"[WARN] baseline unreadable ({e}); treating every hit as new")
        return {}
    return {e["key"]: e for e in data.get("entries", [])}


def scan_history(root: Path) -> list[tuple[str, str, str]]:
    """Grep every blob ever committed. Empty history is trivially clean.

    Returns (key, path, description) so the caller can match against the
    baseline.
    """
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
    hits: list[tuple[str, str, str]] = []
    cur = ""
    for ln in out.stdout.splitlines():
        # Track which file each hunk belongs to, so ALLOW applies to history the
        # same way it applies to the tree. Without this the audit flags its own
        # banned-token list the moment it is committed, and the only way to get a
        # green run is to stop reading the output -- which defeats the gate.
        if ln.startswith("+++ b/"):
            cur = ln[6:].strip()
            continue
        if ln.startswith("+") and cur not in ALLOW:
            if pat.search(ln):
                desc = ln.strip()[:100]
            elif _bare_v_hit(ln):
                desc = f"[bare name] {ln.strip()[:90]}  <- use \"VASCO 2025 sample\""
            else:
                why = _conditional_hit(ln)
                if not why:
                    continue
                desc = f"[uncited] {ln.strip()[:90]}  <- {why}"
            hits.append((_key(cur, ln), cur, desc))
    return hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--tree-only", action="store_true")
    ap.add_argument("--update-baseline", action="store_true",
                    help="record every current history hit as accepted. Only "
                         "for a deliberate, reviewed decision -- it makes the "
                         "gate green by lowering the bar, not by fixing "
                         "anything.")
    args = ap.parse_args()
    root = Path(args.root)

    tree = scan_tree(root)
    hist = [] if args.tree_only else scan_history(root)

    if args.update_baseline:
        if args.tree_only:
            print("[FAIL] --update-baseline needs the history scan")
            return 2
        payload = {
            "note": "Known history hits, accepted on review. Git history cannot "
                    "be edited in place, so these are recorded rather than "
                    "removed. Keys are sha256(path + line)[:16]; the text is "
                    "deliberately not stored, so this file never reproduces "
                    "what it describes. Every entry weakens the gate by one "
                    "line -- keep the list short and justify each addition.",
            "recorded": "2026-10-15",
            # Deduplicated on key. One line touched by several commits
            # produces the same key each time, and listing it twice would
            # overstate the baseline's size while understating how much a
            # single removal reopens.
            "entries": [
                {"key": k, "path": pth, "why": "comment or prose, not data"}
                for k, pth in sorted({(k, pth) for k, pth, _ in hist})
            ],
        }
        BASELINE_PATH.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"[OK]   baseline written: {len(payload['entries'])} entr"
              f"{'y' if len(payload['entries']) == 1 else 'ies'} "
              f"-> {BASELINE_PATH.name}")
        return 0

    baseline = load_baseline()
    # Count distinct excused lines, not hit occurrences: the number that matters
    # is how many lines the gate is no longer checking.
    known = sorted({h[0] for h in hist if h[0] in baseline})
    new_hits = [h for h in hist if h[0] not in baseline]
    seen = {h[0] for h in hist}
    stale = [k for k in baseline if k not in seen]

    if tree:
        print(f"[FAIL] {len(tree)} working-tree hit(s):")
        for f, i, l in tree:
            print(f"   {f}:{i}: {l}")
    else:
        print("[OK]   working tree clean")

    if not args.tree_only:
        if new_hits:
            print(f"\n[FAIL] {len(new_hits)} NEW history hit(s) -- not in the "
                  f"baseline. These cannot be fixed by editing files; the "
                  f"history itself carries them:")
            for _, pth, desc in new_hits[:40]:
                print(f"   {pth}: {desc}")
        elif known:
            print(f"\n[OK]   git history: no new hits "
                  f"({len(known)} known, baselined -- see {BASELINE_PATH.name})")
        else:
            print("[OK]   git history clean")

        if known and new_hits:
            print(f"       ({len(known)} further hit(s) are baselined and not "
                  f"counted above)")
        if stale:
            print(f"\n[WARN] {len(stale)} baseline entr"
                  f"{'y' if len(stale) == 1 else 'ies'} no longer "
                  f"{'matches' if len(stale) == 1 else 'match'} anything in "
                  f"history -- the history was rewritten, or the entry was "
                  f"wrong. Re-run with --update-baseline to prune.")

    ok = not tree and not new_hits
    print("\nINDEPENDENCE AUDIT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
