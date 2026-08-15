#!/usr/bin/env python3
"""PTF veto stage with coverage measured, not assumed.

PTF is a modern-epoch survey, so its semantics are the OPPOSITE of the
SuperCOSMOS stage: a counterpart means the source is still there, hence not a
vanishing candidate, hence DROPPED. MNRAS 2022 applies a quality gate with it,
which this reproduces: a match counts only if COALESCE(ngoodobs,0) > 0.

WHY A COVERAGE COLUMN

PTF is not an all-sky survey. Its footprint over 2009-2012 is wide but uneven,
so "no PTF counterpart" has two very different meanings: PTF looked and saw
nothing, or PTF never looked. Conflating them does not create false removals --
for a drop-matches veto an uncovered row is KEPT, which is the conservative
direction -- but it does make the removal rate uninterpretable, because the
denominator silently includes sky the survey could never have tested.

So each chunk asks one query for two numbers:

    n5_ngood   PTF objects within 5" passing the ngoodobs gate  -> the veto
    n60_any    PTF objects of any kind within 60"               -> coverage

A row with n60_any = 0 is almost certainly outside PTF coverage. The ledger
reports the removal rate both over all rows and over covered rows only; the
second is the scientifically meaningful one and the first is what a naive run
would have quoted.

RESUMABLE, chunked, flags-only. Nothing is dropped here -- use
tools/shrink_stage.py --keep unmatched to apply it.

Usage:
    python3 tools/ptf_stage_coverage.py \\
      --in-csv <run-dir>/stages/stage_S1_SCOS_upload.csv \\
      --out-dir <run-dir>/stages/ptf --chunk-size 2000

How to validate: --limit-chunks 1 first; n5_ngood must be <= n60_any for every
row (asserted), and the covered fraction should be high but not exactly 1.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd

IRSA_SYNC = "https://irsa.ipac.caltech.edu/TAP/sync"
TABLE = "ptf_objects"


# IRSA's ADQL will not accept a conditional aggregate -- SUM(CASE WHEN
# CONTAINS(...) ...) returns "Invalid or unsupported ADQL query string". Both
# halves are fine on their own, so each chunk costs two simple queries instead
# of one clever one. Verified against the service before writing this.
def adql_veto(table: str, r_deg: float) -> str:
    return ("SELECT DISTINCT u.row_id AS row_id "
            f"FROM TAP_UPLOAD.t1 AS u, {table} AS p "
            "WHERE CONTAINS(POINT('ICRS',p.ra,p.dec),"
            f"CIRCLE('ICRS',u.ra,u.dec,{r_deg:.10f}))=1 "
            "AND COALESCE(p.ngoodobs,0)>0")


def adql_coverage(table: str, r_deg: float) -> str:
    return ("SELECT u.row_id AS row_id, COUNT(*) AS n60_any "
            f"FROM TAP_UPLOAD.t1 AS u, {table} AS p "
            "WHERE CONTAINS(POINT('ICRS',p.ra,p.dec),"
            f"CIRCLE('ICRS',u.ra,u.dec,{r_deg:.10f}))=1 "
            "GROUP BY u.row_id")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-csv", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--radius-arcsec", type=float, default=5.0)
    ap.add_argument("--coverage-arcsec", type=float, default=60.0)
    ap.add_argument("--chunk-size", type=int, default=2000)
    ap.add_argument("--limit-chunks", type=int, default=0)
    ap.add_argument("--table", default=TABLE)
    ap.add_argument("--url", default=IRSA_SYNC)
    ap.add_argument("--sleep", type=float, default=2.0)
    args = ap.parse_args(argv)

    stilts, curl = shutil.which("stilts"), shutil.which("curl")
    if not stilts or not curl:
        raise SystemExit("[FAIL] need stilts and curl on PATH")
    out = Path(args.out_dir)
    (out / "parts").mkdir(parents=True, exist_ok=True)

    d = pd.read_csv(args.in_csv)
    n = len(d)
    queries = {"veto": adql_veto(args.table, args.radius_arcsec / 3600.0),
               "cov": adql_coverage(args.table, args.coverage_arcsec / 3600.0)}
    nch = (n + args.chunk_size - 1) // args.chunk_size
    todo = nch if not args.limit_chunks else min(nch, args.limit_chunks)
    print(f"[PTF] {n:,} rows, {nch} chunks of {args.chunk_size}, running {todo}")
    print(f"[PTF] veto {args.radius_arcsec}\" (ngoodobs>0), "
          f"coverage probe {args.coverage_arcsec}\"")
    for k, v in queries.items():
        print(f"[PTF] ADQL[{k}]: {v}")
    print()

    def run_query(kind: str, up_vot: Path, dest: Path, tag: str) -> bool:
        r = subprocess.run(
            [curl, "-sS", "--max-time", "900", "--retry", "4",
             "--retry-delay", "10", "-o", str(dest),
             "-F", "REQUEST=doQuery", "-F", "LANG=ADQL", "-F", "FORMAT=csv",
             "-F", f"QUERY={queries[kind]}",
             "-F", "UPLOAD=t1,param:uploadfile",
             "-F", f"uploadfile=@{up_vot}", args.url],
            text=True, capture_output=True)
        if r.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
            print(f"[PTF] {tag} {kind} FAILED (curl {r.returncode}) {r.stderr[:200]}")
            dest.unlink(missing_ok=True)
            return False
        head = dest.read_text()[:400]
        if head.lstrip().startswith("<") or "QUERY_STATUS" in head:
            print(f"[PTF] {tag} {kind} error document:\n{head[:300]}")
            dest.unlink(missing_ok=True)
            return False
        return True

    t0 = time.time()
    for i in range(todo):
        tag = f"part {i+1}/{todo}"
        pv = out / "parts" / f"veto_{i+1:04d}.csv"
        pc = out / "parts" / f"cov_{i+1:04d}.csv"
        if pv.exists() and pc.exists() and pv.stat().st_size and pc.stat().st_size:
            print(f"[PTF] {tag}: cached")
            continue
        sub = d.iloc[i * args.chunk_size:(i + 1) * args.chunk_size]
        up_csv = out / "parts" / f"_up_{i+1:04d}.csv"
        up_vot = out / "parts" / f"_up_{i+1:04d}.vot"
        sub.rename(columns={"src_id": "row_id"})[["row_id", "ra", "dec"]] \
           .to_csv(up_csv, index=False)
        subprocess.run([stilts, "tcopy", f"in={up_csv}", "ifmt=csv",
                        f"out={up_vot}", "ofmt=votable"], check=True)
        for kind, dest in (("veto", pv), ("cov", pc)):
            if dest.exists() and dest.stat().st_size:
                continue
            if not run_query(kind, up_vot, dest, tag):
                print("[PTF] re-run to resume from this chunk")
                return 1
            time.sleep(args.sleep)
        up_csv.unlink(missing_ok=True)
        up_vot.unlink(missing_ok=True)
        el = time.time() - t0
        print(f"[PTF] {tag}  {len(sub)} rows  {el:.0f}s  "
              f"~{el/(i+1)*(todo-i-1):.0f}s left", flush=True)

    def gather(pattern: str, cols: list[str]) -> pd.DataFrame:
        fs = sorted((out / "parts").glob(pattern))
        got = []
        for p in fs:
            if not p.stat().st_size:
                continue
            try:
                got.append(pd.read_csv(p))
            except pd.errors.EmptyDataError:
                continue
        return pd.concat(got, ignore_index=True) if got else pd.DataFrame(columns=cols)

    vet = gather("veto_*.csv", ["row_id"])
    cov = gather("cov_*.csv", ["row_id", "n60_any"])
    done = min(todo * args.chunk_size, n)
    f = d.iloc[:done][["src_id", "ra", "dec"]].copy()
    f["has_ptf_match_ngood"] = f.src_id.isin(set(vet.row_id)) if len(vet) else False
    f = f.merge(cov.rename(columns={"row_id": "src_id"}), on="src_id", how="left")
    f["n60_any"] = f.n60_any.fillna(0).astype(int)
    f["covered"] = f.n60_any > 0
    # A 5" match inside a 60" circle must also appear in the coverage count.
    bad = int((f.has_ptf_match_ngood & ~f.covered).sum())
    if bad:
        print(f"[WARN] {bad} rows matched at {args.radius_arcsec}\" but show no "
              f"neighbour at {args.coverage_arcsec}\" -- inconsistent, investigate")
    f["n5_ngood"] = f.has_ptf_match_ngood.astype(int)
    f.to_csv(out / "ptf_flags.csv", index=False)
    cov = f[f.covered]
    res = {"rows": int(len(f)), "covered": int(f.covered.sum()),
           "covered_pct": round(100 * f.covered.mean(), 3),
           "matched": int(f.has_ptf_match_ngood.sum()),
           "matched_pct_all": round(100 * f.has_ptf_match_ngood.mean(), 3),
           "matched_pct_covered": round(100 * cov.has_ptf_match_ngood.mean(), 3)
           if len(cov) else None,
           "radius_arcsec": args.radius_arcsec,
           "coverage_arcsec": args.coverage_arcsec, "table": args.table}
    (out / "ptf_ledger.json").write_text(json.dumps(res, indent=2) + "\n")

    print(f"\n[PTF] {len(f):,} rows queried")
    print(f"  inside PTF coverage      : {res['covered']:,} ({res['covered_pct']:.1f}%)")
    print(f"  PTF match (would DROP)   : {res['matched']:,} "
          f"({res['matched_pct_all']:.2f}% of all)")
    print(f"  ... of COVERED rows only : {res['matched_pct_covered']}%   "
          f"<- the interpretable rate")
    print(f"\n[OUT] {out/'ptf_flags.csv'}\n[NOTE] flags only; apply with "
          f"tools/shrink_stage.py --keep unmatched")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
