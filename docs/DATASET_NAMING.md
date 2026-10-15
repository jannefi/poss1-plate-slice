# How this repository names the catalogues it compares against

A short convention, so that the same dataset is called the same thing
everywhere and so that two different datasets are never conflated.

## The rule

**The 107,875-row transient catalogue produced by the VASCO team's detection
pipeline is called the "VASCO 2025 sample".** That is the owning team's own
term for it, and this repository uses it in documentation, code comments,
figure labels and commit messages alike.

**Single-letter labels are not used in this repository.** Letter notation for
these catalogues circulates in the wider literature and discussion, and it is
perfectly legitimate there — this is a house style, not a criticism of it.
Inside this repository it is banned for two practical reasons:

1. **A single letter is not self-describing.** A reader who arrives at one file
   cannot tell what it refers to, and the letters are not used consistently
   across the sources that employ them.
2. **It hides a distinction that matters.** There is more than one catalogue in
   this family and they are not interchangeable — see below. A bare letter
   makes substituting one for the other invisible; a name does not.

When a different catalogue of theirs is meant, name that one specifically
rather than reaching for a letter or a row count.

## The distinction the convention protects

| | rows | where it appears |
|---|---:|---|
| the **VASCO 2025 sample** | 107,875 | the dataset scored in Bruehl et al. (2026) |
| its northern sibling | 106,339 | the PASP 2025 shadow-test sample |
| the published candidate list | 5,399 | Solano et al. (2022), public via SVO |

These are three different sets. A figure quoted against one of them is wrong
against the others, and the row counts are close enough that a mix-up survives
casual reading. Writing the name out is what prevents it.

## Figures and images

The same rule applies to anything rendered into an image — axis labels, legend
entries, annotations, screenshots of a terminal. A label baked into a figure is
easy to miss, because proofreading attention goes to the prose around it and
not to the pixels. Generic labels ("reference catalogue", "target") are
preferred in figures even where a full name would fit.

## Enforcement

`tools/audit_independence.py` enforces this automatically over both the working
tree and the full git history, and exits non-zero on a violation, so it works
as a pre-commit hook and in CI:

```bash
python3 tools/audit_independence.py
```

It distinguishes two kinds of hit. **Hard** tokens — filesystem paths, artifact
filenames — are never excused. **Naming** tokens are excused on a line that is
plainly about photometry, because Johnson V is a real band and "the V band" is
not a naming violation. That exemption is deliberately narrow: it excuses the
whole line, so anything added to it weakens the gate.

This file is on the audit's allow-list, because its job is to name the thing it
bans. It is the only documentation file that is.

### The history baseline

Git history cannot be edited in place. Removing a line from it means rewriting
every commit after it and force-pushing, which breaks every existing clone —
a cure worse than the disease when the line in question is a comment rather
than data. But leaving the gate permanently red is worse still: a check that
always fails stops being read, and then it protects nothing.

So the audit carries `tools/audit_independence_baseline.json`: **eight
pre-existing history hits, reviewed on 2026-10-15 and accepted**. All eight are
comments or prose in tracked source — the full list, with a reason for each, is
in that file. **No row of any reference catalogue has ever been committed to
this repository**, which was verified separately by scanning the history for
path and artifact-filename tokens; the only such hit is a comment in
`scripts/gcp/push_to_gcp.sh` asserting that a file is *not* copied to a VM.

Three properties keep the baseline from becoming a loophole:

- **It stores hashes, not text.** A baseline that quoted the offending lines
  would reintroduce the very tokens the audit exists to exclude, and would then
  need allow-listing itself. The key is `sha256(path + line)[:16]`, so change
  the line and the key changes — an entry cannot be stretched to cover
  something new.
- **It never applies to the working tree.** Anything you can fix by editing a
  file, you must; the baseline covers only what editing cannot reach.
- **It is visible and it shrinks.** Every run prints how many lines are
  excused, and an entry that stops matching anything is reported as stale
  rather than kept silently.

`--update-baseline` re-records the current hits. It makes the gate green by
lowering the bar rather than by fixing anything, so it is for a deliberate,
reviewed decision and not for getting a red run to pass.

## How to validate

```bash
printf 'more rows than V\n' > docs/_gate_test.md && python3 tools/audit_independence.py --tree-only; rm docs/_gate_test.md
```

That must report a `[bare name]` hit and exit non-zero. A gate that has never
been seen to fail has not been tested.
