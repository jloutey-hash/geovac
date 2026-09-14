"""DELTA #2 -- close the seven loci the WIDENED C16 scope exposed.

Widening the six new p60 entries' `files` lists (they had shipped narrower
than the OLDER p60 entries sitting next to them) took C16 from PASS to FAIL
with 7 live loci.  That is the gate working: the narrow scope was the reason
the walls register's "Proposition D" and the single-value floor in
claims_register/INDEX were invisible on 2026-09-12 AND 2026-09-13.

Two classes among the seven:

  (a) LEGITIMATE WITHDRAWAL RECORDS that quote the retired wording in order to
      retire it.  These are correct text; they simply lack the standardized
      `[retracted YYYY-MM-DD: key]` token the gate keys its exemption on.  The
      gate's own marker note says 46 of 69 entries still rely on hand-authored
      exemption vocabulary -- the class that produced three false-clean entries
      on 2026-09-03.  Token added.

  (b) ONE GENUINELY STALE locus: `paper_60.done.md` C8#15 uses "Proposition D"
      as a live name with no re-attribution anywhere near it, while criterion
      22 forty lines later carries the full re-attribution.  Same file,
      inconsistent.

Idempotent.
"""
from __future__ import annotations

import sys

DOD = "docs/qa/paper_60.done.md"
WAL = "docs/walls/register.md"

EDITS = [
    # (a) withdrawal records -- add the standardized token
    (DOD, "frames-withdrawal-token", "[retracted 2026-09-12: p60-frames-completeness]",
     '**The 2026-09-11 frames reading — "overcompleteness is the price of one-centre\n'
     '    completeness" — is WITHDRAWN;** re-asserting it = MATERIAL.',
     '**The 2026-09-11 frames reading — "overcompleteness is the price of one-centre\n'
     '    completeness" — is WITHDRAWN [retracted 2026-09-12: p60-frames-completeness];**\n'
     '    re-asserting it = MATERIAL.'),

    (DOD, "removability-withdrawal-token", "[retracted 2026-09-12: p60-removability-corollary]",
     '**The removability\n    corollary is WITHDRAWN**',
     '**The removability\n    corollary is WITHDRAWN [retracted 2026-09-12: p60-removability-corollary]**'),

    (DOD, "crit22-propd-token", "[retracted 2026-09-12: p60-prop-d-as-new]",
     'Presenting Proposition D as a new proposition = MATERIAL. Backing',
     'Presenting it as a new proposition of ours = MATERIAL\n'
     '    [retracted 2026-09-12: p60-prop-d-as-new — the LABEL is retired; the\n'
     '    `l`-vs-`m` application is what the paper claims]. Backing'),

    (DOD, "crit22-propd-name", "the block-diagonal congruence result — a",
     '22. **The `l`-selection loss is NOT a conditioning effect [SYMBOLIC].** Proposition D — a',
     '22. **The `l`-selection loss is NOT a conditioning effect [SYMBOLIC].** The\n'
     '    block-diagonal congruence result — a'),

    (WAL, "walls-propd-token", "[retracted 2026-09-12: p60-prop-d-as-new]",
     'the label "Proposition D" is retired and only the `l`-vs-`m` application is ours)',
     'the label "Proposition D" is retired [retracted 2026-09-12: p60-prop-d-as-new] '
     'and only the `l`-vs-`m` application is ours)'),

    # (b) the genuinely stale locus
    (DOD, "crit15-propd-stale", "the block-diagonal congruence result (Löwdin / Slater-Koster 1954)",
     'it does NOT recover `l`-selection — Proposition D is untouched.',
     'it does NOT recover `l`-selection — the block-diagonal congruence result\n'
     '    (Löwdin / Slater-Koster 1954) is untouched.'),
]


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name)
            continue
        if t.count(old) != 1:
            missed.append((name, t.count(old)))
            continue
        loaded[path] = t.replace(old, new)
        applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied:
        print(f"  ok    {n}")
    for n in skipped:
        print(f"  skip  {n} (already applied)")
    for n, c in missed:
        print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
