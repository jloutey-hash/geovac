"""group2 Batch-4 remediation, docs: the stale-value loci the completeness-critic
and claims-P58 reviewer found outside the papers.

claim-matrix:289 (FCI-A He hybrid) -- still "0.35% ... FALSE-POSITIVE/NO-TEST"
  though Batch 3 re-measured (0.26%). Update the number + status note.
claim-matrix:483 (Paper 58 QFD) -- claim column says "84-digit" while the row's
  OWN note and the paper (L1296) say 60. Fix the internal inconsistency.
group2.done.md C8 literals -- "Be 0.90% / Li 1.07%" are the superseded pre-ERI-fix
  values ("a .done.md ratifying a retired value"). Update to 0.71% / 1.03%.

No LaTeX backslashes in the edited spans; plain Python edits, idempotent.
"""
from __future__ import annotations
import sys

MTX = "docs/claim_test_matrix.md"
DOD = "docs/qa/group2.done.md"

EDITS = [
    (MTX, "mtx-289-fcia", "He 0.26% (grid hybrid)",
     "| FCI-A | He 0.35% (grid hybrid) / 0.19% (graph-native) | vacuous (E<0) / wrong h1 config | **FALSE-POSITIVE/NO-TEST** | backfill hybrid-config recompute test |",
     "| FCI-A | He 0.26% (grid hybrid) / 0.19% (graph-native) | reproduces this session (paper-convention hybrid direct CI, 0.26%); no committed hybrid-config pin | **NO-TEST (reproduces)** | 0.35%→0.26% (exact-rule re-measure 2026-09-13); backfill hybrid-config recompute test |"),

    (MTX, "mtx-483-60digit", "H₂ 60-digit certified FCI energy",
     "H₂ 84-digit certified FCI energy",
     "H₂ 60-digit certified FCI energy"),

    (DOD, "dod-c8-fcia", "graph-native CI He 0.19% / Be 0.71% / Li 1.03%",
     "graph-native CI He 0.19% / Be 0.90% / Li 1.07%",
     "graph-native CI He 0.19% / Be 0.71% / Li 1.03%"),
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
            skipped.append(name); continue
        if t.count(old) != 1:
            missed.append((name, t.count(old))); continue
        loaded[path] = t.replace(old, new); applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied: print(f"  ok    {n}")
    for n in skipped: print(f"  skip  {n} (already applied)")
    for n, c in missed: print(f"  MISS  {n}: anchor count={c}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
