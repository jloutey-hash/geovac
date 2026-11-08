"""group2 Batch-3 remediation: RED backing-test pins updated to the re-measured
exact-rule values (this session). test_direct_ci.py Li n4/n5 + Be n4 were pinned
to pre-2026-08-29 values (only the He pin was updated at the time), so they have
been RED. test_coupled_one_norm was pinned to 91.65; the exact-rule value is
89.80 (measured this session; matches the table's 80.5 non-identity + identity).

No LaTeX; plain Python edits with per-edit count checks (expect gives the
required occurrence count for a value that appears more than once)."""
from __future__ import annotations
import sys

DCI = "tests/test_direct_ci.py"
CC = "tests/test_coupled_composition.py"

# (path, name, marker_present_means_done, old, new, expected_count)
EDITS = [
    # Li n4: full-precision value appears 3x (docstring xref + assert + msg)
    (DCI, "li4-val", "-7.398718", "-7.395921", "-7.398718", 3),
    (DCI, "li4-doc", "E = -7.39872 Ha (1.06%)",
     "E = -7.39592 Ha (1.10%)", "E = -7.39872 Ha (1.06%)", 1),
    (DCI, "li4-pct", "1.06% above it (matches",
     "1.10% above it (matches", "1.06% above it (matches", 1),

    # Li n5
    (DCI, "li5-val", "-7.400705", "-7.397751", "-7.400705", 3),
    (DCI, "li5-doc", "E = -7.40071 Ha (1.03%)",
     "E = -7.39775 Ha (1.07%)", "E = -7.40071 Ha (1.03%)", 1),
    (DCI, "li5-pct", "1.03% above it (matches",
     "1.07% above it (matches", "1.03% above it (matches", 1),

    # Be n4
    (DCI, "be4-val", "-14.562659", "-14.535460", "-14.562659", 3),
    (DCI, "be4-doc", "E = -14.5627 Ha (0.71%)",
     "E = -14.5355 Ha (0.90%)", "E = -14.5627 Ha (0.71%)", 1),
    (DCI, "be4-pct", "0.71% above it (matches",
     "0.90% above it (matches", "0.71% above it (matches", 1),

    # coupled 1-norm re-pin + docstring
    (CC, "onenorm-pin", "one_norm - 89.80) < 0.3",
     "assert abs(one_norm - 91.65) < 1.0",
     "assert abs(one_norm - 89.80) < 0.3", 1),
    (CC, "onenorm-doc", "= 89.80 Ha\n        (exact-rule 2026-09-13",
     '''"""Coupled 1-norm (all terms incl. identity) ~ 89.0 Ha
        (exact-rule 2026-08-29; was ~85.69 under the retired rule A).
        Non-identity part measured 81.77; identity adds ~9.9."""''',
     '''"""Coupled (CB) 1-norm, all terms incl. identity = 89.80 Ha
        (exact-rule 2026-09-13 recompute; was pinned 91.65, and ~85.69
        under the retired pair-diagonal rule; table 80.5 non-identity)."""''', 1),
]


def main() -> int:
    loaded: dict[str, str] = {}
    applied, skipped, missed = [], [], []
    for path, name, marker, old, new, expect in EDITS:
        if path not in loaded:
            with open(path, encoding="utf-8") as fh:
                loaded[path] = fh.read()
        t = loaded[path]
        if marker in t:
            skipped.append(name); continue
        if t.count(old) != expect:
            missed.append((name, t.count(old), expect)); continue
        loaded[path] = t.replace(old, new); applied.append(name)
    for path, t in loaded.items():
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(t)
    for n in applied: print(f"  ok    {n}")
    for n in skipped: print(f"  skip  {n} (already applied)")
    for n, c, e in missed: print(f"  MISS  {n}: count={c}, expected={e}")
    print(f"applied {len(applied)}, skipped {len(skipped)}, MISSED {len(missed)}")
    return 3 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
