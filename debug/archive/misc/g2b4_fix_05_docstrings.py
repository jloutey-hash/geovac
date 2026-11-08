"""group2 Batch-4: two stale test docstrings the code-P58 reviewer flagged (the
asserts are correct; only the docstrings are stale).
  test_paper58_qfd L92: "84-digit certified value" -> 60-digit (matches the
    assert 1e-40 and the module's own corrected note at L11-13).
  test_paper58_decompactification_front header L27: "t_c = 2.7456" -> 2.664
    (the test body L166/L170 already computes 2.664 and REJECTS 2.7456).
Leave the historical-note "84 digits" at qfd L13 and the rejected-value 2.7456
at decompactification L166/L172 (those are correct as-is).
No backslashes; plain Python edits, idempotent."""
from __future__ import annotations
import sys

QFD = "tests/test_paper58_qfd.py"
DEC = "tests/test_paper58_decompactification_front.py"

EDITS = [
    (QFD, "qfd-doc-60", "at R=1.4 it reproduces the 60-digit certified value.",
     "at R=1.4 it reproduces the 84-digit certified value.",
     "at R=1.4 it reproduces the 60-digit certified value."),
    (DEC, "dec-tc-header", "t < t_c = 2.664, the root of",
     "t < t_c = 2.7456, the root of",
     "t < t_c = 2.664, the root of"),
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
