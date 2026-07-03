#!/usr/bin/env python3
r"""
Deterministic paper<->test reference integrity check  --  the /qa C13 criterion.

Verifies that the test references the papers cite inline stay honest against the
claim->test record (docs/claim_test_matrix.md) and the actual test suite. Two
layers:

  GATE (deterministic, FAILs):
    every inline `tests/test_*.py` reference in a gated paper resolves to a
    REAL, LIVE test file under tests/ (not renamed, not typo'd, not archived).
    A stale reference -- a paper citing a test that was renamed or archived --
    is the real failure mode this catches (the C11/C12 analog for test refs).

  ADVISORY (reported, does NOT fail):
    * inline refs that resolve only under tests/_archive/ (a live claim cited
      to a dead-end / superseded test -- a smell);
    * matrix coverage: claim->test rows in claim_test_matrix.md whose test is
      not cited inline ANYWHERE in the papers (a deliberate gap is allowed --
      not every matrix claim has a clean prose anchor -- so this is advisory).

Why existence is the gate, not "every ref in the matrix": the matrix is a
curated LOAD-BEARING subset; papers legitimately cite more tests than the
matrix tracks (Paper 32 cites ~20, the matrix lists ~5), so requiring every
inline ref to appear in the matrix would over-fail. Existence has no such
false positive.

Scope: --gate <branch-folder> gates a branch (default = the trunk target);
--gate-corpus gates everything.

Exit 0 = no MISSING ref in the gated scope. Exit 1 = >=1 MISSING.

Usage:
  python debug/qa/check_paper_test_refs.py
  python debug/qa/check_paper_test_refs.py --gate group5_qed_gauge
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
PAPERS = ROOT / "papers"
TESTS = ROOT / "tests"
ARCHIVE = TESTS / "_archive"
MATRIX = ROOT / "docs" / "claim_test_matrix.md"

TRUNK = {
    "Paper_0_Geometric_Packing.tex",
    "paper_1_spectrum.tex",
    "Paper_7_Dimensionless_Vacuum.tex",
    "paper_32_spectral_triple.tex",
    "paper_38_su2_propinquity_convergence.tex",
    "group3_foundations_synthesis.tex",
}

# an inline reference: either tests/-prefixed, or bare but .py-suffixed; the
# name may end in a glob (e.g. test_qed_*.py = a family of tests). Underscores
# are LaTeX-escaped (\_) in the .tex.
REF = re.compile(
    r"tests/test(?:\\?_[A-Za-z0-9*]+)+(?:\.py)?"      # tests/-prefixed
    r"|test(?:\\?_[A-Za-z0-9*]+)+\.py"                 # bare but .py-suffixed
)
MATRIX_TEST = re.compile(r"\btest_[A-Za-z0-9_]+\b")
# bare FUNCTION-name citation inside \texttt{...}, no .py suffix and no tests/
# prefix (the group5 1st-cert P33 class: a paper citing a nonexistent test
# FUNCTION, invisible to the file-level REF pattern). Hardened 2026-07-03.
BARE_FN = re.compile(r"\\texttt\{(test(?:\\_[A-Za-z0-9]+)+)\}")


def stem(ref: str) -> str:
    """'tests/test\\_fock\\_projection.py' -> 'test_fock_projection'
       'tests/test\\_qed\\_*.py'          -> 'test_qed_*' (glob kept)."""
    s = ref.replace("\\_", "_").replace("\\", "")
    s = s.split("/", 1)[1] if "/" in s else s
    return s[:-3] if s.endswith(".py") else s


def test_status(name: str) -> str:
    if "*" in name:                                   # glob ref: >=1 match suffices
        if any(TESTS.glob(f"{name}.py")):
            return "LIVE"
        if ARCHIVE.exists() and any(ARCHIVE.rglob(f"{name}.py")):
            return "ARCHIVED"
        return "MISSING"
    if (TESTS / f"{name}.py").exists():
        return "LIVE"
    if ARCHIVE.exists() and any(ARCHIVE.rglob(f"{name}.py")):
        return "ARCHIVED"
    return "MISSING"


def _gate_substr(argv: "list[str]") -> "str | None":
    for i, a in enumerate(argv):
        if a.startswith("--gate="):
            return a.split("=", 1)[1]
        if a == "--gate" and i + 1 < len(argv):
            return argv[i + 1]
    return None


def main() -> int:
    gate_corpus = "--gate-corpus" in sys.argv
    gate_substr = _gate_substr(sys.argv)
    files = sorted(p for p in PAPERS.rglob("*.tex") if "archive" not in p.parts)

    def is_gated(p: pathlib.Path) -> bool:
        if gate_corpus:
            return True
        if gate_substr:
            return gate_substr in str(p).replace("\\", "/")
        return p.name in TRUNK

    scope = ("the whole corpus" if gate_corpus else
             f"papers matching '{gate_substr}'" if gate_substr else "the trunk target")

    cited_stems: set[str] = set()
    gated_missing, audit_missing, archived = [], [], []
    for p in files:
        raw = p.read_text(encoding="utf-8", errors="replace")
        rel = p.relative_to(ROOT)
        for m in REF.finditer(raw):
            name = stem(m.group(0))
            cited_stems.add(name)
            st = test_status(name)
            ln = raw.count("\n", 0, m.start()) + 1
            if st == "MISSING":
                (gated_missing if is_gated(p) else audit_missing).append((rel, ln, name))
            elif st == "ARCHIVED":
                archived.append((rel, ln, name))
        # bare function-name citations: the named function must exist in tests/
        for m in BARE_FN.finditer(raw):
            fn = m.group(1).replace("\\_", "_")
            if fn.endswith(".py") or "*" in fn:
                continue
            span_l, span_r = max(0, m.start() - 40), m.end() + 8
            ctx = raw[span_l:span_r]
            if "tests/" in ctx or ".py" in raw[m.end():m.end() + 4]:
                continue  # the file-level REF pattern owns these
            hit = (ROOT / "tests" / (fn + ".py")).exists() or any(
                ("def " + fn + "(") in q.read_text(encoding="utf-8", errors="replace")
                for q in (ROOT / "tests").glob("*.py"))
            if not hit:
                ln = raw.count("\n", 0, m.start()) + 1
                (gated_missing if is_gated(p) else audit_missing).append(
                    (rel, ln, fn + " (function)"))

    # matrix coverage (advisory): matrix tests never cited inline anywhere
    matrix_tests: set[str] = set()
    if MATRIX.exists():
        for m in MATRIX_TEST.finditer(MATRIX.read_text(encoding="utf-8", errors="replace")):
            matrix_tests.add(m.group(0)[:-3] if m.group(0).endswith(".py") else m.group(0))
    uncovered = sorted(t for t in matrix_tests
                       if test_status(t) == "LIVE" and t not in cited_stems)

    print(f"\nInline test references: {len(cited_stems)} distinct   "
          f"[gated scope: {scope}]")
    print(f"  GATED-MISSING {len(gated_missing)}, AUDIT-MISSING {len(audit_missing)}, "
          f"ARCHIVED {len(archived)}, matrix-uncovered {len(uncovered)}")

    if gated_missing:
        print(f"\n*** GATED MISSING ({len(gated_missing)}) -- inline ref to a "
              f"non-existent test in {scope}: ***")
        for rel, ln, name in gated_missing:
            print(f"  {rel}:{ln}  ->  tests/{name}.py  NOT FOUND")
    if archived:
        print(f"\n--- ARCHIVED (advisory) ({len(archived)}) -- live claim cites a "
              f"test that lives only in tests/_archive/: ---")
        for rel, ln, name in archived:
            print(f"  {rel}:{ln}  ->  {name} (archived)")
    if audit_missing:
        print(f"\n--- AUDIT MISSING (advisory, out-of-scope) ({len(audit_missing)}): ---")
        for rel, ln, name in audit_missing:
            print(f"  {rel}:{ln}  ->  tests/{name}.py  NOT FOUND")
    if uncovered:
        print(f"\n--- MATRIX COVERAGE GAPS (advisory) ({len(uncovered)}) -- "
              f"matrix tests not cited inline in any paper: ---")
        for t in uncovered:
            print(f"  {t}")

    if gated_missing:
        print("\nRESULT: FAIL (paper cites a non-existent test in the gated scope)")
        return 1
    print(f"\nRESULT: PASS (every inline test reference in {scope} resolves to a live test)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
