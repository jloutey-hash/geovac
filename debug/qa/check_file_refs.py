#!/usr/bin/env python3
r"""
Deterministic paper->file reference integrity check  --  the /qa C14 criterion.

The C13 companion (check_paper_test_refs.py) verifies inline `tests/test_*.py`
references resolve. It does NOT cover the OTHER file paths papers cite as
backing -- `geovac/` code modules, `debug/` scripts and memos, `benchmarks/`,
`demo/`. The group3 run-4 defect (Paper 55 `thm:jlo_depth2_reading_A` cited a
nonexistent `geovac/jlo_chi.py`) is exactly that gap: a load-bearing claim
pointing at a code module that does not exist. This check closes it.

GATE (deterministic, FAILs) -- the SERIOUS class:
    every inline file reference into a PERMANENT code/artifact dir
    (geovac/, benchmarks/, demo/) in a gated paper resolves to a REAL file
    under the repo (glob refs need >=1 match). A claim pointing at a code
    module that does not exist -- the group3 run-4 `geovac/jlo_chi.py`
    defect -- is broken BACKING and is the failure mode this gates.

ADVISORY (reported, does NOT fail) -- the hygiene class:
    references into debug/ (sprint memos and scripts). debug/ is the
    transient clean-room directory (CLAUDE.md sec.9) that is pruned over
    time, so a paper citing `debug/foo_memo.md` goes stale BY DESIGN once
    the memo is cleaned up. These are dangling documentation pointers
    (a smell -- papers should cite the permanent record, not transient
    debug/ files), reported for a hygiene sweep but NOT a cert blocker.

    tests/ refs are intentionally EXCLUDED here -- they are C13's domain
    (which also reasons about tests/_archive/ and the claim->test matrix).

Scope: --gate <substr> gates papers whose path contains <substr> (default =
the trunk target); --gate-corpus gates every paper.

Exit 0 = no MISSING ref in the gated scope. Exit 1 = >=1 MISSING.

Usage:
  python debug/qa/check_file_refs.py
  python debug/qa/check_file_refs.py --gate group3_foundations
  python debug/qa/check_file_refs.py --gate-corpus
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
PAPERS = ROOT / "papers"

TRUNK = {
    "Paper_0_Geometric_Packing.tex",
    "paper_1_spectrum.tex",
    "Paper_7_Dimensionless_Vacuum.tex",
    "paper_32_spectral_triple.tex",
    "paper_38_su2_propinquity_convergence.tex",
    "group3_foundations_synthesis.tex",
}

# A file reference: one of the tracked dirs, then a path (LaTeX-escaped
# underscores `\_` allowed, glob `*` allowed), then a tracked extension.
# tests/ is deliberately not a tracked dir here (C13 owns it).
REF = re.compile(
    r"(?:geovac|debug|benchmarks|demo)/[\w*/.\\-]+\.(?:py|md|json|txt|tex|csv)"
)


def norm(ref: str) -> str:
    r"""'geovac/jlo\_chi.py' -> 'geovac/jlo_chi.py' (strip LaTeX escapes)."""
    return ref.replace("\\_", "_").replace("\\%", "%").replace("\\", "")


def resolves(relpath: str) -> bool:
    """True iff relpath (repo-relative) exists; globs need >=1 match."""
    if "*" in relpath:
        return any(ROOT.glob(relpath))
    return (ROOT / relpath).exists()


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

    distinct: set[str] = set()
    gated_missing, audit_missing, debug_hygiene = [], [], []
    for p in files:
        raw = p.read_text(encoding="utf-8", errors="replace")
        rel = p.relative_to(ROOT)
        for m in REF.finditer(raw):
            ref = norm(m.group(0))
            distinct.add(ref)
            if not resolves(ref):
                ln = raw.count("\n", 0, m.start()) + 1
                top = ref.split("/", 1)[0]
                if top == "debug":            # transient clean-room dir: advisory
                    debug_hygiene.append((rel, ln, ref))
                elif is_gated(p):             # serious class (code/artifact), in scope
                    gated_missing.append((rel, ln, ref))
                else:                         # serious class, out of scope
                    audit_missing.append((rel, ln, ref))

    print(f"\nInline geovac/benchmarks/demo (gated) + debug (advisory) file "
          f"references: {len(distinct)} distinct   [gated scope: {scope}]")
    print(f"  GATED-MISSING {len(gated_missing)} (serious: code/artifact dirs), "
          f"AUDIT-MISSING {len(audit_missing)}, debug-hygiene {len(debug_hygiene)}")

    if gated_missing:
        print(f"\n*** GATED MISSING ({len(gated_missing)}) -- paper cites a "
              f"non-existent code/artifact file in {scope}: ***")
        for rel, ln, ref in gated_missing:
            print(f"  {rel}:{ln}  ->  {ref}  NOT FOUND")
    if audit_missing:
        print(f"\n--- AUDIT MISSING (serious, advisory out-of-scope) "
              f"({len(audit_missing)}): ---")
        for rel, ln, ref in audit_missing:
            print(f"  {rel}:{ln}  ->  {ref}  NOT FOUND")
    if debug_hygiene:
        print(f"\n--- DEBUG-HYGIENE (advisory) ({len(debug_hygiene)}) -- dangling "
              f"pointers into the transient debug/ clean-room dir; not a cert "
              f"blocker, candidate for a hygiene sweep: ---")
        for rel, ln, ref in debug_hygiene:
            print(f"  {rel}:{ln}  ->  {ref}  NOT FOUND")

    if gated_missing:
        print("\nRESULT: FAIL (paper cites a non-existent code/artifact file "
              "in the gated scope)")
        return 1
    print(f"\nRESULT: PASS (every gated geovac/benchmarks/demo file reference "
          f"in {scope} resolves; debug/ refs advisory)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
