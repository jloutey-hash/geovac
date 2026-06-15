#!/usr/bin/env python3
r"""
Deterministic K-label hard-prohibition screen  --  the /qa C5 tripwire backstop.

CLAUDE.md §13.5 forbids presenting the Paper 2 combination rule
K = pi(B + F - Delta) as anything stronger than an Observation -- never
"derived", "a theorem", "proven", "conjecture", or "conjectural". An LLM
reviewer can MISS this under a self-contradicting hedge (the /qa run-#4 blind
spot, 2026-06-15: a reviewer SAMPLED the K-appearances and missed a planted
"K is now derived as a theorem" because a later clause said "does not promote
K to a theorem"). This script is the deterministic backstop -- the C5 analog
of the C11 internal-title check.

Unlike the title check (an EXACT string match -> clean PASS/FAIL), the
K-prohibition is SEMANTIC: "K is NOT derived", a "non-selection theorem",
"K WOULD become derived", a future-work "Derive K" goal are all COMPLIANT,
while "K IS derived / a theorem / conjectural" is the violation. So the screen
finds every tier word that sits in the SAME CLAUSE as a K-rule reference and
clears the compliant patterns (negation / hypothetical / future-goal /
non-selection); anything left is SUSPECT.

Scope partition (matches the /qa trunk target):
  * TRUNK papers (0,1,7,32,38 + group3 synthesis) -> a SUSPECT here FAILS the
    gate (exit 1). This set has had the conjecture->Observation sweep, so it
    should be clean.
  * the rest of the corpus -> SUSPECT listed as an ADVISORY AUDIT (does NOT
    fail the gate); this is where the un-swept "conjectural" drift lives, for
    PI decision on a corpus-wide sweep.

It BACKS the claims-reviewer (guarantees exhaustive enumeration, the thing the
LLM failed at); it does not replace adjudication.

Exit 0 = no SUSPECT in the trunk set. Exit 1 = >=1 trunk SUSPECT.

Usage:
  python debug/qa/check_k_label.py            # trunk gate + corpus audit
  python debug/qa/check_k_label.py --all      # also print COMPLIANT enumeration
  python debug/qa/check_k_label.py --gate-corpus   # FAIL on ANY corpus SUSPECT
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
PAPERS = ROOT / "papers"

# the /qa trunk target (where C5 is certified clean)
TRUNK = {
    "Paper_0_Geometric_Packing.tex",
    "paper_1_spectrum.tex",
    "Paper_7_Dimensionless_Vacuum.tex",
    "paper_32_spectral_triple.tex",
    "paper_38_su2_propinquity_convergence.tex",
    "group3_foundations_synthesis.tex",
}

# a reference to the K = pi(B+F-Delta) COMBINATION rule (NOT the modular K=J)
K_ANCHOR = re.compile(
    r"combination\s+(?:rule|formula|law)"
    r"|B\s*\+\s*F\s*-\s*\\?Delta"
    r"|\\?alpha[-\s]?conjecture"
    r"|K\s*(?:\\stackrel\s*\{[^}]*\}\s*\{=\}|=)\s*\\?pi\s*\\?\(?\s*B",
    re.IGNORECASE,
)

# prohibited tier words (positive assertion that K is more than an Observation)
TIER = re.compile(
    r"\b(deriv(?:e|ed|es|ation|able)|conjectur(?:e|ed|es|al)|theorem|"
    r"proven|proof|proves?|predict(?:s|ed|ion)?)\b",
    re.IGNORECASE,
)

NEG = re.compile(r"\b(not|no|never|cannot|can't|without|non|nor|n't)\b", re.IGNORECASE)
HYP = re.compile(
    r"\b(would|could|if|were|may|might|should|hypothetic\w*|future|any|remains?\s+to\s+be|"
    r"open|toward|prospect\w*)\b",
    re.IGNORECASE,
)
NONSEL = re.compile(
    r"non-?selection|single[-\s]mechanism|no\s+morphism|no\s+mechanism|"
    r"impossibilit\w*|cannot\s+be\s+generated|not\s+\w+\s+generate|no-?go|"
    r"eliminat\w*",
    re.IGNORECASE,
)
# imperative future-work goal at clause start: "Derive the combination rule",
# "Prove the combination rule"
IMPERATIVE = re.compile(r"^\s*(?:\\item\s*)?(?:\\\w+\{[^}]*\}\s*)?(derive|prove)\b", re.IGNORECASE)
CLAUSE_BREAK = re.compile(r"[.;:]\s|\n\s*\n|\\item")


def _blank(m: "re.Match[str]") -> str:
    return " " * len(m.group(0))


_XREF = (r"theorem|conjecture|proposition|lemma|corollary|prediction|observation"
         r"|eq|equation|section|figure|table")


def strip_latex(text: str) -> str:
    text = re.sub(r"(?<!\\)%[^\n]*", _blank, text)                  # comments
    text = re.sub(r"\\(?:begin|end)\{[^}]*\}", _blank, text)        # env markers
    text = re.sub(
        r"\\(?:label|ref|eqref|autoref|pageref|cref|Cref|cite[a-z]*)\{[^}]*\}",
        _blank, text,
    )                                                               # refs/labels
    # font wrappers -> spaces (keep content): "\textit{Derive ...}" -> "  Derive ..."
    text = re.sub(r"\\(?:emph|textit|textbf|textsf|texttt|text|mathrm|mathit|mathbf)\{",
                  _blank, text)
    # cross-reference words: "Theorem~\ref" (now "Theorem~  ") and "Theorem 3"
    text = re.sub(rf"\b(?:{_XREF})s?\s*~", _blank, text, flags=re.IGNORECASE)
    text = re.sub(rf"\b(?:{_XREF})s?\s+\d", _blank, text, flags=re.IGNORECASE)
    return text


def line_of(text: str, pos: int) -> int:
    return text.count("\n", 0, pos) + 1


def clause_of(text: str, start: int, end: int) -> str:
    """The clause (between the nearest .;: / blank-line / \\item breaks) holding
    a tier word at [start, end)."""
    left = 0
    for m in CLAUSE_BREAK.finditer(text, 0, start):
        left = m.end()
    rb = CLAUSE_BREAK.search(text, end)
    right = rb.start() if rb else len(text)
    return text[left:right]


def scan_file(path: pathlib.Path) -> "list[tuple[int,str,str,str]]":
    raw = path.read_text(encoding="utf-8", errors="replace")
    text = strip_latex(raw)
    hits = []
    for m in TIER.finditer(text):
        clause = clause_of(text, m.start(), m.end())
        if not K_ANCHOR.search(clause):          # tier word must share a clause
            continue                              # with a K-rule reference
        low = clause.lower()
        # local "before the tier word" inside the clause, for direction
        rel = m.start() - text.rindex(clause) if clause in text else 0
        before = low[max(0, rel - 60):rel]
        neg_ok = NEG.search(before) and "no longer" not in before
        verdict = "SUSPECT"
        if neg_ok or HYP.search(low) or NONSEL.search(low) or IMPERATIVE.search(clause):
            verdict = "COMPLIANT"
        ln = line_of(text, m.start())
        snippet = re.sub(r"\s+", " ", clause.strip())[:150]
        hits.append((ln, m.group(0), verdict, snippet))
    return hits


def _gate_substr(argv: "list[str]") -> "str | None":
    for i, a in enumerate(argv):
        if a.startswith("--gate="):
            return a.split("=", 1)[1]
        if a == "--gate" and i + 1 < len(argv):
            return argv[i + 1]
    return None


def main() -> int:
    show_all = "--all" in sys.argv
    gate_corpus = "--gate-corpus" in sys.argv
    gate_substr = _gate_substr(sys.argv)        # branch sweep: --gate group5_qed_gauge
    files = sorted(p for p in PAPERS.rglob("*.tex") if "archive" not in p.parts)

    def is_gated(p: pathlib.Path) -> bool:
        if gate_corpus:
            return True
        if gate_substr:
            return gate_substr in str(p).replace("\\", "/")
        return p.name in TRUNK                   # default gated scope = the trunk target

    scope = ("the whole corpus" if gate_corpus else
             f"papers matching '{gate_substr}'" if gate_substr else "the trunk target")

    gated_suspect, audit_suspect, compliant = [], [], 0
    for p in files:
        gated = is_gated(p)
        for ln, word, verdict, snip in scan_file(p):
            rel = p.relative_to(ROOT)
            if verdict == "COMPLIANT":
                compliant += 1
                if show_all:
                    print(f"  [COMPLIANT] {rel}:{ln}  ({word})  {snip}")
            elif gated:
                gated_suspect.append((rel, ln, word, snip))
            else:
                audit_suspect.append((rel, ln, word, snip))

    total = compliant + len(gated_suspect) + len(audit_suspect)
    print(f"\nK-rule clauses carrying a tier word: {total}  "
          f"(COMPLIANT {compliant}, GATED-SUSPECT {len(gated_suspect)}, "
          f"AUDIT-SUSPECT {len(audit_suspect)})   [gated scope: {scope}]")

    if gated_suspect:
        print(f"\n*** GATED SUSPECT ({len(gated_suspect)}) -- candidate C5 "
              f"hard-prohibition touch(es) in {scope}: ***")
        for rel, ln, word, snip in gated_suspect:
            print(f"  [{word}] {rel}:{ln}\n      {snip}")

    if audit_suspect:
        print(f"\n--- AUDIT (advisory, out-of-scope; does NOT fail the gate) "
              f"({len(audit_suspect)}) -- K-label drift for the branch sweep: ---")
        for rel, ln, word, snip in audit_suspect:
            print(f"  [{word}] {rel}:{ln}\n      {snip}")

    if gated_suspect:
        print("\nRESULT: FAIL (K-label criterion -- SUSPECT assertion(s) in gated scope)")
        return 1
    print(f"\nRESULT: PASS (no SUSPECT K-tier assertion in {scope})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
