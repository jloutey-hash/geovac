r"""C15 — inline arXiv-ID consistency check (deterministic).

Catches the failure class the LLM citation reviewers systematically miss:\ an
arXiv ID written in the *body prose* that does not match any arXiv ID in the
paper's own bibliography (a transposed/typo'd inline ID while the bibitem is
correct). This was the /qa group1 re-cert-1 calibration gap (an inline
`arXiv:2504.10830` seed slipped past the bibitem-focused citation reviewers).

Rule: an arXiv ID mentioned inline (in the body, before
`\begin{thebibliography}`, outside LaTeX comments) that is NOT in the paper's
bibliography but is a *near-match* (same length, Hamming distance <= 2) to a
bibitem ID is almost certainly a transposed/typo'd version of that bibitem ID
and is FAILed. (Inline IDs far from every bibitem ID are NOT flagged --- a
paper may legitimately mention an arXiv ID inline whose bibitem lists only the
journal reference; that is not this check's target.)

This is a *complement* to the LLM citation-reviewer (C4), not a replacement:\
it catches transposed/typo'd IDs deterministically and cheaply; the reviewer
still owns wrong-title / wrong-venue / wrong-attribution (right ID, wrong
metadata — e.g. the bizi "Spectral action in Lorentzian signature" title slip,
which this check cannot see because the ID itself is correct).

Version suffixes (`v4`) and old-style IDs (`math-ph/0110001`) are normalized.
LaTeX-comment lines (content after an unescaped `%`) are ignored.

Scope: `--gate <substr>` gates papers whose path contains <substr> (default =
the trunk set); `--gate-corpus` gates every paper.

  python debug/qa/check_inline_arxiv.py --gate group1
  python debug/qa/check_inline_arxiv.py --gate-corpus
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

# arXiv IDs: new-style NNNN.NNNNN(vN) and old-style archive/NNNNNNN(vN).
ARXIV = re.compile(
    r"arXiv:\s*((?:[a-z][\w.-]*/)?\d{4}\.\d{4,5}|[a-z-]+(?:\.[A-Z]{2})?/\d{7})"
    r"(?:v\d+)?",
    re.IGNORECASE,
)


def _strip_comment(line: str) -> str:
    """Drop content after an unescaped % (LaTeX comment)."""
    out = []
    i = 0
    while i < len(line):
        c = line[i]
        if c == "\\" and i + 1 < len(line):
            out.append(line[i : i + 2])
            i += 2
            continue
        if c == "%":
            break
        out.append(c)
        i += 1
    return "".join(out)


def _norm(arxiv_id: str) -> str:
    return arxiv_id.strip().lower()


def _near(a: str, b: str) -> bool:
    """True iff a,b are same-length and differ in <= 2 positions (a transposed
    or single/double-digit typo of the same ID), but not identical."""
    if a == b or len(a) != len(b):
        return False
    return sum(1 for x, y in zip(a, b) if x != y) <= 2


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

    scope = (
        "the whole corpus" if gate_corpus
        else f"scope '{gate_substr}'" if gate_substr
        else "the trunk"
    )

    failures: list[str] = []
    n_checked = 0
    for p in files:
        if not is_gated(p):
            continue
        n_checked += 1
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
        # split body vs bibliography
        bib_start = next(
            (i for i, ln in enumerate(lines)
             if "\\begin{thebibliography}" in ln or "\\bibliography{" in ln),
            len(lines),
        )
        declared: set[str] = set()
        for ln in lines[bib_start:]:
            for m in ARXIV.finditer(_strip_comment(ln)):
                declared.add(_norm(m.group(1)))
        for lineno, ln in enumerate(lines[:bib_start], start=1):
            for m in ARXIV.finditer(_strip_comment(ln)):
                aid = _norm(m.group(1))
                if aid in declared:
                    continue
                near = next((d for d in declared if _near(aid, d)), None)
                if near is not None:
                    failures.append(
                        f"  {p.relative_to(ROOT)}:{lineno}  ->  inline "
                        f"arXiv:{m.group(1)}  is a near-match (transposed/typo'd?) "
                        f"of bibitem arXiv:{near} but does not equal it"
                    )

    print(f"inline arXiv-ID check over {n_checked} papers  [{scope}]\n")
    if failures:
        print("\n".join(failures))
        print(
            f"\nRESULT: FAIL -- {len(failures)} inline arXiv ID(s) in {scope} "
            f"are near-match transpositions of a bibitem ID (likely typos)."
        )
        return 1
    print(
        f"RESULT: PASS -- no inline arXiv ID in {scope} is a near-match "
        f"transposition of a bibitem ID."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
