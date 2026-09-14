#!/usr/bin/env python
"""C24 (PROPOSED) -- paragraph-splice detection in the .tex sources.

WHY THIS EXISTS
---------------
On 2026-09-12 an applier inserted a paragraph at an anchor that ended
MID-SENTENCE.  The result stranded a dangling "For" at the end of one
paragraph and left a later paragraph starting "water's $A_1$ block ...",
i.e. mid-clause and lowercase.

**Every gate passed.**  C10 compiled the paper with zero undefined references,
C19 found no eaten escapes, C16/C21/C22/C14 were silent.  They were all
working; none of them reads prose.  The defect was found by re-reading the
rendered seam, which is not a repeatable process.

The class is not new, and the corpus has already recorded an instance of it in
the criteria document itself: `docs/qa/criteria.md` records that C20's OWN
registration block was "spliced into the middle of the C19 sentence -- the same
mid-sentence-splice class the delta run caught", displacing a C18 sentence that
then needed re-homing.  Script-driven paragraph insertion is now the corpus's
standard editing mechanism, which makes the class systematic.

TWO CORRECTIONS THE FIRST DRAFT OF THIS GATE NEEDED
---------------------------------------------------
Both were found by `_prose_gate_probe.py`, which re-plants the REAL defect and
requires the gate to fire.  The first draft reported PASS on all three scopes
with ZERO advisory findings, and would have shipped as a gate that examines
nothing:

  1. **Environment tracking counted `\\begin{document}`.**  That never closes
     until the last line, so `inside_env` was true for the entire body and
     every paragraph was skipped.  Only prose-suppressing environments are
     tracked now, by name.
  2. **The conjunction was the wrong shape.**  The draft required one seam to
     BOTH end unterminated and be followed by a lowercase opener, reasoning
     that this would keep false positives near zero.  But when a whole
     paragraph is spliced into the middle of a sentence, the two halves land at
     OPPOSITE ENDS of the insertion -- the sentence breaks before it and
     resumes after it.  The real defect therefore triggered neither half of the
     conjunction.  The signals are now reported independently.

WHAT IT LOOKS FOR
-----------------
  (a) a prose paragraph whose final non-space character is not terminal
      punctuation (. ? ! : ;) -- a sentence that stops without ending;
  (b) a prose paragraph that opens with a lowercase ordinary word -- a
      sentence that starts already in progress.

Structural openers (\\section, \\item, \\begin, \\label, comments, math) and
anything inside tabular / verbatim / display-math / bibliography environments
are exempt.

Usage:
    python debug/qa/check_prose_continuity.py --gate paper_60
    python debug/qa/check_prose_continuity.py --gate group6 --list
    python debug/qa/check_prose_continuity.py --selftest
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qa_scopes import describe, resolve  # noqa: E402

TERMINAL = set(".?!:;")

# Environments inside which prose judgement is suppressed.  NOTE: `document`
# is deliberately absent -- counting it was the bug that made the first draft
# examine nothing.
SUPPRESS_ENV = {
    "tabular", "tabularx", "tabular*", "array", "verbatim", "lstlisting",
    "equation", "equation*", "align", "align*", "gather", "gather*",
    "multline", "multline*", "eqnarray", "eqnarray*", "displaymath",
    "thebibliography", "figure", "table", "tikzpicture", "matrix",
    "smallmatrix", "pmatrix", "bmatrix", "cases", "split",
}
ENV_RE = re.compile(r"\\(begin|end)\{([^}*]*\*?)\}")

STRUCTURAL_OPEN = re.compile(
    r"^\s*(?:%|\\(?:section|subsection|subsubsection|paragraph|begin|end|item|"
    r"label|caption|bibitem|newcommand|renewcommand|def|input|include|"
    r"bibliography|appendix|maketitle|tableofcontents|noindent|centering|"
    r"hline|toprule|midrule|bottomrule|addcontentsline|setlength|vspace|"
    r"hspace|clearpage|newpage|footnotetext|author|title|date|abstract)\b"
    r"|\s*[\[\$&]|\s*\\\[|\s*\\\\)"
)

LEADING_WRAP = re.compile(
    r"^\s*(?:\\(?:textbf|textit|emph|texttt|textsc|underline|textrm|mbox|"
    r"textnormal)\s*\{|\{)"
)

LOWER_WORD = re.compile(r"^[a-z][a-z'\-]{1,}(?:\s|$|[.,;:)])")

NONTERMINAL_OK = re.compile(
    r"(?:\\\\|\\(?:end|begin)\{[^}]*\}|\$\$|\\\]|\\\[|"
    r"\\(?:label|ref|cite|eqref|citep|citet)\{[^}]*\})\s*$"
)

# A line consisting solely of macro calls is structure, not prose:
# \maketitle, \appendix, \clearpage, \bibliography{refs}, ...
MACRO_ONLY = re.compile(
    r"^\s*(?:\\[A-Za-z@]+\*?(?:\[[^\]]*\])?(?:\{[^{}]*\})*\s*)+$"
)


def _strip_comment(line: str) -> str:
    out, esc = [], False
    for ch in line:
        if esc:
            out.append(ch)
            esc = False
            continue
        if ch == "\\":
            out.append(ch)
            esc = True
            continue
        if ch == "%":
            break
        out.append(ch)
    return "".join(out)


def paragraphs(text: str):
    lines = text.split("\n")
    buf: list[str] = []
    start = 1
    for i, raw in enumerate(lines, start=1):
        if raw.strip() == "":
            if buf:
                yield start, buf
                buf = []
            start = i + 1
        else:
            if not buf:
                start = i
            buf.append(raw)
    if buf:
        yield start, buf


#: closing delimiters that may follow the real terminal punctuation
CLOSERS = ")]}\"'`"

# "We summarize the result as" -> \textbf{Theorem (...).} is a real LaTeX
# idiom: prose introducing a DISPLAYED STATEMENT legitimately ends without
# terminal punctuation.  Exempted by statement word rather than by baseline,
# because the idiom generalizes.  Note this deliberately does NOT exempt a
# run-in tier tag such as \textbf{[MEASURED]}, which is exactly the shape a
# spliced paragraph has -- exempting all bold run-ins would have re-hidden
# the 2026-09-12 defect.
_OPENS_STATEMENT = re.compile(
    r"^\s*\\(?:textbf|emph|textit|textsc)\{\s*\(?(?:Theorem|Lemma|Proposition|"
    r"Corollary|Definition|Remark|Example|Claim|Proof|Conjecture|Axiom|Fact|"
    r"Observation|Assumption|Hypothesis|Notation)\b"
)

_OPENS_SUPPRESSED = re.compile(
    r"^\s*(?:\\\[|\$\$|\\begin\{(?:" + "|".join(
        re.escape(e) for e in sorted(SUPPRESS_ENV)) + r")\})"
)


def _classify(text: str):
    """Split into blocks tagged prose / display / structural.

    A prose paragraph that runs INTO a display equation legitimately ends
    without punctuation ("... we have" then \\begin{equation}), and one that
    RESUMES after a display legitimately starts lowercase ("vanishes (or ...").
    Both are normal LaTeX and neither is a splice, so each judgement is made
    only when the relevant NEIGHBOUR is also prose.  Missing this was the
    difference between 107 corpus-wide findings and 0.
    """
    blocks = []
    depth = 0
    for ln, raw in paragraphs(text):
        joined = "\n".join(_strip_comment(l) for l in raw)
        stripped = joined.strip()
        entered_at = depth
        for kind, name in ENV_RE.findall(joined):
            if name in SUPPRESS_ENV:
                depth += 1 if kind == "begin" else -1
        depth = max(0, depth)
        # A block is display only if it OPENS with a suppressed environment or
        # sits inside one.  Testing "contains a suppressed env" instead was too
        # coarse: an inline \begin{smallmatrix} mid-sentence disqualified the
        # whole prose paragraph, which silently dropped a real detection.
        opens_display = bool(_OPENS_SUPPRESSED.match(stripped))
        if not stripped:
            kind_ = "blank"
        elif entered_at > 0 or depth > 0 or opens_display:
            kind_ = "display"
        elif STRUCTURAL_OPEN.match(_strip_comment(raw[0])) or MACRO_ONLY.match(
                stripped.split("\n")[-1]):
            kind_ = "structural"
        else:
            kind_ = "prose"
        blocks.append({"line": ln, "text": stripped, "kind": kind_})
    return blocks


#: trailing decorations that may follow the real terminal punctuation
TRAILING_MACRO = re.compile(
    r"(?:\\(?:checkmark|qed|square|blacksquare|hfill|medskip|smallskip|"
    r"bigskip|par|noindent)\b\s*)+$"
)


def analyse(path: str):
    """Return (unterminated, lowercase_open) findings for one .tex file."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()

    # The preamble is never prose.  Judging it produced false positives on
    # \newtheorem / \newcommand / \newcolumntype blocks in three papers.
    cut = text.find("\\begin{document}")
    if cut != -1:
        text = "\n" * text[:cut].count("\n") + text[cut:]

    blocks = _classify(text)
    unterminated, lower_open = [], []

    for i, b in enumerate(blocks):
        if b["kind"] != "prose":
            continue
        nxt = blocks[i + 1]["kind"] if i + 1 < len(blocks) else None
        prv = blocks[i - 1]["kind"] if i > 0 else None

        # (b) opens mid-clause -- only meaningful if the PREVIOUS block is prose
        if prv == "prose":
            head = _strip_comment(b["text"].split("\n")[0])
            prev = None
            while prev != head:
                prev = head
                head = LEADING_WRAP.sub("", head, count=1)
            head = head.strip()
            if LOWER_WORD.match(head):
                lower_open.append((b["line"], head[:70]))

        # (a) stops without terminating -- only if the NEXT block is prose and
        # is not a displayed statement the sentence is introducing
        nxt_statement = (i + 1 < len(blocks)
                         and bool(_OPENS_STATEMENT.match(blocks[i + 1]["text"])))
        if nxt == "prose" and not nxt_statement:
            last = b["text"].rstrip()
            last_line = last.split("\n")[-1]
            if NONTERMINAL_OK.search(last) or MACRO_ONLY.match(last_line):
                continue
            trimmed = TRAILING_MACRO.sub("", last).rstrip()
            trimmed = trimmed.rstrip(CLOSERS).rstrip()
            tail_char = trimmed[-1] if trimmed else ""
            if tail_char and tail_char not in TERMINAL and tail_char not in "}]$,":
                unterminated.append((b["line"], last_line[-70:]))
    return unterminated, lower_open


def selftest() -> int:
    import tempfile
    cases = [
        ("clean.tex", "Alpha beta gamma.\n\nDelta epsilon zeta.\n", 0, 0),
        ("dangling.tex", "Alpha beta gamma.  For\n\nDelta epsilon.\n", 1, 0),
        ("resumes.tex", "Alpha beta gamma.\n\nwater's block is here.\n", 0, 1),
        ("real-shape.tex",
         "Alpha beta.  For\n\n\\textbf{[X]} Inserted paragraph here.\n\n"
         "water's block resumes.\n", 1, 1),
        ("math-ok.tex",
         "See the display\n\\begin{equation}\nx = 1\n\\end{equation}\n\n"
         "Next paragraph.\n", 0, 0),
    ]
    ok = True
    with tempfile.TemporaryDirectory() as d:
        for name, body, want_u, want_l in cases:
            p = os.path.join(d, name)
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(body)
            u, l = analyse(p)
            good = (len(u) == want_u and len(l) == want_l)
            ok = ok and good
            print(f"  [{'ok' if good else 'FAIL'}] {name}: "
                  f"unterminated={len(u)} (want {want_u}), "
                  f"lowercase-open={len(l)} (want {want_l})")
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", default="trunk")
    ap.add_argument("--list", action="store_true", help="list every finding")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()

    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
    all_tex = sorted(glob.glob(os.path.join(root, "papers", "**", "*.tex"),
                               recursive=True))
    files, warnings = resolve(args.gate, all_tex)
    for w in warnings:
        print(w)

    tot_u = tot_l = 0
    for f in files:
        u, l = analyse(f)
        tot_u += len(u)
        tot_l += len(l)
        if (u or l) and args.list:
            print(f"\n  {os.path.relpath(f, root)}")
            for ln, txt in u:
                print(f"    UNTERMINATED  L{ln}: ...{txt!r}")
            for ln, txt in l:
                print(f"    LOWERCASE-OPEN L{ln}: {txt!r}")

    n = tot_u + tot_l
    print(f"\n  unterminated paragraphs : {tot_u}")
    print(f"  lowercase-opening paras : {tot_l}")
    if n:
        print(f"\nRESULT: FAIL ({n} prose-continuity finding(s) in "
              f"{describe(args.gate, files)}; run --list to see them)")
        return 1
    print(f"\nRESULT: PASS (no prose-continuity findings in "
          f"{describe(args.gate, files)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
