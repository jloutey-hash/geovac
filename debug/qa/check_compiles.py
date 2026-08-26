#!/usr/bin/env python
"""C10 -- compile integrity, including UNDEFINED REFERENCES.

Why this script exists (2026-08-22 lesson).  C10 was being certified by running
pdflatex and reading its exit code.  That is not enough:

    pdflatex -halt-on-error EXITS 0 on undefined references and citations.

They are warnings, not errors.  So a paper can carry a dangling
``\\ref{sec:does_not_exist}`` -- which renders as a bare ``??`` in the PDF --
through every "C10 green" report ever emitted.  The 2026-08-22 full certifying
run found exactly that: THREE dangling cross-document references in Paper 19
(``sec:relationship``, ``eq:Vpk``, ``sec:conv_w1e_cross_domain_wall`` -- all
labels belonging to OTHER papers, which LaTeX can never resolve) plus an
undefined ``Note1`` endnote in the group3 synthesis, all pre-existing and all
invisible to the exit-code check.

This gate compiles each paper to a fixed point (up to three passes, since
cross-references need two and endnote machinery can need three) and FAILs on:
  - undefined references
  - undefined citations
  - genuine LaTeX errors (nonzero exit)
Font-substitution warnings are advisory: they are a local font-availability
matter, not a document defect.

External-BibTeX papers (2026-08-24 fix).  Papers that use ``\\bibliography{...}``
(an external ``.bib``) instead of an embedded ``thebibliography`` need a bibtex
run, which the pdflatex-only loop never did -- so EVERY citation in such a paper
came back undefined and the gate FAILed the whole paper spuriously (group4's
paper_20 was the only external-bib paper in scope, which is why it surfaced
there).  When the first pass writes a ``\\bibdata`` line to the ``.aux`` we run
the full revtex cycle instead: pdflatex, bibtex, pdflatex, bibtex, pdflatex,
pdflatex.  The SECOND bibtex cycle is required because revtex ``\\bibnote``
annotations can themselves contain ``\\cite`` (paper_20's ``Childs2021`` sits
inside the Trotter-feasibility bibnote), and that nested citation only becomes a
top-level ``\\citation`` after the first ``.bbl`` is processed.  bibtex's own
exit code is ignored on purpose: apsrev4-2 emits journal-macro warnings and
exits 1 on a perfectly good build, so only pdflatex's exit code and the final
log's undefined-ref/cite lines decide the verdict.

Usage:
    python debug/qa/check_compiles.py [--gate <substring>] [--passes N]
"""
from __future__ import annotations

import argparse
import glob
import io
import os
import re
import subprocess
import sys

# NB "Notes.bib" has NO leading dot: revtex4-2 writes <base>Notes.bib,
# so the dotted form never matched and left strays behind (found
# 2026-08-22 by the delta code review).
AUX_SUFFIXES = (".aux", ".log", ".out", "Notes.bib", ".toc", ".bbl",
                ".blg", ".spl", ".fdb_latexmk", ".fls")

UNDEF_REF = re.compile(r"Reference `([^']*)' on page [^ ]* undefined")
UNDEF_CITE = re.compile(r"Citation `([^']*)' on page [^ ]* undefined")


def _clean(d: str, base: str) -> None:
    for suf in AUX_SUFFIXES:
        f = os.path.join(d, base + suf)
        if os.path.exists(f):
            try:
                os.remove(f)
            except OSError:
                pass


def compile_one(path: str, passes: int):
    d, fname = os.path.split(path)
    base = os.path.splitext(fname)[0]
    # Clean BEFORE compiling, not only after.  A leftover .aux from an earlier
    # build can make an endnote/citation resolve that would NOT resolve from a
    # clean checkout, so without this the gate returns different verdicts for
    # the same source depending on what ran before it.  (Found 2026-08-22: the
    # target-scoped run reported P19/P24 ok while the corpus sweep, which
    # started clean, reported undefined `Note1` in both.)
    _clean(d, base)

    def _pdflatex() -> int:
        proc = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", fname],
            cwd=d, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return proc.returncode

    def _bibtex() -> None:
        # exit code intentionally ignored: apsrev4-2 exits 1 on a good build
        subprocess.run(["bibtex", base], cwd=d,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    rc = _pdflatex()  # pass 1: writes .aux (+ <base>Notes.bib for revtex)

    auxp = os.path.join(d, base + ".aux")
    uses_bibtex = False
    if os.path.exists(auxp):
        aux = io.open(auxp, encoding="utf-8", errors="ignore").read()
        uses_bibtex = "\\bibdata{" in aux

    if uses_bibtex:
        # full revtex cycle; the SECOND bibtex cycle resolves \cite nested
        # inside \bibnote annotations (see module docstring)
        _bibtex(); rc = _pdflatex()
        _bibtex(); rc = _pdflatex(); rc = _pdflatex()
    else:
        for _ in range(max(0, passes - 1)):
            rc = _pdflatex()

    logp = os.path.join(d, base + ".log")
    log = ""
    if os.path.exists(logp):
        log = io.open(logp, encoding="utf-8", errors="ignore").read()

    refs = sorted(set(UNDEF_REF.findall(log)))
    cites = sorted(set(UNDEF_CITE.findall(log)))
    fatal = "! LaTeX Error" in log or "! Emergency stop" in log
    fonts = log.count("Font Warning")

    _clean(d, base)

    return dict(base=base, rc=rc, refs=refs, cites=cites, fatal=fatal,
                fonts=fonts)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", default="", help="substring filter on the path")
    ap.add_argument("--passes", type=int, default=3)
    args = ap.parse_args()

    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    files = sorted(glob.glob(os.path.join(root, "papers", "**", "*.tex"),
                             recursive=True))
    files = [f for f in files if "/archive/" not in f.replace("\\", "/")]
    if args.gate:
        files = [f for f in files if args.gate in f.replace("\\", "/")]

    bad = 0
    for f in files:
        r = compile_one(f, args.passes)
        problems = []
        if r["fatal"] or r["rc"] != 0:
            problems.append(f"LaTeX error (exit {r['rc']})")
        if r["refs"]:
            problems.append("undefined ref(s): " + ", ".join(r["refs"]))
        if r["cites"]:
            problems.append("undefined citation(s): " + ", ".join(r["cites"]))
        note = f"  [{r['fonts']} font warn]" if r["fonts"] else ""
        if problems:
            bad += 1
            print(f"  FAIL {r['base']}{note}")
            for p in problems:
                print(f"       - {p}")
        else:
            print(f"  ok   {r['base']}{note}")

    scope = f" in scope '{args.gate}'" if args.gate else ""
    if bad:
        print(f"\nRESULT: FAIL ({bad} document(s) with compile/reference "
              f"defects{scope})")
        return 1
    print(f"\nRESULT: PASS ({len(files)} document(s) compile with zero "
          f"undefined references or citations{scope})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
