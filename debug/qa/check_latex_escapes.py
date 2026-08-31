#!/usr/bin/env python
"""C19 -- eaten-escape corruption gate.

The failure class: a paper edit applied through a Python replacement string in
which a LaTeX control sequence was written WITHOUT a raw-string prefix, so
Python consumed the backslash as an escape before the text ever reached the
file.  ``\\ref`` becomes CR + ``ef``; ``\\textbf`` becomes BS + ``extbf``;
``\\times`` becomes TAB + ``imes``.

NOT detectable here: a swallowed ``\\'`` (as in ``\\'e`` -> ``'e``)
leaves no control character behind, so nothing distinguishes it from
an ordinary apostrophe. That form must be caught by C10 rendering or
by review, not by this gate -- do not claim it.

Why this needs its OWN gate rather than riding on C10 (compiles):

    ``Sec.~\\ref{sec:obstruction}`` corrupted to ``Sec.~<CR>ef{sec:obstruction}``
    COMPILES CLEAN.  LaTeX renders the literal text "ef{sec:obstruction}".
    No error, no warning, no undefined reference -- the cross-reference simply
    silently disappears from the document.

Three instances have now reached the corpus this way (two TAB corruptions in
Paper 59 during the v4.107.0 QA remediation, one of which destroyed an
[OBSERVATION] tag; a backspace pair in the C17 registry regexes that silently
dead-lettered a guard; and the CR corruption in Paper 59's non-classicality
attribution pointer, found by the 2026-08-22 claims re-run).  Each was found by
a human or an LLM reviewer noticing odd rendering -- which is exactly the
low-salience, zero-variance class that a grep should own instead.

The standing cure is upstream (write the script to a file with raw strings,
never inline a LaTeX-bearing replacement in a bash heredoc); this gate is the
backstop for when that discipline slips.

Usage:
    python debug/qa/check_latex_escapes.py [--gate <substring>] [--selftest]
"""
from __future__ import annotations

import argparse
import glob
import io
import os
import re
import sys

# Control characters that a Python escape produces from a LaTeX macro, paired
# with the macro tail they leave behind.
ESCAPE_ARTIFACTS = [
    ("\r", "ef",     r"\ref"),
    ("\r", "aggedright", r"\raggedright"),
    ("\b", "extbf",  r"\textbf"),
    ("\b", "extit",  r"\textit"),
    ("\b", "exttt",  r"\texttt"),
    ("\b", "egin",   r"\begin"),
    ("\b", "ibitem", r"\bibitem"),
    ("\t", "imes",   r"\times"),
    ("\t", "extbf",  r"\textbf"),
    ("\t", "au",     r"\tau"),
    ("\t", "ag",     r"\tag"),
    ("\n", "ewcommand", r"\newcommand"),
    ("\n", "onumber",   r"\nonumber"),
    ("\f", "rac",    r"\frac"),
    ("\v", "space",  r"\vspace"),
    ("\a", "lpha",   r"\alpha"),
]

# Bare control characters that are NEVER legitimate in a .tex source.  A raw
# TAB is excluded here because it is legal (and common) as indentation inside
# tikzpicture/tabular bodies -- a TAB that is actually corruption shows up in
# ESCAPE_ARTIFACTS instead, glued to a macro tail.  A TAB is flagged only when
# it sits mid-line against non-whitespace, which indentation never does.
# A lone CR cannot be legitimate: CRLF is normalised away before
# scanning (see main()), so any surviving \r is a swallowed escape --
# and catching it generically covers \right, \rho, \rangle, \rule,
# \raggedright, ... rather than only the two tails once listed.
# Same reasoning for \f (\frac), \v (\vspace), \a (\alpha), \b.
BARE_CONTROL = re.compile(r"[\r\b\x0b\x0c\x07]")


def scan_text(text: str):
    """Return a list of (line_no, kind, detail, excerpt)."""
    findings = []
    lines = text.split("\n")

    for ctrl, tail, macro in ESCAPE_ARTIFACTS:
        if ctrl == "\n":
            # A swallowed newline leaves the tail at the START of a line.
            for i, line in enumerate(lines, 1):
                if line.startswith(tail + "{") or line.startswith(tail + "["):
                    findings.append((i, "eaten-escape",
                                     f"line begins with '{tail}' -- looks like "
                                     f"a swallowed '{macro}'", line[:90]))
            continue
        needle = ctrl + tail
        start = 0
        while True:
            j = text.find(needle, start)
            if j < 0:
                break
            ln = text[:j].count("\n") + 1
            findings.append((ln, "eaten-escape",
                             f"{ctrl!r} + '{tail}' -- looks like a swallowed "
                             f"'{macro}'", text[max(0, j - 30):j + 50]
                             .replace("\n", " ")))
            start = j + 1

    for m in BARE_CONTROL.finditer(text):
        ln = text[:m.start()].count("\n") + 1
        findings.append((ln, "bare-control",
                         f"raw {m.group()!r} in .tex source",
                         text[max(0, m.start() - 30):m.start() + 50]
                         .replace("\n", " ")))

    # TAB, line-aware: benign as indentation, suspicious mid-line.
    for i, line in enumerate(lines, 1):
        for j, ch in enumerate(line):
            if ch != "\t":
                continue
            if line[:j].strip() == "":
                continue          # pure indentation -- legitimate
            findings.append((i, "bare-control",
                             "raw '\\t' mid-line in .tex source "
                             "(swallowed '\\t...' escape?)",
                             line[max(0, j - 30):j + 50]))

    return findings


def selftest() -> int:
    positives = [
        "See Sec.~\refoo",                      # noqa -- literal CR+ef
        "value \times 2",                        # literal TAB+imes
        "\bextbf{bold}",                         # literal BS+extbf
        "onumber{}\nfoo",                        # swallowed newline form
    ]
    positives[0] = "See Sec.~" + "\r" + "ef{sec:x}"
    positives[1] = "value " + "\t" + "imes 2"
    positives[2] = "\b" + "extbf{bold}"
    positives[3] = "onumber{x}"
    # Forms the 2026-08-22 delta review found MISSED by the tail-list-only
    # version. These are far more common in this corpus's math than \ref,
    # and are now covered generically by the CR entry in BARE_CONTROL.
    positives += [
        "\\left( x " + "\r" + "ight)",        # \right
        "the density " + "\r" + "ho_{ab}",     # \rho
        "inner " + "\r" + "angle x",           # \rangle
        "a " + "\r" + "ule{1pt}{1pt}",         # \rule
        "angle " + "\t" + "heta_{12}",         # \theta, mid-line
        "\\begin{center}" + "\r" + "aggedright",  # \raggedright
    ]
    negatives = [
        r"See Sec.~\ref{sec:x}",
        r"value \times 2",
        r"\textbf{bold}",
        r"Th\'eor\`eme relatif au mouvement d'un point attir\'e",
        r"\begin{equation}\nonumber x \end{equation}",
    ]
    bad = 0
    for t in positives:
        if not scan_text(t):
            print(f"  SELFTEST MISS (should flag): {t!r}")
            bad += 1
    for t in negatives:
        f = scan_text(t)
        if f:
            print(f"  SELFTEST FALSE POSITIVE: {t!r} -> {f}")
            bad += 1
    print(f"selftest: {len(positives)} positive / {len(negatives)} negative "
          f"cases, {bad} failure(s)")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", default="", help="substring filter on the path")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        return selftest()

    root = os.path.join(os.path.dirname(__file__), "..", "..")
    files = sorted(glob.glob(os.path.join(root, "papers", "**", "*.tex"),
                             recursive=True))
    files = [f for f in files
             if "/archive/" not in f.replace("\\", "/")]
    if args.gate:
        # Named scopes resolve to an explicit file list; anything else keeps
        # the historical substring behaviour.  `trunk` NEEDS a named scope:
        # no trunk paper's path contains "trunk", so the substring filter
        # matched nothing and this gate reported PASS on 0 papers -- it had
        # never examined a trunk document.  Scope per docs/qa/trunk.done.md.
        NAMED = {
            "trunk": [
                "papers/group3_foundations/Paper_0_Geometric_Packing.tex",
                "papers/group3_foundations/paper_1_spectrum.tex",
                "papers/group3_foundations/Paper_7_Dimensionless_Vacuum.tex",
                "papers/group1_operator_algebras/paper_32_spectral_triple.tex",
                "papers/group1_operator_algebras/"
                "paper_38_su2_propinquity_convergence.tex",
                "papers/synthesis/group3_foundations_synthesis.tex",
            ],
        }
        if args.gate in NAMED:
            want = set(NAMED[args.gate])
            # glob yields ABSOLUTE paths; the scope list is repo-relative, so
            # match by suffix. Comparing them directly matched nothing.
            files = [f for f in files
                     if any(f.replace("\\", "/").endswith(w) for w in want)]
            hit = {w for w in want
                   if any(f.replace("\\", "/").endswith(w) for f in files)}
            if hit != want:
                print(f"   [scope] WARNING: {len(want - hit)} file(s) in the "
                      f"'{args.gate}' scope were not found: {sorted(want - hit)}")
        else:
            files = [f for f in files if args.gate in f.replace("\\", "/")]

    total = 0
    for f in files:
        text = io.open(f, encoding="utf-8", errors="ignore",
                       newline="").read()
        # CRLF line endings are normal on this platform: neutralise them
        # before scanning so only INTERIOR control characters survive.
        text = text.replace("\r\n", "\n")
        for ln, kind, detail, excerpt in scan_text(text):
            rel = os.path.relpath(f, root).replace("\\", "/")
            print(f"  [{kind}] {rel}:{ln}  {detail}")
            print(f"      ...{excerpt.strip()}...")
            total += 1

    scope = f" in scope '{args.gate}'" if args.gate else ""
    if total:
        print(f"\nRESULT: FAIL ({total} eaten-escape / bare-control "
              f"corruption(s){scope})")
        return 1
    print(f"\nRESULT: PASS (no eaten-escape corruption in {len(files)} "
          f"paper(s){scope})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
