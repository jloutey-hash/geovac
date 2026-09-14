"""DELTA remediation -- the ten paper-level findings (clean anchors).

Supersedes delta_fix_02_paper.py, whose three failing anchors were mis-wrapped;
the attempted repair went through a bash heredoc and was mangled by the
backslash-halving rule that `memory/feedback_no_heredoc_backslashes.md` exists
to prevent.  This file carries every anchor verbatim as read from the target.

See delta_fix_02_paper.py's docstring for the per-finding rationale.  The
headline one: M1 is an OVER-CORRECTION I introduced -- "values are truncated,
not rounded" is false, provably so from the exact 1s^2 value -(27/16)^2 =
-2.84765625, which truncates to -2.8476 while the paper prints -2.8477.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "rounded to four decimals"

EDITS = [
 ("M1-truncation",
  "ladder;\\ values are truncated, not rounded.",
  "ladder;\\ values are rounded to four decimals."),

 ("M2-sector-abstract",
  "while $2\\,^{1}S$ sits $\\gvq{p60_exc_ratio_k452}{1.08}\\times$ above it;\\ at\n$K=105$ the posing cost falls by",
  "while $2\\,^{1}S$ sits $\\gvq{p60_exc_ratio_k452}{1.08}\\times$ above it;\\ and on\nthe $s$-only ladder at $K=105$ the posing cost falls by"),

 ("M2-sector-body",
  "ladder at $K=105$ reads $\\gvq{p60_posing_cost_ground}{4.21}$,",
  "ladder --- $s$-only at this rung, where the $K=452$ pair above is full\n$s\\!+\\!p\\!+\\!d\\!+\\!f$ --- at $K=105$ reads $\\gvq{p60_posing_cost_ground}{4.21}$,"),

 ("M3-scope-abstract",
  "frontier on which the \\emph{conditioning} axis is removable for $M=2$ and\nnon-collinear $M=3$,",
  "frontier on which the \\emph{conditioning} axis is removable on $s$-sector\nshared-scale bases at $M=2$ and non-collinear $M=3$,"),

 ("M3-scope-conclusion",
  "Established for $M=2$ and non-collinear $M=3$;\\ the\ncollinear case is open and is not claimed.",
  "Established on $s$-sector shared-scale bases at $M=2$ and non-collinear\n$M=3$;\\ the collinear case is open and is not claimed."),

 ("M4-perfectly-conditioned",
  "a perfectly-conditioned gerade sector",
  "a gerade sector whose conditioning is flat"),

 ("M5-control-clause",
  "A control confirms the rotation is doing the work\nrather than the preconditioning as such:",
  "A control confirms that it is the rotation that \\emph{bounds} the\ngrowth, not the preconditioning alone:"),

 ("M6-completeness-qualified",
  "soften it rather than escaping it, because they are a\n\\emph{complete} set at one scale",
  "soften it rather than escaping it, because they are a\nset complete at one scale \\emph{in the atomic metric} --- a property the\nmolecular problem does not inherit, as measured below"),

 ("H1-bernstein",
  "exponent, and the Bernstein $\\Theta(\\kappa)$ floor for \\emph{both} problems",
  "exponent, and the $\\Theta(\\kappa)$ floor of Theorem~73 of~\\cite{gslw2019} for\n\\emph{both} problems"),

 ("H2-theorem73",
  "rather than contradicting it:\\ that floor\nconstrains polynomial approximation of $x^{-1/2}$ on $[\\kappa^{-1},1]$ under\nthe QSVT parity and boundedness constraints, and preconditioning changes\nthe operator rather than the polynomial, so $x^{-1/2}$ is never approximated on\nthe bad interval.",
  "rather than contradicting it:\\ that floor\nlower-bounds the number of block-encoding queries for an eigenvalue\ntransformation on the \\emph{promised spectral interval}, and preconditioning\nchanges that interval --- from $[\\kappa^{-1},1]$ to $[1/2.23,1]$ --- rather\nthan the polynomial approximated on it."),
]


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    bad = 0
    for name, old, _ in EDITS:
        if t.count(old) != 1:
            print(f"  {name}: anchor count={t.count(old)}")
            bad += 1
    if bad:
        print("ABORT -- no file written")
        return 2
    for name, old, new in EDITS:
        t = t.replace(old, new)
        print(f"  ok  {name}")
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print(f"applied: {len(EDITS)} paper fixes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
