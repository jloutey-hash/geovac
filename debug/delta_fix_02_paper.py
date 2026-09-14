"""DELTA remediation -- the eight paper-level findings.

Claims delta (7 SMALL) + citation delta (2 SMALL).  Each verified by the PM.

M1 -- THE OVER-CORRECTION, and it is mine.  The fix for the mixed-basis chain
     ADDED "values are truncated, not rounded".  That is FALSE, and provably so
     without any measurement: the 1s^2 rung is exactly -(27/16)^2 =
     -2.84765625, which TRUNCATES to -2.8476 and ROUNDS to -2.8477 -- and the
     paper prints -2.8477.  All four rungs are rounded (s: -2.874468 ->
     -2.8745; +p: -2.894672 -> -2.8947).  I reasoned "truncation" from the old
     THREE-decimal display and carried it over when I moved to four.
M2 -- the s-only sector label reached the conclusion and not the abstract or
     body: the same silent-basis-mix defect the chain fix was written to close,
     surviving one locus away.
M3 -- the third-lever scope omits "s-sector" at the abstract and conclusion;
     the string occurs exactly ONCE in the paper.  C8#15 lists it first.
M4 -- "perfectly-conditioned gerade sector" survived the abstract/conclusion
     sweep at its own locus; tab:resource charges that route d_inv = 18.
M5 -- the control sentence contradicts itself: clause A says preconditioning is
     not doing the work, clause B measures it halving the exponent 1.96 -> 0.98.
     The rotation supplies BOUNDEDNESS, which is the true statement.
M6 -- "complete set at one scale" unqualified inside a MOLECULAR comparison,
     270 lines before the caveat that says completeness there is atomic-metric
     only.
H1 -- "Bernstein" survives uncited at one locus and names a mechanism the cited
     source does not use (Theorem 73 is proved via reflection operators, not a
     Bernstein/Markov polynomial inequality).
H2 -- Theorem 73 constrains the block-encoding QUERY COUNT for eigenvalue
     transformation, not polynomial approximation under parity/boundedness
     constraints.  PM NOTE: I could not fetch the theorem text directly (arXiv
     abstract page and ar5iv both failed to surface it); two independent
     reviewers quote it identically from the PDF.  The change is in the
     CONSERVATIVE direction -- it weakens our characterization of someone
     else's theorem -- so it is safe under either reading.  The escape argument
     survives: preconditioning changes the promised spectral interval.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "rounded to four decimals"

EDITS = [
 ("M1-truncation",
  "is a pure $\\ell_{\\max}$ ladder;\\ values are truncated, not\nrounded.",
  "is a pure $\\ell_{\\max}$ ladder;\\ values are rounded to four decimals."),

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
  "A control confirms the rotation is doing the work rather than the\npreconditioning as such:\\ the naive $\\mathrm{blockdiag}(P,P)$",
  "A control confirms that it is the rotation that \\emph{bounds} the growth,\nnot the preconditioning alone:\\ the naive $\\mathrm{blockdiag}(P,P)$"),

 ("M6-completeness-qualified",
  "the Coulomb--Sturmians soften it rather than escaping it because they are a\n\\emph{complete} set at one scale",
  "the Coulomb--Sturmians soften it rather than escaping it because they are a\nset complete at one scale \\emph{in the atomic metric} (a property the\nmolecular problem does not inherit --- see below)"),

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
