"""REMEDIATION 4/5 -- body content corrections.

/qa paper_60 FULL 2026-09-12.  Seven findings, each verified by the PM against
primary text or by direct measurement before editing.

  A. tau OVER-SCOPE (code-A MATERIAL-2).  The paper gloss "the tau (DST-I)
     algebra --- which is to say Toeplitz minus Hankel, the structure of this
     section" is FALSE: tau membership needs the coefficient sequence to
     TERMINATE, and the chirp symbol's does not.  MEASURED by the PM: the DST-I
     leaves 4.12% / 1.98% / 1.10% off-diagonal on the cross block at n=16/64/160,
     against 2.3e-14 for a genuinely-tau tri(1,2,1).  The claim-matrix row
     already records the hypothesis gap flatly; the paper asserted the opposite.
     The corrected version is STRONGER for the paper, and identifies the residue.
  B. ORTHONORMALITY overstatement (claims-molecular F5) -- written by the PM on
     2026-09-12; as phrased it voids this paper's own molecular result.
  C. "leaves the growth untouched" (claims-molecular F4, code-A MATERIAL-1).
     MEASURED: 2766->42008 over a 16x range in N is exponent 0.98 against the
     raw column's 1.96, so the naive control HALVES the exponent.  "Unbounded"
     is the true word.  The control's own blindness is remediation 5.
  D. SCOPE paragraph drops "non-collinear" (claims-molecular F2a).
  E. the rank-(M-1) rotation stated unqualified for general M (F2b) -- the
     v5.11.0 wording at a locus the v5.11.4 scoping sweep did not reach.
  F. "stronger than either of the two above" unqualified (claims-molecular F7).
     MEASURED: gerade alpha = 1.13 = O(1) vs preconditioned alpha = 125.5 at
     n=160, so end-to-end the gerade lever wins ~111x where it applies.
  G. "the Coulomb-Sturmians evade it" (F6) -- they pay it at N^2 instead of
     N^6, as the same parenthesis says.
  H. "the one leg of this section that is derived" (claims-atomic M6) -- the
     section carries at least four derived legs, three of them tier-labelled.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "does not put our matrices inside it"

EDITS = [
 ("A-tau",
  """and for matrices in the $\\tau$ (DST-I)
algebra --- which is to say Toeplitz minus Hankel, the structure of this
section --- their result is an identity rather than an asymptotic:\\ the
eigenvalues are the ratio symbol sampled on the grid, with the weight
cancelling identically.  So the measurement below confirms a theorem instead of""",
  """and for matrices in the $\\tau$ (DST-I)
algebra their result is an identity rather than an asymptotic:\\ the
eigenvalues are the ratio symbol sampled on the grid, with the weight
cancelling identically.  \\emph{Our Toeplitz-minus-Hankel form does not put our
matrices inside it}:\\ $\\tau$ membership requires the coefficient sequence to
terminate, which the matching polynomial's does and the chirp symbol's does not,
and the DST-I correspondingly leaves $4.1\\%$/$2.0\\%$/$1.1\\%$ off-diagonal
weight on the cross block at $n=16$/$64$/$160$ against $10^{-14}$ for a
genuinely-$\\tau$ matrix.  The identity therefore covers the $\\tau$
\\emph{idealisation} of this problem rather than the problem, and the gap is
exactly the residue reported above:\\ the grid-sampled $\\tau$ prediction
reproduces $\\pi^2/24$ to $9\\times10^{-6}$ while the true $\\sigma_{\\max}$
departs by $0.26\\%$.  So the measurement below confirms a theorem instead of"""),

 ("B-orthonormal",
  """the basis is exactly orthonormal in the metric actually used, the statement this
paper measures independently as the intra-center block being the identity, so""",
  """the basis is exactly orthonormal \\emph{within each center} in the metric
actually used --- the statement this paper measures independently as the
intra-center block being the identity, and not a statement about the two-center
set, whose smallest generalized eigenvalue does collapse --- so"""),

 ("C-control",
  """without it leaves the growth untouched ($2766\\to42008$ over the same range).""",
  """without it leaves the growth \\emph{unbounded} ($2766\\to42008$ over the same
range --- an exponent of $0.98$, halved from the raw column's $1.96$ but still
growing, against the bounded $44$)."""),

 ("D-scope",
  """$s$-sector shared-scale bases at $M=2$ and $M=3$, together with the pricing""",
  """$s$-sector shared-scale bases at $M=2$ and \\emph{non-collinear} $M=3$ (the
collinear case is open and is not claimed), together with the pricing"""),

 ("E-rotation",
  """--- which is why the fixed, geometry-independent rank-$(M-1)$ rotation of
Sec.~\\ref{sec:resource} removes it, and why the coupling's participation ratio
stays small.""",
  """--- which is why the fixed, geometry-independent rank-$(M-1)$ rotation of
Sec.~\\ref{sec:resource} removes it at $M=2$ and for non-collinear $M=3$, and why
the coupling's participation ratio stays small."""),

 ("F-superlative",
  """A third conditioning lever exists, it is stronger than
either of the two above, and---unlike the gerade lever---it does not require
equivalent centers.""",
  """A third conditioning lever exists;\\ it is stronger than
either of the two above \\emph{in coverage and on the conditioning axis}, and
---unlike the gerade lever---it does not require equivalent centers.  It is not
uniformly cheaper end to end:\\ where the gerade lever applies its
subnormalization is $O(1)$ (the block's least eigenvalue is bounded below by the
same $2.555041\\ldots$ constant) against this lever's $O(n)$, so on
$\\mathrm{H}_2^+$ the gerade route remains the cheaper one."""),

 ("G-evade",
  """the Coulomb--Sturmians evade it because they are a""",
  """the Coulomb--Sturmians soften it rather than escaping it, because they are a"""),

 ("H-oneleg",
  """Eq.~\\eqref{eq:T0_closed} is the one leg of this section that is derived rather""",
  """Eq.~\\eqref{eq:T0_closed} is the one leg of this section's \\emph{scaling}
measurements that is derived rather"""),
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
        print("ABORT")
        return 2
    for name, old, new in EDITS:
        t = t.replace(old, new)
        print(f"  ok  {name}")
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print(f"applied: {len(EDITS)} body corrections")
    return 0


if __name__ == "__main__":
    sys.exit(main())
