"""Apply the contraction-seam capture to Paper 60.

Two paragraphs inserted after the chirp paragraph (eq:chirp_decay), before the
composition-wall paragraph:

  1. [SYMBOLIC + MEASURED]  the transcendental tagging owed under CLAUDE.md
     Sec. 4 -- the pi^2 of eq:sigma_law is truncation-side, the (2pi)^-1/2 and
     pi/4 of eq:chirp_decay are continuum-side; both calibration-tier M2.
  2. [MEASURED]  the weight-independence of the law, with its control and its
     honest scope (NOT V_0-independence).

Idempotent: refuses to run twice.
"""
from __future__ import annotations

import sys

PAPER = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"

ANCHOR = (
    "hypothesis covers this symbol class is left open.\n"
)

MARKER = "opposite sides of the compactness boundary"

NEW = r"""
\textbf{[SYMBOLIC + MEASURED]} The two transcendentals of this section sit on
opposite sides of the compactness boundary of Paper~18, and the classification
is not cosmetic:\ it tracks which of the two walls is removable.  The $\pi^2$
of Eq.~\eqref{eq:sigma_law} is \emph{truncation-side}.  It is the
Kac--Murdock--Szeg\H{o} constant $c_1$, which for a quadratic symbol zero is
the first Dirichlet eigenvalue of $-d^2/dx^2$ on the unit interval, and it is
present for every symbol with such a zero whether or not a Bessel function
appears anywhere in it.  In the Fock polar angle $\theta=\pi-\chi$ that is
directly visible:\ minimising $\langle\theta^2\rangle$ over
$\mathrm{span}\{\sin a\chi\}_{a\le n}$ --- a band-limited concentration problem
containing no symbol at all --- gives $n^2\min\langle\theta^2\rangle\to\pi^2$
(Richardson $9.86949$ against $\pi^2=9.86960$), the near-null direction of the
cross block realises that minimum with r.m.s.\ spread $\pi/n$ (Richardson
$3.14158$), and $1-\sigma_{\max}=(kR)^2\langle\theta^2\rangle/24$ holds on that
direction to $0.3\%$ at $n=320$ across $kR=0.5$--$4$.  The $\pi^2$ is thus the
price of the finite basis, not of the second center.  The $(2\pi)^{-1/2}$ and
the $\pi/4$ of Eq.~\eqref{eq:chirp_decay} are \emph{continuum-side}:\ the
normalisation and the branch phase of a Bessel asymptotic at the
$p\to\infty$ pole, which is the decompactified direction.  Under Paper~18's
taxonomy both are calibration-tier, sub-mechanism M2 --- the
$\pi^2\!\cdot\mathbb{Q}$ half of that ring for the first, the
$\sqrt\pi\cdot\mathbb{Q}$ half for the second.  The consequence is the one
Sec.~\ref{sec:resource} draws operationally:\ a truncation-side price is a
property of the matrix, which a preconditioner reaches, while a continuum-side
price is a property of the symbol, which no congruence of the finite section
can touch.

\textbf{[MEASURED]} Reading the symbol as a quotient yields a prediction, and
it survives its own control.  Let the metric be deformed by any smooth positive
radial weight $W$ on the Fock sphere --- the intra-center block becoming the
finite section of $W$ and the cross-center block that of $W\,j_0$.  The
generalised symbol is then $W j_0/W=j_0$, which does not see $W$ at all, so
both the exponent and the constant should survive unchanged.  They do:\ at
$kR=2$ and $n=160$ the collapse $(1-\sigma_{\max})(n/kR)^2$ reads $0.4072$,
$0.4123$, $0.4107$ and $0.4164$ for $W=1$, $1+\tfrac45\cos\chi$,
$2+\sin\chi$ and $e^{-\chi}$ against $\pi^2/24=0.4112$, with fitted exponents
between $-1.98$ and $-2.01$ throughout.  The control carries the claim:\ a
weight that \emph{vanishes} at the degeneracy, $W=1+\cos\chi$, moves the
constant to $0.828$ while leaving the exponent at $-1.97$, so the agreement
above is a measurement and not an insensitivity.  \emph{Scope, because the
natural overstatement is close by.}  This is not independence of the weighting
potential $V_0$:\ a position-space-local $V_0$ acts on momentum space by
convolution rather than multiplication and therefore leaves this class
entirely.  What is shown is narrower and still worth having --- the law is
carried by the translation phase, not by any momentum-space weight multiplying
it.  It also predicts what a metric \emph{outside} the class should do, and the
one measured case agrees:\ the $L^2$ overlap degrades with the same exponent
but a constant larger by $\approx\!1.4$, which is what a non-multiplication
metric is expected to look like.  Whether any position-space $V_0$ preserves
the constant is open.
"""


def main() -> int:
    with open(PAPER, encoding="utf-8") as fh:
        text = fh.read()
    if MARKER in text:
        print("ALREADY APPLIED -- marker present; nothing done.")
        return 1
    if text.count(ANCHOR) != 1:
        print(f"ANCHOR not unique (count={text.count(ANCHOR)}); aborting.")
        return 2
    text = text.replace(ANCHOR, ANCHOR + NEW)
    with open(PAPER, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"applied: inserted {len(NEW)} chars after the chirp paragraph")
    return 0


if __name__ == "__main__":
    sys.exit(main())
