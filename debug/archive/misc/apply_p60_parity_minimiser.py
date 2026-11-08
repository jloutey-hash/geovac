"""Make the minimiser claim precise: it carries the antipodal parity factor.

Paper 60 now asserts "the minimiser is the Dirichlet ground state in the band
index" (added earlier today when the independent-route claim was withdrawn).
That sentence had NO backing test -- a coverage gap under the claim->artifact
rule, created by this session.

Closing it turned up that the claim is true only up to an alternation.  The
basis index is the chi-index and the Dirichlet mode is natural in the
theta-index, with sin(a chi) = (-1)^(a+1) sin(a theta).  Measured:

    corr(band minimiser, sin(pi a/(n+1)))            = 0.0000001  at n = 320
    corr(band minimiser, (-1)^(a+1) sin(pi a/(n+1))) = 0.9999998  at n = 320

i.e. without the parity factor the two are EXACTLY ORTHOGONAL.  The SW
near-null direction matches the alternating mode to 0.9999993 at n = 320.

This is the same antipodal parity that C23 run #2 flagged as the caution on the
contraction reading (n^-1 U_{n-1}(cos(pi - z/n)) -> (-1)^(n+1) j0(z)) -- a
second, independent appearance of the same factor.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "antipodal parity"

OLD = """\\emph{This is a change of representation and not
an independent route:\\ it is the $b\\equiv1$ case of the same theorem, and the
minimiser is the Dirichlet ground state in the band index.}"""

NEW = """\\emph{This is a change of representation and not
an independent route:\\ it is the $b\\equiv1$ case of the same theorem, and the
minimiser is the Dirichlet ground state in the band index.}  It is worth
writing that identification out, because it carries the \\emph{antipodal
parity} and is false without it:\\ the basis index is $\\chi$ while the
Dirichlet mode is natural in $\\theta$, and
$\\sin(a\\chi)=(-1)^{a+1}\\sin(a\\theta)$, so the minimiser is
$c_a\\propto(-1)^{a+1}\\sin\\bigl(\\pi a/(n+1)\\bigr)$.  Measured at $n=320$, that
mode correlates with the band minimiser to $1-2\\times10^{-7}$ and with the
near-null direction of the cross block to $1-7\\times10^{-7}$, while the
unalternated $\\sin(\\pi a/(n+1))$ is \\emph{exactly orthogonal} to both.  The
same factor is the caution attached to the contraction reading of $j_0$, where
the antipodal Mehler--Heine limit converges along parities rather than
outright."""


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    if t.count(OLD) != 1:
        print(f"anchor count={t.count(OLD)}; aborting")
        return 2
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t.replace(OLD, NEW))
    print("applied: minimiser claim made parity-precise")
    return 0


if __name__ == "__main__":
    sys.exit(main())
