"""Correct the transcendental-tagging paragraph after C23 run #3.

Three findings, all verified locally before editing (C23 hard rule):

  T1/T2 PRIOR ART.  The reading of c_1 = pi^2 as a Dirichlet eigenvalue, and
  the band-limited second-moment extremal problem, are BOTH in Boettcher-Widom
  -- the source this paper already cites two paragraphs earlier, whose TITLE
  names the Wirtinger-Sobolev inequality that carries the constant, and which
  itself observes that c_alpha is independent of b.  Verified here: the bibitem
  is at line 1695 and the title reads "From Toeplitz eigenvalues through
  Green's kernels to higher-order Wirtinger-Sobolev inequalities".  So the
  Bessel-free measurement is the b == 1 case of the same theorem -- a change of
  REPRESENTATION, not an independent route, and the paragraph must stop
  presenting it as independent evidence.

  T3 -- THE SERIOUS ONE.  The removability sentence names the wrong mechanism
  and its second half is false.  Verified at line 1084: the preconditioner's
  matching polynomial is chosen to SHARE THE SYMBOL'S ZERO (Serra), so the
  lever is symbol-side, not truncation-side; and preconditioning IS a
  congruence of the finite section, one that changes the symbol to f/g, so
  "no congruence can touch a property of the symbol" is false as written.

Also flags the internal tension the scan surfaced: our own Toeplitz scan memo
recorded both "the conditioning and the non-locality are one fact" and "they
are two independent facts" at different points, and this paragraph's split
adopts the second.

NOT changed: the KMS inline attribution, reported as a C20 defect.  Checked --
the bibitem exists at line 1691 and the sentence cites it; the flagged phrase
is a second mention inside the SAME sentence that already carries the \\cite.

Idempotent.
"""
from __future__ import annotations

import sys

P = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
MARKER = "What this does not license"

OLD_START = "\\textbf{[SYMBOLIC + MEASURED]} The two transcendentals of this section sit on"
OLD_END = ("price is a property of the symbol, which no congruence of the finite section\n"
           "can touch.\n")

NEW = r"""\textbf{[SYMBOLIC + MEASURED]} The two transcendentals of this section are
tagged here against Paper~18, and the tagging is a statement of
\emph{provenance} --- where each constant comes from --- nothing more.  The
$\pi^2$ of Eq.~\eqref{eq:sigma_law} is the Kac--Murdock--Szeg\H{o} constant
$c_1$, and its reading as an eigenvalue is already in the source cited above:\
B\"ottcher and Widom~\cite{bottcher_widom2005} present $c_\alpha$ as the least
eigenvalue of the boundary-value problem $(-1)^\alpha u^{(2\alpha)}=\lambda u$
on $[0,1]$ with $u^{(j)}(0)=u^{(j)}(1)=0$, which at $\alpha=1$ is the Dirichlet
problem for $-d^2/dx^2$ and gives $\pi^2$;\ their title names the extremal
(Wirtinger--Sobolev) inequality that carries it, and the independence of
$c_\alpha$ from $b$ is likewise their observation, not ours.  The constant is
therefore fixed by the order of the symbol's zero together with the finite
section, and carries no Bessel content:\ nothing about $j_0$ beyond
$b(1)=(kR)^2/24$ enters it.  We record the same statement in the basis's own
variable, where it is directly visible --- minimising
$\langle\theta^2\rangle$ over $\mathrm{span}\{\sin a\chi\}_{a\le n}$ in the
Fock polar angle $\theta=\pi-\chi$ gives $n^2\min\langle\theta^2\rangle\to\pi^2$
(Richardson $9.86949$ against $\pi^2=9.86960$), the near-null direction attains
that minimum with r.m.s.\ spread $\pi/n$ (Richardson $3.14158$), and
$1-\sigma_{\max}=(kR)^2\langle\theta^2\rangle/24$ holds on it to $0.3\%$ at
$n=320$ across $kR=0.5$--$4$.  \emph{This is a change of representation and not
an independent route:\ it is the $b\equiv1$ case of the same theorem, and the
minimiser is the Dirichlet ground state in the band index.}  The
$(2\pi)^{-1/2}$ and the $\pi/4$ of Eq.~\eqref{eq:chirp_decay} have a different
provenance --- the normalisation and the branch phase of a Bessel asymptotic at
the $p\to\infty$ pole.  Under Paper~18's taxonomy both constants are
calibration-tier, sub-mechanism M2:\ the $\pi^2\!\cdot\mathbb{Q}$ half of that
ring for the first, the $\sqrt\pi\cdot\mathbb{Q}$ half for the second.

\emph{What this does not license.}  An earlier version of this paragraph read
the provenance as \emph{predicting} removability --- the first price being ``a
property of the matrix, which a preconditioner reaches'' and the second ``a
property of the symbol, which no congruence can touch''.  Both halves are
wrong.  The preconditioner of Sec.~\ref{sec:resource} is built \emph{from the
symbol}:\ its matching polynomial is chosen precisely to share the symbol's
zero~\cite{serra1997}.  And preconditioning \emph{is} a congruence of the
finite section, one that replaces the symbol by $f/g$.  The defensible
statement concerns the symbol on both sides and turns on the \emph{kind} of
feature rather than on the origin of its constant:\ a banded congruence
multiplies the symbol by a trigonometric polynomial, which can cancel a zero of
finite order but cannot alter a decay or smoothness class.  The conditioning
pole is a zero of order two;\ the locality pole is the chirp's $j^{-5/4}$
envelope.  That asymmetry, and not the provenance of $\pi^2$, is why the first
yields and the second does not.  \emph{Open, and load-bearing for the split
adopted here:}\ whether the two poles are genuinely independent features or two
faces of one object is not settled --- our own reading has gone both ways ---
and this paragraph takes the weaker, measured position.
"""

BIB_OLD = """A.~B\\"ottcher and H.~Widom, ``From Toeplitz eigenvalues through Green's
kernels to higher-order Wirtinger--Sobolev inequalities,''
arXiv:math/0412269 (2004).
"""
BIB_NEW = """A.~B\\"ottcher and H.~Widom, ``From Toeplitz eigenvalues through Green's
kernels to higher-order Wirtinger--Sobolev inequalities,'' in
\\textit{The Extended Field of Operator Theory}, Operator Theory: Advances and
Applications, Birkh\\"auser, Basel (2007), pp.~73--87;\\
doi:10.1007/978-3-7643-7980-3\\_4;\\ arXiv:math/0412269 (2004).
"""


def main() -> int:
    with open(P, encoding="utf-8") as fh:
        t = fh.read()
    if MARKER in t:
        print("ALREADY APPLIED")
        return 1
    i = t.find(OLD_START)
    if i < 0:
        print("start anchor missing")
        return 2
    j = t.find(OLD_END, i)
    if j < 0:
        print("end anchor missing")
        return 3
    t = t[:i] + NEW + t[j + len(OLD_END):]
    if t.count(BIB_OLD) != 1:
        print(f"bib anchor count={t.count(BIB_OLD)}; body applied, bib skipped")
    else:
        t = t.replace(BIB_OLD, BIB_NEW)
    with open(P, "w", encoding="utf-8") as fh:
        fh.write(t)
    print("applied: tagging paragraph corrected + bibitem given published coords")
    return 0


if __name__ == "__main__":
    sys.exit(main())
