"""DELTA citation NIT, verified at source by the PM: "no non-analytic
functions" over-states what Tao-McCurdy-Rescigno's basis contains.

Their odd-m DVR functions carry an explicit (xi^2-1)^{1/2} factor, introduced
in their own words "for the non-analytic behavior of the wave function at
xi = 1" (Sec. II.A, Eq. 11).  That is a one-electron angular-momentum boundary
condition at the focal axis -- the prolate-spheroidal analogue of the sin(theta)
behaviour ordinary spherical harmonics have for odd m -- and has nothing to do
with r12 or the electron-electron cusp, which is what the surrounding argument
is about.  The load-bearing claim is unaffected; the wording is not.

Verified independently by the PM against the accepted manuscript text before
applying (the reviewer's quote was not taken on trust).

Write-tool script file per memory rule feedback_no_heredoc_backslashes.
"""
import io
import sys

PATH = "papers/group2_quantum_chemistry/paper_12_algebraic_vee.tex"

OLD = r"""Rescigno~\cite{tao_mccurdy_rescigno2010}, working in these same
prolate spheroidal coordinates with a polynomial angular basis, no
$r_{12}$ factors and no non-analytic functions, reach
$-1.17442$~Ha---$0.05$~mHa from exact---at angular truncation
$l_{\max} = 6$."""

NEW = r"""Rescigno~\cite{tao_mccurdy_rescigno2010}, working in these same
prolate spheroidal coordinates with a polynomial angular basis and no
$r_{12}$ dependence anywhere in it, reach
$-1.17442$~Ha---$0.05$~mHa from exact---at angular truncation
$l_{\max} = 6$.  (Their odd-$m$ radial functions do carry an explicit
$(\xi^2-1)^{1/2}$, which they introduce for the non-analytic behaviour
at $\xi = 1$;\ that is a one-electron boundary condition on the focal
axis, the analogue of the $\sin\theta$ behaviour of ordinary spherical
harmonics at odd $m$, and is unrelated to $r_{12}$.  The point stands
as stated:\ nothing in their basis represents the \emph{interelectronic}
non-analyticity.)"""

with io.open(PATH, encoding="utf-8") as fh:
    text = fh.read()

if OLD not in text:
    print("FAILED TO MATCH")
    sys.exit(1)

with io.open(PATH, "w", encoding="utf-8") as fh:
    fh.write(text.replace(OLD, NEW, 1))

print("Paper 12: TMR basis description tightened to the r12-specific claim")
