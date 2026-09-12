"""Upgrade the -2.90250 withdrawal from "could not be confirmed" to a bound.

Because T' is a matrix of pure numbers, a sub-family's secular matrix is exactly
the principal submatrix of M (verified bit-identically).  Cauchy interlacing then
caps the largest root, and since E = -p^2/2 with p > 0, no sub-family can beat
the family it is drawn from.  So the cited figure cannot have come from the
locked-scale posing at any K or any selection -- a proof rather than a
comparison of numbers.
"""
import io

PAP = "papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex"
REG = "debug/qa/numeric_registry.py"

# ------------------------------------------------------------------ registry
reg = io.open(REG, encoding="utf-8").read()
NEW = '''    "p60_best102_locked": dict(
        value=7.25, convention="constant: mHa above the exact He ground state "
                               "reached by the BEST 102 configurations (ranked by "
                               "ground-state weight) drawn from the K=244 pool, "
                               "locked-scale posing -- the sharpest selection test "
                               "of Avery's '102 optimized configurations'",
        q=None,
        provenance="MEASURED 2026-09-08, debug/p60_avery_102_probe.py. ABOVE the "
                   "pool's own 7.057 mHa, as Cauchy interlacing requires: a "
                   "principal submatrix cannot have a larger top eigenvalue. 200 "
                   "random 102-subsets reach only 903 mHa. The cited Avery figure "
                   "is 1.224 mHa, unreachable in this posing at any K.",
        aliases={7.057: "the full K=244 pool", 903.4: "best of 200 random 102-subsets"}),
'''
anchor = '    "p60_cond_S_converged": dict('
assert anchor in reg and "p60_best102_locked" not in reg
reg = reg.replace(anchor, NEW + anchor, 1)
io.open(REG, "w", encoding="utf-8").write(reg)
print("registry: p60_best102_locked added")

# --------------------------------------------------------------------- paper
tex = io.open(PAP, encoding="utf-8").read()

OLD = r"""\textbf{[OBSERVATION]} A comparison to the $-2.90250$~Ha reached with $102$
optimized Coulomb--Sturmian configurations~\cite{avery2006} is \emph{withdrawn}
as a characterization of the present method:  that figure ($1.2$~mHa short)
matches our scale-optimized ladder ($1.64$~mHa short at $K=100$) rather than our
locked-scale one ($7.70$~mHa at the same $K$), and which posing produced it could
not be confirmed against the primary source.  The machinery is correct, including
the higher-$\ell$ angular correlation (Gaunt / Wigner-$3j$ multipole coupling)."""

NEW_T = r"""\textbf{[INTERNAL THEOREM]} \emph{Selection cannot rescue the locked
posing.}  Because $T'$ is a matrix of pure numbers, its entries do not depend on
which \emph{other} configurations are present, so the secular matrix of any
sub-family $A$ is \emph{exactly} the corresponding principal submatrix of $M$
(verified bit-identically).  Cauchy interlacing then gives
$\lambda_{\max}(M_A)\le\lambda_{\max}(M)$, and since $E=-p_\kappa^{2}/2$ with
$p_\kappa=\lambda_{\max}>0$,
\begin{equation}
  E(A)\;\ge\;E(M)\qquad\text{for every sub-family } A\subseteq M,
  \label{eq:no_selection}
\end{equation}
the families being nested in $n_{\max}$ so that $E$ decreases monotonically along
the ladder.  No choice of configurations, however optimized, beats the family it
is drawn from.  \textbf{[MEASURED]} The best $102$ by ground-state weight out of
a $K=244$ pool returns $\gvq{p60_best102_locked}{7.25}$~mHa against that pool's
own $7.06$~mHa --- \emph{above} it, as Eq.~\eqref{eq:no_selection} requires ---
while $200$ random $102$-subsets reach only $903$~mHa.  Note what supplies the
bound:\ it is the same pure-number property that makes the encoding attractive.
Were $T'$ basis-dependent, sub-families would not be principal submatrices and
interlacing would not apply.

\textbf{[OBSERVATION]} A comparison to the $-2.90250$~Ha reached with $102$
optimized Coulomb--Sturmian configurations~\cite{avery2006} is therefore
\emph{withdrawn}, and on stronger grounds than a mismatch of numbers:\ by
Eq.~\eqref{eq:no_selection} that figure ($1.2$~mHa short) \emph{cannot} have come
from the locked-scale Goscinskian posing, at any $K$ or under any selection ---
our largest computed family stands at $6.82$~mHa and the ladder is monotone.  It
matches instead our scale-optimized ladder ($1.64$~mHa at $K=100$).  Which posing
produced it --- or which basis, the genuine shared-scale Coulomb--Sturmians being
the other candidate --- could not be confirmed against the primary source, and we
flag it as the single result that would overturn this section were it to prove
otherwise.  The machinery is correct, including the higher-$\ell$ angular
correlation (Gaunt / Wigner-$3j$ multipole coupling)."""

assert OLD in tex, "withdrawal locus not found"
tex = tex.replace(OLD, NEW_T, 1)
io.open(PAP, "w", encoding="utf-8").write(tex)
print("paper: interlacing theorem added; withdrawal upgraded to a bound")
