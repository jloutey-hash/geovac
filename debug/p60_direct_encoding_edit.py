"""Fold the direct-encoding construction into Paper 60: the second power of n is
no longer hypothetical."""
from pathlib import Path

P = Path("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex")
src = P.read_text(encoding="utf-8")
E = []

# ---- the table row is no longer a hypothetical -----------------------------
E.append((
    r"preconditioned, $G$ direct    & $125.4$  & $16.1$            & $2.0\times10^{3}$\\",
    r"preconditioned, $G$ direct    & $125.5$  & $24.3$            & $3.1\times10^{3}$\\"))

# ---- the verdict sentence + the construction -------------------------------
E.append((
    r"""the floor it would be $n$.  \textbf{So the lever is worth one power of $n$ as
priced here, and two if a direct block-encoding of $G$ is found.}  That is the
open item, and the size of the prize is visible in the waste:\ $\|G\|=0.372$
while the composed encoding carries $\alpha=3196$, a factor
$8.6\times10^{3}$ paid for nothing but the order of operations.""",
    r"""the floor it would be $n$.  The third row is not hypothetical:\ the direct
encoding exists, and the rest of this subsection constructs it.  \textbf{So the
lever is worth two powers of $n$, taking the metric penalty from $n^{3}$ to
$n$.}

\textbf{[SYMBOLIC + MEASURED]} The construction is available because
$G$'s symbol is a \emph{ratio} of two symbols that vanish to the same order.
Both $1-\sigma$ and $g$ have quadratic zeros at $\chi=\pi$, so the quotient
tends to $(kR)^2/24$ there, and at $\chi\to0$ it tends to $\tfrac14$;\ in the
variable $s=kR\cot(\chi/2)$,
\begin{equation}
\mathrm{ratio}(s)=\bigl(1-j_0(s)\bigr)\frac{s^{2}+(kR)^{2}}{4s^{2}},
\qquad \|\mathrm{ratio}\|_\infty=0.3716=\|G\|
\label{eq:ratio_symbol}
\end{equation}
---the sup coincides with $\|G\|$, so a circulant-embedded encoding of the
Toeplitz-minus-Hankel matrix $B$ built from this symbol's own cosine
coefficients carries $\alpha=0.372$ rather than the composed $3196$.  That
recovers exactly the factor the composition wasted.

$B$ is \emph{not} $G$:\ they differ by the finite-section commutator, measured
at $11.5$--$11.8\%$ in operator norm and \emph{not growing} with $n$.  The
point is that the difference does not matter, because $B$ enters only as a
whitening.  Any $X$ with $X^\dagger SX=I$ preserves the spectrum
(Eq.~\eqref{eq:amplitude_floor}), and taking $X=P^{-1/2}B^{-1/2}$ gives
$X^\dagger(I-C)X=B^{-1/2}GB^{-1/2}$, whose conditioning is
\emph{bounded}---$1.222,1.228,1.232,1.234$ at $n=20,40,80,160$, settling near
$1.23$---so one further $O(1)$-degree transformation absorbs it.  And
$\|X\|$ lands on the amplitude floor:\ $\|P^{-1/2}B^{-1/2}\|/\|(I-C)^{-1/2}\|$
$=1.0097,1.0049,1.0025,1.0013$ over the same range, tightening toward $1$.
\begin{center}
\begin{tabular}{crrrr}
\toprule
$n$ & $\mathrm{cond}(B)$ & $\|G-B\|/\|G\|$ & $\mathrm{cond}(B^{-1/2}GB^{-1/2})$ & $\|X\|/\text{floor}$\\
\midrule
$20$  & $2.198$ & $0.115$ & $1.222$ & $1.0097$\\
$80$  & $2.228$ & $0.118$ & $1.232$ & $1.0025$\\
$160$ & $2.229$ & $0.118$ & $\mathbf{1.234}$ & $\mathbf{1.0013}$\\
\bottomrule
\end{tabular}
\end{center}
So the whole metric factor is three explicitly-known pieces:\ a DST-I, a
diagonal computed from the index, and a degree-$\approx\!24$ QSVT on a
directly-encoded Toeplitz-minus-Hankel matrix whose subnormalization is
$O(1)$---none of it growing with basis size
(\texttt{tests/test\_paper60\_direct\_encoding.py})."""))

# ---- scope ------------------------------------------------------------------
E.append((
    r"""above under Table~\ref{tab:resource}'s cost model.  A direct block-encoding of
$G$ is not constructed here.  And by the block-structure result above, none of
this recovers $\ell$-selection, which no congruence can.""",
    r"""above under Table~\ref{tab:resource}'s cost model, and the direct construction
of Eq.~\eqref{eq:ratio_symbol}.  What is \emph{not} established is the circuit:\
the $O(\log^2N)$ sine transform is cited rather than compiled, and the
circulant embedding of $B$ is standard but not laid out here, so the counts
remain a resource model rather than a gate count.  And by the block-structure
result above, none of this recovers $\ell$-selection, which no congruence
can---the metric is now cheap to \emph{apply}, and the sparsity it destroys
stays destroyed."""))

for old, new in E:
    n = src.count(old)
    assert n == 1, f"anchor matched {n} times:\n{old[:110]}"
    src = src.replace(old, new)

P.write_text(src, encoding="utf-8")
print(f"{len(E)} edits applied, each matched exactly once.")
