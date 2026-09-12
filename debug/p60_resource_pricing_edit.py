"""Fold the end-to-end resource pricing into Paper 60."""
from pathlib import Path

P = Path("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex")
src = P.read_text(encoding="utf-8")

OLD = r"""\emph{Scope.} What is established is the \emph{conditioning} statement, on
$s$-sector shared-scale bases at $M=2$ and $M=3$.  The end-to-end resource claim
additionally requires pricing the sine-transform circuit and a block-encoding of
$G$, neither of which is attempted here;\ and by the block-structure result
above, none of this recovers $\ell$-selection, which no congruence can."""

NEW = r"""\textbf{[RESOURCE MODEL]} Priced end to end the lever is smaller than the
conditioning table suggests, and saying so precisely is more useful than the
raw ratio.  Two costs must be tracked rather than one:\ the QSVT degree
$d_{\rm inv}$, and the \emph{subnormalization} $\alpha$ of the block-encodings,
which multiplies when encodings are composed.

The amplitude is not negotiable.  Any $X$ with $X^\dagger SX=I$ satisfies
$X=S^{-1/2}U$ for some unitary $U$, so
\begin{equation}
\|X\|=\|S^{-1/2}\|\quad\text{exactly, independently of how $X$ is factored}
\label{eq:amplitude_floor}
\end{equation}
(verified to $6\times10^{-13}$ over $n=20$--$160$).  \emph{No factorization can
lower the amplitude floor}, and the untreated route already attains it.  What
preconditioning buys is depth, and only depth.  The factor $P^{-1/2}$ is exactly
implementable and costs no calls to any encoding of the metric:\ the DST-I has an
$O(\log^2N)$ quantum circuit~\cite{klappenecker2001}, and the diagonal
$\lambda_k^{-1/2}=\bigl(2+2\cos\tfrac{k\pi}{n+1}\bigr)^{-1/2}$ is computed from
the index.
\begin{center}
\begin{tabular}{lccc}
\toprule
route & $\alpha$ & $d_{\rm inv}$ & $\alpha\,d_{\rm inv}$\\
\midrule
untreated                     & $125.4$  & $3.1\times10^{5}$ & $3.9\times10^{7}$\\
preconditioned, $G$ composed  & $3196$   & $16.1$            & $5.1\times10^{4}$\\
preconditioned, $G$ direct    & $125.4$  & $16.1$            & $2.0\times10^{3}$\\
\bottomrule
\end{tabular}
\end{center}
\noindent(at $n=160$, $kR=2$, $\epsilon_E=1.6$~mHa).  In exponents, which is the
honest form:\ untreated $\alpha\sim n$ and $d_{\rm inv}\sim n^{2}$, so the
product scales as $n^{3}$;\ preconditioned with $G$ obtained by \emph{composing}
$P^{-1/2}$ with $I-C$, the degree goes flat but $\alpha$ inherits
$\|P^{-1/2}\|^{2}\sim n^{2}$, giving $n^{2}$;\ with a direct encoding of $G$ at
the floor it would be $n$.  \textbf{So the lever is worth one power of $n$ as
priced here, and two if a direct block-encoding of $G$ is found.}  That is the
open item, and the size of the prize is visible in the waste:\ $\|G\|=0.372$
while the composed encoding carries $\alpha=3196$, a factor
$8.6\times10^{3}$ paid for nothing but the order of operations.

\emph{Scope.} What is established is the \emph{conditioning} statement, on
$s$-sector shared-scale bases at $M=2$ and $M=3$, together with the pricing
above under Table~\ref{tab:resource}'s cost model.  A direct block-encoding of
$G$ is not constructed here.  And by the block-structure result above, none of
this recovers $\ell$-selection, which no congruence can."""

assert src.count(OLD) == 1, "scope anchor not found"
src = src.replace(OLD, NEW)

BIB_ANCHOR = r"""\bibitem{serra1997}"""
BIB_NEW = r"""\bibitem{klappenecker2001}
A.~Klappenecker and M.~R\"otteler, ``Discrete cosine transforms on quantum
computers,'' in \textit{Proc.\ IEEE R8-EURASIP Symp.\ Image and Signal
Processing and Analysis (ISPA)}, Pula, Croatia, pp.\ 464--468 (2001);
arXiv:quant-ph/0111038.

\bibitem{serra1997}"""
assert src.count(BIB_ANCHOR) == 1
src = src.replace(BIB_ANCHOR, BIB_NEW)

P.write_text(src, encoding="utf-8")
print("paper 60: resource pricing + klappenecker2001 bibitem")
