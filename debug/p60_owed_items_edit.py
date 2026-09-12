"""Paper 60 owed items (2026-09-12): west_ruedenberg2013 removal + the
preconditioner lever."""
from pathlib import Path

P = Path("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex")
src = P.read_text(encoding="utf-8")
EDITS = []

# --- 1. west_ruedenberg2013 cannot support the principal-angle attribution ---
EDITS.append((
    r"angles~\cite{amos_hall1961,king1967,west_ruedenberg2013} between",
    r"angles~\cite{amos_hall1961,king1967} between"))

EDITS.append((
    r"""\bibitem{west_ruedenberg2013}
A.~C.~West, M.~W.~Schmidt, M.~S.~Gordon, and K.~Ruedenberg, ``A
comprehensive analysis of molecule-intrinsic quasi-atomic, bonding, and
correlating orbitals.\ I.\ Hartree--Fock wave functions,''
\textit{J.\ Chem.\ Phys.}\ \textbf{139}, 234107 (2013).

""", ""))

# --- 2. the preconditioner lever -------------------------------------------
EDITS.append((
    r"""\textbf{[MEASURED]} The large-$R$ lever is quantified alongside:""",
    r"""\textbf{[MEASURED]} A third conditioning lever exists, it is stronger than
either of the two above, and---unlike the gerade lever---it does not require
equivalent centers.  Because the ill-conditioning is a symbol \emph{zero} of
known order and location, it is removable by preconditioning in the sense of
Serra~\cite{serra1997}:\ a trigonometric polynomial sharing the zero bounds the
preconditioned spectrum in a fixed interval for every $n$.  The matching
polynomial here is $g(\chi)=2+2\cos\chi$, whose quadratic zero sits at
$\chi=\pi$ exactly where $1-\sigma$'s does;\ in this basis its Hankel part
vanishes identically ($n+m\ge2$ while $g_j=0$ for $j\ge2$), so the
preconditioner is \emph{exactly} the tridiagonal $\mathrm{tri}(1,2,1)$, which
the discrete sine transform diagonalizes in closed form (verified against the
closed-form eigenpairs to $9\times10^{-15}$).  Writing
$G=P^{-1/2}(I-C)P^{-1/2}$, the ungerade conditioning collapses and
\emph{stops growing}:
\begin{center}
\begin{tabular}{crrr}
\toprule
$n$ & $\mathrm{cond}(I+C)$ & $\mathrm{cond}(I-C)$ & $\mathrm{cond}(G)$\\
\midrule
$10$  & $2.383$ & $81.9$    & $2.096$\\
$40$  & $2.537$ & $1225.7$  & $2.219$\\
$160$ & $2.554$ & $19126.5$ & $\mathbf{2.229}$\\
\bottomrule
\end{tabular}
\end{center}
The whitening is unaffected in substance:\ any $X$ with $X^\dagger SX=I$ leaves
the generalized spectrum invariant, and $X=P^{-1/2}G^{-1/2}$ is such an $X$, so
the QSVT degree is set by $\mathrm{cond}(G)\to2.23$ rather than by
$\mathrm{cond}(S)$.  On the resource model of Table~\ref{tab:resource} that is
$d_{\rm inv}\approx16$, \emph{flat in basis size}, against $3\times10^{5}$ at
$n=160$ untreated.  This escapes the Bernstein $\Theta(\kappa)$ floor quoted
above rather than contradicting it:\ that floor constrains polynomial
approximation of $x^{-1/2}$ on $[\kappa^{-1},1]$, and preconditioning changes
the operator rather than the polynomial, so $x^{-1/2}$ is never approximated on
the bad interval.  The statement that the metric cost is ``mitigated, not
dissolved'' is therefore \emph{too pessimistic on the conditioning axis}, and
should be read as applying to the two axes below.

\textbf{[MEASURED]} What preconditioning does \emph{not} buy is locality, and
the reason is the second pole.  The symbol is pathological at both ends and
only one end is a zero:\ $P^{-1/2}$ cures $\chi=\pi$, while at $\chi\to0$ the
chirp survives untouched.  Measured as the bandwidth needed for a fixed
relative Frobenius accuracy, $S^{-1/2}$ requires a \emph{fixed fraction} of the
matrix ($b/n=0.72$ at $10^{-2}$, flat over $n=32$--$256$) whereas $G^{-1/2}$
requires $b=11\to17$ over the same range;\ but at $10^{-3}$ the advantage
erodes ($b=27\to141$), exactly as algebraic decay predicts.  The two are
distinguished cleanly by the profile exponent, which for $G^{-1/2}$ sits at
$-1.19$---the chirp's own $-5/4$, \emph{independent of $n$}---while for
$S^{-1/2}$ it drifts with $n$ ($-0.90\to-0.78$) as the $\chi=\pi$ singularity
sharpens.  So preconditioning removes one of the two pathologies exactly and
leaves the other exactly where it was.  Together with the block-structure
result above, the honest summary is that the metric's \emph{conditioning} cost
is removable, its \emph{locality} cost is capped by the chirp, and its
\emph{$\ell$-block} cost is not a conditioning effect at all.

\emph{Scope.} Measured on the homonuclear two-center $s$-sector symbol, for
which the parity blocks are $I\pm C$;\ whether the same construction reaches a
polyatomic block with symmetry-inequivalent centers (the
$\mathrm{H}_2\mathrm{O}$ $A_1$ case above) is untested here.  The end-to-end
resource claim additionally requires pricing the sine-transform circuit and a
block-encoding of $G$, neither of which is attempted in this paper;\ what is
established is the conditioning statement.

\textbf{[MEASURED]} The large-$R$ lever is quantified alongside:"""))

# --- 3. bibliography: Serra ------------------------------------------------
EDITS.append((
    r"""\bibitem{kac_murdock_szego1953}""",
    r"""\bibitem{serra1997}
S.~Serra, ``Preconditioning strategies for asymptotically ill-conditioned block
Toeplitz systems,'' \textit{BIT}\ \textbf{34}, 579 (1994);\ and ``On the
extreme eigenvalues of Hermitian (block) Toeplitz matrices,''
\textit{Linear Algebra Appl.}\ \textbf{270}, 109 (1998).

\bibitem{kac_murdock_szego1953}"""))

for old, new in EDITS:
    n = src.count(old)
    assert n == 1, f"anchor matched {n} times, expected 1:\n{old[:90]}"
    src = src.replace(old, new)

P.write_text(src, encoding="utf-8")
print(f"{len(EDITS)} edits applied, each matched exactly once.")
