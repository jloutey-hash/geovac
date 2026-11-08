"""Paper 60 remediation (2026-09-11): KMS attribution + Halmos rescope + the
l-block/conditioning separation.  Written as a script per the no-heredoc-
backslashes rule -- every replacement is checked to fire exactly once."""
from pathlib import Path

P = Path("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex")
src = P.read_text(encoding="utf-8")

EDITS = []

# --- 1. Prior-art attribution for eq:sigma_law -------------------------------
EDITS.append((
    "law read on different windows.  The same singular spectrum carries the",
    r"""law read on different windows.

\textbf{[PRIOR ART]} Eq.~\eqref{eq:sigma_law} is not new \emph{as an
asymptotic law}, and we claim only the identification.  It is the
Kac--Murdock--Szeg\H{o} extreme-eigenvalue
asymptotic~\cite{kac_murdock_szego1953}:\ writing the relevant symbol in the
normal form $|1-t|^{2\alpha}b(t)$ used by B\"ottcher and
Widom~\cite{bottcher_widom2005}, one has
$\lambda_{\min}\sim(c_\alpha/n^{2\alpha})\,b(1)$, with $c_1=\pi^2$ due to Kac,
Murdock and Szeg\H{o}.  Our symbol is the case $\alpha=1$:\ with
$\theta=\pi-\chi$, $1-j_0(kR\cot(\chi/2))=(kR)^2\theta^2/24+O(\theta^4)$, so
$b(1)=(kR)^2/24$ and $c_1b(1)/n^2$ reproduces Eq.~\eqref{eq:sigma_law}
constant and all;\ the conditioning exponent $2$ is the order of the symbol's
maximum and nothing else.  What is ours is the \emph{identification} --- that
the two-center Shibuya--Wulfman metric in the sine basis is such a finite
section, Toeplitz minus Hankel since $\langle n|a|m\rangle=c_{n-m}-c_{n+m}$,
with symbol $j_0(kR\cot(\chi/2))$;\ equivalently that the Shibuya--Wulfman
operator is multiplication by the translation phase
$e^{i\mathbf{p}\cdot\mathbf{R}}$ on the Fock sphere, whose angular average
$j_0(pR)$ becomes trivial at $p=0$.  The physical reading follows:\ a basis
complete enough to carry wavelengths longer than $R$ cannot tell the two
centers apart.  Equivalently, $1-\sigma_{\max}=\tfrac16(R/L_{\max})^2$ with
$L_{\max}=2n/(\pi k)$ the longest wavelength the truncated basis carries, so
the degeneracy switches on exactly when $L_{\max}$ exceeds the bond length ---
which a complete basis must eventually do.

\textbf{[MEASURED]} One residue is worth stating rather than hiding.  Our
symbol does \emph{not} satisfy the smoothness hypothesis under
which~\cite{bottcher_widom2005} proves the constant:\ at the opposite end
$\chi\to0$ the symbol is a chirp (amplitude $\sim\chi$, phase $\sim2kR/\chi$)
whose Fourier coefficients decay only as $|c_j|\sim j^{-5/4}$, so
$\sum_jj|c_j|$ diverges.  The constant nevertheless holds.  Relatedly, the
$\sim\!1\%$ figure quoted above is the asymptotic's own $O(1/n)$ term rather
than scatter --- the relative residue halves under each doubling of $n$,
falling $0.134\to0.010$ across $n=10\to160$ at $kR=2$.  Whether a weaker
hypothesis covers this symbol class is left open.

The same singular spectrum carries the"""))

# --- 2. Halmos: canonical form is his; the norm identity is a consequence ----
EDITS.append((
    r"is $\|[P_A,P_B]\|=\max_k\sigma_k\sqrt{1-\sigma_k^{2}}$~\cite{halmos1969,bottcher_spitkovsky2010}, which \emph{saturates}",
    r"""is $\|[P_A,P_B]\|=\max_k\sigma_k\sqrt{1-\sigma_k^{2}}$ --- a one-line
consequence of the two-subspaces canonical
form~\cite{halmos1969,bottcher_spitkovsky2010}, which resolves the pair into
$2\times2$ blocks at the principal angles $\theta_k$, in which the commutator
has norm $\sin\theta_k\cos\theta_k$ --- which \emph{saturates}"""))

# --- 3. The l-block loss is NOT a functional of the sigma spectrum -----------
EDITS.append((
    r"""severity information in $1-\sigma_{\max}$
(\texttt{tests/test\_paper60\_sigma\_law.py}).""",
    r"""severity information in $1-\sigma_{\max}$
(\texttt{tests/test\_paper60\_sigma\_law.py}).

\textbf{[SYMBOLIC]} A third cost does \emph{not} live in this spectrum, and
separating it strengthens rather than weakens the obstruction.  If $X$ is
invertible and block diagonal with respect to
$\mathcal{H}=\bigoplus_\ell\mathcal{H}_\ell$, and $X^\dagger SX$ is block
diagonal, then $S=X^{-\dagger}(X^\dagger SX)X^{-1}$ is block diagonal too,
being a product of block-diagonal factors.  Contrapositively:\ if $S$ is not
$\ell$-block diagonal, \emph{no} block-diagonal congruence --- L\"owdin,
canonical, or Cholesky --- orthogonalizes it.  The two-center metric couples
$\ell$ while preserving $m$ (axial symmetry), so $m$-selection survives
orthogonalization and within-$m$ $\ell$-selection cannot, at every condition
number $>1$:\ the loss is already present for arbitrarily weak inter-center
coupling and does not relax as $\mathrm{cond}(S)\to1^{+}$.  It is therefore
independent of Eq.~\eqref{eq:sigma_law}, and the sparsity-destroying
orthogonalization excluded below is excluded \emph{structurally}, not on
conditioning grounds."""))

# --- 4. Bibliography --------------------------------------------------------
EDITS.append((
    r"""\bibitem{halmos1969}""",
    r"""\bibitem{kac_murdock_szego1953}
M.~Kac, W.~L.~Murdock, and G.~Szeg\H{o}, ``On the eigen-values of certain
Hermitian forms,'' \textit{J.\ Rational Mech.\ Anal.}\ \textbf{2}, 767 (1953).

\bibitem{bottcher_widom2005}
A.~B\"ottcher and H.~Widom, ``From Toeplitz eigenvalues through Green's
kernels to higher-order Wirtinger--Sobolev inequalities,''
arXiv:math/0412269 (2004).

\bibitem{halmos1969}"""))

for old, new in EDITS:
    n = src.count(old)
    assert n == 1, f"anchor matched {n} times, expected 1:\n{old[:90]}"
    src = src.replace(old, new)

P.write_text(src, encoding="utf-8")
print(f"{len(EDITS)} edits applied, each matched exactly once.")
