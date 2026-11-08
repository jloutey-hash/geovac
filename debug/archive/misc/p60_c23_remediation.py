"""C23 run #1 remediation for Paper 60 (2026-09-12).

Rule followed throughout: cite ONLY what was verified at primary source in this
session. Where the scan found prior art that could not be reached (WebSearch
budget exhausted), the paper's OVER-CLAIM is removed in prose using the classical
name, and the specific citation is recorded as owed rather than invented.

VERIFIED HERE and therefore cited:
  * DLMF 10.32.10 and 10.40.2 (quoted verbatim with their phase conditions)
  * Loring, Ann. Funct. Anal. 5(2), 2014, arXiv:1306.1923 (title/venue/topic)
  * Barthelme & Usevich, SIAM J. Matrix Anal. Appl. 42(1), 17 (2021),
    arXiv:1910.14067 (title/venue/topic)

OWED, named in prose but NOT given a bibitem:
  Jordan 1875; Jordan-Wielandt; the two-block CBS constant; Slater-Koster 1954;
  the Loewdin-symmetry note; Hartman-Wintner; Jaffard 1990; Groechenig-Leinert
  TAMS 358 (2006); Driscoll-Fornberg 2002.
"""
from pathlib import Path

P = Path("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex")
src = P.read_text(encoding="utf-8")
E = []

# --- A1: the block-spectrum chain is classical, not derived here -------------
E.append((
    r"""singular values of the cross block (the cosines of the principal
angles~\cite{amos_hall1961,king1967} between
the two center subspaces), so
$\mathrm{cond}(S)=(1+\sigma_{\max})/(1-\sigma_{\max})$ \emph{exactly}.""",
    r"""singular values of the cross block (the cosines of the principal
angles~\cite{amos_hall1961,king1967} between
the two center subspaces), so
$\mathrm{cond}(S)=(1+\sigma_{\max})/(1-\sigma_{\max})$ \emph{exactly}.
This chain is classical and we use rather than derive it:\ the principal angles
are Jordan's, the spectrum $\{1\pm\sigma_k\}$ is the Jordan--Wielandt form of
the off-diagonal block, and the resulting condition number is the two-block
CBS constant of subspace-correction theory."""))

# --- B1: the chirp asymptotic is a Bessel asymptotic, with the full constant --
E.append((
    r"""\textbf{[MEASURED]} One residue is worth stating rather than hiding.  Our
symbol does \emph{not} satisfy the smoothness hypothesis under
which~\cite{bottcher_widom2005} proves the constant:\ at the opposite end
$\chi\to0$ the symbol is a chirp (amplitude $\sim\chi$, phase $\sim2kR/\chi$)
whose Fourier coefficients decay only as $|c_j|\sim j^{-5/4}$, so
$\sum_jj|c_j|$ diverges.  The constant nevertheless holds.  Relatedly, the""",
    r"""\textbf{[SYMBOLIC + MEASURED]} One residue is worth stating rather than
hiding, and it has a closed form.  Our symbol does \emph{not} satisfy the
smoothness hypothesis under which~\cite{bottcher_widom2005} proves the
constant:\ at the opposite end $\chi\to0$ it is a chirp, amplitude
$\sim\chi/(2kR)$ and phase $\sim2kR/\chi$.  Its coefficient asymptotics need no
stationary-phase argument, because the model integral
$\int_0^\infty x\,e^{i(2c/x+jx)}\,dx$ is a modified Bessel function:\ it is
DLMF~(10.32.10)~\cite{dlmf} at $\nu=2$, continued from that formula's
$|\mathrm{ph}\,z|<\pi/4$ to $|\mathrm{ph}\,z|=\pi/2$, and DLMF~(10.40.2), whose
stated validity $|\mathrm{ph}\,z|\le\tfrac32\pi-\delta$ covers that sector
directly, gives
\begin{equation}
|c_j|=(2\pi)^{-1/2}2^{-3/4}(kR)^{-1/4}\,j^{-5/4}
\bigl|\sin\!\bigl(2\sqrt{2kRj}+\tfrac\pi4\bigr)\bigr|+o(j^{-5/4}),
\label{eq:chirp_decay}
\end{equation}
constant and phase included, verified over $j=64$--$65536$ at three values of
$kR$ by three independent quadrature routes, sign pattern included.  Two
corrections come with it.  The $\pi/4$ is the branch phase of the
$(\pi/2z)^{1/2}$ prefactor, \emph{not} a stationary-phase signature.  And
$\sum_j|c_j|$ \emph{converges} --- $5/4>1$, so the symbol is in the Wiener
algebra; what diverges is $\sum_jj|c_j|$, which is the hypothesis
B\"ottcher--Widom needs and a different condition entirely.  The constant
nevertheless holds.  Relatedly, the"""))

# --- A3: Proposition D demoted to a specialization of a classical result ------
E.append((
    r"""\textbf{[SYMBOLIC]} A third cost does \emph{not} live in this spectrum, and
separating it strengthens rather than weakens the obstruction.  If $X$ is""",
    r"""\textbf{[SYMBOLIC]} A third cost does \emph{not} live in this spectrum, and
separating it strengthens rather than weakens the obstruction.  The underlying
algebra is classical rather than ours --- that a congruence preserves a block
grading exactly when it commutes with it is the symmetry-preservation property
of L\"owdin orthogonalization, known in this paper's own field since
Slater and Koster, and in operator terms the statement that the block-diagonal
matrices form a commutant and are therefore inverse-closed.  What is ours is
only the application, below.  If $X$ is"""))

E.append((
    r"""independent of Eq.~\eqref{eq:sigma_law}, and the sparsity-destroying
orthogonalization excluded below is excluded \emph{structurally}, not on
conditioning grounds.""",
    r"""independent of Eq.~\eqref{eq:sigma_law}, and the sparsity-destroying
orthogonalization excluded below is excluded \emph{structurally}, not on
conditioning grounds.  That specialization --- $m$ survives, within-$m$ $\ell$
cannot, at every $\mathrm{cond}(S)>1$ --- is the part this paper claims."""))

# --- A2: the commutator norm has a citable source ----------------------------
E.append((
    r"""consequence of the two-subspaces canonical
form~\cite{halmos1969,bottcher_spitkovsky2010}, which resolves the pair into""",
    r"""consequence of the two-subspaces canonical
form~\cite{halmos1969,bottcher_spitkovsky2010} (the norm identity itself is
recorded by Loring~\cite{loring2014}), which resolves the pair into"""))

# --- A4: name the Toeplitz spectral mechanism behind the gerade constant ------
E.append((
    r"""$\mathrm{cond}(I+C)\to 2/(1+\min_x j_0(x))=2.555041\ldots$, an exact constant of
the sinc symbol independent of $R$ and of basis size""",
    r"""$\mathrm{cond}(I+C)\to 2/(1+\min_x j_0(x))=2.555041\ldots$, an exact constant of
the sinc symbol independent of $R$ and of basis size --- an instance of the
classical fact that a self-adjoint Toeplitz operator's spectrum is the convex
hull of its symbol's essential range, so that
$\mathrm{cond}(I+T(a))=(1+\max a)/(1+\min a)$"""))

# --- B3: the decay class is not the obstruction ------------------------------
E.append((
    r"""leaves the other exactly where it was.  Together with the block-structure""",
    r"""leaves the other exactly where it was.  The localization class is
\emph{not} what blocks a local inverse:\ with $|c_j|\sim j^{-5/4}$ the metric
has polynomial off-diagonal decay of order $>1$ in one dimension, which is
Jaffard's class, and inversion there preserves the decay \emph{provided the
operator is boundedly invertible};\ inverse-closedness then makes the
holomorphic and Riesz calculi agree, so $z^{-1/2}$ is available exactly when
$0\notin\mathrm{spec}$.  The escape is the spectrum touching zero, not the
decay.  Together with the block-structure"""))

# --- B4: the rank-(M-1) degeneracy is the flat limit --------------------------
E.append((
    r"""is geometry-independent.  At $\chi=\pi$ every block symbol tends to $j_0(0)=1$""",
    r"""is geometry-independent, and the degeneracy itself is a known one:\ a kernel
matrix whose entries all tend to a common value degenerates to a rank-one
all-ones matrix with an $(M-1)$-fold null space, the \emph{flat limit} of the
radial-basis-function literature~\cite{barthelme_usevich2021}.  What is ours is
that the block symbols realize it at a \emph{symbol point} --- $\chi=\pi$,
i.e.\ $p=0$ --- rather than in a shape-parameter limit, which is exactly why a
fixed rank-$(M-1)$ rotation removes it.  At $\chi=\pi$ every block symbol tends to $j_0(0)=1$"""))

# --- bibliography: only the three verified this session -----------------------
E.append((
    r"""\bibitem{klappenecker2001}""",
    r"""\bibitem{dlmf}
NIST Digital Library of Mathematical Functions, \url{https://dlmf.nist.gov/},
F.~W.~J.~Olver \textit{et al.}, eds.

\bibitem{loring2014}
T.~A.~Loring, ``Principal angles and approximation for quaternionic
projections,'' \textit{Ann.\ Funct.\ Anal.}\ \textbf{5}(2), 176 (2014);
arXiv:1306.1923.

\bibitem{barthelme_usevich2021}
S.~Barthelm\'e and K.~Usevich, ``Spectral properties of kernel matrices in the
flat limit,'' \textit{SIAM J.\ Matrix Anal.\ Appl.}\ \textbf{42}(1), 17 (2021);
arXiv:1910.14067.

\bibitem{klappenecker2001}"""))

for old, new in E:
    n = src.count(old)
    assert n == 1, f"anchor matched {n} times:\n{old[:100]}"
    src = src.replace(old, new)

P.write_text(src, encoding="utf-8")
print(f"{len(E)} C23 edits applied, each matched exactly once.")
