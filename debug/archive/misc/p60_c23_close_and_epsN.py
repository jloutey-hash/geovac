"""Items 2 and 3 (2026-09-12): the eps_N one-direction result + closing the C23
owed citations with primaries verified THIS session.

Verified here (search or primary read), and only these are cited:
  Jordan 1875 (Bull. SMF 3, 103)            -- principal angles, NUMDAM-confirmed
  Bjorck-Golub 1973 (Math. Comp. 27, 579)   -- their computation, DOI-confirmed
  Eijkhout-Vassilevski 1991 (SIAM Rev. 33, 405) -- the CBS constant
  Hartman-Wintner 1954 (Amer. J. Math. 76, 867) -- Toeplitz spectrum = conv(ess range)
  Lowdin 1950 (J. Chem. Phys. 18, 365)      -- read in Slater-Koster's own footnote 12
  Slater-Koster 1954 (Phys. Rev. 94, 1498)  -- PRIMARY READ, p.1500 Sec.II
  Rokob-Szabados-Surjan                     -- existence + abstract; PDF cert-blocked
  Jaffard 1990 (Ann. IHP C 7, 461)
  Groechenig-Leinert 2006 (TAMS 358)
  Driscoll-Fornberg 2002 (Comput. Math. Appl. 43, 413)

NOT verified, so named in prose with no bibitem: the "Jordan-Wielandt" label
(Stewart-Sun / Horn-Johnson are books, unopened). And the scan's claim that
Slater-Koster's APPENDIX states the symmetry theorem is NOT relied on -- what
was read is their Sec. II use of Lowdin for the multi-centre non-orthogonality
problem, which is what the citation now carries.
"""
from pathlib import Path

P = Path("papers/group2_quantum_chemistry/paper_60_sturmian_secular_quantum.tex")
src = P.read_text(encoding="utf-8")
E = []

# ---- A1 ---------------------------------------------------------------------
E.append((
    r"""the principal angles
are Jordan's, the spectrum $\{1\pm\sigma_k\}$ is the Jordan--Wielandt form of
the off-diagonal block, and the resulting condition number is the two-block
CBS constant of subspace-correction theory.""",
    r"""the principal angles
are Jordan's~\cite{jordan1875} (their computation
is~\cite{bjorck_golub1973}), the spectrum $\{1\pm\sigma_k\}$ is the
Jordan--Wielandt form of the off-diagonal block, and the resulting condition
number is the two-block CBS constant of subspace-correction
theory~\cite{eijkhout_vassilevski1991}."""))

# ---- the eps_N measurement: the degeneracy is ONE DIRECTION -----------------
E.append((
    r"""\textbf{[SYMBOLIC + MEASURED]} One residue is worth stating rather than
hiding, and it has a closed form.""",
    r"""\textbf{[MEASURED]} One reading of this degeneracy should be resisted, and
the measurement excluding it is cheap.  Because the intra-center block is
exactly the identity, the $A$-set is \emph{orthonormal} in the $S$ metric, so
the deficit of a displaced Sturmian against the one-center span is directly
computable:\ $\varepsilon_N^2=1-\sum_{i\le N}\langle\chi^A_i,\chi^B_1\rangle^2$
by Bessel's inequality.  It does \emph{not} tend to zero.  It plateaus at
$0.380$, $0.696$ and $0.907$ for $kR=1,2,4$, flat from $N=16$ to $N=256$
(\texttt{tests/test\_paper60\_one\_direction.py}).  The one-center shared-scale
set is therefore measurably far from complete \emph{in the molecular metric} ---
between $38\%$ and $91\%$ of a displaced basis function lies outside its span
--- so the degeneracy is \emph{not} the statement that one center's span
already contains the other's.  What $\sigma_{\max}\to1$ asserts is weaker and
more useful:\ there exists a \emph{combination} in the $B$-span that the
$A$-span captures, concentrated where the symbol attains its supremum.  The
near-dependence is \emph{one direction}, not a diffuse property of the basis
--- which is why the fixed, geometry-independent rank-$(M-1)$ rotation of
Sec.~\ref{sec:resource} removes it, and why the coupling's participation ratio
stays small.  It also qualifies the motivation offered below for this basis:\
the Coulomb--Sturmians' completeness at a single scale is completeness in the
\emph{atomic} metric, and the molecular problem is posed in $V_0$, where the
one-center set is not complete.

\textbf{[SYMBOLIC + MEASURED]} One residue is worth stating rather than
hiding, and it has a closed form."""))

# ---- A3 ---------------------------------------------------------------------
E.append((
    r"""the symmetry-preservation property
of L\"owdin orthogonalization, known in this paper's own field since
Slater and Koster, and in operator terms the statement that the block-diagonal
matrices form a commutant and are therefore inverse-closed.""",
    r"""the symmetry-preservation property
of L\"owdin's symmetric orthogonalization~\cite{lowdin1950}, applied to exactly
this multi-center non-orthogonality problem by Slater and
Koster~\cite{slater_koster1954} and since analyzed for the conditions under
which the symmetry actually survives~\cite{rokob2008};\ in operator terms it is
the statement that the block-diagonal matrices form a commutant and are
therefore inverse-closed."""))

# ---- A4 ---------------------------------------------------------------------
E.append((
    r"""classical fact that a self-adjoint Toeplitz operator's spectrum is the convex
hull of its symbol's essential range, so that""",
    r"""classical fact that a self-adjoint Toeplitz operator's spectrum is the convex
hull of its symbol's essential range~\cite{hartman_wintner1954}, so that"""))

# ---- B3 ---------------------------------------------------------------------
E.append((
    r"""has polynomial off-diagonal decay of order $>1$ in one dimension, which is
Jaffard's class, and inversion there preserves the decay""",
    r"""has polynomial off-diagonal decay of order $>1$ in one dimension, which is
Jaffard's class~\cite{jaffard1990}, and inversion there preserves the decay"""))

# ---- B4 ---------------------------------------------------------------------
E.append((
    r"""the \emph{flat limit} of the
radial-basis-function literature~\cite{barthelme_usevich2021}.""",
    r"""the \emph{flat limit} of the
radial-basis-function literature~\cite{driscoll_fornberg2002,barthelme_usevich2021}."""))

# ---- bibliography -----------------------------------------------------------
E.append((
    r"""\bibitem{dlmf}""",
    r"""\bibitem{jordan1875}
C.~Jordan, ``Essai sur la g\'eom\'etrie \`a $n$ dimensions,''
\textit{Bull.\ Soc.\ Math.\ France}\ \textbf{3}, 103 (1875).

\bibitem{bjorck_golub1973}
\AA.~Bj\"orck and G.~H.~Golub, ``Numerical methods for computing angles between
linear subspaces,'' \textit{Math.\ Comp.}\ \textbf{27}, 579 (1973).

\bibitem{eijkhout_vassilevski1991}
V.~Eijkhout and P.~Vassilevski, ``The role of the strengthened
Cauchy--Buniakowskii--Schwarz inequality in multilevel methods,''
\textit{SIAM Rev.}\ \textbf{33}, 405 (1991).

\bibitem{hartman_wintner1954}
P.~Hartman and A.~Wintner, ``The spectra of Toeplitz's matrices,''
\textit{Amer.\ J.\ Math.}\ \textbf{76}, 867 (1954).

\bibitem{lowdin1950}
P.-O.~L\"owdin, \textit{J.\ Chem.\ Phys.}\ \textbf{18}, 365 (1950).

\bibitem{slater_koster1954}
J.~C.~Slater and G.~F.~Koster, ``Simplified LCAO method for the periodic
potential problem,'' \textit{Phys.\ Rev.}\ \textbf{94}, 1498 (1954).

\bibitem{rokob2008}
T.~A.~Rokob, \'A.~Szabados, and P.~R.~Surj\'an, ``A note on the symmetry
properties of L\"owdin's orthogonalization schemes,'' in
\textit{Zahradn\'ik Festschrift} (2008).

\bibitem{jaffard1990}
S.~Jaffard, ``Propri\'et\'es des matrices bien localis\'ees pr\`es de leur
diagonale et quelques applications,'' \textit{Ann.\ Inst.\ H.\ Poincar\'e
Anal.\ Non Lin\'eaire}\ \textbf{7}, 461 (1990).

\bibitem{grochenig_leinert2006}
K.~Gr\"ochenig and M.~Leinert, ``Symmetry and inverse-closedness of matrix
algebras and functional calculus for infinite matrices,''
\textit{Trans.\ Amer.\ Math.\ Soc.}\ \textbf{358}, 2695 (2006).

\bibitem{driscoll_fornberg2002}
T.~A.~Driscoll and B.~Fornberg, ``Interpolation in the limit of increasingly
flat radial basis functions,'' \textit{Comput.\ Math.\ Appl.}\ \textbf{43},
413 (2002).

\bibitem{dlmf}"""))

for old, new in E:
    n = src.count(old)
    assert n == 1, f"anchor matched {n} times:\n{old[:110]}"
    src = src.replace(old, new)

P.write_text(src, encoding="utf-8")
print(f"{len(E)} edits applied, each matched exactly once.")
