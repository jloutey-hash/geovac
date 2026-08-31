r"""RECONSTRUCTION step 5: the N--F hub paragraph (degeneracy caveat) and the
last retired-scaling mentions.

Final measured (post-order-fix):
  hub I(2p-1,2p+1), eigensolver-returned member: N 0.148  O 1.063  F 0.311
  sampled maximum over the multiplet:            N 0.69   O 1.34   F 1.32
  ground degeneracies:                           N 4      O 7      F 6
  exact core closure (every support det carries 1s^2) at all three
The ln 8 / ln 16 "attained ceilings" are RETIRED (they were artifacts of a
diagonal-only single-orbital entropy routine that discarded the spin
coherence <up|rho|dn>).
"""
import io

P = "papers/group6_precision_observations/paper_26_entanglement.tex"
s = io.open(P, encoding="utf-8").read()
n = 0


def rep(old, new, label):
    global s, n
    assert s.count(old) == 1, f"{label}: {s.count(old)} matches"
    s = s.replace(old, new)
    n += 1
    print(f"  ok  {label}")


rep(r"""\medskip\noindent\textbf{Nitrogen--Fluorine ($\mathbf{Z = 7}$--$\mathbf{9}$, 7--9 electrons):}
Hub transitions to $2p(\pm 1)$.  P-shell correlation dominates
the mutual information network.""",
    r"""\medskip\noindent\textbf{Nitrogen--Fluorine ($\mathbf{Z = 7}$--$\mathbf{9}$, 7--9 electrons):}
The hub transitions to the $2p(\pm 1)$ pair, which carries the
mutual information:\ $I(2p_{-1}, 2p_{+1}) = 0.148$ (N), $1.063$ (O),
$0.311$ (F) \emph{for the multiplet member the eigensolver returns}.
These ground states are $4$-, $7$- and $6$-fold degenerate, and the
value is \emph{member-dependent}:\ it ranges from $0$ on
single-determinant members up to a sampled maximum of $0.69$ (N),
$1.34$ (O), $1.32$ (F).  The quoted values are therefore
eigensolver-dependent, not basis constants, and even the qualitative
reading ``the valence network is alive'' is member-dependent.

\emph{Retracted 2026-08-29.}  Earlier versions reported
$1.837/1.785/1.644$ here and an information-theoretic ceiling of
$\ln 8 \approx 2.08$ (N/F) and $\ln 16 \approx 2.77$ (O), ``attained
exactly by explicit members''.  Both the values and the ceilings were
artifacts of a single-orbital entropy routine that built $\rho_i$
diagonal-only, discarding the spin coherence
$\langle{\uparrow}\rvert\rho_i\lvert{\downarrow}\rangle$ --- nonzero
for precisely these sector-mixed multiplet members.  Because a diagonal
Shannon entropy majorizes the true von Neumann entropy, the routine
over-reported $s_i$ (e.g.\ $0.798$ against a true $0.313$) and the
derived ``mutual information'' violated subadditivity on half the
orbital pairs.  No analytic ceiling is asserted in its place; the
maxima above are MEASURED.

What \emph{is} degeneracy-robust is the core-valence statement, and it
holds as a theorem for all three:\ every determinant in each ground
multiplet carries $1s^2$, so $\rho_{1s}$ is pure for any member and
$I_{\text{cv}} = 0$ exactly.""",
    "N--F hub paragraph + ceiling retraction")

rep(r"""  \frac{S_{\text{bond}}}{S_{\text{core}}} \approx 50,""",
    r"""  \frac{S_{\text{bond}}}{S_{\text{core}}} \approx 40,""",
    "core/bond ratio equation")

rep(r"""Hamiltonians with $O(Q^{2.5})$ Pauli terms for composed molecular""",
    r"""Hamiltonians whose Pauli count is exactly linear in the qubit count at
fixed basis, for composed molecular""",
    "intro Pauli scaling")

rep(r"""The $O(Q^{2.5})$ Pauli scaling of the GeoVac composed qubit""",
    r"""The Gaunt-driven Pauli sparsity of the GeoVac composed qubit""",
    "SS III.D Pauli scaling")

rep(r"""GeoVac Hamiltonians achieve $O(Q^{2.5})$ Pauli scaling with""",
    r"""GeoVac Hamiltonians achieve Gaunt-driven Pauli sparsity with""",
    "discussion Pauli scaling")

io.open(P, "w", encoding="utf-8").write(s)
print(f"\n{n} edits applied")
