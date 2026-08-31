r"""Paper 26: the I_cv = 0 'theorem' premise fails at Z=7 under the exact
ERI rule (16 support determinants lack 1s^2, carrying 0.0014% of the
multiplet weight), and the stale tier line still advertises the retracted
ceiling attainment."""
import io

P = "papers/group6_precision_observations/paper_26_entanglement.tex"
s = io.open(P, encoding="utf-8").read()

old = r"""What \emph{is} degeneracy-robust is the core-valence statement, and it
is robust as a \emph{theorem} rather than a sampled observation: every
determinant in each ground multiplet carries $1s^2$, so $\rho_{1s}$ is
pure for \emph{any} member, hence $I_{\text{cv}} = 0$ exactly --- for
every member, without exception.  (An earlier version of this caveat
claimed the valence network exceeds $0.3$ in every member; that was an
artifact of sampling sector-mixed eigenvectors, which cannot reach the
single-determinant members.)  \emph{Tier:} the multiplet structure and
the $I_{\text{cv}} = 0$ statement are derivations from a
machine-verified premise (every support determinant carries $1s^2$,
pinned sampling-free by the backing test); the ceiling attainment is
exact and constructive (explicit members, backed to $10^{-9}$); the
degeneracies are MEASURED; the member-representative values are
eigensolver-dependent by construction."""

new = r"""What is degeneracy-robust is the core-valence statement, though it is
now a \emph{bound} rather than the exact theorem claimed previously.
For O and F every determinant in the ground multiplet carries $1s^2$, so
$\rho_{1s}$ is pure for any member and $I_{\text{cv}} = 0$ exactly.  For
N this premise fails under the exact Coulomb selection rule:\ sixteen
support determinants lack $1s^2$, because the restored $m$-changing
multipoles include the core$\,\to\,$valence double excitation
$\langle 1s\,1s | 2p_{-1}\,2p_{+1}\rangle$.  They carry
$5.8\times10^{-5}$ of the multiplet weight ($0.0014\%$), giving a
measured $I_{\text{cv}} \le 4.3\times10^{-4}$ over sampled members ---
so core closure at N is \emph{approximate at the $10^{-4}$ level}, not
exact.  (An earlier version of this caveat claimed the valence network
exceeds $0.3$ in every member; that was an artifact of sampling
sector-mixed eigenvectors, which cannot reach the single-determinant
members.)  \emph{Tier:} the multiplet structure is MEASURED; the
$I_{\text{cv}} = 0$ statement is a derivation from a machine-verified
premise \emph{for O and F only}, and a MEASURED bound at N; the
degeneracies are MEASURED; the member-representative values are
eigensolver-dependent by construction."""

assert s.count(old) == 1, s.count(old)
io.open(P, "w", encoding="utf-8").write(s.replace(old, new))
print("  ok  I_cv theorem scoped to O/F + N bound; stale tier line fixed")
