r"""RECONSTRUCTION step 3: Paper 26 abstract + provenance, rebuilt from the
preserved text and today's final (post-order-fix) measurements."""
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


rep(r"""off-diagonal electron repulsion contributes $100\%$ of the
entanglement.""",
    r"""off-diagonal electron repulsion contributes effectively $100\%$ of the
entanglement.""",
    "abstract: effectively-100% hedge")

rep(r"""Any unitary rotation of orbital indices, no matter how small, fills
the ERI tensor from $42\%$ to $99\%$ density in a step-function
transition at the identity.""",
    r"""Any generic (angular-mixing) unitary rotation of orbital indices, no
matter how small, fills the ERI tensor from $17\%$ to $100\%$ density
in a rapid but graded transition at the identity.""",
    "abstract: sparsity transition")

rep(r"""maps reveal a hub migration pattern---$1s \to 2s \to 2p$---that
tracks the partially filled shell across the first row of the
periodic table.""",
    r"""maps reveal a hub migration pattern---$1s \to 2s \to 2p$---that
tracks the partially filled shell across the first row of the
periodic table (for $Z = 7$--$9$ both the values and the presence of a
nonzero valence network are multiplet-member-dependent).""",
    "abstract: hub degeneracy caveat")

rep(r"""Fourth, core-valence decoupling: the core-valence mutual
information falls to $\approx 2\times10^{-3}$ by $Z = 4$ and below
$10^{-3}$ for $Z \geq 5$, quantitatively
justifying the composed fiber bundle factorization used in the
GeoVac molecular framework.
In the composed molecular architecture, per-block entanglement is
$R$-independent, with core/bond entropy ratios of $50\times$ at
matched qubit count.""",
    r"""Fourth, core-valence decoupling: the core-valence mutual
information falls to $\approx 4\times10^{-3}$ by $Z = 4$, decays
monotonically through $Z = 6$, and reaches the
$\approx 10^{-15}$ floating-point noise floor for $Z \geq 7$ --- an
exact core-closure property of this basis (the $1s$ occupation is
exactly $2$), so it constrains the basis rather than justifying the
composed fiber bundle factorization used in the GeoVac molecular
framework.
In the composed molecular architecture, per-block entanglement is
$R$-independent, with core/bond entropy ratios of $40\times$ at
matched qubit count.""",
    "abstract: core-valence")

rep(r"""and that the $O(Q^{2.5})$ Pauli scaling of the GeoVac qubit
Hamiltonians is a direct consequence of this basis choice.""",
    r"""and that the Gaunt sparsity of the GeoVac qubit Hamiltonians ---
whose Pauli count is exactly linear in the qubit count at fixed basis
--- is a direct consequence of this basis choice.""",
    "abstract: Pauli scaling")

io.open(P, "w", encoding="utf-8").write(s)
print(f"\n{n} abstract edits applied")
