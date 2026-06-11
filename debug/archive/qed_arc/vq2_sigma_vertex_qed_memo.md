# VQ-2: Pauli sigma vertex for graph-native QED on S^3

**Date:** 2026-05-01
**Status:** NEGATIVE -- sigma vertex is structurally wrong for QED
**Baseline comparison:** 1/8 (scalar Fock), 4/8 (Dirac graph Rule B)
**This investigation:** 2/8 selection rules survive

---

## Investigation summary

VQ-2 tests whether the Pauli spin matrices sigma_x, sigma_y, sigma_z can serve
as a VECTOR vertex coupling for graph-native QED on the finite Fock graph,
potentially recovering some of the 4 continuum QED selection rules that remain
broken on both the scalar Fock graph (1/8) and the Dirac graph (4/8).

The investigation is a clean negative: the sigma vertex is structurally
incapable of serving as a QED vertex because it preserves both n and l
(only changing kappa/j within the same shell), making the coupling graph
disconnected and unable to mediate inter-shell energy exchange.

---

## Step 1: Sigma matrices in the coupled (j, m_j) basis

The Pauli spin matrices sigma_a act on the spin-1/2 part of the Dirac spinor.
In the coupled (j, m_j) basis used by DiracLabel, they are constructed via CG
decomposition:

    sigma_z[a,b] = sum_{m_s} CG(l,m_l,1/2,m_s|j_a,m_j_a) * (2*m_s) * CG(l,m_l,1/2,m_s|j_b,m_j_b)
    sigma_+ = CG product for m_s: -1/2 -> +1/2
    sigma_- = CG product for m_s: +1/2 -> -1/2
    sigma_x = sigma_+ + sigma_-    (REAL matrix)
    sigma_y = -i(sigma_+ - sigma_-)  (PURELY IMAGINARY matrix)

**Pauli algebra verified within each (n,l) block:**
- sigma^2 = sigma_x^2 + sigma_y^2 + sigma_z^2 = 3*I (PASS, all blocks)
- [sigma_x, sigma_y] = 2i*sigma_z (PASS)
- [sigma_y, sigma_z] = 2i*sigma_x (PASS)
- [sigma_z, sigma_x] = 2i*sigma_y (PASS)
- sigma_x, sigma_z: real symmetric (Hermitian)
- sigma_y: purely imaginary anti-symmetric (Hermitian: (A*)^T = A)
- All blocks traceless

**Normalization bug found and fixed:** The initial sigma_+/sigma_- construction
divided by 2 (following the convention sigma_+ = (sigma_x + i*sigma_y)/2), but
the CG product computes the STANDARD raising/lowering operator (without the 1/2),
so sigma_x = sigma_+ + sigma_- without division. The corrected normalization
gives sigma_a^2 = I within each block, as required.

---

## Step 2: Sigma-edge graph topology

The sigma-edge graph has edges wherever any sigma_a has a nonzero off-diagonal
matrix element (i.e., a != b and |sigma_mu[a,b]| > 0 for some mu).

**Key structural property: sigma preserves both n AND l.** It only changes
kappa (equivalently j) within the same (n,l) block. This means:

- No inter-shell edges (unlike scalar Fock graph T+/T- ladders)
- No inter-l edges (unlike Dirac Rule B E1 dipole)
- Graph is **disconnected**: beta_0 = 3 components
  - Component 1: n=1, l=0 (2 nodes, 1 edge) -- the GS pair
  - Component 2: n=2, l=0 (2 nodes, 1 edge)
  - Component 3: n=2, l=1 (6 nodes, 10 edges)

**Graph statistics:**
- V = 10, E = 12
- beta_1 = 5 (independent cycles)
- GS pendant: YES (degree 1 for both GS nodes)
- Degree sequence: [1, 1, 1, 1, 2, 2, 4, 4, 4, 4]

**sigma_z contributes NO edges:** sigma_z is diagonal in the (j, m_j) basis,
so sigma_z[a,b] = 0 for a != b. Only sigma_x and sigma_y create off-diagonal
couplings. sigma_z acts as an on-site (diagonal) vertex weighting only.

---

## Step 3: QED self-energy

Two contraction conventions tested:

**Hermitian contraction** (V . G_gamma . V^dagger):
- Sigma(GS) = I (nonzero) -- structural zero NOT recovered
- Tr(Sigma) = 6.646
- Eigenvalues: {0.224(x2), 0.457(x2), 0.643(x2), 1.000(x4)}
- Positive semidefinite: YES
- Zero eigenvalues: 0

**Transpose contraction** (V . G_gamma . V^T):
- Sigma(GS) = 0 (exactly zero)
- Tr(Sigma) = 0.284

**Mechanism for transpose GS zero:** In the GS 2-node pair, the only edge
carries sigma_x (real symmetric) and sigma_y (imaginary antisymmetric). In the
transpose contraction:

    sigma_x @ sigma_x^T = sigma_x^2 = I
    sigma_y @ sigma_y^T = sigma_y @ (-sigma_y) = -sigma_y^2 = -I

These cancel exactly: I + (-I) = 0. The Hermitian contraction uses sigma_y^dagger
= sigma_y (not sigma_y^T = -sigma_y), so both contribute positively: I + I = 2I.

**This is NOT a physical selection rule.** The transpose zero is a kinematic
artifact of the sigma_y antisymmetry, not a vertex-parity or SO(4) channel count
protection. The physically correct contraction for Hermitian operators is V.G.V^dagger,
which gives Sigma(GS) = I != 0. The pendant-edge mechanism (GS is degree 1)
still operates, as in the scalar graph QED.

---

## Step 4: Selection rule census

| Rule | Status | Scalar (1/8) | Dirac B (4/8) | Sigma (2/8) |
|------|--------|--------------|---------------|-------------|
| Delta(m_j) conservation | SURVIVES | N/A | SURVIVES | SURVIVES |
| Spatial parity (E1) | BROKEN | BROKEN | SURVIVES | BROKEN |
| Gaunt/CG sparsity | SURVIVES | SURVIVES | SURVIVES | SURVIVES |
| Vertex parity | BROKEN | BROKEN | BROKEN | BROKEN |
| SO(4) channel count | BROKEN | BROKEN | BROKEN | BROKEN |
| Charge conjugation | BROKEN | BROKEN | BROKEN | BROKEN |
| Furry's theorem | BROKEN | BROKEN | SURVIVES | BROKEN |
| Ward identity | BROKEN | BROKEN | BROKEN | BROKEN |

**Result: 2/8 rules survive** -- WORSE than the Dirac graph baseline (4/8).

The sigma vertex actually LOSES rules compared to the Dirac graph:
- Spatial parity: sigma preserves l (all couplings have l_a = l_b), so
  l_a + l_b is always even. This is NOT E1 parity (which requires odd).
  The rule is trivially satisfied in a degenerate sense, but provides no
  selection power.
- Furry's theorem: the sigma vertex creates nonzero tadpoles (Tr per
  polarization != 0 in the sigma-edge subgraph).

---

## Step 5: Algebraic characterization

- Self-energy matrix is real (all imaginary parts < 1e-10)
- sigma_x, sigma_z matrix elements are CG products = rational (exact)
- sigma_y matrix elements are i * (CG products) = i * rational (exact)
- Bilinear products sigma_mu[a,c] * sigma_mu[d,b] cancel the i factors,
  giving rational * rational or (i*rational)*(i*rational) = -rational,
  net result: all rational
- **Number field: Q (rationals)** -- simpler than scalar Fock graph
  Q[sqrt(2), sqrt(3), sqrt(6)] or Dirac graph Q[sqrt(2), sqrt(17), ...]
- pi-free: YES (trivially, since all quantities are rational)

---

## Why the sigma vertex fails structurally

The 4 continuum QED selection rules that are broken on ALL finite graph
constructions (vertex parity, SO(4) channel count, Ward identity, charge
conjugation) share a common requirement: the INTERACTION vertex must couple
states in DIFFERENT energy shells via a photon that carries energy-momentum.

The sigma vertex is purely **intra-shell**: sigma_a preserves both n and l,
only changing the spin quantum number (kappa/j) within the same (n,l) block.
This means:

1. **No vertex parity:** sigma preserves n, so n_1 + n_2 is always even.
   The parity rule requires n_1 + n_2 + q to be odd, which needs the photon
   to change the electron's shell.

2. **No SO(4) channel count:** The SO(4) channel function W(n_ext, n_int, q)
   requires coupling between different n-shells. Intra-shell coupling has
   no channel structure to constrain.

3. **No Ward identity:** [D, Sigma] != 0 because the Dirac operator D is
   diagonal in n (proportional to eigenvalues n + 3/2) while Sigma from
   the intra-shell sigma vertex commutes with n but not with the kappa-dependent
   part of D.

4. **No charge conjugation:** The sigma vertex lacks the C-symmetry of the
   vector photon coupling.

The sigma vertex provides VECTOR quantum numbers (3 polarizations, proper
Pauli algebra) but on the WRONG coupling topology (intra-shell). What
graph-native QED needs is the combination of:
- Inter-shell coupling (provided by the existing graph edges)
- Vector quantum numbers (provided by sigma_a)
- Both features simultaneously (neither alone is sufficient)

This suggests the correct next step would be to tensor the sigma vertex
structure with the existing graph adjacency, creating a vertex of the form
A_{ab} * sigma_mu (where A is the adjacency matrix). This is essentially
the sigma-dot-nabla coupling of Dirac QED, but on the discrete graph.

---

## Comparison with VQ-1 (sigma.L)

VQ-1 tested sigma.L (spin-orbit coupling operator) as vertex. VQ-1 found
sigma.L is ALSO intra-shell AND zero on l=0, making it even more restrictive
than bare sigma_a. The VQ-2 sigma vertex at least couples l=0 states (via
sigma_x and sigma_y, which change m_j within the same kappa block), but
still fails the inter-shell requirement.

The pattern from VQ-1 and VQ-2: any operator that acts only on INTERNAL
quantum numbers (spin, angular momentum) without coupling to the graph
adjacency structure cannot serve as a QED vertex. The photon must be
associated with graph EDGES (which encode energy-shell transitions), not
with node-local operators.

---

## Data files

- Computation script: `debug/vq2_sigma_vertex_qed.py`
- Numerical data: `debug/data/vq2_sigma_vertex_qed.json`
- This memo: `debug/vq2_sigma_vertex_qed_memo.md`
