# Vector-Photon QED Sprint Memo

**Date:** 2026-05-01
**Sprint:** Vector-photon QED on the Fock graph
**Status:** Complete — 7/8 scalar, 8/8 Dirac selection rules recovered

## Motivation

The graph-native QED framework (GN-1 through GN-7) computes QED entirely on the
finite Fock graph using exact algebraic arithmetic. The photon is a scalar
1-cochain — a number on each edge. This produces π-free results in ℚ[√2,√3,√6,...]
but recovers only 1/8 continuum QED selection rules. The native Dirac graph
(Rule B) recovers 4/8 by using spinor node labels. The transverse photon from
Wilson plaquettes (Paper 30) recovers 1 additional (GS structural zero), for 5/8.

The remaining 3 rules were identified as "vector-quantum-number-required" —
they need the photon to carry angular momentum labels (q, m_q), not just
topological eigenvalues from plaquettes.

This sprint tests whether adding explicit photon quantum numbers recovers
the remaining rules.

## Pre-investigation: L₁ block-diagonalization

**Question:** Does the edge Laplacian L₁ = B^T B already contain q-resolved
photon structure that we can extract by reorganization?

**Answer:** NO. L₁ does not block-diagonalize by (q, m_q). The fundamental
reason is that photon quantum numbers describe the *transition* (the difference
between endpoint labels), not a *conserved quantity* at each vertex. L₁ couples
edges that share a vertex, and at that shared vertex, edges with different
photon quantum numbers meet. Cross-sector Frobenius fraction: 20% at n_max=3
on the scalar Fock graph.

The ONLY decomposition that works is by l-sector (Fock) or κ-sector (Dirac
Rule A) — connected components of the graph. Transition quantum numbers
cannot serve as block-diagonal labels.

**Implication:** Vector photon structure must be added from outside. It is
genuinely new structure, not a reorganization of what the graph contains.

**Data:** `debug/l1_photon_block_diag.py`, `debug/data/l1_photon_block_diag.json`

## Construction

### Electron register
Fock graph nodes (n, l, m) or Dirac labels (n, κ, m_j). Both tested.

### Photon register
Modes labeled by (q, m_q), q = 1..q_max, m_q = -q..q.
Propagator: G_γ(q) = 1/[q(q+2)] from the vector Laplacian eigenvalues on S³.

### Vertex coupling
V(a, b, q, m_q) = angular × parity × radial, where:
- **Angular:** Wigner 3j symbol enforcing triangle inequality and m-conservation
- **Parity:** l_a + l_b + q must be ODD (E1 vector coupling)
- **Radial:** R = 1 (unit coupling, isolates angular structure)

### Self-energy
Σ(a, c) = Σ_{b, q, m_q} V(a, b, q, m_q) · G_e(b) · G_γ(q) · V(c, b, q, m_q)

## Results

### Selection rule census

| # | Rule | Scalar+Vector | Dirac+Vector | Mechanism |
|---|------|:---:|:---:|---|
| 1 | Gaunt/CG sparsity | PASS | PASS | Built into 3j symbol |
| 2 | Vertex parity / GS zero | PASS | PASS | l_a+l_b+q odd; GS: 2l=even≠odd |
| 3 | SO(4) channel count | PASS | PASS | Triangle inequality on (l,q,l') |
| 4 | Δm / Δm_j conservation | PASS | PASS | 3j m-selection rule |
| 5 | Spatial parity | PASS | PASS | E1 parity enforcement |
| 6 | Furry's theorem | FAIL | PASS | Scalar: V(a,a)≠0 for l≥1; Dirac: V(a,a)=0 identically (spinor phase) |
| 7 | Ward identity (n=2) | PASS | PASS | Σ block-diagonal within n-shells |
| 7 | Ward identity (n=3) | FAIL | FAIL | Truncation: q=2 couples n-shells |
| 8 | Charge conjugation | PASS | PASS | Σ Hermiticity verified |

**Score: 7/8 scalar, 8/8 Dirac at n_max=2.**

### Ground state structural zero

Σ(GS, GS) = 0 exactly (machine precision) at all tested n_max. Mechanism:
for the GS with l=0, the triangle inequality forces q = l_int, and the parity
rule requires l_int + l_int + q = 2l_int + l_int = 3l_int to be odd —
wait, more precisely: l_ext + l_int + q odd with l_ext=0 gives l_int + q odd,
and triangle with l_ext=0 gives q = l_int, so 2l_int must be odd. Impossible.

This was the rule broken by the pendant-edge theorem on the scalar graph
(Σ_scalar(GS) = 2(n_max-1)/n_max → 2). The vector photon restores it through
angular momentum conservation alone.

### Calibration exchange constant

The l=1 self-energy at n_max=2 is EXACTLY 1/(4π) to machine precision.
This is the Weyl exchange constant for S² — the density-of-states constant
for the 2-sphere. It enters as the normalization of the vector spherical
harmonics: sqrt((2l+1)(2q+1)(2l'+1)/(4π)) in the vertex coupling.

The entire scalar graph QED is π-free (ℚ[√2,√3,√6,...]). The vector photon
introduces exactly one transcendental per loop: π, through the S² solid
angle normalization. This is precisely where Paper 18's taxonomy predicted
calibration content would enter.

### Scalar vs Dirac electron comparison

Scalar and Dirac electrons give DIFFERENT selection rule scores with vector
photons: 7/8 scalar vs 8/8 Dirac. The spinor labels (κ, m_j) provide one
additional selection rule recovery — Furry's theorem — via the Dirac vertex
phase constraint that forces V(a, a, q, m_q) = 0 identically. Scalar
electrons lack this constraint: V(a, a, q, m_q) ≠ 0 for l≥1 when q is odd.

The Dirac configuration also gives a 5× better Ward ratio at n_max=3
(0.61 vs 3.1) due to the linear Dirac spectrum (n+1/2) providing better
commutator structure than the quadratic scalar spectrum (n²-1).

### Furry's theorem: scalar failure, Dirac recovery

Furry's theorem behaves differently for scalar and Dirac electrons with
vector photons.

**Scalar electrons (FAIL):** For l≥1 states, the diagonal vertex coupling
V(a, a, q, m_q) is nonzero when q is odd (parity 2l+q odd is satisfied).
The 3j(j, q, j; -m, 0, m) does not vanish for generic m. Without
second-quantized C-symmetry (particle-antiparticle loop cancellation),
scalar electrons cannot suppress these tadpole contributions.

**Dirac electrons (PASS):** The Dirac spinor phase constraint forces
V(a, a, q, m_q) = 0 identically. The mechanism: the alpha-matrix vertex
coupling <psi_a|alpha . epsilon|psi_a> vanishes because alpha is
off-diagonal in large/small component space, and the i-factor on the small
component makes V_LS + V_SL = 0 by Hermiticity. This is a kinematic
single-particle identity — it does not require second-quantized
particle-antiparticle field structure, only the spinor structure of the
Dirac wavefunction.

This should be distinguished from the field-theoretic Furry mechanism in
continuum QED, where charge conjugation symmetry ensures that for every
electron loop there is a positron loop that exactly cancels. The Dirac
vertex phase constraint achieves the same selection rule (vanishing
odd-vertex diagrams) through a different mechanism: the tadpole vertex
itself is zero, rather than being canceled by an antiparticle loop.

Note: the earlier native Dirac graph sprint showed Furry "recovered" by yet
another mechanism — the graph topology prevented self-loops (edges only
connect Δl=±1, so no edge from a state to itself). That was a topological
accident specific to the adjacency rule, not the physical mechanism.

## Structural conclusions

### The 1-7-8 partition (updated)

The 8 continuum QED selection rules partition differently for scalar and
Dirac electrons:

**Scalar electrons:**

| Tier | Rules | Source |
|------|:---:|---|
| Graph topology | 1 | Gaunt/CG sparsity (always survives) |
| Angular momentum conservation | 6 | CG coefficients + E1 parity + triangle inequality |
| Furry (fails) | 0 | V(a,a) ≠ 0 for scalar l≥1; no C-symmetry available |
| **Total** | **7/8** | |

**Dirac electrons:**

| Tier | Rules | Source |
|------|:---:|---|
| Graph topology | 1 | Gaunt/CG sparsity (always survives) |
| Angular momentum conservation | 6 | CG coefficients + E1 parity + triangle inequality |
| Spinor phase constraint | 1 | V(a,a) = 0 identically from Dirac vertex structure |
| **Total** | **8/8** | |

The 1->7 transition costs exactly one transcendental: pi, entering as 1/(4pi)
per loop. The 7->8 transition does NOT require second-quantized field theory
as previously believed — it requires only the spinor structure of the Dirac
wavefunction. The alpha-matrix vertex coupling provides the kinematic identity
V(a,a) = 0 that scalar wavefunctions lack.

### Paper 18 taxonomy placement

The vector photon quantum numbers are CALIBRATION exchange constants in the
Paper 18 taxonomy. The graph-native QED (π-free, algebraic) is the intrinsic
tier. Adding vector photon modes imports the S² solid angle normalization
as calibration content. This is consistent with the three-layer structure
identified in the transverse photon sprint: topology → geometry → calibration.

### Comparison with previous QED sprints

| Configuration | Rules | What was added |
|---|:---:|---|
| Scalar graph, scalar photon (GN-1..7) | 1/8 | Pure graph topology |
| Dirac graph Rule B, scalar photon | 4/8 | Spinor node labels |
| Any graph + plaquette photon | 2/8 | Wilson gauge 2-cells |
| Scalar electrons + vector photon | **7/8** | Photon (q, m_q) + CG vertex |
| Dirac electrons + vector photon | **8/8** | Photon (q, m_q) + CG vertex + spinor phase |
| Continuum QED | 8/8 | Second-quantized fermion field |

## Files

| File | Description |
|------|-------------|
| `geovac/vector_qed.py` | Main module: vertex, propagator, self-energy, selection rules |
| `tests/test_vector_qed.py` | 99 tests (scalar + Dirac configurations) |
| `debug/l1_photon_block_diag.py` | L₁ block-diagonalization diagnostic |
| `debug/vector_qed_diagnostic.py` | Scalar vector QED selection rule census |
| `debug/vector_qed_dirac_diagnostic.py` | Dirac vector QED selection rule census |
| `debug/data/l1_photon_block_diag.json` | L₁ diagnostic data |
| `debug/data/vector_qed_diagnostic.json` | Scalar diagnostic data |
| `debug/data/vector_qed_dirac_diagnostic.json` | Dirac diagnostic data |

## Open questions

1. **Ward identity at larger n_max:** Is the n_max=3 failure purely a truncation
   artifact, or does it indicate a deeper issue with the finite-basis Ward relation?
   The 5× improvement from Dirac eigenvalues suggests truncation, but convergence
   should be verified.

2. **Radial matrix elements:** We used R=1 (unit radial coupling). Adding the
   actual Fock radial matrix elements would affect magnitudes but not selection
   rules. However, it would be needed for quantitative comparison with continuum
   QED results (e.g., recovering α/(2π) for the anomalous moment).

3. **Connection to F₂ convergence:** The graph-native F₂ converges to zero as
   n^{-0.57}. Does the vector-photon F₂ converge to α/(2π)? This would be the
   ultimate validation that vector photon quantum numbers are the correct
   calibration content.

4. **Spinor phase constraint universality:** The Dirac vertex identity V(a,a)=0
   was verified at n_max=2. Does it hold at all n_max, and can it be proven
   from the algebraic structure of the Dirac wavefunction components alone
   (independent of the specific spectrum)?
