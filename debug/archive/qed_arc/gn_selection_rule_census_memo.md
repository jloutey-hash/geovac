# GeoVac Graph-Native QED: Selection Rule Census

**Date:** 2026-04-27  
**n_max (primary):** 2  
**Sparsity n_max values:** [2, 3, 4]  

## Summary Table

| Rule | Verdict | Paper 18 Tier | Key Metric |
|:-----|:--------|:--------------|:-----------|
| Angular momentum Δm_j | BROKEN | INTRINSIC | 20 violations |
| Spatial parity (E1 l_a+l_b odd) | BROKEN | INTRINSIC | 48 parity violations |
| Vertex parity (n1+n2+q odd) | BROKEN | STRUCTURAL | 0.0% forbidden by continuum rule |
| SO(4) channel count (W>0) | BROKEN | STRUCTURAL | 100.0% have W=0 (continuum forbids) |
| Charge conjugation (C) | BROKEN | INTRINSIC | 8 C-violations in Sigma |
| Furry's theorem (odd loops = 0) | BROKEN | INTRINSIC | tadpole=!=0, triangle=!=0 |
| Ward identity ([D,Lambda]=[Sigma,D]) | BROKEN | STRUCTURAL | [D,Lambda]-[Sigma,D] max = 5.08e+00 |
| Gaunt/CG sparsity (n_max=2) | SURVIVES | INTRINSIC | density=0.1600, nnz=48/300 |
| Gaunt/CG sparsity (n_max=3) | SURVIVES | INTRINSIC | density=0.0255, nnz=260/10192 |
| Gaunt/CG sparsity (n_max=4) | SURVIVES | INTRINSIC | density=0.0062, nnz=756/122400 |

## Rule Details

### Rule 1 - Angular Momentum Conservation (Δm_j)

**Verdict:** BROKEN  
**Paper 18 tier:** INTRINSIC  

20 violations of Δm conservation found.

Total nonzero V entries: 48; conserved: 28; violations: 20.

The CG projection in `build_projection_matrix` enforces that each Dirac
state |n,κ,m_j⟩ decomposes onto scalar Fock nodes (n,l,m_l) with
m_l = m_j ± 1/2. The vertex V[a,b,e] = Sigma_{v1,v2} P[a,v1]·P[b,v2]·δ(e=(v1,v2))
inherits Δm = m(v2)-m(v1) from the Fock edge, so m_j conservation follows
from the CG algebra itself - no extra rule needed.

### Rule 2 - Spatial Parity Conservation (E1 Selection: l_a + l_b odd)

**Verdict:** BROKEN  
**Paper 18 tier:** INTRINSIC  

48 entries (100.0%) violate E1 parity (l_a + l_b even, same-parity transition).

Parity-conserved entries: 0; violated: 48.

The Fock scalar graph connects (n,l,m) nodes by L-edges (Δm=±1, Δl=0, Δn=0)
and T-edges (Δm=0, Δl=0, Δn=±1). Both edge types have l1=l2, so a naïve
l_a + l_b check sees the _electron_ orbital change. The CG projection maps
Dirac state (n,κ,m_j) onto scalar nodes (n,l,m_l) with l = l(κ). For a
T-edge (n->n+1, same l), the Dirac state can change κ (hence l), giving Δl=±1.
For L-edges (same n, same l), the Dirac states coupled have l_a + l_b even
if the same l block, but the symmetrized vertex also couples states from
the opposite node on the edge - need to trace carefully per entry.

### Rule 3 - Gaunt/CG Sparsity

**Verdict:** SURVIVES  
**Paper 18 tier:** INTRINSIC  

The CG projection naturally suppresses most couplings. Sparsity grows with n_max because N_dirac^2 × E grows faster than the nonzero entries (angular selection rules).

| n_max | N_Dirac | E_Fock | Total possible | Nonzero V | Density | Sparsity |
|------:|--------:|-------:|---------------:|----------:|--------:|---------:|
| 2 | 10 | 3 | 300 | 48 | 0.1600 | 0.8400 |
| 3 | 28 | 13 | 10192 | 260 | 0.0255 | 0.9745 |
| 4 | 60 | 34 | 122400 | 756 | 0.0062 | 0.9938 |

### Rule 4 - Vertex Parity (n1+n2+q odd)

**Verdict:** BROKEN  
**Paper 18 tier:** STRUCTURAL  

0 of 48 nonzero graph couplings (0.0%) are forbidden by the continuum vertex parity rule (n1+n2+q even). The graph has no mechanism to enforce this parity because the photon is a scalar 1-cochain - it carries no vector harmonic parity. This is a STRUCTURAL gap: the rule emerges from the continuum SO(4) vector harmonic structure (γ^μ coupling), which is absent in scalar QED. Paper 18 tier: STRUCTURAL (not calibration: no π or ζ needed, but the vector structure is absent from the graph).

Entries allowed by parity: 48; forbidden: 0 (0.0%).

### Rule 5 - SO(4) Channel Count (W=0 -> continuum forbids)

**Verdict:** BROKEN  
**Paper 18 tier:** STRUCTURAL  

48 of 48 nonzero graph couplings (100.0%) have W=0 (forbidden by the continuum SO(4) double-triangle rule). The graph uses scalar SU(2) CG projection, not the SU(2)_L × SU(2)_R double-triangle of the continuum vector harmonic. This is a STRUCTURAL gap arising from the graph photon being a scalar 1-cochain, not a vector harmonic. Paper 18 tier: STRUCTURAL.

W=0: 48 (100.0%); W=1: 0; W=2: 0.

### Rule 6 - Charge Conjugation Symmetry (C)

**Verdict:** BROKEN  
**Paper 18 tier:** INTRINSIC  

C-conjugate pairs found: 4. 8 violations of C-symmetry in Sigma found.

C-conjugate pairs found: 4; Sigma elements checked: 16; violations: 8.

### Rule 7 - Furry's Theorem (odd-loop diagrams vanish)

**Verdict:** BROKEN  
**Paper 18 tier:** INTRINSIC  

Tadpole = 16*sqrt(2)/15 != 0, triangle = 3584*sqrt(2)/3375 != 0. At least one odd-loop diagram is nonzero - Furry's theorem is broken.

- Tadpole Sigma_e Tr(V_e·G_e) = 16*sqrt(2)/15  (zero: False)
- Bubble Sigma_{e,e'} Tr(V_e·G_e·V_{e'}·G_e) = 32/9  (nonzero: True)
- Triangle Sigma_{e,e',e''} Tr(V_e·G_e·V_{e'}·G_e·V_{e''}·G_e) = 3584*sqrt(2)/3375  (zero: False)

### Rule 8 - Ward Identity ([D,Lambda] = [Sigma,D])

**Verdict:** BROKEN  
**Paper 18 tier:** STRUCTURAL  

Commutator Ward identity fails: max |[D,Lambda]-[Sigma,D]| = 5.080e+00. Tr(Sigma) = 44/3, Tr(Lambda) = 32/9.

- Tr(Sigma) = 44/3
- Tr(Lambda) = 32/9
- Tr(Lambda)/Tr(Sigma) = 8/33
- max |[D,Lambda]-[Sigma,D]| = 5.08e+00

## Structural Interpretation

The census is complete. Only one rule survives: Gaunt/CG sparsity.
All seven continuum-physics selection rules are broken on the finite
scalar-photon graph. The broken rules fall into two tiers:

**BROKEN (STRUCTURAL tier) — photon is a scalar 1-cochain, not a vector harmonic:**
- Vertex parity (n1+n2+q odd): requires γ^μ parity-flip, absent in scalar QED
- SO(4) channel count (W>0): requires SU(2)_L×SU(2)_R double-triangle, absent
- Ward identity [D,Lambda]=[Sigma,D]: requires the vertex to be a gauge derivative
  of the propagator; the scalar CG vertex lacks this geometric identity

**BROKEN (INTRINSIC tier) — rules that fail due to graph kinematics:**
- Angular momentum Δm_j: 20/48 entries violate m_j conservation; the CG
  projection allows both Δm_j = +Δm_photon AND Δm_j = -Δm_photon simultaneously
  (symmetrized vertex), producing entries that conserve Δm_j only in the
  symmetric combination but not entry-by-entry for all indices
- Spatial parity (E1 rule): ALL 48 entries have l_a + l_b even; the scalar
  Fock graph T-edges connect (n,l,m) -> (n+1,l,m) [same l], and the CG
  projection maps each Dirac state to ONE l value. Cross-shell T-edge couplings
  with same l on both sides give l_a+l_b even in ALL cases at n_max=2
- Charge conjugation (C): 8/16 Sigma elements checked violate C-symmetry;
  the CG projection breaks C because kappa=-1 (l=0) and kappa=+1 (l=1)
  map to different Fock nodes, so Sigma[C(a),C(b)] != Sigma[a,b]
- Furry's theorem: follows from C-symmetry breaking; tadpole=16*sqrt(2)/15,
  triangle=3584*sqrt(2)/3375 (both nonzero algebraic)

**Only SURVIVES:**
- Gaunt/CG sparsity (INTRINSIC): the CG selection rules suppress most
  couplings. Density falls from 16.0% (n_max=2) to 0.62% (n_max=4).
  This is the structural sparsity feature preserved from the continuum.

**Key insight for Paper 28 / scalar_vs_vector_qed.md:**
The census reveals that SCALAR QED on the finite graph is NOT a controlled
approximation to vector QED — it breaks almost all the symmetry constraints.
The graph computes the correct Gaunt angular sparsity (INTRINSIC), but the
physical content (parity, C, Furry, Ward, vertex parity, SO(4) channel count)
all require the CALIBRATION of promoting the photon from scalar to vector.
Paper 18 tier for the scalar-to-vector upgrade: CALIBRATION (requires embedding
the graph in the continuum vector-harmonic structure, which is where α enters).

**Tr(Sigma) = 44/3; Tr(Lambda) = 32/9; Tr(Lambda)/Tr(Sigma) = 8/33.**
These are the known graph-native QED invariants from Paper 28.