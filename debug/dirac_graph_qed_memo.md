# Native Dirac Graph QED: Sprint Memo

**Date:** 2026-04-27
**n_max (primary):** 2
**n_max (extended):** 3, 4 (pendant analysis)

## Summary

This sprint builds QED directly on the Dirac graph with nodes (n_fock, kappa, m_j),
bypassing the CG projection used in the scalar-Fock-graph QED (GN-1 through GN-7).
Two adjacency rules are tested:

- **Rule A (kappa-preserving):** Edges at fixed kappa, lifting scalar T+/- and L+/- ladders.
  Graph splits into kappa-sectors (analogous to l-shell splitting of scalar graph).
- **Rule B (E1 dipole):** Edges enforce Delta l = +/-1 (parity flip), |Delta j| <= 1,
  |Delta m_j| <= 1. Mixes kappa. Single connected component.

## Task 1: Hodge Decomposition

### Comparison Table

| Graph | n_max | V | E | beta_0 | beta_1 | GS pendant | GS degree |
|:------|------:|--:|--:|-------:|-------:|:------------|----------:|
| Scalar Fock | 2 | 5 | 3 | 2 | 0 | YES (deg=1) | 1 |
| Scalar Fock | 3 | 14 | 13 | 3 | 2 | YES (deg=1) | 1 |
| Dirac Rule A | 2 | 10 | 8 | 3 | 1 | NO | 2 |
| Dirac Rule A | 3 | 28 | 29 | 5 | 6 | NO | 2 |
| Dirac Rule B | 2 | 10 | 20 | 1 | 11 | NO | 5 |
| Dirac Rule B | 3 | 28 | 106 | 1 | 79 | NO | 5 |

### Exact L0 and L1 Eigenvalues (n_max=2)

**Rule A:**
- L0 eigenvalues: {0: 3, 2: 4, 2-sqrt(2): 1, 2+sqrt(2): 1, 4: 1}
- L1 eigenvalues: {0: 1, 2: 4, 2-sqrt(2): 1, 2+sqrt(2): 1, 4: 1}
- SVD identity: VERIFIED (nonzero spectra match)
- pi-free: YES (algebraic over Q, contains sqrt(2))

**Rule B:**
- L0 eigenvalues: {0: 1, 7/2-sqrt(17)/2: 1, 7/2+sqrt(17)/2: 1, 4: 3, 5: 2, 11/2-sqrt(41)/2: 1, 11/2+sqrt(41)/2: 1}
- L1 eigenvalues: same nonzero + {0: 11}
- SVD identity: VERIFIED
- pi-free: YES (algebraic over Q, contains sqrt(17), sqrt(41))

### Structural Finding: GS is NOT a Pendant

The ground state (n_fock=1, kappa=-1, m_j=+/-1/2) is NEVER a pendant vertex
(degree 1) on either Dirac graph, at any n_max tested (2, 3, 4).

- **Rule A:** GS degree = 2 (constant). The two GS nodes connect via the
  m_j ladder (m_j=+1/2 <-> m_j=-1/2) AND via the n-ladder (n=1 <-> n=2
  at same kappa=-1). But there are TWO connections per GS node, not one.
- **Rule B:** GS degree = 5 (constant). The dipole selection rule connects
  each 1s_{1/2} state to all accessible p_{1/2} and p_{3/2} states at
  n_fock=1,2 with Delta l = 1, |Delta m_j| <= 1.

This is a **structural difference** from the scalar Fock graph, where the
GS node |1,0,0> has degree 1 (single T+ edge to |2,0,0>). The scalar
pendant-edge theorem Sigma(GS) = 2(n_max-1)/n_max does NOT transfer.

## Task 2: Photon Propagator

The photon propagator G_gamma = L1^+ (Moore-Penrose pseudoinverse of the
edge Laplacian) is computed exactly via sympy at n_max=2.

| Graph | E | rank(G_gamma) | beta_1 (gauge modes) |
|:------|--:|:--------------|---------------------:|
| Scalar Fock (n_max=2) | 3 | 3 | 0 |
| Dirac Rule A (n_max=2) | 8 | 7 | 1 |
| Dirac Rule B (n_max=2) | 20 | 9 | 11 |

Rule B has 11 gauge zero modes (beta_1 = 11), much richer than Rule A (1)
or the scalar graph (0). This means the Dirac graph has substantial
cycle structure in its photon sector.

## Task 3: Electron Propagator

At t=0 (free Dirac propagator), D_Dirac = diag(lambda_i) with
lambda_i = chi * (n_CH + 3/2), where chi = +1 for kappa < 0, -1 for kappa > 0.
G_e = diag(1/lambda_i), purely rational in Q.

The key observation: at t=0, G_e is diagonal, so the bare vertex trace
Tr(V_bare @ G_e) = 0 because V_bare is purely off-diagonal (identity
coupling on edges). This means F2 extraction requires t != 0 or a
different normalization convention. This is structurally different from
the scalar Fock graph where the CG projection creates diagonal entries
in V_bare.

## Task 4: Self-Energy

### Exact Self-Energy (n_max=2)

| Quantity | Scalar Fock | Dirac Rule A | Dirac Rule B |
|:---------|:------------|:-------------|:-------------|
| Tr(Sigma) | 44/3 | **17/2** | **103/20** |
| GS block zero | NO | NO | NO |
| GS block | [[1,1],[1,1]] | **[[5/8, 0], [0, 5/8]]** | **[[103/160, 21/160], [21/160, 103/160]]** |
| is_rational | mixed | YES | YES |
| is_pi_free | YES | YES | YES |
| is_Hermitian | YES | YES | YES |
| is_PSD | YES | YES | YES |

**Key finding: GS block is NOT zero on either Dirac graph.** The structural
zero (Sigma(GS) = 0) from the continuum vertex parity rule is NOT recovered.
However, the structure differs:

- **Scalar Fock:** GS block = [[1,1],[1,1]] — rank-1, degenerate,
  from CG-projected identity coupling
- **Rule A:** GS block = diag(5/8, 5/8) — DIAGONAL, no m_j mixing.
  This is because Rule A preserves kappa AND edges at fixed kappa do
  not mix m_j (the m_j ladder connects m_j+1 to m_j-1, contributing
  symmetrically). The value 5/8 is rational.
- **Rule B:** GS block = [[103/160, 21/160], [21/160, 103/160]] —
  off-diagonal m_j mixing from dipole edges that connect different
  m_j values across kappa sectors.

### Exact Eigenvalues of Sigma (n_max=2)

**Rule A:** {1/2: 4, 3/4: 2, 5/4 - sqrt(2)/2: 2, 5/4 + sqrt(2)/2: 2}
- Algebraic but NOT rational (contains sqrt(2))
- No zero eigenvalues — Sigma is strictly positive definite

**Rule B:** {1/5: 2, 1/4: 3, 49/100: 1, 267/400 - 3*sqrt(881)/400: 1,
             267/400 + 3*sqrt(881)/400: 1, 33/40: 1, 27/20: 1}
- Algebraic but NOT rational (contains sqrt(881))
- No zero eigenvalues — Sigma is strictly positive definite

**Comparison:** The scalar Fock graph Sigma has 5 zero eigenvalues (matching
the L1 kernel). The Dirac graph Sigma has NO zero eigenvalues for either
rule, meaning the photon's gauge structure does not protect any electron
subspace from self-energy corrections.

### Exact Vertex Correction (n_max=2, t=0)

**Rule A:**
- Tr(Lambda) = 44/15
- Eigenvalues: {-1/5: 2, 4/15: 2, 2/5: 2, 1/2 - sqrt(2)/5: 2, 1/2 + sqrt(2)/5: 2}
- Lambda has NEGATIVE eigenvalue (-1/5, multiplicity 2) — Lambda is indefinite
- Contains sqrt(2)

**Rule B:**
- Tr(Lambda) = 662/375
- Eigenvalues: {2/125: 2, 2/15: 3, 9/50: 2, 89/250 - sqrt(881)/250: 1,
                89/250 + sqrt(881)/250: 1, 98/375: 1}
- Lambda is positive semidefinite (no negative eigenvalues)
- Contains sqrt(881)

### n_max=3 Self-Energy (Rule B, numeric)

- Tr(Sigma) = 9.620 (vs 103/20 = 5.15 at n_max=2)
- GS block: [[0.259, 0.018], [0.018, 0.259]] — off-diagonal shrinks
- GS block NOT zero

### Pendant Analysis (n_max=2,3,4)

| Rule | n_max | GS degree | GS pendant |
|:-----|------:|----------:|:-----------|
| A | 2 | 2 | NO |
| A | 3 | 2 | NO |
| A | 4 | 2 | NO |
| B | 2 | 5 | NO |
| B | 3 | 5 | NO |
| B | 4 | 5 | NO |

GS degree is CONSTANT (independent of n_max) for both rules. This is
because the GS connects only to states at n_fock=1 and n_fock=2 via
either the m_j ladder (Rule A) or the dipole E1 rule (Rule B), and
these connections do not change as higher shells are added.

## Task 5: Selection Rule Census

### Comparison Table (n_max=2)

| Rule | Scalar Fock | Dirac A | Dirac B |
|:-----|:------------|:--------|:--------|
| 1. Delta m_j conservation | BROKEN (20 violations) | **SURVIVES** | **SURVIVES** |
| 2. Spatial parity (E1 l_a+l_b odd) | BROKEN (all even) | BROKEN (all even) | **SURVIVES** (all odd) |
| 3. Gaunt/CG sparsity | SURVIVES (84% sparse) | SURVIVES (82% sparse) | SURVIVES (56% sparse) |
| 4. Vertex parity (n1+n2+q odd) | BROKEN (0% forbidden) | BROKEN (mixed) | BROKEN (mixed) |
| 5. SO(4) channel count (W>0) | BROKEN (all W=0) | N/A (structural) | N/A (structural) |
| 6. Charge conjugation (C) | BROKEN (8 violations) | BROKEN | BROKEN |
| 7. Furry (odd-loop = 0) | BROKEN (tadpole != 0) | **SURVIVES** (tadpole=0) | **SURVIVES** (tadpole=0) |
| 8. Ward identity | BROKEN (5.08) | BROKEN | BROKEN |
| **Total surviving** | **1/8** | **3/8** | **4/8** |

### Detailed Rule Analysis

**Rule 1 (Delta m_j): RECOVERED on both Dirac graphs.**
On the Dirac graph, edges directly connect specific (n,kappa,m_j) states.
By construction, both Rule A (|Delta m_j| <= 1) and Rule B (|Delta m_j| <= 1)
enforce angular momentum conservation at every edge. The scalar Fock graph
breaks this because the CG projection creates entries that couple Dirac
states with different m_j through the same Fock edge.

**Rule 2 (Spatial parity): RECOVERED on Rule B only.**
Rule B enforces Delta l = +/-1 at every edge (parity flip). This means
l_a + l_b is always ODD, which is exactly the E1 selection rule. Rule A
preserves kappa (hence l), so l_a = l_b and l_a + l_b is always EVEN —
same-parity transitions, violating E1.

**Rule 7 (Furry's theorem): RECOVERED on both Dirac graphs.**
The tadpole vanishes because G_e(t=0) is diagonal and V_e is off-diagonal:
Tr(sum_e V_e @ G_e) = sum_{(i,j) edge} (G_e[i,j] + G_e[j,i]) = 0.
The triangle also vanishes: odd-loop traces of products of off-diagonal
matrices with a diagonal matrix are zero by parity.
**Mechanism:** The identity vertex on the Dirac graph has a natural
parity structure (odd under graph-transposition) that the CG-projected
scalar vertex lacks. This is NOT the continuum vertex parity (which
requires vector photon quantum numbers) but a graph-topological parity.

**Rule 6 (Charge conjugation): STILL BROKEN on both graphs.**
C-conjugation maps (n, kappa, m_j) -> (n, -kappa, -m_j). The self-energy
Sigma is NOT C-symmetric because the Dirac eigenvalues have opposite
signs for opposite kappa: lambda(kappa<0) > 0, lambda(kappa>0) < 0.
The photon propagator treats edges between different kappa sectors
asymmetrically through the graph topology.

**Rule 4 (Vertex parity): STILL BROKEN on both graphs.**
The Dirac graph has no explicit photon quantum number q. The identity
vertex carries no angular momentum transfer quantum number — it simply
couples the two endpoints of an edge. Vertex parity n1+n2+q odd requires
a vector photon with definite q, which is absent.

## Structural Interpretation

### What Changed: Scalar Fock -> Dirac Graph

The Dirac graph recovers 3 (Rule A) or 4 (Rule B) selection rules out of 8,
compared to 1/8 on the scalar Fock graph. The recovered rules are:

1. **Delta m_j conservation** — because edges directly connect specific
   (kappa, m_j) states instead of going through CG projection
2. **Spatial parity** (Rule B only) — because E1 dipole edges enforce
   Delta l = +/-1 by construction
3. **Furry's theorem** — because the identity vertex on the Dirac graph
   is off-diagonal (odd parity), making odd-loop traces vanish at t=0

### What Did NOT Change

The structural zero (Sigma(GS) = 0) is NOT recovered. The GS is no longer
a pendant vertex (degree 2 or 5 vs degree 1 on scalar graph), so the
pendant-edge theorem does not apply. However, the GS block is still nonzero
because the photon propagator G_gamma connects edges involving the GS node
to themselves (the GS node participates in multiple edges, each of which
contributes to Sigma through V_e @ G_gamma @ V_e'^T).

The vertex parity, SO(4) channel count, charge conjugation, and Ward
identity remain broken. These require:
- **Vertex parity, SO(4):** vector photon with definite angular momentum
  quantum numbers (the Dirac graph photon is still a scalar 1-cochain)
- **C symmetry:** requires the photon propagator to be C-invariant, which
  fails because the graph topology is kappa-asymmetric
- **Ward identity:** requires the vertex to be the gauge derivative of the
  propagator; the identity vertex is not a gauge derivative

### Paper 18 Classification

The Dirac graph QED lives in the INTRINSIC tier:
- All eigenvalues are algebraic (Q[sqrt(2), sqrt(17), sqrt(41), sqrt(881)])
- No transcendentals (pi, zeta) anywhere
- Exact rational traces: Tr(Sigma_A) = 17/2, Tr(Sigma_B) = 103/20,
  Tr(Lambda_A) = 44/15, Tr(Lambda_B) = 662/375

The algebraic irrationals come from the graph topology (Laplacian
characteristic polynomial roots), NOT from CG coefficients (which are
absent). The Dirac graph algebraic content is simpler than the scalar
Fock graph (Q[sqrt(2)] for Rule A vs Q[sqrt(2), sqrt(3), sqrt(6)] for
scalar Fock).

### Significance for Paper 28

**The Dirac graph is an intermediate construction** between the scalar
Fock graph (1/8 rules) and the full continuum vector QED (8/8 rules).
It recovers the rules that depend on spinor quantum numbers (Delta m_j,
spatial parity, Furry) but cannot recover the rules that require a
vector photon (vertex parity, SO(4), Ward).

This cleanly separates the selection rule failures into two categories:
- **Spinor-recoverable (3 rules):** fixed by using Dirac nodes
  instead of scalar Fock nodes
- **Vector-photon-required (4 rules):** require promoting the photon
  from a scalar 1-cochain to a vector harmonic

The charge conjugation failure is in between: it depends on the
electron sector asymmetry (kappa sign), not on the photon structure,
but is also not fixed by simply using Dirac nodes.

## Exact Algebraic Data

### Self-Energy Eigenvalues (n_max=2)

| Rule A | Multiplicity |
|:-------|:------------|
| 1/2 | 4 |
| 3/4 | 2 |
| 5/4 - sqrt(2)/2 | 2 |
| 5/4 + sqrt(2)/2 | 2 |

| Rule B | Multiplicity |
|:-------|:------------|
| 1/5 | 2 |
| 1/4 | 3 |
| 49/100 | 1 |
| 33/40 | 1 |
| 27/20 | 1 |
| 267/400 - 3*sqrt(881)/400 | 1 |
| 267/400 + 3*sqrt(881)/400 | 1 |

### Vertex Correction Eigenvalues (n_max=2, t=0)

| Rule A | Multiplicity |
|:-------|:------------|
| -1/5 | 2 |
| 4/15 | 2 |
| 2/5 | 2 |
| 1/2 - sqrt(2)/5 | 2 |
| 1/2 + sqrt(2)/5 | 2 |

| Rule B | Multiplicity |
|:-------|:------------|
| 2/125 | 2 |
| 2/15 | 3 |
| 9/50 | 2 |
| 98/375 | 1 |
| 89/250 - sqrt(881)/250 | 1 |
| 89/250 + sqrt(881)/250 | 1 |

## Key Files

- `debug/dirac_graph_qed_sprint.py` — Main computation script (all 5 tasks)
- `debug/dirac_graph_qed_exact.py` — Exact sympy verification
- `debug/data/dirac_graph_qed_sprint.json` — Full float64 results
- `debug/data/dirac_graph_qed_exact.json` — Exact algebraic results
- `debug/dirac_graph_qed_memo.md` — This memo
