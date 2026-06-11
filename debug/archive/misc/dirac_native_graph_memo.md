# Native Dirac Graph Analysis -- Direction 4 Investigation

**Date:** 2026-04-24
**Status:** NEGATIVE -- DO NOT BUILD

## Question

Should GeoVac build qubit Hamiltonians directly on a native Dirac graph
with nodes (n, kappa, m_j), replacing the current three-step scalar
pipeline (scalar graph -> perturbative SO -> (kappa, m_j) relabeling)?

## Background

The current relativistic composed pipeline (T3, `composed_qubit_relativistic.py`)
follows a three-step process:

1. Build the scalar graph Laplacian h1 on nodes (n, l, m) with off-diagonal
   kappa = -1/16 inter-shell couplings
2. Add perturbative Breit-Pauli spin-orbit as a diagonal correction in the
   (kappa, m_j) basis
3. Build jj-coupled two-body ERIs using X_k angular coefficients on (kappa, m_j) labels

Direction 4 asks whether a single-step pipeline -- nodes natively labeled
(n, kappa, m_j) with h1 from the Dirac-Coulomb spectrum -- would be
structurally advantageous.

## Results

### 1. Node/qubit count: IDENTICAL

| n_max | Scalar spatial | Scalar Q | Dirac Q | Ratio |
|------:|---------------:|---------:|--------:|------:|
|     1 |              1 |        2 |       2 | 1.000 |
|     2 |              5 |       10 |      10 | 1.000 |
|     3 |             14 |       28 |      28 | 1.000 |
|     4 |             30 |       60 |      60 | 1.000 |
|     5 |             55 |      110 |     110 | 1.000 |

This is a representation-theoretic identity, not a coincidence. For each
(n, l), the scalar basis has 2(2l+1) spin-orbitals; the Dirac basis has
(2l+2) from kappa=-(l+1) plus 2l from kappa=+l (when l >= 1), totaling
4l+2 = 2(2l+1). The qubit register size is identical by construction.

### 2. One-body structure: SCALAR WINS

**Scalar h1:** Graph Laplacian with off-diagonal kappa = -1/16 T-plus/T-minus
couplings between adjacent n-shells within each (l, m) channel, plus
diagonal -Z^2/(2n^2). The off-diagonal topology IS the GeoVac framework
identity -- it encodes the inter-shell transitions that make graph-native
CI work.

**Dirac-Coulomb h1:** EXACTLY diagonal with E_Dirac(n, kappa). No
off-diagonal couplings at all. Fine structure is built into the diagonal,
but the graph topology is completely absent. The one-body sector reduces
to an energy list, not a graph.

He (Z=2) one-body comparison at n_max=2:

| State         | NR E (Ha)  | Dirac E (Ha) | SO shift   |
|:--------------|:-----------|:-------------|:-----------|
| n=1,l=0,j=0.5| -2.000000  | -2.000027    | -2.67e-05  |
| n=2,l=0,j=0.5| -0.500000  | -0.500004    | -4.18e-06  |
| n=2,l=1,j=1.5| -0.500000  | -0.500001    | -5.57e-07  |
| n=2,l=1,j=0.5| -0.500000  | -0.500002    | -1.67e-06  |

The fine-structure corrections are O(alpha^2) ~ 5e-5 relative. For quantum
simulation purposes this is negligible.

**Critical insight:** The scalar graph's kappa = -1/16 off-diagonal couplings
are the DEFINING FEATURE of GeoVac (Paper 0, Paper 7). The graph-native CI
(Sprint 3C, 0.19% He accuracy) works BECAUSE of these couplings -- they
encode the inter-shell mixing that captures correlation. Removing them in
favor of an exactly diagonal h1 removes the framework's identity.

### 3. Two-body ERI density

ERI density comparison from T0 (pair-diagonal convention):

| l_max | d_scalar | d_spinor | ratio  |
|------:|---------:|---------:|-------:|
|     0 |   1.0000 |   0.2500 |  0.250 |
|     1 |   0.0781 |   0.0429 |  0.550 |
|     2 |   0.0276 |   0.0211 |  0.764 |
|     3 |   0.0144 |   0.0123 |  0.856 |
|     4 |   0.0090 |   0.0081 |  0.898 |
|     5 |   0.0062 |   0.0057 |  0.925 |

The spinor ERI density is LOWER than scalar at every l_max (ratio < 1).
However, this density advantage is offset by a LARGER orbital index space
(both p and q run over all spinor indices, not spatial + spin separately),
resulting in more total nonzero Pauli terms.

### 4. Concrete Pauli counts: He at n_max=2

| Pipeline                   | Q  | N_Pauli | Ratio vs scalar |
|:---------------------------|---:|--------:|----------------:|
| Scalar (n,l,m)             | 10 |     111 |            1.00 |
| Dirac labels + NR h1       | 10 |     471 |            4.24 |
| Dirac labels + Dirac h1    | 10 |     471 |            4.24 |

The Pauli count INCREASES 4.24x in the Dirac basis compared to the
scalar basis. This ratio matches the known post-TR T3 ratio exactly.

The h1 choice (NR vs exact Dirac) does not affect the Pauli count at all --
the 4.24x inflation comes entirely from the jj-coupled angular coefficients
X_k having more nonzero entries than the scalar Gaunt coefficients c^k.

### 5. Three interpretations of "native Dirac graph"

**Interpretation 1 (Minimal):** Just swap h1 diagonal to exact Dirac-Coulomb
energies. Same Q, same ERIs, same Pauli count. Gain is exact fine structure
at O(alpha^2), which is a 1-line code change and not worth a new architecture.

**Interpretation 2 (Moderate):** Build a Dirac GRAPH with off-diagonal
couplings. Would need to identify what plays the role of kappa = -1/16
in the Dirac sector. The Dirac operator on S^3 is first-order (eigenvalues
n + 3/2), while the scalar Laplacian is second-order (eigenvalues n^2 - 1).
The natural Dirac coupling (sigma dot r-hat) connects kappa to -kappa WITHIN
each n-shell (intra-shell), not BETWEEN shells (inter-shell) like the scalar
kappa couplings. No known inter-shell Dirac graph coupling constant exists.

**Interpretation 3 (Full):** Use the Ihara zeta Rule B adjacency (E1 dipole,
parity-flip Delta-l = +/-1) from Sprint RH-C. This is a genuinely new Dirac-
specific graph structure -- but finding the coupling constant (the Dirac
kappa-analog) is a research investigation, not an engineering task. The
MEMORY.md note ("pi-free Dirac graph viable; future architecture sprint")
refers to this open question.

## Verdict

**DO NOT BUILD** a native Dirac graph Hamiltonian as an engineering task.

**Reasons:**

1. Same qubit count -- no savings from relabeling
2. Dirac h1 is exactly diagonal -- loses the graph topology that defines GeoVac
3. Pauli counts INCREASE 4.24x from jj-coupled angular coefficients
4. Current scalar + perturbative SO already uses (kappa, m_j) labels for ERIs
5. Fine structure gain is O(alpha^2), negligible for quantum simulation
6. The scalar graph off-diagonal kappa couplings are essential for graph-native CI

**What to do instead:**

- **Keep the current approach:** scalar graph + perturbative SO (production pipeline)
- **If exact fine structure needed:** swap h1 diagonal to E_Dirac (1-line change, no architectural work)
- **The real Dirac graph direction** is finding the Dirac kappa-analog coupling
  constant for the Ihara zeta Rule B adjacency -- but this is a research
  investigation (MEMORY.md "future architecture sprint"), not engineering

## Files

- Analysis script: `debug/dirac_native_graph_analysis.py`
- Results data: `debug/data/dirac_native_graph_analysis.json`
- This memo: `debug/dirac_native_graph_memo.md`
