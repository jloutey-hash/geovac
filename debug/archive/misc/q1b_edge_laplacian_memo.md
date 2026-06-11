# Track Q1-B: Edge Laplacian of the GeoVac S^3 Graph at n_max=3

**Date:** 2026-04-15
**Verdict:** NEGATIVE (clean)

## Graph topology

The n_max=3 GeoVac graph has 14 nodes (quantum states (n,l,m)) and 13 edges
(angular m-ladder + radial n-transitions). The graph is **disconnected** with
3 connected components:

1. **s-chain** (3 nodes): (1,0,0) -- (2,0,0) -- (3,0,0) [path P_3]
2. **p-block** (6 nodes): n=2,3 with l=1 [grid with angular + radial edges]
3. **d-chain** (5 nodes): (3,2,-2) -- ... -- (3,2,2) [path P_5]

This gives first Betti number beta_1 = M - N + c = 13 - 14 + 3 = 2 (two
independent cycles, both in the p-block where the angular and radial edges
form a cycle).

## Edge Laplacian spectrum (exact)

L_edge = B^T B is a 13x13 integer matrix. Its characteristic polynomial factors
completely over Q(sqrt(5)):

    lambda^2 (lambda-5)(lambda-3)^3 (lambda-2)(lambda-1)^2 (lambda^2 - 5*lambda + 5)(lambda^2 - 3*lambda + 1)

Eigenvalues:

| Exact value       | Numerical    | Mult | Origin          |
|:-------------------|:-------------|:----:|:----------------|
| 0                  | 0            | 2    | beta_1 = 2      |
| (3 - sqrt(5))/2   | 0.38197      | 1    | P_5 (d-chain)   |
| 1                  | 1.00000      | 2    | P_3 (s-chain)   |
| (5 - sqrt(5))/2   | 1.38197      | 1    | P_5 (d-chain)   |
| 2                  | 2.00000      | 1    | p-block         |
| (3 + sqrt(5))/2   | 2.61803      | 1    | P_5 (d-chain)   |
| 3                  | 3.00000      | 3    | P_3 + p-block   |
| (5 + sqrt(5))/2   | 3.61803      | 1    | P_5 (d-chain)   |
| 5                  | 5.00000      | 1    | p-block         |

The four golden-ratio eigenvalues (involving sqrt(5)) come from the path
graph P_5 of the d-orbital chain.

## Spectral invariants

| Invariant               | Exact     | Numerical    |
|:------------------------|:----------|:-------------|
| Tr(L_edge)              | 26        | 26           |
| Product of nonzero eigs | 1350      | 1350         |
| zeta_edge(1)            | 77/10     | 7.700        |
| zeta_edge(2)            | 3067/300  | 10.223       |
| zeta_edge(3)            | 185797/9000 | 20.644     |
| zeta_edge(4)            | 13332907/270000 | 49.381 |

**All spectral invariants lie in Q(sqrt(5))** -- they are algebraic numbers
with no transcendental content. This is a structural fact: the edge Laplacian
is an integer matrix, so its eigenvalues are algebraic by the fundamental
theorem of algebra.

## Checks against K ingredients

| Target              | Value     | Nearest invariant  | Rel error |
|:--------------------|:----------|:-------------------|:----------|
| B = 42              | 42        | Trace = 26         | 38%       |
| F = pi^2/6          | 1.6449    | 42/Tr = 21/13      | 1.8%      |
| Delta = 1/40        | 0.025     | none within 5%     | --        |
| B + F - Delta       | 43.620    | none within 5%     | --        |
| K = pi(B+F-Delta)   | 137.036   | none within 5%     | --        |

The sole near-miss (42/Tr = 21/13 = 1.6154 vs F = 1.6449) is a rational
approximation to a transcendental number -- a coincidence, not a structural
identity. 21/13 is a Fibonacci-fraction approximation to the golden ratio,
which is present in the spectrum through the P_5 path graph.

## Node-edge spectral relationship

The 11 nonzero eigenvalues of L_edge = B^T B match the 11 nonzero eigenvalues
of L_node = B B^T to machine precision (max difference 8.9e-16), confirming
the SVD theorem. L_node has 3 zero eigenvalues (one per connected component),
while L_edge has 2 (first Betti number).

## Structural conclusion

The edge Laplacian cannot contain F = pi^2/6 or any transcendental K
ingredient because its spectrum is algebraic by construction. This is
consistent with the Phase 4B-4G structural decomposition:

- **B = 42** lives on the node Laplacian (Casimir-weighted sector trace)
- **F = pi^2/6** enters through the infinite Fock degeneracy Dirichlet
  series D_{n^2}(s=4) = zeta_R(2), which is an analytic continuation
  to infinite shells -- not accessible from any finite graph spectrum
- **Delta = 1/40** is a Dirac mode count at n=3

The edge Laplacian adds no new information beyond the node Laplacian
for this graph: their nonzero spectra are identical. The 2 additional
zero eigenvalues encode the cycle structure (beta_1 = 2) of the p-block,
which is topological but does not connect to K.

**The "photon Laplacian" interpretation does not produce new K ingredients.**
