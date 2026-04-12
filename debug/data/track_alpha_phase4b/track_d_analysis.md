# Track alpha-D: Hopf Projection as Graph Morphism

## Construction

- **S^3 graph**: GeoVac `GeometricLattice(max_n=3, Z=1)`. Nodes are (n,l,m) for 1<=n<=3, 0<=l<n, -l<=m<=l (14 nodes). Edges from Paper 7 selection rules: angular Delta m=+/-1 within (n,l), radial Delta n=+/-1 with same (l,m).
- **S^2 quotient**: collapse 2l+1 m-substates of each (n,l) into one node labeled (n,l); 6 nodes with multiplicities 1,1,3,1,3,5 (sum=14). Intersector adjacency counts S^3 edges crossing sectors (these are exactly the radial Delta n=+/-1 edges). Intra-sector edges are angular L+/- transitions and become the fiber.
- **S^1 fibers**: each (n,l) sector with l>=1 has 2l+1 m-states connected by angular L+/- transitions. Paper 7 selection rules give a PATH graph P_{2l+1}, not a cycle. Cycle alternative computed for completeness.

## S^3 graph at n_max = 3

Nodes: 14, edges: 13

Laplacian eigenvalues (S^3): [0.0, 0.0, 0.0, 0.3819660112501052, 0.9999999999999998, 1.0000000000000007, 1.3819660112501049, 1.999999999999999, 2.6180339887498953, 3.0, 3.0, 3.0000000000000004, 3.6180339887498945, 4.999999999999999]

## S^2 quotient

Intersector adjacency A_S2:
```
     0    1    0    0    0    0
     1    0    0    1    0    0
     0    0    0    0    3    0
     0    1    0    0    0    0
     0    0    3    0    0    0
     0    0    0    0    0    0
```
Laplacian eigenvalues (S^2): [0.0, 0.0, 0.0, 0.9999999999999998, 3.0, 6.0]

## Fiber Laplacian (path-graph) eigenvalue sums

- [2, 1]: eigs=[0.0, 0.9999999999999998, 3.0], logdet'=1.0986122886681096
- [3, 1]: eigs=[0.0, 0.9999999999999998, 3.0], logdet'=1.0986122886681096
- [3, 2]: eigs=[0.0, 0.3819660112501052, 1.3819660112501049, 2.6180339887498953, 3.6180339887498945], logdet'=1.6094379124341005
Sum of fiber log det' (path): 3.8066624897703196

## Spectral invariants comparison

| Invariant | S^3 | S^2 |
|---|---|---|
| n_nodes | 14 | 6 |
| zero_modes | 3 | 3 |
| trace | 26.0 | 10.0 |
| zeta_1 | 7.699999999999999 | 1.5000000000000002 |
| zeta_2 | 10.223333333333331 | 1.1388888888888893 |
| zeta_3 | 20.644111111111105 | 1.0416666666666672 |
| logdet_prime | 7.2078598714324755 | 2.8903717578961645 |
| det_prime | 1350.0 | 17.999999999999996 |
| alg_connectivity | 0.3819660112501052 | 0.9999999999999998 |
| largest_eig | 4.999999999999999 | 6.0 |
| spectral_entropy | 2.2360712934207645 | 0.8979457248567797 |

## Graph exchange constant candidates

- zeta1_diff = 6.199999999999999
- logdet_diff = 4.317488113536311
- vn_entropy_diff = 1.3381255685639848
- fiber_logdet_sum_path = 3.8066624897703196
- fiber_logdet_sum_cycle = 7.613324979540639

## Candidates vs targets (rel < 5%)

| Candidate | Value | Target | Target value | Rel err |
|---|---|---|---|---|
| casimir_weighted_trace | 42 | B_42 | 42 | 0.000e+00 |
| pi_times_(casimir+F-Delta) | 137.036 | pi*(B+F-Delta)_137.0357 | 137.036 | 0.000e+00 |
| pi_times_(casimir+F-Delta) | 137.036 | K_alpha_137.036 | 137.036 | 4.767e-07 |
| casimir_times_pi_plus_F | 137.115 | pi*(B+F-Delta)_137.0357 | 137.036 | 5.731e-04 |
| casimir_times_pi_plus_F | 137.115 | K_alpha_137.036 | 137.036 | 5.736e-04 |
| S2_logdet_plus_42 | 44.8904 | B+F_43.6449 | 43.6449 | 2.854e-02 |
| S2_logdet_plus_42 | 44.8904 | B+F-Delta_43.6199 | 43.6199 | 2.913e-02 |
| S2_logdet_plus_42 | 44.8904 | K/pi_43.6199 | 43.6199 | 2.913e-02 |
| casimir_times_pi | 131.947 | K_alpha_137.036 | 137.036 | 3.714e-02 |
| casimir_weighted_trace | 42 | K/pi_43.6199 | 43.6199 | 3.714e-02 |
| casimir_weighted_trace | 42 | B+F-Delta_43.6199 | 43.6199 | 3.714e-02 |
| casimir_times_pi | 131.947 | pi*(B+F-Delta)_137.0357 | 137.036 | 3.714e-02 |
| casimir_weighted_trace | 42 | B+F_43.6449 | 43.6449 | 3.769e-02 |

## Verdict

AMBIGUOUS / TAUTOLOGICAL. Interpretation of the results:

1. **The raw graph Laplacian spectral invariants of L_{S^3} and L_{S^2} do NOT encode K.** Eigenvalues are small integers and golden ratios (Fibonacci values from the n-chain paths). logdet'(L_{S^3}) = 7.208, logdet'(L_{S^2}) = 2.890, zeta_1(L_{S^3}) = 7.7, zeta_1(L_{S^2}) = 1.5. None of these hit K, K/pi, B=42, or B+F-Delta to better than 3% in any natural combination.

2. **The 'graph exchange constant' (logdet/zeta/entropy differences between S^3 and S^2) does not hit any target.** logdet_diff = 4.317, zeta1_diff = 6.2, VN entropy diff = 1.338, fiber logdet sum = 3.807. None are close to 42, 43.62, pi^2/6, or 137.

3. **The integer B = 42 IS recovered exactly**, but as the degeneracy-weighted Casimir sum Sum_{(n,l)} (2l+1)*l(l+1) truncated at n_max=3, summing contributions 0+0+6+0+6+30 = 42. This is the quotient-space Casimir trace, not a graph Laplacian spectral invariant. It reproduces Paper 2's B = 42 construction exactly and therefore gives K = pi*(B+F-Delta) to 4.8e-7 relative error. **But this is Paper 2's existing formula restated in graph language, not a new derivation.** The tautology is: B=42 comes from SO(3) representation theory of the S^2 base, which is identical to weighting (n,l) sectors by (2l+1)*l(l+1). Calling the sectors 'nodes of the quotient graph' does not change the computation.

4. **Conclusion**: the Hopf-quotient spectral route does not produce an independent derivation of K at n_max=3. The graph morphism is well-defined, but its spectral invariants live in a different algebraic world (small integers, Fibonacci) from the Paper 2 targets (42, pi^2/6, 1/40, 137). The only combination that hits is the one that bypasses the graph Laplacian entirely and uses the SO(3) Casimir directly -- which is Paper 2's original construction.

Recommendation: **CLEAN NEGATIVE for the 'spectral invariants of the graph quotient encode K' hypothesis.** The honest partial result is that the S^2-quotient node set carries the SO(3) Casimir weighting, so B=42 is naturally expressed as a weighted trace on the quotient; this is a re-labeling of Paper 2, not new physics. No independent prediction of K or K/pi from Laplacian spectra was obtained.