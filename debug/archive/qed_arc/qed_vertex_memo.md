# QED Vertex Coupling and Two-Loop Sunset on S^3

## Summary

Module `geovac/qed_vertex.py` implements the QED vertex coupling on S^3 and 
the two-loop sunset diagram with vertex selection rules, demonstrating how
the vertex restriction exposes non-trivial transcendental content from the
even/odd parity split of the Dirac spectrum.

## Key Finding: Catalan G and beta(4), NOT zeta(3)

**The previous agent's hypothesis was incorrect.** The even/odd split of the 
Dirac Dirichlet series D(4) does NOT expose zeta(3). Instead, it exposes 
Catalan's constant G = beta(2) and Dirichlet beta(4).

### Exact decomposition (verified by PSLQ at 80 digits)

```
D_even(4) = pi^2/2 - pi^4/24 - 4G + 4*beta(4)
D_odd(4)  = pi^2/2 - pi^4/24 + 4G - 4*beta(4)
D(4)      = pi^2 - pi^4/12        [G and beta(4) cancel]
```

where:
- G = Catalan's constant = sum_{n>=0} (-1)^n/(2n+1)^2 = 0.9159...
- beta(4) = sum_{n>=0} (-1)^n/(2n+1)^4 = 0.9889...

### Mechanism

The even/odd Dirac modes have eigenvalues lambda_{2k} = 2k + 3/2 and 
lambda_{2k+1} = 2k + 5/2. The step-2 sums involve Hurwitz zeta at 
quarter-integer shifts:

```
D_even(s) = 2^{-s} * [8 * zeta(s-2, 3/4) - (1/2) * zeta(s, 3/4)]
D_odd(s)  = 2^{-s} * [8 * zeta(s-2, 5/4) - (1/2) * zeta(s, 5/4)]
```

The Hurwitz zeta at shift 3/4 decomposes via:
- zeta(s, 1/4) + zeta(s, 3/4) = (4^s - 2^s) * zeta_R(s)  [Riemann zeta]
- zeta(s, 1/4) - zeta(s, 3/4) = 4^s * beta(s)  [Dirichlet beta]

So the half-sum gives Riemann zeta (pi^{even}) and the half-difference gives 
Dirichlet beta. The even/odd Dirac split accesses these individually.

### Why NOT zeta(3)?

The Dirichlet beta function beta(s) = L(s, chi_4) is the L-function for the
unique non-principal character mod 4 (chi_4(n) = (-1)^{(n-1)/2} for odd n,
0 for even n). It is a fundamentally different transcendental family from 
Riemann zeta:

- beta(odd) = rational * pi^{odd}  (e.g., beta(1) = pi/4, beta(3) = pi^3/32)
- beta(even) = new transcendentals  (e.g., beta(2) = G, beta(4) = 0.9889...)

The zeta(3) hypothesis failed because PSLQ cannot find D_even(4) in 
{1, pi^2, pi^4, zeta(3)} -- this was interpreted as "PSLQ needs more terms"
but was actually because D_even(4) is genuinely outside that span. The correct
basis is {1, pi^2, pi^4, G, beta(4)}.

## Module Contents

### Functions (public)

| Function | Description |
|:---------|:------------|
| `reduced_coupling_squared(n1, n2, n_gamma)` | Reduced coupling |V|^2 (approximate 6j normalization) |
| `vertex_allowed_triples(n_max_e, n_max_gamma)` | List all allowed (n1, n2, q) triples |
| `two_loop_sunset_unrestricted(n_max, s1, s2)` | Unrestricted = D(2s1) * D(2s2), pi^{even} |
| `two_loop_sunset_vertex_restricted(n_max, s1, s2)` | With vertex + photon propagator |
| `two_loop_vertex_correction(n_max, s1, s2)` | Pair-existence restriction analysis |
| `two_loop_odd_even_split(n_max, s)` | D_even/D_odd decomposition with PSLQ |
| `two_loop_sunset_parity_weighted(n_max, s1, s2)` | Parity-weighted double sum |
| `two_loop_photon_line(n_max, s_e, s_gamma)` | Mixed Dirac + Hodge-1 propagator |
| `decompose_two_loop_result(value)` | Extended PSLQ into zeta basis |
| `verify_vertex_factorization_failure(n_max)` | Documents factorization breaking |
| `flat_space_vertex_sum(N)` | Flat-space nearest-neighbor analog |
| `two_loop_transcendental_classification()` | Paper 18 taxonomy |

### Internal helpers

| Function | Description |
|:---------|:------------|
| `_dirac_D(s)` | Full Dirac Dirichlet series via Hurwitz |
| `_dirac_D_even(s)` | Even-n sub-sum via Hurwitz at 3/4 |
| `_dirac_D_odd(s)` | Odd-n sub-sum via Hurwitz at 5/4 |
| `_dirichlet_beta(s)` | Dirichlet beta function via Hurwitz |
| `_decompose_with_beta_basis(value)` | PSLQ into {1, pi^2, pi^4, G, beta(4)} |

## Paper 18 Taxonomy Update

The vertex two-loop result adds a new transcendental tier to the operator-order grid:

| Tier | Source | Content | Example |
|:-----|:-------|:--------|:--------|
| Rational | Graph topology, selection rules | Q | Gaunt integrals |
| Calibration pi | Second-order operators (D^2) | pi^{even} | D(4) = pi^2 - pi^4/12 |
| Dirichlet beta | Quarter-integer Hurwitz shifts | G, beta(4) | D_even(4), D_odd(4) |
| Odd Riemann zeta | First-order operators (|D|) | zeta(3), zeta(5) | Dirac sector D3 |

The Dirichlet beta tier is structurally distinct from both calibration pi and 
odd Riemann zeta. It arises from the Dirichlet L-function L(s, chi_4) at the 
character chi_4 mod 4, which is the natural arithmetic object at 
quarter-integer Hurwitz shifts.

## SO(4) Clebsch-Gordan Weighted Vertex (v2)

### The physical vertex weight

The QED vertex on S^3 couples Dirac spinor at level n1 (positive chirality,
SO(4) rep ((n1+1)/2, n1/2)) to level n2 (negative chirality, rep (n2/2, (n2+1)/2))
via a transverse vector at level q, which has two SO(4) components:

    V_A = ((q+1)/2, (q-1)/2))    V_B = ((q-1)/2, (q+1)/2))

After summing over all magnetic quantum numbers using the Wigner-Eckart theorem
on each SU(2) factor, the SU(2) CG orthogonality gives a factor of
(2*j_final+1) per channel. So the ONLY n-dependent vertex structure is the
**channel count** W(n1, n2, q) = number of SO(4) components satisfying the
triangle inequality in BOTH SU(2) factors.

    W(n1, n2, q) in {0, 1, 2}

- W=0: forbidden by selection rule
- W=1: only one SO(4) component contributes (boundary of triangle range)
- W=2: both components contribute (interior of triangle range)

### Exact pair-weight formula (verified n1, n2 <= 12)

    W_total(n1, n2) = sum_{q: allowed} W(n1, n2, q) = 2*min(n1,n2) - 1 - delta_{n1,n2}

    N_q(n1, n2) = count of allowed q = min(n1, n2)

The ratio W_total / N_q -> 2 for large min, with O(1/min) correction.

### Transcendental consequence: depth-2 multiple Hurwitz zeta

The min(n1,n2)-weighted double sum factorizes as:

    S_min = sum_{k=1}^inf T(k)^2

where T(k) = 2*zeta(s-2, k+3/2) - (1/2)*zeta(s, k+3/2) is the Dirac tail.

**PSLQ result at 100 digits with 18-element weight-8 basis:**
S_min is NOT identifiable in the span of:
{1, pi^2, pi^4, pi^6, pi^8, G, beta(4), zeta(3), zeta(5), G^2, G*beta(4),
 beta(4)^2, zeta(3)^2, zeta(3)*zeta(5), pi^2*zeta(3), pi^2*G, pi^2*beta(4)}

This means the CG-weighted vertex introduces a **genuinely new transcendental
object** -- a depth-2 multiple Hurwitz zeta value on the half-integer lattice --
that lies outside all standard transcendental constants considered in the
GeoVac framework.

### Paper 18 taxonomy impact

The weighted vertex adds a fifth tier:

| Tier | Source | Content | Example |
|:-----|:-------|:--------|:--------|
| Rational | Graph topology, selection rules | Q | Gaunt integrals, W(n1,n2,q) |
| Calibration pi | Second-order operators (D^2) | pi^{even} | D(4) |
| Dirichlet beta | Quarter-integer Hurwitz | G, beta(4) | D_even(4), D_odd(4) |
| Odd Riemann zeta | First-order operators (\|D\|) | zeta(3), zeta(5) | D3 |
| Depth-2 Hurwitz | Nested spectral sums (CG vertex) | sum T(k)^2 | S_min |

The depth-2 tier arises from the PRODUCT structure of two-loop diagrams: the
single Dirichlet series D(s) decomposes cleanly into Hurwitz + beta, but
the NESTED sum of squared tails does not. This is the spectral analog of the
standard QFT fact that multi-loop diagrams produce irreducible MZVs.

## Test Results

41/41 non-slow tests pass (9 slow tests deselected). Key verified results:
- Selection rule matches hodge1_s3.vertex_coupling exactly
- D(4)^2 factorization verified to 1e-60
- D_even(4) PSLQ: [-24, 0, 12, -1, -96, 96] = pi^2/2 - pi^4/24 - 4G + 4*beta(4)
- D_odd(4) has opposite G and beta(4) coefficients (cancel to 1e-30)
- D_even(4) + D_odd(4) = D(4) to 1e-60 (exact Hurwitz)
- Vertex-restricted sunset converges and increases with n_max
- All forbidden pairs are (0, m) and (n, 0) (41 at n_max=20)
- Flat-space adjacent-pair sum = pi^2/3 - 3 (verified to 1%)
- SO(4) channel count W(n1,n2,q) range {0,1,2} verified for n1,n2 < 8, q < 10
- Channel count symmetry W(n1,n2,q) = W(n2,n1,q) verified
- Pair weight formula W_total = 2*min-1-delta verified for n1,n2 <= 12
- N_q = min(n1,n2) verified for n1,n2 <= 9
- S_min NOT identified in 18-element weight-8 MZV+beta PSLQ basis (100 digits)

## Bugs Fixed from Previous Agent

1. Test `test_dipole_selection_rule`: contradictory assertions on line 58/60 
   (asserted both True and False for same call). Removed the erroneous line.

2. Test `test_only_00_is_forbidden`: incorrectly claimed only (0,0) is forbidden.
   Actually ALL (0,m) and (n,0) pairs are forbidden (triangle + parity makes
   q = |n1-n2| the only candidate, and 2*max(n1,n2) is always even).

3. Test `test_correction_is_rational`: expected only (0,0) contribution but
   the correction sums all 2*n_max+1 forbidden pairs.

4. Even/odd split computed by direct truncated sum (N=2000) compared against 
   exact Hurwitz with tolerance 1e-30 -- impossible due to O(1/N) truncation 
   error. Fixed by implementing exact Hurwitz evaluation for D_even and D_odd
   via `_dirac_D_even` and `_dirac_D_odd`.

5. PSLQ tests used basis {1, pi^2, pi^4, zeta(3)} which is WRONG for D_even(4).
   The correct basis is {1, pi^2, pi^4, G, beta(4)}. All structural theorem 
   tests rewritten with correct basis.
