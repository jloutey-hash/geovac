# Order-one defect of the true CH window as an edge-crossing sum

Date: 2026-10-31. Driver: `debug/order_one_edge_crossing.py`. Data:
`debug/data/order_one_edge_crossing_1_2_3.json`, `debug/data/order_one_edge_crossing_4.json`,
run logs `debug/data/order_one_edge_crossing_run*.txt`. Nothing under papers/, tests/, geovac/,
docs/ was modified.

**Verdict: GO.** The order-one defect is an exact edge-crossing sum, verified elementwise for
every pair at n_max = 1..4, with two controls. The register's phrasing was right in substance:
the slope weight f'(lambda_i) is replaced by the Dirac gap. The gauge companion turns out to be
the same object, not just the same mechanism. It is the f'-weighted half-trace of the order-one
defect at a = b = chi.

## 1. The quantity measured

`tests/test_true_window_real_structure.py::_stats(n, 1)` computes, for every pair of compressed
band-limited multipliers,

    C(a,b) = [[D_W, T_a], B_b],   T_a = P a P,  D_W = P D P,  B_b = J T_b^* J^{-1},

and the same file asserts B_b = T_b (to 1e-12). So C(a,b) = [[PDP, PaP], PbP]. The test reports
max-abs in the ambient label basis (V C V^dag) and the operator norm.

## 2. Derivation

Write Q = 1 - P. P is a spectral projection of D, so [PDP, PaP] = P[D,a]P. Expand
P[D,a]P b P = P[D,a] b P - P[D,a] Q b P, and the same for the other ordering:

    C(a,b) = P [[D,a], b] P  -  P [D,a] Q b P  +  P b Q [D,a] P.

In the continuum [D, a] = c(da) is Clifford multiplication, and it commutes with the function b.
So [[D,a],b] = 0 and the first term drops. In the Dirac eigenbasis ([D,a]_ik = (lambda_i -
lambda_k) a_ik), for i, l in W:

    C_il = sum_{k not in W} [ -(lambda_i - lambda_k) a_ik b_kl  +  b_ik (lambda_k - lambda_l) a_kl ].

This is a pure edge-crossing sum. Every term is a path W -> exterior -> W: one step out across
the edge and one step back, through a single exterior state k. No interior term survives. The
weight is the Dirac gap across the crossed edge, (lambda_i - lambda_k) on the outgoing leg of
the first term and (lambda_k - lambda_l) on the returning leg of the second. No spectral-function
slope appears, because no f is involved.

Structural corollaries:
- The same sum taken over interior k gives -C, and the sum over all k gives 0. The defect is
  exactly the part of the vanishing continuum row sum that the window cuts off, which is the
  gauge mechanism verbatim.
- A band-N multiplier moves the Dirac level by at most N - 1. Then a_ik with k exterior needs i
  within N_a - 1 of the top level, which is the "Toeplitz boundary term" vanishing that
  `test_residual_is_a_boundary_term` pins. It is now a consequence of the formula.

**Relation to the gauge companion.** Take a = b = chi = i(b - b^dag), which is Hermitian. Then

    C(chi,chi)_ii = -2 sum_{k not in W} (lambda_i - lambda_k) |chi_ik|^2,

and so

    k_pg = -sum_{i in W, k not in W} |chi_ik|^2 f'(lambda_i)(lambda_i - lambda_k)
         = (1/2) Tr_W [ f'(D_W) [[D_W, chi_W], chi_W] ].

The gauge identity [chi,[chi,D]] = 0 is the order-one identity [[D,a],b] = 0 at a = b = chi. The
gauge breaking is the order-one defect of chi with itself, traced against the heat slope. The
crossing matrix elements chi_ik (i in W, k not in W) are the same in both sums. Only the outer
weight differs: the order-one operator carries the bare gap, and the heat coefficient multiplies
the gap by f'(lambda_i). One side effect is that the half-trace form computes k_pg from window
data alone, without exterior states.

## 3. Numerical verification

Ambient: j <= (3 n_max - 1)/2. Every state that a band <= 2 n_max multiplier reaches from W is
present, so the crossing sum is complete, not truncated. The Clebsch-Gordan coefficients are an
in-driver Racah formula, checked against sympy to 1.1e-16. The Dirac eigenbasis is built per
j-block.

| n_max | dim W | dim Q | pairs | max rel. (own basis) | max rel. vs TEST's V C V^dag (elementwise) | test max-abs / op-norm / #fail reproduced |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| 1 | 4 | 24 | 25 | 3.0e-16 | 5.9e-16 | 3/4, 0.75, 4 |
| 2 | 16 | 166 | 900 | 4.6e-15 | 7.4e-15 | 5/3, 1.8634, 513 |
| 3 | 40 | 530 | 8281 | 1.3e-14 | 1.7e-14 | 21/8, 2.9020, 5856 |
| 4 | 80 | 1220 | 41616 | 2.3e-14 | (not run; test table stops at 3) | — |

The "vs TEST" column maps the crossing formula onto the test's own ambient labels and compares
every entry of every pair against `_setup`/`_stats` from the test module, with C built from the
test's Bs. Absolute residuals are at most 2.7e-15, against defect entries of order 1 to 4.

Gauge identity, (1/2) Tr_W f'(D_W) C(chi,chi) against `window_k` of the gauge test, at fixed tau
in {1, 4} and n_max in {2, 3, 4, 5}: max relative difference 1.9e-15. Example values: k_pg =
-2.3544 (tau = 1, n = 2) and -0.9324 (tau = 4, n = 5).

## 4. Controls (non-tautology)

The formula is an algebraic consequence of [D,P] = 0 and one analytic input, [[D,a],b] = 0.
The controls check that the input is what does the work.

1. **The continuum identity holds for D and fails for |D|.** ([[D,a],b])_WW is computed with
   all ambient intermediates and maximised over all pairs. With D it is 4.6e-16, 1.4e-15,
   1.9e-15 and 2.7e-15 (n = 1..4), against gap scale 1.5 to 4.5. With |D| in place of D (the
   sign(D) confusion) it is 0.75, 1.385 and 2.07 (n = 1..3).
2. **The crossing formula fails for |D|.** Repeating the comparison with |D|, the crossing
   formula misses the directly computed |D|-defect by 0.75, 1.385 and 2.07 (n = 1..3). At
   n = 1 the |D|-defect is exactly 0 while the formula gives 0.75. So matching is not automatic.
   It requires the Clifford identity, and therefore the sign structure of D.
3. **Crossing restriction.** Over interior k the same sum equals -C, and over all k it equals 0.
   Dropping the restriction therefore gives 0, not C.
4. **Gauge half-trace with |D|.** Replacing D by |D| inside the double commutator gives
   (1/2) Tr f' C = 0 to 3e-16 at every (tau, n), against |k_pg| between 0.47 and 4.7. The
   identity singles out D.

## 5. Scope and caveats

- This is an operator identity for each pair (a, b), valid for every band-limited pair. The
  test's max-abs and op-norm numbers are norms of this crossing sum. The formula does not by
  itself give the growth law of those norms with n_max. That needs the size of the exterior
  matrix elements a_ik near the edge, times gaps of order n_max.
- Exterior states were enumerated completely, not truncated. The formula needs a_ik only for
  j <= (3 n_max - 1)/2.
- The gauge relation is checked here only for the gauge test's single chi (band 1/2, the YC
  instrument). The derivation is general for any Hermitian band-limited chi and any even or
  odd f, since only [f(D), P] = 0 is used.
- Proposed guard, not written: per the CLAUDE.md §9 guard-writing rule it is a separate step.
  It would assert the elementwise identity at n_max = 2, 3 against `_stats`, and the |D|
  mismatch controls 1 and 2 as the rejected wrong answer.
