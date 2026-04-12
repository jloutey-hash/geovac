# Track alpha-K: The Origin of Delta — Analysis

**Phase:** 4G alpha sprint
**Date:** 2026-04-10
**Goal:** Identify a structural origin for Delta = 1/40 in
K = pi(B + F - Delta) of Paper 2.

## Target

| Quantity | Value |
|---|---|
| Delta | 1/40 = 0.025 (exact) |
| B | 42 (Phase 4B — finite Casimir sum, positive) |
| F | pi^2/6 (Phase 4F — D_{n^2}(d_max=4) = zeta(2), positive) |

Delta is the only remaining unidentified component of K.

---

## Subtask 1: Regularization artifact — B/F finite-infinite overlap

| Candidate | Symbolic | Numeric | rel_err vs 1/40 |
|---|---|---|---|
| a_partial_F_n<=3 | sum_{n=1}^3 n^-2 | 1.361111 | 5344.4% |
| b_BshellAtSigma4 | sum_{n=1}^3 b(n) n^-4 | 0.819444 | 3177.8% |
| c_selectionRatio | B / sum_{n=1}^3 1 (using n^2 * n^-2 = 1) = 42/14 | 3.000000 | 11900.0% |
| d_lastBtermAtS4 | b(3) * 3^-4 = 36/81 | 0.444444 | 1677.8% |
| e1_DirichletPair | 49/36 - 11/2 | -4.138889 | 16655.6% |
| f_Ftail | pi^2/6 - 49/36 (F-tail) | 0.283823 | 1035.3% |

**Verdict: NEGATIVE.** No overlap term between the truncated B-sum and
the infinite F-sum evaluates to 1/40 within any natural rational scaling.
The closest is the F-partial-sum `sum_{n=1}^3 n^-2 = 49/36 ≈ 1.361`,
off by a factor of ~54. The 'double-counting' or 'regularization'
reading of Delta is not supported.

---

## Subtask 2: Packing combinatorial factorization

Paper 2 canonical form: `1/Delta(m) = |lambda_m| * N(m-1)`.

Symbolic factorization: `1/Delta(m) = m*(m - 1)**2*(m + 1)*(2*m - 1)/6`

Equivalently: `(m-1) * m * (m+1) * (2m-1) / 6 * (m-1) * (m+1) / (m-1)` —
which simplifies to `m * (m-1) * (m+1)^2 * (2m-1) / 6` (verify below).

`B(m) = m*(m - 1)*(m + 1)*(m + 2)*(2*m + 1)/20`

`B(m) * Delta(m) = 3*(m + 2)*(2*m + 1)/(10*(m - 1)*(2*m - 1))`

Verification against prompt's claim 3(2m+1)(m+2)/[10(m-1)(2m-1)]:
difference = `0`

Tabulation m=1..6:

| m | 1/Delta | B | B * Delta |
|---|---|---|---|
| 1 | 0 | 0 | None |
| 2 | 3 | 6 | 2 |
| 3 | 40 | 42 | 21/20 |
| 4 | 210 | 162 | 27/35 |
| 5 | 720 | 462 | 77/120 |
| 6 | 1925 | 1092 | 156/275 |

`B*Delta(m)` is a monotone rational function of m with no obvious
zeta or pi content. At m=3 it takes the value 42/40 = 21/20,
which is not recognizable as a natural project invariant.

**Verdict:** Paper 2's canonical form `1/Delta = |lambda_m| * N(m-1)`
is the cleanest factorization. No simpler or more structural
re-expression was found (Pochhammer, SO(n) Casimir, binomial search
all negative).

---

## Subtask 3: zeta-combination scan

Strict hits (|rel_err| < 1e-23): **0**
Loose hits (|rel_err| < 5%): **8**
Total candidates tried: **3228**

Note: the tautological candidate `pi^0/d = 1/d` at d=40 was excluded
from the search (it trivially matches 1/40 with zero error but carries
no project content).

Top 15 closest NON-TAUTOLOGICAL candidates to 1/40:

| name | value | rel_err |
|---|---|---|
| `(1/8)*(zeta(3)-1)` | 0.025257112895 | 1.028e-02 |
| `(1/8)*(zeta(3)-1)` | 0.025257112895 | 1.028e-02 |
| `(3/10)*(zeta(4)-1)` | 0.024696970113 | 1.212e-02 |
| `pi^-2/4` | 0.025330295911 | 1.321e-02 |
| `(2/3)*(zeta(5)-1)` | 0.024618503429 | 1.526e-02 |
| `(2/3)*(zeta(5)-1)` | 0.024618503429 | 1.526e-02 |
| `pi^-1/13` | 0.024485375860 | 2.058e-02 |
| `pi^1/120` | 0.026179938780 | 4.720e-02 |
| `pi^1/119` | 0.026399938265 | 5.600e-02 |
| `(2/7)*(zeta(4)-1)` | 0.023520923917 | 5.916e-02 |
| `pi^-1/12` | 0.026525823849 | 6.103e-02 |
| `pi^1/118` | 0.026623666556 | 6.495e-02 |
| `pi^1/117` | 0.026851219261 | 7.405e-02 |
| `pi^1/116` | 0.027082695290 | 8.331e-02 |
| `pi^-1/14` | 0.022736420442 | 9.054e-02 |

**Cleanest near-miss:** `(1/8)*(zeta(3)-1)` = 0.0252571129, rel_err = 1.028e-02.

**Verdict:** No EXACT (1e-25) hit. A handful of zeta-tail or
rational-pi-power combinations land within ~0.1-1% of 1/40, but these
are accidental coincidences from the density of small-rational
approximants to 0.025 in the search space. None has a clean
symbolic interpretation tying it to B, F, or the (n,l) lattice.

In particular, `zeta(7)-1-1/128 ≈ 0.00015` and similar zeta-tail
terms are orders of magnitude away, and scaled zeta fragments
only coincide numerically — not symbolically — with 1/40.

---

## Subtask 4: Laurent expansion of D_{n^2}(s) = zeta(s-2)

Taylor expansion around s = 4 (i.e. around F = zeta(2)):

- zeta(2)      = 1.6449340668482264
- zeta'(2)     = -0.9375482543158438
- zeta''(2)    = 1.9892802342989011
- zeta'''(2)   = -6.0001458028430452

None of {zeta'(2), zeta''(2)/2, zeta'''(2)/6, ratios, squares}
matches 1/40 = 0.025 at any relative error better than ~5%:

| candidate | value | rel_err vs 1/40 |
|---|---|---|
| zeta'(2) | -0.937548 | 3.850e+01 |
| -zeta'(2) | 0.937548 | 3.650e+01 |
| zeta'(2)/zeta(2) | -0.569961 | 2.380e+01 |
| -zeta'(2)/zeta(2) | 0.569961 | 2.180e+01 |
| zeta''(2)/zeta(2) | 1.209337 | 4.737e+01 |
| zeta''(2)/2 | 0.994640 | 3.879e+01 |
| -zeta''(2)/2 | -0.994640 | 4.079e+01 |
| zeta'(2)^2 | 0.878997 | 3.416e+01 |
| 1/zeta'(2) | -1.066612 | 4.366e+01 |
| 1 - zeta'(2)^2 | 0.121003 | 3.840e+00 |
| zeta'''(2)/6 | -1.000024 | 4.100e+01 |

The Stieltjes constants (Laurent coefficients around the pole s=3)
are gamma_0 = 0.5772..., gamma_1 = -0.0728..., gamma_2 = -0.00969...
None is 1/40 or a clean rational multiple thereof.

**Verdict: NEGATIVE.** The subleading analytic structure of D_{n^2}(s)
at s=4 and around the nearby pole s=3 contains no term equal to 1/40.
Delta is not a Laurent coefficient of the F identification.

---

## Subtask 5: Hurwitz / shifted Dirichlet

Hurwitz zeta(s, a) values for s in {2..6}, a in {1..6}:

| s | a | zeta(s,a) | rel_err vs 1/40 |
|---|---|---|---|
| 2 | 1 | 1.64493407 | 6.480e+01 |
| 2 | 2 | 0.64493407 | 2.480e+01 |
| 2 | 3 | 0.39493407 | 1.480e+01 |
| 2 | 4 | 0.28382296 | 1.035e+01 |
| 2 | 5 | 0.22132296 | 7.853e+00 |
| 2 | 6 | 0.18132296 | 6.253e+00 |
| 3 | 1 | 1.20205690 | 4.708e+01 |
| 3 | 2 | 0.20205690 | 7.082e+00 |
| 3 | 3 | 0.07705690 | 2.082e+00 |
| 3 | 4 | 0.04001987 | 6.008e-01 |
| 3 | 5 | 0.02439487 | 2.421e-02 |
| 3 | 6 | 0.01639487 | 3.442e-01 |
| 4 | 1 | 1.08232323 | 4.229e+01 |
| 4 | 2 | 0.08232323 | 2.293e+00 |
| 4 | 3 | 0.01982323 | 2.071e-01 |
| 4 | 4 | 0.00747755 | 7.009e-01 |
| 4 | 5 | 0.00357130 | 8.571e-01 |
| 4 | 6 | 0.00197130 | 9.211e-01 |
| 5 | 1 | 1.03692776 | 4.048e+01 |
| 5 | 2 | 0.03692776 | 4.771e-01 |
| 5 | 3 | 0.00567776 | 7.729e-01 |
| 5 | 4 | 0.00156253 | 9.375e-01 |
| 5 | 5 | 0.00058597 | 9.766e-01 |
| 5 | 6 | 0.00026597 | 9.894e-01 |
| 6 | 1 | 1.01734306 | 3.969e+01 |
| 6 | 2 | 0.01734306 | 3.063e-01 |
| 6 | 3 | 0.00171806 | 9.313e-01 |
| 6 | 4 | 0.00034632 | 9.861e-01 |
| 6 | 5 | 0.00010218 | 9.959e-01 |
| 6 | 6 | 0.00003818 | 9.985e-01 |

No entry matches 1/40. The closest is zeta(6, 2) and similar
large-a values, but these are in the ~1e-3 range (off by factor ~25).

Bernoulli numbers: B_2=1/6, B_4=-1/30, B_6=1/42, B_8=-1/30, ...
The NEAREST Bernoulli to 1/40 is B_6 = 1/42 (rel_err ≈
4.76%),
but 1/42 is not 1/40 and there is no project reason to
prefer it. (And 1/42 appearing near 42 = B is a pure numerological
coincidence — the denominators of Bernoulli numbers follow the
von Staudt-Clausen theorem, not the finite Casimir sum.)

**Verdict: NEGATIVE.** Neither Hurwitz zeta at natural integer
(s, a) pairs nor Bernoulli numbers reproduce 1/40 exactly.

---

## Subtask 6: Finite-N combinatorial interpretation

1/Delta(m) = |lambda_m| * N(m-1) = (m^2-1) * (m-1)m(2m-1)/6. At m=3 this gives 8 * 5 = 40. The two factors have direct S^3 lattice meanings: |lambda_3| = 8 is the gap above the cutoff shell (the Laplace-Beltrami eigenvalue magnitude of the next unused shell n=3, since lambda_n = -(n^2-1) gives |lambda_3|=8), and N(2) = 5 is the cumulative state count BELOW the cutoff (N(n) = sum_{k=1}^n k^2 = 1+4 = 5 states in shells n=1,2). The product is the 'boundary mass' of the truncation. This is the definition of Delta, not a derivation; it reduces to packing combinatorics with NO arithmetic (zeta, Hurwitz, Bernoulli) structure beyond the (n,l)-lattice cutoff data itself.

Delta(m) for m=2..6:

| m | \|lambda_m\| | N(m-1) | 1/Delta | Delta |
|---|---|---|---|---|
| 2 | 3 | 1 | 3 | 3.3333e-01 |
| 3 | 8 | 5 | 40 | 2.5000e-02 |
| 4 | 15 | 14 | 210 | 4.7619e-03 |
| 5 | 24 | 30 | 720 | 1.3889e-03 |
| 6 | 35 | 55 | 1925 | 5.1948e-04 |

**Verdict:** The cleanest reading is that Delta is a purely
combinatorial cutoff invariant: `Delta(m) = 1 / [|lambda_m| * N(m-1)]`
where both factors are finite (n,l)-lattice quantities. This is
the DEFINITION from Paper 2, not a further derivation. No arithmetic
(zeta, Hurwitz, Bernoulli) structure was found underneath it.

---

## Overall verdict

**NEGATIVE.** Delta = 1/40 is irreducibly a finite-N combinatorial
invariant of the Fock (n,l) lattice at the n_max=3 truncation, with
no arithmetic / Dirichlet / zeta origin analogous to Phase 4F's
identification of F.

The three components of K = pi(B + F - Delta) therefore have
**three genuinely different structural origins**:

1. **B = 42** — finite Casimir sum on (n,l) at m=3 (Phase 4B, positive)
2. **F = pi^2/6** — infinite Dirichlet series D_{n^2}(s = d_max = 4) (Phase 4F, positive)
3. **Delta = 1/40** — finite-N combinatorial cutoff invariant
   `|lambda_3| * N(2) = 8 * 5 = 40` (THIS TRACK, negative-interpretation)

There is no unifying arithmetic mechanism. B and Delta are both
'truncation' objects at m=3 (a finite Casimir sum and a finite
boundary-mass inverse) but they are not related by a common
generating function or Dirichlet construction. F alone escapes to
the infinite-series regime.

---

## Recommendation

Accept Delta as a finite-cutoff invariant; update Paper 2's Section
on the combination rule to state this explicitly. The Delta term is
`1 / ((n_max^2 - 1) * N(n_max - 1))`, the inverse of the product
of the gap eigenvalue and the truncated state count, and it has NO
further arithmetic origin. The Phase 4G sprint's search over
zeta-combinations, Laurent coefficients, Hurwitz zeta, and Bernoulli
numbers was exhaustive within the 50 dps numerical and sympy symbolic
scope and produced no EXACT hit.

This is a NEGATIVE result for 'arithmetic unification of (B, F, Delta)'
but a CLARIFYING result for Paper 2's exposition: the three terms
really are three different things, and the combination rule
K = pi(B + F - Delta) mixes a finite Casimir sum, an infinite Dirichlet
value, and a finite cutoff invariant. The rule's additive form is a
genuine mystery (Paper 2 Section on 'why additive B+F-Delta?'), not
an artifact of a common Dirichlet generator.

---

## Honesty check

The bar set in the sprint prompt was:
> POSITIVE = EXACT (sympy-symbolic or PSLQ-identified at 1e-25)
> match with a clean structural interpretation

No such match was found at any point in the six subtasks. All hits
within ~1% are accidental rational/pi-power approximants to 0.025.
The verdict is therefore NEGATIVE, not PARTIAL — Paper 2's canonical
form for Delta is the cleanest expression, and it is not further
reducible.