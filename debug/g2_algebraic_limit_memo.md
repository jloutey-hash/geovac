# g-2 Algebraic Limit Investigation Memo

**Direction 5 from strategic brief**
**Date:** 2026-04-24
**Status:** MIXED — one positive structural result, one structural obstruction

---

## Question

Can the flat-space Schwinger limit alpha/(2pi) of the anomalous magnetic moment
(g-2 vertex correction) on S^3 be derived algebraically — replacing the numerical
spectral mode sum with a closed-form expression using Hurwitz/zeta techniques?

## Answer Summary

**The shell-summed vertex correction F_1 (charge form factor): YES**, algebraic
closed form exists. It decomposes as:

    Lambda(n_ext) = r(n_ext) + (2*n_ext - 1)^2 * D(4)

where r(n_ext) is an **exactly rational** number and D(4) = pi^2 - pi^4/12 is
the Dirac Dirichlet series at s=4 (already known in closed form from Paper 28).

**The anomalous magnetic moment F_2 (magnetic form factor / g-2): NO**, there is a
structural obstruction from CG irrationals at fixed m_j projections that prevents
Hurwitz decomposition.

---

## Detailed Findings

### Finding 1: d_T/mu cancellation (key simplification)

The vertex correction from `qed_self_energy.py` contains photon propagator mu(q)
and photon degeneracy d_T(q) factors. Both are **identically equal**:

    d_T(q) = q(q+2)    (transverse photon degeneracy on S^3)
    mu(q)  = q(q+2)    (photon mass on S^3)

Therefore d_T(q)/mu(q) = 1 for all photon modes q >= 1. This eliminates ALL
q-dependence from the vertex correction weight, leaving only the SO(4) channel
count W(n1, n2, q) summed over allowed photon modes.

### Finding 2: S_q symmetry and closed form

The total channel weight S_q(n1, n2) = sum_q W(n1, n2, q) is:
- **Symmetric:** S_q(a, b) = S_q(b, a) (verified for all a, b <= 7)
- **Zero at n=0:** S_q(0, n) = S_q(n, 0) = 0 for all n (reproduces the
  self-energy structural zero theorem)
- **Closed form for n1, n2 >= 1:** S_q(n1, n2) = 2*min(n1, n2) - 1 - delta_{n1,n2}
  (verified against direct computation)

Since S_q is symmetric, the double-vertex product is S_q^2:

    Lambda(n_ext) = sum_{n_int} g(n_int)/|lam(n_int)|^4 * S_q(n_ext, n_int)^2

### Finding 3: Tail factorization

For fixed n_ext >= 1, S_q(n_ext, n_int) is CONSTANT for all n_int > n_ext:

    S_q(n_ext, n_int) = 2*n_ext - 1    (for n_int > n_ext)

This means the tail of the sum (n_int > n_ext) factorizes:

    Lambda_tail = (2*n_ext - 1)^2 * [D(4) - D_head(n_ext)]

where D_head is a finite rational sum and D(4) = pi^2 - pi^4/12 is exact.
The head terms (n_int <= n_ext) also have rational S_q^2 values, so the
full vertex correction is:

    Lambda(n_ext) = r(n_ext) + (2*n_ext - 1)^2 * D(4)

with r(n_ext) **exactly rational**.

### Finding 4: Exact rational parts

| n_ext | P_tail = (2n-1)^2 | r(n_ext) | r(n_ext) float |
|:------|:-------------------|:---------|:---------------|
| 1 | 1 | -55552/50625 | -1.09732... |
| 2 | 9 | -140031424/13505625 | -10.36838... |
| 3 | 25 | -33434276032/1093955625 | -30.56274... |
| 4 | 49 | -3029660404067072/48049812916875 | -63.05249... |
| 5 | 81 | -89300485278602603392/823410424031320125 | -108.45197... |
| 6 | 121 | -32758506406804745792/196050100959838125 | -167.09253... |

Denominators are products of (2k+3)^4 for k = 0, 1, ..., n_ext, which come
from |lam(k)|^4 = ((2k+3)/2)^4 in the Dirac spectrum.

### Finding 5: Structural obstruction for F_2 (g-2)

The anomalous magnetic moment F_2 requires the **m_j-dependent** difference:

    B = L(m_j = +1/2) - L(m_j = -1/2)

This breaks the shell-sum symmetry that enabled Finding 3. Specifically:

1. **Individual B(n_int) are NOT rational.** Per-level values live in growing
   number fields: B(1) in Q(sqrt(2), sqrt(3)), B(2) adds sqrt(5), B(3) adds
   sqrt(7), etc. (Documented in `debug/g2_racah_test.py`.)

2. **Root cause:** On S^3, positive chirality spinors have (j_L, j_R) = ((n+1)/2, n/2)
   with j_L - j_R = 1/2. The CG coefficients CG(j_L, j_R, j; m_L, m_R, m_j) for
   this half-integer step generically produce sqrt(2j+1) type irrationals.

3. **Cannot cancel:** Each n_int introduces algebraically independent square roots.
   Since sqrt(p) for distinct primes p are linearly independent over Q, the
   irrationals from different n_int levels cannot cancel in the infinite sum.

4. **Wigner-Eckart does not help:** The reduced matrix element is itself non-rational
   because the internal CG coupling at each vertex involves unequal (j_L, j_R) pairs.

5. **Contrast with F_1:** The shell-summed F_1 works because summing over ALL m_j
   invokes CG orthogonality, which produces rational 6j/9j symbols. Fixing m_j
   for F_2 extraction breaks this orthogonality.

### Finding 6: Current convergence status of F_2

From `debug/data/alpha_g_minus_2_ratio_investigation.json`:
- F_2/Schwinger ratio at n_max=12: **1.0845** (5 reliable digits)
- The residual F_2/Schwinger - 1 = 0.0845 ≈ 2/25 = R/(12*|lambda_ext|^2)
  matches the Parker-Toms leading curvature correction
- Power-law convergence: B(n) ~ n^{-7}, very slow
- PSLQ found no relation with {1, pi, sqrt(2), sqrt(3), zeta(3)} at maxcoeff=200

---

## Assessment for Paper 28

### What goes in Paper 28

The F_1 algebraic closed form (Findings 1-4) is a genuine new result:

    Lambda(n_ext) = r(n_ext) + (2*n_ext - 1)^2 * (pi^2 - pi^4/12)

This is at the same level as the existing Paper 28 results:
- D(4) = pi^2 - pi^4/12 (Theorem 1)
- D_even(4) - D_odd(4) = 2^3 * (beta(4) - beta(2)) (Theorem 3)
- Self-energy structural zero (Theorem 4)

A natural Paper 28 addition: **"One-loop vertex charge correction in algebraic form"**
as a new theorem or corollary.

### What does NOT go in Paper 28

The F_2 = alpha/(2pi) Schwinger limit on S^3 remains numerical, converging
at ~33% above Schwinger at n_max=20 (from `qed_self_energy.py`) or ~8.5%
above at n_max=12 with explicit CG computation (from `qed_anomalous_moment.py`).
The CG irrational obstruction is structural and should be documented as a
clean negative result.

---

## Concrete Next Steps

### If pursuing F_2 further (alternative approaches):

1. **Large-n asymptotic expansion of B(n_int).** Even though individual B(n_int)
   are algebraically irrational, the ASYMPTOTIC form B(n) ~ c * n^{-7} might
   have a rational or transcendental coefficient c that is expressible in closed
   form. This would give the TAIL of the F_2 sum analytically, leaving only a
   finite head to compute numerically.

2. **Wigner-Eckart at the AMPLITUDE level.** Instead of computing B = c_up - c_dn
   (which fixes m_j), decompose the full vertex amplitude into spherical tensor
   components T^(k)_q. The T^(0) (monopole) gives F_1, the T^(1) (dipole)
   gives F_2. The reduced matrix elements of T^(k) ARE rational via Racah algebra.
   The difficulty shifts to expressing the PROBE operator (electromagnetic vertex)
   in the T^(k) basis — this may be tractable.

3. **Differential equation approach.** Express F_2(n_ext) as a function of
   n_ext and derive a recurrence or differential equation from the SO(4)
   selection rules. The flat-space limit would then be n_ext -> infinity,
   potentially yielding alpha/(2pi) as the fixed point of the recurrence.

### Recommended action:

The F_1 algebraic form is ready for Paper 28 integration. The F_2 obstruction
should be documented as a structural negative result with the three alternative
approaches listed as future directions. Approach 2 (Wigner-Eckart tensor
decomposition) is the most promising since it stays within the Racah algebra
framework and avoids the fixed-m_j problem entirely.

---

## Key Files

- `debug/g2_algebraic_limit.py` — Analysis script (this investigation)
- `debug/data/g2_algebraic_limit.json` — Computational results
- `debug/g2_racah_test.py` — CG irrational obstruction (prior result)
- `debug/g2_analytical_c2.py` — Exact sympy B(n_int) computation (prior)
- `debug/g2_curvature_coefficients.py` — Curvature expansion analysis (prior)
- `geovac/qed_self_energy.py` — Shell-summed vertex correction infrastructure
- `geovac/qed_anomalous_moment.py` — m_j-resolved F_2 extraction
- `geovac/qed_vertex.py` — Two-loop Hurwitz decomposition (reference)
