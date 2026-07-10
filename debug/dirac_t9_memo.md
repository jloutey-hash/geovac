# Track T9 Memo: D^2 on S^3 -- Paper 18 Empty-Cell Probe

**Track:** T9 (Dirac-on-S^3 Tier 3 sprint)
**Date:** 2026-04-15
**Status:** COMPLETE -- empty cell filled (degenerate with scalar calibration)

---

## 1. Lichnerowicz Verification

The squared Dirac operator on unit S^3 satisfies the Lichnerowicz formula:

    D^2 = nabla*nabla + R/4

where R = 6 is the scalar curvature of the unit 3-sphere. Since R/4 = 3/2:

    D^2 eigenvalues: mu_n = (n + 3/2)^2 = (2n+3)^2 / 4,  n = 0, 1, 2, ...
    nabla*nabla eigenvalues: mu_n - 3/2 = n^2 + 3n + 3/4

with degeneracies g_n = 2(n+1)(n+2) (full Dirac, both chiralities).

Verification: (n + 3/2)^2 - 3/2 = n^2 + 3n + 9/4 - 3/2 = n^2 + 3n + 3/4. PASS.

---

## 2. Spectral Zeta Closed Forms

### Partial-fraction decomposition

Substituting m = n+1 and using the identity m(m+1) = [(2m+1)^2 - 1]/4:

    zeta_{D^2}(s) = 2^{2s-1} * [lambda(2s-2) - lambda(2s)]

where lambda(a) = (1 - 2^{-a}) * zeta_R(a) is the Dirichlet lambda function.

Since zeta_R(2k) = rational * pi^{2k} (Bernoulli number formula), every term is a rational multiple of pi^{even}. **This is the key structural result.**

### Closed forms at s = 1, 2, 3, 4

| s | zeta_{D^2}(s) | Numerical | Convergent? |
|---|---|---|---|
| 1 | -pi^2/4 | -2.4674 | No (analytic continuation) |
| 2 | pi^2 - pi^4/12 | 1.7522 | Yes (first convergent) |
| 3 | pi^4/3 - pi^6/30 | 0.4234 | Yes |
| 4 | pi^6(168 - 17pi^2)/1260 | 0.1654 | Yes |

**Pattern:** zeta_{D^2}(s) = c_{2s-2} * pi^{2s-2} + c_{2s} * pi^{2s} at every integer s, where c_j are rationals. This is a two-term polynomial in pi^2.

### Structural mechanism

D^2 eigenvalues are (2n+3)^2/4 -- perfect squares of odd half-integers. The spectral zeta sums therefore reduce to

    Sum_{k=3,5,7,...} f(k) * k^{-2s}

which are Dirichlet series over odd integers with exponent 2s (always even). By the Euler product / Bernoulli formula, sums of (odd integer)^{-2s} are always rational multiples of zeta_R(2s) = rational * pi^{2s}. **No odd-zeta content can appear.**

This is in structural contrast with the first-order Dirac spectral zeta:

    zeta_D(s) = Sum g_n * |lambda_n|^{-s} = Sum g_n * (2n+3)^{-s} / 2^s

At odd s, this involves Sum (odd)^{-s} = lambda(s), and when s is odd, lambda(s) contains zeta_R(s) which is an independent transcendental (not a power of pi). This is how Tier 1 D3 found zeta(3): zeta_D(4) involves zeta_R(3) through the degeneracy-weighted sum.

**Squaring the operator converts odd^{-s} to odd^{-2s}, eliminating all odd-zeta content.** This is a structural theorem, not an accident of small s.

---

## 3. Hopf-Equivariant Decomposition

The Hopf U(1) decomposition of the Dirac spectrum has an important structural feature: all charges are half-integer. At level n (CH convention):

    n=0: charges {-1/2, +1/2}, mult 2 each
    n=1: charges {-3/2, -1/2, +1/2, +3/2}, mult (2,4,4,2)
    n=2: charges {-5/2, ..., +5/2}, mult (2,4,6,6,4,2)

Integer charges (including q=0) have zero multiplicity at every level. This reflects the spinor character of the Dirac eigenspaces -- they live on a shifted lattice compared to the scalar (Fock) harmonics whose charges are integers.

Since D^2 shares eigenspaces with |D|, the per-charge spectral zetas zeta_{D^2, q}(s) inherit the same half-integer charge structure. Each per-charge sum still involves only (odd)^{-2s} denominators, so **no per-charge decomposition can produce odd-zeta content either.** The structural argument is universal across all Hopf sectors.

---

## 4. Paper 18 Empty-Cell Verdict

### The 2x2 grid

Paper 18 Section IV organizes transcendental content by operator order and bundle type:

| | Scalar bundle | Spinor bundle |
|---|---|---|
| **1st-order** | N/A (Laplacian is 2nd-order) | zeta(odd) + alpha^2 [Tier 1 D3, Tier 2 T5] |
| **2nd-order** | pi^{even} [Paper 24 calibration] | **pi^{even} [THIS RESULT]** |

### Verdict: FILLED, DEGENERATE WITH SCALAR CALIBRATION

The (2nd-order, spinor-bundle) cell produces pi^{even} -- the same transcendental class as the (2nd-order, scalar) cell. The Lichnerowicz constant shift R/4 = 3/2 modifies rational coefficients but does not introduce new transcendental content.

**The primary discriminant for transcendental content is the operator order, not the bundle type:**
- 1st-order operators (Dirac D) produce zeta(odd) content through sums over odd^{-s} at odd s
- 2nd-order operators (Laplacian, D^2) produce pi^{even} content through sums over integer^{-2s} or odd^{-2s} at integer s

The bundle (scalar vs spinor) affects:
- The degeneracy weights (n+1)^2 vs 2(n+1)(n+2)
- The eigenvalue lattice (integer vs half-integer)
- The Hopf charge lattice (integer vs half-integer)
- The rational coefficients in pi^{2k} expansions

But NOT the transcendental class.

---

## 5. Paper 24 Cross-Reference

Paper 24's Coulomb/HO asymmetry thesis distinguishes:
- 1st-order complex-analytic operators (Bargmann-Segal/Euler) --> linear spectra --> no pi
- 2nd-order Riemannian operators (Laplace-Beltrami) --> quadratic spectra --> calibration pi

T9 extends this to spinors: D^2 is a 2nd-order operator on spinor sections and produces pi^{even}, consistent with the thesis. The fact that D^2 = (1st-order)^2 does not reintroduce the 1st-order (odd-zeta) content -- squaring maps odd denominators to even powers, which is precisely the mechanism that produces pi^{even} from Bernoulli numbers.

This provides the structural explanation for why the HO rigidity theorem (Paper 24 Theorem 3) and the calibration pi phenomenon are operator-order effects: the difference between SUM_n a_n * n^{-s} (odd-zeta possible at odd s) and SUM_n a_n * n^{-2s} (pi^{even} only at integer s) is precisely the 1st-vs-2nd-order distinction.

---

## 6. Recommendation for Papers

**Paper 18 Section IV:** Add a paragraph documenting that the (2nd-order, spinor-bundle) cell is degenerate with the scalar calibration cell. The structural statement: "Operator order is the primary discriminant for transcendental content on round S^3. The Lichnerowicz identity D^2 = nabla*nabla + R/4 shifts eigenvalues by a rational constant, preserving the transcendental class. Per-charge Hopf decomposition confirms this universally."

**Paper 24 Section V:** Add a corollary to the HO rigidity theorem: "The 1st-order/2nd-order operator-order distinction is bundle-independent. D^2 on spinor sections produces pi^{even} at every integer s, degenerate with the scalar Laplacian. The mechanism is algebraic: squaring maps Dirichlet lambda(odd) to lambda(even), and lambda(2k) = rational * pi^{2k} by the Bernoulli formula."

---

## Files

- `debug/tier3_t9_squared_dirac.py` -- computation script
- `debug/data/tier3_t9_squared_dirac.json` -- results
- `debug/dirac_t9_memo.md` -- this memo
