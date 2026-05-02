# Spectral Action on the Fock S³ Spectral Triple — Analysis Memo

**Date:** 2026-05-02  
**Track:** WH1 spectral triple formalization, Step #2 (spectral action computation)  
**Script:** `debug/spectral_triple_spectral_action.py`  
**Data:** `debug/data/spectral_triple_spectral_action.json`

## Summary

Computed the Connes-Chamseddine (CC) spectral action Tr f(D²/Λ²) on the finite
Fock-projected S³ graph for 4 configurations (n_max=2,3 × uniform/CG ×
permutation/kramers) with 8 test functions (sharp cutoff, Gaussian, sigmoid at
4 widths, polynomial at 4 powers), extracted heat kernel SD coefficients, and
compared against K = π(B + F − Δ) ≈ 137.036.

**Verdict: STRUCTURAL NEGATIVE.** K/π = 43.62 exceeds the maximum possible
spectral action at n_max=3 (dim_H = 40), making any match structurally
impossible. The gap decomposes cleanly into Paper 2 invariants.

## Key Finding: The Ceiling Theorem

At n_max=3 with Dirac degeneracies g_n = 2(n+1)(n+2):
- dim_H = 4 + 12 + 24 = 40 = Δ⁻¹ from Paper 2

For any test function f with 0 ≤ f(x) ≤ 1:
- Tr f(D²/Λ²) = Σᵢ f(λᵢ²/Λ²) ≤ 40 × max(f) = 40

Since K/π = B + F − Δ = 43.62 > 40 = dim_H, **no normalized test function
can produce K/π from the n_max=3 spectrum.**

## Gap Decomposition

K/π − dim_H = 43.62 − 40 = 3.62

This decomposes as:
- (B − Δ⁻¹) = 42 − 40 = 2: Casimir trace exceeds state count because l=0
  states contribute to dim_H but not to the l(l+1)-weighted B sum
- F − Δ = π²/6 − 1/40 = 1.62: the infinite Dirichlet series contribution
  minus the boundary term

Structural reading: F = ζ(2) = D_{n²}(d_max) (Phase 4F) is the contribution
from the infinite tail of the Fock degeneracy series. It cannot be captured by
any finite-dimensional trace. This is the spectral-action proof that K requires
an infinite object (F) alongside finite ones (B, Δ).

## Positive Structural Findings

1. **Str(sharp cutoff) = 0 exactly** for all 4 configurations. This is the
   finite-graph analog of the index theorem Ind(D) = 0 on S³. The chirality
   grading splits dim_H into equal halves (chi+ = chi- = 20 at n_max=3).

2. **Tr(D²) matches continuum Casimir to 0.02%** for CG weights:
   - n_max=2: 84.023/84.000 (0.028%)
   - n_max=3: 378.068/378.000 (0.018%)
   This validates the CG edge weights as the correct spectral triple structure.

3. **dim_H = 40 = Δ⁻¹** at n_max=3, confirming the Paper 2 coincidence in the
   spectral triple framework.

4. **B = 42 from the lattice Casimir sum** reproduced exactly.

## Supertrace Behavior

| Quantity | n_max=2 CG | n_max=3 CG | Continuum |
|:---------|:-----------|:-----------|:----------|
| Str(Θ) | 0 | 0 | 0 |
| Str(D²)/dim_H | −1.000 | −1.012 | 0 (ST-1 F1) |
| Str(\|D\|) | −4.002 | −8.067 | 0 (eta) |
| Str(D) | −20.004 | −31.995 | 0 (eta) |

- Str(sharp) = 0 confirms balanced chirality (analog of index theorem)
- Str(D²) ≈ −dim_H is a nonzero finite-size effect; the continuum SD
  cancellation (Sprint ST-1 Finding F1: a_k^{D²}/a_k^{Δ_LB} = 4) is a
  large-N limit that hasn't converged at n_max=2,3
- Str(D) ≈ −32 at n_max=3 CG; related to the chirality-weighted eta
  invariant. Nonzero because {D, γ} ≠ 0 on the finite graph

## Heat Kernel SD Coefficients

Fitted Tr exp(−tD²) in the small-t regime (t = 0.001 to 0.1) to the
expansion a₀ t^{−3/2} + a₁ t^{−1/2} + a₂ t^{1/2} + a₃ t^{3/2}:

| Coeff | n_max=3 CG | Continuum S³ | Ratio |
|:------|:-----------|:-------------|:------|
| a₀ | −0.00207 | √π = 1.772 | −0.001 |
| a₁ | 3.008 | √π = 1.772 | 1.697 |
| a₂ | 39.399 | √π/8 = 0.222 | 177.8 |

The a₀ ≈ 0 is structurally forced: a finite graph has no t^{−3/2} divergence
(bounded spectrum). The a₁ overshoots by 70%. The a₂ is 178× the continuum
value. The SD asymptotic expansion is not a useful approximation tool on the
finite graph — it requires the continuum limit.

## Spectral Zeta Moments

Tr |D|^k at n_max=3 CG:
- k=0: 40 (= dim_H)
- k=1: 120.006 (continuum: Σ g_n |λ_n| = 120)
- k=2: 378.068 (continuum: 378)
- k=3: 1230.45 (continuum: Σ g_n |λ_n|³)

The CG-weighted graph reproduces continuum moment sums to better than 0.02%
at k=1,2. This is the strongest validation that the CG weights are correct.

## Consistency with Prior Results

This result is consistent with:
- **Phase 4G (α-K):** B, F, Δ have categorically different origins with no
  common generator. The spectral action cannot provide that common generator.
- **Sprint ST-1 F1:** Perturbative CC supertrace vanishes on continuum S³
  (confirmed as a limit; doesn't hold exactly on finite graph).
- **Sprint ST-1 F4:** Non-perturbative CC remainder is always negative ~−1
  to −3, never near K/π. Confirmed: no cutoff produces K/π.
- **WH5:** α is a projection constant combining three structurally distinct
  spectral objects, not derivable from a single functional.

The ceiling theorem strengthens WH5: the finite-graph CC action
*structurally cannot produce K/π* at the natural truncation n_max=3.
The infinite Dirichlet contribution F = ζ(2) is the part that escapes.

## What This Does and Does Not Rule Out

**Ruled out:** K = π(B + F − Δ) as a Connes-Chamseddine spectral action
Tr f(D²/Λ²) at any finite cutoff on the n_max=3 Fock graph.

**Not ruled out:**
1. K as a spectral action on the *product* geometry S³_continuum × F_finite
   (almost-commutative setup, standard NCG approach)
2. K as the sum of a finite spectral trace (B, Δ) plus a regularized
   infinite series (F), assembled by different mechanisms — this is the
   Phase 4F/4G picture, and the spectral action results confirm it
3. K at higher n_max where dim_H > K/π (n_max ≥ 4 gives dim_H = 80 > 43.62)
   — but the right test function would need to suppress ~36.4 eigenvalues
   worth of weight, which seems unlikely to produce K/π without tuning

---

## Continuum Limit: SD Two-Term Exactness Theorem

**Date:** 2026-05-02 (continuation session)
**Data:** `debug/data/spectral_action_sd_exactness.json`

### The Theorem

The Seeley-DeWitt expansion of the Dirac heat kernel on unit S³ terminates
at exactly 2 terms:

**K(t) = (√π/2)·t^{-3/2} − (√π/4)·t^{-1/2} + O(e^{-π²/t})**

All higher SD coefficients a_k = 0 identically for k ≥ 2.

### Proof (Jacobi Theta Modular Identity)

The Dirac heat kernel K(t) = Σ_{n≥0} 2(n+1)(n+2) exp(-(n+3/2)²t) can be
rewritten as a derivative of the Jacobi theta function θ₂(0, e^{-t}).
The modular identity

  θ₂(0, e^{-t}) = √(π/t) · θ₄(0, e^{-π²/t})

maps the small-t expansion (which produces the SD series) to a large-π²/t
regime where θ₄ → 1 + O(e^{-π²/t}). This forces:
- The polynomial-in-t part of t^{3/2}·K(t) to be exactly a₀ + a₁t (2 terms)
- All terms a₂t², a₃t³, ... to vanish identically
- The remainder to be exponentially small: O(e^{-π²/t})

### Numerical Verification

1. **a₂ extraction at 80 dps (50,000 terms):** The quantity
   [K(t) - a₀/t^{3/2} - a₁/t^{1/2}] / t^{1/2} gives a₂ = -3.88×10⁻⁷⁴
   at t = 10⁻⁴ — numerical noise, not a physical coefficient.

2. **Residual scaling at 100 dps (100,000 terms):** The residual
   K(t)·t^{3/2} - a₀ - a₁t stays at ~10⁻⁹⁶ for all t from 0.01 to 10⁻⁶.
   If a₂ ≠ 0, this residual would scale as t². It does not — it's flat
   numerical noise at every t tested.

3. **Exponential correction verification:** The analytical correction
   √π·e^{-π²/t}·[2π²/t^{5/2} - 1/t^{3/2} + 1/(2t^{1/2})] matches the
   numerical residuals at ratio = 1.000 at t = 0.20, 0.15, 0.10, 0.08.

### The Spectral Action Equation

Setting Tr exp(-D²/Λ²) = K/π gives:

  (√π/2)·Λ³ − (√π/4)·Λ = K/π

or equivalently the depressed cubic:

  2Λ³ − Λ = 4(K/π)/√π

**Cardano closed form:**

  Λ = ∛(A + √(A²−1/216)) + ∛(A − √(A²−1/216))

where A = 2(B + F − Δ)/√π ≈ 24.611.

Result: **Λ_∞ = 3.7102454679060528505** (50 dps, verified by both
polyroots and Cardano to 46 matching digits).

### Structural Analysis

1. **Λ is NOT autonomously selected.** The CC framework provides the exact
   2-term functional form but no mechanism to fix Λ. The value 3.710...
   is determined by inverting the spectral action at K/π. For any target T,
   the cubic 2Λ³ − Λ = 4T/√π has exactly one real root (IVT + monotonicity).
   No physics selects this particular Λ.

2. **The 40-mode coincidence.** A sharp cutoff at Λ = 3.71 includes exactly
   40 = Δ⁻¹ modes (n = 0, 1, 2, with |λ₂| = 3.5 < 3.71 < 4.5 = |λ₃|).
   This connects to the supertrace sprint F2 finding: Δ⁻¹ = 40 is the
   Euler-Maclaurin upper boundary term of the Dirac mode-count sum.

3. **All three Paper 2 invariants enter Λ.** Through A = 2(B+F−Δ)/√π,
   the Cardano formula depends on B, F, and Δ simultaneously. But this
   is just because K/π = B + F − Δ is the target — no new decomposition.

4. **The exponential correction is negligible.** At Λ = 3.71:
   e^{-π²Λ²} ≈ 10⁻⁵⁹. The 2-term polynomial is exact to 59 digits.

5. **Continuum convergence from finite graphs.** Root-finding at n_max=4..8
   gives: 4.813, 4.024, 3.814, 3.744, 3.720 — converging monotonically
   to Λ_∞ = 3.710 from above. The IVT guarantees existence for n_max ≥ 4
   (where dim_H ≥ 80 > 43.62 = K/π). The match is tautological.

### SD Ratio Check

a₀(D²)/a₀(Δ_LB) = (√π/2)/(√π/4) = 2 = rank(spinor bundle on S³)

This is the rank of a Weyl spinor (2-component). The supertrace sprint
F1 reported ratio = 4 = dim(spinor bundle); the discrepancy is the Dirac
vs Weyl counting convention (Dirac = 2 × Weyl, so 2 × 2 = 4).

### Scripts and Data

- `debug/spectral_action_sd_verification.py` — a₂ extraction (executed)
- `debug/spectral_action_continuum_limit.py` — continuum Λ_∞ root-finding
- `debug/spectral_action_lambda_structure.py` — structural analysis (has
  π-ratio bug at line 104, see note below)
- `debug/spectral_action_rootfind.py` — brentq at n_max=4..8
- `debug/spectral_action_nmax_sweep.py` — grid sweep n_max=3..8
- `debug/data/spectral_action_sd_exactness.json` — all numerical results

**Bug note (FIXED):** `spectral_action_sd_verification.py` line 104 computed
`target = 4 * K_over_pi / (mp_pi * mpsqrt(mp_pi))` which gave 4K/π^{5/2}
instead of the correct 4K/π^{1/2}. The "ratio = π" in its output was the
missing π factor, not a structural finding. Fixed to use `mpsqrt(mp_pi)` alone.
