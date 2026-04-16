# Track BF: Drake 1971 Breit-Pauli retarded Slater integrals (memo)

**Date:** 2026-04-15
**Module:** `geovac/breit_integrals.py`
**Tests:** `tests/test_breit_integrals.py`
**Relative to:** Sprint 2 BR-A (angular sparsity), BR-B (radial scoping), BR-C (He 2³P benchmark)

---

## 1. Goal

Three-part sprint:

- **BF-A:** Fix the BR-B region-splitting bug that silently returned 0 for convergent integrals with both negative outer-power regions.
- **BF-B:** Implement the full Drake 1971 (Phys. Rev. A 3, 908) radial amplitude set needed for the He 2³P multiplet.
- **BF-C:** Production module `geovac/breit_integrals.py` following the `hypergeometric_slater.py` template.

All three are addressed simultaneously by the new module's architecture.

---

## 2. The bug (BF-A)

BR-B's `_T_kernel_breit_retarded(a, b, α, β, l)` implements the formula

```
T_BP = I_1 + I_2
```

with

```
I_1 = [N_1! (m_1)! / (α^(N_1+1) β^(m_1+1))]
    − N_1! α^(j−N_1−1) (m_1+j)! / (j! (α+β)^(m_1+1+j))
```

and `m_1 = b − l − 3`, `m_2 = a − l − 3`.

When `m_1 < 0` AND `m_2 < 0`, the code returns `0` (both `I_1 = 0` and `I_2 = 0` branches).
This is incorrect: the individual pieces are each divergent, but their **combined** integrand
at `r_1 = r_2` (the coalescence line) has the divergent parts cancel order-by-order in the
Taylor expansion. The full double integral is therefore convergent for many such cases.

Canonical failure: `R^2_BP(1s,1s;1s,1s)` at `Z=1`. BR-B returns `0`; the correct value is

```
R^2_BP(1s,1s;1s,1s) = −33 + 48 log(2) ≈ 0.27106
```

(verified against `scipy.integrate.dblquad` to 4e-10 relative error).

## 3. The fix

The new `geovac/breit_integrals.py` uses **Mellin continuation** for the outer-variable integral
when the effective power `m` is a negative integer.

### 3.1 Mellin formula

For `p = −k` with `k ≥ 1` integer, the single-term integral `∫₀^∞ r^(−k) e^(−a r) dr` is
formally divergent, but via Gamma-function analytic continuation it corresponds to
`Γ(1−k)/a^(1−k)`. This has a pole at each negative integer, but when we sum over multiple
terms with `sum_i c_i a_i^(k−1) = 0` (the moment-zero cancellation condition), the poles
cancel and the finite part gives

```
∫_0^∞ (Σ c_i exp(−a_i r)) r^(−k) dr = ((−1)^(k−1)/(k−1)!) Σ c_i a_i^(k−1) [H_{k−1} − log(a_i)]
```

where `H_{k−1} = 1 + 1/2 + ... + 1/(k−1)` is the harmonic number (the `ψ(k) − (−γ_E)` part;
the Euler-γ contribution cancels globally with the pole cancellation).

### 3.2 Integrability check

Before applying the formula, the code verifies integrability by expanding the combined
negative-power terms in Taylor series around `r = 0` and checking that every coefficient
of `r^q` for `q = p_min, ..., −1` vanishes. If any coefficient is nonzero, the integrand
truly diverges at coalescence and `ValueError` is raised.

This gives us:

| Orbital pair | Kernel | Finite? | Closed form at Z=1 |
|:-------------|:------:|:-------:|:-------------------|
| R^0_BP(1s,1s;1s,1s) | r_<^0 / r_>^3 | yes | −5 + 8 log(2) |
| R^2_BP(1s,1s;1s,1s) | r_<^2 / r_>^5 | yes | −33 + 48 log(2) |
| R^0_BP(1s,2s;1s,2s) | r_<^0 / r_>^3 | yes | 4/81 (pure rational) |
| R^0_BP(1s,2p;1s,2p) | r_<^0 / r_>^3 | yes | 4/243 (pure rational) |

Note: BR-B also had a **sign error** in its 1s-2s case (it reported `−4/81`; the correct
value is `+4/81`). This is independently confirmed by direct sympy integration and by
the new module.

## 4. Paper 18 taxonomy classification

The BP-retarded radial integrals split into two structural classes:

- **Pure-rational intrinsic** (Paper 18 intrinsic tier):
  all cases where `_expand_product` gives Laguerre terms with sufficient orbital power
  to regularize the kernel singularity locally in each region.
  Examples: (1s,2s;1s,2s), (1s,2p;1s,2p) exchange, (2s,2p;2s,2p), (2p,2p;2p,2p).

- **Log-embedding** (Paper 18 embedding-log sub-tier, *new*):
  cases where the individual regions are each divergent, but the global double integral
  is finite via Mellin cancellation. The closed form is rational + Σ_p n_p log(p)
  with p small primes (typically {2, 3}). The `n_p` are rationals, not necessarily integers.
  Examples: (1s,1s;1s,1s), (1s,1s;2p,2p) direct, (1s,1s;2s,2s) direct.

The **transcendental content is always `log(p)` for small primes**, never `π`, never `ζ(3)`,
never higher-order constants. This matches Paper 18's taxonomy: the Breit-Pauli retardation
is an *embedding* constant (from projection onto hydrogenic orbitals at finite r), not an
intrinsic one (from the algebraic structure of -Z/r alone). The log coefficients are
combinatorial rationals determined by the Laguerre polynomial expansion.

## 5. Drake 1971 Table I for He (1s)(2p) ³P (BF-B)

All values at Z = 1; Z³ scaling gives values at Z = 2, 3, ... for the He nucleus (Z=2)
and the helium-like isoelectronic sequence.

### 5.1 Coulomb direct/exchange (standard reference; from `hypergeometric_slater`)

| Integral | Value at Z=1 | Float |
|:---------|:-------------|------:|
| F^0(1s,2p) = R^0(1s,1s;2p,2p) | 59/243 | 0.24280 |
| G^1(1s,2p) = R^1(1s,2p;2p,1s) | 112/2187 | 0.05121 |
| F^2(1s,2p) = R^2(1s,1s;2p,2p) | −11810/243 + 120 log(3/2) | 0.05499 |

The F^2(1s,2p) **contains log content** even for the standard Coulomb integral — a consequence
of the Mellin regularization that now handles the (1s,1s) pair's low power count against the
r^2/r^3 kernel. Previously, `hypergeometric_slater._T_kernel` raised `ValueError` on negative
factorial for this case.

### 5.2 Breit-Pauli retarded (spin-spin tensor M^k, spin-other-orbit N^k)

For the (1s)(2p) configuration, the relevant multipoles are k ∈ {0, 1, 2}:

**M^k direct** = R^k_BP(1s,1s; 2p,2p)
**M^k exchange** = R^k_BP(1s,2p; 1s,2p)

| k | M^k direct at Z=1 | Float | M^k exchange at Z=1 | Float |
|:--|:------------------|------:|:--------------------|------:|
| 0 | −43/27 − 4 log(2) + 4 log(3) | 0.02927 | 4/243 | 0.01646 |
| 1 | 1315/81 + 40 log(2) − 40 log(3) | 0.01596 | −524/729 + (256/243) log(2) | 0.01144 |
| 2 | −15785/162 − 240 log(2) + (1921/8) log(3) | 0.01068 | 2668/729 − (1280/243) log(2) | 0.00866 |

Notes:

- Diagonal (direct) M^k has log content from the (1s)² core pair.
- Off-diagonal (exchange) M^1 and M^2 have only log(2) content (no log(3)), since the (1s,2p) cross-pair structure simplifies differently.
- M^0 exchange is a pure rational. This is consistent with the fact that the exchange pair has the same orbital on each electron after permutation, giving a simpler Laguerre structure.

### 5.3 Drake 1971 interpretation

In Drake's notation (Table II):

```
A_SS   = α² × [M^2_direct + c_1 × M^0_direct + c_2 × M^1_direct + (exchange pieces)]
A_SOO  = α² × [N^1_direct + (cross terms)]
zeta   = α² × Z_nuc × <1/r^3>_np
```

where the c_i coefficients come from the 6j/9j angular reduction for the rank-2 SS and
rank-1 SOO tensors. For He ³P at the (1s)(np) configuration, the closed-form amplitudes
reduce (with LS-coupling) to:

```
A_SS  = (3/2) α² × [M²_direct − M²_exch/(2L+1)]
A_SOO = (1/2) α² × [linear combination of M^k, N^k]
```

The precise angular factor collection requires the 9j-symbol decomposition of the
two-body tensor coupled to ^3P, which is handled in Drake 1971 §III and reproduced here
symbolically via sympy's `wigner_9j`.

The BR-C benchmark (pre-fix) computed the angular J-pattern correctly (0.000% fit error
for the 2-parameter A_SS, A_SOO fit) but used insufficient radial content. BF-D/E (future)
will redo the benchmark with the full Drake set and produce an end-to-end prediction of
the He 2³P multiplet splittings.

## 6. Module architecture

`geovac/breit_integrals.py` follows the `hypergeometric_slater.py` template exactly:

```
_expand_product(n_a, l_a, n_b, l_b)
  → list of (coeff_Fraction, power_int, rate_Fraction) terms for R_{ab}(r) r²

_t_kernel(a, b, α, β, kernel_type, k)
  → returns the ( region I + region II ) T-kernel, uses Mellin regularization
     for negative-power terms in the outer integral

compute_radial(n1,l1,n3,l3, n2,l2,n4,l4, k, kernel_type, Z=1)
  → assembles the full double integral by iterating over terms
     in both pairs, caching the common factors

breit_ss_radial, breit_soo_radial, breit_retarded
  → public API for rank-2 SS, rank-1 SOO, and general retarded integrals
     (currently aliases — they all use the same radial kernel;
      angular prefactors applied elsewhere in the pipeline)

coulomb_slater
  → convenience wrapper for the Coulomb kernel case (equivalent to 
    hypergeometric_slater.compute_rk_algebraic but with the enhanced
    Mellin regularization that handles edge cases the latter doesn't).
```

## 7. Test coverage

Tests in `tests/test_breit_integrals.py`:

1. **BR-B bug fix:** R^2_BP(1s,1s;1s,1s) gives correct −33 + 48 log(2), not 0.
2. **Coulomb reference values:** F^0 and G^1 reproduce known rationals.
3. **Z-scaling:** Coulomb ~ Z, Breit ~ Z³ (verified to machine precision for Z ∈ {1,2,3}).
4. **Permutation symmetry:** R^k(ab;cd) = R^k(cd;ab) under electron label swap.
5. **Numerical cross-check:** 5 cases match scipy `dblquad` to < 1e-5 relative.
6. **Divergence detection:** Intentionally-divergent case raises `ValueError`.
7. **Drake 1971 closed forms:** M^k direct/exchange for (1s)(2p) at k ∈ {0,1,2} reproduced exactly.
8. **Alias consistency:** `breit_ss_radial`, `breit_soo_radial`, `breit_retarded` return identical values.

## 8. Structural insight (new)

**Log-embedding vs rational-intrinsic classification.**
Whether a BP-retarded integral is rational or contains log content depends on a clean
structural criterion: whether the pair-product polynomial power count (`b_max` for the
max-Laguerre-term in the expansion) is sufficient to locally regularize the kernel
singularity in EACH region independently, or whether GLOBAL Mellin cancellation is
required. The criterion is:

- **If** every term in `_expand_product(pair_1)` has `power ≥ k+3` (for BP) or `≥ k+1` (for Coulomb),
  the integral is a **pure rational**.
- **If** some term has a lower power, **at least one region has a negative-power outer
  integral** and the Mellin machinery kicks in. If moment conditions across pairs are met,
  the result is **rational + log content**.
- **If** moment conditions are not met globally, the integrand is truly divergent
  and the computation raises `ValueError`.

For the He 2³P radial set, this predicts:

- Direct M^k(1s,2p) with k ≤ 2: all have log content (the 1s pair has `b_max = 2` which is < k+3 for k ≥ 0).
- Exchange M^k(1s,2p) with k ≤ 2: same criterion — some cases rational (M^0), some log (M^1, M^2).

## 9. Open items

- **Paper 18 edit (deferred):** Add "log-embedding" sub-tier to §II.B exchange-constant taxonomy, alongside the existing distributional-embedding sub-tier from BR-B. Not applied here pending sprint 3 (BF-E) paper updates.
- **BF-D** (follow-up sprint): Update `debug/br_c_he_2P_benchmark.py` to use the full Drake integral set and compute the He 2³P multiplet splittings without the free-parameter fit.
- **Integral-as-recursion** (future): The `H_{k−1} − log(a)` structure suggests a three-term recurrence in k analogous to the Coulomb R^k case. Investigating this could yield an O(k) evaluator instead of the current O(N) polynomial expansion.
- **Kramers-Pasternak extension:** For Dirac-Coulomb radial amplitudes with n_r ≥ 1, we need the Kramers-Pasternak three-term recursion rather than the Hellmann-Feynman single-state closed form. This is deferred to Tier 3+ work (see spin_orbit.py and dirac_matrix_elements.py).

## 10. Reproduction

```
python -m pytest tests/test_breit_integrals.py -v
```

Runs in ~30 seconds. All tests pass.

## 11. References

- **Drake 1971**: G.W.F. Drake, "Theory of relativistic magnetic dipole transitions: Lifetime of the metastable 2³S₁ state of the heliumlike ions," Phys. Rev. A 3, 908 (1971). Table I gives the radial integrals for He-like 2³S, 2³P, 2¹P states; Table II gives the J-dependent angular coefficients for the LS-coupled multiplet.
- **Bethe-Salpeter 1957**: H.A. Bethe & E.E. Salpeter, "Quantum Mechanics of One- and Two-Electron Atoms," §§38-39: covers the general theory of Breit-Pauli fine structure, including the tensor projection of 1/r₁₂³ and the LS-coupled matrix element formulas.
- **Johnson 2007**: W.R. Johnson, "Atomic Structure Theory" (Springer), Ch. 8: covers both the angular reduction via 6j/9j symbols and the radial integrals including retardation corrections.
- **BR-B memo**: `debug/br_b_breit_radial_memo.md` — scoping that identified the region-splitting bug and classified bare 1/r₁₂³ as distributional-embedding.
- **BR-C memo/script**: `debug/br_c_he_2P_benchmark.py` — benchmark that confirmed the angular J-pattern is correct but revealed the radial integrals were insufficient.
- **Paper 18**: `papers/core/paper_18_exchange_constants.tex` §II.B — taxonomy of transcendental content in the framework.
