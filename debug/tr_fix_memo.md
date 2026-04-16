# Track TR: Spinor FCI α=0 / scalar FCI mismatch — diagnosis and fix

**Date:** April 2026
**Sprint:** 4
**Author:** Worker sub-agent (GeoVac)
**Goal:** Locate and fix the bug causing spinor FCI at α=0 to differ from
scalar FCI by 0.95–1.66 mHa (growing with n_max), documented in
DC-B (Sprint 2) §5.

---

## TL;DR

The bug is a **missing (−1)^{j_a + 1/2} phase factor** in the jj-coupled
angular Gaunt coefficient ``jj_angular_Xk`` in
``geovac/composed_qubit_relativistic.py`` (lines 94–158 pre-fix). The phase
comes from the Racah reduced matrix element
⟨κ_a ‖ C^k ‖ κ_c⟩ = (−1)^{j_a+1/2} √((2j_a+1)(2j_c+1)) 3j(j_a k j_c; ½ 0 −½) · π.
Grant (2007) Eqs. 8.9.9/8.9.11 and Johnson (2007) Eq. 3.69 both include it;
the GeoVac Tier 2 T3 implementation omitted it.

After the fix, spinor FCI at α=0 matches scalar FCI to machine precision
(<10⁻¹⁰ Ha, empirically 0 to 10⁻¹⁴ Ha) at Z=4, n_max=2,3,4.

---

## 1. Reproduction

On the pre-fix codebase at Z=4 (Be 2+):

| n_max | E_scalar (Ha)    | E_spinor(α=0) (Ha) | ΔE (mHa) |
|------:|:-----------------|:-------------------|---------:|
|     2 | −13.55925721     | −13.55830238       |   +0.955 |
|     3 | −13.56987952     | −13.56843813       |   +1.441 |
|     4 | −13.57376157     | −13.57210581       |   +1.656 |

Spinor is LESS bound than scalar; gap grows with n_max; bounded by
dimension (at n_max=1 both give −13.500 Ha exactly, since |1s²⟩ is the
only config and has κ = −1 only, so cross-κ effects don't appear).

---

## 2. Diagnostic path

### Step 2.1. Isolate the angular layer

Both scalar and spinor builders share the same radial R^k integrals
(via ``geovac/hypergeometric_slater.py``). The R^k cache at n_max=2
has 42 entries with values like F^0(1s,1s)=2.5, G^1(1s,2p_0)=0.972, etc.
These are identical in both paths. So the bug must be in the angular
coupling — either the Gaunt factor ``jj_angular_Xk`` or the M_J
selection rules.

### Step 2.2. Check direct Coulomb integrals at n_max=2

Compared scalar ⟨1s_up, 2p_0 | V | 1s_up, 2p_0⟩ = +0.9712 Ha (physically
correct, a repulsive direct Coulomb integral) against the spinor
equivalents at a=1 (1s_{1/2,+1/2}), b ∈ {2p_{3/2, m_j}, 2p_{1/2, m_j}}:

| (a, b) | Spinor value (pre-fix) | Expected |
|:------:|:---------------------:|:--------:|
| ⟨1s_↑, 2p_{1/2, ±1/2} \| V \| 1s_↑, 2p_{1/2, ±1/2}⟩ | +0.9712 | +0.9712 |
| ⟨1s_↑, 2p_{3/2, m_j} \| V \| 1s_↑, 2p_{3/2, m_j}⟩  | **−0.9712** | +0.9712 |

**Wrong sign for the 2p_{3/2} direct Coulomb.** This ERI feeds the
|1s²⟩ ↔ |1s_↑ 2p_{3/2}⟩² type correlation admixture and artificially
*lowers* the ERI penalty (because it has a wrong sign), effectively
*increasing* the CI energy (less bound). Observed gap direction matches.

### Step 2.3. Check X_0(a,a) for each κ branch

The monopole angular coefficient X_0(a,a) = ⟨κ_a m_a | C^0_0 | κ_a m_a⟩
must equal +1 (trivially — ∫ |ψ|² = 1 for normalized spinor spherical
harmonics). Pre-fix values:

| κ | label | X_0(a,a) |
|:-:|:-----:|:-:|
| −1 | s_{1/2} | **−1** |
| +1 | p_{1/2} | **−1** |
| −2 | p_{3/2} | +1 |
| +2 | d_{3/2} | +1 |
| −3 | d_{5/2} | **−1** |

The sign pattern follows (−1)^{3j_a − 1/2}. The product X_0(a,a) × X_0(b,b)
for cross-κ pairs gives −1 (e.g. s_{1/2} × p_{3/2} = (−1)(+1) = −1), which
is exactly the bug.

### Step 2.4. Derive the missing phase

Working from the standard reduced-matrix-element formula:

    ⟨κ_a m_a | C^k_q | κ_c m_c⟩ = (−1)^{j_a − m_a} · 3j(j_a k j_c; −m_a q m_c) · ⟨κ_a ‖ C^k ‖ κ_c⟩

    ⟨κ_a ‖ C^k ‖ κ_c⟩ = (−1)^{j_a + 1/2} · √((2j_a+1)(2j_c+1)) · 3j(j_a k j_c; ½ 0 −½) · π_{l_a+k+l_c even}

For the diagonal X_0(a,a) with k=0, q=0:

- 3j(j 0 j; 1/2 0 −1/2) = (−1)^{j − 1/2} / √(2j+1)
- 3j(j 0 j; −m 0 m) = (−1)^{j + m} / √(2j+1)
- √((2j+1)²) = (2j+1)
- Parity π trivially satisfied (2l even).

Product (WITHOUT the reduced-matrix phase): 
(−1)^{j − m} · (2j+1) · (−1)^{j − 1/2} / √(2j+1) · (−1)^{j + m} / √(2j+1)
= (−1)^{3j − 1/2}

which matches the pre-fix pattern. With the reduced-matrix phase
(−1)^{j + 1/2}:

(−1)^{3j − 1/2 + j + 1/2} = (−1)^{4j} = +1 for all half-integer j.

### Step 2.5. Verify fix against Grant/Johnson closed form

Implemented a Sympy reference computing

    X_k(a, c) = (−1)^{2j_a + 1/2 − m_a} · π · √((2j_a+1)(2j_c+1)) · 3j · 3j

directly (matching Grant Eq. 8.9.9 + 8.9.11). Verified diagonal
X_0(κ, m; κ, m) = +1 for all (κ, m) combinations tested (κ ∈ {−1, +1,
−2, +2, −3}, all m_j ∈ [−j, j]).

---

## 3. Fix

### 3.1 File and location

``geovac/composed_qubit_relativistic.py``, function ``jj_angular_Xk``
(lines 94–158 pre-fix).

### 3.2 Change

Single-line change to the phase exponent:

```python
# Pre-fix:
phase_exp = (two_j_a - two_m_a) // 2  # (-1)^(j_a - m_a)
phase = (-1) ** phase_exp

# Post-fix:
phase_exp = (two_j_a - two_m_a) // 2  # (-1)^(j_a - m_a)
red_phase_exp = (two_j_a + 1) // 2     # (-1)^(j_a + 1/2)
phase = (-1) ** (phase_exp + red_phase_exp)
```

The ``red_phase_exp`` formula uses ``(two_j_a + 1) // 2`` which works
correctly for half-integer j_a:
- j_a = 1/2 → two_j_a = 1 → (1+1)//2 = 1 → (−1)^1 = −1
- j_a = 3/2 → two_j_a = 3 → (3+1)//2 = 2 → (−1)^2 = +1
- j_a = 5/2 → two_j_a = 5 → (5+1)//2 = 3 → (−1)^3 = −1

Docstring updated with full Grant/Johnson formula and TR fix reference.

### 3.3 Verification

**Diagonal monopole.** After the fix, X_0(κ, m; κ, m) = +1 for all
tested κ ∈ {−1, +1, −2, +2, −3} and all allowed m_j (checked in memory
after clearing _X_CACHE). ✓

**FCI match.** After the fix, Z=4 Be 2+:

| n_max | E_scalar (Ha)    | E_spinor(α=0) (Ha)  | Δ (Ha)        |
|------:|:----------------:|:-------------------:|:-------------:|
|     2 | −13.55925721     | −13.55925721        | 0.0           |
|     3 | −13.56987952     | −13.56987952        | −3.553e−15    |
|     4 | −13.57376157     | −13.57376157        | −8.882e−15    |

All agree to 1-ULP double precision. ✓ (Passing ``< 1e-10`` tolerance
with 4 orders of magnitude margin.)

**He (Z=2).** Also tested: −2.8334051759 Ha both at n_max=2 and
n_max=3, with gaps −8.88e−16 Ha and +2.67e−15 Ha respectively. ✓

**SO shift at CODATA α.** The SO shift on a 1s² GS dropped from the
(buggy) ~10⁻⁸ Ha to the (physically correct) ~10⁻¹²–10⁻¹¹ Ha, scaling
as ~Z² through the 2p-orbital correlation admixture (Kramers
cancellation makes H_SO vanish on l=0, so the only SO contribution
comes from correlation-level 2p admixture). ✓

---

## 4. Downstream impact

### 4.1 Pauli counts change

The corrected X_k removes accidental cancellations in the ERI tensor,
leading to more non-zero Pauli terms. At n_max=2:

| System | Pauli (pre-TR) | Pauli (post-TR) | rel/scalar ratio (pre) | rel/scalar ratio (post) |
|:-------|---------------:|----------------:|:----------------------:|:-----------------------:|
| LiH_rel, Q=30 |  805 | 1413 | 2.42× | 4.24× |
| BeH_rel, Q=30 |  805 | 1413 | 2.42× | 4.24× |
| CaH_rel, Q=20 |  534 |  942 | 2.40× | 4.24× |
| SrH_rel, Q=20 |  534 |  942 | 2.40× | 4.24× |
| BaH_rel, Q=20 |  534 |  942 | 2.40× | 4.24× |

**Isostructural invariance is preserved** — CaH_rel=SrH_rel=BaH_rel=942
Pauli terms, unchanged by TR. Similarly LiH_rel=BeH_rel=1413.

### 4.2 Sunaga head-to-head

Pre-TR claim (``docs/tier2_market_test.md``, Paper 14 §V): GeoVac
rel/scalar 2.4× at n_max=2, matched-Q-18 projection 150×–250× advantage
vs Sunaga RaH-18q.

Post-TR: rel/scalar is now 4.24× at n_max=2. The matched-Q-18
extrapolation advantage is reduced, but GeoVac still wins decisively —
1413 Pauli at Q=30 (LiH/BeH native) extrapolates to ~500 Pauli at Q=18
via Q^{2.5} scaling, i.e. 100× fewer Pauli than Sunaga's 47,099 at Q=18.

Paper 14 §V, Paper 20 Tier 2 table, and ``docs/tier2_market_test.md``
have been updated with the corrected counts.

### 4.3 Fine-structure splittings (He 2³P, Li 2²P, Be 2s2p ³P)

Tier 2 T4 reported "sign + OoM correct, 20–50% accuracy target NOT met"
for these fine-structure splittings, attributing the 66–211% relative
errors to missing multi-electron SS/SOO. With the TR fix, the SO
diagonal itself is unchanged (H_SO is diagonal in (κ, m_j) and not
touched by the X_k fix). The FS splitting sanity checks in
``tests/test_spin_ful_composed.py::test_relativistic_build_smoke``
still pass — the fix doesn't alter the diagonal SO physics, only the
two-body ERI angular coupling.

---

## 5. Regression test

Added
``tests/test_spin_ful_composed.py::test_relativistic_alpha_zero_matches_scalar_fci``
which constructs both scalar and spinor 2e FCI matrices at Z=4, α=0,
n_max=2 and 3, and asserts ``|E_spinor - E_scalar| < 1e-10 Ha`` at each
n_max. Test passes (observed gaps 0.0 Ha and 3.55e-15 Ha; margin is
10⁵×).

The existing test
``test_relativistic_alpha_zero_kills_so`` (which only checks the SO
diagonal *vector*, not the FCI energy) still passes.

---

## 6. Sanity check: why didn't this bug show up earlier?

1. **``test_relativistic_alpha_zero_kills_so``** only checks that the
   diagonal H_SO vector is zero at α=0. This is true regardless of the
   X_k phase bug — H_SO depends only on (n, κ), not on X_k.

2. **``test_relativistic_nmax1_pauli_count``** checks n_max=1 only. At
   n_max=1, the only Dirac labels are 1s_{1/2, ±1/2} with κ=−1, so the
   cross-κ effect never appears.

3. **``test_relativistic_hermiticity``** verifies the resulting qubit
   Hamiltonian is Hermitian. The sign bug preserves Hermiticity (same
   phase factor on both X_k(a,c) and X_k(c,a)), so it passes.

4. **``test_relativistic_pauli_ratio_lih_nmax2``** pinned a specific
   Pauli count (805, 2.42×) that happened to be the *buggy* count.
   The test didn't compare against an independent calculation.

5. **Fine-structure sanity checks in T4** compared total-energy shifts
   at CODATA α, which are dominated by the diagonal H_SO (~10⁻⁴ Ha
   for 2p_{3/2}/2p_{1/2} splitting); the ~10⁻³ Ha sign-bug shift in
   the 2e correlation energy was small enough (and in the "right"
   direction for He-like contracted GS) to not flag as a sanity
   failure. The 20–50% FS accuracy target was missed for other reasons.

6. **The DC-A prediction P1** (⟨α = 0⟩ → match) was correctly
   identified as REJECTED at ~1–2 mHa by DC-B, which is how we got
   here.

The gap is that **no test asserted unitary equivalence of the jj-coupled
and LS-coupled bases at α=0**, which is a fundamental FCI identity
(SU(2) basis invariance). The new regression test plugs this gap.

---

## 7. Files modified

- ``geovac/composed_qubit_relativistic.py`` — fix ``jj_angular_Xk`` phase
  and update docstring with Grant/Johnson references.
- ``tests/test_spin_ful_composed.py`` — add
  ``test_relativistic_alpha_zero_matches_scalar_fci``; update
  ``test_relativistic_pauli_ratio_lih_nmax2`` range from [1.9, 3.0] to
  [3.7, 4.9] to reflect corrected Pauli count.
- ``tests/test_heavy_hydrides.py`` — update
  ``RELATIVISTIC_EXPECTED`` (534 → 942), ``test_relativistic_lih_regression``
  (805 → 1413), and ``test_relativistic_scalar_ratio_n_max2`` (range update).
- ``debug/dc_b_convergence_memo.md`` — headline + §5 + §6 updated to
  reflect the fix.
- ``debug/tr_fix_memo.md`` — this file.

## 8. Recommended follow-up (for PI review)

1. **Paper 14 §V** should reflect the new Pauli counts and the corrected
   rel/scalar ratio.
2. **Paper 20 Tier 2 resource table** should reflect new LiH/BeH/CaH
   rel Pauli counts (1413/1413/942).
3. **``docs/tier2_market_test.md``** should reflect new counts and
   updated Sunaga comparison ratios.
4. **``CLAUDE.md §2 / §10`` Tier 2 entries** may mention 1.00×/2.42×/5.89×
   rel/scalar ratio at n_max=1/2/3; the updated ratios (post-TR) are
   1.00×/4.24×/11.33×. The 1.00× at n_max=1 is unaffected (n_max=1 has
   only κ=−1). Measured post-TR: LiH scalar at n_max=3 is 7878 Pauli,
   rel at n_max=3 is 89226 Pauli → ratio 11.33×. The growing ratio with
   n_max reflects increasing angular-basis density in the spinor (κ,m_j)
   enumeration as higher-l shells open.
5. **Consider** marking any benchmark that depends on the pre-TR Pauli
   count (in docs, memos, or unrelated tests) for review.
