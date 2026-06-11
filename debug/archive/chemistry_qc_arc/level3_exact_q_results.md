# Level 3 Exact Q-Matrix Results

**Date:** 2026-03-28 (Track B, Q-matrix improvement)

---

## 1. Summary

Implemented the exact Q-matrix for the coupled-channel solver by computing dP/dR algebraically from Hellmann-Feynman quantities. The exact Q = P*P + dP/dR replaces the closure approximation Q ~ P*P.

**Key finding:** The exact Q consistently improves over the closure at l_max >= 1, reducing error by 14-19%. At l_max=0 the improvement is marginal because the overcorrection is dominated by the diagonal Q (DBOC), which is identical in both formulations (dP_mm/dR = 0 by antisymmetry).

---

## 2. Derivation

The stored P convention is:
    P_stored[mu, nu] = V_ad[nu, mu] / (mu_mu - mu_nu)

where V_ad = evecs^T @ V_coupling @ evecs is the coupling in the adiabatic basis.

Differentiating:
    dP/dR = [dV_ad[nu,mu]/dR * Delta - V_ad[nu,mu] * dDelta/dR] / Delta^2

where Delta = mu_mu - mu_nu and:
- dDelta/dR = V_ad[mu,mu] - V_ad[nu,nu]  (Hellmann-Feynman)
- dV_ad[nu,mu]/dR = sum_{k!=nu} P_s[nu,k] * V_ad[k,mu] + sum_{k!=mu} P_s[mu,k] * V_ad[nu,k]

All quantities are available from the existing Hellmann-Feynman computation. No finite differences needed. Verified against central FD (delta=1e-7) to < 1e-9 accuracy.

Key properties:
- dP_mm/dR = 0 identically (P_mm = 0 for all R by antisymmetry)
- dP/dR is antisymmetric: dP_mn/dR = -dP_nm/dR
- The correction only affects off-diagonal Q elements
- Diagonal Q (DBOC) is exact in both closure and exact formulations

---

## 3. Results: He Ground State

### l_max comparison (Z=2, n_basis=15, n_channels=3, N_R_radial=2000)

| l_max | q_mode='full' (closure) | q_mode='exact' | Improvement | Single-channel |
|:-----:|:-----------------------:|:--------------:|:-----------:|:--------------:|
| 0     | 1.123%                  | 1.104%         | 1.7%        | 0.164%         |
| 1     | 0.363%                  | 0.310%         | 14.6%       | 0.559%         |
| 2     | 0.293%                  | 0.238%         | 18.8%       | 0.633%         |
| 3     | 0.272%                  | 0.219%         | 19.5%       | 0.651%         |

**Exact He:** -2.903724 Ha

### All q_modes at l_max=0

| q_mode    | Energy (Ha) | Error (%) | Notes |
|:----------|:-----------:|:---------:|:------|
| none      | -2.9064     | 0.094%    | P-coupling only, overshoots below exact |
| diagonal  | -2.8714     | 1.114%    | DBOC on diagonal, overcorrects |
| full      | -2.8711     | 1.124%    | Closure Q = P*P, overcorrects |
| exact     | -2.8717     | 1.104%    | Q = P*P + dP/dR, slight improvement |

### dP/dR magnitude (R=2.0, l_max=0)

- dP/dR_01 peak: 0.143 (comparable to P_01 = -0.332)
- Q_exact_01 peak: 0.157
- Q_closure_01 peak: 0.017
- The dP/dR correction is ~8x larger than the off-diagonal closure Q

---

## 4. Analysis

### Why l_max=0 improvement is small

At l_max=0, the overcorrection is 1.1% (vs 0.16% adiabatic). The DBOC adds +0.035 Ha to the potential, while the P-coupling cancels ~97% of this. The residual overcorrection is almost entirely from the diagonal Q (DBOC), which is exact in both formulations.

The dP/dR correction modifies only the off-diagonal Q. At l_max=0 with 3 channels, the off-diagonal Q is a small perturbation compared to the diagonal DBOC. The exact Q moves the energy by ~0.0006 Ha toward exact, but this is <2% of the total overcorrection.

### Why higher l_max shows more improvement

At l_max >= 1, the additional angular channels create more off-diagonal coupling. The dP/dR correction grows with the number of channels and with the coupling strength between channels of different l. The 14-19% improvement at l_max=1-3 reflects the increasing importance of off-diagonal non-adiabatic coupling in the multichannel problem.

### Convergence direction

With q_mode='exact', the coupled-channel error decreases monotonically from l_max=1 to l_max=3: 0.310% -> 0.238% -> 0.219%. This confirms that the exact Q does not introduce any spurious behavior.

---

## 5. Implications

1. **Exact Q is correct but not transformative at l_max=0:** The 1.1% overcorrection at l_max=0 is intrinsic to including DBOC without the full P-coupling cancellation. The coupled-channel formalism captures the P-coupling, but the balance between diagonal DBOC and off-diagonal P-coupling is imperfect at small channel count.

2. **Useful at higher l_max:** The 19% improvement at l_max=3 (0.27% -> 0.22%) is worthwhile. For production calculations at l_max >= 2, q_mode='exact' is recommended.

3. **Path to sub-0.1%:** The remaining 0.22% error at l_max=3 is likely from basis truncation (n_channels=3) and the finite angular basis (n_basis=15). Increasing n_channels and n_basis should converge toward the 0.05% FD result.

4. **Backlog item resolved:** The exact dP/dR derivation is complete and algebraic. No finite differences needed. The implementation adds ~20% computational overhead (an additional O(N^3) matrix-vector product per R point) which is negligible compared to the radial solve.

---

## 6. Implementation Details

### New code
- `_compute_dPdR()` function in `geovac/algebraic_coupled_channel.py`
- `compute_exact_dPdR` parameter on `compute_algebraic_coupling()`
- `q_mode='exact'` option on `solve_hyperspherical_algebraic_coupled()`

### New tests (5 added to test_algebraic_coupled_channel.py)
| # | Test | Status | What it verifies |
|:-:|:-----|:------:|:-----------------|
| 11 | test_dPdR_matches_finite_difference | PASS | Algebraic dP/dR vs FD to < 1e-6 |
| 12 | test_dPdR_diagonal_zero | PASS | dP_mm/dR = 0 (antisymmetry) |
| 13 | test_q_exact_diagonal_equals_closure | PASS | Q_exact_mm = Q_closure_mm |
| 14 | test_exact_lmax_convergence_monotonic | PASS | Error decreases l_max=1,2,3 |
| 15 | test_exact_improves_over_full | PASS | Exact < full at l_max=2 |

All 15 tests pass (10 existing + 5 new).
