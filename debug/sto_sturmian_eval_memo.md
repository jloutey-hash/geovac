# STO → shared-k Coulomb-Sturmian evaluator for `noci_engine` — scoping memo (2026-10-09)

Diagnostic only. Nothing outside `debug/` was modified.

**Drivers:** `debug/sto_sturmian_eval_convergence.py` (Q1), `debug/sto_sturmian_eval_headtohead.py` (Q2/Q3),
`debug/sto_sturmian_eval_nah.py` (Q4), `debug/sto_sturmian_eval_fitquality.py` (shape audit).
**Data:** `debug/data/sto_sturmian_eval_{convergence,headtohead,nah,fitquality}.{json,log}`.

## Verdict: BORDERLINE (the Gaussian-fit error is real, but swapping in Sturmians is not how to fix it)

- The Sturmian route is **exact to tolerance at manageable term counts.** The worst case is the NaH shell set, at about 50 terms per STO for a 1e-10 residual norm. In the LiH test the integrals reach machine precision (1e-16) at 32 terms.
- The Gaussian-fit error **is material at the precision Paper 58 prints** (NaH Table II, and its E(Na) provenance number). It does **not** change any qualitative conclusion.
- **The material error is not specific to the Sturmian replacement.** Deeper least-squares Gaussian fits, at 8 or more primitives, remove it completely. In every case that matters here the exact per-shell route dominates the infinite shared-k expansion.
- The three-centre wall is untouched.
- The decision gate's GO conditions are met to the letter. But the payoff the gate assumed (removing the fit error) can be had more cheaply another way. What remains specific to the Sturmian route is **decidability and structure, not accuracy**. That is exactly "exact ≠ accurate" (`memory/polyatomic_state_of_play.md` §1).

## Guardrail (CLAUDE.md §3.5, "Sturmian" + "molecule"): outside scope

I loaded Paper 8, Theorem `thm:structural`. It says that in a molecular Sturmian basis with shared momentum scale p₀ (the single-n D-matrix form), (H, S) carry their R-dependence through one orthogonal SO(4) congruence. The eigenvalues are therefore R-independent. That theorem is about the **encoding (variational) basis.**

Here the variational space is unchanged: the atom-centred STOs with per-shell ζ, which is exactly what `noci_engine` uses now. The shared-k set enters only as a **quadrature-free evaluator** of those STOs' integrals. As the number of terms goes to infinity the evaluator returns the STO integrals exactly (shown below to 1e-16). No eigenproblem is ever posed in the Sturmian basis, so no p₀ congruence exists to cancel. The guardrail does not apply. If someone later used the expansion's span as the basis, it would.

## Q1 — Convergence (exact, no quadrature)

Method: the Avery coefficients c_{n,l}(ρ), extended to r^{N-1}e^{-ζr} by applying (−d/dρ)^{N-1-l}. The Sturmian L² Gram matrix is tridiagonal (G_{n,n+1} = −½√((n−l)(n+l+1)/(n(n+1)))). The Sturmian kinetic matrix is k²δ − (k²/2)G. Together these give the L² and H¹ (kinetic) norms of the truncation residual exactly, as tail quadratic forms.

Checks:
- Gram formula vs exact Fraction integrals: 1e-41.
- Derivative coefficients vs exact projection: 2e-40.
- Completeness (‖Σc_nS_n‖² = ‖f‖²): ≤2e-39.

Everything depends only on ρ = ζ/k, through t = |ρ−1|/(ρ+1). The count is symmetric under ρ ↔ 1/ρ to within ±10%. Roughly, N(ε) ≈ ln ε / ln t plus a small l- and N-dependent offset.

Terms needed for a relative L² residual below 1e-6 / 1e-10 (the kinetic-norm counts are within ±10%):

| ρ (t) | 1s | 2s | 2p | 3s | 3d |
|:--|:--:|:--:|:--:|:--:|:--:|
| 1 (0) | 1 / 1 | 2 / 2 | 1 / 1 | 3 / 3 | 1 / 1 |
| 1.25 or 0.8 (0.11) | 8 / 12 | 9 / 14 | 8 / 13 | 11 / 15 | 9 / 13 |
| 2 or 1/2 (0.33) | 16 / 24 | 18 / 27 | 17 / 26 | 20 / 30 | 19 / 28 |
| 4 or 1/4 (0.60) | 35 / 53 | 39 / 59 | 39 / 58 | 44 / 64 | 42 / 63 |
| 8 or 1/8 (0.78) | 72 / 110 | 82 / 122 | 81 / 121 | 90 / 132 | 90 / 131 |
| 16 or 1/16 (0.88) | 147 / 224 | 167 / 248 | 167 / 247 | 185 / 268 | 185 / 268 |

**Realistic shell sets.** k was optimised on a grid to minimise the worst-case count, at a 1e-10 L² residual.

| set | per-centre optimum (terms per STO) | single global k |
|:--|:--|:--|
| NaH (Na 1s 10.63 / 2s 3.3 / 2p 3.44 / 3s 0.836; H 1.0) | Na: k ≈ 2.84 → 1s 50, 2s 12, 2p 12, 3s 53. H: k = ζ_H → 1 | k ≈ 2.84 → max 53 (H 1s needs 36) |
| H₂O (O 1s 7.66 / 2s 2.25 / 2p 2.23; H 1.0) | O: k ≈ 3.94 → 24 / 23 / 22. H: 1 | k ≈ 2.82 → max 36 (O 1s 35, H 36) |
| LiH (Li 1s 2.69 / 2s 0.65; H 1.0) | Li: k ≈ 1.26 → 26 / 26. H: 1 | k ≈ 1.26 → max 26 |

At 1e-6 with per-centre k, NaH needs at most 35 terms, H₂O 15 and LiH 17 (with a global k, H₂O needs 23). The optimum sits near the geometric mean √(ζ_min ζ_max) of each centre's exponents, as the ρ ↔ 1/ρ symmetry predicts. The driver for the term count is the **core/valence exponent ratio on one centre**: 12.7 for Na. Each Sturmian channel set is shared by all STOs of that l on a centre. So the 1e-10 evaluator set is about 53 s-functions plus 12 × 3 p-functions on Na, and 1 function on H.

**Structural caveat that dominates the cost question.** An STO r^{N-1}e^{-ζr}Y_lm is an **exact finite** combination of N−l hydrogenic functions sharing the rate ζ (Z = n'ζ, n' = l+1..N). That is 1 term for 1s/2p/3d, 2 for 2s/3p and 3 for 3s. The infinite shared-k expansion is only needed if you insist on one k per centre. Both routes need the same engine change: one Z per orbital instead of one Z per centre (see Q2). So as an STO evaluator, the per-shell finite route dominates the shared-k one in every case.

## Q2 — Which integrals become exact, and where the walls are

**One-electron, two-centre (diatomic): exact.** Every class (S, T, V_A, V_B) is a two-centre integral of poly × exp functions. It closes through the Mulliken auxiliaries A_m(α), B_m(β) in prolate spheroidal coordinates. This is the same route `geovac/qfd_core.py` already uses for hydrogenic s orbitals.

I implemented it independently for s-type STOs (`two_center` in `sto_sturmian_eval_headtohead.py`, mpmath, no quadrature). Applied to the Sturmian-expanded STOs, the integrals converge geometrically onto the direct exact STO values (LiH, global k = 1.261):

| terms | max\|ΔS\| | max\|Δh\| | ΔE_FCI (Ha) |
|:--:|:--:|:--:|:--:|
| 8 | 7.8e-5 | 4.1e-5 | 3.4e-5 |
| 12 | 2.3e-7 | 2.9e-7 | 6.9e-8 |
| 16 | 2.3e-9 | 2.8e-9 | 1.2e-12 |
| 20 | 1.7e-11 | 1.7e-11 | <1e-11 |
| 32 | 4e-17 | 1e-16 | ~1e-15 |

The energy error is second order (variational), as expected. A deliberately poor k = 1.0 needs about 26 terms for 1e-11. Shared-k structure also makes T and the own-nucleus V trivial (k²δ − (k²/2)G and −Zk/n·δ). With a common k on both centres, every cross pair sits at β = (a−b)R/2 = 0, where B_m(0) is rational.

**Two-electron, diatomic: exact but not cheaper.**
- One-centre densities (AA|·) need no expansion: a product of two STOs on one centre is already a single poly × exp.
- Two-centre densities (the hybrid and exchange classes) must be expanded. The cost multiplies by the product of term counts. For NaH (Na shell × H 1s) | (Na shell × H 1s) exchange with per-centre k: up to 53² ≈ 2,800 Sturmian exchange integrals per STO integral. With a global k (H needs 36 terms): ~3.6e6. The per-shell finite route needs at most 9, and the current 6-Gaussian evaluator needs 6⁴ = 1,296 primitive quartets.

Two Sturmian-specific *structural* effects:
1. **The exchange τ-series terminates.** `qfd_core.exchange_closed_form` notes that the Neumann τ sum terminates exactly when both centres of a density carry the same exponent (q = 0). Under a single global k, every Sturmian quartet satisfies q = 0. So the heteronuclear exchange, which is non-terminating for LiH at 30 certified digits, becomes truncation-free per Sturmian quartet. The truncation simply moves from τ to the Sturmian index n.
2. **Equal rates everywhere land on the NaN branch.** Paper 58 records that the shell-route hybrid evaluator with l_a, l_b > 0 requires Z_B < Z_A strictly, and returns NaN when the rates coincide. A shared-k evaluator puts every pair on that branch, so the removable-limit branch would have to be built first.

Either way, `two_center_eri` / `qfd_assemble` currently assert **one Z per centre** (aabb, hybrid, exchange). Sturmians need Z = nk per orbital. The memory file already scopes that as "plumbing plus re-validation".

**Three-centre (polyatomic): the wall is untouched.** The expansion is per centre, so it never reduces the number of centres.
- The one-electron three-centre class ⟨χ_Y|1/r_X|χ_Z⟩ is already closed (v4.81.0) and is unaffected.
- The two-body T2 = (XY|XZ) stays elliptic **even at a fully shared exponent.** Paper 59's own reference case uses ζ = 1 on all three centres. Its two scales c₁ = s(1−s) and c₂ = t(1−t) are Feynman-parameter functions, degenerate only on s = t, so equal orbital exponents do not reach the degenerate diagonal.

So the wall bites exactly at Paper 58's scope-section BeH₂/H₂O RHF ladder ("Gaussian-fitted Slater shapes in the same engine"). That ladder cannot leave Gaussians by this route. The diatomic cases (H₂, LiH, NaH) are unaffected.

## Q3 — Head-to-head, Gaussian-fit vs exact

**H₂, 1s ζ = 1, R = 1.4.** Everything here is closed form: Mulliken for h; Slater/Ruedenberg/Sugiura formulas for the four ERI classes. The exact FCI is −1.1065566060914, which reproduces Paper 58 §qfd's certified value. Gaussian evaluator error, and the part of it that comes from h (exact g substituted):

| primitives | ΔE (Ha) | from h | from g | max\|Δh\| | max\|Δg\| |
|:--:|:--:|:--:|:--:|:--:|:--:|
| 4 | 2.67e-3 | 2.68e-3 | −1.1e-6 | 1.5e-3 | 4.2e-5 |
| 6 (LSQ) | 3.02e-4 | 3.02e-4 | −3.4e-7 | 1.7e-4 | 1.2e-6 |
| STO-6G (production) | 3.04e-4 | 3.05e-4 | −1.3e-7 | 1.7e-4 | 1.2e-5 |
| 8 | 4.5e-5 | 4.5e-5 | −9e-8 | 2.6e-5 | 1.1e-7 |
| 12 | 1.9e-6 | 1.9e-6 | −2e-9 | 1.1e-6 | 3e-9 |

This reproduces Paper 58's sequence (2.7e-3, 3.0e-4, 4.5e-5, 8.0e-6, 1.9e-6) and its reading. **More than 99.9% of the fit error is in the one-electron block (the nuclear cusp).** The error scales about ×7 per two primitives. The ERIs are already at 1e-6 to 1e-9.

**LiH, s-only STO (Li 1s 2.69, 2s 0.65; H 1.0), R = 3.25, Paper 58's fixed-geometry case.** Reference: exact h, g at 12 primitives. That gives E_FCI = −7.9672764 total, NOCI-3 = −7.9659720.

| evaluator | ΔE_FCI | ΔE_NOCI3 | max\|Δh\| | max\|Δg\| vs 12 |
|:--|:--:|:--:|:--:|:--:|
| LSQ 4 | 2.90e-2 | 2.91e-2 | 1.4e-2 | 5.0e-4 |
| LSQ 6 (test_paper58 `_lih_system`) | 3.30e-3 | 3.32e-3 | 1.5e-3 | 6.1e-5 |
| production (STO-6G 1s + STO-3G 2s) | 2.91e-3 | 2.83e-3 | **6.4e-2** | **2.7e-2** |
| LSQ 8 | 4.8e-4 | 4.9e-4 | 2.3e-4 | 7.2e-6 |
| LSQ 12 | 2.0e-5 | 2.0e-5 | 9.4e-6 | — |

- **The fit error is roughly 3 mHa at the depth Paper 58 uses.** It scales as ζ² through the Li 1s cusp.
- The **production** shapes have element errors 40 times larger (max\|Δh\| 6.4e-2, max\|Δg\| 2.7e-2). STO-3G 2s is a *different function* (⟨fit|STO⟩ = 0.99968, versus 1 − 6e-7 for the 6-Gaussian fits). Its energy error is no bigger only because of variational cancellation.
- The fixed-geometry numbers Paper 58 quotes are differences within one evaluator, so they are robust. NOCI-3 sits 1.30 mHa above FCI with exact h, against 1.32 mHa printed.

**Cost.** In mpmath at 60 digits, the exact Mulliken h for this 3-orbital set takes 0.01 s directly. Expanding in Sturmians costs 1 s at 16 terms and 12 s at 32 terms (it scales as the fourth power of the term count per pair). For comparison, the NaH 7-orbital MD build was measured at 0.27 s with 6 primitives and 2.6 s with 12.

## Q4 — Is the fit error material in any Paper 58 number? Yes, at printed precision; no qualitative change

Setup: the NaH ladder exactly as in `tests/test_paper58_nah_ladder.py` (atom-only ζ; the paper's 0.25 a₀ grid; its 3-point parabola). It was re-evaluated with the paper's production shapes ("prod") and with engine least-squares fits of **all four** shapes at 6 to 12 Gaussians. The L-sequence converges, and L8 already agrees with L12 to 1e-4 eV; L12 stands in for exact STO. "prod" reproduces Table II to every printed digit.

| quantity (Table II / text) | printed (prod) | converged (L12) | shift |
|:--|:--:|:--:|:--:|
| E(Na), atom-only provenance | −161.0845 | **−161.1192** | −34.7 mHa |
| 3-det R_eq (a₀) | 3.736 | 3.744 | +0.008 (+4.8% → +5.0% vs expt) |
| 3-det D_e (eV) | 1.071 | 1.082 | +0.011 (+1.0%) |
| 3-det % FCI | 91.1 | 92.0 | +0.9 pp |
| covalent 2-det R_eq / D_e / % | 3.954 / 0.686 / 58.4 | 3.963 / 0.695 / 59.2 | |
| all 4 det | 3.713 / 1.080 / 91.9 | 3.722 / 1.091 / 92.9 | |
| FCI (91 det) | 3.595 / 1.175 / 100 | 3.595 / 1.175 / 100 | none |
| Löwdin 3-det R_eq / D_e / % | 3.658 / 0.136 / 11.6 | 3.679 / **0.163** / **13.9** | D_e +20% relative |
| fixed R = 3.5: 3-det % / Löwdin 3-det % | 89.3 / 10.8 | 90.0 / 12.9 | |

Readings:
1. **Every bold Table II headline moves in its last printed digit or beyond.** R_eq 3.736 → 3.744, D_e 1.071 → 1.082, 91.1 → 92.0. The Löwdin 3-det D_e moves by 20% relative.
2. **Nothing qualitative changes.** Three determinants still bind and recover over 85% of FCI binding; Löwdin truncation damage is still severe (13.9% against 92.0%); FCI is untouched. The paper's central claims stand.
3. **The E(Na) comparison is entirely fit artifact.** Paper 58 contrasts E(Na) = −161.0845 with "a Clementi minimal-STO reference of ≈ −161.12". The converged-evaluator value is −161.1192. The apparent 35 mHa gap is the STO-6G/STO-3G evaluator, not the basis.
4. **The culprit is the shape choice, not the depth of the 6-Gaussian fits.** L6 (6-Gaussian least-squares fits everywhere) already gives D_e 1.081 and 92.0%. The shift comes from the hardcoded STO-3G 2s (3 primitives, overlap 0.99968) and STO-6G 1s, which `_shapes()` uses for NaH instead of `fit_sto_shape`.
5. **Prose discrepancy to flag.** Paper 58 §nah says "the radial functions are six-Gaussian fits of Slater shapes". `test_paper58_sto_fit_quality` guards ⟨fit|STO⟩ = 1.000000 only for the 2p/3s least-squares fits. The NaH 2s shape actually used is a 3-Gaussian STO-3G contraction at 0.99968. The paper argues the fits are essentially perfect, so the D_e deficit must be basis incompleteness. That conclusion survives, since D_e moves by 1%, but it rests on a misdescription of which fits were used.
6. **Not affected:**
   - The census g-row density corroboration: zeros there are threshold decisions on symmetry-permitted entries, not fit-sensitive values.
   - The fixed-geometry LiH ratios: within-evaluator differences.
   - The H₂O/BeH₂ RHF ladder: reported "as scope", and it is a three-centre case anyway.

## Bottom line for the PI

1. **Cheap, high-value fix (a PI call; paper-number edit).** Evaluate the NaH ladder with `fit_sto_shape` at 8 or more primitives for every shape. Then either re-quote Table II, E(Na) and the Löwdin row (3.744 / 1.082 / 92.0; 0.163 / 13.9; −161.119), or correct the "six-Gaussian fits" sentence. No Sturmian machinery is needed for this.
2. **Sturmian evaluator: BORDERLINE.**
   - It is exact and convergent, and the term counts are fine for diatomics (at most about 53 per STO for NaH at 1e-10).
   - It gives no accuracy beyond what deeper fits give. It is dominated by the per-shell finite hydrogenic route, which needs the same per-orbital-Z plumbing.
   - Its ERI cost multiplies with the term counts.
   - It cannot help polyatomics.
   - Its genuine gains are decidability: rational B_m(0), and an exchange τ-series that terminates for every quartet under one k. Those matter only if the goal is certifying Paper 58 NaH digits Lindemann-style, as §qfd did for H₂/LiH.
   - If that is the goal, build the per-orbital-Z unblock in `two_center_eri`/`qfd_assemble` first, together with the equal-rate removable-limit branch. Then choose per-shell (finite) or shared-k (terminating τ) per class.
