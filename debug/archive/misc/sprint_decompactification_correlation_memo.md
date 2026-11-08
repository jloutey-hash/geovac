# Sprint memo — Does electron–electron interaction move the two-center decompactification front only through Z → Z_eff?

**Date:** 2026-09-06 · **Type:** DIAGNOSTIC (measurement only) · **Systems:** H₂⁺/HeH²⁺ (Rung 0), H₂ and HeH⁺ (Rungs 1–2)
**Driver:** `debug/decompactification_correlation_ladder.py` · **Data:** `debug/data/decompactification_correlation_ladder.json`
**Predecessor:** `debug/sprint_decompactification_R_sweep_memo.md` (H₂⁺: per-shell front R\*(n) ≈ 2.14·ℓ_n, ℓ_n = n/Z)

## Verdict: **GO at the principal-angle level — with a FINDING at the coherence level**

- **GO (the hypothesis as posed).** The decompactification front of the *occupied one-body space* — the R at which the
  occupation-weighted principal angle between the A- and B-centered components of the natural orbitals reaches 45°
  (measure M1, the object of Paper 60's compound-matrix identity) — matches the one-electron 1s–1s law evaluated at
  the *empirical* Z_eff(R) to **−0.7 % (H₂, HF), −1.2 % (H₂, FCI), 0.0 % (HeH⁺, HF), −0.2 % (HeH⁺, FCI)**; the
  pointwise residual |Δcos(R)| ≤ 0.006 over the whole R-range where the second center holds > 1 % of an electron.
  Both rungs, both molecules, three basis sets, two ERI-grid resolutions. Screening moves the front; nothing else does
  at this level, to about 1 %.
- **FINDING ("in some other way").** The *signed cross-center coherence* |γ_AB|/√(γ_AA γ_BB) (measure M2) has the same
  front at HF (identical to M1 for one MO) but at FCI its front moves **inward by 23 % (H₂: 0.987 vs 1.285) and 9 %
  (HeH⁺: 0.690 vs 0.758)** — a shift not expressible as a decay length. Mechanism (isolated below): the antibonding
  natural orbital enters the per-center density with weight amplified by (1+S)/(1−S) relative to the bonding one, so
  even n_u = 0.025 (H₂ at R = 1.4) pulls the coherence from S to a value between S (MO limit) and S² (Heitler–London
  limit). This is occupation redistribution (bond-order collapse), an *occupation-number* effect, not a Z_eff effect.
- Hard control **passed**: ℓ_eff is a measured output (varies 1.46 → 0.83 across R at HF; HF and FCI differ by 12 %
  at R = 4; equals no basis exponent); basis doubling moves R\* by ≤ 0.8 % and ζ_eff by ≤ 0.1 % (H₂) / ≤ 0.8 % (He center).

## Design (measures stated explicitly)

**Substrate (Rungs 1–2).** Per-center even-tempered 1s-STO sets: H {0.5, 0.75, 1.125, 1.688, 2.531}, He {0.9, 1.35,
2.025, 3.038, 4.556}; "double" = 10 interleaved over the same span; "wide" = 8 over 0.35–3.62 / 0.63–6.52. Heavy
atom at A (origin). One-electron integrals closed-form (same center) or prolate Gauss–Legendre (cross center; overlap
vs exact Mulliken form 7e-14); kinetic via the STO identity −½∇²χ_ζ = (ζ/r_c − ζ²/2)χ_ζ, symmetrized. ERIs: multipole
expansion about A on a graded radial grid (the Paper 60 route) — validated: He HF −2.86167 (limit −2.8617), He FCI
−2.87888 (s-limit −2.8790), (AA|BB) vs `two_center_eri.aabb_value` 8e-5, minimal-basis FCI vs the Paper 60 driver
1.1e-4, H atom −0.499997. Canonical orthogonalization at S-eigenvalue 1e-7. RHF and singlet FCI (spatial product
basis, 1-RDM γ = 2CCᵀ). Natural orbitals split into A-block and B-block components c_k^A, c_k^B.

- **ℓ_eff per center** from the dominant NO's component φ_k^X(r) = Σ_{μ∈X} c_μ χ_μ: route A = log-slope of φ over its
  [r₅₀, r₉₅] cumulative-norm window; route B = ζ = 3/(2⟨r⟩) (exact for a pure 1s); route C = the same moment of the
  XX-block density (occupation-weighted; ≡ B at HF).
- **Measured front:** M0 = |cos θ| of the dominant NO's A/B components; **M1** = Σ n_k|cos θ_k| / Σ n_k (occupation-weighted
  principal cosine); **M2** = Tr(P_AB S_BA)/√(Tr(P_AA S_AA)Tr(P_BB S_BB)) (signed coherence; M2 = M1 for a single MO);
  M3 = largest singular value between span{φ_k^A} and span{φ_k^B} over NOs with n_k > 0.01 (occupation-blind).
  R\* = R at which the measure crosses 1/√2 from above (PCHIP + brentq on the R grid).
- **Predicted front:** cos_pred(R) = exact ⟨1s(ζ_A)_A | 1s(ζ_B)_B⟩ at the empirical ζ_X(R) — the Rung-0 one-electron law at
  Z_eff, evaluated exactly rather than through a fitted constant. R\*_pred = its 1/√2 crossing.

## Rung 0(a) — one-electron heteronuclear law (control)

Grid (Z_A,Z_B) ∈ {(1,1),(2,1),(3,1),(2,2)} × n_A,n_B ∈ 1..4, exact overlaps. Two facts a naive "R\* = c(ℓ_A+ℓ_B)" misses:

1. **The absolute crossing |S| = 1/√2 exists only if the united-atom overlap is above it.** For 1s–1s exponents (a,b),
   S(0) = (2√(ab)/(a+b))³, so the front exists only for exponent ratio t < t_c = **2.746**. On the grid only 13/64 entries
   cross: the same-Z, same-n ones (reproducing yesterday: R\*/ℓ = 1.565, 2.194, 2.189, 2.130 for n = 1..4, both Z) and
   a single heteronuclear same-n pair, (Z=2,1; 1s,1s): **R\* = 0.754** — *below* the homonuclear Z=2 value 0.782 and
   half the additive prediction 0.78×(0.5+1.0) = 1.17. The dense 1s–1s ratio scan (a=1, b=t) has R\* = 1.565, 1.278,
   0.901, 0.388 at t = 1, 1.35, 1.81, 2.44 and none beyond t_c: **for mismatched shapes the absolute front is not a
   tail-reach quantity** — the R=0 principal angle already spends most of the 45° budget.
2. **The tail-reach front R_rel (|S| falls to |S|_max/√2 past its peak) combines the two decay lengths as a near-geometric
   mean, not additively and not max-dominated.** 1s–1s over the four charge pairs: R_rel = 0.79·2√(ℓ_Aℓ_B)
   (max dev 2.0 %; free power p = q = 0.49). Dense scan t ∈ [1,8]: geometric max dev 9 %, additive 23 %, max-law 62 %,
   harmonic 47 %; free power p = q = 0.44 (rms_log 0.011). R_rel(t) = 1.565, 1.173, 0.896, 0.701, 0.628 at
   t = 1, 1.81, 3.28, 5.94, 8. Mixed-n pairs are near-orthogonal at all R (|S|_max ≤ 0.18) and mixed-Z n ≥ 2 pairs
   have multi-peak overlaps whose "relative" crossings (up to 36 bohr) are far-tail artifacts — not front-like, excluded.

Since H₂'s components are equal by symmetry and HeH⁺'s have ζ_He/ζ_H ≈ 1.05 at the front (t ≪ t_c), the Rung-1/2
prediction uses the exact 1s–1s overlap directly; no fitted constant is load-bearing.

## Rung 0(b) — HeH²⁺ η-equation l-mixing vs R (b = R(Z_B − Z_A) = −R)

Ground η-eigenvector of −l(l+1) + c²η² + bη (N = 50 Legendre basis; A matches the prolate solver and
`_angular_sep_const` to 1e-13; the homonuclear recomputation matches yesterday's curve to 4e-16), c²(R) from the exact
prolate solver (Z_A=2, Z_B=1; E_elec → −2 − 1/R, He⁺ + proton).

| R | HeH²⁺ part.-deficit | entropy (bits) | ⟨η⟩ | H₂⁺ part.-deficit | entropy |
|--:|--:|--:|--:|--:|--:|
| 0.5 | 0.024 | 0.17 | −0.18 | 0.0001 | 0.002 |
| 1.0 | 0.117 | 0.57 | −0.39 | 0.0014 | 0.015 |
| 2.0 | 0.440 | 1.35 | −0.71 | 0.013 | 0.10 |
| 4.0 | 0.578 | 1.89 | −0.87 | 0.113 | 0.52 |
| 8.0 | 0.712 | 2.38 | −0.94 | 0.475 | 1.15 |
| 12.0 | 0.752 | 2.67 | −0.96 | 0.437 | 1.35 |

Shape difference, two structural features: (i) the heteronuclear mixing turns on at **first order** (bη couples l ↔ l±1;
deficit ∝ R^2.1 at small R) where the homonuclear one is second order (c²η² couples l ↔ l±2; ∝ R^4.0); (ii) the
heteronuclear deficit does **not saturate near 0.5** — it keeps rising (0.75 at R = 12, a shoulder at R ≈ 2.7–3.6) as the
electron localizes on the He focus (⟨η⟩ → −1), i.e. the η-label is lost to *localization*, not to gerade mixing. Both
curves are continuous.

## Rungs 1–2 — H₂ and HeH⁺

R\* = front (bohr); ζ_eff = route-B exponent at R\*; "pred" = R\* of the exact 1s–1s law at the empirical ζ_eff
(routes B / A / C); residual relative to pred(B). Base basis; doubled and wide in brackets.

| molecule · rung | R\*(M1) | R\*(M0) | ζ_A, ζ_B at R\* | R\*_pred (B / A / C) | residual M1 | R\*(M2) | residual M2 |
|:--|--:|--:|--:|--:|--:|--:|--:|
| H₂ · HF  | 1.283 [1.283, 1.282] | 1.283 | 1.21, 1.21 | 1.292 / 1.291 / 1.292 | **−0.7 %** [−0.7, −0.7] | 1.283 | −0.7 % |
| H₂ · FCI | 1.270 [1.269, 1.269] | 1.274 | 1.22, 1.22 | 1.285 / 1.285 / 1.283 | **−1.2 %** [−1.2, −1.2] | 0.987 [0.987, 0.989] | **−23.1 %** |
| HeH⁺ · HF  | 0.757 [0.751, 0.755] | 0.757 | 2.11, 2.02 | 0.757 / 0.750 / 0.757 | **0.0 %** [−1.0, +0.1] | 0.757 | 0.0 % |
| HeH⁺ · FCI | 0.757 [0.751, 0.754] | 0.758 | 2.12, 2.01 | 0.758 / 0.752 / 0.759 | **−0.2 %** [−1.1, −0.2] | 0.690 [0.685, 0.690] | **−9.0 %** |

**Residual Δcos(R) = M1 − cos_pred(R), the (iv) curves.** H₂ HF: −0.004 (R=0.5) → −0.002 (R\*) → 0 (R≈2.4) → +0.008
(R=6). H₂ FCI: −0.004 to −0.005 for R ≤ 4, → −0.0006 at R = 6. HeH⁺ HF: −0.006 (R=0.3) → +0.001 (R=0.9–1.1) → −0.006
(R=2.6–3.0) → −0.013 at R = 5 where the H center holds 0.2 % of an electron. HeH⁺ FCI: −0.008 → 0 → −0.004 → −0.012.
Sign at the front: negative at every rung and molecule (the measured overlap is slightly *smaller* than a single
exponential of the same moment), FCI ≈ 2× HF in magnitude for H₂; all sub-percent in cos.

**The coherence front, and why it moves.** Writing w_k² for the per-center norm of NO k, M2 = Σ n_k w_k² cos θ_k /
Σ n_k w_k². In the minimal-basis form w_g² = 1/(2(1+S)) and w_u² = 1/(2(1−S)), so the antibonding NO is weighted by
(1+S)/(1−S) ≈ 5 at S = 0.67: at R = 1.4 (n_g = 1.969, n_u = 0.025, cos_g = +0.674, cos_u = −0.68) this formula gives
0.591 against the measured 0.596. Limits: MO (n_u = 0) M2 = S; Heitler–London M2 = S². H₂ runs from the first toward
the second as n_u grows (0.008 at R = 0.5, 0.064 at 2.0, 0.51 at 4.0, 0.885 at 6.0; M2 = 0.0047 vs S = 0.047 at R = 6).
HeH⁺ stays near the MO limit (n_u ≤ 0.015) and shifts 9 %. This is the one place correlation proper — not screening —
moves a front, and it does so through occupation numbers, i.e. through the *k-th compound* of the same one-electron
overlap, exactly the object Paper 60 §manyelectron says does not compound *as a metric*. The subspace route M3 is
occupation-blind: at FCI it admits the σ_u component once n_u > 0.01 and its span contains a combination with larger
cross-overlap (H₂ front 3.15), a property of the span, not of the occupied density; it coincides with M1 at HF and for
HeH⁺.

## Hard control (basis does not pin ℓ_eff)

1. **ℓ_eff is an output.** H₂: ζ_eff(HF) = 1.459 (R=0.5) → 1.186 (1.4; the textbook optimal H₂ exponent 1.19) →
   0.828 (6.0); ζ_eff(FCI) = 1.454 → 1.193 → 0.997 (returns to the H atom). HF vs FCI differ by 0.6 % at R = 1.4,
   12 % at R = 4, 20 % at R = 6 (HF's per-center orbital keeps diffusing as HF dissociates to its ionic mixture).
   HeH⁺: He center 2.43 → 1.85 (R=1.46) → 1.62 (R=5; He s-limit single-ζ 1.69); H center 2.33 → 1.64 → 0.68 (population
   0.002). None equals a basis exponent.
2. **Basis doubling** (5 → 10 per center; 6–7 of 20 dropped at 1e-7, S_min = 5e-13): ΔR\*(M1) = −5e-5 (H₂ HF), −6e-4
   (H₂ FCI), −0.006 = 0.8 % (HeH⁺); max|Δζ_A| = 0.09 % (H₂), 0.8 % (He center); |Δζ_B| ≤ 3.5 % only at R ≥ 3 where the
   H center holds < 5 % of an electron; ΔE ≤ 1.2 mHa. **Wide span** (8 per center, 0.35–3.62 / 0.63–6.52): ΔR\* ≤ 0.0024,
   Δζ ≤ 0.1 % / 0.5 % / 2.9 % (same qualifier). Residuals unchanged to 0.1 % across all three sets.
3. **Two routes for ζ** agree to 0.3 % at the front (route A 1.210 vs route B 1.209 at R = 1.3) and ≤ 2 % elsewhere; route C
   ≡ B at HF and within 1 % at FCI. Fit windows [r₅₀, r₉₅] = [1.1, 2.6] bohr at R = 1.3 (H₂), R² ≥ 0.9997 everywhere
   except the 0.2 %-populated H center at R = 5 (0.994). *Asymptotic-regime statement, honest version:* the local
   exponent drifts 1.23 → 1.17 across the window and 1.11 at 2r₉₅ — the component is not one exponential (this drift is
   what the Δcos residual measures), and the r → ∞ tail of a finite even-tempered set is its smallest exponent, a basis
   artifact, so the window is the norm-carrying region, not the true asymptote.
4. **ERI-grid resolution control** (base basis, nr 3000 → 6000, nth 300 → 600, R ∈ {0.6, 0.8, 1.0, 1.3, 1.46, 2.0, 4.0,
   6.0}): ΔM1 ≤ 1.6e-4 (HF, R = 6) and ≤ 5e-5 (FCI); ΔM2 ≤ 2e-5; Δζ ≤ 0.03 %; ΔE_FCI ≤ 0.16 mHa (H₂), ≤ 0.03 mHa
   (HeH⁺). Fronts on the control sub-grid: H₂ HF 1.2824 / FCI 1.2692 / M2 0.9872, HeH⁺ HF 0.7568 / FCI 0.7564 /
   M2 0.6892 — all within 0.001 bohr of the standard-resolution values. Near the front the gates improve 4×
   (moment overlap 4e-4 → 1e-4, (aa|aa) 1.8e-3 → 4.4e-4 at R = 1.3). The large-R (aa|aa) gate does *not* improve
   (2.5e-2 → 4.9e-2 at R = 6): it is the L_max = 24 multipole truncation of the steepest B-centred function's
   self-repulsion far from the expansion center, not grid resolution. **L_max control** (24 → 40 at R ∈ {1.3, 4, 6}
   H₂ and {0.8, 4} HeH⁺): every reported quantity changes by ≤ 2e-5 (ΔE < 1e-6 Ha, Δζ < 1e-3 %), but the gate itself
   does *not* converge — it worsens (5.1e-3 → 1.6e-2 at R = 4; 2.5e-2 → 2.9e-2 at R = 6). So the self-repulsion of the
   steepest B-centred function at R ≥ 4 is a 1–3 % defect of the single-center multipole route (steep density far from
   the expansion center) that neither grid nor L_max cures; it is harmless here only because that function's
   variational weight at large R is negligible (the controls bound its effect at 2e-5 in every front, exponent and
   coherence). Recorded as a limitation of the route, not certified away.

## Cross-checks

- **Energies.** H₂ R = 1.4: HF −1.12837 (s-only HF), FCI −1.15451 (exact −1.17447; `solve_level4_h2_multichannel`
  l_max=2 −1.15763); D_e(basis-consistent) = 0.1545 Ha = 89 % of exact; R = 6: −1.00033 → two H atoms. HeH⁺ R = 1.46:
  FCI −2.94247 (exact −2.9787); against the same basis's He limit −2.87888, D_e = 0.0636 Ha = 85 % of 0.0750; R = 5:
  −2.87912 → He + H⁺. `solve_level4_h2_multichannel(Z_A=2, Z_B=1, l_max=2, n_alpha=80, n_Re=150)` gives −2.76541 at
  R = 1.46 and −2.59481 at R = 2.0 — the adiabatic midpoint-origin setting is not converged for HeH⁺ (Paper 15 needs the
  charge-center origin and l_max = 4); recorded as reference points only, not as a validation of either side.
- **Measured front, two routes:** M0 vs M1 agree to 0.3 % (1.274 vs 1.270); M3 = M1 at HF.
- **Integral gates near the front** (R ≈ 0.7–1.5): L=0-moment overlap vs exact ≤ 4e-4, steepest (aa|aa) vs 5a/8 ≤ 2e-3,
  potential of each A density at B vs closed form ≤ 2e-5, ERI-matrix asymmetry 4e-15, FCI C-matrix asymmetry ≤ 3e-8.
  The moment/(aa|aa) gates degrade to 8e-3/2.5e-2 at R = 6 (angular width 1/(ζR) of the steepest B function); see
  control 4.

## Interpretation

At the level at which the front was defined — the principal angle between the two single-center projectors, now
generalized to the occupied one-body space — many-body physics moves it only by changing the exponential reach of the
per-center component the electrons actually occupy: Z → Z_eff(R), to about 1 %, at HF and at FCI, homonuclear and
heteronuclear. The residual is a sub-percent *shape* term (a contracted multi-ζ component is not one exponential) with a
definite sign, not a screening term. What correlation does move is a different object: the signed coherence, i.e. the
bond-order-weighted overlap, and it moves it through occupation numbers of the antibonding natural orbital — the
compound-matrix structure Paper 60 identifies — not through any length scale. Rung 0 adds two one-electron facts the
prediction needed: mismatched per-center shapes have no absolute front at all beyond exponent ratio 2.75, and the
tail-reach front of two unequal tails is their geometric mean, not their sum.

## Caveats / not computed

- s-only per-center basis (no pσ polarization); the front is an s-shell object and the substrate binds correctly
  (89 % / 85 % of D_e), but a pσ control was not run.
- Two electrons only; the M2 mechanism is written in its two-orbital form. Whether the principal-angle GO survives at
  N > 2 (where several NOs are strongly occupied) is untested.
- HeH⁺ at R ≥ 3 has < 5 % of an electron on H; ζ_B and Δcos there are reported but not load-bearing.
- The M3 subspace measure is reported as an occupation-blind route only; its FCI crossing is a threshold artifact.
- R\* is a grid interpolation (PCHIP); grid spacing 0.1 near the fronts.
- Guardrail (§3.5, Papers 8–9) acknowledged: genuine two-center machinery, R-dependence measured; no single-center /
  shared-p₀ encoding proposed and no binding claim made. No paper, CLAUDE.md, or test edits.


## PM correction (2026-09-06, on capture)

The Rung 0(a) threshold **t_c = 2.746 is a scan-grid artifact**: the driver took the first t on its discrete ratio scan at which the crossing was lost. The closed form stated in the same sentence, S(0) = (2 sqrt(ab)/(a+b))^3 = 1/sqrt2, has root **t_c = 2.664**, and an independent prolate-quadrature route (`tests/test_paper58_decompactification_front.py`) confirms S(R) is monotone decreasing for every t, so the united-atom value is the maximum and the crossing is lost exactly at 2.664 (present at t = 2.65, absent at 2.68). Papers, registry and CHANGELOG carry 2.664. Every other number in this memo reproduces from the data file.
