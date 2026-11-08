# Sprint memo: the prolate separation constant Λ as an exact label (H2+ → H2 → LiH)

**Date:** 2026-09-26. **Verdicts:** A **PASS** (the physics is textbook; see §4) · B **GO** (a diagnostic/compression label, not a conserved one) · C **STOP** (the gate fails on the Route C LiH orbital set). The guardrail was acknowledged and respected: all work stays in two-centre prolate coordinates (Level 2). There is no single-centre or shared-p₀ basis and no graph concatenation.

## 1. Construction

One electron, with Z_A at z = −R/2 (η = −1) and Z_B at z = +R/2. Let Z_± = Z_A ± Z_B.

    H = (F + G)/(ξ²−η²),   F = (2/R²)F̃,  G = (2/R²)G̃
    F̃ = −∂ξ(ξ²−1)∂ξ + m²/(ξ²−1) − R Z₊ ξ
    G̃ = −∂η(1−η²)∂η + m²/(1−η²) + R Z₋ η
    Λ = (ξ² G̃ + η² F̃)/(ξ²−η²)

Λ is dimensionless. On an eigenstate it equals the η-equation separation constant A, including the −p²η² term (p² = −R²E/2).

In a product basis f_a(ξ) g_b(η) e^{imφ}, every operator is a Kronecker combination of 1D blocks:

    W  = r2⊗1 − r0⊗a2                  (the ξ²−η² Jacobian; the metric)
    Hm = (2/R²)[Fx⊗1 + r0⊗Gy]           H c = E W c
    N  = r2⊗Gy + Fx⊗a2                  Λ c = λ W c

Here Fx = Kx − R Z₊ r1 and Gy = diag l(l+1) + R Z₋ a1. **These three matrices are exactly Atkinson's operator determinants Δ₀, Δ₁, Δ₂ of the two-parameter eigenproblem** (F̃ − εξ² + λ)x = 0, (G̃ + εη² − λ)y = 0, with ε = ER²/2. Atkinson's theorem then gives two results **at every finite truncation**, not merely in the limit:

- W⁻¹Hm and W⁻¹N commute exactly.
- Their joint eigenvectors are decomposable, x ⊗ y.

The one-body commutation is therefore exact at every cutoff. That is the same pattern as the SO(4) Casimir on the (n,l,m) labels.

**Blocks.**

- **η:** normalized associated Legendre P_l^m. The kinetic block is diag l(l+1) exactly. The η and η² blocks are closed form (square roots of rationals, the same reduced element as `_az_element` without its n² factor).
- **ξ:** (ξ²−1)^{m/2} e^{−α(ξ−1)} L_a^{(m)}(2α(ξ−1)). The integrands are polynomial × e^{−t}, so a Gauss–Laguerre rule with N+m+6 nodes is **exact**.

**Independent-route cross-check.** The 1D blocks were rebuilt with the production mpf monomial-moment route (`prolate_recondition._radial_blocks/_angular_blocks`, Gegenbauer family). The unit-normalized blocks agree to **2.8×10⁻¹⁴** for m = 0, 1, 2 and α = 0.7, 3 (`debug/prolate_separation_constant_xcheck.py`).

**Transcendentals.** In stage A the matrices contain only √rationals, the Laguerre roots (algebraic), and the free scale α. None of these is named in the Paper 18 sense. The eigenvalues E(R) and λ(R) are transcendental functions of R: the usual spheroidal content.

In stages B and C the imported engines bring their own seeds:

- e^{−c} from ∫₁^∞ ξ^k e^{−cξ};
- i_l(α) from the Legendre expansion of e^{−αη} in the STO cores.

Both are pre-existing and nothing new was introduced.

## 2. Stage A: one electron (H2+, and Z_A = 3, Z_B = 1)

Driver: `debug/prolate_separation_constant_stageA.py`. Log and JSON: `debug/data/prolate_separation_constant_stageA.{log,json}`.

- **Grid:** 90 log-spaced R ∈ [0.01, 80]; m = 0, 1, 2; basis 30 ξ × 40 η.
- **Basis converged:** at R = 80 the 30×40 and 40×70 bases differ by 7×10⁻¹⁴ Ha (H2+) and 1×10⁻⁴ Ha (3,1).
- **Energy checks:** H2+ E(R=2) = −1.1026342 and E(R=1) = −1.4517863, which are the standard values.

### 2.1 Commutation

Relative ‖[Ĥ, Λ̂]‖₂ / (‖Ĥ‖‖Λ̂‖), measured in the W-orthonormal frame:

| | R ≥ ~1 | R = 0.01 (cond W ~10⁶) |
|:--|:--|:--|
| H2+ | 1×10⁻¹⁶ – 2×10⁻¹⁴ | ≤ 9×10⁻¹⁰ |
| (3,1) | 7×10⁻¹⁷ – 2×10⁻¹⁴ | ≤ 5×10⁻¹⁰ |

The small-R residue is float64 rounding amplified by cond(W), not truncation.

**Control: Λ with the wrong sign on Z_A − Z_B.** For (3,1) the commutator becomes 3.7×10⁻⁴ (R ≈ 0.94), 3.9×10⁻³ (R = 4.3) and 3.9×10⁻² (R = 80). For H2+ it is identical, as it must be, because Z₋ = 0.

**Separability certified.** Two further checks hold for every tracked state at every R:

- Joint eigenvectors reshape to rank 1: σ₂/σ₁ ≤ 10⁻⁵ for all states, and ≤ 10⁻¹¹ for the low ones.
- The joint pair (E, λ) is an exact eigenpair of both separated 1D Sturm–Liouville problems, to ≤ 2×10⁻⁷ and typically 10⁻¹¹.

**State labels.** Each state is labelled by its **Klein-oscillation nodal indices (n_ξ, n_η)**, read as eigenvalue indices of the two 1D problems. These labels are R-invariant, so they track states through every crossing without any interpolation choice.

### 2.2 United-atom limit (R → 0)

The labels should satisfy l = n_η + m and n = n_ξ + n_η + m + 1. They do:

- At R = 0.01, λ − l(l+1) is −3×10⁻⁵ (1sσ) or smaller for H2+. Here the correction is O(R²).
- For (3,1), λ − l(l+1) is at most 2×10⁻⁴. The correction is O(R), because the R Z₋ η term is first order.
- E → −Z_UA²/2n² in both cases (Z_UA = 2 for H2+, 4 for the (3,1) system).

**All 75 tracked states (49 H2+, 26 (3,1)) reproduce l(l+1) and the united-atom level.** Λ → L² about the midpoint is confirmed.

### 2.3 Separated-atom limit (R → ∞), checked against the SO(4) Runge–Lenz construction

Asymptotic analysis of the η equation near η = ∓1 (derivation in the stage-A driver docstring and §4) predicts

    λ + p² = R·(Z_other + A_z) + O(1),    A_z = (Z_self/n)·k,

where k is an eigenvalue of the energy-scaled Runge–Lenz A_z on shell n at fixed m. The fit is λ + p² = sR + c₀ + c₁/R on R = 65.4, 72.3, 80. The k spectrum is taken from `so4_generators` (`tests/test_graph_s3_nonconvergence.py`), which builds A_z from Pauli's reduced elements; no integer is typed in by hand.

**(3,1).** All 26 tracked states sit on the A_z spectrum to |Δk| ≤ 0.008 (≤ 0.004 except one n = 6 state), for shells up to n = 6 on the Z = 3 centre plus H(1s). The first-order Stark shift, ΔE = −(3n/2Z_self)·k·Z_other/R², gives k values with **the same sign in every case** and magnitudes within 3–15%; the residue is the polarizability tail at finite R.

**Sign of k.** Positive k means the electron is displaced toward the partner nucleus. So Λ's leading large-R term is exactly the operator that diagonalizes the partner's Stark field: the parabolic basis.

**H2+.**

- Every state that correlates to n ≤ 3 lands within |Δk| ≤ 0.012.
- States correlating to n = 4–5 have not reached the asymptotic regime at R = 80 (|Δk| up to 0.2). This is the known slow n⁶-scaled polarizability approach, not a failure of the label.
- Each g/u partner pair shares one (n, k), as it must.

### 2.4 Correlation table (lowest states, m = 0, 1)

Columns:

- UA and physical parity (for H2+, inversion parity = η-parity × (−1)^m);
- nodal labels (n_ξ, n_η);
- λ at R = 0.01, and at R = 1.91 (near equilibrium);
- the separated-atom (SA) centre and n;
- k from Λ, and the A_z spectrum on that shell.

**H2+**

| UA | (n_ξ,n_η) | λ(0.01) | λ(1.91) | SA | k from Λ | A_z spec |
|:--|:--|:--|:--|:--|:--|:--|
| 1sσg | (0,0) | 0.0000 | −0.749 | H n=1 | −0.0000 | {0} |
| 2pσu | (0,1) | 2.0000 | 1.269 | H n=1 | −0.0000 | {0} |
| 2sσg | (1,0) | 0.0000 | −0.228 | H n=2 | −0.9998 | {±1} |
| 3pσu | (1,1) | 2.0000 | 1.720 | H n=2 | −0.9998 | {±1} |
| 3dσg | (0,2) | 6.0000 | 5.779 | H n=2 | +0.9997 | {±1} |
| 4fσu | (0,3) | 12.000 | 11.883 | H n=2 | +0.9997 | {±1} |
| 3sσg | (2,0) | 0.0000 | −0.110 | H n=3 | −1.9972 | {0,±2} |
| 4dσg | (1,2) | 6.0000 | 5.877 | H n=3 | +0.0002 | {0,±2} |
| 5gσg | (0,4) | 20.000 | 19.926 | H n=3 | +1.9880 | {0,±2} |
| 2pπu | (0,0) | 2.0000 | 1.840 | H n=2 | +0.0000 | {0} |
| 3dπg | (0,1) | 6.0000 | 5.823 | H n=2 | +0.0000 | {0} |
| 3pπu | (1,0) | 2.0000 | 1.926 | H n=3 | −0.9971 | {±1} |
| 4fπu | (0,2) | 12.000 | 11.893 | H n=3 | +0.9953 | {±1} |

**Z_A = 3, Z_B = 1** (A = the Z = 3 centre)

| UA (Z=4) | (n_ξ,n_η) | λ(1.91) | SA | k from Λ | A_z spec |
|:--|:--|:--|:--|:--|:--|
| 1sσ | (0,0) | −7.296 | A n=1 | +0.0004 | {0} |
| 2pσ | (0,1) | 1.179 | A n=2 | +1.0000 | {±1} |
| 2sσ | (1,0) | −3.052 | A n=2 | −1.0000 | {±1} |
| 3pσ | (1,1) | 1.987 | A n=3 | −0.0000 | {0,±2} |
| 3dσ | (0,2) | 5.300 | **B (H) n=1** | −0.0000 | {0} |
| 3sσ | (2,0) | −2.270 | A n=3 | −1.9999 | {0,±2} |
| 4fσ | (0,3) | 11.676 | A n=3 | +1.9999 | {0,±2} |
| 2pπ | (0,0) | 0.522 | A n=2 | +0.0000 | {0} |
| 3dπ | (0,1) | 5.507 | A n=3 | +1.0000 | {±1} |
| 3pπ | (1,0) | 0.960 | A n=3 | −1.0000 | {±1} |

Note the Z = 3 n = 3 level is exactly degenerate with H(1s) at −1/2 Ha. The 3dσ → H(1s) and 4fσ → Li²⁺(n=3, k=+2) states have the same leading slope s = 3. Λ separates them only at O(1) in c₀, and ⟨η⟩ assigns the centre.

### 2.5 Crossings

**All detected same-symmetry crossings are true crossings between states whose Λ differs by O(1–40).** This includes crossings between states that share m and inversion parity, in both systems. Examples:

- H2+ 2sσg × 3dσg at R ≈ 4.06: Λ = −0.88 vs 4.82.
- (3,1) 3pσ × 3dσ at R ≈ 1.16: Λ = 2.05 vs 5.79.

Among the 8 lowest states per block, H2+ has 46 crossings over R ∈ [0.1, 80] and (3,1) has 7. Three R < 0.04 hits are united-atom degeneracies and are excluded. No avoided crossing is possible between tracks, since the labels are Klein indices. So the Wigner–von Neumann non-crossing rule is lifted exactly where Λ differs, which is everywhere.

**GATE A: PASS.** Commutation is exact at every truncation, and both limits are reproduced with a quantitative Runge–Lenz check.

## 3. Stage B: H2 (production Paper 12 CI, (j,l,μ) = (4,4,2), α = 1.2, 975 functions)

Driver: `debug/prolate_separation_constant_stageB.py`. Data: `debug/data/prolate_separation_constant_stageB_442_a1.2.json` and `debug/data/stageB_442.out`.

**Setup.**

- The CI S, H1 and V_ee come unchanged from `geovac.prolate_recondition` (direct engine).
- The one-electron W, Hm and N per μ are assembled from that module's own exact mpf 1D blocks.
- Consistency checks: the product-basis S block reproduces h⁶·cc·W⊗W bit-exactly; the MO-pair basis is orthonormal to 4×10⁻⁹; E₀ = −1.17405771 in the MO basis equals the recondition solve.

**One-body.** Within each μ block, [Hm, Λ] = 0 to 10⁻¹⁵–3×10⁻¹² (Atkinson). The MOs are the joint H2+ eigenfunctions, parity impurity < 10⁻²⁶.

**Two-body.**

- Relative ‖[H, Λ₁]‖ = 9.1×10⁻³. The one-body part alone gives 4.5×10⁻¹¹, so all of it is V_ee.
- ‖[H, Λ₁]ψ₀‖ = 2.86 Ha (in Λ units).
- ⟨Λ₁+Λ₂⟩ = −0.795, with standard deviation 0.651.

For scale, 1σg has Λ = −0.4445, so 1σg² gives −0.889.

**Configuration weights (MO pairs).**

| Configuration | Weight |
|:--|:--|
| 1σg² | **92.27%** |
| 1σg·σg(Λ = −0.129) | 5.1% |
| 1σg·σg(Λ = 0.024) | 0.9% |
| 1σu² (Λ = 1.637) | 0.8% |
| 1πu² | 0.19% |

**Truncation** (ΔE in mHa relative to the full CI of 975 configurations):

| kept | Λ-sum | Λ-max | orbital-energy-sum | random (+HF) |
|:--|:--|:--|:--|:--|
| 50% | 0.210 | 0.134 | 1.953 | 38.4 |
| 30% | 0.537 | 0.317 | 3.996 | 60.0 |
| **20%** | **0.869** | **0.551** | 6.704 | 68.7 |
| 10% | 2.640 | 2.399 | 8.140 | 81.2 |

**Orbital-level truncation** at equal orbital count, compared with **natural orbitals** (the strongest baseline, but it needs the answer first):

| n_orb | Λ-ordered | energy-ordered | natural orbitals |
|:--|:--|:--|:--|
| 10 | 8.4 mHa | 28.7 mHa | 0.98 mHa |
| 30 | 0.55 mHa | 10.7 mHa | 0.08 mHa |

**GATE B: GO, on both criteria.**

- The dominant configuration has 92.3% of the weight (> 90%).
- A Λ-sum-truncated CI loses 0.87 mHa while dropping 80% of configurations (< 1 mHa, ≥ 30% dropped).

**Honest reading.**

1. The 92% weight is the ordinary fact that H2+ MOs are good H2 orbitals. It is not specific to Λ.
2. The substantive finding is that **an a priori, one-electron, exactly computable label ranks configurations 8–10× better than orbital energy** at 20–50% retention. My reading of why: virtual-orbital energies in a finite basis are a poor correlation ranking, while small Λ selects low nodal/angular content in η.
3. Λ loses to natural orbitals by about 10× at orbital level. So Λ is a good *a priori* selector, not the best selector.
4. Λ₁+Λ₂ is not a conserved label for H2 (σ = 0.65).

## 4. Stage C: LiH (Route C orbital set)

Driver: `debug/prolate_separation_constant_stageC.py`. Data: `debug/data/stageC_{build,analyse}.out` and the cached integrals `prolate_separation_constant_stageC_ints.npz`.

**Setup.**

- Orbital set: `lih_orbitals_c4(R=3.015, extended, 1π)`, 10 orbitals, 2025 determinants.
- Engines unchanged. The FCI reproduces the logged C4 energy, **−8.00794 Ha**.
- The Λ(Z_A=3, Z_B=1) matrix is built with the engine's moment helpers. The same pieces re-assemble the engine's h1 exactly (relative deviation 0.0).

**Findings.**

1. **One-body commutation is lost.** ‖[h1, Λ]‖_rel = 2.2×10⁻², because the Route C set (STO cores with η-exponential factors, mixed exponents) is not a product basis, and Atkinson's structure needs one. MOs of the bare Li³⁺+H⁺ h1 carry large Λ leakage: ‖(Λ−⟨Λ⟩)φ‖ = 1.8–4.7 for all σ MOs except the Li 1s core (0.04). The π pair has zero leakage.
2. The leading CSF has 59% of the weight on these bare-nucleus MOs (< 90%).
3. **Truncation is no better than energy ordering.**

| kept | Λ-sum | orbital-energy-sum | random |
|:--|:--|:--|:--|
| 50% | 2.01 mHa | 1.97 mHa | 138 mHa |
| 20% | 7.5 mHa | 10.9 mHa | 275 mHa |
| 10% | 47.6 mHa | 34.0 mHa | 302 mHa |

   At orbital level with n_orb = 8, Λ-ordering gives 12.0 mHa, ε-ordering 37.9 mHa, and natural orbitals 0.41 mHa.

**GATE C: STOP.** Neither criterion is met: the dominant weight is < 90%, and a 20%-kept truncation loses 7.5 mHa.

**Reasons.** There are two:

- The Route C set is too small and not product-structured for Λ to be represented faithfully.
- The bare-nucleus Λ ignores core screening, so the valence MOs are not Λ-eigenfunctions of any relevant one-body operator.

**What would change it** (not attempted; these are engine-level decisions for the PI):

- a product-structured (ξ-Laguerre × η-Legendre) valence set;
- a screened Λ, built from the core-shielded one-body operator.

## 5. New versus textbook (honest)

**Textbook, verified against primary-source abstracts.**

- Λ is the Erikson–Hill constant of motion. H. A. Erikson & E. L. Hill, *Phys. Rev.* **75**, 29 (1949): H, L_z and this operator form a complete commuting set, giving the dynamical meaning of the separation parameter.
- Its derivation from the Runge–Lenz vector, which is the "L² ↔ A_z" reading, is C. A. Coulson & A. Joseph, *Int. J. Quantum Chem.* **1**, 337 (1967), "A constant of the motion for the two-centre Kepler problem".
- The existence of same-symmetry crossings and the united-atom/separated-atom correlation diagrams are standard H2+ lore. I did not check a specific primary source for the correlation table; mark it [unverified attribution].

**So stage A's physics is 1949/1967.**

**What this sprint adds (modest, and not checked for prior statement):**

1. **Finite-truncation exactness via Atkinson's operator determinants.** In any product prolate basis, (W, Hm, N) = (Δ₀, Δ₁, Δ₂), so [H, Λ] = 0 and separability (rank-1 joint eigenvectors) hold at every cutoff, to rounding. This is the prolate analogue of the GeoVac pattern "the Casimir is exact on the labels at every n_max." Multiparameter spectral methods for spheroidal problems exist, so this may already be stated somewhere; I did not search.
2. **Quantitative SO(4) reading of the separated-atom limit.** (λ+p²)/R − Z_other lands on the spectrum of the corpus's own Runge–Lenz A_z (Z/n)k construction to ≤ 8×10⁻³ ((3,1), n ≤ 6), with sign agreement against first-order Stark. The operator identity is Coulson–Joseph; the numerical check against `so4_generators` is ours.
3. **Stage B.** The one-electron Λ is a good *a priori* configuration selector for H2: 8–10× better than orbital-energy ordering, and still about 10× worse than natural orbitals. Λ₁+Λ₂ is not conserved (the V_ee commutator is fully O(1)). This is the only piece beyond textbook. Its scope is one molecule, one basis and one geometry.
4. **Stage C.** The selector does not transfer to the Route C LiH set. Its product structure is where the exactness came from.

## 6. Test

`tests/test_paper11_separation_constant.py`: 14 tests, 3.5 s. It is self-contained and imports `so4_generators` from the S³ test.

Checks:

- **(a)** Commutation < 1e-11 for (3,1) and (1,1), at R = 0.5, 2, 8 and m = 0, 1, plus a non-triviality check that ptp(Λ) > 1.
- **(b)** United-atom λ → l(l+1) and E → −Z²/2n², with (n, l) read from the Klein nodal labels. States with s, p and d character are required.
- **(c)** Separated-atom k within 0.02 of the A_z spectrum for the lowest four (3,1) states, of which at least two have n > 1.

**Fire-tested** with `debug/qa/fire_test.py`; each plant made its guard FAIL:

| Plant | Guard that fired | Wrong answer rejected |
|:--|:--|:--|
| Λ with the wrong sign on the Z_A−Z_B term (`+R·Zm·a1 → −R·Zm·a1`) | (a) | A sign-flipped Λ |
| Λ → 2Λ (still commutes) | (b) | A mis-normalized Λ |
| Λ → Λ + (R³/40)·H (commutes, and passes (a) and (b)) | (c) | A Λ that commutes and has the right united-atom limit but the wrong separated-atom slope |

Control: the third plant stays green on (a) and (b), as designed, so (c) is the only check that rejects it.

**Wrong answer the test does *not* reject:** Λ + f(H) with f vanishing faster than R at both ends.

## 7. Files

`debug/prolate_separation_constant{,_xcheck,_stageA,_stageB,_stageC}.py`; data in `debug/data/` (stageA json/log, stageB_442 json/out, stageC out + npz); `tests/test_paper11_separation_constant.py`. No geovac/, paper, CLAUDE.md, CHANGELOG or QA registry was touched.
