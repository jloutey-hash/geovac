# Scoping memo: what the prolate product basis gains and gives up (2026-09-26)

**Type:** scoping / diagnostic. No engine code was written, and geovac/, the papers, CLAUDE.md, CHANGELOG and the QA registries are untouched.

**Probes:** three cheap float64 scripts (about 1 minute each):

- `debug/product_basis_tradeoff_probe.py`
- `debug/product_basis_tradeoff_probe_zeff.py`
- `debug/product_basis_tradeoff_probe_zeff2.py`

Their logs are in `debug/data/product_basis_tradeoff_probe*.log` and `.out`.

**Guardrails:** everything stays two-centre prolate (Level 2) with per-centre exponents. Papers 8–9 are not triggered: there is no single centre and no shared p₀. FCI-M is not triggered: there is no graph concatenation.

## 0. Two corrections to the framing

1. **Route C's orbitals are all product functions. What they lack is tensor closure.**
   - Every Route C orbital is rank-1 in (ξ, η). Take the Li 1s STO: since r_A = (R/2)(ξ+η),
     e^{−ζ r_A} = e^{−(ζR/2)(ξ−1)} · e^{−(ζR/2)(1+η)}, exactly.
   - The same holds for STO powers, because r_A^k is a polynomial in (ξ+η), and for p-functions, because r_A cos θ_A = (R/2)(ξη+1).
   - Atkinson's exactness needs something stronger than "each function is a product". It needs the Galerkin space to be a **full tensor product V_ξ ⊗ V_η**, so that W, H_m and N are Kronecker combinations. The generalization with a non-orthonormal η-overlap a₀ is W = r₂⊗a₀ − r₀⊗a₂; the determinant identities still hold.
   - The Route C set is a union of 8 σ + 2 π rank-1 functions with distinct ξ-rates and η-factors, not a closed tensor space. That is the whole origin of stage C's ‖[h₁,Λ]‖ = 2.2×10⁻².
   - A corollary: a *triangular* truncation (j + l ≤ N) would also break exactness. Paper 12's rectangular (j, l) truncation keeps it.
2. **In an analytic product basis the tight core is cheap (probe P1).** The "0.36 Ha wall" was not a property of product structure. It came from grid quadrature and grid-orbital MOs (increment 7).

## 1. The ledger

| # | Property | Depends on product (tensor) structure? | Evidence (corpus + this probe) |
|:-:|:--|:--|:--|
| G1 | **Λ commutes with h exactly at every cutoff** | **YES.** Needs a full tensor space *and* a separable one-body operator: bare Coulomb, −Z_A/r_A − Z_B/r_B plus a constant. | Stage A/B: 10⁻¹⁶–10⁻¹². Route C (not closed): 2.2×10⁻². P2 below: any core screening breaks it (10⁻³ relative) even inside a closed tensor space. |
| G2 | **x⊗y decomposable eigenvectors** (rank-1 MOs) | **YES**, for joint eigenfunctions of a separable h | Stage A: σ₂/σ₁ ≤ 10⁻¹¹. P2: screened σ-valence eigenvectors have σ₂/σ₁ = 0.04–0.6. π stays at 0.03. |
| G3 | **Λ as an a-priori CI selector** | **YES** (it inherits G1), and it also needs the separable h's MOs to be good orbitals | H₂: 20% of configurations kept for 0.87 mHa, 8–10× better than ε-ordering, but about 10× worse than NOs. LiH stage C: 7.5 vs 10.9 mHa at 20%, marginal; leading CSF 59%. |
| G4 | **Closed-form / algebraic integrals** | **NO.** This depends on the *form* of the factors (polynomial-η, exp·poly-ξ), not on tensor closure. | Pure polynomial-η product (Paper 12): the Neumann l-sum *terminates*, and the B-seeds are closed-form {E₁, γ, ln}. Route C one-body stays closed-form (incomplete Γ, e^{±β}; C1 exact to 1e-50). Route C core ERIs use the Legendre-η series of e^{−βη}, truncated at L = 24 (i_l(β) seeds, converged to 1e-15/1e-28). So Route C keeps an exact, deterministic evaluation, but *convergent* rather than *terminating*. A tensor-closed Route C set would have the same status. |
| G5 | **Conditioning** | **PARTLY.** Product structure lets conditioning be repaired by exact 1D re-basing (Kronecker-factored). | Paper 12: monomial cond 1.8e16 (3,3) → 4.1e26 (5,5), cured by Laguerre×Legendre re-basing (~10⁶ gain). The same trick on the Route C ladder gives 2.9e9 → 1.4e5, with LiH ≤ 3.6e6. Pure Laguerre×Legendre 30×40: cond(W) = **6.0e3** at R = 2.6–3.45 (P2). Increment 7's mixed grid MOs: 1e6. Paper 12 warns that raw and normalized cond differ by up to 6 orders. |
| G6 | **ERI sparsity / selection rules** | **NO extra.** Only M_L (plus g/u, homonuclear only) holds either way; Λ produces no ERI zeros. | V_ee–Λ commutator is O(1) (H₂: 9.1×10⁻³ relative; ‖[H,Λ₁]ψ₀‖ = 2.86). Route C M = 10 ERI tensor is **47.8% nonzero** (M_L selection only). Within-m Löwdin preserves M_L, so the Gaunt-sparsity dead end (Löwdin 14× inflation) does not bite in prolate: there is no l-rule to lose. |
| G7 | **Qubit Pauli counts** | **NO** (follows from G6) | The corpus has **no prolate Pauli count at all** (grep over papers: none). Predicted to be dense within M_L blocks either way. |
| G8 | **Factored MO ERIs** (a payoff of G2) | **YES**, for rank-1 MOs | Neumann V_ee = Σ_l X_l(ξ₁,ξ₂)Y_l(η₁)Y_l(η₂), with the Jacobian adding four separable terms per l. For rank-1 MOs, (pq\|rs) contracts in 1D with no AO→MO transform. Not measured. It is useless for LiH σ valence (G2 is lost). |

**What product structure costs** (the costs claimed in the task, re-measured):

| # | Claimed cost | Verdict | Numbers |
|:-:|:--|:--|:--|
| C1 | "A tight Li 1s needs many ξη products" | **REFUTED as stated.** It is *one* product. In a single-α Laguerre×Legendre basis the cost is primitives, not rank. | Bare Z = 3 1s (exact −4.5 Ha, R-independent), core e^{−2.6875 r_A}. Errors and core norm outside the span: 10×10, α = R/2: **+1.2e-4 Ha**, outside 2e-6. 20×20: ≤ 1.6e-6 Ha for any α. 30×40: ≤ 2.4e-10 Ha, outside ≤ 1e-12. |
| C2 | "Resolution is R-dependent" | **TRUE only in small bases with a valence-tuned α** | 10×10, α = 1.2: error +0.4 → +7.2 mHa over R = 2.6 → 3.45 (R-spread **6.8 mHa**, 4× chemical accuracy). α ∝ R removes it (10×10, α = R/2: spread 1.2e-5 Ha). The 0.36 Ha increment-7 swing was the grid; C2 of Route C measured the grid ERI swing (26 mHa for (cc\|cc)), which an analytic basis also avoids. |
| C3 | "Basis size for heteronuclear systems" | **TRUE, and this is the real cost** | Primitive tensor 10×10 to 20×20 per m (100–400 functions) vs Route C's 8 σ orbitals. Orbital FCI caps at M ~ 16–20 (sparse FCI, v5.15.19), so the primitives must be **contracted**, and the contraction choice is the crux (§2). |

**Net.** Product structure buys the exact label (G1–G3) and the 1D conditioning repair (G5). It costs *only* primitive count (C3); tight cores and R-flat resolution are affordable analytically (C1, C2). It buys **nothing** for integral algebraicity (G4), sparsity (G6) or Pauli counts (G7).

## 2. The crux: the orbitals that carry the label are not LiH's good valence orbitals

After contraction, only eigenfunctions of a *separable* one-body operator keep an exact Λ, and for LiH those are poor σ-valence orbitals (stage C: bare Li³⁺+H⁺ MOs give a 59% leading CSF). Good orbitals come from a *screened* operator, which is non-separable.

**P2 — quantified.** Setup:

- h_eff = T − 1/r_A − 1/r_B (separable) + δV_H(r_A) + σ|c⟩⟨c|.
- δV_H = −2e^{−2ζ r_A}(ζ + 1/r_A) is the Hartree residual of Li 1s², a short-range, non-separable piece. Its matrix is a sum of four Kronecker products, computed exactly.
- σ|c⟩⟨c| is the core level-shift projector; c is itself a product function.
- Core exchange is omitted. The electron is a one-electron valence model.
- Basis: pure 30×40 tensor.

Results:

- ‖[h_eff, Λ(1,1)]‖_rel rises from 3.7e-15 (separable only) to **1.7e-3**. The level-shift projector adds nothing: its ⟨P⟩ on the valence is 3e-4. **The breaker is core penetration (δV)**, not orthogonality bookkeeping.
- **σ-valence label purity** is max |⟨joint (E,Λ) eigenstate of h(Z_eff,1) | ψ⟩|² over R = 2.6–3.45:

| Z_eff of the label | lowest σ valence | next σ | π lowest (m = 1) |
|:--|:--|:--|:--|
| 1.0 (asymptotic, a-priori) | 0.59–0.61 | 0.30–0.58 | **1.000** |
| 1.5 | 0.89–0.96 | 0.59–0.85 | 0.92–0.94 |
| 2.0 | **0.96** | 0.56–0.79 | 0.72–0.79 |
| 2.5 | 0.85–0.91 | 0.44–0.63 | 0.51–0.62 |

Reading:

1. **π is label-clean.** It does not penetrate the s core, so Λ(1,1) is essentially exact on screened π (purity 1.000, leakage ~0.2).
2. **The lowest σ valence is 96% one Λ(2,1) state**, so the separable part *does* dominate the occupied σ.
3. **The correlating σ virtuals are 56–85% pure**, and those are exactly what a CI selector must rank.
4. The best Z_eff (≈ 1.5–2) is chosen *a posteriori* by purity. It is a label-operator choice, not a Hamiltonian parameter, so it breaks no §13.5 rule, but it is not a-priori.
5. **Mechanism.** Short-range penetration makes the effective charge run from 3 to 1 across the valence, and Λ(Z,1) can hold only one constant charge. A radial node spherical about A is not a prolate coordinate surface.

## 3. Middle paths, and whether each re-enters a dead end

| Path | Keeps an exact label? | Dead-end check |
|:--|:--|:--|
| (a) Product valence set on top of a projected/orthogonalized core | No. Exact only if the core is absent from the valence operator; P2 shows penetration (δV), not the projector, is what breaks it. | Using the projected core **as the Hamiltonian** re-enters the **frozen-core wall** (+5.5% outward and worsening, v5.15.6; increments 2b/6). It is acceptable only as a device to *define orbitals and labels*, with all four electrons kept active in the CI (Route C). |
| (b) Screened Λ(Z_eff, 1) | Approximate: occupied σ 96%, π ~100%, σ virtuals 56–85% | New, and not a dead end. It *is* the "effective-charge label"; do not promote it to a Hamiltonian Z_eff (Paper 17's Z_eff/PK ceiling lives there). |
| (c) Tensor-closed per-centre exponent families (Route C's factors closed into V_ξ⊗V_η) | Exact for **bare** h. Useless for screened orbitals (§2). | Not Papers 8–9 (two-centre, per-centre rates). **Two-block per-block exponent** (Paper 12 H₂) showed that extra ξ-rates buy ~0 accuracy on H₂, but LiH's core rate is load-bearing (the 2nd core exponent gives −17 to −23 mHa). The Löwdin/Gaunt-sparsity row does not apply (G6). |
| (d) Pure Laguerre×Legendre all-electron LiH (no atom-centred STO), with α ∝ R or a 20×20 primitive space | Exact on bare-h MOs | P1 says the one-body core is clean. Untested: the (cc\|cc)-class ERIs through the terminating Neumann at a high polynomial degree, and the contraction problem (C3). Not a known dead end, but it lands on the §2 crux: its label-carrying MOs are the bad ones. |

## 4. The payoff if the label were restored for LiH, honestly

- **Energy.** None expected. The LiH energy is limited by correlation, not orbital selection:
  - the orbital lever is near-exhausted: M = 16 → 20 buys 3 mHa, with a ceiling of −8.032;
  - the gap to exact is the cusp plus three-body terms (VMC −8.047, marriage −8.037, the −8.062 additive estimate).
  A better CI selector can only reach the same ceiling more cheaply.
- **Compression / active space.**
  - The production marriage already truncates by **natural orbitals** (74 determinants).
  - On H₂, NOs beat Λ by ~10× (0.98 vs 8.4 mHa at 10 orbitals).
  - On LiH, Λ beat ε-ordering only marginally (7.5 vs 10.9 mHa) and lost to NOs by ~30× (0.41 vs 12.0 mHa at n_orb = 8).
  - A screened Λ could at best reach H₂-like performance on σ, still behind NOs. A cheap-NO proxy (MP2/CISD NOs) is the fair baseline and is standard practice.
- **Qubits / Pauli.** No ERI zeros (G6) means no Pauli reduction. Qubit count falls only through active-space size, where NOs dominate.
- **Where it could still pay.**
  1. π/δ polarization-shell selection: the label survives screening there, so Λ could order which π/δ product functions to add to the ladder a priori.
  2. Atlas bookkeeping: Λ is **FORCED** (exact at every cutoff) for one-electron two-centre problems and for bare-nucleus orbital sets. Core screening is precisely where it becomes approximate, and that boundary is now measured (§2).

**Verdict on an engine change: not worth it** for LiH energy or qubits. A tensor-closed or pure-product LiH engine would restore an exact label on the *wrong* orbitals.

## 5. Ranked analyses feasible in one session

| Rank | Analysis | Measures | Cost | Decision gate |
|:-:|:--|:--|:--|:--|
| **1** | **Screened-Λ selector on the cached Route C FCI.** Reuse the stage C `analyse` path and `prolate_separation_constant_stageC_ints.npz`. Build MOs from the core-screened Fock-like h (h₁ + 2J_c − K_c from the cached ERIs) and Λ(Z_eff,1) for Z_eff ∈ {1, 1.5, 2} via `lambda_matrix`. Add a cheap-NO baseline (NOs of a CISD in the same 10 orbitals). | Leading-CSF weight, and the truncation table (Λ-sum vs ε vs random vs CISD-NO vs exact NO) | ~1 h to write; minutes to run (10 orbitals, 2025 determinants, integrals cached) | **GO** if 20%-kept loses < 1 mHa AND beats ε by ≥ 3× AND is within 3× of CISD-NOs. **STOP** otherwise. *Prediction from P2: STOP/BORDERLINE*, because σ virtuals are only 56–85% pure. |
| 2 | **Tensor closure of the Route C one-body set.** V_ξ = span of the distinct ξ-factors, V_η = span of the distinct η-polynomials, one-body only (c4 `one_body_general` + `lambda_matrix`, mpf). | ‖[h₁, Λ(3,1)]‖ (theory: back to ~1e-14), dimension growth, normalized cond(S), bare one-body spectrum vs the unclosed set | ~1–2 h; a ~50–100-function mpf one-body build, estimated minutes | Confirms the §0 correction numerically. **GO** (closure is cheap) if dim ≤ 5× and normalized cond < 1e8. It does *not* change the §4 verdict, so it is a bookkeeping result for Paper 11 §separation_label (the "Route C loses it" sentence should say *not tensor-closed*, not *not product*). |
| 3 | **P2 extension: add core exchange −K_c (non-local).** Two-electron valence in h_eff via a stage-B-style product CI. | Whether exchange raises or lowers the σ purity; the 2e leading weight with screened-Λ MOs | ~half session. Needs a Neumann (cv\|cv) exchange on the product basis (C2 `build_Xtab_pair` exists but is mpf-slow at 1200 primitives; a smaller 12×16 tensor is feasible) | Informational. Promote to rank 1's inputs only if purity of the σ virtuals rises above 0.9. |
| 4 | **Pure-product core ERI check.** (cc\|cc) for the contracted Z = 3 1s from a 20×20 Laguerre×Legendre expansion through the terminating Neumann (`prolate_recondition` X-table, high degree), vs the closed form 5ζ/8, at R = 2.6/3.0/3.45. | Whether the atom-centred STO in Route C is *needed*, or only convenient | ~half session (mpf, one pair density of degree ~40 in ξ and η) | GO if the error and R-spread are < 1e-6 Ha. It unlocks option (d) but hits §2, so it is low value. |
| — | **Prolate Pauli count** (fills the G7 gap; coverage, not a product-basis question) | N_Pauli for Route C M = 10 via JW on the cached h₁/ERI | ~1 h | Predicted dense within M_L; PI call for Paper 14 completeness. |

## 6. Recommendation

1. **Do not change the engine.** Product structure buys an exact label and cheap conditioning repair. It does not buy algebraicity, sparsity or Pauli terms. For LiH the label lives on the wrong orbitals: penetration of the Li core makes σ valence and virtual orbitals only 56–96% label-pure, and the energy is limited by correlation, not selection.
2. **Run analysis 1** (≤ 1 h plus minutes). It is the one cheap, decisive test of whether any Λ variant is useful for LiH. P2 predicts STOP. If it lands GO, rank 3 is the follow-on.
3. **Optionally run analysis 2**, and carry its §0 correction into Paper 11 §separation_label. That paragraph attributes the loss to a "non-product basis". The accurate statement is **"not tensor-closed"** (every Route C orbital is itself a product), and the scope sentence should also say that screening, not only basis structure, breaks the label. This is a PM paper edit, not done here.
4. **Record, for the atlas:**
   - The "tight core needs many products / R-dependent resolution" cost is **refuted in an analytic basis** (P1: 10×10 with α = R/2 gives a 1.2e-5 Ha R-spread; 20×20 ≤ 1.6e-6). The 0.36 Ha wall was the grid.
   - The boundary where Λ goes from exact to approximate is **core penetration**, measured at 1.7e-3 relative commutator, σ purity 0.6 (Z_eff = 1) to 0.96 (best), π ≈ 1.

**Caveats.** The P1/P2 probes are one-body. The valence model is a single electron with Hartree screening only (no core exchange, no valence–valence repulsion). The Z_eff scan picks its optimum a posteriori. The stage-C comparison numbers are quoted from `sprint_prolate_separation_constant_memo.md`, not re-run.
