# Sprint L2 — Synthesis Memo (Lorentzian extension closure at finite cutoff)

**Sprint:** L2 (Lorentzian extension of the WH1-PROVEN Riemannian foundation via Krein-space spectral triple at signature (3, 1))
**Dates:** 2026-05-16 (L2-A through L2-E proper) — 2026-05-17 (L2-F/G synthesis + paper edits)
**Verdict:** **CLOSED-AT-FINITE-CUTOFF.** All five sub-sprints (L2-A scoping, L2-B Krein space, L2-C Lorentzian Dirac, L2-D Connes axiom audit, L2-E Krein-level modular Hamiltonian) closed at finite cutoff via bit-exact Riemannian-limit checks and bit-exact load-bearing falsifier passes. WH1 PROVEN survives all load-bearing checks.

**Source memos (in chronological order):**

- `debug/sprint_l2a_scoping_memo.md` — GO-WITH-CAVEATS, 3.5 months scaling
- `debug/sprint_l2_falsifiers.md` — 17 falsifiers catalogued, 4 load-bearing
- `debug/l2_b_krein_construction_memo.md` — Krein space + fundamental symmetry
- `debug/l2_c_lorentzian_dirac_memo.md` — Lorentzian Dirac via vdD Prop 4.1
- `debug/l2_d_connes_axiom_audit_31_memo.md` — BBB axioms at (m, n) = (4, 6)
- `debug/l2_e_modular_hamiltonian_lorentzian_memo.md` — Krein-level modular Hamiltonian
- `debug/l2_paper_edits_round1_memo.md` — L2-B/C/D paper edits applied
- `debug/l2_g_synthesis_memo.md` — this synthesis sprint's summary

---

## §1. Executive summary

Sprint L2 was launched in response to a named literature trigger: Nieuviarts 2025 (arXiv:2502.18105) presents a twisted-spectral-triple morphism Riemannian → pseudo-Riemannian that, in principle, shortens the multi-month Lorentzian-propinquity construction that the May-2026 Lorentz-boost NO-GO scoping memo identified as the dominant cost. Sprint L2-A audit verified GO at a 3.5-month-equivalent scope; the actual execution collapsed to a two-session sprint (2026-05-16 through 2026-05-17) via parallel sub-agent dispatch, with all five sub-sprints landing bit-exact at finite cutoff. **The full Lorentzian-extension arc closes at the operator-system level on the Lorentzian Krein space at finite n_max × finite N_t.** The continuum-limit / Lorentzian-propinquity-rate-level extension (Sprint L3) remains open as a separate frontier; no published Lorentzian propinquity construction exists as of May 2026.

The substantive new structural finding (the headline of the sprint) is that **Paper 42 §7.2 open question O3 (the H_local ≠ D_W finding) is signature-INDEPENDENT at the Riemannian limit**. The spectral-action--vs--modular-Hamiltonian generator distinction is NOT a peculiarity of the round-S³ Riemannian truncation; it is a deeper structural feature of the framework's modular content. At the Riemannian limit on the Lorentzian Krein space, the residual ‖H_local − D_L^W‖_F equals the Riemannian-side residual bit-exact at every tested n_max. At N_t > 1 the gap widens (refined upward by temporal-derivative content), but the Riemannian-limit baseline is preserved.

WH1 PROVEN (Paper 38, qualitative-rate Latrémolière-propinquity GH-convergence theorem on the Riemannian truncated triple) survives all load-bearing checks. The Lorentzian extension is structurally additive on top of the WH1 Riemannian foundation, not a re-test of it.

---

## §2. The L2-A → L2-B → L2-C → L2-D → L2-E cascade

The Lorentzian-extension architecture follows the van den Dungen 2016 Proposition 4.1 lift recipe verbatim, with the Camporesi–Higuchi spatial spinor bundle on S³ in the Peskin–Schroeder chiral basis at signature (3, 1) West-coast:

### §2.1 L2-A scoping (1–2 days, 2026-05-16 morning)

**Verdict:** GO-WITH-CAVEATS. 3.5-month-equivalent scope decomposes into five sub-sprints (B/C/D/E proper + F catalogue + G synthesis). Three named blockers: B1 (Lorentzian propinquity, not addressed by L2), B2 (Krein-space spectral triple at (m, n) = (4, 6), addressed by L2-B + L2-D), B3 (Lorentzian-native heat-trace formalism, not addressed by L2; the Wick-rotation temporal-circle special case is sufficient for Sprint TD Track 4 and L2-E). Three named risks (R1 spinor-bundle embedding ambiguity, R2 temporal BC choice, R3 continuous spectrum) all resolved by structural requirement in the sub-sprints.

### §2.2 L2-B Krein space (2026-05-16 morning)

**Verdict:** CLOSED at every tested cutoff. Construction:

- **Hilbert space carrier:** K_{n_max, N_t} = H_GV^{n_max} ⊗ L²(ℝ_t)_cutoff with bounded uniform temporal grid t ∈ [−T_max, T_max].
- **Fundamental symmetry:** J = γ⁰ ⊗ I_{N_t} where γ⁰ is the Cl(3, 1) time gamma matrix in the chiral basis (a real ±1-entry off-diagonal swap on the chirality-doubled H_GV).
- **Convention lock-in:** West-coast metric η = diag(+, −, −, −), Peskin–Schroeder chiral basis where γ⁵ is diagonal (making FullDiracLabel.chirality directly the γ⁵ eigenvalue label). The Dirac-basis alternative is rejected on Riemannian-limit clarity grounds.
- **NOT compactified:** the temporal slot is the bounded grid, NOT the periodic circle S¹_β. Periodic BC would create CTCs on S³ × S¹ per Geroch's theorem, structurally incompatible with the bounded-wedge thermal-physics reading.

**Bit-exact axioms** (Frobenius residual = 0.0 in float64) at every (n_max, N_t) ∈ {1, 2, 3} × {1, 11, 21}: J² = +I, J*J = I, J* = J, K = K^+ ⊕ K^− with dim K^± = dim K/2.

**Load-bearing falsifier L2B-FALS-4** (Riemannian-limit recovery at N_t = 1): PASS bit-identically. dim K_{N_t=1} = N_Dirac(n_max) = (2/3)·n_max·(n_max+1)·(n_max+2) = 4, 16, 40 at n_max = 1, 2, 3 — matches Paper 32 Def 3.2 (`def:H_GV`) exactly; basis labels via FullDiracLabel.__eq__ identical object-by-object; J matrix Frobenius residual = 0.0. **Confirms the Camporesi-Higuchi spatial spinor bundle is structurally compatible with the Cl(3, 1) gamma-matrix embedding in the chiral basis.** WH1 PROVEN is NOT re-opened.

### §2.3 L2-C Lorentzian Dirac (2026-05-16 midday)

**Verdict:** CLOSED-WITH-STRUCTURAL-FINDING. Construction via van den Dungen 2016 Proposition 4.1 (arXiv:1505.01939) instantiated on (M, g) = (S³ × ℝ, ds²_S³ − dt²) at signature (s, t) = (3, 1):

$$D_L = i \cdot \big[\gamma^0 \otimes \partial_t + D_{\rm GV} \otimes I_{N_t}\big]$$

- **∂_t centered FD + Dirichlet zero BC is structurally forced**, not chosen — the only choice giving anti-Hermitian ∂_t on a finite grid (forward/backward/upwind break Krein-self-adjointness; periodic re-creates CTCs).
- **The sign i^t = +i for t = 1 is derived structurally** from the D_L^× = D_L requirement (sympy-exact via {γ⁰, D_GV} = 0 and ∂_t^† = -∂_t); the alternative i^t = −i would give anti-Krein-self-adjoint D_L.

**Three load-bearing falsifiers pass bit-exact** at every (n_max, N_t) ∈ {1, 2, 3} × {1, 11, 21}:
- (a) Krein-self-adjointness D_L^× = γ⁰ D_L^† γ⁰ = D_L (residual = 0.0);
- (b) Riemannian-limit recovery D_L|_{N_t=1} = i·D_GV bit-identically (load-bearing L2C-FALS-3); spectrum |λ_n| = n+1/2 with multiplicities g_n^Dirac = 2n(n+1) matches D_GV exactly;
- (c) real spectrum on Krein-positive cone K^+ at machine precision (|Im(λ)| ≤ 2.32×10⁻¹⁷).

**Structural finding L2C-FALS-2 (chirality anticommutation):** {γ⁵, D_L} ≠ 0 in general. Decomposition: {γ⁵, D_L} = 2i (γ⁵ D_GV) ⊗ I_{N_t}, with ‖{γ⁵, D_L}‖_F = 2 ‖D_GV‖_F · √N_t. This is NOT a bug — it is a structural consequence of GeoVac's chirality-diagonal D_GV in the chiral basis (truthful CH commutes with γ⁵, not anticommutes). The diagnostic feeds into Sprint L2-D's BBB axiom audit at (m, n) = (4, 6).

### §2.4 L2-D Connes axiom audit at BBB (m, n) = (4, 6) (2026-05-16 afternoon)

**Verdict:** CLOSED-WITH-STRUCTURAL-FINDING. BBB 2018 Table 1 signs at (m, n) = (4, 6) verified directly from PDF (arXiv:1611.07062 v2): ε = +1, ε'' = −1, κ = −1, κ'' = +1. BBB Table 3 (s, t) ↔ (m, n) translation at (s, t) = (3, 1) West-coast gives m = t+s = 4, n = t−s = −2 ≡ 6 (mod 8), confirming (m, n) = (4, 6).

**Lift of the 4-spinor charge conjugation** U_4 = i·γ² = (iσ_y)_chir ⊗ (iσ²)_spin to the Krein space via `lorentzian_J_spatial_matrix(basis)`: J_L = J_{L, spatial} ⊗ K_t with K_t identity-unitary times complex conjugation. Why J_L² = +I (and not −I as for J_GV): the chirality factor i·σ_y in chirality space squares to −I_chir; combined with i·σ² in spin space which squares to −I_spin, their product squares to (−I_chir) ⊗ (−I_spin) = +I.

**Four BBB-predicted-sign axioms pass bit-exact** (Frobenius residual = 0.0 in float64) at every panel cell:
1. (i) J_L² = +I (BBB ε = +1, LOAD-BEARING L2D-FALS-1);
2. (ii) {J_L, γ⁵} = 0 (BBB ε'' = −1);
3. (iii) {J_L, η} = 0 (BBB εκ = −1);
4. (iv) J_L · D_L = +D_L · J_L (BBB universal Sec 5(v)).

**Structural finding (load-bearing scope finding):** the BBB universal axiom χD = −Dχ FAILS on truthful Camporesi-Higuchi D_GV. Residual ‖{γ⁵, D_L}‖_F = 2 ‖D_GV‖_F at N_t = 1 (6.00, 18.33, 38.88 at n_max = 1, 2, 3). This is NOT a basis-convention bug: it is the structural consequence of three mutually-inconsistent ingredients — GeoVac's chirality-as-γ⁵ identification (L2-B convention), chirality-diagonal D_GV (Paper 32 §III), and BBB universal axiom. **Three resolutions catalogued:** R1 (accept the structural finding: truthful is Krein-self-adjoint with BBB-predicted (4, 6) signs on J-relations, not a full BBB indefinite spectral triple); R2 (use offdiag CH: satisfies BBB axiom, breaks Riemannian-limit recovery bit-exactness); R3 (redefine γ⁵ as off-diagonal grading: rebuilds L2-B basis convention). **Recommended path for L2-E:** R1 + R2 in parallel (mirrors the Riemannian R3.5/Paper 42 truthful-vs-offdiag trade).

**M3 trivialization L0 prediction is CONVENTION-DEPENDENT.** Under Paper 28's n_fock-parity convention (the one that accesses Catalan G via D_even - D_odd quarter-integer Hurwitz), the (3, 1) chirality-symmetric truncation gives bit-identical D_even - D_odd to the Riemannian side at every n_max ∈ {1, ..., 5} — L0 prediction FALSIFIED under this reading. Under chirality-pairing convention, the truncation trivially gives 0 by symmetry — L0 prediction tautologically confirmed but does not access Paper 28's M3 content. **Structural refinement:** the BBB sign-flip {J, γ⁵} = 0 is a spectral-triple-axiom statement, NOT a vertex-parity-sum statement. M3's mechanism is sectional in n_fock-parity, NOT in spacetime signature. Refines the master Mellin engine domain partition cleanly (Paper 18 §III.7).

WH1 PROVEN is NOT re-opened — the load-bearing L2D-FALS-1 (J_L² = +I) passes bit-exact.

### §2.5 L2-E Krein-level modular Hamiltonian (2026-05-16 evening through 2026-05-17 morning)

**Verdict:** CLOSED-AT-FINITE-CUTOFF. Krein-level Paper 42 redo on the hemispheric wedge of S³ × ℝ at signature (3, 1).

**Construction:**

- **Wedge:** W_L = P_W^{spatial} ⊗ P_{t ≥ 0}, where P_W^{spatial} = (1/2)(I + R_polar) is the Paper 42 hemispheric wedge on H_GV (the m_j-reflection involution) and P_{t ≥ 0} is the t ≥ 0 half-line projector on ℂ^{N_t}. At N_t = 1 the temporal projector is the identity I_1.
- **BW choice of local Hamiltonian:** H_local := K_L^{α, W} / β at β = 2π, with K_L^{α, W} the wedge-restricted BW-α generator inheriting integer eigenvalues two_m_j from Paper 42 Def 5.1.
- **Wedge KMS state:** ρ_W^L = e^{-K_L^{α, W}}/Z, **β-independent at the algebra-action level** under the BW choice.

**Four LOAD-BEARING falsifiers pass bit-exact** at every tested (n_max, N_t) ∈ {1, 2, 3} × {1, 11, 21}:

| Falsifier | Statement | Max residual |
|:----------|:----------|:-------------|
| L2E-FALS-1 | BW-α period closure σ_{2π}^{L, α}(O) = O | ≤ 4×10⁻¹⁶ |
| L2E-FALS-2 | BW-γ Tomita period closure σ_{2π}^{L, TT}(a) = a | ≤ 4×10⁻¹⁶ |
| L2E-FALS-3 | Flow conjugacy σ_t^{L, TT}(a) = σ_{−t}^{L, α}(a) at general t | ≤ 4×10⁻¹⁶ |
| Riemannian-limit recovery | K_L^{α, W}\|_{N_t=1} = K_α^W, K_L^{TT}\|_{N_t=1} = K_{TT} bit-identically | 0.0 |

**Six-witness collapse** at the Krein level: all six instantiations (BW, HH_{M=1}, HH_{M=2}, Sew_{M=1}, Unruh_{a=1}, Unruh_{a=2}) give bit-identical Δ_L and K_L^{TT} — cross-witness consistency residual exactly 0.0 at every cell.

**The four-witness Wick-rotation theorem** (Hartle-Hawking + Sewell + Bisognano-Wichmann + Unruh) **lifts from "structural correspondence at the metric-functional level"** (Sprint TD Track 4, Sprint Unruh-pendant 2026-05-09/10) **to "literal identification at the operator-system level (Lorentzian, finite cutoff)"** at every Krein cutoff. **This is the framework's first operator-system-level literal identification of a Lorentzian QFT theorem at finite cutoff.**

---

## §3. The substantive structural findings

The five sub-sprints produced four substantive structural findings — three secondary (L2-B/C/D-level) and one headline (L2-E):

### §3.1 L2-B: Cl(3, 1) chiral-basis embedding is structurally compatible with Camporesi–Higuchi at the Hilbert-space level

The L2-A audit §5.7 named a structural worry: the Camporesi–Higuchi spatial spinor bundle (used in Paper 38 for WH1 PROVEN) might be incompatible with the Cl(3, 1) gamma-matrix embedding at signature (3, 1) West-coast. The bit-identical Riemannian-limit check passing at n_max ∈ {1, 2, 3} — with basis labels matching object-by-object via FullDiracLabel.__eq__ — resolves the worry at the structural level. The mechanism is concrete: GeoVac's chirality doubling (`full_dirac_basis` Weyl + anti-Weyl) is exactly the chirality structure that γ⁵ encodes in the chiral basis, and γ⁰ as the off-diagonal swap acts as the natural lift of the chirality-flip operation to H_GV. **WH1 PROVEN is NOT re-opened.**

### §3.2 L2-C: i^t = +i sign and ∂_t centered FD + Dirichlet zero BC both derived structurally

Two non-trivial convention choices that the L2-A audit Risk R2 named as "free parameters" both reduce to structural requirements rather than free choices:

- **Sign of i^t for (s, t) = (3, 1):** the van den Dungen Prop 4.1 prescription is "D_L = i^t · D̸_{g_r}". The L2-A audit named this as a sign ambiguity to verify. L2-C shows the sign +i is structurally forced by the Krein-self-adjointness requirement D_L^× = D_L. The alternative −i gives anti-Krein-self-adjoint D_L. Direct sympy verification via {γ⁰, D_GV} = 0 and ∂_t^† = −∂_t.
- **Temporal boundary conditions for ∂_t on the finite grid:** L2-A audit Risk R2 flagged four candidates (Dirichlet zero, Neumann, periodic, free-truncation). Only centered FD + Dirichlet zero gives anti-Hermitian ∂_t on a finite grid (forward/backward/upwind break Krein-self-adjointness; periodic re-creates CTCs on S³ × S¹). The choice is therefore structurally locked, not a free parameter.

### §3.3 L2-D: BBB universal axiom χD = −Dχ structurally fails on GeoVac's chirality-diagonal D_GV

The L2-C structural finding ({γ⁵, D_L} ≠ 0) was diagnosed and resolved at the L2-D BBB-classification level: the BBB universal axiom χD = −Dχ (BBB Sec 5(v), signature-independent) fails on truthful Camporesi–Higuchi D_GV. The mechanism is structural: GeoVac's chirality-as-γ⁵ identification (L2-B basis convention) + chirality-diagonal D_GV (Paper 32 §III definition) + BBB universal axiom are a mutually-inconsistent triple. Pick any two; the third must go. **Not a basis-convention bug** — no basis change can flip {γ⁵, D_GV} from zero (commutator) to nonzero (anticommutator) when both operators are diagonal in the chirality index.

The structural reading is that GeoVac's chirality grading is structurally aligned with the spectral content of D_GV (both diagonal in the same basis), which is structurally distinct from BBB's universal requirement that χ act as an anticommuting Z_2-grading on D. The L2-D structural finding is a load-bearing scope finding for the Lorentzian construction at (3, 1) West-coast: with truthful D_GV, the Krein construction is a "Krein-self-adjoint pair with BBB-predicted (4, 6) signs on the J-relations," not a full BBB indefinite spectral triple in the strict Sec 5(v) sense. **Three resolutions catalogued (R1/R2/R3)**, with R1+R2 in parallel recommended for L2-E.

**M3 trivialization L0 prediction is convention-dependent.** Under Paper 28's n_fock-parity convention, the L0 prediction is FALSIFIED; under chirality-pairing convention, it is tautologically confirmed but doesn't access M3 content. The structural refinement is that BBB sign-flip {J, γ⁵} = 0 is a spectral-triple-axiom statement, not a vertex-parity-sum statement — M3's mechanism is sectional in n_fock-parity, NOT in spacetime signature.

### §3.4 L2-E (the headline): Paper 42 §7.2 H_local verdict is signature-INDEPENDENT at the Riemannian limit

Paper 42 §7.2 names a load-bearing scope finding: the framework's intrinsic Camporesi-Higuchi Dirac D_W is NOT the right local Hamiltonian for the BW vacuum at β = 2π. The right local Hamiltonian is H_local := K_α^W / β, NOT D_W. This is the spectral-action--vs--modular-Hamiltonian generator distinction, listed as Paper 42 §10 open question O3.

Sprint L2-E's substantive new structural finding: **this distinction is signature-INDEPENDENT at the Riemannian limit.** Concretely:

$$\big\| H_{\mathrm{local}}\big|_{(3, 1), N_t = 1} - D_L^W\big|_{N_t = 1}\big\|_F \;=\; \big\| H_{\mathrm{local}}\big|_{(3, 0)} - D_{\mathrm{GV}}^W\big\|_F$$

bit-exact at every tested n_max. Values: 2.1332, 6.5275, 13.854 at n_max = 1, 2, 3 respectively.

**Mechanism of the bit-exact match:** at N_t = 1, D_L = i·D_GV exactly (the temporal derivative ∂_t vanishes on the singleton grid). Wedge restriction gives D_L^W = i·D_GV^W. The norm ‖H_local − i·D_GV^W‖_F equals ‖H_local − D_GV^W‖_F because H_local has real diagonal entries (K_L^{α, W}/β has real integer-rational diagonal).

At N_t > 1 the Lorentzian-side residual is **refined upward** by the temporal-derivative content of D_L = i(γ⁰ ⊗ ∂_t + D_GV ⊗ I): the spectral-action object has more content than D_GV, so the gap between H_local and D_L^W widens monotonically with temporal cutoff.

**Refinement of Paper 42 §10 O3:** from "open question about signature-dependence" to **"signature-INDEPENDENT structural distinction at the Riemannian limit, refined upward by temporal-derivative content at N_t > 1."** The distinction is NOT a peculiarity of the round-S³ Riemannian truncation — it is a deeper structural feature of the framework's modular content, surfacing across signatures. The **deeper** question (why does the spectral-action Dirac fail to be the modular generator on the BW vacuum, structurally?) remains the residual O3 open question.

---

## §4. What this means for WH1 / Paper 38 / Paper 40 / Paper 42 framing

- **WH1 PROVEN survives.** Sprint L2 is structurally additive (a new (3, 1) layer on top of the Riemannian foundation), not a re-test. The load-bearing Riemannian-limit checks at N_t = 1 (L2-B Hilbert space, L2-C Dirac operator, L2-D Connes axioms with `J_GV² = -I` on the Riemannian side preserved bit-exactly, L2-E modular Hamiltonian) all pass bit-identically, confirming that the Lorentzian construction sits on top of WH1's Riemannian foundation without disturbing it.

- **Paper 38 unchanged.** The qualitative-rate Latrémolière-propinquity GH-convergence theorem on the Riemannian truncated triple is the WH1 PROVEN result; Sprint L2 does not touch it.

- **Paper 40 unchanged.** The unified compact-Lie-group propinquity convergence theorem with universal 4/π rate constant is independent of signature.

- **Paper 42 §10 O1 (Lorentzian extension) closed at finite cutoff.** Paper 42 §10 O3 (spectral-action vs modular-Hamiltonian generator distinction) refined from "open question about signature-dependence" to "signature-INDEPENDENT structural distinction at the Riemannian limit." Paper 42 §11 new section "Lorentzian closure at finite cutoff" documents the closure; the math.OA-style standalone writeup is outlined as Paper 43 (companion to Papers 38, 39, 40, 42).

- **Paper 32 §VIII.E gains a fifth subsection** (§VIII.E.E "Krein-level unified-strong four-witness theorem") closing the Lorentzian-extension arc at the operator-system level on the Krein space.

- **Paper 34 §III.27 Wick rotation projection** lifts from "structural correspondence at the metric-functional level" to "literal identification at the Krein operator-system level (finite cutoff)" — the Lorentzian-side analog of Sprint L1's Riemannian-side closure of §III.27.

- **Paper 31 §8 Sig/Op partition** empirically verified at 28/28 projections by the combination of Sprints L1 + L2 (the partition was an a priori prediction at Sprint L0).

---

## §5. What's NOT closed by Sprint L2

The honest scope discipline of Sprint L2 names three distinct frontiers that remain open as separate sprints:

### §5.1 Continuum Lorentzian propinquity (Sprint L3)

Sprint L2-E closes the operator-system-level Wick-rotation theorem **at finite n_max and finite N_t**. The continuum-limit / Lorentzian-propinquity-rate-level extension remains open and is the named Sprint L3 target. **No published Lorentzian propinquity construction exists as of May 2026.** Latrémolière 2017/2026, Toyota 2023, Hekkelman–McDonald 2024 a/b, Leimbach–van Suijlekom 2024, and Farsi–Latrémolière 2024/2025 cover only Riemannian / Hilbert-space structures.

Constructing a Lorentzian propinquity is **original NCG-mathematics work at the 6-12 month scale**. The Nieuviarts 2025 twist-morphism (arXiv:2502.18105) is a candidate shortcut — but requires a separate L2-Nieuviarts-scoping pass to verify applicability to GeoVac S³ = SU(2) (odd-dim caveat that Sprint L0 audit flagged).

Sprint L3 is named open. The framework's Lorentzian-side claims sit at finite cutoff only.

### §5.2 Cross-manifold W2b ($\mathcal{T}_{S^3} \otimes \mathcal{T}_{\mathrm{Hardy}(S^5)}$)

The cross-manifold W2b wall (named in the multi-focal sprint outcome of 2026-05-08) is **structurally distinct** from the (3, 1) extension closed by Sprint L2. W2b targets the tensor product of the round-S³ Camporesi-Higuchi triple with a holomorphic-sector S⁵ Bargmann-Segal triple (Paper 24 §V); it is blocked at the NCG-framework level by the Coulomb/HO category mismatch (Riemannian-vs-Hardy-sector asymmetry).

Sprint L2 closes the Lorentzian extension within $\mathcal{T}_{S^3}$ at signature (3, 1); it does NOT address the cross-manifold question. Paper 32 §VIII.D "frontier-of-field framing" preserves W2b as an open structural extension requiring NCG-framework-level innovation, distinct from the L2 architecture.

### §5.3 Calibration data (W3, second packing axiom)

Sprint L2 is structurally additive at the spectral-triple-machinery level. It does NOT address the inner-factor-input-data question (W3 in CLAUDE.md §1.7 multi-focal-wall taxonomy):

> What generates calibration data (Yukawas, gauge couplings, mass spectra) in the framework?

W3 remains open. Sprint L2's Lorentzian extension reproduces the existing Riemannian framework's calibration content at signature (3, 1); it does NOT generate new calibration data.

The Sprint W3 spectral-zeta candidate (2026-05-08) was tested and falsified the same evening by three independent follow-up tracks; the W3 calibration-data question is back to the 14-vague-speculations / 0-concrete-proposals state. Sprint L2 does not change this.

---

## §6. Recommended next-direction options

Four options surfaced from the synthesis pass. Each is well-scoped at the 1–6 month range. The PI should select based on strategic priorities.

### Option A: Draft Paper 43 (math.OA standalone Lorentzian-extension paper)

**Scope:** ~2 weeks of writing + 1 week of review. Sibling to Papers 38, 39, 40, 42 in `papers/standalone/`. Natural Zenodo deposit when ready.

**Rationale:** Sprint L2 produced a clean operator-system-level closure that warrants a math.OA-style standalone writeup. The closure is non-trivial (first operator-system literal identification of a Lorentzian QFT theorem in the framework); the construction is clean (van den Dungen 2016 Prop 4.1 + BBB 2018 at (m, n) = (4, 6) + Paper 42 verbatim with temporal slot); the H_local signature-independence finding is a substantive new structural result.

**Risk:** thin — Sprint L2 closes at finite cutoff only, not in the GH limit. The Paper 43 writeup would explicitly NOT supersede Paper 42 — Paper 42 is the Riemannian closure, Paper 43 is the Lorentzian extension, and together they constitute the four-witness Wick-rotation arc at the operator-system level.

**Outline drafted:** `papers/standalone/paper_43_lorentzian_extension_outline.md`.

### Option B: Open Sprint L3 (Lorentzian propinquity construction)

**Scope:** multi-month (6-12 months as audited at L2-A); possibly shortened via the Nieuviarts 2025 twist-morphism if L2-Nieuviarts-scoping verifies applicability to S³ = SU(2) (odd-dim caveat).

**Rationale:** the natural mathematical follow-on to Sprint L2's finite-cutoff closure. Would close the framework's Lorentzian-side claims at the qualitative-rate / propinquity level, matching the Riemannian-side WH1 PROVEN status.

**Risk:** highest of the four options. Original NCG-mathematics work. Sprint L0 audit explicitly noted L3 "may be skippable" — the operator-system-level finite-cutoff identification is sufficient for the four-witness Wick-rotation theorem.

### Option C: Return to state-side dictionary direction

**Scope:** 1-2 months. Sprint 3+ targets in the state-side complement of Paper 34's projection dictionary (§III.28 apparatus identity + the candidate state-side dictionary enumeration in Paper 34 §VIII open questions: mutual information, conditional entropy, fidelity, trace distance, relative entropy, Wasserstein-Kantorovich).

**Rationale:** the state-side complement was opened by §III.28 (Sprint TD Track 2 + Track 5) but only partially explored. The dictionary-completion arc (Sprint 1+2+3, 2026-05-15) closed at 28 projections on the spectral side; the state-side has only one entry so far (§III.28). Natural follow-on to fill out the state-side panel.

**Risk:** low. The state-side mechanism is already established (apparatus identity, PSLQ-disjoint from master Mellin engine). Filling out the panel is mostly mechanical.

### Option D: Return to physics-side precision catalogue

**Scope:** open-ended. Multi-focal-composition wall, §V.D convention exposures, §V.C autopsies — the physics-side arc that was active before the L0 audit was triggered.

**Rationale:** the physics-side arc has accumulated substantial results (9 systems in the precision catalogue, 6 §V.D convention exposures, 6 §V.C autopsies). Continuing the arc would deepen the framework's empirical reach.

**Risk:** low. Established arc with known target classes.

---

## §7. Memory and CLAUDE.md status snapshot

**Files modified:**

- `papers/standalone/paper_42_modular_hamiltonian_four_witness.tex` — new §11 "Lorentzian closure at finite cutoff" (~250 lines); §10 O3 refinement paragraph
- `papers/synthesis/paper_32_spectral_triple.tex` — new §VIII.E.E "Krein-level unified-strong theorem" (~200 lines); new bibitem `paper42`
- `papers/observations/paper_34_projection_taxonomy.tex` — §III.27 status flip + L2-E closure paragraph; §V.E new L2-E closure paragraph
- `papers/core/paper_31_universal_coulomb_partition.tex` — §8 new "Empirical verification of the partition by Sprints L1 + L2" subsection
- `CLAUDE.md` — §2 new L2-E sprint bullet (~600 words); new Sprint L2 closure synthesis bullet (~500 words); §1.7 WH1 PROVEN entry extended with L2-D/L2-E sentences
- `memory/MEMORY.md` — new index entry pointing to `sprint_l2_closure.md`

**Files created:**

- `debug/sprint_l2_synthesis_memo.md` — this synthesis memo (~3500 words)
- `papers/standalone/paper_43_lorentzian_extension_outline.md` — Paper 43 outline (~2000 words)
- `memory/sprint_l2_closure.md` — session-restoration memory file

**Compile status (three-pass):** Paper 42 SUCCESS; Paper 32 SUCCESS; Paper 34 SUCCESS; Paper 31 SUCCESS. Three pre-existing warnings on Paper 34 unrelated to Sprint L2 edits (`sec:matches`, `sec:curvature_coefficients`, `tab:catalog_off`).

**WH1 PROVEN status:** unchanged. Sprint L2 is structurally additive on top of the Riemannian foundation; all load-bearing falsifiers pass.

End of memo.
