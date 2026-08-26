# Sprint memo — Quantum cost model for a non-Hermitian (xTC) ground-state operator

**Date:** 2026-08-23
**Scope:** Grounded cost model only (literature verification + linear-algebra reasoning). No paper or `geovac/` code touched. Feeds the decision on whether the xTC accuracy win (see `debug/sprint_xtc_poc_li_memo.md`) can survive when the effective operator is carried onto a quantum computer.

**The operator.** The xTC effective Hamiltonian H̃ = H + D + K + L3(contracted) is **non-Hermitian but diagonalizable** with a **real** ground eigenvalue and a real ground eigenvector (PoC uses `scipy.linalg.eig`, takes the real ground state). Write H̃ = V Λ V⁻¹. The new axis this operator introduces, absent from Paper 14's Hermitian cost analysis (Pauli count, 1-norm λ, real-time Trotter), is the **eigenvector condition number κ_V = κ₂(V)** — the Bauer–Fike / departure-from-normality factor. Bauer–Fike is the reason κ_V must appear in *any* spectral method on H̃: a perturbation of size δ to the matrix moves an eigenvalue by up to κ_V·δ, so eigenvalue estimation to precision ε on a non-normal operator cannot cost less than ~κ_V/ε. Everything below is organized around where each route pays κ_V.

Paper 14's Hermitian machinery does NOT transfer directly:
- **VQE Rayleigh quotient** ⟨ψ|H|ψ⟩/⟨ψ|ψ⟩ is not stationary at, and does not bound, an eigenvalue of a non-Hermitian H̃.
- **Real-time Trotter / textbook QPE** require H̃ Hermitian so that e^{-iH̃t} is unitary; e^{-iH̃t} is not unitary here.
- **1-norm λ** (block-encoding subnormalization) still makes sense — H̃ is still a sum of Paulis (complex coefficients) — and is **measured comparable to plain** (PoC: 1-norm ~0.90× plain; not densified). λ is therefore NOT the load-bearing risk.

---

## The four routes

### Route 1 — QEVE / QEVT (Low & Su, arXiv:2401.06240, SIAM J. Comput. 2024/2026)

**Primary source verified (full-text HTML).** This is the purpose-built tool: eigenvalue processing of **non-normal** block-encoded operators, explicitly tabulating measures of non-normality (Jordan condition number, numerical range, pseudospectrum).

- **QEVE (estimation, for operators with real spectra) — Theorem 3.** Query complexity to the block encoding:
  **O( α · κ_S / ε · log(1/p_fail) )**
  where κ_S = condition number of the diagonalizing (Jordan/eigenvector) matrix — this is exactly κ_V — α = block-encoding normalization (= our λ), ε = precision, p_fail = failure probability.
  - **κ_V enters LINEARLY.** This is the Bauer–Fike floor; it is optimal, not an artifact — no eigenvalue method can beat linear-in-κ_V on a non-normal operator.
  - α (=λ) enters linearly (as in the Hermitian qubitization cost).
  - ε: 1/ε Heisenberg-limited (the diagonalizable-input Heisenberg scaling is the paper's headline).
  - Spectral gap Δ: no explicit query multiplier, BUT to *identify* the ground eigenvalue distinctly you must estimate to ε < Δ; effective cost O(α·κ_V/min(ε,Δ)). Gap enters as a resolution requirement, same as Hermitian QPE.
- **QEVT (transformation, e.g. a spectral-projector polynomial for ground-state prep) — Theorem 4.** Worst-case **κ_V² · ‖p‖ / ‖p(A/α)|ψ⟩‖ · log(...)** — quadratic in κ_V, plus a polynomial-norm and an initial-overlap factor. QEVT is the tool if you need the eigenSTATE; QEVE is the tool if you need the eigenVALUE.

**Verdict-lean: (i) MILD in κ_V** — linear, and linear is the best achievable. QEVE is the natural fit because the xTC operator has a real spectrum and a real ground state, satisfying QEVE's hypothesis. GO if κ_V stays O(1–10).

**Directly-analogous primary source — arXiv:2511.21867 (2026), "Accuracy and resource advantages of quantum eigenvalue estimation with non-Hermitian transcorrelated electronic Hamiltonians"** (verified, full-text HTML). This paper does *exactly this route on exactly this operator class* (QEVE on TC/xTC electronic Hamiltonians), and it both confirms the scaling and supplies the missing empirical number:
- Reports the QEVE T-gate count as **∝ α·κ_S/ε** (concrete prefactor ≈ (12K+4μ)·368800·√3·π), i.e. linear in κ_S — matches Theorem 3.
- **Measured Jordan condition numbers κ_S ≈ 1–3 for full TC Hamiltonians (second-row atoms, STO-6G); xTC reduces κ_S further ("quite significantly in almost all instances").** This is the load-bearing empirical fact: κ_V is MILD for transcorrelated electronic Hamiltonians, and xTC (GeoVac's route) makes it milder.
- TC/xTC 1-norm ≈ comparable-to-slightly-reduced vs the Hermitian STO-6G one-norm (~8–82 Ha), corroborating GeoVac's own "not densified, ~0.90×" measurement.
- Their honest cap: TC's *accuracy* advantage survives only for small systems (STO-6G TC error grows for O/F/Ne); "turning to the transcorrelated picture does not always lead to better quantum gate counts." Crucially, that erosion is a **basis-accuracy** effect, NOT a κ_V blow-up — the non-Hermitian penalty itself stays mild. Qubit-count advantage is universal.
- **Note for PI:** this paper (Nov-2025/2026) largely pre-empts the QEVE-on-xTC framing; GeoVac's distinct angle would be the Coulomb-Sturmian basis + the measured cusp/sparsity edge, not the QEVE-on-TC idea itself.

### Route 2 — Non-Hermitian QITE (Motta et al., arXiv:1901.07653; McArdle et al. VarQITE arXiv:1804.03023; TC test arXiv:2407.10523)

Imaginary-time evolution e^{-βH̃}|ψ₀⟩/‖·‖ projects onto the eigenvector of smallest Re(eigenvalue) — this **holds for non-Hermitian H̃** (verified: arXiv:2407.10523 states "ITE also holds for non-Hermitian Hamiltonians, hence VarQITE is applicable for TC systems," via the similarity-transform relation).

Cost drivers:
- **β to converge ~ (1/Δ)·log(1/ε)** — ground-state error ~ e^{-2βΔ}. This is a **poly(1/Δ)** dependence, and it is *shared with every Hermitian ground-state method* (QPE gap resolution, adiabatic, cooling). Δ, not κ_V, dominates.
- **Correlation-domain / k-locality cost.** Base (non-variational) QITE (Motta) approximates each non-unitary step by a unitary supported on a domain of D qubits set by the correlation length; the measurements to determine the step scale **exponentially in the domain size (~3^D–4^D)**. This is a k-locality penalty, present even in the Hermitian case. For xTC the effective operator is **2-body (k=2), same locality as the Coulomb ERI** — so this penalty is NOT made worse by TC. Variational QITE (McArdle) trades the domain-exponential for a fixed ansatz + an ill-conditioned McLachlan/QFI linear system M ẋ = V solved each step (its own conditioning issue).
- **Non-Hermitian extras:** (a) energy/expectation values become complex → measure real AND imaginary parts (constant ~2× overhead; verified 2407.10523); (b) the norm ‖e^{-βH̃}‖ is non-monotone with **transient growth governed by the pseudospectrum ~κ_V**, so the post-selection / renormalization sampling picks up a factor that scales with κ_V during the transient; (c) McLachlan matrix loses Hermiticity/PSD structure. Empirically, 2407.10523 report stable convergence on TC (real part monotone-decreasing, imaginary part → 0) — no κ_V instability observed at their scale.

**Verdict-lean: MILD in κ_V, but (ii) poly(1/Δ) in the gap.** The load-bearing penalty for this route is Δ and the correlation-domain exponential (both shared with Hermitian ground-state finding), not κ_V. κ_V shows up only as a transient-norm sampling factor. Route 2 is a reasonable NISQ fallback; its risk is the generic gap/locality risk, not non-Hermiticity.

### Route 3 — Bivariational / non-Hermitian VQE and NOQE (arXiv:2507.04783; arXiv:2205.09039; arXiv:2608.12830)

Left/right (biorthogonal) ansatz, generalized Rayleigh quotient **L(θ,φ) = ⟨φ|H̃|ψ⟩ / ⟨φ|ψ⟩** (verified, 2507.04783).

Cost driver — **measurement-variance amplification by conditioning** (verified, 2507.04783): the estimator variance is magnified by the **reciprocal of the biorthogonal overlap squared, 1/|⟨φ|ψ⟩|²**, and a small left–right overlap is "a signature of extreme non-normality." For a near-defective eigenvector |⟨L|R⟩| ~ 1/κ_V, so:
- **Base bivariational: measurement cost ~ κ_V² / ε²** — QUADRATIC in κ_V, and /ε² (no Heisenberg limit; near-term estimator).
- **NOQE (2205.09039):** solves the GEVP H c = E S c from non-orthogonal references; base measurement cost **O(M³)** (M = #references), eigenvalue error amplified by κ(S) (overlap-matrix conditioning — a cousin of κ_V).
- **Improved NOQE (2608.12830, verified):** with overlap thresholding, per-element shot count drops to **O(M)**, and "eigenvalue sensitivity is controlled by the condition number of the *retained* overlap matrix" — i.e. **LINEAR in κ (of the thresholded block)**, and in practice sub-linear on structured (H-chain/ring) instances.

**Verdict-lean: base (ii) quadratic in κ_V (fatal if κ_V large); improved variant (i) linear-in-κ but still /ε² measurement.** For the measured κ_V ~ 1–3 regime this is acceptable as a NISQ route; it degrades fast if κ_V ever climbs.

### Route 4 — Hermitian folding: solve H̃†H̃ or (H̃+H̃†)/2

Attempt to dodge non-Hermiticity by building a Hermitian surrogate. Both variants are worse, and the first is **biased**:

- **H̃†H̃ (PSD; ground = smallest singular value σ_min):**
  - **Eigenvalues of H̃†H̃ = singular values² of H̃, which for a NON-normal H̃ are NOT |eigenvalues|.** Its ground state is the min right-singular vector, ≠ the ground eigenstate of H̃. So folding **computes the wrong object** unless H̃ is normal. To recover an eigenvalue you must root-search E over σ_min(H̃ − E·I) = 0; near an ill-conditioned eigenvalue the pseudospectrum is fat, so σ_min(H̃−E·I) stays tiny over a window of width **~κ_V·ε**, reintroducing κ_V into the E-resolution.
  - **κ squares:** κ(H̃†H̃) = κ(H̃)².
  - **λ squares:** block-encoding the product multiplies subnormalizations → ~λ²; Pauli count multiplies (product of two operators; the ~5% Gaunt/xTC sparsity is destroyed).
  - **gap squares/shrinks:** relevant gap becomes σ₂² − σ₁² ≈ (σ₂−σ₁)(σ₂+σ₁).
- **(H̃+H̃†)/2 (Hermitian):** λ_min gives the bottom of the numerical range (Bendixson) — a **bound**, equal to Re(E₀) only when H̃ is normal. **Biased for non-normal H̃**, not the eigenvalue.

**Verdict-lean: (ii) fatal AND biased** — squares κ, squares λ, shrinks the gap, destroys sparsity, and (H̃†H̃) computes singular values / (symmetrized) computes a numerical-range bound rather than the eigenvalue. Correct only in the normal limit — precisely where folding is unnecessary. Worst of the four.

---

## Summary table

| Route | Cost scaling | Driving property | Verdict-lean |
|:------|:-------------|:-----------------|:-------------|
| 1. QEVE (estimation, real spectrum) | O(α·κ_V/ε·log 1/p) | κ_V **linear** (Bauer–Fike floor); α=λ linear; ε Heisenberg | **(i) MILD** — GO if κ_V=O(1–10). κ_V~1–3 measured for TC (2511.21867), xTC lowers it |
| 1'. QEVT (state prep) | ~κ_V²·‖p‖/‖p·ψ‖·log | κ_V **quadratic** | (ii) worse; use only if the eigenSTATE is needed |
| 2. Non-Herm. QITE | β~(1/Δ)log(1/ε) × (3^D–4^D domain); +2× complex; +κ_V transient sampling | **Δ (poly 1/Δ)** + k-locality; κ_V only a transient factor | (i) mild in κ_V / (ii) poly(1/Δ). xTC k=2 ⇒ locality not worse than Coulomb |
| 3a. Bivariational VQE | ~κ_V²/ε² measurements | 1/|⟨L|R⟩|² ~ **κ_V² (quadratic)** | (ii) quadratic — fatal if κ_V large; OK at κ_V~1–3 |
| 3b. Improved NOQE (thresholded) | O(M) shots/elt; ~κ/ε² | κ of retained overlap **linear** | (i) linear-in-κ, but /ε² near-term |
| 4. H̃†H̃ or (H̃+H̃†)/2 | κ², λ², gap², sparsity lost | κ **squared** + **spectrum changed** | (ii) FATAL + biased (computes singular values / numerical range, not eigenvalues) |

## The single load-bearing property

**κ_V = κ₂(V), the eigenvector condition number (departure-from-normality / Bauer–Fike factor) of H̃.**
It is the only genuinely NEW axis the TC transform introduces: λ is measured comparable to plain (~0.90×, not densified) and Δ is the generic gap that every ground-state method — Hermitian or not — must pay. The best route (QEVE) pays κ_V *linearly*, which is the optimal, unavoidable Bauer–Fike price; the fold route (4) pays it *squared* and also computes the wrong object; the bivariational route (3a) pays it *squared* in measurements. So the entire quantum-cost GO/STOP for xTC hinges on **whether κ_V(H̃_xTC) stays O(1–10)**.

**Prognosis: favorable but unmeasured for GeoVac's operator.** The directly-analogous primary source (arXiv:2511.21867) measures κ_S ≈ 1–3 for full TC electronic Hamiltonians and shows xTC reduces it further — i.e. mild, and xTC is on the good side. GeoVac has NOT yet measured κ_V for its own **Coulomb-Sturmian** xTC operator. The eigenvectors are already computed inside the PoC's `scipy.linalg.eig` call — the one missing number is a two-line addition: `kappa_V = numpy.linalg.cond(V)` on the right-eigenvector matrix at production γ. That measurement converts this cost model from "prognosis" to "decision."

---

## Verified references (arXiv ID → what it actually says → extracted cost)

1. **arXiv:2401.06240** — Low & Su, *Quantum eigenvalue processing*, SIAM J. Comput. (2024/2026). Full-text verified. **QEVE Thm 3: O(α·κ_S/ε·log(1/p))**, κ_S = eigenvector/Jordan condition number, **linear**; Heisenberg 1/ε; handles complex spectra and ill-conditioned Jordan bases. **QEVT Thm 4:** worst-case **κ_S²**. Tabulates non-normality measures (Jordan cond, numerical range, pseudospectrum).
2. **arXiv:2511.21867** (2026) — *Accuracy and resource advantages of QEVE with non-Hermitian transcorrelated electronic Hamiltonians*. Full-text verified. QEVE T-count **∝ α·κ_S/ε** (linear, matches Thm 3). **Measured κ_S ≈ 1–3 for full TC (2nd-row atoms, STO-6G); xTC reduces it "quite significantly."** TC/xTC 1-norm comparable to Hermitian (~8–82 Ha). Verdict: TC accuracy survives for small systems; gate-count erosion for O/F/Ne is a *basis-accuracy* effect, not κ_V. Qubit advantage universal. (This paper does the QEVE-on-TC route directly — largely pre-empts that framing.)
3. **arXiv:1901.07653** — Motta et al., *Determining eigenstates and thermal states … via QITE*, Nat. Phys. 16, 205 (2020). Abstract verified; standard results: e^{-βH} projects to ground state, error ~e^{-2βΔ} ⇒ **β~(1/Δ)log(1/ε)**; per-step unitary supported on a correlation domain of D qubits with measurement cost **exponential in D (~3^D–4^D)**.
4. **arXiv:1804.03023** — McArdle et al., *Variational ansatz-based QITE (VarQITE)*. Fixed ansatz; solves McLachlan/QFI linear system each step (conditioning of that Gram matrix is the variational cost issue). O(m) state preps per step for the QFI matrix.
5. **arXiv:2407.10523** — *VarQITE for MPS ansatz with tests on transcorrelated Hamiltonians*. Full-text verified. "ITE holds for non-Hermitian Hamiltonians ⇒ VarQITE applicable to TC." Energy is complex-valued (measure real+imag parts, ~2×); observed **stable** convergence (real part monotone ↓, imag → 0). Does NOT give a rigorous gap/κ_V cost.
6. **arXiv:2507.04783** — non-Hermitian variational eigensolver. Full-text verified. Biorthogonal loss **⟨φ|H|ψ⟩/⟨φ|ψ⟩**; variance amplified by **1/|⟨φ|ψ⟩|² (~κ_V²)**; small left–right overlap = signature of extreme non-normality.
7. **arXiv:2205.09039** — Baek et al., NOQE, PRX Quantum 4, 030307 (2023). GEVP Hc=ESc from non-orthogonal refs; base measurement **O(M³)**, error amplified by κ(S).
8. **arXiv:2608.12830** — *Improved Measurement Cost Scaling in NOQE*. Abstract verified. Overlap thresholding ⇒ per-element shots **O(M)** (from O(M³)); sensitivity **controlled by κ of the retained overlap (linear-in-κ)**; sub-linear in practice on H-chains/rings.

(Supporting, from prior-art scan `debug/sturmian_quantum_priorart.md`, not re-verified here: arXiv:2112.02554 fault-tolerant GEVP with explicit Õ(κ_B) metric cost; arXiv:2506.13534 outer-scan-over-E + inner QSVE, Õ(√(NK)).)
