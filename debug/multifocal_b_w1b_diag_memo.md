# Phase B Wall Diagnostic — W1b: Magnetization-Distribution Cross-Register

**Date:** 2026-05-07
**Sub-sprint:** B-W1b-diag (scoping, not closure)
**Author:** Sub-agent (PM dispatch; no production code modified)
**Source walls:** Phase A synthesis §1 (W1b row); HF-4 Zemach failure memo; multi-focal-wall pattern
**Frame:** GeoVac is an almost-commutative spectral triple in the Marcolli–van Suijlekom lineage (WH1 PROVEN). The multi-focal-composition wall has six refined sub-walls; W1b is the magnetization-distribution focal-length wall, physically distinct from W1a (recoil = color-electric coordinate operator) because the proton magnetization is a color-magnetic distribution and the Zemach radius $r_Z$ is a categorically different focal length from the charge radius $R_p$.

---

## 1. Diagnostic question

The HF-4 sprint exhibited W1b empirically: the residual to the 21 cm line could not be closed by any GeoVac-internal manipulation, and substituting the externally measured Zemach radius $r_Z = 1.045$ fm as a Layer-2 focal length closed the residual to $+18$ ppm — within Eides Tab. 7.3's expected multi-loop + nuclear-polarizability budget. The diagnostic question is therefore not "can GeoVac compute $r_Z$?" (it cannot, by Wilson-SU(3) matter-element CG obstruction) but rather: **what would a magnetization-distribution operator on the proton register actually be, what input data would it need, what calibration neighborhood is it in, and is it reducible to W1a?**

The five sub-questions Q1–Q5 each get a short section. Q5 (the verdict) is then the load-bearing part for the gating decision into Phase C.

---

## 2. Q1 — QCD-internal input: scalar $r_Z$ vs full $\rho_M(r)$ vs moment list

The Zemach correction in Eides §7.2 is

$$\frac{\Delta \nu_Z}{\nu_F} = -2 Z \alpha m_e r_Z, \qquad r_Z = \int d^3 r \int d^3 r' \, \rho(r')\, m(|r - r'|)$$

where $\rho(r)$ is the charge distribution and $m(r)$ is the magnetization-density distribution, both normalized to unit total. The Zemach radius is the *first moment* of the charge–magnetization convolution. Higher-order corrections (Friar 1979; Karshenboim 2005 review; Eides–Grotch–Shelyuto 2007 monograph) bring in additional moments: the Friar moment $\langle r^3 \rangle_{(2)}$ for the third-Zemach moment, polarizability terms $\Delta_\text{pol}$ that depend on the off-shell two-photon-exchange amplitude, and recoil-corrected Zemach $r_Z(1 + \mathcal{O}(m_e/m_p))$.

For the leading hyperfine correction at parts-per-million on hydrogen, **the input is the scalar $r_Z$**. For sub-ppm precision (the Eides 2024 paper [doi:10.1016/j.physletb.2024.139049] addresses exactly this regime, halving the polarizability uncertainty), the input upgrades to a *finite list of moments* — typically $\{r_Z, \langle r^3 \rangle_{(2)}, \Delta_\text{pol}\}$ — each independently extracted from elastic and inelastic electron-proton scattering form-factor data combined with spin-structure-function input. The full $\rho_M(r)$ as a function is *not* the natural input to atomic-physics calculations; the moments are. This is because the convolution $\rho \star m$ in the hyperfine matrix element sees only finitely many radial moments at any finite order in $\alpha$.

The Eides 2024 paper [Phys. Lett. B (2024)] specifically uses *new spin-structure data* to reduce the input uncertainty on $\Delta_\text{pol}$ to $\pm 0.6$ ppm; their structural form is unchanged from Eides–Grotch–Shelyuto 2007: a scalar Zemach radius times the Fermi contact frequency, plus polarizability and recoil corrections each with its own scalar input. The *input data structure* for GeoVac is therefore: $\{r_Z, \Delta_\text{pol}, \mathcal{O}(\alpha)\text{-recoil terms}\}$, all scalars from QCD-internal data.

**Q1 answer:** the natural input is a **finite moment list**, dominated at parts-per-million precision by the single scalar $r_Z$. Full $\rho_M(r)$ is available from form-factor parametrizations (Bernauer–Mainz 2014; Lin–Hammer–Meißner 2021) but is not what the hyperfine calculation consumes — only its finite moments are. GeoVac's natural Layer-2 input is therefore *one to three scalar focal lengths*, each tagged as a calibration constant from QCD/scattering data, not the function $\rho_M(r)$.

This is structurally analogous to how the Sommerfeld fine-structure constant $\alpha$ enters Paper 28 calculations: a single scalar from CODATA, not a function. The categorical distinction from the charge radius $R_p$ is that $R_p$ enters $\Delta E_{1s}$ via the Foldy correction ($\rho_E$ moment), while $r_Z$ enters $\Delta\nu_\text{hf} / \nu_F$ via the convolution moment ($\rho_E \star \rho_M$). They are different physical quantities, computed from different form factors at different $Q^2$ ranges, and one cannot be derived from the other without additional QCD input.

---

## 3. Q2 — Eides 2024 framework: structural form and translation cost

The Eides–Grotch–Shelyuto (2007 monograph; refreshed in the 2024 PLB paper) framework is **classical-r-substitution into a smeared Fermi contact** *as far as the atomic-physics-side Hamiltonian is concerned*. The proton is treated as an external classical source whose finite extent is encoded in a single scalar $r_Z$ that multiplies the point-like Fermi contact term:

$$H_\text{hf}^\text{Zemach} = (1 - 2 Z \alpha m_e r_Z) \, A_\text{hf}^\text{point}\, \mathbf{I}\cdot\mathbf{S}$$

This is structurally option (a) in the diagnostic question: classical-r-substitution. The proton's spatial wavefunction does *not* appear as an operator on a quantum register; the proton is a point with a scalar size correction.

Crucially: the "operator-level prescription on a multi-particle Hamiltonian" (option b) is how this is written in NRQED and pNRQED *but only as a matching theory* — at every step in the matching chain, one operator-valued atomic-physics computation matches against an external proton-side scattering-data input. The proton register is never quantized in the Eides program. The closer-to-(b) treatments (Pachucki–Patkóš–Yerokhin 2023 PRL recoil; Khabibullin–Karshenboim 2023 polarizability) all use Foldy–Wouthuysen reduction to *eliminate* the nuclear coordinate operator and absorb its content into scalar shifts on the electronic Hamiltonian — they go *toward* the structural form GeoVac already has, not toward operator-valued nuclear coordinates.

There is a more sophisticated treatment (option c): Carlson–Vanderhaeghen 2008 PRA two-photon-exchange amplitude, where the proton's intermediate states are summed over in a dispersion integral and the proton is *not* treated as a static potential. But even there, the *output* of the dispersion integral is a finite scalar correction $\Delta_\text{pol}$ to the atomic-physics Hamiltonian. The "operator" lives at the proton vertex of the two-photon-exchange diagram; on the atomic-physics Hilbert space it is again a scalar.

**Q2 answer:** the Eides 2024 framework is structurally option (a) — classical-r-substitution — *at the level of the atomic-physics Hamiltonian*. Translation to GeoVac's nuclear register is therefore *trivial as a calibration prescription*: substitute $r_Z$ as a Layer-2 input. This is exactly what HF-4 did at +18 ppm closure. **The Eides framework does not in fact build a magnetization-distribution operator on a quantum nuclear register; it presupposes that the proton is a classical scalar source.**

This is informative: Eides is the calibration neighborhood for the *answer*, not for an *operator construction*. If GeoVac wants a magnetization-distribution operator on the proton register, the Eides framework provides the calibration ($r_Z$ as a number) but does *not* provide a precedent for the operator structure. The operator construction (if pursued) is a GeoVac-internal mathematical exercise calibrated against the Eides $r_Z$ output, not a port from Eides.

---

## 4. Q3 — Operator construction sketch

A magnetization-distribution operator on the proton register would need to do three things: (i) act on the proton's internal quantum numbers, (ii) carry a length scale (the Zemach radius), (iii) couple to the electron 1s contact density via convolution. The natural construction has three layers, increasing in cost:

### 4a. Multipole-expansion operator (cheapest)

Expand the magnetization density as $\hat{m}(\hat{\mathbf{R}}_p) = \sum_{L,M} \hat{m}_{LM}\hat{O}_{LM}(\hat{\mathbf{R}}_p)$ with operator-valued multipole amplitudes $\hat{m}_{LM}$ and tensor operators $\hat{O}_{LM}$ on the proton register (only $L=0,1$ contribute for spin-1/2 proton by parity). The convolution with the electron contact density reduces — at leading order in $r_Z$ — to a single scalar matrix element scaled by $r_Z$. Pauli weight is small ($\sim L^2$ per multipole). The catch: this operator encodes $r_Z$ as a *coefficient*, not as something derived from the proton wavefunction — structurally a Wigner-D-rotation-like multiplicative operator with $r_Z$ pulled in from outside. It implements W1b but does not close it without external $r_Z$.

### 4b. Operator-valued $\hat{R}_n$ on the HO register (medium cost)

The Track NI proton register encodes proton spatial structure in a 3D harmonic oscillator basis $|n_r, l, m_l\rangle$ at $\hbar\omega \sim$ MeV. The proton spatial coordinate operator $\hat{\mathbf{R}}_p$ on this register is well-defined: its matrix elements come from Moshinsky–Talmi brackets (already in `geovac/nuclear/moshinsky.py`). A magnetization-density operator can then be written as

$$\hat{m}(\mathbf{r}) = \mu_p\, \rho_M(\mathbf{r} - \hat{\mathbf{R}}_p)$$

where $\rho_M$ is a fixed radial profile (e.g. dipole form factor, Gaussian, or the lattice-QCD-fit Bernauer parametrization). The hyperfine operator is then the convolution of $\hat{m}$ with the electron contact density $|\psi_{1s}(0)|^2$ — but with $\hat{\mathbf{R}}_p$ now operator-valued on the nuclear register.

The cost: a multipole expansion of $|r_e - \hat{R}_p|^{-3}$ where one argument is operator-valued. This is structurally identical to the W1a operator construction (Pachucki–Patkóš–Yerokhin recoil) — the *same* infrastructure (Moshinsky–Talmi brackets, multipole expansion across a quantum nucleus) is needed. The Pauli weight grows as $O(Q_\text{nuc} \cdot Q_\text{elec})$ per multipole where $Q_\text{nuc}$ is the proton register size; for Track NI's 16-qubit nuclear register, this is $\sim O(160)$ Pauli strings per multipole, which the framework can handle.

The remaining input is $\rho_M(r)$ as a function (the radial profile). At parts-per-million precision, a Gaussian or dipole form-factor profile fits the data to sub-percent and only its finite moments matter — so the input is again a finite list of QCD-internal scalars.

### 4c. Spectral-triple promotion of the proton register (expensive)

The cleanest construction upgrades the Track NI proton register to a *full spectral triple* $(\mathcal{A}_p, \mathcal{H}_p, D_p)$. Then the Sprint H1 AC machinery applies: $\mathcal{T}_\text{atom} = \mathcal{T}_\text{elec} \otimes \mathcal{T}_p$, $D = D_e \otimes I + \gamma_e \otimes D_p$. The hyperfine coupling is an inner fluctuation $\omega \in \Omega_D^1$, and the magnetization operator is one component of the off-diagonal $D_F$-like block. KO-dim arithmetic and J construction generalize directly. **This is the same construction as the proposed A8'-for-W1a Pachucki port** (Phase A Candidate 1) — once built, W1a (color-electric) and W1b (color-magnetic) close as different inner-fluctuation components on the same composed triple. Layer 4c is no more expensive than 4b at the qubit level — the spectral-triple structure is organizational. But it generalizes cleanly to higher moments ($\langle r^3 \rangle_{(2)}$, $\Delta_\text{pol}$) as additional inner-fluctuation components.

**Q3 answer:** the natural operator class is a **multipole expansion of the proton magnetization density across registers, with the proton coordinate $\hat{\mathbf{R}}_p$ promoted to an operator** (layer 4b), ideally inside a spectral-triple-promoted proton register (layer 4c). Pauli weight is modest ($\sim 100$–$200$ strings per multipole at Track NI's 16-qubit nuclear register); no new qubits are needed. The radial profile $\rho_M(r)$ enters as Layer-2 calibration data — structurally a *finite list of moments* per Q1 — not as a fundamental input.

The construction is *very nearly the same* as the W1a operator construction. The difference is the cross-coupling: W1a wants $-Z/|r_e - \hat{R}_p|$ (color-electric, single moment $R_p^2$); W1b wants the convolution $\rho_E \star \rho_M$ in the contact term (color-magnetic, single moment $r_Z$). The infrastructure is shared.

---

## 5. Q4 — Reducible to W1a in disguise?

This is the load-bearing question for Phase C scope. The argument both ways:

### Argument FOR reducibility

Once the operator-valued $\hat{\mathbf{R}}_p$ exists on the Track NI proton register (the W1a closure object), the magnetization-density convolution $\rho_E(r_e) \star \rho_M(\hat{R}_p)$ is a moment expansion of that operator: $H_\text{hf}^\text{full} = A_\text{hf}^\text{point} (1 - 2 Z \alpha m_e \langle \hat{r}_Z \rangle + \mathcal{O}(\alpha^2 r_Z^2))$ where $\langle \hat{r}_Z \rangle$ is a scalar matrix element of an operator built from $\hat{\mathbf{R}}_p$. If the W1a sprint closes the cross-register two-body Coulomb via multipole expansion, *the same multipole machinery applies to the convolution* — only the radial weight (Gaussian/dipole vs Coulomb) and the angular structure (couples to electron spin vs orbital) differ. Wigner-3j, Moshinsky–Talmi, cross-register tensor product (Tool 8) are all unchanged. W1b is then a *calibration moment integral* on top of W1a infrastructure: the operator structure is identical, only the input data differs.

### Argument AGAINST reducibility

The magnetization distribution $\rho_M$ is *physically distinct* from the charge distribution $\rho_E$: $\rho_E$ is the quark color-electric density; $\rho_M$ is the quark color-magnetic-moment density. They are independently parameterized in form-factor extractions ($G_E^p$ vs $G_M^p$), have different $Q^2$-dependence, and $r_Z$ is a *convolution* not derivable from either alone. Three categorical distinctions matter: (1) two independent QCD inputs — W1a needs charge-distribution moments only; W1b additionally needs magnetization-distribution moments. (2) The full $\rho_M$ profile is not redundant with $\rho_E$; lattice QCD computes $G_M^p(Q^2)$ and $G_E^p(Q^2)$ separately, with separate renormalization and polarizability contributions. (3) Higher moments differ: at sub-ppm the Friar moment $\langle r^3 \rangle_{(2)}$ depends on the *correlation* between $\rho_E$ and $\rho_M$ beyond their individual radii. W1a closure is therefore *necessary but not sufficient* for W1b: the same infrastructure suffices, but additional Layer-2 calibration is required.

### Reconciliation

**Operator infrastructure is shared.** Both walls require the same A8'-class spectral-triple promotion of the proton register, the same Moshinsky–Talmi multipole expansion across the operator-valued $\hat{R}_p$, the same cross-register tensor product, the same Pauli-string assembly. Once this infrastructure exists, both close architecturally together.

**Input data is distinct.** W1a needs $\langle r^2 \rangle_E$; W1b needs $r_Z$ (charge × magnetization convolution moment) plus the radial profiles for higher moments. W1b is "W1a + one extra Layer-2 calibration scalar."

Strict answer: **structurally yes (operator infrastructure shared), calibration-wise no (additional QCD-input scalar required).** Phase C verdict: closing W1a *does* close W1b architecturally, but the Phase C-W1a deliverable should explicitly include the magnetization-density operator as a second inner-fluctuation component on the same composed triple. This is exactly how the Connes SM construction handles multi-channel inner fluctuations ($\omega_\text{gauge}$ and $\omega_\text{Higgs}$ as separate components on the same triple, classified cleanly in Sprint H1). W1a's recoil and W1b's Zemach are analogous: two inner-fluctuation components on the same Track-NI-promoted composed atomic spectral triple.

---

## 6. Q5 — Verdict

The four-option verdict (a/b/c/d) and the reasoning:

- **(a) Tooling-addressable, Wigner-D + multipole + tensor-product machinery sufficient:** *True for the operator construction* (layer 4a or 4b, §4), but the input data $r_Z$ is not generated by the framework. Pure tooling answer is incomplete.

- **(b) Downstream of W1a — closing W1a closes W1b for free, modulo a calibration moment integral:** **This is the correct verdict.** Per Q4 reconciliation, the operator infrastructure is shared and the Phase C-W1a sprint can deliver both walls' architectural closure simultaneously, *provided* the sprint scope explicitly includes the magnetization-density inner-fluctuation component and adds $r_Z$ to the Layer-2 input list. The closure is not free in absolute terms (an extra QCD-input scalar plus an extra inner-fluctuation component), but it is free *in incremental sprint cost beyond W1a* — the additional engineering is small relative to the W1a closure itself.

- **(c) Genuinely independent of W1a, requires distinct infrastructure (e.g. field-theoretic input from QCD):** *False.* The QCD input ($r_Z$) is the same kind of Layer-2 calibration scalar as the proton charge radius (already in `R_PROTON_BOHR`). No new infrastructure beyond the W1a operator promotion is needed.

- **(d) Calibration-side, the framework will always need $r_Z$ or $\rho_M$ as Layer-2 input:** *Partially true* — the *value* of $r_Z$ is and will remain a Layer-2 calibration input, since the framework's Wilson-SU(3) construction (Paper 30 + Sprint ST-SU3) admits the gauge content but cannot autonomously compute matter elements (Sprint 5 Track S5 CG obstruction). However, the *operator structure* that consumes $r_Z$ is not calibration-side; it is a missing GeoVac-internal architectural piece that closes structurally with W1a.

**Final verdict: (b) downstream of W1a, with explicit Phase C scope addition to declare W1b's Layer-2 input ($r_Z$, optionally $\Delta_\text{pol}$ at sub-ppm precision) and include the magnetization-density inner-fluctuation component in the W1a operator construction.**

---

## 7. Phase C-W1a sprint scope must include W1b

Given the (b) verdict, the Phase C-W1a sprint's scope statement should explicitly include the W1b architectural closure. The augmented scope is:

**Phase C-W1a (with W1b inclusion):**

1. **Promote the Track NI proton register to a spectral triple** $(\mathcal{A}_p, \mathcal{H}_p, D_p)$ with KO-dim 3 (Riemannian-spinor, matching the GeoVac sector). Compose with the electronic GeoVac sector via $\mathcal{T}_\text{atom} = \mathcal{T}_\text{GV} \otimes \mathcal{T}_p$. KO-dim arithmetic: $3 + 3 = 6 \equiv 6 \pmod 8$, so $J^2 = +I$ on the composed triple — a textbook Riemannian tensor-product result.
2. **Build the cross-register two-body coordinate operator** $V_{eN}(r_e, \hat{R}_p) = -Z/|r_e - \hat{R}_p|$ via multipole expansion across the operator-valued $\hat{R}_p$ on the proton HO register. Verify multipole termination at $L_\text{max} = 2 \max(l_e, l_N)$ across mismatched exponents (the Phase B-W1a-diag check, also gating).
3. **Build the magnetization-density operator** as a *second inner-fluctuation component* on the same composed triple. The radial profile $\rho_M(r)$ enters as a fixed parameterization (dipole form factor at first pass; Bernauer–Mainz 2014 fit at second pass). The convolution with the electron contact density gives the Zemach correction with $r_Z$ as the leading scalar Layer-2 input.
4. **Calibration target.** W1a closes against Pachucki–Patkóš–Yerokhin 2023 PRL recoil structure (residual at hydrogen 1S). W1b closes against Eides Tab. 7.3 / Karshenboim 2005 Zemach structure (residual at 21 cm). The single Phase C-W1a sprint deliverable validates both at once.
5. **Layer-2 input list.** $\{R_p, r_Z, [\Delta_\text{pol}], [m_p/m_e]\}$ — all CODATA / Eides values, no fits.
6. **Verdict expectation.** W1a + W1b both close at architectural level; residual on the 21 cm line drops from +18 ppm (HF-4) to a budget consistent with multi-loop QED + nuclear polarizability ($\sim$ +12 to +18 ppm per Eides Tab. 7.3) — i.e. W2a-territory residuals, not W1-territory residuals.

The incremental cost of including W1b on top of W1a is small: one additional inner-fluctuation component, one additional Layer-2 scalar in the input list, and one additional calibration target (21 cm hyperfine shift) on top of the W1a hydrogen-1S shift. The Pauli-string-budget addition is a single multipole expansion of comparable weight to the recoil V_ne, $\sim$ 100–200 additional cross-register Pauli terms.

**Phase B-W1b-diag's gating decision:** Phase C-W1b should NOT be a separate sprint. It should be folded into Phase C-W1a as an explicit deliverable. The synthesis memo's "B-W1b-diag" → "Phase C-W1b" gating arrow in §7 of the Phase A synthesis is therefore corrected: the gate is to Phase C-W1a-augmented, not to a separate Phase C-W1b.

---

## 8. Honest scope and uncertainty

Covered: the Track NI cross-register architecture (Paper 23 §VI; `geovac/nuclear/nuclear_electronic.py` + `form_factor.py`); the Eides–Grotch–Shelyuto + Karshenboim hyperfine framework as the calibration neighborhood; operator-construction sketch via multipole expansion + spectral-triple promotion; the reducibility analysis to W1a.

Not covered deeply: (1) sub-ppm precision targets ($\Delta_\text{pol}$, recoil-corrected Zemach, Friar moment $\langle r^3 \rangle_{(2)}$ each enter as additional inner-fluctuation components — downstream, out of scope at ppm-level closure); (2) the exact form-factor parameterization for $\rho_M$ — choice among dipole / Galster / Bernauer / Lin–Hammer–Meißner shifts $r_Z$ within current uncertainty ($\sim 1$%); (3) muonic-hydrogen vs hydrogen $r_Z$ extraction (sub-percent shift, sub-ppm effect on 21 cm); (4) explicit KO-dim arithmetic on the spectral-triple-promoted proton register — quoted textbook $3+3 = 6 \pmod 8$ but verifying this for the specific Camporesi–Higuchi spinor + HO-basis Riemannian Dirac on the proton is a 1–2 week Phase C check; if KO-dim forces antimatter doubling à la Sprint H1's $J_F^2 = -I$, construction generalizes with larger Hilbert space.

Confidence: **high** on the (b) verdict (operator infrastructure shared, separate calibration scalar). **Medium-high** on the Phase C scope (operator construction works; KO-dim arithmetic on promoted proton register needs explicit verification). **Medium** on calibration-target detail (form-factor parameterization shifts $r_Z$ within uncertainty at sub-ppm but not ppm).

---

## 9. Summary for synthesis

**Verdict: (b) — W1b is downstream of W1a.**

The operator infrastructure is shared between W1a and W1b: both walls require promoting the Track NI proton register to a full spectral triple, building the cross-register two-body operator via multipole expansion across an operator-valued $\hat{R}_p$, and assembling the result as cross-register Pauli strings. Once that machinery exists, W1b adds a single magnetization-density inner-fluctuation component (small incremental cost) and one Layer-2 calibration scalar ($r_Z$ from QCD/scattering data, Eides 2024 value $r_Z = 1.045$ fm).

The Phase C-W1a sprint scope should be augmented to explicitly include W1b's magnetization-density inner-fluctuation component as a deliverable, with $r_Z$ declared in the Layer-2 input list alongside $R_p$. No separate Phase C-W1b sprint is needed. The single augmented sprint closes both walls architecturally and lands the 21 cm residual in the multi-loop QED + nuclear polarizability budget ($\sim$ +12 to +18 ppm per Eides Tab. 7.3) — exactly where the residual then properly belongs to W2a (multi-loop renormalization), not to W1.

The calibration-side input ($r_Z$) is and will remain a Layer-2 scalar from QCD/scattering data, since the Wilson-SU(3) construction (Paper 30 + Sprint ST-SU3) admits the gauge content but cannot autonomously compute matter elements (Sprint 5 Track S5 CG obstruction). This is consistent with GeoVac's structural-skeleton scope statement and with the Paper 18 sixth tier (inner-factor / Layer-2 input data).

The 16th projection in Paper 34's living catalogue would naturally be the magnetization-distribution Layer-2 focal length, sibling to but categorically distinct from the rest-mass / charge-radius projections (14th, 15th). Updating Paper 34 is a small additional task, optional and gated on the Phase C-W1a-augmented sprint outcome.

---

**End of Phase B W1b diagnostic memo.**
