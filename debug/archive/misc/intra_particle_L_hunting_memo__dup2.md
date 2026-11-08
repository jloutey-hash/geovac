# Intra-particle / intra-nuclear length-scale hunting memo

**Date:** 2026-05-09
**Source:** post-§IV.5 audit (`debug/paper34_v_base_unit_audit_memo.md`); follow-up to §F.1 magnetization-density 17th-projection candidate.
**Question.** §IV.5 currently shows [L] injected by four projections (Fock $a_0$, stereographic $r$, Wigner D $\vec{R}_{AB}$, mol-frame $R$). The audit flagged magnetization-density / Zemach $r_Z$ as a fifth [L] candidate. **Are there other intra-particle / intra-nuclear length scales that deserve projection status, and do they all fit one mechanism or multiple?**

This memo is a research catalogue. No paper edits. Tone: working-hypothesis register on structural readings; sub-percent literature values cited where current.

---

## 1. Catalogue of candidate intra-particle / intra-nuclear [L] scales

| # | Scale | Symbol | Current value (literature) | Physics role | Hamiltonian term modified | Mechanism | GeoVac treatment as of 2026-05-09 | Candidate slot |
|:-:|:------|:-------|:--------------------------:|:-------------|:--------------------------|:----------|:---------------------------------|:---------------|
| 1 | Proton charge radius | $r_p$ | 0.84075(64) fm (CODATA 2022) | Charge-distribution correction to electronic Coulomb potential | $V_\text{ne}(r) \to \int d^3r' \rho_E(r')/|\vec{r}-\vec{r}'|$; Foldy/Friar correction $+\tfrac{2\pi}{3}Z\alpha\langle r^2\rangle_E\,\delta^3(r)$ | Charge-density convolution — modifies $V_{ne}$ | External scalar via `R_PROTON_BOHR` (point-source); no operator | (a) charge-density projection |
| 2 | Proton Zemach radius | $r_Z$ | 1.045(20) fm Eides 2024; 1.013(10)(12) fm lattice 2024 | Magnetization-distribution correction to hyperfine | $A_\text{hf}^\text{point} \to A_\text{hf}^\text{point}(1 - 2Z\alpha m_e r_Z)$ | Magnetization-density convolution — modifies hyperfine contact | `geovac/magnetization_density.py` ✓ operator-level (HF-4, μH) | (b) magnetization-density projection ← §F.1 candidate |
| 3 | Proton magnetic radius | $r_M$ | 0.776–0.851 fm (Mainz/world) | Magnetic form-factor slope; controls Zemach moment via $\rho_M$ profile | Same hyperfine term as $r_Z$, but at sub-leading order | Profile shape parameter — distinguishes Gaussian vs exponential vs dipole $\rho_M(r)$ | Sub-feature of (b): different $r_Z$-calibrated profile families | sub-feature of (b) |
| 4 | Proton Friar moment | $\langle r^3\rangle_{(2)}$ | $\sim 2.85\,\text{fm}^3$ (Friar–Payne) | Sub-leading recoil-corrected Zemach in muonic atoms; controls $O((r_Z/a_0)^2)$ correction | Sub-leading hyperfine + Lamb-shift contact term | Higher moment of $\rho_E \star \rho_M$ — independent of $r_Z$ in general | $M_2$ moment computed in `magnetization_density.py` but only as profile-dependent diagnostic | sub-feature of (b)+(a) |
| 5 | Deuteron charge radius | $r_d$ | 2.1415(45) fm (atomic D); 2.1256 fm (μD) | Same as $r_p$ but for deuteron | Same as #1 | Same as #1 | External scalar via `R_DEUTERON_BOHR` (point source) | (a), one species at a time |
| 6 | Deuteron Zemach radius | $r_Z(D)$ | 2.593(16) fm (Friar–Payne 2005) | Same as $r_Z$ but for deuteron | Same as #2 | Same as #2 | `magnetization_density.py` already accepts $r_Z(D)$ as input (D HFS row) | (b), one species at a time |
| 7 | Deuteron quadrupole moment | $Q_d$ | 0.285699(15)(18) fm² | Tensorial coupling of nuclear charge distribution to electronic field gradient | $H_Q = -e\,Q_d\,T^{(2)}_{ij}\,(\partial_i E_j)/6$, an $L\geq 1$ multipole; couples to electronic $l\geq 1$ gradients | Tensor (rank-2) charge-distribution coupling — categorically different from $r_p$ scalar | None; flagged in Paper 23 as "sub-ppm for s-states" | (c) tensor-distribution projection |
| 8 | NN scattering length (singlet) | $a_s$ | $-23.7153(43)$ fm | Asymptotic singlet $^1S_0$ NN potential calibration | NN potential $V_\text{NN}(r)$ as $r\to\infty$ | Asymptotic-potential parameter — sets effective interaction range, not a distribution moment | None; nuclear sector uses Minnesota Gaussian (Paper 23) which is an explicit potential | (d) asymptotic-potential / effective-range tier (not a distribution length) |
| 9 | NN scattering length (triplet) | $a_t$ | 5.4114(15) fm | Same as #8 but $^3S_1$ | Same as #8 | Same as #8 | None | (d) |
| 10 | NN effective range (s/t) | $r_s, r_t$ | 2.665(56) fm; 1.7468(19) fm | Effective-range expansion correction to scattering length | Same as #8 | Same as #8 (next-order coefficient) | None | (d) |
| 11 | Pion Compton wavelength | $1/m_\pi$ | $\sim 1.41$ fm | Natural [L] of strong interaction; Yukawa range | Implicit in any Yukawa-tail NN potential | Compton wavelength — a length set by an external mass scale (chiral EFT) | None; framework has no native $\pi$-meson-exchange operator | (e) Compton-wavelength tier (length from an external mass projection) |
| 12 | Nuclear polarizability | $\Delta_\text{pol}$ | $\sim$ +1–2 ppm hydrogen, $\sim$ +44 ppm deuteron HFS | Inelastic nuclear-excitation contribution to hyperfine | Sub-leading addition to hyperfine, indexed by an inelastic structure function $S_1(\nu, Q^2)$ | Inelastic-channel contribution — has dimensional analog of length but lives in a non-local inelastic kernel | None; framework has no excited-nucleon spectrum | (f) inelastic-channel tier (NOT a length scale; included for completeness) |
| 13 | Neutron skin / charge radius | $r_n^2$ | $-0.1161(22)$ fm² (negative!) | Asymmetry between proton and neutron charge distributions in nuclei | Neutron contribution to nuclear charge distribution | Charge-density of neutron component | None; first-row + frozen-core treatments are point-nucleus | sub-feature of (a) for non-trivial $A>1$ nuclei |

**Twelve scales catalogued.** Items 1–7 are intra-nuclear [L] scales of distribution type (charge / magnetization / tensor). Items 8–11 are interaction-range or external-mass-derived lengths. Item 12 is included to flag that nuclear polarizability is *not* a length scale despite often being budgeted alongside Zemach contributions. Item 13 is a sub-feature for heavier nuclei.

---

## 2. Structural mechanism analysis

The catalogue partitions cleanly into four mechanism classes plus one special case.

**Class (a): scalar charge-density convolution.** $r_p$, $r_d$, $r_n$. Modifies $V_{ne}$ via $\int d^3r' \rho_E(r')/|\vec{r}-\vec{r}'|$. At leading order the only effect on s-state energies is the Foldy/Friar contact shift $\propto \langle r^2\rangle_E\,|\psi(0)|^2$. The Hamiltonian term modified is the **electron–nucleus Coulomb operator**.

**Class (b): scalar magnetization-density convolution.** $r_Z$, $r_Z(D)$, $r_M$ (sub-feature), Friar moment $\langle r^3\rangle_{(2)}$ (sub-feature). Modifies $A_\text{hf}^\text{point}$ via $\int d^3r\,d^3r'\,\rho_E(r')\,\rho_M(|\vec{r}-\vec{r}'|)$. The Hamiltonian term modified is the **Fermi-contact spin–spin coupling**. This class has the operator-level construction in `magnetization_density.py`.

**Class (c): tensor (rank-≥2) distribution.** $Q_d$, magnetic dipole moments at higher order, octupole moments for spin-≥3/2 nuclei. Modifies *gradient* couplings; couples only to electronic $l \geq 1$ for $Q$-class operators. The Hamiltonian term modified is a **multipole expansion of the nuclear potential at $L \geq 2$**. Categorically distinct from (a) and (b): not a scalar distribution moment, but the leading coefficient of a tensor-coupled multipole.

**Class (d): asymptotic-potential parameter.** $a_s$, $a_t$, $r_s$, $r_t$. NOT a distribution moment of any nuclear density. Modifies the **NN interaction itself**, not an electronic Hamiltonian term. Lives in the nuclear sector (Paper 23). The Coulomb-projection foundational geometry of GeoVac (Fock + S^3) does not host this — it is a property of the strong potential, which the framework currently embeds via Minnesota Gaussians at a fixed scale.

**Class (e): Compton-wavelength length from a mass scale.** $1/m_\pi$ and analogs ($1/m_K$, $1/m_W$, $\hbar/(m_e c)$ which is already implicit in $a_0$). Length derived directly from a particle mass via the Compton relation. *Not an independent projection*: it is the rest-mass projection (§III.14) reread in length units. The framework already has this — it's the same axis as $m$.

**Special case (f): inelastic-channel tier.** $\Delta_\text{pol}$. Not a length; an inelastic structure function. Listed for completeness because it competes with $r_Z$ in precision budgets. Belongs in Paper 18's inner-factor input tier (§IV.6), not §III as a projection.

### The structural test (one mechanism vs many)

The audit memo's question — "are these all the same mechanism with sub-features, or categorically different?" — has a clean answer:

- **Within class (a)**: all sub-features of one mechanism (different species, same convolution). A *charge-density projection* with one variable $r_E$ per nuclear species suffices.
- **Within class (b)**: all sub-features of one mechanism. Currently realized as `magnetization_density.py`. A *magnetization-density projection* with one variable $r_Z$ per species (and $r_M$, Friar moment as profile-shape sub-features) suffices.
- **Across (a) and (b)**: **categorically distinct**. They modify different Hamiltonian terms (Coulomb vs Fermi-contact); they couple to different electronic operators (density at origin vs spin-density at origin); they have different selection rules (Foldy is l=0 only via $\delta^3(r)$, Zemach is l=0 only via the same delta but on the *spin-spin* sector). Empirical anchor: hydrogen 21cm depends almost entirely on $r_Z$ and trivially on $r_p$, while electronic Lamb-shift in hydrogen depends primarily on $r_p$ and trivially on $r_Z$. The mechanisms cannot be unified by a single profile-broadening prescription.
- **Class (c) vs (a),(b)**: **categorically distinct**. Tensor coupling to a gradient is structurally different from scalar coupling to a contact density. Different selection rules (e.g., $Q_d$ requires $l\geq 1$ on the electronic side).
- **Class (d)**: **categorically distinct from all electronic-side classes**. Nuclear-internal interaction parameter, not an electronic Hamiltonian modification.

So the framework needs *at least three* §III entries on the electronic-side intra-nuclear axis — (a), (b), (c) — not one with sub-features.

---

## 3. Structural relations among the scales

QCD provides several relations the framework should respect or test. Each relation is a *structural prediction* that, if it holds in the framework, would compress the parameter count.

**(R1) $r_p^2 \neq r_Z^2$ but they share a chiral-perturbation expansion.** $r_p^2 = -6\,\text{d}G_E/\text{d}Q^2|_0$; $r_M^2 = -(6/G_M(0))\,\text{d}G_M/\text{d}Q^2|_0$; $r_Z = \int d^3r\,d^3r'\,\rho_E(r')\rho_M(|\vec{r}-\vec{r}'|)$. The Zemach radius is a *convolution moment* of the two distributions, not a fundamental input. The naive Friar identity $r_Z \approx r_p + r_M\cdot\text{shape factor}$ is only approximate. **The framework currently treats $r_p$, $r_Z$, $r_M$ as three independent scalar inputs.** A genuine QCD-internal model would express $r_Z$ as a derived quantity; chiral lattice (Phys. Rev. D 110 L011503, 2024) does this, but their result $r_Z = 1.013(10)(12)$ fm is itself an input from QCD lattice computation. Within GeoVac, these stay as Layer-2 calibration constants — three inputs, not two.

**(R2) Tensor moments and scalar moments are independent.** The Wigner–Eckart theorem says they cannot mix at leading order: $\langle r^2\rangle_E$ (scalar) and $Q_d$ (tensor rank-2) have orthogonal selection rules. So adding $Q_d$ as a class-(c) projection does not duplicate the class-(a) input.

**(R3) Asymptotic-potential parameters live in the strong sector.** $a_s$ vs $r_p$ are not related by any chiral relation that the framework can tap into without a Yang–Mills-with-quarks input. Class (d) sits structurally outside the electronic-side projection family.

**(R4) Compton-wavelength scales fold into rest-mass.** $1/m_\pi$ is just rest-mass projection (§III.14) restricted to the pion species, then expressed in length units via $\hbar c$. No new projection needed; inherited.

**(R5) Nuclear polarizability is non-local in the radial coordinate.** $\Delta_\text{pol}$ requires an inelastic structure function $S_1(\nu, Q^2)$, which is a function of two variables not a length. It does *not* fit into any §III projection as currently structured; it lives in Paper 18 §IV.6's inner-factor input tier.

---

## 4. Recommendation

**Recommended count: three §III electronic-side intra-nuclear projections (a), (b), (c), distinct from each other; one §IV.6 inner-factor inelastic entry (f); class (d) is nuclear-sector and out of §III scope.**

| Slot | Class | Variables | Dim | Mechanism | Hamiltonian term |
|:-----|:------|:---------:|:---:|:----------|:-----------------|
| §III.17 candidate | (a) charge-density | $r_E$ per nucleus species | $[L]$ | Convolution $V_{ne}\to V_{ne}\star\rho_E$ | $V_{ne}$ |
| §III.18 candidate | (b) magnetization-density | $r_Z$ per nucleus species | $[L]$ | Convolution $A_\text{hf}^\text{contact}\to(\rho_E\star\rho_M)$ | Fermi-contact $A_\text{hf}\,\hat{I}\cdot\hat{S}\,\delta^3(r)$ |
| §III.19 candidate | (c) tensor-distribution | $Q_N$ per nucleus species (rank-2 leading; higher rank as needed) | $[L]^2$ | Multipole coupling to $\partial_i E_j$ / gradients | Quadrupole $H_Q$, octupole $H_O$, ... |
| Paper 18 §IV.6 inner factor | (f) inelastic | $\Delta_\text{pol}$ per nucleus | dimensionless ppm-class | Inelastic structure function | adds to total hyperfine |

This is a *minimal* taxonomy. It is justified by the empirical separation:

- HF-4 anchor ($r_Z$) lives in (b), not (a) or (c).
- The Foldy / Friar correction to the electronic Lamb shift (currently external in Paper 36 LS-1..LS-4) lives in (a). When the framework adds an operator-level $V_{ne}\to\int d^3r' \rho_E(r')/|\vec{r}-\vec{r}'|$ (the natural sibling of `magnetization_density.py`), it will need a (a) projection slot.
- The deuteron quadrupole $Q_d$ (currently flagged as "sub-ppm for s-states" in Paper 23 §VIII) will need (c) when the catalogue extends to molecular HFS observables (HD/D2 rotational HFS, where electronic gradients at the nucleus are nonzero and Q_d couples at percent level — already measured at sub-percent precision per Komasa–Pachucki 2020).

### Falsifier statement

A clean negative for "three projections" would be: every intra-nuclear [L] mechanism on the catalogue folds into one universal "nuclear-distribution focal length" projection, where the Hamiltonian term modified depends on a tensor index and the framework computes it once. **This is structurally false** because the Hamiltonian terms are different operators: $V_{ne}$ and $A_\text{hf}\,\hat{I}\cdot\hat{S}$ are not related by any tensor-rank rotation; they are independent inner fluctuations of the spectral triple per the Connes-Chamseddine $\omega = \omega_\text{gauge} + \omega_\text{Higgs}$ analog Paper 23 §VII.4 invokes for $\omega_\text{recoil} + \omega_\text{magn}$.

A clean positive — and the recommended action — is exactly what the catalogue shows: three electronic-side classes with categorically distinct Hamiltonian targets, one nuclear-side class outside §III scope, and one inelastic class in the inner-factor input tier.

---

## 5. Honest scope: what the framework cannot represent

- **NN scattering length / effective range (class d).** The framework's natural geometry on the nuclear side is the Bargmann–Segal lattice (Paper 24) or the HO + Minnesota Gaussian (Paper 23). Neither hosts a parameter of "asymptotic NN potential" type. The Minnesota potential's two Gaussian widths are external scalar inputs in `geovac/nuclear/`, not derived. An explicit chiral-EFT extension of Paper 23 — replacing Minnesota with a scattering-length-calibrated low-energy effective potential — would add this class as a fifth Layer-2 input tier in Paper 18 §IV. Two-body scattering observables themselves require continuum-spectrum operators (asymptotic plane-wave states) that the Fock-projection geometry does not natively host.

- **Pion-exchange Yukawa range (class e).** $1/m_\pi$ is downstream of pion rest-mass projection. Requires a pion-meson species in the framework's particle content. None exists. Adding one would be the same architectural extension as the chiral-EFT direction above.

- **Inelastic structure functions (class f).** $\Delta_\text{pol}$ is a non-local inelastic kernel, not a length. Lives in Paper 18 §IV.6 inner-factor tier.

- **Lattice-QCD-derived radii.** $r_Z = 1.013(10)(12)$ fm from Phys. Rev. D 110 L011503 (2024) is a *result* the framework cannot reproduce; it is consumed as Layer-2 input.

These are not framework defects — they are precise statements of where the GeoVac structural skeleton ends and external QCD/strong-interaction calibration begins. The same pattern as the multi-focal-composition wall (CLAUDE.md §3): the framework cleanly handles selection rules, transcendental signatures, scaling laws; it does not autonomously generate calibration data.

### Final recommendation

Promote three projections — charge-density (a), magnetization-density (b), tensor-distribution (c) — as §III.17–§III.19 candidates. Do *not* unify them into one projection. The audit memo §F.1's single "magnetization-density projection" candidate handles only (b); class (a) is a distinct sibling with its own operator-level analog needed when the framework eventually constructs a `geovac/charge_density.py` module to handle Foldy / Friar contact corrections to $V_{ne}$ at the operator level. Class (c) is a third sibling for tensor couplings, currently a forward-looking slot.

**Files referenced:** `debug/paper34_v_base_unit_audit_memo.md` §F.1; `geovac/magnetization_density.py`; `papers/group6_precision_observations/paper_34_projection_taxonomy.tex` §III.16 (Breit retardation as the most recent precedent for adding a projection); `papers/group4_quantum_computing/paper_23_nuclear_shell.tex` §VII.4–VIII (operator-level Zemach + deuterium HFS).

**Sources for current values:**
- Proton charge radius: [CODATA 2022 (NIST)](https://physics.nist.gov/cgi-bin/cuu/Value?rp), [CODATA 2022 arxiv](https://arxiv.org/html/2409.03787v1) — $r_p = 0.84075(64)$ fm
- Proton Zemach radius: [Lattice QCD 2024](https://arxiv.org/abs/2309.17232) — $r_Z = 1.013(10)(12)$ fm; Eides 2024 atomic value 1.045 fm
- Proton magnetic radius: [Lin–Hammer–Meissner 2020 (arXiv:2002.05167)](https://arxiv.org/abs/2002.05167); [Mainz A1 2010](https://arxiv.org/abs/1505.01489) — $r_M \in [0.776, 0.851]$ fm
- Deuteron quadrupole moment: [Komasa–Pachucki PRL 2020 (arXiv:2010.06888)](https://arxiv.org/abs/2010.06888) — $Q_d = 0.285699(15)(18)$ fm²
- Deuteron charge radius: [Pohl et al. Metrologia 2017 (arXiv:1607.03165)](https://arxiv.org/abs/1607.03165) — $r_d = 2.1415(45)$ fm
- np scattering parameters: [Babenko–Petrov 2007 (arXiv:0704.1024)](https://arxiv.org/abs/0704.1024) — $a_t = 5.4114(15)$ fm, $a_s = -23.7153(43)$ fm
- Deuteron Zemach radius: [Friar–Payne 2005](https://www.sciencedirect.com/science/article/pii/S0370269303017416) — $r_Z(D) = 2.593(16)$ fm
- Nuclear polarizability ($\Delta_\text{pol}$): [Pachucki et al. PRA 107 (2023)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.107.052802); [Pachucki 2025 review](https://arxiv.org/html/2506.08879)
