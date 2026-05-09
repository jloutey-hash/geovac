# External uniform fields in Paper 34: do Zeeman and Stark warrant new §III projections?

**Date:** 2026-05-09
**Context:** post-base-unit-audit hunting pass, per the audit memo's closing section. The audit identified three loose threads; this memo follows up on thread (3): **external uniform fields as candidate new projections**, focused on uniform magnetic field $B$ (Zeeman regime) and uniform electric field $E$ (Stark regime). The check is structural, not a sprint: do existing §III projections absorb each, or is a new entry required?

## 1. Dimensional analysis in atomic units

Atomic units: $\hbar = m_e = e = a_0 = 1$, energy in Hartree (Ry), and $c = 1/\alpha \approx 137$. Charge $[Q]$ is the integer label, not a continuous scale.

**Magnetic field $B$.** The natural objects associated with $B$ are:

- *Cyclotron frequency:* $\omega_c = eB/(m_e c)$, dimension $[E]$ in atomic units. A uniform $B$-field induces Landau quantization $E_n = \hbar\omega_c (n + \tfrac{1}{2})$.
- *Magnetic length:* $\ell_B = \sqrt{\hbar/(eB)}$, dimension $[L]$. Sets the scale of cyclotron orbits.
- *Bohr-magneton energy:* $\mu_B B = (e\hbar/2m_e c) B$, dimension $[E]$. Linear Zeeman coupling to spin/orbital angular momentum.

The two scales are conjugate: $m_e \omega_c \ell_B^2 = \hbar$ exactly, so $\omega_c \ell_B^2 = 1$ in atomic units. They are not independent — specifying $B$ fixes both. **Either $\omega_c$ ($[E]$) or $\ell_B$ ($[L]$) can serve as the projection's primary variable; the other is its conjugate.**

This is structurally parallel to the rest-mass projection (variable $m$, dimension $[M]$, but at $c = 1$ the inverse Compton length is the conjugate $[L]$ scale), and to the observation/temporal-window projection (variable $\beta$ in $[T]$, equivalent to $T = 1/\beta$ in $[E]$). In each case, a single physical variable provides one base-unit injection with a fixed conjugate.

**Electric field $E$.** The natural object is the Stark Hamiltonian $H_\text{Stark} = -e\,\vec{\mathcal{E}}\cdot\vec{r}$ (with $\vec{\mathcal{E}}$ the field strength to avoid notation collision with energy). In atomic units $|\vec{\mathcal{E}}|$ has dimension $[E]/[L]$ (force per unit charge); the natural energy scale is $|\vec{\mathcal{E}}| \cdot a_0$, giving an $[E]$ injection. There is no characteristic *length* attached to a uniform $\vec{\mathcal{E}}$ in the same way $\ell_B$ is attached to $B$ — the Stark Hamiltonian is linear in $r$ everywhere, with no quantization.

**Verdict on dimensional analysis:** Both $B$ and $\vec{\mathcal{E}}$ inject a new $[E]$ scale (via $\omega_c$ or $|\vec{\mathcal{E}}|a_0$ respectively). $B$ additionally injects a conjugate $[L]$ scale ($\ell_B$); $\vec{\mathcal{E}}$ does not. By the §IV.5 base-unit accounting, $B$ and $\vec{\mathcal{E}}$ would each be a fifth $[E]$ injection (joining Fock, spectral action, Drake-Swainson, Breit retardation), with $B$ also being a fifth $[L]$ injection (joining Fock, stereographic, Wigner $D$, mol-frame).

## 2. Mechanism analysis against existing §III projections

### 2.1 Spectral action (#6) — does it absorb $B$?

Connes-Chamseddine spectral action evaluates $\mathrm{Tr}\, f(D^2/\Lambda^2)$ where $\Lambda$ is a UV cutoff. Coupling an external gauge field is well-defined inside this framework: replace $D \to D_A = D + A$ where $A$ is a connection 1-form. For a uniform external $B$, the relevant $A$ is the symmetric-gauge or Landau-gauge vector potential.

But this is **distinct from the §III.6 projection as currently named**. Spectral action injects the variable $\Lambda$ (UV cutoff scale) and $\alpha$ (gauge coupling, internal to the spectral triple). $\Lambda$ is an IR-to-UV regulator for *vacuum fluctuations*; $\omega_c$ is an IR scale for *external matter response*. They live at opposite ends of the energy axis. A uniform external $B$ is not absorbed by tuning $\Lambda$.

Could one read $B$ as entering through the *connection $A$* on the spectral triple, with no new projection needed? Yes formally — adding a U(1) gauge connection is the standard inner-fluctuation mechanism in NCG (Connes-Marcolli). But in §III.6 the gauge content is *internal* (via Wilson plaquettes, Paper 25/30), parameterized by $\beta_{\rm Wilson} = 1/g^2$ and dimensionless. The external uniform $B$ is a *fixed background field*, not a fluctuating gauge variable being summed over. **§III.6 does not naturally absorb external fields** — its variables are $\Lambda, \alpha$, neither of which is the strength of an external classical field.

Verdict: spectral action does NOT absorb $B$ or $\vec{\mathcal{E}}$.

### 2.2 Rest-mass projection (#14) — does it absorb $B$?

Rest-mass projection is $\omega_n^2 \to \omega_n^2 + m^2$ (KG-1). The Landau-level Hamiltonian on a uniform $B$ background is $H = \frac{1}{2m_e}(\vec{p} - e\vec{A}/c)^2$, whose spectrum on a 2D plane is $E_n = \hbar\omega_c(n+\tfrac{1}{2})$ — *not* a shift $\omega_n^2 \to \omega_n^2 + (eB)^2$. The mathematical structure is qualitatively different: rest-mass is an additive shift in $\omega^2$ of a free spectrum; Landau is a complete *re-quantization* of the kinetic operator with a new ladder structure ($n \in \mathbb{Z}_{\geq 0}$, eigenvalue spacing $\hbar\omega_c$).

The closer rest-mass analog might be the diamagnetic shift $\Delta E_\text{diam} = (e^2 B^2/8m_e c^2)\langle r_\perp^2\rangle$, which is quadratic in $B$ and acts as a perturbative additive shift. But this is the *non-resonant* contribution, valid only when $\omega_c \ll$ atomic energies. For strong $B$ (Landau regime), it fails completely.

Verdict: rest-mass does NOT absorb $B$ at the structural level. Stark linear shift in atoms is also categorically distinct from rest-mass: $-\vec{\mathcal{E}}\cdot\vec{r}$ has zero diagonal element on $\ell$-states (parity), so the leading Stark effect is *quadratic* and non-additive in the spectrum.

### 2.3 Observation / temporal-window (#15) — does it absorb $B$?

This is the most structurally interesting comparison. Both Matsubara and Landau quantize a continuous direction to give a discrete tower with a characteristic spacing:

| Projection | Compactified direction | Quantum number | Spacing |
|---|---|---|---|
| Observation/temporal-window (#15) | Time $t \to S^1_\beta$ | $k \in \mathbb{Z}$ | $2\pi/\beta$ |
| Cyclotron orbit (Zeeman?) | Cyclotron orbit angle $\phi \in [0, 2\pi)$ | $n \in \mathbb{Z}_{\geq 0}$ | $\hbar\omega_c$ |

Superficially parallel. But there are three structural differences that make them *not* the same mechanism:

1. **No spatial direction is compactified by a uniform $B$.** A uniform magnetic field on $\mathbb{R}^3$ leaves all three spatial directions non-compact. What it does is *confine* perpendicular motion to circular orbits of radius $\ell_B$ — but the underlying configuration space is unchanged. The Matsubara mechanism, by contrast, *literally identifies* $t \sim t + \beta$.

2. **The π appearance is different.** Matsubara modes are $\omega^t_k = 2\pi k/\beta$ — π is in the *eigenvalue spacing*. Landau levels are $E_n = \hbar\omega_c(n + \tfrac{1}{2})$ — **no π in the spacing**. The integration measure $d^2 r_\perp$ in the Landau problem produces a $2\pi$ factor in the *density of states* (degeneracy per Landau level $= eB/(2\pi\hbar c)$ per unit area), but this is a phase-space π, not a spectral π. The first π-bearing Landau eigenvalue does not exist as Matsubara's $(n=0, k=1)$ does.

3. **Bosonic vs fermionic boundary condition has no Landau analog.** Periodic vs antiperiodic time encodes statistics in #15. Landau levels are the same for any spinless particle; spin enters only through the additional Zeeman coupling $g_s\mu_B \vec{B}\cdot\vec{S}$, which is structurally an *added* term, not a modification of the orbital quantization.

Verdict: observation/temporal-window does NOT absorb $B$. The compactification analogy fails at the level of the spatial direction, the π appearance, and the statistical content. **The Landau spectrum looks like a harmonic oscillator (Bargmann-Segal sibling), not like a Matsubara tower.**

Stark is even further from #15: $-\vec{\mathcal{E}}\cdot\vec{r}$ does not compactify any direction, does not produce a discrete tower, and does not introduce π.

## 3. Compactification-template comparison

Paper 35's structural prediction is: **π enters a GeoVac observable iff the evaluation includes a continuous integration over a temporal/spectral parameter promoted from the discrete graph spectrum.**

Does Zeeman fit this template?

The Landau spectrum $E_n = \hbar\omega_c(n+\tfrac{1}{2})$ is structurally *identical* to a 2D harmonic oscillator at frequency $\omega_c$ — i.e., the Bargmann-Segal lattice (Paper 24) at radial dimension $d = 2$. Paper 24 proves the BS lattice is **bit-exactly π-free**. Therefore *the bare Landau spectrum is π-free*, just as the bare KG spectrum on $S^3 \times \mathbb{R}$ is.

Where does π enter Zeeman observables? In the *integration measure*: the Landau-level degeneracy per unit area is $eB/(2\pi\hbar c)$, the magnetic flux quantum is $\Phi_0 = 2\pi\hbar c/e$, etc. These are phase-space measure constants, structurally akin to the $1/(4\pi)$ Hopf-base measure of vector-photon promotion (#11) or the $\sqrt{\pi}$ Seeley-DeWitt coefficient of spectral action (#6). They are calibration-π, not Matsubara-2π.

**Verdict on compactification template:** Zeeman is **categorically different from Matsubara**. The bare Landau spectrum is HO-class and π-free; the π in Zeeman observables enters through the 2D phase-space measure (a calibration $\pi$ that lives in the master Mellin engine M1 mechanism per CLAUDE.md §1.7 WH1 / Sprint TS). Stark is even further removed: no discrete tower, no compactification, no inherent π.

## 4. Structural sketch if a new §III projection is required

Even though Zeeman and Stark are not absorbed by existing projections, the question is whether they *structurally* warrant a new §III entry, or whether they are better treated as Layer-2 inputs (external classical-field calibration, like $r_Z$ or Penin-Pivovarov $\alpha^4$).

**Candidate §III.17a: External uniform magnetic field (Zeeman/Landau)**
- Source: bare Fock-projected (or Bargmann-Segal) graph.
- Target: graph with covariant derivative $D \to D - i e\vec{A}/(\hbar c)$ for $\vec{A}$ symmetric gauge of uniform $\vec{B}$; spectrum acquires a Landau ladder for free electrons or a Zeeman shift $g_s\mu_B B m_s$ for bound electrons.
- Variable: $\omega_c = eB/(m_e c)$ (or equivalently $\ell_B$).
- Dimension: $[E]$ (with conjugate $[L]$ pinned).
- Transcendental signature: HO-class, π-free at the spectrum level (Bargmann-Segal sibling); calibration-π enters through 2D phase-space measure $1/(2\pi)$ per Landau level when integrating density of states. Lives in M1 ring of master Mellin engine for measure factors.
- Role per §IV.5: scale-extension (new $[E]$ instance, joining Fock / spectral action / Drake-Swainson / Breit retardation as the fifth $[E]$ injection).

**Candidate §III.17b: External uniform electric field (Stark)**
- Source: bare Fock-projected graph.
- Target: graph with diagonal perturbation $-\vec{\mathcal{E}}\cdot\vec{r}$; spectrum acquires linear Stark shift for degenerate states (hydrogen $n \geq 2$) or quadratic Stark shift $\propto \alpha_{\rm pol} |\vec{\mathcal{E}}|^2$ for non-degenerate states.
- Variable: $\mathcal{E} = |\vec{\mathcal{E}}|$.
- Dimension: $[E]/[L]$, equivalently $[E]$ when contracted with $a_0$.
- Transcendental signature: ring-preserving for rational $\mathcal{E}$ (in atomic units); polarizability $\alpha_{\rm pol}$ for hydrogen $1s$ is $9a_0^3/2 \cdot 4\pi\epsilon_0$, a pure rational in atomic units — no π in the spectrum.
- Role per §IV.5: scale-extension (new $[E]$ instance, parallel to but distinct from Zeeman).

**A more conservative alternative: bundle them as one "external classical background field" projection** with sub-sectors for $B$ (Zeeman-Landau) and $\mathcal{E}$ (Stark). This parallels how spectral action (#6) bundles UV-cutoff with internal gauge $\alpha$, or how Wilson plaquette (#10) bundles $\beta_{\rm Wilson}$ with the underlying SU(2) Haar measure. **Either decomposition is defensible**; the conservative version emphasizes the shared structural role (background-field coupling), the split version emphasizes the dimensional asymmetry ($B$ has conjugate $[L]$, $\mathcal{E}$ does not).

## 5. Empirical anchors that would test the framework

The §III projection question is structural; the empirical question is whether GeoVac at finite $n_{\max}$ produces correct numerical values for precision Zeeman/Stark observables.

- **Hanneke et al. (2008, 2011) electron $g$-factor in Penning trap.** $g/2 - 1 = a_e = 1.159\,652\,180\,73(28) \times 10^{-3}$. The framework computes $a_e$ on Dirac-S³ via the Schwinger asymptote (HF-2, sprint memo). Penning-trap $a_e$ is measured *in the presence of a uniform $B$* — the cyclotron and anomaly frequencies are extracted by fitting Landau spectra. For framework-internal consistency, the way GeoVac handles uniform $B$ would feed into the systematic of $a_e$ extraction, but at the precision of Hanneke (±$28 \times 10^{-13}$) any Zeeman framework error is absorbed in the Layer-2 multi-loop QED budget, not in Zeeman geometry.

- **Hydrogen Zeeman tensor structure (Lamb-Retherford-style).** In low $B$, hydrogen $1s$ has a linear Zeeman splitting $\Delta E = g_s\mu_B B m_s$. In strong $B$ (Paschen-Back), the LS-coupled multiplet decouples and individual $(m_\ell, m_s)$ states become eigenstates. The transition between regimes is the Paschen-Back point. The framework could in principle compute the full Zeeman tensor on the Dirac-S³ basis; this would be a clean test of whether external-$B$ adds any structure beyond the existing spinor-lift transcendental signature. Predicted to remain in $\alpha^2 \cdot \mathbb{Q}$ ring (Paper 18 spinor-intrinsic tier), with $\omega_c$ as Layer-2 input.

- **Stark shifts in Rydberg atoms** (e.g., Gallagher et al., precision measurements in cesium Rydberg states). The framework's hydrogenic-with-Z_eff Rydberg states should reproduce the linear Stark shift for hydrogen $n=2$ ($\Delta E = 3 e a_0 \mathcal{E}$ exact rational coefficient) and the quadratic Stark shift $\Delta E = -\tfrac{1}{2}\alpha_{\rm pol} \mathcal{E}^2$ where $\alpha_{\rm pol}^{\rm 1s} = 9/2$ exact rational in atomic units. *These are Layer-1 graph computations; no new transcendental enters.* That is the strongest argument that **Stark is Layer-2 input only**, not a new projection — every Stark coefficient at every order in $\mathcal{E}$ is in the bare-graph algebraic ring, with $\mathcal{E}$ as a parameter.

- **Penning-trap g-2 differential between B-field configurations.** If GeoVac's Zeeman handling is structurally complete, then the difference of $a_e$ measurements at different cyclotron frequencies should agree with the framework to within multi-loop budget (~LS-8a wall). This is a *consistency* test, not a new prediction.

## 6. Honest verdict

**Stark is Layer-2 input only — not a new §III projection.**
- Stark Hamiltonian $-\vec{\mathcal{E}}\cdot\vec{r}$ is a diagonal perturbation in the Fock-projected basis with rational matrix elements (computable via the same machinery as $V_\text{ne}$, indeed it is a uniform-external version of $V_\text{ne}$ at infinite source distance). All Stark shifts at all orders are rational in atomic units. Adding a new projection for "uniform external classical background that couples linearly to position" would be cosmetic — the mechanism is just a Hamiltonian perturbation with a numerical parameter, indistinguishable structurally from the Layer-2-input role of $r_Z$ or external g-factors.
- **Recommendation: treat $\mathcal{E}$ as a Layer-2 calibration parameter** in any Stark precision row of the Paper 34 catalogue. No new §III entry.

**Zeeman is borderline; the cleanest reading is "sub-feature of the Bargmann-Segal projection (#3) when $B$ is strong, sub-feature of the spinor lift (#7) when $B$ is weak."**
- In the **strong-$B$ Landau regime**, the cyclotron motion produces an HO-class spectrum; this is structurally a Bargmann-Segal lattice (Paper 24) at radial dimension 2 (perpendicular plane only), with the Landau ladder as the SU(2) (N, 0) tower at $d=2$. **Zeeman in the Landau regime is a Bargmann-Segal projection on the cyclotron plane, with $\omega_c$ playing the role of $\hbar\omega$.**
- In the **weak-$B$ Zeeman regime** (atomic-physics limit), the linear Zeeman shift $g_s\mu_B B m_s + g_\ell\mu_B B m_\ell$ acts on the existing spinor-lift basis with no new geometry. This is a sub-feature of the Camporesi-Higuchi spinor lift (#7), parameterized by the external dimensionless ratio $\omega_c/\text{Ry} \sim \alpha^2 B/B_0$ where $B_0 = m_e^2 c^3/(e\hbar) \approx 4.7 \times 10^9$ T is the atomic-units $B$ scale.
- The two regimes are connected by the Paschen-Back transition; in both, the framework has the right machinery without a new projection.
- **Recommendation: treat $\omega_c$ (or $B$) as a Layer-2 calibration parameter** that engages either the Bargmann-Segal projection (Landau regime) or the spinor-lift projection (Zeeman regime). The transition between regimes is a sub-feature of the existing dictionary, not a new projection.

**No new §III entries are warranted by external uniform fields.** Both $B$ and $\mathcal{E}$ are external classical parameters (Layer-2) that engage existing projections (Bargmann-Segal, spinor lift, Fock) when their value enters a precision observable. The base-unit accounting in §IV.5 is preserved without expansion.

**Finer point worth recording:** the audit's question "is $[E]$ injected by a fifth projection?" is answered NO at the projection level. $[E]$ remains four-instance (Fock, spectral action, Drake-Swainson, Breit retardation). External $B$ and $\mathcal{E}$ are Layer-2 *parameters* that engage the foundational $[E]$ scale (Ry) by external comparison ($\omega_c/$Ry, $\mathcal{E} a_0/$Ry); they don't *inject* a new $[E]$ scale in the §III sense. This is consistent with §IV.5's "scale-extension" classification being reserved for projections that *introduce a new variable* with its own dimension; external classical fields are external numerical parameters, not new variables in the sense of $\Lambda$ or $K$.

**The structural distinction that matters:** internal (graph-derived) variables vs external (classical-background) parameters. §III projections inject internal variables. External classical fields are Layer-2. The audit's hunting question pinned this distinction, which had been implicit.

**Open follow-up if the framework ever runs a precision Penning-trap or atomic-Stark catalogue row:** verify empirically that $a_e^{\rm Penning}$ extraction works with $\omega_c$ as Layer-2 parameter on the Bargmann-Segal-as-Landau interpretation, and that hydrogen Stark $\alpha_{\rm pol}^{1s} = 9/2$ a.u. exactly reproduces from graph-native computation. Both are predicted to land cleanly without invoking any new projection. If either fails — for example, if the Landau-regime-Bargmann-Segal correspondence requires a structural modification not currently in §III.3 — that would be the empirical signal to revisit the §III.17 question.
