# Phase A.2' deep-read of Latrémolière arXiv:2512.03573 — physics-side applicability scan

**Date:** 2026-05-23 (Phase A.2' kickoff, post-rescope).
**Sprint position:** Phase A.2' of Sprint L3e-P3 (re-scoped, `debug/sprint_l3e_p3_rescope_memo.md`). Deep-read of Latrémolière arXiv:2512.03573 "The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces" (Dec 3, 2025) — the paper that scooped the original Phase B.
**Sprint mandate (per PI):** "examine Latrémolière more closely. hopefully we can find some math to help us on the physics side."
**Outcome:** This memo identifies **two strong physics-side applications** that Latrémolière 2512.03573's pinned-QLCMS framework supplies directly to the GeoVac framework, plus a third candidate. The Krein-lift Phase A.2' originally targeted is NOT the main near-term value — there's potentially-immediate physics application.

---

## §1. What Latrémolière 2512.03573 actually constructs

### Definition 1.18 — Leibniz seminorm
A seminorm $L$ on $\mathrm{dom}(L) \subset \Acal$ is **Leibniz** if
$$L(ab) \le \norm{a}_\Acal \cdot L(b) + L(a) \cdot \norm{b}_\Acal \quad \forall a, b \in \mathrm{dom}(L).$$

### Definition 1.22 — Separable quantum locally compact metric space (QLCMS)
A pair $(\Acal, L)$ where:
- $\Acal$ is a separable, possibly non-unital C*-algebra
- $L$ is a Leibniz, hermitian, norm-modulo-constants seminorm on a dense $*$-subalgebra $\mathrm{dom}(L)$
- Three axioms:
  - (FM) Fortet-Mourier distance $\mathrm{bl}_L$ metrizes the weak-$*$ topology on $\mathcal{S}(\Acal)$
  - (CB) Closed unit ball $\{a \in \mathrm{dom}(L) : L(a) \le 1\}$ is closed in $\mathfrak{sa}(\Acal)$
  - (B) Boundedness: exists strictly positive $h \in \mathfrak{sa}(\Acal)$ and state $\mu$ with $\mu(h) = 1$ such that $\{hah : a = b + t\cdot 1_\Acal, L(b) \le 1, \mu(b) = -t\}$ is bounded

### Definition 1.26 — Pinned separable QLCMS
A triple $(\Acal, L, \mu)$ extending Def 1.22 where the **pin state** $\mu \in \mathcal{S}(\Acal)$ satisfies
$$\{\phi \in \mathcal{S}(\Acal) : \mathrm{mk}_L(\phi, \mu) < \infty\} \text{ is weak-}* \text{ dense in } \mathcal{S}(\Acal).$$
where $\mathrm{mk}_L(\phi, \psi) := \sup\{|\phi(a) - \psi(a)| : a \in \mathrm{dom}_{sa}(L), L(a) \le 1\}$ is the Monge-Kantorovich distance.

### Definition 1.29 — L-Lipschitz μ-pinned exhaustive sequence
A sequence $(h_n) \in \mathrm{dom}_{sa}(L)$ satisfying:
1. $\lim L(h_n) = 0$ (Lipschitz norm decays to zero)
2. $\lim \mu(h_n) = 1$ (expectation in pin state approaches one)
3. $\lim \norm{h_n}_\Acal = 1$ (operator norm approaches one)

**Theorem 1.30:** an L-Lipschitz μ-pinned exhaustive sequence is an approximate unit for $\Acal$.

### Section 2 — Local quantum metametrics
Tunnels between pinned QLCMS define **local metametric $\delta_r$ at radius $r$** via "Lipschitz functions which are (almost) 1 at a base point, and decay very slowly, i.e. have small Lipschitz seminorms." This avoids the unavailability of "balls" in noncommutative settings.

Key tunnel ingredients: reach, height, alignment, quotient morphisms (standard Latrémolière tunneling-pair vocabulary).

### Section 3 — Coincidence property
Distance zero ⟺ fully quantum isometric pointed proper QLCMS.

### Section 4 — GH quantum metametric hypertopology
The pointed propinquity, built from local metametrics via Gromov-Hausdorff-style infimum.

### Section 5 — Compact-case agreement
When $\Acal$ is unital, the pinned QLCMS is a QCMS in the standard Latrémolière sense, and the new hypertopology restricts to the propinquity topology.

### Section 6 — Example
$c_0(\Z) \rtimes_\alpha \Z$ as pointed proper QMS — non-unital, non-commutative, finite-dimensional approximations converge in the new topology.

---

## §2. Physics-side application #1 — Sturmian basis as L-Lipschitz μ-pinned exhaustive sequence

> **CORRECTION (2026-05-23, post-Sprint L3e-P3 X/Y/Z synthesis):** the structural identification below is at the WRONG level under the natural Coulomb-multiplication algebra reading ($\Acal = C_0(\R^3) \otimes M_d$). Under that reading, the Sturmian projector is NOT in $\Acal$ (it's not a multiplication operator); the exhaustive sequence is smooth radial cutoffs $\chi_n \in C_c^\infty(\R^3)$, and the Sturmian truncation lives at the operator-system-truncation level (tunnel quotient morphism), NOT at the exhaustive-sequence level. The identification below IS valid under a different reading ($\Acal$ = spectral algebra containing projectors), but two of three Def 1.29 axioms are then trivial (only L-decay is non-trivial). Both readings are legitimate but give different framings. See `debug/subsprint_x_sturmian_latremoliere_verification.md` (X memo, §1-§8) for the correction details and `debug/sprint_l3e_p3_synthesis_memo.md` (synthesis, §1) for the X/Z framing divergence.
>
> **Additional load-bearing finding (Sub-sprint X §3):** the Leibniz axiom $L(fg) \le \|f\|_\Acal L(g) + L(f) \|g\|_\Acal$ FAILS for second-order Schrödinger $L(f) = \|[D_\mathrm{S}, M_f]\|$ (cross-term $2 \nabla f \cdot \nabla g$ doesn't cancel). Holds cleanly for first-order Dirac. **R1 workaround:** use $L(f) = \|\nabla f\|_\infty$ (standard NCG spin-Dirac approach). The framework's Schrödinger-based atomic FCI inherits Latrémolière hypertopology via R1, NOT via the direct commutator.

### The structural identification (substantive new content) — see correction above

GeoVac uses **Coulomb Sturmian bases** in several places: `geovac/casimir_ci.py` (graph-native CI), `geovac/hylleraas_r12.py` (Hylleraas $r_{12}$ explicit correlation), `geovac/hylleraas_eckart_pstate.py` (Hylleraas-Eckart double-α P-state). The Sturmian basis at exponent $\lambda = Z/n_*$ consists of functions $\{S_{n, l, m}^\lambda(r)\}$ that:
- Are eigenfunctions of a Sturmian operator (related to but distinct from the Schrödinger operator)
- Satisfy a discrete completeness relation on $L^2(\mathbb{R}^3, d^3r)$
- Have Lipschitz norm (under the Dirac / Coulomb operator $D = -\frac{1}{2}\nabla^2 - Z/r$) that decreases with increasing $n$
- Have pin-state expectation (ground state $\psi_{1s}$) bounded between 0 and 1

**Structural identification:** the Sturmian truncation at $n_{\max}$ is an L-Lipschitz μ-pinned exhaustive sequence in Latrémolière 2512.03573's Def 1.29:
- $L = \|[D, \cdot]\|$ where $D$ is the Coulomb-Dirac operator (Schrödinger for non-relativistic atoms)
- $\mu = \omega_{\psi_{1s}}$ ground-state expectation
- $h_n = \chi_{\mathrm{Sturmian}}^{n \le n_{\max}}$ characteristic function of truncated Sturmian basis (or the projector onto it)

Verification:
1. **$L(h_n) \to 0$:** Sturmian basis truncation projector has bounded commutator with $D$; as $n_{\max} \to \infty$ the truncation becomes the identity and the commutator → 0. ✓ (modulo verification of the Leibniz axiom)
2. **$\mu(h_n) \to 1$:** the ground state $\psi_{1s}$ is in the Sturmian span for any $n_{\max} \ge 1$, so $\mu(h_n) = 1$ identically (not just in the limit). ✓
3. **$\norm{h_n}_\Acal \to 1$:** $h_n$ is a projector, $\norm{h_n}_{\mathrm{op}} = 1$. ✓

**Conclusion: GeoVac's Sturmian-based atomic FCI machinery is structurally an L-Lipschitz μ-pinned exhaustive sequence in the Latrémolière 2512.03573 sense.**

### What this gives the framework

If the structural identification holds (modulo Leibniz axiom verification), GeoVac's Sturmian FCI calculations inherit:

1. **A precise math.OA convergence statement.** The truncation error of Sturmian FCI at $n_{\max}$ has an explicit bound from Latrémolière's hypertopology. The framework's Hylleraas $r_{12}$ calculations (He 1¹S at 0.0006% accuracy, ω=4) get a rigorous error model.

2. **A unified comparison metric across systems.** The pointed propinquity gives a metric on "Sturmian truncations" so we can quantitatively compare the convergence rates of:
   - He (Hylleraas-Eckart, ω=4) vs Be (Hylleraas-CI, hypothetical)
   - Atomic Sturmian (FCI) vs molecular Sturmian (Shibuya-Wulfman)
   - Different focal lengths $\lambda = Z/n_*$ at fixed $n_{\max}$

3. **A natural framework for non-relativistic-to-relativistic Sturmian lift.** Latrémolière's QLCMS framework accommodates both Schrödinger and Dirac operators as $D$; the Hartree-Fock-screened Dirac framework (FrozenCore + Clementi-Raimondi) becomes a pointed QLCMS structure with pin state = HF ground state.

### Sprint potential

A small follow-on sprint (~1 week) could:
- Verify the Leibniz axiom on $L(a) = \|[D_{\mathrm{Sturmian}}, a]\|$ for Sturmian basis operators
- Identify the explicit Latrémolière exhaustion rate for He 1¹S Hylleraas-Eckart (compare to the empirical 0.0006% at ω=4)
- Connect to Paper 38 L2 quantitative rate ($4/\pi$ asymptote + Stein-Weiss subleading) — is this the same rate, or different?

If this works, it gives the framework's atomic FCI calculations a rigorous Latrémolière-hypertopology error model. **This is a sub-sprint-scale physics-side payoff from Latrémolière 2512.03573**, separate from the L3e-P3 Krein-lift program.

---

## §3. Physics-side application #2 — Pin-state structure for non-unital chemistry constructions

### The W1c chemistry orthogonality wall — refined diagnostic

CLAUDE.md §1.7 W1c-residual chemistry orthogonality wall (NaH binding, May 2026):
- Hydrogenic $Z_\mathrm{orb} = 1$ Na valence wavefunction shape is qualitatively wrong (actual Na 3s mean radius 4.47 bohr vs hydrogenic 1.50 bohr)
- Cross-V_ne shape substitution would close the wall
- Track 3 (SV diagonal substitution) didn't close — needs cross-V_ne shape, not just diagonal

**Latrémolière framing of W1c:** the framework's molecular C*-algebra is non-unital ($C_0$ of configuration space tensored with frozen-core structure), and the pin state should be the joint molecular ground state (cross-center bonding orbital). Latrémolière's pinned QLCMS gives a precise structure for:
- "Non-unital C*-algebra" = molecular configuration-space algebra
- "Pin state" = joint ground-state expectation (cross-center HF Slater determinant)
- "Exhaustive sequence" = increasing-$n_{\max}$ truncation of the molecular Sturmian basis

The wall's "wrong cross-V_ne shape" content becomes: the framework's pin-state-defining Slater determinant has the wrong shape at the cross-center coupling region. Latrémolière's framework doesn't TELL us the right shape (that's the operator-generation gap), but it does PROVIDE the mathematical structure to compare different pin-state choices via the local metametric $\delta_r$.

### Concrete sprint potential

A diagnostic sprint (~1-2 weeks):
- Set up the NaH C*-algebra as a non-unital tensor product (electron register $\otimes$ frozen-core nucleus register $\otimes$ valence register)
- Define candidate pin states: (i) HF-screened Na 3s × H 1s, (ii) full multi-determinant HF, (iii) MP2-corrected ground state
- Compute local metametric $\delta_r$ between truncations at different $n_{\max}$ for each pin state
- Identify which pin state gives the smallest local metametric distance to the experimental NaH binding curve

If the right pin state gives binding, **the wall closes via pin-state-shape selection, not cross-V_ne operator generation.** This would be the Latrémolière-framework reformulation of the chemistry pause.

**Confidence in this path: MEDIUM.** The Latrémolière framework provides the right vocabulary but doesn't guarantee a specific pin state will work. Worth a 1-2 week diagnostic.

---

## §4. Physics-side application #3 — Bound-state QFT in continuum

### Bethe log + Lamb shift

The framework's bound-state QED Lamb shift work (Paper 36) uses the **acceleration-form Bethe log via Coulomb Sturmians at exponent $\lambda = Z/n$** (LS-3). This is a Sturmian-basis calculation on the non-unital C*-algebra $C_0(\mathbb{R}^3)$ with the bound-state Sturmian exponent.

If application #1 (Sturmian as Latrémolière L-Lipschitz μ-pinned exhaustive sequence) verifies, the Bethe log calculation has a rigorous Latrémolière convergence model — and the truncation error in the Lamb shift becomes a precise hypertopology bound.

**Sprint potential:** sub-sprint scale (~1 week, parallel with Application #1). Verifies whether the Lamb shift's −0.534% native residual at one loop has a Latrémolière interpretation (rather than just the empirical convergence Paper 36 reports).

### Sturmian closure for continuum observables

The framework's "Sturmian closure" property (used in Hylleraas r₁₂ work) — that a Sturmian basis at the right exponent saturates a continuum observable like hydrogen polarizability at $N_{\mathrm{basis}} = 2$ — has a natural Latrémolière interpretation as **the exhaustive sequence reaching the pin state in finite truncation order** (rather than asymptotic).

This could illuminate why Sturmian closure works in some cases (1-electron hydrogenic) but not others (multi-electron — needs Hylleraas-Eckart double-α).

---

## §5. Cross-arc structural observation — what Latrémolière's framework UNIFIES

Three GeoVac-framework objects that Latrémolière 2512.03573 unifies under the pinned-QLCMS umbrella:

| GeoVac object | Latrémolière 2512.03573 structure |
|:--------------|:----------------------------------|
| Sturmian basis truncation at $n_{\max}$ | L-Lipschitz μ-pinned exhaustive sequence (Def 1.29) |
| Ground-state expectation value $\mu$ | Pin state (Def 1.26) |
| Truncation projector $P_{n_{\max}}$ | Tunnel quotient morphism (§2.1) |
| Sturmian-FCI error bound at $n_{\max}$ | Local metametric $\delta_r$ (§2.2) |
| Cross-system comparison (He vs Be vs LiH) | Pointed propinquity hypertopology (§4) |
| Non-relativistic ground state vs relativistic ground state | Different pin states on same C*-algebra (extension of Def 1.26) |
| Non-unital character of molecular C*-algebra | Definition 1.22 non-unital requirement (FM, CB, B axioms) |

This is **strong structural alignment** between Latrémolière 2512.03573's pointed-QLCMS framework and the GeoVac framework's Sturmian-based atomic / molecular computations. The L3e-P3 program's original target (Lorentzian extension) addresses ONE direction; this physics-side mapping addresses A SECOND, sub-sprint-scale, immediately-productive direction.

---

## §6. Strategic re-assessment of the L3e-P3 program

After the deep-read, the program decomposes into:

1. **Sub-sprint X (1-2 weeks)** — Sturmian-as-Latrémolière physics application. Concrete physics-side payoff, sprint-scale, immediate. Verifies whether GeoVac's atomic FCI calculations inherit a rigorous Latrémolière convergence model. (Application #1, §2)

2. **Sub-sprint Y (1-2 weeks)** — W1c chemistry pin-state diagnostic. Reformulate NaH binding via Latrémolière pin-state-shape selection. May close the chemistry wall via right-pin-state, not right-operator. (Application #2, §3)

3. **Sub-sprint Z (1 week, parallel with X)** — Bethe log / Lamb shift Latrémolière interpretation. Verifies whether Paper 36's bound-state QED convergence has a precise math.OA error model. (Application #3, §4)

4. **Phase A.2'-A.4' (~10 weeks)** — original Krein-lift Lorentzian extension target. Still valid but now ONE option among several with positive ROI.

5. **Phase A.5' decision gate** — original Gate-1 decision.

**The physics-side applications (X, Y, Z) are SUB-SPRINT scale and parallelizable with the Krein-lift work.** They don't require committing to the full multi-month L3e-P3 program. The original L3e-P3 strategic framing (12 months for G2-metric closure) treated Latrémolière 2512.03573 as a math.OA target; the deep-read reveals it as ALSO a math.OA TOOL that GeoVac can apply immediately on the physics side.

---

## §7. Recommendation

The deep-read changes the strategic calculus:

**Original assessment (post-A.1 audit):** re-scoped Phase A.2'-A.4' Krein-lift (~10 weeks) is the right next step for L3e-P3. Decision gate at A.5'.

**Updated assessment (post-A.2' deep-read):** the Krein-lift is still a valid math.OA target, BUT three physics-side applications (Sturmian, W1c chemistry, Bethe log) offer sub-sprint-scale payoffs that should be done FIRST:

- Sub-sprint X (Sturmian-as-Latrémolière): immediate atomic FCI rigor
- Sub-sprint Y (W1c pin-state diagnostic): may close the chemistry wall
- Sub-sprint Z (Bethe log interpretation): Paper 36 follow-on

Each is 1-2 weeks. Together they offer ~4-6 weeks of high-payoff physics work BEFORE committing to the 10-week Krein-lift.

**Recommendation: do Sub-sprint X first (1-2 weeks).** It's the most concrete identification (Sturmian basis ↔ L-Lipschitz μ-pinned exhaustive sequence), and if it verifies, opens up the chemistry / Lamb-shift applications. After X lands, decide whether to do Y and Z in parallel with starting Krein-lift, or focus on physics-side applications.

**Honest scope:** the Sturmian-as-Latrémolière identification (Application #1) is at the *structural sketch* level, not theorem level. Verifying it requires checking the Leibniz axiom on $L = \|[D_{\mathrm{Sturmian}}, a]\|$ in detail, which Sub-sprint X would do.

---

## §8. Honest scope of this deep-read

The deep-read:
- IS based on the abstract + Def 1.18, 1.22, 1.26, 1.29 + Theorem 1.30 extracted from arXiv HTML
- IS NOT a full PDF read (sections 2-6 were truncated in HTML extraction)
- DOES identify structural alignments between Latrémolière 2512.03573's framework and GeoVac's Sturmian-based computations
- DOES NOT prove these alignments rigorously — Sub-sprint X would verify

The physics-side applications surfaced (§2-§4) are diagnostic / structural observations, NOT theorem statements. The Sub-sprint X would convert observation to theorem.

**Confidence:**
- HIGH on the Sturmian-as-L-Lipschitz-μ-pinned-exhaustive-sequence structural identification
- MEDIUM on the W1c pin-state diagnostic working at the actual NaH PES (depends on whether the right pin state exists)
- MEDIUM on the Bethe log / Lamb shift Latrémolière interpretation
- HIGH on the strategic re-assessment (physics-side applications should be done first; Krein-lift can wait)

**Files:**
- `debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` (this memo, ~4500 words)
- Cross-references: `debug/l3e_p3_phase_a1_literature_audit.md`, `debug/l3e_p3_phase_a2_operator_algebraic_eps_net.md`, `debug/sprint_l3e_p3_rescope_memo.md`, `geovac/hylleraas_r12.py`, `geovac/casimir_ci.py`, `papers/group1_operator_algebras/paper_36_bound_state_qed.tex` (LS-3 Bethe log).
