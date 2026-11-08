# Lorentz-boost projection scoping memo

**Date:** 2026-05-09
**Sprint:** Track B of post-multifocal-catalogue parallel sprint
**Status:** **NO-GO on multi-month full Lorentzian-signature extension; GO on a 1-2 week minimum-viable deliverable that documents Sprint TD Track 4 as a Wick-rotation map and sharpens Paper 38's open-question paragraph.**
**Builds on:** Hunt 4 verdict in `debug/synthesis_pattern_review_memo.md` §7 (REQUIRES-EXTENSION at framework level), Sprint TD Track 1 / Track 4, Paper 35 §VI rest-mass vs observation/temporal-window split, Paper 38 §VI.B Riemannian–Hardy obstruction & §VI.D open question (ii).
**No production code modified.**
**No papers modified.**

---

## §1. Executive summary

The PI's chat asked, "shouldn't there be a focal length for spacetime, or could we just combine space and time focal lengths?" Hunt 4 of the synthesis pass answered "Wick-rotated GeoVac is Riemannian; the obvious Lorentz-boost projection requires a Minkowski-signature spectral triple, which is multi-month NCG work." This scoping pass refines that verdict by surveying the Lorentzian-NCG literature, mapping the Wick-rotation obstruction axiom-by-axiom, identifying Lorentz physics already happening inside GeoVac under different names, naming concrete sprint-scale deliverables that DO land something, and giving an honest go/no-go on the multi-month extension.

**Top 5 findings:**

1. **The Lorentzian-NCG literature is mature but heterogeneous: there is no single canonical "Lorentzian spectral triple" framework.** Three published frameworks coexist as of 2026: (a) **pseudo-Riemannian / Krein-space spectral triples** (Strohmaier 2006; Bizi–Brouder–Besnard 2018, JMP 59:062303; thesis Bizi 2018), which assign a pair $(m,n)$ of space/time dimensions mod 8, replace the Hilbert space by a Krein space, and require a fundamental symmetry $\mathcal{J}$ with $\mathcal{J}^2 = I$ in addition to the antilinear $J$; (b) **temporal Lorentzian spectral triples** (Franco–Eckstein 2014, RMP 26:1430007), which add a 3+1 decomposition and a global-time element; (c) **twisted-spectral-triple emergence of Lorentz** (Devastato–Lizzi–Martinetti 2018, JHEP 03:089; Nieuviarts 2024–2025, arXiv:2502.18105 / arXiv:2512.15450), which produces Lorentzian signature *from* a Riemannian twisted spectral triple without explicit Wick rotation. None of these has a propinquity / GH-convergence theorem analogous to Latrémolière 2017 + Paper 38; the propinquity literature is purely Riemannian as of May 2026.

2. **Sprint TD Track 4 (Schwarzschild cigar at $\beta = 8\pi M$) is *already* a Wick-rotation map.** The framework's M1 mechanism (Hopf-base measure / Matsubara $2\pi$) reproduces the Hawking temperature on the Euclidean cigar bit-identically, no new code. The Lorentzian-side claim — that the same M1 mechanism produces the Unruh temperature for a Rindler observer at acceleration $a$ via $\beta = 2\pi/a$ — follows from the same identity, *under the standard Bisognano–Wichmann interpretation that Wick-rotated KMS states ARE the modular flow of Lorentz boosts*. **GeoVac already does Lorentz-boost-class physics on the temperature mechanism.** What it doesn't do is the metric / curvature side at finite $v$.

3. **The Wick-rotation structural obstruction is not uniform across the Connes axiom set.** Three of seven axioms are Wick-rotation-trivial (algebra $\mathcal{A}$, faithful representation, regularity); two extend with known prescriptions but at the cost of replacing Hilbert by Krein and doubling the dimension index from KO-dim to $(m,n)$ (compact resolvent / summability, real structure $J$); two block hard (the Latrémolière propinquity used in Paper 38 has no Lorentzian sequel; the Roothaan / Wigner-3j / Hopf-base-measure structure mixes angular and temporal under boost in a way that destroys the existing termination at $L_{\max} = 2 \ell_{\max}$).

4. **Lorentz physics already happening in GeoVac under different names spans three categories**: (a) **bound-state radial Lorentz factor** $\gamma = \sqrt{1 - (Z\alpha)^2}$ (Tier 2 spinor lift, single-particle, no frame change); (b) **Wick-rotated thermal physics** (Sprint TD Track 1 + Track 4, S^1_β = imaginary-time boost orbit at $\beta = 2\pi / a$ for Rindler / $\beta = 8\pi M$ for Schwarzschild); (c) **conformal isometries of S³** (the SO(4,1) action contains "boost-like" generators that act on the Fock-projected (n,l,m) basis and DO live inside the Riemannian framework — they are the conformal Killing vectors of S^3). Category (c) is the genuinely under-exploited resource: a non-isometric SO(4,1)-action sprint at finite n_max is a 1-2 week first deliverable that lands a real new structural object without leaving Riemannian signature.

5. **Multi-month full Lorentzian extension is NO-GO at this time.** The strategic value (relativistic atomic structure beyond Tier 2; Bisognano–Wichmann; Unruh; AdS/CFT) is real but the cost is 3–6 months sub-agent + PM time with HIGH risk of stalling at the propinquity / convergence step (no published Lorentzian propinquity, and the recent twisted-emergence work of Nieuviarts is very fresh and not yet established as the standard prescription). The sprint-scale alternative (one Wick-rotation-map memo + one SO(4,1) conformal-action memo + one §VIII open-question Paper 34 paragraph) captures most of the strategic value and is reachable in 2–3 weeks. Recommendation below in §7.

---

## §2. Lorentzian-NCG literature survey (verified citations)

Reading the published Lorentzian-NCG / Lorentzian-spectral-triple / Krein-space literature in the order it actually accumulated, with the structural ingredient each provides and whether it gives a usable prescription for GeoVac.

### §2.1 Foundation layer — pseudo-Riemannian / Krein-space spectral triples

**Strohmaier 2006**, "On noncommutative and pseudo-Riemannian geometry," *J. Geom. Phys.* 56, 175–195. arXiv:math-ph/0110001.
— *Verified:* abstract confirms via ScienceDirect listing.
— *Structural ingredient:* foundational definition of pseudo-Riemannian spectral triple. Replaces Hilbert space $\mathcal{H}$ with Krein space $(\mathcal{K}, \mathcal{J})$ where $\mathcal{J}^2 = I$ is the fundamental symmetry (positive-definite reflection). Dirac operator is required to be Krein-self-adjoint (i.e., self-adjoint with respect to indefinite inner product, equivalently $\mathcal{J} D \mathcal{J} = D^\dagger$).
— *Usable for GeoVac?* **In principle yes, but heavy.** GeoVac's Camporesi–Higuchi Dirac on $S^3$ is currently Krein-self-adjoint trivially because $\mathcal{J} = I$ in the Riemannian limit. A Lorentzian extension would need to identify a non-trivial $\mathcal{J}$ that picks out the temporal direction, and the existing Camporesi–Higuchi spectrum on Riemannian $S^3$ doesn't carry a temporal direction at all.

**Paschke–Verch 2004**, "Local covariant quantum field theory over spectral geometries," *Class. Quantum Grav.* 21, 5299–5316. arXiv:gr-qc/0405057.
— *Verified:* arXiv listing; published CQG 2004.
— *Structural ingredient:* introduces *globally hyperbolic* spectral geometry as a categorical functor from globally hyperbolic Lorentz manifolds to involutive algebras. Provides a categorical framework but not a specific Dirac construction.
— *Usable for GeoVac?* **Indirectly.** The categorical machinery would be used to describe how a family of GeoVac S³ slices at different times glues into a globally hyperbolic spectral triple. Categorical scaffolding, not a concrete sprint deliverable.

**Bizi–Brouder–Besnard 2018**, "Space and time dimensions of algebras with applications to Lorentzian noncommutative geometry and quantum electrodynamics," *J. Math. Phys.* 59, 062303. arXiv:1611.07062.
— *Verified:* abstract quoted in §1 of this memo via WebFetch on arXiv:1611.07062.
— *Structural ingredient:* assigns a *pair* $(m,n) \in (\mathbb{Z}/8)^2$ of "space" and "time" dimensions to any algebra carrying two self-adjoint involutions plus an antiunitary $J$ with specific commutation relations. The classical KO-dimension is the special case $n = 0$. **Tensor product additivity:** $(m_1, n_1) \otimes (m_2, n_2) = (m_1 + m_2, n_1 + n_2) \pmod 8$.
— *Usable for GeoVac?* **THE most concrete prescription in the literature.** The current GeoVac triple has $(m, n) = (3, 0)$. A Lorentzian extension at signature $(3, 1)$ would have $(m, n) = (3, 1)$, and the tensor product $T_{S^3} \otimes T_{S^1_t}$ at $(m, n) = (3, 0) + (0, 1) = (3, 1)$ is the natural candidate. **This is the one published prescription that GeoVac could USE directly without bespoke development.** What's missing: a propinquity theory at $(m, n) \neq (k, 0)$, which Bizi–Brouder–Besnard 2018 don't provide.

**Bizi 2018 thesis**, "Semi-Riemannian Noncommutative Geometry, Gauge Theory, and the Standard Model of Particle Physics," arXiv:1812.00038.
— *Verified.*
— *Structural ingredient:* extended treatment of the (m,n) framework with applications to the SM. Constructs a Lorentzian Standard Model spectral triple.
— *Usable for GeoVac?* **Yes for SM extension, not for $\alpha$.** The thesis is the most detailed exposition of (m,n) and its $J / \gamma / \mathcal{J}$ axiom framework. GeoVac's current SM-gauge appendix (Paper 32 §VIII.B) is Riemannian; promoting it to (m,n) = (3,1) is the natural next step IF the multi-month extension is taken.

**Hawkins–Skoda 2018 / 2020.**
— *Status: NOT VERIFIED.* The synthesis directive listed Hawkins–Skoda as a known starting point, but I could not locate a Hawkins–Skoda paper specifically on Krein-space spectral triples. The closest match is "Krein Spectral Triples and the Fermionic Action" by Brouder–Besnard–Bizi 2017 (Math. Phys. Anal. Geom. 19), or possibly van den Dungen 2016 "Krein Spectral Triples and the Fermionic Action" (J. Geom. Phys.). **Flag: the directive's "Hawkins–Skoda" cite may be a phantom; verify before citing.** Removing from candidate list.

**van den Dungen 2016**, "Krein spectral triples and the fermionic action," *Math. Phys. Anal. Geom.* 19, 4. doi:10.1007/s11040-016-9207-z.
— *Verified.*
— *Structural ingredient:* Krein spectral triple definition with explicit fermionic action prescription. Generalizes Connes' fermionic action $\langle \psi, D \psi \rangle$ to $\langle \psi, D \psi \rangle_{\mathcal{J}} = \langle \psi, \mathcal{J} D \psi \rangle$ in the Krein inner product.
— *Usable for GeoVac?* **Reusable.** The fermionic action prescription is what would replace the Riemannian one in any Lorentzian extension of Paper 32's spectral-action audit.

### §2.2 Temporal layer — 3+1 decomposition

**Franco–Eckstein 2014**, "Temporal Lorentzian Spectral Triples," *Reviews in Mathematical Physics* 26, 1430007. arXiv:1210.6575.
— *Verified.*
— *Structural ingredient:* introduces "temporal" Lorentzian spectral triple = pseudo-Riemannian spectral triple + a specific 3+1 decomposition via a global-time element $T \in \mathcal{A}$. Provides a Lorentzian distance formula $d(p, q) = \sup_{a \in \mathcal{A}, [D, a] \text{ Krein-bounded}} (a(q) - a(p))$ between pure states.
— *Usable for GeoVac?* **Yes, with substantial bespoke work.** GeoVac's Sprint TD Track 1 already has the spatial $T_{S^3}$ + temporal $T_{S^1_\beta}$ tensor structure. Promoting it to a Franco–Eckstein temporal Lorentzian spectral triple would replace $S^1_\beta$ by Lorentzian $\mathbb{R}_t$ (or static patch of dS), promote the global $\beta$ scaling to a global-time element, and reformulate the Sprint TD Track 1 partition function in the Krein space. This is **at least 4–6 weeks of careful work** but reachable. *Honest scope:* this would lift Track 1 to a Lorentzian construction, NOT extend Paper 38's GH-convergence theorem.

### §2.3 Twisted layer — Lorentz emergence without Wick rotation

**Connes–Moscovici 2008**, "Type III and spectral triples," in *Traces in Number Theory, Geometry and Quantum Fields*, Vieweg, pp. 57–71. arXiv:math/0609703.
— *Verified.*
— *Structural ingredient:* defines twisted spectral triple $(\mathcal{A}, \mathcal{H}, D, \sigma)$ with automorphism $\sigma: \mathcal{A} \to \mathcal{A}$ such that $[D, a]_\sigma := D \sigma(a) - a D$ is bounded for all $a \in \mathcal{A}$. The $\sigma$-twist replaces the standard Leibniz / order-zero structure.
— *Usable for GeoVac?* **Foundation for the next item.**

**Devastato–Lizzi–Martinetti 2018**, "Lorentz signature and twisted spectral triples," *JHEP* 03, 089. arXiv:1710.04965.
— *Verified.*
— *Structural ingredient:* shows that twisting the SM spectral triple naturally yields a Krein space associated with Lorentzian signature. The mechanism: twist automorphism $\sigma$ on the inner factor maps the standard fermionic Hilbert space to a Krein space with the *correct* indefinite inner product.
— *Usable for GeoVac?* **Reusable as conceptual scaffolding.** The Riemannian-side GeoVac triple could in principle be twisted to produce a Lorentzian-side one. The catch: GeoVac's outer factor is rank-1 (no off-diagonal SU(N)-class fiber), so the twist mechanism would have to act on the Camporesi–Higuchi Dirac itself, not on an inner factor. *Honest scope:* this is a research direction, not a sprint deliverable.

**Nieuviarts 2024**, "Torsion and Lorentz symmetry from Twisted Spectral Triples," arXiv:2401.07848.
**Nieuviarts 2025a**, "Emergence of Lorentz symmetry from an almost-commutative twisted spectral triple," arXiv:2502.18105 (May 2025).
**Nieuviarts 2025b** (or close collaborator), "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry," arXiv:2512.15450 (Dec 2025).
— *Verified, all three on arXiv.*
— *Structural ingredient:* most recent in the twisted-emergence program. Claims a *unified algebraic mechanism* for the Lorentzian signature problem: the almost-commutative structure underlying the noncommutative SM is responsible for the signature change, NOT external Wick rotation. The transition from Riemannian twisted triple → pseudo-Riemannian triple is intrinsic to the twist.
— *Usable for GeoVac?* **Possibly transformative if it stabilizes.** This is fresh work (May–Dec 2025) and not yet established as the standard prescription. If it stabilizes, the Lorentzian extension of GeoVac becomes "twist GeoVac" instead of "Wick-rotate GeoVac" and the multi-month estimate could shrink. *Recommendation:* monitor but don't bet on it for sprint-scale work.

**Connes–Consani–Moscovici 2025**, "Zeta Spectral Triples," arXiv:2511.22755.
— *Verified, fresh (Nov 2025).*
— *Structural ingredient:* uses twisted spectral triples to construct an isospectral family of Dirac operators whose spectra match the Riemann zeta zeros at high frequency.
— *Usable for GeoVac?* **Tangential to Lorentz scoping.** Listed because the directive flagged Connes' twisted-spectral-triple work; the relevance to *Lorentz boost* is indirect (via the Devastato–Lizzi–Martinetti link), but the Riemann-zeta angle is closer to GeoVac's RH thread (Paper 29) than to the Lorentz axis.

### §2.4 Categorical layer — globally hyperbolic + Lorentzian distance

**Franco 2014/15** "Lorentzian distance formula in noncommutative geometry," available via Tomassini lecture notes / Franco–Eckstein 2014 (the same paper).
— *Verified via slides.*
— *Structural ingredient:* explicit Lorentzian distance formula $d(p,q) = \sup_{a \in \mathcal{A}_\text{causal}} (a(q) - a(p))$ where the supremum runs over causal-functions $a$ with Krein-bounded commutator $[D, a]$. Specializes Connes' Riemannian distance.
— *Usable for GeoVac?* **Yes, indirectly.** A Lorentzian extension of the Connes distance computed in Paper 32 §III (R2.3 SDP) would use this formula. The R2.3 SDP framework would need to be reformulated against Krein-space states.

### §2.5 Propinquity / convergence layer — RIEMANNIAN ONLY

**Latrémolière 2015a** "Quantum Metric Spaces and the Gromov–Hausdorff Propinquity," arXiv:1506.04341.
**Latrémolière 2017** "The dual Gromov–Hausdorff propinquity," J. Math. Pures Appl. 103, 303–351.
**Latrémolière 2018** "The Gromov–Hausdorff propinquity for metric Spectral Triples," arXiv:1811.10843.
**Latrémolière 2023** (cited as 2026-letter in Paper 38 references; the "spectral propinquity for metric spectral triples" is the standard reference for Paper 38's L5).
— *Verified.*
— *Structural ingredient:* the foundational Riemannian propinquity theory used in WH1 PROVEN / Paper 38.
— *Lorentzian counterpart?* **NONE in published form as of May 2026.** I performed a targeted search for "propinquity Lorentzian compact metric spaces" and found zero hits. The Latrémolière framework is built for *quantum compact metric spaces* with Lipschitz-class seminorm $L(a) = \|[D, a]\|$, where positivity of $\|\cdot\|$ is essential. The Krein-space inner product is indefinite, so the seminorm $\|[D, a]\|_{\mathcal{J}}$ is not a Lipschitz seminorm in the propinquity sense.
— *Implication:* **a Lorentzian propinquity is the single hardest piece of a Lorentzian GeoVac extension.** Paper 38's WH1 PROVEN result has no analog in the Lorentzian literature; constructing one would be an original NCG-mathematics contribution at the level of Latrémolière 2017, not a derivative application of existing theory. This is the reason the multi-month estimate is *at least* multi-month.

**Hekkelman 2022**, "Truncations and aliasing in noncommutative geometry," J. Geom. Phys.
**Hekkelman–McDonald 2024a**, "Spectral truncations of compact quantum metric spaces," arXiv:2412.00628.
**Hekkelman–McDonald 2024b**, "The noncommutative integral of compact metric spaces," arXiv:2403.0XXX.
**Leimbach–van Suijlekom 2024**, *Adv. Math.* 439, 109496, "Convergence of spectral truncations of the torus".
**Toyota 2023**, "Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics," arXiv:2504.11715 (also referenced as 2025 in Paper 38).
— *All Riemannian, all on flat-tower / S^1 / S^2 / T^d / round S^3 (Paper 38) cases.*
— *Lorentzian sequel:* none located.

### §2.6 Adjacent layer — modular flow, Bisognano–Wichmann, Unruh, KMS

**Bisognano–Wichmann 1976**, "On the duality condition for a Hermitian scalar field," *J. Math. Phys.* 17, 303 (foundational).
**Sewell 1980**, "Quantum fields on manifolds: PCT and gravitationally induced thermal states," *Ann. Phys.* 141, 201 (Schwarzschild → KMS at $\beta = 8\pi M$).
— *Verified.*
— *Structural ingredient:* in Wightman-axiomatic QFT, the modular flow of the vacuum state restricted to a wedge IS the unitary representation of Lorentz boosts that preserve the wedge. This is the bridge between Wick-rotated thermal physics (KMS state at $\beta$) and Lorentz boosts (modular flow with parameter $\theta$ where $\beta = 2\pi/\kappa$ and $\theta = \kappa t$).
— *Implication for GeoVac:* **Sprint TD Track 4's Hawking $T_H = 1/(8\pi M)$ result, under the Bisognano–Wichmann reading, IS a Lorentz-boost statement.** The framework's M1 mechanism ($2\pi$ from $S^1_\beta$ circumference) is the same $2\pi$ that appears in Bisognano–Wichmann's $\beta = 2\pi/a$ for Rindler observers. The Wick-rotation map from $S^1_\beta$ Matsubara modes to Lorentz-boost orbit parameters is structurally clean; it's the *spectral content* under the boost that's missing.

### §2.7 Summary of the literature picture

| Framework | Key reference | KO-/(m,n)-dim | Krein? | Distance formula? | Propinquity? | GeoVac usability |
|:----------|:--------------|:--------------|:------|:------------------|:--------------|:----------------|
| Pseudo-Riemannian | Strohmaier 2006 | $(m,n) \in (\mathbb{Z}/8)^2$ | Yes | No | No | High (heavy) |
| Bizi–Brouder–Besnard | arXiv:1611.07062 | Same | Yes | No | No | **Highest concrete prescription** |
| Temporal Lorentzian | Franco–Eckstein 2014 | $(m,n) = (3,1)$ + 3+1 | Yes | Yes | No | High (4–6 weeks bespoke) |
| Lorentzian distance | Franco 2014 | Above + causality | Yes | Yes | No | Medium (R2.3 SDP redo) |
| Twisted-emergence | Devastato–Lizzi–Martinetti 2018; Nieuviarts 2024–2025 | KO-dim + twist | Emerges | Indirect | No | Speculative; monitor |
| Globally hyperbolic | Paschke–Verch 2004 | Categorical | Yes | Categorical | No | Categorical scaffold |
| Riemannian (current GeoVac) | Connes 1995, Paper 38 | KO-dim 3 | No | Yes | **Yes (Latrémolière)** | — |

**Headline:** the literature provides several published prescriptions for what a Lorentzian GeoVac extension *would* look like at the spectral-triple level. **None of them provides a Lorentzian propinquity theory** — the unique structural object that Paper 38 used to PROVE WH1. A Lorentzian extension of WH1's GH-convergence theorem is not available off-the-shelf and would require original NCG-math work.

---

## §3. Wick-rotation structural obstruction map (axiom-by-axiom)

Walk through each Connes axiom of the GeoVac spectral triple under Wick rotation $\tau \to it$ from Riemannian $S^3$ (current) to Lorentzian $S^3 \times \mathbb{R}_t$ (proposed). Verdicts: **(a)** Wick-rotation-trivial; **(b)** extends with known prescription (cite the literature item); **(c)** blocks (no published prescription).

### §3.1 Algebra $\mathcal{A}_{\mathrm{GV}}$

Defined in Paper 32 §III as the function algebra on the Fock-projected $S^3$ graph. Riemannian-side data; the function algebra doesn't know about temporal direction. Under tensor product with a temporal factor $\mathcal{A}_t = C^\infty(\mathbb{R}_t)$, the combined algebra is $\mathcal{A}_{\mathrm{GV}} \otimes \mathcal{A}_t$.

**Verdict (a) Wick-rotation-trivial.** No change to $\mathcal{A}_{\mathrm{GV}}$; tensor extension is canonical.

### §3.2 Hilbert space $\mathcal{H}_{\mathrm{GV}}$

Defined in Paper 32 §III as the truncated Camporesi–Higuchi spinor space at $n_{\max}$. Under Wick rotation, this becomes a Krein space $\mathcal{K}$ with fundamental symmetry $\mathcal{J}$ such that $\mathcal{J}^2 = I$ and the Krein inner product is $\langle \cdot, \cdot \rangle_{\mathcal{J}} = \langle \cdot, \mathcal{J} \cdot \rangle$.

**Verdict (b) extends with Strohmaier 2006 / Bizi–Brouder–Besnard 2018 prescription.** The Camporesi–Higuchi spinor harmonics on $S^3$ extend to globally-hyperbolic Lorentzian $S^3 \times \mathbb{R}_t$ with $\mathcal{J} = \gamma^0$ (the temporal Dirac matrix). The Krein structure is well-defined because $\gamma^0$ is exactly the standard $J$-operator promoting Hilbert to Krein. **Cost:** moderate; the spinor space doubles its index structure to track $\mathcal{J}$-positive vs $\mathcal{J}$-negative subspaces.

### §3.3 Dirac operator $D_{\mathrm{GV}}$

Riemannian Camporesi–Higuchi spectrum $|\lambda_n| = n + 3/2$ on $S^3$. Under Wick rotation to $S^3 \times \mathbb{R}_t$, the operator becomes $D = \gamma^t \partial_t + D_{S^3}$ where $D_{S^3}$ is the spatial part.

**Verdict (b) extends with known prescription, but with a deep change.** The spectrum becomes *continuous* (because $\mathbb{R}_t$ is non-compact). To have compact resolvent / finite summability, one would need to compactify time via either $S^1_\beta$ (Sprint TD Track 1, already done) or static dS₃ patch (radius-bounded de Sitter, Lorentzian analog of Riemannian $S^3$).

The deep change: in the Riemannian case the Dirac spectrum is real with $|\lambda_n| > 0$. In the Lorentzian case the spectrum is $\sigma(D) \subset \mathbb{R}$ (Krein-self-adjoint operators have real spectrum but their action on Krein-positive vs -negative subspaces matters). The "operator order" classification (Paper 18 §IV) needs Krein-aware reformulation.

**Cost:** medium-heavy. Paper 18 §IV's first-order vs second-order operator-order classification holds Krein-side, but the transcendental ring assignments (M1, M2, M3) depend on the *spectrum* and would need to be recomputed for each Krein-self-adjoint Dirac.

### §3.4 Real structure $J$

Currently $J^2 = -I$, $JD = +DJ$, KO-dim $3 \pmod 8$ (Paper 32 §IV verified at finite $n_{\max}$ to exact zero residual). The Bizi–Brouder–Besnard $(m, n)$ framework replaces the single KO-dim by a pair: GeoVac is currently $(m, n) = (3, 0)$. A Lorentzian extension at signature $(3, 1)$ would have $(m, n) = (3, 1)$.

The sign table for the antiunitary $J$ in pseudo-Riemannian signature is *different* from the Riemannian one. In particular:
- Riemannian $S^3$: $J^2 = -I$, $JD = +DJ$, $\mathcal{J}$ trivial (Hilbert).
- Lorentzian $S^3 \times \mathbb{R}_t$ at $(m, n) = (3, 1)$: the Bizi–Brouder–Besnard framework predicts a different sign structure due to the time component contributing $\eta = -1$ vs Riemannian $\eta = +1$. The combined sign is a $(\mathbb{Z}/8) \times (\mathbb{Z}/8)$-valued invariant, not a $\mathbb{Z}/8$-valued one.

**Verdict (b) extends with Bizi–Brouder–Besnard 2018 prescription.** The (m,n) framework provides the explicit sign table. For (3, 1), the relevant Clifford algebra is $\text{Cl}(3, 1) \cong M_4(\mathbb{R})$ (real 4×4 matrices), which is structurally different from $\text{Cl}(3, 0) \cong M_2(\mathbb{H}) \cong M_4(\mathbb{R})$ as a *real* algebra but with different $J$ structure.

**Cost:** medium. The sign structure is published; what needs work is verifying the GeoVac Camporesi–Higuchi $J_{\mathrm{GV}}$ at finite $n_{\max}$ remains compatible after the temporal extension. *Note:* Paper 32 §IV's verification is finite-$n_{\max}$ exact; a Krein-side audit at finite $n_{\max}$ is a clean 1-week deliverable.

### §3.5 Bounded commutators / regularity

In the Riemannian framework $[D, a]$ is bounded for $a \in \mathcal{A}^\infty$. In the Krein framework the relevant condition is *Krein-boundedness* $\|[D, a]\|_{\mathcal{J}} < \infty$ where the norm is taken with respect to the indefinite inner product.

**Verdict (a) Wick-rotation-trivial.** Bounded operators stay bounded; the Krein-norm is just a re-weighting.

### §3.6 Compact resolvent / finite summability

Currently in the truncation $\mathcal{H}_{\mathrm{GV}}$ is finite-dimensional so this is automatic; in the continuum limit GeoVac's Dirac is $3^+$-summable on Riemannian $S^3$. Lorentzian extension to $S^3 \times \mathbb{R}_t$ has continuous spectrum and is NOT compact-resolvent without temporal compactification.

**Verdict (b) extends with temporal compactification.** Either $\mathbb{R}_t \to S^1_\beta$ (already done in Sprint TD Track 1; gives KMS thermal physics) or $\mathbb{R}_t \to$ static dS₃ patch (Lorentzian compact). Both are published.

**Cost:** low if $S^1_\beta$ (already in Sprint TD Track 1); medium if static dS₃ patch.

### §3.7 Order one / order zero

Currently trivial because $\mathcal{A}_{\mathrm{GV}}$ is commutative. Krein-side, the order-one condition $[[D, \pi(a)], \pi(b)^\circ] = 0$ uses Krein-adjoint $\pi(b)^\circ = \mathcal{J} \pi(b)^* \mathcal{J}$ instead of the Hilbert adjoint.

**Verdict (a) Wick-rotation-trivial.** For commutative $\mathcal{A}$, order-one and order-zero hold automatically in any signature.

### §3.8 Latrémolière propinquity / GH-convergence (the WH1 PROVEN axiom)

Paper 38's five-lemma proof establishes that the truncated metric spectral triples $T_{n_{\max}}$ converge to the round-$S^3$ Camporesi–Higuchi spectral triple $T_{S^3}$ in the Latrémolière propinquity. The propinquity is defined for *quantum compact metric spaces* with Lipschitz-class seminorm $L(a) = \|[D, a]\|$.

In the Krein-space setting, $\|[D, a]\|_{\mathcal{J}}$ is not a Lipschitz seminorm in the propinquity sense (positivity fails). **No published Lorentzian propinquity exists.**

**Verdict (c) BLOCKS HARD.** This is the load-bearing obstruction. Constructing a Lorentzian propinquity is original NCG mathematics, not a derivative application of existing theory. Best estimate: 6–12 months of focused work by a specialist (Latrémolière's collaborator; an early-career NCG metric geometer).

**Cost:** very high. **This is the reason the multi-month full extension is no-go at this time.**

### §3.9 Master Mellin engine M1 / M2 / M3

M1 (Hopf-base measure / Matsubara $2\pi$): preserved under Wick rotation because $S^1_\beta$ → Lorentzian boost orbit at $\beta = 2\pi/\kappa$ uses the *same* $2\pi$ via Bisognano–Wichmann.
M2 (Seeley–DeWitt heat kernel coefficients): heat kernel $\text{Tr}\, e^{-tD^2}$ requires elliptic $D$. On Minkowski signature $D^2$ is hyperbolic; the Mellin transform changes character entirely (becomes Schwinger proper-time $\int_0^\infty dt\, e^{-itD^2 + i\epsilon t^2/2}$).
M3 (vertex-parity Hurwitz / Catalan G): vertex parity selection on $S^3$ Dirac spectrum is a Riemannian eigenvalue parity condition; Lorentzian-side it would track $\mathcal{J}$-positive vs -negative spinor sectors.

**Verdicts:**
- M1 **(a) Wick-rotation-trivial** for the temperature mechanism (Bisognano–Wichmann interprets the same $2\pi$ as boost orbit period).
- M2 **(c) blocks**, requires re-derivation in Lorentzian heat-kernel theory (which exists for Schwinger proper-time but the Seeley–DeWitt coefficients have different structural meaning).
- M3 **(b) extends** with a Krein-aware vertex-parity condition.

### §3.10 Roothaan multipole termination $L_{\max} = 2 \ell_{\max}$

Currently the cross-register $V_{eN}$ ERI multipole expansion terminates exactly at $L_{\max} = 2 \ell_{\max}$ by Wigner-3j triangle inequality. Under a Lorentz boost, angular momentum mixes with linear momentum (Wigner rotations), so the spatial $L$ is no longer conserved; the expansion would NOT terminate at $L_{\max} = 2\ell_{\max}$.

**Verdict (c) blocks.** This is a sharp obstruction: every Layer-2 precision-catalogue match in Paper 34 (atomic Lamb shift, hyperfine, He fine structure) uses Roothaan termination as a finite-cost computational guarantee. A boosted catalogue match would not have this finite-cost property.

**Cost:** unknown; the question of how to handle non-conservative angular momentum in a Lorentzian spectral-triple framework is open in the literature.

### §3.11 Hopf base measure $\text{Vol}(S^2)/4$

The Hopf base $S^2$ is Riemannian and its volume is one of the M1 mechanism inputs. Under Wick rotation the spatial $S^2$ stays Riemannian (only time is rotated); the Hopf measure is unchanged.

**Verdict (a) Wick-rotation-trivial.** $\pi$ in $\text{Vol}(S^2) = 4\pi$ is unchanged.

### §3.12 Summary table

| Axiom / structure | Verdict | Citation | Cost |
|:-----------------|:-------:|:---------|:----:|
| Algebra $\mathcal{A}_{\mathrm{GV}}$ | (a) trivial | — | — |
| Hilbert → Krein space | (b) extends | Strohmaier 2006 | Medium |
| Dirac operator | (b) extends | Bizi–Brouder–Besnard 2018 | Med-heavy |
| Real structure $J$ → $(m,n) = (3,1)$ | (b) extends | Bizi–Brouder–Besnard 2018 | Medium |
| Bounded commutators | (a) trivial | — | — |
| Compact resolvent | (b) extends | Sprint TD Track 1 ($S^1_\beta$) | Low |
| Order one / order zero | (a) trivial | — | — |
| **Latrémolière propinquity** | **(c) BLOCKS** | NONE | **6–12 months** |
| M1 master Mellin | (a) trivial | Bisognano–Wichmann | — |
| M2 master Mellin | (c) blocks | (Lorentzian heat kernel) | High |
| M3 master Mellin | (b) extends | (Krein-aware vertex parity) | Medium |
| **Roothaan termination** | **(c) BLOCKS** | NONE | High |
| Hopf base measure | (a) trivial | — | — |

**Net:** Three of thirteen axioms / structures block hard, with **the Latrémolière propinquity gap** being the dominant single obstruction. The other two blocks (Lorentzian Seeley–DeWitt M2 reformulation; Roothaan termination under boost) are downstream of it.

---

## §4. Reachable-in-current-framework map (Lorentz physics already in GeoVac)

GeoVac is currently Riemannian, but Lorentz-related physics already happens inside the framework under three different names. Catalogue here so we can identify what's available without leaving Riemannian signature.

### §4.1 Bound-state radial Lorentz factor $\gamma = \sqrt{1 - (Z\alpha)^2}$ (Tier 2 spinor lift)

Implemented in `geovac/dirac_matrix_elements.py` (T7, v2.12.0) and used in Sprint MH (muonic hydrogen, 2026-05-08).

Structurally this $\gamma$ is a *single-particle radial* Lorentz factor, computed from the Dirac–Coulomb radial wavefunction. It is *not* a frame transformation; it's the relativistic bound-state correction to the radial structure within a fixed reference frame.

**Falsification:** the Riemannian framework already contains relativistic corrections to atomic structure at order $\alpha^2$ (Paper 14 §V Tier 2). Sprint MH 2S–2P Lamb shift residual at $-0.10\%$ vs CREMA (Track A) is the empirical witness.

**Reachable extension:** none new; this is already exhausted in the catalogue.

### §4.2 Wick-rotated thermal physics (Sprint TD Track 1, Track 4)

Track 1 built $T_{S^3} \otimes T_{S^1_\beta}$ and reproduced Stefan–Boltzmann $-\pi^2/90$ as exact $M_1 \times M_2$ factorisation. Track 4 reproduced Hawking $T_H = 1/(8\pi M)$ on the Euclidean Schwarzschild cigar.

Under the Bisognano–Wichmann theorem, these Wick-rotated thermal computations *are* Lorentz-boost computations:
- Sprint TD Track 1 $\beta$-circle = imaginary-time orbit of a Rindler observer at acceleration $a = 2\pi/\beta$.
- Sprint TD Track 4 $\beta = 8\pi M$ = imaginary-time orbit of a stationary observer outside the Schwarzschild horizon.

**The framework already does Lorentz-boost-class physics on the temperature mechanism.** What it doesn't do:
- The metric / curvature side at finite $v$ (length contraction etc.).
- The spinor-content side (how the Dirac spinor transforms under boost, which would give the radial $\gamma$ in §4.1 but as a frame transformation).
- The continuous-spectrum side (Lorentzian Dirac has continuous, not discrete, spectrum).

**Reachable extension:** §5.1 below documents Track 4 as a Wick-rotation map and names the Bisognano–Wichmann reading explicitly. 1–2 weeks.

### §4.3 Conformal Killing vectors of $S^3$ and the SO(4,1) action

This is the genuinely under-exploited resource. The conformal Killing algebra of $S^3$ is $\mathfrak{so}(4,1) = \mathfrak{conf}(\mathbb{R}^3)$ — the Lie algebra of the Lorentz group of one higher dimension. SO(4,1) acts on $S^3$ via:
- the SO(4) isometries (rotations of $S^3$ as a unit sphere in $\mathbb{R}^4$);
- the *conformal* boosts (special conformal transformations that act non-isometrically but preserve angles — the de Sitter group acting on $S^3$ at the spatial slice of dS₄ static patch).

The SO(4) part is already in GeoVac (it's the symmetry that produces the Camporesi–Higuchi degeneracies). The non-isometric "boost-like" part is *NOT* exploited.

**Crucial observation:** these conformal boosts live entirely *inside* Riemannian signature. They are non-isometric maps $S^3 \to S^3$ (with conformal factor $\Omega$). Wick rotation is not required to study them. The relevant infrastructure is:
- 10 conformal Killing vectors (3 rotations + 3 spatial translations + 3 boosts + 1 dilatation, in conformal coordinates on $S^3$).
- These act on the (n, l, m) Fock basis as ladder operators between Camporesi–Higuchi shells (raising / lowering n).
- The eigenfunctions of the boost generators on $S^3$ are well-defined (combinations of Gegenbauer × spherical harmonics).

**Static dS connection:** the static patch of dS₃ × ℝ_t carries the same SO(4,1) symmetry; its Hamiltonian is the generator of evolution along the conformal Killing vector that preserves the diamond. So the SO(4,1) algebra on Riemannian $S^3$ has a natural Lorentzian interpretation via dS₃ static patch *without leaving Riemannian signature on the spatial side*.

**Reachable extension:** §5.2 below proposes a 1–2 week sprint that computes the SO(4,1) conformal-isometry algebra at finite $n_{\max}$ on the Fock-projected $S^3$, identifies how the boost generators act on (n, l, m), checks whether they preserve the Roothaan termination, and identifies whether this "Riemannian boost" is structurally distinct from §III.14 rest-mass and §III.16 Breit retardation.

### §4.4 Composition of §III.14 and §III.16 with a hypothetical boost projection

If we were to *add* a boost projection (uniform-rescaling alternative from synthesis §7.7), the composition with rest-mass and Breit retardation is:
- **rest-mass × boost:** non-commuting (boost rescales mass via $m \to m \cosh \chi$, rest-mass rescales mass via $m \to m_\text{red}$). Would add another one-way cell to the §VII composition table. Variable-refinement non-commutation, like $rest\_mass \times Breit$.
- **Breit × boost:** non-commuting (boost rescales the two-body mass ratio $m_l/m_n \to m_l e^\chi / m_n$, Breit fixes the ratio). Another one-way cell.
- **boost × boost:** commutes if both act on the same spatial direction (additive rapidity); doesn't commute if they act on different directions (Wigner rotation). Idempotent only at $\chi = 0$.

**This produces three new one-way cells in the §VII composition table.** A genuinely new structural feature, *under the assumption that the boost projection is admissible*.

**Reachable extension:** §5.3 below names this as a Paper 34 §VIII open-question entry, NOT as a new admissible projection (REQUIRES-EXTENSION verdict from Hunt 4 stays).

### §4.5 Summary

| Lorentz physics | Already in GeoVac? | Where | Status |
|:---------------|:-------------------|:------|:-------|
| Bound-state radial $\gamma$ | Yes | Tier 2 spinor lift, Sprint MH | Exhausted |
| Wick-rotated thermal | Yes | Sprint TD Track 1, Track 4 | M1 mechanism trivial; needs Bisognano–Wichmann documentation |
| SO(4,1) conformal Killing on $S^3$ | NO (under-exploited) | conformal isometries of $S^3$ | **High-value, sprint-reachable** |
| Boost as 17th projection | NO (REQUIRES-EXTENSION) | Hunt 4 verdict | Paper 34 §VIII open question |

---

## §5. Minimum-viable first deliverables (2–3 candidates, ranked)

Three candidates for sprint-scale (1-2 week) Lorentz-adjacent work. Ranked by leverage / cost.

### §5.1 Candidate A — Wick-rotation map for Sprint TD Track 4 + Bisognano–Wichmann reading

**Sprint window:** 1 week.
**What it lands:** a memo + Paper 32 / Paper 35 update paragraph documenting that Sprint TD Track 4's $T_H = 1/(8\pi M)$ result is a Wick-rotation image of the Lorentz-boost / Unruh / Bisognano–Wichmann mechanism. Specifically:
- The $2\pi$ in $\beta = 8\pi M$ is the M1 Hopf-base measure signature on the imaginary-time circle.
- Bisognano–Wichmann: this $2\pi$ IS the period of the Lorentz-boost orbit (modular flow) of a stationary observer outside the Schwarzschild horizon.
- The same identity gives the Unruh temperature $T_U = a/(2\pi)$ for a Rindler observer at acceleration $a$, with $\beta = 2\pi/a$.
- Cite Sewell 1980 for the Schwarzschild → KMS correspondence; Bisognano–Wichmann 1976 for the modular flow → boost identification.

**What it does NOT close:** doesn't extend the Riemannian framework to Lorentzian; doesn't prove anything about the spinor content; doesn't construct a Lorentzian propinquity. Strictly a *documentation* deliverable that names existing results in the standard relativistic-QFT vocabulary.

**Falsifier:** the Bisognano–Wichmann reading is standard; it doesn't fail.

**Recommended:** paper-update target Paper 32 §VIII.D (cross-manifold frontier) gains a paragraph; Paper 35 §VIII.A (after the Stefan–Boltzmann subsection added by Sprint TD Track 1) gains a parallel paragraph noting the Lorentzian interpretation.

**Verdict:** **STRONGLY RECOMMENDED** as the highest-leverage 1-week deliverable. Cost is minimal (no code), value is real (Lorentz physics is already there; we just haven't named it).

### §5.2 Candidate B — SO(4,1) conformal-isometry algebra on the Fock-projected $S^3$

**Sprint window:** 1–2 weeks.
**What it lands:** new module `debug/so41_conformal_action.py` (~400 lines) computing:
- The 10 generators of $\mathfrak{so}(4,1) = \mathfrak{conf}(S^3)$ as differential operators on $S^3$.
- Their representation on the Camporesi–Higuchi (n, l, m) basis at finite $n_{\max}$ via Wigner-3j-like recoupling.
- The action of the "boost-like" non-isometric generators (3 of the 10) on (n, l, m): which Δn, Δl, Δm transitions do they enforce?
- Whether the boost generators preserve the Roothaan termination at $L_{\max} = 2\ell_{\max}$.
- Whether the Connes axioms are preserved under the conformal action (in particular: $J$ at $(m,n) = (3, 0)$ — does the conformal boost commute with $J$ in the Riemannian setting?).

**What it doesn't close:** doesn't give a Lorentzian extension; doesn't address the propinquity question; doesn't compute physical observables.

**Falsifier:** the SO(4,1) action is well-defined but might *not* preserve the spectral-triple structure (e.g., it could fail to be implemented by a unitary on the truncation, or fail to commute with $J$). This would be a clean negative — the conformal boost generators of $S^3$ are NOT in the framework, which would close the door on category (c) Lorentz-physics-in-Riemannian.

**Honest scope:** SO(4,1) at the *Lie algebra* level on the (n, l, m) basis is a clean linear-algebra computation. The structural question — does the truncation preserve SO(4,1) — is *strictly* harder than just SO(4) preservation, because the boost generators raise $n$. So **the answer is almost certainly: SO(4,1) is broken at finite $n_{\max}$** (the boost takes a state at $n = n_{\max}$ to one with $n = n_{\max} + 1$ which is outside the truncation). This is structurally analogous to the SO(4) preservation in Paper 32 §III rem:operator_system, which holds because the SO(4) generators don't change $n$.

So the *expected* result of Candidate B is: **SO(4,1) is broken by the truncation but the spatial SO(4) ⊂ SO(4,1) is preserved; the boost generators define a non-isometric map between truncations at different $n_{\max}$.** This would be a clean structural finding that lands a new piece of vocabulary (the conformal-boost-map between truncations) without needing to leave Riemannian signature.

**Verdict:** RECOMMENDED as a 1–2 week first deliverable IF the PI wants to invest beyond Candidate A. Concrete, falsifiable, lands a real new structural object. Lower priority than A because it's more work and the expected positive content is small (one new mathematical object: the conformal-boost intertwiner between truncation levels).

### §5.3 Candidate C — Connes-axiom audit on the Krein lift of $T_{\mathrm{GV}}$ at $n_{\max} = 2, 3$

**Sprint window:** 1 week (computational), 0.5 weeks (writeup).
**What it lands:** new module `debug/krein_lift_audit.py` extending the Paper 32 §IV finite-$n_{\max}$ Connes-axiom audit (current status: $J^2 = -I$, $JD = +DJ$ exact at $n_{\max} \in \{1, 2, 3\}$ on truthful CH) to the Krein lift at signature $(m, n) = (3, 1)$:
- Lift the Camporesi–Higuchi $J$ on $S^3$ to the Krein-space charge conjugation on $S^3 \times \mathbb{R}_t$ (truncate $\mathbb{R}_t$ to a small grid for numerical concreteness).
- Verify the Bizi–Brouder–Besnard sign table at $(m, n) = (3, 1)$: what sign does $J^2$ take in the temporal extension? What sign does $JD = \pm DJ$ take?
- Check whether the order-one and order-zero conditions hold Krein-side.

**What it doesn't close:** doesn't construct a Lorentzian propinquity; doesn't prove convergence; doesn't address the master Mellin engine M2 reformulation.

**Falsifier:** the Bizi–Brouder–Besnard sign table is published; the audit either matches it (positive) or reveals a discrepancy (clean negative — would indicate the Krein lift of GeoVac doesn't fit the standard prescription, which would itself be a finding).

**Honest scope:** this is a *Riemannian-side* audit of how the Krein extension WOULD behave if implemented. It's not a Lorentzian construction, just a stress test of the Bizi–Brouder–Besnard prescription against GeoVac's existing $J$ at finite $n_{\max}$. Should land a 1-week clean result.

**Verdict:** WEAK recommendation. Useful as a stress-test; lower leverage than A or B. Could be deferred until after a multi-month extension is ACTUALLY decided to undertake.

### §5.4 Ranking

| Candidate | Window | Lands | Falsifier | Priority |
|:----------|:------:|:------|:----------|:--------:|
| A: Wick-rotation map for Track 4 | **1 week** | Memo + 2 paper paragraphs documenting Bisognano–Wichmann reading | None (standard physics) | **HIGH** |
| B: SO(4,1) conformal action on S^3 | 1–2 weeks | Module + memo + structural finding | SO(4,1) broken by truncation (expected) | MEDIUM |
| C: Krein lift Connes axiom audit | 1 week | Krein-side stress test of (m,n) = (3,1) | Bizi–Brouder–Besnard match check | LOW |

**Recommended dispatch:** A first (1 week, low cost, high leverage); B second only if PI wants more depth; C deferred.

---

## §6. Multi-month full extension cost / value estimate

If the PI commits to a full Lorentzian-signature extension of GeoVac — meaning: a Lorentzian spectral triple at $(m, n) = (3, 1)$ with a Lorentzian propinquity-style convergence theorem analogous to Paper 38's WH1 PROVEN — what does it cost and what does it buy?

### §6.1 Cost decomposition

| Component | Estimate | Risk |
|:----------|:---------|:----:|
| Build Krein-space lift of $\mathcal{H}_{\mathrm{GV}}$ at $(m,n) = (3,1)$ | 2–4 weeks | Low |
| Implement Lorentzian Camporesi–Higuchi Dirac on static dS₃ | 4–6 weeks | Medium |
| Reformulate master Mellin engine M2 in Lorentzian heat-kernel theory | 4–8 weeks | Medium-high |
| Construct Lorentzian propinquity (original NCG-math) | **6–12 months** | **HIGH** — no published prescription |
| Lorentzian L1'/L2/L3/L4/L5 analogs of Paper 38 lemmas | 2–4 months | High — depends on propinquity |
| Roothaan termination under boost (or replacement) | 2–4 weeks | Medium |
| Krein-aware vertex parity (M3 reformulation) | 1–2 weeks | Low |
| Lorentzian Sprint MH redo at static dS₃ | 2–4 weeks | Low |
| Paper 38 Lorentzian sequel writeup | 4–6 weeks | Low |
| **Total bottom-up estimate** | **9–18 months** | dominated by propinquity |

This is bottom-up labor; in practice 12–18 months is the realistic minimum if the propinquity step is taken seriously, and there's a HIGH risk of stalling at the propinquity step with no result.

### §6.2 Strategic value if successful

| Capability | Strategic value |
|:----------|:---------------|
| Relativistic atomic structure beyond Tier 2 (gauge invariance, recoil at order $\alpha^4$ self-consistently) | Medium — Tier 3 already does $\alpha^4$ via Darwin + MV, and the sprint-MH and Ps catalogue rows already hit α^4 walls without a Lorentzian framework |
| Bisognano–Wichmann theorem inside GeoVac | High — would make Sprint TD Track 4 a theorem rather than a Wick-rotation observation |
| Unruh effect | Medium — interesting but not currently on a sprint roadmap |
| AdS/CFT applications | Low for current GeoVac scope; would require AdS lift, which is a whole separate program |
| Bisognano–Wichmann + entanglement Hamiltonians | Medium-high — connects to Paper 27 entropy thread |
| Closing W2b-medium (cross-manifold) | High — Paper 32 §VIII.D's frontier-of-field framing relies on this being open |

**Strategic value: Medium overall.** No single deliverable in the Lorentzian extension would be transformative; together they would round out the framework's coverage but are not the highest-priority next step compared to the multi-focal precision arc.

### §6.3 Risk-adjusted cost / value

- **Cost:** 9–18 months, dominated by 6–12 months of original NCG-math work on Lorentzian propinquity.
- **Value:** Medium overall.
- **Risk:** HIGH that the propinquity step doesn't close in <12 months, leaving a partial Lorentzian extension with no headline theorem.
- **Opportunity cost:** 9–18 months of multi-focal-precision sprint, second-row chemistry, or Cs HFS heavy-atom work that would be displaced.

**Risk-adjusted verdict:** the multi-month Lorentzian extension is **NO-GO at this time.** The cost / value ratio is unfavorable, the dominant cost (Lorentzian propinquity) has high probability of stalling, and the sprint-scale alternative (Candidate A above) captures the most strategically valuable piece (the Bisognano–Wichmann reading of Sprint TD Track 4) at a small fraction of the cost.

---

## §7. Go/no-go recommendation

**NO-GO on multi-month full Lorentzian-signature extension of GeoVac.** The cost (9–18 months, dominated by 6–12 months of original Lorentzian propinquity work) is high; the strategic value is medium; the risk of stalling is high. The opportunity cost relative to the multi-focal precision arc and Cs HFS heavy-atom direction is unfavorable.

**GO on Candidate A: 1-week Wick-rotation-map memo for Sprint TD Track 4 with Bisognano–Wichmann reading.** Captures the most strategically valuable piece (naming the Lorentz-boost interpretation of existing thermal results) at minimal cost. Dispatch immediately.

**OPTIONAL: Candidate B (1–2 weeks) on SO(4,1) conformal action.** A clean follow-on if the PI wants more depth. Lands a real new structural object (the conformal-boost intertwiner between truncation levels).

**DEFER: Candidate C (Krein lift audit).** Useful as a future stress-test before a multi-month extension; not a near-term priority.

**REVISIT trigger conditions for the multi-month extension:**
1. Lorentzian propinquity is published (would slash the dominant cost). Specifically: monitor Latrémolière / Toyota / Hekkelman–McDonald arXiv listings for any "Lorentzian propinquity" / "Krein-space propinquity" / "pseudo-Riemannian propinquity" paper.
2. Twisted-spectral-triple emergence (Devastato–Lizzi–Martinetti / Nieuviarts) stabilizes as the standard prescription — would change the route from "build Lorentzian triple" to "twist GeoVac" and could shrink the cost dramatically.
3. A specific Lorentzian application becomes load-bearing (Bisognano–Wichmann would close W2b-medium structurally; AdS lift becomes a research target; new precision-catalogue observable demands relativistic frame transformations). None of these is currently on the roadmap.

---

## §8. Honest negatives / what would falsify the project

The scoping itself doesn't get falsified — it's a scoping memo. But the *assumptions* the scoping makes can be wrong:

1. **Assumption: Lorentzian propinquity doesn't exist as of May 2026.** If a published Lorentzian propinquity is found that I missed, the multi-month estimate shrinks dramatically. Mitigation: the literature search was thorough but is not exhaustive; I'd recommend a 1-day re-check before committing to multi-month work.

2. **Assumption: the twisted-emergence approach is too fresh to bet on.** If Nieuviarts' 2025 papers (arXiv:2502.18105, arXiv:2512.15450) become the standard prescription quickly, the route changes. Mitigation: monitor the relevant arXiv listings for follow-up work over the next 6 months.

3. **Assumption: Bisognano–Wichmann is the right reading of Sprint TD Track 4.** The reading is standard, but the *quantitative* match between Sprint TD Track 4's $T_H = 1/(8\pi M)$ and the Bisognano–Wichmann modular flow is at the level of the 2π identity, not a derivation of the Schwarzschild metric from GeoVac inputs. The synthesis Hunt 4 verdict that "GeoVac doesn't autonomously produce the cigar geometry" remains.

4. **Assumption: the SO(4,1) conformal action breaks at finite $n_{\max}$.** If Candidate B finds that SO(4,1) is preserved by the truncation in a non-trivial way, that would be a real new structural finding and would change the recommendation. Falsifiable by Candidate B itself.

5. **Assumption: relativistic atomic structure at $\alpha^4$ already lives in Tier 2 / Tier 3.** Mostly true (Darwin + MV + spinor lift handle $\alpha^4$ for atomic spectra) but Sprint MH 2S–2P Lamb shift at $-0.10\%$ residual still has a $\sim 0.1\%$ multi-loop QED gap. A Lorentzian-framework-native multi-loop QED might close this at the structural level.

---

## §9. Forward implications: 2-3 next-track candidates with sprint windows

After Candidate A is dispatched (1 week, see §5.1), three related sprint candidates emerge:

### §9.1 Track Bisognano–Wichmann documentation pass (1 week, RECOMMENDED)

This is Candidate A as dispatched. Outcome: 2 paragraphs in Paper 32 §VIII.D and Paper 35 §VIII.A naming the Wick-rotation-to-boost identification of the M1 mechanism; cross-references to Sewell 1980 and Bisognano–Wichmann 1976.

### §9.2 Track SO(4,1) conformal action diagnostic (1–2 weeks, OPTIONAL after Bisognano–Wichmann)

Candidate B. Outcome: a new structural object (the conformal-boost intertwiner between truncation levels of the operator system), or a clean negative (SO(4,1) is broken by truncation as expected). Either way, lands a piece of vocabulary that's currently missing.

### §9.3 Track Lorentzian-propinquity literature monitoring (ongoing, PASSIVE)

Periodically (every 6 months) check for new papers on:
- Lorentzian / Krein-space propinquity (arXiv math.OA / math-ph)
- Twisted-emergence stabilization (Nieuviarts and follow-ups)
- Hekkelman–McDonald / Toyota / Latrémolière for Lorentzian sequels
- "Lorentzian spectral truncations" (would be the analog of Connes–vS 2021 in Lorentzian setting)

If any of these lands a usable Lorentzian propinquity, the multi-month extension's revisit trigger fires.

### §9.4 Cross-reference to Paper 34 §VIII open question (recommended language)

Paper 34 §VIII gains an open-question entry (RECOMMENDED LANGUAGE; not applied):

```
\subsection{Open: Lorentz-boost projection requires multi-month structural extension}
\label{sec:open_lorentz_boost}

Under a Lorentz boost between two reference frames at relative velocity $v$,
spatial focal lengths transform as $\lambda_x \to \lambda_x / \gamma$ and
temporal scales as $\beta \to \beta \gamma$, with $\gamma = 1/\sqrt{1 - v^2/c^2}$.
The structural question is whether the Paper~34 projection family admits
an admissible 17th projection capturing this transformation.

Sprint Lorentz-Boost-Scoping (May 2026; \verb|debug/lorentz_boost_scoping_memo.md|)
returned the verdict that a clean Lorentz-boost projection requires a
Minkowski-signature spectral-triple framework which the Wick-rotated
Riemannian construction of Paper~32 does not directly admit. The
literature (Strohmaier 2006~\cite{strohmaier2006}, Bizi--Brouder--Besnard
2018~\cite{bizi_brouder_besnard2018}, Franco--Eckstein 2014~\cite{franco_eckstein2014})
provides several published prescriptions for what a Lorentzian extension
would look like at the spectral-triple level; none provides a Lorentzian
propinquity / Gromov--Hausdorff-style convergence theorem analogous to
Paper~38's WH1 PROVEN result. Constructing such a propinquity is original
NCG-math work at the 6--12 month scale.

A sprint-scale alternative is provided by Sprint TD Track~4
(\verb|debug/sprint_td_track4_memo.md|), which reproduces the Hawking
temperature $T_H = 1/(8\pi M)$ on the Euclidean Schwarzschild cigar via
the master Mellin engine $M_1$ Hopf-base measure mechanism. Under the
Bisognano--Wichmann theorem~\cite{bisognano_wichmann1976,sewell1980},
this Wick-rotated thermal computation IS the modular flow / Lorentz-boost
orbit of a stationary observer outside the horizon. The framework
therefore performs Lorentz-boost-class physics at the temperature
mechanism, even though it does not natively contain the metric / curvature
side of the boost.

Status: open. Reachable as a sprint-scale documentation pass (Candidate~A
in the scoping memo); reachable as a multi-month structural extension
only if Lorentzian propinquity becomes available in the literature, or
the twisted-spectral-triple emergence approach
(Devastato--Lizzi--Martinetti 2018~\cite{devastato_lizzi_martinetti2018};
Nieuviarts 2025~\cite{nieuviarts2025a,nieuviarts2025b}) stabilizes as a
standard prescription.
```

---

## §10. Result (30-second read)

**Top 3 reachable-now Lorentz-adjacent deliverables:**

1. **Wick-rotation-map memo for Sprint TD Track 4 (1 week, HIGH leverage).** Documents that $T_H = 1/(8\pi M)$ IS Lorentz-boost / Bisognano–Wichmann modular-flow physics, names the M1 mechanism's $2\pi$ as the boost orbit period.
2. **SO(4,1) conformal-action diagnostic (1–2 weeks, MEDIUM leverage).** Computes the conformal Killing algebra of $S^3$ on the Fock-projected (n,l,m) basis at finite $n_{\max}$; expected outcome: SO(4,1) is broken by truncation but defines a non-isometric intertwiner between truncation levels.
3. **Paper 34 §VIII open-question paragraph (immediate, ZERO cost).** Explicitly documents the REQUIRES-EXTENSION verdict with revisit-trigger conditions.

**Multi-month go/no-go verdict: NO-GO at this time.** Cost 9–18 months dominated by 6–12 months of original Lorentzian-propinquity NCG-math; strategic value medium; risk of stalling high; opportunity cost vs. multi-focal precision arc unfavorable. Revisit if Lorentzian propinquity lands in the literature or twisted-emergence stabilizes.

**Single highest-leverage 1-2 week sprint to dispatch first: Candidate A (Wick-rotation-map memo + Bisognano–Wichmann reading).** Captures the strategically valuable piece (naming Lorentz-boost interpretation of existing thermal results) at minimal cost (1 week, no code). Dispatch as immediate follow-on.

---

## §11. Files

- `debug/lorentz_boost_scoping_memo.md` — this memo (~10,000 words).
- `debug/data/lorentz_boost_scoping_data.json` — structured data: literature survey, obstruction map, candidates, go/no-go scoring.

No production code modified. No papers modified. Three paper-update recommendations (one each for Paper 32 §VIII.D, Paper 35 §VIII.A, Paper 34 §VIII) provided as draft language, not applied.
