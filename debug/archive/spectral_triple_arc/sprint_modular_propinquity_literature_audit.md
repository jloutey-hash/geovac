# Sprint Modular Propinquity — Literature Audit

**Date:** 2026-05-23.
**Sprint position:** Track 0 of the dual modular propinquity parallel sprint (Tracks M-X, M-Y, M-Z, M-H1 are downstream applications).
**Mandate:** locate Latrémolière's modular and dual-modular propinquity papers; extract canonical definitions and theorems; assess substantive vs surface-level overlap with the framework's Roothaan autopsy methodology; recommend which paper anchors each application track.

---

## Executive summary

Latrémolière's **modular propinquity** (arXiv:1608.04881, *Dissertationes Math.* 544, 2019) and **dual modular propinquity** (arXiv:1811.04534, *J. Noncommut. Geom.* 15, 2021) are noncommutative Gromov–Hausdorff distances on the class of **metrized quantum vector bundles** (MQVB) — Hilbert C*-modules $\mathcal{M}$ over quasi-Leibniz QCMS $(\mathfrak{A}, \mathsf{L}_{\mathfrak{A}})$ equipped with a **D-norm** (the module analog of a Lipschitz seminorm; the D-norm generalizes the operator norm associated with a metric connection on a hermitian vector bundle). The non-dual ("bridge") version uses *modular bridges* and is shown to be a pseudo-metric up to full modular quantum isometry. The dual version uses *modular tunnels* (which are explicit isometric embeddings into a common MQVB), drops one Leibniz-type axiom relative to the bridge version, and is **complete** — the central reason for its introduction. The covariant modular propinquity from these two papers feeds directly into Latrémolière's spectral propinquity on metric spectral triples (arXiv:1811.10843, *Adv. Math.* 404, 2022), which is the natural setting for Sprint H1 (Higgs as inner fluctuation) and the GeoVac math.OA arc more broadly.

**The "dual" in "dual modular propinquity" refers to dual modes of construction (tunnels = explicit isometric embeddings vs bridges = pivot-based numerical quantification) — it does NOT refer to dual Hilbert C*-modules, Morita equivalence, or Roothaan-style additive decomposition of observables.** The PI's flagged vocabulary similarity to the Roothaan autopsy methodology is **surface-level only**: both decompose into pieces, but the math.OA "dual" is tunnel-vs-bridge construction duality, while the GeoVac "autopsy" is additive component decomposition across Paper 34 projections. The substantive overlap is at the framework level — Hilbert-module pin-state pairings — not at the decomposition-structure level.

**Recommendations for the four parallel application tracks:**
- **M-X (Sturmian-FCI Latrémolière error model):** anchor on Latrémolière arXiv:2512.03573 (already deep-read) for the pinned-QLCMS skeleton, plus arXiv:1608.04881 for D-norm vocabulary in case the Sturmian truncation is reformulated as a metrized quantum vector bundle (module-level rather than algebra-level).
- **M-Y (NaH pin-state W1c diagnostic):** anchor on arXiv:2512.03573 (pin-state framework) plus arXiv:1811.04534 (dual modular propinquity for completeness) since the W1c diagnostic is comparing distinct pin-state choices on the same MQVB.
- **M-Z (Bethe log Latrémolière error model):** anchor on arXiv:2512.03573 plus arXiv:1811.10843 (spectral propinquity for metric spectral triples) — Paper 36's bound-state QED uses Dirac S³ data, which is exactly the spectral-triple regime.
- **M-H1 (Higgs as inner fluctuation, POSITIVE-THIN):** anchor on arXiv:1608.04881 + arXiv:1811.04534 + arXiv:1811.10843 — the Higgs sits in the inner-fluctuation $\omega = a\,dD\,b$ slot of the spectral triple, and the metrized-quantum-vector-bundle / D-norm framework is the natural NCG vocabulary for it.

---

## §1. Literature search results

Six relevant Latrémolière papers identified in the modular / dual modular line, plus one independent application paper:

| arXiv | Authors | Title | Year | Venue | Status |
|:------|:--------|:------|:----:|:------|:------:|
| 1608.04881 | Latrémolière | The Modular Gromov-Hausdorff Propinquity | 2016 (v5 Feb 2018) | *Dissertationes Math.* 544 (2019), 70 pp | Published ✓ |
| 1703.07073 | Latrémolière | Heisenberg Modules over Quantum 2-tori are metrized quantum vector bundles | 2017 | *Canad. J. Math.* | Published ✓ |
| 1803.06601 | Latrémolière | Convergence of Heisenberg Modules over Quantum 2-tori for the Modular Gromov-Hausdorff Propinquity | 2018 | *J. Operator Theory* 84 (2020), 211–237 | Published ✓ |
| 1811.04534 | Latrémolière | The Dual Modular Gromov-Hausdorff Propinquity and Completeness | 2018 (v3 Jul 2020) | *J. Noncommut. Geom.* 15 (2021), 347–398 | Published ✓ |
| 1811.10843 | Latrémolière | The Gromov-Hausdorff propinquity for metric Spectral Triples | 2018 | *Adv. Math.* 404 (2022), 108393, 56 pp | Published ✓ |
| 2112.11000 | Latrémolière | Continuity of the Spectrum of Dirac Operators of Spectral Triples for the Spectral Propinquity | 2021 | *Math. Ann.* (2023) | Published ✓ |
| 2211.11107 | Aguilar, Yu | The Fell topology and the modular Gromov-Hausdorff propinquity | 2022 | (preprint) | Preprint |
| 2512.03573 | Latrémolière | The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces | Dec 2025 | (preprint) | Preprint — covered by Phase A.2' deep-read |
| 2504.11715 | Farsi, Latrémolière | Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics | Apr 2025 | (preprint) | Preprint |

**Note on 2512.03573:** the Dec 2025 paper introduces a *pointed* (not modular) quantum-metric framework — that paper is about C*-algebras with a pin state, not Hilbert modules with a D-norm. It is structurally complementary to the modular line: 2512.03573 handles non-compact algebras with a distinguished state; the modular line handles modules over algebras with a Lipschitz seminorm. The two could be combined ("pointed modular propinquity for proper Hilbert C*-modules") but no published paper does so as of May 2026.

**One-paragraph summaries:**

- **arXiv:1608.04881 (Modular propinquity).** Introduces the **metrized quantum vector bundle** $(\mathcal{M}, \langle\cdot,\cdot\rangle_{\mathcal{M}}, \mathsf{D}_{\mathcal{M}}, \mathfrak{A}, \mathsf{L}_{\mathfrak{A}})$ — a left Hilbert $\mathfrak{A}$-module with a densely defined **D-norm** $\mathsf{D}_{\mathcal{M}}$ generalizing the operator norm given by a metric connection on a hermitian vector bundle. Defines the **modular bridge** (a unital C*-algebra $\mathfrak{D}$, a "pivot" $x \in \mathfrak{D}$ with $\|x\|_{\mathfrak{D}} = 1$, two unital $*$-monomorphisms $\pi_{\mathfrak{A}}, \pi_{\mathfrak{B}}$, plus families of anchors $(\omega_j)$ in the unit ball of $\mathsf{D}_{\mathcal{M}}$ and co-anchors $(\eta_j)$ in the unit ball of $\mathsf{D}_{\mathcal{N}}$). Builds the **modular propinquity** $\Lambda^{\mathrm{mod}}_{\mathcal{B}}(\Omega_{\mathfrak{A}}, \Omega_{\mathfrak{B}})$ as the infimum of trek lengths over compatible modular bridges. Main result: the modular propinquity is a pseudo-metric, with distance zero exactly between fully modular quantum isometric MQVBs. Application: continuity of Heisenberg modules over quantum 2-tori for the modular propinquity.

- **arXiv:1811.04534 (Dual modular propinquity).** Introduces the **dual modular propinquity** $\Lambda^{*\mathrm{mod}}_{\mathcal{T}}(\Omega_{\mathfrak{A}}, \Omega_{\mathfrak{B}})$ as the infimum of $\chi(\tau)$ over **modular tunnels** $\tau$ — pairs of full modular quantum isometries from a common MQVB onto $\Omega_{\mathfrak{A}}$ and $\Omega_{\mathfrak{B}}$. Modular bridges from arXiv:1608.04881 are shown to give modular tunnels, so the dual modular propinquity is dominated by the modular propinquity. **The dual version is weaker but COMPLETE — this is the central reason for its introduction.** The paper also extends MQVBs to **metrical** quantum vector bundles (a QCMS acting on a Hilbert module over a possibly DIFFERENT QCMS), which is the framework used in arXiv:1811.10843 for spectral triples.

- **arXiv:1811.10843 (Spectral propinquity).** Defines the **spectral propinquity** on metric spectral triples $(\mathfrak{A}, \mathcal{H}, D)$ using the **covariant modular propinquity** (a group-action-equivariant version of the modular propinquity) as the load-bearing component. The Dirac operator $D$ is encoded as a metric structure on $\mathcal{H}$ via the seminorm $\mathsf{L}(a) = \|[D, a]\|$; the spectral triple is read as a metrical quantum vector bundle ($\mathcal{H}$ is a Hilbert module over $\mathbb{C}$, with $\mathfrak{A}$ acting on it). Main theorem: the spectral propinquity is null iff the Dirac operators are unitarily equivalent AND the underlying C*-algebras are $*$-isomorphic.

- **arXiv:2112.11000 (Spectrum continuity).** Proves that if a sequence of metric spectral triples converges in the spectral propinquity, then the spectra of the Dirac operators converge. This is the downstream functoriality theorem.

- **arXiv:1703.07073 + 1803.06601 (Heisenberg modules application).** Verifies that Heisenberg modules over quantum 2-tori are MQVBs and proves their continuity in the modular propinquity. This is the canonical worked example.

---

## §2. Canonical definitions

### 2.1 Quasi-Leibniz quantum compact metric space (background — Latrémolière 1608.04881 Def 2.4)

A **quasi-Leibniz quantum compact metric space** $(\mathfrak{A}, \mathsf{L})$ is a unital C*-algebra $\mathfrak{A}$ with a seminorm $\mathsf{L}$ on a dense Jordan–Lie subalgebra $\mathrm{dom}(\mathsf{L}) \subset \mathfrak{sa}(\mathfrak{A})$ satisfying:
1. $\{a \in \mathrm{dom}(\mathsf{L}) : \mathsf{L}(a) = 0\} = \mathbb{R} \cdot 1_{\mathfrak{A}}$;
2. The Monge–Kantorovich metric $\mathrm{mk}_{\mathsf{L}}(\varphi, \psi) = \sup\{|\varphi(a) - \psi(a)| : \mathsf{L}(a) \le 1\}$ metrizes the weak-* topology on $\mathcal{S}(\mathfrak{A})$;
3. $\mathsf{L}$ is lower semicontinuous;
4. For some admissible $F$: $\max\{\mathsf{L}(a \circ b), \mathsf{L}(\{a, b\})\} \le F(\|a\|, \|b\|, \mathsf{L}(a), \mathsf{L}(b))$ where $\circ$ is the Jordan product and $\{\cdot,\cdot\}$ is the Lie product (this is the *quasi-Leibniz* axiom).

When $F(x, y, l_x, l_y) = x l_y + y l_x$ this reduces to the **Leibniz** inequality. The QCMS substrate is what GeoVac inherits from Paper 38 (Camporesi-Higuchi spectral triple on $S^3$).

### 2.2 Hilbert C*-module and module morphism (Latrémolière 1608.04881 Defs 3.1, 3.4, 3.5)

A **left Hilbert module** $(\mathcal{M}, \langle\cdot,\cdot\rangle_{\mathcal{M}})$ over a C*-algebra $\mathfrak{A}$ is a left $\mathfrak{A}$-module equipped with a sesquilinear map $\langle\cdot,\cdot\rangle_{\mathcal{M}} : \mathcal{M} \times \mathcal{M} \to \mathfrak{A}$ satisfying:
1. $\langle a\omega, \eta \rangle_{\mathcal{M}} = a \langle\omega, \eta\rangle_{\mathcal{M}}$;
2. $\langle\omega, \eta\rangle_{\mathcal{M}}^* = \langle\eta, \omega\rangle_{\mathcal{M}}$;
3. $\langle\omega, \omega\rangle_{\mathcal{M}} \ge 0$;
4. $\langle\omega, \omega\rangle_{\mathcal{M}} = 0$ iff $\omega = 0$;
5. (Hilbert) $\mathcal{M}$ is complete for $\|\omega\|_{\mathcal{M}} = \sqrt{\|\langle\omega, \omega\rangle_{\mathcal{M}}\|_{\mathfrak{A}}}$.

A **module morphism** $(\Theta, \theta)$ from $(\mathcal{M}, \langle\cdot,\cdot\rangle_{\mathcal{M}})$ to $(\mathcal{N}, \langle\cdot,\cdot\rangle_{\mathcal{N}})$ is a pair where $\theta : \mathfrak{A} \to \mathfrak{B}$ is a $*$-morphism and $\Theta : \mathcal{M} \to \mathcal{N}$ is $\mathbb{C}$-linear with $\Theta(a\omega) = \theta(a)\Theta(\omega)$ and $\langle\Theta(\omega), \Theta(\eta)\rangle_{\mathcal{N}} = \theta(\langle\omega, \eta\rangle_{\mathcal{M}})$ (the inner-product condition appears in the dual modular paper, Def 2.6).

### 2.3 D-norm: the canonical module Lipschitz seminorm (Latrémolière 1608.04881 Def 3.8)

**Definition (D-norm / metrized quantum vector bundle).** A $(F, G, H)$-**metrized quantum vector bundle** $(\mathcal{M}, \langle\cdot,\cdot\rangle_{\mathcal{M}}, \mathsf{D}_{\mathcal{M}}, \mathfrak{A}, \mathsf{L}_{\mathfrak{A}})$ is given by:
1. An $F$-quasi-Leibniz QCMS $(\mathfrak{A}, \mathsf{L}_{\mathfrak{A}})$ called the **base quantum space**;
2. A left Hilbert $\mathfrak{A}$-module $(\mathcal{M}, \langle\cdot,\cdot\rangle_{\mathcal{M}})$;
3. A norm $\mathsf{D}_{\mathcal{M}}$ defined on a dense $\mathbb{C}$-subspace $\mathrm{dom}(\mathsf{D}_{\mathcal{M}}) \subset \mathcal{M}$ such that:
   - (a) $\|\cdot\|_{\mathcal{M}} \le \mathsf{D}_{\mathcal{M}}$ (D-norm dominates Hilbert norm);
   - (b) the set $\{\omega \in \mathcal{M} : \mathsf{D}_{\mathcal{M}}(\omega) \le 1\}$ is compact for $\|\cdot\|_{\mathcal{M}}$ (closed unit-ball compactness);
   - (c) **inner quasi-Leibniz inequality:** $\mathsf{D}_{\mathcal{M}}(a\omega) \le G(\|a\|_{\mathfrak{A}}, \mathsf{L}_{\mathfrak{A}}(a), \mathsf{D}_{\mathcal{M}}(\omega))$;
   - (d) **modular quasi-Leibniz inequality:** $\max\{\mathsf{L}_{\mathfrak{A}}(\Re\langle\omega,\eta\rangle_{\mathcal{M}}), \mathsf{L}_{\mathfrak{A}}(\Im\langle\omega,\eta\rangle_{\mathcal{M}})\} \le H(\mathsf{D}_{\mathcal{M}}(\omega), \mathsf{D}_{\mathcal{M}}(\eta))$.

The norm $\mathsf{D}_{\mathcal{M}}$ is the **D-norm**. The (b) compactness and (a) domination conditions are the module analogs of Rieffel's QCMS axioms; (c) and (d) are the module analogs of the Leibniz inequality.

**Classical prototype (motivation from the paper introduction):** if $V$ is a hermitian complex vector bundle on a compact connected Riemannian manifold $M$ with metric connection $\nabla$, then the $C(M)$-module $\Gamma$ of continuous sections of $V$ is a left Hilbert module, and the D-norm is the operator norm $\mathsf{D}_{\mathcal{M}}(\omega) := \|\omega\|_{\mathcal{M}} + \|\nabla \omega\|_{\mathcal{M}}$ (or analog) associated with the metric connection. The D-norm is the abstract NCG generalization of "connection + module norm."

### 2.4 Dual modular paper's RELAXED MQVB definition (Latrémolière 1811.04534 Def 2.8)

The dual modular paper **drops axiom (d)** (the modular quasi-Leibniz inequality) — it shows this axiom is not needed for the dual modular propinquity construction. The relaxed definition uses only:
- $(\mathfrak{A}, \mathsf{L}_{\mathfrak{A}})$ an $F$-quasi-Leibniz QCMS;
- $(\mathcal{M}, \langle\cdot,\cdot\rangle_{\mathcal{M}})$ a left Hilbert $\mathfrak{A}$-module;
- $\mathsf{D}_{\mathcal{M}}$ satisfying (a) domination, (b) closed unit-ball compactness, (c) inner quasi-Leibniz only.

This relaxation matters because modular tunnels obtained from modular bridges do not automatically satisfy axiom (d). The relaxed definition is later extended (Def 2.12) to **metrical quantum vector bundles** where the base QCMS $\mathfrak{B}$ acts on the module $\mathcal{M}$ via a $*$-morphism into the adjoinable-operator algebra — this allows $\mathcal{M}$ to be a Hilbert module over a DIFFERENT C*-algebra (e.g., $\mathcal{M}$ a Hilbert module over $\mathbb{C}$, with $\mathfrak{A}$ acting on it) and is the framework used for spectral triples.

### 2.5 Modular bridge (Latrémolière 1608.04881 Def 4.4)

**Definition (modular bridge).** Let $\Omega_{\mathfrak{A}} = (\mathcal{M}, \langle\cdot,\cdot\rangle_{\mathcal{M}}, \mathsf{D}_{\mathcal{M}}, \mathfrak{A}, \mathsf{L}_{\mathfrak{A}})$ and $\Omega_{\mathfrak{B}} = (\mathcal{N}, \langle\cdot,\cdot\rangle_{\mathcal{N}}, \mathsf{D}_{\mathcal{N}}, \mathfrak{B}, \mathsf{L}_{\mathfrak{B}})$ be two MQVBs. A **modular bridge** $\gamma = (\Omega_{\mathfrak{A}}, \Omega_{\mathfrak{B}}, \mathfrak{D}, x, \pi_{\mathfrak{A}}, \pi_{\mathfrak{B}}, (\omega_j)_{j \in J}, (\eta_j)_{j \in J})$ from $\Omega_{\mathfrak{A}}$ to $\Omega_{\mathfrak{B}}$ consists of:
1. A unital C*-algebra $\mathfrak{D}$;
2. A **pivot** $x \in \mathfrak{D}$ with $\mathcal{S}_1(\mathfrak{D}|x) \ne \emptyset$ and $\|x\|_{\mathfrak{D}} = 1$;
3. Two unital $*$-monomorphisms $\pi_{\mathfrak{A}} : \mathfrak{A} \hookrightarrow \mathfrak{D}$ and $\pi_{\mathfrak{B}} : \mathfrak{B} \hookrightarrow \mathfrak{D}$;
4. A nonempty set $J$;
5. **Anchors** $(\omega_j)_{j \in J} \subset \mathcal{D}_1(\Omega_{\mathfrak{A}})$ (unit ball of $\mathsf{D}_{\mathcal{M}}$);
6. **Co-anchors** $(\eta_j)_{j \in J} \subset \mathcal{D}_1(\Omega_{\mathfrak{B}})$.

The bridge has a **length** $\lambda(\gamma) = \max\{\varsigma(\gamma), \varrho(\gamma)\}$ where $\varsigma(\gamma)$ is the height (Hausdorff distance between state-space images) and $\varrho(\gamma) = \max\{\varrho_{\flat}(\gamma), \varrho^{\sharp}(\gamma) + \varpi(\gamma)\}$ is the reach (combination of basic bridge reach, modular reach measuring anchor/co-anchor pairing distortion via the **deck seminorm** $\mathsf{dn}_\gamma$, and imprint measuring how densely anchors fill the unit ball).

### 2.6 Modular tunnel (Latrémolière 1811.04534 Def 2.20 referencing earlier work)

**Definition (F-tunnel).** Let $(\mathfrak{A}_1, \mathsf{L}_1)$ and $(\mathfrak{A}_2, \mathsf{L}_2)$ be two $F$-quasi-Leibniz QCMS. An $F$-**tunnel** $\tau = (\mathfrak{D}, \mathsf{L}_{\mathfrak{D}}, \pi_1, \pi_2)$ from $(\mathfrak{A}_1, \mathsf{L}_1)$ to $(\mathfrak{A}_2, \mathsf{L}_2)$ is:
1. An $F$-quasi-Leibniz QCMS $(\mathfrak{D}, \mathsf{L}_{\mathfrak{D}})$;
2. For $j \in \{1, 2\}$, a $*$-morphism $\pi_j : (\mathfrak{D}, \mathsf{L}_{\mathfrak{D}}) \twoheadrightarrow (\mathfrak{A}_j, \mathsf{L}_j)$ which is a **quantum isometry** (i.e., $\pi_j$ is a $*$-epimorphism with $\mathsf{L}_j(b) = \inf\{\mathsf{L}_{\mathfrak{D}}(a) : \pi_j(a) = b\}$).

A **modular tunnel** between MQVBs is the same construction lifted to the module level: pairs of full modular quantum isometries from a common MQVB onto $\Omega_{\mathfrak{A}}$ and $\Omega_{\mathfrak{B}}$. The **extent** $\chi(\tau)$ is the Hausdorff distance (in the Monge-Kantorovich metric) between the state-space images of $\pi_1, \pi_2$.

**The contrast between bridge and tunnel:** a bridge embeds via two $*$-monomorphisms $\pi_{\mathfrak{A}}, \pi_{\mathfrak{B}}$ into a *containing* algebra $\mathfrak{D}$; a tunnel projects via two $*$-epimorphisms $\pi_1, \pi_2$ from a *common* algebra $\mathfrak{D}$ onto the two targets. Bridges quantify "embedding distortion" via the pivot + anchors construction; tunnels quantify "projection distortion" via the Hausdorff distance between state-space images. **The two are dual notions of isometric embedding** — embeddings INTO vs surjections ONTO — which is what "dual" refers to in "dual modular propinquity."

### 2.7 Connection on a module (Latrémolière 1608.04881 introduction)

The classical prototype the D-norm generalizes is the **metric connection** $\nabla$ on a hermitian vector bundle over a compact Riemannian manifold $M$. A connection gives a covariant derivative $\nabla : \Gamma(V) \to \Gamma(V \otimes T^* M)$, and the D-norm is then $\mathsf{D}_{\mathcal{M}}(\omega) := \|\omega\|_{\mathcal{M}} + \|\nabla \omega\|_{\Gamma(V \otimes T^*M)}$ (or appropriate analog). Latrémolière's framework does NOT require the full strength of a connection (parallel transport, curvature, etc.) — only the operator norm $\mathsf{D}_{\mathcal{M}}$ derived from it. This is the sense in which "the D-norm encodes some aspects of the connection."

In the noncommutative setting, "metric connections on Hilbert modules" exist under generous conditions (per 1608.04881 introduction), but uniqueness is unclear. The framework adopts the perspective that the metric information needed to work with vector bundles is: (i) the inner product on the module of sections, (ii) the D-norm encoding connection-derived operator-norm data. This is genuinely more general than asking for a connection.

### 2.8 Morita equivalence in the propinquity sense

**Not directly addressed in the modular propinquity papers.** The papers treat Hilbert C*-modules as the framework's objects but do not formalize Morita equivalence as a propinquity-respecting equivalence. The implicit notion is **full modular quantum isometry** (the Hilbert-module isomorphism that preserves both the base-algebra QCMS structure and the D-norm) — this is the "isomorphism" in the metric category. Morita equivalence (which preserves the algebra category up to module category equivalence) is a strictly weaker notion than full modular quantum isometry.

---

## §3. Main theorem statements

### 3.1 Modular propinquity main theorem (Latrémolière 1608.04881 §§5-6)

**Definition 5.6 (modular Gromov-Hausdorff $\mathcal{B}$-propinquity).** Let $\mathcal{C}$ be a nonempty class of $(F, G, H)$-MQVBs and $\mathcal{B}$ a class of modular bridges compatible with $\mathcal{C}$. The modular Gromov-Hausdorff $\mathcal{B}$-propinquity between $\Omega_{\mathfrak{A}}, \Omega_{\mathfrak{B}} \in \mathcal{C}$ is:
$$\Lambda^{\mathrm{mod}}_{\mathcal{B}}(\Omega_{\mathfrak{A}}, \Omega_{\mathfrak{B}}) = \inf\{\lambda(\Gamma) : \Gamma \in \mathrm{Treks}[\Omega_{\mathfrak{A}} \stackrel{\mathcal{B}}{\longrightarrow} \Omega_{\mathfrak{B}}]\}.$$

**Proposition 5.10 (finiteness + dominates quantum propinquity).** $\Lambda^{\mathrm{mod}}_{\mathcal{B}}(\Omega_{\mathfrak{A}}, \Omega_{\mathfrak{B}}) < \infty$; moreover the modular propinquity dominates the underlying quantum propinquity between the base QCMS.

**Proposition 5.15 (triangle inequality + symmetry).** $\Lambda^{\mathrm{mod}}_{\mathcal{B}}$ is symmetric and satisfies the triangle inequality.

**Proposition 5.16 (full modular quantum isometry $\Rightarrow$ distance zero).** If there exists a full modular quantum isometry $(\theta, \Theta)$ from $\Omega_{\mathfrak{A}}$ to $\Omega_{\mathfrak{B}}$, then $\Lambda^{\mathrm{mod}}_{\mathcal{B}}(\Omega_{\mathfrak{A}}, \Omega_{\mathfrak{B}}) = 0$.

**Main coincidence theorem (§6 of arXiv:1608.04881).** The modular propinquity is null between two MQVBs iff they are fully modular quantum isometric. (This is the converse to Proposition 5.16; the full coincidence theorem is proved in §6 over multiple pages.)

**Combined:** $\Lambda^{\mathrm{mod}}_{\mathcal{B}}$ is a **pseudo-metric** on the class of MQVBs, and is a **metric up to full modular quantum isometry** — the natural NCG analog of Gromov-Hausdorff distance on classes of metric spaces.

### 3.2 Dual modular propinquity completeness theorem (Latrémolière 1811.04534)

**Theorem-Definition 2.19 (dual GH propinquity for QCMS).** For $\mathcal{T}$ an appropriate class of $F$-tunnels, $\Lambda^*_{\mathcal{T}}((\mathfrak{A}_1, \mathsf{L}_1), (\mathfrak{A}_2, \mathsf{L}_2)) := \inf\{\chi(\tau) : \tau \in \mathrm{Tunnels}[(\mathfrak{A}_1, \mathsf{L}_1) \stackrel{\mathcal{T}}{\longrightarrow} (\mathfrak{A}_2, \mathsf{L}_2)]\}$ is a complete metric on $F$-QCMS up to full quantum isometry.

**Main theorem of arXiv:1811.04534 (extension to modules):** the **dual modular propinquity** $\Lambda^{*\mathrm{mod}}_{\mathcal{T}}(\Omega_{\mathfrak{A}}, \Omega_{\mathfrak{B}})$ defined via modular tunnels is:
1. A pseudo-metric on the class of $(F, H)$-MQVBs;
2. Null iff $\Omega_{\mathfrak{A}}$ and $\Omega_{\mathfrak{B}}$ are fully modular quantum isometric;
3. Dominated by the modular propinquity $\Lambda^{\mathrm{mod}}_{\mathcal{B}}$ (modular bridges give modular tunnels);
4. **Complete** — Cauchy sequences of MQVBs have a limit MQVB (up to full modular quantum isometry).

**The completeness in (4) is the central reason for introducing the dual modular propinquity.** The (non-dual) modular propinquity from arXiv:1608.04881 is not known to be complete.

### 3.3 Compatibility with spectral triples (Latrémolière 1811.10843)

The metrical-quantum-vector-bundle extension of MQVBs (arXiv:1811.04534 Def 2.12) accommodates spectral triples: a metric spectral triple $(\mathfrak{A}, \mathcal{H}, D)$ defines a metrical quantum vector bundle where:
- $\mathcal{H}$ is a Hilbert module over $\mathbb{C}$ (the inner product is $\mathbb{C}$-valued, not $\mathfrak{A}$-valued);
- $\mathfrak{A}$ acts on $\mathcal{H}$ via a $*$-morphism into bounded operators;
- The Lipschitz seminorm $\mathsf{L}(a) = \|[D, a]\|_{\mathcal{H}}$ is the QCMS seminorm on $\mathfrak{A}$;
- The D-norm on $\mathcal{H}$ encodes the Dirac operator data.

The **spectral propinquity** (arXiv:1811.10843) is built using the **covariant modular propinquity** (group-equivariant version) on these metrical-quantum-vector-bundle representations of spectral triples. Main result (Theorem of arXiv:1811.10843): the spectral propinquity is null iff the Dirac operators are unitarily equivalent AND the underlying C*-algebras are $*$-isomorphic. Downstream functoriality (arXiv:2112.11000): convergence in the spectral propinquity implies convergence of Dirac spectra.

---

## §4. Scope of applicability

| Aspect | Modular propinquity (1608.04881) | Dual modular propinquity (1811.04534) | Spectral propinquity (1811.10843) | Pointed-QLCMS (2512.03573) |
|:-------|:--------------------------------:|:-------------------------------------:|:---------------------------------:|:--------------------------:|
| Compact algebra | required | required | required | NO (non-compact OK) |
| Unital | required | required | required | non-unital OK |
| Hilbert vs C*-module | Hilbert C*-module over base QCMS | Hilbert C*-module (relaxed axioms) | metrical QVB (module over ℂ + algebra action) | C*-algebra only, no module structure |
| Pin state required | NO | NO | NO | YES (Def 1.26) |
| Tracial requirement on base | NO | NO | NO | NO (but Fortet-Mourier compatibility) |
| Completeness of metric | unknown | YES (load-bearing) | inherited from dual modular | YES (metametric) |
| Lipschitz axiom (Leibniz) | quasi-Leibniz | quasi-Leibniz (relaxed) | inherited | Leibniz hermitian norm |
| Examples worked | Heisenberg over quantum 2-tori | extends modular examples | metric spectral triples on quantum tori, Sierpinski gasket | $c_0(\mathbb{Z}) \rtimes \mathbb{Z}$ + classical $C_0(X)$ |

**For GeoVac applications:**
- M-X (Sturmian-FCI): non-compact algebra ($C_0(\R^3)$ or relatives) — requires pointed-QLCMS framework or non-compact extension. The modular line as stated doesn't directly apply, but the spectral propinquity framework via metrical quantum vector bundles does (the framework is compact for the Sturmian-truncated algebras).
- M-Y (NaH pin state): the framework needs a pin state to distinguish candidate Slater determinants. arXiv:2512.03573 provides this. The dual modular propinquity provides completeness for the metric on pin-state-distinguished MQVBs.
- M-Z (Bethe log Lamb shift): bound-state QED operates on the metric spectral triple regime — arXiv:1811.10843 directly applies.
- M-H1 (Higgs as inner fluctuation): inner fluctuations $\omega = a \, dD \, b$ are bimodule endomorphisms on the spectral triple's Hilbert module, so the modular / dual modular framework is directly relevant.

---

## §5. Connection to GeoVac vocabulary

| GeoVac structure | Latrémolière modular / dual-modular structure | Anchor paper |
|:-----------------|:--------------------------------------------|:------------:|
| Camporesi-Higuchi spinor bundle on $S^3$ | left Hilbert module over $C(S^3)$ (or its truncation) | 1608.04881 Def 3.1 |
| Truncated CH spectral triple $T_{n_{\max}}$ (Paper 38) | metrical quantum vector bundle (Hilbert module + algebra action + D-norm derived from CH Dirac) | 1811.04534 Def 2.12 |
| Sturmian basis at exponent $\lambda = Z/n_*$ | basis of the Hilbert module $\mathcal{M}$ at a specific D-norm choice (focal length encoded in D-norm) | 1608.04881 + 2512.03573 |
| Berezin reconstruction (Paper 38 L4) | candidate modular bridge / tunnel morphism between truncated and full MQVB | 1608.04881 Def 4.4 + 1811.04534 Def 2.20 |
| Spinor Dirac $D_{CH}$ on $S^3$ | first-order operator on module; gives D-norm $\mathsf{D}_{\mathcal{M}}(\omega) = \|\omega\|_{\mathcal{M}} + \|D_{CH} \omega\|_{\mathcal{M}}$ | 1608.04881 §1 intro |
| Schrödinger Coulomb $-\frac{1}{2}\nabla^2 - Z/r$ | Bochner Laplacian of metric connection $\nabla$ on trivial line bundle over $\R^3$; D-norm $\mathsf{D}_{\mathcal{M}}(f) = \|f\|_{L^2} + \|\nabla f\|_{L^2}$ | 1608.04881 §1 intro |
| Inner fluctuation $\omega = a \, dD \, b$ (Sprint H1 Higgs candidate) | bimodule endomorphism on metric spectral triple; modulates D-norm | 1811.10843 + 1608.04881 |
| Paper 34 §III.20 Phillips-Kleinman projection | candidate modular morphism / tunnel construction restricting MQVB to valence sector | 1811.04534 Def 2.9 |
| Roothaan multipole termination | sparsity property of D-norm at $L_{\max} = 2 \ell_{\max}$ — NOT a propinquity-vocabulary object | (no direct mapping) |
| Coulomb $S^3$ vs Bargmann $S^5$ four-layer asymmetry | Hilbert modules over different base QCMS — the framework treats them as different MQVBs with no canonical comparison | (no direct mapping) |

**Substantive overlap:** the Camporesi-Higuchi spectral triple (Paper 38) IS a metric spectral triple in the sense of arXiv:1811.10843, hence a metrical quantum vector bundle. Paper 38's propinquity convergence theorem is exactly the (covariant) modular / spectral propinquity machinery applied to GeoVac's $S^3$ structure. The four-witness Wick-rotation theorem (Paper 42) operates on the same MQVB substrate. **Sprint H1's Higgs-as-inner-fluctuation question lives natively in this vocabulary** — inner fluctuations are exactly the bimodule endomorphisms the dual modular framework handles.

---

## §6. Roothaan-autopsy substantive overlap check (HONEST)

**The vocabulary similarity is at the framework level, not the decomposition-structure level.**

**Roothaan autopsy structure (Paper 34 §V.C):** an observable $\mathcal{O}$ (e.g., H 21cm hyperfine) is decomposed additively as
$$\mathcal{O} = \mathcal{O}_{\mathrm{BF}} + \mathcal{O}_{a_e} + \mathcal{O}_{\mathrm{red.mass}} + \mathcal{O}_{\mathrm{Zemach}} + \cdots$$
where each component sits at a specific Paper 34 projection / focal length and the decomposition is additive across components. The decomposition exposes literature convention mismatches (§V.D) and identifies which components are framework-native vs Layer-2 inputs.

**Dual modular propinquity structure (Latrémolière 1811.04534):** a metric $\Lambda^{*\mathrm{mod}}$ on the class of MQVBs is defined via modular tunnels (instead of modular bridges). The "dual" refers to the **bridge ↔ tunnel duality** — bridges embed INTO a containing algebra, tunnels project FROM a common algebra. There is no additive decomposition structure; the metric is defined via an infimum over tunnels, not via a sum over components.

**Surface vs substantive overlap analysis:**

| Question | Roothaan autopsy | Dual modular propinquity |
|:---------|:----------------|:------------------------|
| What is "dualized"? | Nothing — autopsy is decomposition, not duality | Bridges → tunnels (embedding ↔ surjection duality) |
| Additive structure? | Yes — sum of components | No — infimum over tunnels |
| Component-by-component meaning? | Each component is a Paper 34 projection | Components are not introduced; tunnels are single objects |
| Bimodule structure? | No — observables are scalars at a focal length | Implicitly — modular morphisms on Hilbert modules |
| Convention-mismatch detection? | Yes (§V.D living catalogue) | No — propinquity is convention-free |
| Output type | Decomposed observable + residual structure | A single distance value |

**Verdict: surface-level vocabulary similarity only.** The word "dual" refers to:
- In Roothaan autopsy context (not actually invoked in Paper 34): would refer to dual observables (transition amplitudes vs energy components) — but Paper 34 does NOT use "dual" as a structural term.
- In dual modular propinquity context: refers to the dual notions of isometric embedding (tunnels project FROM a common space, bridges embed INTO a common space — the two are dual constructions).

The PI flagged this because the SUM-OF-PIECES quality of both objects looks superficially similar. **But the math.OA "dual" is a tunnel-vs-bridge duality of CONSTRUCTION METHODS, while the GeoVac "autopsy" is an ADDITIVE DECOMPOSITION of a single observable across projection-tagged components.** These are structurally distinct: a tunnel-vs-bridge construction does NOT produce additive decompositions; an additive decomposition does NOT live in a metric-construction-duality framework.

**Where the substantive overlap actually IS:** at the framework level. Both Latrémolière 2512.03573 and GeoVac's Roothaan-autopsy methodology operate on Hilbert-module pin-state pairings (sub-sprint Y's reading of the W1c chemistry pin-state). The pin state is the structural component shared between the two languages. But the decomposition / propinquity layer above the pin state is structurally different.

**Recommendation:** do not claim the Roothaan-autopsy structure has a math.OA equivalent in the dual modular propinquity vocabulary. The right framing is: the Roothaan autopsy operates at the Layer-2 itemization level (Paper 34 §V.C-D); the propinquity / modular / dual modular framework operates at the Layer-1 / spectral-triple substrate level. They are complementary, not equivalent.

---

## §7. Net assessment — paper recommendations per application track

### Track M-X (Sturmian-FCI Latrémolière error model)

**Anchor papers:**
- Latrémolière arXiv:2512.03573 (already deep-read) — pinned-QLCMS framework, non-unital extension, pin-state structure
- Latrémolière arXiv:1608.04881 — modular propinquity & D-norm vocabulary for the module-level interpretation (Hilbert spinor bundle as left Hilbert module over $C(S^3)$ truncation)
- Latrémolière arXiv:1811.04534 — relaxed MQVB definition (Def 2.8) is the right substrate if the Sturmian basis is interpreted at the module level

**Sub-question to answer in M-X:** is the Sturmian truncation at $n_{\max}$ a sequence of MQVB substructures (module-level), an L-Lipschitz μ-pinned exhaustive sequence (algebra-level, per 2512.03573), or both? The X/Z framing divergence from `sprint_l3e_p3_synthesis_memo.md` §1 has not yet committed to a single answer; M-X should resolve this with the module-level vocabulary on the table.

### Track M-Y (NaH pin-state W1c diagnostic)

**Anchor papers:**
- Latrémolière arXiv:2512.03573 — pin-state structure (Def 1.26) is the load-bearing component
- Latrémolière arXiv:1811.04534 — dual modular propinquity for COMPLETENESS — the W1c diagnostic compares candidate pin states (i) hydrogenic, (ii) SV-corrected, (iii) bonding-orbital, (iv) MP2-perturbed; the completeness of $\Lambda^{*\mathrm{mod}}$ guarantees Cauchy sequences of candidate pin states have a limit pin state (the "right answer" of the diagnostic)

**Sub-question to answer in M-Y:** compute the explicit $\Lambda^{*\mathrm{mod}}$-distance between candidate pin states; identify which candidate gives the smallest distance to the experimental NaH binding curve.

### Track M-Z (Bethe log Lamb shift Latrémolière error model)

**Anchor papers:**
- Latrémolière arXiv:1811.10843 — spectral propinquity for metric spectral triples (the LS-3 acceleration form on Dirac-S³ is a metric-spectral-triple computation)
- Latrémolière arXiv:2112.11000 — spectrum continuity in the spectral propinquity (the Lamb shift IS a Dirac-spectrum derived observable)
- Latrémolière arXiv:2512.03573 — pinned-QLCMS framework for the LS-7/LS-8a residual decomposition (non-loop physics envelope as operator-system-content beyond the pin-state local metametric)

**Sub-question to answer in M-Z:** identify the explicit spectral propinquity bound for the Sturmian truncation in LS-3 at $N = 12, 16, 20$. Compare to the empirical $\sim 1/N^4$ rate.

### Track M-H1 (Higgs as inner fluctuation in POSITIVE-THIN verdict)

**Anchor papers:**
- Latrémolière arXiv:1608.04881 — modular propinquity for MQVBs (the AC factor $\mathcal{A}_{GV} \otimes \mathcal{A}_F$ acts on a Hilbert module via inner fluctuations)
- Latrémolière arXiv:1811.04534 — dual modular propinquity for relaxed MQVB definition (the Higgs candidate inner fluctuation $\omega_{\mathrm{Higgs}}$ may not satisfy modular quasi-Leibniz, so the relaxed framework is needed)
- Latrémolière arXiv:1811.10843 — spectral propinquity bookkeeping (the H1 falsifier — "every Hermitian $D_F$ derived from GeoVac structure produces gauge-1-forms only" — is a spectral-propinquity-level statement)

**Sub-question to answer in M-H1:** formalize the H1 falsifier in the modular propinquity vocabulary. Does the GeoVac POSITIVE-THIN verdict (admits Higgs structurally but Yukawa is empirical input) translate to a propinquity-level statement about which inner fluctuations are reachable by modular morphisms?

### Cross-track recommendation

All four tracks should adopt a **shared notation convention** from the Latrémolière literature:
- $\mathsf{L}$ for the Lipschitz seminorm on the algebra (Latrémolière convention; matches GeoVac's existing $L = \|[D, \cdot]\|$ usage);
- $\mathsf{D}_{\mathcal{M}}$ for the D-norm on the module;
- $\Lambda^{\mathrm{mod}}$ for the modular propinquity, $\Lambda^{*\mathrm{mod}}$ for the dual modular propinquity;
- $\mathrm{mk}_{\mathsf{L}}$ for the Monge-Kantorovich distance on the state space.

This avoids notation clashes with GeoVac's existing math.OA papers (38, 39, 40, 42, 43, 44, 45, 46, 47) which use $\Lambda$ for the (single, non-modular) propinquity.

---

## §8. Honest scope and confidence

**This audit:**
- IS based on direct PDF extraction of arXiv:1608.04881 (pages 1-10, 23-40) and arXiv:1811.04534 (pages 1-12) — the canonical references
- DOES include direct paraphrase of Definitions 2.4, 3.1, 3.4, 3.5, 3.8, 4.4, 4.13–4.18, 5.6, 5.10–5.16 from arXiv:1608.04881 and Definitions 2.4, 2.7, 2.8, 2.9, 2.11, 2.12, 2.14, 2.15, 2.20, 2.22 + Theorem-Definitions 2.16, 2.19, 2.21 from arXiv:1811.04534
- DOES include web-search summaries of arXiv:1811.10843 (spectral propinquity) and arXiv:2112.11000 (spectrum continuity), but did NOT extract full content from these (web-summary level only)
- DOES NOT include full proof-level extraction (the §6 coincidence theorem of arXiv:1608.04881 is multi-page and would require a full read)

**Confidence:**
- HIGH on Definitions 2.4, 3.8, 4.4 of arXiv:1608.04881 (direct extraction)
- HIGH on Definitions 2.8, 2.12, 2.20 of arXiv:1811.04534 (direct extraction)
- HIGH on the bridge-vs-tunnel duality reading (extracted from both papers)
- HIGH on the Roothaan-autopsy substantive-vs-surface overlap analysis (§6) — clean structural comparison, not in dispute
- MEDIUM on the spectral propinquity content (arXiv:1811.10843) — summary level
- MEDIUM on the specific recommended-paper-per-track mapping (§7) — based on framework matching, not on having read every paper end-to-end

**What this audit does NOT do:**
- It does not produce the application tracks' actual computations — those are M-X, M-Y, M-Z, M-H1's tasks
- It does not adjudicate the X/Z framing divergence (algebra vs spectral choice of $\mathcal{A}$) — this is M-X's task
- It does not propose a Lorentzian extension — this remains the Krein-lift program from Phase A.2'
- It does not claim a substantive Roothaan-autopsy ↔ dual modular propinquity equivalence — the honest answer is that they are structurally distinct (§6 verdict)

---

## §9. Files and sources

**Direct PDF extractions (sources of primary definitions):**
- `webfetch-1779557971069-pktqnz.pdf` — Latrémolière arXiv:1608.04881 v5 (Feb 2018), 64 pages, 710 KB. Pages 1-10 and 23-40 extracted.
- `webfetch-1779557974582-90d7wk.pdf` — Latrémolière arXiv:1811.04534 v3 (Jul 2020), 52 pages, 553 KB. Pages 1-12 extracted.

**Web-search summaries (secondary):**
- arXiv:1803.06601 (Heisenberg modules application)
- arXiv:1811.10843 (spectral propinquity)
- arXiv:2112.11000 (spectrum continuity)
- arXiv:2211.11107 (Fell topology connection)
- arXiv:2504.11715 (Riemannian metric path continuity)

**Cross-references:**
- `debug/sprint_l3e_p3_synthesis_memo.md` — three-sub-sprint X+Y+Z synthesis
- `debug/l3e_p3_phase_a1_literature_audit.md` — May-23 Phase A.1 literature audit (Latrémolière 2512.03573 finding)
- `debug/l3e_p3_phase_a2prime_latremoliere_deep_read.md` — Phase A.2' deep-read of 2512.03573

**Output:** `debug/sprint_modular_propinquity_literature_audit.md` (this file, ~4200 words, 9 sections).
