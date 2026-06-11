# Sprint G6 scoping — graviton from CC-style metric variation on the GeoVac substrate

**Date:** 2026-05-28
**Type:** Scoping memo (no code, no production change). Identifies structural obstacles, candidate paths, minimum-cost first test.
**Verdict:** G6 in full is multi-month with a real structural obstacle (gravitons require explicit gamma-matrix / tetrad structure that GeoVac's discrete CH Dirac doesn't manifest preserve). The minimum-cost diagnostic — linearized perturbation theory on a truncated substrate — is 2-3 weeks and would distinguish "GeoVac hosts gravitons natively" from "GeoVac is structurally non-graviton" with a clean named falsifier. **Recommendation: do G6-Diag (2-3 weeks) before committing to full G6 (multi-month).**

## 1. Goal

Derive linearized graviton dynamics on the GeoVac discrete substrate via Connes-Chamseddine spectral-action variation. If successful, this would extend the Sprint G1/G2 Einstein-Hilbert + cosmological constant structure (which lives in the spinor sector at the level of the BACKGROUND action) to PERTURBATIONS — the spin-2 modes that propagate as gravitational waves.

The CC continuum approach:
1. Background: spectral triple $(\mathcal{A}, \mathcal{H}, D_0)$ with $D_0 = \gamma^\mu D_\mu$ on a Riemannian manifold
2. Perturb metric: $g \to g + h$ with $h$ TT symmetric rank-2
3. Tetrad/vierbein varies: $e^\mu_a \to e^\mu_a + \delta e^\mu_a$
4. Dirac varies: $\delta D = i\,\delta e^\mu_a \gamma^a \partial_\mu + ...$
5. Spectral action varies: $S[D + \delta D] = S[D] + \delta S + \delta^2 S + ...$
6. Second-order $\delta^2 S$ at the extremum gives the Fierz-Pauli graviton kinetic action

The structural-skeleton question for GeoVac: do steps 2-6 transfer to the discrete substrate?

## 2. The structural obstacle

GeoVac's substrate is the discrete Camporesi-Higuchi spectrum on $S^3$: eigenvalues $|\lambda_n| = (n+3/2)/R$, multiplicities $g_n = 2(n+1)(n+2)$, Hilbert space $\mathcal{H}_{n_{\max}}$ truncated at level $n_{\max}$.

This substrate preserves the **spectral content** of the continuum CH Dirac (eigenvalues + multiplicities + spinor labels $(n, \kappa, m_j)$). But it does NOT manifestly preserve the **explicit Lorentz / gamma-matrix structure** $D = \gamma^\mu D_\mu$ that CC's continuum gravity requires.

Consequence: the metric perturbation operation $\delta g \to \delta D$ doesn't have a clean discrete analog. The standard CC graviton construction relies on the tetrad $e^\mu_a$ in $\gamma^\mu = e^\mu_a \gamma^a$ as the natural "carrier" of spin-2 metric perturbations. GeoVac's discrete substrate has implicit (not manifest) gamma structure, so this carrier isn't directly available.

**Three candidate ways around this obstacle:**

## 3. Three paths

### Path P1: Explicit gamma-matrix re-derivation

Re-derive the CH Dirac from explicit gamma matrices on $S^3$ at finite cutoff, keeping the full Lorentz / tetrad structure as manifest data on top of the discrete spectrum. Then perturb the metric naturally via $\delta e^\mu_a$.

**Effort:** 2-4 months. Substantial new code: angular-momentum-projected gamma matrices, spin connection, tetrad fluctuations on the discrete substrate.

**Structural risk:** The discretization may not preserve all symmetries of the continuum gamma structure. There's no guarantee that the resulting "discrete gamma" data captures the same physics as the eigenvalue-only CH spectrum that GeoVac actually uses.

**Cleanness:** Most direct CC analog. If it works, the result is "GeoVac graviton = CC graviton on discrete substrate."

### Path P2: Quadratic / bilinear extension

Define gravitons not as linear perturbations of $D$ but as quadratic forms (bilinears of spinor matrix elements) on $\mathcal{H}_{n_{\max}}$. The spin-2 structure emerges from spinor $\otimes$ spinor decomposition (spin-1/2 $\otimes$ spin-1/2 = spin-0 $\oplus$ spin-1, but products of spin-1 give spin-2).

**Effort:** 3-6 months conceptual + technical. Requires extending CC spectral action to "second-quantized" perturbations.

**Structural risk:** This is a NEW conceptual extension of spectral action, not standard CC. Higher risk of non-uniqueness (multiple inequivalent extensions). May not yield Fierz-Pauli graviton.

**Cleanness:** Conceptually novel. Could be its own paper if it works.

### Path P3: External / hybrid continuum graviton

Treat $h_{\mu\nu}$ as a continuum field on $S^3$ that couples to the discrete substrate via the Fock projection. The graviton is "external" to the discrete framework; the discrete substrate acts as a "detector" or "matter" coupled to the gravitational field.

**Effort:** 1-3 months. Requires derivation of coupling between continuum $h_{\mu\nu}$ and discrete CH spectrum (the "spectral coupling" to gravity).

**Structural risk:** Not natively discrete. Loses the structural skeleton of GeoVac. Result would be "GeoVac substrate as matter in CC gravity background" rather than "GeoVac substrate hosts gravity."

**Cleanness:** Standard, well-defined. But less ambitious — doesn't tell us about discreteness of gravity.

## 4. Minimum-cost diagnostic test (G6-Diag)

Before committing to any of Paths P1-P3 (all multi-month), the load-bearing question is: **does the GeoVac substrate even have eigenmodes with spin-2 angular momentum content?**

This question can be answered by a sprint-scale-medium (2-3 weeks) diagnostic computation.

### 4.1 The math

Background: $D_0$ = CH Dirac on truncated $\mathcal{H}_{n_{\max}}$. Perturb $D = D_0 + \epsilon V$ where $V \in \text{Herm}(\mathcal{H}_{n_{\max}})$.

Spectral action expansion to second order:
$$S[D_0 + \epsilon V] = S[D_0] + \epsilon S^{(1)}[V] + \tfrac{1}{2}\epsilon^2 S^{(2)}[V, V] + O(\epsilon^3)$$

For Gaussian cutoff $f(x) = e^{-x}$, the quadratic form is:
$$S^{(2)}[V, V] = \tfrac{1}{2\Lambda^4}\,\mathrm{Tr}\bigl[e^{-D_0^2/\Lambda^2}\{D_0, V\}^2\bigr] - \tfrac{1}{\Lambda^2}\,\mathrm{Tr}\bigl[e^{-D_0^2/\Lambda^2} V^2\bigr] + (\text{Duhamel corrections})$$

This is a quadratic form on the real vector space $\text{Herm}(\mathcal{H}_{n_{\max}})$, which has dimension $(\dim \mathcal{H})^2$.

Diagonalize: find eigenmodes $V_\alpha$ with eigenvalues $\kappa_\alpha$. These are the "propagating modes" of perturbation theory at this background.

### 4.2 Angular momentum classification

The CH Hilbert space at $n_{\max}$ decomposes into $SO(4) = SU(2)_L \times SU(2)_R$ irreps. At level $n$, spinor harmonics transform as $((n+1)/2, n/2) \oplus (n/2, (n+1)/2)$ (the two chiralities).

Hermitian operators $V$ on $\mathcal{H}$ decompose into $SO(4)$ irreps via tensor products of these representations. Specifically:
- $(0, 0)$ scalar singlet
- $(1/2, 1/2)$ vector
- $(1, 0) \oplus (0, 1)$ chiral spin-1 (gauge bosons)
- $(1, 1)$ spin-2 (graviton)

**The diagnostic**: for each eigenmode $V_\alpha$ of $S^{(2)}$, decompose into $SO(4)$ irreps. Check for nonzero $(1, 1)$ content.

### 4.3 Concrete first test ($n_{\max} = 2$)

- $\dim \mathcal{H}_{n_{\max}=2} = 4 + 12 + 24 = 40$
- $\dim \text{Herm}(\mathcal{H}) = 40^2 = 1600$ (real dimension)
- Quadratic form $K$ = $1600 \times 1600$ real symmetric matrix
- Diagonalize: standard linear algebra (1-2 seconds)
- $SO(4)$ decomposition: known formulas for product irreps, classify each eigenmode

Fully tractable. Larger $n_{\max}$ available if needed for convergence checks.

### 4.4 Named falsifiers

**POSITIVE (G6 viable):** At least one eigenmode $V_\alpha$ of $S^{(2)}$ has nonzero $(j_L, j_R) = (1, 1)$ content under $SO(4)$ decomposition, with nonzero eigenvalue $\kappa_\alpha$ (nondegenerate, propagating mode). This is a GeoVac graviton candidate. Full G6 (multi-month, Path P1 most likely) is justified.

**NEGATIVE (G6 structurally blocked):** All eigenmodes $V_\alpha$ carry only $(j_L, j_R)$ irreps with $j_L + j_R \leq 1$ (i.e., scalar, spin-1/2, or spin-1 content only). The GeoVac substrate doesn't natively support spin-2 perturbations. Pivot to Path P3 (hybrid) or accept that gravity is structurally external to GeoVac.

**AMBIGUOUS:** $(1, 1)$ content is present but only in "spurious" modes (zero eigenvalue, gauge-like, etc.). Requires more careful analysis — extend to $n_{\max} = 3, 4$ for convergence; check whether the modes propagate or are pure gauge.

### 4.5 Effort estimate for G6-Diag

| Component | Effort |
|---|---|
| Implementation of $S^{(2)}$ as quadratic form | 1 week |
| $SO(4)$ angular momentum decomposition | 3-5 days |
| Diagonalization + mode classification | 3-5 days |
| Result interpretation + memo | 1 week |
| **Total** | **2-3 weeks** |

This is sprint-scale-medium. Comparable to typical week-long sprints.

## 5. The deep question (worth thinking about)

Even if G6-Diag finds $(1, 1)$ eigenmodes, there's a deeper question: **are these modes physical gravitons, or just "spin-2 bilinears of spinor data"?**

In CC's continuum, the graviton arises BECAUSE the Dirac operator has an explicit metric dependence via the tetrad. The metric IS the physical input; the graviton IS the spin-2 fluctuation of that input. Spectral action is a functional of the metric, and varying it gives gravity equations.

In GeoVac's discrete substrate, there's no "physical metric" input that the substrate varies over. The CH spectrum is FIXED (per Fock rigidity, Paper 23). Perturbations $V$ are just Hermitian matrices — they don't have an a priori interpretation as "metric perturbations."

So even if $(1, 1)$ eigenmodes exist, calling them "gravitons" requires an additional argument: that these modes couple to matter (or other modes) the way gravitons should, that their kinetic action is Fierz-Pauli, etc. This is content beyond just the existence of spin-2 perturbations.

**The cleanest reading:** G6-Diag's $(1, 1)$ content (if positive) would be a NECESSARY condition for graviton dynamics, not sufficient. Sufficient conditions require building out Path P1 or P2.

**The cleanest negative reading:** Absence of $(1, 1)$ content (if negative) is sufficient to rule out gravitons at the substrate level. Forces Path P3 (hybrid) or "GeoVac is non-gravitational."

## 6. Alternative forward paths in the gravity arc

Instead of (or alongside) G6-Diag, two other multi-month paths are well-defined:

### G4: Cigar geometry parameterized by black-hole mass M

Compute the spectral action on the Euclidean Schwarzschild-like cigar (Sprint TD Track 4 already did $T_H = 1/(8\pi M)$). Chase Bekenstein-Hawking $S = A/4$ as the concrete target.

**Effort:** Multi-month. Tractable: well-defined within standard CC, clear target, doesn't require new conceptual extensions.

**Cleanness:** Most concrete gravity-side target. Result either confirms BH entropy from spectral action on discrete substrate, or identifies where the construction fails.

### G5: Decompactified $S^3 \times \mathbb{R}_\tau$

Extend G2 from compact $S^1_\beta$ to non-compact $\mathbb{R}_\tau$ to make the $\beta \to \infty$ minimum (zero-temperature de Sitter vacuum) explicit. Sprint-scale, $\sim 1$ week.

**Effort:** Sprint-scale (1-2 weeks).

**Cleanness:** Mechanical extension of G2, doesn't address gravity directly but cleans up the gravity-side interpretation from G2.

## 7. Recommendation

Three options, in order of effort:

1. **G5 first (1-2 weeks).** Sprint-scale cleanup of G2's $\beta \to \infty$ minimum. Doesn't touch G6.
2. **G6-Diag (2-3 weeks).** The load-bearing diagnostic. Either:
   - POSITIVE result $\to$ commit to G6 full (multi-month, Path P1)
   - NEGATIVE result $\to$ pivot to G4 or accept GeoVac is structurally non-graviton
3. **G4 (multi-month).** Bekenstein-Hawking from cigar spectral action. Concrete gravity-side target. Doesn't require G6.

**My recommendation: G6-Diag.** It's the cheapest test of the structural question "does GeoVac substrate host gravitons natively?" Either result is informative. The cost (2-3 weeks) is small compared to the multi-month commitments of G4 or full G6, and the answer materially affects which of those paths is worth pursuing.

If POSITIVE: G6 full becomes the natural next sprint, with the diagnostic providing strong evidence the construction will succeed.

If NEGATIVE: G6 is structurally closed. The gravity arc's future is then G4 (Bekenstein-Hawking) as a concrete spectral-action gravity result without graviton dynamics. The structural-skeleton-scope reading of GeoVac sharpens: framework reproduces Einstein-Hilbert + cosmological constant at the BACKGROUND level (G1, G2 on Dirac sector) but does not host gravitational waves.

## 8. Files

- `debug/g6_scoping_memo.md` — this scoping memo (no code, no production change)

No CLAUDE.md / Paper / CHANGELOG updates from this scoping pass. Updates would come either from running G6-Diag (next sprint candidate) or from PI decision to defer.

## 9. Cross-references

- **Sprint G1** (`sprint_g1_path1_spectral_action_S3_radius.md`) — background action's two-term exactness on Dirac sector
- **Sprint G2** (`sprint_g2_spectral_action_S3_x_S1.md`) — propagation to 4D thermal product
- **Sprint G3** (`sprint_g3_scalar_TT_S3.md`) — spinor-bundle specificity of two-term exactness; structural distinction between Dirac and tensor sectors
- **Paper 23 Fock rigidity theorem** — the $S^3$ metric is rigid under the framework's Coulomb projection; structural foundation of the "Path P3 hybrid only" obstacle
- **memory/geovac_structural_skeleton_scope_pattern** — gravity coupling and cosmological constant are Class 1 calibration data; G6-Diag negative result would extend this to "graviton dynamics is also outside framework scope"
- **CC spectral action literature** (Chamseddine-Connes 1997 onward, Connes-Marcolli 2008) — standard continuum derivation that the discrete framework attempts to mirror
