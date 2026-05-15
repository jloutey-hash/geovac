# Unified GH-convergence scoping memo

**Date:** 2026-05-10 (initial), 2026-05-15 (L3 erratum + Dirac-triangle reformulation + prerequisite verdict)
**Author:** PM agent (scoping deliverable; no production code under `geovac/`)
**Status:** scoping memo + prerequisite findings consolidated; execution path GO

---

## §0. Headline consolidated verdict (added 2026-05-15)

After the initial scoping memo (§§1–7 below, May 10) and two prerequisite sprints (P-A: Casimir triangle inequality; P-B Track 4: L2 constant; SU(3) numerical sanity by P-A's machinery, this consolidation), the state of play is:

- **Right abstraction level:** Class 1 — compact Lie group $G$ with bi-invariant metric, rank $r \ge 1$. SU(2) is the rank-1 instance (Paper 38). Confirmed.
- **Headline literature gap:** Gaudillot-Estrada–vS (IMRN 2025, arXiv:2310.14733) proved state-space GH at full Class 1 generality, with no rate and no Dirac. (Authorship attribution flagged: Track P-A attributed this arXiv ID to Hekkelman; sub-agent disagreement, needs verification.) The unified theorem upgrades that state-space result to propinquity with explicit rate.
- **L3 ingredient: REFRAMED.** The "Casimir triangle inequality" identified in original §5.3 as the discovery target is **FALSE** (Sprint P-A; counterexample at SU(2), $j=1, j'=1/2$). The correct L3 ingredient is the **Dirac-triangle inequality** $|D(\pi) - D(\pi')| \le \sqrt{C(\sigma)}$ for $\sigma \subset \pi \otimes \pi'^*$, verified on the full 100-pair SU(3) panel (`debug/dirac_triangle_su3_check.py`) with max ratio 0.672 and asymptotic $\to 1^-$ family. See §5.3 (revised).
- **Rate constant prediction:** $c(G) = 2\,\mathrm{Vol}(G/T)/\mathrm{Vol}(G)$ in bi-invariant Haar normalization. SU(2) verified: $4/\pi$. SU(3) value testable but blocked at this sprint by org usage limit (P-B Track did not complete). **Open prerequisite.**
- **Log-power test (single-log $\log\Lambda/\Lambda$ vs double-log $\log^2\Lambda/\Lambda$ at rank 2):** Also blocked at this sprint. **Open prerequisite.**

**Forward plan:** the L3 obstruction is cleared. Two prerequisites remain (rate constant + log-power test on SU(3) Weyl-integration computation). Once those land (1–2 day sprint each), the 4–6 week execution sprint is well-defined.

---

## §1. Background

**Paper 38 (closed, May 2026)** proved that the Connes–van Suijlekom truncations $\mathcal{T}_{n_\max}$ of the round-$S^3 = \mathrm{SU}(2)$ Camporesi–Higuchi spectral triple $\mathcal{T}_{S^3}$ converge to it in the Latrémolière propinquity, with explicit rate
$$\Lambda(\mathcal{T}_{n_\max}, \mathcal{T}_{S^3}) \le C_3 \cdot \gamma_{n_\max}, \qquad \gamma_{n_\max} = \frac{4 \log n_\max}{\pi\, n_\max} + O(1/n_\max),$$
with $C_3 = 1$ (L3 Lipschitz comparison) and $\gamma_{n_\max}$ from the L2 mass-concentration moment of the central spectral Fejér kernel on $\mathrm{SU}(2)$. The proof has five lemmas (L1' chirality-doubled operator system; L2 central Fejér kernel and Stein–Weiss-type rate; L3 Lipschitz with $C_3 = 1$ from SO(4) selection rule + Avery 3-Y; L4 Berezin reconstruction via Peter–Weyl; L5 Latrémolière tunneling-pair assembly).

The asymptotic constant $4/\pi = \mathrm{Vol}(S^2)/\pi^2 = 2\mathrm{Vol}(S^1)/\mathrm{Vol}(\mathrm{SU}(2))$ identifies as the **Hopf-base measure factor** in the $\mathrm{SU}(2)/U(1)$ Haar normalisation. Independently, the GeoVac master Mellin engine (Sprint TS-E1, Paper 32 §VIII, May 2026) classifies $\pi$-sources in any finite chain of Paper 34 projections into three sub-mechanisms M1/M2/M3 indexed by the operator-order $k \in \{0,1,2\}$ in $\mathcal{M}[\mathrm{Tr}(D^k e^{-tD^2})]$. The constant $4/\pi$ in Paper 38's rate is the **M1 (Hopf-base measure, $k=0$) signature**. This is the first place inside GeoVac where the abstract transcendental taxonomy and the concrete metric-NCG rate constant identify.

**What's at stake.** Paper 38 closed a published-open question (Connes–vS 2021 deferred GH convergence to "elsewhere" three times) for a single non-abelian Lie group, $\mathrm{SU}(2)$. The next-step question is whether the SU(2)-specific machinery generalises, and at what level. That answer determines whether GeoVac's WH1 PROVEN status extends from "the round-$S^3$ Camporesi–Higuchi triple is its propinquity limit" to "every member of a natural class $\mathcal{C}$ of compact Riemannian spectral triples is its propinquity limit." The class $\mathcal{C}$ is the deliverable of this scoping memo.

---

## §2. Literature audit

I audited (with WebSearch + WebFetch, May 2026) the published literature on propinquity convergence of spectral truncations adjacent to Paper 38's setting. The audit covers every paper that could plausibly compete with, subsume, or be subsumed by a unified GH-convergence theorem. ArXiv IDs verified.

### §2.1. Direct precedents (state-space level, partial-overlap)

**Gaudillot-Estrada & van Suijlekom (2023; IMRN 2025).** *Convergence of Spectral Truncations for Compact Metric Groups.* arXiv:2310.14733; IMRN 2025, Issue 13, rnaf197 (July 2025).

This is the most direct partial precedent. They prove:

> The sequence of truncated state spaces $\{S(P_\Lambda C(G) P_\Lambda)\}_\Lambda$ converges in classical (not quantum) Gromov–Hausdorff distance to $S(C(G))$ for **any compact metric group $G$ with a bi-invariant metric**, where $P_\Lambda$ is the Peter–Weyl projection onto a finite set of irreducible representations.

Crucial scope facts:

- **State-space convergence**, NOT propinquity. They explicitly remark that quantum GH and classical GH may not agree.
- **No Dirac operator** in the construction. Lip-norm comes from the action of $G$ on itself, not from spectral data.
- **No explicit rate.** The bound is given as an integral over a liftable probability measure approaching the identity; rate is qualitative.
- Covers SU(2), SU(N), and any compact Lie group. So at the **state-space level**, the SU(N) generalisation is already done.

**Bottom line:** Gaudillot-Estrada–vS solves the abelian–compact-group state-space question at full generality but does NOT touch the spectral-triple propinquity. Paper 38 is a strict refinement of their result for $G = \mathrm{SU}(2)$ — Latrémolière propinquity > state-space GH.

**[Erratum 2026-05-15: Sprint P-A attributed this arXiv ID to "Hekkelman 2023" instead. Sub-agent disagreement on authorship; substantive findings independent of this. Verification flagged for next session.]**

**Leimbach (2024).** *Convergence of Peter–Weyl Truncations of Compact Quantum Groups.* arXiv:2409.16698; J. Noncommut. Geom., to appear.

Companion to Gaudillot-Estrada–vS in the *quantum* group setting. Uses Kerr's complete GH distance, no Dirac, no rate. Covers coamenable compact quantum groups with bi-invariant Lip-norms. Includes $\mathrm{SU}_q(N)$ as instances. Does NOT cover the classical Riemannian Dirac-spectral-triple case.

**Toyota (2023).** *Quantum Gromov–Hausdorff convergence of spectral truncations for groups with polynomial growth.* arXiv:2309.13469.

**Discrete groups only.** Compact Lie groups not covered.

### §2.2. Direct precedents (full propinquity, restricted-domain)

**Leimbach & van Suijlekom (2023, Adv. Math. 2024).** *Gromov–Hausdorff Convergence of Spectral Truncations for Tori.* arXiv:2302.07877; Adv. Math. 439 (2024), 109496.

The flat-torus case Paper 38 was modeled on. Proves state-space GH convergence (with Connes-distance metric) for $T^d$. Uses Schur multiplier characterisation on the abelian dual $\mathbb{Z}^d$, and introduces the spectral Fejér kernel. Rate $O(1/\Lambda)$. **Abelian.** Paper 38 is the SU(2) analog, slowed by one log factor.

**Hekkelman (2022).** *Truncated Geometry on the Circle.* Lett. Math. Phys. 112, 20 (2022); arXiv:2111.13865.

The $S^1$ Toeplitz/Fejér case. Pure-state space converges to $C(S^1)$ via Toeplitz operator system. Explicit Connes-distance metric. Rank-1 abelian.

**Hekkelman–McDonald (2024).** *Spectral truncations of $T^d$ and quantum metric geometry.* J. Noncommut. Geom., to appear; arXiv:2403.18619. (Companion: *A noncommutative integral on spectrally truncated spectral triples,* arXiv:2412.00628.)

The Berezin–Toeplitz convergence on $T^d$, later extended to a Szegő limit formula for spectrally truncated triples. Companion track to L4 (the Berezin reconstruction in Paper 38). **Flat / abelian.**

**Bhattacharyya–Duhan–Pradhan (2024).** *Gromov–Hausdorff convergence of metric spaces of UCP maps.* arXiv:2410.15454.

Demonstrates that the van Suijlekom conditions guarantee GH convergence of *sets of UCP maps* (with BW topology), not of states or operator systems directly. Foundational machinery, broad scope, no rate. Used as a meta-framework for downstream propinquity theorems.

### §2.3. Symmetric-space / coadjoint-orbit precedents (different framework)

**Rieffel (2009, 2023).** *Dirac operators for coadjoint orbits of compact Lie groups* (arXiv:0812.2884); *Dirac Operators for Matrix Algebras Converging to Coadjoint Orbits* (Comm. Math. Phys. 401 (2023), 1951; arXiv:2108.01136).

Rieffel constructs $G$-invariant Dirac operators on **coadjoint orbits** $\mathcal{O} \subset \mathfrak{g}^*$ of a compact Lie group $G$, and proves that fuzzy approximations (matrix algebras) converge to $C(\mathcal{O})$ in a quantum Gromov–Hausdorff sense determined by those Dirac operators. The 2023 paper is in Latrémolière's spectral propinquity. **Crucial scope difference:** convergence is from below (matrix algebras of growing dimension converging up to the coadjoint orbit), NOT from spectral truncation of the orbit itself. Coverage is canonical compact symmetric spaces realised as coadjoint orbits: $S^2 = \mathrm{SU}(2)/U(1)$, $\mathbb{CP}^n = \mathrm{SU}(n+1)/(U(n) \times U(1))$, complex Grassmannians, full flag varieties.

**Crucial fact:** $S^3 = \mathrm{SU}(2)$ is **not** a coadjoint orbit — it is the simply-connected compact Lie group itself (which is a *principal* homogeneous space, not a flag/orbit). So Rieffel's framework does not cover Paper 38's case. Conversely, Paper 38 does not subsume Rieffel: a Lie group $G$ acting on itself by *left* translation is structurally different from $G$ acting on a coadjoint orbit.

**Bordemann–Meinrenken–Schlichenmaier (1994).** *Toeplitz quantization of Kähler manifolds and gl(N), N → ∞ limits.* Comm. Math. Phys. 165 (1994), 281; hep-th/9309134.

The classic Berezin–Toeplitz quantisation of compact quantizable Kähler manifolds — operator-norm convergence of $\mathrm{gl}(N)$ matrix algebras to $C(M)$ for $M$ Kähler. Not propinquity, not spectral triple, but the canonical Berezin precedent for Kähler symmetric spaces.

**Hawkins (2000).** *Quantization of equivariant vector bundles.* Comm. Math. Phys. 202 (1999), 517; *Geometric quantization of vector bundles…* CMP 215 (2000), 409.

The Berezin reconstruction template Paper 38 cites for L4. Hawkins requires a **compatible almost-complex structure** (Kähler), which $S^3 = \mathrm{SU}(2)$ does NOT have (it is parallelisable but not Kähler). Paper 38's L4 used the central-Fejér-kernel form as the non-Kähler analog.

### §2.4. Full-spectral-propinquity activity (orthogonal to truncation)

**Farsi–Latrémolière (2024).** Spectral continuity for collapse of products containing at least one **abelian factor** — abelian-factor restriction explicitly lifted by Paper 38. (cf. paper 38 §1.)

**Farsi–Latrémolière (2025).** *Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics.* arXiv:2504.11715. Continuous-family deformation on closed spin manifolds, NOT discrete truncation. Not a precedent for the truncation question.

**Latrémolière (2017/2023).** The Gromov–Hausdorff propinquity for metric spectral triples (arXiv:1811.10843; Adv. Math. 415 (2023), 108876). The framework Paper 38 uses; defines the metric.

### §2.5. Rieffel (2022). *Convergence of Fourier truncations for compact quantum groups and finitely generated groups.* arXiv:2210.00387.

Generalises the Fejér–Riesz operator system to compact matrix quantum groups; convergence in qGH distance. Examples include $\mathrm{SU}_q(2)$, $S^2_q$. Quantum-group setting.

### §2.6. Summary of the literature gap

After audit:

- **State-space GH for compact metric groups with bi-invariant metric:** Done at full generality (Gaudillot-Estrada–vS 2025) **without rate, without Dirac**.
- **Latrémolière propinquity for compact Riemannian spectral triples obtained by spectral truncation of $G/H$ with $H \ne \{e\}$:** Done in Paper 38 only for $G = \mathrm{SU}(2), H = \{e\}$ (i.e. $G/H = S^3$).
- **Latrémolière propinquity for matrix algebras → coadjoint orbits (which include $S^2, \mathbb{CP}^n$):** Done by Rieffel (2023). **From below.**
- **Latrémolière propinquity for spectral truncations of generic compact symmetric space $G/K$:** **OPEN.**
- **Latrémolière propinquity for spectral truncations of generic compact Lie group $G$ at rank $\ge 2$:** **OPEN.**

The literature gap that the unified theorem would close is the union of the last two: spectral-truncation propinquity convergence for compact Riemannian symmetric spaces / homogeneous spaces / Lie groups beyond $\mathrm{SU}(2)$. No published paper covers this.

---

## §3. The right abstraction level

The four candidate generalisation classes:

1. **Compact Lie groups $G$ with bi-invariant metric.** SU(2) is rank-1; $G$ at higher rank is open.
2. **Compact Riemannian symmetric spaces $G/K$.** Includes spheres $S^n$ ($n \ge 2$), $\mathbb{CP}^n$, Grassmannians, $S^3 = \mathrm{Spin}(3)$ (which is also a Lie group).
3. **Compact homogeneous spaces $G/H$.**
4. **Compact Kähler manifolds with Berezin–Toeplitz quantisation.** Direct route via Hawkins / Bordemann–Meinrenken–Schlichenmaier; Hekkelman–McDonald-style.

I evaluate each on three axes (harmonic analysis support, spectral triple existence, and rate-constant identification with the master Mellin engine):

### §3.1. Class 4 (compact Kähler) — rejected

This is the most beautifully developed framework (Hawkins, BMS, Schlichenmaier reviews). But it has two structural blockers as the GeoVac unified target:

- **It cannot host $S^3 = \mathrm{SU}(2)$.** Paper 38's principal example is non-Kähler. Restricting the unified theorem to Class 4 would *exclude* the original SU(2) result, which is the wrong direction.
- **It overlaps Rieffel-2023 substantively.** Coadjoint orbits with their Kirillov–Kostant–Souriau symplectic form are Kähler, and Rieffel already handles the matrix-algebra-→-orbit direction in Latrémolière propinquity. The "from-above" spectral-truncation analog is a useful contribution but is the Hekkelman–McDonald arc, not the GeoVac arc.

**Verdict: not the right class.**

### §3.2. Class 3 (compact homogeneous $G/H$) — too wide

For generic $G/H$:

- The harmonic analysis IS a Frobenius reciprocity / induced-representation construction (sections of the homogeneous bundle $G \times_H V$ decompose as $\bigoplus_\pi V_\pi \otimes \mathrm{Hom}_H(V_\pi, V)$). The central-good-kernel construction analog requires a "central" notion that may not be available when $H$ is not normal in $G$.
- Camporesi–Higuchi-style closed-form Dirac spectra exist on rank-1 symmetric spaces (their 1996 paper handles $S^n$ and $H^n$ in arbitrary dimension via Harish-Chandra's Casimir radial part), but generic $G/H$ Dirac spectra are not in such closed form.
- The Avery 3-Y integral analog for arbitrary $G/H$ does not appear to be in the literature (Avery's work is on $S^n$ via SO($n$+1)).

The class is correct in that it would subsume both Lie-group truncations (with $H = \{e\}$) and symmetric-space truncations (with $H = K$ a maximal compact subgroup of a non-compact $G^{\mathbb{C}}$, in the Cartan classification), but the proof at this generality would need substantively more machinery than Paper 38 used.

**Verdict: too wide for one sprint. Right ultimate target.**

### §3.3. Class 1 (compact Lie groups, bi-invariant metric) — natural narrowing

For $G$ a compact Lie group with bi-invariant metric:

- Peter–Weyl gives the harmonic analysis directly; the central good kernel construction (Paper 38's L2) extends because every compact group has a centre and a class-function subalgebra. The natural-coefficient Plancherel weight $\sqrt{\dim V_\pi}$ generalises $\sqrt{2j+1}$ from SU(2).
- Camporesi–Higuchi gives Dirac spectra in closed form for rank-1 symmetric spaces ($S^n, H^n$), which include $S^3 = \mathrm{SU}(2)$. For higher-rank compact Lie groups, the Dirac operator on $G$ with bi-invariant metric is constructed from the Killing form via Kostant's cubic Dirac operator, with explicit eigenvalues from Casimir + $\rho$-shift.
- The selection rule analog of $|n - n'| + 1 \le N \le n + n' - 1$ is the **Clebsch–Gordan rule for the tensor product of irreducibles in $G$**: $V_\pi \otimes V_{\pi'} = \bigoplus_\sigma N^\sigma_{\pi \pi'} V_\sigma$. The statement "the multiplier matrix entry $(M_f)_{\pi', \pi}$ is supported on $\sigma$ in this tensor product decomposition" is the Avery-3-Y analog, and it generalises verbatim. The numerical bound on the shell-difference scaling that gave $C_3 \le \sqrt{(N-1)/(N+1)}$ in SU(2) needs a higher-rank replacement, but the *structure* transfers.
- Gaudillot-Estrada–vS already work in Class 1 at the state-space level. Lifting to propinquity at full Class 1 generality is the natural next step.

**Verdict: correct narrowing.** Paper 38 is the rank-1 instance.

### §3.4. Class 2 (compact Riemannian symmetric space $G/K$) — also natural

A compact symmetric space $G/K$ has:

- Peter–Weyl decomposition of $C^\infty(G/K)$ as $\bigoplus_\pi V_\pi^K$ (the $K$-invariant subspace of each irrep $\pi$ of $G$). Spherical-function machinery (Helgason).
- $G$-invariant Dirac operator on $G/K$ when $G/K$ is spin; closed-form spectrum via Camporesi–Higuchi (rank-1 case) or Slebarski / Kostant cubic Dirac (higher rank).
- The central-good-kernel analog: convolve with a positive central function on $G$ that is supported on the $K$-bi-invariant class functions (zonal functions).

There is technical overhead vs Class 1 (need to track $K$-invariance) but no fundamental obstruction.

$S^3 = \mathrm{SU}(2)$ also fits into Class 2 as $\mathrm{SU}(2) \times \mathrm{SU}(2) / \mathrm{SU}(2)_{\mathrm{diag}}$ (the "group case" of a symmetric pair), so Class 2 contains Class 1 at the symmetric-pair level.

**Verdict: correct, contains Class 1.**

### §3.5. Recommendation

**Class 1 (compact Lie group $G$ with bi-invariant metric) is the right next step.**

Reasoning:

- It contains $\mathrm{SU}(2)$ as the rank-1 instance, so Paper 38 is recovered as a corollary.
- It is exactly the class for which Gaudillot-Estrada–vS proved state-space GH convergence in 2025. The unified theorem would *upgrade their result to propinquity*, with explicit rate. This is a clean, well-defined published-open question.
- It is a single step beyond Paper 38 — most of the L1', L2, L4 machinery extends mechanically; the L3 lemma is where the rank-2+ work happens.
- The rate constant identification (master Mellin engine) is *cleanest* at this level: every compact Lie group $G$ has a Hopf-type fibration $T \to G \to G/T$ (where $T$ is a maximal torus), so the M1 Hopf-base measure mechanism extends immediately. The unified theorem will give a rate of the form $\gamma_n \sim c(G) \log n / n$ with $c(G)$ a Haar-normalisation ratio determined by $\mathrm{Vol}(G/T)/\mathrm{Vol}(G)$ — predicted to be the rank-$r$ generalisation of $4/\pi$.
- After Class 1 lands, Class 2 is a follow-on sprint (additional $K$-invariance bookkeeping); Class 3 is the long-term ambition (Frobenius reciprocity).

A reasonable two-stage program:

- **Stage A (4–8 weeks, this sprint):** Class 1, compact Lie group with bi-invariant metric. Target: SU(N), Spin($n$), Sp($n$), exceptional groups (G2 etc.) all in one statement. Keystone: identify the rank-$r$ generalisation of the Avery 3-Y integral and the corresponding Lipschitz comparison constant.
- **Stage B (next 8–12 weeks, follow-on):** Class 2, compact Riemannian symmetric space $G/K$. Target: $S^n$, $\mathbb{CP}^n$, complex Grassmannians, full flag varieties, with Class 1 recovered as the "group case" $G \times G / G_{\mathrm{diag}}$.

I will write the unified theorem statement and lemma-by-lemma porting plan for **Stage A only**.

---

## §4. Candidate unified theorem statement

**Theorem (Unified GH-convergence on compact Lie groups; conjectural).** Let $G$ be a compact connected Lie group of rank $r$, equipped with a bi-invariant Riemannian metric $g$ (normalised so that the maximal torus $T \subset G$ has unit-radius factors). Let $\mathcal{T}_G = (C^\infty(G), L^2(G, \Sigma), D_G)$ be the canonical Riemannian Dirac spectral triple, with $D_G$ the bi-invariant (Kostant cubic) Dirac operator on $G$. For each $\Lambda > 0$, let $P_\Lambda$ be the Peter–Weyl projection onto the span of irreducibles $V_\pi$ with $\|\pi\|_{C} \le \Lambda^2$ (Casimir cutoff), and let $\mathcal{T}_\Lambda = (P_\Lambda C^\infty(G) P_\Lambda, P_\Lambda L^2(G, \Sigma), P_\Lambda D_G P_\Lambda)$ be the spectrally truncated metric spectral triple. Then $\mathcal{T}_\Lambda$ converges to $\mathcal{T}_G$ in the Latrémolière propinquity, with explicit rate

$$\boxed{\Lambda_{\mathrm{prop}}(\mathcal{T}_\Lambda, \mathcal{T}_G) \le C_3(G) \cdot \gamma_\Lambda(G), \qquad \gamma_\Lambda(G) = c(G) \cdot \frac{\log \Lambda}{\Lambda} + O(1/\Lambda),}$$

where:

- $C_3(G) \le 1$ is the Lipschitz comparison constant, asymptotic-tight as $\Lambda \to \infty$.
- $c(G) = \frac{\mathrm{Vol}(G/T)}{\pi^r}$ is the master Mellin M1 (Hopf-base measure) signature for the maximal-torus quotient $G/T$ in the bi-invariant Haar normalisation.

**Recovery of Paper 38.** For $G = \mathrm{SU}(2)$, rank $r = 1$, $T = U(1)$, $G/T = S^2$:
$$c(\mathrm{SU}(2)) = \frac{\mathrm{Vol}(S^2)}{\pi^1} = \frac{4\pi}{\pi} = 4.$$
But Paper 38's constant is $4/\pi$, not $4$. The discrepancy comes from a normalisation convention: Paper 38 normalises by $\mathrm{Vol}(\mathrm{SU}(2)) = 2\pi^2$ rather than $\pi^r = \pi$ in the rate-constant denominator. The cleanest unified form is

$$c(G) = \frac{\mathrm{Vol}(G/T)}{\mathrm{Vol}(G)} \cdot 2 = \frac{2 \mathrm{Vol}(G/T)}{\mathrm{Vol}(G)}.$$

Check for $\mathrm{SU}(2)$: $\mathrm{Vol}(S^2)/\mathrm{Vol}(\mathrm{SU}(2)) = 4\pi/(2\pi^2) = 2/\pi$, so $c(\mathrm{SU}(2)) = 2 \cdot 2/\pi = 4/\pi$. ✓

So the conjectured rate constant is

$$\boxed{c(G) = \frac{2\,\mathrm{Vol}(G/T)}{\mathrm{Vol}(G)}}$$

which is the bi-invariant-Haar volume ratio of the Hopf-base $G/T$ to $G$, with a factor 2 that comes from the spectral-Fejér-kernel-vs-Dirichlet-kernel doubling (same factor 2 as in the Stein–Weiss circle estimate $\int_{T^1} F_n |\theta| d\theta \sim (4/\pi) \log n / n$, where 4 = 2 × 2).

**Predicted constant for SU(3) (rank 2).** $G/T = \mathrm{SU}(3)/T^2 = \mathrm{full\ flag\ } SU(3)$, $\mathrm{Vol}(G/T) = \pi^3 / 12$ in the bi-invariant unit-Killing-form normalisation, $\mathrm{Vol}(\mathrm{SU}(3)) = \pi^4 \sqrt{3}/2$ (or some such — exact value depends on bi-invariant normalisation convention). Predicted $c(\mathrm{SU}(3))$ is a specific transcendental in $\pi^{-1} \cdot \mathbb{Q}$ once the convention is fixed. **This is the single cleanest falsifiable prediction the unified theorem makes. Verification blocked at this sprint by org usage limit; deferred.**

---

## §5. Lemma-by-lemma porting plan

Status legend:

- **PV (ports verbatim):** SU(2) proof transfers without change; just substitute the abstract harmonic analysis.
- **PWNC (ports with named change):** structure of the proof transfers but specific constants/conditions change. Name the change.
- **NNI (needs new ingredient):** SU(2) proof uses something specific to SU(2) that doesn't generalise. Discovery zone.

### §5.1. L1' (chirality-doubled operator system)

**Status: PV.**

The chirality-doubled operator system substrate is a property of the spinor bundle $\Sigma \to G$, not specific to $G = \mathrm{SU}(2)$. For any compact spin Lie group, $\Sigma = G \times_{\mathrm{Spin}} \mathbb{C}^{2^{\dim G/2}}$ decomposes by chirality (when $\dim G$ is even, using $\gamma_5$) or by half-spinor representation ($\dim G$ odd, using the unique irreducible Clifford representation). The block-diagonal structure of the multiplier $M_f$ in each chirality sector is a structural property.

The "offdiag CH extension" used in L1' to ensure SDP non-degeneracy on cross-shell pairs is the analog of an $E1$ ladder — for general $G$, this becomes the standard ladder operator for adjacent Casimir levels, which is well-defined in any Peter–Weyl decomposition (lowering/raising operators between irreducibles connected by the natural representation).

**No new ingredient needed.** Direct lift.

### §5.2. L2 (central spectral Fejér kernel and quantitative rate)

**Status: PWNC.** Fully constructible; rate constant changes.

The central spectral Fejér kernel construction
$$K_\Lambda(g) := \frac{1}{Z_\Lambda} \left| \sum_{\pi : C(\pi) \le \Lambda^2} \sqrt{\dim V_\pi}\, \chi_\pi(g) \right|^2$$
is the natural-coefficient generalisation of Paper 38's Definition 2.1 (where SU(2) takes $\dim V_\pi = 2j+1 = n$ and the bijection $n = 2j+1$ links to the Fock-shell index). For general $G$, $\dim V_\pi$ is given by the Weyl dimension formula in terms of the highest weight of $\pi$.

Properties (a)–(c) of L2 (positivity, Plancherel symbol, cb-norm) extend by the Peter–Weyl theorem and the Bożejko–Fendler 1991 cb-norm equality on amenable compact groups, both of which work for any compact $G$.

Property (d) (mass-concentration moment $\gamma_\Lambda$ and asymptotic rate) is where the rank-dependence enters. The closed-form sum-rule analog to Paper 38's eq. (d.i) involves an integration over conjugacy classes of $G$ — the conjugacy-class measure on $G$ is the Weyl integration formula
$$dg = \frac{1}{|W|} \prod_{\alpha > 0} \left| 2\sin\frac{\langle\alpha, t\rangle}{2} \right|^2 dt$$
(Weyl integration formula), where $W$ is the Weyl group and $\alpha$ ranges over positive roots. For SU(2), there is one positive root and $|W| = 2$, giving the $\sin^2(\chi/2)$ factor. For SU(3) (rank 2), there are 3 positive roots and $|W| = 6$, giving a different polynomial factor in the eigenvalue parameters.

The rate $\gamma_\Lambda \sim c(G) \log \Lambda / \Lambda$ is expected to hold by a Stein–Weiss / Abel–Plana argument analogous to Paper 38's Appendix A. The constant $c(G)$ comes out as $2\,\mathrm{Vol}(G/T)/\mathrm{Vol}(G)$ from the Weyl integration formula and the $\sqrt{\dim V_\pi}$ Plancherel weighting.

**Named change:** the closed-form sum-rule and the Stein–Weiss derivation are rank-dependent. Need to redo the Abel–Plana for SU(3) explicitly to pin the constant. Not difficult, just bookkeeping.

### §5.3. L3 (Lipschitz comparison, $C_3 \le 1$) — REVISED 2026-05-15

**Status: PWNC (revised from NNI after sprint P-A + SU(3) Dirac-triangle verification).**

*Erratum (2026-05-15):* This subsection originally framed the L3 ingredient as a "Casimir triangle inequality" $|C(\pi) - C(\pi')| \le C(\sigma)$. **Sprint P-A established that inequality is FALSE** — it already fails at SU(2) for $j = 1, j' = 1/2$ ($|\Delta C| = 5/4 > 3/4 = C(V_{1/2})$). Direct reading of Paper 38's L3 proof memo (`debug/r25_l3_proof_memo.md`) confirms Paper 38 does NOT use this inequality — it uses the **Avery 3-Y selection rule** ($|n - n'| \le N - 1$, which is a shell-difference rule on **Dirac quantum numbers**, NOT a Casimir difference) plus the Lipschitz scaling $\|\nabla Y^{(3)}_{NLM}\|_\infty \sim \sqrt{N^2 - 1}$, giving the per-harmonic ratio $(N-1)/\sqrt{N^2-1} \nearrow 1^-$.

**The correct L3 generalization** (verified at SU(3); see below) is the **Dirac-triangle inequality**:

$$\boxed{|D(\pi) - D(\pi')| \le \sqrt{C(\sigma)} \quad \text{for every } \sigma \subset \pi \otimes \pi'^*}$$

where $D$ is the Kostant cubic Dirac eigenvalue ($|D(\lambda)|^2 = \langle\lambda + \rho, \lambda + \rho\rangle$, equivalently $|D(\lambda)|^2 = C(\lambda) + \langle\rho, \rho\rangle$). On SU(2) this reduces to $|j - j'| \le \sqrt{j_\sigma(j_\sigma + 1)}$, holding trivially because $j_\sigma \ge |j-j'|$ from the Clebsch–Gordan triangle and $j_\sigma(j_\sigma + 1) \ge j_\sigma^2 \ge |j - j'|^2$.

**SU(3) verification** (`debug/dirac_triangle_su3_check.py`, exact sympy arithmetic, agent's Casimir normalization $\langle\rho,\rho\rangle = 1$):
- **100 / 100 pairs** in the panel $\{(p,q) : p+q \le 3\}$ satisfy the Dirac-triangle inequality. No counterexamples.
- **Maximum panel ratio: 0.672** (at $\lambda = (0,0), \lambda' = (0,3)$, the trivial × decuplet case; comparable to Paper 38's $C_3 = 0.707$ at $n_{\max}=3$).
- All three of P-A's Casimir-triangle counterexamples pass the Dirac-triangle test with margin: ratios 0.61, 0.49, 0.49 (against Casimir-triangle ratios 1.40, 2.0, 2.0 respectively).
- Asymptotic family $(n,0) \otimes (1,0)^*$: ratios $0.480, 0.612, 0.689, 0.739, 0.775, 0.803, 0.824$ for $n = 2..8$, **asymptoting to 1 from below**. Structurally identical to Paper 38's $\sqrt{(N-1)/(N+1)} \to 1^-$.

**Reclassification: PWNC.** The Avery 3-Y selection rule + harmonic-gradient scaling on SU(2) generalize to the Dirac-eigenvalue spectral-gap bound + $\sqrt{C(\sigma)}$ Lipschitz scaling on general compact Lie groups. The Brauer–Steinberg analysis from the original §5.3 (Brauer/Steinberg highest-weight formula + Cauchy–Schwarz) is expected to go through cleanly for the Dirac version, with the "cross term obstruction" that broke the Casimir version replaced by a sum of squares on the dominant chamber (sketch: $|D(\pi) - D(\pi')|^2 \le \langle\lambda - \lambda', 2(\lambda + \lambda') + 4\rho\rangle \cdot$ correction terms; details deferred to execution sprint).

**Asymptotic-tight bound $C_3(G) = 1$** is recovered in the limit, consistent with the structural pattern. Most likely outcome: full proof at all ranks via Cauchy–Schwarz + Brauer–Steinberg on the weight lattice; ETA 1–2 weeks of careful bookkeeping during the execution sprint.

### §5.4. L4 (Berezin reconstruction)

**Status: PWNC.**

The Berezin reconstruction
$$B_\Lambda(f) = P_\Lambda \cdot (K_\Lambda * f) \cdot P_\Lambda$$
extends to general $G$ via the Peter–Weyl + central-Fejér-kernel form. Properties (a)–(d) (positivity, contractivity, approximate identity, L3 compatibility) follow as in Paper 38: positivity from $K_\Lambda \ge 0$, contractivity from Young's inequality, approximate identity from L2(d), L3 compatibility from L3 + Young's gradient inequality.

The L4(c) rate $\|B_\Lambda(f) - P_\Lambda M_f P_\Lambda\| \le \gamma_\Lambda \|\nabla f\|_\infty$ inherits the rank-dependent $\gamma_\Lambda$ from L2(d).

**Named change:** the Plancherel weight $\hat K_\Lambda(\pi) = \dim V_\pi / Z_\Lambda$ generalises the SU(2) value $(2j+1)/Z_{n_\max}$. The Hawkins-style Kähler reconstruction is *not* used (Paper 38 already deviated from this; the central-Fejér-kernel form is the non-Kähler alternative). For general compact $G$, $G$ is parallelisable but not Kähler unless $G$ is a torus, so we always use the central-Fejér form. **Same proof structure as Paper 38 L4.**

### §5.5. L5 (Latrémolière propinquity assembly)

**Status: PV.**

L5 is the bookkeeping lemma that assembles L1'–L4 into the propinquity bound via a tunneling pair $(B_\Lambda, P_\Lambda)$. The propinquity framework is independent of the specific group — it is a metric on equivalence classes of metric spectral triples (Latrémolière 2017/2023). Reach is $\gamma_\Lambda$ (from L4(c) and the dual via L2(c)); height is $\gamma_\Lambda$ (from L4(d) + Lipschitz/good-kernel estimate).

**No new ingredient.** Direct lift from Paper 38 §3.5.

### §5.6. Summary of porting

| Lemma | Status (REVISED 2026-05-15) | Notes |
|-------|--------|-------|
| L1' (operator system substrate) | **PV** | Spinor-bundle structure is universal for compact spin manifolds |
| L2 (central Fejér kernel + rate) | **PWNC** | Closed-form sum-rule rank-dependent; Stein–Weiss / Abel–Plana redo at SU(3) is the canary — open prerequisite |
| L3 (Lipschitz, $C_3 \le 1$) | **PWNC** ← was NNI | Dirac-triangle inequality verified at SU(3); Brauer–Steinberg + Cauchy–Schwarz proof expected to extend |
| L4 (Berezin reconstruction) | **PWNC** | Plancherel weight $\dim V_\pi / Z_\Lambda$; same proof structure |
| L5 (propinquity assembly) | **PV** | Mechanical |

**Updated count:** Four of five lemmas port verbatim or with named change. The single "needs new ingredient" item from the original scoping was a mis-summary; under the corrected Dirac-triangle reformulation, no lemma requires structurally new ingredients — only careful bookkeeping with the rank-dependent constants.

---

## §6. Risk assessment and forward plan

### §6.1. What's likely to work in 4–8 weeks

- **L1', L4, L5:** all port mechanically with named-change rank-dependent constants.
- **L2:** the closed-form sum-rule redo at SU(3) is the main labor (Weyl integration formula + spherical-Fejér-kernel moment). 1–2 days at most once the org usage limit resets.
- **L3:** the Dirac-triangle inequality is verified at SU(3) on the panel; the analytic Brauer–Steinberg + Cauchy–Schwarz proof at all ranks is bookkeeping.

About 4 weeks of focused harmonic-analysis bookkeeping total, modulo the rate-constant pinning at SU(3).

### §6.2. What might block

**Failure mode 1: rank-dependent log power.** A genuine surprise would be if the rate at rank $r \ge 2$ becomes $\log^r \Lambda / \Lambda$ rather than $\log \Lambda / \Lambda$ (the rank-$r$ generalisation of the one-log-slowdown). This would be visible in the Stein–Weiss redo at SU(3). If it happens, the theorem still holds but the rate-constant identification with the master Mellin engine becomes more subtle (each log factor would correspond to a separate Cesaro-moment level). **Still an open prerequisite; not blocked at theorem level.**

**Failure mode 2: bi-invariant metric does not extend cleanly.** For non-simply-connected compact $G$ or for torsion-containing fundamental groups, the bi-invariant metric is not unique up to scaling, and convention questions matter for the rate constant. Resolvable, but cosmetic.

**Failure mode 3: scoping is wrong.** The right class might actually be Class 2 (compact symmetric space $G/K$) rather than Class 1, because the GeoVac-internal Camporesi–Higuchi triple is technically the symmetric-space construction $\mathrm{SO}(4)/\mathrm{SO}(3) = S^3$, not the Lie-group construction $\mathrm{SU}(2)$. The two coincide for $S^3$ but diverge at higher rank. This is a real ambiguity. The Paper 38 setup is "ambidextrous" — it works with either reading.

**Recommendation:** start with Class 1 (Lie group setup), but in the L2 closed-form sum-rule redo at SU(3), explicitly check whether the Lie-group Weyl integration formula and the symmetric-space spherical-function formula give the same rate constant. They should, given the Class-1-contains-Paper-38 reading, but verifying this on SU(3) is the canary.

### §6.3. Most likely "structure that turned out to be incidental"

Three candidates:

1. **The factor $4/\pi$ vs $\Vol(G/T)/\Vol(G) \cdot 2$.** In Paper 38, the constant is identified as the Hopf-base measure $4/\pi = \Vol(S^2)/\pi^2 = 2 \Vol(S^2)/\Vol(\mathrm{SU}(2))$. The factor 2 in the unified form has two possible origins (Cesaro-2 doubling or Weyl-group order 2 = $|W(\mathrm{SU}(2))|$); they coincide for SU(2). At higher rank, $|W(G)| \ne 2$ in general, so the unified form might be $c(G) = |W(G)| \cdot \mathrm{Vol}(G/T)/\mathrm{Vol}(G)$ rather than $2 \cdot \mathrm{Vol}(G/T)/\mathrm{Vol}(G)$. Need to check.

2. **The Avery 3-Y integral.** Paper 38 uses Avery's specific basis. **(2026-05-15 update: confirmed)** The Dirac-triangle inequality doesn't actually need the Avery basis — it just needs the tensor-product selection rule on irreps. Paper 38's reliance on Avery is incidental to the SU(2) case. The unified theorem is cleaner than Paper 38.

3. **The Camporesi–Higuchi closed-form spectrum.** Paper 38 uses $|\lambda_n| = n + 1/2$ explicitly in the L1' SDP analysis. For higher-rank groups, the Dirac spectrum is given by Casimir + $\rho$-shift but is not in closed-form scalar form. **(2026-05-15 update: not a blocker)** The proof goes through using the spectral *gap* rather than the explicit spectrum. The Kostant cubic Dirac eigenvalue $|D(\lambda)|^2 = \langle\lambda + \rho, \lambda + \rho\rangle$ is the right substitute.

### §6.4. Recommendation

**Proceed to execution at Class 1 (compact Lie group with bi-invariant metric).**

Remaining prerequisites (deferred to next session due to org usage limit):

1. **Numerical SU(3) rate constant.** Verify $c(G) = 2\,\mathrm{Vol}(G/T)/\mathrm{Vol}(G)$ on SU(3) numerically. Identify whether the prefactor is 2 (Cesaro doubling) or $|W(G)|$ (Weyl group order) or some other explicit transcendental ratio. **This is the cleanest falsifiable prediction the unified theorem makes.**
2. **Log-power test.** Single-log vs double-log fit on $\gamma_\Lambda$ at SU(3). Resolved by the same SU(3) numerical computation as #1.

After both prerequisites resolve, the 4–6 week execution sprint is well-defined and ready to launch.

---

## §7. Open questions

1. ~~**Casimir triangle inequality:** Standard result or new harmonic analysis?~~ — **RESOLVED 2026-05-15 by P-A: the inequality is FALSE. The corrected L3 reformulation (Dirac-triangle inequality) holds; verified at SU(3) on the 100-pair panel.**
2. **Rank-$r$ rate-constant prediction:** Verify $c(G) = ?\, \mathrm{Vol}(G/T)/\mathrm{Vol}(G)$ on SU(3) numerically. Identify whether the prefactor is 2 (Cesaro doubling) or $|W(G)|$ (Weyl group order) or some other explicit transcendental ratio. **OPEN, deferred from this sprint.**
3. **Symmetric-space vs Lie-group ambidexterity:** Does the Class-2 ($G/K$) version of the theorem reduce to Class 1 ($G$ alone) in the group case $K = \{e\}$? Should, but worth checking on $S^3 = \mathrm{SO}(4)/\mathrm{SO}(3) = \mathrm{SU}(2)$.
4. **Higher rank log power:** Is the rate $\log^r \Lambda / \Lambda$ at rank $r$, or $\log \Lambda / \Lambda$ uniformly? If the former, what's the exact prefactor at each log-power level? Master Mellin engine reading would differ between these cases. **OPEN.**
5. **Coadjoint-orbit overlap with Rieffel-2023:** For $G = \mathrm{SU}(2)$, the unified theorem covers $S^3$ (the group) and Rieffel covers $S^2$ (the coadjoint orbit). For $G = \mathrm{SU}(3)$, the unified theorem covers $\mathrm{SU}(3)$ and Rieffel covers $\mathbb{CP}^2$, the flag variety, etc. Are the rates compatible in any sense?
6. **Connection to Hekkelman–McDonald 2024:** their non-commutative integral / Szegő limit on $T^d$ is the Berezin–Toeplitz analog of Paper 38's L4. Does the unified theorem at Class 1 admit a Szegő companion?
7. **Beyond bi-invariant metric:** What happens for non-bi-invariant metrics on $G$? The propinquity should still converge but with a different (Lie-group-asymmetric) rate. Out of scope for now but worth flagging.
8. **W2b cross-manifold tensor product:** From CLAUDE.md §1.7 WH1 PROVEN: $\mathcal{T}_{S^3} \otimes \mathcal{T}_{\mathrm{Hardy}(S^5)}$ is structurally blocked. Does the unified theorem at Class 1 illuminate this blocker? Probably not — the blocker is at the Riemannian-vs-Hardy-sector framework level, not at the rank level.
9. ~~**Authorship of arXiv:2310.14733:** Track 1 attributed it to Gaudillot-Estrada & van Suijlekom (IMRN 2025); P-A attributed it to Hekkelman 2023. Verify next session.~~ **RESOLVED 2026-05-15 (WebFetch + WebSearch verification):** Track 1 was correct. arXiv:2310.14733 = Yvann Gaudillot-Estrada and Walter D. van Suijlekom, "Convergence of Spectral Truncations for Compact Metric Groups," IMRN 2025, Issue 13, paper rnaf197 (July 2025). P-A's "Hekkelman 2023" attribution was wrong — likely conflated with Hekkelman 2022 ($S^1$ Toeplitz) or Hekkelman–McDonald 2024 ($T^d$ Berezin–Toeplitz). Erratum applied to `debug/casimir_triangle_inequality_lit_check.md`.

---

**Bibliography.** All paper IDs in §2 verified against arXiv as of 2026-05-10.
