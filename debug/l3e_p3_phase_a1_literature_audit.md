# Sprint L3e-P3-A.1 — Literature Audit Memo

**Date:** 2026-05-23
**Sprint:** L3e-P3-A.1 (literature deep-dive for Pointed Latrémolière propinquity via Mondino-Sämann synthetic Lorentzian GH bridge)
**Predecessor:** `debug/sprint_l3e_p3_detailed_scoping_memo.md` (multi-month scoping), `debug/l3_literature_audit_memo.md` (May-17 baseline audit)
**Method:** WebSearch + WebFetch across arXiv abstracts, HTML, and PDF where accessible. ~10 fetches + 8 web searches.
**Honest scope:** abstract-level scan with HTML extraction of definitions/theorems where available. PDF extraction was unreliable (text streams compressed); HTML rendered better. Verdicts at the structural-claim level, not full proof level.

---

## §1. Audit purpose and scope

Phase A.1 of Sprint L3e-P3 is the literature deep-dive that determines whether the planned bridge construction (Phases A.2–A.4: define operator-algebraic ε-net analog, prove correspondence to Mondino-Sämann causal-diamond ε-net, exhibit GeoVac wedge as ε-net) is well-targeted given the current published landscape. The scoping memo's §6 named four concurrent-work risks, with Mondino-Sämann moving into operator algebras at "MEDIUM probability" and Latrémolière extending to Lorentzian at "LOW probability." This audit re-verifies those probability estimates against publications through May 2026.

**Three target streams** (per the dispatch):
1. **Mondino-Sämann synthetic Lorentzian GH lineage** — verify current state and extract precise definitions of causal-diamond ε-net, time-separation function τ, pre-compactness theorem, and unbounded-case extensions.
2. **Latrémolière + Farsi-Latrémolière propinquity 2024-2025 updates** — re-verify the May-17 audit's verdict that "no published Lorentzian propinquity exists" and look for any pointed/locally-compact extension that might pre-empt Sprint L3e-P3-B.
3. **Hekkelman-McDonald non-compact assessment** — assess whether their Tauberian/spectral-asymptotics framework complements or competes with pointed propinquity.

**Net headline finding (preview, full discussion §3 and §5):** Latrémolière himself published on December 3, 2025, **"The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces"** (arXiv:2512.03573), introducing exactly the operator-algebraic pointed-quantum-metric framework that Sprint L3e-P3 Phase B was scoped to construct from scratch over 6 months. This is a **partial scoop**, not a total one: Latrémolière 2512.03573 is Riemannian and unitless of Krein-wedge content, so the GeoVac-specific Lorentzian closure (Phase C, Paper 49) is unscooped, but Phase B's standalone Paper 48 ("Pointed Latrémolière propinquity for non-compact quantum metric spaces") would substantially duplicate Latrémolière's December 2025 contribution. The scoping memo's risk assessment of "Latrémolière scoop-risk = LOW" (May 17) is now obsolete — it became HIGH between May 17 and December 3, and we are post-publication on the operator-algebraic side. This re-prices the entire Sprint L3e-P3 program; recommendations in §6.

---

## §2. Stream 1 — Mondino-Sämann synthetic Lorentzian GH lineage

### §2.1 Canonical references and current versions

| arXiv | Authors | Title | Date / latest version | Status |
|:------|:--------|:------|:----------------------|:-------|
| 2209.14384 | Müller et al. | Lorentzian metric spaces & their GH-convergence | 2022, v3 2025 | Published / canonical |
| 2209.12736 | (multi-author) | GH metrics & dimensions of Lorentzian length spaces | 2022 | Published |
| 2412.04311 | Bykov, Minguzzi, Suhr | Lorentzian metric spaces & GH-convergence: unbounded case | Dec 2024, rev May 2025 | Published |
| 2504.10380 | Mondino, Sämann | Lorentzian GH convergence and pre-compactness | April 2025, **v4 Dec 9, 2025** | Active; key reference |
| 2506.10852 | Braun, Sämann | Gromov's reconstruction theorem & measured Lorentzian GH | June 2025 | Published |
| 2510.13069 | Che, Perales, Sormani | Synthetic Lorentzian GH convergence (Sakovich-Sormani lineage) | Oct 2025 | Published |
| 2510.24423 | Minguzzi | Results on Lorentzian metric spaces | Oct 2025, *Gen. Rel. Grav.* 58 (2026) | Published |
| 2511.18389 | (multi-author) | Intrinsic Timed Hausdorff Convergence | Nov 2025 | Published |
| 2605.09101 | Kubota | Lorentzian coarea inequality | May 9, 2026 | Just released; uses Lorentzian Hausdorff measure (McCann-Sämann) |
| 2605.11271 | (multi-author) | Convergence of Lorentzian spaces and curvature bounds for generalized cones | May 2026 | Just released |
| 2605.03172 | (multi-author) | Stability of Synthetic Timelike Ricci Bounds under C⁰-Limits | May 2026 | Just released |

The lineage has accelerated through 2025 and the first half of 2026. Eight papers in the May 2025 – May 2026 window. The community is mature and actively publishing.

### §2.2 Precise definitions extracted from Mondino-Sämann arXiv:2504.10380 v4

From the HTML rendering of v4 (December 9, 2025), I extracted the following definitions and theorems:

**Definition 2.1 (chronological / causal relations on a Lorentzian pre-length space).** Given an extended time-separation function $\ell$ satisfying the reverse triangle inequality:
- Timelike: $\ll \;:=\; \ell^{-1}((0, \infty])$
- Causal: $\le \;:=\; \ell^{-1}([0, \infty])$

**Definition 2.3 (Lorentzian pre-length space).** A Lorentzian pre-length space carries the extended time-separation $\ell$ and the time-separation function $\tau := \max(0, \ell)$.

**Definition 3.2 (ε-net).** "An ε-net $S$ for $A$ is a collection of causal diamonds $S = (J_i)_{i \in \Omega}$ satisfying:
- (i) $\tau(J_i) \le \varepsilon$ for all $i \in \Omega$;
- (ii) $A \subseteq \bigcup_{i \in \Omega} J_i$."

(Causal diamonds $J(p, q) := \{x : p \le x \le q\}$ for $p, q$ in the space with $p \le q$, and $\tau(J_i) := \tau(p_i, q_i)$.)

**Definition 3.6 (LGH-convergence of subsets).** $A_n \xrightarrow{LGH} A$ requires four ingredients:
- (i) equal-cardinality finite ε-nets;
- (ii) correspondences between vertex sets with distortion $\to 0$;
- (iii) extension property preserving distortion bounds;
- (iv) forward density: vertices of $A_n$ approach points of $A$ via monotone sequences.

**Definition 3.8 (covered Lorentzian pre-length space).** A tuple $(X, \ell, o, \mathcal{U})$ where $o \in X$ is a base point and $\mathcal{U} = (U_k)_{k \in \mathbb{N}}$ is a countable family with:
- $U_k \subseteq U_{k+1}$,
- $o \in U_k$ for all $k$,
- $\sup_{x, y \in U_k} \tau(x, y) < \infty$ (bounded timelike diameter on each $U_k$).

**Definition 3.12 (pLGH-convergence — pointed Lorentzian Gromov-Hausdorff).** Applies the LGH framework to **covered** Lorentzian pre-length spaces, requiring $U_{k,n} \xrightarrow{LGH} U_{k,\infty}$ as $n \to \infty$ for each cover level $k$.

**Theorem 6.2 (Pre-compactness I — pointed version).** "Any sequence in $\mathfrak{X}$ has a converging subsequence; i.e., for any sequence $((X_n, \ell_n, o_n, \mathcal{U}_n))_{n \in \mathbb{N}} \subset \mathfrak{X}$ there exists a subsequence $(n_j)_j \subset \mathbb{N}$ and a covered Lorentzian pre-length space $(X, \ell, o, \mathcal{U})$ such that $(X_{n_j}, \ell_{n_j}, o_{n_j}, \mathcal{U}_{n_j}) \xrightarrow{pLGH} (X, \ell, o, \mathcal{U})$ as $j \to \infty$."

**Remark 3.9 (why pointed and covered are needed).** "unlike metric spaces lacking natural exhaustions, Lorentzian spaces require explicit covers since time-separation level-sets are typically non-compact."

**Critical correction to the May-17 audit.** The May-17 memo treated Mondino-Sämann 2504.10380 as compact-only and named the non-compact extension as a separate problem. **In v4 (December 9, 2025) the paper already covers the non-compact case via the pLGH-convergence on covered Lorentzian pre-length spaces.** The extension to non-compact carriers on the synthetic side is published. The scoping memo §2 Phase A.1 was correct to list 2412.04311 (Bykov-Minguzzi-Suhr) as the non-compact extension, but Mondino-Sämann v4 has already incorporated a parallel construction (covers + base points).

### §2.3 Bykov-Minguzzi-Suhr unbounded extension (arXiv:2412.04311)

Bykov-Minguzzi-Suhr drop the boundedness condition entirely. Their characterization of unbounded Lorentzian metric spaces uses three minimal axioms:
- (a) the reverse triangle inequality for chronologically related events;
- (b) Lorentzian distance continuity and relative compactness of chronological diamonds;
- (c) a distinguishing condition via the Lorentzian distance function.

A "countably generating condition" is added to ensure the Polish property. They prove (pre)length spaces are GH-stable under this extension.

**Mechanism for non-compactness on the synthetic side:** rather than restrict to bounded carriers, they characterize the natural quasi-uniformity (resp. quasi-metric) on sequenced Lorentzian metric spaces, and GH-convergence is defined directly on this quasi-uniform structure.

**Comparison with Mondino-Sämann pLGH:** Bykov-Minguzzi-Suhr work in the *unpointed* unbounded setting; Mondino-Sämann v4 work in the *pointed* unbounded setting via covered structures. Both are published. The two approaches are complementary: pointed convergence (Mondino-Sämann) gives a basepoint-aware limit (good for cosmological models with a distinguished event), unpointed (Bykov-Minguzzi-Suhr) gives basepoint-free limit.

### §2.4 Recent extensions and 2026 follow-ups

- **Braun-Sämann 2506.10852** (June 2025): "Gromov's reconstruction theorem & measured Lorentzian GH" — adds a measure structure. Defines normalized bounded Lorentzian metric measure spaces with three notions of convergence. Bounded only.
- **Che-Perales-Sormani 2510.13069** (Oct 2025): "Synthetic Lorentzian Gromov-Hausdorff convergence" — Sakovich-Sormani lineage; introduces "timed metric spaces" and "intrinsic timed-Hausdorff convergence." Distinct from Mondino-Sämann's pre-length-space approach (uses time functions rather than time-separation functions). Two parallel synthetic frameworks now exist.
- **Minguzzi 2510.24423** (Oct 2025, *Gen. Rel. Grav.* 58, 2026): "Results on Lorentzian metric spaces" — proves every Lorentzian metric space admits a Cauchy time function; constructive proof, novel even for smooth spacetimes.
- **2511.18389**: "Intrinsic Timed Hausdorff Convergence" — extends Sakovich-Sormani distance notions; shows intrinsic timed-Hausdorff convergence implies (timeless) GH-convergence and big-bang convergence.
- **2605.09101** (May 9, 2026): Kubota, Lorentzian coarea inequality — uses Lorentzian Hausdorff measure (McCann-Sämann), establishes coarea formula. Recent.
- **2605.11271** (May 2026): Convergence of Lorentzian spaces and curvature bounds for generalized cones. Recent.
- **2605.03172** (May 2026): Stability of synthetic timelike Ricci bounds under C⁰-limits, applications to impulsive gravitational waves. Recent.

### §2.5 Definitions table — synthetic-side objects relevant to bridge construction

| Object | Synthetic-side definition (Mondino-Sämann v4) | Phase A.2 target operator-algebraic analog |
|:-------|:------------------------------------------------|:--------------------------------------------|
| Time-separation function $\tau$ | $\tau := \max(0, \ell)$ with reverse triangle inequality | KMS-correlation function $\tau_L(\omega_1, \omega_2)$ on Krein-positive state pairs (per scoping memo §2 A.2) |
| Causal diamond $J(p, q)$ | $\{x : p \le x \le q\}$ for $p \le q$ | Modular-flow time-band on wedge $W_L$ at modular parameter $\beta = \varepsilon^{-1}$ (per scoping memo §2 A.2) |
| ε-net | Collection of causal diamonds with $\tau$-diameter $\le \varepsilon$ covering $A$ | Truncated operator system $O^L_{n_{\max}, N_t}$ at cutoff $(n_{\max}, N_t)$ with $\varepsilon = O(\max(1/n_{\max}, 1/N_t))$ |
| Cover $\mathcal{U} = (U_k)$ | $U_k \subseteq U_{k+1}$, $o \in U_k$, bounded $\tau$-diameter on each | Increasing family of truncations at varied cutoffs |
| Base point $o$ | Distinguished event in $X$ | BW vacuum state $\omega_W^L = e^{-K_\alpha^W}/Z$ (Paper 43 §4.2) |
| Pre-compactness | Theorem 6.2: subsequences converge in pLGH | Open: pre-compactness of Krein-pointed operator-system sequences |

The Phase A.2 target operator-algebraic ε-net construction is well-defined as a research program. Each row in the table is a concrete bridge target. The structural mismatch the scoping memo §2 Phase A flagged ("synthetic uses causal diamonds, operator-algebraic uses Hilbert inner product") is genuine but not insurmountable — the BW vacuum gives both sides a canonical base point, and the modular flow gives both sides a notion of "time-separation diameter."

---

## §3. Stream 2 — Latrémolière + Farsi-Latrémolière trajectory 2024–2026

### §3.1 Headline finding: Latrémolière arXiv:2512.03573 (December 3, 2025)

**This is the dominant finding of the audit and re-prices Sprint L3e-P3-B in particular.**

| Field | Value |
|:------|:------|
| Title | The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces |
| Author | Frédéric Latrémolière (solo) |
| arXiv | 2512.03573 |
| Submitted | 3 December 2025 |
| Category | math.OA (primary), MSC 46L89, 46L30, 58B34 |
| Length / structure | 6 sections + appendix (estimated 30+ pages from PDF size) |
| Predecessor | arXiv:1406.0233 (Latrémolière 2014, "Topographic Gromov-Hausdorff quantum Hypertopology for Quantum Proper Metric Spaces") |

**Abstract (paraphrased from extraction):** Introduces a hypertopology on **pointed proper quantum metric spaces** — separable, possibly non-unital C*-algebras endowed with a Lipschitz seminorm analog, a **distinguished state**, and a particular type of approximate units. Provides an analog of the GH distance on proper metric spaces. **When restricted to quantum compact metric spaces, the new topology is compatible with the topology of the Gromov-Hausdorff propinquity.** Includes new examples of noncompact, noncommutative pointed proper quantum metric spaces as limits of finite-dimensional QCMS.

**Definitions extracted from HTML rendering of 2512.03573:**

**Definition 1.22 (Separable quantum locally compact metric space).** A pair $(\mathfrak{A}, \mathsf{L})$ where $\mathfrak{A}$ is a C*-algebra and $\mathsf{L}$ is a **Leibniz hermitian norm-modulo-constants** satisfying:
- (i) the Fortet-Mourier distance metrizes the weak* topology on the state space,
- (ii) closed unit-ball closure property,
- (iii) bounded Monge-Kantorovich tightness condition.

**Definition 1.26 (Pinned separable quantum locally compact metric space).** A triple $(\mathfrak{A}, \mathsf{L}, \mu)$ where $(\mathfrak{A}, \mathsf{L})$ is a separable QLCMS and $\mu \in \mathcal{S}(\mathfrak{A})$ is a **distinguished state** such that the states with finite Monge-Kantorovich distance from $\mu$ form a weak*-dense subset. The state $\mu$ is called **the pin**.

**Definition 1.29 (L-Lipschitz μ-pinned exhaustive sequence).** Used to characterize proper quantum metric spaces; the precise definition involves exhaustion via Monge-Kantorovich distance.

**Section structure of 2512.03573:**
1. Pointed Proper Quantum Metric Spaces (1.1 locally compact, 1.2 proper)
2. Local Quantum Metametrics (2.1 tunnels, 2.2 metametrics)
3. Coincidence property of the local metametrics (3.1 target sets, 3.2 property)
4. Gromov-Hausdorff quantum metametric (4.1 metametric, 4.2 hypertopology)
5. Classical and quantum compact cases
6. Noncompact noncommutative example: $c_0(\mathbb{Z}) \rtimes_\alpha \mathbb{Z}$ as limit of finite-dim matrix algebras

The construction is a **metametric** (relaxed triangle inequality), explicitly described by the author as "an object beyond metrics." It satisfies distance-zero ⟹ "full isometric" equivalence. The full GH quantum metametric (§4) induces the hypertopology.

**Compact-case agreement (Lemma 1.23 and §5):** when $\mathfrak{A}$ is unital, the pinned separable QLCMS is a QCMS in the Latrémolière propinquity sense, and the new topology restricts to the propinquity topology. (This was the Phase B.3 compact-case-agreement theorem the scoping memo planned to prove.)

**Non-unital examples (Example 1.12, §6):** $C_0(X)$ for locally compact $X$; $c_0(\mathbb{Z}) \rtimes_\alpha \mathbb{Z}$ crossed product as limit of finite-dim matrix algebras. The non-compact extension via base-state restriction (Phase B.4 in the scoping memo) is explicitly worked out.

**Lorentzian / Krein content: NONE.** The paper is entirely Riemannian / C*-algebra in scope. No reference to Lorentzian structures, Krein spaces, or BBB signature classification.

### §3.2 What this means for Sprint L3e-P3-B

The scoping memo §3 Phase B planned four sub-deliverables:
- B.1 Pointed QCMS definition (4–6 weeks)
- B.2 Pointed tunneling pair (4–6 weeks)
- B.3 Pointed propinquity metric + compact-case agreement (4–6 weeks)
- B.4 Non-compact extension via base-state restriction (4–6 weeks)
- B.5 Paper 48 draft (2–4 weeks)

**Latrémolière 2512.03573 implements B.1–B.4 already.** The "pinned separable QLCMS" is the pointed QCMS (B.1). The "local quantum metametric" via tunnels is the pointed tunneling pair (B.2). The §4 GH quantum metametric is the pointed propinquity (B.3). The §6 $c_0(\mathbb{Z}) \rtimes_\alpha \mathbb{Z}$ example is the non-compact extension (B.4). The compact-case agreement is Lemma 1.23 (B.3).

**Paper 48 in its scoped form (a standalone math.OA paper introducing pointed Latrémolière propinquity) would now substantially duplicate Latrémolière 2512.03573.** This is not an irrelevant overlap — it's an essentially-complete scoop on the operator-algebraic side. Sprint L3e-P3-B as currently scoped is no longer the right next move.

### §3.3 Related 2024–2026 Latrémolière / Farsi-Latrémolière work (still Riemannian)

The May-17 audit's verdict that Latrémolière 2024–2026 stays Riemannian is **confirmed** for these papers:
- **arXiv:2404.00240** (Farsi-Latrémolière, v3 Oct 2025): Collapse in NCG and spectral continuity, Riemannian spin manifolds, U(1) principal bundles. No Lorentzian content.
- **arXiv:2504.11715** (Farsi-Latrémolière, April 2025): Continuity for spectral propinquity of Dirac operators on analytic path of Riemannian metrics. Title is explicit.
- **arXiv:2602.23080** (early 2026): "Noncommutative coarse metric geometry" — develops bridge between Latrémolière proper QCMS and W*-metric approach. Coarse-geometric, not Lorentzian.

**What's specifically new in 2025–2026 relative to the May-17 baseline:**
- arXiv:2512.03573 (Latrémolière, Dec 2025): the pointed-QCMS hypertopology — described above; this is the major new entry.
- arXiv:2604.02117 (early 2026): "The Bures metric and the quantum metric on the density space of a C*-algebra: the non-unital case" — non-unital Bures-quantum-metric framework; related to Latrémolière's QLCMS via the "quantum Lipschitz triple" construction.
- arXiv:2602.23080: noncommutative coarse metric geometry, bridges proper QMS to coarse and W*-metric approaches.

**No Lorentzian Latrémolière publication.** The author's trajectory remains firmly Riemannian. The May-17 LOW-scoop-risk verdict on Lorentzian-side Latrémolière contributions is unchanged. **But the May-17 verdict did not anticipate that Latrémolière would publish a pointed QMS hypertopology in December 2025, scooping Phase B's operator-algebraic content.**

### §3.4 Critical re-pricing of Sprint L3e-P3-B scoop risk

The scoping memo §6 concurrent-work table priced "Independent pointed-propinquity definition" as "Always possible, mitigation: pre-submit Paper 48 to arXiv ASAP." The mitigation window has closed: Latrémolière pre-submitted his version on December 3, 2025, six months before this audit.

**Re-pricing:**
- Mondino-Sämann moving into operator algebras: probability LOW (no signal; pure synthetic trajectory).
- Latrémolière moving into Lorentzian: probability LOW (no signal in 2024–2026 publications).
- **Independent pointed-propinquity definition: PROBABILITY = 1, ALREADY HAPPENED.**
- Hekkelman-McDonald non-compact extension: probability MEDIUM in 2026–2027 in their Tauberian direction; assessment in §4 below.

---

## §4. Stream 3 — Hekkelman-McDonald non-compact assessment

### §4.1 Trajectory and current papers

| arXiv | Authors | Title | Date | Framework |
|:------|:--------|:------|:-----|:----------|
| 2412.00628 | Hekkelman, McDonald | A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity | Dec 2024, accepted *J. Funct. Anal.* Jul 2025 | Tauberian / Szegő limit / quantum ergodicity |
| 2106.02235 | (earlier work) | Semiclassical Weyl law and exact spectral asymptotics in NCG | 2021 | Tauberian |
| 2604.15008 | Ponge (monograph) | Noncommutative Geometry, Spectral Asymptotics, and Semiclassical Analysis | 2026 | Generalizes McDonald-Sukochev-Zanin to Condition (W) |

**Framework characterization (from extraction).** Hekkelman-McDonald 2412.00628 develops semiclassical Weyl laws and a noncommutative integral on **spectrally truncated spectral triples** (Connes-van Suijlekom paradigm). Key tools: Szegő limit formula, quantum ergodicity, vacuum-state uniqueness for C*-dynamical systems. **This is Tauberian / spectral-asymptotics analysis, NOT Gromov-Hausdorff-style metric convergence.**

The framework treats the **truncated → full continuum limit** at the level of spectral integrals (Connes-style noncommutative integration), not at the level of state-space distance functions. It is structurally complementary to propinquity-style metric convergence.

**Compact vs. non-compact.** The abstract specifies "compact spectral triples." The framework is positioned within the compact case. Non-compact extensions appear to live in the Ponge 2026 monograph (arXiv:2604.15008), which generalizes McDonald-Sukochev-Zanin to "Condition (W)" and explicitly covers "open manifolds with conformally cusp metrics of finite volume."

**No pointed structure.** Neither Hekkelman-McDonald nor Ponge introduce basepoint structures. The framework is integral-theoretic and asymptotic, working at the level of trace formulae rather than state-space metrics.

### §4.2 Complementary or competing with pointed propinquity?

**Complementary.** Hekkelman-McDonald-Ponge work in a different category than pointed propinquity:
- Tauberian / integral-theoretic: spectral functions, traces, Connes' integral, Weyl laws.
- Pointed propinquity (Latrémolière 2512.03573): state-space metric, Monge-Kantorovich distance, hypertopology.

These two frameworks could in principle coexist on a single non-compact spectral triple. For the GeoVac G2-metric question, the Hekkelman-McDonald framework gives a *spectral* statement about cutoff removal (already addressed in Paper 47 via norm-resolvent convergence), and the Latrémolière 2512.03573 framework gives a *metric* statement (the target of Sprint L3e-P3 Phase C).

**No bridge published.** No paper in the Hekkelman-McDonald-Ponge lineage incorporates propinquity, and no paper in the Latrémolière lineage incorporates Tauberian Weyl laws. The two threads run in parallel.

**For Sprint L3e-P3 specifically.** The Hekkelman-McDonald framework does not provide an alternative or competing scaffold for the bridge construction. It also doesn't pose a scoop risk for Sprint L3e-P3-C (the GeoVac G2-metric closure at the pointed-propinquity level). The Ponge 2026 monograph's non-compact spectral-triple extension is at the integral / asymptotic level, complementary to the metric-propinquity direction.

### §4.3 What about an explicit Hekkelman-McDonald 2026 extension?

The May-17 audit's note on the Hekkelman-McDonald 2024 paper as a Tauberian residue of the master Mellin engine still applies. No subsequent paper from Hekkelman or McDonald has extended the Tauberian framework to Lorentzian or to propinquity-style convergence as of May 2026. The monograph form (Ponge 2604.15008) consolidates the asymptotic / integral side without crossing into metric convergence.

---

## §5. Cross-stream synthesis

### §5.1 Where the published bridge would sit

The Phase A bridge construction (synthetic Lorentzian ε-net ↔ operator-algebraic ε-net) would connect:
- **Synthetic side (PUBLISHED):** Mondino-Sämann pLGH-convergence on covered Lorentzian pre-length spaces (2504.10380 v4, Theorem 6.2 pre-compactness), or Bykov-Minguzzi-Suhr unbounded LGH (2412.04311).
- **Operator-algebraic side (PUBLISHED):** Latrémolière pointed proper QMS hypertopology (2512.03573, GH quantum metametric §4).

**No paper bridges the two as of May 2026.** The Mondino-Sämann lineage has not crossed into C*-algebras / spectral triples. The Latrémolière 2512.03573 paper has no Lorentzian content. Neither author cites the other. The two communities operate in disjoint subject areas (math.DG / math.MG for synthetic Lorentzian; math.OA for noncommutative metric).

**Strategic position:** the bridge IS unfilled, but the two pillars are both now mature and PUBLISHED. Sprint L3e-P3-A.2 (define operator-algebraic ε-net analog) can now leverage Latrémolière 2512.03573 directly — the pointed-QCMS framework is no longer Phase B output but Phase B INPUT. The bridge construction becomes much cleaner because the operator-algebraic side already has the pointed-pin structure formalized.

### §5.2 What's been published vs what's missing

| Component | Status | Source |
|:----------|:-------|:-------|
| Synthetic Lorentzian GH (compact) | PUBLISHED | Müller et al. 2209.14384; Mondino-Sämann 2504.10380 |
| Synthetic Lorentzian GH (unbounded, unpointed) | PUBLISHED | Bykov-Minguzzi-Suhr 2412.04311 |
| Synthetic Lorentzian GH (pointed, covered) | PUBLISHED | Mondino-Sämann 2504.10380 v4 §3.8, §6.2 |
| Lorentzian Hausdorff measure | PUBLISHED | McCann-Sämann (cited by Kubota 2605.09101) |
| Causal-diamond ε-net definition | PUBLISHED | Mondino-Sämann 2504.10380 v4 Def 3.2 |
| Operator-algebraic Lorentzian spectral triple at fixed cutoff | PUBLISHED | Bizi-Brouder-Besnard 2018, van den Dungen 2016, GeoVac Paper 43 |
| Operator-algebraic pointed propinquity (compact, Riemannian) | PUBLISHED | Latrémolière 2512.03573 |
| Operator-algebraic pointed propinquity (non-compact, Riemannian) | PUBLISHED | Latrémolière 2512.03573 §4–§6 |
| Operator-algebraic Lorentzian propinquity (any flavor) | **NOT PUBLISHED** | Open |
| Bridge: synthetic Lorentzian ε-net ↔ Latrémolière pointed QMS | **NOT PUBLISHED** | Open — was Sprint L3e-P3 Phase A |
| Operator-algebraic Lorentzian pointed propinquity | **NOT PUBLISHED** | Open — was Sprint L3e-P3 Phase B |
| GeoVac G2-metric closure | **NOT PUBLISHED** | Open — was Sprint L3e-P3 Phase C |

**Two layers of openness remain:**
1. **The Lorentzian extension of Latrémolière 2512.03573.** Taking the pointed-QMS framework and extending it from Riemannian to Lorentzian / Krein-signature. No published precedent.
2. **The bridge construction.** Connecting Mondino-Sämann pLGH (synthetic side) to a hypothetical Lorentzian extension of Latrémolière 2512.03573 (operator-algebraic side). No published precedent.

Sprint L3e-P3 originally bundled (1) + (2) + (Phase C application) into one 12-month program. With Latrémolière 2512.03573 published, the program redivides: the (Riemannian) pointed-QMS framework is now an INPUT, and the original (1) + (2) become the remaining novelty.

---

## §6. Recommendations for Phase A.2 design

### §6.1 Recommendation: refocus Phase A on the Lorentzian extension of Latrémolière 2512.03573

The original Phase A.2 (define operator-algebraic ε-net analog) and Phase A.3 (correspondence theorem) remain valid, but the operator-algebraic substrate should now be **Latrémolière's pointed proper QMS framework (2512.03573) lifted to the Krein-signature / Lorentzian setting**, not a from-scratch construction.

**Re-scoped Phase A deliverables:**

**A.1' (this memo, ~1 week, COMPLETE).** Literature audit confirming Latrémolière 2512.03573 as the operator-algebraic input.

**A.2' — Krein-lift of Latrémolière 2512.03573 pointed QMS (3–6 weeks).** Take Latrémolière's pinned separable quantum locally compact metric space definition (Def 1.22, 1.26) and lift it to the Krein-signature setting: replace the Hilbert-space-valued Lipschitz seminorm with a Krein-self-adjoint version, replace the distinguished state μ with the BW vacuum $\omega_W^L$, replace the Monge-Kantorovich distance with a Krein-positive Wasserstein-Kantorovich distance (using Paper 44's Krein-positive state space). Check whether Latrémolière's three axioms (Fortet-Mourier metrization, closed-ball closure, Monge-Kantorovich tightness) transport to the Krein setting.

**A.3' — Bridge to Mondino-Sämann pLGH (3–4 weeks).** With the Krein-lifted pointed QMS framework in hand, prove the correspondence: every truncated Lorentzian Krein spectral triple $\mathcal{T}^L_{n_{\max}, N_t}$ corresponds to a covered Lorentzian pre-length space $(X, \ell, o, \mathcal{U})$ (Mondino-Sämann Def 3.8) on $S^3 \times \mathbb{R}_t$, with $o$ the BW vacuum event and $\mathcal{U}$ the modular-flow cover at scales $\beta_k = k$.

**A.4' — GeoVac wedge as Krein-pointed QMS (2–3 weeks).** Apply A.2' + A.3' to the specific Paper 43 / Paper 44 wedge construction. Verify the BW vacuum is the canonical pin state.

**A.5' — Phase A sprint memo (1 week).** Decision gate, same as original.

**Net Phase A re-scope: ~9–14 weeks total (vs. original 8 weeks).** Slightly longer due to the careful Krein-lift step. But the deliverable changes character substantially: instead of building pointed propinquity from scratch, Phase A delivers a *Lorentzian extension* of an already-published operator-algebraic framework. This is sharper, more publishable, and well-targeted.

### §6.2 Phase B (original Paper 48) is scooped — should NOT be written as scoped

**Paper 48 as scoped (pointed Latrémolière propinquity, math.OA standalone) would substantially duplicate Latrémolière 2512.03573.** The scoping memo §3 wrote: "the Phase B output (Paper 48) is the most original-NCG-math contribution of the program." That statement is now false.

**Alternative Phase B targets that ARE original:**
- **Paper 48' — Krein-lifted pointed propinquity (Lorentzian extension of Latrémolière 2512.03573).** The pointed-QMS framework in Krein signature. Novel: nobody has lifted Latrémolière's machinery to Krein.
- **Paper 48'' — The bridge theorem (Mondino-Sämann pLGH ↔ Latrémolière-Krein pointed QMS).** Novel: nobody bridges the two communities.
- **Or merged: a single paper combining the Krein lift + bridge + GeoVac G2-metric closure.** This would be Paper 48-49 merged, written as one work with the Lorentzian extension at its center.

The PI's prior judgment that Paper 48 should be split from Paper 49 was contingent on Paper 48 being "the pointed propinquity definition for the broader NCG community." With Latrémolière having already done that for the Riemannian case, the natural Paper 48 is now "the Lorentzian / Krein extension of Latrémolière 2512.03573, with GeoVac wedge as application." Single paper, ~25-30 pages.

### §6.3 Phase C (GeoVac G2-metric closure) remains unscooped and valid

The Sprint L3e-P3-C objective — pointed-propinquity convergence on the truncated Lorentzian Krein wedge to the genuine non-compact $S^3 \times \mathbb{R}_t$ — is unaffected by Latrémolière 2512.03573 (which is non-Lorentzian). Phase C remains a clean GeoVac contribution.

With the re-scoped Phase A and merged Phase B-48/C-49, the total program is now **~6-9 months** instead of 12. Two papers (Paper 48 = Krein lift + bridge + GeoVac application) or one paper.

### §6.4 Decision recommendation

**Continue Sprint L3e-P3** but with the re-scoped Phase A as described in §6.1. Do not write Paper 48 as originally scoped. The Phase A.5' decision gate (~10 weeks) should determine whether the Krein-lift is structurally clean (in which case proceed to merged Paper 48/49) or whether Krein structure obstructs the Latrémolière framework (in which case redesign at a deeper level).

**Alternative: pause and reassess.** The Latrémolière 2512.03573 publication is recent (Dec 2025) and the proof techniques are not yet fully digested by the community. There is a case for waiting 6–12 months to see what extensions Latrémolière (or others) publish before committing to a Krein-lift. This is the diagnostic-before-engineering option — and matches the audit's purpose.

**My recommendation (writing as audit author, not PI):** the re-scoped Phase A is contained at ~10 weeks with explicit decision gate. The pause-and-reassess option costs the 2026 calendar window. The original "60% Phase A success" probability now updates conditional on the Krein-lift specifically — I'd estimate 50–70% for clean Krein-lift of Latrémolière 2512.03573, with the 30–50% obstruction case being a useful negative result (Krein-incompatibility of Latrémolière's Leibniz hermitian norm). Either outcome is publishable.

---

## §7. Concurrent-work risk assessment

### §7.1 Re-priced scoop matrix (May 2026)

| Risk | Source | Current state | Recommendation |
|:-----|:-------|:--------------|:---------------|
| Independent Lorentzian pointed propinquity | Latrémolière | Already published 2512.03573 (Riemannian only) — Lorentzian extension would require either Latrémolière himself or a collaborator | Pre-submit Paper 48 ASAP if proceeding |
| Mondino-Sämann move into operator algebras | Mondino-Sämann | No signal in 2024–2026 trajectory (entirely synthetic) | LOW risk; revisit at A.5' decision gate |
| Independent Krein-lift of Latrémolière 2512.03573 | Anyone in NCG | High potential interest given Dec 2025 publication; ~6–12 months for an active group to spot the gap | MEDIUM risk; pre-submit ASAP |
| Bridge synthetic ↔ operator-algebraic | Anyone in NCG metric geometry | Open and visible; ~12–24 months for someone to write | LOW-MEDIUM risk |
| Hekkelman-McDonald extension to Lorentzian | Hekkelman, McDonald | Different framework (Tauberian); MEDIUM risk in their Tauberian direction but not the metric direction | LOW risk for metric propinquity |
| Ponge or McDonald-Sukochev-Zanin lineage extension to Krein | Ponge / MSZ lineage | Non-compact spectral triples already addressed (cusp manifolds), but at integral-theoretic level | LOW risk for metric propinquity |

### §7.2 Specific scoop-vector analysis

**The most concerning vector is "independent Krein-lift of Latrémolière 2512.03573."** Latrémolière's December 2025 paper is timely and visible. Active NCG groups (Hawkins, Connes-vS lineage, Marcolli-vS lineage, Sitarz, Krajewski, Connes-Suijlekom) could spot the Lorentzian extension as a natural next step. The 6–12 month timeline is plausible if anyone is currently working in that direction.

**Mitigation:** if Sprint L3e-P3 proceeds with the re-scoped Phase A.2'–A.4', the Phase A deliverables should be written as arXiv-deposit-ready memos at each step. The Phase A.5' decision gate at ~10 weeks would coincide with arXiv-submission of a Phase A report. This minimizes window-of-exposure.

### §7.3 What about scoop on the GeoVac-specific Phase C?

**The GeoVac-specific Phase C is harder to scoop** because it requires:
- Knowledge of GeoVac's Camporesi-Higuchi spectral triple (Paper 32 §III)
- The four-witness Wick-rotation closure (Paper 42)
- The Lorentzian Krein wedge construction (Paper 43)
- The pointed-propinquity machinery (now Latrémolière 2512.03573)

The first three are GeoVac-internal. An outside group could in principle build the analogous Lorentzian / Krein wedge construction on flat Minkowski space, but it would not be GeoVac's $S^3 \times \mathbb{R}_t$ pointed-propinquity closure. **Phase C is essentially unscoopable in its GeoVac-specific form.** This argues for moving directly to the merged Paper 48/49 (Krein lift + GeoVac application) once Phase A returns positive.

---

## §8. Verdict and recommendation

**The Phase A.2-A.4 bridge construction is no longer well-targeted as originally scoped, but a re-scoped Phase A targeting the Lorentzian extension of Latrémolière 2512.03573 IS well-targeted.**

The original scoping memo's framing — "Phase B Paper 48 builds pointed propinquity from scratch, Phase C applies it to GeoVac" — has been displaced by Latrémolière 2512.03573 (December 3, 2025), which builds pointed propinquity from scratch on the Riemannian side. The remaining novelty is:
1. Lift Latrémolière's framework to Krein signature (~3–6 weeks; Phase A.2');
2. Bridge to Mondino-Sämann pLGH (~3–4 weeks; Phase A.3');
3. Apply to GeoVac G2-metric closure (~3–4 months; Phase B'/C', combined).

The total scope shrinks from 12 months to 6–9 months. The deliverable changes from two standalone papers (Paper 48 + Paper 49) to one or two papers (a merged Krein-lift / GeoVac application, possibly with Paper 49 split off as the bound-state QFT applications follow-up).

**Probability of Phase A success** (re-scoped, focused on Krein-lift): I estimate 50–70%. The Krein structure could either accommodate Latrémolière's Leibniz hermitian norm cleanly (positive outcome) or surface a structural incompatibility (negative outcome — also publishable as a documented obstruction). Either way, ~10 weeks to a decision gate.

**Recommended next action:** confirm with PI whether to proceed with the re-scoped Phase A (10 weeks, leading to a merged-paper decision at A.5') or to pause-and-reassess for 6 months to let the Latrémolière 2512.03573 framework settle in the literature. My audit-author recommendation is to proceed with the re-scoped Phase A: the contained 10-week sprint is the same diagnostic-before-engineering pattern that succeeded for earlier multi-month commitments, and the deliverable (Krein-lift attempt) is publishable regardless of outcome.

**Open follow-up not in this audit's scope:** the Phase C bound-state QFT applications (Hawking, Unruh, Lamb shift in curved background) are unaffected by the scoop and remain valid GeoVac downstream sprints once the merged Paper 48' lands.

---

End of memo.

## Sources

Synthetic Lorentzian GH (verified May 2026):
- arXiv:2209.14384 — Müller et al., "Lorentzian metric spaces and their Gromov-Hausdorff convergence" (2022, v3 2025)
- arXiv:2209.12736 — "GH metrics and dimensions of Lorentzian length spaces" (2022)
- arXiv:2412.04311 — Bykov, Minguzzi, Suhr, "Lorentzian metric spaces & GH-convergence: unbounded case" (Dec 2024)
- arXiv:2504.10380 v4 — Mondino, Sämann, "Lorentzian GH convergence and pre-compactness" (April 2025, v4 Dec 9, 2025)
- arXiv:2506.10852 — Braun, Sämann, "Gromov's reconstruction theorem & measured Lorentzian GH" (June 2025)
- arXiv:2510.13069 — Che, Perales, Sormani, "Synthetic Lorentzian GH convergence / Gromov compactness for intrinsic timed-Hausdorff distance" (Oct 2025)
- arXiv:2510.24423 — Minguzzi, "Results on Lorentzian metric spaces," Gen. Rel. Grav. 58 (2026)
- arXiv:2511.18389 — "Intrinsic Timed Hausdorff Convergence and Its Implications" (Nov 2025)
- arXiv:2605.09101 — Kubota, "Lorentzian coarea inequality" (May 2026)
- arXiv:2605.11271 — "Convergence of Lorentzian spaces and curvature bounds for generalized cones" (May 2026)
- arXiv:2605.03172 — "Stability of synthetic timelike Ricci bounds under C⁰-limits" (May 2026)

Latrémolière propinquity (verified May 2026):
- arXiv:1302.4058 — Latrémolière 2013, Quantum GH Propinquity
- arXiv:1406.0233 — Latrémolière 2014, Topographic GH Hypertopology for Proper QMS (predecessor of 2512.03573)
- arXiv:1811.10843 — Latrémolière 2018, GH propinquity for metric spectral triples (Adv. Math. 2022)
- arXiv:2404.00240 — Farsi-Latrémolière, Collapse in NCG (v3 Oct 2025)
- arXiv:2504.11715 — Farsi-Latrémolière, Continuity on Riemannian path (April 2025)
- **arXiv:2512.03573 — Latrémolière, "The quantum Gromov-Hausdorff Hypertopology on the class of pointed Proper Quantum Metric Spaces" (Dec 3, 2025) — KEY FINDING**
- arXiv:2602.23080 — "Noncommutative coarse metric geometry" (early 2026)
- arXiv:2604.02117 — "The Bures metric and the quantum metric on the density space of a C*-algebra: the non-unital case" (early 2026)

Hekkelman-McDonald / Tauberian non-compact (verified May 2026):
- arXiv:2412.00628 — Hekkelman, McDonald, "Noncommutative integral on spectrally truncated spectral triples, link with quantum ergodicity" (Dec 2024, accepted J. Funct. Anal. July 2025)
- arXiv:2106.02235 — "Semiclassical Weyl law and exact spectral asymptotics in NCG" (2021)
- arXiv:2604.15008 — Ponge, "Noncommutative Geometry, Spectral Asymptotics, and Semiclassical Analysis" (2026 monograph)

Bridging context:
- arXiv:2005.08544 — "Gromov-Hausdorff convergence of state spaces for spectral truncations" (2020)
- arXiv:2302.09117 — Bassi, "Isometry groups of inductive limits of metric spectral triples and GH convergence" (J. London Math. Soc. 2023)
- arXiv:2410.15454 — "UCP maps and GH convergence"
- arXiv:1505.01939 — van den Dungen, Krein spectral triples (2016)
- arXiv:1611.07062 — Bizi-Brouder-Besnard (m,n) classification (2018)

Pre-existing baseline:
- `debug/l3_literature_audit_memo.md` (May 17, 2026 baseline)
- `debug/sprint_l3e_p3_detailed_scoping_memo.md` (full scoping)
