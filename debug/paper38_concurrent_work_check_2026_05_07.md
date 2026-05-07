# Paper 38 Concurrent-Work Freshness Check (2026-05-07)

**Scope:** post-2024-12 work in math.OA / math-ph / math.DG potentially overlapping Paper 38's claim — qualitative-rate Latrémolière propinquity convergence of finite spectral truncations $\mathcal{T}_{n_{\max}}$ to the round-$S^3$ Camporesi–Higuchi spectral triple $\mathcal{T}_{S^3}$, via a five-lemma proof using SU(2) Peter–Weyl harmonic analysis.

**Verdict: CITE_RECOMMENDED.** No direct competitor identified. Paper 38's specific result — propinquity convergence to a *non-abelian compact Lie group* spectral triple — remains novel as of 2026-05-07. Several near-overlaps in adjacent frameworks should be cited for completeness and to position Paper 38 against the surrounding literature; six specific recommendations below.

---

## 1. The specific concern

Has anyone proved or sketched GH/propinquity convergence of truncated spectral triples to a round-$S^3$ / SU(2) / non-abelian compact Lie group spectral triple in 2025–2026?

**Answer: No.** Every candidate identified either (a) treats abelian/flat structures only ($S^1$, $T^d$, $S^2$ Berezin–Toeplitz, abelian factors in collapse), (b) treats discrete groups with polynomial growth (explicitly excludes compact Lie groups), (c) treats inductive limits of AF algebras (different category), (d) treats analytic paths of metrics (continuous family, not discrete truncation), or (e) develops a different convergence object (UCP-map metric spaces, noncommutative integral, quantum ergodicity).

---

## 2. Direct competitors

**None found.** A focused search across:
- arXiv math.OA listings 2025 (sampled first 50 papers)
- Author searches: van Suijlekom, Latrémolière, Farsi, Hekkelman, McDonald, Leimbach, Aguilar, Marcolli, Connes, Pérez-Sánchez
- Keyword searches: "spectral truncation propinquity SU(2) sphere", "operator system Connes distance compact Lie group convergence", "Camporesi-Higuchi propinquity", "fuzzy sphere SU(2) spectral truncation convergence"
- Recent journal scans: Comm. Math. Phys., Adv. Math., J. Geom. Phys., J. Math. Phys., J. Funct. Anal.

returned no paper proving propinquity convergence of finite spectral truncations to a non-abelian compact Lie group spectral triple. The Connes–van Suijlekom 2021 framework (arXiv:2004.14115, the load-bearing reference for Paper 38) explicitly defers GH convergence to "elsewhere" three times, and to date, follow-ups have addressed only flat / polynomial-growth / inductive-limit categories.

---

## 3. Near-overlaps — citation recommendations

Six papers worth citing in Paper 38, ranked by relevance.

### 3.1 Bhattacharyya–Duhan–Pradhan, arXiv:2410.15454 (Oct 2024, rev. Feb 2025)
**"Gromov–Hausdorff convergence of metric spaces of UCP maps."**
Proves GH convergence of sequences of sets of unital completely positive maps under van Suijlekom's operator-system framework. The right *framework* (van Suijlekom truncations) but applied to UCP-map metric spaces with the BW-topology, not the Latrémolière metric-spectral-triple propinquity. Does not specify SU(2)/$S^3$ (abstract is general). **Cite** as the closest 2024–2025 framework analog; explicitly distinguish from Paper 38's claim, which works in the metric-spectral-triple propinquity (Latrémolière 2017/2023, arXiv:1811.10843).

### 3.2 Farsi–Latrémolière, arXiv:2504.11715 (April 2025)
**"Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics."**
Proves that a polynomial path of Riemannian metrics on a closed spin manifold induces a continuous field in the spectral propinquity of metric spectral triples. Different convergence setup (continuous family of metrics on a fixed manifold vs. Paper 38's discrete sequence of truncated triples converging to a continuous limit). **Cite** as the most recent (April 2025) result in the spectral propinquity framework; situate Paper 38 as the discrete-truncation analog.

### 3.3 Farsi–Latrémolière, arXiv:2404.00240 (March 2024, rev. October 2025)
**"Collapse in Noncommutative Geometry and Spectral Continuity."**
Proves spectral continuity for collapse of products of spectral triples *with at least one abelian factor* — covers U(1) principal bundles over Riemannian spin manifolds. Explicitly does **not** cover non-abelian like SU(2). **Cite** as a *scope-contrast* reference: the abelian-factor restriction is exactly what Paper 38 lifts.

### 3.4 Hekkelman–McDonald, arXiv:2412.00628 (Dec 2024, rev. Aug 2025; J. Funct. Anal., to appear, DOI: 10.1016/j.jfa.2025.111154)
**"A noncommutative integral on spectrally truncated spectral triples, and a link with quantum ergodicity."**
Develops the noncommutative integral and a Szegő limit formula for spectrally truncated spectral triples in the Connes–van Suijlekom paradigm; defines geodesic-flow ergodicity for compact spectral triples. Does *not* prove GH or propinquity convergence. **Cite** as a recent (2025) companion direction in the spectral-truncation program, complementary to Paper 38's metric-side result.

### 3.5 Toyota, arXiv:2309.13469 (Sep 2023)
**"Quantum Gromov–Hausdorff convergence of spectral truncations for groups with polynomial growth."**
Proves quantum-GH convergence of spectral truncations in the Lip-norm-from-high-derivatives framework, restricted to discrete groups with polynomial growth. **Compact Lie groups like SU(2) are not in scope** (SU(2) is compact, not polynomial-growth-discrete). **Cite** as a *scope-contrast* reference: Paper 38 occupies the orthogonal compact-Lie-group corner of the same general program.

### 3.6 Aguilar (recent, January 2026)
**"Spectral triples on Effros–Shen AF algebras."** (arXiv listing under K. Aguilar, exact identifier not retrieved in this scan.)
Constructs spectral triples on AF algebras using paths-category presentation (Christensen–Ivan style) and proves quantum-GH propinquity convergence for sequences of Effros–Shen and UHF algebras. Inductive-limit / AF setting; not non-abelian Lie group. **Optional cite** for completeness if Paper 38 wants to position against the broader 2025–2026 propinquity activity.

---

## 4. Already-cited references (sanity check)

Per the Paper 38 outreach plan, the load-bearing references that should already be in the bibliography:
- **Connes–van Suijlekom 2021** (arXiv:2004.14115, CMP) — operator-system spectral truncation framework. ✓ Should be in Paper 38.
- **Latrémolière 2017/2023** (arXiv:1811.10843, Adv. Math. 415, 108876) — metric-spectral-triple Gromov–Hausdorff propinquity. ✓ Pinned in pre-submission edit (Track 1, edit #1).
- **Leimbach–van Suijlekom 2024** (Adv. Math. 439, 109496; arXiv:2302.07877) — torus result that Paper 38 generalizes to SU(2). ✓ Should be in Paper 38.
- **Camporesi–Higuchi 1996** (J. Geom. Phys. 20, 1–18) — Dirac eigenfunctions on $S^3$. ✓ Standard reference.
- **Krajewski 1998 / Paschke–Sitarz 2000** — finite spectral triple classification (relevant to the AC-extension cousin work, less directly to Paper 38).

---

## 5. Other recent activity worth knowing about (low-priority)

- **Connes–Consani–Moscovici, arXiv:2511.22755 (Nov 2025), "Zeta Spectral Triples."** Active spectral-triple work; not a direct overlap but indicates the field remains live.
- **arXiv:2512.15450 (Dec 2025), "Emergence of Time from a Twisted Spectral Triple in Almost-Commutative Geometry."** Recent AC-style work.
- **arXiv:2502.18105 (May 2025), "Emergence of Lorentz symmetry."** Spectral-action / Yang–Mills–Higgs activity.

None of these is a competitor; they confirm the spectral-triple program is active in 2025–2026 but address different questions.

---

## 6. Recommendation

**Submission strategy unchanged.** Paper 38's claim — qualitative-rate Latrémolière propinquity convergence of finite truncations to the round-$S^3$ Camporesi–Higuchi spectral triple via SU(2) Peter–Weyl, with quantitative-rate $4/\pi$ asymptote on the L2 lemma — is novel as of 2026-05-07. No paper found between the December 2024 reference cutoff in the original literature survey and today proves the same theorem or a near-superset.

**Bibliography additions for Paper 38 before arXiv push:**
1. Bhattacharyya–Duhan–Pradhan 2024/2025 (UCP-map GH convergence) — distinguish from metric-spectral-triple propinquity.
2. Farsi–Latrémolière 2025 (analytic-path continuity) — situate against discrete truncation.
3. Farsi–Latrémolière 2024/2025 (collapse, abelian-factor restriction) — scope contrast.
4. Hekkelman–McDonald 2024/2025 (NC integral, quantum ergodicity) — companion direction.
5. Toyota 2023 (polynomial-growth groups) — scope contrast.
6. (Optional) Aguilar 2026 (Effros–Shen AF) — broader propinquity activity.

These are five-to-seven-line additions per bibitem, no rewrite of the proof or framing required. The paper's positioning as the first non-abelian compact-Lie-group case in the Connes–van Suijlekom program is preserved and sharpened by this round of citations.

**No re-check needed before submission unless held >2 weeks.**

---

## Sources searched

- arXiv abstracts: 2504.11715, 2412.00628, 2410.15454, 2309.13469, 2404.00240
- arXiv listing math.OA 2025 (first 50 papers)
- WebSearch queries on author names (Hekkelman, McDonald, van Suijlekom, Latrémolière, Aguilar, Pérez-Sánchez, Leimbach) and keywords (spectral triple propinquity SU(2), spectral truncation Gromov-Hausdorff sphere, Camporesi-Higuchi propinquity, fuzzy sphere SU(2) convergence)
- Recent journals via search (J. Geom. Phys., Comm. Math. Phys., Adv. Math., J. Funct. Anal., J. Math. Phys.)
