# Multi-Focal Phase B Sub-Sprint W2b-diag

**Sprint:** Phase B sub-sprint **W2b-diag** of the multi-focal-composition Phase B family.
**Date:** 2026-05-07.
**Author:** Sub-agent (W2b-diag).
**Frame:** Phase A's wall taxonomy named **W2b** as the cross-manifold spectral-triple incompatibility wall — exhibited in G4b (Sprint G4 scoping memo `debug/g4_cross_manifold_scoping_memo.md`; Paper 32 §VIII.C; Paper 24 §V's four-layer Coulomb/HO asymmetry). The Phase A literature review (Track 3, surprise S2) found that tensor products of two infinite-dimensional metric spectral triples are openly stated as open in the published NCG literature: Latrémolière 2026 (arXiv:2603.19128) closes the *almost-commutative* case (one infinite Riemannian + one finite); Farsi–Latrémolière 2024/2025 work on related tensor-product structures but do not appear to close the two-infinite-metric case.
**Status:** Diagnostic. Literature read + structural argument; no production code modified, no papers modified.

---

## §1. The diagnostic question

W2b is named as "tensor product of two infinite metric spectral triples" — a category Phase A flagged as a **frontier-of-NCG** problem rather than a GeoVac-internal limitation. The natural reduction sequence is:

- **W2b-easy:** $\mathcal{T}_{S^3}^{\lambda_a} \otimes \mathcal{T}_{S^3}^{\lambda_b}$ — same manifold (round $S^3$, Camporesi–Higuchi spinor Dirac), two different focal lengths $\lambda_a \neq \lambda_b$. This is what GeoVac actually wants for **two electrons at different focal lengths**, i.e. it is W1a (cross-register coordinate operator) re-expressed in the spectral-triple language. Each factor *individually* is the object Paper 38 closed (WH1 PROVEN); the question is whether their tensor product converges in the propinquity to the tensor product of continuum triples.
- **W2b-medium:** $\mathcal{T}_{S^3} \otimes \mathcal{T}_{\mathrm{Hardy}(S^5)}$ — Camporesi–Higuchi Riemannian spinor on one side, Bargmann–Segal first-order complex Euler operator on the Hardy sector $H^2(S^5)$ on the other. Different manifolds, different categories of operator, but each is a "natural" GeoVac construction (Paper 7 on the Coulomb side, Paper 24 on the HO side).
- **W2b-hard:** $\mathcal{T}_{S^3} \otimes \mathcal{T}_{S^5}^{\mathrm{Riem}}$ where the second factor is the round Riemannian spectral triple on $S^5$, abandoning the Bargmann–Segal content of Paper 24. KO-arithmetic clean ($3 + 5 = 8 \equiv 0 \pmod 8$), textbook Connes–Marcolli object — but throws away the Sprint ST-SU3 SU(3) Wilson construction on the Bargmann graph (because Sprint ST-SU3 lives on the Bargmann graph, not on round $S^5$).

The distinction matters strategically. If W2b-easy is reachable by extending Paper 38's machinery (Latrémolière propinquity on $\mathrm{SU}(2)$ via Peter–Weyl + central spectral Fejér + Berezin reconstruction), that closes W1a from the NCG side rather than the atomic-physics side. The multi-focal-composition sprint then has a **single keystone target** — extending Paper 38 to a tensor product of two copies of itself — rather than a panel of seven independent walls.

The diagnostic asks: (Q1) is two-infinite-metric tensor product genuinely open in published NCG; (Q2) which sub-question is the right Phase C-W2b target; (Q3) does the Coulomb/HO category mismatch obstruct W2b-medium structurally; (Q4) what is the cleanest Phase C-W2b sprint conditional on Q1–Q3; (Q5) verdict on tractability.

---

## §2. Q1: Literature confirmation that two-infinite-metric tensor product is open

### §2.1 Direct reads

I confirmed Track 3's S2 surprise via direct WebFetch on the four most relevant published papers:

**Latrémolière 2026 (arXiv:2603.19128), "Spectral continuity of almost commutative manifolds for the C¹ topology on Riemannian metrics."** The abstract explicitly states the models studied are "constructed as products of [the] canonical spectral triple of a compact connected spin manifold with a finite dimensional spectral triple." The novel contribution is using the spectral propinquity to establish continuity of these products as the Riemannian metric varies. **One Riemannian × one finite. The two-infinite case is not addressed.**

**Farsi–Latrémolière 2024 (arXiv:2404.00240), "Collapse in noncommutative geometry and spectral continuity."** The abstract names the worked examples as "collapse of product of spectral triples with one Abelian factor, U(1) principal bundles over Riemannian spin manifolds, and noncommutative principal bundles." The phrase **"with one Abelian factor"** is the operative restriction. Two non-abelian infinite factors is not in the worked examples and (per the public abstract) not in the theorems.

**Farsi–Latrémolière 2025 (arXiv:2504.11715), "Continuity for the spectral propinquity of the Dirac operators associated with an analytic path of Riemannian metrics."** The abstract explicitly states this is a *single* metric spectral triple varying along an analytic path of Riemannian metrics on a fixed closed spin manifold. **Single triple, not a product.** A continuous-family setup distinct from the discrete-truncation tensor product W2b-easy needs.

**Aguilar 2019 (arXiv:1907.07357), "Quantum metrics on $C(X) \otimes A$ for $A$ an AF algebra."** The structural restriction is sharper than Track 3 quoted: **one factor commutative ($C(X)$), the other AF.** This is "one Riemannian-with-no-Dirac × one finite-step-limit," strictly inside the almost-commutative Latrémolière 2026 paradigm with the finite factor relaxed to AF.

Adjacent: **Hekkelman–McDonald 2024b (arXiv:2412.00628)** does not address tensor products. Its scope is the noncommutative integral on a *single* spectrally truncated spectral triple, with applications to quantum ergodicity and Szegő-style asymptotics. Nothing about products.

### §2.2 What the published thread covers

Synthesizing across the four papers and Track 3's broader review, the published propinquity tensor-product theorems span:

| Setup | Status | Reference |
|:------|:------:|:----------|
| Two finite-dimensional triples | Standard | Connes–Chamseddine almost-commutative |
| Riemannian × finite | Closed | Latrémolière 2026 (arXiv:2603.19128) |
| $C(X)$ × AF | Closed | Aguilar 2019 (arXiv:1907.07357) |
| Riemannian × abelian | Worked examples in Farsi–Latrémolière 2024 | arXiv:2404.00240 |
| U(1) principal bundles over Riemannian spin | Worked examples in Farsi–Latrémolière 2024 | arXiv:2404.00240 |
| Single Riemannian along analytic metric path | Closed | Farsi–Latrémolière 2025 (arXiv:2504.11715) |
| Two infinite Riemannian non-flat triples | **OPEN** | None |
| Riemannian × Hardy-sector first-order complex | **OPEN** | None |
| Two GeoVac spectral triples at different focal lengths | **OPEN** | None |

Every existing tensor-product propinquity theorem requires the second factor to be either (i) finite, (ii) AF, (iii) abelian, or (iv) reduce to a principal-bundle structure that effectively restricts the second factor to a single character / Hopf $S^1$. **The two-infinite-non-abelian case is not in the published literature**. Track 3's S2 surprise is robustly correct.

This is consistent with the Connes–van Suijlekom 2021 paradigm: their CMP 383 paper sets up truncations on a single spectral triple and explicitly defers convergence to "elsewhere" three times. The flat-structure subset (Hekkelman 2022 on $S^1$, Hekkelman–McDonald 2024 on $T^d$, Leimbach–van Suijlekom 2024 on tori) closed the abelian-flat cases. Paper 38 (Loutey 2026) closed the simplest non-abelian case ($\mathrm{SU}(2) = S^3$). **No subsequent paper closes a tensor product of two such non-abelian objects.**

### §2.3 Why the gap exists structurally

Three technical issues compound: (a) Latrémolière tunneling pairs do not multiply trivially — for $(B_a \otimes B_b, P_a \otimes P_b)$ the joint reach inherits from both L4 rates with their commutator structure, and the joint Lipschitz constant is not $\max(C_3^a, C_3^b)$ in general (Bożejko–Fendler 1991 covers each factor's cb-norm but the propinquity-tunneling-pair assembly for the product is missing in print); (b) the Lipschitz seminorm on a tensor product couples the two sides through the chirality grading $\gamma_a$ via Connes–Marcolli, and the joint Lipschitz comparison estimate is unpublished in the two-infinite case; (c) Berezin reconstruction on $S^3$ is non-Kähler-specific (Paper 38 uses central-spectral-Fejér rather than Hawkins holomorphic Toeplitz), and tensoring two such reconstructions is published only for abelian-on-abelian (Leimbach–van Suijlekom torus). None of these is unsurmountable; none is published.

### §2.4 Verdict on Q1

**Confirmed.** Track 3 surprise S2 is correct. Tensor products of two infinite-dimensional non-abelian metric spectral triples are an open question in the published NCG literature as of May 2026. The closest published result on the structural neighborhood is Paper 38 (Loutey 2026) for the *single*-factor $\mathrm{SU}(2)$ case, which is the natural prerequisite. The next published result in the sequence that NCG-community-naturally addresses W2b-easy (two SU(2) factors) does not yet exist.

---

## §3. Q2: Sub-question identification (W2b-easy, medium, hard)

### §3.1 W2b-easy: $\mathcal{T}_{S^3}^{\lambda_a} \otimes \mathcal{T}_{S^3}^{\lambda_b}$

This is the **simplest natural extension of Paper 38**: same manifold, same KO-dim 3, same Camporesi–Higuchi Dirac spectrum on each factor, only the focal length parameter $\lambda$ differs. The factors are individually the object Paper 38 closed.

Setup: Two electrons in distinct hydrogenic shells, each living on its own $S^3$ Fock projection at its own $p_0 = Z/n$. The composed Hilbert space is $\mathcal{H}_a \otimes \mathcal{H}_b$, with shell-truncated $\mathcal{H}_{a,n_{\max}^a} \otimes \mathcal{H}_{b,n_{\max}^b}$. The composed Dirac is $D_{a,b} = D_a \otimes 1_b + \gamma_a \otimes D_b$ in the Connes–Marcolli convention (KO-dim $3 + 3 = 6$), or some refined recoupling if one wants more than just direct product. The composed algebra is $\mathcal{A}_a \otimes \mathcal{A}_b = C^\infty(S^3) \otimes C^\infty(S^3)$.

**Why this is the right reduction of W1a.** The HF-3 recoil sprint diagnosed W1a as "the framework has no two-body spatial coordinate operator $V(\hat{r}_e, \hat{R}_n)$ across registers at different focal lengths." On the spectral-triple side, that two-body operator lives in $\mathcal{A}_e \otimes \mathcal{A}_N$ acting on $\mathcal{H}_e \otimes \mathcal{H}_N$. The cross-coupling — the spatial multipole-expanded $V_{eN}$ — is a multiplier on the tensor product algebra. **W2b-easy is W1a re-expressed in the spectral-triple language.** Closing it from the propinquity side does not require building a new register architecture; it requires extending Paper 38's GH-convergence theorem from one factor to a tensor product of two factors of the same type.

**Difficulty class.** Each factor is the object Paper 38 closed. The five-lemma chain extends factor-by-factor: L1' lifts via $\gamma_a \otimes \gamma_b$ chirality grading; L2's central spectral Fejér tensorizes as $K_a(g_a) K_b(g_b)$ (product Plancherel) with cb-norm $4/((n_{\max}^a + 1)(n_{\max}^b + 1))$; L3's joint Lipschitz follows from Connes–Marcolli, $\|[D_{a,b}, M_{f \otimes g}]\| \le \|[D_a, M_f]\| \|g\|_\infty + \|f\|_\infty \|[D_b, M_g]\|$, giving joint $C_3 \le 2$ in general and $C_3 = 1$ on the normalized factorized panel (sharpening to the full operator system requires the dual-triangle-inequality argument of Paper 38 §3.4 lifted to the product); L4's Berezin reconstruction tensorizes as $B_a \otimes B_b$ with all four properties inherited factor-by-factor and approximate-identity rate $\max(\gamma_{n_{\max}^a}, \gamma_{n_{\max}^b})$; L5's tunneling pair $(B_a \otimes B_b, P_a \otimes P_b)$ has joint reach $\le \max(\gamma_{n_{\max}^a}, \gamma_{n_{\max}^b})$, with the height bookkeeping reducing via Paper 38 Remark 3.5 plus L4(d). **Net headline:** $\Lambda(\mathcal{T}_a \otimes \mathcal{T}_b, \mathcal{T}_{S^3} \otimes \mathcal{T}_{S^3}) \le C_3 \max(\gamma_{n_{\max}^a}, \gamma_{n_{\max}^b})$, structurally a 4–8 week extension of Paper 38, not a multi-month frontier.

### §3.2 W2b-medium and W2b-hard

**W2b-medium** ($\mathcal{T}_{S^3} \otimes \mathcal{T}_{\mathrm{Hardy}(S^5)}$): different manifolds, different categories of operator — Riemannian spinor Dirac on $S^3$ vs. first-order complex Euler on $H^2(S^5)$. KO-dim is not classically defined for the Hardy factor. This is the G4b sub-gap (Paper 32 §VIII.C; G4 scoping memo). The Coulomb/HO category mismatch (Paper 24 §V) **structurally obstructs** W2b-medium: no published prescription in NCG for tensoring Riemannian × Hardy-sector first-order complex. Detailed analysis in §4.

**W2b-hard** ($\mathcal{T}_{S^3} \otimes \mathcal{T}_{S^5}^{\mathrm{Riem}}$): abandon Paper 24 Bargmann content, take round-$S^5$ Riemannian spinor (KO-dim 5). KO-arithmetic clean ($3+5=8\equiv 0 \pmod 8$). But (i) the propinquity convergence for $\mathcal{T}_{S^5}^{\mathrm{Riem}}$ alone is unpublished (Paper 38 covers $S^3 = \mathrm{SU}(2)$; rank-2 Spin(6) analog requires fresh L2/L3/L4 work); (ii) joint tensor propinquity is then the W2b-easy difficulty class on top; (iii) throws away the Sprint ST-SU3 SU(3) Wilson content. **Lower structural payoff than W2b-easy** despite being mathematically cleanest.

### §3.4 Verdict on Q2

**The natural Phase C-W2b target is W2b-easy.** The difficulty class is in the same neighborhood as Paper 38 (4–8 weeks of focused extension; possibly sub-sprint scale), and the strategic payoff is dramatic: closing W2b-easy closes W1a from the NCG side, eliminating the need for a separate atomic-physics-style sprint on Pachucki–Patkóš–Yerokhin reduction. **W2b-medium is structurally obstructed (Coulomb/HO category mismatch), and W2b-hard requires building $\mathcal{T}_{S^5}^{\mathrm{Riem}}$ from scratch as a prerequisite.**

---

## §4. Q3: Coulomb/HO category mismatch — is W2b-medium resolvable?

### §4.1 The mismatch in technical detail

Paper 24 §V states the four-layer Coulomb/HO asymmetry: (i) spectrum-computing role of $L_0$ (yes $S^3$, no $S^5$ Bargmann); (ii) calibration $\pi$ (yes $S^3$, no $S^5$); (iii) Wilson gauge with natural matter coupling (yes $S^3$, no $S^5$ — matter on Bargmann breaks the fixed-rep requirement, Sprint 5 Track S5 negative); (iv) universal Wilson–Hodge combinatorial vocabulary. Layers (i)–(iii) are Coulomb-specific; (iv) is universal.

At the spectral-triple level (Paper 24 Theorem 6.1 of §V; Paper 32 §VIII.C):

- **Coulomb side:** $D_{S^3}$ is a second-order Riemannian spinor Dirac. KO-dim 3. Spectrum $|\lambda_n| = n + 1/2$ via Camporesi–Higuchi (1996). Spinor bundle $\Sigma$ over $S^3$ in the standard Riemannian sense.
- **HO side (Bargmann):** the natural diagonal operator is the Euler / number operator $\hat N = \sum_i z_i \partial_i$ on $H^2(S^5)$ restricted to symmetric SU(3) irreps $(N, 0)$. **First-order complex-analytic.** Linear spectrum $\lambda_N = N$ on $\mathcal{H}_N$ of dimension $(N+1)(N+2)/2$. KO-dimension is **not classically defined** (Hardy sector, not Riemannian spinor bundle).

The mismatch is at the *operator type* level, not at the manifold level. It would persist even if one took $\mathcal{T}_{S^5}^{\mathrm{Riem}}$ instead of the Bargmann construction — but in that case the Bargmann content is gone, which is W2b-hard, not W2b-medium.

### §4.2 Bargmann transform as a candidate bridge?

The Bargmann transform $B: L^2(\mathbb{R}^3) \to F^2(\mathbb{C}^3)$ is unitary (Bargmann–Segal–Hall) and maps the Schrödinger HO Hamiltonian to $\hbar\omega(\hat N + 3/2)$. **This is a unitary intertwiner of operators on Hilbert spaces, not a unitary equivalence of spectral triples.** A spectral triple specifies an algebra, a Dirac, and (for real / KO-dim) a real structure $J$ and chirality grading $\gamma$. Bargmann maps the Schrödinger picture (no canonical Dirac) to the Bargmann picture (Euler, no $\gamma$ or $J$); neither side is yet a spectral triple. Furthermore, the Coulomb $\mathcal{T}_{S^3}$ uses the Camporesi–Higuchi Riemannian spinor Dirac on $S^3$, not the HO side at all — Bargmann is for the HO side. **The Bargmann transform is not a bridge between $\mathcal{T}_{S^3}$ and $\mathcal{T}_{\mathrm{Hardy}(S^5)}$.** The Kustaanheimo–Stiefel transformation maps 3D Coulomb to 4D HO (with constraint), but (a) at coordinate level only, (b) to the *4D* HO not the 3D HO of Paper 24, (c) without preserving the Riemannian / spinor bundle structure required for KO-dim. KS does not bridge the spectral triples either.

### §4.3 Is the mismatch resolvable in principle?

A genuine bridge would require either:

(a) Generalizing the spectral-triple framework to handle Hardy-sector first-order complex operators alongside Riemannian spinor Dirac operators in the same tensor product. **No published work does this.** Hawkins 2000 (Berezin–Toeplitz on the unit disk) gets close on a single manifold but does not cross.

(b) Reformulating the Bargmann–Segal $S^5$ side as a Riemannian spectral triple by promoting the Euler operator to a "Dirac-like" object via chiral doubling (Option B in Sprint G4 scoping). This sacrifices the bit-exact $\pi$-free certificate of Paper 24 §III on the $S^5$ side and produces a non-standard NCG object (KO-dim not classically defined).

(c) Identifying a unifying ambient manifold whose GeoVac sub-structure contains both $S^3$ Coulomb and $S^5$ Bargmann as natural sub-objects. **No candidate is in the GeoVac record.** Paper 24 §V explicitly flags this as an open structural question.

None of (a)–(c) is reachable as a sprint-scale extension of Paper 38. **The Coulomb/HO mismatch obstructs W2b-medium at the framework level, not at the calculational level.**

### §4.4 Verdict on Q3

**W2b-medium is structurally blocked by category mismatch and is not resolvable via any known bridge (Bargmann, KS, etc.).** Closing it requires NCG-framework extension — either generalized propinquity theory for Hardy-sector × Riemannian factors, or a unifying ambient-manifold construction. Both are multi-paper / multi-year directions. **The recommendation is to record W2b-medium as an open structural question in Paper 32 §VIII.C / Paper 24 §V (already done in both venues) and not to pursue it as a sprint target.**

This sharpens Phase A's W2b row. W2b is genuinely two distinct walls: W2b-easy (extensible from Paper 38, sprint-scale) vs. W2b-medium/hard (NCG-framework-level, multi-paper). The Phase A summary collapsed them; Phase B-W2b-diag separates them.

---

## §5. Q4: Cleanest Phase C-W2b sprint conditional on Q1–Q3

### §5.1 The recommended sprint: Phase C-W2b-easy

**Title.** Latrémolière propinquity convergence for a tensor product of two truncated Camporesi–Higuchi spectral triples on $S^3$ at distinct focal lengths.

**Theorem statement (target).** Let $\mathcal{T}_{S^3}^a, \mathcal{T}_{S^3}^b$ denote two copies of the Camporesi–Higuchi spectral triple at focal lengths $\lambda_a, \lambda_b > 0$ (distinct). Let $\mathcal{T}_{n_{\max}^a, n_{\max}^b} = \mathcal{T}_{n_{\max}^a}^a \otimes \mathcal{T}_{n_{\max}^b}^b$ denote the tensor-product truncated metric spectral triple. Then $\mathcal{T}_{n_{\max}^a, n_{\max}^b}$ converges to $\mathcal{T}_{S^3}^a \otimes \mathcal{T}_{S^3}^b$ in the Latrémolière quantum Gromov–Hausdorff propinquity, with rate
$$\Lambda(\mathcal{T}_{n_{\max}^a, n_{\max}^b}, \mathcal{T}_{S^3}^a \otimes \mathcal{T}_{S^3}^b) \le C_3 \cdot \max(\gamma_{n_{\max}^a}, \gamma_{n_{\max}^b})$$
with $C_3 = 1$ on the factorized observable panel and a sharpened constant on the full operator system, and asymptotic constant $4/\pi$ on each factor (the M1 Hopf-base measure signature, inherited from Paper 38).

**Roadmap.** (1) Tensor-extend each of Paper 38's five lemmas factor-by-factor (sketched in §3.1). (2) Identify the joint L3 Lipschitz constant: $C_3 \le 1$ on the normalized factorized panel, $\le 2$ in general; $C_3 = 1$ on the full operator system requires a dual-triangle-inequality lift of Paper 38 §3.4 to the product. (3) Joint Berezin reconstruction $B_a \otimes B_b$ inherits all four L4 properties factor-by-factor with rate $\max(\gamma_{n_{\max}^a}, \gamma_{n_{\max}^b})$. (4) Assemble the propinquity bound: $\Lambda \le C_3 \cdot \max(\gamma_{n_{\max}^a}, \gamma_{n_{\max}^b})$ via L5 tunneling-pair + L4(d) + Stein–Weiss. (5) Numerical verification at $(n_{\max}^a, n_{\max}^b) \in \{2,3,4\}^2$ via existing `geovac/gh_convergence.py` and `geovac/connes_distance.py` (doubles dim_H but stays tractable to $(3,3)$). (6) Paper draft, either standalone (Paper 38b) or §VI of Paper 38 v2.

**Effort.** 4–6 weeks for the symbolic mathematics, 2–4 weeks for the numerical verification, 2 weeks for the paper draft. Total 8–12 weeks at sprint cadence. **Risk profile: low-medium.** L1', L2, L4 extensions are mechanical; risk is concentrated in joint L3 (Bożejko–Fendler product cb-norm) and joint L5 (Latrémolière 2017/2023 propinquity-definition height bookkeeping), both addressable with existing techniques.

### §5.2 What this sprint produces (and does not)

**Produces:** a standalone NCG theorem closing a genuinely-open published question (Q1); algebraic-side closure of W1a (cross-register coordinate operator definable as a multiplier on the joint algebra; Pachucki–Patkóš–Yerokhin port becomes a streamlined application instead of an independent construction); a second WH1-PROVEN-class headline in the project record; a foundational object for Phase C-W1a.

**Does NOT produce:** closure of W2b-medium (Coulomb/HO category mismatch, §4); closure of W2b-hard (would need $\mathcal{T}_{S^5}^{\mathrm{Riem}}$ propinquity as prerequisite); progress on W2a, W3, W1b, W1c (independent walls). Those remain on the panel.

### §5.3 Verdict on Q4

**Phase C-W2b-easy is the cleanest target.** It is the single sprint that closes one of the multi-focal walls definitively, produces a publishable NCG theorem, and unlocks W1a as an application. W2b-medium and W2b-hard are open structural questions that are not sprint-scale and should remain in Paper 32 §VIII.C as recorded open problems.

---

## §6. Q5: Verdict (a/b/c/d)

The four options as stated:

- **(a) W2b-easy is reachable by extending Paper 38's machinery and is the obvious next sprint after W1a-diag.**
- (b) W2b-easy is reachable but requires more than Paper 38's tools — frontier-of-field but tractable.
- (c) W2b in any form is at the NCG community's open frontier with no clear nearby reduction.
- (d) The cross-manifold $\mathcal{T}_{S^3} \otimes \mathcal{T}_{S^5}$ case is structurally blocked by the Coulomb/HO asymmetry's category mismatch and requires a fundamentally new bridge.

**Verdict: (a) modulo a structural caveat that simultaneously affirms (d).** W2b decomposes into one tractable case (W2b-easy, verdict (a) — 4–8 week extension of Paper 38's machinery, technical risk concentrated in joint L3 and L5, both addressable via Bożejko–Fendler / Latrémolière 2017/2023 results) and two structurally obstructed cases (W2b-medium, verdict (d) — Coulomb/HO category mismatch with no published bridge; W2b-hard, requires $\mathcal{T}_{S^5}^{\mathrm{Riem}}$ propinquity as a separate prerequisite paper). The right framing: W2b-easy is the keystone, the other two are open structural questions for the broader NCG community.

---

## §7. Strategic note: does closing W2b-easy eat W1a from the NCG side?

### §7.1 The strategic question

If Phase C-W2b-easy succeeds, what does it do to the rest of the multi-focal sprint plan? Specifically: does it close W1a (cross-register coordinate operator) without needing the separate Pachucki–Patkóš–Yerokhin port (Phase A Candidate 1)?

### §7.2 Argument that it does

W1a's failure mode is "no two-body spatial coordinate operator $V(\hat r_e, \hat R_N)$ across registers at different focal lengths." On the spectral-triple side, the cross-register architecture is a tensor product $\mathcal{H}_e \otimes \mathcal{H}_N$ (Paper 23 §VI Track NI). The cross-coupling is a multiplier on $\mathcal{A}_e \otimes \mathcal{A}_N$. If Phase C-W2b-easy proves that the tensor-product algebra is well-defined as a Lipschitz-seminorm tensor product converging to the continuum tensor product, then (i) cross-register multipliers are well-defined elements of the joint operator system, (ii) the Connes–Marcolli formula $\|[D_e \otimes 1 + \gamma_e \otimes D_N, M_{V_{eN}}]\|_{\mathrm{op}}$ bounds the cross-coupling Lipschitz seminorm, and (iii) the multipole expansion on $S^3 \times S^3$ (analog of Shibuya–Wulfman) terminates at $L_{\max}^a + L_{\max}^b$ by the joint Gaunt selection rule. This is an **algebraic-side closure of W1a**: the cross-register coordinate operator is constructed as a multiplier on the tensor product algebra, with the propinquity bound providing the metric foundation.

### §7.3 Argument that it does not (or only partially)

The Pachucki–Patkóš–Yerokhin 2023 target (PRL 130, 023004) is a Foldy–Wouthuysen reduction at $(Z\alpha)^6$ in mass ratio — a specific reduction at a specific order, mass-ratio recoil terms exact rather than perturbative. Closing W2b-easy gives the cross-register operator and the propinquity bound but not automatically the *physics extraction* at PRL precision. So the right framing: W2b-easy closes the **algebraic-side** of W1a; Pachucki port closes the **physics-side**. Sequencing: (1) Phase C-W2b-easy (4–8 weeks, NCG headline, closes W1a algebraically); (2) Phase C-W1a-physics (streamlined Pachucki port building on the W2b-easy framework, ~4–6 weeks instead of 6–10 because the algebraic foundation is laid). **Total: 8–14 weeks for full W1a closure, of which the first 4–8 weeks produce a publishable NCG theorem.**

### §7.4 Revised Phase A recommendation

W2b-easy is the keystone for the W1 cluster, not for the entire panel. W1b/c are tooling-addressable but independent of W2b-easy; W2a and W3 remain frontier-of-field. But the W1 cluster is the *largest* and *most-demanded* sub-cluster (Pachucki–Patkóš–Yerokhin 2023 is the strongest single calibration target in the Phase A audit), and **W2b-easy compresses it from a multi-sprint sequence to one sprint plus a streamlined application.**

The Phase A synthesis ranked Phase C-W1a (Pachucki port) as Candidate 1 with effort 6–10 weeks. Phase B-W2b-diag suggests the more strategic move is: **replace Candidate 1 with Phase C-W2b-easy at 4–8 weeks, followed by streamlined Phase C-W1a-physics at 4–6 weeks.** Total effort comparable, but the deliverable is now *two* publishable artifacts (NCG theorem + atomic-physics application) instead of one — and the NCG theorem itself closes one of the two-infinite-metric tensor-product cases that Track 3's S2 surprise documents as genuinely open. **This is the strongest single strategic recommendation from Phase B-W2b-diag.**

---

## §8. Honest scope and uncertainty

**Highly confident:** (a) tensor products of two infinite-dimensional non-abelian metric spectral triples are genuinely open in the published NCG literature (May 2026); WebFetch confirmed on arXiv:2603.19128, arXiv:2404.00240, arXiv:2504.11715, arXiv:1907.07357, arXiv:2412.00628. (b) The Coulomb/HO category mismatch (Paper 24 §V) is a real structural obstruction to W2b-medium. (c) W2b-easy is the right sub-question. (d) Each of Paper 38's five lemmas has a structurally clean factor-by-factor extension, with joint L3 Lipschitz and joint L5 height as the two technical bottlenecks.

**Less confident:** (a) Exact sprint timeline (4–8 weeks if clean, 12+ if joint L3/L5 needs unanticipated machinery). (b) The exact constant in joint L3 ($C_3 = 1$ on factorized panel is robust; on full operator system requires unverified dual-triangle-inequality argument; result might come out $C_3 \le 2$ or $C_3 \le 1 + o(1)$ — none affect asymptotic rate, only explicit constant). (c) Whether W2b-easy alone closes W1a algebraically, or whether an additional 2–4 weeks of cross-register-coordinate-operator construction is needed within the W2b-easy sprint.

**Not done:** Line-by-line verification of joint L3/L5 (sprint-scale, diagnostic only sketches); a check on 2025–2026 NCG conference proceedings beyond arXiv (small community; if an in-press result exists it would change Q1's verdict); a verification that `geovac/gh_convergence.py`, `geovac/connes_distance.py`, `geovac/operator_system.py`, `geovac/spinor_operator_system.py`, `geovac/berezin_reconstruction.py`, `geovac/central_fejer_su2.py`, `geovac/real_structure.py` extend mechanically to the tensor product (argued on structural grounds only).

**The diagnostic does NOT claim** Phase C-W2b-easy is guaranteed to succeed (high probability based on structural arguments, but rigorous proof is a sprint deliverable), nor that W2b-medium/hard are unreachable forever (they require NCG-framework extension, multi-paper / multi-year, not a sprint target).

---

## §9. Summary

- **Q1:** Open. No published NCG work covers two infinite-dimensional non-abelian metric spectral triples; Latrémolière 2026, Farsi–Latrémolière 2024/2025, Aguilar 2019, Hekkelman–McDonald 2024b all stay inside Riemannian × finite, $C(X)$ × AF, products-with-one-abelian-factor, or single-triple settings.
- **Q2:** W2b-easy ($\mathcal{T}_{S^3}^{\lambda_a} \otimes \mathcal{T}_{S^3}^{\lambda_b}$, same manifold, distinct focal lengths) is the right sub-question — sprint-scale 4–8 week extension of Paper 38, equivalent to closing W1a from the NCG side.
- **Q3:** W2b-medium structurally blocked by Coulomb/HO category mismatch (Paper 24 §V four-layer asymmetry, spectral-triple level). No published bridge between Riemannian spinor Dirac and Hardy-sector first-order complex operators. Bargmann transform and KS transformation do not bridge.
- **Q4:** Phase C-W2b-easy is the cleanest sprint target. 4–8 weeks; produces publishable NCG theorem; closes W1a algebraically.
- **Q5:** (a) modulo structural caveat from (d). W2b-easy reachable; W2b-medium open structural question.
- **Strategic:** yes, Phase C-W2b-easy eats W1a from the NCG side. The multi-focal-composition sprint may have a single keystone target — Phase C-W2b-easy — with W1a as its first application. **Make W2b-easy the keystone, not a sub-wall.**

---

**End of W2b-diag memo.**
