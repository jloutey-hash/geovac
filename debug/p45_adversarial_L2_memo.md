# Adversarial referee review — Paper 45 Lemma L2 (joint cb-norm) and Paper 38 antecedent

**Date:** 2026-06-09. **Mandate:** refute. **Scope:** Paper 45 §2.4 + §4.1 (Lemma L2) + Prop. reach_height; Paper 38 §L2/L4/L5; `debug/l3b_2_sub_sprint_B_cb_norm_memo.md`; `debug/r25_l2_proof_memo.md`; `geovac/central_fejer_su2.py`. Diagnostic only — no papers or production code modified.

**Verdict summary**

| Vector | Verdict |
|:--|:--|
| V1 object identification | **WRONG** (counterexample, bit-exact) |
| V2 citations (BF 1991, Pisier Ch. 8) | **WRONG citation** (BF), **UNVERIFIED-IMPORT/likely-wrong** (Pisier Thm 8.10); the needed fact is elementary |
| V3 downstream reach_P | **SERIOUS-GAP** (proof invalid; statement plausibly recoverable; rate not automatic) |
| V4 U(1) factor / N_t-independence | **FIXABLE-GAP** (displayed symbol ≠ displayed kernel; sup = 1 conclusion holds) |
| V5 Paper 38 antecedent | **WRONG** in L2(b),(c) and eq:B_def display; L4(a,c,d), L2(d), L3, reach_B, height_B **HOLD** under the convolution reading |
| Lemma L2 overall | **WRONG as stated** — proves a true fact about a different object; the true cb-norm of the map the proof architecture uses is exactly **1** |
| reach_P dependence on L2 | **SERIOUS-GAP** — the invoked invariant is irrelevant to a partial inverse; a corrected proof needs a Leimbach–vS-style antiderivative/transference lemma not present in either paper |

---

## V1. Object identification — WRONG

Three distinct objects are in play:

1. **Fourier/multiplier symbol of the kernel** $m_K(\pi) = \frac{1}{d_\pi}\int_G K(g)\overline{\chi_\pi(g)}\,dg$ — the scalar by which convolution $S_K f = K*f$ acts on the $\pi$-isotypic block. This is what a cb-norm computation of $S_K$ needs.
2. **Plancherel mass distribution of the kernel amplitude** $h = \sum_{j\le j_{\max}}\sqrt{2j+1}\,\chi_j$: the fraction $|\hat h(j)|^2/\|h\|_2^2 = (2j+1)/Z_{n_{\max}}$ of the $L^2$ mass of $h$ in shell $j$.
3. The Leimbach–vS **antiderivative Schur-multiplier norm** (their Lemma 3.4/3.5) — the genuinely small $O(1/n)$ quantity that drives the *inverse* direction in the torus proof.

**Papers 38/40/45 label object (2) as object (1).** Paper 38 L2(b) (paper_38...tex:524–528) claims $T_K$ acts on the weight-$j$ isotypic component by $(2j+1)/Z$, vanishing for $j > j_{\max}$; Paper 45 eq:K_SU2_plancherel (paper_45...tex:541–544) repeats it.

**Bit-exact counterexample** ($n_{\max}=2$, $Z=3$). $\chi_{1/2}^2 = \chi_0 + \chi_1$ gives
$K = |\chi_0 + \sqrt2\chi_{1/2}|^2/3 = \chi_0 + \tfrac{2\sqrt2}{3}\chi_{1/2} + \tfrac23\chi_1$, so the true symbol is

| $N=2j+1$ | true $m_K$ | claimed $N/Z$ |
|:--:|:--:|:--:|
| 1 | **1** | 1/3 |
| 2 | $\sqrt2/3 \approx 0.471$ | 2/3 |
| 3 | $2/9$ | 0 |

Verified two independent ways (exact Clebsch–Gordan expansion in sympy; numerical Weyl-integration quadrature on the maximal torus, agreement to 6+ digits at $n_{\max}\in\{2,3\}$). General pattern (checked $n_{\max}\le 4$ exact, $\le 80$ float): $m_K(1) = 1$ always; $m_K$ **decreasing**; support $N \le 2n_{\max}-1$, **not** $N \le n_{\max}$.

Consequences:

- **(a) contradicts (b) inside Lemma L2 itself.** $\widehat K(\text{trivial}) = \int_G K\,dg$ under any multiplier convention. Paper 38 L2(a) proves $\int K = 1$ (and `verify_normalization_symbolic` verifies it in sympy); L2(b) claims $\widehat K(0) = 1/Z \ne 1$. One of (a),(b) must be false; (a) is true.
- **The true cb-norm of $S_K$ is exactly 1**, for every $n_{\max}$, $N_t$: $K \ge 0$ with $\int K = 1$ makes $S_K$ (and the compressed Berezin $B = P\,M_{K*\cdot}\,P$) **unital completely positive**, and a UCP map has $\|\cdot\|_{\mathrm{cb}} = \|\varphi(1)\| = 1$. A claimed cb-norm $2/(n_{\max}+1) < 1$ is incompatible with unitality.
- **Paper 45 self-contradicts at line 1106–1108**: "both $B^{\mathrm{joint}}$ and $P^{\mathrm{joint}}$ are UCP (with the Berezin cb-norm bound $\|B\|_{\mathrm{cb}} \le 2/(n_{\max}+1) \le 1$…)". UCP ⟹ cb-norm $=1$; both halves of the sentence cannot hold for $n_{\max} \ge 2$. Same contradiction at L4(b) (lines 983–985, 1016–1024): the Young's-inequality proof gives constant 1; the "sharper" $2/(n_{\max}+1)$ bound would give $\|B(1)\| \le 2/(n_{\max}+1) < 1 = \|B(1)\|$.
- **The two factors of the same joint kernel use two different conventions.** The U(1) "Plancherel symbol" (eq:K_U1_plancherel) is multiplier-type (1 at trivial character, decreasing); the SU(2) one is mass-type (1/Z at trivial, increasing). Under the consistent multiplier convention the joint cb-norm is $1\cdot 1 = 1$; under the consistent mass convention the U(1) mass is uniform $1/N_t$ and the joint sup is $\frac{2}{n_{\max}+1}\cdot\frac{1}{N_t}$ — **$N_t$-dependent**. The headline value $2/(n_{\max}+1)$ exists only by mixing conventions across the tensor factors.
- **eq:B_def / eq:joint_berezin are wrong displays of the Berezin map.** Paper 38 def:berezin (718–728) and Paper 45 eq:joint_berezin (962–967) define $B$ as the weighted sum with weights $N/Z$ and claim equivalence with the convolution form $P M_{K*f} P$ ("the two definitions agree… by joint Plancherel", paper_45:970–975). They do not agree — at $f = 1$ the sum form gives $B(1) = (1/Z)\,I \to 0$ while the convolution form gives $B(1) = I$. Under the sum-form reading, L4(c) fails maximally at $f=1$ (deficit $1 - 1/Z \to 1$ vs RHS $0$) and the tunneling pair is not UCP; under the convolution reading all four L4 properties go through (the proofs of L4(a),(c),(d) in both papers in fact only use the convolution form). The convolution form is the only viable definition; note its multiplier support $N \le 2n_{\max}-1$ matches Paper 45's own envelope correction rem:envelope_v2 (842–859), which the sum form contradicts.

**What $2/(n_{\max}+1)$ actually is:** the sup of the normalized Plancherel window density $\dim V_\pi / Z$ — a correct computation about object (2), i.e. about the block-scalar map "multiply the $j$-block by the window weight." No step of the tunneling-pair assembly uses that map.

## V2. Citation check — WRONG citation + inapplicable-deep-theorem-for-elementary-fact

- **`bozejko_fendler1991` bibitem exists** in both papers (P45:1667–1671, P38:1338–1342) and is the Arch. Math. **57** (1991) 290–298 paper, *"Herz–Schur multipliers and uniformly bounded representations of **discrete** groups."* Web-verified: its content is that coefficients of uniformly bounded representations of **discrete** groups are Herz–Schur multipliers with $\|\varphi\|_{M_0A(G)} \le \|\pi\|^2\|\xi\|\|\eta\|$. It contains **no** statement of the form "central Fourier multiplier on a compact group has Schur cb-norm = $\ell^\infty$ norm of symbol." Hypotheses (discrete $G$; a u.b. representation) are not satisfied by, and the conclusion is not about, $\mathrm{SU}(2)\times U(1)$. The famous BF *multiplier* paper is the **1984** Boll. UMIA note ($B_2(G) = M_0A(G)$, cb Fourier multipliers = Herz–Schur multipliers); the sub-sprint memo (§10) asserts the 1991 paper is "the extended version with full proofs and amenability hypothesis" of the 1984 one — false; they are different results. `geovac/central_fejer_su2.py:646` cites "Bozejko-Fendler **1984**" while the papers cite **1991**: the project is internally inconsistent about which paper it means.
- **`pisier2001` bibitem exists** (P45:1800–1803; LNM 1618, 2nd ed. 2001). Web-checked TOC of the 2nd edition: the Herz–Schur/multiplier chapter is *"Hankelian Schur multipliers. Herz–Schur multipliers"* (Ch. 6 in 0-indexed, Ch. 7 in 1-indexed numbering); Chapter 8 in either convention is *"The similarity problem for cyclic homomorphisms on a C\*-algebra"* or *"Completely bounded maps in the Banach space setting"* — not multipliers on amenable groups. "Ch. 8 Thm 8.10" as a central-multiplier cb-equality: **UNVERIFIED-IMPORT** (full text not accessible), with strong structural evidence the cited statement is not there. The origin is `debug/r25_l2_proof_memo.md` "Theorem 5.1," a from-memory transcription (stated there for $T_m$ on $L^p(G)$ — a category error for cb-norms) that propagated verbatim through the sub-sprint memo into both papers.
- **The fact actually needed is elementary and needs no citation**: a central Fourier multiplier acts block-scalar on $\bigoplus_\pi M_{d_\pi}$ (the Plancherel/dual side), so its cb-norm there is $\sup_\pi|m(\pi)|$ — finite-dimensional linear algebra; amenability is irrelevant. Caveat the papers also miss: on the *function* side ($Z(C(G))$ with sup norm, where Lemma L2 places $S_K$), "cb-norm = $\sup|m|$" is **not** a general identity (the multiplier norm there is a $B(G)$-type norm $\ge \sup|m|$); it holds for our $K$ only because $K \ge 0$, $\int K = 1$ forces the norm to be attained at the trivial representation — i.e., the correct value is 1, again. Citing an inapplicable deep theorem for an elementary fact — and then evaluating the elementary fact on the wrong symbol — is the full anatomy of the error.
- Takesaki IV.4.14 (commutative tensor factorization): statement true and standard (commutative ⟹ nuclear); theorem number UNVERIFIED-IMPORT; harmless.

**Hypothesis checklist**

| Import | Statement as used in P38/45 | Actual content | Hypotheses met by our objects? |
|:--|:--|:--|:--|
| BF 1991 Arch. Math. 57 | central multiplier cb = $\sup$ symbol, amenable compact $G$ | u.b.-rep coefficients ⊂ $B_2(G)$, **discrete** $G$ | **No** (not discrete; no u.b. rep; conclusion absent) |
| Pisier 2001 "Ch. 8 Thm 8.10" | same | Ch. 8 = cyclic homomorphisms / Banach-space cb maps; Herz–Schur is Ch. 6/7 | **No at chapter level**; number unverifiable |
| Block-scalar cb fact | applied to symbol $(2j+1)/Z$ | true on dual side for the *actual* symbol $m_K$ | Math yes; **applied to wrong symbol & wrong side** |
| Takesaki IV.4.14 | $Z(C(G_1\times G_2)) \cong Z\hat\otimes Z$ | standard | Yes (number unverified) |
| Paper 38 L2(c) "verbatim" | cb $= 2/(n_{\max}+1)$ | L2(b) false; L2(a) contradicts L2(b) at $j=0$ | **No** |

## V3. Downstream use (Prop. reach_height) — SERIOUS-GAP

Paper 45:1170–1182: "*Lemma L2 gives $\|S_K\|_{\mathrm{cb}} = 2/(n_{\max}+1)$ finite and non-zero on the central subalgebra, so $\sigma$ exists there and is bounded by the same scale.*" Paper 38:868–873 is the same move ("via the L2(c) cb-norm equality… by symmetry of the convolution").

This is a non-argument three times over: (i) the constant is not the cb-norm of the map being inverted; (ii) a forward cb-norm bounds no partial inverse — a partial inverse of a multiplier is bounded by $1/\min|\text{symbol on support}|$, not by $\sup$; (iii) "bounded by the same scale" is backwards — a *small* forward norm makes inverses *large*.

Computed inverse scales (this review):

| Symbol used | min on support | $\|\sigma\|$ scale |
|:--|:--|:--|
| claimed mass symbol, window $N \le n_{\max}$ | $1/Z$ at $N=1$ | $Z = n_{\max}(n_{\max}+1)/2 = O(n_{\max}^2)$ |
| true symbol $m_K$, window $N \le n_{\max}$ | $\to \approx 0.318$ (empirical, $n\le 80$; $\approx 0.471, 0.470, 0.396$ at $n=2,3,4$) | $\le \approx 3.15$, **uniformly bounded** |
| true symbol $m_K$, full achievable envelope $N \le 2n_{\max}-1$ | $m_K(2n_{\max}-1) = \frac{2}{(2n_{\max}-1)(n_{\max}+1)}$ | $O(n_{\max}^2)$ |

So neither claimed scale is right, and the corrected argument faces a real fork: if the lift $\sigma$ only needs the window $N \le n_{\max}$, the symbol-inverse is $O(1)$ and a reach_P $\le C'\gamma$ bound with $C' = O(1)$ is plausible; if it must cover the natural-substrate envelope $N \le 2n_{\max}-1$ (which Paper 45's own rem:envelope_v2 says the multiplier algebra reaches), the naive symbol-inverse costs $O(n_{\max}^2)$ and **destroys the rate**. Moreover the inverse-multiplier bound controls operator-side norms, not the Lipschitz seminorm of the lifted *function* — the actual content of reach_P. The correct mechanism in the literature is Leimbach–vS Lemma 3.4/3.5 (web-verified): an **antiderivative Schur multiplier** $w(n) = (1-m(n))\,n_\mu/\|n\|^2$ plus Bożejko–Fendler-type *transference* (Schur ⟷ Fourier), which uses $[D,\cdot]$ to gain the decay instead of inverting the kernel tail — giving exactly $\|\mathrm{id} - \sigma_\Lambda\rho_\Lambda\| \le \gamma_\Lambda\|[D,\cdot]\|$. No SU(2) analog of that lemma is constructed in Paper 38 or 45. (The `central_fejer_su2.py:647` docstring even remembers this: "the symbol-side estimate that drives Lemma 3.4's antiderivative trick" — the $2/(n_{\max}+1)$ value was apparently *intended* as that small quantity, then mis-derived as a cb-norm of the Fejér multiplier.)

**Does reach_P $\le \gamma$ survive?** Plausibly yes in substance: Gaudillot-Estrada & van Suijlekom (arXiv:2310.14733, IMRN) prove GH convergence of truncated state spaces for **all compact metric groups** — independent corroboration of the qualitative conclusion for the scalar $C(\mathrm{SU}(2))$ sector (cited by Paper 40, not by Papers 38/45). But: their result is GH-of-state-spaces for $C(G)$ truncations, not Latrémolière propinquity for the chirality-doubled CH operator system; transporting it, or building the SU(2) antiderivative lemma, is genuine missing work. If the corrected constant $C' > 1$, the headline becomes $\Lambda \le \max(C', C_3)\,\gamma^{\mathrm{joint}}$ — rate unchanged, constant degraded; the §6 numerical panel (which "quotes $\gamma^{\mathrm{joint}}$ directly," P45:853–854) would then understate the proven bound by $C'$.

## V4. U(1) factor and N_t-independence — FIXABLE-GAP

The displayed kernel $K^{U(1)}_{N_t,T} = \frac{1}{N_t}\bigl|\sum_{|k|\le K_{\max}}e^{i2\pi k\theta/T}\bigr|^2$ (eq:K_U1_def) has true Fourier symbol $(N_t-|k|)/N_t$ on $|k|\le N_t-1$ (numerically verified at $N_t = 3, 5$), **not** the displayed $\max(0,(N_t+1-2|k|)/(N_t+1))$ (eq:K_U1_plancherel), which is the symbol of the *one-sided order-$K_{\max}$* Fejér kernel $\frac{1}{K_{\max}+1}|\sum_{k=0}^{K_{\max}}e^{ik\theta}|^2$ — a second, independent kernel/symbol mismatch. Both candidate symbols are multiplier-type with $\sup = 1$ at $k=0$, so the factor-value 1 and its $N_t$-independence hold under either; fix is to make kernel and symbol consistent. But Remark rem:cb_Nt_indep's headline ("joint cb-norm $2/(n_{\max}+1)$ independent of $N_t$") survives only as the trivial corrected statement "joint cb-norm $= 1$": as noted in V1, under the mass convention used on the SU(2) factor the U(1) mass is uniform $1/N_t$ and the product is $N_t$-dependent.

## V5. Paper 38 antecedent — same confusion, plus what survives

Same object confusion, same wrong citation, and Paper 45 inherits "L2(c) verbatim" (P45:753–759; sub-sprint memo §4 transcribes it as Lemma 4.1). Additional findings:

- **Code is circular, not confirmatory.** `geovac/central_fejer_su2.py:340–358` returns $(2j+1)/Z$ *by formula* (docstring asserts it is "the diagonal matrix element on $V_j \otimes V_j^*$" — the precise locus of the error); `central_multiplier_cb_norm` (630–651) returns `Rational(2, n+1)` *by formula*; the audit (line ~717) compares the function to the same formula. The sub-sprint memo's "numerical verification panel" (§8) verifies arithmetic of the formula against itself.
- **Contagion**: Paper 40 L2(b),(c) (paper_40:599–658) repeats the false symbol claim for all compact $G$ ($\widehat K(\pi) = \dim V_\pi/Z$ — the proof at 639–645 states the correct symbol definition and then asserts the false evaluation; the squared-modulus cross-terms via the Clebsch–Gordan series are dropped), uses it at reach_P (line 1605), and at 586–591 claims the $\sqrt{\dim V_\pi}$ weight is *determined by* the cb-norm requirement — a rationale that collapses since any normalized non-negative central kernel has smoothing cb-norm 1. Paper 46 inherits via Paper 45 (not audited here). The memos `r25_l2_proof_memo.md` (part (g), "Theorem 5.1") and `l3b_2_sub_sprint_B_cb_norm_memo.md` (Theorem 3.1, Lemma 4.1, Theorem 7.1) are the source documents.
- **What survives unconditionally** (these are findings too): L2(a) normalization; L2(d) $\gamma_{n_{\max}}$ closed form, $4/\pi$ asymptote and uniform bound (defined directly as a kernel moment — no dependence on the symbol claims); L3 ($C_3$); L4(a),(c),(d) and contractivity-with-constant-1, *under the convolution-form definition of $B$* (which is the form all four proofs actually use); reach_B and height_B in L5 (cite only L3/L4 + Stein–Weiss); height_P $= 0$. The false constant $2/(n_{\max}+1)$ does **not** enter the headline numerically — it is load-bearing only as the (invalid) existence/boundedness argument for $\sigma$ in reach_P.

## Corrected statements (for a future fix sprint; not applied)

1. **L2(b)′:** the convolution symbol of $K_{n_{\max}}$ is $m_K(N) = \frac{1}{N Z}\sum_{a,b \le n_{\max}}^{\mathrm{CG}}\sqrt{ab}$, with $m_K(1)=1$, decreasing, support $N \le 2n_{\max}-1$.
2. **L2(c)′:** $S_K$ and $B = P M_{K*\cdot}P$ are unital CP; $\|S_K\|_{\mathrm{cb}} = \|B\|_{\mathrm{cb}} = 1$ for all $n_{\max}, N_t$ (joint included; trivially $N_t$-independent). The quantity $2/(n_{\max}+1)$ is $\sup_j$ of the Plancherel window density $\dim/Z$ and plays no role in the tunnel bounds.
3. **B-def′:** define $B$ by the convolution form only; delete the $N/Z$-weighted sum display or restate it with weights $m_K(N)$ over $N \le 2n_{\max}-1$.
4. **reach_P′:** requires an SU(2) analog of Leimbach–vS Lemmas 3.4–3.5 (antiderivative Schur multiplier + transference), or adaptation of Gaudillot-Estrada–van Suijlekom 2310.14733; expected $\mathrm{reach}_P \le C'\gamma$ with $C' = O(1)$, fork on window-vs-envelope to be resolved. Until then the main theorems of Papers 38/45 hold with reach_P as a **named gap**.
5. Fix the U(1) kernel/symbol pair to one consistent convention; replace the BF-1991/Pisier-Ch.8 citations with either no citation (elementary block-scalar fact, dual side) or BF 1984 / Pisier Ch. 6 for the transference form actually needed in reach_P′.

## Overall verdicts

**Lemma L2 (P45 lem:L2 and P38 L2(b),(c)): WRONG as written.** The "Plancherel symbol" is the mass distribution of the kernel amplitude, not the multiplier symbol; the computed sup is the cb-norm of a map the proof never uses; the true cb-norm of the map it does use is exactly 1 (unital CP). Internal contradictions: L2(a) vs L2(b) at the trivial rep; "UCP with cb-norm $< 1$" at P45:1106; sum-form vs convolution-form Berezin.

**reach_P (P45 prop:reach_height eq:reach_P; P38 L5): SERIOUS-GAP.** The cited invariant is irrelevant to a partial inverse (and the "same scale" claim is backwards); no valid proof of reach_P $\le \gamma^{\mathrm{joint}}$ exists in either paper. The conclusion is plausibly true (Gaudillot-Estrada–vS for the scalar sector; uniformly bounded window-inverse $\approx 3.15$ computed here), the rate is plausibly unchanged, and the headline constant changes at most by an $O(1)$ factor — but closing it requires a genuinely new lemma (SU(2) antiderivative transference), not an edit.

**Files:** paper_45_lorentzian_propinquity.tex:530–568, 728–806, 783–792, 962–1024, 1106–1108, 1170–1182; paper_38_su2_propinquity_convergence.tex:498–530, 557–563, 718–728, 745–800, 863–873, 994–1010; paper_40_unified_propinquity_convergence.tex:586–658, 1605; geovac/central_fejer_su2.py:340–358, 630–651; debug/r25_l2_proof_memo.md §5.1; debug/l3b_2_sub_sprint_B_cb_norm_memo.md §3–§7.

**Web sources:** [BF 1991 context (Steenstrup, arXiv:1001.0326)](https://arxiv.org/abs/1001.0326); [u.b. reps & cb multipliers survey (arXiv:1707.08329)](https://arxiv.org/pdf/1707.08329); [Pisier LNM 1618 2nd ed. TOC (Springer)](https://link.springer.com/book/10.1007/b55674); [Leimbach–vS tori (arXiv:2302.07877)](https://arxiv.org/html/2302.07877v2); [Gaudillot-Estrada–vS compact metric groups (arXiv:2310.14733)](https://arxiv.org/abs/2310.14733).
