# Paper 45 adversarial referee review — L5 assembly (§6) and framework validity

**Sprint:** Phase 1b of `docs/corpus_accessibility_plan.md` — adversarial review, mandate REFUTE.
**Date:** 2026-06-09. **Scope:** `papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex` (P45), with Paper 46, Paper 38, the L3b-2 sub-sprint D memo, the Phase 1B-A closure memo, and the two cited Latrémolière papers (both fetched and read from the actual PDFs).
**Overall verdict: thm:main does NOT survive as stated.** Two independent fatal defects (V1 degeneracy, V2 fabricated framework citations), plus a per-constituent counterexample (height_B). Several components survive and are catalogued in §7.

---

## 1. V1 — The degeneracy question. Verdict: **WRONG (counterexample)**

### 1.1 The seminorm kernel strictly contains the scalars (full-space level)

P45's operator system (eq:op_system_def, P45:498–519) is spanned by pure tensors $M^{\rm spat}_{N,L,M}\otimes M^{\rm temp}_p$ with the temporal multipliers **momentum-polynomial diagonal**, $M^{\rm temp}_p=\mathrm{diag}(\omega_k^p)$ (eq:temp_multipliers, P45:512–514); Vandermonde gives the equivalent basis of momentum indicators $\mathrm{diag}(\delta_{k,q})$ (P45:515–519). The spatial label set includes $N=1$, whose Avery 3-Y multiplier is a multiple of the identity (and L1′(a), P45:649–650, guarantees unitality regardless).

**Counterexample.** For $N_t\ge 2$ take $a = \mathbf 1_{\rm spat}\otimes\mathrm{diag}(\delta_{k,q})\in\mathcal O^L$, $a\notin\mathbb C\cdot\mathbf 1$. By the paper's own structural identity (eq:L3_struct_id, P45:867–871), $[D_L,a]=i[D_{\rm GV},\mathbf 1]\otimes a_t=0$. Hence the full-space Lipschitz seminorm $L_{\rm op}(a)=\|[D_L,a]\|$ has $\ker L_{\rm op}\supseteq \mathbf 1\otimes\mathcal A_{\rm temp}$, dimension $N_t>1$.

By Latrémolière's actual framework this is disqualifying: **Prop. 1.10 of arXiv:1811.10843** (pp. 5–6 of the fetched PDF) proves a spectral triple is *metric* iff $(\mathfrak A,\mathsf L_D)$ is a Leibniz quantum compact metric space, and derives $\{a:\mathsf L_D(a)=0\}=\mathbb R 1_{\mathfrak A}$ as part of that definition (the MK metric must metrize the weak-* topology). With a non-scalar kernel element, states differing on $\mathbf 1\otimes\mathcal A_{\rm temp}$ are at MK distance $+\infty$; the object is not a quantum compact metric space and the propinquity is undefined on it.

### 1.2 Worse: the K⁺-compressed seminorm of def:weak_form_propinquity is **identically zero**

P45's Definition def:weak_form_propinquity (P45:589–604) takes $\mathcal T^+=(P_+\mathcal O^LP_+,\;\mathcal K^+,\;P_+D_LP_+)$ — the **compressed** Dirac. Krein-self-adjointness of $D_L$ with $J_L=\gamma^0\otimes I$ *forces* $\{\gamma^0,D_{\rm GV}\}=0$ (P45:485–492), hence $P_+D_{\rm GV}P_+=\tfrac14(D_{\rm GV}+\{\gamma^0,D_{\rm GV}\}+\gamma^0D_{\rm GV}\gamma^0)=0$. So $P_+D_LP_+\cong i\,I_{\rm Weyl}\otimes\partial_t$ — **the spatial Dirac is annihilated by the K⁺ compression**. Since every compressed multiplier $W\otimes m$ has momentum-diagonal $m$, $[P_+D_LP_+,P_+aP_+]=W\otimes i[\partial_t,m]=0$ for **every** $a\in\mathcal O^L$. The K⁺-restricted Lipschitz seminorm is $\equiv 0$ on the entire operator system.

**Paper 46 already documents this.** P46 eq:LPwk_def (P46:523–525) defines exactly this seminorm and states (P46:529–547): "$P_+D_L^{\rm off}P_+=0$ … $P_+D_LP_+$ *keeps only the temporal piece*, and $L^{\rm wk}(a)$ depends only on the temporal commutator $[\partial_t,M^{\rm temp}]$. On the natural substrate, $M^{\rm temp}$ is momentum-diagonal and $[\partial_t,M^{\rm temp}]=0$…" P46's Lemma 3.2 (lem:temporal_invisibility, P46:591–645) proves $[D^{\rm diag},a]\equiv 0$ for all $a\in\mathcal O^L$. P46 frames its own $L_{\rm op}$ as "a strictly larger seminorm than $L^{\rm wk}$" (P46:542–547) — i.e., the corpus itself certifies that Paper 45's weak-form seminorm carries no metric content at all. Neither paper anywhere addresses the kernel condition $L(a)=0\iff a\in\mathbb C\mathbf 1$.

### 1.3 The continuum reference triple is degenerate too — in a *different* direction

def:tunneling_pair (P45:1089–1104) restricts the continuum triple the same way. On the continuum, multiplication operators do not commute with $\partial_t$, so $P_+D^L_\infty P_+ = i\partial_t$ gives $L^{\rm wk}_\infty(f)=\|\partial_t f\|_\infty$ — kernel $=$ all $f$ constant in $\tau$ ($\cong C^\infty(S^3)$, infinite-dimensional). So def:weak_form_propinquity feeds the Latrémolière propinquity two inadmissible arguments with *mismatched* degeneracies (truncated: everything invisible; continuum: spatial directions invisible). There is no reading under which $\Lambda(\mathcal T^+_{\rm trunc},\mathcal T^+_\infty)$ is a defined, finite quantity in the cited framework; naively extended, both state spaces have infinite diameter in the degenerate directions and the GH-type distance is $+\infty$, contradicting the bound $2.0746$.

### 1.4 What the panel numbers actually are

`geovac/lorentzian_propinquity_compact_temporal.py:563–592` computes `propinquity_bound = max(reach_B_theoretical, gamma_su2, height_B_theoretical, 0)` — closed-form **rate formulas**, never a distance between the constructed objects. The Table (P45:1358–1382) "$\Lambda^L$ bound" column equals $\gamma^{SU(2)}_{n_{\max}}$ at every cell — Paper 38's spatial values, independent of $N_t$. The sub-sprint D memo concedes the structure: "the operator-system Lipschitz seminorm is governed by the spatial Dirac, not by the temporal Fourier truncation" (`debug/l3b_2_sub_sprint_D_propinquity_memo.md` §7).

The $N_t=1$ "bit-exact Riemannian recovery" (prop:riemannian_limit, P45:1387–1413) is therefore recovery of a **formula**, not of the defined object: at $N_t=1$, $\partial_t=0$, so the literal $\mathcal T^+$ has *zero* Dirac and zero seminorm, while Paper 38's triple has the full CH commutator seminorm — they cannot be equal as metric objects. The bit-exact match is automatic because the bound formula never depended on the temporal factor.

### 1.5 Root cause: the temporal substrate is not the temporal analog of the spatial substrate

Paper 38's spatial multipliers are **compressed multiplication operators** $P\,M_{Y_{NLM}}P$ (Toeplitz-type, Avery 3-Y). The honest temporal analog is the Connes–van Suijlekom circle truncation $P\,M_{e^{iq\tau}}P$ — a truncated **shift** in the momentum basis. P45 instead uses momentum-**diagonal** multipliers, i.e., functions *of the momentum operator* (functions of the Dirac's temporal part), which trivially commute with $\partial_t$. Consequences:

- The two definitions of $B^{\rm joint}$ given in P45 are **inequivalent**: the sum form (eq:joint_berezin, P45:962–967) sends $e^{iq\tau}\mapsto \hat K(q)\,\mathrm{diag}(\delta_{k,q})$ (rank-1 diagonal), while the compression form $P^{\rm joint}(K*f)P^{\rm joint}$ (P45:970–975) sends it to $\hat K(q)\,S_q$ (truncated shift). The claimed agreement "on the central subalgebra by joint Plancherel" (P45:974–975) is false for $q\ne 0$ on the U(1) factor.
- **L4(c)** (approximate identity, P45:986–993) fails for the temporal factor under the only reading in which it is meaningful ($P M_f P$ = compressed multiplication): at $N_t=3,T=2\pi,q=1$, $\|B(e^{i\tau})-PM_{e^{i\tau}}P\|=\|\tfrac12\mathrm{diag}(\delta_{k,1})-S_1\|\ge 1 > \gamma^{U(1)}\!\cdot\!\|\partial_\tau e^{i\tau}\| = 0.722$. The L4(c) proof (P45:1026–1050) actually bounds $\|PM_{K*f}P-PM_fP\|$ and silently identifies $B(f)$ with $PM_{K*f}P$, which fails on the temporal factor.
- **L3's headline identity** (eq:L3_struct_id) is true *only because* of this substitution; with the honest Toeplitz substrate it fails ($[\partial_t,PM_{e^{iq\tau}}P]\ne 0$).

The degeneracy (V1) and the Berezin inconsistency are two symptoms of one substrate-level conflation: the temporal function algebra was replaced by its Fourier-dual.

### 1.6 Direct counterexample to prop:reach_height eq:height_B

Under P45's own height definition (P45:1186–1196): take $f=\mathbf 1\otimes f_t$ with $\|\partial_t f_t\|_\infty=1$ (admissible for $N_t\ge3$). Then $\|f\|_{\rm Lip}=1$, while $B^{\rm joint}(f)=\mathbf 1\otimes\mathrm{diag}(\hat K(k)c_k)$ commutes with $D_L$, so $\|B(f)\|_{\rm Lip^{op}}=0$ and $\mathrm{height}_B\ge 1$ for **all** $(n_{\max},N_t\ge 3)$. Since the claimed bound $C_3\gamma^{\rm joint}\to 0$ in the joint limit, eq:height_B (P45:1135–1136) is **false asymptotically** (at the small panel cells $\gamma^{\rm joint}>1$ masks the violation). The "Stein–Weiss … the U(1) factor inherits the standard Fejér-on-the-circle gradient control" step (P45:1201–1207) is exactly where this fails: Fejér gradient control holds for compressed multiplication operators, not for momentum-diagonal images whose commutator seminorm is identically zero.

---

## 2. V2 — Citation fidelity. Verdict: **WRONG — fabricated attributions; actual device is van Suijlekom's, scenario (iii)**

Both candidate source papers were fetched and read.

### 2.1 Claimed vs. actual, arXiv:1811.10843 = Adv. Math. 404 (2022) 108393 (the cited `latremoliere_metric_st_2017`)

| P45 claim (location) | Actual content of cited paper (verified, page refs to v6 PDF) |
|---|---|
| "Thm 5.5 — the weak-form metric-spectral-triple propinquity bound for a direct UCP tunneling pair"; "$\Lambda$ is the infimum over all UCP tunneling pairs of the four-term maximum" (thm:main proof, P45:1234–1245) | **No Theorem 5.5 exists.** §5 ("examples") contains Thm 5.1 (Sierpiński gasket), Thm 5.2 (quantum/fuzzy tori), Remark 5.3 — last numbered item (pp. 36–43); references begin p. 44. The propinquity is $\Lambda^{\rm spec}=\inf\{\chi(\tau):\tau$ a metrical **tunnel**$\}$ (Defs 2.19/2.21/2.22, pp. 16–17), tunnels built from **quantum isometries = *-epimorphisms** (Def 1.12 p. 7, Def 1.15), never UCP pairs. No "weak-form" anywhere. |
| "reach$_B$ … (Latrémolière 2018/2023 Def. 3.4)" (P45:1147–1150) | Item 3.4 is a **Remark** about $\||\beta^g\||_{\mathscr M}\le 1$ for module endomorphisms in covariant modular systems (p. 19). Nothing about reach or Berezin maps. |
| "height$_B$ … (ibid. Def. 3.5)" (P45:1150–1153) | Item 3.5 is a **Notation**: $G[r]$ = closed balls in a proper metric monoid (p. 20). Nothing about height. |

### 2.2 Same check against the *other* cited Latrémolière paper (Trans. AMS 368 (2016) 365–411, arXiv:1302.4058), since the internal memo conflates the two ("Latrémolière, arXiv:1811.10843, Trans. AMS 368 (2016) 365–411, §4 weak-form", sub-sprint D memo §1.2)

Lemma 3.4 = 1-level state characterization; Remark 3.5 = pivot norm $\ge 1$; Def 3.6 = **bridge** $(\mathfrak D,\omega,\pi_{\mathfrak A},\pi_{\mathfrak B})$ with unital ***-monomorphisms** and pivot; item 5.5 is a **Definition** (itinerary for a trek). A bridge has ONE height (Def 3.16) and ONE reach, both defined via the pivot — not four constituents attached to two UCP maps. Standing hypotheses: Leibniz QCMS over unital C*-algebras with $\{L=0\}=\mathbb C 1$ (Prop 2.17, Def 2.19).

**Conclusion:** the citation is fabricated in both candidate sources — the numbered items exist but say unrelated things, and the central object ("direct UCP tunneling pair" with $\max(\mathrm{reach}_B,\mathrm{reach}_P,\mathrm{height}_B,\mathrm{height}_P)$) exists in neither.

### 2.3 What the device actually is — scenario (iii)

The $(B,P)$ pair with (i) approximate identity $\|B(f)-PM_fP\|\le\gamma\,{\rm Lip}(f)$, (ii) Lipschitz contraction, (iii) partial inverse on a dense subalgebra is the **C¹-approximate order isomorphism** technique of van Suijlekom, *"Gromov–Hausdorff convergence of state spaces for spectral truncations"*, J. Geom. Phys. 158 (2020) 103866, arXiv:2005.08544 (operator system spectral triples; UCP maps are the native morphisms; conclusion is GH convergence of state spaces under the Connes/MK metric), as used for tori by Leimbach–vS, Adv. Math. 439 (2024) 109496, arXiv:2302.07877 — both already in P45's bibliography for other purposes. Even there, the framework presupposes that the MK metric on each side metrizes the weak-* topology (finite distances), which V1's degeneracy destroys. So a vocabulary swap alone does not rescue thm:main; it rescues only the spatial-quotient statement (§6 below).

### 2.4 Provenance note

The fabricated numbers were inserted by the Phase 1B-A "L5 closure" sprint (2026-05-24): its memo records as deliberate edits "Connects each constituent to the corresponding Latrémolière 2017 §3.4 / §3.5 definition" and "Opens with explicit citation of Latrémolière 2018/2023 §5 Theorem 5.5" (`debug/sprint_phase1b_a_paper45_l5_closure_memo.md` §1, Edits 2–3). The pre-existing sub-sprint D memo had cited only "Latrémolière 2017 §4 (or Paper 38 §L5 / Leimbach–vS 2024 §4)". I.e., the "closure" sprint *increased* the specificity of citations without checking the source — a verify-the-verifier failure mode. (Paper 38's main theorem is its §4; "Thm 5.5" matches no GeoVac-internal numbering either.)

---

## 3. V3 — Hypothesis checklist (Latrémolière 1811.10843 membership requirements)

| Hypothesis (source) | Truncated $\mathcal T^+_{n_{\max},N_t,T}$ | Continuum $\mathcal T^+_{\mathcal M}$ (def:tunneling_pair) |
|---|---|---|
| Unital **C*-algebra** (Def 1.1, p. 2) | **FAILS** — $\mathcal O^L$ is a *-closed linear subspace, "not multiplicatively closed" (P45:521–525) | full algebra OK; restricted version VERIFIABLE-BUT-UNSTATED ($P_+C^\infty P_+$ is an algebra since scalar $f$ commutes with $\gamma^0$) |
| **Leibniz** Lip-norm (Prop 1.10, pp. 5–6; needs algebra structure) | **FAILS** — not statable on an operator system | VERIFIABLE for the full triple; vacuous/odd for the restricted seminorm |
| **Kernel = scalars** (Prop 1.10 proof, p. 5) | **FAILS** — $L^{\rm wk}\equiv 0$ (§1.2); even full-space $L_{\rm op}$ has $\ker\supseteq\mathbf 1\otimes\mathcal A_{\rm temp}$ | **FAILS** for restricted ($\ker = $ time-independent functions); holds for the full (unrestricted) triple |
| MK metric **metrizes weak-*** (Def 1.8) | **FAILS** (distances $+\infty$) | **FAILS** (restricted); holds unrestricted |
| Lower semicontinuity of L | OK (finite-dimensional) | OK |
| Self-adjoint D, dense domain | OK ($P_+D_LP_+$ Hermitian) | OK |
| Tunnel morphisms = **quantum isometries / *-epimorphisms** (Defs 1.12, 1.15) | **FAILS** — $(B,P)$ are UCP maps; no tunnel/bridge is ever constructed; no quasi-Leibniz middle space exhibited | — |
| $\Lambda^{\rm spec}=\inf$ extent over tunnels (Def 2.22) | **NOT USED** — replaced by an uncited four-max over a UCP pair | — |

Note Latrémolière's §5.1 example (pp. 36–38) shows what an actual convergence proof in his framework looks like: an explicit Leibniz tunnel on $\mathfrak A\oplus\mathfrak A$ with seminorm $\mathsf S(a,b)=\max\{\mathsf L(a),\mathsf L_T(b),C\|a-b\|\}$ and coordinate *-epimorphisms. Nothing of this kind appears in P45 (or P38).

---

## 4. V4 — K⁺-preservation bookkeeping. Verdict: **SERIOUS-GAP (and the asserted bounds are false under the restricted seminorm)**

What is verified-in-paper: $[J_L,B^{\rm joint}(f)]=0$ and $[J_L,P^{\rm joint}]=0$ at the operator level (rem:Kplus_tunnel_valid, P45:688–713; L4(e), P45:1061–1065). This makes the compressed pair well-defined as maps. It does **not** make it a tunneling pair "between the K⁺-restricted triples" in any metric sense, because:

1. **Unit-ball swap.** Every constituent in prop:reach_height (P45:1156–1207) is computed as a supremum over $\{f:\|\nabla^{\rm joint}f\|_\infty\le 1\}$ — the classical joint gradient ball, an *extraneous third seminorm*. Under the K⁺-restricted triples' own Lipschitz seminorms (the only ones Definition def:weak_form_propinquity provides), the continuum-side unit ball is $\{f:\|\partial_tf\|_\infty\le1\}$, which contains functions of arbitrarily large spatial gradient; the truncated-side unit ball is all of $P_+\mathcal O^LP_+$ ($L^{\rm wk}\equiv0$). The suprema of the deficits over those balls are $+\infty$ (e.g. reach$_B$: $f=\sin(Kx_{\rm spat})$, deficit $\sim\gamma_{SU}K\to\infty$). The phrase "the K⁺-compression of the deficit is bounded by the unrestricted deficit" (P45:1164–1168) is correct per fixed $f$ (compression contracts operator norm) but the ball-level inequality runs the other way — larger ball, larger sup.
2. **State-space step missing.** In any honest framework (Latrémolière bridges/tunnels, or vS state-space GH) the reach/height are state-space quantities; nothing in P45 shows that states of the restricted systems extend/restrict compatibly, and with degenerate seminorms they provably do not behave (infinite MK distances).
3. height$_B\le C_3\gamma$ is false outright (§1.6); height$_P=0$ "no Lipschitz distortion" (P45:1209–1211) is asserted for a projection without specifying with respect to which pair of seminorms — under the restricted seminorms it is vacuous-true only because everything has seminorm 0.

---

## 5. V5 — Paper 38 antecedent. Verdict: **FIXABLE-GAP (citation framing), with one unproven robustness assertion**

Paper 38 §L5 (P38:818–928) uses the same four-constituent $(B,P)$ device but is more honest about provenance: rem:height_definition (P38:915–928) calls its height "the spectral-triple-level **analog** of the height in [Latrémolière §4]" — an admission that the definitions are home-made. P38 contains **no** "Thm 5.5 / Def 3.4 / Def 3.5" attribution; those entered the corpus in P45 via the Phase 1B-A sprint (§2.4). However:

- P38's main theorem (P38:934–949) still claims convergence "in the Latrémolière quantum Gromov–Hausdorff propinquity," which is the wrong framework for an operator-system substrate (no algebra, no Leibniz, tunnels-not-UCP-pairs). The honest home is van Suijlekom's state-space GH distance (§2.3), where the content plausibly survives.
- P38 is aware of its own kernel-condition wobble: L1′ (P38:431–485) introduces the offdiag CH Dirac "purely SDP-bounding: it ensures that the operator-norm constraint … is non-degenerate on the full operator system" while L3/L4 use the *truthful* CH Dirac, and bridges the two with "related by a bounded compact perturbation, and the propinquity bound … is robust under this choice" (P38:472–474). That robustness is asserted, not proved: Latrémolière's own bounded-perturbation lemma (1811.10843 §5.1, Eq. (5.1)) requires $2r\||T\||<1$ (perturbation small relative to inverse diameter), and the E1 offdiag coupling is $O(1)$. Whether $\ker\|[D^{\rm truthful}_{\rm CH},\cdot]\|\cap\mathcal O_{n_{\max}}=\mathbb C\mathbf 1$ is UNKNOWN in-paper (the WH1 R3.5 corpus result claims it for the *offdiag* operator). So P38 needs (i) reframing in vS distance, (ii) a one-page kernel-condition verification for the seminorm actually used. Its spatial convergence content is very likely true.

---

## 6. Strongest true statement, and the cost

**As stated, thm:main (P45:1216–1232) is not a theorem about any defined metric quantity.** The strongest true statements extractable, in decreasing strength:

- **(S1) What the present construction actually proves.** *The spatial quotient converges.* On $\mathcal O^L/(\mathbf 1\otimes\mathcal A_{\rm temp})$ the full-space seminorm $L_{\rm op}$ descends to Paper 38's spatial seminorm, and (modulo P38's own vS-reframing + kernel check, §5) the truncations of $C(S^3)$ converge in van Suijlekom's state-space GH distance at rate $C_3\gamma^{SU(2)}_{n_{\max}}$; the temporal factor is a **metrically invisible spectator** carried $J_L$-equivariantly ($[J_L,\cdot]=0$ verified bit-exact for all maps). This is exactly what the panel numbers ($\Lambda$ column $=\gamma_{SU}$, $N_t$-independent) measure. It is a Riemannian SU(2) statement plus a Krein-equivariance decoration; it has **no Lorentzian metric content**.
- **(S2) What a repaired construction could prove.** Replace the temporal substrate by the Connes–vS circle truncation (Toeplitz compressions $PM_{e^{iq\tau}}P$). Then the kernel condition plausibly holds, L4(c) becomes true with the Fejér rate, and combining Paper 38's SU(2) side with Leimbach–vS's circle/torus technique gives an honest **vS-GH convergence theorem for spectral truncations of $S^3\times S^1_T$, equivariant under the BBB (4,6) Krein structure** — at the cost that eq:L3_struct_id is false there ($C_3^{\rm joint}$ and the rate acquire genuine temporal contributions; the "temporal contributes nothing" headline disappears). This would still be a Wick-rotated *Riemannian product* theorem; "Lorentzian" survives only as the equivariance bookkeeping of $(J_L, D_L,$ BBB signs$)$.
- **(S3)** The verified algebra: L2's $\|S_{K^{\rm joint}}\|_{\rm cb}=2/(n_{\max}+1)$ (correct, seminorm-independent); eq:L3_struct_id as a bit-exact algebraic identity *of this substrate*; L4(a),(b),(e); the closed-form rate panel as formula consistency with P38.

**Cost to the novelty claim.** "First Lorentzian propinquity convergence theorem on truncated Krein spectral triples in the published math.OA literature" (abstract, P45:174–176) must be withdrawn in its current form: (a) the quantity bounded is not the Latrémolière propinquity (no tunnel, wrong category, fabricated theorem citation); (b) the objects are not quantum compact metric spaces (degenerate seminorms — identically zero on the truncated K⁺ side); (c) the Lorentzian/temporal direction carries zero metric content even in the most charitable reading, so the claim is vacuous *as a Lorentzian claim* — consistent with the bit-exact $N_t=1$ recovery being automatic. Downstream exposure: Paper 46 ("free upgrade" of a zero seminorm; same kernel defect at $L_{\rm op}$ level), Paper 47 §7 / Paper 48 / Paper 49 inherit the P45 panel and the "$\Lambda^{\rm P45}$" quantity; the Mondino–Sämann bridge papers transport the panel "verbatim" and inherit the defect.

## 7. What survives (findings)

(i) L2 joint cb-norm theorem — clean. (ii) eq:L3_struct_id / P46 Lemma 3.2 as bit-exact algebra — clean, but it is the *diagnosis* of the degeneracy, not a feature. (iii) $[J_L,\cdot]=0$ equivariance of the whole tunneling apparatus — clean. (iv) Paper 38's spatial convergence content — likely true after vS-reframing + kernel check. (v) The honest-scope discipline (Q1/G1, P45:1483–1498) correctly flags that no Krein-signature metric exists in the literature — but mislocates the consequence: the missing ingredient invalidates the *present* weak form too, since the K⁺ trade destroys the Dirac's metric content rather than preserving it.

---

### Sources verified
Latrémolière, arXiv:1811.10843v6 = Adv. Math. 404 (2022) 108393 (full PDF read: Def 1.1, Prop 1.10, Defs 1.12/1.15, Defs 2.19–2.22, Rem 3.4, Not 3.5, §5 = Thms 5.1, 5.2, Rem 5.3, refs from p. 44). Latrémolière, arXiv:1302.4058 = Trans. AMS 368 (2016) 365–411 (Lem 3.4, Rem 3.5, Def 3.6, Def 3.10, Def 5.5, bridge height Def 3.16). van Suijlekom, J. Geom. Phys. 158 (2020) 103866, arXiv:2005.08544; Leimbach–vS, Adv. Math. 439 (2024) 109496, arXiv:2302.07877 (web-verified).
