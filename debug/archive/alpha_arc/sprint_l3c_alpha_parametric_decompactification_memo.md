# Sprint L3c-α — Parametric de-compactification on the natural substrate

**Date:** 2026-05-23
**Sprint goal:** Verify that the Paper 45 K⁺-weak-form / Paper 46 strong-form (main theorem, natural substrate) propinquity convergence bound survives the **joint parametric limit** $(n_{\max}, N_t, T) \to (\infty, \infty, \infty)$ along an explicit coupled scaling $T = T(N_t)$ satisfying $T \to \infty$ and $T/N_t \to 0$.
**Sprint outcome:** **CLOSED at the analytical level on the natural substrate** in one session. The result follows from Paper 45 / Paper 46 main theorems with no new analytical content — every $T$-dependent ingredient in the propinquity bound has a $T$-independent structural form, with the sole $T$-dependence entering through L4's $\gamma^{U(1)} = O(T/N_t)$ Fejér-on-circle rate, which vanishes along any coupled scaling $T(N_t)/N_t \to 0$.
**Scope:** **Parametric** de-compactification along a sequence of compact-temporal triples with $T_n \to \infty$, NOT the genuine non-compact limit $\mathcal{T}^L_{S^3 \times \mathbb{R}_t}$. The latter is L3c-γ (inductive limit, separate sub-sprint).
**No production code changes. No paper edits applied** (queued for PI decision per §6 below).

---

## §1. Setup and notation

Adopt Paper 45 §3 notation verbatim:
- $\Krein_{n_{\max}, N_t, T} = \HGV^{n_{\max}} \otimes \C^{N_t}$, Krein space at compact temporal radius $T$.
- $\JL = J_{\mathrm{spatial}} \otimes I_{N_t}$ at BBB $(m, n) = (4, 6)$ Peskin-Schroeder chiral basis.
- $\DL = i(\gamma^0 \otimes \partial_t + \DGV \otimes I_{N_t})$ Lorentzian Dirac (Paper 45 §3, vdD 2016 Prop 4.1).
- Natural substrate $\OpL_{n_{\max}, N_t}$: chirality-doubled scalar spatial multipliers $M^{\spat} = \mathrm{blkdiag}(W, W)$ tensored with momentum-diagonal temporal multipliers $M^{\temp}$.
- Joint Berezin map $\Bjoint = B^{\SU(2)} \otimes B^{U(1)_T}$ (PURE_TENSOR).

For L3c-α, introduce a **scaling sequence** $T : \mathbb{N} \to \R_{>0}$ assigning to each $N_t$ a temporal radius $T(N_t)$. The L3c-α theorem will hold for any scaling sequence satisfying

$$\boxed{\,T(N_t) \to \infty \quad \text{and} \quad T(N_t)/N_t \to 0 \quad \text{as } N_t \to \infty\,}$$

(call this **admissible scaling**). Canonical examples:
- **Sub-linear**: $T(N_t) = N_t / \log N_t$ — gives $T/N_t = 1/\log N_t \to 0$ slowly; $T \to \infty$ near-linearly.
- **Square-root**: $T(N_t) = \sqrt{N_t}$ — gives $T/N_t = 1/\sqrt{N_t} \to 0$ at $N_t^{-1/2}$ rate; $T \to \infty$ as $N_t^{1/2}$.
- **General**: $T(N_t) = N_t^{1-\epsilon}$ for any $\epsilon \in (0, 1)$.

All three are admissible; the choice trades de-compactification rate against propinquity decay rate.

---

## §2. Theorem (L3c-α, parametric de-compactification)

**Theorem (L3c-α).** Let $T(\cdot)$ be any admissible scaling sequence (§1). On the natural substrate of Paper 45 / Paper 46 main theorem, the K⁺-weak-form propinquity bound
$$\Lprop\bigl(\mathcal{T}^L_{n_{\max}, N_t, T(N_t)},\ \mathcal{T}^L_{S^3 \times S^1_{T(N_t)}}\bigr) \;\le\; C_3^{\mathrm{joint}}(n_{\max}) \cdot \gamma^{\mathrm{joint}}\bigl(n_{\max}, N_t, T(N_t)\bigr)$$
holds at every $(n_{\max}, N_t)$, with $C_3^{\mathrm{joint}}(n_{\max}) = \sqrt{1 - 1/n_{\max}} \to 1^{-}$ (Paper 46 §4.3 envelope-aware form, $T$- and $N_t$-independent) and
$$\gamma^{\mathrm{joint}}\bigl(n_{\max}, N_t, T(N_t)\bigr) \;=\; \max\bigl(\gamma^{\SU(2)}(n_{\max}),\ \gamma^{U(1)}\bigl(N_t, T(N_t)\bigr)\bigr) \;\to\; 0$$
as $(n_{\max}, N_t) \to \infty$, with $\gamma^{\SU(2)} = O(\log n_{\max}/n_{\max})$ and $\gamma^{U(1)} = O(T(N_t)/N_t) \to 0$ by admissibility.

The strong-form analog (Paper 46 main theorem, natural chirality-doubled scalar substrate) gives the same bound $\Lambda^{\mathrm{strong}}(\nmax, N_t, T(N_t)) = \Lprop(\nmax, N_t, T(N_t))$ bit-exact under the "free upgrade" reading.

**Honest scope.** The limit triple at each finite stage $(n_{\max}, N_t)$ is $\mathcal{T}^L_{S^3 \times S^1_{T(N_t)}}$ — a **compact-temporal** triple at radius $T(N_t)$. As $N_t \to \infty$ the radius $T(N_t) \to \infty$, but at every finite stage the target is still compact-temporal. The genuine non-compact target $\mathcal{T}^L_{S^3 \times \mathbb{R}_t}$ does not appear as a literal limit object in this theorem; reaching it requires the inductive-limit construction of Sprint L3c-γ.

---

## §3. Proof outline (per-lemma transport)

We show that each of L1', L2, L3, L4, L5 transports from the fixed-$T$ regime to the coupled-scaling regime without modification. The key observation: every $T$-dependent quantity in Paper 45 lives in **L4 only** (the Berezin reach), and that dependence is the standard Fejér-on-circle rate $O(T/N_t)$.

### L1' (operator-system substrate)

Paper 45 Lemma L1' (`lem:L1prime`) constructs $\OpL_{n_{\max}, N_t}$ as the (chirality-doubled scalar spatial) $\otimes$ (momentum-diagonal temporal) multiplier algebra. The construction is **stated independently of $T$** — the temporal multipliers are diagonal in the Fourier basis of $S^1_T$ at radius $T$, but the algebraic content (commutative diagonal subalgebra of $\Mat_{N_t}(\C)$) is $T$-independent. Propagation number $= 2$, Krein-positivity preservation, Riemannian-limit recovery at $N_t = 1$ all hold uniformly in $T$.

**Transport status: TRIVIAL** — L1' is structurally $T$-blind.

### L2 (joint cb-norm)

Paper 45 Lemma L2 (`lem:L2`, eq:L2_main):
$$\cbnorm{S_{\Kjoint_{n_{\max}, N_t, T}}} = \cbnorm{S_{\KSU_{n_{\max}}}} \cdot \cbnorm{S_{\KU_{N_t, T}}} = \frac{2}{n_{\max} + 1} \cdot 1 = \frac{2}{n_{\max} + 1}.$$

The U(1)_T cb-norm is **exactly $1$** at every $(N_t, T)$ — this is the headline Bożejko-Fendler simplification on the amenable compact group $\SU(2) \times U(1)_T$ (Paper 45 §5 Lemma L2 proof). The product $\cbnorm = 2/(n_{\max} + 1)$ depends only on $n_{\max}$.

**Transport status: TRIVIAL** — L2 cb-norm is $T$- and $N_t$-independent.

### L3 (joint Lichnerowicz)

Paper 45 Lemma L3 (`lem:L3`, eq:L3_struct_id):
$$[\DL, a_s \otimes a_t] = i [\DGV, a_s] \otimes a_t \quad \text{(bit-exact)}.$$

This structural identity has two vanishing-cross-term contributions:
1. $[\gamma^0, a_s] = 0$ because natural-substrate $a_s = \mathrm{blkdiag}(W, W)$ has identical chirality blocks.
2. $[\partial_t, a_t] = 0$ because temporal multiplier $a_t$ and $\partial_t = i \cdot \mathrm{diag}(\omega_k)$ are both diagonal in the momentum basis.

**Neither (1) nor (2) depends on $T$.** The Fourier momenta $\omega_k = 2\pi k/T$ scale with $T$, but the *commutator* $[\partial_t, a_t]$ vanishes structurally regardless of the specific momentum values. The joint Lichnerowicz constant
$$C_3^{\mathrm{joint}}(n_{\max}, N_t) \le C_3^{\SU(2)}(n_{\max}) = \sqrt{1 - 1/n_{\max}} \to 1^{-}$$
is **$T$- and $N_t$-independent**.

**Transport status: TRIVIAL** — L3 is the load-bearing $T$-blind ingredient. This is the central reason L3c-α works on the natural substrate (in contrast to the enlarged substrate, where the chirality-flip term $[\gamma^0, M^{\flip}] = 2 \mathrm{diag}(W, -W) \neq 0$ produces the $O(T)$ obstruction documented in β-L3 §1.3).

### L4 (joint Berezin reconstruction)

Paper 45 Lemma L4 (`lem:L4`, eq:L4_approx_id):
$$\opnorm{\Bjoint(f) - \Pjoint M_f \Pjoint} \le \gamma^{\mathrm{joint}}_{n_{\max}, N_t, T} \cdot \norm{\nabla^{\mathrm{joint}} f}_{\infty},$$
with $\gamma^{\mathrm{joint}} = \max(\gamma^{\SU(2)}, \gamma^{U(1)})$, $\gamma^{\SU(2)} = O(\log n_{\max}/n_{\max})$, $\gamma^{U(1)} = O(T/N_t)$.

**This is the load-bearing $T$-dependent ingredient.** Under admissible scaling $T(N_t)/N_t \to 0$,
$$\gamma^{U(1)}\bigl(N_t, T(N_t)\bigr) = O\bigl(T(N_t)/N_t\bigr) = o(1) \to 0.$$
The proof structure (Paper 45 §5 L4 proof, lines 965-989) is the standard Katznelson Ch. I Fejér-on-circle approximation: $\|K^{U(1)}_{N_t, T} \ast f_t - f_t\|_\infty \le C_{\mathrm{Fej}} \cdot (T/N_t) \cdot \|\partial_t f_t\|_\infty$ where $C_{\mathrm{Fej}}$ is the standard Fejér constant (numerical value $\sim 1$, independent of $T$). The implicit constant is $T$-independent because the Fejér kernel on $S^1_T$ is the standard Fejér kernel rescaled to circumference $T$, and the rescaling preserves the rate constant.

**Transport status: WORKS UNDER ADMISSIBLE SCALING.** $\gamma^{U(1)} \to 0$ iff $T(N_t)/N_t \to 0$. All other L4 properties (positivity, contractivity, L3 compatibility, Krein-positivity preservation) are $T$-blind: positivity and contractivity follow from convolution-form / UCP / Plancherel structure (independent of $T$); L3 compatibility uses L3 which is $T$-blind; Krein-positivity uses $[\JL, \Bjoint(f)] = 0$ which is structural.

### L5 (Latrémolière propinquity assembly)

Paper 45 §6 (eq:propinquity_bound):
$$\Lprop \le \max(\mathrm{reach}_B, \mathrm{reach}_P, \mathrm{height}_B, \mathrm{height}_P) = \max(\gamma^{\mathrm{joint}}, \gamma^{\mathrm{joint}}, C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}, 0) = C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}.$$

The four contributions in the max are:
- $\mathrm{reach}_B \le \gamma^{\mathrm{joint}}$ — uses L4(c) approximate identity.
- $\mathrm{reach}_P \le \gamma^{\mathrm{joint}}$ — uses L4(c) symmetric argument.
- $\mathrm{height}_B \le C_3^{\mathrm{joint}} \cdot \gamma^{\mathrm{joint}}$ — uses L3 compatibility + L4(c).
- $\mathrm{height}_P = 0$ — $P_{\mathrm{joint}}$ is a projection.

Every contribution decomposes as (a $T$-independent constant) × $\gamma^{\mathrm{joint}}$. The constants are L2 ($T$-independent), L3 ($T$-independent), or $0$. Therefore the propinquity bound at each $(n_{\max}, N_t, T(N_t))$ is structurally identical to the fixed-$T$ bound, with $\gamma^{\mathrm{joint}}$ evaluated at the coupled cell.

**Transport status: TRIVIAL given L4** — L5 is the propinquity bookkeeping; once L4 gives $\gamma^{\mathrm{joint}} \to 0$, L5 gives $\Lprop \to 0$ verbatim.

---

## §4. The load-bearing observation

The reason L3c-α closes essentially trivially on the natural substrate is captured in one sentence:

> **Every $T$-dependent quantity in the Paper 45 propinquity bound is concentrated in L4's $\gamma^{U(1)} = O(T/N_t)$ Fejér-on-circle rate; all other constants ($C_3^{\mathrm{joint}}$, cb-norm, propagation number, reach/height structure) are $T$-independent.**

This is the structural reason the natural-substrate path is the **right** opening move for G2 — it isolates the $T$-dependence in a single, controllable ingredient.

By contrast, the enlarged substrate (Paper 46 Appendix B) has an additional $O(T)$ contribution in the **gradient norm** itself (the $2 \|D_t\|_{\mathrm{op}} \|a^{\mathrm{flip}}\|_{\mathrm{op}}$ term, β-L3 §1.3), which is not paid off by the $\gamma^{U(1)}$ rate. That's why β-L3 found "rate does NOT survive joint $T \to \infty$" on the enlarged substrate — but only on the enlarged substrate. The natural substrate is structurally clean for parametric de-compactification.

---

## §5. Numerical-verification plan (optional sub-sprint L3c-α.2)

The theorem of §2 is analytical and does not require numerical verification beyond the panel cells already in Paper 45 §7 Table~1 ($(n_{\max}, N_t) \in \{(2,3), (3,5), (4,7)\}$ at $T = 2\pi$). However, a small numerical panel would confirm the rate scales correctly under coupled scaling.

**Proposed panel (if PI wants empirical confirmation):**

| $N_t$ | $T(N_t)$ at $T = N_t/\log N_t$ | $T(N_t)$ at $T = \sqrt{N_t}$ | $\gamma^{U(1)}(N_t, T)$ scaling |
|:------|:--------------------------------|:------------------------------|:--------------------------------|
| $3$ | $3/\log 3 \approx 2.73$ | $\sqrt{3} \approx 1.73$ | baseline (compare to fixed $T = 2\pi$) |
| $5$ | $5/\log 5 \approx 3.11$ | $\sqrt{5} \approx 2.24$ | $\gamma^{U(1)} \propto 1/\log N_t$ vs $1/\sqrt{N_t}$ |
| $7$ | $7/\log 7 \approx 3.60$ | $\sqrt{7} \approx 2.65$ | |
| $11$ | $11/\log 11 \approx 4.59$ | $\sqrt{11} \approx 3.32$ | |
| $21$ | $21/\log 21 \approx 6.90$ | $\sqrt{21} \approx 4.58$ | |

The numerical computation would verify $\Lprop(\mathcal{T}^L_{n_{\max}, N_t, T(N_t)}, \mathcal{T}^L_{S^3 \times S^1_{T(N_t)}})$ for $n_{\max} \in \{2, 3, 4\}$ and these $N_t$ values, using the existing `geovac/gh_convergence_tensor.py` infrastructure (verified at Paper 45 panel) with the temporal-radius parameter swept along the coupled sequence.

**Status: DEFERRED to optional sub-sprint L3c-α.2.** The analytical statement (§2) is the load-bearing deliverable; the numerical sweep is confirmation that would not change the verdict.

---

## §6. Decision gates and downstream sprints

### Decision gate 1 (this sprint outcome)

L3c-α closes **POSITIVE** at the analytical level on the natural substrate (Paper 45 / Paper 46 main theorem). The result is a one-sentence theorem (§2) with a per-lemma transport audit (§3) that adds no new analytical content — Paper 45's proof of the main theorem applies verbatim at each coupled cell.

### Decision gate 2 (open L3c-β?)

L3c-β would extend the same parametric statement to the **enlarged substrate** (Paper 46 Appendix B). The β-L3 finding documents that the rate fails on the enlarged substrate under joint $T \to \infty$ at the current gradient form. L3c-β would require either:

- **Option (i):** A refined gradient norm that doesn't pick up the $O(T)$ term from the chirality-flip contribution. Candidates: replace $\norm{a^{\flip}}_{\mathrm{op}}$ with a $T$-modulated version, e.g., $T^{-1/2} \norm{a^{\flip}}_{\mathrm{op}}$ scaled so the Lichnerowicz bound becomes $T$-independent.
- **Option (ii):** A different propinquity construction that absorbs the $O(T)$ term into a height contribution rather than the rate (β-L3 §5.4 option (ii)).

**Estimated cost for L3c-β: 6–10 weeks** (scoping memo §5 estimate). Risk MEDIUM — original work, may terminate as partial negative.

**Recommendation: DEFER L3c-β.** L3c-α closes the natural-substrate case and is sufficient for the "stepping stone" role toward L3c-γ. L3c-β can be deferred until either (a) PI explicitly wants the enlarged-substrate parametric statement or (b) L3c-γ surfaces a need for it.

### Decision gate 3 (open L3c-γ?)

L3c-γ would construct the **inductive limit** from the sequence of compact-temporal triples $\{\mathcal{T}^L_{S^3 \times S^1_{T_n}}\}_{n \in \mathbb{N}}$ (with $T_n \to \infty$ along an admissible scaling) to the genuine non-compact triple $\mathcal{T}^L_{S^3 \times \mathbb{R}_t}$.

**Construction outline** (per scoping memo §4 Path P3):
1. Verify that the spectral triples $\{\mathcal{T}^L_{S^3 \times S^1_{T_n}}\}$ form an inductive system: refinement maps $\mathcal{T}^L_{S^3 \times S^1_{T_n}} \hookrightarrow \mathcal{T}^L_{S^3 \times S^1_{T_{n+1}}}$ via Fourier-mode extension when $T_{n+1} > T_n$.
2. Verify that the inductive limit (in the operator-system / spectral-triple category) is the right Krein triple on $S^3 \times \mathbb{R}$ — natural candidate is Paper 43's bounded-interval Krein construction $\mathcal{T}^L_{S^3 \times [-T_{\max}, T_{\max}]}$ with $T_{\max} \to \infty$.
3. Bound the embedding propinquity $\Lprop(\mathcal{T}^L_{S^3 \times S^1_{T_n}}, \mathcal{T}^L_{S^3 \times \mathbb{R}_t})$ via tail estimates as $T_n \to \infty$.

L3c-α provides the building block for step 1 (the parametric stability statement). Steps 2–3 are the new analytical content.

**Estimated cost for L3c-γ: 6–10 weeks** (scoping memo §5). Risk MEDIUM-HIGH — multiple delicate steps, no direct precedent for inductive-limit propinquity in the Lorentzian / Krein setting.

**Recommendation: OPEN L3c-γ as the natural next sprint.** This is where the substantive new content for G2 sits — L3c-α is essentially a corollary of Paper 45.

---

## §7. Honest scope summary

| Question | Status |
|:---------|:-------|
| Does the Paper 45 main-theorem rate survive joint $(n_{\max}, N_t, T) \to \infty$ on the natural substrate with admissible scaling $T(N_t)/N_t \to 0$? | **YES** (§2 Theorem, §3 per-lemma transport) |
| Does this give the genuine non-compact propinquity $\Lprop(\mathcal{T}^L_{n_{\max}, N_t, T(N_t)}, \mathcal{T}^L_{S^3 \times \mathbb{R}_t}) \to 0$? | **NO** — gives only the parametric statement with compact-temporal target at each finite stage. The genuine non-compact limit needs L3c-γ. |
| Does this close G2 (Paper 45 §1.4)? | **PARTIAL.** G2 asks for the non-compact $\R_t$ limit; L3c-α establishes parametric stability under $T \to \infty$, which is a necessary but not sufficient condition. |
| Does L3c-α produce a publishable new theorem? | **PROBABLY NOT STANDALONE.** The theorem follows from Paper 45 main theorem essentially as a corollary. The right venue is either (a) a corollary/remark in Paper 45 ("the convergence is uniform along admissible scaling sequences"), or (b) a stepping-stone lemma in the eventual Paper 47 / 48 covering L3c-γ. |
| Does the L3c-α result open the L3c-γ sprint? | **YES.** L3c-α is the analytical foundation for the parametric stability needed in L3c-γ's inductive limit. |

---

## §8. Recommended paper edits (NOT applied; PI decision)

If the PI accepts L3c-α as closed at the analytical level, three paper edits are natural:

1. **Paper 45 §6.2 or §7.2** — new Corollary or Remark stating the parametric stability:
   > *Corollary (parametric de-compactification).* Let $T(\cdot)$ satisfy $T(N_t) \to \infty$ and $T(N_t)/N_t \to 0$ as $N_t \to \infty$. Then under Theorem~\ref{thm:main}, $\Lprop(\mathcal{T}^L_{n_{\max}, N_t, T(N_t)}, \mathcal{T}^L_{S^3 \times S^1_{T(N_t)}}) \to 0$ as $(n_{\max}, N_t) \to (\infty, \infty)$, with $\gamma^{\mathrm{joint}} = \max(O(\log n_{\max}/n_{\max}), O(T(N_t)/N_t)) \to 0$ from L4 and $C_3^{\mathrm{joint}}$ $T$-independent from L3.

2. **Paper 45 §1.4 G2 statement** — refine to clarify that G2 is the *inductive-limit* question (L3c-γ), not the parametric statement (L3c-α). One-sentence update.

3. **Paper 46 Appendix B / §7.2** — cross-reference: "the L3c-α parametric result extends bit-exact to the strong-form via the free-upgrade reading; on the enlarged substrate of Appendix B, the result does NOT extend (β-L3 §1.3 obstruction)."

These edits are scope-honest (don't claim more than the theorem gives) and naturally fit the existing paper structure. **They are NOT applied in this memo** — queued for PI decision per §9 below.

---

## §9. Open questions for PI

1. **Accept L3c-α as closed?** The analytical statement (§2) follows from Paper 45 with no new content; the per-lemma audit (§3) is bookkeeping. Is this sufficient to call L3c-α CLOSED, or does the PI want the numerical-verification panel (§5) before declaring closure?
2. **Apply paper edits (§8)?** Three edits are scope-honest and natural fits to Paper 45 / Paper 46. Apply, or defer until L3c-γ also lands?
3. **Open L3c-γ as next sprint?** This is the substantive new content for G2 (inductive limit to genuine $\mathbb{R}_t$). 6–10 weeks estimate. Alternative: arXiv submission of Papers 45/46 first, then L3c-γ.
4. **L3c-β positioning.** Defer indefinitely, or queue as follow-on after L3c-γ lands?

---

## §10. Sprint verdict

**Sprint L3c-α: CLOSED at the analytical level on the natural substrate (Paper 45 / Paper 46 main theorem).**

- Theorem statement (§2): parametric propinquity convergence under admissible scaling.
- Per-lemma transport audit (§3): every $T$-dependent quantity isolated in L4's $\gamma^{U(1)} = O(T/N_t)$, vanishes under $T/N_t \to 0$.
- Load-bearing observation (§4): natural-substrate $C_3^{\mathrm{joint}}$ and L2 cb-norm are structurally $T$-independent because the time-direction Lipschitz invisibility ($\{\gamma^0, \DGV\} \neq 0$ but $[\gamma^0, M^{\spat}] = 0$) leaves no $T$-residue in the bound apart from the Fejér rate.
- Numerical verification (§5): deferred to optional L3c-α.2, would confirm but not change the verdict.

**Honest scope (§7):** L3c-α gives the parametric statement, not the genuine non-compact limit. G2 closure requires L3c-γ.

**Next sprint recommendation (§6 Decision gate 3):** **OPEN L3c-γ** (inductive limit to genuine $\mathbb{R}_t$) as the substantive next sprint for G2. 6–10 weeks estimate, MEDIUM-HIGH risk, leverages L3c-α as analytical foundation.

**Confidence:** HIGH on the analytical statement (§2). HIGH on the per-lemma transport audit (§3) — every step is a direct application of Paper 45 proofs with $T$ varying. MEDIUM on whether this should be presented as a standalone result or absorbed as a corollary in Paper 45 / a stepping stone in Paper 47 (open question for PI).
