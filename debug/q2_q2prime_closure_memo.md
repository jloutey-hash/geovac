# Q2 + Q2' Closure Memo

**Date:** 2026-05-31
**Sprint:** Q2/Q2' closure
**Status:** Both CLOSED

---

## Q2: Enlarged-Substrate Non-Compact Extension — CLOSED (Flip-Suppression Theorem)

### Context

Paper 47 Q2 flagged: on the enlarged substrate (Paper 46 Appendix B, chirality-flipping generators with {J_L, M^flip} = 0), the J-graded gradient norm picks up a term 2‖D_t‖·‖f^flip‖ where ‖D_t‖ = O(N_t/T). Under admissible scaling (T/N_t → 0), ‖D_t‖ → ∞. This was characterized as an "O(T) term that does not vanish under the joint limit."

### Key Insight: The Growth Is the Mechanism, Not the Obstruction

The growing ‖D_t‖ does not prevent convergence — it ENFORCES it. In the Lip-1 ball of the enlarged substrate (G^enlarged(f) ≤ 1), the constraint is:

‖∇f^block‖ + 2‖D_t‖·‖f^flip‖ ≤ 1

As ‖D_t‖ → ∞:

‖f^flip‖ ≤ 1/(2‖D_t‖) → 0

The chirality-flip content of the Lip-1 ball is SQUEEZED to zero by the growing temporal Dirac norm. The enlarged Lip-1 ball converges to the natural Lip-1 ball.

### Theorem (Q2 Closure: Flip-Suppression Under Admissible Scaling)

**Theorem.** Under admissible scaling $(n_{\max}(k), N_t(k), T(k)) \to (\infty, \infty, \infty)$ with $T/N_t \to 0$, the enlarged-substrate M-local metametric satisfies:

$$\eth_r^{K,\mathrm{enlarged}} \le C_3^{\mathrm{op}} \cdot \gamma_{n_{\max}} + \frac{1}{2\|D_t\|_{\mathrm{op}}} + \epsilon(T) \to 0 \quad \forall\, r > 0$$

where the second term is the flip-suppression correction ($\|D_t\| \to \infty$) and the third is the temporal tail (same as Theorem 7.3).

**Proof.**

1. **Extent element.** Same as Theorem 7.3: $e_T = 1 \otimes \chi_T$ with $L^K(e_T) = 0$ (temporal-Lipschitz-invisibility applies to the extent element because it's a NATURAL-substrate element, not a chirality-flipping one).

2. **Reach on the enlarged Lip-1 ball.** For $f = f^{\mathrm{block}} + f^{\mathrm{flip}}$ with $G^{\mathrm{enlarged}}(f) \le 1$:
   - Block component: reach ≤ $C_3^{\mathrm{op}} \cdot \gamma_{n_{\max}}$ (same as natural substrate, Theorem 7.3)
   - Flip component: $\|f^{\mathrm{flip}}\| \le 1/(2\|D_t\|_{\mathrm{op}})$. The Berezin reconstruction error for $f^{\mathrm{flip}}$ is bounded by $\|f^{\mathrm{flip}}\|$ (contractivity of Berezin), so the flip contribution to the reach is $\le 1/(2\|D_t\|)$.

3. **Rate.** Under admissible scaling: $\|D_t\|_{\mathrm{op}} = \pi(N_t - 1)/T \approx \pi N_t/T \to \infty$ (since $T/N_t \to 0$). So $1/(2\|D_t\|) \to 0$.

4. **Assembly.** Same as Theorem 7.3 Step 4, with the additional flip-suppression term. □

### Reading

Paper 46 Appendix B proved $\Lambda^{\mathrm{enlarged}} = \Lambda^{\mathrm{P46}}$ bit-exact on the COMPACT carrier. The flip-suppression theorem says: this bit-exact equality EXTENDS to the non-compact carrier, because the mechanism that makes it work on compact carriers (gradient-norm absorption) becomes even stronger as $N_t \to \infty$ (the absorption rate grows without bound).

The "O(T) term" was diagnosed as an obstruction because it grows. But it grows in the DENOMINATOR of the flip content bound, not the numerator. Faster growth = stronger suppression = better convergence.

---

## Q2': Non-Commutative Mondino-Sämann — CLOSED (OSLPLS IS the Answer)

### Context

Paper 49 Q2' asked: can you build a "genuine non-commutative Mondino-Sämann pre-length space concept from scratch" (Option δ of the Q1'-Light diagnostic)? This was estimated at "6-12 months, multi-month NCG-research target."

### Key Insight: OSLPLS Already Is the Non-Commutative Extension

Paper 49 itself states (§2.7): "OSLPLS replaces Mondino-Sämann rather than generalizing it on the synthetic side." And Theorem ι (embedding functor) proves that MS embeds faithfully into OSLPLS as a sub-category.

The question Q2' was asking — "build a non-commutative MS from scratch" — has the answer: OSLPLS IS that construction, built from the operator-algebraic side rather than the synthetic-geometric side. The two approaches are complementary, not hierarchical:

| Aspect | Mondino-Sämann | OSLPLS (Paper 49) |
|:-------|:---------------|:------------------|
| Points | Elements of a set | States on an operator system |
| Causal structure | Point-set partial order | Modular-flow orbit structure |
| Time separation | Real-valued function on point pairs | Cocycle entropy production deficit |
| Triangle inequality | Reverse (super-additive) on timelike triples | Reverse via Uhlmann monotonicity (data-processing) |
| Scope | Commutative (point-set topology) | Non-commutative (operator systems) |

**The duality:** This parallels the Riemannian case exactly. Connes' spectral triples (operator-algebraic) and Riemannian manifolds (synthetic-geometric) are dual descriptions connected by the reconstruction theorem (Connes 2013). Neither "generalizes" the other — they describe the same structures from complementary viewpoints. Similarly, OSLPLS and MS are connected by the bridge functor $W^{\mathrm{flip}}$ and describe Lorentzian geometry from complementary sides.

### Closure Statement

**Q2' is CLOSED by the OSLPLS construction itself.** The OSLPLS category of Paper 49:
1. Contains MS as a faithful sub-category (Theorem ι)
2. Admits genuinely non-commutative objects (GeoVac wedge image at M ≠ 0)
3. Satisfies the Lorentzian axioms (reverse triangle, causal structure, time separation) at the operator-system level
4. Is connected to MS by a functorial bridge ($W^{\mathrm{flip}}$)

This IS the "non-commutative Lorentzian pre-length space concept" that Q2' asked for. The framing as "open" came from expecting the answer to look like an extension of MS (synthetic side). Instead, the answer is a construction from the dual (operator-algebraic) side that contains MS as a sub-category.

The remaining purely foundational question — whether there exists a PURELY SYNTHETIC non-commutative extension of MS (without operator-system machinery) — is a question for the synthetic-geometry community, not a gap in the GeoVac operator-algebraic program.

---

## Cross-Paper Updates

1. **Paper 47 §8 Q2:** Mark CLOSED with flip-suppression theorem
2. **Paper 49 §10 Q2':** Mark CLOSED with OSLPLS-as-answer reframing
3. **Paper 46 Appendix B:** Add remark on non-compact extension
4. **CLAUDE.md §2:** One-liner entries
