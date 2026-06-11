# Casimir triangle inequality: literature audit and verdict

**Date:** May 2026
**Sprint:** Unified compact-Lie-group propinquity scoping (follow-on to Paper 38)
**Status:** **VERDICT: FALSE.** Counterexample found at SU(2) and confirmed at SU(3).

**Erratum 2026-05-15:** Throughout this memo, references to "Hekkelman 2023" / "Hekkelman, *Convergence of spectral truncations for compact metric groups*" at arXiv:2310.14733 are **incorrect attribution**. The correct authors of arXiv:2310.14733 (verified against arXiv abstract page + IMRN published version) are **Yvann Gaudillot-Estrada and Walter D. van Suijlekom**, published as "Convergence of Spectral Truncations for Compact Metric Groups," International Mathematics Research Notices, Volume 2025, Issue 13, July 2025, paper ID rnaf197, DOI 10.1093/imrn/rnaf197. The misattribution was likely a conflation with Hekkelman 2022 (the $S^1$ Toeplitz paper, arXiv:2111.13865) or Hekkelman–McDonald 2024 (the $T^d$ Berezin–Toeplitz paper, arXiv:2403.18619 / 2412.00628), which are different papers on different spaces. The substantive findings of this memo (the Casimir-triangle inequality being false, the corrected Dirac-triangle reformulation being right) are independent of the attribution.

## 1. Statement

Let $G$ be a compact connected Lie group, $\pi, \pi'$ irreps with highest weights $\lambda, \lambda'$, and $\sigma$ any irrep appearing in $\pi \otimes \pi'^*$. Question:

$$|C(\pi) - C(\pi')| \le C(\sigma)? \tag{*}$$

with $C(\pi) = \langle \lambda + \rho, \lambda + \rho \rangle - \langle \rho, \rho\rangle$.

## 2. The prompt-claimed SU(2) verification is incorrect

The prompt asserts: *"For SU(2), this inequality reduces to $|j(j+1) - j'(j'+1)| \le J(J+1)$ for $J \in [|j-j'|, j+j']$, which is the Avery 3-Y selection rule and is easily verified."*

This conflates two different inequalities. The Avery 3-Y rule is the **triangle rule on $J$**: $|j-j'| \le J \le j+j'$. That is unconditionally true (it defines the range). The Casimir inequality (*) is a separate claim, and it fails.

**SU(2) counterexample.** Take $j = 1$, $j' = 1/2$.
- $C(j) = 1 \cdot 2 = 2$
- $C(j') = (1/2)(3/2) = 3/4$
- $|\Delta C| = 5/4$
- Tensor product: $V_1 \otimes V_{1/2} = V_{1/2} \oplus V_{3/2}$
- For $\sigma = V_{1/2}$: $C(\sigma) = 3/4$
- $|\Delta C| = 5/4 > 3/4 = C(\sigma)$. **Inequality fails.**

**Algebraic explanation in SU(2).** With $j > j' \ge 0$:
$$C(j) - C(j') = j(j+1) - j'(j'+1) = (j - j')(j + j' + 1).$$
The smallest Casimir in $V_j \otimes V_{j'}$ is at $J_{\min} = j - j'$:
$$C(J_{\min}) = (j - j')(j - j' + 1).$$
The inequality $C(j) - C(j') \le C(J_{\min})$ requires
$$(j - j')(j + j' + 1) \le (j - j')(j - j' + 1) \iff j + j' + 1 \le j - j' + 1 \iff 2 j' \le 0,$$
which fails for any $j' > 0$.

**The inequality holds only when $j'(j'+1) = 0$**, i.e., $j' = 0$ (trivial rep), where $C(\pi') = 0$ and $V_{j'}^* \otimes V_j = V_j$ contains $\sigma = V_j$ with $C(\sigma) = C(j) \ge |\Delta C|$ trivially.

## 3. SU(3) confirms the failure broadly

Numerical verification (`debug/casimir_triangle_su3.py`, sympy exact arithmetic) over all pairs $(\lambda, \lambda')$ with $|\lambda|, |\lambda'| \le 3$ in Dynkin labels:

| Total tested | Pass | Fail |
|---|---|---|
| 100 pairs | 68 | **32** |

**Headline counterexample.** $\lambda = (3,0)$ (decuplet), $\lambda' = (1,0)$ (fundamental).
- $C(3,0) = 6$, $C(1,0) = 4/3$, $|\Delta C| = 14/3 \approx 4.67$
- $\lambda'^* = (0,1)$, decomposition $(3,0) \otimes (0,1) = (3,1) \oplus (2,0)$ (the standard $\mathbf{10} \otimes \bar{\mathbf{3}} = \bar{\mathbf{6}} \oplus \mathbf{24}$, dimension $30 = 10 \cdot 3$ ✓)
- For $\sigma = (2,0)$: $C(\sigma) = 10/3 \approx 3.33$
- $|\Delta C|/C(\sigma) = 14/10 = 1.40$. **Violated by 40%.**

**Pattern of failures.** Every pair where one of $\lambda, \lambda'$ is far above the other in the dominance order produces a violation at the smallest-Casimir component of the decomposition. The largest violation ratio in the panel is **2.0**, e.g.:

- $\lambda = (0,2), \lambda' = (0,3)$: $|\Delta C| = 8/3$, $\sigma = (1,0)$ with $C(\sigma) = 4/3$, ratio 2.0
- $\lambda = (3,0), \lambda' = (2,0)$: $|\Delta C| = 8/3$, $\sigma = (1,0)$ with $C(\sigma) = 4/3$, ratio 2.0

The 8-of-8 sanity (8 ⊗ 8 = 1+8+8+10+$\bar{10}$+27 ✓) and 3 ⊗ $\bar 3$ = 1+8 sanity passed, so the tensor product code is correct.

## 4. Literature audit

Searches performed:

1. "Casimir eigenvalue inequality tensor product representation compact Lie group"
2. "Casimir tensor product highest weight inequality bound"
3. "Brauer-Steinberg formula tensor product highest weights compact Lie group"
4. "Casimir difference Cauchy-Schwarz weight lattice Lie algebra"
5. "Casimir tensor product spectral truncation Connes spectral triple"
6. "Lipschitz constant compact Lie group Connes distance representation Casimir"
7. "Rieffel coadjoint orbit Casimir Lipschitz seminorm compact Lie group"
8. "|C(λ) - C(μ)| Casimir highest weight bound"
9. "Hawkins Berezin reconstruction compact Lie group spectral truncation"
10. "Steinberg formula highest weight tensor product convex hull dominant chamber"
11. "length function compact Lie group Casimir Lip-norm Connes spectral length"

Specific papers fetched:

- **Connes–van Suijlekom, *Spectral truncations in NCG and operator systems*** (arXiv:2004.14115): no Casimir inequality used.
- **Rieffel, *Dirac operators for coadjoint orbits of compact Lie groups*** (arXiv:0812.2884): no Casimir-difference inequality. Uses Levi-Civita connection.
- **van Suijlekom et al., *GH convergence of state spaces for spectral truncations*** (arXiv:2005.08544): focused on circle/torus, no Casimir inequality.
- **Leimbach–van Suijlekom, *GH convergence of spectral truncations for tori*** (arXiv:2302.07877, *Adv. Math.* 2024): proof goes through Schur–Fourier multipliers and a spectral Fejér kernel; no Casimir inequality.
- **Hekkelman, *Convergence of spectral truncations for compact metric groups*** (arXiv:2310.14733): generalizes Leimbach–vS to general compact metric groups via Lip-norms induced by group action; **no Casimir-difference inequality used** in the proof. The Lip-norm framework operates without ever bounding $|C(\pi) - C(\pi')|$ by $C(\sigma)$.

**Verdict on the literature audit.** The inequality (*) does not appear as a standard or named result in compact Lie group representation theory, in Rieffel's coadjoint-orbit Dirac operator work, in the Connes–van Suijlekom spectral truncation framework, or in the GH-convergence sequence (Leimbach–vS, Hekkelman). The closest related result is the elementary fact that $C(\sigma) \ge 0$ for every $\sigma$, which is trivial.

The inequality (*) **could** have been a candidate for the L3 Lipschitz step in the unified theorem, but it is not present in the extant proofs because **it is false**. The published proofs (Leimbach–vS, Hekkelman) bypass it entirely by working with Schur multipliers / matrix coefficient norms, which give a completely different route.

## 5. The correct (weaker) bound that *does* hold

The substitute that **does** hold is:

$$\boxed{|C(\pi) - C(\pi')| \le C(\pi) + C(\pi') = C(\pi) + C(\pi'^*)}$$

This is the SU(2) Avery rule (the **triangle inequality on $J$**, not on Casimirs) plus the elementary observation that $C(\pi) \le C(\pi \otimes \pi'^*)$ in some averaged sense, but even that needs care.

A cleaner correct statement: **in $V_\pi \otimes V_{\pi'}^*$, the largest Casimir component is bounded above by $C(\pi) + C(\pi'^*) + 2 \langle \lambda, \lambda' \rangle$**. The largest highest weight is $\lambda + \lambda'^*$, giving $C(\sigma_{\max}) = \langle \lambda + \lambda'^* + \rho, \lambda + \lambda'^* + \rho\rangle - \langle \rho, \rho\rangle$. By Cauchy–Schwarz, $|C(\pi) - C(\pi')| \le C(\sigma_{\max})$ when $\sigma_{\max} = $ Cartan product is contained.

So the inequality **holds for $\sigma$ = the Cartan-product (largest) summand**, but **fails for the smallest summand**, which is what (*) requires. In the SU(2) case the failure is at $\sigma = V_{|j-j'|}$; in SU(3) at the smallest-Casimir $\sigma$ of the decomposition.

## 6. Why the proof sketch breaks down

Following the prompt's proof attempt:

1. **Brauer–Steinberg.** A summand $\sigma$ of $\pi \otimes \pi'^*$ has highest weight of the form $\mu = \lambda - w \lambda'$ for some $w \in W$ (Weyl group). Specifically, $\mu \in \lambda + \mathrm{wts}(V_{\pi'^*})$ shifted into the dominant chamber.

2. **Casimir difference identity.**
$$C(\pi) - C(\pi') = \langle \lambda - \lambda', \lambda + \lambda' + 2\rho\rangle.$$

3. **Casimir of $\sigma$.** $C(\sigma) = \langle \mu, \mu + 2\rho\rangle$.

4. **Where the proof fails.** In the SU(2) case, the smallest summand has $\mu = (j - j') \cdot \omega$ in fundamental-weight notation. With $\rho = \omega$ for SU(2),
   $$C(\sigma) = (j-j')(j-j'+1)$$
   while
   $$C(\pi) - C(\pi') = (j-j')(j+j'+1)$$
   and the difference $(j-j')(j+j'+1) - (j-j')(j-j'+1) = (j-j') \cdot 2 j' \ge 0$ with equality only when $j j' = 0$.

In Cauchy–Schwarz language: $\langle \lambda - \lambda', \lambda + \lambda' + 2\rho\rangle$ has TWO contributions, $\langle \lambda - \lambda', \lambda - \lambda'\rangle + \langle \lambda - \lambda', 2 \lambda' + 2\rho\rangle$. The first is the size of $\lambda - \lambda'$, but the **second** can be large and positive (when $\lambda - \lambda'$ aligns with $\rho$), and there is no Cauchy–Schwarz wrinkle that bounds it by $C(\sigma)$.

Concretely: if $\lambda$ and $\lambda'$ are colinear in the dominant chamber (e.g., $\lambda = 3\omega_1, \lambda' = \omega_1$ in SU(3)) the inner-product cross-terms add constructively, blowing up $|\Delta C|$ relative to the smallest $C(\sigma)$.

## 7. Implications for the unified theorem L3 lemma

Paper 38's L3 Lipschitz comparison constant $C_3 = 1$ in SU(2) was obtained by a route that doesn't actually use (*). It uses the **Avery 3-Y selection rule** ($|j-j'| \le J \le j+j'$) plus an explicit harmonic gradient bound on $S^3$ (the natural-panel calculation). The "Casimir triangle inequality" framing in the scoping memo is a **mis-summary** of what L3 actually proved.

For the unified theorem, the correct statement of L3 is:

> $\|[D, B_{n_{\max}}(f)]\|_{op} \le C_3(G) \cdot \|\nabla f\|_\infty$

where $C_3(G)$ depends on the **gradient operator on $G$** and the **harmonic-analysis dual of the Berezin reconstruction**, not on a Casimir inequality. The right object is the matrix coefficient $\langle \nabla \chi_\pi, \nabla \chi_{\pi'} \rangle$ on $G$, which is bounded by Cauchy–Schwarz in the standard way.

**Net effect on the scoping memo:** the L3 step in the unified theorem does **not** require (*) and was wrongly identified as the load-bearing item. The actual load-bearing items are:

- L2 (central spectral Fejér kernel): proven for SU(2), the SU(N) generalization needs the fundamental-domain volume and the abelianizing assumption (central kernel).
- L4 (Berezin reconstruction): the Hawkins 2000 formalism on Kähler-quantizable coadjoint orbits gives this for all coadjoint-orbit cases. For non-Kähler $G$ (most compact Lie groups), Hawkins-style reconstruction needs the Cartan-subgroup decomposition, which does NOT generalize trivially.

The cleanest path forward is to follow Hekkelman 2023 (compact metric groups) directly, since that proof goes through Schur-multiplier estimates and Lip-norms induced by the group action, and **does not require any Casimir-difference inequality**.

## 8. Recommendation

1. **Update the unified GH scoping memo** (`debug/unified_gh_scoping_memo.md`): retract the "Casimir triangle inequality" framing as the load-bearing L3 lemma. The real L3 obstacle (if any) is the harmonic-gradient bound on $G$, which Hekkelman 2023 already handles via Lip-norms induced by group translation.

2. **Use Hekkelman 2023 (arXiv:2310.14733) as the master reference** for the unified theorem rather than re-proving via Berezin / central-Fejér transport. That paper already proves the convergence theorem for **all compact metric groups**, which is strictly stronger than what's needed here.

3. **The rate constant $c(G) = 2 \mathrm{Vol}(G/T)/\mathrm{Vol}(G)$ as the M1 Hopf-base measure signature** is a separate question from the Casimir inequality; it does NOT depend on (*). The constant is a pure geometric invariant of the group; verifying it is the M1 signature requires running the master-Mellin-engine analog of MR-A/MR-B from the SU(2) sprint at higher rank.

4. **No revision needed to Paper 38.** Paper 38's L3 proof is for SU(2) only and uses an explicit harmonic-gradient bound, not the (false) general inequality. The SU(2) panel calculation $C_3 = 1$ stands.

## 9. Confidence level

**HIGH** that the inequality is FALSE.
- Algebraic argument in SU(2) is one line: $(j-j')(j+j'+1) > (j-j')(j-j'+1)$ for $j' > 0$.
- Numerical scan in SU(3) found 32 / 100 violations, with maximum ratio 2.0.
- Literature audit confirms (*) is not a known theorem.

**HIGH** that Paper 38's L3 proof for SU(2) does not actually depend on the false inequality. The natural-panel calculation in `geovac/central_fejer_su2.py` (L3 lemma module) computes $C_3$ via direct evaluation of the SU(2) gradient on harmonic polynomials, not via (*).

**MODERATE** that the unified theorem is reachable via the Hekkelman 2023 route, which sidesteps the Casimir-difference question entirely.

## 10. Files

- `debug/casimir_triangle_su3.py` — sympy exact-arithmetic SU(3) verifier (Freudenthal weight enumeration + Brauer–Klimyk tensor product decomposition + Casimir inequality check; sanity verified on 3⊗3̄, 3⊗3, 8⊗8).
- `debug/casimir_triangle_inequality_lit_check.md` — this memo.

## 11. Sources

- Hekkelman, *Convergence of spectral truncations for compact metric groups*, arXiv:2310.14733 (2023).
- Leimbach–van Suijlekom, *Gromov-Hausdorff convergence of spectral truncations for tori*, *Adv. Math.* 439 (2024), arXiv:2302.07877.
- Connes–van Suijlekom, *Spectral truncations in noncommutative geometry and operator systems*, *Comm. Math. Phys.* 383 (2021), arXiv:2004.14115.
- Rieffel, *Dirac operators for coadjoint orbits of compact Lie groups*, arXiv:0812.2884; expanded as *Dirac operators for matrix algebras converging to coadjoint orbits*, *Comm. Math. Phys.* 401 (2023).
- van Suijlekom et al., *Gromov-Hausdorff convergence of state spaces for spectral truncations*, *J. Geom. Phys.* 162 (2021), arXiv:2005.08544.
- Hawkins, *Berezin–Toeplitz quantization on Lie groups*, *J. Funct. Anal.* 2009, arXiv:0806.3063.
