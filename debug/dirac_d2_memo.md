# Track D2 — Dirac analog of B = 42

**Date:** 2026-04-15
**Sector sprint:** Dirac-on-S³ Tier 1
**Author:** Worker agent (D2)
**Status:** Complete
**Verdict for D5:** **B-doesn't-lift** (nuanced — see §4)

---

## 1. Summary

Paper 2's truncated Casimir trace on the scalar Laplace–Beltrami spectrum of S³ at the m=3 cutoff is
$$
B = \sum_{n=1}^{3}\sum_{l=0}^{n-1} (2l+1)\,l(l+1) = 42 \qquad \Bigl(=\tfrac{m(m-1)(m+1)(m+2)(2m+1)}{20}\Bigr|_{m=3}\Bigr).
$$

**The question:** Does any analogous finite truncated trace on the Dirac spectrum of S³ (Camporesi–Higuchi: $|\lambda_n| = n + 3/2$, $g_n^{\text{Dirac}} = 2(n+1)(n+2)$, $g_n^{\text{Weyl}} = (n+1)(n+2)$) reproduce $B$ as a closed-form identity, promoting $B$ to a sector-independent $S^3$ invariant?

**The answer:** No, but three partial hits are recorded symbolically.

All computations were pure sympy (zero floating point). Algebraic content: the Dirac sector is in $\mathbb{Q}$ (Track D1's π-free certificate), so every quantity below is an exact rational — no exchange-constant content in Paper 18's taxonomy (intrinsic/calibration/embedding/flow) enters this track.

---

## 2. What was tested (pure symbolic)

Using the CH cutoff $n_{\text{CH}} \in \{0,1,2\}$ to match Paper 2's $m=3$ (i.e.\ $n_{\text{Fock}} \le 3$):

### Candidate (a) — truncated Casimir traces

| Weight | Dirac sum | Weyl sum | Match to B=42 |
|---|---|---|---|
| $\sum g_n$ | 40 | 20 | neither |
| $\sum |\lambda_n| \cdot g_n$ | 120 | 60 | neither (60 = 10B/7) |
| $\sum |\lambda_n|^2 \cdot g_n$ | **378 = 9·B** | 189 = 9B/2 | Dirac hits 9·B |
| $\sum (2n+3)^2 \cdot g_n$ | **1512 = 36·B** | 756 = 18·B | integer multiples |

### Candidate (b) — single-level coincidences at n=5
$g_5^{\text{Weyl}} = (5+1)(5+2) = 42 = B$ and $g_5^{\text{Dirac}} = 84 = 2B$. **Outside** Paper 2's cutoff window; no independent argument selects $n=5$ from packing data $(d_{\max}, N_{\text{init}}, l_{\max}) = (4,2,2)$.

### Candidate (c) — cumulative sums
Closed forms:
$$\sum_{n=0}^{N} g_n^{\text{Dirac}} = \tfrac{(N+1)(N+2)(2N+3)}{3},\qquad \sum_{n=0}^{N} g_n^{\text{Weyl}} = \tfrac{(N+1)(N+2)(2N+3)}{6}.$$
Dirac values at $N=0,1,2,3,\dots$: 4, 16, **40**, 80, 140, 224, … Weyl: 2, 8, 20, 40, 70, **112**, 168, 240, … — neither cumulative sum hits 42 at any integer cutoff. At $N=2$, the Dirac cumulative = $\boxed{40 = \Delta^{-1}}$ (Phase 4H SM-D, already documented in CLAUDE.md §2).

### Candidate (d) — variant weightings and small-integer combinations
No $(c_0, c_1, c_2)$ with $|c_i| \le 3$ gives $c_0 g_0^{\text{Dirac}} + c_1 g_1^{\text{Dirac}} + c_2 g_2^{\text{Dirac}} = 42$ (Dirac degeneracies 4, 12, 24 are all $\equiv 0 \pmod 4$, but 42 is not). Five combinations work for the Weyl (degeneracies 2, 6, 12), but none with coefficients $(c_n) = $ something structurally meaningful (e.g., matching Paper 2's $(2l+1) l(l+1)$).

### Candidate (e) — closed-form comparison (the headline)

**The single clean single-level identity** (at the *top* of the m-window, $n_{\text{CH}} = m-1$):
$$
|\lambda_{m-1}|\cdot g_{m-1}^{\text{Weyl}} \;=\; \tfrac{(2m+1)\,m(m+1)}{2}.
$$
Compare to $B(m) = \tfrac{m(m-1)(m+1)(m+2)(2m+1)}{20}$:
$$
\frac{B(m)}{|\lambda_{m-1}|\cdot g_{m-1}^{\text{Weyl}}} \;=\; \frac{(m-1)(m+2)}{10}.
$$
This equals 1 iff $(m-1)(m+2) = 10$, solved uniquely by $\boxed{m = 3}$.

At $m=3$ exactly: $|\lambda_2| \cdot g_2^{\text{Weyl}} = (7/2)\cdot 12 = 42 = B$. A genuine symbolic equality at Paper 2's cutoff, and at no other integer $m$.

Similarly, $\sum |\lambda_n|^2 g_n^{\text{Dirac}} = (N+1)(N+2)(N+3)(2N+3)(2N+5)/10$ at $N=m-1$ has ratio $2(2m+3)/(m-1)$ to $B(m)$, equal to 9 uniquely at $m=3$.

---

## 3. Why the partial hits are coincidences, not identities

Paper 2's $B(m) = m(m-1)(m+1)(m+2)(2m+1)/20$ carries a factor $(m-1)(m+2)$. The Dirac/Weyl cumulative traces $(N+1)(N+2)(2N+3)/3$ etc.\ carry **no such factor**. The mismatch is structural:

- The scalar Laplacian on $S^3$ has a zero mode at $n=1, l=0$ (constant function). Its Casimir weight $l(l+1)$ vanishes, giving $B(1)=0$ and the $(m-1)$ factor.
- The Dirac operator on $S^3$ has **no zero mode** ($|\lambda_n| = n+3/2 \ge 3/2$). Every level contributes to any Dirac trace, and the ground-level contribution is nonzero: $g_0^{\text{Weyl}} = 2$, $|\lambda_0| = 3/2$.

This is the direct spectral-theoretic reason no universal Dirac identity reproduces $B$: the Dirac spectrum is gapped, the scalar spectrum is not, and $B(m)$'s polynomial signature encodes that gap.

The single hit at $m=3$ — where $(m-1)(m+2) = 2\cdot 5 = 10$ — exploits the Paper 0 packing integers $(N_{\text{init}}, l_{\max}+2) = (2, 4)$… **not quite**. $N_{\text{init}}=2$ and $l_{\max}+2 = 4$ multiply to 8, not 10. The numerical coincidence $(m-1)(m+2) = 10$ at $m=3$ has no clean Paper 0 interpretation.

---

## 4. Verdict for D5

**B-doesn't-lift (nuanced).**

- **Does not lift**: no closed-form Dirac truncation reproduces $B(m)$ for general $m$. Every apparent hit is a single-point rational coincidence whose ratio-as-a-function-of-$m$ is non-constant.
- **Nuance worth flagging to D5**: at Paper 2's specific cutoff $m=3$ and ONLY at that cutoff, the top-of-window Weyl weighted degeneracy hits $B$ exactly: $|\lambda_2|\cdot g_2^{\text{Weyl}} = 42$. Similarly $\sum |\lambda|^2 g^{\text{Dirac}} = 9B$ at $m=3$ only. These are artifacts of the elementary identity $(m-1)(m+2) = 10$ at $m=3$, not spectral theorems.
- **Already known**: the Dirac sector reproduces $\Delta^{-1} = 40$ via $g_3^{\text{Dirac}}$ (Phase 4H SM-D, unchanged). That identity is structural (the CH formula evaluated at a single level); the B-style truncated-trace identities are not.

Combined with D3 and D4 outputs, this contributes to the "three-tier coincidence" end of the D5 verdict matrix. If D3 also fails to lift $F = \pi^2/6$ to the Dirac Dirichlet series, and D4's Hopf decomposition gives no clean Paper-2 link, the sprint-level verdict is "K = π(B + F − Δ) remains a three-tier coincidence, now formally documented."

---

## 5. Structural observation (why this matters even as a negative)

Phase 4H SM-D identified $\Delta^{-1} = 40$ as a **single-level** Dirac degeneracy at $n_{\text{CH}}=3$. Paper 2's $B = 42$ is a **truncated cumulative trace** over scalar shells with Casimir weighting. These are structurally different object types:

| Ingredient | Origin | Weighting |
|---|---|---|
| $B = 42$ | cumulative truncated trace, scalar $(-\Delta_{\text{LB}})$ spectrum on $S^3$ | Casimir $l(l+1)$ |
| $F = \pi^2/6$ | infinite Dirichlet series, scalar degeneracies at $s = d_{\max}=4$ | $n^{-4}$ |
| $\Delta^{-1} = 40$ | single-level spin-Dirac degeneracy at $n=3$ | unweighted |

After D2, $\Delta$ still sits cleanly in the Dirac sector as a single-level invariant, but $B$ stays in the scalar sector as a cumulative invariant. The "common-generator" hope — that Dirac-on-$S^3$ produces all three of $B, F, \Delta$ from one spectral construction — is not supported by the Casimir-trace branch tested here. The remaining avenue (D4, Hopf-equivariant) is the only route left that could restore a unified structure.

---

## 6. Deliverables

- `debug/dirac_d2_casimir_trace.py` — pure-sympy script, candidates (a)–(e).
- `debug/data/dirac_d2_casimir_trace.json` — full structured output, including closed forms, ratios, and the exact solve $(m-1)(m+2)=10 \Rightarrow m=3$.
- `debug/dirac_d2_memo.md` — this memo.

No modifications to `geovac/dirac_s3.py`; its API (dirac_eigenvalue_abs, dirac_degeneracy) was sufficient. No modifications to Paper 2 (deferred to D5 per sprint plan).
