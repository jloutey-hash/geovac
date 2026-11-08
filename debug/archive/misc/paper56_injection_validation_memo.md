# Paper 56 `thm:injection_g4` — genuine period-map validation

**Date:** 2026-06-16
**Scope:** Replace the tautological/hardcoded backing of Paper 56's headline
theorem `thm:injection_g4` (cosmic-Galois injection
`U*_GV ↪ G_4`) with a genuine period-map test, or report exactly what
genuine backing requires.
**Deliverable test:** `tests/test_paper56_injection_g4_periodmap.py`
**Paper edits:** none (per task; this is a validation/finding memo).

---

## VERDICT

**GENUINELY-REFUTED (C4, the substantive content) + NEEDS-IMPLEMENTATION
(a richer period map) for any positive backing.**

The C4 closed-immersion / depth-1 M3-column injectivity — which Paper 56's
own honest-scope remark (`rem:injection_honest_scope`) names as *"the only
place where the injection has genuinely new structural content"* — is
**refuted** under the only concrete per-sector M3 period evaluation the
corpus actually supplies. Computing the Gram matrix of the real M3 period
vectors (instead of hardcoding `eye(n)`) gives **rank 1, det 0** at both
n_max = 2 and n_max = 3. Distinct M3 generators map to *collinear* period
values, not Q-linearly independent ones.

C1 and C3 are genuinely backable and pass; C2 is definitional (not
falsifiable in any non-trivial way). C4 — the load-bearing claim — fails
its own honest test.

---

## What the paper actually specifies (the construction)

1. **Period map** (`def:period_map`, paper line ~1174): sends the M3
   primitive generator `x_{(n,l),1}` to `m_3^{(n,l)} ∈ M_3 ⊂ MT(Z[i,1/2])`
   via the master Mellin engine `M[Tr(D e^{-tD²})]` on the
   Camporesi–Higuchi S³ spectrum.

2. **Master Mellin engine** (Paper 55 `eq:mellin_extract`, `eq:ch_spectrum`):
   the k=1 (M3) output is the **global** Dirac spectral zeta
   `D(s) = 2 ζ(s−2, 3/2) − ½ ζ(s, 3/2)` — a single function summed over
   the *whole* spectrum (all shells n). **It is not sector-resolved on its
   face.** Every M3 output on CH-S³ is a Q-linear combination of the
   *shared* Hurwitz quarter-integer values {ζ(s,1/4), ζ(s,3/4), ζ(s,5/4)}
   (Paper 55 §sec:m3, `thm:m3_cyclotomic_mixed_tate`).

3. **The only per-sector M3 quantity in the corpus** (Paper 55 lines
   ~2205–2236): the chirality-symmetrized per-sector eta value
   `η_{(n,l)} = (2l+1)(2n+1)` for l<n, `n(2n+1)` for l=n — explicitly
   identified there as *"the leading t⁰ η-density supertrace coefficient =
   M3"* per sector. It is a **single integer per sector**, n_max-independent.

So under the construction the paper gives, the M3 generator of sector
(n,l) maps to `η_{(n,l)} · b` where `b` is the *shared* Hurwitz building
block. The period-coefficient **vector** of each generator is therefore
`η_{(n,l)}` times a fixed vector — i.e. all generators are collinear.

## What the current (tautological) test does

`tests/test_paper56_injection_g4.py`:
- **C1**: `simplify(s_a*s_b − s_a*s_b) == 0` — `lhs` and the compared
  `product` are the *same expression*; the residual is identically 0 by
  construction. No period map, no evaluation.
- **C4**: `gram = eye(n); assert gram.det() == 1` — the Gram is **hardcoded
  to the identity**, which *is* the claim under test (non-degeneracy)
  asserted as a premise. Plus `test_goncharov_deligne_faithfulness_consequence`
  is literally `assert True`.
- The remaining tests are count-bookkeeping (`len(gens)`, closed-form panel
  totals) — fine as bookkeeping, but they verify arithmetic of the residual
  *counts*, not the *theorem*.

## What I built and the result

`tests/test_paper56_injection_g4_periodmap.py` constructs the period map
from the paper's stated construction and computes — does not assume — the
Gram matrix:

```
n_max=2: η_(n,l) = [3, 3, 5, 15, 10]   (sectors (1,0)(1,1)(2,0)(2,1)(2,2))
         distinct values 4/5  →  (1,0) and (1,1) COLLIDE at 3
         GENUINE Gram rank = 1, det = 0   →  C4 injectivity REFUTED
n_max=3: η_(n,l) = [3, 3, 5, 15, 10, 7, 21, 35, 21]
         distinct values 7/9  →  (1,0)≡(1,1)=3 and (3,1)≡(3,3)=21
         GENUINE Gram rank = 1, det = 0   →  C4 injectivity REFUTED
```

Test run: **14 passed, 2 xfailed** (the 2 strict-xfail are the honest
`rank == N` C4 assertion, recorded as a *live falsifier*: if a richer
period map ever makes injectivity hold, the strict xfail turns into an
XPASS failure and forces a re-read). The passing tests positively assert
the refutation: genuine Gram is rank 1 (not `eye(N)`), det 0, and there
exist explicit period-value collisions.

### Per-compatibility assessment

| Comp. | Current test | Genuine status |
|:--|:--|:--|
| **C1** multiplicativity | tautology `x*y − x*y` | **GENUINELY-BACKED** — evaluate the M3 generators on their real scalar period content, check `π(x·y)=π(x)·π(y)` (passes). Honest but weak: multiplicativity of a Q-algebra map is essentially the universal property of `Sym_Q(V)`. |
| **C2** depth-1 coproduct | `assert coproduct_form is not None` | **DEFINITIONAL** — `Δx = x⊗1 + 1⊗x` holds by *definition* of the primitive Hopf substrate; there is nothing falsifiable to compute. Not tautological in a misleading way, but not load-bearing either. |
| **C3** SL₂→Levi via χ₄ | `det(g)=1` on 5-element panel | **GENUINELY-BACKED** — det of explicit non-diagonal rational matrices; real computation, passes. (This one was already fine.) |
| **C4** closed-immersion / M3 injectivity | `gram = eye(n)` hardcoded | **GENUINELY-REFUTED** — real Gram is rank 1, det 0. This is the substantive content of the theorem per the paper's own remark. |

## Why C4 fails / what a positive backing would require (NEEDS-IMPLEMENTATION scope)

The refutation is structural, not a precision artifact: the per-sector M3
content is a *single scalar* (one Hurwitz building block × η), so the
period vectors are 1-dimensional and collinear. The paper invokes
Goncharov–Deligne 2005 faithfulness to assert Q-linear independence of
distinct generators' periods, but **the construction never produces a
period vector of dimension > 1 per sector**, so faithfulness has nothing
to act on at the sector-resolution level. The hardcoded `eye(n)` papered
over exactly this.

A genuine *positive* backing of C4 would require, concretely, one of:

1. **A sector-resolved M3 period map of rank ≥ 2.** Define
   `π(x_{(n,l),1})` as an honest vector in the 3-dim quarter-integer
   Hurwitz space {ζ(s,1/4), ζ(s,3/4), ζ(s,5/4)} (or a depth/weight-graded
   Glanois `B⁴` basis), derived from a **per-sector** parity-decomposition
   of `D(s)` — i.e. split `D_even/D_odd` *within* each shell (n,l), not on
   global Fock parity. The corpus does **not** currently contain such a
   per-sector decomposition; the parity split is on the global Fock quantum
   number n (Paper 55 §sec:m3), giving one vector for the whole spectrum,
   not N(n_max) independent ones. This is the missing object. Implementing
   it means deriving, from the CH spinor decomposition, the per-(n,l)
   coefficients of the three Hurwitz blocks — a real derivation sprint, not
   a test-writing task, and it may well still collapse (the honest-scope
   remark's "M1/M2 collapse" reasoning applies: the period side is
   genuinely low-dimensional per weight).

2. **A weight/multi-s evaluation.** Evaluate each generator's period at a
   *grid* of integer arguments `s ∈ {2,3,4,…}` so each generator becomes a
   vector `(η_{(n,l)}·D-block(s))_s`. But because the per-sector content is
   `η_{(n,l)} × (shared s-dependent block)`, the resulting matrix is still
   rank 1 (outer product of the η-vector with the block-vs-s vector). Same
   collapse. This does **not** rescue C4.

3. **Honest restatement of the theorem.** If Reading A (abelianization) is
   the intended geometry, the truthful claim is *not* "distinct M3
   generators → independent periods" but rather "the M3 column factors
   through a 1-dimensional (rank-1) period quotient per weight," i.e. the
   injection is of `U*_GV^{ab}` onto a **rank-1** image, with the
   sector-multiplicity collapsing exactly as M1/M2 do. Under that reading
   C4's "closed immersion" survives only as a statement about the
   abelianized/collapsed image, and the current `eye(n)` panel is simply
   the wrong object to verify it.

**Recommendation (not applied — papers not edited per task):** the C4
panel in Paper 56 §sec:injection_panel and the `eye(n)` claim in
`test_paper56_injection_g4.py` overstate what the construction supports.
Either (a) supply the rank-≥2 per-sector period map of route 1 (derivation
sprint, uncertain outcome), or (b) restate C4 honestly per route 3
(rank-1/abelianized image). This should go to the PI — it touches the
status of a paper's headline theorem (§9 disposition: "a load-bearing
claim with no/weak/false-positive backing" → raise to PI).

## Files

- **Created:** `tests/test_paper56_injection_g4_periodmap.py` — genuine
  period-map test (14 passed, 2 strict-xfail = live C4 falsifier).
- **Read (no edits):** `papers/group3_foundations/paper_56_tannakian_substrate.tex`
  (§sec:injection_g4, `def:period_map`, `thm:injection_g4`,
  `sec:injection_panel`); `papers/group3_foundations/paper_55_periods_of_geovac.tex`
  (§sec:m3, per-sector η); `geovac/pro_system.py`; `geovac/tannakian.py`;
  `tests/test_paper56_injection_g4.py`; `debug/sprint_injection_g4_memo.md`.
