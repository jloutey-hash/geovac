# Sprint NA-1 — Hidden non-abelian structure probe

**Date:** 2026-06-06
**Sprint:** NA-1 of Q5'-Tannakian-Closure (after TC-2 series + cosmic-Galois comparison sub-agent reports)
**Driver:** `debug/compute_q5p_na1_non_abelian_probe.py`
**Module:** `geovac/tannakian.py`, `geovac/pro_system.py` (re-used substrate)
**Data:** `debug/data/sprint_q5p_na1_non_abelian_probe.json`
**Wall time:** 0.25 s
**Discipline:** bit-exact `sympy.Rational` / `sympy.Integer`.

---

## 1. TL;DR

**Verdict: STRUCTURAL — abelian-by-construction.**

The four sub-agent reports identified a sharp structural gap between GeoVac's $U^*_{GV} = \mathbb{G}_a^\infty \rtimes SL_2$ (abelian unipotent radical) and the classical motivic Galois groups (Brown's $\mathcal{G}_{MT(\mathbb{Z})}$, Connes-Marcolli's $U^*_{CM}$ — both have **free non-abelian** unipotent radicals). NA-1 probes whether GeoVac has hidden non-abelian content not yet detected, or whether the abelianness is genuinely structural.

**Answer:** the abelianness is *structurally forced* by the choice of primitive Hopf algebra $\mathcal{H}_{GV} = \mathrm{Sym}_\mathbb{Q}(V_{n_{\max}})$. By the Cartier-Milnor-Moore theorem (characteristic 0), any primitive Hopf algebra has abelian Tannakian dual. There is no hidden non-abelian content at the substrate level.

**For non-abelian content, the substrate would need to be enriched** to the cofree-cocommutative shuffle Hopf algebra $T(V)$ with deconcatenation coproduct — this is the motivic-style structure whose Tannakian dual is the free non-abelian unipotent group. Whether GeoVac's *physical content* (Mellin moments, spectral data, master Mellin engine) naturally respects the primitive coproduct (abelian dual) or the shuffle coproduct (non-abelian dual) is **a named multi-year open question** with one sprint-scale concrete test identified.

---

## 2. Two possible structural readings (both honest)

### Reading A — GeoVac is the abelianization of cosmic Galois.

If GeoVac's physical content respects the primitive coproduct (which is what the current substrate uses), then $U^*_{GV}$ is precisely the **abelianization** of Connes-Marcolli's $U^*_{CM}$:\ the quotient by the commutator ideal of the free Lie algebra. Under this reading, GeoVac is a *weight-only* geometric realization that drops the depth/iteration content of motivic Galois. The "geometric seed" picture sharpens to:\ GeoVac sees the *bare grading structure* of cosmic Galois (one generator per Mellin slot, one $\mathbb{G}_a$ per generator) plus the $SL_2$ enhancement from $S^3$ geometry, but does NOT see the bracket structure that makes ζ(3), ζ(5), ζ(7), ... algebraically independent in the motivic sense.

### Reading B — GeoVac's substrate needs shuffle enrichment.

If GeoVac's physical content *actually* respects the shuffle coproduct (which is what one would discover by computing depth-2 Mellin moments and finding they factorise asymmetrically rather than symmetrically), then the primitive Hopf algebra we built was a *first-order approximation* and the correct substrate is the shuffle Hopf algebra $T(V)$. Under this reading, GeoVac IS a candidate geometric realization of (a sub-quotient of) the full motivic Galois group, with $U^*_{GV}$ upgrading to a free non-abelian unipotent group whose Lie algebra is free on the Mellin generators.

**Both readings are consistent with the current evidence.** Distinguishing them requires the sprint-scale test in §6.

---

## 3. Bit-exact panel

### Part A — Commutativity on witness panel ($n_{\max} = 2$, 15 generators, 15² = 225 pair-on-rep tests)

For every pair $(g_1, g_2)$ of primitive generators and every witness rep $V_{g_{\mathrm{ref}}}$, the commutator $[X_{g_1}^{V_{g_{\mathrm{ref}}}}, X_{g_2}^{V_{g_{\mathrm{ref}}}}]$ vanishes bit-exactly. This is trivial by construction (each $V_g$ has only $X_g$ non-zero), but the verification confirms no implementation slip.

**3{,}375 / 3{,}375 bit-exact zero residuals** ($15 \times 15$ pairs $\times$ 15 ref reps).

### Part B — Forcing argument: commuting nilpotents on dim 2 are proportional

Symbolic computation via `sympy.solve`:\ given $X_1 = E_{12}$ canonical, find all $X_2 \in M_2(\mathbb{Q})$ satisfying $X_2^2 = 0$ (nilpotent) and $[X_1, X_2] = 0$ (commuting with $X_1$). The solution variety is

$$X_2 = b_p \cdot E_{12} = b_p \cdot X_1, \qquad b_p \in \mathbb{Q} \text{ free}.$$

One-dimensional, proportional to $X_1$. **Conclusion:\ even if we tried to construct a single 2-dim rep where two primitive generators both act non-trivially, the abelian Hopf-algebra structure (which requires all generator actions to commute) would force them to act *proportionally*.** There is no non-abelian arrangement on a 2-dim rep within the primitive Hopf-algebra category.

### Part C — Primitive vs shuffle coproduct on depth-2 element

For a depth-2 element $x_1 \cdot x_2$ (viewing $x_1, x_2$ as generators of $V$):

**Primitive coproduct** on $\mathrm{Sym}(V)$:
$$\Delta_{\mathrm{prim}}(x_1 \cdot x_2) = x_1 x_2 \otimes 1 + x_1 \otimes x_2 + x_2 \otimes x_1 + 1 \otimes x_1 x_2.$$

Four terms, **symmetric** in $x_1 \leftrightarrow x_2$ (the cross-terms $x_1 \otimes x_2$ and $x_2 \otimes x_1$ are both present).

**Deconcatenation coproduct** on $T(V)$:
$$\Delta_{\mathrm{dec}}(x_1 \otimes x_2) = (x_1 \otimes x_2) \otimes 1 + x_1 \otimes x_2 + 1 \otimes (x_1 \otimes x_2).$$

Three terms, **asymmetric** in ordering (only the $x_1 \otimes x_2$ partition appears, not the $x_2 \otimes x_1$).

The structural difference is what distinguishes abelian-dual primitive from non-abelian-dual shuffle. At depth 1 (single generator) both coproducts agree:\ $\Delta(x_g) = x_g \otimes 1 + 1 \otimes x_g$. The non-abelian content of motivic Hopf algebras lives *exclusively* at depth $\geq 2$ via the ordering asymmetry.

### Part D — Multi-year open question

Sprint-scale concrete test:\ take two specific Mellin moments at depth 1 (an M2 Seeley-DeWitt coefficient and an M3 vertex-parity Hurwitz coefficient). Compute their *joint* Mellin transform at depth 2. Two possible outcomes:

- **Symmetric** factorisation (M2$_a$ · M3$_b$ = M3$_b$ · M2$_a$, primitive product):\ GeoVac respects primitive coproduct, Reading A is correct, $U^*_{GV}$ is the abelianization.
- **Asymmetric** factorisation (deconcatenation pair, with ordering content):\ GeoVac respects shuffle coproduct, Reading B is correct, substrate needs enrichment.

The test is sprint-scale (1-2 weeks) and concrete. It is the first natural follow-on flagged by this sprint.

---

## 4. Honest scope

**Closed by NA-1:**
- The current substrate $\mathcal{H}_{GV} = \mathrm{Sym}(V)$ with primitive coproduct gives abelian Tannakian dual by construction (Cartier-Milnor-Moore).
- No hidden non-abelian content exists at the witness-panel level (bit-exact verification).
- The natural enrichment direction (cofree-cocommutative shuffle Hopf on $T(V)$) is identified explicitly with reference to the motivic standard.

**Open (named, sprint-scale):**
- Does GeoVac's depth-2 physical content respect primitive or shuffle? Concrete test in §3 Part D.

**Open (multi-year):**
- If shuffle wins:\ build the enriched Hopf algebra structure on GeoVac's substrate and verify it transports the per-cutoff Tannakian closure (TC-2a/b/c/d) to a non-abelian closure with free Lie unipotent.
- If primitive wins:\ formalise the abelianization claim — $U^*_{GV} = $ abelianization of $U^*_{CM}$ — as a clean structural statement.

---

## 5. Discipline checks

**Curve-fit audit (`feedback_audit_numerical_claims`):** Zero free parameters. The claim "primitive Hopf → abelian dual" is a classical theorem (Cartier-Milnor-Moore, characteristic 0). The bit-exact Part B verification has zero free parameters (single solve of a polynomial system).

**Discrete-for-skeleton:** `sympy.Rational` throughout.

**Tag transcendentals:** None appear.

**Diagnostic-before-engineering:** This IS a diagnostic-before-engineering sprint — we recognised after the cosmic-Galois sub-agent reports that engaging with the abelian-vs-non-abelian question structurally is the right pre-requisite for any future probe.

**WH1 PROVEN unaffected.** **Hard prohibitions check clean.**

---

## 6. Sprint-scale follow-ons (named, not done here)

1. **Depth-2 Mellin moment test** (§3 Part D). 1-2 sprints. Probes Reading A vs Reading B.
2. **Hodge-theoretic SL_2 probe.** Identify whether GeoVac's $SL_2$ corresponds to the Mumford-Tate $SL_2$ of a polarized variation of mixed Hodge structures on the $S^3$ substrate. Track A sub-agent flagged this as the natural alternative if our $SL_2$ is not simply the geometric $\mathrm{Spin}(3)$. 2-4 sprints.
3. **Explicit injection $U^*_{GV} \hookrightarrow \mathcal{G}_4$ at the period-ring level.** From Track C:\ provable today using Eskandari et al. 2025 + Deligne 2010. 1-2 sprints.
4. **Hunt for depth-2 level-4 cyclotomic MZV in a GeoVac observable.** If we find one, GeoVac's exhaustion of $\mathcal{G}_4$-periods gains evidence. 2-3 sprints exploratory.

---

## 7. Files

### Produced
- `debug/compute_q5p_na1_non_abelian_probe.py` — driver (~330 lines, 0.25 s).
- `debug/data/sprint_q5p_na1_non_abelian_probe.json` — data dump.
- `debug/sprint_q5p_na1_non_abelian_probe_memo.md` — this memo.

### References
- Cartier, P. ``A primer of Hopf algebras'' in Frontiers in Number Theory, Physics, and Geometry II (Springer, 2007).
- Brown, F. ``Mixed Tate motives over $\mathbb{Z}$'' Ann. Math. 175 (2012) 949–976.
- Connes, A.; Marcolli, M. ``Noncommutative Geometry, Quantum Fields and Motives'' AMS (2008), chapter 17.

---

## 8. One-line verdict

**STRUCTURAL — abelian-by-construction.** $\mathcal{H}_{GV} = \mathrm{Sym}(V)$ with primitive coproduct has abelian Tannakian dual by Cartier-Milnor-Moore;\ no hidden non-abelian content;\ the natural enrichment to cofree-cocommutative shuffle Hopf algebra $T(V)$ is identified;\ whether GeoVac's depth-2 physical content respects primitive or shuffle is a concrete sprint-scale follow-on. Two structural readings of the cosmic-Galois comparison remain honest:\ Reading A (GeoVac is the abelianization of $U^*_{CM}$) and Reading B (substrate needs shuffle enrichment).
