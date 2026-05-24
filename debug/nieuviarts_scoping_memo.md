# Nieuviarts Twist-Morphism Scoping for GeoVac S³ = SU(2) Camporesi–Higuchi

**Date:** 2026-05-16
**Sprint window:** 1–2 day scoping pass
**Trigger:** L1 literature update flagged trigger (b) FIRED — Nieuviarts twisted-spectral-triple stabilization (arXiv:2502.18105v3, May 2025; proceedings 2512.15450, May 2026). The May 2026 Track B scoping memo had explicitly flagged this exact contingency ("If Nieuviarts' 2025 papers become the standard prescription quickly, the route changes") and committed to re-checking. This is that re-check.

## Verdict

**NO-GO. The Nieuviarts twist morphism does NOT apply to the GeoVac S³ = SU(2) Camporesi–Higuchi spectral triple. The obstruction is structural: GeoVac sits at KO-dim 3 (mod 8), which is odd, and the Nieuviarts construction explicitly requires even-dimensional manifolds at all three published stages of the program. The BBB direct lift via (m,n) = (3,0) → (3,1) Krein-space embedding remains the only published path, with the 3–6 month sub-budget intact and the 6–12 month Lorentzian-propinquity blocker unchanged.**

## Evidence

### arXiv:2502.18105v3 (Nieuviarts, May 2025) — verified real

**Abstract (verbatim, extracted from arXiv abstract page):**

> "This article demonstrates how the transition from a (Riemannian) twisted spectral triple to a pseudo-Riemannian spectral triple arises within an almost-commutative spectral triple. This opens a new perspective on the Lorentzian signature problem, showing that the almost-commutative structure at the heart of the noncommutative standard model of particle physics could be the origin of the emergence of a Lorentzian spectral triple, starting from self-adjoint Dirac operators and conventional inner product structures, in the framework of twisted spectral triples. We present an alternative to Wick rotation, acting on the metric and the Christoffel symbols in a way that does not introduce any complex numbers."

**Definitive odd-dim restriction (extracted from arXiv HTML, Section 2.1, Definition 2.2):**

> "The fact that no such procedure for twisting odd-dimensional manifold's spectral triple have been found will be one of the justifications to focus on the study of even-dimensional manifolds."

The restriction is structural, not editorial. The twist-by-grading mechanism (Def. 2.2) requires decomposing the Hilbert space via grading projections $p_\pm := (1 \pm \Gamma)/2$, which only makes sense when $\Gamma$ exists as a $\mathbb{Z}_2$ grading commuting with the algebra and anti-commuting with the Dirac — i.e. on **even-dim KO-cycles** in the Connes sense. The "even-dimensional" requirement appears 12+ times in the paper as structural, not optional. No worked examples on specific manifolds are given (Section 5.3 references "the 4-dimensional case" without computing it).

### arXiv:2402.05839 (Nieuviarts, 2024) — verified real, foundational predecessor

**Abstract (verbatim):**

> "We present a connection between twisted spectral triples and pseudo-Riemannian spectral triples, rooted in the fundamental interplay between twists and Krein products. A concept of morphism of spectral triples is introduced, transforming one spectral triple into its dual. **In the case of even-dimensional manifolds**, we demonstrate how this construction implements a local signature change via the parity operator induced by the twist..." (emphasis added)

Even-dim restriction is **explicit in the abstract**. Submission history shows revisions through March 2026 (v6); the restriction was not lifted in any revision (it is inherited by 2502.18105v3 explicitly).

### arXiv:2512.15450v2 (Nieuviarts, May 2026 proceedings) — verified real

**Abstract (verbatim):**

> "This proceeding presents a synthesis of recent results on the emergence of pseudo-Riemannian structures from twisted spectral triples within the almost-commutative framework. It provides a unified algebraic mechanism for addressing the Lorentzian signature problem, demonstrating how the almost-commutative structure underlying the noncommutative Standard Model of particle physics may give rise to Lorentzian spectral triple from a purely Riemannian setting. This notably offers an alternative to Wick rotation, provided by a notion of morphism connecting twisted and pseudo-Riemannian spectral triples."

Proceedings synthesize 2402.05839 + 2502.18105 — no new manifolds, no extension to odd-dim. S³, SU(2), and spheres of any dimension are not mentioned.

### Phantom check

All three arXiv IDs verified real. No new phantoms. (Phantom-check rule honored from L1 update — "Hawkins–Skoda" remained the only caught phantom.)

## Applicability test against GeoVac structure

| GeoVac feature | Nieuviarts requirement | Verdict |
|:---|:---|:---|
| $\mathcal{A}_{GV} = \mathbb{C}^{V_\text{Fock}}$ commutative | Any involutive algebra OK | PASS |
| $\mathcal{H}_{GV}$ = Camporesi–Higuchi spinor bundle on $S^3$ | Hilbert space with $\mathbb{Z}_2$ grading $\Gamma$ from KO-dim even | **FAIL** |
| $D_{GV}$ diagonal CH, eigenvalues $\pm(n+3/2)$ | Self-adjoint Dirac with compact resolvent | PASS |
| $J$ quaternionic, $J^2 = -I$, $JD = +DJ$, KO-dim 3 | KO-dim must be **even** (Table 1 of 2502.18105v3 lists cases $\{0,2,4,6\} \bmod 8$) | **FAIL** |
| Real structure verified at $n_\max \in \{1,2,3\}$ on truthful CH | Twisted real structure $J := K\hat{J}$ requires twist $K$ and grading $\Gamma$ | **FAIL** (no $\Gamma$) |
| $S^3 = SU(2)$ is 3-dim (odd) | Even-dim only (explicit, Def. 2.2 + abstract of 2402.05839) | **FAIL** |
| Where Lorentz emerges in Nieuviarts | **Inner / AC factor**, not outer geometric factor | Even if construction applied, GeoVac's outer factor is the part that needs Lorentz signature |

The last row is the second structural problem even if dim parity is set aside: Nieuviarts's Lorentz emergence happens **on the almost-commutative inner factor** (the abstract is explicit — "the almost-commutative structure at the heart of the noncommutative standard model... could be the origin of the emergence"). GeoVac needs Lorentz signature on the **outer** $S^3$ factor (the geometric Dirac side), because the physics question is about modular flow on the round-$S^3$ Bisognano–Wichmann reading, not about an electroweak signature in the inner factor. The construction is positioned in the wrong layer of the AC tensor product for GeoVac's use case.

## Why this NO-GO is reliable

1. **The odd-dim restriction is structural, not editorial.** Twist-by-grading mechanism requires $\Gamma$; $\Gamma$ exists for even KO-dim only; this is Connes' axiomatic framework, not Nieuviarts's preference. Lifting the restriction would require new mathematics, which the v3 (May 2025) paper explicitly says has not been found.

2. **The restriction is consistent across all three Nieuviarts papers.** 2402.05839 (foundational, with v6 March 2026), 2502.18105v3 (May 2025, the construction paper), 2512.15450 (May 2026 proceedings synthesis) — none extends to odd-dim.

3. **The Layer in which Lorentz emerges is wrong for GeoVac.** Even if the odd-dim block were removed, GeoVac would need outer-factor Lorentz signature on $S^3 \to S^3 \times \mathbb{R}$-like extension, while Nieuviarts injects the signature on the inner AC factor sitting alongside a Riemannian outer factor.

4. **The May 2026 Track B scoping memo predicted this contingency.** §9 of `debug/lorentz_boost_scoping_memo.md` named "If Nieuviarts' 2025 papers become standard prescription quickly, the route changes" as an assumption to monitor; the re-check landed on "no, the restriction is structural." Track B's NO-GO on the multi-month extension remains the standing verdict.

## Cost impact

**Zero shortening.** The BBB direct lift via (3,0) → (3,1) Krein-space embedding (Bizi–Brouder–Besnard 2018, arXiv:1611.07062) remains the only published path. The May 2026 Track B cost estimate stands:

| Item | Estimate (unchanged) |
|:---|:---|
| Krein lift of $\mathcal{H}_{GV}$ to $(m,n)=(3,1)$ | 2–4 weeks (mechanical) |
| Lorentzian Connes axioms on Krein side (BBB Table 2: $J^2 = +I$, $JD = -DJ$, $\gamma \to \gamma^5$) | 2–4 weeks (axiom-side) |
| **Lorentzian propinquity (original NCG-math)** | **6–12 months, HIGH stall risk** |
| Lorentzian L1'/L2/L3/L4/L5 analogs of Paper 38 lemmas | 2–4 months (gated on propinquity) |
| **Bottom-up total** | **9–18 months, dominated by propinquity** |

The 6–12 month Lorentzian-propinquity blocker is the load-bearing cost. Nieuviarts does not address propinquity; it addresses signature change. Even if (hypothetically) the odd-dim restriction were lifted, GH-convergence in the Lorentzian setting would still be unaddressed in the literature. So the cost-shortening hope was bounded above by the 2–4 weeks of axiom-side mechanical work even in the best case, not by the dominant 6–12 month propinquity cost.

## Recommendation: no sprint open

**Confirm BBB direct lift as the named path. Defer until propinquity-side load-bearing requirement materializes.** The recommendation from May 2026 Track B is unchanged:

- **Status of multi-month Lorentzian extension: NO-GO**, 9–18 months, opportunity cost dominated by multi-focal precision arc + second-row chemistry + Cs HFS heavy-atom direction.
- **Status of sprint-scale Lorentzian-adjacent work: GO** for Candidate A (Wick-rotation map + SO(4,1) conformal-action memo + Paper 34 §VIII open-question paragraph) per Track B §7. Already partly executed by Sprint Unruh+BW (May 9–10).
- **Re-trigger condition:** A specific Lorentzian application becomes load-bearing for a GeoVac observable (Track B trigger (c)); or a Lorentzian propinquity is published (Track B trigger (a)); or the odd-dim restriction is lifted in published work (Track B trigger (b) — re-fires only if explicitly lifted, not just because new Nieuviarts work appears).

The May 9–10 Unruh + BW sprints landed the **structural-correspondence** verdict on the M1 Hopf-base / Vol($S^1$) = $2\pi$ unification (Paper 32 §VIII Remark `rem:bisognano_wichmann_reading`; Paper 35 §VIII). That structural reading does NOT require a Lorentzian spectral triple; it operates at the metric-functional-value level (M1 mechanism). Promotion to literal identification at the operator-system level would require the multi-month extension, and the cost-benefit remains unfavorable.

## Honest scope of this scoping pass

- **What I read fully:** abstracts of all three Nieuviarts papers (verbatim quoted above); key Definition 2.2 + Section 2.1 text of 2502.18105v3 via arXiv HTML extraction (sufficient to nail the odd-dim restriction definitively); existing Paper 32 §VIII Lorentzian references; existing `debug/lorentz_boost_scoping_memo.md` for cost-estimate cross-check.
- **What I could NOT read:** the full PDF of 2502.18105v3 (PDF text extraction failed — binary handling issue, common with arXiv). The HTML extraction of Section 2.1 + Definition 2.2 + Table 1 + Section 5.3 references is sufficient for the verdict because the odd-dim restriction is **explicit in a definition** and **inherited from a published abstract** (2402.05839). No deeper reading would reverse it.
- **What I am NOT claiming:** I have not verified every theorem of 2502.18105v3 or assessed whether the twist construction has bugs. The verdict is purely on applicability: even if the construction is sound, it does not apply to KO-dim 3.
- **Where math.OA specialist would add value:** confirming that the obstruction is structural rather than removable (it is — the $\Gamma$ grading is part of the Connes axiomatic framework, not Nieuviarts's choice — but a specialist would back this with the cleanest reference). Not needed for the NO-GO at sprint resolution.

## Cross-references

- `debug/lorentzian_literature_update_2026_05_16.md` — L1 update that flagged trigger (b)
- `debug/lorentz_boost_scoping_memo.md` (May 2026 Track B) — original NO-GO with cost estimates this memo confirms unchanged
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` §VIII — case-exhaustion theorem and GH-convergence (Riemannian) where Lorentzian extension is named as open question, not load-bearing
- Bizi–Brouder–Besnard 2018, arXiv:1611.07062 — verified real, the path that remains
- Nieuviarts arXiv:2502.18105v3, 2402.05839, 2512.15450v2 — all verified real, all carrying the even-dim restriction
