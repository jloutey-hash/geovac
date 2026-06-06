# Sprint Q5'-HardParts — joint closure of the two genuinely hard Stage-1/Stage-2 follow-ons

Date: 2026-06-05 (third Q5' arc of the day, after v3.59.0 morning [Stage1-Followon] and v3.60.0 afternoon [Stage1-Prosystem])
Scope: two parallel sub-sprints (Track A: Stage 2 Hopf algebra of Mellin moments; Track B: strict-strong-form cochain-morphism functoriality) dispatched as the genuinely "hard parts" of the cosmic-Galois $U^*$ program. Track A returned first; Track B was re-dispatched after a token-quota interruption and returned with a sharper-than-expected verdict.

Sub-sprint memos (detail):
- `debug/sprint_q5p_stage2_hopf_memo.md` (Track A: Stage 2 candidate Hopf algebra $\mathcal{H}_{\mathrm{GV}}$)
- `debug/sprint_q5p_strict_strong_memo.md` (Track B: strict-strong-form cochain-morphism functoriality)

## 1. TL;DR

Two parallel sub-sprints closed the two named "hard" follow-ons of the cosmic-Galois $U^*$ program. Track A returned the first scoping step of Stage 2 with a concrete bit-exact finite-cutoff candidate Hopf algebra $\mathcal{H}_{\mathrm{GV}}^{(n_{\max})}$ + its candidate motivic Galois group $U^{*(n_{\max})}_{\mathrm{GeoVac}}$; verdict POSITIVE. Track B sharpened the v3.60.0 honest-scope caveat into a bit-exact dichotomy: pro-system functoriality is strict at the class level (closed v3.60.0) and bit-exactly non-strict with $n_{\max} = 3$ fixed-point structural drift at the cochain-morphism level (closed today); verdict POSITIVE-STRUCTURAL-DRIFT.

**Joint structural finding (visible only because both tracks landed):** the candidate Stage-2 Hopf algebra is *abelian primitive* (Track A) precisely because cochain-morphism non-functoriality on commutative $\mathcal{A}$ (Track B) blocks the natural non-abelian nesting that would correspond to a non-primitive coproduct. Stated as a forced reading:

> The abelian primitive-cospan shape of $\mathcal{H}_{\mathrm{GV}}$ at the v3.60.0 pro-system substrate is *forced* by the bit-exact strict/non-strict dichotomy between class-level and cochain-morphism-level functoriality. Non-abelian Stage-2 content requires a richer substrate than commutative $\mathcal{A}$ + class-level pro-system data — exactly the three multi-year enrichment ingredients Track A flagged (nested-Hopf-tower $J^*(S^3)$, cross-shell off-diagonal Dirac, JLO/CM bicomplex Cuntz extension).

Both verdicts are bit-exact `sympy.Rational` (Track A: 437 axiom + truncation cells; Track B: degree-3 closure residual panel + pull-back identity at $n_{\max} \in \{2, 3, 4\}$ for both JLO and CM-$\eta$). The Mellin slot $k \in \{0, 1, 2\}$ is preserved by $\Delta$ (Track A), and the strict/non-strict dichotomy holds for both JLO and CM-$\eta$ towers (Track B retires v3.59.0's "CM is the cochain-clean alternative" reading at depth > 1).

## 2. Joint structural finding (umbrella content)

The two tracks together give a sharper Stage-1 closure than either alone. At the cocycle-class level the v3.60.0 pro-system gives strict inverse-system data on the same $K_0$ sector decomposition $\{[e_{(n, l)}]\}$; Track B today shows that this strictness does NOT lift to cochain-morphism level on commutative $\mathcal{A}$; Track A today shows that the natural Stage-2 substrate built from the strict-class-level data is abelian primitive Hopf with $U^{*(n_{\max})}_{\mathrm{GeoVac}} = \mathbb{G}_a^{3 N(n_{\max})}$.

The forced-reading synthesis: **non-abelian pro-unipotent content in the Stage-2 motivic Galois group is structurally blocked by the cochain-morphism non-functoriality of the underlying JLO bicomplex on commutative $\mathcal{A}$**. Any future non-abelian Hopf enrichment must come from a substrate enrichment, not from a clever choice of coproduct on the current substrate — Track A's three flagged ingredients (sub-sector nesting via $J^*(S^3)$; cross-shell off-diagonal Dirac; JLO/CM Cuntz extension) are categorically the right directions because each one *breaks* one of the structural conditions that forces the primitive coproduct: respectively, the sector-disjointness of CH idempotents; the sector-locality of the Dirac structure; the commutativity of the underlying algebra.

The Mellin-slot factorisation $\mathcal{H}_{\mathrm{GV}} = \mathcal{H}_{\mathrm{GV}}^{[0]} \otimes \mathcal{H}_{\mathrm{GV}}^{[1]} \otimes \mathcal{H}_{\mathrm{GV}}^{[2]}$ from Track A is the Hopf-algebraic shadow of the case-exhaustion theorem's M1/M2/M3 partition (Paper 32 §VIII Theorem `thm:pi_source_case_exhaustion`): Stage 2's motivic Galois group acts diagonally on the three mechanisms, *enforcing* the partition structurally at the cosmic-Galois level rather than as an empirical observation. This is the substantive new content beyond the case-exhaustion theorem at the Stage-2 substrate level.

## 3. Sub-sprint summaries

### Track A — Stage 2 candidate Hopf algebra

The candidate $\mathcal{H}_{\mathrm{GV}}^{(n_{\max})} = \mathrm{Sym}_{\mathbb{Q}}(V_{n_{\max}})$ on $3 N(n_{\max})$-dimensional space of primitive generators $x_{(n, l), k}$ ($k \in \{0, 1, 2\}$ corresponding to M1, M3, M2). All five Hopf axioms (coassociativity, counit-left/right, antipode-left/right, bialgebra compatibility) bit-exact at $n_{\max} \in \{2, 3\}$ — 269 zero residuals on the axiom panel + 42 generator-level $k$-preservation checks + 126 truncation Hopf-hom residuals = 437 bit-exact zero residuals total. Mellin-slot $k$-grading preserved by $\Delta$ exactly, giving the tensor decomposition $\mathcal{H}_{\mathrm{GV}} = \mathcal{H}_{\mathrm{GV}}^{[0]} \otimes \mathcal{H}_{\mathrm{GV}}^{[1]} \otimes \mathcal{H}_{\mathrm{GV}}^{[2]}$. Pro-system truncation $P_{n+1 \to n}$ from v3.60.0 lifts bit-exactly to a Hopf-algebra homomorphism. Candidate motivic Galois group: $U^{*(n_{\max})}_{\mathrm{GeoVac}} = \mathbb{G}_a^{3 N(n_{\max})}$ (dimensions 15, 27 at $n_{\max} = 2, 3$), with the M1/M3/M2-respecting product structure $\prod_k \mathbb{G}_a^{N(n_{\max})}$.

Structural shape: $\mathcal{H}_{\mathrm{GV}}$ is the **abelian primitive-cospan** of the Connes–Kreimer Hopf algebra of Feynman graphs (Connes–Kreimer 1998, arXiv:hep-th/9808042). Same generator-and-grading shape, same pro-system truncation, but the CK sub-graph coproduct is replaced by the primitive coproduct here because CH Fock sector idempotents are pairwise orthogonal (no proper sub-sector nesting). Multi-year continuation toward full Stage-2 Tannakian construction (Connes–Marcolli 2007 arXiv:math/0409306) requires substrate enrichment via one of three flagged ingredients.

Driver: `debug/compute_q5p_stage2_hopf.py`. Data: `debug/data/sprint_q5p_stage2_hopf.json`. Memo: `debug/sprint_q5p_stage2_hopf_memo.md`.

### Track B — strict-strong-form cochain-morphism drift

Adapted Sub-Sprint 2c's JLO bicomplex driver from $n_{\max} = 2$ to $n_{\max} \in \{3, 4\}$ on a justified diagnostic sub-panel (old + new palindromic 4-tuples + degree-1 sanity at full panel). Degree-3 closure $(b\phi_2 + B\phi_4)$ residuals at $t^0$:

| Cocycle | 4-tuple | $n_{\max} = 2$ | $n_{\max} = 3$ | $n_{\max} = 4$ |
|:----:|:----|:----:|:----:|:----:|
| JLO | $(e_{(2,1)}, e_{(2,2)}, e_{(2,2)}, e_{(2,1)})$ palindrome | $+1/196608$ | $-1/65536$ | $-1/65536$ (FIXED) |
| JLO | $(e_{(2,1)}, e_{(2,1)}, e_{(2,2)}, e_{(2,2)})$ | $-1/196608$ | $+1/65536$ | $+1/65536$ (FIXED) |
| JLO | $(e_{(1,0)}, e_{(1,1)}, e_{(1,1)}, e_{(1,0)})$ palindrome | $0$ | $0$ | $0$ |
| CM-$\eta$ | $(e_{(2,1)}, e_{(2,2)}, e_{(2,2)}, e_{(2,1)})$ palindrome | $-59/786432$ | $+1/8192$ | $+1/8192$ (FIXED) |
| CM-$\eta$ | $(e_{(1,0)}, e_{(1,1)}, e_{(1,1)}, e_{(1,0)})$ palindrome | $+3/65536$ | $+3/65536$ | $+3/65536$ |

Bit-exact pull-back closure identity: $\frac{1}{196608} + \left(-\frac{1}{49152}\right) = -\frac{1}{65536}$ with cochain-level pull-back increment $\Delta(B\phi_4)_{3 \to 2} = -1/(3 \cdot 2^{14}) = -1/49152$. The pro-system cochain-level pull-back $P^*_{3 \to 2}[\phi_4^{(3)}]$ does NOT equal $\phi_4^{(2)}$ on common inputs; the increment is bit-exact and explains the drift from $n_{\max} = 2$'s "boundary cutoff" $\pm 1/196608$ to $n_{\max} \ge 3$'s "interior cutoff" $\pm 1/65536$ fixed point.

Three independent witnesses of strict-strong-form failure:
1. Cutoff-stable failure $n_{\max} \ge 3$ (bit-exact fixed point).
2. Bit-exact pull-back closure identity.
3. CM-$\eta$ also fails at degree 3 (retiring v3.59.0's degree-1-only-clean cochain-clean reading).

Empirical confirmation of the published-open Connes 1994 Ch. IV / Loday §1.4 gap on entire-cyclic functoriality under algebra truncations.

Driver: `debug/compute_q5p_strict_strong.py`. Data: `debug/data/sprint_q5p_strict_strong.json`. Memo: `debug/sprint_q5p_strict_strong_memo.md`.

## 4. Paper edits applied

- **Paper 32 §VIII:** two new remarks applied sequentially in canonical order:
  - `rem:q5p_stage2_hopf_substrate` (Track A) inserted after `rem:q5p_prosystem_functoriality`.
  - `rem:q5p_strict_strong_drift` (Track B) inserted after `rem:q5p_stage2_hopf_substrate`.
  - One-sentence cross-reference update to the closing of `rem:q5p_prosystem_functoriality` (the v3.60.0 honest-scope caveat now points to the new `rem:q5p_strict_strong_drift` for the bit-exact dichotomy).
  - Pages: 67 → 68 (three-pass clean, zero undefined refs, zero errors).
- **Paper 55 §subsec:open_m2_m3:** two new \emph{emph}-prefixed paragraphs applied sequentially after the v3.60.0 pro-system paragraph and before the "Honest scope" paragraph:
  - Track A's Stage-2 Hopf paragraph.
  - Track B's strict-strong-form drift paragraph.
  - One-sentence cross-reference update to the closing of the v3.60.0 pro-system paragraph (points to the new strict-strong-form drift paragraph for the bit-exact dichotomy).
  - Pages: 24 → 25 (three-pass clean).
- **Paper 18:** no edit (Mellin engine §III.7 is upstream of both the Stage-2 Hopf substrate question and the cochain-morphism functoriality question).

## 5. Verification gate compliance (§13.4)

- **Test gate ✓** — No production code touched. Bit-exact rational arithmetic throughout both tracks. 437 axiom + truncation cells (Track A) + degree-3 closure residual panel at $n_{\max} \in \{2, 3, 4\}$ for both JLO and CM-$\eta$ (Track B) + the bit-exact pull-back closure identity verified independently.
- **Dead-end gate ✓** — Neither track matches any entry in §3 failed approaches. Both are genuinely new cohomological work.
- **Prime directive gate ✓** — Neither track modifies discrete structure (quantum number labeling, selection rules, channel structure, Gaunt couplings). All work at the cohomological-class and cochain-symbol level.
- **Consistency gate ✓** — Both tracks consistent with v3.60.0 class-level closure (extend rather than contradict). Track A's abelian primitive structure is forced by the sector-disjointness of CH idempotents (commutative algebra). Track B's bit-exact closure identity (cf. § 3) is a structural consistency check between the $n_{\max} = 2$ and $n_{\max} \ge 3$ residuals.
- **Equation gate ✓** — Track A: all five Hopf axioms verified bit-exactly at $n_{\max} \in \{2, 3\}$ (52 + 104 + 104 + 9 + 42 + 126 = 437 cells). Track B: closure residuals bit-exact at $n_{\max} \in \{2, 3, 4\}$ + pull-back closure identity. No new equations going into papers without bit-exact verification.

## 6. Honest scope

**Closed at theorem grade (bit-exact at finite cutoff):**

- Track A: Hopf axioms (coassociativity / counit / antipode / bialgebra compatibility) at $n_{\max} \in \{2, 3\}$; Mellin-slot $k$-grading preservation by $\Delta$; truncation $P_{n+1 \to n}$ as Hopf-algebra homomorphism.
- Track A: explicit $U^{*(n_{\max})}_{\mathrm{GeoVac}} = \mathbb{G}_a^{3 N(n_{\max})}$ with M1/M3/M2-respecting product structure.
- Track B: bit-exact strict/non-strict dichotomy at the cochain-morphism level for both JLO and CM-$\eta$ towers at $n_{\max} \in \{2, 3, 4\}$ on the diagnostic sub-panel.
- Track B: bit-exact pull-back closure identity $\frac{1}{196608} + \left(-\frac{1}{49152}\right) = -\frac{1}{65536}$ explaining the drift mechanism.
- Joint: $\mathcal{H}_{\mathrm{GV}}^{[k]}$ tensor factorisation is structurally forced by the M1/M2/M3 case-exhaustion partition (Paper 32 §VIII Theorem `thm:pi_source_case_exhaustion`) at the Stage-2 substrate level.

**Structural sketch (not yet theorem):**

- The forced-reading synthesis that non-abelian Stage-2 content requires substrate enrichment via one of Track A's three flagged ingredients (nested-Hopf-tower / cross-shell off-diagonal Dirac / JLO-CM Cuntz extension). This is a sharpened reading, not a no-go theorem.
- Tauberian rate uniformity in a neighborhood of $s = 3/2$ remains open (Karamata 1962 / Korevaar 2004 precedent, named in v3.59.0).

**Named open follow-ons (carrying forward):**

- M3 continuum residue (analog of v3.59.0 Track 2's M1/M2 closure on the CM-$\eta$ side; quarter-integer Hurwitz shifts; Paper 28 T9 precedent). Sprint-scale.
- Tauberian rate uniformity. Sprint-scale.
- Multi-year strict-strong-form cochain-level pro-functoriality at general cohomology degree (Connes 1994 Ch. IV / Loday §1.4 published gap).
- Multi-year Stage 2 full Tannakian construction with non-abelian pro-unipotent content via substrate enrichment.

**Hard prohibitions check (§13.5):** No changes to natural geometry hierarchy. No fitted/empirical parameters introduced. No deletion of negative results from §3. No removal of "conjectural" label from Paper 2 combination rule.

**Curve-fit audit (`feedback_audit_numerical_claims`):** Track A's Hopf axioms are bit-exact zero on derivation (not curve-fit); the abelian-primitive shape is *forced* by commutativity + sector-disjointness, not chosen. Track B's residuals are derived bit-exact computations against Sub-Sprint 2c ground truth; the closure identity $\frac{1}{196608} + \left(-\frac{1}{49152}\right) = -\frac{1}{65536}$ is a bit-exact algebraic verification, not curve-fit. Selection bias: the verdict gates were articulated BEFORE running computations; both tracks landed at the strongest gate option (POSITIVE / POSITIVE-STRUCTURAL-DRIFT), not BORDERLINE / STOP.

**Discrete-for-skeleton (`feedback_discrete_for_skeleton`):** Both tracks bit-exact `sympy.Rational` throughout; zero floats, zero PSLQ, zero transcendentals introduced at finite cutoff. Continuum-Mellin transcendentals (Track 2 of v3.59.0's $\sqrt{\pi}/2$) untouched.

**Tag transcendentals (`feedback_tag_transcendentals`):** Zero transcendentals appear in this sprint at finite cutoff. The Stage-2 $U^*_{\mathrm{GeoVac}}$ structure is intrinsically integer-graded / $\mathbb{Q}$-rational at the substrate level; transcendental content enters only at the continuum periods level via the case-exhaustion theorem (Paper 32 §VIII).

## 7. Files

### Memos
- `debug/sprint_q5p_hard_parts_2026_06_05_memo.md` (this umbrella)
- `debug/sprint_q5p_stage2_hopf_memo.md` (Track A detail)
- `debug/sprint_q5p_strict_strong_memo.md` (Track B detail)

### Drivers (new this sprint)
- `debug/compute_q5p_stage2_hopf.py` (Track A; ~650 lines, 0.02 s wall, bit-exact)
- `debug/compute_q5p_strict_strong.py` (Track B; ~520 lines, ~1.1 s wall for $n_{\max} \in \{2, 3, 4\}$ on diagnostic sub-panel)

### Data (new this sprint)
- `debug/data/sprint_q5p_stage2_hopf.json` (Track A; exact rational axiom + truncation panel)
- `debug/data/sprint_q5p_strict_strong.json` (Track B; degree-3 closure residual panel + pull-back identity)

### Paper files modified
- `papers/group1_operator_algebras/paper_32_spectral_triple.tex` (+2 remarks, +1 cross-reference update, 67 → 68 pages)
- `papers/group3_foundations/paper_55_periods_of_geovac.tex` (+2 paragraphs, +1 cross-reference update, 24 → 25 pages)

## 8. Stage 1 / Stage 2 status update

**Stage 1 status:** bit-exactly constructed across SEVEN diagnostic axes:
1. Symbol level (v3.58.0 Sub-Sprint 1: $(16, 84, 36)$ at $n_{\max} = 2$).
2. Polynomial closed forms across cutoffs (v3.58.0 Sub-Sprint 2a: $M_1, M_2, M_3$ as degree 3, 5, 4 polynomials).
3. Mixed JLO + CM-residue cochain framework (v3.58.0 Sub-Sprint 1 + v3.59.0 Track 1).
4. Mellin extraction localization (v3.58.0 Sub-Sprint 2b + v3.59.0 Track 2: M1 at $s = 3/2$ pole, M2 at integer $s$, M3 in $\eta$-pairing).
5. Bicomplex cocycle class (v3.58.0 Sub-Sprint 2c JLO HP-even $(+2, -2, +2, +2, -4)$ + v3.59.0 Track 1 CM-$\eta$ $(3, 3, 5, 15, 10)$).
6. Cutoff pro-system functoriality at the CLASS level (v3.60.0 Pro-system: sector-local closed forms $\chi_{(n, l)}, \eta_{(n, l)}$; bit-exact pull-back).
7. **Cochain-morphism level structural drift / dichotomy (THIS SPRINT Track B):** strict at class level, bit-exactly non-strict with $n_{\max} = 3$ fixed-point drift at cochain-morphism level.

**Stage 2 status:** multi-year, now with concrete bit-exact substrate (THIS SPRINT Track A):
- Candidate Hopf algebra $\mathcal{H}_{\mathrm{GV}}^{(n_{\max})}$ + candidate motivic Galois group $U^{*(n_{\max})}_{\mathrm{GeoVac}} = \mathbb{G}_a^{3 N(n_{\max})}$ explicit.
- Abelian primitive at current substrate; substrate enrichment is the multi-year continuation path.
- Three candidate enrichment ingredients flagged structurally.

## 9. One-line verdict

Two parallel sub-sprints closed the two genuinely "hard" follow-ons of the cosmic-Galois $U^*$ program: Track A delivered the first scoping step of Stage 2 with a concrete bit-exact candidate Hopf algebra + motivic Galois group ($U^{*(n_{\max})}_{\mathrm{GeoVac}} = \mathbb{G}_a^{3 N(n_{\max})}$, M1/M3/M2-tensor-decomposed) and Track B sharpened the v3.60.0 honest-scope caveat into a bit-exact strict/non-strict dichotomy at the cochain-morphism level ($\pm 1/196608 \to \pm 1/65536$ fixed point from $n_{\max} = 3$, with bit-exact pull-back closure identity); together, the abelian-primitive substrate is forced by the cochain-morphism non-functoriality, with non-abelian Stage-2 content requiring substrate enrichment via one of three flagged structural ingredients. Stage 1 of the cosmic-Galois $U^*$ bridge now bit-exactly constructed on seven diagnostic axes; Stage 2 has a concrete bit-exact starting object.
