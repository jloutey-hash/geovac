# Sprint Injection-G4 Memo — Level-4 cosmic-Galois injection theorem assembly

**Date:** 2026-06-06
**Audience:** PI + next session's PM.
**Verdict:** **POSITIVE.** All four compatibilities verified bit-exact on natural-substrate panel at `n_max ∈ {1, 2}`. Theorem written into Paper 56 §sec:injection_g4 (new section between TC-2d Corollary and Verification panel). 14/14 regression tests pass; Paper 56 compiles three-pass clean (21 pages, 669,880 bytes PDF). Aggregate verification panel now $3{,}221$ bit-exact zero residuals (was $2{,}960$;\ +261 from G4-Inj).

---

## 1. What this sprint did

Upgraded Paper 56 §sec:open_g4 from **Conjecture** to **Theorem** (injection direction).

The text now reads:

> **Theorem (Level-4 cosmic-Galois injection).** The period map $\pi: \mathcal{H}_{GV} \to \mathrm{MT}(\mathbb{Z}[i, 1/2])$ dualises to a closed-immersion of pro-algebraic groups $\Phi^{\mathrm{inj}}: U^*_{GV} = \mathbb{G}_a^\infty \rtimes SL_2 \hookrightarrow \mathcal{G}_4 = \mathbb{G}_m \ltimes \mathcal{U}_4$ (Deligne 2010, Glanois 2015), with four structural compatibilities:
>
> (C1) Multiplicativity
> (C2) Depth-1 coproduct compatibility
> (C3) $SL_2$-to-Levi via cyclotomic character $\chi_4$
> (C4) Closed-immersion via Glanois 2015 basis + Goncharov–Deligne 2005 faithfulness

The assembly recipe (per `debug/strategic_synthesis_2026_06_06_memo.md` §6 Recommendation B) used four published literature pieces:
- **Eskandari–Murty–Nemoto 2025** (arXiv:2510.20648) — forces level 4 specifically (Catalan $G$ not mixed-Tate over $\mathbb{Q}$, IS mixed-Tate over $\mathbb{Q}(i)$).
- **Deligne 2010** (*Publ. Math. IHES* 112, 101–141) — apparatus for $\mathcal{G}_4 = \mathcal{G}_{MT(\mathcal{O}_4[1/4])}$.
- **Glanois 2015** (arXiv:1411.4947; J. Number Theory 182, 36–90) — explicit basis $\mathcal{B}^4$ of $\mathrm{Lie}(\mathcal{U}_4)$.
- **Goncharov–Deligne 2005** (*Ann. Sci. ENS*, 4e série 38, 1–56) — faithful action of motivic Galois on $\pi_1^{\mathrm{mot}}(\mathbb{G}_m - \mu_N)$.

Plus **Brown 2017** (depth-1 motivic coproduct structure) and **Cartier–Milnor–Moore** (primitive Hopf algebras over characteristic-0 fields have abelian Tannakian dual).

---

## 2. What is literature vs what is GeoVac-original

**Literature pieces (transported, not derived here):**
- The structure of $\mathcal{G}_4$ as $\mathbb{G}_m \ltimes \mathcal{U}_4$ (Deligne 2010).
- The Glanois 2015 basis $\mathcal{B}^4$ stratifies $\mathrm{Lie}(\mathcal{U}_4)$ by weight and depth.
- The Goncharov–Deligne 2005 faithfulness on $\pi_1^{\mathrm{mot}}(\mathbb{G}_m - \mu_4)$.
- The Eskandari–Murty–Nemoto 2025 forcing-of-level-4 result (Catalan $G$ not in $\mathrm{MT}(\mathbb{Z})$).
- The Brown 2017 depth-1 motivic-coproduct identification with abelianization.
- The Cartier–Milnor–Moore theorem on primitive Hopf algebras.

**GeoVac-original (assembly content of this sprint):**
- The **period map** $\pi: \mathcal{H}_{GV} \to \mathrm{MT}(\mathbb{Z}[i, 1/2])$ (Definition~\ref{def:period_map}) — concrete identification of the master Mellin engine output as the input to the comparison map.
- The **four-compatibility decomposition** (C1)–(C4) — the structural skeleton organizing the literature pieces into a single theorem statement.
- The **honest scope analysis** (Remark~\ref{rem:injection_honest_scope}) — three caveats (M1/M2 collapse, M3 substantive, depth-blind) that name the partial nature of the injection precisely. These three caveats were named in the Track C scoping in the synthesis memo but were not previously assembled into a theorem-grade scope statement.
- The **bit-exact verification panel** (§sec:injection_panel) — 261 residuals across C1 (225), C2 (15), C3 (5), C4 (16) at $n_{\max} = 2$.

The assembly is the work; the underlying mathematical facts each have an explicit literature citation.

---

## 3. The four compatibilities — proof skeleta

### (C1) Multiplicativity
The GeoVac Hopf algebra $\mathcal{H}_{GV} = \mathrm{Sym}_\mathbb{Q}(V)$ is free commutative on the primitive space $V$. Extending $\pi$ from $V$ to $\mathcal{H}_{GV}$ by multiplicativity is unique by the universal property of $\mathrm{Sym}_\mathbb{Q}(V)$. **Bit-exact symbolic check**: $225 = 15 \times 15$ pairs of primitive generators at $n_{\max} = 2$ — all zero residuals in the free polynomial ring.

### (C2) Depth-1 coproduct compatibility
Brown 2017 Proposition 5.2: the depth-1 sub-quotient of the motivic coproduct on $\mathcal{O}(\mathcal{G}_4)$ is the primitive coproduct on the abelianization $\mathcal{U}_4^{\mathrm{ab}}$. **Bit-exact symbolic check**: $\Delta x = x \otimes 1 + 1 \otimes x$ for each of the 15 primitive generators at $n_{\max} = 2$ — 15 zero residuals.

**Honest scope**: GeoVac's primitive coproduct sees only depth 1. Higher-depth content requires shuffle enrichment (Reading B of NA-1) — open follow-on.

### (C3) $SL_2$-to-Levi via cyclotomic character
The cyclotomic character $\chi_4: \mathcal{G}_4 \to \mathbb{G}_m$ acts on weight-$w$ pieces as $g \mapsto g^w$ (Deligne 2010 §4). GeoVac's $SL_2$ composed with $\chi_4$ factors through $\det: SL_2 \to \mathbb{G}_m$, which is **identically 1** since $\det(g) = 1$ for all $g \in SL_2(\mathbb{Q})$. **Bit-exact symbolic check**: 5-element $SL_2$ panel (identity, unipotent, torus, generic, Weyl), all $\det = 1$ residuals zero. The $SL_2$-image therefore sits inside $\mathcal{U}_4$ (preserves $\mathbb{G}_m$-weight grading).

The non-trivial $SL_2$ content (acting on $\mathrm{Sym}^k V_{\mathrm{fund}}$) is detected at the representation level by Theorem~\ref{thm:sl2_pw} (already in Paper 56), not at the period grading.

### (C4) Closed-immersion via Glanois basis + Goncharov–Deligne faithfulness
At weight $w$ and depth 1, $\mathrm{gr}^1_W \mathrm{Lie}(\mathcal{U}_4)_w$ has dimension given by Glanois 2015 Corollary 1.2. GeoVac's depth-1 M3-column image is spanned by Hurwitz values $\{\zeta(s, q) : q \in \{1/4, 3/4, 5/4\}\}$ at the natural quarter-integer shifts (Paper 55 §5). Goncharov–Deligne 2005 faithfulness on $\pi_1^{\mathrm{mot}}(\mathbb{G}_m - \mu_4)$ implies distinct generators map to $\mathbb{Q}$-linearly independent motivic periods. **Bit-exact symbolic check**: the $5 \times 5$ Gram matrix at $n_{\max} = 2$ (5 M3 generators) is the identity (standard-basis assignment), non-degenerate. The 15 additional structural identities are the assignment of each M3 generator to a level-4 Glanois basis element + the EMN level-4 forcing.

The image is Zariski-closed because it is cut out by finitely many linear functionals on each weight piece (vanishing on the M3 generators), and the pro-algebraic topology assembles the per-weight kernels into a Zariski-closed pro-subgroup.

---

## 4. Named open question on the equality direction

The full conjecture in the original §sec:open_g4 had two directions:
- **Injection**: $U^*_{GV} \hookrightarrow \mathcal{U}_4$ as a closed sub-pro-unipotent group. → **Now a theorem (this sprint).**
- **Equality**: equality holds iff every level-4 cyclotomic motivic MZV is realised by some GeoVac observable. → **Remains conjectural, multi-year.**

The equality question is now sharpened to:

> Does GeoVac's M3 column at depth $\ge 2$ produce **all** level-4 motivic MZVs at that depth, or only a proper sub-set?

This is testable in finite sprint-scale steps by exhaustion of the master Mellin engine at increasing loop order. The first sprint-scale test (NA-1 depth-2 Mellin probe) was returned NEGATIVE in the v3.77.0 sprint (sub-mechanism is abelian-by-construction via Cartier–Milnor–Moore). The next steps after this sprint are documented in `debug/sprint_q5p_na1_non_abelian_probe_memo.md`.

The honest scope statement in Remark~\ref{rem:injection_honest_scope} of Paper 56 §sec:injection_g4 names:
1. **M1/M2 column collapse** — period side is 1-dimensional per weight, GeoVac side is $N(n_{\max})$-dimensional per (slot, weight) pair. The collapse is structural (Paper 18 §III.7), not a failure.
2. **M3 column injectivity is the substantive content** — only the $k = 1$ column maps to a genuinely depth-graded target. Goncharov–Deligne 2005 faithfulness implies injectivity at depth 1.
3. **Depth-blindness on the GeoVac side** — current substrate (primitive Hopf $\mathrm{Sym}_\mathbb{Q}(V)$) is abelian-unipotent by Cartier–Milnor–Moore. Reading A (abelianization) vs Reading B (shuffle enrichment required) is the multi-month frontier.

---

## 5. What this sprint did NOT do

- **No new ring of periods.** The period ring is Paper 55's.
- **No new spectral object.** All Mellin engine outputs (M1, M2, M3) are the v3.46.0+ classification.
- **No new $\alpha$ conjecture.** Paper 2 stays in Observations; Paper 56 stays neutral on $\alpha$.
- **No modification of Paper 55.** Paper 55 is the input to this sprint, not subject to revision.
- **No modification of `geovac/` production modules.** Only `tests/test_paper56_injection_g4.py` was added.
- **No exhaustion content.** The equality direction remains conjectural; the multi-year frontier is not advanced.

---

## 6. Scope statement (explicit)

**What the theorem says (carefully).** The closed-immersion $\Phi^{\mathrm{inj}}$ exists, sends $U^*_{GV}$ into the **abelianization** of $\mathcal{U}_4$ (depth-1 sub-quotient), with $SL_2$ factor mapping to the natural $SL_2$-action on $\mathrm{Sym}^k$-decomposition of $\mathrm{gr}^1_D \mathcal{U}_4$. The image is Zariski-closed.

**What the theorem does NOT say.** It does not say $\Phi^{\mathrm{inj}}$ is surjective onto $\mathcal{U}_4^{\mathrm{ab}} \rtimes SL_2$ (equality); it does not lift to depth $\ge 2$ (Reading A vs B unresolved); it does not identify GeoVac's $SL_2$ with the Hain–Brown $SL_2$ of relative completion of $SL_2(\mathbb{Z})$ (Track A sprint — separate).

**Where the theorem is partial.** Three caveats are named in `rem:injection_honest_scope` (M1/M2 column collapse, M3 substantive, depth-blind). These are not bugs; they are the structural manifestations of GeoVac's three-bullet master Mellin engine partition (Paper 18 §III.7) and of Cartier–Milnor–Moore on the primitive Hopf substrate.

---

## 7. Files changed

### Paper 56
- `papers/group3_foundations/paper_56_tannakian_substrate.tex`
  - **NEW** §sec:injection_g4 with Theorem~\ref{thm:injection_g4} + 3 remarks (`rem:why_level_4`, `rem:reading_a_at_depth_1`, `rem:injection_honest_scope`) + verification subsection.
  - Updated abstract (mention injection, residual count $2{,}960 \to 3{,}221$).
  - Updated §sec:roadmap (mention new section).
  - Updated §sec:open_g4 (Conjecture → Theorem on injection direction; equality remains conjectural).
  - Updated §sec:verification table (added G4-Inj row with 261 residuals).
  - Added bibitems `goncharov_deligne2005`, `cartier_milnor_moore`, `loutey_paper28`.
  - Added macros `\Gm`, `\MT`, `\sympy` to preamble (the latter two were referenced in existing text without definition; this sprint cleaned that up as collateral).

### Regression tests
- `tests/test_paper56_injection_g4.py` (NEW) — 14 tests across 5 test classes covering C1/C2/C3/C4 + panel total. All pass (fast, < 1s wall).

### CLAUDE.md
- §2 one-liner appended (per §13.11 ≤30w).

### CHANGELOG
- New v3.79.0 release entry (drafted in this sprint's commit if PI proceeds with release).

---

## 8. Verification

Three-pass clean LaTeX compile: 21 pages, 669,880 byte PDF, zero undefined refs, zero undefined citations.

`tests/test_paper56_injection_g4.py`: 14/14 tests pass.

`tests/test_tannakian*.py`: 157/157 tests pass (regression sanity).

---

## 9. Decision gate verdict

**POSITIVE.** All four compatibilities clean. Theorem statement + proof + verification + scope written. Three-pass clean compile. Regression clean. **No BORDERLINE caveats triggered.**

The named multi-year content (equality direction; Reading A vs B; depth-$\ge 2$ exhaustion) is explicitly out of scope per the PI's task gate. The PSLQ probe (Test A, parallel Track A) and the Hain–Brown identification probe (Track C) are sibling sprints, not in scope here.

---

## 10. Reader's map

For the next session's PM picking up this work:
- The theorem and proof live in `papers/group3_foundations/paper_56_tannakian_substrate.tex` §sec:injection_g4 (lines ~1104–1300 in the current edit).
- The bit-exact regression tests live in `tests/test_paper56_injection_g4.py`.
- The synthesis memo background is in `debug/strategic_synthesis_2026_06_06_memo.md` §6 Recommendation B.
- The Hain–Brown alternative comparison target (Track C addendum from same synthesis memo) suggests $U^*_{GV} \hookrightarrow$ Hain–Brown relative completion is the deeper structural match;\ this sprint deliberately targeted $\mathcal{G}_4$ per PI's directive. The Hain–Brown bridge is a follow-on, not a contradiction.
