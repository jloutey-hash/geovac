# Sprint Hodge-SL₂ — canonical memo

**Date:** 2026-06-14. **Supersedes** `debug/sprint_hodge_sl2_scoping_memo.md`
(the in-flight plan).
**Verdict:** **BORDERLINE (CM)** in the decision-gate sense — but a
POSITIVE-grade *result*: it is more explanatory than the "SL₂ = MT"
outcome the sprint originally chased.

**Drivers:** `debug/sprint_hodge_sl2_step1_complex_structure.py`,
`..._step23_mumford_tate.py`, `..._audit_qi_triangulation.py`.
**Data:** `debug/data/sprint_hodge_sl2_{step1,step23,audit_qi}.json`.
**Discipline:** all linear algebra exact `sympy.Rational`; period checks
`mpmath` @ 40 dps; curve-fit audit applied to the Q(i) triangulation.

---

## 1. Headline

GeoVac's $SL_2$ (the reductive Levi factor of
$U^* = \mathbb{G}_a^\infty \rtimes SL_2$, forced by $S^3 = SU(2)$) is
**NOT** the Mumford–Tate group of GeoVac's canonical weight-1 Hodge
structure. GeoVac realizes the **complex-multiplication point** with CM
field **$\mathbb{Q}(i)$**: the MT group is the 1-dimensional CM torus
inside $SL_2$, not the full 3-dimensional $SL_2$.

Three payoffs make this the rich outcome:
1. It **structurally explains the Hain–Brown period-negative**: Hain–Brown
   is the generic (full-$SL_2$, non-CM, modular) story; GeoVac sits at the
   CM point where modular forms degenerate, so GeoVac's periods are
   CM-pure over $\mathbb{Q}(i)$ (Catalan, $\beta$-values), not modular.
   The June period-level negative *is* the CM degeneration.
2. It **sharpens Reading A**: GeoVac doesn't merely drop iteration — it
   sits at the most arithmetically special point, where the symmetry
   group collapses from $SL_2$ to a CM torus.
3. The **$\mathbb{Q}(i)$ triangulation** (Hodge CM field = level-4 field
   = period field) is certified by audit as **structurally forced by
   half-integer spin**, not coincidence.

Plus a **paper-correction**: Paper 56 §sec:open_g4_hodge's claim that the
"$J^2=-1$, $JD=+DJ$" KO-dim-3 real spectral triple is *verified* at finite
cutoff is not bit-exact (see §2 step 1).

## 2. Steps and results (all bit-exact unless noted)

**Step 1 — complex structure exists; KO-dim-3 partial.** The Kramers J
($J^2 = -I$, the SU(2) quaternionic structure) exists bit-exact at
$n_{\max} \in \{2,3\}$. With CG-weighted adjacency the off-diagonal part
of $D$ anticommutes with $J$ **exactly** ($\{J, A_{CG}\} = 0$) — closing
the codebase's flagged "future work." But the **diagonal** CH part
*commutes* with $J$ ($[J,\Lambda]=0$, since $\Lambda$ is $m_j$-blind),
so $D = \Lambda + \kappa A$ satisfies **neither** $JD=+DJ$ nor $JD=-DJ$;
the residual is exactly $2\Lambda J$. The diagonal and CG-off-diagonal
pieces sit in **opposite $J$-parity sectors**. Hence: the complex
structure exists, but the clean KO-dim-3 real spectral triple does not
hold on the natural finite substrate. (Verdict: PARTIAL; informative.)

**Steps 2–3 — polarization + Mumford–Tate.** Restricting Kramers J to
the $n=1$, $\kappa=-1$, $j=1/2$ doublet gives the canonical complex
structure $J = \begin{psmallmatrix}0&-1\\1&0\end{psmallmatrix}$ (rational,
$J^2=-I$, $J$-invariant block). With the symplectic polarization
$Q = \begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix}$: $J$ is symplectic
and the Riemann form $QJ = I$ is positive-definite → a **genuine
polarized weight-1 Hodge structure**. Mumford–Tate group:
$\mathrm{minpoly}(J) = t^2+1 \Rightarrow \mathbb{Q}[J] \cong \mathbb{Q}(i)$;
the Hodge circle $\cos t\,I + \sin t\,J$ generates the norm-1 torus of
$\mathbb{Q}(i)$ — a **1-dimensional** CM torus (vs $\dim SL_2 = 3$).
Sanity: a $3/4/5$ Pythagorean rotation is a genuine $SL_2(\mathbb{Q})$
point of this torus acting via the standard rep. **Verdict: the MT group
is the CM torus by $\mathbb{Q}(i)$, NOT $SL_2$.** (Robust across modeling
choice: any $j=1/2$ doublet gives a rational $J$ with $J^2=-I$, hence
$\mathbb{Q}(i)$; $SL_2$ would require a transcendental $J$, which
algebraic GeoVac data cannot produce.)

**Audit — Q(i) triangulation is structurally forced.** The claim "Hodge
CM field = level-4 field = period field = $\mathbb{Q}(i)$, for one
reason" passes a structural (not numerical) audit with a falsifiable
null:
- (A) Hodge spin-gating: $\varepsilon = (-1)^{2j} = -1$ for all 16 spinor
  states (complex structure exists); scalar $(-1)^{2l} = +1$ (no complex
  structure).
- (B) the two 4's: $\mathrm{disc}\,\mathbb{Q}(i) = -4$ and
  $\mathrm{cond}\,\chi_{-4} = 4$ — same $4 = \mathbb{Q}(\mu_4)$.
- (C) period spin-gating: $\beta(2) = $ Catalan (conductor 4,
  $\mathbb{Q}(i)$, from the Dirac vertex parity, Paper 28) vs scalar
  $\zeta(2) = \pi^2/6$ (conductor 1, pure-Tate $\mathbb{Q}$).
The scalar null removes $\mathbb{Q}(i)$ from **both** sides
simultaneously → **common cause = half-integer spin**, not coincidence.
Honest scope: level-4 $= \mathbb{Q}(\mu_4)$ is literature (Deligne 2010 +
EMN 2025); the audit verifies the common spin cause + same-field
identity, it does not re-derive the motivic level-4 result.

## 3. Paper changes (Paper 56 §sec:open_g4_hodge)

1. **Correction:** the "$J^2=-1$, $JD=+DJ$ verified" KO-dim-3 claim →
   the honest finite-substrate statement (complex structure exists;
   clean JD relation does not; opposite-parity diagonal vs CG-offdiagonal).
2. **Upgrade convention → identification:** new Proposition — GeoVac's
   canonical weight-1 Hodge structure is CM-type with CM field
   $\mathbb{Q}(i)$; its MT group is the CM torus inside $SL_2$;
   structurally explains the Hain–Brown negative; ties the Hodge
   realization to the level-4 field $\mathbb{Q}(i)$ via the
   spin-forced triangulation.

## 4. Honest scope / what this does NOT claim

- Does NOT claim $SL_2$ = MT (it isn't; the verdict is CM torus).
- Does NOT re-derive the motivic level-4 identification (literature).
- The modeling choice ($V_{\mathrm{fund}}$ = $j=1/2$ Dirac doublet,
  complex structure = restricted Kramers J) is canonical but a choice;
  the alternative framing (abstract $SL_2$'s $SO(2)$ maximal compact)
  gives the same rational $J$ → same $\mathbb{Q}(i)$ CM verdict.
- Equality direction of the injection theorem (multi-year) untouched.
  WH1 PROVEN, WH6, Paper 2 combination rule — all unaffected.

## 5. Files

- Drivers + data as above.
- `tests/test_paper56_hodge_sl2.py` — verification panel.
- `papers/group3_foundations/paper_56_tannakian_substrate.tex`
  §sec:open_g4_hodge — correction + CM-identification proposition.
