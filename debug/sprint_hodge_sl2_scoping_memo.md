# Sprint Hodge-SL₂ — scoping memo + execution plan

**Date:** 2026-06-14
**Status:** SCOPED, execution begun (step 1 driver in flight).
**Charter:** Test whether GeoVac's SL₂ (the reductive Levi factor of
$U^* = \mathbb{G}_a^\infty \rtimes SL_2$, forced by $S^3 = SU(2)$,
acting via standard + Sym² on the Peter–Weyl panel) is the
**Mumford–Tate group of a polarized weight-1 Hodge structure** on
$V_{\mathrm{fund}} = \mathbb{Q}^2$. Goal: upgrade Paper 56
§sec:open_g4_hodge from an adopted *convention* (Hodge / Hodge-de-Rham
realization) to an actual *identification* — or return an informative
NEGATIVE that weakens the convention.

**Provenance:** named + deferred three times on 2026-06-06 (strategic
synthesis §5 item 3 "Recommendation C"; §2 Track B; NA-1 probe memo),
never executed (no prior sprint memo exists — confirmed by keyword
sweep). **Excludes** the multi-year equality/exhaustion wall (PI
direction 2026-06-14).

---

## 1. The question, concretely

GeoVac's SL₂ is currently a purely *algebraic* symmetry — it shuffles
representations (V, Sym²V, …). The **Mumford–Tate group** is the
symmetry group of an actual *Hodge structure* — the kind of object that
is literally the $H^1$ of an elliptic curve. The sprint asks whether
our algebraic SL₂ is secretly *that*.

Standard dictionary: a polarized weight-1 Hodge structure on
$\mathbb{Q}^2$ ⟺ a Hodge circle $h: U(1) \to SL_2(\mathbb{R})$, i.e.

- a **complex structure** $J_h$ on $\mathbb{R}^2$ with $J_h^2 = -1$, plus
- an $SL_2$-invariant **symplectic polarization** $Q$ with the Riemann
  positivity $Q(x, J_h x) > 0$.

The **Mumford–Tate group** is the smallest $\mathbb{Q}$-algebraic
subgroup of $SL_2$ whose real points contain $\mathrm{im}(h)$. Generic
$h \Rightarrow MT = SL_2$; degenerate $\Rightarrow$ a proper torus
(CM type).

## 2. What it operates on (grounded in code, this session)

- **SL₂ action** — `geovac/tannakian.py`: `sl2_standard_action`
  ($\rho(g) = g$ on $\mathbb{Q}^2$), `_sl2_sym2_action`
  ($\mathrm{Sym}^2(g)$ on $\mathbb{Q}^3$). Standard + Sym² is *exactly*
  the rep structure of the MT group of a weight-1 rank-2 HS.
- **Real structure** — `geovac/spectral_triple.py`, two J's:
  - `permutation` (default): $J^2 = +I$, $JD = DJ$. All Connes axioms
    pass — but $J^2 = +1$ is the **wrong square** for a complex
    structure.
  - `kramers`: $J^2 = -I$ (the SU(2) = S³ quaternionic structure — the
    **right square**), but $JD = -DJ$ fails on `uniform` adjacency; the
    docstring flags **`adjacency_weights='cg'`** (Wigner-3j-weighted,
    exact/π-free) as the fix, marked "future work."

**Load-bearing finding.** The candidate complex structure is the
Kramers J ($J^2 = -1$), and its D-compatibility ($JD = -DJ$, the KO-dim-3
real-spectral-triple condition) is an open item the codebase itself
never closed. **Step 1 resolves it.**

**Paper-correction flag.** Paper 56 §sec:open_g4_hodge cites "$J^2 = -1$,
$JD = +DJ$" as *verified* Connes axioms supplying the real-Hodge
prerequisite. But the $J^2 = -1$ (Kramers) structure does **not**
satisfy a clean JD relation on the natural uniform substrate — that is
precisely what step 1 must *establish*, not assume. The §sec:open_g4_hodge
text should be tightened regardless of the sprint outcome.

## 3. The computation — 5 steps (all small + exact)

1. **Complex structure (LOAD-BEARING).** Kramers J, confirm $J^2 = -I$;
   test $JD = -DJ$ on CG-weighted adjacency, exact. PASS ⇒ KO-dim-3 real
   structure confirmed on the natural CG substrate; a D-compatible
   complex structure exists. *(This driver: step 1.)*
2. **Polarization.** $SL_2$-invariant symplectic $Q$ on $V_{\mathrm{fund}}$;
   Riemann positivity $Q(x, J_h x) > 0$, exact / controlled precision.
3. **Mumford–Tate group.** $\mathbb{Q}$-Zariski closure of the Hodge
   circle in $SL_2$; decide $= SL_2$ (generic) vs proper torus (CM) vs
   $GL_2$.
4. **Match.** MT action on $V_{\mathrm{fund}}$, $\mathrm{Sym}^2 V_{\mathrm{fund}}$
   vs GeoVac $\rho$, bit-exact.
5. **Honest VHS identification.** Confirm the standard weight-1 VHS
   ($H^1$ of a generic elliptic curve); confirm **non-extension** to
   Hain–Brown (no pro-unipotent radical) — which structurally explains
   the period-level Hain–Brown negative and is the Hodge-side
   restatement of Reading A.

## 4. Decision gate

- **POSITIVE:** $SL_2 = MT(\text{standard polarized weight-1 HS})$,
  bit-exact. Capstone of the cosmic-Galois map.
- **BORDERLINE:** HS exists, $MT$ = proper subgroup (torus ⇒ CM-type).
- **NEGATIVE:** no canonical complex structure / positivity fails ⇒
  weakens the adopted Hodge convention. Informative either way.

## 5. Expected value (honest)

**Consolidation, not frontier.** Three facts bound the upside:
(i) Hain–Brown tested NEGATIVE at the period level (depth 1 & 2) — the
rich elliptic-enrichment door is already shut; (ii) the MT group of a
generic non-CM elliptic curve is "just" $SL_2$ with no unipotent
radical; (iii) a POSITIVE is the Hodge-side restatement of Reading A
(we carry the reductive symmetry, not the iteration). Most likely
outcome: theorem-grade *"$SL_2 = MT$ of the standard VHS; Hain–Brown
excluded because we have the reductive factor but no unipotent
radical."* That completes the structural identification on **both**
Levi factors ($\mathbb{G}_a$ = abelianized motivic unipotent via the
injection theorem; $SL_2$ = MT of the standard VHS via this sprint).
The natural capstone for the map / outreach, not a new expedition.

## 6. Effort / files / risk

- **Effort:** ~2 weeks. Core computation small + exact; the second half
  is honest VHS identification + literature grounding (Hain 2014;
  standard Mumford–Tate theory).
- **Risk (conceptual, not compute):** step 1 may fail — the Kramers J
  may not become cleanly D-compatible even with CG weights. That is the
  one outcome that *changes* the picture (weakens the Hodge convention).
  Low-probability, high-information.
- **Hygiene:** no production `geovac/` changes (uses existing
  `adjacency_weights='cg'`); no §3 dead-end re-derived; no guardrail
  paper triggered.

## 7. Day-by-day

- **D1–2:** Step 1 (Kramers J + CG, $J^2=-I$, $JD=-DJ$). *Checkpoint:
  does the HS prerequisite exist?*
- **D3–4:** Step 2 polarization + Riemann positivity.
- **D5–6:** Step 3 MT group. *Checkpoint: $SL_2$ / torus / $GL_2$?*
- **D7–8:** Step 5 honest VHS identification + Hain–Brown non-extension.
- **D9–10:** Paper 56 §sec:open_g4_hodge edit (convention → identification),
  `tests/test_paper56_hodge_sl2.py`, canonical sprint memo (replaces this
  scoping memo at close), CHANGELOG, CLAUDE.md §2 one-liner.

## 8. Files

- `debug/sprint_hodge_sl2_step{1..N}_*.py` — step drivers.
- `debug/data/sprint_hodge_sl2_*.json` — results.
- `debug/sprint_hodge_sl2_memo.md` — canonical memo at close.
- `tests/test_paper56_hodge_sl2.py` — verification.
- `papers/group3_foundations/paper_56_tannakian_substrate.tex`
  §sec:open_g4_hodge — convention → identification (or negative).
