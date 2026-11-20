# Dirac-on-S³ Tier 1 Sprint Plan

**PI approved:** 2026-04-15 (Dirac-on-S³ pivot committed; Leader + Explorer briefs synthesized).
**Target:** 3–4 weeks (likely much less if the algebraic-first philosophy pays off as it has in Phases 4B–4H and in the algebraic registry precedents).
**Framing:** curiosity-driven — "what is Coulomb, structurally, that no other force is, and can the Dirac-on-S³ sector turn K = π(B + F − Δ) from a three-tier coincidence into a theorem?" Publishable in every one of the four verdict outcomes.

---

## Governing principle for this sprint (and all future sprints)

**Algebraic-first.** If a track prompt is about to commit a worker to a long numerical computation, the PM must pause and ask: is there an algebraic structure? Before ANY long numerical run is dispatched, the Decomposer agent (`agents/DECOMPOSER.md`) is consulted or the question "what does Paper 18 say about this transcendental" is answered. Precedent: Paper 12 Neumann, Gaunt selection rules, Level 2 Laguerre, Level 3 S/K recurrence, Paper 24 Bargmann-Segal, `hypergeometric_slater.py`. Every one of those collapsed weeks of numerics to a symbolic one-liner once the right structure was written down. **Reference: Paper 18 (exchange constants taxonomy).**

The PI's expectation is that this sprint's 3–4 week estimate will collapse substantially if the algebraic route is followed. Track completions should be reported in days when possible, and a track that is taking more than a week of worker time is a signal to stop and invoke the Decomposer.

---

## PI decisions (recorded)

1. **D4 reshape: YES.** Drop "Fock rigidity theorem for Dirac" as a Tier 1 deliverable. Replace with Orduz Hopf-equivariant decomposition (Explorer Finding 3.2). Fock rigidity is deferred to Tier 2 or beyond per Explorer Gap #1.
2. **D1 builds both Weyl and full 4-component Dirac.** Let D2/D3/D4 pick which sector they need.
3. **Paper 2 framing update: defer until D5.** No edits to Paper 2 during the sprint.
4. **Tier 1b "Fock rigidity for Dirac" follow-on: deferred to Tier 2.** Explorer Gap #1 confirms no off-the-shelf theorem exists to extend; this is a multi-month research problem, not a sprint track.

---

## Proposed tracks (final)

### Track D1 — Dirac-on-S³ infrastructure module

**Goal:** Build `geovac/dirac_s3.py` with Camporesi–Higuchi spectrum and spinor spherical harmonics on S³, exposed through a (n, l, m, σ) labeling compatible with existing Fock graph nodes. Both single-chirality (Weyl, g_n = (n+1)(n+2)) and full Dirac (g_n = 2(n+1)(n+2)) sectors. All degeneracies and eigenvalues verified symbolically as integers / half-integers in exact Fraction arithmetic (π-free certificate analogous to Paper 24).

**Algebraic-first constraint:** Do NOT use floating-point numerical diagonalization to verify spectra. Camporesi-Higuchi gives closed forms — verify symbolically using sympy. If a matrix element looks like it needs integration, stop and check whether Harish-Chandra's Casimir formula or SU(2)×SU(2) decomposition gives it closed-form first.

**Success:** Module + tests, π-free in exact rational arithmetic, |λ_n| and g_n match Camporesi-Higuchi symbolically for n = 0..6 in both sectors. g_3^Dirac = 40 reproduced as exact Fraction. Spinor harmonics on S³ indexed by (n, l, m, σ) with explicit map to existing (n, l, m) Fock nodes.

**Failure:** Spinor harmonics don't align cleanly with (n, l, m) — forces parallel graph rather than extension. (Low risk per Explorer Finding 1.1/1.2: SU(2)×SU(2) decomposition is clean.)

**Deliverable:** `geovac/dirac_s3.py`; `tests/test_dirac_s3.py`; `docs/dirac_s3_design_memo.md`.

**Depends on:** Nothing. Start here.

**Guardrails to cite in the worker prompt:**
- Section 3 "Local discrete Dirac on Fock graph" (Ginsparg-Wilson — Explorer Finding 4.2). If the worker tries to write a local nearest-neighbor discrete Dirac, GW theorem forbids chirality preservation without doublers. The Dirac operator in D1 lives on the spinor bundle over S³, not on a simplicial complex — keep this distinction explicit.
- Paper 24 § III (π-free certification pattern). Mirror that certification style.

**References to provide:**
- Camporesi-Higuchi 1996 (arXiv:gr-qc/9505009)
- Bär 1996 (*J. Math. Soc. Japan* 48)
- Paper 18 §IV (where relativistic transcendentals should live in the taxonomy)

---

### Track D2 — Dirac analog of B = 42

**Goal:** Compute the finite truncated Dirac Casimir trace on S³ at the m=3 cutoff that Paper 2 uses for scalars. Does Σ |λ_n|²·g_n (or Σ |λ_n|·g_n, or Σ g_n) up to n=2 equal 42, a simple rational multiple of 42, or something unrelated?

**Algebraic-first constraint:** This is pure symbolic arithmetic over integers and half-integers. There should be zero floating point. If the worker needs more than a sympy notebook, something is wrong.

**Success:** Clean rational value ∈ {42, 42/k, k·42 for small k}, with a one-line spectral-theoretic reason. Promotes B to a Dirac-spectral invariant and is one quarter of the common-generator theorem candidate.

**Failure:** Unrelated rational (e.g., 30, 56) with no tie to 42. Confirms B lives only on the scalar Laplace–Beltrami sector. (This is a clean verdict, not a dead end.)

**Deliverable:** `debug/track_dirac_s3/D2_casimir_trace.py` + notebook-style memo with symbolic derivation.

**Depends on:** D1.

**Expected runtime:** hours, not days. If this takes more than a day, stop and check the setup.

---

### Track D3 — Dirac analog of F = π²/6

**Goal:** Compute D_{g_n^Dirac}(s) = Σ_{n≥1} g_n · n^{−s} at s = d_max = 4 (the packing exponent that gave F = ζ(2) for scalars in Phase 4F, α-J). Also test close variants: s ∈ {3, 5, 6}, Weyl vs. full, degeneracy-weighted vs. uniform. Does any variant give a closed form in ζ-values that links structurally to F = π²/6?

**Algebraic-first constraint:** Use sympy's `summation` and `zeta` symbolically. Target: exact equalities, not numerical coincidences. A "match to 10⁻⁸" is not a result; an exact sympy equality IS a result. If the series doesn't admit closed form, try regularized differences of ζ(s−k) as in Phase 4F.

**Success:** Closed-form symbolic equality at some natural s, especially if it equals F, 2F, F/2, or π²/3. Promotes F to a Dirac-spectral invariant.

**Failure:** No variant gives clean closure; series are genuine mixtures of multiple ζ-values with no tie to F. Keeps F on the scalar sector.

**Deliverable:** `debug/track_dirac_s3/D3_dirichlet_series.py` + memo extending Phase 4F's data artifact.

**Depends on:** D1.

**Expected runtime:** hours. Phase 4F's F = ζ(2) identification was a one-liner once the right sum was written down.

---

### Track D4 — Orduz Hopf-equivariant Dirac decomposition (RESHAPED per Explorer Finding 3.2)

**Goal:** Using the Orduz decomposition of the S³ Dirac operator as a tower of U(1)-charge-q-twisted Dirac operators on S² (base) with line-bundle-twist q (fiber), compute the charge-q decomposition of the 40 states at n=3 for the full Dirac case. Does the decomposition over q ∈ {−3, −2, −1, 0, 1, 2, 3} yield a structural link to Paper 2's B = 42 (via partial sums), F = π²/6 (via base/fiber Dirichlet series), or a new invariant that unifies with the scalar Hopf construction from Paper 2?

**Algebraic-first constraint:** Twisted Dirac spectra on S² are off-the-shelf (Baum–Friedrich, Bär). Use symbolic Fraction arithmetic throughout. The question is combinatorial: how does 40 = 2·4·5 split over fiber charges, and does the resulting partition connect to 42 or ζ(2)?

**Success:** Structural identity linking the Dirac ν=3 Hopf decomposition to Paper 2's B/F/Δ in exact arithmetic. If this works, it is the first genuinely new α structural result since Phase 4F.

**Failure:** Decomposition gives an uninteresting combinatorial split (e.g., q-uniform or q-monotone) with no closed-form tie. Rules out the Hopf-Dirac mechanism; combined with D2/D3 gives the three-tier coincidence verdict.

**Deliverable:** `debug/track_dirac_s3/D4_hopf_equivariant.py` + 3–5 page memo. If positive, becomes the basis for a new section in Paper 2 (after D5 review).

**Depends on:** D1.

**Guardrails to cite:**
- Section 3 "Hopf graph morphism (Phase 4B, α-D)" — the scalar Hopf graph quotient was negative. D4 tests the Dirac analog, which is categorically different (first-order operator, nontrivial charge-q twists). The distinction must be stated in the worker prompt.
- Section 3 "S⁵ spectral geometry (Phase 4E, α-I)" — do NOT extend D4 to higher spheres; Explorer Finding 3.1 confirms S⁵/S⁷ are uninteresting for α (140, 672 miss all targets).

**References to provide:**
- Orduz, "S¹-equivariant Dirac operators on the Hopf fibration" (notes and papers 2018-2020)
- Trautman, "The Hopf fibration—seven times in physics," *Int. J. Theor. Phys.* (2002)
- Bär & Dahl, "Small eigenvalues of the Dirac operator on circle bundles"
- Paper 2 §III (Hopf base/fiber construction as currently written)

**Expected runtime:** 2–5 days. Highest-variance track, highest-payoff track. If it's dragging past a week, the Decomposer gets invoked.

---

### Track D5 — Theorem-vs-coincidence verdict

**Goal:** Synthesize D2/D3/D4 outputs into a single decisive statement about K = π(B + F − Δ). Four outcomes mapped onto four paper-level conclusions:

| Outcome | Conclusion | Writeup |
|:---|:---|:---|
| D2 + D3 + D4 all positive | Common-generator theorem candidate | Open Tier-1b proof sprint; Paper 2 §IV major update |
| Only D4 positive | Dirac-Hopf gives one structural lift; B and F remain scalar | Paper 2 §IV substantial update with Hopf-Dirac link |
| D2 + D3 positive, D4 negative | B, F both lift to Dirac spectrum; Hopf structure not the link | Paper 2 §IV update; Tier-1b "find the link" sprint |
| All negative | Three-tier coincidence formally documented | Paper 2 §IV closure statement; α remains conjectural with structural-mystery framing |

**Algebraic-first constraint:** D5 is synthesis; no new computation. If a D5 worker wants to re-run anything numerically, something has gone wrong upstream.

**Success:** Unambiguous verdict in one of the four states; Paper 2 §IV rewrite ready for PI review; explicit recommendation on whether to open Tier-1b or move directly to Tier 2 (spin-ful composed qubit encoding).

**Deliverable:** `docs/dirac_s3_verdict.md`; proposed Paper 2 §IV rewrite as a drop-in replacement (not a splinter file — full text ready for in-place edit).

**Depends on:** D2, D3, D4 complete.

---

### Track D6 — Dirac analog of the 18-symbolic-proof integrity suite (DEFERRED)

**Goal:** Extend `tests/test_fock_projection.py` with a Dirac analog under the Biedenharn-Johnson-Lippmann algebra. 4–8 new symbolic proofs.

**Status:** **DEFERRED to after D5 verdict.** Explorer Gap #1 confirms there's no off-the-shelf Dirac Fock projection theorem to verify against. D6 prematurely commits to infrastructure for a theorem we haven't established exists.

**Revisit when:** D5 verdict lands and the post-verdict plan calls for either Tier-1b proof work or Tier 2 qubit work that needs symbolic integrity gates.

---

## Sequencing

```
Day 0:         PI approves sprint plan.
Days 1–3:      Track D1 (infrastructure). Blocks D2/D3/D4.
Days 3–5:      Tracks D2, D3, D4 in parallel. D2/D3 are cheap; D4 is the long pole.
Days 5–10:     Track D4 completion (if clean algebraic structure); otherwise invoke Decomposer.
Days 6–11:     Track D5 synthesis. Paper 2 §IV draft.
Day 12:        PI review. Verdict statement. Post-verdict plan (Tier-1b or Tier 2).
```

**Realistic expectation:** this compresses to 1–2 weeks if the algebraic structures fall out cleanly. The Leader's 3–4 week estimate assumed numerical-first workflows; with algebraic-first, the timeline shortens substantially.

---

## PM prompt preamble (for the worker-dispatching session)

Every sub-agent dispatch in this sprint must include the following preamble verbatim before the track-specific task:

```
GOVERNING PHILOSOPHY (non-negotiable):
GeoVac has repeatedly found that continuum problems admit algebraic solutions
— Paper 12 Neumann expansion, Gaunt integrals, Level 2/3 Laguerre recurrence,
hypergeometric_slater.py, Paper 24 Bargmann-Segal — and the algebraic route
is consistently faster, more accurate, and more structurally informative than
numerical methods. BEFORE committing to any numerical computation taking more
than a few minutes, you MUST:
  1. Identify the transcendental content per Paper 18's exchange-constant
     taxonomy (intrinsic / calibration / embedding / flow).
  2. Check whether quantum-number recursion, group-theoretic selection rules,
     or exact Fraction arithmetic gives a closed form.
  3. If the algebraic route is unclear, stop and report back rather than
     running a long numerical job.

Numerical output in floating point is not a result. A sympy symbolic equality
IS a result. A "match to 10⁻⁸" requires an exact-arithmetic follow-up or it
doesn't count.

Reference: Paper 18 (exchange constants taxonomy) is the governing document
for where transcendentals can and cannot enter.
```

---

## Papers affected

- **Paper 2** (conjectural, α combination rule) — target of D5 §IV rewrite. Framing locked until D5 verdict.
- **Paper 23** (nuclear shell, Fock rigidity theorem) — NOT updated in Tier 1. Dirac Fock rigidity is deferred.
- **Paper 24** (Bargmann-Segal, HO rigidity) — referenced for π-free certification pattern in D1; not modified.
- **Paper 18** (exchange constants taxonomy) — governing reference throughout; not modified unless D5 verdict requires a new classification tier for Dirac-spectral invariants.

---

## Memory references (for future PM sessions reading this plan)

- `feedback_algebraic_first.md` — the governing philosophy
- `project_dirac_s3_pivot.md` — the three-tier research direction
- `feedback_tc_correction.md` — always use particle-number-projected FCI (relevant if D4 touches multi-electron physics)
- `feedback_approach.md` — PI prefers curiosity-driven framing over product/publication
