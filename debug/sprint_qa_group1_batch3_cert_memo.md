# /qa group1 Batch 3 — first FULL cert (P29/39/40/50/52 + synth) FAIL→remediated (2026-06-23, v4.44.0)

PI invoked `/qa group1 batch3` — the first all-5-dimension cert of Batch 3 = Papers **29** (Hopf graphs
Ramanujan / Ihara), **39** (tensor k-fold state-space GH), **40** (universal 4/π rate), **50** (CFT₃
F-theorem), **52** (type-III correspondence) + the group1 synthesis. Prior Batch-3 runs were code-only
validation (v4.43.5, code-side zombies) + the v4.34.0 coverage backfill.

Base HEAD 1477df1 (v4.43.7). Worktree `../geovac-qa-seed-g1b3`. 6 seeds (one per gating dimension, covering
all 5), controls M1–M5. `debug/qa/group1_batch3_seed_key.json`. 16 LLM reviewers (code/claims/citation ×5 +
synthesis ×1), opus, path-pinned, blind; + 5 deterministic gates.

## Calibration: PASS — 6/6 sensitivity, clean specificity
- S-code-50 (F-theorem bit-exact `<10^(-40)`→`<10^(40)` flipped exponent) → code-50 ✓ (MATERIAL; transitively backstopped)
- S-code-29 (crossing discrimination `dev>0`→`dev>-1e9`) → code-29 ✓ (NIT; l.132 sibling backstops)
- S-claims-C7 (p39 l.88 "state-space GH"→"Latrémolière propinquity") → claims-39 ✓ (+ code-39 cross)
- S-claims-C8 (p40 l.262 "proof sketch/mechanism"→"complete symbolic proof for all ranks") → claims-40 ✓ (+ code-40 cross)
- S-citation-29 (Marcolli-vS arXiv 1301.3480→1031.3480) → citation-29 ✓ (404, WRONG-ID)
- S-synthesis-C9 (p50 "bit-exact reproduction"→"first-principles derivation") → synthesis ✓ (LARGE)
- M1–M5 all respected. Worktree removed + leak-scan CLEAN (6 seeds confined to worktree).

## Verdict: FAIL — calibrated panel + genuine non-seed MATERIAL defects in ALL FOUR LLM dimensions

### Per-dimension scorecard
| Dimension | Exercised | Calibrated | Clean | Carrying defects |
|---|---|---|---|---|
| Code (C1/C2) | ✅ ×5 | ✅ | ❌ | code-40 LARGE (p40 general-G backing), code-39 (incomplete zombie sweep) |
| Claims (C3/C5/C6/C8) | ✅ ×5 | ✅ | ❌ | claims-29 (p29 title + Obs over-scope) |
| Citation (C4) | ✅ ×5 | ✅ | ❌ | citation-39 (aguilar venue), citation-40 (Vinberg ref + Leimbach-vS) |
| Synthesis (C9) | ✅ ×1 | ✅ | ❌ | fabricated P29 theorem cites + Obs→theorem + p40 universality |
| Deterministic (C5/C11/C13/C14/C15) | ✅ | n/a | ✅ PASS | — |

### Genuine MATERIAL defects (verified vs primary text/code) — REMEDIATED this pass
1. **Synthesis (C9, LARGE): fabricated Paper-29 theorem citations + Observation→theorem promotion + over-scope.** The synthesis cited "Paper 29, Theorems 5.1–5.2 and Corollary 5.3" (l.213/452) and wrapped the Ramanujan result in `\begin{theorem}` claiming "for every n_max∈{2,3,4,5} satisfies strictly" — but Paper 29 has **NO `\begin{theorem}`** (confirmed grep), and code-29 measured S³ **crosses** at n=5 (dev=+0.198). Fixed: demoted to a remark-tier Observation, removed the fabricated cites (→ Paper 29 Observation [Ramanujan] + §7), corrected to the finite-size statement (sub-Ramanujan at small sizes, 3/4 families cross at modest sizes).
2. **Synthesis (C9): Paper-40 4/π universality stated as fully proven** — added the inline tier caveat (rigorous at rank 1; proof-sketch mechanism for general G; numerics-pinned ranks 1–2; rank-uniform proof a named gap; per-group extractions fit-sensitive).
3. **claims-29 (C8, MATERIAL): Paper 29 Observation 1 over-scope** ("for every graph tested" contradicted by its own §7) + §5 "leave as open question" stale vs §7's resolution. Fixed: scoped Observation 1 to small cutoffs (n_max≤4 S³, N_max≤3 S⁵) + finite-size note; §5 now forward-refs §7's crossing resolution.
4. **code-39 (C1/C2): the v4.43.5 code-side-zombie sweep was INCOMPLETE in `gh_convergence_tensor.py`.** Residual withdrawn-claim zombies (module docstring "Latrémolière propinquity"/C₃≤2; `FiveLemmaStatusTensor.L3_T` withdrawn-Pythagorean "C₃<1, =1+o(1)"; `TensorPropinquityBound`/`c_lipschitz_full_pythagorean` docstring; `epsilon_cross_bound` false-Pythagorean-identity mechanism; comment residuals) **+ the sub-flavor: `test_five_lemma_status_l3_done` ASSERTED `"Pythagorean" in L3_T`** (pinning the zombie). Fixed: all → state-space GH + triangle bound C₃≥1→√2 + WITHDRAWN flags; test → assert "triangle bound"/"WITHDRAWN" + guard. Strings/docstrings only (no logic change); affected tests green.
5. **citation-40 (C4): Leimbach-vS mischaracterized** in p40 as proving "Latrémolière propinquity convergence" (they prove state-space GH — same C7 class as the seed). Fixed → state-space GH. **citation-39 (C4): aguilar2019 wrong venue** ("Banach J. Math. Anal., to appear 2019" → Banach Center Publ. 120 (2020), 9–36) + paraphrased title. Fixed.

### RAISED to PI — named re-cert preconditions (NOT fixed this pass)
- **A (LARGE) — p40 general-G universality test-backfill (the confirmed v4.43.5 Flavor-B follow-on):** code-40 ran the machinery and found the general-G C₃=1 / L3 verifier (`verify_dirac_triangle`/`run_panel`) is imported for constructors but **NEVER CALLED** → the general-rank keystone has no live assertion; the rank-2 Weyl-integration + Reading-A-vs-B tests are `@pytest.mark.slow` (default-skipped → the cert run exercises only the SU(2) rank-1 anchor); the per-group rate-constant table values (SU(3) 1.243/Sp(2) 1.087/G2 1.177/SU(4) 0.900) have **no test** and the live G2 panel gives **1.695** (≠ paper 1.177, fit-sensitive). FIX = wire the general-G verifier into a default-collected test + migrate the `debug/qa/_resurrected/` drivers to a permanent home + add a fast reduced rank-2 assertion + reconcile/caveat the table values. Dedicated sprint. The p40 paper prose is otherwise honestly hedged (rigorous-at-rank-1 + named gap); the synthesis caveat is now added.
- **B — p29 title** "The GeoVac Hopf Graphs Are Ramanujan" (unqualified) overclaims at the most-read line though the abstract/§7/conclusion qualify it (finite-size). Editorial judgment for PI.
- **C — p29 "algebraic integer" math point:** claims-29 flagged Cor. int_alg ("every Ihara zero is an algebraic integer") justified via a NON-monic Ihara–Bass polynomial (algebraic over ℚ ≠ algebraic integer); likely the reciprocal-zeros / Hashimoto eigenvalues (monic char-poly) are the algebraic integers. Needs primary-math verification.
- **D — citation NITs:** p40 "Vinberg's lemma" load-bearing proof step has no bibitem (add a Vinberg 1990 reference); terras2010 "Thm 2.3" unverifiable; cosmetic key-year labels.

### Clean dimensions/papers (controls held)
- code-50: the F-theorem backing is BACKED-SOUND (only the seed); M3/M4 respected.
- claims-50, claims-52, citation-50, citation-52: clean — the honest "bit-exact MATCH not derivation" (p50) and structural-correspondence/anti-ontological (p52) framings confirmed sound (M4 respected); code-29 confirmed the Ramanujan backing genuine (M2); code-52 confirmed p52 conceptual; the "61+ digits" is a true conservative claim.

**Status:** v4.44.0 remediates defects 1–5 (synthesis + p29 + gh_convergence_tensor zombie sweep + 2 citations). Papers 29/39/40/synthesis compile clean (errors=0); C5/C11/C13/C14/C15 PASS; topo proofs 18/18 + the changed L3_T test green. **Batch 3 still needs a full 5-dimension clean re-cert AFTER the named p40 backfill (A) + PI rulings on B/C/D.**

---

## Re-cert preconditions RESOLVED (2026-06-23, v4.45.0)

PI directed "let's get started fixing those." All four named preconditions from the FULL-cert FAIL are now addressed.

**A (LARGE) — p40 general-G universality test-backfill — DONE.**
- Migrated the two backing drivers (`dirac_triangle_extended_verify.py`, `sp2_g2_rate_constant.py`) out of the prunable `debug/qa/_resurrected/` → permanent `tests/rank2_rate_support/` (+ README). The S⁵/l2 provenance files stay (Paper-50 backing is in production `qed_two_loop.py`, not those drivers).
- **Wired the previously-uncalled `verify_dirac_triangle`/`run_panel`** into a DEFAULT test `test_general_G_dirac_triangle_C3_leq_1` (reduced panels SU(3)/Sp(2)/G2, fail_count=0, 0<sup ratio≤1) — the general-G C₃≤1 keystone now has a live default assertion. **Un-slowed the Haar check** (`test_rank2_weyl_integration_haar_normalized` default at n_quad=40). Kept big-panel + Haar-80 + G2 A-over-B as `@slow`. Removed the `skipif` (drivers permanent).
- **Honest-scope nuance surfaced & documented:** the all-σ Dirac triangle is PANEL-BOUNDED — a naive uniform sum≤4 panel FAILS (G2 (1,0)v(0,4): max_ratio 2.42) because a small-Casimir σ exceeds the bound, while the PRV/max-σ bound (what C₃ needs) still holds. The validated `@slow` test uses the original group-specific Casimir bounds (SU(3) p+q≤5, Sp(2) a+b≤3, G2 a+b≤2 — all fail_count=0). Paper-40 prose is honestly hedged already (rank-1-rigorous + asymptotic-PRV + named gap); claim_test_matrix rows 112/113 updated (113 NO-TEST → BACKED-SOUND). Default p40 tests 3 pass; slow 3 pass.

**C — p29 "algebraic integer" math point — DONE (verified, MATERIAL).** Confirmed computationally: the Ihara *zeros* (roots of the non-monic Bass det `det(I−sA+s²Q)`, leading coeff `det Q = ∏(deg−1)≠±1`) are algebraic over ℚ but NOT algebraic integers — the paper's own closed-form factor `4s²+1` has zeros ±i/2 (minimal poly `4s²+1`, non-monic). The **reciprocal** zeros = Hashimoto T-eigenvalues (monic integer char poly `det(λI−T)`) ARE algebraic integers (±2i, `s²+4`). Corrected Cor `int_alg` (retitled "reciprocal Ihara zeros"; proof via the monic Hashimoto char poly; explicit `4s²+1`/`s²+4` example) + both synthesis spots (Headline 1 + the Observation block).

**B — p29 title — DONE (PI chose "reframe to finite-size finding").** "The GeoVac Hopf Graphs Are Ramanujan" → **"Finite-Size Graph-RH for the GeoVac Hopf Graphs: Ihara Zeta, Bound-Crossing, and a Scope Boundary on Selberg-on-Hydrogen."** Corpus-wide bibitem sweep (9 citing papers across groups 1/3/4/5 + the synthesis paraphrase + p29 self-comment); **C11 PASS corpus-wide (ALL groups)**.

**D — p40 Vinberg missing bibitem — DONE.** The "Vinberg's lemma" (two usages) is the standard Weyl-dominance fact (dominant orbit rep maximises ⟨·,ρ⟩); re-attributed to the already-grounded `bourbaki_lie_8` rather than guess an unverified "Vinberg 1990" (which would risk a new C4 defect).

**Status:** Papers 29/40/synthesis compile clean (errors=0); C5/C11(all)/C13/C14/C15 PASS; topo proofs 18/18 + p40 default 3/3 + slow 3/3 green. Batch 3 is now ready for a **clean full 5-dimension re-cert** (the LARGE backfill that blocked it is closed).
