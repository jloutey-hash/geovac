# Sprint L2 — Sequential Sub-Sprint Plan

**Date:** 2026-05-16
**Sprint:** L2 (Lorentzian extension of GeoVac to signature (3, 1))
**Companion:** `debug/sprint_l2a_scoping_memo.md`, `debug/data/sprint_l2a_bibliography.json`
**Total wall-clock:** ~14 weeks (~3.5 months) for L2-B through L2-E + synthesis. Add 1-2 weeks buffer → 4 months.
**L3 deferred** (Lorentzian propinquity, 6-12 months original NCG-math; not required for L2 closure).

---

## Sub-sprint dispatch order

```
[1] L2-F        falsifier-cataloguing            1 wk, concurrent with L2-B   no prereqs
[2] L2-B        Krein space construction         3 wk                          no prereqs
[3] L2-C        Lorentzian Dirac operator        4 wk                          requires L2-B
[4] L2-D        Connes axiom audit at (3,1)      2 wk                          requires L2-B + L2-C
[5] L2-E        Krein-level Paper 42 redo        4 wk                          requires L2-B + L2-C + L2-D
[6] L2-G        synthesis memo + paper edits     1 wk                          requires all above
```

L3 (Lorentzian propinquity) is the named multi-month follow-up, deferred.

---

## L2-B brief — Krein space construction on S³ × ℝ at signature (3, 1)

```
CONTEXT FILES (read once):
  - papers/standalone/paper_42_modular_hamiltonian_four_witness.tex  (Riemannian Krein analog J_TT, BW-α)
  - papers/synthesis/paper_32_spectral_triple.tex §III, §IV          (current GeoVac Riemannian triple)
  - debug/lorentzian_l0_audit_memo.md                                (28-projection partition)
  - debug/sprint_l2a_scoping_memo.md §2 (BBB tables), §4.1-4.2       (Krein space construction)
  - debug/data/sprint_l2a_bibliography.json                          (reference list)
  - geovac/modular_hamiltonian.py                                    (existing TomitaConjugation, R-side reference)

TASK:
  Construct the Krein space K_{n_max} = H_GV^{n_max} ⊗ L²(R_t)_{cutoff}
  with fundamental symmetry J = γ⁰ (temporal Dirac matrix) at finite n_max ∈ {1, 2, 3},
  using the Cℓ(3, 1) gamma matrix algebra in West-coast convention.

  The temporal direction is truncated to a small uniform grid t ∈ [-T_max, T_max]
  with N_t grid points (parameter; default N_t = 21). NO compactification to S¹_β
  — this would introduce CTCs per the Sprint L2-A audit §3.8.

CONSTRAINTS:
  - West-coast Lorentzian convention throughout (gamma matrix algebra Cℓ(1, 3)
    with metric signature (+, -, -, -)).
  - Spatial gamma matrices must reduce to GeoVac Riemannian Camporesi-Higuchi
    spinor structure in the t → 0 limit (Paper 32 §III is the bit-exact target).
  - DO NOT modify any production geovac/ module that has Riemannian usage;
    new module `geovac/krein_space_construction.py` is fresh.
  - Failed approaches to avoid:
      * S³ × S¹_β temporal compactification → CTC-violating per Geroch.
      * Nieuviarts twist-morphism → odd-dim NO-GO confirmed in Paper 42 §3.2.
      * "Just complexify time" without explicit spinor-bundle Wick rotation
        → fails to give Krein-self-adjoint Dirac.

SUCCESS CRITERIA (falsifier conditions):
  - Krein axioms verified at finite n_max bit-exactly:
      J² = +I; J* J = I; indefinite inner product is Hermitian.
  - K = K⁺ ⊕ K⁻ positive/negative-definite split is explicit.
  - Riemannian-limit check: the t → 0 / N_t = 1 reduction recovers
    Paper 32 §III H_GV (the existing Riemannian Hilbert space) bit-identically.
  - Tests: ~30 new tests in tests/test_krein_space_construction.py.

OUTPUT FORMAT:
  - Memo: debug/l2_b_krein_construction_memo.md (~2500-3500 words; structure
    matching Paper 42 memos in style and detail).
  - Module: geovac/krein_space_construction.py (~250 lines).
  - Tests: tests/test_krein_space_construction.py.
  - Data: debug/data/l2_b_krein_construction.json (Krein axiom residuals at n_max=1,2,3).
  - Report Krein axiom residuals at machine-precision scale (≤ 10^{-14}) at every tested n_max.

PAPER UPDATES (recommended, apply at L2-G synthesis):
  - Paper 32 §VIII.E new subsection "Krein space construction (Sprint L2-B)".
  - Paper 34 §V.E Lorentzian transfer audit table refinement: §III.27 status
    moves from "structural correspondence" to "Krein-construction-exists" at finite cutoff.

RISKS:
  - R1: spatial gamma matrix embedding into Cℓ(3, 1) not unique. Resolution:
    fix West-coast convention upfront; document explicitly.
  - Riemannian-limit check fails: STOP and escalate. This would mean the
    Camporesi-Higuchi spatial spinor bundle is not compatible with (3, 1)
    structurally, which re-opens WH1 PROVEN at the spinor-lift level.
```

---

## L2-C brief — Lorentzian Dirac operator on S³ × ℝ

```
CONTEXT FILES (read once):
  - papers/synthesis/paper_32_spectral_triple.tex §III                (D_GV definition)
  - papers/standalone/paper_42_modular_hamiltonian_four_witness.tex §4, §6
  - debug/l2_b_krein_construction_memo.md                              (Krein space from L2-B)
  - geovac/krein_space_construction.py                                 (Krein space module from L2-B)
  - geovac/full_dirac_operator_system.py                               (existing Camporesi-Higuchi spatial D_GV)
  - debug/sprint_l2a_scoping_memo.md §4.3                             (Lorentzian Dirac construction)
  - debug/data/sprint_l2a_bibliography.json                            (van den Dungen 2016 Prop 4.1 reference)

TASK:
  Implement the van den Dungen 2016 Proposition 4.1 lift D_L = i^t D̸_g_r on
  (S³ × R, ds² - dt²) at finite n_max. The Wick-rotated Riemannian metric is
  g_r = ds²_S³ + dt². D̸_g_r is the standard Dirac on (S³ × R, g_r) Riemannian.

  For t = 1 (one time direction): D_L = i D̸_g_r, where D̸_g_r decomposes as
  D̸_g_r = γ⁰ ∂_t + D̸_S³ via the obvious spatial-temporal splitting.

  Verify Krein-self-adjointness: D_L^× = D_L, where the Krein-adjoint is
  defined via T^× = η T† η with η = γ⁰ as the fundamental symmetry.
  Verify chirality anticommutation: γ_5 D_L = -D_L γ_5.

  At finite n_max truncate the spatial direction (Camporesi-Higuchi at n_max,
  as in Paper 32 §III) and temporal direction (small grid on R_t from L2-B).

CONSTRAINTS:
  - West-coast convention (same as L2-B).
  - DO NOT compactify time to S¹.
  - The "i^t" factor in van den Dungen Prop 4.1 is structural; verify the sign.
  - Failed approaches to avoid:
      * Just naive "γ⁰ ∂_t + D̸_S³" without the Wick-rotation Riemannian-side
        intermediate step — this can fail Krein-self-adjointness depending
        on convention.

SUCCESS CRITERIA:
  - D_L^× = D_L: Krein-self-adjoint, machine precision residual ≤ 10^{-12}.
  - γ_5 D_L = -D_L γ_5: chirality anticommutation, machine precision.
  - Riemannian-limit check: setting N_t = 1 and Wick-rotating back, D_L → D_GV
    (Paper 32 §III) bit-identically.
  - Spectrum is real (Krein-self-adjoint operators have real spectrum
    on Krein-positive subspaces).
  - Tests: ~30 new tests in tests/test_lorentzian_dirac.py.

OUTPUT FORMAT:
  - Memo: debug/l2_c_lorentzian_dirac_memo.md (~2500-3500 words).
  - Module: geovac/lorentzian_dirac.py (~300 lines).
  - Tests: tests/test_lorentzian_dirac.py.
  - Data: debug/data/l2_c_lorentzian_dirac.json (Krein-adjoint residuals,
    chirality anticommutator residuals, spectrum at finite n_max=1,2,3).

PAPER UPDATES (apply at L2-G synthesis):
  - Paper 32 §VIII.E continued with Lorentzian Dirac subsection.
  - Reference to van den Dungen 2016 Prop 4.1.

RISKS:
  - R2: temporal-direction truncation is non-canonical. Test multiple boundary
    conditions (Dirichlet, periodic-but-with-care-for-CTC, hard-wall) and
    document the choice.
  - R3: continuous spectrum in R_t direction means finite-cutoff approximations
    are unavoidable. Ensure tests are formulated so they don't depend on
    fine details of temporal discretization.
```

---

## L2-D brief — Connes axiom audit at signature (3, 1)

```
CONTEXT FILES (read once):
  - papers/synthesis/paper_32_spectral_triple.tex §IV                  (existing R-side audit at KO-dim 3)
  - debug/l2_b_krein_construction_memo.md, debug/l2_c_lorentzian_dirac_memo.md
  - geovac/krein_space_construction.py, geovac/lorentzian_dirac.py
  - geovac/real_structure.py                                            (existing R-side J at KO-dim 3)
  - debug/sprint_l2a_scoping_memo.md §2 (BBB sign table), §4.4         (Connes axiom audit)
  - papers/observations/paper_34_projection_taxonomy.tex §V.E           (L0 audit table predictions)

TASK:
  Verify Paper 32 §IV Connes axiom audit at signature (3, 1) using BBB Table 1
  sign conventions. At (m, n) = (4, 6) (which corresponds to s=3, t=1 in
  West-coast Lorentzian convention), BBB gives:
    ε = +1, ε'' = -1, κ = +1, κ'' = -1.
  These translate to:
    J² = +I        (different from R-side J² = -I at KO-dim 3)
    Jγ_5 = -γ_5 J  (different sign from KO-dim 3)
    Jη = +ηJ
    ηγ_5 = -γ_5 η  (combined sign).

  Construct the Lorentzian J_L at signature (3, 1) on the Krein space K_{n_max}
  built in L2-B, with the Lorentzian D_L from L2-C, and verify each Connes
  axiom bit-exactly at finite n_max ∈ {1, 2, 3}:
    (i)   J_L² = +I
    (ii)  J_L D_L = +D_L J_L     (this is the BBB sign at (4,6) — verify)
    (iii) J_L γ_5 = -γ_5 J_L
    (iv)  γ_5 D_L = -D_L γ_5     (same as R-side)
    (v)   Order-zero: [a, J_L b J_L^{-1}] = 0 for a, b ∈ A_GV
    (vi)  Order-one: [[D_L, a], J_L b J_L^{-1}] = 0 for a, b ∈ A_GV

  PREDICTION TO VERIFY: M3 sub-mechanism (vertex-parity Hurwitz / Catalan G)
  trivializes on chirality-symmetric Dirac spectrum at (3, 1). Specifically:
  the (3, 1) chirality grading γ_5 has a different sign structure than (3, 0),
  which has the effect of making the M3 vertex-parity sum identically zero
  on the chirality-symmetric truncation of the Camporesi-Higuchi spectrum.
  This is the L0 audit prediction (debug/lorentzian_l0_audit_memo.md §4 trivially_closed_at_3_1).

CONSTRAINTS:
  - West-coast convention (consistent with L2-B/C).
  - Use BBB Table 1 sign conventions verbatim.
  - Failed approaches to avoid: assuming KO-dim 3 signs at (3, 1); they
    differ via BBB Table 1.

SUCCESS CRITERIA:
  - All six Connes axioms verified at finite n_max=1,2,3 with residuals ≤ 10^{-12}.
  - M3 trivialization on chirality-symmetric truncation:
    compute vertex-parity Hurwitz sum at (3, 1) at finite n_max and verify
    it vanishes (machine precision); compare with the Riemannian-side
    non-vanishing value to show clean structural distinction.
  - Tests: ~25 new tests in tests/test_connes_axiom_audit_31.py.

OUTPUT FORMAT:
  - Memo: debug/l2_d_connes_axiom_audit_31_memo.md (~2500 words).
  - Module: geovac/connes_axiom_audit_31.py (~200 lines).
  - Tests: tests/test_connes_axiom_audit_31.py.
  - Data: debug/data/l2_d_connes_axiom_audit_31.json (residual table per axiom per n_max).

PAPER UPDATES:
  - Paper 32 §IV: extend the existing Connes axiom audit table from KO-dim 3
    only to KO-dim 3 + (m, n) = (4, 6) Lorentzian.
  - Paper 32 §VIII.E.D new subsection on the Connes axiom audit at (3, 1).
  - Paper 32 §VIII.E.D.M3-prediction: cross-reference to L0 audit and verify.
  - Paper 34 §V.E Lorentzian transfer audit: update M3 row from "predicted"
    to "verified" (or note negative if prediction fails).
  - CLAUDE.md §1.7 WH1 entry: brief note on Riemannian J² = -I vs Lorentzian
    J² = +I as the two natural cases.

RISKS:
  - R4: sign convention bookkeeping is error-prone. Cross-verify against
    BBB Table 1 + 2 + 3 (relating (p,q), (s,t), (m,n)) and the West-coast
    physics convention.
  - If M3 trivialization prediction FAILS: this would be a clean structural
    negative, refining the master Mellin engine domain partition (Sprint MR-A/B/C)
    in unexpected ways. Document carefully and update Paper 34 §V.E and
    Paper 32 §VIII case-exhaustion theorem prediction.
```

---

## L2-E brief — Operator-system-level Wick-rotation literal identification at Krein level

```
CONTEXT FILES (read once):
  - papers/standalone/paper_42_modular_hamiltonian_four_witness.tex   (Riemannian unified-strong proof; THE reference)
  - papers/synthesis/paper_32_spectral_triple.tex §VIII rem:bisognano_wichmann_reading
  - debug/l2_b_krein_construction_memo.md, debug/l2_c_lorentzian_dirac_memo.md, debug/l2_d_connes_axiom_audit_31_memo.md
  - geovac/krein_space_construction.py, geovac/lorentzian_dirac.py, geovac/connes_axiom_audit_31.py
  - geovac/modular_hamiltonian.py                                     (Riemannian TomitaConjugation, BW-α, BW-γ)
  - debug/sprint_l2a_scoping_memo.md §4.7, §5.4

TASK:
  Re-derive the Paper 42 unified-strong four-witness theorem (Hartle-Hawking +
  Sewell + Bisognano-Wichmann + Unruh) at the Krein level on a hemispheric
  wedge of S³ × R (Lorentzian).

  Riemannian side gives two constructions:
    BW-α (geometric):     K_α = J_polar (integer spectrum, σ_{2π}^α(O) = O bit-exact)
    BW-γ (Tomita-Takesaki): K_TT = -log Δ via polar decomposition on GNS HS-space
                            (also σ_{2π}^TT(O) = O bit-exact at finite n_max)
  These are conjugate at operator-action level: σ_t^TT(O) = σ_{-t}^α(O).

  Krein side (TASK):
    1. Define the hemispheric wedge W_L on S³ × R as the analog of the
       hemispheric wedge P_W on S³ used in Paper 42 §4.
    2. Construct the wedge KMS state ρ_W on the Krein algebra at β = 2π.
    3. Build the Lorentzian modular Hamiltonian K_L_α (geometric) and
       K_L_TT (Tomita-Takesaki, via polar decomposition on a Krein GNS space).
    4. Verify σ_{2π}^L(O) = O bit-exact at finite n_max ∈ {1, 2, 3}.
    5. Verify operator-action conjugacy σ_t^L_TT(O) = σ_{-t}^L_α(O) bit-exact.

  This mirrors Paper 42 Riemannian-side structure. The Lorentzian distinction
  is that the modular flow generators K_L_α, K_L_TT live on a Krein-space
  algebra and the polar decomposition for K_L_TT must respect the indefinite
  inner product.

  ALSO: verify the load-bearing scope finding of Paper 42 §7.2 — that
  H_local = K_α^W / β is NOT the Camporesi-Higuchi Dirac D_W — in the
  Lorentzian case. Does the same finding hold? Does it change? This is
  one of the open questions O3 of Paper 42; Sprint L2-E should land a
  verdict.

CONSTRAINTS:
  - Mirror Paper 42 structure throughout — sections, naming, table layout.
  - Use the West-coast convention (consistent with L2-B/C/D).
  - The hemispheric wedge construction must reduce to the Paper 42
    Riemannian P_W in the t = 0 / N_t = 1 reduction limit.
  - Bit-exact verification at finite n_max is the falsifier standard
    (matching Paper 42's bit-exact Riemannian results).

SUCCESS CRITERIA:
  - σ_{2π}^L_α(O) = O at finite n_max=1,2,3, residual ≤ machine precision.
  - σ_{2π}^L_TT(O) = O at finite n_max=1,2,3, residual ≤ machine precision.
  - Conjugacy σ_t^L_TT(O) = σ_{-t}^L_α(O) at general t (e.g., t = 1, π, 2π),
    residual ≤ machine precision.
  - Riemannian-limit recovery: in the t → 0 / N_t = 1 reduction limit, the
    Lorentzian K_L_α, K_L_TT reduce to Paper 42's Riemannian K_α, K_TT
    bit-identically.
  - Tests: ~50 new tests in tests/test_modular_hamiltonian_lorentzian.py.

OUTPUT FORMAT:
  - Memo: debug/l2_e_modular_hamiltonian_lorentzian_memo.md (~3500 words,
    mirroring Paper 42 §5-§8 structure).
  - Module: geovac/modular_hamiltonian_lorentzian.py (~400 lines).
  - Tests: tests/test_modular_hamiltonian_lorentzian.py.
  - Data: debug/data/l2_e_modular_hamiltonian_lorentzian.json.

PAPER UPDATES:
  - Paper 42 §10 open question O1 (Lorentzian extension to (3,1)): CLOSE
    at finite cutoff level. Update Paper 42 with new §11 "Lorentzian
    closure at finite cutoff" if PI authorizes (3-5 pages).
  - Paper 32 §VIII.E continued with Krein-level unified-strong theorem.
  - Paper 32 §VIII.D Lorentzian-side cross-manifold W2b closure:
    update from "frontier-of-field" to "Lorentzian-side closed at finite cutoff".
  - Paper 34 §V.E Lorentzian transfer audit §III.27 entry:
    "structural correspondence" → "literal identification at Krein level (finite cutoff)".

RISKS:
  - R5: H_local choice may differ at (3, 1). Document carefully if so;
    refine Paper 42 §7.2 load-bearing scope finding.
  - R6: GNS Krein space construction is non-trivial — Krein-positive
    cone vs Krein-negative cone may give two GNS spaces, only one
    of which is suitable for the polar decomposition. Resolution: follow
    van den Dungen 2016 §2 (Krein-positive completion) and document.
```

---

## L2-F brief — Falsifier catalogue (parallel, lightweight)

```
CONTEXT FILES:
  - debug/sprint_l2a_scoping_memo.md §5.5 (falsifier list)
  - debug/lorentzian_l0_audit_memo.md
  - papers/standalone/paper_42_modular_hamiltonian_four_witness.tex (Riemannian falsifier model)

TASK:
  Codify the falsifiers from L2-B/C/D/E into a single tracking document.
  Each falsifier:
    - statement (what would falsify the closure)
    - test (specific numerical or symbolic check)
    - threshold (what value of the test triggers "falsified")
    - dependency (which sprint actually runs the test)
    - escalation (what to do if falsified)

  Output: a single memo + table cross-referencing each falsifier to the
  sprint that tests it.

OUTPUT FORMAT:
  - Memo: debug/sprint_l2_falsifiers.md (~1500 words).
  - Falsifier list (Markdown table) + per-falsifier prose paragraph.

ESTIMATED WALL-CLOCK: 1 week, concurrent with L2-B.
```

---

## L2-G brief — Synthesis memo + paper edits

```
CONTEXT FILES:
  - All L2-B/C/D/E memos and data.
  - debug/sprint_l2a_scoping_memo.md
  - debug/sprint_l2_falsifiers.md (L2-F output)
  - papers/standalone/paper_42_modular_hamiltonian_four_witness.tex
  - papers/synthesis/paper_32_spectral_triple.tex
  - papers/observations/paper_34_projection_taxonomy.tex
  - papers/core/paper_31_universal_coulomb_partition.tex

TASK:
  Write the Sprint L2 synthesis memo summarizing all closures and apply
  paper-edit recommendations from L2-B/C/D/E. Specifically:

  1. Apply Paper 32 §VIII.E paper edits (Krein space; Lorentzian Dirac;
     Connes axiom audit at (3,1); operator-system-level Wick-rotation
     literal identification).
  2. Apply Paper 32 §IV extension of Connes axiom audit table.
  3. Apply Paper 34 §V.E updates (28-row partition table refined).
  4. Apply Paper 42 §11 (if PI authorizes; alternatively §10 O1 update).
  5. Apply Paper 31 §8 Sig/Op partition refinement.
  6. CLAUDE.md §1.7 WH1 entry update; CLAUDE.md §2 sprint bullet.
  7. Memory file `memory_sprint_l2.md` for PM session restoration.

  Also: outline the standalone Lorentzian-extension paper (which would be
  drafted after the L3 question is closed, but the outline is documentation).

DELIVERABLE: synthesis memo, paper edits, memory file, paper outline.

ESTIMATED WALL-CLOCK: 1 week.
```

---

## Cumulative timing

| Week | Activity |
|:----:|:---------|
| 1 | L2-B begins; L2-F concurrent |
| 2 | L2-B continues; L2-F closes |
| 3 | L2-B verification + memo |
| 4 | L2-C begins |
| 5 | L2-C continues |
| 6 | L2-C continues |
| 7 | L2-C verification + memo |
| 8 | L2-D |
| 9 | L2-D verification + memo |
| 10 | L2-E begins |
| 11 | L2-E continues |
| 12 | L2-E continues |
| 13 | L2-E verification + memo |
| 14 | L2-G synthesis |

Total: **14 weeks** (~3.5 months) without buffer; **16 weeks** (~4 months) with 2 weeks buffer for unforeseen complications. Sprint L3 (Lorentzian propinquity) is deferred and would add 6-12 months.

End of plan.
