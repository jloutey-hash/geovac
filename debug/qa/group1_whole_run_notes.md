# /qa group1 — WHOLE-GROUP run (chunked panel) — running notes (2026-06-23)

Deterministic layer (whole group, first): **ALL PASS** — C5, C10(compile 14/14 errors=0), C11, C13, C14(436 debug/ advisory), C15, C16(0 live zombies).

Seeds: 11 planted in worktree (3 code + 4 claims + 3 citation + 1 synth). Key: debug/qa/group1_whole_seed_key.json.
Controls (must-not-flag): M1 p42 −I/+I closure; M2 p43 spatial D_W Π_W-even; M3 p53 Cesàro-positive/heat-fails; M4 p29 reciprocal-zeros-alg-int; M5 p40 panel-bounded C₃; M6 p50 S⁷ DONE.

## CODE dimension (1/paper)

### Wave 1
- **code-29**: SEED CODE-29 (line132 vacuous count) → CAUGHT as N1, graded NIT (sibling per-family `dev>0` asserts backstop) ✓ correct materiality. Genuine: **M2 p29 Cor int_alg (reciprocal=alg-integer/monic vs Ihara-zeros-not) has NO test — LARGE coverage gap** (math true, 3-line backfill); M1 Dirac B-3 tests hang ~40min un-@slow (SMALL hygiene); NITs: prop-reflection partial, chirality test absent in worktree, no inline tier tags. VERIFY M2.
- **code-39**: no code seed (unseeded code agent). Verdict **C1/C2 PASS, 0 material**. Confirmed Pythagorean withdrawal COMPLETE in gh_convergence_tensor.py (~40 occ all WITHDRAWN-flagged; theorem-string guards zombie) — M-control spirit respected ✓. NITs: SVD triangle-tight test @slow (skipped default); claim_test_matrix only 1 p39 row; stale "pythagorean" alias naming.

### Wave 2
- **code-43**: SEED CODE-43 (line82 `<1e12` vacuous HS-orthog) → CAUGHT as M1 MATERIAL ✓. Else BACKED-SOUND (mixed-parity proof tight, 1/π² closed form tight). 
- **code-53**: SEED CODE-53 (line113 `<1e6` vacuous heat-fails) → CAUGHT as MATERIAL-2 ✓. **GENUINE non-seed MATERIAL-1**: plane-test docstring zombie (test_paper53_plane_bochner_riesz.py L4-5,108 attribute the −0.13 positivity-fail to Cesàro/Markov-Cesàro — the withdrawn conflation; live = Cesàro PRESERVES, heat fails). VERIFY+FIX+add C16 entry.
- **code-45**: no code seed. BACKED-SOUND keystone (K⁺ annihilation 20/20). **GENUINE non-seed MATERIAL (SMALL)**: `lorentzian_theorem_statement()` (geovac/lorentzian_propinquity_compact_temporal.py:780-812) returns RETRACTED convergence thm + C₃→1⁻ unflagged + `TestTheoremStatement` (test_lorentzian_propinquity.py:347-365) PINS the zombie (asserts keywords present). v4.43.5 flavor surviving in p45 module. VERIFY+FIX.
- **code-46**: PASS 0 material (M-control C₃=1 + withdrawn-√ refuted-and-tagged ✓). NIT S1 C₃=1 inherited-not-directly-tested.
- **code-49**: BACKED-SOUND 0 material (M-control Datta-not-Umegaki + "generically" ✓). NIT-1 96/96 backing in debug/ (Flavor-B durability).
- **code-50**: SEED CLAIMS-D (paper50:122 "derivation as genuine 3d CFTs") CROSS-CAUGHT by code as M1 ✓. M-control S⁷ Erratum + wrong-spectrum control ✓.
- **code-analytic-trio{47,48,52}**: SEED CLAIMS-C (paper47:135 "hold at the metric level") CROSS-CAUGHT as M1 LARGE ✓. Else descopes correctly tagged. NIT debug_bridge_correction.py misfiled in tests/.

CODE dim calibration: CODE-29/43/53 all caught (graded right); CLAIMS-C/D cross-caught. SENS so far 5/5 of the seeds reachable by code.
GENUINE non-seed CODE MATERIAL: (1) p53 plane-test docstring Cesàro-zombie; (2) p45 lorentzian_theorem_statement() + pin. 
COVERAGE GAPS (log): p29 int_alg no-test (LARGE, recurring); p29 Dirac-B3 40min un-@slow; Flavor-B debug/-path backing (p49 96/96, p46 panel/F2/F3, p50 S⁷ ladder, p47/48 panels).

## CLAIMS dimension (chunked) — Wave 3 dispatched: A[29,39,40] B[42,43,44] C[45,46,47] D[48,49,50,52,53]
## SYNTHESIS — Wave 3 dispatched
## CITATION dimension (chunked) — Wave 4
- **cite-A[29,39,40,42,43]**: SEED CITE-A (p42 connes_vs2021 arXiv:2007.09988) CAUGHT ✓ (+ gave correct 2004.14115). Else GROUNDED; NITs only (cosmetic key/year drift, p43 connes_rovelli off-by-one page).
- **cite-B[44,45,46,47,48]**: SEED CITE-B (p46 vdD Prop 6.3) CAUGHT ✓ (+ correct 4.1; §6 is electro-weak). Else GROUNDED; fabricated-Latrémolière refs correctly quarantined (positive control). NITs cosmetic.
- **cite-C[49,50,52,53]**: SEED CITE-C (p53 Muckenhoupt & G.Weiss) CAUGHT ✓ (correct = Stein). KPS F-theorem source verified bit-exact. Else GROUNDED; NITs cosmetic.
Citation dim: 3/3 seeds caught, clean beyond seeds (NITs only).

## SCORECARD
Sensitivity 10/11 seeds caught (CODE-29/43/53, CLAIMS-A/B/C/D, CITE-A/B/C). MISS: **SYNTH** — synthesis reviewer read past the "established bound" lead-in (reassured by the retraction 12 lines below) ⇒ synthesis dimension UNDER-CALIBRATED.
Specificity: clean — all M1–M6 controls respected across the panel; no false-positive on a known-good.

## GENUINE non-seed MATERIAL (the FAIL drivers) — ALL REMEDIATED
1–3. **√(1−1/nmax) "C3^op" zombie class** (v4.40/4.41 fixed p45/p46, never transported): p47 thm:inner (l.581) + thm:g2_metric (l.1345) + synthesis eq:p46_main (l.1358) cited the WITHDRAWN operator-norm envelope as live → all collapsed to C₃=1 (Paper 38 L3); p47 g2_metric γ=2/(nmax+1) mass-max→O(log/n) Berezin rate. FIXED.
4. **p53 plane-test docstring** attributed −0.13 positivity-fail to Cesàro → re-attributed to HEAT weight (Cesàro preserves; obstruction=rate). FIXED.
5. **p45 lorentzian_theorem_statement()** returned the RETRACTED weak-form convergence theorem + pinning test asserted the zombie keywords → rewrote fn to the K⁺-annihilation headline + rewrote test to guard. FIXED.
+ **C16 maintenance**: added `withdrawn-c3op-envelope-sqrt` entry (the √(1−1/nmax) form C16 didn't previously catch). C16 re-run PASS (35 occ all flagged).

## VERDICT: FAIL → remediated (NOT a certification).
- Deterministic ✓ clean · Code ✓ calibrated · Claims ✓ calibrated · Citation ✓ calibrated · **Synthesis UNDER-CALIBRATED (missed SYNTH plant)**.
- 5 genuine MATERIAL defects (1 class spanning 3 sites + 2 code zombies) found & fixed; all compile errors=0; C16 green; affected tests green.
- Re-cert precondition: a clean confirmation re-run with a SHARPENED synthesis reviewer (the "live word near an intact retraction still counts" lesson — analog of the run-#4 K-sentence fix).

## COVERAGE GAPS / NITs to raise (not fixed this run)
- p29 Cor int_alg (reciprocal=alg-integer) NO test — LARGE recurring gap, math true (3-line backfill).
- p29 Dirac-B3 tests hang ~40min un-@slow — quick @slow follow-up.
- Flavor-B debug/-path load-bearing backing (p49 96/96, p46 panel, p50 S⁷ ladder, p47/48 panels) — migration backfill (cf. tests/rank2_rate_support/ precedent).
- claims NITs: stale .tex top-comments (p40 l.7, p39 l.2 propinquity), stale `pythagorean_op_norm` eq label; analytic-trio: debug_bridge_correction.py misfiled in tests/.
