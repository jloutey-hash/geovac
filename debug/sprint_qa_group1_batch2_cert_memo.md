# Sprint memo — `/qa group1` Batch 2 FULL certifying run (2026-06-22, v4.43.0)

**Target:** Group 1 Batch 2 = Papers 42, 43, 44, 53 + group1 synthesis (the Lorentzian-arc
FOUNDATION papers). **First all-5-dimension cert** of this batch — the prior two Batch-2 runs
(rc2 v4.28.0, rerun2 v4.32.0) were claims+citation only and both noted "next clean re-run should
certify." This run shows that prediction was optimistic: the **code dimension**, exercised in full
for the first time on Batch 2, carried a FAIL.

## 1. Method
- Base HEAD 7017802 (v4.42.1). Throwaway worktree `../geovac-qa-seed-group1-batch2` (branch
  `qa-seed-group1-batch2`), removed + leak-scanned at end (real corpus confirmed clean on all 5 seed sites).
- Deterministic gates first: C5/C11/C13/C14/C15 `--gate group1` all PASS.
- 5 blind calibration seeds (one per gating dimension) + 5 must-not-flag controls (M1–M5).
  Key in `debug/qa/group1_batch2_seed_key.json`.
- 13 reviewers, opus, path-pinned to the worktree, forbidden the real corpus: code-reviewer ×4
  (42/43/44/53), claims-reviewer ×4 (42/43/44/53), citation-reviewer ×4 (42/43/44/53), claims-reviewer ×1 (synthesis C9).

## 2. Calibration — PANEL CALIBRATED
**Sensitivity 5/5** (every seed caught):

| Seed | Class | Catcher |
|---|---|---|
| S-code-C2 (test l.90 `<1e9` vacuous tolerance, R²=I) | C2 | code-42 |
| S-claims-C7 (p44 l.1282 "Latrémolière propinquity" on Paper 38) | C7 | claims-44 |
| S-claims-C14-p53 (thm:interior "all (a)-(d) on disk … metric convergence") | C14 | claims-53 **+** code-53 |
| S-citation-C4 (p53 latremoliere2025 → arXiv:1811.10843 wrong-ID) | C4 | citation-53 |
| S-synthesis-C9 (synth "finite disk converges despite its Dirichlet boundary") | C9 | synthesis |

**Specificity clean** — M1 (p42 four-witness), M2 (p43 Krein closure), M3 (p44 l.1484 own-object
propinquity), M4 (p53 disk-obstructed/plane), M5 (synth faithful) all respected.
*Caveat on M1:* see §3 — a genuine code finding overturned the M1 *protection* (the control was
authored around a paper claim that is itself false). This is not a specificity failure (no
verified-correct control was flagged); it is a **mis-authored control** — the reviewer was right.

## 3. Verdict: FAIL → remediated
The code dimension (under-exercised in rc2/rerun2) carried the verdict. Genuine (non-seed) findings:

### MATERIAL (fixed this run)
1. **p42 §5.5(II) false closure justification (code-42, LARGE).** The paper (l.1261) said using
   $H_{\mathrm{local}}=D_W$ "the bit-exact period closure would not hold." This is **false**: the
   period closure is *defined* as the conjugation $\sigma_{2\pi}(O)=U O U^{-1}$ with $U=e^{2\pi i\beta D_W}=-I$
   (D_W diagonal, half-integer spectrum), and a global scalar cancels: $(-I)O(-I)=O$. code-42 verified
   $\sigma_{2\pi}(O)=O$ with $D_W$ at residual $1.9\times10^{-14}$. **Reconcile:** the conclusion
   ($H_{\mathrm{local}}=K_\alpha^W$, not $D_W$) survives; the genuine distinction is *operator-level*
   ($e^{2\pi i\beta D_W}=-I$ spinor double-cover vs $e^{2\pi i K_\alpha^W}=+I$), not flow-closure.
   **Fix:** rewrote §5.5(II) (l.1256–1278) + the sibling echo at l.698 to the operator-level lift;
   tied it to the WH7 compact-boost / $e^{2\pi i K}=+I$ reading. *This overturned my own M1 control —
   I had encoded the paper's imprecise claim as "correct."*
2. **p42 verch2001 descriptor (citation-42, SMALL).** In-text "nuclearity-modular" attached to Verch's
   *spin-statistics* paper (nuclearity is a different Verch work). **Fix:** →"generally covariant
   spin-statistics and modular structure" (×2, l.315/1559).
3. **p44 Connes–vS proposition numbers (citation-44, SMALL).** Cited "Proposition 4.2/4.3"; web confirms
   §4 covers Toeplitz operator systems but not the exact numbers, and one source contradicted 4.2
   (suggested 4.9). **Fix:** softened both to "§4" (l.923, l.1424) — confirmed pointer, no asserted prop number.
4. **p43 deleted-debug citation (code-43, MATERIAL).** §5.2 cited a **deleted** file
   `debug/h_local_orthogonality_formal_proof.py` as the "Source" for a load-bearing-adjacent
   closed-form — doubly wrong (the proof is *in-paper*, Lemma 5.5). **Fix:** neutralized → CHANGELOG +
   the self-contained in-paper proof.

### Coverage gaps (NO-TEST — logged, queued for backfill)
Recorded in `docs/claim_test_matrix.md` (Batch-2 section), raised to PI:
- **p43 §5.2 Pythagorean** HS-orthogonality + closed form $r^2=\kappa_g^2 S(n)/4\pi^2+D(n)$ + 18-cell
  panel: in-paper proof is self-contained, but the numerical verification has no `tests/` regression.
- **p53 disk-obstruction negatives** (min g≈−0.13; interior rate $\Lambda^{+0.07}$) **+ prop:prop2**:
  only `debug/archive/` drivers + JSON. The negative-result half + the one numbered Proposition are
  unbacked in the permanent suite. → resurrect the `debug/archive/` drivers into `tests/`.
- **p53 cited test `@slow`/skipped-by-default** — fast-gate or de-slow.
- **p42 cross-witness collapse test BACKED-WEAK (code-42, SMALL)** — collapse is by-construction
  ($\beta$-independent $\rho_W$); test compares bit-identical paths. Test-hardening (instantiate each
  witness from its own $\beta$-dependent $H_{\mathrm{local}}$).

### NITs (drained / logged)
- p42 l.495/1650 "propinquity convergence" on Paper 38 → "state-space GH convergence" (FIXED — the
  rc2 sweep missed these two in p42).
- Cosmetic citation slips (uncited orphans connes1995/strohmaier2006/takesaki1970/marcolli_vs2014/
  franco_eckstein2014 in p43; cheeger1983 in p53; key-year labels; bizi title plural) — part of the
  corpus-wide cosmetic-citation debt.

## 4. Citation apparatus — grounded
All three citation reviewers (42/43/44) + 53: **no fabricated arXiv IDs, no wrong-ID resolving to a
different paper** (the canonical failure mode), beyond the planted seed (caught). che_perales = DGA 103
(2026) confirmed grounded (not flagged). The math.OA/physics apparatus (Bisognano–Wichmann, Tomita–
Takesaki, Connes–Rovelli, Connes–vS, Camporesi–Higuchi, van den Dungen, BBB, Nieuviarts, Latrémolière,
Stein–Weiss, Stempak, Colzani) all resolve correctly.

## 5. Files changed
- `papers/group1_operator_algebras/paper_42_*.tex` — §5.5(II) + l.698 closure-prose fix; verch2001
  descriptor ×2; l.495/1650 propinquity→state-space GH; `\mathbb{1}`→`I`.
- `papers/group1_operator_algebras/paper_43_*.tex` — deleted-debug citation neutralized.
- `papers/group1_operator_algebras/paper_44_*.tex` — Connes–vS prop numbers → §4 ×2.
- `docs/claim_test_matrix.md` — Batch-2 section (load-bearing claims + 3 coverage gaps).
- `docs/qa/group1.done.md` — Batch-2 run log.
- `debug/qa/group1_batch2_seed_key.json` — seed key (created during the run).

All three edited papers compile errors=0 / undef refs+cites=0 (6 p42 "undefined" = pre-existing
font-shape warnings, not refs).

## 6. Honest scope
- **Closed this run:** the prose/citation MATERIAL tail of Batch 2 (the false-closure justification,
  the deleted-debug citation, the descriptor + prop-number citation slips, the two C7 NITs).
- **Theorem-grade survivors confirmed by code review:** p42 BW-α/BW-γ period closure + flow conjugacy
  + J_TT²=+I (GNS); p43 Lorentzian four-witness + N_t=1 recovery + Connes-axiom BBB signs; p44 prop=2
  (genuinely computed) + Riemannian limit + Krein-positivity-trivial; p53 plane Bochner–Riesz convergence.
- **Open / queued:** the 3 coverage gaps (p43 Pythagorean numeric panel; p53 disk-negatives + prop2;
  p53 slow-gate) + p42 collapse-test hardening — a dedicated **test-backfill sprint** (resurrect-from-
  `debug/archive/` per [[feedback_resurrect_pruned_artifacts]]) before re-cert.
- **Not yet certified.** This is FAIL→remediated. A clean full re-run (cert) is still required; the
  coverage-gap backfill is the named precondition.

## 7. Code-dimension re-validation pass (v4.43.1, 2026-06-22)
PI-directed narrow re-run of ONLY C1/C2 (Papers 42/43/44/53), same shape as the Batch-1 rem3 pass.
Worktree `../geovac-qa-seed-g1b2-code` (base HEAD 59c6ea4 = v4.43.0), removed + leak-scanned (corpus clean).
Seed key `debug/qa/group1_batch2_code_seed_key.json`. 4 code-reviewers, opus, blind, path-pinned.

**Validation pass, NOT a cert** (1 of 5 gating dimensions ⇒ INCONCLUSIVE-for-cert).

**Calibration — sensitivity 4/4.** Fresh vacuous-tolerance seeds (distinct sites from the full run), all caught:
| Seed | Site | Catcher |
|---|---|---|
| S-code-p42 | `test_modular_hamiltonian.py:124` (PQP `<1e9`) | code-42 |
| S-code-p43 | `test_modular_hamiltonian_lorentzian.py:180` (K_L=K_α `<1e9`) | code-43 (+code-42) |
| S-code-p44 | `test_operator_system_lorentzian.py:180` (max_residual `<1e9`) | code-44 |
| S-code-p53 | `test_paper53_plane_bochner_riesz.py:96` (rate `p>-1e9`) | code-53 |

**Fix-validation (the point of the pass):**
- **p42 §5.5(II) v4.43.0 fix HOLDS** — code-42 *independently recomputed*: $e^{i2\pi D_W}=-I$ exactly,
  $\sigma_{2\pi}^{D_W}(O)=O$ closes at residual ~1e-16, "the prose MATCHES the computation." The substance
  of the v4.43.0 correction is confirmed by fresh computation.
- **p44 prop=2 genuinely computed** (k-fold product ranks; dim O¹=14 < 64=dim O²; not hardcoded) — confirmed.
- **p43 deleted-debug citation gone** — confirmed; the §5.2 in-paper Lemma 5.5 is self-contained.

**Zero new MATERIAL code defects** beyond the 4 seeds and the already-logged coverage gaps
(p43 §5.2 Pythagorean numeric panel; p53 disk-negatives + prop:prop2; p53 `@slow`-gated test). These
were re-found (expected — still open, backfill pending), consistent with the logged state.

**One NEW NIT drained:** the v4.43.0 §5.5(II) formula carried a spurious β — $e^{i2\pi\beta D_W}$ at
$\beta=2\pi$ is $e^{i4\pi^2(n+1/2)}\ne-I$; the operator the code uses is $e^{i2\pi D_W}=-I$. Corrected
(symmetric with $e^{i2\pi K_\alpha^W}=+I$); substance unaffected; p42 compiles errors=0/undef=0.

**Verdict:** the code dimension is **MATERIAL-clean beyond the known coverage-gap tail**, and the v4.43.0
fixes are confirmed holding. The coverage-gap backfill (resurrect p53 `debug/archive/` drivers into `tests/`;
add the p43 Pythagorean numeric regression; de-slow / fast-gate the p53 plane test) remains the named
precondition for a full Batch-2 re-cert.

## 8. Coverage-gap test-backfill (v4.43.2, 2026-06-22)
PI-directed (the named re-cert precondition). All 4 gaps closed; recompute-from-framework, validated
against the live machinery (not hardcode-and-assert; the v4.34.0 group3 discipline).

| Gap | Backfill | Assertion |
|---|---|---|
| p53 disk-obstruction negative | `tests/test_paper53_disk_obstruction.py::test_disk_markov_berezin_positivity_fails` | positive f=(ρ/R)² → min g < −0.10 at Λ=200 (positivity fails) + log-log slope > 0 (rate no-decay, Λ^{+0.07}) |
| p53 prop:prop2 | `..._disk_obstruction.py::test_disk_operator_system_prop2` (4 cells) | O² fills $M_N(\C)$ (prop=2), N∈{8,12,12,18}, via live `operator_system_dim` |
| p43 §5.2 Pythagorean | `tests/test_paper43_pythagorean_hs_orthogonality.py` (7 cells) | ⟨H_local,D_W⟩_HS<1e-12; chirality pairing max diff/sum<1e-14; ‖H−D‖²=‖H‖²+‖D‖² (<1e-10); Riemannian + Lorentzian |
| p53 plane @slow | `test_paper53_plane_bochner_riesz.py` de-slowed | both markers removed; runs by default (~1.5s) |
| p42 collapse BACKED-WEAK | `test_modular_hamiltonian.py` +2 tests | bit-identical ρ/Δ/K_TT across distinct-β witnesses (<1e-13) **+ negative control** (β-dependent ρ without 1/β, ≈4.7e-4 vs <1e-13 = 9-order separation) |

Paper citations repointed debug/→tests/ (p53 prop:prop2 + rem:interior_correction; p43 §5.2). claim_test_matrix
Batch-2 rows → BACKED-SOUND. C13 PASS; p43/p53 compile errors=0/undef=0; 80 passed/4 slow-skip; no production
code changed. **The code-dimension precondition for the Batch-2 re-cert is satisfied.**

## 9. RE-CERT — full 5-dimension, post-backfill (v4.43.3, 2026-06-22)
PI-directed full re-cert from HEAD 81874e1 (v4.43.2). Worktree `../geovac-qa-seed-g1b2-cert`, removed +
leak-scanned (corpus clean). Seed key `debug/qa/group1_batch2_cert_seed_key.json`. 13 reviewers (code/
claims/citation ×4 + synthesis ×1), opus, blind, path-pinned.

**Calibration — sensitivity 5/5** (fresh seed sites): S-code-p43 (`max_alpha_residual<1e9`)→code-43;
S-code-p44 (`residual>-1e9`)→code-44; S-claims-p44-C14 (OPEN→"established")→claims-44 + code-44;
S-citation-p42 (bizi 1611.07026 wrong-ID)→citation-42; S-synthesis-C9 (disk-converges Λ^{-1.30})→synthesis.
**Specificity clean** — M1 (§5.5(II) corrected), M2 (p53 backfill), M3 (p43 §5.2), M4 (p42 collapse), M5
(p44 prop=2/§4), M6 (che_perales) all respected; the v4.43.2 backfill confirmed holding.

**VERDICT: FAIL → remediated.** Genuine (non-seed) material defects, all fixed:
| # | Defect | Dim | Grade |
|---|---|---|---|
| 1 | synthesis §6.1 "$e^{i2\pi n}=1$ would FAIL with $D_{CH}$" (v4.43.0 §5.5(II) sweep missed the synthesis) | C9 | LARGE |
| 2 | p53 l.656/680 descoped Paper 45 as compact-carrier precedent → Paper 40 | C7/C14 | SMALL |
| 3 | p42 van den Dungen inline-cited, no `\bibitem` → added | C4 | SMALL |
| 4 | p42 l.462 Connes–vS "Prop. 4.2" (unconfirmable) → "§4" | C4 | precision |

NITs drained: p43 thm residual 4e-16→7e-15 (matched its own tables); p42 §5.4 propinquity→state-space GH;
**code-42 §5.5(II) coverage gap CLOSED** (`test_operator_level_period_lift_DW_vs_Kalpha` — verifies
$e^{i2\pi D_W}=-I$, $e^{i2\pi K_\alpha}=+I$, lifts differ, conjugation flow closes with $D_W$; cited inline).

Compiles clean; C5/C11/C13/C14/C15 PASS; 81 passed/4 slow-skip. **Honest cap:** still FAIL→remediated,
NOT a clean PASS — this run's genuine findings were a synthesis-sweep miss + citation/cross-ref consistency
(no code-dimension regression); a fresh full re-cert is the confirmation step toward a trustworthy PASS.
