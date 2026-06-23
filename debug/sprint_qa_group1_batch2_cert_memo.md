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

## 10. FRESH re-cert — full 5-dimension (v4.43.4, 2026-06-22)
PI-directed confirmation re-cert from HEAD 83bd369 (v4.43.3). Worktree `../geovac-qa-seed-g1b2-rc`,
removed + leak-scanned (corpus clean). Seed key `debug/qa/group1_batch2_rc_seed_key.json`. 13 reviewers,
opus, blind, path-pinned.

**Calibration — sensitivity 5/5** (fresh seed sites): S-code-p43 (`abs(hs_ip)<1e9`)→code-43; S-code-p53
(`min(mins)<1e9`)→code-53; S-claims-p43-C7 (state-space GH→Latrémolière)→claims-43; S-citation-p44
(hekkelman 13865→13685 wrong-ID)→citation-44; S-synthesis-C9 (finite cutoff→continuum limit)→synthesis.
**Specificity clean** — M1–M7 respected (the v4.43.x fixes confirmed holding).

**VERDICT: FAIL → remediated.** 2 genuine material defects, **both recurring classes in new locations**:
| Class | This run's locations | Fix |
|---|---|---|
| 1. False-closure framing ("$D_{CH}$ would not produce the closure") | **p42 abstract** (l.127-129) + **Paper 32** (l.6720-6724) | corpus-wide grep sweep → operator-level distinction; re-grep now CLEAN |
| 2. §5.2 closed-form ($r^2=\kappa_g^2 S(n)/4\pi^2+D(n)$, $1/\pi^2$, PSLQ) untested | p43 §5.2 (only ⟨H,D⟩=0/Pythagorean were backed by v4.43.2) | new `test_pythagorean_residual_closed_form` ($\|H\|^2 4\pi^2=S(n)$, $\|D\|^2=D(n)$); §5.2 prose cites it |

Class-1 history: §5.5(II) body (v4.43.0) → synthesis (v4.43.3) → p42 abstract + Paper 32 (this run). The
location-by-location fixing is why fresh adversaries kept finding the class; **the corpus-wide grep sweep
is the convergent fix.** NITs drained: connes_rovelli end-page 2918→2917; witness-pair ~38% magnitude
pinned (`>1e-6`→`0.30<r<0.45`); cross-witness named-test caveat (operator-level test is load-bearing).

Compiles clean; C5/C11/C13/C14/C15 PASS; 45 passed/1 slow-skip (p43/p44 tests) + cross-witness 8 passed.
**Honest cap:** still FAIL→remediated. But both recurring classes are now swept corpus-wide (false-closure
grep-clean; §5.2 closed-form backed), so the remaining tail should be thin — a further confirmation re-cert
is the certification step.

---

## §11. rc2 confirmation re-cert (2026-06-23, v4.43.6) — FAIL→remediated; CODE dimension caught two genuine MATERIAL defects

PI invoked `/qa group1 batch2` hoping the v4.43.4 corpus-wide sweeps left only a thin tail (a PASS).
Full 5-dimension certifying run from HEAD 6c8e2c5 (v4.43.5). Worktree `../geovac-qa-seed-g1b2-rc2`,
5 fresh seeds (one per gating dimension incl. a fresh C5 K-prohibition seed) + 5 controls M1–M5
(`debug/qa/group1_batch2_rc2_seed_key.json`). 13 reviewers (claims×4, synthesis×1, citation×4, code×4),
all opus, path-pinned to the worktree, forbidden the real corpus.

**Calibration: PASS — 5/5 sensitivity, clean specificity.**
- S-claims-p42-C5 (K "is a derived consequence with a first-principles derivation") → claims-42 ✓ (MATERIAL/LARGE C5)
- S-synthesis-C9 ("σ_2π(O)=O in the continuum limit") → synthesis ✓ (MATERIAL/LARGE)
- S-code-p42 (test l.762 `<1e9`) → code-42 ✓ (vacuous tol; SMALL, n=3 sibling backstops)
- S-code-p44 (test l.180 `<1e9`) → code-44 ✓ (NIT, L176 `assert ok` backstops)
- S-citation-p53 (latremoliere2015 vol 103→113) → citation-53 ✓ (SMALL)
- M1–M5 all respected (no control false-flagged). Worktree removed + leak-scan CLEAN.

**Verdict: FAIL** — calibrated panel + 2 verified non-seed MATERIAL defects, BOTH caught by the CODE
dimension (running the operators / running both weights), both accepted by claims+citation+synthesis.
This is the recurring arc pattern: the code dimension catches what prose review accepts.

### MATERIAL-1 (code-43, LARGE) — Paper 43 §5.2 formal proof: Input I2 false for the spatial Dirac
Theorem 5.7's proof routes through Lemma `lem:parity_HS` (A even × B odd) with B = the whole wedge Dirac
D_W^L claimed Π_W-**odd** (I2). The backing test (`test_paper43_pythagorean_hs_orthogonality.py` l.46) and
a direct check show the **spatial** wedge Dirac D_GV^W is **diagonal** (eigenvalue χ·(n+½)), i.e.
Π_W-**EVEN** ([Π_W,D_W]=0, {Π_W,D_W}≠0) — contradicting I2 AND the paper's own §5.3 BBB-axiom finding
({χ,D}=0 fails ⇒ D_GV chirality-even). The orthogonality is true via the **chirality sign-pairing** the
test actually exercises (χ=+1/−1 partners, equal H, opposite D). The temporal γ⁰ part IS Π_W-odd, so the
Lemma correctly handles only that factor.
**Fix (v4.43.6):** corrected I2 to mixed parity (I2a temporal-odd via `lem:parity_HS`; I2b spatial-even-
sign-paired via a NEW `lem:pairing_HS`); split Theorem 5.7's proof accordingly; fixed honest-scope (i) and
the "implied universality" passage. Same fix transported to the synthesis (I3 + the two-lemma proof). The
conclusion + closed form stand (test unchanged, 11/11).

### MATERIAL-2 (code-53, MATERIAL) — Paper 53 §4 positivity leg: −0.13 is the heat weight, not Cesàro
The paper claimed the **Cesàro** weight (Def 4.1, s≥2) fails positivity, min g ≈ −0.13. Direct recompute
(`scratchpad/check_cesaro.py`): the Cesàro weight at s=2,4,6 is **positivity-preserving** (min g > 0 at
every Λ); the −0.13 is the truncated **heat** weight e^{−tλ} (Gibbs), which the paper itself concedes is
trivially non-positive. The backing test computed the heat weight while calling itself Markov–Cesàro. The
genuine finite-disk obstruction is the **non-decaying approximate-identity rate** (e1 plateaus ~0.3,
weight-robust: heat slope +0.10, Cesàro +0.055 — verified `scratchpad/check_cesaro_rate.py`), which forces
the plane pivot. Positivity is NOT the obstruction.
**Fix (v4.43.6):** corrected all Paper-53 positivity statements (header comment, abstract ×2, §Results,
§Honest-scope, thm:interior header + part (b), rem:interior_correction, proof structure, rem:numerics) —
Cesàro preserves positivity on the disk; heat fails; obstruction = non-decaying rate. Rewrote the test:
`test_disk_cesaro_positivity_holds_heat_fails` (Cesàro min g>0 ∧ heat min g<−0.10) +
`test_disk_interior_rate_does_not_decay` (e1 plateaus >0.15 both weights). Updated claim_test_matrix rows
53 + 43. Tests 9/9.

### NITs (logged, non-blocking)
- citation-43 BBB Table-1 signs at (4,6): reviewer (ar5iv, self-hedged) read ε=−1/κ″=−1; but ε=+1 is
  FORCED by the code-verified J²=+I and the reviewer's own "relations are BBB-consistent." Almost certainly
  an ar5iv parse artifact; no load-bearing impact (the four boxed relations are bit-exact-verified). **No
  edit made** (a blind edit from an unreliable source would risk introducing the real error); flagged for a
  dedicated BBB-PDF check.
- recurring C7 "propinquity" wording in source comments / §1.2(a) thread heading / L5 lemma name /
  thm:plane_propinquity header (state-space-GH-neutralized); Cesàro-order s≥2(disk) vs s≥1/2(plane) is
  correct (different objects); code-44 "exactly 0.0" prose vs ≤1e-14 gate; citation key-year labels.

**Status:** v4.43.6 fixes both MATERIAL defects + compiles clean (43/53/synthesis errors=0), C5/C11/C13/C14/
C15 PASS, topo proofs 18/18 + 43/53 tests green. Both defects were in the paper-BODY proof/claim layer that
the prior Batch-2 runs (claims/citation-heavy) never exercised against the live machinery. A further
confirmation re-cert is the certification step.

---

## §12. CODE-DIMENSION-only validation (2026-06-23, v4.43.7) — v4.43.6 fixes HELD, code MATERIAL-clean

PI invoked `/qa group1 batch2 code dimension only` — a fresh adversarial pass on the two v4.43.6
code-touching fixes (the §11 honest cap: "the new corrections have not themselves been through a fresh
adversarial pass"). **Validation pass, NOT a cert** (1 of 5 gating dimensions ⇒ INCONCLUSIVE-for-cert at
the target level). Base HEAD d28a373 (v4.43.6). Worktree `../geovac-qa-seed-g1b2-code`, 4 fresh code seeds
(one per Batch-2 paper-with-tests, distinct from rc2's l.762/l.180), controls M1–M5 (the v4.43.6 fixes +
established-sound backing). `debug/qa/group1_batch2_code_dimension_seed_key.json`. 4 code-reviewers, opus,
path-pinned, blind, RAN the tests.

**Calibration: PASS — 4/4 sensitivity, clean specificity.**
- S-code-42 (l.425 cross-witness ρ residual `<1e9`) → code-42 ✓ (vacuous tol; NIT — Δ/K_TT siblings l.429/432 backstop at <1e-13)
- S-code-43 (l.230 BW-γ Tomita `<1e9`) → code-43 ✓ (dead gate; NIT — battery <1e-12 + verdict backstop)
- S-code-44 (l.222 `assert dims == dims`) → code-44 ✓ (tautology; NIT — `prop==2` l.219 + sibling dims==192/400 backstop)
- S-code-53 (l.108 Cesàro positivity `>-1e9`) → code-53 ✓ (graded SMALL/FALSE-POSITIVE — the only gate on the corrected positivity leg, least backstopped)
- **M1–M5 all respected.** code-43 affirmed the v4.43.6 mixed-parity §5.2 proof sound; code-53 confirmed by direct recompute the Cesàro-vs-heat framing correct; code-44 prop=2 genuinely computed; code-42 collapse+β-control genuine; closed-form+operator-lift sound. Zero controls false-flagged.

**Verdict (code dimension): VALIDATION PASS — calibrated + MATERIAL-clean.** Zero genuine non-seed MATERIAL
defects. Every reviewer "fix" recommendation targeted a SEED; the real corpus already has the correct tight
tolerances at all four sites (verified post-removal: `<1e-13` / `<1e-12` / `==[14,64]` / `>0.0`). **The two
v4.43.6 fixes HELD under fresh adversarial review** (the run's purpose): the §5.2 mixed-parity proof and the
p53 Cesàro-positivity re-attribution are both confirmed sound.

**Roll-up: INCONCLUSIVE-for-cert** — only the code dimension exercised; claims/citation/synthesis/deterministic
not re-run this pass (a Batch-2 cert needs all 5). This validates the v4.43.6 remediation, it does not certify
Batch 2.

### NITs logged (pre-existing, non-blocking, all backstopped or non-load-bearing)
- code-42: `verify_witness` STRONG_IDENTIFICATION verdict gate `<1e-10` (production) vs true ~1e-16, periodicity tests assert on the verdict string — backstopped bit-exactly by `test_K_geometric_integer_spectrum` (integer to 1e-14 + odd-parity exact). NIT.
- code-43: `test_six_witness_collapse_lorentzian` checks only scalar periodicity residuals, not Δ_L/K_TT directly — reviewer verified the stronger claim is TRUE (ρ/Δ_L/K_TT bit-identical across β). Optional hardening. NIT.
- code-43: §5.2 "8.9e-16 / 18-cell six-witness panel" cited from an archived driver (only the BW κ_g=1 cell has a live test); reviewer reproduced the number + it's a proven κ_g-independent corollary. Coverage NIT.
- code-43: M3 Reading-A test tautological (identical arithmetic twice) — mirrors a weak, honestly-labeled convention OBSERVATION. NIT.
- code-44: wedge Prop 9.2 numeric panel pinned only in a slow/skipped test (fast test checks bound well-formedness). Non-load-bearing connective material. NIT.

No code/paper edits this run (MATERIAL-clean). The Batch-2 cert still needs a full 5-dimension clean re-cert.
