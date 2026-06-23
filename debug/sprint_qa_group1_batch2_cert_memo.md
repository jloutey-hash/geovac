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
