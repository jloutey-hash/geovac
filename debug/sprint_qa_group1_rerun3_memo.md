# Sprint memo — `/qa group1` certifying re-run Batch 3 (rr3)

**Date:** 2026-06-21
**Version:** v4.33.0
**Target:** group1 Batch 3 — Papers 29, 39, 40, 50, 52 + group1 synthesis
**Verdict:** FAIL (calibrated 4/4 sensitivity, 5/5 specificity) → remediated
**Seed key:** `debug/qa/group1_rr3_seed_key.json`

Completes the third and final batch of the group1 certifying re-runs (rr1/rr1b on Batch 1, rr2 on Batch 2, rr3 here). Same combined method as rr2: a full calibrated panel with an **exhaustive per-bibitem citation verify folded in**, so one run both calibrates the verdict and drains the batch's C4 tail.

## 1. Panel

15 agents, all path-pinned to the worktree `../geovac-qa-seed-group1-rr3` (branch `qa-seed-group1-rr3`, base HEAD e04b3e1) and forbidden from reading the real corpus:
- **claims ×5** (opus): papers 29, 39, 40, 50, 52
- **citation ×5** (opus, EXHAUSTIVE — every bibitem): papers 29, 39, 40, 50, 52
- **code ×4**: papers 29, 39, 40, 50 (two — 29 & 39 — died on a transient 529 Overloaded with 0 tokens on the first dispatch; re-dispatched, both returned)
- **synthesis ×1** (opus): group1 synthesis

## 2. Calibration

**4 planted seeds, all caught (sensitivity 4/4):**
| seed | site | flip | catcher |
|---|---|---|---|
| S-claims-C8 | p29 abstract l.43 | finite-size → "asymptotic statement, establishing the graph RH" | claims-29 ✓ |
| S-citation-C4 | p50 KPS bibitem | arXiv:1105.4598 → 1105.5498 | citation-50 ✓ |
| S-code-C2 | test_paper40 l.44 | `<0.05` → `<1e9` (vacuous) | code-40 ✓ |
| S-synthesis-C9 | synthesis l.136-138 | "reproduces … bit-exactly" → "derives … from a CFT calculation" | synthesis ✓ |

**5 controls, none false-flagged (specificity 5/5):** p29 finite-size body (correct); p39 √2-triangle (SOUND); p40 rank-1 4/π + named-gap honest scope; p50 KPS bit-exact arithmetic match (not derivation); p52 Category III at state-space-GH level. The v4.29.0 rc3 fixes held.

Panel CALIBRATED → the FAIL verdict is trustworthy.

## 3. Genuine MATERIAL defects remediated (real corpus)

- **p40 C4 LARGE — vinberg1990 misattribution.** The "dominance maximises inner product" Weyl-chamber lemma (l.1176, l.1414), load-bearing for the PRV-summand inequality → C₃(G)=1 asymptotic tightness, was cited to Vinberg's *Some commutative subalgebras of a universal enveloping algebra* (shift-of-argument / Mishchenko–Fomenko; does NOT contain it). The lemma is the classical dominance-order fact. Reattributed to Humphreys §13 (added `humphreys1972`, removed `vinberg1990`). Inequality true; referent corrected. **Reconciled against primary text** (read l.1168-1188 + bibitem l.2356-2367 before acting).
- **p40 claims — abstract/theorem tier.** Abstract asserted 4/π universality flatly ("class-wide theorem") without the §L2 fit-sensitivity/named-gap caveat. Ported tiered language into abstract; added inline Tier note to `thm:universal_constant_analytical`. `\ref{sec:L2}` resolves.
- **p52 C14 LARGE — Q4 Lorentzian zombie.** "Lorentzian extension of the state-space GH convergence (45,46,47) extends Category III to non-compact carriers" presented descoped 45/46 (K⁺ annihilation degeneracy theorem) + partial 47 as established. Restated: metric-level Lorentzian extension does not currently close; only Paper 47 norm-resolvent arrow survives. + C7 ×2 (p38 "propinquity proof"→state-space GH; Q2 p40 all-groups hedged to match p52's own §body).
- **p39 code C2 — symbolic-bound bug.** `c3_full_pythagorean_bound_symbolic` scanned only 3 corners → returned 7/6 at (10,4), under-reporting the true interior max √(121/87) at (8,5) and contradicting the paper's "Asymmetric supremum" Remark. Fixed to exact full-grid scan (matches numeric `c3_full_triangle_bound`); added `test_symbolic_asymmetric_interior_attainment`; corrected stale withdrawn-Pythagorean docstring. Convergence theorem + ≤√2 cap unaffected (they ride the already-correct numeric bound). + C7 NIT: two "joint propinquity bound" prose labels → "state-space-GH bound".
- **p50 C7 + citation.** "Propinquity-rate / propinquity-bound" → state-space GH (inherited from P38). Squashed-S³ partition-function claim mis-cited JKPS (round-S³ only) → repointed to Hama–Hosomichi–Lee (`hama_hosomichi_lee2011`, arXiv:1102.4716). "Weyl Dirac" / "Free Dirac (Weyl, CH)" → corrected (no Weyl projection in odd dim).
- **p29 citation + title.** "Kotani–Sunada's Theorem 1.3 supplies the Ξ-completion (product over degree strata)" — specific theorem number attached to a statement the source does not appear to contain → restated as the natural generalisation / our proposal (KS kept as the general irregular-Ihara–Bass reference). bander1966 title "Group Theory of the Hydrogen Atom" → "and the Hydrogen Atom (I)".

## 4. Coverage gaps logged (verified-true, weak/no test — raise to PI)

- **p29 Observation 2 block-decomposition** (det M_sym/antisym = P22/P12, + 22+24 Dirac-B): NO test, only debug memo. The P12·P22 *factorization* is tested (`test_ihara_zeta.py` slow); the block-*origin* identity is not.
- **p29 Alon–Boppana bound-crossing**: JSON smoke test that doesn't recompute spectra (already logged rc3).
- **p50 wedge-KMS entropy + CHM identification (Thm 4.2)**: load-bearing Sec-3/Thm-4.2 headlines with no `tests/` backing, only debug memos. (The KPS F-value match itself IS tested, BACKED-SOUND.)

Disposition: log here + matrix + CHANGELOG; not remediated this round (consistent with rc3's disposition of similar verified-true coverage gaps). Candidate for a dedicated test-backfill sprint.

## 5. Reconcile catches (over-flags, not acted)

- p52 "propinquity-style metric" (l.302) = generic family-language for the abstract metric class; GeoVac instance correctly pinned to P38 state-space GH.
- p50 KPS match correctly stated as reproduction-not-derivation; `test_paper50_f_theorem.py` BACKED-SOUND (factor-4 S⁵ scalar-multiplicity fix verified, pre-fix value provably excluded).
- p39 k-fold √k extension correctly stamped "sketch; only k=2 proved".

## 6. Honest scope

- **Theorem grade (verified this round):** p40 SU(2) rank-1 4/π convergence (code-40 reran: value 0.0078 at n=512, genuinely converges; the `<0.05` real-corpus tolerance is meaningful); p39 √2-triangle tightness + diagonal closed form (PANEL-VERIFIED on real CH matrices); p50 KPS S³/S⁵ bit-exact arithmetic match (~1e-82); p29 Ramanujan verdicts on the Hashimoto spectrum + Kramers/chirality obstruction + S⁵ N_max=3 closed-form factorization.
- **Structural sketch / named gap:** p40 rank-uniform 4/π analytical proof (named gap; robust content = Reading A over Reading B at rank ≥2); p39 k-fold √k extension (k=2 proved only).
- **Numerical observation:** the rank-2 A-over-B discrimination (G2 cleanest, c_can 1.6951 closer to A than B).
- **Named open follow-ons:** the four coverage gaps in §4 (test backfill); a CLEAN re-run of each of the three batches to convert FAIL→remediated into a certified PASS.

## 7. Verification

- All 5 papers compile errors=0/undef=0 (new `\ref{sec:L2}`, `humphreys1972`, `hama_hosomichi_lee2011` resolve).
- `geovac/gh_convergence_tensor.py` full-grid fix verified: (10,4)→√(121/87), (3,3)→√(6/5) regression; symbolic == numeric.
- New `test_symbolic_asymmetric_interior_attainment` passes; `test_gh_convergence_tensor.py` non-slow green.
- Deterministic gates C5/C11/C13/C14/C15 group1 PASS.
- Worktree removed; seed-leak scan of real corpus CLEAN (no p29/p50/p40-test/synthesis seed present).

## 8. Status

All three group1 certifying-re-run batches (1/2/3) now FAIL→remediated with their C4 tails drained. **Next:** one CLEAN re-run per batch (convergence) to earn a certified PASS. No genuine defect is expected to survive a clean run, but the gate requires the calibrated clean pass to certify.
