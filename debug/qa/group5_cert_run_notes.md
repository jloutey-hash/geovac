# /qa group5 — certifying FULL run notes (2026-07-03/04)

**Shape:** FULL certifying pass (the delta-2 remediation-diff verification folded in
per PI direction). Worktree `../geovac-qa-cert-group5` @ v4.65.1 (6f9228b); seeds
committed under an innocuous subject (the delta-2 lesson). Verdict below.

## Calibration

- **25 seeds** (16 code / 4 claims / 4 citation / 1 synthesis; 2 per Sonnet-tier
  agent, 1 per Opus chunk in non-first papers), catalog classes S1–S9.
  **Sensitivity 25/25** — 23 caught by original assignees; 2 initial misses
  (S-C15 zeta float-under-SYMBOLIC-PROOF by Sonnet code-P51; S-CL2 census-cost
  inversion by Opus claims-P28/P33) recovered by fresh re-dispatches that caught
  them (run-6 protocol). Multiple cross-catches (P30 tier seed caught by its code
  reviewer via git archaeology; citation seeds cross-caught by claims).
- **Specificity 8/8** controls clean (incl. the delta-2 remediation texts: the CP²
  38% floor pipeline independently re-derived; Drake–Swainson Table I re-fetched
  from the primary PDF, both values bit-for-bit).
- One agent died at the spend limit (code-P28) and was re-dispatched fresh (2/2 on
  its seeds).

## Panel: 15 original + 2 recovery + rescan + critic + gap-closure (≈19 agents,
## ≈3.9M subagent tokens)

Deterministic first (all 7 scripts + C10 compiles 9/9): green at HEAD, re-green
after every remediation wave.

## Genuine MATERIAL findings — ALL remediated in-run (2 loop-until-dry cycles)

Papers (second-locus zombies of already-corrected readings + tier calibration):
P51 Q2 carried the Wald-rejected scalar-vs-Dirac factor-of-2 reading; the
"Connection to G8" clause carried the G4-5d-rejected φ(2)/φ(1) S_BH prefactor;
the Möbius "continuum mechanism remains open" clause predated B4; the
spinor-cone attribution sentence had Cheeger on both sides (magnitude→scalar
lineage / sign→anti-periodic extraction / spinor→Fursaev–Miele now);
thm:scalar_ak SYMBOLIC PROOF tag exceeded its backing (→ INTERNAL THEOREM,
classical θ₃/Poisson); the "by Hermiticity" Furry mechanism harmonized to
m_l-conservation at 5 P28/P33 loci (matches rem:mechanism + the symbolic test's
Im M ≡ 0 ∧ M ≢ 0 structure); P2 EM boundary factor-of-2 (g_D(3)/2 = 20 = Δ⁻¹/2,
disclosed) + nine→twelve mechanisms; P41 "equivalence…established" (seed) region
also had 44-vs-60 dual-vertex (monopole sites per the driver), 3σ-vs-2σ, and
rank-scope wording; P36 7→68 ppm, muonic Antognini column 202.53→202.16,
fifteen→twenty-eight projections; P30 "~20%"→2–35% (recomputed 1.8/34/23/5.3%),
39/39→49/49, Theorem→Proposition; synthesis carried the withdrawn "L₁ = the
kinetic term" reading (second-locus miss of the v4.65.0 P30 correction) — fixed +
C16 entry.

Coverage gaps closed with new pins (~25 new tests, all green): replica S_BH
algebra (tests/test_paper51_replica_sbh.py — the Wald-remark premise was
untested); P25 golden-ratio L₁ spectrum (sympy-exact, was memo-cited only); P30
full-graph U(1) reduction + Table-1 n_max=4 row + MC-vs-LO 2–35% band; P33 census
rules 3/8 at n_max=3; P2 Hermiticity-premise roots + mp.dps fixture; P41 W6/W7
numeric-field re-derivations + §6.5/§6.7 v5 JSON pins (σ_comb −0.0803±0.0017
z=−48.3; c_σ_ens 1.196 vs c_μ_comb 1.242); P28 Sommerfeld D₄–D₆ value pins
(D₄/D₅ independently cross-verified by direct Euler-sum partial summation) +
three-loop convergence-table pins (printed rows were stale → recomputed
1.062/3.572/5.908/7.872/14.798).

Production: the §3-documented-broken EM-tail estimator
(geovac/qed_three_loop.py, log-modulated p=1 regime) now carries a docstring
warning + RuntimeWarning + a documented-failure pin replacing its vacuous
tail>0 test (measured: 67× under at n_max=20, 3.7× over at n_max=50).

Durability: P2 p-value driver → tests/paper2_support/ (archived-measured; full
1.92e9 sweep not regression-gated — reduced-slice pin is a named follow-up);
the three pruned muonic drivers resurrected from 56d29dd^ →
tests/paper36_precision_support/ (UNMODIFIED, not re-validated — named
follow-up); papers repointed; hodge1 labeling-convention note added at P28
eq:hodge1_spectrum (the PI-named follow-up stays deferred, now disclosed inline).

Registry: 3 new C16 entries (su2-kinetic-equals-L1, wald-factor2-cone-coefficient,
sbh-phi2-prefactor), pattern-self-tested to CATCH the retired phrasings.

## Convergence

Cycle 1 = panel findings → remediation (~50 loci). Cycle 2 = adversarial rescan of
the remediation diff (3 SMALL: comment sign, mis-pointed test cite — fixed; the
"645" count dismissed by direct collection) + completeness-critic → the one
uncovered load-bearing region (P28 MZV/S_min middle third) reviewed by a focused
gap agent (its S_min/w10 core verdict: "genuinely and unusually well-proven") →
its 3 MATERIALs fixed. Final state: gates 7/7 + C10 9/9 + full backing suite
660+ green; no unremediated MATERIAL. DRY.

## Honest ceiling (named, per the coverage-honesty rule)

Archived-measured (preserved, not regression-gated): the P2 p-value 5.2e-9
number; the P36 muonic-H/muonium chains (drivers resurrected, unvalidated); P51's
five L6 production-substrate numbers (dead debug drivers — standing C14-advisory
debt); Sommerfeld D₄–D₆ Euler-sum derivations (value pins only); the LS-2
velocity form behind the "3.3×" ratio. Known open: F₂-evenness-in-t untested
(matrix); hodge1/P28 eigenvalue labeling (PI-named, now disclosed inline); the 15
slow TKappa tests not re-audited (pre-existing, green at v4.64.0). The gap-closure
and recovery agents ran unseeded (their findings are self-evidencing; their
"clean" rests on the panel-level 25/25 calibration). Seed-commit subject was
innocuous this run; two code agents still identified the seed pattern from
criteria.md and reported correctly regardless.

## Verdict: PASS → group5 CERTIFIED (5th certified branch)

Per the group4-8th-cert standard: fully calibrated panel + all dimensions
exercised + in-run loop-until-dry converged to zero unremediated MATERIAL +
deterministic surface green. No group5 §C8 headline was ever wrong across the
arc (the defects were second-locus zombies, tier calibration, stale secondary
numbers, and coverage gaps — now pinned).
