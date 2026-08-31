# DELTA-4 verification run (2026-08-29) — run record

Fired by PI ("do the delta") after the post-cert-3 remediation, the ERI
correction arc, and the Paper 26 reconstruction.

## Shape

**DELTA-verification run**, not a certifying run. Verdict vocabulary is
**CLEAN-DELTA / DEFECTS** — a delta can never return PASS. A clean delta is
the precondition for firing the next FULL certifying run.

## Scope (diff since the last calibrated run = cert-3)

- **Production code:** `casimir_ci` (M_L delta + Condon–Shortley order),
  `lattice_index` + `composed_qubit` (q sign, then the c^k(d,b) order),
  `sturmian_solver`, `tc_integrals`, `nuclear/{harmonic_shell,potential_sparsity}`,
  guard comment in `sturmian_secular`; `debug/archive/misc/entanglement_first_row`
  (diagonal-only entropy → full RDM).
- **Papers:** 14, 20 (re-priced twice — exact rule, then the order fix),
  26 (RECONSTRUCTED from a preserved PDF after a truncation incident),
  24 + 27 (Minnesota magnitude caveat), 32 + 34 (duration language),
  group4 + group6 syntheses.
- **Gates:** C16 (bare-graph entry), C17 (two new families + DoD scope),
  C18 (LaTeX-tie + symbolic escapes + runtime exemption), criteria.md.
- **Tests:** ~20 files re-pinned.

Note: the working tree also carries pre-existing uncommitted drift from
earlier sessions (P23, P35, several docs). That is OUT of this delta's scope.

## Dimensions

| Dimension | Affected? | Dispatch |
|---|---|---|
| Deterministic | yes | run whole-target, all 5 branches — **all PASS** |
| Code / test-backing | yes | 1 agent (sonnet tier → **2 seeds**) |
| Paper claims / prose | yes | 2 agents (P14+P20; P26+P24+P27), Opus |
| Synthesis | yes | 1 agent (group4 synthesis), Opus |
| External citations | **no** | no new references this delta — not dispatched |

## Deterministic layer (run BEFORE dispatch)

group1 / group3 / group4 / group6 / synthesis — C10, C11, C13, C14, C16,
C17, C18, C19 all **PASS**. Papers 14 + 20 + 26 compile three-pass clean
(0 errors, 0 undefined refs).

## Seeds (7, sandbox `C:\Users\jlout\Desktop\geovac-qa-delta4` ONLY)

Deliberately drawn from classes cert-3's calibration did NOT cover (its 8
code seeds were all one class, tolerance-loosening) plus this arc's own
recurring failure modes:

| id | class | catcher |
|---|---|---|
| D1a | witness tolerance neutered (1e-12 → 1e-1 on c^k = −√6/5) | code |
| D1b | exact-closure assertion made vacuous | code |
| D2a | retired exponent restored (propagation class) | claims-A |
| D2b | linearity claim inverted vs its own equation | claims-A |
| D3a | honest-scope inversion + self-contradiction | claims-B |
| D3b | retracted claim re-asserted in new clothing | claims-B |
| D4  | half-propagated correction (d-block ordering) | synthesis |

Placement-guarded (no comment-line anchors). Key:
`debug/qa/delta4_seed_key.json`.

## Controls (5, known-good — INDEPENDENTLY established, not same-day belief)

K1 order-independence of counts (re-measured after the order fix) ·
K2 He −2.896298 / 0.2558% (measured twice, stable to 1.3e-5) ·
K3 the 4-fold ERI symmetry group (correct for complex harmonics) ·
K4 Paper 22's 6.0608% (independent exact-sympy enumeration) ·
K5 five honest-scope sentences that must NOT be flagged as hedging.

*Cert-3 lesson applied:* controls must be independently established. Cert-3
designated two same-day, un-re-verified results as controls and both turned
out WRONG, voiding that run's specificity measurement.

## Process commitments for this run

- No calibration headline until **every** dispatched dimension has reported.
- Seeds are removed from the sandbox and the sandbox deleted at the end;
  the real corpus is swept for every seed pattern before the verdict.
- Verdict is CLEAN-DELTA only if zero verified-genuine MATERIAL findings
  AND the panel caught its seeds.

*(Results appended below when the agents report.)*

---

# RESULTS — verdict: **DEFECTS**

Three of four dispatched dimensions have reported (code pending at time of
writing). The verdict does not depend on the fourth: the three that reported
returned a large, PM-verified genuine ledger.

## Calibration so far

| dimension | seeds | caught | notes |
|---|---|---|---|
| synthesis | D4 | **1/1** | named 9.23 as the retired d-block coefficient |
| claims-B (P26/P24/P27) | D3a, D3b | **2/2** | caught the honest-scope inversion AND the re-asserted ceiling |
| claims-A (P14/P20) | D2a, D2b | **1/2** | caught D2a (retired exponent in the correction note); **MISSED D2b** ("sub-linearly" contradicting its own equation two lines above) |
| code | D1a, D1b | pending | — |

**D2b is a genuine miss and it matters**: it is an internal-contradiction
defect of exactly the class the enumeration mandate targets, in a sentence
adjacent to the equation it contradicts. Claims-A found ~25 other defects,
so it is not an uncalibrated agent — but its miss is on the class the
sharpened rule-7 prompt is supposed to make unmissable.

All 5 controls held: no reviewer flagged a control as a defect.

## Genuine findings — the headline

**The correction arc was applied to the abstracts, the `sec:eri_rule`
sections and ~6 tables, and NOT to the rest of the two papers or the
synthesis.** Both papers and the synthesis currently publish retired
pair-diagonal values as live headline numbers in their conclusions,
discussion sections and roughly half their tables.

PM-verified on the REAL corpus (seeds are sandbox-only and confirmed absent):

- **P14 `sec:111`** derives the retired coefficient *analytically*: "111
  non-identity Pauli terms per s+p block", "65 nonzero ERI quartets",
  "N_Pauli/Q = 111/10 = 11.1". Exact-rule values: 279, 107, 27.90. This is
  the paper's only analytic derivation of the linearity coefficient and it
  derives the retired number.
- **P14 conclusion** still lists "Pauli term count: O(Q^{3.15})" while the
  same conclusion says Q^{3.8} — self-contradicting.
- **P14 tab:multi_center** — every row is exactly 11.10·Q + 1 (LiF 778,
  CO/N2/F2 1,111, NaCl 556), cited three lines later as evidence for
  27.90·Q. 6 live occurrences of "1,111".
- **P20 tab:scaling** — the entire table (2.24/2.18/2.22/2.21/2.20/2.19,
  mean 2.21) is retired and is cross-referenced as if it showed 3.17;
  8 live occurrences. A Q=100 extrapolation ("≈6,000×") rides on it.
- **P26** — the reconstruction left four LARGE zombies: a stale
  `I(2s,3s) = 1.324` that the paper's own restored correction note claims
  to have removed; a Paper-14-retired "51×–1712×"; a Conclusions paragraph
  still stating the pre-correction "42% → 99% step-function"; and Li
  `I_cv = 0.213` contradicting its own Table II (0.228).
- **group4 synthesis** — six MATERIAL-LARGE, all the same shape: patched at
  the loci naming the correction, unpatched at the loci depending on it.
  **Remediated in this run** (4 passes, 15 edits): atomic exponent, the
  "scaling unchanged / constant factor" discipline (6 loci), the
  pair-diagonal-as-legitimate-approximation framing, eq:linearity's label
  and parenthetical, the market-test numbers, the balanced count, the
  Paper-38 GH-rate denial (×2), and the Trenev provenance over-claim.

## Two-way verdict (upgrades found)

- The equality of the within-molecule exponent across all 8 systems is
  currently framed as an empirical fit; it is **forced** by the exact
  factorization N = c(n_max)·Q with rational c. Under-stated.
- Isostructural invariance survived a corpus-wide re-pricing **bit-exactly**
  — a correction-invariant result stated flatly.
- The correction has **two** independent arbiters (4-D quadrature from the
  definition, and agreement with the precision-physics path), stated more
  weakly than the evidence supports.

## Deterministic layer

All gates PASS on all five branches (C10, C11, C13, C14, C16, C17, C18,
C19). **This is the point of the exercise:** the deterministic layer is
green while ~25 MATERIAL claim defects are live. C17's registries gate the
specific number-families registered; they cannot gate "an analytic
derivation of the retired coefficient" or "a table whose every row encodes
the old law". Registry gates gate *numbers*, not *reasoning*.

## Owed before the next delta

1. Complete the P14/P20 re-pricing at ALL dependent loci (not just the ones
   naming the correction) — the ledger above is the worklist.
2. Fix the four P26 reconstruction zombies.
3. Re-run the delta. A delta may never return PASS; a CLEAN delta remains
   the precondition for a certifying run.

---

# CODE DIMENSION — reported (completing the record)

## Calibration: 1/2, and the miss was caused by MY error

| seed | verdict |
|---|---|
| D1a (witness tolerance 1e-12 -> 1e-1) | **CAUGHT**, but classified NIT rather than MATERIAL. The reviewer reasoned correctly that the *historical* bug zeroes the coefficient exactly, so 1e-1 still catches it — while noting a ~20% window on a claimed EXACT closed form would pass a differently-broken value. Under-classified, not missed. |
| D1b (exact-closure 1e-12 -> 1e-1) | **MISSED** — but `test_paper26_entanglement.py` was **0/7 runnable** in the sandbox (a driver dependency was not copied), so the agent could not exercise the assertion. |

**PROCESS FAILURE (mine): I deleted the sandbox while the code agent was
still running.** The agent reported the directory vanishing mid-session and
correctly declined to keep probing it, which cut short one planned empirical
check (whether the Araki–Lieb guard would have failed under the OLD
diagonal-only entropy routine — still open, item 6 in its report). Two
lessons: (i) never tear down a review sandbox until every dispatched agent
has reported; (ii) a sandbox must be dependency-complete, or tests silently
do not run and their dimension is uncalibrated.

## The finding that justified the whole run

**`potential_sparsity.compute_eri_tensor` carried the SAME Condon–Shortley
order bug** — `ck_bd_list = ac_k_map.get((b_idx, d_idx))`. I had fixed that
module's `ck_coefficient` q-sign and never looked at its assembly.

It is invisible to `test_paper22_density::test_potential_independence_lmax2`,
the only test exercising it, because that test compares the nonzero **mask**
and the **density** across five potentials — both invariant under a sign
flip, since |c^k(b,d)| == |c^k(d,b)|. The reviewer ran it with the bug
present and it passed.

**Sweeping after that find turned up five MORE instances**:
`nuclear/harmonic_shell.py:413`, `sturmian_solver.py` (three sites),
`tc_integrals.py:498`. **All six fixed; a final sweep confirms zero remaining
assembly-order instances in production.** Total for the arc: the same
one-line defect in **13 places** across two distinct expressions (the q sign
and the factor order).

Also fixed: an unguarded `np.savez` to `geovac/cache/` (fails on any fresh
checkout with no cache dir). Repairing that introduced, and then caught, a
substring-replace indentation break — `"        np.savez"` is a substring of
`"            np.savez"`. Syntax-checked and re-verified.

`test_paper22_density` 20 passed; the P14/P20/P22/P26 set 33 passed.

## What this says about the gate

The reviewer's summary is the honest one: *"'all fixed, all tested green' is
not yet true end-to-end."* Both headline fixes are genuinely correct and
correctly wired **where their tests can run** — but a mask-and-density test
cannot see a sign flip, and a sandbox missing a driver silently runs nothing.
Green is not evidence unless the test can both execute and discriminate.

---

# POST-DELTA SWEEP — a 15th site, of the OTHER expression

Re-running the consumers of the five newly-fixed modules (the background
suites had been killed by an environment reset, so the fixes were unverified)
gave 339 passed / 4 failed. Three were straightforward value moves. The
fourth was a signal:

`test_he_max_n2_angular_adds_integrals` failed as **"265 vs 265"** — the TC
angular gradient no longer adds any integrals over radial-only. Both numbers
are 265, which is precisely the **no-M_L count** for He n_max=2 (physical:
107).

**Cause: `geovac/tc_integrals.py`'s ERI assembly never imposes the Coulomb
selection rule `m_a + m_b = m_c + m_d`.** This is the SAME missing-delta
defect `casimir_ci` carried — a fifteenth site, and the first of that
expression found outside `casimir_ci`. Its factor order was corrected in this
run; its selection rule was never there.

**Deliberately NOT fixed.** Adding the rule changes the TC operator and needs
its own verification pass. This session has already demonstrated the cost of
a half-verified physics change (the order fix was applied, reverted on an
unexplained factor of 2, then re-applied once the factor turned out to be the
k_orb/Z convention). The three dependent tests are `xfail(strict=True)` with
the precise reason rather than re-pinned to numbers that would move again.

Running tally for the arc: the same conceptual defect in **15 sites** across
two expressions —
- **q-sign** (`q = mc - ma`): 7 modules, all fixed (1 was a false positive:
  `sturmian_secular`, whose q is negated at the call site).
- **factor order** (`c^k(b,d)` vs Condon-Shortley `c^k(d,b)`): 7 sites, all
  fixed.
- **missing M_L delta**: `casimir_ci` (fixed) and `tc_integrals` (**OPEN**).

## Verification state at close

- Syntax-checked: all 7 modified production modules parse.
- Consumers of the 5 newly-fixed modules: **339 passed, 0 failed, 3 xfailed**
  after the quarantine (was 4 failed).
- P22/P14/P26/P20 set: 33 passed.
- Deterministic gates: PASS on all 5 branches.
- **NOT re-run since the last six production fixes:** the full bounded suite
  (killed by the environment reset). The consumer set above is a superset of
  the affected modules, but it is not the whole suite.

---

## Post-delta addendum (2026-08-30): the TC M_L rule, fixed and split

The "deliberately NOT fixed" item above was worked, and it split cleanly in
two. `tc_integrals` assembles the transcorrelated ERI in two passes -- a
radial-gradient pass and an angular-gradient pass -- and only the first was
salvageable.

**Radial pass: FIXED, and it lands exactly.** With the Coulomb selection rule
`m_a + m_b = m_c + m_d` added, He n_max=2 gives **107 nonzeros with 0
M_L-violating entries** -- the physical ERI support, matching the corrected
`lattice_index`/`casimir_ci` count exactly. The retired 65 came from the
pair-diagonal sign error; the 265 seen mid-arc came from the missing rule.
Both are now accounted for. LiH composed radial-only re-pins 562 -> 1354
Pauli.

**Angular pass: a SIXTEENTH site, and a sharper finding than expected.** The
angular-gradient assembly contributes 26 entries beyond radial-only at He
n_max=2, and **all 26 violate L_z conservation** -- none conserve it. For a
rotationally invariant correlator u(r12) the transcorrelated Hamiltonian must
commute with total L_z, so every one of those entries is a bookkeeping
artifact. They are structured, not noise: magnitude ~5.3e-3 against a
conserving max of 1.25 (ratio ~1%), e.g. the (1,0,0)(1,0,0)|(2,1,-1)(2,1,-1)
entry carrying M_L 0 -> -2 at +5.290349e-03.

**No filter was added.** Dropping the violating entries would hide the error
and leave the survivors unverified -- the same masking failure the arc has
been correcting elsewhere. Instead the test now *asserts* the finding: 26
added, all unphysical, radial-only all clean.

This sharpens an existing dead-end rather than opening a new one. CLAUDE.md
section 3 already records "TC angular gradient for l>0 orbitals: adds 2.66x
Pauli terms for 0.01 pp accuracy; radial-only is optimal." On this basis the
angular term buys not 0.01 pp but *nothing* -- its entire contribution to the
ERI tensor is outside the physical sector.

Arc tally, final: same conceptual defect in **16 sites** across three
expressions -- q-sign (7 modules, 1 false positive), factor order (7 sites),
missing M_L delta (`casimir_ci`, `tc_integrals` radial; `tc_integrals`
angular is a distinct assembly bug, recorded not fixed).

Tests: `tests/test_tc_angular.py` 13 passed, 1 skipped, 0 xfail.

---

## Remediation pass 2 (2026-08-30): the dependent-loci ledger

Worked the P14/P20/P26 re-pricing the delta left open, plus everything the
sweep turned up on the way.

### Paper 14

**`sec:111` rewritten, not renumbered** (label -> `sec:block_pauli_count`;
zero inbound refs, and the old name was itself a zombie once 111 retired).
Every number in that subsection was *derived* from the pair-diagonal rule,
so swapping the numerals would have left a false derivation carrying true
values. The replacement is measured and closes:

- 107 ERIs = 25 direct + 82 non-direct (retired 65 = 25 + 40).
- 279 non-identity Pauli = 55 + 224 (retired 111 = 55 + 56). The **55 direct
  survive verbatim** -- the monopole argument c^0(l,m;l,m)=1 involves no sign
  convention and no m-rule -- and are exactly the all-Z strings; all 224
  others carry an X or Y. Measured, not assumed.
- A closed derivation of 107 replaces the retired 7^2+4^2: with P_k the
  k-allowed orbital pairs (|P_0|=7, |P_1|=12, |P_2|=9), the union over k of
  the M_L-constrained products |P_k x P_k| is **107 exactly**. The naive
  sum of squares is 274; it over-counts because it ignores the M_L coupling
  between the two pairs and double-counts quartets admitted by several k.
- Parity forces |dl_ac| = |dl_bd|: 59 quartets at 0, 48 at 1.
- The **pure-l-shell collapse is withdrawn**. It claimed m-conservation
  forces a=c, killing all exchange (55 Pauli, 5.5 per qubit, "explains the
  reduced transition-metal coefficient"). The exact rule constrains only the
  sum, which ordinary exchange <m1 m2|m2 m1> satisfies identically: a pure-d
  shell carries **85 nonzero ERIs, not 25**, 60 of them off-diagonal.

**Other P14 loci:** `tab:multi_center` re-priced from the shipping library
(1,953/2,790/2,790/2,790/1,395, every row exactly 27.90 x Q, switched to the
non-identity convention so no stray +1); three live `1{,}111` zombies in
prose; seven live `Q^{3.15}` zombies (the QWC comparison at L603 also
*reverses* direction -- 3.36 is now below the term exponent, not above);
and L1142, which asserted 778 / 18,383 Pauli for H2O while the table four
lines above already carried the corrected 1,954 / 99,247.

### Paper 20

`tab:scaling` re-measured by the paper's own method (two-point fit at
n_max = 1, 2). **The correction strengthens the structural claim while
cutting the headline extrapolation.** All six molecules now give
alpha = **2.816 identically** -- not "nearly constant" with +-0.02 spread,
but equal to the digits shown, and forced rather than empirical: the Pauli
count is exactly linear in Q at each fixed basis (1.5 Q at n_max=1,
27.90 Q at n_max=2) and Q doubles-and-a-half by the same factor 5 for every
molecule, so alpha = 1 + log_5(18.6) = 2.8163... independent of the
molecule. The Q=100 extrapolated Gaussian ratio falls from ~6,000x to
**~370x**.

`tab:multicenter` (a second multi-center table, in Paper 20) was **caught by
the new C17 entry, not by the manual sweep** -- re-priced with measured
Pauli and greedy QWC (1,953/2,790/1,395; QWC 21 -> 64).

### Paper 26

Four zombies, all measured at n_max=2 (the basis the section states it
uses) rather than patched:

- The **Lithium bullet was wrong in every particular.** It quoted
  I(2s,3s)=1.324 -- an n_max=3 value for an orbital that does not exist at
  n_max=2, which the section's *own* correction note flags before leaving
  the bullet unchanged -- and claimed "the 1s orbital drops to 9% of the
  total mutual information". Measured at n_max=2: the dominant pair is
  I(1s,2s)=0.218 and the **1s share is 89.8%**, with 2s at 91.2%. Lithium
  has no hub shift; 1s and 2s are co-hubs of the same pair. Rewritten.
- He I(1s,2s) 0.877 -> **0.748** (1s share 99.0%, sole hub, unchanged
  qualitatively). Be: the first genuine shift -- 1s collapses to **0.5%**,
  2s becomes the hub at **65.1%** (paper said 43%), strongest single pair
  already p-p at 0.202.
- The conclusion's hub-migration sentence now carries the Li caveat: the
  migration begins at beryllium.
- Li I_cv 0.213 -> **0.228**. Decidable without new measurement: Table II
  already said 0.228 and `test_paper26_entanglement.py` pins 0.2276; the
  prose carried the retired 2*S_core value.
- Pauli advantage 51x--1712x -> **~50x--320x** (P14's corrected canonical
  range); conclusion step-function densities 42%->99% -> **17.1%->100%**.

### SEVENTEENTH SITE (found, quantified, NOT applied)

`geovac/composed_qubit_relativistic.py::_build_spinor_eri_block` carries the
same Condon-Shortley factor-order defect corrected at seven other sites: it
binds the dict key as `(b, d)` and multiplies `X_k(a,c) * X_k(b,d)` where the
rule requires `X_k(d,b)`. Its q-sign and its M_L rule are already correct.
(Three other `mc - ma` hits in the tree were re-audited and are genuine
non-defects: `sturmian_secular` negates q at the call site, and the two
`composed_qubit` counting loops test only `|q| <= k`, which is sign-blind.)

Measured impact, LiH_rel n_max=2: ERI support identical (2676 both) and
sum|ERI| identical -- **magnitudes unchanged, only signs: 1080/2676 entries
(40.4%) flip**. Sign flips change JW cancellation, so downstream counts do
move: N_Pauli 1413 -> 1501 (LiH, BeH), 942 -> 998 (CaH); lambda 40.594 ->
39.533, 143.96 -> 142.52, 18.679 -> 17.914. The n_max=1 rows are
bit-identical, which is the correct control: all-s spinors make the two
orderings coincide, so the fix must be a no-op there, and is.

**Not applied.** Both orderings satisfy the required 4-fold permutational
symmetry -- particle exchange and hermiticity hold either way, as the
M_L-constrained phase algebra predicts and a direct test confirmed (0
violations both ways) -- so **no internal discriminator exists**. The scalar
analogue was adopted only after verification against an independent
reference value; the relativistic path has no such check in hand, and a
half-verified change to it is precisely the failure this arc has been
correcting (it already cost one revert-then-reapply cycle earlier in the
run). Recorded with numbers; owed its own verified pass.

Separately, `tab:spinor_resource`'s "scalar" column is **not comparable**:
it reports the relativistic spec evaluated by the spinless builder, which
returns the spinor-dimension count (1413/89226/942/59484) rather than the
scalar molecule's count (LiH is 837 at Q=30, 42138 at Q=84). So the
"rel/scalar 4.24x/11.33x" ratios are not the advantage they appear to be.
Both facts are now stated inline in the table's caption as a vintage note.

### Gates

C17 gained `composed-retired-scaling-and-counts` for the second-locus forms
the existing entry does not match (bare `1712\times` without the thin-space
comma, `1{,}111`, `Q^{3.15}`, `2.21 +- 0.02`, `6{,}000\times`, "exponent of
~2.2"). Discrimination proven both ways per the hard rule: **fires on all
six retired wordings, silent on all six corrected wordings**. It earned its
place immediately by catching the Paper 20 multi-center table.

Deterministic layer green on group4 and group6: C10 (three-pass, and
checked for undefined refs/citations rather than trusting exit code -- 0 of
each; P14 29pp, P20 12pp, P26 7pp), C16, C17, C18, C19, C11.

### Test-suite side effect of the exact rule (measured, not hidden)

`tests/test_ecosystem_export.py` appeared to hang. It does not -- it got
expensive, and by exactly the amount the correction predicts. Measured build
times for `ecosystem_export.hamiltonian('LiH')`:

| max_n | Q   | non-identity Pauli | build   |
|------:|----:|-------------------:|--------:|
| 2     |  30 |                837 |   0.3 s |
| 3     |  84 |             42,138 |  25.3 s |
| 4     | 180 |            750,774 | 830.8 s |

Under the retired rule LiH at max_n=3 was 7,878 terms; the exact rule makes
it 42,138 and the growth compounds. The implied local slope,
log(750774/42138)/log(180/84) = **3.78**, independently reproduces the
documented "local slope rises toward ~3.8 by n_max=4" -- a free cross-check
on the re-pricing from a completely different direction.

Three propinquity tests sweep max_n in {2,3,4} and so now cost ~14 min each.
Marked `@pytest.mark.slow` per CLAUDE.md SS14 (>10 s), with the measurement
in the marker comment so the skip is documented rather than silent. The
bound is still exercised in the default run by
`test_propinquity_bound_finite_and_positive` (max_n=2) and the
metadata/caching/universality tests. **The 830.8 s completion was measured
before marking** -- a slow marker on a genuinely hung test would be a
cover-up; this one is a cost statement.

File now runs **28 passed, 79 skipped in 16 s**.

### Verification at close

- Tests: test_paper14_eri_rule + test_paper20_balanced_lambda +
  test_tc_angular + test_h2_bond_pair_qubit = 35 passed, 14 skipped.
  test_paper26_entanglement = 7 passed. test_heavy_hydrides + the 18
  symbolic S^3 proofs (test_fock_projection + test_fock_laplacian) = 45
  passed. test_ecosystem_export = 28 passed, 79 skipped.
- Papers compile three-pass with **0** undefined references, **0** undefined
  citations, **0** undefined control sequences (checked in the logs, not by
  trusting pdflatex's exit code -- it exits 0 on dangling refs):
  P14 29 pp, P20 12 pp, P26 7 pp.
- Deterministic gates PASS on group4 and group6: C10, C11, C13, C14, C16,
  C17, C18, C19.
- Self-audit caught one of my own errors before it shipped: I wrote the
  closed form as 2.8157 when 1 + log_5(18.6) = **2.8163**. Corrected in the
  paper, these notes, and the CHANGELOG.

### Not done / owed

- The 17th site (relativistic factor order) needs an external verification
  route; `tab:spinor_resource` re-pricing is blocked behind it.
- `tc_integrals`' angular-gradient L_z violation is recorded, not fixed.
- CLAUDE.md SS1.6 and SS13.4a still cite the retired `O(Q^2.5)`. Both are
  PM-locked sections (SS13.5), so they are flagged to the PI rather than
  edited.

---

## Addendum 2 (2026-08-30, PI direction): record both orderings

PI call on the 17th site: "maybe we record both tests for now?" Done --
`tests/test_paper14_spinor_eri_ordering.py` (6 tests, 4.0 s, default run).

The point is to make the undecided state a **measurement rather than a
recollection**. It pins, for LiH_rel:

- what the two orderings SHARE: identical support (2676) and identical
  |ERI| entrywise -- the orderings permute signs only, which is exactly why
  every count-based figure in the re-pricing was safe while the values
  were not;
- the size of the open question: **1080/2676 = 40.4%** sign flips;
- the **control**: at n_max=1 (all-s spinors) the orderings must coincide,
  and are bit-identical -- which is what makes the n_max=2 divergence
  trustworthy rather than an artifact of the source substitution;
- the downstream numbers **both ways**: shipped 1413 Pauli / lambda 40.594,
  corrected 1501 / 39.533;
- and, most usefully, **the discriminators that do NOT discriminate**:
  particle exchange and hermiticity hold for BOTH orderings, parametrised
  over both so the negative is explicit. The module docstring carries the
  reason (R^k is pair-swap symmetric; the X_k(c,a) phase is evened out by
  M_J conservation), so the symmetry battery is blind *by construction* --
  a future attempt will not waste a cycle re-deriving that.

The corrected variant is built by **source substitution on the shipped
module**, not by forking it, so the two cannot drift apart in any respect
except the one under study.

**Tripwire, and it was tested.** If someone applies the correction, the
source anchor stops matching and the file fails with a message naming the
obligation (re-price tab:spinor_resource, rewrite this test to pin the
single surviving ordering, see task #13). Verified by simulating the fix:
all six tests fail with that message, and the module was restored
bit-for-bit (`git status` clean). An untested guard is worse than none.

Paper 14's vintage note now cites the test, so the open state is reachable
from the table itself and not only from this memo (C13 PASS).

---

## Addendum 3 (2026-08-30, PI-authorized): the DoD was steering reviewers wrong

Before firing the delta, the group4 **definition-of-done** -- the frozen
criteria the reviewers work against -- was found to encode two now-false
statements, both in the unsafe direction.

**(a) It exempted the exact class that turned out to be wrong.** Item 2
classed scaling exponents as *ROBUST under CF-1*, reasoning that "the
pair-diagonal vs global-M_L re-pricing is a constant factor within valence
class, so the log-log slope and the linearity are unchanged", and licensed
them as MEASURED. Every item in that exempted list moved: atomic
3.15 -> 3.8, composed 2.5 -> 3.17, N_Pauli 11.10 -> 27.90 Q, P20 two-point
2.21 -> 2.816. The "constant factor" premise is the same single-point
artifact already retired in CLAUDE.md. **INVERTED:** one tier, not two --
both exponents and multipliers are rule-sensitive. The 1-norm O(Q^1.69) and
QWC O(Q^3.36) exponents were never re-derived and are now marked
**UNVERIFIED**, not robust. Added reviewer instruction: a clean log-log fit
is not evidence -- this class read O(Q^2.5) for months *because* the fit was
clean.

**(b) It recorded a superseded decision.** Item 3 held "DECIDED A (disclose
the pair-diagonal rule); B deliberately left on the table". There was never
a choice: the rule was a wrong-sign Gaunt argument. **Replaced with the
dissolution**, and the cert criterion re-aimed -- what is MATERIAL is no
longer "failing to name which of two live rules is in use" (there is one)
but *any pair-diagonal-era number presented as live*, plus any surviving
framing that implies a choice exists. The known unremediated relativistic
instance is named so reviewers confirm the disclosure rather than re-flag
the numbers.

Three further loci repeated the same two retired premises (carried-decision
item 1, the branch-risk paragraph, the Paper 20 watch-note's "334 Pauli,
13x fewer QWC") and were corrected with them; the C8 keystone watch-note --
the enumerated headline list -- still carried the entire pair-diagonal-era
set and was re-priced.

**This is a correction toward reality, not a relaxation: the bar gets
stricter**, since the exponents move from the exempt tier into the
sensitive one.

### The registry earned its keep three times in one pass

1. C17 fired on **my own** new DoD text: I had written "(was X)", which is
   not in the exemption vocabulary. Conformed to "retired:" rather than
   widening the pattern -- "was" is far too common a word to exempt without
   blunting the gate. The gate was right; my wording was off-register.
2. Adding `docs/qa/group4.done.md` to the new entry's file list (it had been
   ungated -- a **coverage gap in a document that steers the reviewers**)
   immediately surfaced the stale C8 keystone list.
3. A new guarded family `composed-lih-market-test-retired` (334/333 Pauli,
   "constant 2.51x", "13x fewer QWC"; discrimination proven 4/4 fire, 4/4
   silent) immediately caught **two live zombies in Paper 14 itself** --
   lines 3097 and 3268 -- that the delta and two manual sweeps had all
   walked past. Re-priced: LiH composed 334 -> **838** Pauli, lambda
   32.6 -> **34.0** Ha, QWC 21 -> **64**; and the headline multipliers with
   them, "2.7x fewer Pauli / 13x fewer QWC" -> **near parity (1.08x) and
   4.3x fewer QWC**. A decisive-looking win the correction removes.

### Flagged, not reconciled

CLAUDE.md's canonical line gives the composed LiH 1-norm as **34.5 Ha**
(ratio 1.007 vs STO-3G's 34.3). The ecosystem path measures **34.0** (ratio
0.99), including the identity, on the same path and convention that
reproduces the 838 count exactly. A ~1.5% gap. The paper takes the measured
value; the discrepancy is raised for the delta reviewers rather than papered
over, because this session has already been bitten once by assuming a
convention (the k_orb factor of 2).
