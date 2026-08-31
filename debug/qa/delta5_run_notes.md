# Delta-5 verification run (2026-08-30)

**Verdict: DEFECTS.** The full certifying run stays blocked.

## Calibration: 4/4 seeds caught, both dimensions

| seed | class | planted in | catcher | result |
|:--|:--|:--|:--|:--|
| S1 | reversed conclusion (d-block "sparser") | P14 sec:block_pauli_count | claims | CAUGHT (M1) |
| S2 | optimality overclaim (kappa class) | P20 after eq:alpha_exact | claims | CAUGHT (M13) |
| S3 | tautological assertion (`not X or True`) | test_tc_angular | code | CAUGHT |
| S4 | retired pair-diagonal rule reinstated | tc_integrals M_L guard | code | CAUGHT |

Both dimensions calibrated 2/2, so the panel's other findings carry weight.
The code agent did not merely spot S4 by reading: it **independently
re-derived 107 via `casimir_ci.two_electron_integral`** and reported that
production gave 65 -- the retired count -- which is the strongest form the
catch could take. Seeds lived only in the worktree; the real corpus was
verified clean of all four before and after, and the worktree is deleted.
(One grep hit on the S4 pattern is `tests/test_paper22_density.py`'s
`_pair_diagonal_nonzero`, a deliberately-named comparison helper with an
accurate docstring -- not a defect.)

## The shape of the failure: propagation, not arithmetic

Every number that was recomputed checks out -- 279 = 55+224, 107 = 59+48,
85, 3.172, 2.8163, 27.90*Q, 30.03, 76x. The claims reviewer independently
re-derived `eq:eri_union` from the paper's own |P_k| values and got the same
59 + 48 = 107 split, and got 85 for the pure-d shell as a sum of squares of
m-multiplicities. **The arithmetic is sound; the sweep did not reach far
enough.** Retired figures survive at loci neither the delta nor two manual
passes touched, and in two cases they *reverse* the direction of the very
correction they were meant to record.

## Verified before reporting (PM re-checked against primary text)

- **M3 CONFIRMED.** P14 l.322 correctly reads "1.8x *more* ... (the retired
  pair-diagonal count, 120, had claimed a 1.3x reduction)"; l.826 still
  asserts "GeoVac achieves 1.3x fewer Pauli terms and 3.8x lower 1-norm".
  The retired claim is live with its direction inverted, in Limitations --
  exactly where the honest disclosure belongs.
- **M14 CONFIRMED.** 107382/1954 = 55.0x; the paper's 138x = 107382/778
  exactly, i.e. the retired count. Same paragraph re-priced LiH correctly
  (76x); H2O was missed.
- **M11 CONFIRMED, against the PM's own work.** 11.10*30 = 333 and
  11.10*20 = 222 exactly, so `tab:spinor_resource`'s scalar column is the
  retired *coefficient*, not the spinor-dimension artifact the new vintage
  note claims. Reason (ii) of that note is wrong and must be rewritten; the
  note's conclusion (read the table as pair-diagonal-era) survives, its
  stated mechanism does not.

## Ledger

Removing the two seeds (M1, M13), **21 genuine material-tier findings** plus
~9 NITs and 3 upgrade candidates.

Reviewer-recommended PI escalation: **M2** (live 51x-1712x), **M3** (He
direction inverted, 2 loci), **M4** (live composed 2.5), **M14** (138x),
**M20** (P26 `eq:hub` is the pre-correction reading), **M21/M22** (P26
"step-function" + unqualified "unique" -- both qualifiers the 2026-08-24
cert installed, stripped; the numbers were re-priced while the
characterization they were meant to correct was not).

PM-fixable: M5, M6, M7, M9, M10, M12, M15, M16, M17, M18, M19, M23 + NITs.
Needs adjudication: **M8** (the lambda / QWC family is live at 8 loci while
the corrected DoD marks it UNVERIFIED -- either re-measure or vintage-mark
uniformly) and **M11**.

## Two-way verdicts (the panel moved both directions)

- **U1:** `eq:eri_union` is a genuine *derivation*, and the paper calls it
  merely "reproducing the measured count". Under-stated.
- **U2:** the 55 direct terms are *provably* convention-independent (the
  monopole argument uses no sign convention and no m-rule) -- a stronger and
  more useful statement than "unaffected by the correction", especially
  given the whole re-pricing turned on a convention.
- **U3:** P20's abstract reports only the fitted 3.17 and omits its own
  stronger result, the forced alpha = 2.8163 identical across six molecules.

## Genuine finding against the PM's own remediation (not a seed)

The code reviewer showed that blanket-marking three propinquity tests
`@pytest.mark.slow` dropped the max_n=2 **value-correctness** check against
the closed form out of the default run; the surviving
`test_propinquity_bound_finite_and_positive` asserts only finite-and-positive.
The claim "still exercised in the default run" was true only in a weak
existence sense. **Fixed same session:** the two tests whose max_n cases are
independent are now parametrized so max_n=2 stays fast and only max_n>=3 is
slow; the monotone test is irreducibly three-point and stays slow as a unit.
File now runs 30 passed / 81 skipped, and the closed-form value check is back.
It also spot-checked the timing claim in that marker comment and found it
honest (17.9s vs the claimed 25.3s, term count 42,139 vs 42,138).

---

## Remediation (2026-08-30, PI-authorized): full ledger worked

All 21 genuine findings, the NITs, and the three upgrades. Every re-priced
number was MEASURED, not inferred.

### New measurements taken to close the ledger

| quantity | measured | retired |
|:--|--:|--:|
| He n_max=2, Q=10 | 287 Pauli, lambda 11.175, QWC 67 | 120 / 11.29 / 25 |
| He n_max=3, Q=28 | 14,078 Pauli, lambda 74.207, QWC 5,569 | 2,659 / 78.36 / 791 |
| isolated s-only block (M=3) | 100% density, 117 Pauli, /Q 19.50 | 118 (identity incl.) |
| isolated s+p block (M=9) | 15.7%, 2,967 Pauli, /Q 164.83 | 8.9% / 976 / 54.22 |
| isolated d-only block (M=5) | 13.6%, 343 Pauli, /Q 34.30 | 4.0% / 56 / 5.60 |
| hub shares at n_max=2 | B(Z=5) 1s 0.1% 2s 19.0%; C(Z=6) 1s 0.0% 2s 3.7% | -- |

The block method was validated before use: it reproduces the paper's own
s-only row exactly (117 non-identity + identity = the tabulated 118).

### The largest corrections

- **The d-block reverses twice over.** `tab:d_block_eri` said the d-only
  block is "the sparsest" at 4.0% density and Pauli/Q = 5.60. Measured:
  13.6% and **34.30** -- denser per qubit than the s+p block's 27.90, and
  consistent with the composed d-block coefficient (30.03 vs 27.90). The
  section lead ("Gaunt selection rules become MORE restrictive at higher
  angular momentum", 20% vs 38.5%) was the same artifact and is rewritten.
- **The exponent gap mostly closed.** P14 claimed a "roughly 1-1.5 in
  exponent ... robust regardless of the per-integral counting" gap. With
  3.8 against 3.9-4.3 it is **0.1-0.5**, and it is precisely NOT robust to
  the counting -- the retired 3.15 is what made 1-1.5 true. Both the intro
  claim and the conclusion bullet are corrected, and the honest surviving
  statement is named: the exact Q-linearity across molecules, not the
  within-molecule exponent gap.
- **He equal-qubit, direction inverted at two loci.** Line 826 asserted
  "1.3x fewer Pauli terms" where the same paper's l.322 already says
  **1.8x MORE**; the Limitations summary said "1.3x and 8.1x fewer". Both
  now carry the measured 1.8x-more / 1.5x-fewer split. The 1-norm
  advantage survives and at Q=28 **improves** (6.8x -> 7.2x) -- one of the
  few comparisons the correction makes better.
- **P20's 370x extrapolation withdrawn, not re-priced.** The paper states
  two lines earlier that the exponent is not constant in n_max; using the
  local ~3.8 in the same formula gives ~4x instead of 370x. And Q=100 is
  reached by ADDING BLOCKS, where the count is exactly linear (exponent 1,
  measured). Neither within-molecule exponent is the right instrument, so
  the extrapolation is dropped and the linearity stated instead.
- **P26's withdrawn qualifiers restored.** "Step-function" and unqualified
  "unique"/"any unitary rotation" were reinstated at two loci while the
  numbers around them were re-priced -- the abstract and eq:step's own
  lead-in both say "rapid but graded" and "essentially unique up to
  (l,m)-label-preserving rotations". The sentences exist to BOUND the
  uniqueness claim; as written they removed the bound.
- **eq:hub was the pre-correction reading.** It asserted the 1s->2s shift
  at Z=3 (denied) and 2s freezing at Z=4 (the opposite of measured).
  Rewritten as 1s -> {1s,2s} -> 2s -> 2p with the measured shares, and the
  "at each transition the previous shell freezes" generalization scoped:
  it holds at every transition EXCEPT the first, since 1s is at 89.8% at
  Z=3.

### M8 adjudicated: marked, not faked

The lambda/QWC/Trotter exponent family could not be honestly re-derived --
SS13.4a wants >= 4 points and the He series needs Q=110, which is out of
reach (the Q=60 QWC grouping did not finish). So the family is
vintage-marked at its **definition site**, where the scope governs all ~18
downstream mentions, plus the abstract (read standalone) and the two tables
that tabulate it. What WAS re-measured -- the lambda values themselves --
is corrected everywhere it appears. `tab:qwc`'s n_max=4 row is left marked
because re-grouping now hits the same O(N^2) wall the caption already
documents for n_max=5; the groups/terms ratio rises at both re-measured
points, so the reported grouping efficiency is if anything optimistic.

### Criteria corrected (they were mis-directing reviewers)

`group6.done.md`'s C8 block for Paper 26 asserted the **retracted** ln 8 /
ln 16 entropy ceilings, an I_cv of ~2e-3 at Z=4 (measured 3.65e-3),
"exactly 0 for Z>=5" (measured 1.71e-3 at Z=5; the paper says Z>=7), a 50x
core/bond ratio (paper: 40x), and a tie to the retired O(Q^2.5). Also
corrected: the N/O/F degeneracies, where the criteria's combinatorial
"4/6/4 = C(4, n_p-2)" loses to the measured **4/7/6** that
`test_paper26_entanglement.py` pins and passes.

### Upgrades applied (verdicts moved both ways)

- **U1** `eq:eri_union` restated as a **derivation**, not a check: its
  right-hand side is built only from the selection rules and the M_L
  constraint, with no appeal to the assembled tensor.
- **U2** the 55 direct terms restated as **provably convention-independent**
  -- the monopole argument invokes neither a sign convention nor an m-rule,
  so given that the whole re-pricing turned on a convention, this is the one
  part of the Pauli count that could not have been at risk.
- **U3** P20's abstract now carries the forced alpha = 2.8163.

### Verification at close

- Deterministic gates: **7/7 PASS on group4 and group6**.
- Papers compile three-pass with **0** undefined refs/citations/control
  sequences (checked in the logs, not via exit code): P14 30 pp, P20 12 pp,
  P26 7 pp.
- Tests: 53 passed, 2 skipped -- including the 18 symbolic S^3 proofs, the
  C17 registry self-test, and the both-orderings spinor record.
- Two framing zombies fixed in section headings ("the pair-diagonal ERI
  *approximation*" -> "the retired pair-diagonal ERI rule: a sign error,
  not an approximation").

### Still owed

- The lambda / QWC exponents want a >= 4-point re-fit once He at Q=110 is
  reachable.
- Task #13 (the 17th site) still needs an external discriminator.
- A fresh delta over THIS remediation before any certifying run.

---

## Addendum: the long-running job finished, and M8 flipped from marked to measured

The lambda/QWC re-measurement completed after the remediation was written.
It reached Q = 60 (n_max=4) after a long greedy-grouping pass, then died on
n_max=5 (Q=110, 519,585 nonzero ERIs) before writing its fits. Three points,
not the four SS13.4a wants -- but three measured points beat a vintage mark,
so M8 is re-dispositioned from **marked** to **measured-but-underpowered**.

| metric | exact rule (Q=10/28/60) | R^2 | max abs log resid | retired |
|:--|--:|--:|--:|--:|
| N_Pauli | **3.779** | 1.0000 | 0.0012 | 3.15 |
| lambda  | **1.792** | 0.9997 | 0.034  | 1.69 |
| QWC     | **4.013** | 0.9976 | 0.202  | 3.36 |

Two things fall out, and the second is a genuine loss.

**A free cross-check.** The Pauli exponent 3.779 was obtained here by a
route independent of the 3.8 already written into the papers, and agrees.
That is the [[feedback_independent_route_crosscheck]] discipline paying off
without being asked for.

**The measurement-group advantage does not survive.** QWC moves 3.36 ->
4.013, which puts it ABOVE the Pauli term exponent (3.779) and INSIDE the
Gaussian O(Q^4-5) band rather than well below it. Under the retired rule it
sat below both. So the shot-budget advantage claimed for VQE is withdrawn,
not re-priced -- and the groups/terms ratio rises at all three points
(0.208->0.233, 0.297->0.396, 0.329->0.344), i.e. the exact rule makes
grouping *less* efficient, consistently.

**This also falsified one of my own delta-5 fixes.** Earlier in this
remediation I had corrected the QWC-vs-Pauli sentence to read "slightly
below the Pauli term exponent of 3.8 but well below the Gaussian
O(Q^4-5) regime". With the measurement in hand that is wrong in *both*
clauses. Corrected again to "slightly above ... and inside rather than
below". A fix applied from the corrected-but-unmeasured side of a question
is still a guess.

Also re-priced from the same run: `tab:qwc`'s n_max=4 row is now MEASURED
(250,402 terms / 86,224 groups / 0.344, retired 31,039/10,199/0.329) rather
than marked unmeasurable; `tab:scaling_summary` carries the three-point
exact-rule exponents with its reduced range stated (the retired fits ran to
Q=110, which the exact rule puts out of reach); and 13 downstream
`Q^{1.69}` mentions moved to `Q^{1.79}` -- the Trotter-step exponents by the
paper's own r ∝ lambda derivation, not by assumption.

Gates: 7/7 on group4 and group6. P14 compiles clean at 30 pp.
Data: `debug/data/exact_rule_lambda_qwc.json`.

**Standing debt unchanged in kind, reduced in size:** the four-point fit
still wants He at Q=110.

---

## Debt paydown (2026-08-30): the 17th site is RESOLVED

The blocker was never effort -- it was that **no internal check could
distinguish the two orderings**. Particle exchange and hermiticity hold for
both, which the M_J-constrained phase algebra predicts and a direct test
confirmed. An external reference was required.

**The discriminator.** The jj-coupled 2p shell (kappa=+1, j=1/2; kappa=-2,
j=3/2 -- six spinors) spans the *same* six-dimensional space as the scalar
2p spin-orbitals (m_l in {-1,0,1}) x (m_s in {+-1/2}). The two bases are
related by the Clebsch-Gordan unitary. The Coulomb operator is
spin-independent, so in the scalar basis its angular part is the
Condon-Shortley product c^k(a,c)*c^k(d,b) -- whose ordering is ALREADY
verified against an independent reference value by the seven-implementation
witness sweep. Rotating that tensor into the jj basis gives exactly the
object the relativistic module must reproduce. Restricting to one shell
makes every radial integral identical, so the comparison is purely angular,
and the decisive statistic is the sign pattern, which no rescaling can
alter.

**Verdict, at k=2 (168 comparable entries):**

| ordering | sign agreement | best-fit scale | residual after scaling |
|:--|--:|--:|--:|
| corrected `X_k(a,c) X_k(d,b)` | **168/168 (100%)** | 1.0000 | **1.7e-16** |
| shipped `X_k(a,c) X_k(b,d)` | 96/168 (57.1%) | 0.2000 | 0.90 |

The CG rotation is unitary to 1.1e-16, and **k=0 is an exact control**: both
orderings match the reference there, as they must, since a scalar monopole
makes them coincide. So the discriminator is calibrated by construction --
it is silent exactly where it should be silent, and decisive where it
should be decisive.

**Applied.** `geovac/composed_qubit_relativistic.py` now uses the corrected
order, with the justification cited in-line. The proof was written into
`tests/test_paper14_spinor_eri_ordering.py` BEFORE the change it justifies;
that file no longer records an open question but proves the resolution, and
it retains the excluded ordering as a live test (its exclusion is the
load-bearing evidence, so deleting it would discard the reason). It also now
asserts that the production source actually carries the corrected loop --
the fix cannot silently regress.

Re-pinned: LiH_rel 1413 -> **1501**, lambda 40.59 -> **39.53**, QWC 121 ->
**142**; CaH/SrH/BaH_rel 942 -> **998** (isostructural, 7 pins).
`tests/test_heavy_hydrides.py` + `test_spin_ful_composed.py`: 41 passed.

**The tripwire worked exactly as designed.** Applying the fix made the
earlier both-orderings test fail with a message naming the obligation
(re-price tab:spinor_resource, rewrite the test). That is what it was for.

---

## Debt 2 paid: the four-point fits

The blocker was diagnosed rather than brute-forced. The n_max=5 run had died
in `count_qwc_groups`, which is O(N^2) greedy -- **not** in the Hamiltonian
build or the 1-norm. Dropping only the QWC leg reached Q=110 and gave the
fourth point for the two exponents that carry weight: the term count and the
fault-tolerant-relevant 1-norm.

He, exact rule:

| Q | N_Pauli | lambda | QWC |
|--:|--:|--:|--:|
| 10 | 287 | 11.175 | 67 |
| 28 | 14,078 | 74.207 | 5,569 |
| 60 | 250,402 | 275.718 | 86,224 |
| 110 | **2,434,441** | **790.007** | -- (2.4M terms; O(N^2) infeasible) |

| metric | alpha | R^2 | max abs log resid | points | retired |
|:--|--:|--:|--:|--:|--:|
| N_Pauli | **3.773** | 1.0000 | 0.0064 | 4 | 3.15 |
| lambda | **1.774** | 0.9997 | 0.0421 | 4 | 1.69 |
| QWC | 4.013 | 0.9976 | 0.2024 | 3 | 3.36 |

**The fourth point confirms rather than moves the fit.** Three points gave
3.779 and 1.792; adding Q=110 shifts them by 0.006 and 0.018. So the
correction's exponents were already stable at three points -- worth knowing,
because it means the SS13.4a four-point bar was protecting against a risk
that had not in fact materialised here.

Both now meet that bar. QWC stays at three points for a stated structural
reason (2.4 million terms against an O(N^2) grouping), not for want of
effort.

Applied to Paper 14: `eq:exact_exponents` and `tab:scaling_summary` carry
the four-point values with the point counts distinguished; 14 downstream
`Q^{1.79}` mentions moved to `Q^{1.77}`; the QWC-vs-Pauli comparison now
cites 3.77. The two-significant-figure headline `Q^{3.8}` is left alone --
3.773 rounds to it, and churning it would be noise.

Gates 7/7 on group4 and group6. P14 compiles clean at 30 pp.
Data: `debug/data/exact_rule_lambda_qwc.json`.

**Both standing debts from the delta-5 remediation are now closed.** What
remains open is a fresh delta over this work before any certifying run.
