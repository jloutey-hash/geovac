# Sprint memo — two production ERI evaluator defects (2026-08-29)

**Origin.** Cert-3 item 20 (Paper 26 §III counts an ERI tensor that omits the
Coulomb selection rule `m_a + m_b = m_c + m_d`). Remediating it uncovered a
second, larger defect in a *different* evaluator, and the pair together move
numbers in a CERTIFIED group4 paper.

---

## 1. The two defects

### (a) `geovac/casimir_ci.py::two_electron_integral` — two bugs

1. **Missing the M_L delta.** The multipole expansion of `1/r₁₂` carries a
   **shared** index `q` between the two Gaunt factors. The code called
   `_gaunt_ck(la,ma,lc,mc,k)` (which internally sets `q = ma−mc`) and
   `_gaunt_ck(lb,mb,ld,md,k)` (`q = mb−md`), letting each factor pick its own
   `q`. Mismatched-`q` products then survive as spurious nonzeros. The shared-`q`
   requirement *is* `m_a+m_b = m_c+m_d`.
   Witness: `⟨2p₊₁ 2p₊₁|2p₋₁ 2p₋₁⟩` (M_L +2 → −2) returned **0.021094**;
   physically it is exactly zero.

2. **Wrong Gaunt argument order (a SIGN error).** Condon–Shortley is
   `c^k(a,c)·c^k(d,b)`; the code used `c^k(a,c)·c^k(b,d)`. Since
   `c^k(1,2) = (−1)^(m₁−m₂) c^k(2,1)`, the two differ by `(−1)^(m_b+m_d)` —
   a sign flip on every integral with `m_b+m_d` odd (**1840 of 3800** at
   n_max=3).

   *This is not an absorbable phase convention.* A diagonal orbital rephasing
   gives `s*_a s*_b s_c s_d`; on `⟨2p₊₁ 2p₀|2p₀ 2p₊₁⟩` that is `+1` while the
   required ratio is `−1`. Contradiction ⇒ genuine error.

### (b) `_ck_coefficient` — dropped every m-changing multipole

Present **identically in two modules**: `geovac/lattice_index.py` (atomic) and
`geovac/composed_qubit.py` (molecular).

Shipped `q = mc − ma`. The 3j bottom row must sum to zero,
`(−m_a) + q + m_c = 0 ⇒ q = m_a − m_c`. With the sign flipped the 3j vanishes
identically unless `m_a = m_c`, so **every m-changing multipole was silently
dropped**: 88 of 126 nonzero `c^k` at l≤2 (**69.8%**).

> This is precisely the *"possible `_ck_coefficient` m-changing-multipole drop
> (needs verification)"* flagged in CLAUDE.md §2 at **v5.0.4** and never
> verified. It is now verified: real, and it reaches production.

**Reach.** `_ck_coefficient` is only invoked under `vee_method='slater_full'` —
which is the **default** in `frozen_core.py` and `locked_shell.py`, and is used
by `ecosystem_export.py::_build_he` (the 37-system `hamiltonian()` library) and
`vqe_benchmark.py`. `JordanWignerEncoder` consumes `li._eri` directly and
reports `eri_density` from it.

---

## 2. Arbitration — three independent confirmations

Neither evaluator got the benefit of the doubt.

1. **4-D quadrature from the definition.** Computed
   `A_k = ∫∫ Y*_a(1) Y*_b(2) P_k(cos γ) Y_c(1) Y_d(2)` by Gauss–Legendre ×
   trapezoid — no addition theorem, no 3j, no Gaunt. Verdict:
   **Condon–Shortley 4, casimir_ci 0** (1 tie on an m-diagonal case).

2. **Two independent evaluators converge — ON THE SPARSITY PATTERN.**
   *(Corrected 2026-08-29, later the same day: this originally read "agree
   exactly", which overstated it. What was verified is the nonzero COUNT and
   support; a follow-on B check found the two still disagree on VALUES for
   m-changing elements — see §7. Counts stand; values are open.)*
   After their (different) fixes, `casimir_ci` and `lattice_index` agree on
   the count:

   | | n_max=2 (M=5) | n_max=3 (M=14) |
   |---|---|---|
   | `lattice_index` as-shipped | 65 | 1,492 |
   | **physical truth** | **107** | **3,800** |
   | `casimir_ci` as-shipped | 265 | 15,364 |

   The two shipped evaluators **bracket** the truth — one over, one under.

3. **Variational physics.** He exact = −2.903724 Ha.

   | test | published | corrected | error vs exact |
   |---|---|---|---|
   | grid-hybrid n_max=5 | −2.893582 | **−2.896311** | 0.349% → **0.255%** |
   | direct-CI n_max=5 | −2.844830 | **−2.847551** | — |

   Restoring dropped couplings *must* lower a variational energy; both moved
   down toward exact and both **stay above** it. Bound respected.

---

## 3. Blast radius

### AFFECTED

- **Paper 26 §III** (group6, the cert-3 item): density **42.40% → 17.12%**,
  265 → 107 nonzero; n_max=4 318,720 → **57,700**. The *"essentially flat …
  does not degrade with basis size"* sub-claim **REVERSES** — corrected
  densities *improve*: 17.12% → 9.89% → 7.12% (M^−0.49).
- **Paper 26 core-closure theorem** — premise now fails at Z=7: 16
  ground-multiplet support determinants lack 1s². **But the violating weight is
  only 0.0014%** (5.79e−5 of 4.0; max amplitude 2.8e−3), so the theorem
  degrades from *exact* to *approximate*, not dead. Mechanism: the restored
  m-changing multipoles include `⟨1s 1s|2p₋₁ 2p₊₁⟩` — a core→valence double
  excitation, i.e. exactly the witness integral.
- **Paper 26 degeneracy table** {7:4, 8:6, 9:4} → measured {7:4, **8:3**, 9:4}.
  Z=8 moved; needs re-derivation.
- **Paper 26 Be I_cv** 0.0022 → 3.65e−3.
- **Paper 14 (CERTIFIED group4)** — see §4. PI-level.
- **CLAUDE.md §2** "Atomic Pauli O(Q^3.15), 1.3x–8.1x vs cc-pVDZ/cc-pVTZ".

### NOT AFFECTED (verified, two-way verdict)

- **Paper 26 §II / Table I** and **all entropy-locus work**: every CI/entropy
  caller filters to `m_i+m_j = 0`, so any matrix element joins two
  m_total=0 configs and M_L holds **automatically**. Verified by exhaustion:
  0 violations among 49 (n_max=2) and 961 (n_max=3) matrix elements.
- **Paper 22 (angular sparsity theorem, 1.44% at l_max=3): headline STANDS.**
  `tests/test_paper22_density.py` imports the `_gaunt_ck` *primitive* (correct)
  and applies global M_L itself. Only its provenance sentence is overstated —
  it says it ties the claim "to the exact routine that builds the physical
  `two_electron_integral`", but that routine did not impose the rule the test
  counts.

---

## 4. PI-LEVEL: Paper 14's atomic advantage inverts

Paper 14 states GeoVac ERI density "drops monotonically from **10.4% at M = 5**"
— exactly the defective 65/625.

| quantity | published | corrected |
|---|---|---|
| ERI density, M=5 | 10.4% | **17.12%** |
| He Q=10 Pauli | 120 | **288** |
| He Q=28 Pauli | 2,659 | **14,211** |
| density scaling | ~M^−0.85…−1 | **M^−0.49** |
| Pauli scaling | O(Q^3.15) | **Q^3.78** |

Head-to-head, against the paper's own Gaussian baselines:

- **Q=10: 288 vs 156 (cc-pVDZ) — the claimed 1.3× advantage becomes 1.8× WORSE.**
- Q=28: 14,211 vs 21,607 (cc-pVTZ) — 8.1× becomes **1.5×**.

The atomic advantage does not vanish, but it **inverts at small Q** and shrinks
5× at large Q. Held for PI decision rather than rewritten unilaterally: this is
a certified paper's central atomic selling point (§9 disposition rule — "raise
large issues to the PI"). Note the *composed* molecular results (O(Q^2.5), the
51×–1712× line) go through a different builder and are **not** implicated.

---

## 4b. The "pair-diagonal convention A" IS this bug

`tests/test_paper14_eri_rule.py` documents an A/B convention axis (CF-1):
rule **A** = pair-diagonal (the product, "makes the sparsity numbers what they
are"), rule **B** = exact global-M_L. Its own comment names the mechanism:

> `# Rule A (q=mc-ma) -> 3j vanishes; rule B keeps it.`

So "convention A" is not a free counting choice — **it is precisely the
wrong-sign `q`.** The docstrings in both modules state the *standard* formula
and then implement the flipped sign; `lattice_index` additionally emits a
UserWarning promising "proper angular coupling via Wigner 3j symbols". A
convention that only exists as a sign error, is described as exactness in its
own docstring, and was flagged at v5.0.4 as a suspected bug, is a bug that was
retroactively read as a convention.

**Disposition (deliberate split):**

- **Atomic path (`lattice_index`) — FIXED.** Rule A there was producing wrong
  *energies*, not merely a different count; the corrected He results are closer
  to exact, stay variationally bounded, and now agree with the independent
  graph-native h1 to 1.6e-4 Ha.
- **Composed path (`composed_qubit`) — NOT TOUCHED.** It carries the flagship
  molecular numbers (O(Q^2.5), 51×–1712×, the 37-system library,
  N_Pauli = 11.10 × Q) and the papers disclose A. Switching it is a PI call.
  Current state pinned by `test_composed_path_still_realizes_rule_a` so the
  two paths cannot drift apart silently while the decision is pending.

**Mitigation worth knowing:** the existing
`test_global_rule_reprices_constant_factor` measures A→B on the composed path
as a **constant factor** (2.51× main-group, 3.25× d-block) with *scaling
unchanged* — so the composed **O(Q^2.5) exponent likely survives** a switch;
only absolute counts and the Gaussian ratios move. **Caveat:** that test
measures the ratio at a *single* basis size (max_n=2) and infers "scaling
unchanged" from one point. On the atomic path the same A→B switch did **not**
preserve the exponent (Q^3.15 → Q^3.78). The composed claim needs a sweep
before it is relied on.

---

## 4c. Two further defects, both PRE-EXISTING (exposed, not caused)

In `debug/archive/misc/entanglement_first_row.py`, which backs Paper 26
Table II / Headline 4:

**(c) The "I_cv = 0" test was TAUTOLOGICAL.**
`compute_mutual_information_matrix` loops only over
`significant = np.where(s_single > 1e-8)[0]`, so a zero-entropy orbital's MI
row is zero **by construction**. 1s had exactly zero entropy, so the assertion
"I_cv == 0 EXACTLY for every member (the theorem's conclusion)" tested the
filter, not the physics — it could never have failed.

**(d) The entropy routines are mutually inconsistent.**
`compute_single_orbital_entropies` disagrees with `compute_subsystem_entropy`
on the *same* single-orbital entropy (Z=7 orbital (2,1,1): **0.7978 vs
0.3125**), and the resulting "MI" **violates subadditivity on 5 of 10 pairs**
(e.g. S_i=0.000186, S_j=0.797813, S_ij=0.312571, which requires
S_ij ≥ |S_i−S_j| = 0.797627).

**Attribution is clean:** both reproduce on *random* CI vectors having nothing
to do with the Hamiltonian, so they are pre-existing. The ERI fix only exposed
(c) by lifting S_1s above the 1e-8 threshold.

*Disposition:* the ceiling / supremum / spread legs are **quarantined**, not
re-pinned — pinning numbers from a routine in this state would manufacture
backing. The sound legs (degeneracy; core-valence leakage computed directly
from CI amplitudes; determinant residual; the bipartite 2·S_core I_cv) still
run.

---

## 4d. Paper 26 Headline 4 — what actually changed

Corrected first-row core-valence MI (bipartite 2·S_core, the sound route):

| | Be | B | C | N | O | F |
|---|---|---|---|---|---|---|
| I_cv | 3.65e−3 | 1.51e−3 | 1.05e−3 | 3.72e−4 | ~0 | ~0 |
| 1s occ | 1.999717 | 1.999887 | 1.999926 | 1.999975 | 2.000000 | 2.000000 |

- **The exact-closure threshold moves from B to O.** What replaces the step is
  a cleaner physical statement: I_cv decays **monotonically** as the core
  contracts, reaching the numerical floor at O. (Arguably a better result than
  the retired step — a trend, not a cliff.)
- **"All 4/6/4 determinants are exact ground eigenstates" now holds only at F.**
  Measured residuals: N (deg 4, 28 support dets) 8.86e−2…1.000; O (deg 3, 4
  dets) 8.93e−14…7.07e−1; F (deg 4, 4 dets) ~1e−14. Determinant-eigenstate-ness
  and core closure are **independent** — O has exact closure yet a genuinely
  superposed multiplet.
- **Z=8 ground degeneracy 6 → 3.**
- Li I_cv 0.2135 → 0.2271; Be 0.0022 → 0.00365.

---

## 5. PI decision + the A→B sweep (2026-08-29, second half)

**PI direction:** *exact rule everywhere, re-price everything honestly* — under
the newly adopted mission statement (CLAUDE.md §1: the forced/free atlas; a
correct exponent is the asset, a flattering wrong one is a liability).

### 5a. The decisive sweep (the paper's own protocol, `composed_lih_scaling_sweep`)

| axis | rule A | rule B |
|---|---|---|
| LiH (Q, N_Pauli), max_n=1..3 | (6,9) (30,333) (84,7878) | (6,9) (30,837) (84,42534) |
| fitted exponent | 2.539 | **3.172** |
| B/A ratio by max_n | — | 1.00× / 2.51× / **5.40×** |

- **The Q^2.5 headline does NOT survive** the exact rule on the basis axis.
  The old "constant factor 2.51×, scaling unchanged" mitigation was a
  single-point artifact: at max_n=2 the ratio is uniformly 2.51× across
  molecules, but along the basis axis it grows.
- **Honest cap on rule A itself:** its local exponent grows with basis
  (2.24 → 3.07 → 3.24 using the max_n=4 point Q=180, N=92,890) — "2.5" was
  always a small-basis fit, under either rule.
- Rule B's 3.17 is still **below the Gaussian 3.9–4.3**: the advantage
  survives re-priced, not inverted.

### 5b. The universality law SURVIVES, re-priced

N_Pauli/Q = **27.900 exactly** on every sampled library system (10 systems,
three periodic-table rows: LiH, BeH₂, H₂O, NaH, MgH₂, CH₄, NH₃, SiH₄, KH,
CaH₂). The structural linearity claim is intact; the coefficient moves
11.10 → 27.90. Re-priced counts at max_n=2: LiH 333→837, BeH₂ 555→1395,
H₂O 777→1953 (Gaussian-comparison ratios shrink by 2.51× at this basis:
51×–1712× → ≈20×–682×).

### 5c. The seven-copy map (final)

| module | verdict | backs |
|---|---|---|
| `lattice_index` | **FIXED** | atomic CI / He / P26 |
| `casimir_ci` (2 distinct bugs) | **FIXED** | P26 §III |
| `composed_qubit` | **FIXED** (PI-approved) | flagship composed, P14/P20 |
| `sturmian_secular` | **correct all along** — its q is *negated at the call site*; guard comment added so a mechanical sweep never "fixes" it into a bug | Paper 60 ✓ clean |
| `sturmian_solver` | FIXED | §3 dead-end BU-1 only (verdict there rests on Löwdin 1-norm, unaffected) |
| `nuclear/harmonic_shell` | FIXED | P22 HO verification path |
| `nuclear/potential_sparsity` | FIXED — **dissolves CF-1**: production enumerator, composed pipeline, and P22's headline global-M_L density are now one convention | P22 realized densities |
| `tc_integrals` | FIXED | old TC debug chain; the *tracked* xTC engine is separate and was correct |

### 5d. Third in-corpus corroboration, found in the audit

`geovac/xtc_angular_sparsity.py` (the tracked engine promoted at v5.0.9) is an
independent, **correct** implementation (`M_Λ = m_a − m_c` at the vertex) —
which is why the xTC memos measured **107/625**: the corpus already contained
the exact corrected count, in a certified artifact, without the discrepancy
against `casimir_ci`'s 265 ever being noticed. Independent-route cross-checks
are now a standing memory rule (`feedback_independent_route_crosscheck.md`).

---

### 5e. Full re-pricing record (second half of 2026-08-29)

**Library census** (`debug/data/library_census_exact_rule.json`, 33 composed
systems): N = 27.90×Q main-group / 30.03×Q d-block, exact; QWC 69 (main) /
109 (d) uniformly (pk-partitioned adds 1 → 70). Balanced builder (12 cells
re-measured via the pinning test's own fixture): isostructural invariance
survives exactly (NaH = KH = 575, HCl = HBr = 9,824, …); λ_ni within ~10% of
the retired values while counts grew 2.4–3.4×.

**Composed 1-norms** (with PK, per basis): LiH 10.15/39.41/202.38, BeH₂
119.51/263.72/641.15, H₂O 9,517.76/17,468.87/23,755.34 (max_n = 1/2/3);
H₂ 6.66/39.51/144.63 (max_n = 2/3/4). Notably H₂O's PK-inflated λ *dropped*
28,055 → 17,469 at max_n=2 (sign-order fix induces coefficient cancellation).

**Papers re-priced in place** (all compile clean, C10/C17/C19 green):
- **Paper 14** (6 passes): abstract, disclosure section rewritten as a
  correction record, tab:composed_pauli / tab:equal_qubit / tab:balanced /
  tab:second_row regenerated, d-block passage reversed, atomic section
  (head-to-head inversion honest per the benchmarking rule: Q=10 is 1.8×
  WORSE, Q=28 is 1.5× better), fitted-exponent display → the single 3.17,
  Track-CB and TC tables vintage-marked.
- **Paper 20** (3 passes): abstract, code listings, sec:eri_rule rewritten,
  tab:resources / tab:molecules / tab:transition_metals regenerated, 190× →
  76×, conclusion re-pointed at the linear system-size axis.
- **group4 synthesis**: 4 loci.
- **CLAUDE.md**: §1.5 positioning + §2 QC-status + best-results rows +
  Phase-4 Paper-22 line (1.44% → 6.06% global).
- **CHANGELOG v5.1.12** written; version bumped (patch — the PI may prefer
  to call this a minor: it is a correction that moves published numbers).

**Tests**: `test_paper14_eri_rule.py` rewritten around the exact rule
(7-implementation witness sweep + repricing pins + anti-regression guards);
`test_paper20_balanced_lambda.py` re-pinned + rule-A-resurfacing guard;
Paper 22's CF-1 regression test flipped to the resolved state (enumerator ==
global D, strictly-above discrimination).

**C17**: family `composed-rule-a-retired-figures` added; on first run it
caught **3 live loci my manual grep sweep had missed**; discrimination proven
both ways (fires on a planted 11.10×Q, silent on the corrected corpus).

**C19 catch**: a heredoc slipped past discipline once more during this work
and TAB-corrupted `\text` in P14 line 355; `check_latex_escapes.py` caught it
immediately. Three heredoc-mangling incidents this session — the standing
rule (`feedback_no_heredoc_backslashes.md`) is confirmed load-bearing; all
LaTeX-bearing edits in this sprint went through Write-tool scripts except the
three slips, each caught.

**Still open at memo time**: He n_max=5 atomic Pauli count (JW on the 519,585
corrected ERIs — heavy, running) → final tab:pauli/tab:eri_density rows +
regeneration of the two hardcoded Paper 14 figures
(`benchmarks/qubit_encoding/generate_paper14_figures.py`); the full-suite
re-pin sweep (suite runs much slower under the denser corrected evaluators —
the ab-initio PK class alone is 5m15s, passing); balanced QWC re-measurement
(marked "---" in the tables); the TC-composed pipeline re-measurement
(vintage-marked in P14); `O(Q^{1.69})` simulation-cost exponent re-fit
(vintage-marked).

---

### 5f. The symmetry discovery: 8-fold was an artifact; the true group is 4-fold

The exact rule revealed that the composed/balanced ERI tensors' **8-fold
permutational symmetry was itself an artifact of the bug**. Complex spherical
harmonics genuinely carry only the 4-fold group — verified bit-exact:
`⟨ab|cd⟩ = ⟨ba|dc⟩` (particle exchange, dev 0.0) and `⟨ab|cd⟩ = ⟨cd|ab⟩`
(hermiticity, dev 0.0), while the single-swap symmetries are broken (max dev
8.8e−2). Under rule A every multipole was m-diagonal, giving accidental
real-orbital symmetry. Two concrete consequences, both **named follow-ons**:

1. **`DirectCI4e` computes wrong energies on exact-rule tensors** — its
   closed-form same-spin block assumes 8-fold. Quarantined
   (`test_balanced_direct_ci`); fix = complex-orbital 4-fold treatment or a
   real-spherical-harmonic transform of the integrals.
2. **FCIDUMP export is lossy** — the format stores 8-fold-unique entries
   (assumes real orbitals), so the round-trip reproduces the 8-fold
   *projection*, not the tensor. Documented in
   `test_ecosystem_export` (round-trip now asserted against the projection,
   with a discrimination guard: the raw tensor must NOT round-trip). Fix =
   real-harmonic transform before export.

Both tests now pin the 4-fold structure as the invariant and treat a
*returning* 8-fold as a wrong-sign-q regression signal.

### 5g. Downstream re-pins (all MEASURED)

- Coupled composition (Track CB): cross/within ERIs 130/195 → **214/321 — the
  ratio is *exactly 2/3* under both rules** (Gaunt-structural, rule-invariant);
  coupled Pauli 854 → 2,702; ratio vs composed 2.56× → 3.23× (structural
  verdict unchanged); 1-norm 85.69 → 91.65.
- Nested (Track DF): all seven rule-A-derived gates re-derived (density 15%→25%
  on measured 17.1%; Pauli/Q 15→32; two-center 150→360 on 288; hetero
  2000→6000 on 4,743). **The Löwdin-inflation dead-end verdict survives:
  16.5× (was 14.3×)**, now pinned as a ratio.
- P26 molecular: S_bond 0.3033139 → **0.3298895, still exactly R-independent**
  (bit-identical blocks mechanism intact); S_core 0.0061155 → 0.0082183;
  bond/core ratio 49.6 → 40.1.
- H₂ bond-pair 111/2,627 → 279/14,178; ecosystem library dict (35 systems) and
  heavy-hydride pins re-set from the census; balanced screened 878 → 2,726;
  Breit-kwarg scalar 334 → 838; general builder 333/555/777 → 837/1395/1953.
- `test_composed_beh2` algebraic-overlap 5.18% mismatch: **PRE-EXISTING and
  independent of this arc** (its chain imports none of the corrected modules)
  — xfail-quarantined with attribution; own follow-up owed.

### 5h. C18 closure (the cert-3 owed gate item, done in passing)

Extended C18 with the LaTeX-tie (`$12$~months`) and symbolic (`$N$-month`)
escapes + a computational-runtime exemption (cross-line window; runtime
wall-clock figures are legitimate measurements, not project chronology);
selftest extended both ways. The extended gate immediately found **6 genuine
loci** the old one passed (5 in P34 — including the known L9022 — and 1 in
P32), all fixed; all 7 branches now PASS.

---

## 6. Honest scope

- **Fixed at theorem grade:** all defects, arbitrated by quadrature-from-
  definition + two-evaluator convergence + variational improvement + the
  xTC-engine corroboration (§5d).
- **MEASURED:** every corrected number above (incl. the A→B sweep and the
  27.900×Q universality).
- **Done:** all 7 production copies audited (6 fixed, 1 proven correct);
  P26 tests+§III corrected; He pins updated; mission statement adopted (§1);
  independent-route cross-check memory rule written.
- **Done (second wave):** the full composed-constellation re-pin (~60 tests
  across 16 files, all green on targeted verification); Papers 14/20/22/24/26/27
  re-priced/caveated and compiling; CLAUDE.md §1.5/§2/§3 re-priced; C16+C17
  registrations with two-way discrimination proofs (the C17 group4 family and
  the extended C18 each caught real loci the manual sweeps missed); the
  cert-3 owed gate items (C16 bare-graph entry, C17 DoD scope, C18 fix,
  criteria hygiene, arithmetic-audit dimension, seed plan, Minnesota
  cross-branch filing at CLAUDE.md §3 + P24 + P27) all closed.
- **In flight:** the final full-suite convergence pass (slow: the corrected
  evaluators make CI builds 3-6x denser; the ab-initio PK fixture alone is
  ~5 min).
- **Process note.** This class — a wrong sign inside a 3j argument, silently
  producing zeros rather than errors — is invisible to every gate we have. It
  was found only because cert-3 forced an independent re-derivation of a
  counted number. The v5.0.4 flag sat unverified because "needs verification"
  was written into a §2 bullet instead of a test. See
  [[feedback_deferral_is_churn]].

---

## 7. Follow-on B (2026-08-29, later): a THIRD evaluator discrepancy — OPEN

Follow-on B set out to fix the two 4-fold consumers (`DirectCI4e`, FCIDUMP)
via a real-spherical-harmonic transform. The transform works (unitary,
produces a real tensor, sparsity cost 1.02x) but does **not** restore 8-fold
symmetry — and chasing that exposed something larger.

**What was established (rigorously):**

1. **The correct Condon-Shortley assembly is `c^k(a,c) · c^k(d,b)`.**
   Rebuilding the He n_max=2 tensor with that order reproduces the
   twice-fixed `casimir_ci` evaluator to **0.0 exactly**; the order actually
   used by `composed_qubit` and `lattice_index` (`c^k(b,d)`) deviates by
   5.9e-2. Both orders preserve particle exchange and hermiticity, which is
   why no symmetry test ever caught it.
2. **`lattice_index` and `casimir_ci` disagree on VALUES**, not just order.
   On the witness `<1s 1s|2p-1 2p+1>`: `lattice_index` gives **+0.03415**,
   `casimir_ci` **-0.01707** — opposite sign *and* a factor of 2. Applying
   the (d,b) order to `lattice_index` fixed the SIGN (-0.03415) but left the
   factor 2. That factor is **unexplained**.
3. **The physicist-notation symmetry audit** (composed LiH): particle
   exchange, hermiticity and full reversal hold to 0.0; a<->c and b<->d break
   by 7.8e-1. For complex spherical harmonics a<->c is genuinely not a
   symmetry, so the 4-fold reading in §5f stands — but the real-harmonic
   transform does not restore a<->c either, which is *not* consistent with a
   genuine ERI tensor and points at the same assembly problem.

**RESOLVED (same session, PI-directed next step).** The factor of 2 was the
**orbital-exponent convention, not a bug**: `casimir_ci.two_electron_integral`
defaults to `k_orb=1.0` while `lattice_index` builds He at `Z=2`, and R^k
scales linearly with the exponent. At matched exponent (`k_orb=2.0`) the two
agree: witness -0.0341411 vs -0.0341496, and across the WHOLE tensor
max|dev| = 1.9e-4 with **107/107 sign agreement** -- 1.9e-4 being exactly the
grid-vs-exact R^k residual visible in the m-diagonal control (1.25000 vs
1.25019). So the only real defect was the Gaunt ORDER, the fix was correct,
and it has been **re-applied** to both modules.

*(Superseded note: the fix was briefly reverted while the factor 2 was
unexplained -- the right call on the information available, since a
half-applied fix in production is worse than either endpoint, but it cost a
cycle. The lesson is narrower than "be cautious": check the unit/normalization
convention of both routes BEFORE concluding two evaluators disagree.)*

**Re-priced after the order fix (COUNTS ALL UNCHANGED -- support is
order-independent, so every Pauli count, density, exponent and the
27.90xQ / 30.03xQ laws stand):**

- **QWC universality survives at new values: 64 (25 main-group) / 111
  (10 d-block)**, was 69/109. Greedy grouping is order-sensitive because the
  sign flips change which JW terms cancel.
- 1-norms moved 2-4%: ecosystem LiH 34.54 -> **34.04**, BeH2 68.81, H2O 372.46,
  NaH 171.99; composed lam_ni LiH 27.7 / BeH2 54.3 / H2O 186.8; balanced
  LiH 75.2 / BeH2 289.6 / H2O 1439.8; all 12 balanced tab:molecules cells.
- **He energies are STABLE**: grid-hybrid -2.896298 (was -2.896311), err
  0.2558% -- the m-changing terms are negligible in a 1s^2-dominated ground
  state, so the 0.255% headline holds.
- Papers 14 + 20 tables re-priced again (15 loci) and compile clean;
  tab:composed_onenorm vintage-marked (its lambda convention differs from the
  live builder's and awaits a convention-matched re-measurement).

**Superseded status line:** Applying the (d,b) order without
resolving the factor 2 would have left production in a state neither tested
nor verified; the modules are back to the state in which all of today's
measurements were taken, so every recorded number remains consistent with
the code that produced it.

**What this does and does not touch:**

- **UNAFFECTED — all counts.** The Gaunt order and the factor 2 do not change
  which entries are nonzero (same support). Every Pauli count, ERI density,
  QWC count, scaling exponent and the 27.90xQ / 30.03xQ laws re-priced today
  stand.
- **OPEN — all values.** 1-norms (34.5 Ha etc.), CI energies (including the
  He 0.349% -> 0.255% improvement), and entropies computed through
  `lattice_index`/`composed_qubit` are provisional until the factor 2 is
  explained. The `casimir_ci` path is the one with both fixes verified
  against quadrature.

**Named next step:** derive the composed/atomic assembly from the
Slater-Condon rule symbolically for one m-changing element and compare
term-by-term against `casimir_ci` — the factor 2 is most likely either a
double-counted (b,d)/(d,b) contribution or a convention difference in the
R^k radial normalization, and the two are distinguishable by inspection of
a single k-term.

**Process note.** This was found only because the independent-route
cross-check rule written this morning was actually applied: two routes
agreeing on counts is not agreement on values, and I had recorded the
stronger claim. Corrected in §2 and the CHANGELOG.
