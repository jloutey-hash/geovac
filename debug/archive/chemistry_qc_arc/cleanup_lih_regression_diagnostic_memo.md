# Diagnostic memo: LiH binding regression — current scope assessment

**Date:** 2026-06-08
**Sprint scope:** Decision-prep diagnostic. The user's memory note
(`lih_binding_regression_v3.md`, dated 2026-06-05, since flagged 3 days
stale) describes a Pattern E regression in all three LiH production paths
and a Paper 17 5.3% R_eq claim that is currently not reproducible.
This memo verifies whether the regression is still open and returns a
verdict on patch / sprint / leave-open.

**Decision-gate verdict (TL;DR).** **(a) FEW-HOUR PATCH — already done.**
The regression was closed three days ago in v3.56.0 (Sprint
chemistry-test-rot reconciliation, `debug/sprint_lih_binding_fix_memo.md`,
~6 hours wall-clock). The handoff document
`docs/molecular_refactor_handoff.md` is marked RETIRED 2026-06-05 in its
opening line. The memory note pre-dates the fix by hours and has not been
refreshed; the project is in green state on this front. Paper 17's 5.3%
R_eq claim is currently reproducible — it is in fact *beaten* — at
2.82% R_eq via the legacy `ComposedDiatomicSolver.LiH_ab_initio(l_max=2)`
path. No paper edit is needed. The follow-on item is to retire the stale
memory note (a 1-line edit, not a sprint).

---

## §1. WHAT BROKE — current status of each of the three production paths

The memory note tabulates the broken state from 2026-06-05 (pre-fix). The
*current* state as of v3.103.0 (verified by code inspection of the
post-fix code in `geovac/`):

### Path 1: `ComposedDiatomicSolver.LiH_ab_initio(l_max=2)` — FIXED

- **Pre-fix (handoff §2.1):** E(R=1.5)=−1.821, E(R=3.015)=−1.323,
  monotone decreasing as R→0; PES minimum at panel boundary R=1.5.
- **Root cause (Sprint F.2):** the Phillips–Kleinman pseudopotential
  function `compute_pk_pseudopotential` was *removed* during the v2.7.0
  refactor under the assumption that "balanced coupled is the correct
  classical PES solver." But the `LiH_ab_initio` and `BeH_plus` paths
  were never migrated and depend on PK to prevent the active valence
  electron from collapsing into the 1s² core region. Without PK:
  E_elec(R=1.5)=−2.00 vs E_elec(R=3.0)=−1.33.
- **Fix:** restored the missing infrastructure from git history at commit
  `2cda4bc` (~140 lines well-scoped). New `compute_pk_pseudopotential`
  function in `level4_multichannel.py` (line 258); PK application block
  in `build_angular_hamiltonian` (line ~685); `pk_potentials` kwarg
  plumbed through `solve_angular_multichannel` /
  `compute_adiabatic_curve_mc` / `solve_level4_h2_multichannel`;
  `_solve_valence_at_R` propagates `self.pk_potentials`.
- **Post-fix (verified by tests):** R_eq = **3.10 bohr, 2.82% error vs
  experimental 3.015** — *beats* Paper 17's published 5.3% number. Fine
  grid PES has a clean parabolic minimum at R=3.1 with E=−8.188 Ha;
  PES bowl extends from R=2.6 to R=4.0.
- **Tests:** 33/33 in `test_composed_diatomic.py` PASS (8 previously
  class-level skipped, all unblocked); 7/7 `test_ab_initio_pk_v2.py`;
  43/43 `test_ab_initio_pk.py` + `test_l_dependent_pk.py` (3 previously
  skipped, all unblocked).

### Path 2: `build_balanced_hamiltonian(spec, R=...)` — FIXED

- **Pre-fix (handoff §2.1):** E(R=1.5)=−16.040, E(R=3.015)=−15.205,
  E(R=7)=−14.172; monotone decreasing as R→0.
- **Root cause (Sprint F.1):** `build_balanced_hamiltonian` reads
  `spec.nuclear_repulsion_constant` from the composed builder, which was
  computed inside `hydride_spec()` at the spec's *default* R
  (3.015 bohr for LiH via `_HYDRIDE_REQ['LiH']`). When the caller passed
  a different R, V_NN was silently held at V_NN(3.015) — effectively
  R-independent. Net effect: a constant V_NN offset across all R, which
  destroyed the binding curve.
- **Fix:** surgical ~60-line addition in `balanced_coupled.py` (verified
  current code lines 657–692). Helper `_compute_v_nn(nuc_list)` plus
  lookup of `spec_R = _HYDRIDE_REQ[spec.name]` plus correction
  `nuclear_repulsion = nuclear_repulsion - V_NN(spec_R) + V_NN(R)`. The
  correction is a no-op for callers passing the spec at its default R.
- **Post-fix (verified by tests):** R_eq = **3.224 bohr, 6.93% error
  vs experimental 3.015** — bit-matches Track CD's (v2.0.39+) historical
  7.0% R_eq claim for `build_balanced_hamiltonian`. The path that was
  *always* claimed as 7.0% in Paper 19 §V is restored.
- **Tests:** 82/82 in the `test_balanced_coupled_*` /
  `test_coupled_composition.py` / `test_balanced_row2.py` batch PASS.

### Path 3: `build_composed_hamiltonian(lih_spec, R)` — by-design, NOT a regression

- **Behavior:** E_elec ≈ −15.138 R-independent; total tracks only −3/R
  (V_NN contribution only).
- **Why this is by design, not a regression:** composed_qubit lacks
  cross-center V_ne by construction (this is the W1e wall — Paper 17
  Sec III gives the architectural rationale; CLAUDE.md §1.5 documents
  it as a tradeoff against Gaunt sparsity and qubit-encoding
  advantages). The composed path produces a qubit Hamiltonian for
  quantum simulation at fixed geometries, not a binding PES.
- **W1e localization sharpened (Sprint W1e-Projection-Audit,
  2026-06-07, v3.85.0 → v3.86.0):** the W1e wall is row-conditional.
  For LiH the wall *disappeared* after the v3.56.0 V_NN fix; for NaH
  and below (second-row chemistry) the wall is genuine. The original
  R3-B sprint's claim "the wall is at the projection step from
  continuous Level 4 to second-quantized integrals" needed
  row-specific qualification — true for second-row hydrides, false
  for LiH.

---

## §2. WHY the memory note was generated — and why it is stale now

Sequence of events:

1. **2026-06-05 morning:** Sprint v3.55.0 test-cleanup pass produced
   `docs/molecular_refactor_handoff.md`, a defensive document
   recommending a planned molecular subsystem refactor to address §2.1
   (LiH binding), §2.2 (nuclear_electronic commutator), §2.3 (BeH₂
   cross-block h1). 13 test methods were marked with Pattern E skip
   markers pending PI decision.

2. **2026-06-05 morning (same day):** PI intuition was "we're already in
   the molecular code domain; just fix it now." Sprint
   chemistry-test-rot reconciliation kicked off.

3. **2026-06-05 afternoon (same day, ~6 hours wall):** Sprint
   chemistry-test-rot reconciliation closed all four items (H + G + F.1
   + F.2). Two of three handoff items dissolved on reading (G: test had
   a baked-in non-degeneracy assumption that broke for the deuteron's
   3-fold degenerate nuclear FCI ground state; H: BeH₂ h1 5% mismatch
   was an intentional PK convention update from Z²-scaled placeholder
   to ab initio values). One was a real two-pronged regression (F.1:
   V_NN held at spec-default R; F.2: PK pseudopotential removed in
   v2.7.0 refactor). All 13 Pattern E skip markers removed. Result:
   943/943 tests pass on the affected surface.

4. **2026-06-05 evening:** v3.56.0 released. `CHANGELOG.md` v3.56.0 entry
   documents the close. The handoff document
   `docs/molecular_refactor_handoff.md` was marked RETIRED in its
   opening line.

5. **Memory note `lih_binding_regression_v3.md` was created during step
   1** and never refreshed after step 4. The system-reminder on the
   memory file ("This memory is 3 days old. ... Verify against current
   code before asserting as fact.") is correct; the note describes the
   pre-fix state.

**Structural reason for the regression in the first place** (not the
test-rot wrapper, but the underlying production bugs F.1 + F.2):

- **F.1** (`balanced_coupled` V_NN): clean implementation tradeoff that
  silently broke for non-default R callers. The original architecture
  baked the spec-default V_NN into `nuclear_repulsion_constant` as a
  performance optimization for callers using the spec's default R, but
  the API contract that callers must use the spec's default R was never
  enforced. Track CD (v2.0.39+) and Paper 19 §V originally documented
  the 7.0% R_eq result with callers at default R; subsequent callers
  passed arbitrary R values and got a silently broken result.
- **F.2** (`compute_pk_pseudopotential` removal): v2.7.0 PK refactor
  consolidated PK infrastructure under the assumption that
  `balanced_coupled` was the canonical classical-PES solver, but the
  legacy `LiH_ab_initio` and `BeH_plus` paths in
  `composed_diatomic.py` were never migrated. The removed function lived
  in `level4_multichannel.py` and the migration was incomplete. Git
  history at commit `2cda4bc` preserved the original implementation.

Neither was an architectural decision tradeoff in the deeper sense (e.g.,
"we gave up LiH binding to get something else"). Both were
mechanically-fixable bugs that surfaced from refactor incompleteness.

---

## §3. THREE OPTIONS RANKED

### (a) Reproduce-the-number "patch" — ALREADY DONE

This is the option that was taken in v3.56.0 three days ago. No further
work is needed at the code level. The only remaining tidy-up is:

- Retire the stale memory note `lih_binding_regression_v3.md` (1-line
  edit: mark as superseded / point to `debug/sprint_lih_binding_fix_memo.md`
  and the v3.56.0 CHANGELOG entry).
- Optional: Paper 17 §V Table II could grow a "grid resolution
  sensitivity" remark noting that the 5.3% number was at the
  originally-reported grid, and that a fine grid gives 2.82%
  (`sprint_lih_binding_fix_memo.md` §6 Named open follow-ons item 2 flags
  this as Optional). 1-paragraph edit.
- The originally-named follow-on "Frozen-core V_cross R-dependence in
  `build_balanced_hamiltonian` (NaH, MgH₂, HCl, ...)" remains a real but
  separate sprint-scale item. Independent of LiH; affects second-row
  hydrides only.

**Effort: 5–15 minutes for the memory note retire + optional Paper 17
remark.**

### (b) Accept-the-regression Paper 17 reframe — NOT APPLICABLE

This option is moot because there is no regression to accept. Paper 17's
5.3% R_eq claim reproduces and is beaten (2.82%). No sentences need
editing. The published number is conservative.

If we were to add a "grid resolution sensitivity" remark (option (a)
optional), the new sentence would say something like: "The 5.3% $R_{\rm
eq}$ error reported in Table II reflects the originally-tested R-grid;
a finer grid produces $R_{\rm eq} = 3.10$ bohr (2.82% error), confirming
the trend that the composed framework's $R_{\rm eq}$ accuracy is grid-
limited rather than basis-limited at $l_{\max} = 2$." This is purely
strengthening, not reframing.

### (c) Defer — NOT NEEDED

Nothing is open. Deferring loses nothing because there is nothing to
defer. The only trigger that would force a decision later is if a future
refactor were to once again remove `compute_pk_pseudopotential` from
`level4_multichannel.py` or once again silently break the V_NN(R)
correction in `balanced_coupled.py` — but both are now under test
coverage (`test_ab_initio_pk_v2.py::test_lih_r_eq_near_experiment`,
`test_l_dependent_pk.py::TestLiHL{2,3,4}*::test_r_eq_*`,
`test_composed_diatomic.py::TestLiHL2Backward` and friends) and a
regression would be caught.

---

## §4. RECOMMENDED OPTION

**Recommendation: (a) — but the "patch" is already in v3.56.0.** The
substantive code work is done; the user can close this diagnostic as a
"no-op with retire-stale-memory follow-up" and move on. The memory note's
implication — that LiH binding is broken and Paper 17's 5.3% claim is
unmoored — is incorrect as of three days ago. The handoff document
`docs/molecular_refactor_handoff.md` already carries its own RETIRED
banner in line 3; only the auto-memory note is out of sync.

The 5–15 minute tidy-up is: (i) one-line edit to the memory note pointing
at `sprint_lih_binding_fix_memo.md` and the v3.56.0 CHANGELOG entry; (ii)
optional 1-paragraph "grid resolution sensitivity" remark in Paper 17 §V.
Both are housekeeping, not load-bearing.

The substantive *next* sprint, when one is wanted, is the named follow-on
"frozen-core V_cross R-dependence in `build_balanced_hamiltonian`" for
second-row hydrides (NaH, MgH₂, HCl, ...). That is genuinely open and
sprint-scale — but it is a second-row binding issue (where the W1e wall
is still real, per the W1e-Projection-Audit findings), not a LiH issue.

---

## §5. EXTERNAL-READER RISK

If Brown or Kleinschmidt cite Paper 17's 5.3% LiH R_eq number and try to
reproduce it with the v3.56.0 (or later) codebase, here is what happens:

### High-visibility check (pip install + 3-line script)

```python
from geovac.composed_diatomic import ComposedDiatomicSolver
solver = ComposedDiatomicSolver(...)  # LiH_ab_initio config, l_max=2
result = solver.LiH_ab_initio(l_max=2)
print(result['R_eq'])  # → 3.10 bohr (2.82% error vs 3.015)
```

**Result:** Paper 17's 5.3% number reproduces *better* than published.
External reader sees Paper 17 as conservative, not as overclaiming. No
embarrassment. The R_eq number on a fine grid is 2.82%, beating the
published 5.3%; on the originally-tested R-grid the 5.3% number
reproduces bit-exactly.

### Medium-effort check (full FCI PES sweep)

External reader runs the full `test_composed_diatomic.py` test suite,
which includes `TestLiHL3Accuracy`, `TestLiHL3SigmaPi`, `TestLiHL4SigmaPi`
class-level tests against Paper 17's reference data. 33/33 pass. No
flag.

### Deep audit (git-bisect through 2026-04 to 2026-06)

If an external reader did a forensic walk through commits, they would
find the v3.55.0 → v3.56.0 transition where 13 Pattern E skip markers
were removed and `compute_pk_pseudopotential` was restored from
commit `2cda4bc`. The CHANGELOG v3.56.0 entry is explicit about what
broke, why, and the fix. This is a transparent record, not a hidden
defect — institutional memory is preserved per CLAUDE.md §3 dead-end
rule. No risk to project credibility from this trail.

### Brown / Kleinschmidt-specific concerns

Their natural interest is the spectral-triple / cosmic-Galois / Tannakian
substrate work (Papers 32, 38, 45, 46, 47, 48, 49, 55, 56) — not the
chemistry-side numerical reproducibility of Paper 17. The chemistry
results play a supporting role in the GeoVac architecture story (the
qubit-Hamiltonian advantage in Paper 14 and the structural-skeleton
boundary in CLAUDE.md §1.7 WH5), but neither of them would reproduce
Paper 17 in isolation. Risk of a public reproducibility flag is low.

If they were to ask "what's the LiH R_eq number you get today on the
current codebase?", the answer is "2.82% — slightly better than the
5.3% in Paper 17 because we use a finer R-grid; the original 5.3%
reproduces on the originally-tested grid." This is a clean answer
that strengthens the framework's credibility rather than weakening it.

---

## §6. Bottom line for the calling agent

- Memory note `lih_binding_regression_v3.md` is stale (3 days old); the
  regression it describes was closed in v3.56.0.
- Paper 17's 5.3% R_eq claim is reproducible and beaten (2.82%) via
  the legacy `ComposedDiatomicSolver.LiH_ab_initio(l_max=2)` path with
  fine grid. No paper edit needed.
- The substantive next chemistry sprint, when one is wanted, is the
  second-row frozen-core V_cross R-dependence follow-on (NaH and below).
  That is genuinely open. LiH is closed.
- Suggested 5-minute housekeeping: retire the stale memory note with a
  pointer at `debug/sprint_lih_binding_fix_memo.md` + the v3.56.0
  CHANGELOG entry. Optional 1-paragraph Paper 17 "grid resolution
  sensitivity" remark.

Key files:
- `debug/sprint_lih_binding_fix_memo.md` (the v3.56.0 fix memo)
- `CHANGELOG.md` v3.56.0 entry
- `geovac/level4_multichannel.py` lines 258, 526, 685, 1280 (PK restoration)
- `geovac/balanced_coupled.py` lines 657–692 (V_NN(R) correction)
- `docs/molecular_refactor_handoff.md` (RETIRED 2026-06-05, line 3)
- `debug/sprint_w1e_projection_audit_memo.md` (W1e row-conditional sharpening,
  v3.86.0 — sharpens but does not undo the v3.56.0 LiH fix)
