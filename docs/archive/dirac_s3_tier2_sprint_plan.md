# Dirac-on-S³ Tier 2 Sprint Plan

**PI approved:** 2026-04-15. Reshaped from Leader proposal by Explorer findings (Szmytkowski angular tables, Martínez-y-Romero radial recursions, unpublished d(l_max) for spinors, OpenFermion-Dirac interoperability).
**Target:** 1-3 weeks under algebraic-first discipline (Tier 1 precedent: hours of compute for "weeks of work").
**Framing:** curiosity-driven — "how does the composed-sparsity mechanism extend when the basis carries spin, and where does relativistic physics enter the exchange-constant taxonomy?" Publishable across all outcomes.

---

## Governing philosophy

Algebraic-first. All worker dispatches carry the Tier 1 preamble verbatim. Reference: Paper 18 taxonomy.

**Honest framing constraint (Explorer Gap #5):** Tier 2 builds spin-orbit corrections on the existing S³ scalar graph using Szmytkowski angular tables and Martínez-y-Romero Dirac-Coulomb radial recursions. This is NOT a new Dirac Fock projection onto S³ (BJL Z₂-supersymmetric algebra, Tier 1 Explorer Finding 2.1). Any paper output must state this clearly.

---

## PI decisions (recorded)

1. **Labeling:** (κ, m_j) per Explorer Finding T2-1.1. Conversion to/from D1's (l, σ, m_j) is mechanical.
2. **Scope:** Three molecules — **LiH, BeH, SrH** — matching Sunaga 2025 for point-by-point market comparison.
3. **Sector:** Weyl (2-component). Full Dirac positronic sector deferred past Tier 3.
4. **SO coupling:** direct (not perturbative X2C). Martínez-y-Romero + Szmytkowski make direct algebraic.
5. **V_ee:** Coulomb only. Breit deferred to Tier 3.
6. **Paper 18 extension:** new subtier "spinor-intrinsic content" for α² fine-structure factors and √(1-(Zα)²) coefficient ring.
7. **d(l_max) deliverable:** folded into Paper 22 as a new section, not standalone.
8. **Market baseline:** Sunaga et al. 2025 (PRA 111, 022817) via OpenFermion-Dirac. DIRAC for accuracy cross-check only.

---

## Proposed tracks (final)

### Track T0 — Spinor-block ERI density d(l_max) (EARLY-WIN, standalone deliverable)

**Goal:** Compute the Paper 22 analog for spinor blocks. Spinor ERI density depends only on l_max (or κ_max), not on V(r), and can be computed from Szmytkowski angular tables alone — no radial integrals, no composed builder needed.

**Algebraic-first constraint:** 6j-symbol triangle conditions (Dyall §9) give sparsity in closed form. Use sympy Wigner 6j symbols throughout; no numerical tensors.

**Predicted outcome:** d_spinor(l_max) ≤ d_scalar(l_max) at every l_max, with equality only at l_max=0. The extra triangle rules (two 6j's) can add zeros but never remove them.

**Success:** d(l_max) computed for l_max = 0, 1, 2, 3, 4 (matching Paper 22's scalar 1.44%/0.90%/0.62% curve). Presented alongside scalar values in a comparative table.

**Failure:** Sparsity density is qualitatively worse than scalar at any l_max. Publishable as a structural finding — "spin doubling destroys angular sparsity."

**Deliverable:** `debug/tier2_t0_spinor_density.py` + `debug/dirac_t0_memo.md` + preliminary data for Paper 22 extension.

**Depends on:** D1 (done). Szmytkowski angular table (Tier 2 Explorer Finding T2-1.1).

**Expected runtime:** 2-3 days.

---

### Track T1 — Spinor matrix elements in (κ, m_j) basis

**Goal:** Build `geovac/dirac_matrix_elements.py` exposing ⟨n', κ', m_j' | r^k | n, κ, m_j⟩, ⟨1/r⟩, ⟨∇⟩ via Szmytkowski angular + Martínez-y-Romero radial recursions. All closed form, exact Fraction arithmetic where possible, with the single transcendental seed being α² (or √(1-(Zα)²), depending on which factorization the Decomposer prefers).

**Algebraic-first constraint:** Do NOT numerically diagonalize radial Dirac operators. Martínez-y-Romero gives three-term recursions. If a matrix element looks like it needs integration, stop and derive the recursion relation first.

**Labeling:** (κ, m_j). Update `SpinorHarmonicLabel` in `geovac/dirac_s3.py` to carry a `j` field (or convert to κ-native labels). D1 §8 flagged this decision.

**Success:** All three matrix element families available in closed form. Machine-precision match to scalar limit when α → 0 and σ → 0. Labeling bridge to D1 verified by regression tests.

**Failure:** Martínez-y-Romero recurrences don't close for some class (e.g., off-diagonal in κ with gradient). Fall-back: identify that class as numerical-required per Paper 18 taxonomy and document.

**Deliverable:** `geovac/dirac_matrix_elements.py`; `tests/test_dirac_matrix_elements.py`; extension to `geovac/dirac_s3.py` for (κ, m_j) labeling; memo `debug/dirac_t1_memo.md`.

**Depends on:** D1, T0 (both for labeling decision and to confirm Szmytkowski angular infrastructure works).

**Guardrails:**
- Do NOT attempt a Dirac Fock projection theorem (Tier 1 Explorer Gap #1).
- Do NOT combine with TC (CUSP-3, decisively dead).

**Expected runtime:** 3-5 days.

---

### Track T2 — Spin-orbit coupling

**Goal:** Compute ⟨ξ(r) L·S⟩ matrix elements where ξ(r) = (1/r) dV/dr = Z/r³ for Coulomb. Use the closed-form ⟨1/r³⟩ hydrogenic result and the L·S diagonal eigenvalue [j(j+1) − l(l+1) − 3/4]/2. Spin-orbit operator in (κ, m_j) is diagonal in (n, κ), so this is purely a block-diagonal additive term.

**Algebraic-first constraint:** ⟨1/r³⟩ is a standard closed form. No numerical integration.

**Success:** SO Hamiltonian builds symbolically; matrix elements are Q(α²) rationals modulo ⟨1/r³⟩ seed; Z-scaling ⟨ξ⟩ ∝ Z⁴ verified across Z={3, 4, 38}.

**Failure:** SO operator doesn't commute with the composed block structure (cross-block couplings appear). Publishable as structural limitation.

**Deliverable:** `geovac/spin_orbit.py`; `tests/test_spin_orbit.py`; `debug/dirac_t2_memo.md`.

**Depends on:** T1.

**Expected runtime:** 1-2 days.

---

### Track T3 — Spin-ful composed pipeline integration

**Goal:** Extend `MolecularSpec` / `build_composed_hamiltonian` with `relativistic=True` kwarg. Wire T1 matrix elements + T2 SO coupling. Build LiH, BeH, SrH spin-ful composed Hamiltonians. Extend Gaunt (3j) to Dyall §9 6j recoupling for the two-electron blocks.

**Algebraic-first constraint:** 6j selection rules are symbolic (sympy `wigner_6j`). No numerical ERI tensor pruning.

**BeH / SrH spec bootstrap:** BeH = BeH₂ spec with one H removed (or a new `beh_spec()` derived from the `beh2_spec` template). SrH = new spec via Z=38 [Ar]3d¹⁰ frozen core — same pattern as the Z=19,20,31-36 frozen-core work in v2.2.0.

**Success:** Three composed Hamiltonians build; Pauli count, QWC groups, 1-norm computed; all existing spinless tests unchanged (regression-free).

**Failure:** 6j recoupling destroys block-diagonal ERI structure (Pauli count inflates > 4× vs spinless). Publishable as structural finding.

**Deliverable:** Extensions to `geovac/composed_qubit.py`, `geovac/molecular_spec.py`, `geovac/ecosystem_export.py`; new specs `beh_spec()`, `srh_spec()`; regression tests.

**Depends on:** T1, T2.

**Expected runtime:** 3-5 days. Largest engineering track.

---

### Track T4 — Market test + fine-structure sanity check

**Goal:** Point-by-point comparison against Sunaga et al. 2025 (PRA 111, 022817) for LiH, BeH, SrH. Report Pauli count, QWC groups, 1-norm, qubit count at equal n_max and at matched-qubit points. Use Senjean's Openfermion-Dirac for the Sunaga baseline.

Secondary: compute fine-structure splittings for He (2³P), Li (2²P doublet), Be (2s2p ³P) at n_max=2,3 and report sign + order-of-magnitude accuracy. Target 20-50% relative, NOT kHz absolute (Explorer T2-3 reality check).

**Algebraic-first constraint:** Analysis only. All computation upstream of T4 is symbolic; T4 tabulates.

**Success:** GeoVac spinor-composed Pauli count ≤ 0.5× Sunaga baseline at matched Q for at least one molecule. Fine-structure splittings have correct sign and correct order of magnitude across all three atoms.

**Failure:** GeoVac Pauli count exceeds Sunaga — document as structural limitation. Splittings wrong sign — indicates labeling or SO sign error in T1/T2.

**Deliverable:** `benchmarks/relativistic_comparison.py`; `docs/tier2_market_test.md`; tables for Paper 14 §V.

**Depends on:** T3.

**Expected runtime:** 2-3 days.

**Guardrail:** Openfermion-Dirac requires compiling DIRAC with Senjean's patch. If that's infeasible in the sprint timeframe, substitute published Sunaga numbers (PRA 111, 022817 Tables) and note that direct runtime comparison is deferred.

---

### Track T5 — π-free certificate extension with new Paper 18 subtier

**Goal:** Extend `verify_pi_free` from Tier 1 D1 to certify that T3's spin-ful composed Hamiltonian coefficients decompose as:
  Q(α²) × {algebraic: √(1-(Zα)²)} × {seed: ⟨r^k⟩, ⟨1/r³⟩}
with no hidden transcendentals. Classify the α² and √(1-(Zα)²) content as a new Paper 18 subtier "spinor-intrinsic content" — first-order-operator analog of Paper 24's "second-order-operator calibration π" result.

**Algebraic-first constraint:** Taxonomic classification, no computation beyond T3's output.

**Success:** Every T3 Hamiltonian coefficient traces to Paper 18 taxonomy. New subtier cleanly added.

**Failure:** Some coefficient carries an unexpected π or ζ(3) that isn't in the taxonomy. Flag and extend the taxonomy further (or identify a bug in T1/T2).

**Deliverable:** Extension to `tests/test_dirac_s3.py`; proposed Paper 18 §IV update (drop-in file, not in-place edit); memo `debug/dirac_t5_memo.md`.

**Depends on:** T3.

**Expected runtime:** 1-2 days.

---

### Track T6 — Verdict memo + paper extensions

**Goal:** Synthesis memo. Drop-in proposals (NOT in-place edits):
- Paper 14 §V "Spinor Composed Encoding" — Pauli scaling law, 1-norm, QWC groups for LiH/BeH/SrH; head-to-head Sunaga comparison.
- Paper 20 new resource-estimation table — same three molecules.
- Paper 22 new section — d_spinor(l_max) alongside the existing scalar d(l_max) theorem.
- Paper 18 §IV — new spinor-intrinsic subtier from T5.
- CLAUDE.md — mechanical edits (version bump v2.10.0 → v2.11.0, §2 Phase 4J summary, §3 any negatives found, §10 new benchmark rows, §11 topic mappings).

**Success:** All four paper proposals + CLAUDE.md edits written as drop-ins, ready for PI review.

**Failure:** — this is synthesis, no failure mode.

**Deliverable:** `docs/tier2_verdict.md`; `docs/paper14_section5_proposal.tex`; `docs/paper22_spinor_section_proposal.tex`; `docs/paper18_spinor_subtier_proposal.tex`; `docs/paper20_tier2_table_proposal.tex`; `docs/claude_md_tier2_updates.md`.

**Depends on:** T4, T5.

**Expected runtime:** 1-2 days.

---

## Sequencing

```
Day 0:    PI approved; sprint launched.
Days 1-3: T0 (early-win standalone deliverable).
Days 3-7: T1 (spinor matrix elements).
Days 5-7: T2 (SO coupling; parallel with T1 tail after labeling is fixed).
Days 7-12: T3 (composed integration, largest track).
Days 11-14: T4 + T5 (parallel).
Days 13-15: T6 (synthesis).
Day 16: PI review. Commit. Tier 3 scoping.
```

**Realistic expectation:** 1-2 weeks if algebraic structures fall cleanly (T0 and T2 are essentially mechanical; T3 is the unknown pole).

---

## Guardrails (sprint-wide)

From Tier 1 Explorer + Tier 2 Explorer:

- **No Dirac Fock projection claims.** BJL Z₂-supersymmetric algebra doesn't lift the Schrödinger Fock map. All Tier 2 work is SO on scalar S³ graph.
- **No S⁵/S⁷ work.** Closed by Phase 4E/4G and Tier 1 Explorer Finding 3.1.
- **No discrete graph Dirac operator.** Ginsparg-Wilson theorem forbids local chirality-preserving graph Dirac.
- **No TC combinations.** CUSP-3 (Tier 1 adjacent) confirmed TC dead on composed basis at all n_max.
- **No KRCI architecture import.** KRCI is a string-based full CI driver, not composed-architecture-compatible (Tier 2 Explorer T2-2.1 Section-3-adjacency flag).
- **No absolute NIST-level fine-structure accuracy targets.** Tier 2 Explorer T2-3 reality check: target sign + order of magnitude, not sub-kHz matching.
- **No multi-center molecules beyond LiH/BeH/SrH.** Further multi-center extension is Tier 3+.

---

## Papers affected

- **Paper 14** — §V "Spinor Composed Encoding" (new section, proposal-only until PI review).
- **Paper 18** — §IV new spinor-intrinsic subtier (proposal-only).
- **Paper 20** — new resource table for LiH/BeH/SrH (proposal-only).
- **Paper 22** — new section on d_spinor(l_max) (proposal-only).
- **Paper 2** — NOT modified (Tier 1 closed that story).

---

## Memory references

- `feedback_algebraic_first.md` — governing philosophy (universal).
- `project_dirac_s3_pivot.md` — three-tier research direction (Tier 2 is the middle tier).
- `feedback_tc_correction.md` — particle-number-projected FCI.
- `feedback_approach.md` — curiosity-driven framing.

---

## PM prompt preamble (verbatim required in all worker dispatches)

```
GOVERNING PHILOSOPHY (non-negotiable):
GeoVac has repeatedly found that continuum problems admit algebraic solutions
— Paper 12 Neumann expansion, Gaunt integrals, Level 2/3 Laguerre recurrence,
hypergeometric_slater.py, Paper 24 Bargmann-Segal, Tier 1 Dirac-on-S³ sprint
(hours of compute for the "3-4 week" sprint) — and the algebraic route is
consistently faster, more accurate, and more structurally informative than
numerical methods. BEFORE committing to any numerical computation taking more
than a few minutes, you MUST:
  1. Identify the transcendental content per Paper 18's exchange-constant
     taxonomy (intrinsic / calibration / embedding / flow / spinor-intrinsic).
  2. Check whether quantum-number recursion, group-theoretic selection rules,
     or exact Fraction arithmetic gives a closed form.
  3. If the algebraic route is unclear, stop and report back rather than
     running a long numerical job.

Numerical output in floating point is not a result. A sympy symbolic equality
IS a result. A "match to 10⁻⁸" requires an exact-arithmetic follow-up or it
doesn't count.

Reference: Paper 18 (exchange constants taxonomy).

TIER 2 FRAMING: Dirac-on-S³ provides spinor labels (D1). Tier 2 builds
spin-orbit corrections on the existing scalar S³ graph using Szmytkowski
angular tables and Martínez-y-Romero Dirac-Coulomb radial recursions. This is
NOT a new Dirac Fock projection — BJL algebra (Tier 1 Explorer Finding 2.1)
does not lift the Schrödinger Fock map.
```
