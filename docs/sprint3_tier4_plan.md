# Sprint 3: Breit Completion + Heavy-Atom Cores + Sunaga Matched-Q

**Version target:** v2.14.0
**Date:** April 2026
**Tracks:** Three parallel tracks (BF is independent of HA/SU; HA and SU are sequential)

---

## Track BF: Breit Completion (Sprint 2 BR-C Follow-Up)

**Goal:** Close the He 2³P fine-structure benchmark to <20% error (from the T8 baseline of 66%) by implementing the full Drake 1971 Breit radial amplitude set and productionizing the module.

**Principle:** Algebraic Deconstruction (every Breit integral should be closed-form)

### Background

Sprint 2 Track BR-C confirmed the angular J-pattern of the He 2³P inversion (2-parameter fit reproduces NIST splittings to 0.000%), but the radial amplitude formulas from ad-hoc Bethe-Salpeter §39 were wrong or incomplete (>600% error). The full treatment requires Drake 1971 (Phys. Rev. A 3, 908) — approximately 10 radial integrals including exchange/cross terms.

Additionally, Sprint 2 found a bug in the BR-B `compute_rk_breit_retarded_algebraic` function: region-splitting skips the m1<0 AND m2<0 case, returning 0 for convergent integrals like (1s,1s;1s,1s) l=2 (published: 83/640).

### Sub-tracks

**BF-A: Fix the BR-B region-splitting bug**

1. Read `debug/br_b_breit_radial.py`, specifically the `_T_kernel_breit_retarded` function and its caller.
2. Identify the branches where m1<0 AND m2<0 are incorrectly skipped.
3. Derive the correct regularization for these cases (should combine into a single convergent integral — the divergent pieces cancel by the transverse tensor structure).
4. Verify against Drake 1971 Table I: (1s,1s;1s,1s) l=2 should give 83/640, etc.
5. Update the test suite (debug-level tests in the BR-B script).

**BF-B: Drake 1971 full radial amplitude set**

1. Read Drake 1971 Section III: the 10-integral set for He-like 2³P states includes:
   - Direct Coulomb F⁰(1s, 2p) and G¹(1s, 2p)
   - Breit SS radial: M²(1s, 2p; 1s, 2p) direct + exchange
   - Breit SOO radial: N¹(1s, 2p) and N¹(2p, 1s)
   - Retarded pieces: M²_ret direct + exchange, N¹_ret

2. For each integral, derive the exact closed form (Bethe-Salpeter rational for same-n, rational + log for cross-n per Sprint 2 BR-B finding).

3. All integrals should be algebraic — if any requires numerical quadrature, flag it as an unexpected structural result.

**BF-C: Production module `geovac/breit_integrals.py`**

Create a clean module with:
- `breit_ss_radial(n_a, l_a, n_b, l_b, Z, k)` — rank-2 spin-spin radial
- `breit_soo_radial(n_a, l_a, n_b, l_b, Z, k)` — rank-1 spin-other-orbit radial
- `breit_retarded(n_a, l_a, n_b, l_b, Z, k)` — retardation corrections
- Exact Fraction arithmetic where possible; sympy Rational + log otherwise
- Unit tests in `tests/test_breit_integrals.py`

**BF-D: He 2³P benchmark with full Breit**

1. Port BR-C's angular J-coefficient machinery to use the new `breit_integrals` module.
2. Compute the three 2³P_J energies via Breit-Pauli perturbation theory on the He 1s2p multiplet.
3. Compare to NIST:
   - 2³P₀−2³P₁: 29616.951 MHz
   - 2³P₁−2³P₂: 2291.178 MHz
   - Full span: 31908.129 MHz
4. Target: <20% error on span (currently 66% with SO only, >600% with ad-hoc BS).

**BF-E: Extend to Li 2²P and Be 2s2p ³P**

If BF-D hits <20% for He:
1. Li 2²P doublet (currently 211% error per T4): 2²P₃/₂ − 2²P₁/₂ = 0.034 cm⁻¹
2. Be 2s2p ³P triplet (currently 78% error per T4): multiplet span ≈ 3.15 cm⁻¹

Both leverage the same `breit_integrals` module.

### Success criteria
- BR-B bug fixed; Drake 1971 Table I values reproduced exactly
- `geovac/breit_integrals.py` with full radial amplitude set and passing tests
- He 2³P <20% error on multiplet span (target met)
- Li and Be benchmarks (bonus if BF-D succeeds)
- Papers 14 §V and 20 Tier-2 table updated with working Breit rows

### Failed approaches to avoid
- Trying to benchmark without fixing BR-B first (will silently get zeros)
- Ad-hoc Bethe-Salpeter §39 minimal form (Sprint 2 confirmed this is under-scoped)

### Files to read first
- `debug/br_b_breit_radial.py` (has the bug)
- `debug/br_c_he_2P_benchmark.py` (angular J-pattern machinery)
- `debug/br_c_he_fine_structure_memo.md` (Sprint 2 scoping diagnosis)
- Drake 1971 if accessible; otherwise Johnson "Atomic Structure Theory" Ch. 8
- `geovac/spin_orbit.py` (closed-form pattern to emulate)

---

## Track HA: Heavy-Atom [Kr]/[Xe] Frozen Cores

**Goal:** Extend `geovac/neon_core.py` FrozenCore to [Kr] (Z=37-54) and [Xe] (Z=55-86) so that heavy-atom hydrides SrH, BaH, RaH can be built.

**Principle:** Natural Geometry Search (incremental) + Algebraic-First (Clementi-Raimondi exponents tabulated, not computed)

### Background

Current FrozenCore coverage:
- [Ne] (Z=11-18): 10-electron core
- [Ar] (Z=19-20, 31-36): 18-electron core + [Ar]3d¹⁰ 28-electron core (v2.2.0)

Missing for Sunaga 2025 comparison:
- [Kr] for SrH (Z=38), YH (Z=39), ..., CdH (Z=48)
- [Kr]4d¹⁰ for InH-XeH (Z=49-54)
- [Xe] for BaH (Z=56), CsH-LuH (Z=57-71)
- [Xe]4f¹⁴ for HfH-HgH (Z=72-80)
- [Xe]4f¹⁴5d¹⁰ for TlH-RnH (Z=81-86), including RaH (Z=88, but Ra uses [Rn] core — out of scope for Sprint 3)

For Sunaga's three published molecules (SrH, BaH, RaH at Q=18): SrH requires [Kr], BaH requires [Xe], RaH requires [Rn]. Sprint 3 targets SrH and BaH; RaH is deferred.

### Sub-tracks

**HA-A: Clementi-Raimondi exponents for [Kr] and [Xe]**

1. Read `geovac/neon_core.py` — understand the existing FrozenCore pattern (tabulated CR exponents, analytical Z_eff(r), density normalization check).

2. Tabulate Clementi-Raimondi Slater orbital exponents for:
   - Kr core: 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p (36 electrons total for Kr itself; 28 for [Kr] core = [Ar]3d¹⁰ shifted + 4s²4p⁶, wait — check: [Kr] = 1s²2s²2p⁶3s²3p⁶3d¹⁰4s²4p⁶ = 36 electrons)
   - Xe core: [Kr]4d¹⁰5s²5p⁶ = 54 electrons

3. Source: Clementi & Raimondi 1963 (J. Chem. Phys. 38, 2686) + Clementi, Raimondi, Reinhardt 1967 (J. Chem. Phys. 47, 1300) for 3rd-4th row.

4. Verify density normalization: integral of n_core(r) over r should equal the core electron count (36 for Kr, 54 for Xe) to within 1%.

**HA-B: Extend atomic_classifier.py for Z=37-54 and Z=55-86**

1. Read `geovac/atomic_classifier.py` — understand the Z=1-36 classification pattern (structure types A/B/C/D/E/F).

2. Add entries for Z=37-54 (fifth row):
   - Z=37-38 (Rb, Sr): alkali/alkaline earth, [Kr] + s valence
   - Z=39-48 (Y through Cd): transition metals, d-block (structure F)
   - Z=49-54 (In-Xe): [Kr]4d¹⁰ + p valence, structure type E

3. Add entries for Z=55-71 (sixth row start):
   - Z=55-56 (Cs, Ba): alkali/alkaline earth, [Xe] + s valence
   - Z=57-71 (La-Lu): f-block + d-block (complex — may defer individual lanthanides)

4. For this sprint: focus on Z=37-38 (Rb, Sr) and Z=55-56 (Cs, Ba) — the two cells needed for SrH and BaH.

**HA-C: SrH and BaH molecular specs**

1. Read `geovac/composed_qubit.py` — see how KH/CaH specs are built (v2.3.0 pattern).

2. Add `srh_spec()` and `bah_spec()` factory functions. Both are single-bond hydrides analogous to KH and CaH₂ (reduced to CaH).

3. Ecosystem export entries for both molecules.

4. Build the composed Hamiltonians and compute:
   - Q (qubit count) at n_max=2
   - N_Pauli
   - 1-norm (total and non-identity)
   - Bond length estimate from PES (if feasible)

**HA-D: Relativistic composed versions**

1. Tier 2 T3 pattern: `srh_spec_relativistic()`, `bah_spec_relativistic()`.
2. Compute Pauli counts with the spinor basis (expect ~2.4× more than scalar per Tier 2 regression).

### Success criteria
- Clementi-Raimondi tabulations for [Kr] and [Xe] with density normalization verified
- atomic_classifier.py entries for Z=37-38 and Z=55-56
- SrH and BaH composed specs with clean Pauli counts and 1-norms
- Relativistic versions of both
- Tests passing; no regression in existing molecule counts

### Failed approaches to avoid
- Trying to solve the 36-electron Kr core quantum-mechanically (infeasible per SCOPE_BOUNDARY)
- Attempting individual lanthanides (f-block open shells) in this sprint

### Files to read first
- `geovac/neon_core.py` (FrozenCore template)
- `geovac/atomic_classifier.py` (classification pattern)
- `geovac/composed_qubit.py` (kh_spec, cah2_spec patterns)
- `geovac/molecular_spec.py` (MolecularSpec, OrbitalBlock)
- SCOPE_BOUNDARY.md (what's currently supported)
- Clementi-Raimondi 1963/1967 papers if accessible

---

## Track SU: Sunaga 2025 Matched-Q=18 Comparison

**Goal:** Execute the head-to-head resource comparison against Sunaga 2025 (PRA 111, 022817) SrH/BaH/RaH at Q=18 that was projected in Tier 2 T4. Gate for promoting Paper 20's Tier-2 table from "native-Q" ratios to "matched-Q" ratios.

**Principle:** Transcendental Cataloging (compare where our transcendental content is cheaper than theirs)

**Dependency:** Requires HA-C complete (SrH and BaH specs must exist).

### Background

Tier 2 Track T4 projected:
- GeoVac LiH rel/Sunaga RaH-18q ratio: 0.017× (native Q=30 for LiH, but not matched)
- GeoVac CaH rel/Sunaga RaH-18q ratio: 0.011× (native Q=30, not matched)

Sunaga's published values (from PRA 111, 022817 main paper Table II): only the RaH-18q cell (47,099 Pauli terms) is public; SrH/BaH/RaH at 18q are in SI Tables S1-S3, flagged DEFERRED.

**For matched-Q=18:** GeoVac's native Q for SrH scalar is ~30 (1 core block + 1 bond pair at n_max=2). To get to Q=18, we would need to truncate the basis — which is not a standard GeoVac operating point. This track assesses whether a matched-Q comparison is even meaningful.

### Sub-tracks

**SU-A: Obtain Sunaga SI data**

1. Check whether Sunaga 2025 SI is accessible (PRA supplemental materials).
2. If accessible: extract SrH/BaH/RaH Q=18 Pauli counts, 1-norm, QWC groups from Tables S1-S3.
3. If not accessible: proceed with SU-B and flag the comparison as "published RaH-18q only"

**SU-B: Native-Q GeoVac resource table for SrH/BaH**

1. At native n_max=2 (Q=30 for single-bond hydrides):
   - N_Pauli scalar, relativistic
   - 1-norm total and non-identity
   - QWC group count
2. Compare to Sunaga RaH-18q = 47,099 Pauli as a single reference point (Tier 2 T4 pattern).

**SU-C: Matched-Q feasibility assessment**

Can GeoVac produce a Q=18 Hamiltonian for SrH? Options:
1. **Truncate basis:** Drop some orbitals from the default block structure. Clean but may not represent the molecule faithfully.
2. **Increase block density:** Use n_max=1 for a minimal basis (Q=10, too small).
3. **Custom spec:** Build a SrH spec with explicitly 9 spatial orbitals (18 spin-orbitals) — requires justification for which 9 orbitals.

Document the structural reason GeoVac operates at native-Q rather than matched-Q. This is itself a valuable contribution to the quantum resource comparison discussion.

**SU-D: Paper 20 Tier-2 table update**

1. Add SrH and BaH rows to Paper 20 §V Tier-2 resource table.
2. Update Sunaga comparison section with:
   - Native-Q ratios (definitive)
   - Matched-Q feasibility note (structural)
   - If SU-A succeeded: actual matched-Q ratios for SrH/BaH

3. Gate on Paper 14 §V updates if resource wins are clear.

### Success criteria
- SrH and BaH native-Q resource tables (Pauli, 1-norm, QWC)
- Matched-Q feasibility assessment (structural note)
- Paper 20 Tier-2 table updated with at least native-Q rows

### Failed approaches to avoid
- Artificially truncating the basis to hit Q=18 without physical justification
- Over-claiming the Sunaga comparison beyond what the data supports

### Files to read first
- `geovac/composed_qubit_relativistic.py` (Tier 2 T3 builder)
- `geovac/ecosystem_export.py` (export pipeline)
- `docs/tier2_market_test.md` (Tier 2 T4 Sunaga comparison methodology)
- `papers/applications/paper_20_resource_benchmarks.tex` (Tier-2 table)

---

## Sprint 3 PM Prompt

```
Read CLAUDE.md, docs/sprint3_tier4_plan.md, and docs/sprint2_final_summary.md
(for context on Sprint 2 results).

Dispatch Track BF and Track HA as parallel sub-agents. Track SU depends on
HA-C, so it waits for the SrH/BaH specs to exist.

Track BF (Breit Completion):
  Sub-agent 1 (BF-A+B+C): Fix BR-B region-splitting bug, implement Drake 1971
    full radial amplitude set, create geovac/breit_integrals.py production
    module with tests. Algebraic-first: every integral must be closed form.

  Sub-agent 2 (BF-D+E): Re-run He 2³P benchmark with fixed radial amplitudes;
    if <20% target met, extend to Li and Be. Update Papers 14 §V and 20
    Tier-2 table with working Breit rows.

Track HA (Heavy Atoms):
  Sub-agent 3 (HA-A+B): Tabulate Clementi-Raimondi exponents for [Kr] and
    [Xe], extend atomic_classifier.py for Z=37-38 and Z=55-56, extend
    FrozenCore in neon_core.py (or create kr_core.py / xe_core.py).

  Sub-agent 4 (HA-C+D): Build srh_spec, bah_spec, and their relativistic
    counterparts. Wire into ecosystem_export. Add tests.

Track SU (Sunaga comparison): waits for HA-C completion.
  Sub-agent 5 (SU-A+B+C+D): Attempt to obtain Sunaga SI data; compute
    native-Q GeoVac resources for SrH and BaH; assess matched-Q feasibility;
    update Paper 20.

Algebraic-first: BF integrals should be exact sympy rationals or rational+log.
HA exponents should be tabulated values (not computed). SU comparison should
be exact counts, no projections beyond stated data.

Exit criteria:
- BF: He 2³P error <20% (or honest-negative if Drake 1971 implementation
  doesn't close the gap); Papers 14/20 updated
- HA: SrH and BaH specs in place with clean Pauli/1-norm; tests passing
- SU: Paper 20 Tier-2 table has SrH/BaH rows; matched-Q feasibility
  documented
```

---

## Sprint 4 Dependencies

Sprint 4 candidates (from the Leader brief, not yet selected):
- QED-on-S³ (Hopf bundle as lattice gauge structure) — paper-worthy observation but not a new computation
- Tier-2 T3 regression fix (spinor/scalar α=0 mismatch found in DC-B §5)
- f-block (lanthanide) atomic classification if Sprint 3 opens the pattern
- Sunaga RaH matched benchmark if [Rn] core implementation becomes feasible

---

## What Sprint 3 Resolves

1. **Fine-structure honest gap closed** (BF target): He 2³P <20% error via Drake 1971. Paper 14 §V table gains a working Breit row.

2. **Heavy-atom pipeline unblocked** (HA target): [Kr] and [Xe] frozen cores tabulated. SrH and BaH enter the molecule library. Sunaga 2025 matched comparison becomes at least partially possible.

3. **Sunaga comparison executed** (SU target): Paper 20 Tier-2 table upgrades from "projected" to "measured" for at least native-Q. If Sunaga SI is accessible, matched-Q too.

All three tracks are in the algebraic-first / tabulated-constants regime with clear pass/fail criteria.
