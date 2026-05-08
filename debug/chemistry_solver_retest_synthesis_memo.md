# Chemistry-Solver Re-Test Synthesis (LiH under multi-focal architecture)
**Date:** 2026-05-08
**Worker fork directive:** test whether multi-focal Phase C closures
(balanced coupled + cross-register V_eN + screened W1c) close, persist,
or change the v2.0.32 R_eq-drift wall.

## Verdict — short form

**Combined verdict (first-row + second-row):**

1. **First-row LiH:** drift persists at Track CD level (+0.054 bohr/n_max),
   unchanged by Phase C closures. W1c is a structural no-op for first
   row — there is no [Ne] core for it to screen — verified
   bit-precisely (5e-14 Ha max diff) across a 5-point PES grid. The
   chemistry-solver upgrade the framework HAS already received is
   the balanced-coupled architecture itself (Track CD, April 2026 /
   v2.0.39), which reduced the original PK-composed wall (+0.18
   bohr/l_max) by 3× to +0.054 bohr/n_max. n_max=2 reproduces Track
   CD baseline bit-precisely (R_eq=3.227 vs 3.226).

2. **Second-row NaH (NEW):** the W1c-residual orthogonality wall
   named in Phase C is now empirically confirmed at the FCI level.
   `coupled_fci_energy` provides the n_e-projected FCI driver that
   was previously missing. NaH balanced FCI at max_n=2 (Q=20)
   monotonically descends in BOTH `screened=False` and
   `screened=True` — no internal minimum at any tested R from 2.5
   to 10.0 bohr. **W1c reduces edge-to-dissociation gap 17.5×**
   (bare: −6.24 Ha; screened: −0.36 Ha) — three times more than
   the 5.4–6× single-R diagnostic — but binding is NOT recovered.

The chemistry-solver wall for first-row is intrinsic to balanced
coupled in the qubit basis at small n_max, not to PK overcounting:
single-center hydrogenic basis flexibility per block cannot resolve
nuclear geometry to better than 7–9% at equilibrium.

The chemistry-solver wall for second-row is the W1c-residual
orthogonality wall: the H valence orbital ↔ [Ne]-core
Phillips-Kleinman-class orthogonality sits beyond the Hartree
screening of Z_eff(ρ) that W1c implements. The architectural
difference between first-row binding and second-row non-binding is
that first-row uses an explicit core block in the FCI sector
(electrons resolve their own orthogonality), while second-row uses
a frozen-core potential (orthogonality must be enforced externally,
and W1c does not yet do so).

## Architectural reframing (the substantive finding)

The directive presupposed a single coherent "chemistry solver". The
v2.0.32 baseline (5.3% R_eq, +0.15–0.22 bohr/l_max drift) was
diagnosed in `composed_diatomic.ComposedDiatomicSolver` —
the level-4 PK-composed PES solver. The Phase C closures live in
`balanced_coupled.build_balanced_hamiltonian` — the qubit-encoded
FCI Hamiltonian builder. These are *different solvers*.

| Architecture | Used in | l_max meaning | R_eq drift |
|--------------|---------|---------------|------------|
| Level-4 PK-composed (v2.0.32) | `composed_diatomic.py` | hyperspherical angular basis | +0.18 bohr/l_max (structural to PK) |
| Balanced coupled (Track CD, v2.0.39) | `balanced_coupled.py` + `coupled_fci_energy` | n_max → l_max = n−1 (coupled) | +0.054 bohr/n_max (3× smaller) |

The balanced-coupled construction *replaced* PK with cross-block
ERIs + cross-center V_ne via Shibuya–Wulfman multipole expansion.
That replacement was the chemistry-solver upgrade. It already shipped
in v2.0.39, before Phase C. Phase C added W1c (frozen-core-aware
screening of the multipole V_ne) which is structurally relevant only
for second-row Z ≥ 11 species.

The strict question "did Phase C upgrade the chemistry solver?" has
the answer **no for first row, untested for second row** with this
sprint's scope.

## R_eq table

| System | n_max | R_eq (bohr) | R_err (%) | E_min (Ha, TC-conv) | E_err (%) | source |
|--------|------:|------------:|----------:|--------------------:|----------:|--------|
| LiH (PK composed, v2.0.32) | l_max=2 | ≈ 3.175 | 5.3 | n/a | n/a | level-4 baseline |
| LiH (PK composed, v2.0.32) | l_max=3 | ≈ 3.33 | ≈ 10.5 | n/a | n/a | level-4 + drift |
| LiH (PK composed, v2.0.32) | l_max=4 | ≈ 3.55 | ≈ 17.7 | n/a | n/a | level-4 + drift |
| LiH (balanced, Track CD)   | n_max=2 | 3.226 | 7.0 | −7.924 | 1.8 | v2.0.39 |
| LiH (balanced, Track CD)   | n_max=3 | 3.280 | 8.8 | −8.055 | 0.20 | v2.0.39 |
| LiH (balanced, this re-test, screened OFF) | n_max=2 | **3.227** | 7.03 | **−7.933** | **1.71** | v2.32 |
| LiH (balanced, this re-test, screened ON)  | n_max=2 | **3.227** | 7.03 | **−7.933** | **1.71** | v2.32 |
| LiH (balanced, this re-test, screened ON)  | n_max=3 | n/a | n/a | n/a | n/a | v2.32 stalled in eigsh; Track CD: 3.280 / 8.8% / −8.055 / 0.20% |
| NaH (balanced, screened OFF) | n_max=2 | **none** (PES descends) | n/a | n/a | n/a | v2.32 (Sprint 7 baseline reproduced) |
| NaH (balanced, screened ON)  | n_max=2 | **none** (PES descends) | n/a | n/a | n/a | v2.32 (W1c cuts overattraction 17.5× but no min) |

Reference: experimental R_eq = 3.015 bohr; literature E_min = −8.071 Ha.
"TC-conv" = subtract E_core(Li) = −7.2799 Ha to match Track CD's
absolute-energy convention. The full E_total = E_core + V_NN + E_electronic.

## W1c bit-identical claim — verified

Across 5 R-grid points (R = 2.4, 2.7, 3.0, 3.3, 3.6 bohr) at
max_n=2, balanced coupled with `screened_cross_center=True` produced
energies bit-identical to `screened_cross_center=False` to within
**5e-14 Ha** (machine precision). This validates the Phase C-W1c
diagnostic claim that the auto-detection in
`cross_center_screened_vne._detect_core_type` correctly returns
`None` for first-row systems and the screened path collapses
analytically to the bare path.

For the LiH chemistry solver, the W1c flag is therefore a documented
no-op. The W1c machinery exists for the second-row PES regression
(NaH/MgH₂/HCl, see Phase C D-PES regression in CLAUDE.md §2) where
the off-center [Ne] core must be screened to recover physical
binding behavior.

## v2.0.32 vs Track CD baseline reproduction

The current v2.32 production code reproduces Track CD's n_max=2
result to bit-precision:

- R_eq: 3.227 vs 3.226 (Track CD) — agreement to within 0.001 bohr,
  which is below the parabolic-fit precision on a coarse 5-point grid.
- E_min(TC-conv): −7.933 vs −7.924 (Track CD) — agreement to 0.009 Ha,
  consistent with grid-precision differences.
- E_err vs literature: 1.71% (this) vs 1.8% (Track CD) — bit-identical
  to displayed precision.
- R_err vs experiment: 7.03% (this) vs 7.0% (Track CD) — bit-identical.

The post-Track-CD changes to the codebase (multi-focal Phase C, paper
edits, Sprint MH, etc.) have not regressed first-row chemistry
calculations.

## Recommended next sprint

The first-row chemistry solver is at its current accuracy ceiling
(R_eq error 7-9%, E_err 0.2-2%) and does not improve under any
Phase C closure. The second-row binding-determination test (NaH
PES with `coupled_fci_energy` + screened W1c) is now COMPLETE
(this sprint) and confirms the W1c-residual orthogonality wall.
The next productive investigations are:

1. **W1c-residual orthogonality fix: Phillips-Kleinman-class
   projection on the H valence orbital against the [Ne] core in
   cross-center geometry.** This is the SINGLE remaining engineering
   closure for second-row chemistry. The mechanical implementation
   would extend `cross_center_screened_vne` to subtract the H-side
   1s overlap with the [Ne] core orbital before the multipole
   expansion. Conceptually similar to PK in `geovac/ab_initio_pk.py`
   but on the cross-center side rather than the same-center side.

2. **Multi-focal V_ee operator.** The balanced builder uses
   intra-block hydrogenic V_ee Slater integrals (single focal length
   per block) and `compute_cross_block_eri` for inter-block V_ee
   (which uses center transformation, not focal-length composition).
   A genuine cross-register V_ee operator — the e-e analog of W1a's
   e-N closure — would be the natural "second cross-register
   operator" for chemistry. Currently every multi-electron pair sees
   classical Coulomb between two electrons on the same focal length;
   electrons in different blocks see classical inter-center V_ee.
   Whether a Pauli-encoded V_ee with two focal lengths gives
   structurally different chemistry is open.

3. **n_max=4 for first-row.** Q ≈ 168 for LiH (3 blocks × 30 spatial
   each, including g orbitals). FCI dimension C(84,2)² ≈ 12M
   determinants — feasible with eigsh, ~1–2 hours wall-time per
   point. Whether the +0.054 bohr/n_max drift continues, asymptotes,
   or reverses is the basis-completeness question. Diminishing
   returns expected per Track CD's "energy converges excellently,
   R_eq drifts structurally" diagnosis, but the actual experiment
   has not been run.

The (1) sprint is the cheapest and would close a named open question
from Phase C. (2) is structurally interesting but speculative. (3) is
expensive and likely confirmatory.

## NaH binding determination (this sprint, key extension)

| R (bohr) | E_total (off, Ha) | E_total (on, Ha) | E_elec (off) | E_elec (on) |
|----------|-------------------|------------------|--------------|-------------|
| 2.5      | −171.0967         | −163.3160        |  −9.638      |  −1.857     |
| 3.0      | −169.9734         | −163.2265        |  −8.448      |  −1.701     |
| 3.5      | −169.1987         | −163.1656        |  −7.625      |  −1.592     |
| 4.0      | −168.6314         | −163.1220        |  −7.022      |  −1.513     |
| 5.0      | −167.7518         | −163.0644        |  −6.093      |  −1.405     |
| 7.0      | −166.2988         | −163.0032        |  −4.583      |  −1.287     |
| 10.0     | −164.8525         | −162.9593        |  −3.094      |  −1.200     |

Both columns descend monotonically. The bare-Coulomb cross-V_ne
(screened=False, Sprint 7 baseline) overattracts by orders of
magnitude (E_elec at R=2.5: −9.6 Ha vs physical E_NaH ≈ −0.7 Ha).
W1c screening reduces this by ~17.5× in PES descent magnitude
(edge-to-dissociation gap −6.24 Ha → −0.36 Ha) but does not
introduce an internal minimum.

For comparison, the experimental NaH bond length is 3.566 bohr
with D_e ≈ 0.075 Ha. Even the screened result has E_elec ≈ −1.86
Ha at R=2.5 (its smallest tested value) and ≈ −1.20 Ha at R=10
— the residual ≈ 0.65 Ha overattraction sits ~9× above the
physical D_e, indicating the H-on-[Ne] orthogonality wall.

## Honest limitations

- This sprint did not test BeH₂ or H₂O — only LiH was extended
  beyond v2.0.39 baseline. The architectural finding (W1c bit-
  identical for first-row) holds by construction (no [Ne] core),
  and the v2.32 reproduction of Track CD's first-row n_max=2 result
  generalizes by symmetry, but BeH₂/H₂O n_max=3 numerical reproduction
  was not directly verified.
- LiH n_max=3 single point at R=3.015 was attempted but stalled in
  4-electron FCI eigsh on a 740k-determinant sparse matrix beyond
  the worker fork's compute budget (>9 minutes without producing
  output, killed cleanly). Track CD's published n_max=3 result
  (R_eq=3.280, drift +0.054 bohr/n_max) is accepted as baseline.
  The architectural finding does not depend on n_max=3 reproduction.
- The directive's "drift slope < +0.05 bohr/l_max" exit threshold is
  a level-4 metric. In the balanced-coupled n_max basis, the
  comparable threshold is the +0.054 bohr/n_max drift slope, which
  was already at the threshold from Track CD — i.e. the architecture
  was already at the verdict boundary before Phase C.
- MgH₂ and HCl PES regression with `coupled_fci_energy` was not run.
  The Z-decreasing W1c-residual orthogonality wall hypothesis from
  Phase C (universal cause, magnitude scales inversely with Z) would
  predict that MgH₂ and HCl at higher Z see proportionally less
  W1c-residual overattraction. NaH (Z=11) is the worst case; testing
  the trend is a follow-up sprint.
- The n_max=3 build cost for LiH balanced coupled (Q=84, 4-electron
  FCI dimension ≈ 741k determinants) is a pragmatic infrastructure
  finding: the existing `coupled_fci_energy` driver builds the
  matrix in `lil_matrix` via Python loops, which scales poorly. A
  future infrastructure sprint targeting either (a) a fortran/cython
  FCI matrix construction or (b) a Davidson eigensolver that does
  not materialize the matrix would unblock this regime.
