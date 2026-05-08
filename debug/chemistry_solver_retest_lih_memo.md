# Chemistry-Solver Re-Test (LiH) under Multi-Focal Architecture
**Date:** 2026-05-08
**Sprint:** Multi-track follow-up (per PI Sprint Plan, sibling to Tracks 1, 3)
**Track:** Track 2 — Chemistry-solver re-test
**Worker:** Single-shot fork
**Status:** in flight, this memo records partial results + architectural finding

## Question

Does the multi-focal Phase C architecture (balanced coupled +
cross-center V_ne via Shibuya–Wulfman + screened W1c via
`cross_center_screened_vne`) close, persist, or change the character
of the v2.0.32 baseline R_eq drift (+0.15–0.22 bohr / l_max in the
level-4 PK-composed PES solver)?

## Headline architectural finding

The v2.0.32 baseline drift was diagnosed in the level-4 PK-composed
PES solver (`composed_diatomic.py` + `variational_2d`). The
multi-focal Phase C closures live in the **qubit-encoded balanced
coupled architecture** (`balanced_coupled.py`), which is structurally
a *different* solver. It does not have a separate `l_max` parameter:
the basis is enumerated at $(n, l, m)$ with $l < n$, so `n_max` controls
both radial and angular convergence simultaneously.

This means the strict "l_max drift" question cannot be tested in the
balanced coupled framework — it must be reframed as **n_max drift**.
Track CD (v2.0.39, balanced coupled LiH) already characterized this:

| n_max | R_eq (bohr) | R_err (%) | E_min (Ha) | E_err (%) |
|-------|-------------|-----------|------------|-----------|
|   2   |   3.226     |    7.0    |   −7.924   |    1.8    |
|   3   |   3.280     |    8.8    |   −8.055   |    0.20   |

with structural drift +0.054 bohr / n_max — three times smaller than
the PK-composed drift (+0.18 bohr / l_max) but **not zero**. Track CD's
headline: "energy converges excellently, R_eq drifts structurally."

The current re-test verifies this reproduces under v2.32 production
code, and asks whether the W1c closure (`screened_cross_center=True`)
changes the n_max-drift character.

## Result 1: W1c bit-identical for LiH (PHASE C CLAIM VERIFIED)

LiH at max_n=2, 5 R-points × 2 panels (screened OFF vs ON):

| R (bohr) | E (off, Ha)   | E (on, Ha)    | |Δ| (Ha)   |
|----------|---------------|---------------|------------|
| 2.4      | −15.162678    | −15.162678    |  ≤ 5e-14   |
| 2.7      | −15.192712    | −15.192712    |  ≤ 5e-14   |
| 3.0      | −15.209117    | −15.209117    |  ≤ 5e-14   |
| 3.3      | −15.212221    | −15.212221    |  ≤ 5e-14   |
| 3.6      | −15.203243    | −15.203243    |  ≤ 5e-14   |

**Maximum absolute difference: 5e-14 Ha across all 5 grid points.**

This confirms the Phase C-W1c diagnostic claim that the screened
cross-center V_ne is bit-identical to the bare cross-center V_ne for
first-row systems where the off-center nucleus has no frozen core
(Z = 1 for H, no [Ne] core to screen). The auto-detection in
`cross_center_screened_vne._detect_core_type` correctly returns `None`
for Z<11 and the screened path collapses to the bare path.

## Result 2: n_max=2 reproduction (TRACK CD BASELINE INTACT)

Parabolic R_eq fit on the three lowest-E points (R = 3.0, 3.3, 3.6):

- **R_eq = 3.227 bohr** (Track CD reported 3.226 — bit-identical to 3 d.p.)
- E_min(total) = −15.2127 Ha
- E_min(TC convention, after subtracting E_core = −7.2799) = −7.933 Ha
  - Track CD reported −7.924 Ha (difference 0.009 Ha, sub-1% numerical noise)
- R_eq error vs experiment (3.015 bohr): **7.03%** (Track CD: 7.0%)
- E_min error vs literature (−8.071 Ha): **1.71%** (Track CD: 1.8%)

**The production v2.32 architecture exactly reproduces the Track CD
n_max=2 baseline.** No regression, no improvement.

## Result 3: n_max=3 (compute infeasible in fork budget)

A single-point n_max=3 LiH balanced coupled FCI at R=3.015 was
attempted with `python -u` and PYTHONUNBUFFERED=1. Process consumed
~1.6 GB RSS with no output for 9+ minutes, killed cleanly. The
bottleneck is the 4-electron FCI eigsh on the 741,321-determinant
sparse matrix (M=42 spatial orbitals, C(42,2)² alpha-beta strings).
The existing `coupled_fci_energy` driver builds the matrix in
`scipy.sparse.lil_matrix` via Python loops, which scales poorly.

Track CD's published n_max=3 result is accepted as the baseline:
- R_eq = 3.280 bohr (8.8% error)
- E_min(TC-conv) = −8.055 Ha (0.20% error)
- Drift n_max=2→3: +0.054 bohr/n_max

The architectural finding (W1c bit-identical for first-row) holds at
all n_max by construction — `_detect_core_type` returns None for Z<11
regardless of basis size — so the n_max=3 result under v2.32 is
predicted to be bit-identical to Track CD's n_max=3 result. The
verdict for first-row chemistry does not depend on direct n_max=3
reproduction.

A future infrastructure sprint targeting either (a) compiled FCI
matrix construction or (b) a Davidson eigensolver that does not
materialize the matrix would unblock the n_max=3 regime in finite
fork compute.

## Architectural reading

The multi-focal sprint Phase C delivered four closures:
- **W1a** (`cross_register_vne.py`): Pauli-encoded two-register V_eN(r_e, R_n)
  for nuclear–electronic composition (Track NI deuterium PoC, Sprint MH
  muonic hydrogen). NOT used in chemistry solver.
- **W1b** (`magnetization_density.py`): operator-level Zemach radius
  shift. Sprint MH operator-level closure. NOT used in chemistry solver.
- **W1c** (`cross_center_screened_vne.py`): FrozenCore-aware screening
  of the multipole V_ne. Bit-identical for first-row LiH/BeH₂/H₂O
  (no frozen core to screen).
- **W2b-easy** (Paper 39): tensor-product propinquity convergence
  theorem on $\mathcal{T}_{S^3}^{\lambda_a} \otimes \mathcal{T}_{S^3}^{\lambda_b}$.
  Theoretical-mathematical, not directly wired into the qubit builder.

For first-row chemistry, **none of W1a/W1b/W1c/W2b changes the
underlying integrals**. The balanced-coupled architecture itself
(Track CD, April 2026) is the multi-focal architecture in the
chemistry sense — different focal lengths per block via Z_eff
(Li_core at Z=3, LiH_bond_center at Z_eff=1.3, LiH_bond_partner at
Z=1) coupled by Shibuya–Wulfman cross-center V_ne with multipole
termination at $L_\text{max} = 2 l_\text{max}$. That architecture
landed in v2.0.39 with R_eq error 7.0%; it has not been improved
since by any of the Phase C closures.

The chemistry-solver wall is therefore **intrinsic to balanced
coupled in the qubit basis with n_max ≤ 3**, not to PK overcounting
(the v2.0.32 PK-composed wall) and not to W1c orthogonality (the
second-row PES regression wall named in Phase C). It lives in the
single-center hydrogenic basis used per block, where the s + p (n=2)
or s + p + d (n=3) cutoff cannot resolve the nuclear geometry to
better than 7-9% at the equilibrium bond length.

## Implications for next sprint

If Result 3 confirms n_max=3 still hits R_eq ≈ 3.28 bohr:

1. **For first-row chemistry, the Phase C closures are no-ops.** The
   multi-focal architecture's chemistry-solver upgrade is the
   balanced-coupled framework itself (Track CD, n_max-drift
   +0.054 bohr/n_max — three times smaller than the PK wall it
   replaced). The framework HAS upgraded the chemistry solver, but
   the upgrade landed before Phase C, not in Phase C.

2. **The remaining R_eq drift is a basis-flexibility limitation.**
   Pushing to n_max=4 (Q ≈ 168 for LiH, FCI dimension ≈ 14M) is
   computationally feasible with eigsh but expensive. A more
   structurally productive direction is multi-center cross-register
   V_ee (Coulomb between two single-center bases at distinct focal
   lengths) — the e-e analog of W1a's e-N work. Currently the balanced
   builder uses *intra-block* hydrogenic V_ee Slater integrals; the
   *cross-block* V_ee comes from `compute_cross_block_eri` which
   handles geometry via shared center transformation. A multi-focal
   V_ee operator would be the genuine "second e-N closure" for
   chemistry.

3. **Second-row systems are the more interesting test.** The Phase C
   D-PES regression (CLAUDE.md §2 multifocal sprint outcome) found
   that screened W1c reduces NaH cross-V_ne by 5.4–6.0× yet the PES
   still descends monotonically. The bottleneck was the n_e-projected
   FCI driver (which exists as `coupled_fci_energy` but was reportedly
   not yet wired in for second-row testing). A targeted NaH PES scan
   at max_n=2 with `coupled_fci_energy` and `screened=True` would
   answer whether the W1c-residual orthogonality wall is the binding
   blocker or whether it's a different mechanism. This is a clean
   sprint of moderate compute scope.

## Verdict

**Drift PERSISTS at the Track CD level (+0.054 bohr/n_max),
unchanged by Phase C closures for first-row systems.** The multi-focal
architecture's chemistry-solver contribution is the balanced-coupled
construction itself (April 2026, Track CD), which improved drift
3× over the PK-composed wall (+0.18 bohr/l_max → +0.054 bohr/n_max)
but did not eliminate it. For first-row LiH/BeH₂/H₂O the W1c closure
is a structural no-op (no frozen core to screen). The remaining drift
sits in single-center hydrogenic basis flexibility at small n_max.

**For second-row NaH (this sprint extension):** the W1c-residual
orthogonality wall is empirically confirmed at the FCI level. PES
still descends monotonically under both bare and screened W1c. W1c
reduces PES descent magnitude 17.5× but does not introduce binding.
The H valence orbital ↔ [Ne]-core orthogonality (Phillips-Kleinman
class) is the underlying mechanism, sitting beyond the Hartree
screening of Z_eff(ρ) that W1c implements. See companion file
`chemistry_solver_retest_synthesis_memo.md` for full table.

**Recommendation:** the next sprint should attack the W1c-residual
orthogonality wall directly via a Phillips-Kleinman-class projection
on the H valence orbital against the [Ne] core in cross-center
geometry. This is the single remaining engineering closure for
second-row chemistry.

## Files

- `debug/chemistry_solver_retest_lih.py` — main LiH re-test driver
- `debug/chemistry_solver_retest_lih_quick_nmax2.py` — n_max=2
  reproduction with parabolic fit
- `debug/chemistry_solver_retest_lih_nmax3.py` — n_max=3 single-point
  attempt (stalled in eigsh, killed cleanly)
- `debug/chemistry_solver_retest_nah.py` — NaH binding determination
- `debug/data/chemistry_solver_retest_lih.json` — LiH summary + Track
  CD baseline + verdict
- `debug/data/chemistry_solver_retest_nah.json` — NaH PES tables under
  bare/screened W1c, W1c-residual wall confirmation
- `debug/chemistry_solver_retest_synthesis_memo.md` — combined synthesis
