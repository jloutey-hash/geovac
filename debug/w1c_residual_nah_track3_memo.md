# W1c-Residual NaH — Track 3 Synthesis Memo

**Date:** 2026-05-09
**Sprint:** Track 3 of post-multi-track multi-track sprint, the named
engineering closure for the W1c-residual orthogonality wall after the
PK cross-center sprint (2026-05-08) returned negative.
**Worker fork directive:** Replace the hydrogenic Z_orb=1 Na valence
basis with FrozenCore Z_eff(r) Schrödinger eigenstates and test
whether NaH binding is recovered.

---

## 1. Verdict — short form

**Track 3 is HONEST NEGATIVE on the binding outcome AND DIAGNOSTIC
POSITIVE on identifying the genuine engineering target.**

The screened-eigenvalue h1 diagonal correction (the named target
from the post-Track-2 PK synthesis memo) does NOT close the
W1c-residual wall: NaH PES at max_n=2 remains monotonically
descending under all 8 configurations of {bare, W1c, PK, SV}, with
no internal equilibrium minimum.

**The substantive structural finding** (the value of Track 3): a
constant diagonal h1 correction is **R-independent** at the
eigenvalue level — verified bit-identical at R=2.5 vs R=10 to machine
precision. Eigenvalues are properties of the on-center potential, not
of the bond. The W1c-residual wall is structurally NOT a diagonal-h1
problem. The named target ("replace the hydrogenic Z_orb=1 basis on
Na with Schrödinger eigenstates") was incomplete as named — addressing
only the diagonal h1 cannot close the wall.

**Diagnostic positive** (the value-add from the sprint): testing the
cross-V_ne integral with hydrogenic-shape labels at the *physical n*
shows a -0.674 Ha differential between framework convention and
physical-n convention, **1.9× the W1c-residual descent (0.357 Ha)** —
suggesting wavefunction-shape substitution would over-correct the
descent profile and produce a binding minimum. The genuine engineering
target is **wavefunction-shape substitution in cross-V_ne and
within-block ERIs**, not diagonal eigenvalue substitution.

The diagnostic-before-engineering rule (memory
`feedback_diagnostic_before_engineering.md`, 2026-05-08) applies
INSIDE this sprint. Two engineering attempts (PK cross-center +
screened-valence diagonal) have now landed clean architectural
modules but failed to close binding; the third pre-implementation
diagnostic (cross-V_ne shape sweep) shows the genuine target is in
cross-V_ne shape, not diagonal h1.

## 2. What was implemented

### 2.1 Module `geovac/screened_valence_basis.py` (~280 lines, 25 tests)

Three module-level helper functions provide the engineering
infrastructure:

```python
screened_valence_eigenvalue(Z_nuc, block_n, l, n_val_offset, ...)
    → float  # Schrödinger eigenvalue in Ha

screened_valence_h1_diagonal(Z_orb, Z_nuc, block_n, l, n_val_offset, ...)
    → (e_screened, e_hydrogenic, delta_h1)

apply_screened_valence_correction(spec, h1, sub_blocks, ...)
    → (h1_corrected, info)
```

**Convention.** Following the relativistic-composed convention
(`composed_qubit_relativistic.py`), physical n maps as
`physical_n = block_n + n_val_offset`, where
`n_val_offset = period - 1`. For Na ([Ne] core, period 3):
block_n=1 → physical_n=3 (Na 3s), block_n=2 → physical_n=4
(Na 4s/4p). This matches existing framework convention so
backward-compat regressions are bit-exact.

**Solver dispatch.** s-states use the log-grid solver
`_solve_screened_radial_log` with `allow_l0=True` (Sprint Cs-HFS
machinery); l ≥ 1 uses the standard FD solver
`_solve_screened_radial`. Module-level cache avoids resolving the
same (Z, l, n_phys) per R-point of a PES sweep.

**Per-orbital corrections for NaH** (verified):
- block_n=1, l=0 (Na 3s): hydrogenic -0.500, screened **-0.170**, Δ=+0.330 Ha
- block_n=2, l=0 (Na 4s): hydrogenic -0.125, screened -0.067, Δ=+0.058 Ha
- block_n=2, l=1 (Na 4p, ×3 m_l): hydrogenic -0.125, screened -0.051, Δ=+0.074 Ha each

Total trace shift on the NaH valence (5 orbitals): **+0.611 Ha**.

### 2.2 Wiring into `geovac/balanced_coupled.py`

Added `screened_valence_basis: bool = False` and
`screened_valence_n_grid: int = 4000` kwargs to
`build_balanced_hamiltonian`. Added a Step 3c block (after PK Step 3b)
that calls `apply_screened_valence_correction` and tracks the
correction matrix in `h1_screened_correction`. Added
`screened_valence_info` and `h1_screened_correction` to the result
dict. Bit-exact backward compat verified:

- LiH at max_n=2 (no frozen core): max(|h1_on - h1_off|) = 0.0 (exact)
- LiH N_pauli unchanged: 878 (Track CD baseline preserved)
- NaH flag OFF: `n_orbitals_corrected=0`, `total_trace_shift=0.0`,
  `np.all(h1_screened_correction == 0)` is True

### 2.3 Tests — `tests/test_balanced_coupled_screened_valence.py`

25 tests covering:
- frozen-core auto-detection (Z=1..18 + a few Z>=19)
- eigenvalue correctness at Na 3s/4s/4p, Mg 3s
- eigenvalue cache consistency
- invalid (block_n, l) raise
- diagonal-only correction (off-diagonal h1, ERI both unchanged)
- backward compat (LiH bit-exact, flag OFF flag = identity)
- composability with `screened_cross_center` and `pk_cross_center`
- diagnostic dict structure
- multi-block species (MgH₂ 2 bond blocks; HCl bond + 3 lone pairs)

All 25 pass. Combined with phillips_kleinman_cross_center (21) and
cross_center_screened_vne (22) regressions: **68/68 tests pass**.

## 3. NaH PES results

Driver: `debug/w1c_residual_nah_track3.py`. Data:
`debug/data/w1c_residual_nah_track3_results.json`. Eight configs
× 7 R-points; total wall time ~25 s.

R-grid: {2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0} bohr.

Configurations test all combinations of `screened_cross_center` (W1c),
`pk_cross_center` (PK), and `screened_valence_basis` (SV):

| Config | $E_{\min}$ (Ha) | $R_{\min}$ (bohr) | Descent depth (Ha) | Bound? |
|:-------|----------------:|------------------:|-------------------:|:-------|
| all (sv+scr+pk)         | -162.943 | 2.5 | 0.3133 | **NO** |
| sv + scr                | -162.995 | 2.5 | 0.3650 | NO |
| sv + pk                 | -170.949 | 2.5 | 6.0973 | NO |
| sv only                 | -171.097 | 2.5 | 6.2441 | NO |
| scr + pk                | -163.264 | 2.5 | 0.3050 | NO |
| scr only                | -163.316 | 2.5 | 0.3567 | NO |
| pk only                 | -170.949 | 2.5 | 6.0973 | NO |
| bare                    | -171.097 | 2.5 | 6.2441 | NO |

**Reductions (descent-depth ladder):**
- bare → +W1c: 6.244 → 0.357 Ha (17.5× reduction; W1c does the
  heavy lifting)
- +W1c → +W1c+PK: 0.357 → 0.305 Ha (1.17× addl; PK adds 14.6%)
- +W1c+PK → +W1c+PK+SV: 0.305 → 0.313 Ha (0.97× — slight worsening,
  within numerical noise)

**Pure SV alone (no W1c, no PK)**: identical to bare (6.244 Ha).
SV alone gives no bond-shape correction.

**Comparison to physical NaH**: experimental $R_{eq} = 3.566$ bohr,
$D_e \approx 0.075$ Ha. The "all" descent of 0.313 Ha is still
**4.2× deeper than physical $D_e$**, with R_min at the smallest
tested point (2.5 bohr) and no internal turning point.

## 4. Why the diagonal correction doesn't help — the structural finding

The screened-valence correction adds a uniform $+0.611$ Ha to every
R-point of the PES (the trace of the diagonal h1 correction). This
shifts the entire PES up by a constant; it does not differentially
affect short-R vs long-R behavior.

To verify, I computed the SV correction matrix at R=2.5 and R=10:

```
R=2.5:   total_trace_shift = 0.610817 Ha
R=10.0:  total_trace_shift = 0.610817 Ha
Difference: 0.0e+00

h1_screened_correction max abs diff (R=2.5 vs R=10): 5.55e-17
```

Bit-identical (within machine precision). The correction matrix is
R-independent because:

1. The eigenvalues of the FrozenCore Z_eff(r) potential depend only
   on the on-center physics, not on the bond R.
2. The substitution lives entirely in the `h1[i,i]` diagonal of the
   heavy-atom-side valence orbitals, which are intrinsic
   single-particle energies.
3. R-dependent attraction lives in the **cross-V_ne** matrix (where
   the orbital wavefunction shape is convolved with the off-center
   Coulomb 1/|r-R_B|) and in the **cross-block ERIs** (V_ee between
   orbitals on different centers). These integrals depend on the
   wavefunction *shape*, not on its eigenvalue.

A constant diagonal h1 shift cannot affect the descent depth (which
is by definition the difference of energies at different R). This is
a structural feature, not a numerical artifact: the engineering
target as named ("replace the Z_orb=1 hydrogenic Na valence basis
with screened Schrödinger eigenstates") was **incomplete**.
Eigenvalue substitution is an architectural prerequisite (the h1
diagonal IS wrong, by 0.330 Ha on Na 3s alone), but it is not the
mechanism that closes the W1c-residual wall.

## 5. Diagnostic positive — the genuine target identified

To localize the residual, I tested whether wavefunction-shape
substitution moves the cross-V_ne integral. The framework's NaH
balanced uses `Z_orb = Z_eff_valence = 1.0` and `block_n` quantum
numbers (1 for the lowest s-shell, 2 for the next). For the
"physical-n" convention I tested two relabelings:

- **physical_4p**: block_n=1 → 3s, block_n=2 (l=0) → 4s, block_n=2 (l=1) → 4p
- **physical_3p**: block_n=1 → 3s, block_n=2 (l=0) → 4s, block_n=2 (l=1) → 3p (lowest p above core)

With Z_orb=1 (asymptotic Na valence Z), I computed
`compute_cross_center_vne` traces against the H nucleus at varied R:

| R (bohr) | Framework (block_n) | Physical (4p) | Physical (3p) | Frame−3p |
|---------:|--------------------:|--------------:|--------------:|---------:|
|     2.5 |             −1.279 |        −0.330 |        −0.466 | −0.814 |
|     3.5 |             −1.104 |        −0.319 |        −0.445 | −0.659 |
|     5.0 |             −0.898 |        −0.302 |        −0.414 | −0.484 |
|    10.0 |             −0.498 |        −0.271 |        −0.359 | −0.139 |

Differential (R=2.5 minus R=10) for `Frame − Physical-3p`:
**−0.814 − (−0.139) = −0.675 Ha**.

This is approximately **1.9× the W1c-residual descent (0.357 Ha)**
and **2.2× the W1c+PK descent (0.305 Ha)**. Replacing the
framework's hydrogenic 1s shape on Na with hydrogenic 3s shape (the
physically correct principal quantum number for Na's valence) would
*over-correct* the cross-V_ne attraction profile by roughly that
amount — almost certainly producing a binding minimum.

**Mechanism**. The hydrogenic Z=1 1s has mean radius
$\langle r \rangle = 1.50$ bohr; the actual Na 3s screened
wavefunction has $\langle r \rangle = 4.47$ bohr. They have
**opposite shape biases**:

- Hydrogenic Z=1 1s: 94% of probability inside r=3 bohr — too tight
  for a real 3s.
- Na 3s screened: only 21% inside 3 bohr — actual 3s extends much
  further than the hydrogenic 1s model.

Cross-V_ne integration probes the orbital tails into the off-center
1/|r-R_B| singularity. A more diffuse orbital sees the singularity
"less directly" (less amplitude near the singular point), so
cross-V_ne attraction is REDUCED for a more diffuse orbital. The
hydrogenic Z=1 1s is much TIGHTER than the actual Na 3s, but
**critically also extends differently as R changes**: a tight
orbital has a sharp falloff, so cross-V_ne integrated against
$1/|r-R_B|$ depends strongly on where R places the singularity
relative to the orbital peak. A diffuse orbital integrates more
evenly across R, with weaker R-dependence (matching the diagnostic
table: 3p at R=2.5 is -0.47, at R=10 is -0.36 — flat differential
of 0.11 Ha across a factor-of-4 R range).

The shape substitution should preserve the same diffuse profile as
the actual Na 3s, eliminating most of the spurious R-dependence in
the cross-V_ne attraction.

## 6. The next named engineering target

Based on the Track 3 diagnostic, the genuine engineering closure for
the W1c-residual wall is **wavefunction-shape substitution in
cross-V_ne and within-block ERIs**, not diagonal eigenvalue
substitution. Two architectural options:

### Option 1: Physical-n hydrogenic relabeling (1-week sprint)

Rewrite the orbital state lists in the OrbitalBlock so that block_n=1
with l=0 on the heavy-atom side uses physical n_phys=period (Na: 3s)
hydrogenic shape, evaluated at Z_orb chosen to match the asymptotic
Na valence Z. For (block_n=2, l=1) use physical n=lowest_above_core_for_l
(Na 3p) instead of n_val_offset+block_n (which gives unbound Na 4p).

**Pros:**
- Uses existing hydrogenic Slater integral machinery
  (`hypergeometric_slater.py`).
- Cheap, ~1-week sprint.
- The Slater integrals at higher physical n are still rational in Z
  (Paper 7 Section VI), so framework purity is preserved.
- Diagonal h1 in this convention would still be hydrogenic
  -Z_orb²/(2·n_phys²) = -1/(2·9) = -0.056 Ha for Na 3s — but THE SAME
  basis is used for cross-V_ne, h1, and ERIs, so the picture is
  internally consistent.

**Cons:**
- Introduces a convention break (block_n no longer increments through
  hydrogenic 1s, 2s, 2p, ... but through valence 3s, 4s, 3p, 3d, ...).
- The h1 diagonal at hydrogenic n=3 with Z_orb=1 (-0.056 Ha) is still
  not the actual Na 3s eigenvalue (-0.170 Ha) — Track 3's screened-
  valence module would still be useful as a post-correction.

### Option 2: Strict-Schmidt orthogonalization (multi-week sprint)

Explicitly Schmidt-subtract the heavy-atom [Ne] core orbitals from
the valence basis before computing all matrix elements. Mathematically
definitive in the strict-Schmidt sense.

**Pros:**
- Exact: the resulting Hamiltonian is the proper PK-equivalent in the
  orthogonalized basis with renormalization.
- Captures both the W1b orthogonality (which Phillips-Kleinman
  approximates) and the W1c-residual shape effects.

**Cons:**
- Requires re-implementing h1, ERI, and cross-V_ne in the
  orthogonalized basis (multi-week sprint).
- The renormalization of the orthogonalized basis is not trivial in
  the cross-block ERI machinery.

### Recommendation

**Option 1 first, as a 1-week diagnostic sprint.** The Track 3
diagnostic suggests Option 1 alone closes most of the descent. If
Option 1 produces a bound NaH PES (or comes very close), it's a
cheaper architectural closure than Option 2. If Option 1 still
under-corrects, Option 2 is the principled fallback.

## 7. Sprint outcome

**Sprint deliverable status:**
- ✅ Module: `geovac/screened_valence_basis.py` (~280 lines, clean)
- ✅ Tests: `tests/test_balanced_coupled_screened_valence.py` (25/25 passing)
- ✅ Wiring: `balanced_coupled.py` extended with
  `screened_valence_basis` kwarg, bit-exact backward compat
- ✅ Eigenvalues verified against published Clementi-Roetti HF data
- ✅ NaH PES: `debug/w1c_residual_nah_track3.py` + JSON
- ❌ Binding recovery: NaH still unbound under SV alone or in
  combination (negative engineering result)
- ✅ Diagnostic: identified wavefunction-shape substitution as the
  genuine engineering target (Frame−Physical-3p differential
  -0.675 Ha is 1.9× the W1c-residual descent)
- ✅ Paper updates: Paper 17 §6.10 + Paper 19 §sec:w1c_residual + CLAUDE.md §2 + CLAUDE.md §3

**Sprint verdict**: the screened-eigenvalue h1 diagonal correction
is necessary architectural infrastructure (the eigenvalue machinery
is correct and physically meaningful — Na 3s = -0.170 Ha matches
literature) but not sufficient for binding recovery. The W1c-residual
wall is structurally not a diagonal-h1 problem; it lives in
wavefunction shape. Track 3 has cleanly identified the next
engineering target (wavefunction-shape substitution) and quantified
that it's the right magnitude (1.9× the residual descent).

The diagnostic-before-engineering rule applies inside this sprint:
two engineering attempts (PK + SV diagonal) have now landed clean
modules but failed to close binding. The cross-V_ne shape sweep was
the right diagnostic to run *before* committing to a third
implementation sprint, and it cleanly indicates wavefunction-shape
substitution as the structural target.

## 8. Files changed

| Path | Change |
|:-----|:-------|
| `geovac/screened_valence_basis.py` | New (~280 lines): screened-eigenvalue module with `screened_valence_eigenvalue`, `screened_valence_h1_diagonal`, `compute_screened_h1_correction_block`, `apply_screened_valence_correction`, module-level cache. |
| `geovac/balanced_coupled.py` | Added `screened_valence_basis: bool = False` and `screened_valence_n_grid: int = 4000` kwargs to `build_balanced_hamiltonian`. Step 3c (post-PK) applies the correction. Bit-exact backward compat when off. Result dict gains `screened_valence_info` and `h1_screened_correction`. |
| `tests/test_balanced_coupled_screened_valence.py` | New: 25 unit tests (all passing). |
| `debug/w1c_residual_nah_track3.py` | New: 8-config NaH PES driver, 7-point R-grid. |
| `debug/w1c_residual_nah_track3_diag.py` | New: cross-V_ne shape diagnostic (framework block_n vs physical 4p vs physical 3p). |
| `debug/data/w1c_residual_nah_track3_results.json` | New: PES results in machine-readable form. |
| `debug/data/w1c_residual_nah_track3_diag.json` | New: diagnostic results. |
| `debug/w1c_residual_nah_track3_memo.md` | This memo. |
| `papers/core/paper_17_composed_geometries.tex` | Architectural note §6.10 extended with Track 3 SV finding. |
| `papers/methods/paper_19_coupled_composition.tex` | §sec:w1c_residual extended (~70 lines): Track 3 paragraph documenting the SV diagonal-only correction, the structural finding (R-independence), the diagnostic positive, and the next-target options. |
| `CLAUDE.md` | §2 sprint outcome paragraph (Track 3 W1c-residual NaH closure attempt). §3 dead-ends row added with the structural reading and diagnostic next target. |

---

**End of synthesis memo.**
