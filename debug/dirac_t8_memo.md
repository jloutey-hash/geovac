# Track T8: Darwin + Mass-Velocity alpha^4 Corrections

**Sprint:** Dirac-on-S3 Tier 3, Track T8.
**Deliverables:** `geovac/fine_structure.py`, `tests/test_fine_structure.py`.
**Consumes:** T2 spin-orbit (`geovac/spin_orbit.py`), T1 closed forms (`geovac/dirac_matrix_elements.py`).

## Three-term decomposition

The full alpha^4 (Breit-Pauli) fine-structure correction has three one-body diagonal terms:

| Term | Formula | Nonzero for |
|:-----|:--------|:------------|
| Spin-orbit (T2) | E_SO = -Z^4 alpha^2 (kappa+1) / [4 n^3 l(l+1/2)(l+1)] | l >= 1 |
| Darwin (T8) | E_D = Z^4 alpha^2 / (2 n^3) | l = 0 only |
| Mass-velocity (T8) | E_MV = -(Z^4 alpha^2)/(2 n^4) [n/(l+1/2) - 3/4] | all (n, l) |

All three are exact rational expressions in (Z, n, l, alpha) -- no pi, no transcendentals beyond alpha. Algebraic in the Paper 18 taxonomy.

## Dirac formula verification

The combined result:

    E_SO + E_D + E_MV = -(Z^4 alpha^2)/(2 n^4) [n/(j+1/2) - 3/4]

Verified as an **exact sympy symbolic equality** for all 16 states with n <= 4 at Z=1, all 9 states with n <= 3 at Z=3, and at symbolic Z. Zero failures across all tested states.

## Dirac accidental degeneracy

The combined fine-structure depends on (n, j) only, not on l separately. Verified degeneracies (exact symbolic):

| Degenerate pair | j | Verified |
|:----------------|:--|:---------|
| 2s_{1/2} = 2p_{1/2} | 1/2 | yes |
| 3s_{1/2} = 3p_{1/2} | 1/2 | yes |
| 3p_{3/2} = 3d_{3/2} | 3/2 | yes |
| 4s_{1/2} = 4p_{1/2} | 1/2 | yes |
| 4p_{3/2} = 4d_{3/2} | 3/2 | yes |
| 4d_{5/2} = 4f_{5/2} | 5/2 | yes |

## He/Li/Be fine-structure comparison: unchanged

The 2p doublet splitting (j=3/2 minus j=1/2) is **unchanged** by Darwin + mass-velocity:

- Darwin is zero for both 2p states (l=1).
- Mass-velocity depends on l, not j. Both 2p_{3/2} (kappa=-2) and 2p_{1/2} (kappa=+1) have l=1, so MV is identical and cancels in the splitting.

The 2p splitting remains alpha^2 Z^4 / 32 (pure SO). The 66-211% errors vs NIST for He/Li/Be are from **multi-electron spin-spin and spin-other-orbit** terms (Direction 3, deferred), not from missing single-particle operators.

However, the individual **level energies** are shifted by Darwin + MV. For hydrogen:

| State | E_SO (Ha) | E_D (Ha) | E_MV (Ha) | E_FS total (Ha) |
|:------|:---------:|:--------:|:---------:|:---------------:|
| 1s_{1/2} | 0 | +alpha^2/2 | -5 alpha^2/8 | -alpha^2/8 |
| 2s_{1/2} | 0 | +alpha^2/16 | -13 alpha^2/128 | -5 alpha^2/128 |
| 2p_{3/2} | +alpha^2/96 | 0 | -7 alpha^2/384 | -alpha^2/128 |
| 2p_{1/2} | -alpha^2/48 | 0 | -7 alpha^2/384 | -5 alpha^2/128 |

## Hydrogen 2p splitting: alpha^2/32 confirmed

The 2p fine-structure splitting is alpha^2/32 whether computed from SO only or from the full three-term sum. This is because Darwin and MV contribute identically (or zero) to both j-branches at fixed l.

The design memo's note about "alpha^2/16" is clarified: alpha^2/16 is the *n=2 level fine-structure energy* E_FS(n=2) = -(alpha^2)/(2*16) * [2/1 - 3/4] = -5 alpha^2/128 at j=1/2, not the splitting. The splitting between j=3/2 and j=1/2 at n=2 is alpha^2/128 - 5 alpha^2/128... wait: E(2p_{3/2}) - E(2s_{1/2}) = -alpha^2/128 - (-5 alpha^2/128) = 4 alpha^2/128 = alpha^2/32. Confirmed.

## Paper 18 classification

All three terms (SO, Darwin, MV) are:
- Diagonal in (n, kappa, m_j)
- Closed-form rationals in (Z, n, l) times alpha^2
- No pi, no transcendentals beyond alpha itself
- **Algebraic** (spinor-intrinsic) in the Paper 18 taxonomy

## Test summary

43 new tests in `tests/test_fine_structure.py`, all passing:
- 8 Darwin tests
- 8 mass-velocity tests  
- 5 Dirac formula verification tests (n_max=2 to 4, Z=1 to 38)
- 8 Dirac degeneracy tests
- 4 fine-structure splitting value tests
- 1 non-relativistic limit test
- 5 input validation tests
- 4 regression tests (T2 SO unchanged)

Plus: 22 T2 + 66 T1 + 51 D1 = 139 upstream tests still passing.

Total: 43 + 139 = 182 tests passing.

## What this means for Papers 14 and 20

The He/Li/Be fine-structure comparison table does **not** improve from T8 alone. The single-particle fine-structure (SO + Darwin + MV) gives the correct doublet splitting structure for hydrogenic atoms, but multi-electron effects (spin-spin, spin-other-orbit, configuration interaction) are needed to match the NIST multiplet splittings for He 2^3P and Be 2s2p 3P. Li 2^2P is the cleanest test (single valence electron), and there the error is dominated by core screening (Z_eff approximation), not by missing relativistic terms.

The T8 deliverable completes the **single-particle alpha^4 fine-structure ladder** as a verified algebraic framework. Its value is:
1. Structural: the Dirac formula is verified as an exact identity, confirming the implementation is correct.
2. Foundation: provides the one-body diagonal terms that any future multi-electron fine-structure calculation (Direction 3) would start from.
3. Taxonomy: all three terms confirmed algebraic in Paper 18's classification.
