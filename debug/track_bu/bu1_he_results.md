# Track BU-1: Sturmian He — Theory + Validation Results

**Date:** 2026-04-02
**Version:** v2.0.33 (pre-release)
**Module:** `geovac/sturmian_solver.py`

## Method

Coulomb Sturmian basis: each orbital (n,l,m) uses Z_eff = n*k where k is a
common scaling parameter. All Sturmians share exp(-kr) decay rate and energy
-k^2/2. Basis is Lowdin-orthogonalized before standard FCI via Slater-Condon
rules. Full R^k Slater integrals + Gaunt angular coefficients.

Comparison: standard hydrogenic FCI with uniform Z_eff optimized variationally.

## Results

### He Ground State Energy (exact: -2.903724 Ha)

| max_n | N_spatial | N_SD | Method       | Z_eff/k_opt | Energy (Ha)  | Error (%) | dE vs std (mHa) |
|-------|-----------|------|--------------|-------------|-------------|-----------|------------------|
| 2     | 5         | 45   | Standard FCI | 1.820       | -2.807809   | 3.3032    | —                |
| 2     | 5         | 45   | Sturmian CI  | 1.812       | -2.811609   | 3.1723    | -3.8             |
| 3     | 14        | 378  | Standard FCI | 1.870       | -2.817209   | 2.9795    | —                |
| 3     | 14        | 378  | Sturmian CI  | 2.180       | -2.844020   | 2.0561    | -26.8            |

### Hydrogen Validation (exact: -0.500000 Ha)

| max_n | Method       | Energy (Ha)  |
|-------|--------------|-------------|
| 1     | Standard FCI | -0.500000   |
| 1     | Sturmian CI  | -0.500000   |
| 2     | Standard FCI | -0.500000   |
| 2     | Sturmian CI  | -0.500000   |

### Variational Bound

All energies above exact: PASS for all cases.

### V_ee Sparsity

| max_n | ERI nonzero (Sturmian) | ERI nonzero (Standard) | Identical? |
|-------|------------------------|------------------------|------------|
| 2     | 65                     | 65                     | YES        |
| 3     | 1492                   | 1492                   | YES        |

Gaunt selection rule sparsity is PRESERVED. The angular coupling coefficients
(Wigner 3j / Gaunt integrals) are identical in both bases because they depend
only on (l,m) quantum numbers, not on radial functions.

### Orbital Exponents (Sturmian, max_n=3, k_opt=2.180)

| n | l | Z_eff = n*k |
|---|---|-------------|
| 1 | 0 | 2.180       |
| 2 | 0 | 4.361       |
| 2 | 1 | 4.361       |
| 3 | 0 | 6.541       |
| 3 | 1 | 6.541       |
| 3 | 2 | 6.541       |

### Overlap Matrix Properties (max_n=3, k=2.180)

- Max off-diagonal: 0.50
- Condition number: 5.83
- Positive definite: YES

## Key Findings

1. **Sturmian CI improves on standard FCI at same max_n.** The advantage grows
   with basis size: 3.8 mHa at max_n=2, 26.8 mHa at max_n=3. This is because
   the n-dependent Z_eff provides automatic radial flexibility that the uniform
   Z_eff basis lacks.

2. **Gaunt selection rule sparsity is exactly preserved.** The ERI nonzero counts
   are identical in both bases. This is expected: Gaunt coefficients are purely
   angular and independent of the radial functions.

3. **Variational bound is respected** in all cases. The Lowdin orthogonalization
   + standard FCI ensures variational correctness.

4. **k_opt shifts upward with max_n** (1.81 -> 2.18). At max_n=3, the 1s orbital
   has Z_eff=2.18 (close to Z=2), while the 3s has Z_eff=6.54 (very contracted).
   The CI mixing compensates for the contracted higher-n Sturmians.

5. **The Sturmian approach does NOT differentiate core from valence automatically.**
   All shells have Z_eff proportional to n: core orbitals have SMALLER Z_eff,
   valence orbitals have LARGER Z_eff. This is the opposite of physical
   core-valence screening. For He (no core-valence distinction), this is not
   a problem. For Li (BU-3), this will be the critical test.

## BU-1 Exit Assessment

**PASS (MARGINAL).** The Sturmian CI improves He energy by 0.13 pp at max_n=2
and 0.92 pp at max_n=3. The advantage is real and grows with basis size.
Sparsity is preserved. Variational bound respected. The improvement is modest
for He (2 electrons, no core-valence), but the growing advantage with max_n
suggests the Sturmian basis captures correlation effects more efficiently.

**Proceed to BU-2** (qubit encoding comparison).

## Wall Times

- max_n=2: ~16s total (both methods)
- max_n=3: ~267s total (both methods, dominated by R^k integrals)
