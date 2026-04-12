# Track alpha-A: Hopf Twisting Spectral Comparison (S^3 vs S^1 x S^2)

## Summary

mpmath precision: 30 dps. All sums truncated when terms drop below convergence tolerance.

### Targets (from Paper 2)

- `K` = 137.036064414482
- `K_over_pi` = 43.6199340668482
- `B` = 42.0
- `B_plus_F` = 43.6449340668482
- `B_plus_F_minus_Delta` = 43.6199340668482
- `F` = 1.64493406684823
- `Delta` = 0.025
- `pi` = 3.14159265358979
- `pi_sq_over_6` = 1.64493406684823
- `alpha_inv_CODATA` = 137.035999084

### Spectral zeta values

| s | zeta_S3(s) | zeta_{S^1 x S^2}(s) | diff (S3 - S1xS2) | ratio |
|---|-----------|---------------------|-------------------|-------|
| 2 | 0.88496503342211321157 | 5.4127042986599804143 | -4.5277392652378672027 | 0.16349776093277502956 |
| 3 | 0.17436675835602830189 | 2.7612420900689163284 | -2.5868753317128880265 | 0.063147943088060191603 |
| 4 | 0.05201327503587812166 | 2.2856377583516032563 | -2.2336244833157251346 | 0.022756569734562828148 |

### Truncated Casimir traces

- `B_S3` (Paper 2, n_max=3) = 42.0  (target B=42)
- `sum_n n^2(n^2-1)` (S^3 Laplacian trace through n=3) = 84.0

S^1 x S^2 truncated Laplacian traces:

- `k=0_lmax=2`: trace(Lap) = 36.0, B-analog (degeneracy-weighted l(l+1)) = 36.0
- `kmax=1_lmax=2`: trace(Lap) = 126.0, B-analog (degeneracy-weighted l(l+1)) = 108.0
- `kmax=2_lmax=2`: trace(Lap) = 270.0, B-analog (degeneracy-weighted l(l+1)) = 180.0
- `kmax=3_lmax=2`: trace(Lap) = 504.0, B-analog (degeneracy-weighted l(l+1)) = 252.0
- `k=0_lmax=3`: trace(Lap) = 120.0, B-analog (degeneracy-weighted l(l+1)) = 120.0
- `kmax=2_lmax=3`: trace(Lap) = 760.0, B-analog (degeneracy-weighted l(l+1)) = 600.0

### Heat kernel coefficients

- S^3: vol = 2pi^2 = 19.7392088, R=6, a1 = 2pi^2 = 19.7392088
- S^1 x S^2: vol = 8pi^2 = 78.95683521, R=2, a1 = 8pi^2/3 = 26.31894507
- vol ratio (S^1 x S^2 / S^3) = 4.0 = 4 (exactly)

### Near-miss candidates (within 5%)

Total: 45 matches. Top 20 by precision:

- B_S3 Paper2 = 42.0 ~ `B` = 42.0, rel_err = 0.0
- B_S3 + zeta_S3(4) = 42.0520132750359 ~ `B` = 42.0, rel_err = 0.00123841
- B_S3 - zeta_S3(4) = 41.9479867249641 ~ `B` = 42.0, rel_err = 0.00123841
- B_S3 + zeta_S3(3) = 42.174366758356 ~ `B` = 42.0, rel_err = 0.00415159
- B_S3 - zeta_S3(3) = 41.825633241644 ~ `B` = 42.0, rel_err = 0.00415159
- zeta_S1xS2(s=4) / zeta_S3(s=4) = 43.9433540144319 ~ `B_plus_F` = 43.6449340668482, rel_err = 0.00683745
- zeta_S1xS2(s=4) / zeta_S3(s=4) = 43.9433540144319 ~ `K_over_pi` = 43.6199340668482, rel_err = 0.0074145
- zeta_S1xS2(s=4) / zeta_S3(s=4) = 43.9433540144319 ~ `B_plus_F_minus_Delta` = 43.6199340668482, rel_err = 0.0074145
- B_S3 + zeta_S1xS2(4) = 44.2856377583516 ~ `B_plus_F` = 43.6449340668482, rel_err = 0.0146799
- pi*(B + zeta_S1xS2(4) - 1/40) = 139.048894424836 ~ `K` = 137.036064414482, rel_err = 0.0146883
- pi*(B + zeta_S1xS2(4) - 1/40) = 139.048894424836 ~ `alpha_inv_CODATA` = 137.035999084, rel_err = 0.0146888
- B_S3 + zeta_S1xS2(4) = 44.2856377583516 ~ `K_over_pi` = 43.6199340668482, rel_err = 0.0152615
- B_S3 + zeta_S1xS2(4) = 44.2856377583516 ~ `B_plus_F_minus_Delta` = 43.6199340668482, rel_err = 0.0152615
- B_S3 + zeta_S3(2) = 42.8849650334221 ~ `K_over_pi` = 43.6199340668482, rel_err = 0.0168494
- B_S3 + zeta_S3(2) = 42.8849650334221 ~ `B_plus_F_minus_Delta` = 43.6199340668482, rel_err = 0.0168494
- B_S3 + zeta_S3(2) = 42.8849650334221 ~ `B_plus_F` = 43.6449340668482, rel_err = 0.0174125
- pi*(B + zeta_S3(2) - 1/40) = 134.648551282114 ~ `alpha_inv_CODATA` = 137.035999084, rel_err = 0.017422
- pi*(B + zeta_S3(2) - 1/40) = 134.648551282114 ~ `K` = 137.036064414482, rel_err = 0.0174225
- B_S3 + zeta_S3(2) = 42.8849650334221 ~ `B` = 42.0, rel_err = 0.0210706
- B_S3 - zeta_S3(2) = 41.1150349665779 ~ `B` = 42.0, rel_err = 0.0210706

### Cleanest NON-TRIVIAL match to K / K/pi / B+F / B+F-Delta

(The trivial match `B_S3 Paper2 = 42 ~ B = 42` is excluded because B is constructed to equal 42; it tests nothing.)

- source: `zeta_S1xS2(s=4) / zeta_S3(s=4)` = 43.9433540144319
- target: `B_plus_F` = 43.6449340668482
- rel_err: 0.00683745

## Interpretation

The comparison tests whether the Hopf twisting (S^3 is a non-trivial S^1 bundle over S^2) leaves a spectral signature in the difference S^3 minus S^1 x S^2. Both manifolds have the same S^2 base Laplacian eigenvalues; the twist only reshuffles how the S^1 fiber momentum combines with them. Because the fiber of S^3 is a Hopf circle (total length 2pi after normalization), the spectra coincide in the k=0 sector (l(l+1) eigenvalues of S^2) but differ for k != 0 where the twist locks k to the Hopf charge of each S^2 irrep.

The truncated Casimir trace B=42 is a discrete, combinatorial cutoff at n_max=3 (Paper 2). S^1 x S^2 has no natural n_max because there is no SO(4) structure mixing k and l; any cutoff is artificial. Comparable cutoffs (l_max=2, with k_max varying) do not reproduce B=42.

## Recommendation

**CLEAN NEGATIVE**: No difference, ratio, or combination of S^3 vs S^1 x S^2 spectral quantities reproduces K, K/pi, B+F, or B+F-Delta at better than 1e-3 precision, let alone the 8.8e-8 precision of the Paper 2 identity. The Hopf twisting does not manifest as a simple spectral-invariant signature distinguishing the two manifolds in the quantities we probed.
