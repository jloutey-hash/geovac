# LiH Overlap-Dependent Kinetic Correction Results

**Date:** 2026-03-11
**Version:** v0.9.37
**Scripts:** `debug/validate_overlap_correction.py`, `debug/validate_overlap_correction_fast.py`
**Configuration:** exact+True, nmax=3

## Method

Added a diagonal kinetic correction to the molecular H1:

```
h1_diag[a] += lambda * sum_b S^2(a,b)
```

where S(a,b) is the STO overlap integral between orbital a on one atom and
orbital b on the other. Only s-orbital (l=0) pairs are included (l>0 overlap
integrals not implemented).

This approximates the Pauli kinetic repulsion from orthogonalizing overlapping
atomic orbitals. The graph Laplacian's intra-atom kinetic energy is R-independent
(confirmed in `debug/KINETIC_DIAGNOSTIC.md`), so this correction supplies the
missing R-dependent kinetic cost.

## Atomic References

- E(Li) = -7.392086 Ha
- E(H) = -0.500000 Ha
- E_sep = -7.892086 Ha
- BSSE = -0.115 Ha (R-independent, lambda-independent)

## Lambda Scan Results

| lambda | Status | R_eq (bohr) | D_cp (Ha) | Notes |
|--------|--------|-------------|-----------|-------|
| 0.00 | NO EQ | --- | --- | Monotonically decreasing |
| 0.05 | NO EQ | --- | --- | Monotonically decreasing |
| 0.10 | NO EQ | --- | --- | Nearly flat at R=3-4 |
| **0.20** | **NEAR-EQ** | **~3.5** | **0.042** | **Shallow interior max** |
| 0.50 | NO EQ | --- | --- | Unbound (D_cp < 0) |
| 1.00 | NO EQ | --- | --- | Unbound |
| 2.00 | NO EQ | --- | --- | Deeply unbound |
| 5.00 | NO EQ | --- | --- | Deeply unbound |
| 10.0 | NO EQ | --- | --- | Deeply unbound |
| 20.0 | NO EQ | --- | --- | Deeply unbound |

**Experimental:** R_eq = 3.015 bohr, D_cp = 0.092 Ha
**Paper (v0.9.11):** R_eq ~ 2.5 bohr, D_cp = 0.093 Ha

## Key Finding: lambda=0.2 Near-Equilibrium

The lambda=0.2 PES shows a shallow interior maximum:

| R (bohr) | E_mol (Ha) | D_cp (Ha) |
|----------|------------|-----------|
| 1.500 | -8.19092 | 0.1840 |
| 2.000 | -8.08573 | 0.0788 |
| 2.500 | -8.03870 | 0.0317 |
| **3.000** | **-8.04742** | **0.0405** |
| **3.500** | **-8.04875** | **0.0418** |
| 4.000 | -8.04405 | 0.0371 |
| 5.000 | -8.02573 | 0.0188 |

The maximum at R~3.5 is barely above the neighboring points (delta D_cp ~ 0.01 Ha).
This is a marginal inflection, not a robust equilibrium:

- **D_cp = 0.042 Ha** vs experiment 0.092 Ha (54% too weak)
- **R_eq ~ 3.5 bohr** vs experiment 3.015 bohr (16% too long)
- The "well" is ~0.01 Ha deep (experiment: ~0.05 Ha deep)

## Baseline PES (lambda=0.0)

| R (bohr) | E_mol (Ha) | D_raw (Ha) | D_cp (Ha) | BSSE (Ha) |
|----------|------------|------------|-----------|-----------|
| 1.500 | -8.47115 | 0.5791 | 0.4642 | -0.1149 |
| 2.000 | -8.28726 | 0.3952 | 0.2803 | -0.1149 |
| 2.500 | -8.17761 | 0.2855 | 0.1706 | -0.1149 |
| 3.000 | -8.11880 | 0.2267 | 0.1118 | -0.1149 |
| 3.500 | -8.09220 | 0.2001 | 0.0852 | -0.1149 |
| 4.000 | -8.07939 | 0.1873 | 0.0724 | -0.1149 |
| 5.000 | -8.04860 | 0.1565 | 0.0416 | -0.1149 |

## Analysis

### Why the Overlap Correction Fails

The correction transitions directly from "bound but monotonically decreasing" (lambda<0.2)
to "unbound" (lambda>0.5) with only a marginal near-equilibrium at lambda~0.2. The
fundamental problem:

1. **Diffuse orbital overlaps don't decay fast enough.** The 2s and 3s orbitals of Li
   (Z=3) and H (Z=1) have large spatial extent. At R=5.0 bohr (well past any equilibrium),
   the total sum(S^2) is still ~0.44 — only 3x smaller than at R=1.5 (sum S^2 ~ 1.32).
   A useful kinetic repulsion needs to be ~10-30x stronger at short R than at long R.

2. **The correction raises E_sep too.** Because the overlaps remain large at all R, the
   correction adds a nearly R-independent energy penalty. This shifts the entire PES
   upward rather than creating differential repulsion at short R.

3. **Only s-orbitals contribute.** The l>0 orbitals (which would provide better spatial
   discrimination) are not included because `compute_overlap_element()` only handles l=0.
   However, the LiH ground state is s-dominated (1s^2 2s^2), so l>0 would contribute
   minimally even if implemented.

### What Would Work Instead

The overlap correction correctly identifies the MECHANISM (orthogonalization cost) but
uses too crude an implementation. More refined approaches:

1. **Weight by 1/n^2 or orbital occupancy.** Core orbitals (n=1) are compact and their
   overlap decays rapidly with R. Valence/diffuse orbitals (n=2,3) should contribute
   less to the kinetic penalty. Weighting T_corr by -Z^2/(2n^2) (the orbital energy)
   would naturally suppress diffuse contributions.

2. **Restrict to n=1 core overlaps only.** The 1s-1s overlap between Li(Z=3) and H(Z=1)
   decays from 0.13 at R=1.5 to ~0.001 at R=5.0 — exactly the R-discrimination needed.
   But it only corrects the core, not the valence.

3. **Off-diagonal kinetic corrections.** The real kinetic repulsion in standard QC comes
   from off-diagonal kinetic energy matrix elements T_AB between AOs on different centers,
   not from diagonal shifts. This would require modifying the H1 matrix structure.

4. **Bond sphere approach (Paper 8).** Use the SO(4) Wigner D-matrix to couple atoms on
   a single S^3, which builds in the correct R-dependent cross-nuclear physics by
   construction rather than as an ad hoc correction.

5. **Attenuate cross-nuclear attraction.** The R_eq problem may be better addressed by
   weakening the Fourier diagonal cross-nuclear term (which is known to overbind at short R)
   rather than adding kinetic repulsion.

## Output Files

- `debug/data/lih_overlap_correction_scan.txt` -- raw data for all lambda values
- `debug/plots/lih_overlap_correction.png` -- D_cp(R) curves
- `docs/OVERLAP_CORRECTION_RESULTS.md` -- this report
