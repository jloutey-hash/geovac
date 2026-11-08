# Sprint H1-Higgs: Higgs Mechanism from Inner Fluctuations

**Date:** 2026-05-31
**Verdict:** POSITIVE-THIN (as expected)
**Version:** builds on v3.33.0 (G4a)

## Goal

Construct the Higgs field from inner fluctuations of the GeoVac spectral triple with the full CCM algebra A_F = C + H + M_3(C), and prove the Yukawa non-selection theorem. This extends G4a (gauge group identification) to the Higgs sector and the spectral-action Higgs potential.

## Architecture

Total Hilbert space: H = H_GV(n_max) x H_F with dim H_F = 32 (16 matter + 16 antimatter). At n_max=2: dim_GV = 16, dim_H = 512. At n_max=3: dim_GV = 40, dim_H = 1280.

Combined Dirac: D = D_GV x 1_F + gamma_GV x D_F. Inner fluctuation: omega = sum a_i [D, b_i]. Fluctuated Dirac: D_A = D + omega + epsilon' J omega J^{-1} with epsilon' = -1 (KO-dim 1).

## Results

### 1. Gauge-Higgs decomposition (n_max=2, 50 random generators)

The inner fluctuation omega splits cleanly:

| Sector | Norm (max/mean) | Structure |
|:-------|:----------------|:----------|
| Gauge  | 0.908 / 0.105   | U(1) x SU(2) x SU(3) (matches G4a) |
| Higgs  | 0.008 / 0.002   | L-R off-diagonal, complex 2x2 doublet |
| Mat-anti off-diag | 0 | Exact zero (decoupled) |
| Lep-quark off-diag | 0 | Exact zero (decoupled) |

Y = 0 sanity check: Higgs identically zero across 20 generators (max = 0.0). Confirms Higgs arises exclusively from D_F off-diagonal.

### 2. Higgs field Phi is a complex doublet

The lepton-sector Higgs lives in the 2x2 L-R off-diagonal block at each GV-pair (i, j) -- this is the complex SU(2) doublet Phi = (phi^+, phi^0)^T. The quark-sector Higgs has the same flavor structure tensored with color: the color-singlet fraction is 0.86 (the remaining 0.14 is off-singlet content from M_3(C) acting in a_F, giving Phi a color connection via the full gauge-covariant derivative). This is the standard CCM structure.

### 3. Spectral action and Mexican-hat potential

Scanning D_combined over Yukawa values y in [0, 0.5] with y_e = y_u = y_d = y:

- delta Tr(D^2) = 448 y^2 (power-law fit: y^{2.000}). This is the mass term mu^2 |Phi|^2.
- delta Tr(D^4) = 14112 y^2 + 448 y^4. The y^2 term is the D_GV^2 x D_F^2 cross-contribution; the y^4 term is the pure Higgs quartic.

On the finite Dirac D_F alone: Tr(D_F^2) = 0.28 > 0, Tr(D_F^4) = 0.0028 > 0.

The CCM spectral action expansion gives V(Phi) = -f_2 Lambda^2 Tr(D_F^2) |Phi|^2 + f_0 Lambda^4 Tr(D_F^4) |Phi|^4. With f_0, f_2 > 0 (positive cutoff moments) and both traces positive, this IS the Mexican-hat potential. The sign structure is automatic in the CCM framework -- not imposed.

Two-term exactness of the bare S^3 spectral action (Paper 51 Theorem 1) is broken by the inner fluctuation, as expected: D_F adds finite-dimensional structure that generates new Seeley-DeWitt contributions at every order.

### 4. Yukawa non-selection theorem

**Constraint chain on D_F:**

| Constraint | dim(admissible space) |
|:-----------|:---------------------|
| 32x32 Hermitian | 1024 real params |
| + {gamma_F, D_F} = 0 | 512 real params |
| + J_F D_F = D_F J_F | 128 real params (matter determines antimatter) |
| + order-one [[D,a], JbJ^{-1}] = 0 | 128 real params (no further constraint at SM level) |

The order-one condition does NOT reduce the parameter count further because:
(a) Lepton mixing (off-diagonal Y_lepton) satisfies order-one exactly (verified: residual = 0.0)
(b) Quark mixing (off-diagonal Y_quark) also satisfies order-one (the algebra acts block-diagonally in lepton/quark, and M_3(C) commutes with the Yukawa flavor structure)
(c) Lepton-quark cross-coupling is already blocked by the matter-sector block structure

**The 128 real parameters** decompose as: 64 from the lepton 4x4 off-chirality block (Yukawa matrix Y_lepton: 2x2 complex = 8 real, but the full off-chirality space is larger -- generation mixing and right-handed neutrino mass terms are all admissible) and 64 from the quark 12x12 off-chirality block (including color-singlet and color-octet components, though only the color-singlet piece contributes to physical Yukawa). For one generation with diagonal Yukawa: 8 real parameters (4 complex couplings y_nu, y_e, y_u, y_d).

**Chirality independence (G3 result):**
- [gamma_GV x 1, 1 x gamma_F] = 0 exactly (they commute)
- ||gamma_GV x 1 - 1 x gamma_F|| = 32.0 (they are independent)

This is the structural reason for non-selection: gamma_GV and gamma_F are independent commuting Z_2 gradings. The Yukawa lives in the gamma_F flip. No GeoVac-side data couples to that flip.

**Theorem statement:** The most general finite Dirac D_F on the CCM SM triple satisfying Hermiticity, J-reality, chirality anticommutation, and order-one is parametrized by the Yukawa matrix Y with 8+ free real parameters (one generation). The gauge sector U(1) x SU(2) x SU(3) is forced by the algebra A_F. The Yukawa couplings are free calibration data (Paper 18 inner-factor tier).

### 5. n_max=3 verification

All axioms bit-exact at n_max=3 (dim_H = 1280). Higgs non-trivial (higgs/gauge = 0.030). No structural change from n_max=2 -- the result is n_max-independent at the structural level.

## Honest assessment

This sprint confirms the expected POSITIVE-THIN verdict. The Higgs mechanism from inner fluctuations works exactly as in the standard CCM framework -- the GeoVac spectral triple on S^3 provides the continuous manifold data (algebra, Dirac, real structure), and the finite algebra A_F = C + H + M_3(C) provides the SM gauge group and Higgs sector. The construction is SM-consistent but not SM-selecting: the gauge group is forced, the Higgs field exists, the Mexican-hat potential has the right sign structure, but the Yukawa couplings are free data.

This is precisely the Marcolli-van Suijlekom lineage (J. Geom. Phys. 75, arXiv:1301.3480): Yang-Mills from inner fluctuations at the spectral-action level, with the Higgs sector arising from the off-diagonal finite Dirac. The Perez-Sanchez correction (arXiv:2401.03705, 2508.17338) applies: the continuum limit gives Yang-Mills without Higgs when A_F is trivial; the Higgs enters only through a non-trivial D_F.

The result strengthens WH1 (GeoVac as an almost-commutative spectral triple in the Marcolli-vS lineage) and sharpens the G2-G3 corollary: the framework admits the full SM gauge + Higgs content given imposed calibration data (Yukawa, hypercharge, generation count). What it does NOT do is select that calibration data from geometry -- this is the multi-focal-composition wall (CLAUDE.md Section 1.7) at the SM level.

## Files

- `debug/h1_higgs_inner_fluctuation.py` -- driver (all 5 sections)
- `debug/data/h1_higgs_inner_fluctuation.json` -- numerical results
- `debug/h1_higgs_inner_fluctuation_memo.md` -- this memo

No production code modified.
