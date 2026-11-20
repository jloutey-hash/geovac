# LiH nmax Convergence Validation Report

**Date:** 2026-03-11
**Version:** v0.9.37
**Status:** PAPER CORRECTION REQUIRED
**Script:** `debug/validate_lih_nmax2_convergence.py`

## Purpose

Validate the convergence trend claimed in Section V.C of the FCI paper:
- D_e_CP(nmax=2) = 0.270 Ha (192% error)
- D_e_CP(nmax=3) = 0.093 Ha (1.0% error)
- BSSE decreases with increasing nmax
- R_eq converges outward toward experiment with increasing nmax

## Atomic Reference Energies

| Atom | nmax=2 (Ha) | nmax=3 (Ha) |
|------|-------------|-------------|
| Li   | -7.101226 | -7.392086 |
| H    | -0.500000 | -0.500000 |

## Convergence Comparison

| Quantity | nmax=2 | nmax=3 | Expt | Trend |
|----------|--------|--------|------|-------|
| D_e_CP at R=3.015 (Ha) | 0.289 | 0.093 | 0.092 | Converging |
| D_e_CP at R_eq (Ha) | 0.644 | 0.464 | 0.092 | Converging |
| BSSE at R=3.015 (Ha) | -0.218 | -0.115 | 0.000 | Converging |
| R_eq (bohr) | 1.50 | 1.50 | 3.015 | WRONG DIRECTION |
| Error in D_e_CP (%) | 214 | 1 | --- | --- |

## Verification of Paper Section V.C Claims

| Claim | Paper Value | Computed Value | Match |
|-------|-------------|----------------|-------|
| D_e_CP(nmax=2) | 0.270 Ha | 0.289 Ha | NO |
| Error(nmax=2) | 192% | 214% | NO |

## WARNING

**PAPER CORRECTION REQUIRED:** Section V.C states D_e_CP(nmax=2) = 0.270 Ha
but computed value is 0.289 Ha.

## Convergence Diagnostics

1. **BSSE convergence:** |BSSE(nmax=2)| = 0.218 Ha > |BSSE(nmax=3)| = 0.115 Ha — PASS
2. **R_eq convergence:** R_eq(nmax=2) = 1.50 > R_eq(nmax=3) = 1.50 bohr — UNEXPECTED
3. **D_e_CP convergence:** D_e_CP(nmax=2) = 0.644 > D_e_CP(nmax=3) = 0.464 Ha — PASS (overbinding decreases)

## Output Files

- `debug/data/lih_nmax_convergence.txt` — numerical data
- `debug/plots/lih_pes_convergence.png` — PES comparison plot

## Next Steps

- Run nmax=4 when computational resources allow (8.2M SDs, estimated >30 min)
- Verify R_eq continues to shift outward toward 3.015 bohr
- Verify BSSE continues to decrease
