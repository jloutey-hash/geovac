# LiH R_eq Investigation Report

**Date:** 2026-03-11
**Version:** v0.9.37
**Status:** PARTIAL — see findings
**Script:** `debug/validate_req_investigation.py`
**Configuration:** `cross_atom_vee='s_only'` (v0.9.11 baseline)

## Purpose

Determine whether R_eq = 2.5 bohr is the raw PES minimum, the CP-corrected
PES minimum, or neither, using the v0.9.11-era `s_only` cross-atom V_ee
configuration that produced the paper's benchmark values.

## Is R_eq = 2.5 bohr the raw or CP-corrected minimum?

| PES type | R_eq (nmax=3) | D_e at R_eq (Ha) |
|----------|---------------|------------------|
| Raw      | 2.90 bohr | 0.7781 |
| CP-corrected | 2.90 bohr | 0.6632 |
| E_mol minimum | 2.90 bohr | --- |
| Experiment | 3.015 bohr | 0.092 |

BSSE R-dependence: std = 0.000000 Ha, range = 0.000000 Ha over
R = 1.5 to 4.0 bohr.

BSSE is effectively R-independent. The raw and CP-corrected PES curves are parallel (identical shape, constant vertical offset). Therefore R_eq(raw) = R_eq(CP).

## Convergence Comparison (s_only, apples-to-apples)

| Quantity | nmax=2 | nmax=3 | Paper claim | Match? |
|----------|--------|--------|-------------|--------|
| R_eq raw (bohr) | 2.89 | 2.90 | ~2.0/~2.5 | --- |
| R_eq CP (bohr) | 2.89 | 2.90 | ~2.5 | NO |
| D_cp at R=3.015 (Ha) | 0.8028 | 0.6614 | 0.093 | NO |
| BSSE at R=3.015 (Ha) | -0.2184 | -0.1149 | -0.115 | YES |
| D_cp(nmax=2) at R=3.015 | 0.8028 | --- | 0.270 | NO |

## Discrepancies

- WARNING: D_cp at R_eq (Ha) nmax=3 = 0.6632, paper = 0.0930, diff = 0.5702 Ha (> 0.005 Ha threshold)
- WARNING: D_cp at R=3.015 (Ha) nmax=3 = 0.6614, paper = 0.0930, diff = 0.5684 Ha (> 0.005 Ha threshold)
- WARNING: D_raw at R=3.015 (Ha) nmax=3 = 0.7763, paper = 0.2050, diff = 0.5713 Ha (> 0.005 Ha threshold)
- WARNING: D_cp(nmax=2) @3.015 nmax=2 = 0.8028, paper = 0.2700, diff = 0.5328 Ha (> 0.005 Ha threshold)

## Output Files

- `debug/data/lih_req_investigation.txt` — full PES data
- `debug/plots/lih_raw_vs_cp_pes.png` — nmax=3 raw vs CP comparison
- `debug/plots/lih_req_investigation.png` — three-panel convergence figure

## Next Steps

- If R_eq discrepancy persists, extend scan below R=1.5 bohr
- Verify paper's R_eq=2.5 claim against original v0.9.11 git checkout
- Run nmax=4 when computational budget allows
