# Sprint 7: Balanced Coupled FCI for Second-Row Molecules (NaH, MgH2)

**Date:** 2026-04-17
**Track:** Sprint 7 balanced second-row FCI
**Status:** COMPLETE -- both molecules show overattraction (no equilibrium)

## Setup

Computed balanced-coupled FCI PES for NaH and MgH2 at n_max=2, following the
Track CE methodology (LiH balanced FCI reference).

- **NaH:** [Ne] frozen core, 1 bond block, 2 valence electrons, Q=20, 239 balanced Pauli terms
- **MgH2:** [Ne] frozen core, 2 bond blocks, 4 valence electrons, Q=40, 1501 balanced Pauli terms

Nuclear repulsion: V_NN = sum(Z_i * Z_j / R_ij) over nuclear pairs (full charges), overriding the spec's nuclear_repulsion_constant to avoid double-counting with E_core and V_cross.

- NaH: V_NN = 11/R
- MgH2: V_NN = 2*12/R + 1/(2R) = 24.5/R (linear H-Mg-H, H at +/-R)

## Results

### NaH PES (2e FCI, Q=20)

| R (bohr) | E_FCI (Ha) | V_NN (Ha) | E_elec (Ha) |
|----------|-----------|----------|-------------|
| 2.000 | -5.8625 | 5.5000 | -11.3625 |
| 2.500 | -5.2377 | 4.4000 | -9.6377 |
| 3.000 | -4.7811 | 3.6667 | -8.4477 |
| 3.300 | -4.5870 | 3.3333 | -7.9204 |
| 3.566 | -4.4514 | 3.0847 | -7.5361 |
| 3.800 | -4.3506 | 2.8947 | -7.2453 |
| 4.000 | -4.2724 | 2.7500 | -7.0224 |
| 4.500 | -4.0871 | 2.4444 | -6.5315 |
| 5.000 | -3.8928 | 2.2000 | -6.0928 |
| 6.000 | -3.4586 | 1.8333 | -5.2919 |
| 8.000 | -2.6049 | 1.3750 | -3.9799 |
| 10.000 | -1.9935 | 1.1000 | -3.0935 |
| 15.000 | -1.3483 | 0.7333 | -2.0816 |
| 20.000 | -1.1003 | 0.5500 | -1.6503 |

**No equilibrium.** PES is monotonically increasing (energy becomes less negative
as R increases). Minimum at R=2.0 bohr (smallest R scanned). Overattraction --
consistent with Track CN finding that NaH n_max=2 balanced coupled does not
produce an equilibrium.

### MgH2 PES (4e FCI, Q=40)

| R (bohr) | E_FCI (Ha) | V_NN (Ha) | E_elec (Ha) |
|----------|-----------|----------|-------------|
| 2.000 | -11.3856 | 12.2500 | -23.6356 |
| 2.500 | -10.2983 | 9.8000 | -20.0983 |
| 2.800 | -9.8704 | 8.7500 | -18.6204 |
| 3.000 | -9.6569 | 8.1667 | -17.8236 |
| 3.261 | -9.4354 | 7.5130 | -16.9485 |
| 3.500 | -9.2656 | 7.0000 | -16.2656 |
| 3.800 | -9.0687 | 6.4474 | -15.5161 |
| 4.000 | -8.9363 | 6.1250 | -15.0613 |
| 4.500 | -8.5718 | 5.4444 | -14.0162 |
| 5.000 | -8.1441 | 4.9000 | -13.0441 |
| 6.000 | -7.1517 | 4.0833 | -11.2350 |
| 8.000 | -5.4345 | 3.0625 | -8.4970 |
| 10.000 | -4.4950 | 2.4500 | -6.9450 |
| 15.000 | -3.7655 | 1.6333 | -5.3988 |

**No equilibrium.** Same overattraction pattern as NaH. PES monotonically
increasing. Minimum at R=2.0 bohr (smallest R scanned).

### Comparison

| Molecule | Q | N_el | N_pauli | Has min | R_eq err | D_e err | Bound |
|----------|---|------|---------|---------|----------|---------|-------|
| LiH (ref) | 30 | 4 | 878 | Yes | 7.0% | 60% | Yes |
| NaH | 20 | 2 | 239 | No | N/A | N/A | No |
| MgH2 | 40 | 4 | 1501 | No | N/A | N/A | No |

## Diagnosis

Both NaH and MgH2 show the same overattraction behavior at n_max=2 that was
previously documented for NaH in Track CN. The electronic energy becomes more
negative monotonically as R decreases, overwhelming the V_NN repulsion wall.

**Root cause (from Track CN):** At n_max=2, the balanced coupled Hamiltonian
overattracts because the cross-center V_ne (Na nucleus attracting the H-side
orbital, and vice versa) is too strong relative to the V_ee repulsion in the
limited basis. The key difference from LiH: NaH and MgH2 have frozen [Ne]
cores (10 electrons), so the valence sector sees a much stronger cross-center
V_ne from Z=11 or Z=12 nuclei without the core electrons screening it
adequately. LiH has an explicit core block that provides some of this
screening.

The n_max=3 convergence for NaH (Track CN) showed energy improvement (+0.40 Ha
over n_max=2) but PES overattraction persisted (no equilibrium). This suggests
the overattraction is structural to the balanced + frozen-core architecture at
these Z values, not just a basis truncation effect.

**Implication for quantum simulation:** The balanced coupled Hamiltonians for
NaH (239 Pauli terms) and MgH2 (1501 Pauli terms) are still valid qubit
Hamiltonians for single-point quantum simulation at fixed geometries. The
overattraction only affects PES scans and equilibrium geometry determination.
For quantum resource estimation at a given R (e.g., experimental R_eq), the
Hamiltonians are well-defined.

## Files

- Script: `debug/sprint7_balanced_second_row.py`
- Data: `debug/data/sprint7_balanced_second_row.json`
- Memo: `debug/sprint7_balanced_second_row_memo.md`

## Timing

- NaH: ~0.2-0.3s per R-point (2e FCI, dim = C(10,1)^2 = 100)
- MgH2: ~135-141s per R-point (4e FCI, dim = C(20,2)^2 = 36,100)
- Total wall time: ~32 min
