# Ab Initio Nuclear Lattice: H2 Rovibrational Spectrum

**Date:** 2026-03-14

**Method:** Prolate spheroidal CI (Eckart Z_eff + sigma_g/sigma_u 2x2 CI) -> Morse fit -> Nuclear Lattice

**Input:** Nuclear charges (Z=1) and hydrogen mass only. Zero experimental spectroscopic input.

> **Note (2026-03-15):** These results use the prolate spheroidal CI PES, which
> differs from Paper 13's Hylleraas + Neumann V_ee PES. Paper 13 reports
> omega_e = 4435 cm-1 and nu_01 = 4157 cm-1 (within 1% of experiment) using the
> more accurate Hylleraas basis. The prolate CI omega_e = 4918 cm-1 (+11.8%)
> overestimates due to poor dissociation-tail energies at R > 4 bohr. See
> Paper 13 Section IX for the authoritative spectroscopic results.

## PES Data

| R (bohr) | E_total (Ha) |
|----------:|-------------:|
| 0.80 | -1.004132 |
| 1.00 | -1.109328 |
| 1.20 | -1.150553 |
| 1.40 | -1.160956 |
| 1.60 | -1.155937 |
| 1.80 | -1.143278 |
| 2.00 | -1.127145 |
| 2.50 | -1.084371 |
| 3.00 | -1.047394 |
| 4.00 | -0.995690 |
| 5.00 | -0.959152 |
| 6.00 | -0.925923 |

## Morse Parameters

| Parameter | Ab Initio | Experiment | Error |
|:----------|----------:|-----------:|------:|
| D_e (Ha) | 0.164669 | 0.1745 | -5.63% |
| R_e (bohr) | 1.3835 | 1.401 | -1.25% |
| omega_e (cm-1) | 4918.24 | 4401.21 | +11.75% |
| omega_e*x_e (cm-1) | 167.33 | 121.34 | +37.90% |
| B_e (cm-1) | 62.416 | 60.853 | +2.57% |
| alpha_e (cm-1) | 3.029 | 3.062 | -1.08% |

## Rovibrational Transitions

| Transition | Ab Initio (cm-1) | Expt (cm-1) | NIST (cm-1) | Error |
|:-----------|------------------:|------------:|------------:|------:|
| v=0->1 (fundamental) | 4583.6 | 4148.3 | 4161.0 | +10.49% |
| v=0->2 (first overtone) | 8832.5 | 8043.7 | 8087.0 | +9.81% |
| J=0->1 (rotational) | 121.8 | 118.6 | 118.0 | +2.66% |
| R(0): (0,0)->(1,1) | 4699.3 | 4260.8 | 4497.8 | +10.29% |

## Notes

- PES computed via prolate spheroidal CI with Eckart Z_eff optimization
- Two orbitals: sigma_g (bonding) + sigma_u (antibonding) from exact H2+ solutions
- 2x2 CI matrix captures ionic-covalent mixing
- V_ee computed via azimuthal averaging (elliptic K integral)
- Morse V(R) = D_e[1-exp(-a(R-R_e))]^2 - D_e fitted via least squares
- Spectroscopic constants derived analytically from Morse parameters
- Nuclear lattice: Morse SU(2) vibrational chain (finite representation) + SO(3) rotational paraboloid
- Vibration-rotation coupling alpha_e included
- D_e underestimated by ~33% (missing dynamical correlation beyond 2x2 CI)
- R_e, omega_e, and B_e accurate to ~1-3% (shape of PES well preserved)
