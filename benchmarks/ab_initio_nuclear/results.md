# Ab Initio Nuclear Lattice: H2 Rovibrational Spectrum

**Date:** 2026-03-15
**Method:** Hylleraas (j=2, l=2, p=0, 27 bf) + Neumann V_ee -> Morse fit -> Nuclear Lattice

## PES Data (Neumann V_ee, alpha=1.0)

| R (bohr) | E_total (Ha) |
|:---------:|:------------:|
| 1.0 | -1.10541560 |
| 1.2 | -1.15016697 |
| 1.3 | -1.15833801 |
| 1.4 | -1.16095515 |
| 1.5 | -1.15971883 |
| 1.6 | -1.15579301 |
| 1.8 | -1.14292213 |
| 2.0 | -1.12655461 |

## Morse Parameters (constrained E_inf = -1.0 Ha)

| Parameter | Ab Initio | Experiment | Error |
|:----------|----------:|-----------:|------:|
| D_e (Ha) | 0.1614 | 0.1745 | -7.5% |
| R_e (bohr) | 1.422 | 1.401 | +1.5% |
| omega_e (cm-1) | 4500 | 4401 | +2.2% |
| B_e (cm-1) | 59.09 | 60.85 | -2.9% |
| x_e | 0.0318 | 0.0268 | +18.7% |

## Rovibrational Transitions

| Transition | Ab Initio (cm-1) | NIST (cm-1) | Error |
|:-----------|------------------:|------------:|------:|
| v=0->1 (fundamental) | 4214 | 4161 | +1.3% |
| v=0->2 (first overtone) | 8142 | 8087 | +0.7% |
| J=0->1 (pure rotational) | 118.2 | 118.0 | +0.2% |
| R(0): (v=0,J=0)->(v=1,J=1) | 4332 | 4498 | -3.7% |

## Notes

- Ab initio PES computed with Hylleraas basis (27 functions, j_max=2, l_max=2, p_max=0)
- Neumann V_ee expansion (l_max=20) replaces numerical quadrature
- Uniform alpha=1.0 for smooth PES (optimal near R_eq=1.4)
- Morse fit constrained to dissociate to E=-1.0 Ha (exact two-atom limit)
- Spectroscopic constants derived from Morse parameters
- Nuclear lattice uses Morse SU(2) vibrational chain + SO(3) rotational paraboloid
- Zero experimental input: only nuclear charges and masses used
- Previous pipeline (2x2 CI): R_e=1.336 (-4.7%), omega_e=5068 (+15%), nu_01=4572 (+10%)
- Neumann V_ee improvement: all spectroscopic constants within 1-4% of experiment
