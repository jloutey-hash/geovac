# Ab Initio Nuclear Lattice: H2 Rovibrational Spectrum

**Date:** 2026-03-15
**Method:** Hylleraas (j=2, l=2, p=0, 27 bf) + Neumann V_ee -> Morse fit -> Nuclear Lattice

## PES Data (Neumann V_ee, alpha optimized per R)

| R (bohr) | E_total (Ha) | D_e % | alpha_opt |
|:---------:|:------------:|:-----:|:---------:|
| 0.80 | -1.004132 | 2.4 | 0.6339 |
| 1.00 | -1.109328 | 62.7 | 0.7693 |
| 1.20 | -1.150553 | 86.3 | 0.8921 |
| 1.40 | -1.160956 | 92.2 | 1.0095 |
| 1.60 | -1.155937 | 89.4 | 1.1205 |
| 1.80 | -1.143278 | 82.1 | 1.2302 |
| 2.00 | -1.127145 | 72.9 | 1.3361 |
| 2.50 | -1.084371 | 48.4 | 1.5934 |
| 3.00 | -1.047394 | 27.2 | 1.8445 |
| 4.00 | -0.995690 | -2.5 | 2.3413 |
| 5.00 | -0.959168 | -23.4 | 2.8392 |
| 6.00 | -0.925965 | -42.4 | 2.9970 |

## Morse Parameters (constrained E_inf = -1.0 Ha, fit to R=1.0-2.5)

| Parameter | Ab Initio | Experiment | Error |
|:----------|----------:|-----------:|------:|
| D_e (Ha) | 0.1615 | 0.1745 | -7.5% |
| R_e (bohr) | 1.418 | 1.401 | +1.2% |
| omega_e (cm-1) | 4435 | 4401 | +0.8% |
| B_e (cm-1) | 59.49 | 60.85 | -2.2% |
| x_e | 0.0318 | 0.0268 | +18.7% |

## Rovibrational Transitions

| Transition | Ab Initio (cm-1) | NIST (cm-1) | Error |
|:-----------|------------------:|------------:|------:|
| v=0->1 (fundamental) | 4157 | 4161 | -0.1% |
| v=0->2 (first overtone) | 8037 | 8087 | -0.6% |
| J=0->1 (pure rotational) | 119.0 | 118.0 | +0.8% |
| R(0): (v=0,J=0)->(v=1,J=1) | 4276 | 4498 | -4.9% |

## Notes

- Ab initio PES computed with Hylleraas basis (27 functions, j_max=2, l_max=2, p_max=0)
- Neumann V_ee expansion (l_max=20) replaces numerical quadrature
- Alpha optimized at each R via golden-section minimization
- Morse fit constrained to dissociate to E=-1.0 Ha (exact two-atom limit)
- Spectroscopic constants derived from Morse parameters
- Nuclear lattice uses Morse SU(2) vibrational chain + SO(3) rotational paraboloid
- Zero experimental input: only nuclear charges and masses used
- Previous pipeline (2x2 CI): R_e=1.336 (-4.7%), omega_e=5068 (+15%), nu_01=4572 (+10%)
- Neumann V_ee improvement: fundamental nu_01 within 0.1% of experiment
