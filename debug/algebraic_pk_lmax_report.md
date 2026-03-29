# Algebraic PK l_max Convergence Report

**Date:** 2026-03-27 16:06
**n_alpha:** 100
**n_Re:** 300
**R_grid:** 18 points over [2.0, 7.0] bohr
**Reference R_eq:** 3.015 bohr (experiment)

## 1. Full Data Table

| PK Mode | l_max | n_ch | R_eq (bohr) | R_eq error (%) | E_min (Ha) | Wall time (s) |
|:--------|:-----:|:----:|:-----------:|:--------------:|:----------:|:-------------:|
| channel_blind | 2 | 9 | 3.2216 | +6.9% | -8.184791 | 466 |
| l_dependent | 2 | 9 | 3.1761 | +5.3% | -8.202393 | 479 |
| algebraic | 2 | 9 | 1.7931 | -40.5% | -8.058932 | 473 |
| channel_blind | 3 | 16 | 3.5832 | +18.8% | -8.206271 | 1063 |
| l_dependent | 3 | 16 | 3.4878 | +15.7% | -8.232460 | 611 |
| algebraic | 3 | 16 | 4.7150 | +56.4% | -8.084130 | 485 |
| channel_blind | 4 | 25 | 3.9557 | +31.2% | -8.255525 | 1121 |
| l_dependent | 4 | 25 | 3.7829 | +25.5% | -8.282225 | 1141 |
| algebraic | 4 | 25 | 5.1870 | +72.0% | -8.134306 | 1287 |

## 2. Drift Rate Comparison

Linear fit: R_eq = slope x l_max + intercept

| PK Mode | Slope (bohr/l_max) | R^2 | l_max range | n_points |
|:--------|:------------------:|:---:|:-----------:|:--------:|
| channel_blind | +0.367 | 1.000 | 2-4 | 3 |
| l_dependent | +0.303 | 1.000 | 2-4 | 3 |
| algebraic | +1.697 | 0.852 | 2-4 | 3 |

**Reference drift rates (Paper 17 / diagnostic arc):**
- channel_blind: +0.23 bohr/l_max (l_max 2-5, diagnostic arc)
- l_dependent: +0.168 bohr/l_max (l_max 2-7, Paper 17)
- r_dependent: +0.160 bohr/l_max (l_max 2-7, Paper 17)

## 3. Assessment

**DRAMATICALLY WORSE.** The algebraic PK projector does not eliminate the l_max
divergence — it amplifies it by 5.6x relative to l_dependent (drift +1.697 vs
+0.303 bohr/l_max). The behavior is qualitatively different from Gaussian PK:

### 3.1. PES Shape Anomaly

At **l_max=2**, the algebraic PK produces a **monotonically decreasing** PES
(no minimum within [2.0, 7.0] bohr). The Morse fit finds R_eq = 1.79 bohr by
extrapolation. This means the projector barrier is far too weak at l_max=2 —
the valence electrons are partially collapsing into the core.

At **l_max=3** and **l_max=4**, a shallow minimum appears but at very large R
(4.72 and 5.19 bohr respectively). The PES well depth is only ~0.005 Ha
(l_max=3) and ~0.006 Ha (l_max=4) — about 10x shallower than the Gaussian PK
PES (~0.07 Ha).

### 3.2. Energy Comparison

The algebraic PK E_min values (-8.06 to -8.13 Ha) are ~0.1–0.15 Ha above
the Gaussian PK values (-8.18 to -8.28 Ha), indicating the rank-1 projector
provides much less total repulsion than the Gaussian barrier.

### 3.3. Root Cause Analysis

The algebraic projector `V_PK = E_shift × |core⟩⟨core|` acts only in the
(0,0) channel subspace (by construction, since the core is mapped to
(l1,l2)=(0,0)). This is correct in principle. However:

1. **Barrier strength**: The Gaussian PK `V(r) = A exp(-Br²)/r²` produces a
   barrier that extends over the full alpha grid and has large matrix elements
   at small alpha (where core density peaks). The rank-1 projector has only
   one nonzero eigenvalue (E_shift ≈ 6.93 Ha), but its effective barrier
   depends on the overlap of the valence wavefunction with the core state. If
   the valence state has low overlap with the core (due to being delocalized
   at molecular geometries), the repulsion is negligible.

2. **R_e scaling**: The projector is injected as `H += R_e × E_shift × |c⟩⟨c|`.
   The R_e factor means the barrier grows linearly with hyperradius, but the
   core state |c⟩ is fixed at a single representative R. At large R_e, the
   core wavefunction may not represent the actual core structure, weakening
   the projector's effectiveness.

3. **Non-monotonic behavior**: The jump from R_eq = 1.79 (l_max=2) to 4.72
   (l_max=3) suggests the projector interacts pathologically with the l_max=2
   basis. At l_max=2 (9 channels), the (0,0) channel weight is large enough
   that the weak projector still provides some repulsion. At l_max=3 (16
   channels), the added channels dilute the (0,0) weight, and the projector
   becomes ineffective.

### 3.4. Conclusion

The algebraic rank-1 PK projector in its current form is **not viable** for
the LiH composed geometry. The atomic-limit approximation (mapping the Level 3
core state into only the (0,0) channel of the Level 4 basis) produces a
projector that is too weak to prevent core collapse. This should be added to
the dead-ends table in CLAUDE.md.

Possible fixes to investigate:
- Scale the energy_shift to match the Gaussian PK's integrated barrier strength
- Use a self-consistent Level 4 core solve (option 3 from the audit) instead
  of the atomic-limit approximation
- Add a finite-width regularization to the rank-1 projector

## 4. PES Data

### channel_blind, l_max=2

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.112076 |
| 2.250 | -8.141540 |
| 2.500 | -8.162555 |
| 2.700 | -8.172341 |
| 2.844 | -8.179321 |
| 2.989 | -8.183565 |
| 3.133 | -8.184270 |
| 3.278 | -8.185283 |
| 3.422 | -8.183847 |
| 3.567 | -8.181422 |
| 3.711 | -8.177216 |
| 3.856 | -8.171429 |
| 4.000 | -8.166241 |
| 4.500 | -8.142358 |
| 5.125 | -8.109468 |
| 5.750 | -8.074304 |
| 6.375 | -8.040715 |
| 7.000 | -8.008980 |

### l_dependent, l_max=2

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.126168 |
| 2.250 | -8.157216 |
| 2.500 | -8.179985 |
| 2.700 | -8.190854 |
| 2.844 | -8.197908 |
| 2.989 | -8.201961 |
| 3.133 | -8.202484 |
| 3.278 | -8.202702 |
| 3.422 | -8.200522 |
| 3.567 | -8.197152 |
| 3.711 | -8.192090 |
| 3.856 | -8.185661 |
| 4.000 | -8.179401 |
| 4.500 | -8.152712 |
| 5.125 | -8.116517 |
| 5.750 | -8.079466 |
| 6.375 | -8.044413 |
| 7.000 | -8.011630 |

### algebraic, l_max=2

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.056214 |
| 2.250 | -8.048290 |
| 2.500 | -8.038122 |
| 2.700 | -8.029670 |
| 2.844 | -8.023783 |
| 2.989 | -8.018348 |
| 3.133 | -8.012865 |
| 3.278 | -8.007853 |
| 3.422 | -8.003504 |
| 3.567 | -7.999212 |
| 3.711 | -7.995338 |
| 3.856 | -7.991843 |
| 4.000 | -7.988229 |
| 4.500 | -7.976750 |
| 5.125 | -7.961998 |
| 5.750 | -7.945789 |
| 6.375 | -7.928820 |
| 7.000 | -7.911130 |

### channel_blind, l_max=3

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.114376 |
| 2.250 | -8.145301 |
| 2.500 | -8.168617 |
| 2.700 | -8.180866 |
| 2.844 | -8.190033 |
| 2.989 | -8.196788 |
| 3.133 | -8.200179 |
| 3.278 | -8.204171 |
| 3.422 | -8.205945 |
| 3.567 | -8.206838 |
| 3.711 | -8.206026 |
| 3.856 | -8.203627 |
| 4.000 | -8.201881 |
| 4.500 | -8.189380 |
| 5.125 | -8.168682 |
| 5.750 | -8.142664 |
| 6.375 | -8.115782 |
| 7.000 | -8.088688 |

### l_dependent, l_max=3

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.133047 |
| 2.250 | -8.166645 |
| 2.500 | -8.193048 |
| 2.700 | -8.207655 |
| 2.844 | -8.217559 |
| 2.989 | -8.224602 |
| 3.133 | -8.228381 |
| 3.278 | -8.231704 |
| 3.422 | -8.232734 |
| 3.567 | -8.232534 |
| 3.711 | -8.230662 |
| 3.856 | -8.227533 |
| 4.000 | -8.224221 |
| 4.500 | -8.207323 |
| 5.125 | -8.180851 |
| 5.750 | -8.151514 |
| 6.375 | -8.121998 |
| 7.000 | -8.093078 |

### algebraic, l_max=3

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.060372 |
| 2.250 | -8.055837 |
| 2.500 | -8.051614 |
| 2.700 | -8.050325 |
| 2.844 | -8.050998 |
| 2.989 | -8.053141 |
| 3.133 | -8.056162 |
| 3.278 | -8.059875 |
| 3.422 | -8.064131 |
| 3.567 | -8.068260 |
| 3.711 | -8.072170 |
| 3.856 | -8.075791 |
| 4.000 | -8.078624 |
| 4.500 | -8.084262 |
| 5.125 | -8.081718 |
| 5.750 | -8.071803 |
| 6.375 | -8.057137 |
| 7.000 | -8.039671 |

### channel_blind, l_max=4

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.126789 |
| 2.250 | -8.161676 |
| 2.500 | -8.189910 |
| 2.700 | -8.206576 |
| 2.844 | -8.219204 |
| 2.989 | -8.229506 |
| 3.133 | -8.236260 |
| 3.278 | -8.243623 |
| 3.422 | -8.248681 |
| 3.567 | -8.252633 |
| 3.711 | -8.254643 |
| 3.856 | -8.254848 |
| 4.000 | -8.255523 |
| 4.500 | -8.250081 |
| 5.125 | -8.235845 |
| 5.750 | -8.214342 |
| 6.375 | -8.190782 |
| 7.000 | -8.166003 |

### l_dependent, l_max=4

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.146910 |
| 2.250 | -8.185534 |
| 2.500 | -8.217897 |
| 2.700 | -8.237515 |
| 2.844 | -8.251063 |
| 2.989 | -8.261677 |
| 3.133 | -8.268755 |
| 3.278 | -8.275213 |
| 3.422 | -8.279222 |
| 3.567 | -8.281740 |
| 3.711 | -8.282369 |
| 3.856 | -8.281533 |
| 4.000 | -8.280334 |
| 4.500 | -8.269630 |
| 5.125 | -8.248972 |
| 5.750 | -8.223839 |
| 6.375 | -8.197442 |
| 7.000 | -8.170719 |

### algebraic, l_max=4

| R (bohr) | E_composed (Ha) |
|:--------:|:---------------:|
| 2.000 | -8.074616 |
| 2.250 | -8.072540 |
| 2.500 | -8.070745 |
| 2.700 | -8.071595 |
| 2.844 | -8.073837 |
| 2.989 | -8.077655 |
| 3.133 | -8.082305 |
| 3.278 | -8.087786 |
| 3.422 | -8.093955 |
| 3.567 | -8.100032 |
| 3.711 | -8.105920 |
| 3.856 | -8.111604 |
| 4.000 | -8.116446 |
| 4.500 | -8.128894 |
| 5.125 | -8.134365 |
| 5.750 | -8.131081 |
| 6.375 | -8.121748 |
| 7.000 | -8.108398 |
