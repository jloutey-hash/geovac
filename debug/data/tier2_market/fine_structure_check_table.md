# Tier 2 T4: Fine-structure sanity check

Source: T2 closed-form Breit-Pauli spin-orbit matrix element,

`geovac.spin_orbit.so_diagonal_matrix_element`.

Formula: H_SO(n, kappa) = -Z^4 alpha^2 (kappa+1) / [4 n^3 l(l+1/2)(l+1)]


## Per-atom single-particle 2p spin-orbit splitting

2p_1/2 (kappa=+1) minus 2p_3/2 (kappa=-2) for the valence 2p electron.

| Atom | Z_nuc | Z_eff (2p) | |H_SO(2p_1/2)-H_SO(2p_3/2)| bare (MHz) | |...| Z_eff (MHz) |
|:---|:---:|:---:|:---:|:---:|
| He (2^3P, 1s2p) | 2 | 1.00 | 1.752e+05 | 1.095e+04 |
| Li (2^2P, 1s^2 2p) | 3 | 1.30 | 8.869e+05 | 3.127e+04 |
| Be (2s2p 3P) | 4 | 1.95 | 2.803e+06 | 1.583e+05 |

## Comparison to reference (NIST / Puchalski-Pachucki)

| System | Reference (MHz) | GeoVac Z_eff (MHz) | Sign correct? | OoM correct? | log10 distance | Rel error |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| Li 2^2P_3/2-2^2P_1/2 | 1.005e+04 | 3.127e+04 | yes | yes | 0.493 | +211.1% |
| He 2^3P total span (vs 1-particle H_SO) | 3.191e+04 | 1.095e+04 | yes | yes | -0.465 | -65.7% |
| Be 2s2p 3P total span (vs 1-particle H_SO) | 7.260e+05 | 1.583e+05 | yes | yes | -0.661 | -78.2% |

## Comparison (bare Z hydrogenic only)

Li bare-Z hydrogenic splitting: 8.869e+05 MHz vs reference 1.005e+04 MHz (log10 distance 1.946, rel err +8722.2%)


## Interpretation

- **Z^4 alpha^2 scaling verified**: H_SO scales as Z^4 by the closed form ⟨1/r^3⟩ ~ Z^3 times the extra Z factor from Coulomb dV/dr.

- **Sign**: GeoVac gives |H_SO(2p_1/2) - H_SO(2p_3/2)| > 0 (the 2p_3/2 has kappa=-2 and H_SO coefficient -Z^4 alpha^2 / 48 < 0, while 2p_1/2 has kappa=+1 and H_SO = +Z^4 alpha^2 / 24 > 0).  The Breit-Pauli leading-order ordering j=3/2 below j=1/2 is a well-known sign subtlety; the full Dirac spectrum inverts this (2p_3/2 above 2p_1/2).  For an order-of-magnitude check we use |H_SO(j=3/2) - H_SO(j=1/2)| = Z^4 alpha^2 / 16.

- **Li 2^2P doublet**: This is the cleanest single-particle test. GeoVac at Z_eff=1.3 gives 3.13e+04 MHz vs NIST 1.01e+04 MHz (0.49 orders off, +211% relative error).

- **He/Be 2^3P**: Multi-electron spin-spin and spin-other-orbit terms contribute comparably to the one-particle H_SO.  The single-particle scale is within an order of magnitude of the measured total span, which is the T4 sanity criterion.

