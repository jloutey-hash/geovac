# BH Phase 0: Entanglement Entropy on the S^3 Hemispheric Wedge

**Date:** 2026-05-31  
**Sprint:** BH-Phase0 extended analysis  
**Files:** `debug/bh_phase0_extended_analysis.py`, `debug/data/bh_phase0_entanglement_entropy.json`

## Setup

The Bisognano-Wichmann (BW) theorem identifies the reduced density matrix of the vacuum state traced over the complement of a Rindler wedge with the thermal Gibbs state at the Unruh temperature. On the truncated Camporesi-Higuchi spectral triple at n_max, the hemispheric wedge of S^3 is the +1 eigenspace of the equatorial reflection R_polar (m_j -> -m_j), and the BW modular Hamiltonian K_alpha is diagonal with eigenvalues two_m_j (positive odd integers).

The entanglement entropy is:

    S_ent = -Tr(rho_W log rho_W)

where rho_W = exp(-K_alpha_W) / Z is the BW canonical KMS state at beta = 2*pi.

## Results

| n_max | dim_W | S_ent     | S / log(n) | S / n^2   |
|------:|------:|----------:|-----------:|----------:|
|     2 |     8 |  1.922212 |     2.7732 |  0.480553 |
|     3 |    20 |  2.691039 |     2.4495 |  0.299004 |
|     4 |    40 |  3.250103 |     2.3445 |  0.203131 |
|     5 |    70 |  3.688641 |     2.2919 |  0.147546 |
|     6 |   112 |  4.049083 |     2.2598 |  0.112475 |
|     7 |   168 |  4.354896 |     2.2380 |  0.088875 |
|     8 |   240 |  4.620386 |     2.2219 |  0.072194 |
|     9 |   330 |  4.854911 |     2.2096 |  0.059937 |
|    10 |   440 |  5.064916 |     2.1997 |  0.050649 |
|    11 |   572 |  5.255031 |     2.1915 |  0.043430 |
|    12 |   728 |  5.428689 |     2.1847 |  0.037699 |

Analytical and numerical values agree to machine precision (max diff 1.8e-15).

## Scaling

Best fit: **S = 1.963 * log(n_max) + 0.540, R^2 = 0.99992**.

The coefficient S / log(n_max) converges monotonically toward 2 from above (2.77 at n_max=2, 2.18 at n_max=12). The asymptotic formula is:

    S_ent -> 2 * log(n_max)   as n_max -> infinity

Competing fits and their R^2 values:

| Model              | R^2        | Verdict     |
|:-------------------|:-----------|:------------|
| S = a * log(n) + b | 0.99992    | **BEST**    |
| S = a * log(n)^2 + b * log(n) + c | 0.999995 | log-quadratic correction |
| S = c * n^alpha     | power-law alpha = 0.557 | REJECTED (not linear or quadratic) |
| S = a * n^2 + b    | 0.828      | **REJECTED** (area law) |

The area-law fit S ~ n_max^2 has R^2 = 0.83 -- strongly rejected. The entropy grows LOGARITHMICALLY with the cutoff, not quadratically.

## Area of the Boundary

The equatorial S^2 "area" in the truncated setting:

- **n_equator** (count of states with two_m_j = 1, the equator-nearest shell) = n_max * (n_max + 1), scaling as n_max^1.78.
- **dim_W** (total wedge dimension) scales as n_max^2.43.

Both grow as power laws in n_max, while S grows only logarithmically. There is no proportionality between S and any area proxy.

## Mechanism

The BW density matrix rho_W has Boltzmann weights exp(-two_m_j). Because two_m_j takes integer values 1, 3, 5, ..., the weights are EXPONENTIALLY suppressed: the k=1 (equator) shell carries 89-96% of the total probability weight across all n_max tested.

The ground shell (two_m_j = 1) has degeneracy n_max * (n_max + 1) ~ n_max^2. A nearly-degenerate distribution over D states has entropy approximately log(D). Therefore:

    S_ent ~ log(n_max^2) = 2 * log(n_max)

The finite-difference ratios delta_S / (2/n) converge to 1 from above (1.15 at n=3, 1.04 at n=12), confirming the derivative d(S)/d(n) -> 2/n.

## Honest Assessment

**Does the entanglement entropy obey a Bekenstein-Hawking area law?**

**No.** The entropy scales as S ~ 2 * log(n_max), not S ~ n_max^2. This is a LOGARITHMIC scaling, not an area law.

**Why not?** The BW modular Hamiltonian K_alpha has a linearly growing spectrum (two_m_j = 1, 3, 5, ..., 2*n_max - 1). The exponential Boltzmann suppression exp(-two_m_j) concentrates the thermal weight on the lowest shell. This means the entropy counts the LOG of the ground-shell degeneracy, not the total number of modes.

In continuum QFT, the area law S ~ A * Lambda^{d-2} arises from summing over ALL UV modes below cutoff Lambda, with each mode contributing O(1) to the entropy. On the truncated spectral triple, the exponential Boltzmann factor kills this mode-counting mechanism.

**Consistency with Paper 50.** The result S ~ 2 * log(n_max) matches the Paper 50 finding ("S(rho_W) ~ 2 * log n_max -- wedge boundary dimension via degeneracy log on lowest K_alpha shell"). This diagnostic confirms and extends that result with an 11-point panel to n_max = 12.

**Relation to Paper 51 gravity.** The BH area law S = A/(4G) in Paper 51 is derived via the spectral action on a conical defect (Sommerfeld-Cheeger geometry), NOT from the BW wedge entanglement entropy. The two computations are structurally different: Paper 51 uses the replica trick on the heat kernel, while this computation uses the BW thermal density matrix directly. The log scaling here is the finite-dimensional artifact; Paper 51's area law comes from the continuum spectral-action integration.
