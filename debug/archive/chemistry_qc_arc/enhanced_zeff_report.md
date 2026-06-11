# Enhanced Z_eff Experiment Report

**Date:** 2026-03-28
**Version:** v2.0.5
**Experiment:** PK-free composed geometry with enhanced Z_eff

## 1. Objective

Test whether encoding core exclusion within Z_eff(r) — rather than as a separate PK barrier — preserves angular channel symmetry and reduces the l_max divergence in composed LiH.

The hypothesis: a modified Z_eff(r) that becomes repulsive near the core radius creates an exclusion zone through the same nuclear coupling machinery (Gaunt integrals, split-region expansion), avoiding the angular symmetry breaking caused by the PK barrier.

## 2. Method

### 2.1 Enhanced Z_eff construction

```
Z_eff_enhanced(r) = Z_eff_standard(r) - A_rep * exp(-r^2 / (2*sigma^2))
```

where:
- `sigma = r_core` = 0.2712 bohr (from ab initio PK `<1/r^2>`-weighted radius)
- `A_rep` determined by energy matching (see below)
- `Z_eff_standard(r)` = CoreScreening Z_eff from Li 1s^2 (goes 3 → 1)

### 2.2 Energy matching

A_rep is set so the integrated repulsive potential energy matches the PK barrier for a hydrogenic 2s valence orbital:

```
integral [A_rep * exp(-r^2/(2s^2)) / r] |R_2s|^2 r^2 dr
  = integral [V_PK(r)] |R_2s|^2 r^2 dr
```

**Result:** A_rep = 40.88 (energy-matched), giving Z_eff(0) = -37.88.

This extreme value reflects a fundamental shape mismatch: the PK barrier (A*exp(-Br^2)/r^2) has 1/r^2 singularity concentrating repulsion at the nucleus, while Z_eff modification produces only 1/r repulsion. Much larger amplitude is needed to match the integrated energy.

A moderate value A_rep = 4.0 (Z_eff(0) = -1) was also tested.

### 2.3 Enhanced Z_eff curves

| r (bohr) | Z_eff_standard | Z_eff_moderate | Z_eff_energy_matched |
|:--------:|:--------------:|:--------------:|:--------------------:|
| 0.00 | 3.000 | -1.000 | -37.884 |
| 0.10 | 2.963 | -0.775 | -35.235 |
| 0.20 | 2.804 | -0.244 | -28.349 |
| 0.30 | 2.553 | +0.383 | -19.626 |
| 0.50 | 2.002 | +1.271 | -5.474 |
| 1.00 | 1.236 | +1.231 | +1.190 |
| 2.00 | 1.007 | +1.007 | +1.007 |
| 5.00 | 1.000 | +1.000 | +1.000 |

Both modes converge to Z_eff = 1 beyond r ~ 1.5 bohr.

### 2.4 Solver configuration

- Level 4 multichannel solver, adiabatic approximation (n_coupled=1)
- Z_A_func = enhanced Z_eff, Z_B = 1.0, no PK
- Z_A_bare = Z_eff_enhanced(0) → **heteronuclear** (all l1+l2 parities included)
- n_alpha = 80, n_Re = 250, n_theta = 48 (quadrature)
- R grid: 8 points over [2.0, 7.0] bohr
- E_composed = E_core + E_elec(Z_A_func) + V_NN(3,1,R) + V_cross_nuc(R)

### 2.5 Key structural difference from standard composed LiH

The standard composed solver uses constant Z_A = 1 (fully screened) + PK barrier. The enhanced Z_eff uses the full screening function Z_eff(r) (3 → 1) modified with a repulsive dip. This means:

1. **Penetration effect:** Valence electrons "see through" core screening at short range (Z_eff ~ 2-3 at r < 0.5 bohr), making them more tightly bound
2. **Heteronuclear basis:** Z_A_bare ≠ Z_B → all (l1,l2) parities included, doubling the channel count at each l_max
3. **Energy scale shift:** E_elec is ~0.1-0.3 Ha more negative than standard composed, due to the penetration effect

## 3. Results

### 3.1 PES scans — Moderate mode (A_rep = 4, Z_eff(0) = -1)

| l_max | n_ch | R_eq | E_min (Ha) | Equilibrium? |
|:-----:|:----:|:----:|:----------:|:------------:|
| 2 | 9 | 2.000* | -8.515 | No (monotonic) |
| 4 | 25 | 2.000* | -8.536 | No (monotonic) |

*Boundary value — PES monotonically decreasing at all R.

**Conclusion:** A_rep = 4 provides insufficient core repulsion. The repulsive dip reaches only Z_eff = -1 at the nucleus, but the penetration effect (Z_eff ~ 2.5 at r = 0.2) more than compensates, making the electrons overbind at all distances.

### 3.2 PES scans — Energy-matched mode (A_rep = 40.88, Z_eff(0) = -37.9)

| l_max | n_ch | R_eq (bohr) | R_eq error | E_min (Ha) | Equilibrium? |
|:-----:|:----:|:-----------:|:----------:|:----------:|:------------:|
| 2 | 9 | 2.000* | -33.7% | -8.113 | No (monotonic) |
| 3 | 16 | ~2.5 | -17.1% | -8.185 | Yes (shallow) |
| 4 | 25 | 2.280 | -24.4% | -8.190 | Yes (shallow) |

*Boundary value.

**PES data (energy-matched, l_max=4):**

| R (bohr) | E_elec (Ha) | E_composed (Ha) |
|:--------:|:-----------:|:---------------:|
| 2.000 | -1.3590 | -8.182456 |
| 2.500 | -1.2661 | **-8.189657** |
| 3.000 | -1.1811 | -8.171344 |
| 3.500 | -1.1089 | -8.146766 |
| 4.000 | -1.0486 | -8.122136 |
| 4.500 | -0.9977 | -8.098979 |
| 5.500 | -0.9148 | -8.056567 |
| 7.000 | -0.8198 | -8.000481 |

The minimum is very shallow: D_e ~ 0.007 Ha from R=2.5 to R=3.0 (vs ~0.07 Ha for standard PK).

### 3.3 Drift rate comparison

| Mode | l=2 | l=3 | l=4 | Drift (bohr/l_max) |
|:-----|:---:|:---:|:---:|:------------------:|
| Enhanced Z_eff (moderate) | 2.000* | — | 2.000* | 0.000 (trivial) |
| Enhanced Z_eff (energy-matched) | 2.000* | ~2.5 | 2.280 | —† |
| l-dependent PK (ref) | 3.176 | 3.488 | 3.783 | +0.303 |
| Channel-blind PK (ref) | 3.222 | 3.583 | 3.956 | +0.367 |
| HeH+ bare (ref) | 1.464 | 1.651 | 1.989 | +0.262 |

*Monotonic PES, minimum at boundary.
†Drift not well-defined: l_max=2 is monotonic (no equilibrium), equilibrium appears only at l_max≥3. At l_max=3→4, R_eq moves from ~2.5 to ~2.3 (negative drift), which is opposite to the PK behavior.

### 3.4 Channel weights at l_max=4

| Channel | Moderate (R=2.0) | Energy-matched (R=2.28) | HeH+ bare (R=1.99, ref) | Composed PK (ref) |
|:--------|:----------------:|:-----------------------:|:-----------------------:|:-----------------:|
| (0,0) | 92.3% | 66.5% | 68.2% | — |
| (1,0) | — | 14.0% | 9.1% | — |
| (0,1) | — | 14.0% | 9.1% | — |
| (0,3)+(3,0) | — | 2.6% | 1.6% | — |
| (1,1) | 1.3% | 1.2% | — | — |
| (2,0)+(0,2) | 5.8% | 1.0% | 11.0% | — |
| **Odd-l total** | **0.5%** | **30.7%** | **20%** | ~0%‡ |

‡Standard composed LiH uses homonuclear channels (Z_A=Z_B=1), so odd-l channels are excluded by construction. The PK breaks symmetry within the even-parity subspace.

## 4. Analysis

### 4.1 Why the enhanced Z_eff approach fails

The experiment reveals a **fundamental incompatibility** between encoding both screening and exclusion in a single Z_eff(r) function:

1. **Shape mismatch:** The PK barrier (1/r^2) concentrates repulsion at the nucleus far more effectively than a Z_eff modification (1/r). To match the PK's energy, A_rep must be ~41, making Z_eff(0) = -38. This is not a parameter-tuning issue — it reflects the different singularity orders of the two approaches.

2. **Penetration-exclusion coupling:** The standard Z_eff(r) encodes core screening (3 → 1). Adding a repulsive dip creates a function that must simultaneously describe:
   - Nuclear attraction at r > 1 bohr (Z_eff ≈ 1)
   - Enhanced attraction at r ≈ 0.3-1.0 bohr (Z_eff ≈ 1.5-2.5, penetration)
   - Strong repulsion at r < 0.3 bohr (Z_eff < 0, exclusion)

   The penetration zone (Z_eff > 1) overbinds the electrons, partially negating the repulsive dip. The standard composed approach avoids this by using constant Z_A = 1 (no penetration) + PK (clean barrier).

3. **Symmetry destruction:** Any non-trivial Z_eff modification makes Z_A_bare ≠ Z_B, switching from homonuclear to heteronuclear channel structure. This introduces odd-l channels (30.7% population), creating MORE angular asymmetry than the PK approach (which operates within the even-parity subspace). The hypothesis that Z_eff would preserve channel symmetry is **refuted** — it produces worse symmetry breaking than PK.

4. **Shallow well:** The energy-matched mode produces a very shallow equilibrium (D_e ~ 0.007 Ha, 10× shallower than PK), and only at l_max ≥ 3. The minimum is not robust — it depends on the balance between penetration attraction and repulsive dip, which shifts with l_max.

### 4.2 Comparison with the decoupled architecture

The standard composed architecture (constant Z_A + PK) succeeds because it **decouples** two distinct physical effects:

| Effect | Standard approach | Enhanced Z_eff approach |
|:-------|:-----------------|:-----------------------|
| Core screening | Constant Z_A_eff = Z - N_core | Smooth Z_eff(r) |
| Core exclusion | PK barrier (Gaussian/r^2) | Repulsive dip in Z_eff |
| Penetration | Not included | Included (unintentionally) |
| Channel symmetry | Preserved (homonuclear) | Destroyed (heteronuclear) |
| Barrier shape | 1/r^2 (concentrated) | 1/r (diffuse) |

The decoupled architecture's advantage is that each piece handles one physical effect cleanly. The Z_eff approach mixes screening, penetration, and exclusion in a single function, losing control over each.

### 4.3 Implications for the l_max divergence problem

This experiment confirms that the l_max divergence cannot be fixed by modifying the nuclear potential alone. The PK pseudopotential operates in a fundamentally different mathematical space (the angular channel subspace) from the nuclear potential. Encoding PK-like physics in Z_eff(r) changes the wrong thing — it modifies the radial nuclear potential at every channel simultaneously, rather than selectively acting on the channels that overlap the core.

The l-dependent PK remains the most promising approach: it acts only on l=0 channels (which actually overlap the 1s^2 core), leaving higher-l channels unaffected. The R-dependent PK scaling (w_PK(R)) further reduces the asymmetry at large R.

## 5. Conclusion

**The enhanced Z_eff approach is not viable for replacing the PK pseudopotential.**

Root causes:
1. 1/r^2 vs 1/r shape mismatch requires extreme A_rep values
2. Penetration effect (Z_eff > 1) overwhelms the repulsive dip
3. Any non-trivial Z_eff modification destroys homonuclear channel symmetry
4. Odd-l channel population (30.7%) exceeds even the bare HeH+ case (20%)

**This should be added to the dead-ends table in CLAUDE.md.**

## 6. Files

- Script: `debug/enhanced_zeff_quick.py` (focused experiment)
- Script: `debug/enhanced_zeff_experiment.py` (full version, not completed)
- Data: `debug/data/enhanced_zeff_experiment.json`
- This report: `debug/enhanced_zeff_report.md`
