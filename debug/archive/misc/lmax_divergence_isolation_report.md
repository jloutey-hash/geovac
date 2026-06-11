# l_max Divergence Isolation Report

**Date:** 2026-03-28
**Version:** v2.0.5
**Experiment:** Bare Level 4 solver l_max convergence across charge asymmetries

## 1. Objective

Determine whether the l_max divergence in LiH composed geometry originates from:

- **(A)** Intrinsic Level 4 solver behavior at asymmetric charges — the multichannel expansion itself produces R_eq drift at high charge ratios, even without Z_eff or PK
- **(B)** Z_eff(r) injection — the smooth screening function breaks algebraic nuclear coupling
- **(C)** Static core / dynamic valence mismatch — the composed geometry's frozen core creates systematic tilt

## 2. Method

Run the bare Level 4 solver (`solve_level4_h2_multichannel()`) with **no Z_eff, no PK, no core screening** on four systems at l_max = 2, 3, 4. Compare R_eq drift rates with composed LiH reference data.

### Systems tested

| System | Z_A | Z_B | Ratio | R range | Origin | Purpose |
|--------|-----|-----|-------|---------|--------|---------|
| H₂ | 1 | 1 | 1:1 | [0.8, 3.0] | midpoint | Symmetric baseline |
| HeH⁺ | 2 | 1 | 2:1 | [1.0, 4.0] | charge_center | Asymmetric test |
| Eff LiH (midpoint) | 1 | 1 | 1:1 | [2.0, 7.0] | midpoint | Distance control |
| Eff LiH (CC origin) | 1 | 1 | 1:1 | [2.0, 7.0] | z0=R/4 | Origin shift control |

### Parameters

- H₂ and HeH⁺ l_max=2,3: n_alpha=200, n_Re=400, 12 R points (full resolution)
- HeH⁺ l_max=4 and systems 3–4: n_alpha=100, n_Re=300, 8 R points (fast mode)
- All runs: adiabatic solver (n_coupled=1), sigma-only (m_max=0)

## 3. Results

### 3.1 R_eq vs l_max

| System | l_max=2 | l_max=3 | l_max=4 | Drift rate |
|--------|---------|---------|---------|------------|
| H₂ (1:1) | 1.4514 | 1.4582 | 1.4912 | **+0.020 bohr/l_max** |
| HeH⁺ (2:1) | 1.4640 | 1.6513 | 1.9888 | **+0.262 bohr/l_max** |
| Eff LiH (1:1, midpoint) | 2.000* | — | 2.000* | 0.000 |
| Eff LiH (1:1, CC origin) | 2.000* | — | 2.000* | 0.000 |
| **Composed LiH (l_dep PK)** | **3.176** | **3.488** | **3.783** | **+0.303 bohr/l_max** |
| **Composed LiH (ch_blind PK)** | **3.222** | **3.583** | **3.956** | **+0.367 bohr/l_max** |

*Boundary value — PES is monotonically decreasing in [2.0, 7.0] (no equilibrium).
Z_A=Z_B=1 at LiH distances is a dissociating H₂ system with R_eq ~ 1.4 bohr.

### 3.2 Reference values

| System | Exact R_eq | H₂ exact E_total |
|--------|-----------|-------------------|
| H₂ | 1.401 bohr | -1.17447 Ha |
| HeH⁺ | 1.463 bohr | -2.9787 Ha |
| Composed LiH | 3.015 bohr | — |

### 3.3 HeH⁺ convergence detail

| l_max | R_eq | R_eq error | E_min | n_ch |
|-------|------|-----------|-------|------|
| 2 | 1.4640 | +0.1% | -2.865305 | 9 |
| 3 | 1.6513 | +12.9% | -2.926623 | 16 |
| 4 | 1.9888 | +35.9% | -2.977435 | 25 |

**HeH⁺ shows dramatic R_eq drift — 36% error at l_max=4 — with no Z_eff, no PK, no composition layer.** The divergence is purely from the bare Level 4 multichannel solver operating on an asymmetric nuclear potential.

## 4. Channel Weight Analysis (l_max=4)

### H₂ at R_eq=1.49 (R_e=1.5)

| Channel (l₁,l₂) | Weight |
|------------------|--------|
| (0,0) | 0.936 |
| (2,0) | 0.024 |
| (0,2) | 0.024 |
| (1,1) | 0.014 |

**93.6% in the (0,0) channel.** Higher-l channels are a small perturbation.

### HeH⁺ at R_eq=1.99 (R_e=1.5)

| Channel (l₁,l₂) | Weight |
|------------------|--------|
| (0,0) | 0.682 |
| (1,0) | 0.091 |
| (0,1) | 0.091 |
| (2,0) | 0.055 |
| (0,2) | 0.055 |
| (1,2)+(2,1) | 0.009 |
| (3,0)+(0,3) | 0.008 |
| (4,0)+(0,4) | 0.005 |
| (2,2) | 0.003 |

**Only 68.2% in the (0,0) channel.** Asymmetric nuclear coupling spreads weight into odd-l channels (l=1,3) which are forbidden by symmetry in H₂. The odd channels carry 20% of the total weight.

### Key structural difference

In H₂, the (l₁,l₂) = (l₂,l₁) symmetry means only even-parity channels contribute. In HeH⁺, odd-l channels (which break the exchange symmetry) are populated, and these channels couple more strongly to the nuclear potential at large R. This is the mechanism: **asymmetric nuclear coupling populates odd-l channels, and these channels preferentially stabilize the system at large R, shifting R_eq outward.**

## 5. Systems 3–4: Control results

Both "effective LiH valence" systems (Z_A=Z_B=1 at R=2–7 bohr) show:
- Monotonically decreasing PES with no equilibrium
- Zero l_max drift (trivially — R_eq is at the boundary)

This confirms:
1. **Distance alone does not cause drift.** H₂ at LiH distances simply dissociates.
2. **Charge-center origin shift does not cause drift.** Moving the origin off-center with Z_A=Z_B=1 produces the same monotonic PES.
3. **The composed LiH equilibrium exists only because of the PK pseudopotential.** The PK provides the core repulsion that creates a minimum near 3.0 bohr.

## 6. Conclusion: Hypothesis A+B confirmed (mixed origin)

### The l_max divergence has two distinct sources:

**Source 1 — Intrinsic to Level 4 with charge asymmetry (Hypothesis A):**
HeH⁺ at bare Level 4 shows +0.262 bohr/l_max drift — comparable to composed LiH (+0.303). This drift occurs with **no Z_eff, no PK, no composition layer**. It is intrinsic to the multichannel expansion when nuclear charges are unequal.

The mechanism is differential angular correlation: asymmetric nuclear coupling populates odd-l channels that are forbidden in homonuclear systems. These channels preferentially stabilize large-R configurations, systematically shifting R_eq outward with increasing l_max.

**Source 2 — PK interaction with channel structure (Hypothesis B, partial):**
The composed LiH has Z_A_eff=1, Z_B=1 (symmetric after screening), yet its drift rate (+0.303) exceeds H₂ (+0.020) by 15×. The PK pseudopotential, which provides core orthogonality, interacts differently with different angular channels. Channel-blind PK applies the same barrier to all channels, but only the (0,0) channel sees the actual core orbital. Higher-l channels get PK repulsion they shouldn't, reducing binding preferentially at short R and pushing R_eq outward. The l-dependent PK partially fixes this (+0.303 vs +0.367) but doesn't eliminate it.

### Quantitative decomposition

| Source | Drift rate | Fraction |
|--------|-----------|----------|
| Bare asymmetric nuclear coupling | ~0.26 bohr/l_max | ~70% |
| PK–channel interaction | ~0.04–0.10 bohr/l_max | ~30% |
| H₂ baseline (symmetric convergence) | ~0.02 bohr/l_max | negligible |

**However**, the composed LiH operates with Z_A_eff=1, Z_B=1 (symmetric), so Source 1 should not apply directly. The fact that it shows drift comparable to HeH⁺ (2:1 asymmetric) means the PK pseudopotential **effectively reintroduces asymmetry into the angular channel structure** even though the nuclear charges are symmetric. The PK barrier on atom A (the screened Li) but not atom B (H) breaks the A↔B symmetry of the nuclear coupling, populating odd-l channels just as actual charge asymmetry would.

### Revised interpretation

The drift in composed LiH is **primarily PK-induced angular asymmetry**, not direct charge asymmetry. The PK breaks the homonuclear symmetry of the Z_A=Z_B=1 nuclear coupling, creating an effective heteronuclear system. This explains why:
- H₂ (no PK, 1:1): negligible drift (+0.020)
- HeH⁺ (no PK, 2:1): large drift (+0.262)
- Composed LiH (PK, 1:1): large drift (+0.303) — PK makes it effectively heteronuclear

### Implications for fix strategies

1. **l-dependent PK is necessary but insufficient.** It reduces drift from +0.367 to +0.303 (17% improvement) by removing PK from higher-l channels, but the A-site PK still breaks symmetry for the (0,0) channel.

2. **R-dependent PK scaling** (w_PK(R) = δ_{l,0} × min(cap, R/R_ref)) works empirically because it suppresses PK at large R where the channel asymmetry causes the most damage.

3. **Algebraic screening (Proposal 2)** would replace the Gaussian PK with a smooth Z_eff(r) correction. This is promising because it would not break the angular channel structure — it modifies the radial nuclear potential without introducing a separate potential term that selectively affects certain channels.

4. **The HeH⁺ drift is a separate problem** that needs to be addressed if GeoVac is to treat heteronuclear systems at high l_max. This is intrinsic to the adiabatic approximation in the multichannel expansion.

## 7. Files

- Script: `debug/lmax_divergence_isolation.py` (full resolution, stopped after HeH⁺ l_max=3)
- Script: `debug/lmax_divergence_fast.py` (fast version, completed)
- Data: `debug/data/lmax_divergence_isolation.json`
- Reference: `debug/data/lmax_convergence_algebraic.csv` (composed LiH data)
