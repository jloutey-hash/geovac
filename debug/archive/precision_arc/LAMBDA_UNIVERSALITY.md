# Lambda Universality Test: H2 vs LiH

**Date:** March 12, 2026
**Version:** v0.9.38+
**Script:** `debug/validate_lambda_universality.py`
**Data:** `debug/data/lambda_universality_h2.txt`, `debug/data/coupled_en_lih.txt`

---

## 1. Goal

Test whether lambda_derived ≈ 0.14 (from LiH force-constant matching) is
universal across diatomic molecules. If yes, lambda is a topological constant
of the GeoVac framework. If no, the force-constant matching is a per-molecule
consistency condition, not a first-principles derivation.

---

## 2. Method

For each molecule, compute:
- k_Morse = 2a²D_e (from experimental spectroscopy)
- k_bare = d²E_bare/dR² at R_eq (from LCAO-FCI PES, lambda=0)
- k_correction/lambda = d²(Delta_E)/dR² at R_eq (from Fock-weighted correction)
- lambda_derived = (k_Morse - k_bare) / (k_correction/lambda)

Both molecules use MolecularLatticeIndex at nmax=3, vee='slater_full',
cross_nuclear_method='exact', cross_atom_vee=True.

---

## 3. Results

### 3a. Morse Parameters

| Parameter | H2 | LiH |
|:----------|:---|:----|
| r_e (bohr) | 1.401 | 3.015 |
| D_e (Ha) | 0.175 | 0.092 |
| a (bohr⁻¹) | 1.029 | 0.598 |
| k_Morse (Ha/bohr²) | 0.369 | 0.066 |
| mu (m_e) | 918.6 | 1606.0 |
| omega_e (cm⁻¹) | 4401 | 1406 |

### 3b. Force Constants

| Quantity | H2 | LiH |
|:---------|:---|:----|
| k_Morse (target) | 0.369 | 0.066 |
| k_bare (lambda=0) | **+0.493** | **-0.127** |
| k_corr/lam (probe=0.05) | +0.203 | +1.380 |
| **lambda_derived** | **-0.61** | **+0.14** |
| omega_bare (cm⁻¹) | 5084 (too high) | — (no min) |

### 3c. H2 Bare PES

The H2 bare LCAO-FCI PES already has a well-defined minimum with force constant
k_bare = 0.493 Ha/bohr² — **larger** than the experimental k_Morse = 0.369.
The Fock-weighted kinetic correction adds positive curvature (+0.203/lam),
making the PES even stiffer.

To match k_Morse, we would need **negative** lambda = -0.61, which is unphysical.

### 3d. LiH Bare PES (from previous analysis)

The LiH bare LCAO-FCI PES has no minimum at R_eq (k_bare = -0.127 < 0). The
Fock-weighted correction adds the missing stiffness, giving lambda_derived ≈ +0.14.

---

## 4. Conclusion: Lambda is NOT Universal

| Molecule | lambda_derived | Sign | Physical? |
|:---------|:---------------|:-----|:----------|
| H2 | -0.61 | Negative | No |
| LiH | +0.14 | Positive | Yes |

**Lambda_derived differs in both magnitude and sign between H2 and LiH.**
It is not a universal topological constant.

### Root Cause

The bare LCAO-FCI PES has qualitatively different errors for H2 vs LiH:

- **H2 (homonuclear):** BSSE creates an overbinding artifact, but the PES
  shape is reasonable — it has a minimum near the correct R, and the force
  constant is within 30% of experiment. The kinetic correction pushes the
  wrong direction.

- **LiH (heteronuclear):** The cross-nuclear attraction (Fourier diagonal)
  overbinds at short R, destroying the PES minimum entirely. The kinetic
  correction partially compensates, creating a minimum.

The force-constant matching formula lambda = (k_Morse - k_bare) / k_correction
is a per-molecule consistency equation that absorbs all the molecule-specific
errors of the LCAO-FCI approximation. It is NOT a derivation of a universal
constant.

### Implications

1. The force-constant matching is a **consistency condition**, not a derivation
2. Lambda cannot replace a proper treatment of BSSE and cross-nuclear attraction
3. The path to first-principles equilibrium geometry must go through fixing the
   underlying LCAO-FCI errors (BSSE, cross-nuclear), not through a universal lambda
4. The LiH result of lambda ≈ 0.14 was coincidental — it happened to compensate
   the specific pattern of errors in LiH's heteronuclear LCAO-FCI

---

## 5. What This Means for the Project

The coupled electron-nuclear framework (Paper 10) is still valuable for:
- Nuclear graph reproduces exact Morse/rovibrational spectra
- Coupling via E_total(v) = E_elec(R(v)) + E_nuc(v) is correct
- The framework correctly identifies v=0 as the ground state

But the kinetic correction parameter lambda is NOT a universal constant.
The force-constant matching is useful as a diagnostic (it quantifies how
much the bare PES curvature differs from experiment), but not as a derivation.

---

## 6. Raw Data

### H2 PES (nmax=3)

| R (bohr) | E(lam=0) | E(lam=0.02) | E(lam=0.05) | E(lam=0.10) |
|:----------|:---------|:------------|:------------|:------------|
| 0.701 | -1.186 | -1.126 | -1.037 | -0.892 |
| 0.841 | -1.269 | -1.221 | -1.151 | -1.035 |
| 0.981 | -1.322 | -1.284 | -1.228 | -1.135 |
| 1.121 | -1.351 | -1.321 | -1.277 | -1.202 |
| 1.261 | -1.363 | -1.340 | -1.305 | -1.247 |
| 1.331 | -1.366 | -1.345 | -1.313 | -1.260 |
| 1.401 | -1.367 | -1.348 | -1.319 | -1.271 |
| 1.471 | -1.365 | -1.348 | -1.322 | -1.278 |
| 1.541 | -1.362 | -1.347 | -1.323 | -1.283 |
| 1.681 | -1.350 | -1.338 | -1.319 | -1.287 |
| 1.961 | -1.308 | -1.300 | -1.288 | -1.267 |
| 2.242 | -1.253 | -1.248 | -1.240 | -1.227 |
| 2.802 | -1.121 | -1.119 | -1.117 | -1.104 |

### LiH PES (nmax=3, from coupled_en_lih.txt)

See `debug/data/coupled_en_lih.txt` for full data.
