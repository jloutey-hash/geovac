# Equilibrium Geometry from Coupled Electron-Nuclear Graph

**Date:** March 12, 2026
**Version:** v0.9.38+
**Script:** `debug/validate_coupled_en.py`
**Module:** `geovac/coupled_en_lattice.py`
**Tests:** `tests/test_coupled_en.py` (22/22 passing)
**Data:** `debug/data/coupled_en_lih.txt`

---

## 1. Goal

Determine whether coupling the electron LCAO-FCI graph to the nuclear Morse
vibrational graph produces equilibrium geometry without fitted parameters.
Specifically: can the Fock-weighted kinetic correction parameter lambda be
**derived** from the nuclear graph's spectroscopic constants rather than fitted?

---

## 2. Method

### 2a. Morse expectation values

For each vibrational state v of LiH, the expectation value of the bond length
is computed analytically (Dahl & Springborg 1988):

    <R>_v = r_e - (1/a) * [psi(2*lambda - 2*v) - ln(2*lambda)]

where psi is the digamma function, a = 0.5984 bohr^-1, lambda = 28.73.

### 2b. Electronic energy at R(v)

LiH LCAO-FCI (nmax=3, 367k SDs) computed at each R = <R>_v for:
- lambda = 0.00 (bare, no kinetic correction)
- lambda = 0.02, 0.05, 0.10 (Fock-weighted kinetic correction)

### 2c. Coupled total energy

    E_total(v) = E_elec(R(v)) + E_nuc(v)

where E_nuc(v) = omega_e(v+1/2) - omega_e_xe(v+1/2)^2 is the Morse spectrum.

### 2d. Force constant matching

The Morse force constant k_Morse = 2*a^2*D_e = 0.0659 Ha/bohr^2 is the
target PES curvature at equilibrium. The corrected PES:

    E(R; lambda) = E_bare(R) + lambda * Delta_E(R)

has force constant k(lambda) = k_bare + lambda * k_correction.
Setting k(lambda) = k_Morse gives lambda_derived.

---

## 3. Results

### 3a. Morse Parameters (LiH)

| Parameter | Value | Unit |
|:----------|:------|:-----|
| r_e | 3.015 | bohr |
| D_e | 0.092 | Ha |
| a | 0.5984 | bohr^-1 |
| lambda | 28.73 | — |
| j (SU(2)) | 29.79 | — |
| v_max | 29 | — |
| k_Morse | 0.0659 | Ha/bohr^2 |

### 3b. <R>_v Expectation Values

| v | <R>_v (bohr) | E_nuc (Ha) | r_inner | r_outer |
|:--|:-------------|:-----------|:--------|:--------|
| 0 | 3.030 | 0.00318 | 2.730 | 3.359 |
| 1 | 3.089 | 0.00937 | 2.552 | 3.657 |
| 2 | 3.151 | 0.01535 | 2.443 | 3.892 |
| 5 | 3.352 | 0.03203 | 2.240 | 4.505 |
| 10 | 3.752 | 0.05559 | 2.054 | 5.525 |
| 20 | 5.054 | 0.08687 | 1.881 | 8.974 |

Key: <R>_0 = 3.030 bohr (only +0.5% from r_e = 3.015). Bond stretches
monotonically with v as expected.

### 3c. Coupled E_total(v) — All Lambda Values

**lambda = 0 (bare):**

| v | R(v) | E_nuc | E_elec | E_total | D_eff |
|:--|:-----|:------|:-------|:--------|:------|
| 0 | 3.030 | 0.003 | -8.116 | -8.113 | 0.224 |
| 1 | 3.089 | 0.009 | -8.111 | -8.102 | 0.219 |
| 5 | 3.352 | 0.032 | -8.096 | -8.064 | 0.204 |
| 10 | 3.752 | 0.056 | -8.084 | -8.029 | 0.192 |

**Minimum at v=0** — E_total monotonically increases with v.
This is expected: E_elec becomes less negative (weaker binding) as R increases,
and E_nuc increases. Both push E_total upward.

**lambda = 0.10 (Fock-weighted):**

| v | R(v) | E_nuc | E_elec | E_total | D_eff |
|:--|:-----|:------|:-------|:--------|:------|
| 0 | 3.030 | 0.003 | -8.041 | -8.037 | 0.148 |
| 1 | 3.089 | 0.009 | -8.042 | -8.033 | 0.150 |
| 5 | 3.352 | 0.032 | -8.051 | -8.019 | 0.159 |
| 7 | 3.501 | 0.042 | -8.053 | -8.011 | 0.161 |
| 10 | 3.752 | 0.056 | -8.052 | -7.997 | 0.160 |

**Minimum at v=0** — but note that E_elec now *increases* (becomes more negative)
from v=0 to v=7, showing the kinetic correction creates a repulsive wall at
short R. The total E_total still has minimum at v=0 because E_nuc(v) rises faster
than E_elec falls.

### 3d. Force Constant Matching — Lambda Derivation

| Source lambda | k_correction/lam | lambda_derived | k_total (check) |
|:-------------|:-----------------|:---------------|:----------------|
| 0.02 | 1.144 | **0.168** | 0.0659 |
| 0.05 | 1.380 | **0.139** | 0.0659 |
| 0.10 | 1.638 | **0.118** | 0.0659 |

**Key result:** k_bare = -0.127 Ha/bohr^2 (negative = no minimum in bare PES).
The Fock-weighted correction adds positive curvature that can be tuned to match
k_Morse = 0.066 Ha/bohr^2.

**lambda_derived ranges from 0.12 to 0.17** depending on which lambda's correction
shape is used. The correction is mildly nonlinear in lambda (higher lambda has
stiffer correction per unit lambda), so the derived value depends on the operating
point. The geometric mean is **lambda ~ 0.14**.

### 3e. Vibrational Frequency Check

| lambda | k (Ha/bohr^2) | omega (cm^-1) | Target |
|:-------|:-------------|:-------------|:-------|
| 0.00 | -0.127 | — (no min) | 1405.65 |
| 0.02 | -0.104 | — (no min) | 1405.65 |
| 0.05 | -0.058 | — (no min) | 1405.65 |
| 0.10 | 0.037 | 1055 | 1405.65 |

lambda = 0.10 gives omega = 1055 cm^-1 (75% of target).
The derived lambda ~ 0.14 would give omega closer to 1406 cm^-1 (by construction,
since it's tuned to match k_Morse).

---

## 4. Key Finding: Lambda CAN Be Derived

**The force constant matching formula works:**

    lambda = (k_Morse - k_bare) / k_correction_per_lambda

This **derives** lambda from:
- k_Morse = 2*a^2*D_e (nuclear spectroscopy → a, D_e)
- k_bare: computed from bare electron PES
- k_correction_per_lambda: computed from the Fock-weighted correction shape

All three quantities are computable without any fitted parameters.
The nuclear graph provides k_Morse through its spectroscopic constants.
The electron graph provides k_bare and k_correction_per_lambda through
the LCAO-FCI at multiple R values.

**lambda_derived ~ 0.14** (vs. fitted lambda = 0.10 from the Fock-weighted
correction study). The 40% discrepancy comes from:
1. The correction is nonlinear in lambda (correction shape depends on lambda)
2. The polynomial fit for k_bare is approximate near a non-existent minimum
3. nmax=3 basis truncation

---

## 5. The Circularity Question

**Is this circular?** Partially.

The nuclear graph uses experimental spectroscopic constants (omega_e, D_e, r_e)
from NIST, which already encode the equilibrium. The electron graph does NOT
know about r_e — it computes E_elec(R) from first principles at any R.

The coupling procedure says: "The electron graph provides the PES shape, the
nuclear graph provides the curvature constraint. Together they determine lambda."

This is **not** a prediction of R_eq from nothing — it's a self-consistency
condition. The true prediction of R_eq would require:
1. An R-independent electron graph (what we have)
2. A way to compute lambda without nuclear input

The force constant matching is a *consistency condition* between the two graphs,
not a first-principles derivation of lambda.

---

## 6. Path Forward

### What would make this non-circular:

1. **Derive a from the electron graph alone.** If the overlap decay rate
   dS/dR at R_eq can be related to the Morse range parameter a, then a (and
   hence k_Morse) would come from the electron graph, not nuclear spectroscopy.

2. **Multiple-molecule universality.** If the derived lambda is the SAME for
   H2, LiH, HCl, CO — i.e., a universal constant — then it's a property of the
   graph topology, not fitted per molecule. This would be the strongest result.

3. **Dimensional analysis.** The Fock-weighted correction has units of
   [energy] × [overlap]^2 × [p0]^4. If lambda = f(KINETIC_SCALE) where
   KINETIC_SCALE = -1/16 is the universal topological constant, that would
   connect to Paper 0's universal packing.

### Immediate next steps:

- Run the same analysis for H2 (simpler, faster) to check if lambda_derived
  is similar to LiH
- Check if lambda_derived * k_correction matches k_Morse for a range of nmax
  values (basis convergence)
- Investigate whether lambda is proportional to 1/(4*KINETIC_SCALE)^2 = 4

---

## 7. Summary

| Question | Answer |
|:---------|:-------|
| Does E_total(v) have minimum at v=0? | **Yes** (for all lambda values) |
| Does v=0 correspond to R_eq ~ 3.015? | **Yes** (<R>_0 = 3.030 bohr) |
| Can lambda be derived from nuclear graph? | **Yes**, via force constant matching |
| lambda_derived | **~0.14** (vs fitted 0.10) |
| Is this first-principles? | **Partially** — still uses experimental omega_e, D_e |
| What's needed for full first-principles? | Universal lambda or electron-graph-only a |
