# LiH Kinetic Repulsion Diagnostic

**Date:** 2026-03-11
**Version:** v0.9.37
**Script:** `debug/diagnose_kinetic_repulsion.py`
**Configuration:** exact+True, nmax=3

## Problem Statement

The balanced PES (exact+True and fourier+s_only) shows monotonically attractive
D_cp(R) with no equilibrium geometry. In real molecules, short-range kinetic
confinement energy creates a repulsive wall. This diagnostic identifies why
the LCAO framework fails to produce this wall.

## Diagnostic 1: Energy Decomposition vs R

| R (bohr) | T (Ha) | V_total (Ha) | η = -V/(2T) | E_total (Ha) | V_cross (Ha) | V_ee (Ha) | V_NN (Ha) |
|----------|--------|-------------|-------------|-------------|-------------|-----------|----------|
| 1.500 | -0.0915 | -8.3796 | -45.7710 | -8.47115 | -4.6894 | 4.0978 | 2.0000 |
| 1.800 | -0.1017 | -8.2469 | -40.5290 | -8.34864 | -4.0239 | 3.8799 | 1.6667 |
| 2.000 | -0.1080 | -8.1792 | -37.8535 | -8.28726 | -3.6615 | 3.7405 | 1.5000 |
| 2.200 | -0.1142 | -8.1199 | -35.5456 | -8.23412 | -3.3442 | 3.6068 | 1.3636 |
| 2.500 | -0.1223 | -8.0554 | -32.9460 | -8.17761 | -2.9536 | 3.4260 | 1.2000 |
| 2.800 | -0.1290 | -8.0095 | -31.0500 | -8.13852 | -2.6410 | 3.2694 | 1.0714 |
| 3.000 | -0.1327 | -7.9861 | -30.0986 | -8.11880 | -2.4673 | 3.1790 | 1.0000 |
| 3.015 | -0.1329 | -7.9846 | -30.0363 | -8.11752 | -2.4553 | 3.1727 | 0.9950 |
| 3.500 | -0.1689 | -7.9233 | -23.4488 | -8.09220 | -2.0550 | 2.7558 | 0.8571 |
| 4.000 | -0.1716 | -7.9078 | -23.0405 | -8.07939 | -1.8496 | 2.6651 | 0.7500 |
| 5.000 | -0.1814 | -7.8672 | -21.6846 | -8.04860 | -1.5156 | 2.5167 | 0.6000 |


### Key Findings

**Does T(R) increase as R decreases?**
**NO — T(R) does NOT increase at short R. The kinetic energy fails to create a repulsive wall.**

- T(R=1.5) = -0.0915 Ha
- T(R=3.015) = -0.1329 Ha
- T(R=5.0) = -0.1814 Ha
- ΔT(1.5→5.0) = 0.0899 Ha

**Virial ratio η = -V/(2T):**
η ranges from -45.771 to -21.685

## Diagnostic 2: Bridge Edge Weights vs R

| R (bohr) | max|H1_AB| | sum|H1_AB| | n_nonzero | bridge_deg_A | bridge_deg_B |
|----------|-----------|-----------|-----------|-------------|-------------|
| 1.500 | 0.027896 | 0.360821 | 14 | 0.4124 | 0.4124 |
| 1.800 | 0.026387 | 0.340140 | 14 | 0.3887 | 0.3887 |
| 2.000 | 0.025016 | 0.321843 | 14 | 0.3678 | 0.3678 |
| 2.200 | 0.023466 | 0.301392 | 14 | 0.3444 | 0.3444 |
| 2.500 | 0.020969 | 0.268728 | 14 | 0.3071 | 0.3071 |
| 2.800 | 0.018430 | 0.235745 | 14 | 0.2694 | 0.2694 |
| 3.000 | 0.016779 | 0.214394 | 14 | 0.2450 | 0.2450 |
| 3.015 | 0.016658 | 0.212826 | 14 | 0.2432 | 0.2432 |
| 3.500 | 0.012974 | 0.165380 | 14 | 0.1890 | 0.1890 |
| 4.000 | 0.009769 | 0.124284 | 14 | 0.1420 | 0.1420 |
| 5.000 | 0.005217 | 0.066176 | 14 | 0.0756 | 0.0756 |


### Key Findings

- Bridge coupling **decays monotonically** with R (STO overlap × conformal factors)
- At R=1.5: max|H1_AB| = 0.027896 Ha
- At R=5.0: max|H1_AB| = 0.005217 Ha
- Bridge edges add kinetic cost via the degree matrix D, but this is **tiny**
  compared to intra-atom kinetic contributions

## Diagnostic 3: H₂ Comparison

**H₂ HAS an equilibrium** near R ≈ 1.2 bohr. This means the kinetic repulsion issue is **heteronuclear-specific** or at least more severe for LiH.

## Diagnostic 4: Intra-Atom Kinetic Independence

The intra-atom off-diagonal H1 (kinetic hopping within Li or H lattice) is
**completely R-independent**. The graph Laplacian edges within each atom are
fixed regardless of inter-nuclear distance.

This means:
- The kinetic energy of each atom's electrons is **frozen** at its isolated-atom value
- When atoms approach, there is no kinetic penalty for orbital compression
- The only R-dependent kinetic contribution comes from bridge edges (tiny)

## Diagnostic 5: Missing Physics

### What creates kinetic repulsion in real molecules

1. **Orthogonalization kinetic energy:** When AOs overlap, Löwdin or Schmidt
   orthogonalization introduces additional nodes in the MOs, raising kinetic energy.
   This goes as ~S²(R) where S is the overlap integral.

2. **Pauli kinetic repulsion:** Antisymmetry forces electrons into higher-momentum
   states when core orbitals overlap. For Li-H, the Li 1s core overlaps with H 1s
   at short R, pushing electrons into antibonding states with higher T.

3. **Kinetic integral T_AB:** In standard QC, the kinetic energy matrix element
   between AOs on different centers, <φ_A|-½∇²|φ_B>, has the correct R-dependence
   to create repulsion. Our bridge edges decay monotonically instead.

### What our framework has

- **Fixed intra-atom graph Laplacian:** Kinetic hopping within each atom is
  R-independent. No orbital compression effect.
- **Monotonically decaying bridges:** Bridge coupling S(R)·Ω_A·Ω_B → 0 as R→0
  would approach a constant, but the conformal factors actually INCREASE at small R,
  partially compensating. Net effect: bridges don't provide repulsion.
- **Cross-nuclear attraction:** This IS R-dependent and creates binding. But without
  compensating kinetic repulsion, it produces monotonic attraction.

## Root Cause

**The graph Laplacian kinetic energy is R-independent for each atom's internal
structure.** In real quantum mechanics, bringing atoms together forces their
electrons to occupy orthogonal states with higher kinetic energy. The discrete
graph Laplacian does not capture this overlap-dependent kinetic cost because:

1. Each atom's adjacency matrix (connectivity) is **fixed** at construction time
2. The degree matrix D (diagonal kinetic contribution) only counts edges, which
   don't change as R decreases
3. Bridge edges contribute to D but are too weak to create sufficient repulsion
4. There is no overlap matrix S_AB that would modify the kinetic operator

## Proposed Fix Paths

### Path A: Overlap-dependent kinetic correction (perturbative)

Add a correction term to the molecular Hamiltonian:

```
H1_corrected = H1 + T_correction(R)

T_correction = lam * Sum_{a in A, b in B} S^2(a,b) * |a><a| * (positive scale)
```

where S(a,b) is the STO overlap between orbital a on atom A and orbital b on
atom B. This would add a positive (repulsive) diagonal correction that grows
as R decreases (overlap increases).

**Pros:** Simple, perturbative, preserves existing architecture
**Cons:** Requires calibration parameter λ, not derived from first principles

### Path B: R-dependent adjacency (modify graph topology)

Make the intra-atom adjacency matrix depend on R:

```
A_intra(R) = A_intra(inf) + dA(R)
```

where δA(R) represents additional "virtual edges" induced by orbital overlap.
When orbitals on different atoms overlap, this effectively adds new paths in
the graph, increasing the degree and thus the kinetic energy.

**Pros:** More physically motivated (overlap creates new connectivity)
**Cons:** Breaks the isolated-atom eigenvalue structure

### Path C: Lowdin kinetic energy (explicit orthogonalization)

Compute the overlap matrix S_AB between atoms, perform Lowdin orthogonalization
S^(-1/2), and add the resulting kinetic energy shift:

```
T_Lowdin = Tr[rho * S^(-1/2) * T * S^(-1/2)] - Tr[rho * T]
```

**Pros:** Exact treatment, standard QC approach
**Cons:** Expensive, may conflict with graph topology

### Path D: Bond sphere approach (Paper 8)

The bond sphere theory puts both atoms on a single S3, which naturally includes
the kinetic coupling between centers. The SO(4) Wigner D-matrix elements already
encode the correct R-dependent kinetic contribution. However, v0.9.16–v0.9.18
tests showed this approach has its own issues (Fourier diagonal overbinding).

## Output Files

- `debug/KINETIC_DIAGNOSTIC.md` — this report
- `debug/plots/lih_kinetic_vs_R.png` — T(R), V(R), virial ratio plots
- `debug/data/lih_kinetic_diagnostic.txt` — raw data
