# Track BD: H2O Composed 1-Norm Analysis

## Summary

The H2O composed 1-norm is **catastrophically inflated** by the Z^2-scaled PK pseudopotential at Z_eff=6.
The PK barrier contributes 98.7% of the total 1-norm, making the raw numbers non-competitive.
Without PK, the H2O electronic 1-norm is 360.81 Ha at Q=70, which is comparable to BeH2's 354.89 Ha (with PK).

## Computed Results

### H2O Composed Qubit Hamiltonian

| n_max | Q  | N_pauli | 1-norm (Ha) | 1-norm w/o PK (Ha) | PK fraction |
|:-----:|:--:|:-------:|:-----------:|:-------------------:|:-----------:|
| 1     | 14 | 22      | 18,913.04   | 251.49              | 98.7%       |
| 2     | 70 | 778     | 28,053.33   | 360.81              | 98.7%       |

### PK Barrier Analysis

The PK pseudopotential for H2O uses Z^2-scaled parameters from Li2+ (Paper 17):
- A = 6.93 * (8/3)^2 = 49.28 Ha
- B = 7.00 * (8/3)^2 = 49.78 Ha

The 1s PK diagonal element is **2,386.77 Ha** -- this single matrix element dominates
the entire 1-norm. For comparison, LiH's PK diagonal is ~7 Ha (A=6.93).

The Z^2 scaling of PK parameters means H2O's PK barrier is (8/3)^2 / (3/3)^2 = 7.1x
larger than LiH's, but the effect on 1-norm is multiplicative across all 4 O-side
valence blocks, producing the ~28,000 Ha total.

### Cross-System Comparison

| System | n_max | Q  | N_pauli | 1-norm (Ha) | Scaling |
|:------:|:-----:|:--:|:-------:|:-----------:|:-------:|
| LiH    | 1     | 6  | 10      | 16.01       |         |
| LiH    | 2     | 30 | 334     | 37.33       | Q^0.53  |
| BeH2   | 1     | 10 | 16      | 198.20      |         |
| BeH2   | 2     | 50 | 556     | 354.89      | Q^0.36  |
| H2O    | 1     | 14 | 22      | 18,913.04   |         |
| H2O    | 2     | 70 | 778     | 28,053.33   | Q^0.24  |

### PK Decomposition (n_max=2)

| System | With PK (Ha) | Without PK (Ha) | PK contrib (Ha) | PK % |
|:------:|:------------:|:----------------:|:----------------:|:----:|
| LiH    | 37.33        | 33.26            | 4.08             | 10.9% |
| H2O    | 28,053.33    | 360.81           | 27,692.52        | 98.7% |

### Verification

- BeH2 n_max=1 (Q=10): 198.20 Ha -- MATCHES expected 198.20 Ha
- BeH2 n_max=2 (Q=50): 354.89 Ha -- MATCHES expected 354.89 Ha
- LiH n_max=2 (Q=30): 37.33 Ha -- MATCHES expected 37.33 Ha

## Analysis

### Why H2O 1-norm is inflated

The root cause is the Z^2-scaled PK pseudopotential at Z_eff=6:

1. **PK parameters scale as Z^2**: A_PK = 6.93 * (Z/3)^2. For Z=8 (oxygen),
   this gives A=49.28, which is 7.1x larger than LiH's A=6.93.

2. **PK diagonal on 1s orbital scales as A * <1s|exp(-Br^2)|1s>**: For Z_eff=6,
   the 1s orbital is very compact (r ~ 1/(2*6) = 0.083 bohr), making the Gaussian
   PK overlap with 1s extremely large. The result is a 2,387 Ha diagonal element.

3. **4 O-side valence blocks**: H2O has 4 blocks with PK (2 bond-O + 2 lone pair),
   each contributing ~2,387 Ha to 1-norm from the 1s diagonal alone.

4. **Physical interpretation**: The PK barrier is preventing valence electrons from
   occupying core-like orbitals. At Z_eff=6, the 1s orbital has binding energy
   -Z_eff^2/2 = -18 Ha, so a 2,387 Ha barrier is ~130x the binding energy.
   This is dramatically over-repulsive.

### Scaling Exponents

The 1-norm scaling exponents are anomalously low (Q^0.24 for H2O, Q^0.36 for BeH2,
Q^0.53 for LiH). This is because the PK contribution is nearly constant (dominated
by the 1s diagonal), while the qubit count grows as n_max^2. Adding more orbitals
dilutes the PK fraction but doesn't add much to the electronic 1-norm.

### Electronic 1-norm (without PK)

Without PK, the H2O electronic 1-norm at Q=70 is 360.81 Ha, comparable to BeH2 at
354.89 Ha (with PK). This is the relevant number for quantum simulation cost if PK
can be handled differently (e.g., as a classical energy shift, or with a reformulated
pseudopotential that doesn't inflate the 1-norm).

### Recommendation for Paper 14 Table 7

The H2O 1-norm should be reported with caveats:
- Report both with-PK and without-PK values
- Note that the Z^2-scaled PK is the dominant cost driver
- The electronic 1-norm (without PK) is competitive
- The PK barrier represents a known limitation of the current Z^2 scaling approach

## Files

- `debug/track_bd/compute_h2o_1norm.py` -- computation script
- `debug/track_bd/h2o_1norm_data.json` -- machine-readable results
- `debug/track_bd/h2o_1norm_results.md` -- this file
