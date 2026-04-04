# GeoVac Quantum Encoding Benchmark Reference

Generated: 2026-04-01 (extracted from existing data, no new computations)
Sources: Paper 14, CHANGELOG.md, benchmarks/qubit_scaling_data.json, debug/track_as/


## 1. Single-Geometry GeoVac Atomic Encodings (He, Z=2, 2 electrons)

From Paper 14, Table 1 (Pauli terms), Table 3 (QWC), Table 4 (1-norm/Trotter), Table 5 (commutator).
Data also in benchmarks/qubit_scaling_data.json.

| System       | n_max | M  |  Q  | Pauli Terms | ERI Density | QWC Groups | 1-Norm (Ha) | lambda/Q | r(1e-3) | r(1e-6) | r_comm(1e-3) | f_anti (%) | E0 (Ha) | Error (%) |
|:-------------|:-----:|:--:|:---:|------------:|:-----------:|:----------:|:-----------:|:--------:|--------:|--------:|:------------:|:----------:|--------:|:---------:|
| GeoVac He    |   2   |  5 |  10 |         120 |    10.4%    |     25     |      11.29  |   1.129  |     253 |   7,987 |       60     |   28.0     | -2.8876 |    0.55   |
| GeoVac He    |   3   | 14 |  28 |       2,659 |     3.9%    |    791     |      78.36  |   2.799  |   1,753 |  55,410 |      303     |   21.1     | -2.8920 |    0.39   |
| GeoVac He    |   4   | 30 |  60 |      31,039 |     2.1%    | 10,199     |     261.57  |   4.360  |   5,849 | 184,959 |      821     |   12.2     | -2.8961 |    0.25   |
| GeoVac He    |   5   | 55 | 110 |     227,338 |     1.3%    |    --      |     657.07  |   5.973  |  14,693 | 464,616 |       --     |    --      | -2.8982 |    0.18   |

Scaling exponents (He, 4-point fit unless noted):

| Metric                     | Exponent | R^2    |
|:---------------------------|:--------:|:------:|
| Pauli terms                |   3.147  | 0.9995 |
| QWC groups (3 pts)         |   3.355  | 1.0000 |
| 1-norm lambda              |   1.694  | 0.9972 |
| Trotter r(1e-3)            |   1.693  | 0.9972 |
| Trotter r(1e-6)            |   1.694  | 0.9972 |
| Commutator r_comm (3 pts)  |   1.47   |   --   |


## 2. Single-Geometry GeoVac Atomic Encodings (H, Z=1, 1 electron)

From benchmarks/qubit_scaling_data.json.

| System       | n_max | M  |  Q  | Pauli Terms | ERI Density | QWC Groups | 1-Norm (Ha) | lambda/Q | r(1e-3) | r(1e-6) |
|:-------------|:-----:|:--:|:---:|------------:|:-----------:|:----------:|:-----------:|:--------:|--------:|--------:|
| GeoVac H     |   2   |  5 |  10 |         120 |    10.4%    |     25     |       7.83  |   0.783  |     176 |   5,540 |
| GeoVac H     |   3   | 14 |  28 |       2,659 |     3.9%    |    789     |      43.00  |   1.535  |     962 |  30,401 |
| GeoVac H     |   4   | 30 |  60 |      31,035 |     2.1%    | 10,199     |     136.91  |   2.282  |   3,062 |  96,811 |
| GeoVac H     |   5   | 55 | 110 |     226,462 |     1.3%    |    --      |     337.91  |   3.072  |   7,556 | 238,937 |

Scaling exponents (H, 4-point fit):

| Metric          | Exponent | R^2    |
|:----------------|:--------:|:------:|
| Pauli terms     |   3.145  | 0.9996 |
| QWC groups      |   3.355  | 1.0000 |
| 1-norm lambda   |   1.570  | 0.9993 |


## 3. Single-Geometry GeoVac Molecular Encoding (H2)

From Paper 14, Table 1. Only one data point available (LCAO encoding, not Level-4 natural geometry).

| System      | n_max/atom | M  |  Q  | Pauli Terms | ERI Density |
|:------------|:----------:|:--:|:---:|------------:|:-----------:|
| GeoVac H2   |     2      | 10 |  20 |         391 |     1.9%    |

**Note:** No H2 Level-4 (natural geometry) qubit encoding has been computed.
The H2 Level-4 solver (Paper 15) operates in mol-frame hyperspherical coordinates
and produces classical PES data (96.0% D_e at l_max=6), but has not been encoded
as a qubit Hamiltonian via JW transformation.

### H2 Classical Accuracy vs Quantum Cost Convergence (Level 4)

**MISSING:** H2 Level-4 qubit encoding at varying l_max does not exist. The Level 4
solver is a coordinate-space Hamiltonian (not second-quantized), so JW encoding
would require the "full N-electron" approach (Track AS), which was shown to be
categorically denser than composed encoding. No l_max -> qubit -> Pauli terms
-> D_e convergence table can be extracted from existing data.

Classical-only convergence from Paper 15 / CHANGELOG (no qubit data):

| l_max | Channels (sigma+pi) | D_e (%) | Source     |
|:-----:|:-------------------:|--------:|:-----------|
|   0   |        1            |  47.0   | Paper 15   |
|   2   |        7            |  84.7   | Paper 15   |
|   4   |       25            |  94.2   | Paper 15   |
|   6   |       61            |  96.0   | Paper 15   |
|  CBS  |       --            |  ~97    | Extrapolated |


## 4. Gaussian Baseline Encodings (Computed Integrals)

From Paper 14 Tables 1, 3, 4, and geovac/gaussian_reference.py.

| System          | Basis    | M  |  Q  | Pauli Terms | ERI Density | QWC Groups | 1-Norm (Ha) | FCI E (Ha) |
|:----------------|:---------|:--:|:---:|------------:|:-----------:|:----------:|:-----------:|:----------:|
| Gauss He        | STO-3G   |  1 |   2 |           4 |    100%     |      1     |       3.38  |   -2.848   |
| Gauss H2        | STO-3G   |  2 |   4 |          15 |     50%     |      5     |       1.98  |   -1.137   |
| Gauss 6-31G*    |          |  4 |   8 |         201 |     50%     |     --     |        --   |     --     |
| Gauss He        | cc-pVDZ  |  5 |  10 |         156 |    100%     |     --     |      42.95  |  -2.8877   |
| Gauss cc-pVDZ*  |          | 10 |  20 |      23,189 |    100%     |     --     |        --   |     --     |
| Gauss He        | cc-pVTZ  | 14 |  28 |      21,607 |    100%     |     --     |     530.47  |  -2.9003   |

*Entries marked * use synthetic integrals with representative sparsity (scaling-fit only, not validated).

### Equal-Qubit Head-to-Head Comparisons (He)

| Q  | GeoVac Pauli | Gaussian Pauli | Ratio | GeoVac lambda | Gaussian lambda | lambda Ratio |
|:--:|:------------:|:--------------:|:-----:|:-------------:|:---------------:|:------------:|
| 10 |     120      |    156 (pVDZ)  | 1.3x  |    11.29      |    42.95 (pVDZ) |     3.8x     |
| 28 |   2,659      | 21,607 (pVTZ)  | 8.1x  |    78.36      |   530.47 (pVTZ) |     6.8x     |


## 5. Gaussian Baseline Encodings (Published, Trenev et al. 2025)

From Paper 14, Table 6, and geovac/composed_qubit.py (GAUSSIAN_LIH_PUBLISHED, GAUSSIAN_H2O_PUBLISHED).
Jordan-Wigner with 2-qubit reduction.

### LiH (Trenev et al.)

| Basis    |  Q  | Pauli Terms |
|:---------|:---:|:------------|
| STO-3G   |  10 |         276 |
| 6-31G    |  20 |       5,851 |
| cc-pVDZ  |  36 |      63,519 |

Fitted exponent: alpha = 4.25, R^2 = 0.999

### H2O (Trenev et al.)

| Basis    |  Q  | Pauli Terms |
|:---------|:---:|:------------|
| STO-3G   |  12 |         551 |
| 6-31G    |  24 |       8,921 |
| cc-pVDZ  |  46 |     107,382 |

Fitted exponent: alpha = 3.92, R^2 > 0.999

**Note:** Trenev et al. do not publish BeH2 data. No published Gaussian BeH2 Pauli counts are available.


## 6. Composed GeoVac Encodings (LiH, BeH2, H2O)

From Paper 14, Tables 6-7, and geovac/composed_qubit.py.

### Pauli Term Counts (Paper 14, Table 6)

| System | Blocks | n_max | M  |  Q  | Pauli Terms | alpha (fit) | R^2   |
|:-------|:------:|:-----:|:--:|:---:|------------:|:-----------:|:-----:|
| LiH    |   3    |   1   |  3 |   6 |          10 |             |       |
| LiH    |   3    |   2   | 15 |  30 |         334 |             |       |
| LiH    |   3    |   3   | 42 |  84 |       7,879 |   2.50      | 0.991 |
| BeH2   |   5    |   1   |  5 |  10 |          16 |             |       |
| BeH2   |   5    |   2   | 25 |  50 |         556 |             |       |
| BeH2   |   5    |   3   | 70 | 140 |      13,131 |   2.51      | 0.991 |
| H2O    |   7    |   1   |  7 |  14 |          22 |             |       |
| H2O    |   7    |   2   | 35 |  70 |         778 |             |       |
| H2O    |   7    |   3   | 98 | 196 |      18,383 |   2.52      | 0.992 |

Universal composed scaling: alpha ~ 2.5 (exponent spread 0.02)

### Additional data point from debug/data/paper14_composed_section_notes.md:

| System | n_max | M  |  Q  | Pauli Terms | N_h1 | N_ERI  |
|:-------|:-----:|:--:|:---:|------------:|-----:|-------:|
| LiH    |   4   | 90 | 180 |      92,899 |  261 | 92,638 |

### 1-Norm and Trotter Steps (Paper 14, Table 7)

| System | n_max |  Q  | 1-Norm (Ha) | r(1e-3) | r(1e-4) |
|:-------|:-----:|:---:|:-----------:|--------:|--------:|
| LiH    |   1   |   6 |      16.01  |     359 |   1,133 |
| LiH    |   2   |  30 |      37.33  |     835 |   2,640 |
| LiH    |   3   |  84 |     202.49  |   4,528 |  14,319 |
| BeH2   |   1   |  10 |     198.20  |   4,432 |  14,015 |
| BeH2   |   2   |  50 |     354.89  |   7,936 |  25,095 |
| BeH2   |   3   | 140 |     735.60  |  16,449 |  52,015 |
| H2O    |   --  |  -- |   MISSING   |    --   |    --   |

**H2O 1-norm:** MISSING from Paper 14 Table 7. The paper caption notes H2O 1-norms
"are inflated by the Z_eff=6 Phillips-Kleinman barrier and should be interpreted
cautiously." The Pauli term counts and scaling exponent are unaffected. Referenced
in Paper 14 Sec IV.D but numerical values not tabulated.

### QWC Groups (Composed Systems)

| System | n_max |  Q  | QWC Groups | Source           |
|:-------|:-----:|:---:|:----------:|:-----------------|
| LiH    |   1   |   6 |      1     | Track AS JSON    |
| LiH    |   2   |  30 |     21     | Track AS JSON    |
| LiH    |   3   |  84 |   MISSING  | Not computed     |
| BeH2   |  all  |  -- |   MISSING  | Not computed     |
| H2O    |  all  |  -- |   MISSING  | Not computed     |

**Note:** Composed system QWC groups are only recorded for LiH in debug/track_as/track_as_results.json.
BeH2 and H2O QWC group counts are not tabulated in Paper 14 or any debug file.

### ERI Density (Composed LiH)

From debug/data/paper14_composed_section_notes.md:

| n_max | M  |  Q  | ERI Density |
|:-----:|:--:|:---:|:-----------:|
|   1   |  3 |   6 |     3.70%   |
|   2   | 15 |  30 |     0.39%   |
|   3   | 42 |  84 |     0.14%   |
|   4   | 90 | 180 |     0.08%   |


## 7. Equal-Qubit GeoVac vs Gaussian Comparison (Composed Systems)

### H2O (Paper 14, Table 8)

Gaussian values interpolated from published Q^3.92 power law.

| n_max |  Q  | GeoVac Pauli | Gaussian (interp.) | Ratio    |
|:-----:|:---:|:------------:|:-------------------:|:--------:|
|   1   |  14 |           22 |           1,126     |    51x   |
|   2   |  70 |          778 |         580,688     |   746x   |
|   3   | 196 |       18,383 |      31,457,102     | 1,712x   |

### LiH (from debug/data/paper14_composed_section_notes.md)

Gaussian values interpolated from published Q^4.25 power law.

| Q  | GeoVac Pauli | Gaussian (interp.) | Ratio    |
|:--:|:------------:|:-------------------:|:--------:|
|  6 |           10 |              33     |   3.3x   |
| 30 |          334 |          30,456     |  91.2x   |
| 84 |        7,879 |       2,423,128     | 307.5x   |
|180 |       92,899 |      61,844,756     | 665.7x   |


## 8. Full N-Electron vs Composed Comparison (LiH, Track AS)

From debug/track_as/track_as_results.json and comparison_table.txt.

| Encoding              |  Q  | dim  | N_Pauli  | N_QWC |  1-Norm (Ha) | Density  |
|:----------------------|:---:|:----:|:--------:|:-----:|:------------:|:--------:|
| Full l=1, g=4         |   7 |  128 |       45 |     4 |        98.93 |  0.27%   |
| Full l=1, g=5         |   8 |  250 |    4,069 |   725 |     1,010.36 |  6.2%    |
| Full l=1, g=6         |   9 |  432 |    7,064 |   830 |     1,756.66 |  2.7%    |
| Full l=1, g=7         |  10 |  686 |   51,907 | 8,415 |     3,762.14 |  5.0%    |
| Full l=2, g=4         |  10 |  768 |    3,288 |    48 |       738.51 |  0.31%   |
| Composed (1,1)        |   6 |  N/A |       10 |     1 |        16.01 |  0.24%   |
| Composed (1,2)        |  22 |  N/A |      226 |    21 |        28.97 |   --     |
| Composed (2,2)        |  30 |  N/A |      334 |    21 |        37.33 |   --     |
| Gaussian STO-3G       |  10 |  N/A |      276 |   -- |          --  |   --     |

Key ratios:
- At comparable Q: Full Q=10 (3,288 terms) vs Composed Q=30 (334 terms) = 10x denser at 3x fewer qubits
- 1-norm: Full (738 Ha) vs Composed (37 Ha) = 20x larger
- Verdict: Composed encoding is categorically sparser


## 9. Scaling Exponent Summary

| System/Method            | Metric       | Exponent | R^2    | Source        |
|:-------------------------|:-------------|:--------:|:------:|:-------------|
| GeoVac He (single)       | Pauli terms  |   3.15   | 0.9995 | Paper 14     |
| GeoVac He (single)       | QWC groups   |   3.36   | 1.0000 | Paper 14     |
| GeoVac He (single)       | 1-norm       |   1.69   | 0.9972 | Paper 14     |
| GeoVac He (single)       | Trotter r    |   1.69   | 0.9972 | Paper 14     |
| GeoVac He (single)       | Comm. bound  |   1.47   |   --   | Paper 14     |
| GeoVac H (single)        | Pauli terms  |   3.15   | 0.9996 | Scaling sweep|
| GeoVac H (single)        | 1-norm       |   1.57   | 0.9993 | Scaling sweep|
| GeoVac LiH (composed)    | Pauli terms  |   2.50   | 0.991  | Paper 14     |
| GeoVac BeH2 (composed)   | Pauli terms  |   2.51   | 0.991  | Paper 14     |
| GeoVac H2O (composed)    | Pauli terms  |   2.52   | 0.992  | Paper 14     |
| Gaussian LiH (Trenev)    | Pauli terms  |   4.25   | 0.999  | Trenev 2025  |
| Gaussian H2O (Trenev)    | Pauli terms  |   3.92   | >0.999 | Trenev 2025  |


## 10. Missing Data Summary

| Item                                      | Referenced In       | Status    |
|:------------------------------------------|:--------------------|:----------|
| H2 Level-4 qubit encoding (any l_max)     | Not referenced      | NOT COMPUTED - Level 4 is coordinate-space, not second-quantized |
| H2O composed 1-norm / Trotter steps       | Paper 14 Sec IV.D   | MISSING - omitted from Table 7 (inflated by Z_eff=6 PK) |
| BeH2 / H2O composed QWC groups            | Not referenced      | NOT COMPUTED |
| LiH composed QWC groups at n_max=3        | Not referenced      | NOT COMPUTED |
| He n_max=5 QWC groups                     | Paper 14 Table 3    | OMITTED (227k terms, O(N^2) grouping too slow) |
| Gaussian BeH2 published Pauli counts      | Not referenced      | NOT AVAILABLE (Trenev et al. do not report BeH2) |
| H2 LCAO encoding 1-norm / QWC             | Not referenced      | NOT COMPUTED |
| Composed LiH/BeH2/H2O commutator bounds   | Paper 14 Sec IV.D   | MENTIONED ("4-5x tightening") but no numerical table |
