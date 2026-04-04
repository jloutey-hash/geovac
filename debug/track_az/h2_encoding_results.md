# Track AZ: H2 Bond-Pair Qubit Encoding Results

## Scaling Sweep (R = 1.4 bohr)

| max_n | M (spatial) | Q (qubits) | N_pauli | 1-norm (Ha) | QWC groups | ERI density | Wall time (s) |
|:-----:|:-----------:|:----------:|:-------:|:-----------:|:----------:|:-----------:|:-------------:|
| 2     | 5           | 10         | 112     | 8.1739      | 21         | 10.40%      | 0.4           |
| 3     | 14          | 28         | 2,627   | 42.0822     | 790        | 3.88%       | 3.6           |
| 4     | 30          | 60         | 30,955  | 133.3748    | 10,195     | 2.10%       | 29.5          |

**Scaling exponent:** N_pauli ~ 0.08 * Q^3.13

This is consistent with the He atomic scaling (Paper 14: O(Q^3.15)) as expected,
since the H2 bond-pair encoding is structurally identical to the He encoding at Z=1.

## R Sweep (max_n = 2, Q = 10)

| R (bohr) | N_pauli | 1-norm (Ha) | V_NN (Ha) |
|:---------:|:-------:|:-----------:|:---------:|
| 0.5       | 112     | 9.4596      | 2.0000    |
| 1.0       | 112     | 8.4596      | 1.0000    |
| 1.4       | 112     | 8.1739      | 0.7143    |
| 2.0       | 112     | 7.9596      | 0.5000    |
| 3.0       | 112     | 7.7930      | 0.3333    |

**Pauli term count is R-independent** (112 terms at all R values). This confirms
that the sparsity pattern is entirely determined by Gaunt integral selection rules,
and only the coefficients (1-norm) vary with R through the nuclear repulsion V_NN = 1/R.

## Encoding Details

- **Basis:** Single-center hydrogenic (n,l,m) at Z_eff = 1
- **h1:** Diagonal, -Z_eff^2 / (2n^2) = -1/(2n^2)
- **ERIs:** Gaunt integral selection rules enforce block-diagonal sparsity
- **Nuclear repulsion:** V_NN = 1/R (constant energy offset)
- **Approximations:** Two-center effects (cross-nuclear attraction, prolate
  spheroidal splitting) not included; this is the composed pipeline's
  single-block approximation
