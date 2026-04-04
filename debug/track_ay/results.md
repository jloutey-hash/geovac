# Track AY: VQE Validation Results

## H2 STO-3G Single Point (R = 1.4 bohr)

| Metric | Value |
|:-------|------:|
| Qubits | 4 |
| Pauli terms | 15 |
| Exact diag energy | -1.137285 Ha |
| VQE energy (scipy) | -1.137254 Ha |
| VQE error (scipy) | 0.031 mHa |
| Converged < 1 mHa | YES |
| VQE energy (qiskit) | -1.133959 Ha |
| VQE error (qiskit) | 3.326 mHa |

## H2 PES Scan (STO-3G, analytical integrals)

| R (bohr) | E_exact (Ha) | E_VQE (Ha) | Error (mHa) | Converged |
|:---------|:-------------|:-----------|:------------|:----------|
| 0.50 | -0.410788 | -0.410508 | 0.281 | YES |
| 0.70 | -0.846326 | -0.845433 | 0.893 | YES |
| 1.00 | -1.078970 | -1.078665 | 0.304 | YES |
| 1.20 | -1.126699 | -1.125830 | 0.869 | YES |
| 1.40 | -1.137276 | -1.137069 | 0.207 | YES |
| 1.80 | -1.110846 | -1.110787 | 0.059 | YES |
| 2.20 | -1.064739 | -1.064678 | 0.062 | YES |
| 2.80 | -1.001125 | -1.001122 | 0.003 | YES |
| 3.50 | -0.957675 | -0.952987 | 4.688 | NO |
| 5.00 | -0.934889 | -0.933419 | 1.471 | NO |

**8/10 points converged within 1 mHa.** The two failures (R=3.5, R=5.0 bohr) are at stretched geometries where the wavefunction becomes multireference and COBYLA has difficulty finding the global minimum. This is a known VQE optimizer limitation, not a Hamiltonian issue -- the exact diagonalization energies are correct at all geometries.

**Integral engine validation:** Our analytical STO-3G integral engine matches the published Szabo & Ostlund FCI energy at R=1.4 to 0.024 mHa, confirming correct implementation.

## LiH Feasibility

| Metric | Value |
|:-------|------:|
| Qubits | 30 |
| Pauli terms | 334 |
| 1-norm | 37.3342 |
| Statevector memory | 17.2 GB |
| Feasible | NO |

## Method

- Ansatz: EfficientSU2 (Ry-Rz layers + linear CX entanglement), reps=3
- Optimizer: COBYLA (scipy.optimize.minimize), maxiter=1000
- Simulator: Statevector (exact, no noise)
- Multiple random restarts: best of 10 for single point, best of 10 for PES
- Integral engine: Analytical STO-3G Gaussian integrals for H2 PES scan (validated at R=1.4 to 0.024 mHa vs Szabo & Ostlund)
- Qubit encoding: Jordan-Wigner via OpenFermion
- Export: GeoVac ecosystem_export -> Qiskit SparsePauliOp
- Also tested: qiskit-algorithms VQE (single point only) -- converges less reliably than scipy multi-restart