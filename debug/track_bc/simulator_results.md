# Track BC: IBM Quantum VQE Demo -- Simulator Results

Date: 2026-04-02

## Statevector Mode (default)

```
python demo/ibm_quantum_demo.py --restarts 5 --maxiter 1000
```

| Metric          | Value          |
|:----------------|:---------------|
| System          | H2 bond-pair   |
| Qubits          | 10             |
| Pauli terms     | 112            |
| 1-norm          | 8.1739 Ha      |
| QWC groups      | 21             |
| Ansatz          | EfficientSU2 (reps=3, 80 params) |
| VQE energy      | 0.215618 Ha    |
| Exact energy    | 0.202142 Ha    |
| VQE error       | 13.476 mHa    |
| Wall time       | 54.8 s         |
| Chemical acc.   | NO             |

## Gaussian STO-3G Comparison

```
python demo/ibm_quantum_demo.py --compare-gaussian --restarts 5 --maxiter 1000
```

| Metric                | GeoVac         | Gaussian STO-3G |
|:----------------------|:---------------|:----------------|
| Qubits                | 10             | 4               |
| Pauli terms           | 112            | 15              |
| 1-norm (Ha)           | 8.1739         | 1.9841          |
| QWC groups            | 21             | 5               |
| VQE energy (Ha)       | 0.215618       | -1.137254       |
| Exact energy (Ha)     | 0.202142       | -1.137285       |
| VQE error (mHa)       | 13.476         | 0.031           |
| Wall time (s)         | 54.4           | 22.8            |

## Shot-Based Mode

```
python demo/ibm_quantum_demo.py --shots 10000 --restarts 3 --maxiter 500
```

| Metric          | Value          |
|:----------------|:---------------|
| Shots           | 10,000         |
| Noise std       | 0.0817 Ha      |
| VQE energy      | 1.287892 Ha    |
| VQE error       | 1085.750 mHa  |

## Analysis

1. **Statevector VQE on GeoVac H2 (10 qubits, 80 parameters):**
   The 10-qubit bond-pair encoding creates a harder optimization landscape
   than the 4-qubit Gaussian STO-3G. With 5 COBYLA restarts, VQE achieves
   ~13 mHa error. This is consistent with known VQE scaling challenges
   for larger qubit counts. Sub-mHa convergence would require many more
   restarts (likely 20+) or a better optimizer (L-BFGS-B, SPSA, etc.).

2. **Gaussian STO-3G VQE (4 qubits, 32 parameters):**
   Converges to 0.031 mHa, confirming Track AY results. The smaller
   system has a smoother optimization landscape.

3. **Shot-based simulation:**
   With 10k shots and 1-norm = 8.17, the noise standard deviation is
   ~0.08 Ha, which overwhelms the optimization. Higher shot counts
   (100k+) would be needed, or advanced measurement techniques
   (variance reduction, importance sampling).

4. **Structural comparison:**
   GeoVac provides 112 Pauli terms at 10 qubits vs Gaussian's 15 terms
   at 4 qubits. The GeoVac advantage is at larger system sizes (LiH:
   334 terms at 30 qubits vs ~63,500 for Gaussian cc-pVDZ, a 190x
   reduction). For H2, the minimal Gaussian basis is more compact.

## Tests

5 fast Hamiltonian builder tests: all pass (10.5s).
Slow VQE tests: require `--slow` flag.
