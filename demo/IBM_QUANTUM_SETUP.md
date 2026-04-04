# IBM Quantum VQE Demo -- Setup Guide

## Quick Start (Local Simulator)

No IBM account required. Runs entirely on your machine.

```bash
# From the project root:
python demo/ibm_quantum_demo.py
```

This runs VQE on the GeoVac H2 bond-pair Hamiltonian (10 qubits, 112 Pauli
terms) using statevector simulation with COBYLA optimization.

## Requirements

### Required packages

```bash
pip install qiskit>=2.0 qiskit-aer>=0.17 openfermion>=1.6
```

### Optional (for IBM Quantum hardware)

```bash
pip install qiskit-ibm-runtime>=0.40
```

## Usage

### Statevector simulation (default)

```bash
python demo/ibm_quantum_demo.py
```

Exact expectation values, no sampling noise. Fastest mode.

### Shot-based simulation

```bash
python demo/ibm_quantum_demo.py --shots 10000
```

Simulates finite-sampling noise with the specified number of shots.
Higher shot counts give more accurate results but are slower.

### Side-by-side Gaussian comparison

```bash
python demo/ibm_quantum_demo.py --compare-gaussian
python demo/ibm_quantum_demo.py --compare-gaussian --shots 10000
```

Runs VQE on both:
- GeoVac H2 bond-pair encoding (10 qubits, 112 Pauli terms)
- Gaussian STO-3G baseline (4 qubits, 15 Pauli terms)

### IBM Quantum hardware

```bash
# Via command-line token:
python demo/ibm_quantum_demo.py --token YOUR_IBM_QUANTUM_TOKEN

# Via environment variable:
export QISKIT_IBM_TOKEN=YOUR_IBM_QUANTUM_TOKEN
python demo/ibm_quantum_demo.py
```

**Important notes for hardware mode:**
- Requires a free IBM Quantum account (https://quantum.ibm.com)
- The GeoVac H2 encoding needs 10 qubits (available on free-tier backends)
- Each VQE iteration submits a remote job; expect significant wall time
- Use `--restarts 1 --maxiter 50` for a quick test run
- Hardware noise will degrade accuracy compared to simulator

### Tuning parameters

```bash
python demo/ibm_quantum_demo.py --restarts 10 --maxiter 2000  # more optimization
python demo/ibm_quantum_demo.py --restarts 1 --maxiter 100    # quick test
```

## What the demo does

1. **Builds H2 Hamiltonian** via `geovac.ecosystem_export.hamiltonian('H2')`
   - Bond-pair encoding at Z_eff=1, max_n=2
   - Jordan-Wigner transformation to 10-qubit Pauli Hamiltonian
   - 112 Pauli terms, 21 QWC measurement groups

2. **Exports to Qiskit** SparsePauliOp format

3. **Computes exact ground state** via full matrix diagonalization

4. **Runs VQE** with EfficientSU2 ansatz (reps=3, linear entanglement)
   - COBYLA optimizer with multiple random restarts
   - Reports energy, error (mHa), and wall time

## Expected results (statevector)

| Metric          | GeoVac H2    | Gaussian STO-3G |
|:----------------|:-------------|:----------------|
| Qubits          | 10           | 4               |
| Pauli terms     | 112          | 15              |
| QWC groups      | 21           | 5               |
| VQE error       | ~10-15 mHa (5 restarts) | < 1 mHa |
| Wall time       | ~30-60 s     | ~5-20 s         |

Note: the 10-qubit GeoVac H2 system has 80 variational parameters,
making COBYLA optimization harder than the 4-qubit Gaussian system.
Sub-mHa convergence requires 20+ restarts (`--restarts 20`).

## Troubleshooting

**ImportError: openfermion not found**
```bash
pip install openfermion
```

**ImportError: qiskit-ibm-runtime not found**
```bash
pip install qiskit-ibm-runtime
```
Only needed for hardware mode; not required for local simulator.

**VQE not converging (error > 10 mHa)**
- Increase restarts: `--restarts 10`
- Increase iterations: `--maxiter 2000`
- The 10-qubit GeoVac system has a harder optimization landscape than 4-qubit STO-3G

**Hardware jobs timing out**
- IBM Quantum free tier has queue times up to several hours
- Use `--restarts 1 --maxiter 50` for testing
- Consider using an IBM Quantum premium plan for faster access
