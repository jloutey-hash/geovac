# geovac-hamiltonians

Qubit Hamiltonians from the GeoVac library, exportable to OpenFermion, Qiskit
and PennyLane.

This package is a thin front door to the main `geovac` package: its
`hamiltonian()` *is* `geovac.ecosystem_export.hamiltonian`, so it always
returns what the library and the papers measure. (Before 0.3.0 it carried its
own copy of the builder, which had fallen behind onto a retired integral rule.
Do not use 0.2.0 or earlier.)

## Install

The package is not on PyPI; install it straight from the repository (this also
installs `geovac`):

```bash
pip install "geovac-hamiltonians @ git+https://github.com/jloutey-hash/geovac.git#subdirectory=geovac-hamiltonians"
```

Optional exports: add `[qiskit]`, `[pennylane]` or `[all]` after the package
name.

## Quick start

```python
from geovac_hamiltonians import hamiltonian

H = hamiltonian('LiH')          # 30 qubits
print(H.n_qubits, H.n_terms, H.one_norm)

of_op = H.to_openfermion()
qk_op = H.to_qiskit()           # requires qiskit
pl_op = H.to_pennylane()        # requires pennylane
```

## Systems and numbers

`hamiltonian()` covers the 37-system library: main-group hydrides of the first
four rows, LiF, CO, N2, F2, NaCl, the first-row transition-metal hydrides,
SrH, BaH, He and H2. Parameters: `R` (bond length in bohr; default the
experimental equilibrium), `max_n` (default 2), `tapered`.

Pauli counts, 1-norms and the comparison against Gaussian bases are reported,
with their conventions and caveats, in Papers 14 and 20 of the
[GeoVac repository](https://github.com/jloutey-hash/geovac). Those
comparisons are at matched qubit count, not matched accuracy.

## License

MIT
