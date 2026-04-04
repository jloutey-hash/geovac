# Track BB: geovac-hamiltonians PyPI Package

## Dependency Chain Analysis

Starting from `ecosystem_export.py`, the import tree is:

```
ecosystem_export.py
  -> composed_qubit.py -> qubit_encoding.py -> [openfermion]
  -> gaussian_reference.py -> qubit_encoding.py -> [openfermion]
  -> lattice_index.py -> lattice.py
  -> qubit_encoding.py -> [openfermion]
```

**Total: 6 geovac modules, ~9,500 lines.**

External dependencies: numpy, scipy, openfermion (required); qiskit, pennylane (optional).

The chain is clean: no module pulls in the research solvers (hyperspherical, level4, etc.).

## Decision: Standalone Bundle

The 6-module chain is small enough to bundle directly into the package with
adjusted internal imports. This makes `pip install geovac-hamiltonians` fully
standalone -- no need for the research repo.

## Bundled Modules

| Package module | Source module | Lines | Purpose |
|---------------|--------------|-------|---------|
| `_ecosystem_export.py` | `geovac/ecosystem_export.py` | 479 | Public API: `hamiltonian()`, `GeoVacHamiltonian` |
| `_composed_qubit.py` | `geovac/composed_qubit.py` | 2,763 | LiH/BeH2/H2O/H2 composed Hamiltonians |
| `_qubit_encoding.py` | `geovac/qubit_encoding.py` | 341 | JW encoding, `build_fermion_op_from_integrals()` |
| `_gaussian_reference.py` | `geovac/gaussian_reference.py` | 420 | STO-3G H2 reference |
| `_lattice_index.py` | `geovac/lattice_index.py` | 4,905 | He atomic FCI (N-electron lattice) |
| `_lattice.py` | `geovac/lattice.py` | 612 | `GeometricLattice` (atomic graph) |

## Import Adjustments

All `from geovac.*` imports were changed to `from geovac_hamiltonians._*`:

- `_composed_qubit.py`: `from geovac.qubit_encoding` -> `from geovac_hamiltonians._qubit_encoding`
- `_gaussian_reference.py`: same
- `_lattice_index.py`: `from .lattice` -> `from geovac_hamiltonians._lattice`
- `_ecosystem_export.py`: all 6 lazy imports adjusted

## Bug Found and Fixed

`_build_h2o()` in `ecosystem_export.py` passed `R=R` to `build_composed_h2o()`,
but the function signature uses `R_OH` as the parameter name. Fixed in the
bundled copy to `R_OH=R`. The original `geovac/ecosystem_export.py` still has
this bug (not modified per task constraints).

## Filesystem Caching

`_lattice_index.py` caches Slater integrals to `{__file__}/../cache/`. In the
installed package this creates a `cache/` directory inside the installed
`geovac_hamiltonians/` package. This works but is not ideal for all deployment
scenarios. A future version could use `platformdirs` or `appdirs` for
user-writable cache locations.

## Build Artifacts

```
dist/
  geovac_hamiltonians-0.1.0-py3-none-any.whl   (92 KB)
  geovac_hamiltonians-0.1.0.tar.gz              (92 KB)
```

## Test Results

29/29 tests pass:
- Import and version: 2
- H2 (bond-pair): 6
- He (atomic): 4
- LiH (composed): 6
- BeH2 (composed): 3
- H2O (composed): 3
- Qiskit export: 1
- PennyLane export: 1
- Error handling: 3

Published Pauli term counts verified:
- H2: 112
- LiH: 334
- BeH2: 556
- H2O: 778

## Future: PyPI Upload

```bash
pip install twine
twine upload dist/*
```

Requires PyPI credentials. For TestPyPI first:
```bash
twine upload --repository testpypi dist/*
pip install --index-url https://test.pypi.org/simple/ geovac-hamiltonians
```
