# Tests Directory

## Unit Tests (pytest)

| File | What it validates |
|------|-------------------|
| `test_fock_projection.py` | 10/10 symbolic proofs: Fock stereographic projection S3 -> R3 |
| `test_fock_laplacian.py` | 8/8 symbolic proofs: Laplace-Beltrami on unit S3 |
| `test_h2_energy_decomposition.py` | H2 Full CI energy decomposition and convergence |
| `test_lih_validation.py` | LiH molecular validation |
| `test_molecular_bugs.py` | Regression tests for molecular Hamiltonian bugs |
| `test_ov_scaling_rigorous.py` | O(V) scaling exponent validation |
| `test_universal_constant_origin.py` | Universal constant -1/16 derivation and origin |

## Benchmark Suites (with assertions)

These files contain validation benchmarks with `assert` statements.
Run as `pytest tests/<file>.py -v` or standalone `python tests/<file>.py`.

| File | What it validates |
|------|-------------------|
| `production_suite.py` | H2, He, H-, Li+, Be2+ energies (golden set) |
| `advanced_benchmarks.py` | Muonic hydrogen, holography, fine structure (AdS/CFT) |
| `heavy_metals.py` | Au78+/Hg79+ solver stability, Li+/Be2+ backward compat |
| `rabi_oscillation.py` | Norm conservation, Rabi period, off-resonance suppression |
| `bulk_physics_puzzles.py` | Geometric g-2, MOND potential falloff |

## Shared Fixtures

`conftest.py` provides shared pytest fixtures (hydrogen/helium solvers).

## Debug & Validation Scripts

Debug scripts live in `debug/` (not here). See `debug/` for:
- `debug_*.py` — standalone analysis scripts
- `heavy_metal_probe.py` — torsion metric exploration
- `universal_bonding.py` — bridge parameter sweeps
- `debug/data/` — generated data files
- `debug/plots/` — generated figures
