# GeoVac Method Stack

Recommended method configurations for multi-electron FCI calculations.

## Current Recommended Stack

### Helium (Z=2, 2 electrons)

```python
from geovac.lattice_index import LatticeIndex

idx = LatticeIndex(
    n_electrons=2,
    max_n=5,              # max_n=4 sufficient for < 0.5%
    nuclear_charge=2,
    vee_method='slater_full',
    h1_method='exact',    # diagonal -Z^2/(2n^2): helium's physical h1
)
H = idx.assemble_hamiltonian()
```

**Performance (slater_full + exact h1), MEASURED 2026-09-25:**

| max_n | n_sd  | Error vs exact (-2.903724) |
|-------|-------|----------------------------|
| 2     | 45    | 2.43%                      |
| 3     | 378   | 2.10%                      |
| 4     | 1,770 | 1.99%                      |
| 5     | 5,995 | 1.94%                      |

Convergence is monotonic from above (variational) and slow: a fixed
k = Z hydrogenic basis has no continuum, so the residual is basis
incompleteness, not a solver defect.

**Corrected 2026-09-25.** This section previously recommended
`h1_method='hybrid'` for He and quoted 0.56% / 0.45% / 0.38% / 0.35%.
Those figures were NOT helium: `hybrid` adds a non-physical kappa*(-A) =
+1/16 hop on every lattice edge to the one-body operator, and in the k = Z
hydrogenic basis T + V_ne is exactly diagonal, so the hop is not part of the
atom's Hamiltonian.  It lowers He by ~48 mHa -- more than helium's entire
correlation energy -- which is what made the numbers look good.
`ecosystem_export.hamiltonian('He')` and `vqe_benchmark.build_geovac_he`
now use `exact`; `hybrid` survives only as an explicit legacy opt-in.

### Lithium (Z=3, 3 electrons)

```python
idx = LatticeIndex(
    n_electrons=3,
    max_n=4,              # max_n=3 sufficient for < 1.2%
    nuclear_charge=3,
    vee_method='slater_full',
    h1_method='exact',    # diagonal -Z^2/(2n^2) only
)
H = idx.assemble_hamiltonian()
```

**Performance (slater_full + exact h1):**

| max_n | n_sd   | Error vs exact | Assembly time |
|-------|--------|---------------|---------------|
| 2     | 120    | 5.04%         | 0.2s          |
| 3     | 3,276  | 1.15%         | 3.4s          |
| 4     | 34,220 | 1.10%         | 338s          |

Convergence is monotonic from above (variational). Exact h1 is
required for Li because the hybrid h1 graph off-diagonal terms
introduce non-physical correlation that causes variational
overshoot (energy dips below exact at max_n >= 3).

## Method Options

### vee_method (two-electron repulsion)

| Method        | Description                                    | Accuracy |
|---------------|------------------------------------------------|----------|
| `slater_full` | Full Slater R^k integrals + Slater-Condon rules | Best     |
| `slater`      | F0 direct Coulomb only (no exchange)           | ~5%      |
| `chordal`     | S3 chordal distance model                      | ~5-18%   |

### h1_method (one-electron Hamiltonian)

| Method   | Description                                     | Use case     |
|----------|-------------------------------------------------|--------------|
| `exact`  | Diagonal hydrogen eigenbasis only (physical)    | He, Li, ...  |
| `hybrid` | Exact diagonal + non-physical +1/16 graph hop   | LEGACY model only; not an atomic Hamiltonian |
| `graph`  | Full graph Laplacian (D-A) + node weights       | Single atoms |

## Reference Energies

| System | E_exact (Ha) | Source |
|--------|-------------|--------|
| He     | -2.9037     | NIST   |
| Li     | -7.4781     | NIST   |
