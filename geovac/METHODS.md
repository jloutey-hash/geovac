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
    h1_method='hybrid',   # exact diagonal + graph off-diagonal
)
H = idx.assemble_hamiltonian()
```

**Performance (slater_full + hybrid h1):**

| max_n | n_sd  | Error vs exact | Assembly time |
|-------|-------|---------------|---------------|
| 2     | 45    | 0.56%         | 0.2s          |
| 3     | 378   | 0.45%         | 0.2s          |
| 4     | 1,770 | 0.38%         | 2.3s          |
| 5     | 5,995 | 0.35%         | 22.5s         |

Convergence is monotonic from above (variational). Hybrid h1 works
well for He because the graph off-diagonal terms provide physical
CI mixing without variational overshoot.

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
| `hybrid` | Exact diagonal (-Z^2/2n^2) + graph off-diagonal | He (2e)      |
| `exact`  | Diagonal hydrogen eigenbasis only               | Li (3e)      |
| `graph`  | Full graph Laplacian (D-A) + node weights       | Single atoms |

## Reference Energies

| System | E_exact (Ha) | Source |
|--------|-------------|--------|
| He     | -2.9037     | NIST   |
| Li     | -7.4781     | NIST   |
