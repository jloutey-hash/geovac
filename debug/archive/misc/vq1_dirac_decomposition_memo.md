# VQ-1: Directional Decomposition of the Camporesi-Higuchi Dirac Operator on S^3

## Summary

Built and verified the three directional component matrices Gamma_mu (mu = x, y, z) of the spin-orbit coupling operator sigma.L in the DiracLabel (n_fock, kappa, m_j) basis at n_max = 2 (10 states).

## Construction

The Dirac operator on unit S^3 acting on spinor harmonics contains the angular operator:

    sigma.L = 2(S_x L_x + S_y L_y + S_z L_z)

where sigma_i = 2 S_i (Pauli matrices in terms of spin-1/2 operators) and L_i are orbital angular momentum components. The three directional components are:

    Gamma_x = 2 S_x L_x  ,  Gamma_y = 2 S_y L_y  ,  Gamma_z = 2 S_z L_z

These are operator products (matrix multiplication, NOT element-wise), and their sum is sigma.L.

### Method

1. **J operators**: Built J_z (diagonal, eigenvalue m_j) and J_+/- (standard ladder operators) in the 10-state DiracLabel basis.

2. **S operators**: Built S_z, S_+/- via Clebsch-Gordan decomposition. Each |j, m_j> state decomposes as sum_ms CG(l, m_j-m_s; 1/2, m_s | j, m_j) |l, m_l> |1/2, m_s>. S acts on the spin part; the CG coefficients translate between the coupled (j, m_j) and uncoupled (m_l, m_s) bases.

3. **L = J - S**: The orbital angular momentum operators are the difference.

4. **Gamma_mu = 2 S_mu @ L_mu**: Each directional component is a matrix product.

## Verification

| Check | Result |
|-------|--------|
| sigma.L = -(kappa+1) diagonal | PASSED (residual 1.67e-16) |
| Gamma_x + Gamma_y + Gamma_z = sigma.L | EXACT (residual 0.00e+00) |
| All Gamma_mu Hermitian | YES (residual < 6e-17) |
| All Gamma_mu real | YES (max |imag| = 0) |
| Tr(S_mu^2) isotropic | YES (all = 5/2) |
| Tr(L_mu^2) isotropic | YES (all = 4) |

### IMPORTANT: sigma.L + 3/2 vs full Dirac eigenvalues

sigma.L + 3/2 gives eigenvalues {-1/2, -1/2, 3/2, 3/2, 3/2, 3/2, 5/2, 5/2, 5/2, 5/2}, which are n-INDEPENDENT. They depend only on kappa.

The full Camporesi-Higuchi Dirac eigenvalues chi*(n_fock+3/2) are n-DEPENDENT: {-7/2, -7/2, 5/2, 5/2, 7/2, 7/2, 7/2, 7/2, 7/2, 7/2}.

This means sigma.L is the ANGULAR part of the Dirac operator. The full Dirac operator on S^3 involves left-invariant vector fields X_i on SU(2), not just the internal angular momentum L_i. The directional decomposition Gamma_x + Gamma_y + Gamma_z decomposes the angular/spin-orbit piece, not the full S^3 differential operator.

## Results

### Sparsity

| Matrix | Nonzero entries | Sparsity |
|--------|----------------|----------|
| Gamma_x | 16/100 | 84.0% |
| Gamma_y | 16/100 | 84.0% |
| Gamma_z | 10/100 | 90.0% |
| sigma.L | 10/100 (diagonal only) | 90.0% |

The first 4 states (n_fock=1, all have l=0/kappa=-1) contribute ZERO to all Gamma_mu. The s-states (l=0) are inert under sigma.L because L|l=0> = 0. All non-trivial structure lives in the n_fock=2 block (states 4-9, l=1 with kappa=-2 and kappa=+1).

### Eigenvalues

All three Gamma_mu have identical eigenvalue spectrum:
    {-1, -1, 0, 0, 0, 0, 0, 0, +1, +1}

The six zero eigenvalues come from: 4 from the l=0 null space (n_fock=1: 2 states, n_fock=2: 2 states with kappa=-1), plus 2 from the kernel within the l=1 block.

### Number field

All matrix entries of all three Gamma_mu are elements of Q(sqrt(2), sqrt(3)).

Distinct nonzero entry values and their algebraic forms:

| Float value | Algebraic form | Squared value |
|-------------|---------------|---------------|
| 0.2357022604 | 1/(3*sqrt(2)) = sqrt(1/18) | 1/18 |
| 0.3333333333 | 1/3 | 1/9 |
| 0.4082482905 | 1/sqrt(6) = sqrt(1/6) | 1/6 |
| 0.4714045208 | sqrt(2)/3 = sqrt(2/9) | 2/9 |
| 0.5773502692 | 1/sqrt(3) = sqrt(1/3) | 1/3 |
| 0.6666666667 | 2/3 | 4/9 |
| 1.0000000000 | 1 | 1 |

The irrational entries are exactly the CG coefficients for l=1, s=1/2 coupling. This is expected: the CG coefficients for coupling angular momentum l with spin-1/2 involve sqrt(rational) by construction.

**PI-free certificate**: All entries are algebraic (elements of Q(sqrt(2), sqrt(3))). No transcendental numbers (pi, zeta values, logarithms) appear. This is consistent with the GeoVac pi-free graph principle.

### Commutators

The three directional components do NOT commute:

| Commutator | max |entry| |
|------------|--------------|
| [Gamma_x, Gamma_y] | 0.8165 ~ sqrt(2/3) |
| [Gamma_x, Gamma_z] | 0.7071 ~ 1/sqrt(2) |
| [Gamma_y, Gamma_z] | 0.7071 ~ 1/sqrt(2) |

The non-commutativity is physical: sigma_x L_x and sigma_y L_y involve different angular momentum components that don't commute with each other.

### Traces

All Tr(Gamma_mu) = 0. This is expected since sigma.L is traceless (sum of eigenvalues -(kappa+1) over all states sums to zero).

### Block structure

The matrices are block-diagonal with respect to n_fock (no inter-shell coupling), but Gamma_x and Gamma_y have off-diagonal elements within the n_fock=2 block connecting states with different kappa (specifically kappa=-2 <-> kappa=+1, i.e. p_{3/2} <-> p_{1/2}). This is the physical j-mixing from the x and y components of spin-orbit coupling.

Gamma_z preserves m_j (as expected for a z-component), and couples kappa=-2 to kappa=+1 only when they share the same m_j. Gamma_x and Gamma_y change m_j by +/-2 (i.e. Delta m_j = +/-1) within the same kappa block, AND couple different kappa values with Delta m_j = 0 or +/-1.

## Structural observations

1. **s-wave null space**: The l=0 sector (4 states at n_max=2) is an exact null space of ALL three Gamma_mu. The Dirac decomposition acts entirely within the l >= 1 sector. This is the same Kramers cancellation that makes H_SO vanish at l=0 (Paper 14, spin_orbit.py).

2. **Isotropic spectrum**: All three Gamma_mu have the same eigenvalue spectrum {-1, 0, 1} with multiplicities {2, 6, 2}. This rotational isotropy is expected from SO(3) symmetry.

3. **Gamma_z is sparser**: 10 vs 16 nonzero entries. The z-component preserves m_j (diagonal in m_j), while x and y components mix m_j values. This makes Gamma_z the natural "easy axis" for computational purposes.

4. **CG number field**: The number field Q(sqrt(2), sqrt(3)) is the same field that appears in graph-native QED (see CLAUDE.md: "F_2_full in Q(sqrt(2))"). The CG coefficients for half-integer spin coupling generate exactly these algebraic numbers.

## Files

- `debug/vq1_dirac_decomposition.py` -- computation script
- `debug/data/vq1_dirac_decomposition.json` -- full numerical results
- `debug/vq1_dirac_decomposition_memo.md` -- this memo
