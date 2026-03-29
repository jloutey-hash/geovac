"""
Tests for composed H₂O qubit Hamiltonian (geovac/composed_qubit.py).

Validates orbital counts, symmetries, selection rules, block structure,
and comparison with BeH₂ and published Gaussian data.
"""

import numpy as np
import pytest

from geovac.composed_qubit import (
    build_composed_h2o,
    build_composed_beh2,
    _enumerate_states,
    GAUSSIAN_H2O_PUBLISHED,
)


# ---------------------------------------------------------------------------
# Fixture: build at small basis for fast tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def h2o_result():
    """Build composed H₂O at max_n_core=2, max_n_val=2 (default)."""
    return build_composed_h2o(max_n_core=2, max_n_val=2, verbose=False)


@pytest.fixture(scope='module')
def beh2_result_for_comparison():
    """Build composed BeH₂ at max_n_core=2, max_n_val=2 for comparison."""
    return build_composed_beh2(max_n_core=2, max_n_val=2, verbose=False)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _get_block(idx: int, blocks) -> str:
    """Identify which block an orbital index belongs to."""
    for name, rng in blocks:
        if idx in rng:
            return name
    return 'unknown'


# ---------------------------------------------------------------------------
# 1. test_h2o_orbital_count
# ---------------------------------------------------------------------------

def test_h2o_orbital_count(h2o_result):
    """M = M_core + 2×M_bond + 2×M_lone,
    where M_bond = M_bond_o + M_bond_h, M_lone = M_lone_o."""
    M = h2o_result['M']
    M_core = h2o_result['M_core']
    M_bond = h2o_result['M_bond']
    M_lone = h2o_result['M_lone']
    assert M == M_core + 2 * M_bond + 2 * M_lone
    # max_n=2: each center has 5 orbitals (1s, 2s, 2p-1, 2p0, 2p+1)
    assert M_core == 5
    assert h2o_result['M_bond_o'] == 5
    assert h2o_result['M_bond_h'] == 5
    assert M_bond == 10  # 5 + 5 per bond
    assert M_lone == 5   # O-side only
    assert M == 5 + 2 * 10 + 2 * 5  # = 35


# ---------------------------------------------------------------------------
# 2. test_h2o_qubit_count
# ---------------------------------------------------------------------------

def test_h2o_qubit_count(h2o_result):
    """Q = 2M (Jordan-Wigner spin-orbital encoding)."""
    assert h2o_result['Q'] == 2 * h2o_result['M']
    # max_n=2: Q = 70
    assert h2o_result['Q'] == 70


# ---------------------------------------------------------------------------
# 3. test_h2o_h1_hermitian
# ---------------------------------------------------------------------------

def test_h2o_h1_hermitian(h2o_result):
    """h1 matrix must be Hermitian (real symmetric for real basis)."""
    h1 = h2o_result['h1']
    assert np.allclose(h1, h1.T, atol=1e-14)


# ---------------------------------------------------------------------------
# 4. test_h2o_cross_eri_zero — 7 independent blocks
# ---------------------------------------------------------------------------

def test_h2o_cross_eri_zero(h2o_result):
    """No ERIs should mix any of the 7 blocks."""
    eri = h2o_result['eri']
    M_core = h2o_result['M_core']
    M_bond = h2o_result['M_bond']
    M_lone = h2o_result['M_lone']
    M = h2o_result['M']

    M_val_o = h2o_result['M_bond_o']

    # Define 7 block ranges
    blocks = [
        ('core', range(0, M_core)),
        ('bond1_o', range(M_core, M_core + M_val_o)),
        ('bond1_h', range(M_core + M_val_o, M_core + M_bond)),
        ('bond2_o', range(M_core + M_bond, M_core + M_bond + M_val_o)),
        ('bond2_h', range(M_core + M_bond + M_val_o, M_core + 2 * M_bond)),
        ('lone1', range(M_core + 2 * M_bond, M_core + 2 * M_bond + M_lone)),
        ('lone2', range(M_core + 2 * M_bond + M_lone, M)),
    ]

    violations = 0
    for p in range(M):
        block_p = _get_block(p, blocks)
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    if abs(eri[p, q, r, s]) > 1e-15:
                        block_q = _get_block(q, blocks)
                        block_r = _get_block(r, blocks)
                        block_s = _get_block(s, blocks)
                        if len({block_p, block_q, block_r, block_s}) > 1:
                            violations += 1

    assert violations == 0, f"Found {violations} cross-block ERI violations"


# ---------------------------------------------------------------------------
# 5. test_h2o_pauli_positive
# ---------------------------------------------------------------------------

def test_h2o_pauli_positive(h2o_result):
    """N_pauli must be a positive integer."""
    N = h2o_result['N_pauli']
    assert isinstance(N, int)
    assert N > 0


# ---------------------------------------------------------------------------
# 6. test_h2o_bond_symmetry — bond 1 and bond 2 blocks identical
# ---------------------------------------------------------------------------

def test_h2o_bond_symmetry(h2o_result):
    """Bond 1 and bond 2 ERI sub-blocks should be identical (both O-H)."""
    eri = h2o_result['eri']
    M_core = h2o_result['M_core']
    M_bond = h2o_result['M_bond']

    b1_start = M_core
    b1_end = M_core + M_bond
    eri_bond1 = eri[b1_start:b1_end, b1_start:b1_end,
                    b1_start:b1_end, b1_start:b1_end]

    b2_start = M_core + M_bond
    b2_end = M_core + 2 * M_bond
    eri_bond2 = eri[b2_start:b2_end, b2_start:b2_end,
                    b2_start:b2_end, b2_start:b2_end]

    assert np.allclose(eri_bond1, eri_bond2, atol=1e-12), \
        f"Bond 1 and Bond 2 ERI blocks differ: max diff = " \
        f"{np.max(np.abs(eri_bond1 - eri_bond2))}"


# ---------------------------------------------------------------------------
# 7. test_h2o_lone_pair_symmetry — lone pair 1 and 2 blocks identical
# ---------------------------------------------------------------------------

def test_h2o_lone_pair_symmetry(h2o_result):
    """Lone pair 1 and lone pair 2 ERI sub-blocks should be identical."""
    eri = h2o_result['eri']
    M_core = h2o_result['M_core']
    M_bond = h2o_result['M_bond']
    M_lone = h2o_result['M_lone']

    lp1_start = M_core + 2 * M_bond
    lp1_end = lp1_start + M_lone
    eri_lp1 = eri[lp1_start:lp1_end, lp1_start:lp1_end,
                  lp1_start:lp1_end, lp1_start:lp1_end]

    lp2_start = lp1_end
    lp2_end = lp2_start + M_lone
    eri_lp2 = eri[lp2_start:lp2_end, lp2_start:lp2_end,
                  lp2_start:lp2_end, lp2_start:lp2_end]

    assert np.allclose(eri_lp1, eri_lp2, atol=1e-12), \
        f"Lone pair 1 and 2 ERI blocks differ: max diff = " \
        f"{np.max(np.abs(eri_lp1 - eri_lp2))}"


# ---------------------------------------------------------------------------
# 8. test_h2o_o_blocks_identical_eri — all 4 O-side valence blocks same ERI
# ---------------------------------------------------------------------------

def test_h2o_o_blocks_identical_eri(h2o_result):
    """All 4 O-side valence blocks (bond1-O, bond2-O, lone1, lone2) should
    have identical ERIs since they all use Z_eff=6 with the same orbital set."""
    eri = h2o_result['eri']
    M_core = h2o_result['M_core']
    M_bond = h2o_result['M_bond']
    M_lone = h2o_result['M_lone']
    M_val_o = h2o_result['M_bond_o']

    # Extract O-side sub-blocks from each of the 4 positions
    offsets = [
        M_core,                          # bond1-O
        M_core + M_bond,                 # bond2-O
        M_core + 2 * M_bond,             # lone1
        M_core + 2 * M_bond + M_lone,    # lone2
    ]

    blocks = []
    for off in offsets:
        block = eri[off:off + M_val_o, off:off + M_val_o,
                    off:off + M_val_o, off:off + M_val_o]
        blocks.append(block)

    # All should match block 0
    for i in range(1, 4):
        assert np.allclose(blocks[0], blocks[i], atol=1e-12), \
            f"O-side block 0 and {i} differ: max diff = " \
            f"{np.max(np.abs(blocks[0] - blocks[i]))}"


# ---------------------------------------------------------------------------
# 9. test_h2o_more_pauli_than_beh2 — at same max_n
# ---------------------------------------------------------------------------

def test_h2o_more_pauli_than_beh2(h2o_result, beh2_result_for_comparison):
    """H₂O should have more Pauli terms than BeH₂ at same max_n
    (more blocks -> more orbitals -> more terms)."""
    assert h2o_result['N_pauli'] > beh2_result_for_comparison['N_pauli'], \
        (f"H₂O N_pauli ({h2o_result['N_pauli']}) should exceed "
         f"BeH₂ N_pauli ({beh2_result_for_comparison['N_pauli']})")


# ---------------------------------------------------------------------------
# 10. test_h2o_eri_density_lower_than_beh2 — 7 blocks vs 5
# ---------------------------------------------------------------------------

def test_h2o_eri_density_lower_than_beh2(h2o_result, beh2_result_for_comparison):
    """Total ERI density should be lower for H₂O than BeH₂ at same max_n
    because H₂O has 7 blocks vs BeH₂'s 5, creating more structural zeros."""
    h2o_density = h2o_result['ERI_density_total']
    beh2_density = beh2_result_for_comparison['ERI_density_total']
    assert h2o_density < beh2_density, \
        (f"H₂O ERI density ({h2o_density:.4%}) should be less than "
         f"BeH₂ ({beh2_density:.4%})")


# ---------------------------------------------------------------------------
# 11. test_h2o_vs_gaussian_sto3g — Q≈12 comparison
# ---------------------------------------------------------------------------

def test_h2o_vs_gaussian_sto3g():
    """At max_n=1 (Q=14, close to Gaussian STO-3G Q=12), compare Pauli counts.
    GeoVac with structural sparsity should have fewer terms."""
    result = build_composed_h2o(max_n_core=1, max_n_val=1, verbose=False)
    Q_gv = result['Q']
    N_gv = result['N_pauli']
    Q_gauss = GAUSSIAN_H2O_PUBLISHED['sto-3g']['Q']
    N_gauss = GAUSSIAN_H2O_PUBLISHED['sto-3g']['N_pauli']

    # At max_n=1: M=7 (1 core + 2x(1+1) bonds + 2x1 lone), Q=14
    assert Q_gv == 14, f"Expected Q=14 at max_n=1, got Q={Q_gv}"

    # N_pauli should be positive and less than Gaussian at similar Q
    assert N_gv > 0, f"N_pauli should be positive, got {N_gv}"
    # GeoVac at Q=14 should have far fewer than Gaussian at Q=12 (551)
    # because of block-diagonal structure. Loose bound: < 551.
    assert N_gv < N_gauss, \
        (f"GeoVac H₂O at Q={Q_gv} has {N_gv} Pauli terms, "
         f"expected fewer than Gaussian STO-3G ({N_gauss} at Q={Q_gauss})")


# ---------------------------------------------------------------------------
# Additional validation tests
# ---------------------------------------------------------------------------

def test_h2o_h1_core_eigenvalues(h2o_result):
    """Core diagonal entries should be -Z^2/(2n^2) for Z=8."""
    h1 = h2o_result['h1']
    states_core = h2o_result['states_core']
    Z = 8
    for i, (n, l, m) in enumerate(states_core):
        expected = -Z**2 / (2.0 * n**2)
        assert abs(h1[i, i] - expected) < 1e-12, \
            f"Core h1[{i},{i}] = {h1[i,i]}, expected {expected}"


def test_h2o_eri_particle_exchange_symmetry(h2o_result):
    """Check chemist notation symmetry: (pq|rs) = (rs|pq)."""
    eri = h2o_result['eri']
    M = h2o_result['M']
    # Spot-check first 10 indices (don't iterate all M^4 for M=35)
    n_checks = min(M, 10)
    for p in range(n_checks):
        for q in range(n_checks):
            for r in range(n_checks):
                for s in range(n_checks):
                    if abs(eri[p, q, r, s]) > 1e-15:
                        assert abs(eri[p, q, r, s] - eri[r, s, p, q]) < 1e-12


def test_h2o_geometry():
    """Verify R_HH computation from R_OH and angle."""
    R_OH = 1.809
    angle = 104.5
    R_HH_expected = 2.0 * R_OH * np.sin(np.radians(angle / 2.0))
    # R_HH ≈ 2.861 bohr
    assert abs(R_HH_expected - 2.861) < 0.01, \
        f"R_HH = {R_HH_expected:.4f}, expected ~2.861 bohr"


def test_h2o_max_n1():
    """max_n=1 should produce M=7 (1 core + 2x(1+1) bonds + 2x1 lone), Q=14."""
    result = build_composed_h2o(max_n_core=1, max_n_val=1, verbose=False)
    assert result['M_core'] == 1
    assert result['M_bond'] == 2   # 1 O-side + 1 H per bond
    assert result['M_lone'] == 1
    assert result['M'] == 7   # 1 + 2*2 + 2*1
    assert result['Q'] == 14
    assert result['N_pauli'] > 0


def test_h2o_no_pk():
    """All diagonals should be negative when PK is disabled."""
    result = build_composed_h2o(max_n_core=2, max_n_val=2,
                                 include_pk=False, verbose=False)
    diag = np.diag(result['h1'])
    assert np.all(diag < 0), f"Found non-negative diagonal without PK: {diag}"


def test_h2o_7_blocks():
    """Verify n_blocks = 7."""
    result = build_composed_h2o(max_n_core=1, max_n_val=1, verbose=False)
    assert result['n_blocks'] == 7


# ---------------------------------------------------------------------------
# Run module directly for quick diagnostics
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("Building composed H₂O for testing...\n")
    result = build_composed_h2o(verbose=True)

    print(f"\n--- Quick validation ---")
    print(f"M = {result['M']}, Q = {result['Q']}, N_pauli = {result['N_pauli']}")
    print(f"h1 Hermitian: {np.allclose(result['h1'], result['h1'].T)}")
    print(f"ERI density: {result['ERI_density_total']:.2%}")
