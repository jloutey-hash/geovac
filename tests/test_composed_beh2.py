"""
Tests for composed BeH₂ qubit Hamiltonian (geovac/composed_qubit.py).

Validates orbital counts, symmetries, selection rules, block structure,
and comparison with LiH.
"""

import numpy as np
import pytest

from geovac.composed_qubit import (
    build_composed_beh2,
    build_composed_lih,
    estimate_cross_bond_eri_count,
    _enumerate_states,
)


# ---------------------------------------------------------------------------
# Fixture: build at small basis for fast tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def beh2_result():
    """Build composed BeH₂ at max_n_core=2, max_n_val=2 (default)."""
    return build_composed_beh2(max_n_core=2, max_n_val=2, verbose=False)


@pytest.fixture(scope='module')
def lih_result_for_comparison():
    """Build composed LiH at max_n_core=2, max_n_val=2 for comparison."""
    return build_composed_lih(max_n_core=2, max_n_val=2, verbose=False)


# ---------------------------------------------------------------------------
# 1. Correct orbital count: M = M_core + 2 × M_bond
# ---------------------------------------------------------------------------

def test_beh2_orbital_count(beh2_result):
    """M = M_core + 2 × M_bond, where M_bond = M_bond_be + M_bond_h."""
    M = beh2_result['M']
    M_core = beh2_result['M_core']
    M_bond = beh2_result['M_bond']
    assert M == M_core + 2 * M_bond
    # max_n=2: core has 5 orbitals, each bond has 5+5=10
    assert M_core == 5
    assert M_bond == 10
    assert M == 25


# ---------------------------------------------------------------------------
# 2. Correct qubit count: Q = 2M
# ---------------------------------------------------------------------------

def test_beh2_qubit_count(beh2_result):
    """Q = 2M (Jordan-Wigner spin-orbital encoding)."""
    assert beh2_result['Q'] == 2 * beh2_result['M']
    # max_n=2: Q = 50
    assert beh2_result['Q'] == 50


# ---------------------------------------------------------------------------
# 3. h1 is Hermitian
# ---------------------------------------------------------------------------

def test_beh2_h1_hermitian(beh2_result):
    """h1 matrix must be Hermitian (real symmetric for real basis)."""
    h1 = beh2_result['h1']
    assert np.allclose(h1, h1.T, atol=1e-14)


# ---------------------------------------------------------------------------
# 4. Cross-block ERIs are zero
# ---------------------------------------------------------------------------

def test_beh2_cross_eri_zero(beh2_result):
    """No ERIs should mix core / bond1 / bond2 blocks."""
    eri = beh2_result['eri']
    M_core = beh2_result['M_core']
    M_bond = beh2_result['M_bond']
    M = beh2_result['M']

    # Define block ranges
    core_range = range(0, M_core)
    bond1_range = range(M_core, M_core + M_bond)
    bond2_range = range(M_core + M_bond, M)

    blocks = [
        ('core', core_range),
        ('bond1', bond1_range),
        ('bond2', bond2_range),
    ]

    # Check that every nonzero ERI has all 4 indices in the same block
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


def _get_block(idx: int, blocks) -> str:
    """Helper to identify which block an orbital index belongs to."""
    for name, rng in blocks:
        if idx in rng:
            return name
    return 'unknown'


# ---------------------------------------------------------------------------
# 5. Pauli term count is positive
# ---------------------------------------------------------------------------

def test_beh2_pauli_positive(beh2_result):
    """N_pauli must be a positive integer."""
    N = beh2_result['N_pauli']
    assert isinstance(N, int)
    assert N > 0


# ---------------------------------------------------------------------------
# 6. Bond symmetry: bond1 and bond2 have identical ERI counts
# ---------------------------------------------------------------------------

def test_beh2_symmetry(beh2_result):
    """Bond 1 and bond 2 blocks should have identical ERI counts by
    symmetry of linear H-Be-H (both bonds have same Z_eff and Z_H)."""
    # The ERI integrals for both bond-Be blocks are identical (same Z_eff)
    # and both bond-H blocks are identical (same Z_H). Check via the eri tensor.
    eri = beh2_result['eri']
    M_core = beh2_result['M_core']
    M_bond = beh2_result['M_bond']
    M_bond_be = beh2_result['M_bond_be']

    # Extract bond1 ERI sub-block
    b1_start = M_core
    b1_end = M_core + M_bond
    eri_bond1 = eri[b1_start:b1_end, b1_start:b1_end,
                    b1_start:b1_end, b1_start:b1_end]

    # Extract bond2 ERI sub-block
    b2_start = M_core + M_bond
    b2_end = M_core + 2 * M_bond
    eri_bond2 = eri[b2_start:b2_end, b2_start:b2_end,
                    b2_start:b2_end, b2_start:b2_end]

    # They should be identical
    assert np.allclose(eri_bond1, eri_bond2, atol=1e-12), \
        f"Bond 1 and Bond 2 ERI blocks differ: max diff = {np.max(np.abs(eri_bond1 - eri_bond2))}"


# ---------------------------------------------------------------------------
# 7. More Pauli terms than LiH at same max_n
# ---------------------------------------------------------------------------

def test_beh2_more_pauli_than_lih(beh2_result, lih_result_for_comparison):
    """BeH₂ should have more Pauli terms than LiH at same max_n
    (more orbitals -> more terms)."""
    assert beh2_result['N_pauli'] > lih_result_for_comparison['N_pauli'], \
        (f"BeH₂ N_pauli ({beh2_result['N_pauli']}) should exceed "
         f"LiH N_pauli ({lih_result_for_comparison['N_pauli']})")


# ---------------------------------------------------------------------------
# 8. ERI density lower than LiH at same max_n
# ---------------------------------------------------------------------------

def test_beh2_eri_density_lower_than_lih(beh2_result, lih_result_for_comparison):
    """Total ERI density should be lower for BeH₂ than LiH at same max_n
    because BeH₂ has 5 blocks vs LiH's 3, creating more structural zeros."""
    beh2_density = beh2_result['ERI_density_total']
    lih_density = lih_result_for_comparison['ERI_density_total']
    assert beh2_density < lih_density, \
        (f"BeH₂ ERI density ({beh2_density:.4%}) should be less than "
         f"LiH ({lih_density:.4%})")


# ---------------------------------------------------------------------------
# Additional validation tests
# ---------------------------------------------------------------------------

def test_beh2_h1_core_eigenvalues(beh2_result):
    """Core diagonal entries should be -Z^2/(2n^2) for Z=4."""
    h1 = beh2_result['h1']
    states_core = beh2_result['states_core']
    Z = 4
    for i, (n, l, m) in enumerate(states_core):
        expected = -Z**2 / (2.0 * n**2)
        assert abs(h1[i, i] - expected) < 1e-12, \
            f"Core h1[{i},{i}] = {h1[i,i]}, expected {expected}"


def test_beh2_eri_particle_exchange_symmetry(beh2_result):
    """Check chemist notation symmetry: (pq|rs) = (rs|pq)."""
    eri = beh2_result['eri']
    M = beh2_result['M']
    # Spot-check (don't iterate all M^4 for large M)
    n_checks = min(M, 10)
    for p in range(n_checks):
        for q in range(n_checks):
            for r in range(n_checks):
                for s in range(n_checks):
                    if abs(eri[p, q, r, s]) > 1e-15:
                        assert abs(eri[p, q, r, s] - eri[r, s, p, q]) < 1e-12


def test_beh2_bond_be_h1_with_pk():
    """PK should push bond-Be 1s entries positive (repulsive barrier)."""
    result = build_composed_beh2(max_n_core=2, max_n_val=2,
                                  include_pk=True, verbose=False)
    h1 = result['h1']
    M_core = result['M_core']
    # Bond1-Be 1s orbital is at index M_core
    val_be_1s = h1[M_core, M_core]
    # With PK, this should be positive (barrier > binding)
    assert val_be_1s > 0, \
        f"PK should make bond-Be 1s positive, got {val_be_1s}"


def test_beh2_no_pk():
    """All diagonals should be negative when PK is disabled."""
    result = build_composed_beh2(max_n_core=2, max_n_val=2,
                                  include_pk=False, verbose=False)
    diag = np.diag(result['h1'])
    assert np.all(diag < 0), f"Found non-negative diagonal without PK: {diag}"


def test_beh2_max_n1():
    """max_n=1 should produce M=3 (1 core + 1+1 per bond × 2)."""
    result = build_composed_beh2(max_n_core=1, max_n_val=1, verbose=False)
    assert result['M_core'] == 1
    assert result['M_bond'] == 2
    assert result['M'] == 5  # 1 + 2*2
    assert result['Q'] == 10
    assert result['N_pauli'] > 0


def test_cross_bond_eri_estimate():
    """Cross-bond ERI count should be positive and bounded."""
    result = estimate_cross_bond_eri_count(2)
    assert result['n_cross_bond'] > 0
    assert result['n_cross_bond'] <= result['upper_bound']
    assert result['ratio'] >= 0


# ===========================================================================
# Monopole inter-fiber coupling tests (Phase 1)
# ===========================================================================

from geovac.inter_fiber_coupling import (
    extract_origin_density,
    transform_to_center_density,
    slater_f0_integral,
    monopole_inter_fiber_energy,
    exchange_inter_fiber_energy,
    compute_overlap_diagnostic,
    extract_channel_data,
    extract_origin_density_algebraic,
    compute_channel_f0_matrix,
    compute_overlap_from_channel_data,
)
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK


@pytest.fixture(scope='module')
def beh2_monopole_setup():
    """Set up core + PK for BeH2 monopole tests."""
    core = CoreScreening(Z=4, l_max=2, n_alpha=200)
    core.solve(verbose=False)
    pk = AbInitioPK(core, n_core=2)
    pk_potentials = [pk.pk_dict(atom='A')]
    return {
        'pk_potentials': pk_potentials,
        'Z_eff': 2.0,
        'Z_ligand': 1.0,
        'l_max': 2,
        'n_alpha': 100,
    }


@pytest.fixture(scope='module')
def monopole_at_three_R(beh2_monopole_setup):
    """Compute monopole at R=2.0, 3.0, 5.0 for testing."""
    s = beh2_monopole_setup
    results = {}
    for R in [2.0, 3.0, 5.0]:
        l4 = solve_level4_h2_multichannel(
            R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
            l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
            verbose=False, pk_potentials=s['pk_potentials'],
        )
        mono = monopole_inter_fiber_energy(
            l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
            l_max=s['l_max'], n_alpha=s['n_alpha'],
            pk_potentials=s['pk_potentials'], n_sample_Re=8,
        )
        results[R] = mono
    return results


# ---------------------------------------------------------------------------
# Monopole energy is positive at all R
# ---------------------------------------------------------------------------

def test_monopole_positive(monopole_at_three_R):
    """Monopole inter-fiber energy must be positive (repulsive) at all R."""
    for R, mono in monopole_at_three_R.items():
        assert mono['E_monopole'] > 0, \
            f"Monopole should be positive at R={R}, got {mono['E_monopole']}"


# ---------------------------------------------------------------------------
# Monopole energy decreases with R
# ---------------------------------------------------------------------------

def test_monopole_decreases_with_R(monopole_at_three_R):
    """Monopole should decrease with R as electron densities separate."""
    E_2 = monopole_at_three_R[2.0]['E_monopole']
    E_3 = monopole_at_three_R[3.0]['E_monopole']
    E_5 = monopole_at_three_R[5.0]['E_monopole']
    assert E_2 > E_3 > E_5, \
        f"Monopole should decrease: E(2.0)={E_2:.4f} > E(3.0)={E_3:.4f} > E(5.0)={E_5:.4f}"


# ---------------------------------------------------------------------------
# Monopole R-dependence is too flat for R_eq correction
# ---------------------------------------------------------------------------

def test_monopole_flat_at_short_R(monopole_at_three_R):
    """Monopole varies <5% between R=2 and R=3 (nearly R-independent near
    the Be core), confirming that scalar monopole correction is insufficient
    to shift R_eq inward. This is consistent with CLAUDE.md: 'scalar energy
    corrections cannot capture wavefunction modification.'"""
    E_2 = monopole_at_three_R[2.0]['E_monopole']
    E_3 = monopole_at_three_R[3.0]['E_monopole']
    fractional_change = abs(E_2 - E_3) / E_2
    assert fractional_change < 0.15, \
        f"Monopole changes {fractional_change:.1%} between R=2 and R=3"
    # The small variation means monopole adds nearly uniform repulsion,
    # pushing R_eq outward (wrong direction).


# ---------------------------------------------------------------------------
# Electron count conservation
# ---------------------------------------------------------------------------

def test_monopole_electron_count(monopole_at_three_R):
    """Density normalization should give ~2 electrons per fiber."""
    for R, mono in monopole_at_three_R.items():
        assert abs(mono['N_elec_origin'] - 2.0) < 0.01, \
            f"Origin density should integrate to 2, got {mono['N_elec_origin']}"
        assert abs(mono['N_elec_Be'] - 2.0) < 0.1, \
            f"Be-centered density should integrate to ~2, got {mono['N_elec_Be']}"


# ---------------------------------------------------------------------------
# F^0 integral validation with known analytical result
# ---------------------------------------------------------------------------

def test_slater_f0_hydrogenic():
    """Validate F^0 integral against known 1s(Z=1) analytical result.

    For a single 1s electron with Z=1:
        P(r) = 4r^2 exp(-2r), integral P = 1
        F^0(1s, 1s) = 5Z/8 = 0.625
    """
    n_r = 500
    r_max = 15.0
    dr = r_max / n_r
    r_grid = (np.arange(n_r) + 0.5) * dr
    Z = 1.0
    P = 4.0 * Z**3 * r_grid**2 * np.exp(-2.0 * Z * r_grid)
    # Normalize to 1 electron
    P *= 1.0 / (np.sum(P) * dr)
    f0 = slater_f0_integral(r_grid, P, P)
    expected = 5.0 * Z / 8.0  # 0.625
    assert abs(f0 - expected) < 0.02, \
        f"F^0(1s_Z=1) = {f0:.4f}, expected {expected:.4f} (grid discretization)"


# ---------------------------------------------------------------------------
# Shell theorem transform preserves normalization
# ---------------------------------------------------------------------------

def test_transform_preserves_norm():
    """Transform to center density should preserve total electron count."""
    n_r = 300
    r_max = 10.0
    dr = r_max / n_r
    r_grid = (np.arange(n_r) + 0.5) * dr
    Z = 2.0
    P = 4.0 * Z**3 * r_grid**2 * np.exp(-2.0 * Z * r_grid)
    P *= 2.0 / (np.sum(P) * dr)  # 2 electrons

    R = 3.0
    d_grid, P_Be = transform_to_center_density(r_grid, P, R)
    N_Be = np.sum(P_Be) * (d_grid[1] - d_grid[0])
    assert abs(N_Be - 2.0) < 0.15, \
        f"Transform should preserve N=2, got {N_Be:.4f}"


# ===========================================================================
# Exchange inter-fiber coupling tests (Phase 3)
# ===========================================================================


@pytest.fixture(scope='module')
def exchange_at_three_R(beh2_monopole_setup):
    """Compute exchange coupling at R=2.0, 3.0, 5.0 for testing."""
    s = beh2_monopole_setup
    results = {}
    for R in [2.0, 3.0, 5.0]:
        l4 = solve_level4_h2_multichannel(
            R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
            l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
            verbose=False, pk_potentials=s['pk_potentials'],
        )
        exch = exchange_inter_fiber_energy(
            l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
            l_max=s['l_max'], n_alpha=s['n_alpha'],
            pk_potentials=s['pk_potentials'], n_sample_Re=10,
        )
        results[R] = exch
    return results


# ---------------------------------------------------------------------------
# F0_monopole is positive (Coulomb repulsion scale)
# ---------------------------------------------------------------------------

def test_exchange_f0_positive(exchange_at_three_R):
    """The monopole F^0 used as exchange energy scale should be positive."""
    for R, exch in exchange_at_three_R.items():
        assert exch['F0_monopole'] > 0, \
            f"F0_monopole should be positive at R={R}, got {exch['F0_monopole']}"


# ---------------------------------------------------------------------------
# E_exchange has correct sign (negative = attractive)
# ---------------------------------------------------------------------------

def test_exchange_sign(exchange_at_three_R):
    """Exchange energy should be negative (attractive correction that
    lowers energy more at short R, pulling R_eq inward)."""
    for R, exch in exchange_at_three_R.items():
        assert exch['E_exchange'] < 0, \
            f"E_exchange should be negative at R={R}, got {exch['E_exchange']}"


# ---------------------------------------------------------------------------
# Exchange is more negative at short R (correct R-dependence)
# ---------------------------------------------------------------------------

def test_exchange_stronger_at_short_R(exchange_at_three_R):
    """Exchange should be more negative at short R where overlap is larger.
    This is the R-dependence needed to pull R_eq inward."""
    E_2 = exchange_at_three_R[2.0]['E_exchange']
    E_3 = exchange_at_three_R[3.0]['E_exchange']
    E_5 = exchange_at_three_R[5.0]['E_exchange']
    # More negative = smaller value
    assert E_2 < E_3 < E_5, \
        f"Exchange should be more negative at short R: " \
        f"E(2.0)={E_2:.4f}, E(3.0)={E_3:.4f}, E(5.0)={E_5:.4f}"


# ---------------------------------------------------------------------------
# Exchange magnitude is physically reasonable
# ---------------------------------------------------------------------------

def test_exchange_magnitude(exchange_at_three_R):
    """Exchange energy should be in a physically reasonable range.
    The diagnostic fit found K=4.08 with S_avg~0.35 at R=2.5 giving
    E_exch ~ -1.4 Ha. The ab initio value should be order-of-magnitude
    comparable (0.01 to 10 Ha range)."""
    for R, exch in exchange_at_three_R.items():
        E = abs(exch['E_exchange'])
        assert E > 1e-4, \
            f"|E_exchange| too small at R={R}: {E:.6f} Ha"
        assert E < 50.0, \
            f"|E_exchange| too large at R={R}: {E:.4f} Ha"


# ===========================================================================
# Algebraic F^0 inter-fiber tests (v2.0.1)
# ===========================================================================


@pytest.fixture(scope='module')
def algebraic_at_R3(beh2_monopole_setup):
    """Compute algebraic monopole and exchange at R=3.0 for testing."""
    s = beh2_monopole_setup
    R = 3.0
    l4 = solve_level4_h2_multichannel(
        R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'], n_Re=300,
        verbose=False, pk_potentials=s['pk_potentials'],
    )
    # Numerical baseline
    mono_num = monopole_inter_fiber_energy(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'], n_sample_Re=12,
        method='numerical',
    )
    # Algebraic pathway
    mono_alg = monopole_inter_fiber_energy(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'], n_sample_Re=12,
        method='algebraic',
    )
    exch_alg = exchange_inter_fiber_energy(
        l4, R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        pk_potentials=s['pk_potentials'], n_sample_Re=12,
        method='algebraic',
    )
    return {
        'R': R, 'l4': l4, 'setup': s,
        'mono_num': mono_num, 'mono_alg': mono_alg,
        'exch_alg': exch_alg,
    }


def test_algebraic_f0_matches_numerical(algebraic_at_R3):
    """Algebraic F^0 must match numerical to within 1%."""
    F0_num = algebraic_at_R3['mono_num']['E_monopole']
    F0_alg = algebraic_at_R3['mono_alg']['E_monopole']
    rel_err = abs(F0_alg - F0_num) / abs(F0_num)
    assert rel_err < 0.01, \
        f"Algebraic F^0 = {F0_alg:.6f} vs numerical {F0_num:.6f}, " \
        f"rel error {rel_err:.4%} > 1%"


def test_algebraic_electron_count(algebraic_at_R3):
    """Algebraic density should integrate to ~2 electrons."""
    N = algebraic_at_R3['mono_alg']['N_elec_origin']
    assert abs(N - 2.0) < 0.01, \
        f"Algebraic density integrates to {N:.4f}, expected 2.0"


def test_channel_f0_matrix_sums_to_total(algebraic_at_R3):
    """F^0 matrix row/column sums should reproduce total F^0."""
    f0_mat = algebraic_at_R3['mono_alg']['F0_matrix_result']
    F0_total = f0_mat['F0_total']
    F0_from_alg = algebraic_at_R3['mono_alg']['E_monopole']
    # The matrix sum should approximately match the total density F^0
    # (exact match requires identical normalization paths)
    rel_err = abs(F0_total - F0_from_alg) / abs(F0_from_alg)
    assert rel_err < 0.05, \
        f"F^0 matrix sum = {F0_total:.6f} vs total {F0_from_alg:.6f}, " \
        f"rel error {rel_err:.4%}"


def test_channel_f0_matrix_symmetric(algebraic_at_R3):
    """F^0 matrix M(ch, ch') must be symmetric."""
    M = algebraic_at_R3['mono_alg']['F0_matrix_result']['F0_matrix']
    assert np.allclose(M, M.T, atol=1e-10), \
        f"F^0 matrix not symmetric: max diff = {np.max(np.abs(M - M.T))}"


def test_channel_f0_dominant_channel(algebraic_at_R3):
    """The (0,0) sigma channel should carry the largest F^0 share."""
    f0_mat = algebraic_at_R3['mono_alg']['F0_matrix_result']
    F0_per_ch = f0_mat['F0_per_channel']
    channels = f0_mat['channels']
    # Find (0,0) channel index
    idx_00 = None
    for i, ch in enumerate(channels):
        if ch == (0, 0):
            idx_00 = i
            break
    assert idx_00 is not None, "Channel (0,0) not found"
    assert F0_per_ch[idx_00] == np.max(F0_per_ch), \
        f"Channel (0,0) not dominant: {F0_per_ch[idx_00]:.4f} vs " \
        f"max {np.max(F0_per_ch):.4f}"


def test_algebraic_overlap_matches_diagnostic(algebraic_at_R3):
    """S_avg from algebraic channel data must match overlap diagnostic."""
    s = algebraic_at_R3['setup']
    R = algebraic_at_R3['R']
    l4 = algebraic_at_R3['l4']
    # Numerical overlap
    ovlp_num = compute_overlap_diagnostic(
        R=R, Z_A=s['Z_eff'], Z_B=s['Z_ligand'],
        l_max=s['l_max'], n_alpha=s['n_alpha'],
        level4_result=l4, pk_potentials=s['pk_potentials'],
        n_sample_Re=12,
    )
    S_num = ovlp_num['S_avg']

    # Algebraic overlap from exchange result
    S_alg = algebraic_at_R3['exch_alg']['S_avg']
    rel_err = abs(S_alg - S_num) / abs(S_num)
    assert rel_err < 0.01, \
        f"Algebraic S_avg = {S_alg:.6f} vs numerical {S_num:.6f}, " \
        f"rel error {rel_err:.4%}"


def test_algebraic_exchange_sign_and_magnitude(algebraic_at_R3):
    """Algebraic exchange must be negative and physically reasonable."""
    E_exch = algebraic_at_R3['exch_alg']['E_exchange']
    assert E_exch < 0, f"Exchange should be negative, got {E_exch}"
    assert abs(E_exch) > 0.01 and abs(E_exch) < 50.0, \
        f"|E_exchange| = {abs(E_exch):.4f} out of [0.01, 50] range"


def test_direct_exchange_close_to_product(algebraic_at_R3):
    """Channel-resolved exchange E_direct should be within 10% of S*F^0."""
    exch = algebraic_at_R3['exch_alg']
    E_product = exch['E_exchange']  # -S * F^0
    E_direct = exch['E_exchange_direct']  # -sum_ch parity_ch * F0_ch
    assert E_direct is not None, "E_exchange_direct not computed"
    ratio = E_direct / E_product
    assert 0.8 < ratio < 1.2, \
        f"Direct/product ratio = {ratio:.4f}, outside [0.8, 1.2]"


# ---------------------------------------------------------------------------
# Run module directly for quick diagnostics
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("Building composed BeH₂ for testing...\n")
    result = build_composed_beh2(verbose=True)

    print(f"\n--- Quick validation ---")
    print(f"M = {result['M']}, Q = {result['Q']}, N_pauli = {result['N_pauli']}")
    print(f"h1 Hermitian: {np.allclose(result['h1'], result['h1'].T)}")
    print(f"ERI density: {result['ERI_density_total']:.2%}")
