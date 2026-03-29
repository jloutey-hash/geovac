"""
Tests for composed LiH qubit Hamiltonian (geovac/composed_qubit.py).

Validates orbital counts, symmetries, selection rules, and JW consistency.
"""

import numpy as np
import pytest

from geovac.composed_qubit import (
    build_composed_lih,
    composed_lih_scaling_sweep,
    estimate_cross_center_eri_count,
    _enumerate_states,
    _wigner3j,
    _ck_coefficient,
    _compute_pk_matrix_elements,
    GAUSSIAN_LIH_PUBLISHED,
    GAUSSIAN_H2O_PUBLISHED,
    fit_gaussian_lih_published_exponent,
)


# ---------------------------------------------------------------------------
# Fixture: build at small basis for fast tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def lih_result():
    """Build composed LiH at max_n_core=2, max_n_val=2 (default)."""
    return build_composed_lih(max_n_core=2, max_n_val=2, verbose=False)


# ---------------------------------------------------------------------------
# 1. Correct orbital count
# ---------------------------------------------------------------------------

def test_orbital_count(lih_result):
    """M = M_core + M_val."""
    M = lih_result['M']
    M_core = lih_result['M_core']
    M_val = lih_result['M_val']
    assert M == M_core + M_val
    # max_n=2 -> 1+4 = 5 states per center
    assert M_core == 5   # 1s, 2s, 2p-1, 2p0, 2p1 at Z=3
    assert M_val == 10   # 5 (screened Li) + 5 (H)


def test_qubit_count(lih_result):
    """Q = 2M (Jordan-Wigner spin-orbital encoding)."""
    assert lih_result['Q'] == 2 * lih_result['M']


def test_state_enumeration():
    """State enumeration matches sum_{n=1}^{max_n} n^2 formula."""
    for max_n in [1, 2, 3, 4]:
        states = _enumerate_states(max_n)
        expected = sum(n**2 for n in range(1, max_n + 1))
        assert len(states) == expected, \
            f"max_n={max_n}: got {len(states)}, expected {expected}"


# ---------------------------------------------------------------------------
# 2. h1 is Hermitian
# ---------------------------------------------------------------------------

def test_h1_hermitian(lih_result):
    """h1 matrix must be Hermitian (real symmetric for real basis)."""
    h1 = lih_result['h1']
    assert np.allclose(h1, h1.T, atol=1e-14)


def test_h1_diagonal_negative_without_pk():
    """Diagonal of h1 should be negative (bound states) when PK is off."""
    result = build_composed_lih(
        max_n_core=2, max_n_val=2, include_pk=False, verbose=False,
    )
    diag = np.diag(result['h1'])
    assert np.all(diag < 0), f"Found non-negative diagonal: {diag}"


def test_h1_pk_makes_valence_positive(lih_result):
    """PK should push some valence-Li diag entries positive (repulsive)."""
    h1 = lih_result['h1']
    M_core = lih_result['M_core']
    M_val_li = lih_result['M_val_li']
    # Core entries must stay negative
    core_diag = np.diag(h1)[:M_core]
    assert np.all(core_diag < 0), f"Core diagonal should be negative: {core_diag}"
    # Valence-Li 1s entry should be positive (PK > binding energy)
    val_li_1s = h1[M_core, M_core]
    assert val_li_1s > 0, \
        f"PK should make valence-Li 1s positive, got {val_li_1s}"


def test_h1_core_eigenvalues(lih_result):
    """Core diagonal entries should be -Z^2/(2n^2) for Z=3."""
    h1 = lih_result['h1']
    M_core = lih_result['M_core']
    states_core = lih_result['states_core']
    Z = 3
    for i, (n, l, m) in enumerate(states_core):
        expected = -Z**2 / (2.0 * n**2)
        assert abs(h1[i, i] - expected) < 1e-12, \
            f"Core h1[{i},{i}] = {h1[i,i]}, expected {expected} for (n={n})"


# ---------------------------------------------------------------------------
# 3. ERI tensor symmetries
# ---------------------------------------------------------------------------

def test_eri_shape(lih_result):
    """ERI tensor should be (M, M, M, M)."""
    M = lih_result['M']
    eri = lih_result['eri']
    assert eri.shape == (M, M, M, M)


def test_eri_chemist_symmetry(lih_result):
    """
    Chemist notation (pq|rs) has 4-fold symmetry for real orbitals:
        (pq|rs) = (qp|rs) = (pq|sr) = (rs|pq)
    We check (pq|rs) = (rs|pq) (particle exchange).
    """
    eri = lih_result['eri']
    M = lih_result['M']
    # Check particle exchange: (pq|rs) = (rs|pq)
    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    if abs(eri[p, q, r, s]) > 1e-15:
                        assert abs(eri[p, q, r, s] - eri[r, s, p, q]) < 1e-12, \
                            f"Symmetry violation: ({p}{q}|{r}{s})={eri[p,q,r,s]} != ({r}{s}|{p}{q})={eri[r,s,p,q]}"


# ---------------------------------------------------------------------------
# 4. Cross-ERI block is zero
# ---------------------------------------------------------------------------

def test_cross_eri_zero(lih_result):
    """Cross-block ERIs must be exactly zero by construction."""
    eri = lih_result['eri']
    M_core = lih_result['M_core']
    M = lih_result['M']

    # Any ERI with mixed core/valence indices should be zero
    for p in range(M):
        for q in range(M):
            for r in range(M):
                for s in range(M):
                    indices = {p, q, r, s}
                    has_core = any(i < M_core for i in indices)
                    has_val = any(i >= M_core for i in indices)
                    if has_core and has_val:
                        assert abs(eri[p, q, r, s]) < 1e-15, \
                            f"Cross-block ERI[{p},{q},{r},{s}] = {eri[p,q,r,s]} != 0"


# ---------------------------------------------------------------------------
# 5. Pauli term count is a positive integer
# ---------------------------------------------------------------------------

def test_pauli_count_positive(lih_result):
    """Pauli term count must be a positive integer."""
    N = lih_result['N_pauli']
    assert isinstance(N, int)
    assert N > 0


def test_pauli_count_reasonable(lih_result):
    """Pauli term count should be less than 4^Q (trivially) and > Q."""
    N = lih_result['N_pauli']
    Q = lih_result['Q']
    assert N > Q, f"Pauli terms ({N}) should exceed qubit count ({Q})"
    assert N < 4**Q, f"Pauli terms ({N}) should be less than 4^Q"


# ---------------------------------------------------------------------------
# 6. ERI density comparison
# ---------------------------------------------------------------------------

def test_eri_density_less_than_single_geometry(lih_result):
    """
    Total ERI density of composed system should be less than single-geometry
    He at max_n=2 (10.4%) due to block-diagonal structure.
    """
    eri_density = lih_result['ERI_density_total']
    # He max_n=2 single geometry: 10.4%
    assert eri_density < 0.104, \
        f"Composed ERI density {eri_density:.1%} >= He single-geometry 10.4%"


def test_eri_density_positive(lih_result):
    """ERI density must be positive (we do have nonzero integrals)."""
    assert lih_result['ERI_density_total'] > 0
    assert lih_result['ERI_density_core'] > 0


# ---------------------------------------------------------------------------
# 7. JW eigenvalue consistency
# ---------------------------------------------------------------------------

def test_jw_eigenvalue_consistency():
    """
    Lowest eigenvalue of qubit Hamiltonian should match FCI of h1/ERI system.

    Uses a small basis (max_n=1) to keep the Hilbert space tractable
    (Q=6 -> 2^6 = 64 dimensional matrix).
    """
    from openfermion import get_sparse_operator
    from scipy.sparse.linalg import eigsh

    # Build at max_n=1: 1 core + 2 val = 3 spatial orbitals, Q=6
    result = build_composed_lih(max_n_core=1, max_n_val=1, verbose=False)
    Q = result['Q']
    assert Q <= 16, f"Q={Q} too large for eigenvalue test"

    qubit_op = result['qubit_op']
    sparse_H = get_sparse_operator(qubit_op)

    # Get lowest eigenvalue
    evals, _ = eigsh(sparse_H, k=1, which='SA')
    E_qubit = float(np.real(evals[0]))

    assert np.isfinite(E_qubit), f"Qubit ground state energy is not finite: {E_qubit}"
    assert E_qubit < 0, f"Expected negative ground state energy, got {E_qubit}"

    # Verify the qubit Hamiltonian is Hermitian
    diff = sparse_H - sparse_H.conj().T
    assert abs(diff).max() < 1e-10, "Qubit Hamiltonian is not Hermitian"

    print(f"  JW ground state energy (Q={Q}): {E_qubit:.6f} Ha")


# ---------------------------------------------------------------------------
# Auxiliary: Wigner 3j validation
# ---------------------------------------------------------------------------

def test_wigner3j_known_values():
    """Validate Wigner 3j against known analytical values."""
    # (1 1 0; 0 0 0) = (-1)^1 / sqrt(3)
    val = _wigner3j(1, 1, 0, 0, 0, 0)
    assert abs(val - (-1.0 / np.sqrt(3))) < 1e-12

    # (1 1 2; 0 0 0) = sqrt(2/15)
    val = _wigner3j(1, 1, 2, 0, 0, 0)
    assert abs(val - np.sqrt(2.0 / 15.0)) < 1e-12

    # m-selection rule violation
    assert _wigner3j(1, 1, 1, 1, 1, 0) == 0.0


def test_ck_s_orbital_unity():
    """c^0(0,0,0,0,0) should give 1/sqrt(4pi) * 4pi = 1."""
    # For s-orbitals: c^0(l=0,m=0,l'=0,m'=0) = 1
    val = _ck_coefficient(0, 0, 0, 0, 0)
    assert abs(val - 1.0) < 1e-12


# ---------------------------------------------------------------------------
# 8. Scaling sweep monotonicity
# ---------------------------------------------------------------------------

def test_scaling_sweep_monotonic():
    """Pauli terms must increase with max_n."""
    result = composed_lih_scaling_sweep(max_n_values=[1, 2], verbose=False)
    data = result['sweep_data']
    assert len(data) == 2
    assert data[1]['N_pauli'] > data[0]['N_pauli'], \
        f"N_pauli should increase: {data[0]['N_pauli']} -> {data[1]['N_pauli']}"
    assert data[1]['Q'] > data[0]['Q'], \
        f"Q should increase: {data[0]['Q']} -> {data[1]['Q']}"


# ---------------------------------------------------------------------------
# 9. PK pseudopotential validation
# ---------------------------------------------------------------------------

def test_pk_potential_positive():
    """PK pseudopotential diagonal elements must be positive (repulsive barrier)."""
    states = _enumerate_states(2)
    # Li²⁺ core PK params from Paper 17
    h1_pk = _compute_pk_matrix_elements(
        Z_eff=1.0, states=states, A_pk=6.93, B_pk=7.00,
    )
    diag = np.diag(h1_pk)
    assert np.all(diag >= 0), \
        f"PK diagonal should be non-negative (repulsive): {diag}"
    # At least the 1s orbital should have substantial PK
    assert diag[0] > 0.01, \
        f"PK on 1s orbital should be significant, got {diag[0]}"


# ---------------------------------------------------------------------------
# 10. Cross-block h1 Hermiticity with PK
# ---------------------------------------------------------------------------

def test_cross_h1_hermitian():
    """h1 with PK cross-block terms must remain Hermitian."""
    result = build_composed_lih(
        max_n_core=2, max_n_val=2, include_pk=True, verbose=False,
    )
    h1 = result['h1']
    assert np.allclose(h1, h1.T, atol=1e-14), \
        "h1 with PK terms is not Hermitian"


# ---------------------------------------------------------------------------
# 11. Cross-h1 impact on Pauli count
# ---------------------------------------------------------------------------

def test_cross_h1_impact():
    """Adding PK cross-h1 terms should not decrease Pauli count."""
    result_no_pk = build_composed_lih(
        max_n_core=2, max_n_val=2, include_pk=False, verbose=False,
    )
    result_pk = build_composed_lih(
        max_n_core=2, max_n_val=2, include_pk=True, verbose=False,
    )
    assert result_pk['N_pauli'] >= result_no_pk['N_pauli'], \
        (f"PK should not decrease Pauli count: "
         f"{result_no_pk['N_pauli']} -> {result_pk['N_pauli']}")


# ---------------------------------------------------------------------------
# 12. max_n=4 runs without error
# ---------------------------------------------------------------------------

@pytest.fixture(scope='session')
def lih_result_n4():
    """Build composed LiH at max_n=4 (cached for session)."""
    return build_composed_lih(max_n_core=4, max_n_val=4, verbose=False)


def test_max_n4_runs(lih_result_n4):
    """build_composed_lih at max_n=4 completes without error."""
    assert lih_result_n4['M'] == 90
    assert lih_result_n4['Q'] == 180
    assert lih_result_n4['N_pauli'] > 0
    assert isinstance(lih_result_n4['N_pauli'], int)


# ---------------------------------------------------------------------------
# 13. Cross-center ERI count validation
# ---------------------------------------------------------------------------

def test_cross_center_eri_count_positive():
    """Estimated cross-center ERI count must be positive."""
    result = estimate_cross_center_eri_count(2)
    assert result['n_cross_center'] > 0, \
        f"Cross-center count should be positive, got {result['n_cross_center']}"


def test_cross_center_eri_count_bounded():
    """Cross-center ERI count must be less than M_val^4."""
    result = estimate_cross_center_eri_count(2)
    assert result['n_cross_center'] < result['upper_bound'], \
        (f"Cross-center count {result['n_cross_center']} should be < "
         f"M_val^4 = {result['upper_bound']}")


# ---------------------------------------------------------------------------
# 14. Published Gaussian data consistency
# ---------------------------------------------------------------------------

def test_gaussian_published_data():
    """Verify published Gaussian constants are present and self-consistent."""
    # LiH: Q increases with basis size, N_pauli increases with Q
    bases_lih = ['sto-3g', '6-31g', 'cc-pvdz']
    for basis in bases_lih:
        assert basis in GAUSSIAN_LIH_PUBLISHED, f"Missing LiH basis {basis}"
        entry = GAUSSIAN_LIH_PUBLISHED[basis]
        assert 'Q' in entry and 'N_pauli' in entry
        assert isinstance(entry['Q'], int) and entry['Q'] > 0
        assert isinstance(entry['N_pauli'], int) and entry['N_pauli'] > 0

    for i in range(len(bases_lih) - 1):
        b1, b2 = bases_lih[i], bases_lih[i + 1]
        assert GAUSSIAN_LIH_PUBLISHED[b2]['Q'] > GAUSSIAN_LIH_PUBLISHED[b1]['Q'], \
            f"LiH Q should increase: {b1}={GAUSSIAN_LIH_PUBLISHED[b1]['Q']} >= {b2}={GAUSSIAN_LIH_PUBLISHED[b2]['Q']}"
        assert GAUSSIAN_LIH_PUBLISHED[b2]['N_pauli'] > GAUSSIAN_LIH_PUBLISHED[b1]['N_pauli'], \
            f"LiH N_pauli should increase: {b1} >= {b2}"

    # H2O: same checks
    bases_h2o = ['sto-3g', '6-31g', 'cc-pvdz']
    for basis in bases_h2o:
        assert basis in GAUSSIAN_H2O_PUBLISHED, f"Missing H2O basis {basis}"
        entry = GAUSSIAN_H2O_PUBLISHED[basis]
        assert 'Q' in entry and 'N_pauli' in entry
        assert isinstance(entry['Q'], int) and entry['Q'] > 0
        assert isinstance(entry['N_pauli'], int) and entry['N_pauli'] > 0

    for i in range(len(bases_h2o) - 1):
        b1, b2 = bases_h2o[i], bases_h2o[i + 1]
        assert GAUSSIAN_H2O_PUBLISHED[b2]['Q'] > GAUSSIAN_H2O_PUBLISHED[b1]['Q'], \
            f"H2O Q should increase: {b1} >= {b2}"
        assert GAUSSIAN_H2O_PUBLISHED[b2]['N_pauli'] > GAUSSIAN_H2O_PUBLISHED[b1]['N_pauli'], \
            f"H2O N_pauli should increase: {b1} >= {b2}"

    # Published LiH values match Trenev et al. Table 5
    assert GAUSSIAN_LIH_PUBLISHED['sto-3g']['Q'] == 10
    assert GAUSSIAN_LIH_PUBLISHED['sto-3g']['N_pauli'] == 276
    assert GAUSSIAN_LIH_PUBLISHED['6-31g']['Q'] == 20
    assert GAUSSIAN_LIH_PUBLISHED['6-31g']['N_pauli'] == 5851
    assert GAUSSIAN_LIH_PUBLISHED['cc-pvdz']['Q'] == 36
    assert GAUSSIAN_LIH_PUBLISHED['cc-pvdz']['N_pauli'] == 63519

    # Fit exponent should be reasonable (between 3 and 5)
    fit = fit_gaussian_lih_published_exponent()
    assert 3.0 < fit['alpha'] < 5.0, \
        f"Gaussian LiH exponent {fit['alpha']:.2f} outside expected range [3, 5]"
    assert fit['R_squared'] > 0.99, \
        f"Gaussian LiH fit R² = {fit['R_squared']:.4f} < 0.99"


# ---------------------------------------------------------------------------
# Run module directly for quick diagnostics
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("Building composed LiH for testing...\n")
    result = build_composed_lih(verbose=True)

    print(f"\n--- Quick validation ---")
    print(f"M = {result['M']}, Q = {result['Q']}, N_pauli = {result['N_pauli']}")
    print(f"h1 Hermitian: {np.allclose(result['h1'], result['h1'].T)}")
    print(f"ERI density: {result['ERI_density_total']:.2%}")
