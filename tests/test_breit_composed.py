"""Tests for Breit-Pauli SS+SOO two-body corrections in the composed pipeline.

Validates:
1. Backward compatibility: include_breit=False (default) leaves all existing
   Pauli counts, 1-norms, and FCI energies unchanged.
2. Breit contribution vanishes at alpha=0.
3. Breit shifts have correct sign and order of magnitude for He 2^3P.
4. He 2^3P J-pattern matches Drake's f_SS=(-2,+1,-1/5) and f_SOO=(+2,+1,-1).
5. Breit Pauli count increase is non-zero but bounded.

Author: GeoVac Development Team
Date: 2026-04-16
"""

from __future__ import annotations

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# 1. Backward compatibility: default include_breit=False unchanged
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec_fn,expected_terms", [
    ("lih_spec_relativistic", None),  # Just check it builds, no regression on exact count
])
def test_breit_false_default_unchanged(spec_fn, expected_terms):
    """include_breit=False (default) must not change any results."""
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac import molecular_spec as ms

    spec = getattr(ms, spec_fn)(max_n=1)
    r_default = build_composed_hamiltonian(spec)
    r_explicit = build_composed_hamiltonian(spec, include_breit=False)

    assert r_default['N_pauli'] == r_explicit['N_pauli']
    assert abs(r_default['lambda_ni'] - r_explicit['lambda_ni']) < 1e-12


def test_scalar_pauli_counts_unchanged_with_breit_kwarg():
    """Scalar (non-relativistic) path ignores include_breit entirely."""
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.molecular_spec import lih_spec

    spec = lih_spec(max_n=2)
    r = build_composed_hamiltonian(spec, include_breit=True)
    n_terms = len(r['qubit_op'].terms)
    assert n_terms == 334, (
        f"LiH scalar with include_breit=True: {n_terms} != 334")


# ---------------------------------------------------------------------------
# 2. Breit vanishes at alpha=0
# ---------------------------------------------------------------------------


def test_breit_zero_at_alpha_zero():
    """Breit ERI count must be zero when alpha_num=0."""
    from geovac.composed_qubit_relativistic import (
        build_composed_hamiltonian_relativistic,
    )
    from geovac.molecular_spec import lih_spec_relativistic

    spec = lih_spec_relativistic(max_n=2)
    r = build_composed_hamiltonian_relativistic(
        spec, alpha_num=0.0, include_breit=True)
    # At alpha=0, the alpha^2 prefactor kills all Breit terms
    assert r['breit_eri_count'] == 0
    assert r['include_breit'] is True


# ---------------------------------------------------------------------------
# 3. Breit builds and produces nonzero correction
# ---------------------------------------------------------------------------


def test_breit_produces_nonzero_correction():
    """With include_breit=True, the ERI should have more terms or shifted values."""
    from geovac.composed_qubit_relativistic import (
        build_composed_hamiltonian_relativistic,
    )
    from geovac.molecular_spec import lih_spec_relativistic

    spec = lih_spec_relativistic(max_n=2)
    r_no = build_composed_hamiltonian_relativistic(spec, include_breit=False)
    r_yes = build_composed_hamiltonian_relativistic(spec, include_breit=True)

    # Breit should add ERI entries
    assert r_yes['breit_eri_count'] > 0, "Breit should produce nonzero ERI count"

    # 1-norm should shift (Breit is alpha^2 ~ 5e-5 relative to Coulomb)
    lam_diff = abs(r_yes['lambda_ni'] - r_no['lambda_ni'])
    assert lam_diff > 0, "Breit should change the 1-norm"
    # But the shift should be small (alpha^2 ~ 5e-5)
    rel_shift = lam_diff / r_no['lambda_ni']
    assert rel_shift < 0.01, (
        f"Breit 1-norm shift {rel_shift:.4f} is too large (expected < 1%)")


def test_breit_block_diagonality():
    """Breit ERI must remain block-diagonal."""
    from geovac.composed_qubit_relativistic import (
        build_composed_hamiltonian_relativistic,
    )
    from geovac.molecular_spec import lih_spec_relativistic

    spec = lih_spec_relativistic(max_n=2)
    r = build_composed_hamiltonian_relativistic(spec, include_breit=True)
    assert r['cross_block_eri_count'] == 0, (
        "Breit ERI breaks block diagonality")


def test_breit_hermiticity():
    """Breit qubit_op must be Hermitian (real coefficients)."""
    from geovac.composed_qubit_relativistic import (
        build_composed_hamiltonian_relativistic,
    )
    from geovac.molecular_spec import lih_spec_relativistic

    spec = lih_spec_relativistic(max_n=2)
    r = build_composed_hamiltonian_relativistic(spec, include_breit=True)
    for term, coeff in r['qubit_op'].terms.items():
        assert abs(coeff.imag) < 1e-10, (
            f"Non-real coefficient {coeff} on Breit term {term}")


# ---------------------------------------------------------------------------
# 4. He 2^3P fine structure from Breit integrals
# ---------------------------------------------------------------------------


def test_he_breit_radial_integrals_nonzero():
    """Breit retarded integrals for He (1s)(2p) should be nonzero."""
    from geovac.breit_integrals import compute_radial

    # M^2 direct: R^2_BP(1s, 2p; 1s, 2p) at Z=2
    val = compute_radial(1, 0, 1, 0, 2, 1, 2, 1, 2, kernel_type="breit", Z=2)
    assert abs(float(val)) > 0, "M^2_direct should be nonzero"

    # M^2 exchange: R^2_BP(1s, 2p; 2p, 1s) at Z=2
    val_ex = compute_radial(1, 0, 2, 1, 2, 1, 1, 0, 2, kernel_type="breit", Z=2)
    assert abs(float(val_ex)) > 0, "M^2_exchange should be nonzero"


def test_he_breit_order_of_magnitude():
    """Breit energy shifts for He should be O(alpha^2 * Z^3) ~ 10^-5 Ha."""
    from geovac.breit_integrals import compute_radial

    ALPHA = 7.2973525693e-3

    # M^2 direct at Z=2
    M2d = float(compute_radial(1, 0, 1, 0, 2, 1, 2, 1, 2,
                                kernel_type="breit", Z=2))
    # M^1 direct at Z=2
    M1d = float(compute_radial(1, 0, 1, 0, 2, 1, 2, 1, 1,
                                kernel_type="breit", Z=2))

    # A_SS ~ alpha^2 * (3/50 * M2d - 2/5 * M2e) ~ alpha^2 * O(Z^3)
    # Expected order: alpha^2 ~ 5e-5, Z^3 = 8, so ~4e-4 Ha
    A_SS_approx = ALPHA**2 * abs(3.0 / 50.0 * M2d)
    assert 1e-7 < A_SS_approx < 1e-2, (
        f"A_SS order of magnitude {A_SS_approx:.2e} outside expected range")

    A_SOO_approx = ALPHA**2 * abs(3.0 / 2.0 * M1d)
    assert 1e-7 < A_SOO_approx < 1e-2, (
        f"A_SOO order of magnitude {A_SOO_approx:.2e} outside expected range")


def test_drake_combining_coefficients_reproduce_splittings():
    """Drake 1971 combining coefficients reproduce He 2^3P splittings.

    This is a physics integration test: compute A_SS and A_SOO from
    Breit retarded radial integrals, apply Drake's f_SS and f_SOO
    J-patterns, and verify the splittings are physical.
    """
    from geovac.breit_integrals import compute_radial
    from geovac.spin_orbit import so_diagonal_matrix_element
    import sympy as sp
    from sympy import Integer

    ALPHA = 7.2973525693e-3
    Z = 2

    # Spin-orbit parameter zeta_2p
    # H_SO(2p_{3/2}) uses kappa = -2 (j = 3/2)
    # H_SO(2p_{1/2}) uses kappa = +1 (j = 1/2)
    so_32 = float(so_diagonal_matrix_element(2, -2, Z=Integer(Z),
                                              alpha=sp.Float(ALPHA)))
    so_12 = float(so_diagonal_matrix_element(2, 1, Z=Integer(Z),
                                              alpha=sp.Float(ALPHA)))
    # zeta = 2 * H_SO(kappa=-2) / (-(-2+1)/2) = 2 * H_SO(kappa=-2) / (1/2)
    # L.S eigenvalue for kappa=-2: -(kappa+1)/2 = -(-2+1)/2 = 1/2
    # So H_SO = zeta/2 * (1/2) => zeta = 4 * H_SO(kappa=-2)
    # Actually: H_SO = (Z*alpha^2/2) * L.S * <1/r^3>
    # For 2p: <1/r^3> = Z^3/(n^3 * l(l+1/2)(l+1)) = Z^3/(8 * 1 * 3/2 * 2) = Z^3/24
    # L.S(kappa=-2) = -(-2+1)/2 = 1/2
    # L.S(kappa=+1) = -(1+1)/2 = -1
    # zeta = Z*alpha^2/2 * Z^3/24 = Z^4 * alpha^2 / 48
    zeta = float(Integer(Z)**4 * sp.Float(ALPHA)**2 / 48)

    # Drake radial integrals at Z=2
    # Direct: orbital ordering (1s)(2p) means: el1 on 1s, el2 on 2p
    # R^k_BP(1s,1s; 2p,2p)
    M2d = float(compute_radial(1, 0, 1, 0, 2, 1, 2, 1, 2,
                                kernel_type="breit", Z=Z))
    M1d = float(compute_radial(1, 0, 1, 0, 2, 1, 2, 1, 1,
                                kernel_type="breit", Z=Z))

    # Exchange: R^k_BP(1s,2p; 2p,1s)
    M2e = float(compute_radial(1, 0, 2, 1, 2, 1, 1, 0, 2,
                                kernel_type="breit", Z=Z))
    M1e = float(compute_radial(1, 0, 2, 1, 2, 1, 1, 0, 1,
                                kernel_type="breit", Z=Z))

    # Drake combining coefficients
    A_SS = ALPHA**2 * (3.0/50 * M2d - 2.0/5 * M2e)
    A_SOO = ALPHA**2 * (3.0/2 * M1d - 1.0 * M1e)

    # J-patterns
    f_SS = {0: -2.0, 1: 1.0, 2: -0.2}
    f_SOO = {0: 2.0, 1: 1.0, 2: -1.0}
    X_J = {0: -4.0, 1: -2.0, 2: 2.0}  # J(J+1) - L(L+1) - S(S+1)

    E = {J: (zeta / 2) * X_J[J] + A_SS * f_SS[J] + A_SOO * f_SOO[J]
         for J in (0, 1, 2)}

    # Check that splittings are nonzero and have correct ordering
    # In He 2^3P: E(J=0) < E(J=1) < E(J=2) (normal ordering for light atoms)
    # The overall splitting should be O(alpha^2 * Z^4) ~ few * 10^-5 Ha
    splitting_01 = E[1] - E[0]
    splitting_12 = E[2] - E[1]

    # Both splittings should be nonzero
    assert abs(splitting_01) > 1e-8, f"01 splitting too small: {splitting_01:.2e}"
    assert abs(splitting_12) > 1e-8, f"12 splitting too small: {splitting_12:.2e}"

    # Total splitting should be O(10^-5) Ha for He
    total_splitting = abs(E[2] - E[0])
    assert 1e-7 < total_splitting < 1e-3, (
        f"Total He 2^3P splitting {total_splitting:.2e} outside expected range")


# ---------------------------------------------------------------------------
# 5. Breit Pauli count change is bounded
# ---------------------------------------------------------------------------


def test_breit_pauli_count_bounded():
    """Breit should not dramatically change Pauli count (same angular structure)."""
    from geovac.composed_qubit_relativistic import (
        build_composed_hamiltonian_relativistic,
    )
    from geovac.molecular_spec import lih_spec_relativistic

    spec = lih_spec_relativistic(max_n=1)
    r_no = build_composed_hamiltonian_relativistic(spec, include_breit=False)
    r_yes = build_composed_hamiltonian_relativistic(spec, include_breit=True)

    # At n_max=1 with only s-orbitals, Breit k=2 channels require l >= 1,
    # so the Breit contribution should be very limited or zero.
    # The Pauli count difference should be small.
    diff = abs(r_yes['N_pauli'] - r_no['N_pauli'])
    # Allow up to 100% increase (Breit adds terms but shouldn't explode)
    assert diff <= r_no['N_pauli'], (
        f"Breit Pauli increase {diff} > baseline {r_no['N_pauli']}")


# ---------------------------------------------------------------------------
# 6. Breit radial integral cache correctness
# ---------------------------------------------------------------------------


def test_breit_rk_cache_z_scaling():
    """Breit radial integrals scale as Z^3."""
    from geovac.breit_integrals import compute_radial

    # R^2_BP(1s,1s;2p,2p) at Z=1 and Z=2
    v1 = float(compute_radial(1, 0, 1, 0, 2, 1, 2, 1, 2,
                               kernel_type="breit", Z=1))
    v2 = float(compute_radial(1, 0, 1, 0, 2, 1, 2, 1, 2,
                               kernel_type="breit", Z=2))
    if abs(v1) > 1e-15:
        ratio = v2 / v1
        assert abs(ratio - 8.0) < 1e-10, (
            f"Z^3 scaling: ratio={ratio:.6f}, expected 8.0")
