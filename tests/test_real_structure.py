"""Tests for geovac/real_structure.py.

Sprint WH1-Connes Step 2: real structure J on truncated operator system.

Test panel
==========

(A) Basic operator-level checks:
    - J as antilinear operator: J(alpha psi) = conj(alpha) J(psi).
    - J^2 = -I (KO-dim 3 expected sign).
    - U is unitary (J as antilinear isometry).

(B) J-D relation:
    - J D = +D J (KO-dim 3 expected sign) on the WEYL sector with the
      truthful CH Dirac.
    - J D = +D J on the FULL DIRAC sector with the truthful CH Dirac.
    - J D = ? on the FULL DIRAC sector with the off-diag CH Dirac
      (this can FAIL because the off-diag perturbation is not J-symmetric;
      we test and report).

(C) J preserves the operator system:
    - Each generator a in O has J a J^{-1} also in O.

(D) Connes' order conditions on operator system O (a substitute for
    the algebra A in the original axioms):
    - Order zero: [a, J b J^{-1}] = 0 for all a, b in O.
    - Order one: [[D, a], J b J^{-1}] = 0 for all a, b in O.

(E) Negative controls:
    - A wrong-sign convention (J^2 = +I) should fail expectations.
    - A J that flips chirality on full_dirac SHOULD give J D = -D J,
      flagging the convention issue.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.real_structure import (
    RealStructure,
    audit_J,
    build_J_full_dirac,
    build_J_weyl,
    verify_J_D_relation,
    verify_J_preserves_O,
    verify_J_squared,
    verify_J_unitary,
    verify_order_one,
    verify_order_zero,
    _full_dirac_J_target,
    _spinor_J_phase,
    _spinor_J_target,
)
from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    camporesi_higuchi_offdiag_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.spinor_operator_system import (
    SpinorLabel,
    SpinorTruncatedOperatorSystem,
    camporesi_higuchi_dirac_matrix,
    spinor_basis,
    spinor_dim,
)


# ---------------------------------------------------------------------------
# Phase formula checks
# ---------------------------------------------------------------------------


def test_spinor_J_phase_returns_unit_complex():
    """sigma(label) is a unit complex number."""
    for n in range(1, 4):
        for l in range(n):
            for two_m_j in range(-(2 * l + 1), 2 * l + 2, 2):
                label = SpinorLabel(n_fock=n, l=l, two_m_j=two_m_j)
                phase = _spinor_J_phase(label)
                assert abs(abs(phase) - 1.0) < 1e-14


def test_spinor_J_phase_squares_to_minus_one_on_pair():
    """sigma(label) * conj(sigma(J_label)) = -1, as required for J^2 = -I."""
    for n in range(1, 4):
        for l in range(n):
            for two_m_j in range(-(2 * l + 1), 2 * l + 2, 2):
                label = SpinorLabel(n_fock=n, l=l, two_m_j=two_m_j)
                target = _spinor_J_target(label)
                phase = _spinor_J_phase(label)
                target_phase = _spinor_J_phase(target)
                # We need phase * conj(target_phase) = -1
                # (this is what makes J^2 = U conj(U) at index = (target, label)
                # produce -delta on diagonal)
                product = phase * np.conj(target_phase)
                assert abs(product - (-1.0)) < 1e-14, (
                    f"phase product = {product}, expected -1, "
                    f"label = {label}, target = {target}"
                )


# ---------------------------------------------------------------------------
# Permutation involutivity
# ---------------------------------------------------------------------------


def test_J_target_involutive_weyl():
    """_spinor_J_target is an involution: J(J(label)) = label."""
    for n in range(1, 4):
        for l in range(n):
            for two_m_j in range(-(2 * l + 1), 2 * l + 2, 2):
                label = SpinorLabel(n_fock=n, l=l, two_m_j=two_m_j)
                target = _spinor_J_target(label)
                target_target = _spinor_J_target(target)
                assert label == target_target, (
                    f"label {label} -> J -> {target} -> J -> {target_target}, "
                    f"expected {label}"
                )


def test_J_target_involutive_full_dirac():
    """_full_dirac_J_target is an involution."""
    for n in range(1, 4):
        for l in range(n):
            for two_m_j in range(-(2 * l + 1), 2 * l + 2, 2):
                for chi in (+1, -1):
                    label = FullDiracLabel(
                        n_fock=n, l=l, two_m_j=two_m_j, chirality=chi
                    )
                    target = _full_dirac_J_target(label)
                    target_target = _full_dirac_J_target(target)
                    assert label == target_target


# ---------------------------------------------------------------------------
# RealStructure unitarity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_J_unitary_weyl(n_max):
    """U is unitary on the Weyl sector."""
    J = build_J_weyl(n_max)
    ok, err = verify_J_unitary(J)
    assert ok, f"U @ U^dagger != I, max err {err:.3e}"


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_J_unitary_full_dirac(n_max):
    """U is unitary on the full Dirac sector."""
    J = build_J_full_dirac(n_max)
    ok, err = verify_J_unitary(J)
    assert ok, f"U @ U^dagger != I, max err {err:.3e}"


# ---------------------------------------------------------------------------
# J^2 = -I (the KO-dim 3 sign)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_J_squared_minus_I_weyl(n_max):
    """J^2 = -I on the Weyl sector at every n_max."""
    J = build_J_weyl(n_max)
    ok, err = verify_J_squared(J, expected_sign=-1)
    assert ok, f"J^2 != -I, max err {err:.3e}"


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_J_squared_minus_I_full_dirac(n_max):
    """J^2 = -I on the full Dirac sector at every n_max."""
    J = build_J_full_dirac(n_max)
    ok, err = verify_J_squared(J, expected_sign=-1)
    assert ok, f"J^2 != -I, max err {err:.3e}"


# ---------------------------------------------------------------------------
# Antilinearity sanity
# ---------------------------------------------------------------------------


def test_J_antilinear_weyl():
    """J(alpha psi) = conj(alpha) J(psi)."""
    n_max = 2
    J = build_J_weyl(n_max)
    psi = np.random.randn(J.dim) + 1j * np.random.randn(J.dim)
    alpha = 1.7 - 2.3j
    LHS = J.apply(alpha * psi)
    RHS = np.conj(alpha) * J.apply(psi)
    assert np.max(np.abs(LHS - RHS)) < 1e-12


def test_J_antilinear_full_dirac():
    """J(alpha psi) = conj(alpha) J(psi) on full Dirac."""
    n_max = 2
    J = build_J_full_dirac(n_max)
    psi = np.random.randn(J.dim) + 1j * np.random.randn(J.dim)
    alpha = -3.2 + 1.1j
    LHS = J.apply(alpha * psi)
    RHS = np.conj(alpha) * J.apply(psi)
    assert np.max(np.abs(LHS - RHS)) < 1e-12


# ---------------------------------------------------------------------------
# J-D relation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_J_commutes_D_weyl_truthful(n_max):
    """J D = +D J on Weyl sector with truthful CH Dirac (D|n,l,m_j> =
    (n + 1/2)|n,l,m_j>, diagonal in n_fock).

    Since truthful CH is diagonal and J is m_j-flipping (which is
    independent of n_fock), this should hold trivially with sign +1.
    """
    basis = spinor_basis(n_max)
    D = camporesi_higuchi_dirac_matrix(basis)
    J = build_J_weyl(n_max)
    ok, err = verify_J_D_relation(J, D, expected_sign=+1)
    assert ok, f"J D != +D J, max err {err:.3e}"


@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_J_commutes_D_full_dirac_truthful(n_max):
    """J D = +D J on full Dirac sector with truthful CH Dirac.

    Truthful: D|n, l, m_j, chi> = chi*(n + 1/2)|n, l, m_j, chi>.
    J flips m_j (block-diagonal in chi). Since D depends on (n, chi)
    only, and J preserves both, J D = D J trivially.
    """
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    J = build_J_full_dirac(n_max)
    ok, err = verify_J_D_relation(J, D, expected_sign=+1)
    assert ok, f"J D != +D J, max err {err:.3e}"


# Note: J D = +D J on the OFFDIAG full Dirac is NOT automatically true.
# We test it and report (HOLDS or APPROX or FAILS).
@pytest.mark.parametrize("n_max", [1, 2, 3])
def test_J_commutes_D_full_dirac_offdiag_status(n_max):
    """J D vs +D J on full Dirac sector with OFFDIAG CH perturbation.

    The offdiag Dirac adds Hermitian E1-style hopping H[i, j] = alpha for
    |Δn|=1, |Δl|=1, |Δm_j|<=1 (within and across chirality). This is
    NOT manifestly J-symmetric unless the perturbation is symmetric under
    m_j -> -m_j on bra/ket simultaneously.

    The off-diag term H[i, j] = alpha is constant along (i, j) pairs that
    satisfy the rule. The rule depends on |Δm_j|, which IS symmetric
    under simultaneous m_j -> -m_j. So we *expect* the offdiag piece to
    commute with J too.

    But the diagonal-lifters l_lift * l + m_lift * 2 m_j include a
    LINEAR term in m_j, which is ANTI-symmetric under m_j -> -m_j. Under
    J, the diagonal entry at index i becomes the entry at index pi(i)
    where pi sends m_j -> -m_j; the m_j-lifter sign flips.

    Therefore, with default lifters (n_lift=1, l_lift=0.1, m_lift=0.005),
    J D != +D J on the offdiag. We test and report rather than assert
    pass/fail.

    EXCEPT: U conj(D) vs +D U:
      The condition is U conj(D) = +D U.
      D real-valued (with default lifters) -> conj(D) = D.
      So we need U D = D U, i.e. D is invariant under U-conjugation.
      The m_j-linear term in diag(D) becomes its OPPOSITE sign at the
      pi-image, so U D != D U in general.
    """
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_offdiag_dirac_matrix(basis)
    J = build_J_full_dirac(n_max)
    ok, err = verify_J_D_relation(J, D, expected_sign=+1, tol=1e-3)
    # We don't ASSERT the result -- we record it and let the test
    # display the actual residual via pytest output.
    print(f"\n  [n_max={n_max}, full_dirac, offdiag]"
          f" J D = +D J: passes={ok}, max err = {err:.3e}")


@pytest.mark.parametrize("n_max", [2, 3])
def test_J_commutes_D_offdiag_no_m_lift_obstructed(n_max):
    """Even with m_lift = 0 in offdiag Dirac, J D = +D J does NOT hold.

    STRUCTURAL OBSTRUCTION (May-2026 finding):

    The offdiag Dirac is a SCALAR (multiplicative) hopping perturbation,
    NOT the genuine spinor-coupled Camporesi-Higuchi Dirac. Charge
    conjugation J on Camporesi-Higuchi spinors involves gamma-matrix
    multiplication on the spin index that is realized via a single
    permutation-phase structure on the (n, l, m_j) basis ONLY for the
    truthful (diagonal) CH Dirac, where the spin information is encoded
    in the eigenvalue.

    For the offdiag CH (a hand-crafted finite-n_max Hermitian
    perturbation), the cross-shell hoppings D[(n,l,m), (n+1,l+/-1,m')]
    are J-symmetric in their support pattern (the offdiag selection rule
    |Δl|=1, |Δm_j|<=1 is invariant under m_j -> -m_j). But the natural
    permutation-phase J = U K with sigma(l, m_j) = i^(2m_j)*(-1)^l
    does NOT have sigma_a = sigma_b for every offdiag-connected pair
    (a, b), so [U, D] != 0 at order 1.

    This test ASSERTS the failure, documenting the obstruction. The
    correct interpretation: the offdiag CH Dirac is NOT in the natural
    domain of the Camporesi-Higuchi charge conjugation; it is a model
    Dirac for SDP-bounded Connes distance, not the physical spinor
    Dirac.

    Resolution path: to get a physical Dirac that anticommutes/commutes
    with J cleanly, one must either:
      (a) Use the truthful CH (diagonal, eigenvalue +/- (n+1/2));
          this works (test above: test_J_commutes_D_full_dirac_truthful).
      (b) Implement the genuine 4-component Dirac as a 2x2 block
          off-diagonal in chirality with the SU(2) angular gradient
          in the off-diagonal block, plus the proper spinor charge
          conjugation matrix C = i*sigma_2 on the spin index.

    Path (b) is beyond Track 2's scope (it requires new spinor
    infrastructure beyond `full_dirac_operator_system.py`).
    """
    basis = full_dirac_basis(n_max)
    # Build offdiag Dirac with m_lift = 0
    D = camporesi_higuchi_offdiag_dirac_matrix(
        basis,
        diag_lifters=(1.0, 0.1, 0.0),  # n_lift, l_lift, m_lift=0
        offdiag_alpha=1.0,
        chirality_coupling=1.0,
    )
    J = build_J_full_dirac(n_max)
    ok, err = verify_J_D_relation(J, D, expected_sign=+1, tol=1e-9)
    assert not ok, (
        f"Expected J D != +D J on offdiag (structural obstruction); "
        f"got ok={ok}, max err {err:.3e}"
    )
    # Document the actual residual size
    print(f"\n  [n_max={n_max}, full_dirac, offdiag, m_lift=0] "
          f"|JD - DJ|_max = {err:.3e} (expected ~ 2 from offdiag rule)")


# ---------------------------------------------------------------------------
# J preserves O (operator system invariance)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2])
def test_J_preserves_O_weyl(n_max):
    """For a in O, J a J^{-1} should also be in O.

    This is a structural property of the Connes-vS truncation: J acts
    on the truncated operator system by an antilinear order-isomorphism,
    so it should map O to O.
    """
    op_sys = SpinorTruncatedOperatorSystem(n_max)
    J = build_J_weyl(n_max)
    ok, failures = verify_J_preserves_O(J, op_sys, tol=1e-9)
    if not ok:
        # Report first few failures for debugging
        msg = (
            f"J does NOT preserve O on Weyl sector at n_max={n_max}. "
            f"{len(failures)} of {len(op_sys.multiplier_matrices)} generators failed.\n"
        )
        for i, res in failures[:3]:
            msg += f"  generator {i} (label={op_sys.multiplier_labels[i]}) "
            msg += f"residual {res:.3e}\n"
        raise AssertionError(f"[n_max={n_max}, weyl] " + msg)


@pytest.mark.parametrize("n_max", [1, 2])
def test_J_preserves_O_full_dirac(n_max):
    """For a in O (full Dirac), J a J^{-1} should also be in O."""
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    J = build_J_full_dirac(n_max)
    ok, failures = verify_J_preserves_O(J, op_sys, tol=1e-9)
    if not ok:
        msg = (
            f"J does NOT preserve O on full Dirac sector at n_max={n_max}. "
            f"{len(failures)} of {len(op_sys.multiplier_matrices)} generators failed.\n"
        )
        for i, res in failures[:3]:
            msg += f"  generator {i} (label={op_sys.multiplier_labels[i]}) "
            msg += f"residual {res:.3e}\n"
        raise AssertionError(f"[n_max={n_max}, full_dirac] " + msg)


# ---------------------------------------------------------------------------
# Order-zero condition: [a, J b J^{-1}] = 0 for all a, b in O
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2])
def test_order_zero_weyl(n_max):
    """Test [a, J b J^{-1}] = 0 on Weyl sector.

    KNOWN RESULT (audit n_max=1: HOLDS trivially with 1 generator;
    n_max=2: FAILS on 124/196 pairs with max residual ~ 5.5e-2;
    n_max=3: FAILS on 2594/3025 pairs with max residual ~ 7.9e-2).

    The condition [a, J b J^{-1}] = 0 is the abelian-classical
    "multiplications commute" identity, which holds in the continuum
    for f, g in C^infty(M) since pointwise multiplication is commutative
    and J b J^{-1} is multiplication by conj(g) in the spinor bundle.

    AT FINITE n_max, the truncation P_{n_max} f P_{n_max} is no longer
    multiplication by a classical function -- it is an OPERATOR SYSTEM
    element (Connes-vS Sec 2). Two operator-system elements a, b in O
    with [a, J b J^{-1}] != 0 give a finite-resolution violation of
    order zero. This is structurally analogous to Connes-vS's failure
    of multiplicative closure (a*b NOT in O).

    The residual scales as the SO(4) Wigner-3j coupling strengths
    between Fock shells (which are O(1) numbers, not O(epsilon) small).
    The order-zero failure is therefore a finite-resolution artifact
    of the truncation -- it should disappear in the n_max -> infinity
    limit by the GH-convergence theorem (which IS the recovery of
    classical commuting multiplications).
    """
    op_sys = SpinorTruncatedOperatorSystem(n_max)
    J = build_J_weyl(n_max)
    ok, failures = verify_order_zero(J, op_sys.multiplier_matrices, tol=1e-9)
    n_total = len(op_sys.multiplier_matrices) ** 2
    if n_max == 1:
        # Trivial case: only identity multiplier
        assert ok, f"n_max=1 should pass trivially; got failures {len(failures)}"
    else:
        # Structural failure expected; document it
        assert not ok, (
            f"n_max={n_max} order-zero is expected to FAIL "
            f"(finite-resolution artifact); got ok={ok}"
        )
        n_fail = len(failures)
        max_res = max(r for _, _, r in failures)
        print(f"\n[n_max={n_max}, weyl] Order-zero: "
              f"{n_fail}/{n_total} pairs fail, max res {max_res:.3e}")


@pytest.mark.parametrize("n_max", [1, 2])
def test_order_zero_full_dirac(n_max):
    """Test [a, J b J^{-1}] = 0 on full Dirac sector.

    Same finite-resolution mechanism as the Weyl test above; expected
    to fail at n_max>=2 with residuals comparable to the Weyl case
    (since the full-Dirac multipliers are block-diagonal in chirality
    with two equal Weyl-sector blocks, the chirality doubling does not
    introduce new operator-system commutators).
    """
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    J = build_J_full_dirac(n_max)
    ok, failures = verify_order_zero(J, op_sys.multiplier_matrices, tol=1e-9)
    n_total = len(op_sys.multiplier_matrices) ** 2
    if n_max == 1:
        assert ok
    else:
        assert not ok, (
            f"n_max={n_max} order-zero is expected to FAIL "
            f"(finite-resolution artifact); got ok={ok}"
        )
        n_fail = len(failures)
        max_res = max(r for _, _, r in failures)
        print(f"\n[n_max={n_max}, full_dirac] Order-zero: "
              f"{n_fail}/{n_total} pairs fail, max res {max_res:.3e}")


# ---------------------------------------------------------------------------
# Order-one condition: [[D, a], J b J^{-1}] = 0
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2])
def test_order_one_full_dirac_truthful(n_max):
    """Test [[D, a], J b J^{-1}] = 0 on full Dirac with truthful CH.

    KNOWN RESULT: HOLDS at n_max=1 (trivially); FAILS at n_max=2 with
    43/196 = 22% pairs failing, max res 1.0e-1; FAILS at n_max=3 with
    1337/3025 = 44% pairs failing, max res 2.0e-1.

    Same finite-resolution mechanism as order-zero: order-one is
    automatic in the continuum classical limit, fails at finite n_max
    because the truncated multipliers and truncated [D, a] no longer
    correspond to multiplication by classical functions and their
    derivatives.
    """
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    D = camporesi_higuchi_full_dirac_matrix(op_sys.basis)
    J = build_J_full_dirac(n_max)
    ok, failures = verify_order_one(J, op_sys.multiplier_matrices, D, tol=1e-9)
    n_total = len(op_sys.multiplier_matrices) ** 2
    if n_max == 1:
        assert ok
    else:
        assert not ok, (
            f"n_max={n_max} order-one is expected to FAIL "
            f"(finite-resolution artifact); got ok={ok}"
        )
        n_fail = len(failures)
        max_res = max(r for _, _, r in failures)
        print(f"\n[n_max={n_max}, full_dirac, truthful] Order-one: "
              f"{n_fail}/{n_total} pairs fail, max res {max_res:.3e}")


# ---------------------------------------------------------------------------
# High-level audit
# ---------------------------------------------------------------------------


def test_audit_J_runs_weyl_n2():
    """Smoke test for high-level audit_J on Weyl sector."""
    result = audit_J(n_max=2, sector="weyl", verbose=False)
    assert result["dim_H"] == spinor_dim(2)
    assert result["U_unitary"][0]
    assert result["J_squared_minus_I"][0]
    assert result["J_commutes_D"][0]


def test_audit_J_runs_full_dirac_n2():
    """Smoke test for high-level audit_J on full Dirac sector."""
    result = audit_J(n_max=2, sector="full_dirac", D_mode="truthful", verbose=False)
    assert result["dim_H"] == full_dirac_dim(2)
    assert result["U_unitary"][0]
    assert result["J_squared_minus_I"][0]
    assert result["J_commutes_D"][0]


# ---------------------------------------------------------------------------
# Negative controls
# ---------------------------------------------------------------------------


def test_wrong_sign_check_J_squared_plus_I():
    """If we check J^2 = +I (wrong sign for KO-dim 3), it should FAIL."""
    J = build_J_weyl(2)
    ok, err = verify_J_squared(J, expected_sign=+1)
    assert not ok, "J^2 = +I should NOT hold for KO-dim 3"
    # The error should be at the order of 2 (since J^2 = -I means
    # |J^2 - (+I)| = |-I - I| = 2)
    assert err > 1.0, f"expected residual ~ 2 for wrong sign, got {err:.3e}"


def test_wrong_sign_check_J_anticommutes_D():
    """If we check J D = -D J (wrong sign for KO-dim 3), it should FAIL on
    truthful CH Dirac."""
    basis = full_dirac_basis(2)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    J = build_J_full_dirac(2)
    ok, err = verify_J_D_relation(J, D, expected_sign=-1)
    assert not ok, "J D = -D J should NOT hold for KO-dim 3"


# ---------------------------------------------------------------------------
# Slow audit at n_max=3 (informational)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_audit_J_full_dirac_n3():
    """Run the full audit at n_max=3 on full Dirac (truthful)."""
    result = audit_J(n_max=3, sector="full_dirac", D_mode="truthful",
                     verbose=True)
    assert result["dim_H"] == full_dirac_dim(3)
    assert result["U_unitary"][0]
    assert result["J_squared_minus_I"][0]
    assert result["J_commutes_D"][0]
