"""Tests for geovac.spinor_certificate (Track T5, Dirac-on-S³ Tier 2 sprint).

Certifies that T2's symbolic spin-orbit Hamiltonian block lives in the
Paper 18 §IV spinor-intrinsic taxonomic ring

    Q(α²) × {γ = √(1 − (Zα)²)} × {hydrogenic radial seeds}

and that the certifier rejects deliberately contaminated coefficients.

All tests use exact sympy arithmetic. No floating-point comparisons.
"""

from __future__ import annotations

import pytest
import sympy as sp
from sympy import Rational, Symbol

from geovac.dirac_matrix_elements import (
    DiracLabel,
    Z_sym,
    alpha_sym,
    gamma_rel,
)
from geovac.spin_orbit import (
    build_so_hamiltonian_block,
    so_diagonal_matrix_element,
)
from geovac.spinor_certificate import (
    SpinorTaxonomyError,
    certify_so_block,
    default_allowed_symbols,
    verify_spinor_pi_free,
)


# ---------------------------------------------------------------------------
# Positive controls: T2 symbolic Hamiltonian certifies cleanly.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_max", [1, 2, 3, 4])
def test_so_block_symbolic_certifies(n_max):
    """T2 ``build_so_hamiltonian_block`` with symbolic (Z, α) passes."""
    block = build_so_hamiltonian_block(n_max=n_max, Z=Z_sym, alpha=alpha_sym)
    assert verify_spinor_pi_free(block) is True


def test_certify_so_block_convenience_wrapper():
    """The one-liner ``certify_so_block`` returns True at n_max=3."""
    assert certify_so_block(n_max=3) is True


def test_so_block_diag_individual_elements():
    """Every diagonal element of the T2 block individually passes."""
    block = build_so_hamiltonian_block(n_max=3, Z=Z_sym, alpha=alpha_sym)
    for label, coeff in block.diag.items():
        assert verify_spinor_pi_free(coeff) is True, (
            f"coefficient for {label} failed: {coeff!r}"
        )


def test_single_matrix_element_2p32():
    """⟨2p_{3/2} | H_SO | 2p_{3/2}⟩ = +Z^4 α² / 96 classifies."""
    coeff = so_diagonal_matrix_element(n=2, kappa=-2, Z=Z_sym, alpha=alpha_sym)
    assert verify_spinor_pi_free(coeff) is True


def test_single_matrix_element_2p12():
    """⟨2p_{1/2} | H_SO | 2p_{1/2}⟩ = -Z^4 α² / 48 classifies."""
    coeff = so_diagonal_matrix_element(n=2, kappa=+1, Z=Z_sym, alpha=alpha_sym)
    assert verify_spinor_pi_free(coeff) is True


def test_s_state_zero_certifies():
    """Kramers l=0 zero (sympy Integer(0)) certifies trivially."""
    coeff = so_diagonal_matrix_element(n=1, kappa=-1, Z=Z_sym, alpha=alpha_sym)
    assert coeff == 0
    assert verify_spinor_pi_free(coeff) is True


def test_gamma_rel_symbol_accepted():
    """γ = √(1 − (Zα)²) is in the allowed ring at the symbol level."""
    expr = gamma_rel * alpha_sym**2 * Rational(7, 16)
    assert verify_spinor_pi_free({"k": expr}) is True


def test_gamma_rel_explicit_sqrt_accepted():
    """The explicit sqrt form of γ also certifies (sqrt is Pow(x, 1/2))."""
    gamma_expanded = sp.sqrt(1 - (Z_sym * alpha_sym) ** 2)
    assert verify_spinor_pi_free({"k": gamma_expanded * Rational(1, 3)}) is True


# ---------------------------------------------------------------------------
# Negative controls: deliberate contamination MUST be rejected.
# ---------------------------------------------------------------------------


def test_rejects_bare_pi():
    """π alone in a coefficient triggers SpinorTaxonomyError."""
    with pytest.raises(SpinorTaxonomyError, match="(?i)π|pi"):
        verify_spinor_pi_free({"bad": sp.pi})


def test_rejects_pi_squared():
    """π² alone triggers rejection."""
    with pytest.raises(SpinorTaxonomyError):
        verify_spinor_pi_free({"bad": sp.pi ** 2})


def test_rejects_pi_in_so_block():
    """Injecting π into a T2-shaped dict triggers rejection.

    The contamination scenario specified in the T5 task: take the
    symbolic SO block and replace one coefficient with something
    π-tainted; confirm the certifier catches it.
    """
    block = build_so_hamiltonian_block(n_max=2, Z=Z_sym, alpha=alpha_sym)
    diag = dict(block.diag)
    # Pick any label at l>0 to overwrite with a π-contaminated value.
    target = next(lab for lab in diag if lab.kappa != -1)
    diag[target] = sp.pi * alpha_sym ** 2 / 96
    with pytest.raises(SpinorTaxonomyError):
        verify_spinor_pi_free(diag)


def test_rejects_zeta():
    """ζ(3) in a coefficient is rejected (odd-zeta tier, not spinor-intrinsic)."""
    with pytest.raises(SpinorTaxonomyError, match="(?i)zeta"):
        verify_spinor_pi_free({"bad": sp.zeta(3) * alpha_sym ** 2})


def test_rejects_log():
    """log(Z) is rejected."""
    with pytest.raises(SpinorTaxonomyError, match="(?i)log"):
        verify_spinor_pi_free({"bad": sp.log(Z_sym)})


def test_rejects_expint():
    """E_1 (expint) is rejected — that's a Level-2 embedding constant."""
    with pytest.raises(SpinorTaxonomyError):
        verify_spinor_pi_free({"bad": sp.expint(1, alpha_sym)})


def test_rejects_unknown_symbol():
    """A Symbol not in the registry is rejected."""
    xyz = Symbol("xyz")
    with pytest.raises(SpinorTaxonomyError, match="(?i)unregistered"):
        verify_spinor_pi_free({"bad": xyz * alpha_sym})


def test_extra_symbols_accepted():
    """User-declared extra symbols can be passed through."""
    xyz = Symbol("xyz")
    assert verify_spinor_pi_free(
        {"ok": xyz * alpha_sym ** 2 * Rational(1, 2)},
        extra_symbols=[xyz],
    ) is True


def test_raise_on_violation_false_returns_False():
    """raise_on_violation=False returns False without raising."""
    ok = verify_spinor_pi_free(
        {"bad": sp.pi * alpha_sym ** 2},
        raise_on_violation=False,
    )
    assert ok is False


# ---------------------------------------------------------------------------
# D1 π-free tests must still pass (regression check).
# ---------------------------------------------------------------------------


def test_d1_verify_pi_free_still_passes():
    """The scalar D1 verify_pi_free on the spinor label generator still
    works — the new module does NOT modify dirac_s3.py.
    """
    from geovac.dirac_s3 import verify_pi_free
    assert verify_pi_free(6, sector="dirac") is True
    assert verify_pi_free(6, sector="weyl") is True


# ---------------------------------------------------------------------------
# Input shape robustness.
# ---------------------------------------------------------------------------


def test_accepts_single_expr():
    """A bare sympy Expr is accepted as input."""
    assert verify_spinor_pi_free(alpha_sym ** 2 * Rational(1, 16)) is True


def test_accepts_iterable():
    """An iterable of Exprs is accepted."""
    exprs = [alpha_sym ** 2 * Rational(1, n) for n in range(1, 5)]
    assert verify_spinor_pi_free(exprs) is True


def test_accepts_block_via_diag_attr():
    """SOHamiltonianBlock is unpacked through its .diag attribute."""
    block = build_so_hamiltonian_block(n_max=2, Z=Z_sym, alpha=alpha_sym)
    assert verify_spinor_pi_free(block) is True


def test_default_registry_contents():
    """The default allowed-symbol registry is exactly {Z, α, γ}."""
    allowed = default_allowed_symbols()
    names = {s.name for s in allowed}
    assert names == {"Z", "alpha", "gamma_rel"}
