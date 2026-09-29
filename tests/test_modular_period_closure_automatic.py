"""The finite-cutoff 2pi period closure is AUTOMATIC for any integer-spectrum K.

Scope note (geovac/modular_hamiltonian.py, 2026-09-28, /qa trunk FULL #9) and
Paper 32 / Paper 42 (finite four-witness "literal identification" demoted):
K_geometric = 2 J_z is POSITED, not derived from a state; sigma_{2pi}(O) = O
holds for ANY operator O whenever K has integer spectrum, because
e^{2 pi i (k_a - k_b)} = 1.  The closure therefore cannot discriminate a
Bisognano-Wichmann (BW) modular flow from any other integer-spectrum
generator, and the verdict label STRONG_IDENTIFICATION means only "the
automatic closure held".

WRONG ANSWER REJECTED: "the period closure discriminates BW" -- i.e. that
only the BW rotation generator (eigenvalue two_m_j) passes.  These tests
feed a generator with no rotation meaning, K = 3*two_m_j + 2*l, and require
that it close the SAME witness battery to the same machine-precision
residual (the residuals are asserted, NOT the module's verdict label, which
would lock in a misleading name); a non-integer control K = sqrt(2)*two_m_j
must fail, so the closure check is live.

Fire-test plant (debug/qa/fire_test.py --plant-in geovac/modular_hamiltonian.py):
verify_modular_periodicity made BW-discriminating (adds the deviation of K
from the BW generator diag(two_m_j) to sigma_{2pi}(O)) -> FIRES.  A pure
relabel of the verdict string stays GREEN (by design).
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.full_dirac_operator_system import (
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)
from geovac.modular_hamiltonian import (
    HemisphericWedge,
    ModularHamiltonian,
    for_bisognano_wichmann,
)


def _mh_with_K(n_max: int, diag_K: np.ndarray) -> ModularHamiltonian:
    basis = full_dirac_basis(n_max)
    return ModularHamiltonian(
        wedge=HemisphericWedge(axis="hopf"),
        basis=basis,
        D=camporesi_higuchi_full_dirac_matrix(basis),
        kappa_g=1.0,
        K_geometric=np.diag(diag_K.astype(np.complex128)),
    )


def _non_bw_integer_K(n_max: int) -> np.ndarray:
    basis = full_dirac_basis(n_max)
    return np.array([3.0 * b.two_m_j + 2.0 * b.l for b in basis])


@pytest.mark.parametrize("n_max", [2, 3, 4])
def test_period_closure_is_automatic_for_any_integer_K(n_max):
    """A non-BW integer-spectrum K passes the witness battery with the same
    verdict as the BW generator: the closure does not discriminate BW."""
    K = _non_bw_integer_K(n_max)
    bw_diag = np.array([float(b.two_m_j) for b in full_dirac_basis(n_max)])
    assert not np.allclose(K, bw_diag)  # genuinely not the rotation generator
    other = _mh_with_K(n_max, K).verify_witness()
    bw = for_bisognano_wichmann(n_max).verify_witness()
    # Assert the closure RESIDUALS, not the module's verdict label: the label
    # STRONG_IDENTIFICATION is misleading (it means only "the automatic closure
    # held") and an honest relabel must not break this guard.
    for name, res in (("BW", bw), ("non-BW", other)):
        per_gen = res["periodicity_per_generator"]
        assert len(per_gen) > 0, f"n_max={n_max}: {name} battery ran no generators"
        assert all(g["residual"] < 1e-10 for g in per_gen), (
            f"n_max={n_max}: {name} closure residuals "
            f"{[g['residual'] for g in per_gen]}")
        assert res["max_periodicity_residual"] < 1e-10
    # same generator battery on both sides
    assert len(other["periodicity_per_generator"]) == len(bw["periodicity_per_generator"])


@pytest.mark.parametrize("n_max", [2, 3])
def test_period_closure_holds_on_arbitrary_operators(n_max):
    """sigma_{2pi}(O) = O on random dense complex O for the non-BW integer K,
    while a non-integer K = sqrt(2)*two_m_j fails (the check is live)."""
    rng = np.random.default_rng(n_max)
    dim = len(full_dirac_basis(n_max))
    mh_int = _mh_with_K(n_max, _non_bw_integer_K(n_max))
    two_mj = np.array([float(b.two_m_j) for b in full_dirac_basis(n_max)])
    mh_irr = _mh_with_K(n_max, np.sqrt(2.0) * two_mj)
    for _ in range(3):
        O = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
        ok, res = mh_int.verify_modular_periodicity(O)
        assert ok and res < 1e-10, f"integer non-BW K: residual {res:.2e}"
        ok_irr, res_irr = mh_irr.verify_modular_periodicity(O)
        assert not ok_irr and res_irr > 1e-2, f"non-integer K closed: {res_irr:.2e}"
