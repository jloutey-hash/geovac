"""GUARD: helium's physical Hamiltonian vs the retired graph-hop model (2026-09-25).

The "graph-native He CI 0.216% (n_max=7)" headline came from
``build_graph_native_fci``'s default h1, which adds a (1/16)*adjacency hop
that is not part of He's Hamiltonian (in the k=Z hydrogenic basis T+V_ne is
exactly diagonal).  The hop lowers He by ~49 mHa -- more than He's whole
correlation energy (42 mHa) -- so the two are different Hamiltonians and must
never be conflated.

What wrong answers this guard rejects:
  * the physical path (h1_method='exact') silently acquiring the hop
    -> its pinned energies move by ~50 mHa;
  * the physical path being wired back to the hybrid default, or the default
    silently losing its hop (both paths coincide) -> the >40 mHa separation
    assertion fails;
  * a "physical" h1 with any off-diagonal element.
Fire-tested 2026-09-25 by re-enabling the hop inside the include_hop=False
branch of casimir_ci._build_graph_h1 (see CHANGELOG / debug memo).
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.casimir_ci import build_graph_native_fci, _build_graph_h1

E_HE = -2.903724377034

# Physical (h1 = -Z^2/2n^2 only) singlet ground-state energies, measured
# 2026-09-25 by debug/he_graph_regression.py and cross-checked against
# build_fci_matrix(k=Z) (1e-15) and LatticeIndex h1_method='exact' (<0.2 mHa).
E_PHYS = {2: -2.833405, 3: -2.842908}


def _e0(**kw) -> float:
    return float(np.linalg.eigvalsh(build_graph_native_fci(Z=2, **kw))[0])


@pytest.mark.parametrize("n_max", [2, 3])
def test_physical_path_energy_pinned(n_max):
    E = _e0(n_max=n_max, h1_method='exact')
    assert abs(E - E_PHYS[n_max]) < 1e-6, (
        f"physical He energy moved: {E:.6f} vs {E_PHYS[n_max]:.6f} "
        "(did the graph hop leak into h1_method='exact'?)")


@pytest.mark.parametrize("n_max", [2, 3])
def test_hop_model_is_not_he(n_max):
    """The hop moves the energy by > 40 mHa: model and He are different."""
    shift = _e0(n_max=n_max) - _e0(n_max=n_max, h1_method='exact')
    assert shift < -0.040, f"hop shift only {1e3*shift:+.1f} mHa -- paths conflated"


def test_physical_h1_has_no_off_diagonal():
    h1, _ = _build_graph_h1(2, 3, include_hop=False)
    assert not np.any(h1 - np.diag(np.diag(h1)))
