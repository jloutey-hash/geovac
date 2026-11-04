"""Paper 19 owed guards (/qa group2 FULL 2026-11-03).

Claims (paper_19_coupled_composition.tex):

  (A) Sec. "Four-electron FCI": "in the other three schemes the one-body and
      electron-repulsion integrals do not change with R (identical at R = 2
      and 8 bohr), so their energies are a constant plus 3/R and cannot bind";
      the balanced configuration is "the only one of the four coupling
      schemes with an R-dependent electronic Hamiltonian".  Its ERIs are
      also R-independent ("the balanced electron-repulsion integrals are also
      R-independent").
      Schemes (Table tab:fci_balanced, LiH n_max = 2):
        decoupled  = composed, no PK, no cross-block
        composed   = composed with PK
        coupled CB = no PK + cross-block V_ee (build_coupled_hamiltonian)
        balanced   = build_balanced_hamiltonian

  (B) Reproducibility note: "at the fixed geometry R = 3.015 the default
      build gives ... -8.0545 Ha at n_max = 3 (frozen-core offset
      -7.2799 Ha)" on the corrected same-spin phase.  Here: the DirectCI4e
      (faithful=False) ground state E of the n_max = 3 balanced build obeys
      E + 7.2799 = -8.054515 +- 5e-5 Ha.  The historical (faithful=True)
      phase sits ~4 mHa lower and must be rejected.

The literals are the paper's printed numbers, not code output.
"""
from __future__ import annotations

import numpy as np
import pytest

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.coupled_composition import build_coupled_hamiltonian
from geovac.molecular_spec import lih_spec

FROZEN_CORE_OFFSET = -7.2799     # Ha, paper's stated offset
PAPER_E_NMAX3 = -8.054515        # Ha, paper (-8.0545), corrected phase


def _build(scheme: str, R: float, max_n: int = 2) -> dict:
    spec = lih_spec(R=R, max_n=max_n)
    if scheme == 'decoupled':
        return build_composed_hamiltonian(spec, pk_in_hamiltonian=False,
                                          verbose=False)
    if scheme == 'composed_pk':
        return build_composed_hamiltonian(spec, pk_in_hamiltonian=True,
                                          verbose=False)
    if scheme == 'coupled_cb':
        return build_coupled_hamiltonian(spec)
    if scheme == 'balanced':
        return build_balanced_hamiltonian(spec, R=R, n_grid_vne=8000,
                                          L_max=4, screened_cross_center=False,
                                          verbose=False)
    raise ValueError(scheme)


@pytest.mark.slow
@pytest.mark.parametrize('scheme', ['decoupled', 'composed_pk', 'coupled_cb'])
def test_paper19_nonbalanced_schemes_r_independent(scheme: str) -> None:
    """(A) h1 and ERIs bit-identical at R = 2 and R = 8 for the three
    non-balanced schemes; only V_NN moves (a nonempty-control on R)."""
    a, b = _build(scheme, 2.0), _build(scheme, 8.0)
    assert a['h1'].shape == b['h1'].shape
    # control: the R argument actually reached the build
    assert abs(a['nuclear_repulsion'] - b['nuclear_repulsion']) > 1.0
    dh = float(np.abs(a['h1'] - b['h1']).max())
    de = float(np.abs(a['eri'] - b['eri']).max())
    assert dh == 0.0, f"{scheme}: h1 changes with R (max |dh1| = {dh:.3e})"
    assert de == 0.0, f"{scheme}: ERIs change with R (max |deri| = {de:.3e})"


@pytest.mark.slow
def test_paper19_balanced_h1_r_dependent_eri_not() -> None:
    """(A) the balanced scheme is the one with an R-dependent h1 (bond-scale
    change, not noise); its ERIs are R-independent."""
    a, b = _build('balanced', 2.0), _build('balanced', 8.0)
    dh = float(np.abs(a['h1'] - b['h1']).max())
    de = float(np.abs(a['eri'] - b['eri']).max())
    assert dh > 0.1, f"balanced h1 barely moves with R (max |dh1| = {dh:.3e})"
    assert de == 0.0, f"balanced ERIs change with R (max |deri| = {de:.3e})"


@pytest.mark.slow
def test_paper19_balanced_nmax3_energy_corrected_phase() -> None:
    """(B) n_max = 3 balanced LiH at R = 3.015, corrected phase:
    E + 7.2799 = -8.054515 +- 5e-5 Ha.  (~12 min: the n_max = 3 build.)"""
    from geovac.balanced_direct_ci import DirectCI4e

    R = 3.015
    ham = _build('balanced', R, max_n=3)
    out = DirectCI4e(ham['h1'], ham['eri'], ham['nuclear_repulsion'],
                     faithful=False).ground_state(tol=1e-10, verbose=False)
    assert out['converged']
    E_paper = out['E'] - FROZEN_CORE_OFFSET
    assert abs(E_paper - PAPER_E_NMAX3) < 5e-5, (
        f"E + 7.2799 = {E_paper:.6f} Ha, paper -8.0545 (+-5e-5)")
