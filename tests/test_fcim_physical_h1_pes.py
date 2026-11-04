"""Guard: FCI-M Table I / Sec. III -- the physical-h1 n_max = 2 LiH PES.

Paper: papers/group2_quantum_chemistry/paper_fci_molecules.tex,
Table ``tab:lih_energies``, Sec. ``sec:energy`` / ``sec:contraction`` and the
Fock-weighted scan Table ``tab:fock_scan``.

The paper's energetics are stated for the PHYSICAL one-electron Hamiltonian:
hydrogenic diagonal + cross-nuclear attraction, with the inter-atomic STO
overlap bridges as the only off-diagonal one-electron terms.  The archived
``MolecularLatticeIndex`` additionally carries the intra-atom graph hop
kappa*A (kappa = -1/16) inside each atomic block; that hop is the
non-physical term retired corpus-wide on 2026-09-25.  This module removes the
intra-atom block of the combined adjacency (bridges kept) and pins:

  * E(R) strictly monotone (rising with R) on R = 1.5 .. 8.0, no interior minimum;
  * E(8) < E_sep - 0.2 Ha, with E_sep = own-basis Li + own-basis H;
  * Boys-Bernardi BSSE = -0.0884 +- 1e-3 Ha, the same at every R tested;
  * lambda = 0.1 Fock-weighted correction: an interior minimum, three-point
    parabola vertex in [3.5, 4.0] bohr (paper: R_eq ~ 3.7).

Expected values are the paper's printed numbers, not values read back from the
code under test.  Fire-test: set ``_REMOVE_INTRA_HOP = False`` (the retired
hop-on build); the BSSE pin, the E(8) table pin and the R_eq window must fail.
"""
from __future__ import annotations

import contextlib
import io
import warnings
from typing import Dict, List

import numpy as np
import pytest
from scipy.sparse import block_diag

import geovac._archive.superseded.molecular_lattice_index as M

pytestmark = [pytest.mark.slow, pytest.mark.filterwarnings("ignore")]

# Physical one-electron Hamiltonian: drop the intra-atom kappa*A hop.
_REMOVE_INTRA_HOP = True

RS: List[float] = [1.5, 2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0]
LAM_RS: List[float] = [1.5, 3.0, 3.5, 4.0, 4.5, 8.0]
BSSE_RS: List[float] = [1.5, 3.015, 8.0]

# Paper Table I (n_max = 2, exact+True, physical h1), printed to 1e-4 Ha.
TABLE_I: Dict[float, float] = {
    1.5: -8.3793, 2.5: -8.0592, 3.015: -7.9899,
    4.0: -7.9312, 5.0: -7.8877, 8.0: -7.8258,
}
E_SEP_PAPER = -7.601          # paper: separated-atom reference, own basis
BSSE_PAPER = -0.0884          # paper: Delta E_BSSE = 0.088 Ha (lowering)


def _quiet(f, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f(*a, **k)


@pytest.fixture(scope="module")
def physical_build(request):
    MLI = M.MolecularLatticeIndex
    orig = MLI._build_combined_adjacency

    def nohop(self) -> None:
        orig(self)
        a_intra = block_diag([self._li_A.lattice.adjacency,
                              self._li_B.lattice.adjacency], format="csr")
        self._adjacency_combined = (self._adjacency_combined - a_intra).tocsr()
        self._adjacency_combined.eliminate_zeros()

    if _REMOVE_INTRA_HOP:
        MLI._build_combined_adjacency = nohop

    def energy(R: float, **kw) -> float:
        m = _quiet(MLI, Z_A=3, Z_B=1, nmax_A=2, nmax_B=2, R=R, n_electrons=4,
                   n_bridges=10, vee_method="slater_full",
                   fci_method="matrix", cross_atom_vee=True, **kw)
        return float(_quiet(m.compute_ground_state, n_states=1)[0][0])

    try:
        E0 = {R: energy(R) for R in RS}
        E01 = {R: energy(R, t_corr_lambda=0.1, t_corr_fock_weighted=True)
               for R in LAM_RS}
        bsse = {R: _quiet(M.compute_bsse_correction, Z_A=3, Z_B=1, nmax_A=2,
                          nmax_B=2, R=R, n_electrons_A=3, n_electrons_B=1,
                          vee_method="slater_full", fci_method="matrix")
                for R in BSSE_RS}
    finally:
        MLI._build_combined_adjacency = orig
    return E0, E01, bsse


def test_table_i_energies(physical_build):
    E0, _, _ = physical_build
    for R, ref in TABLE_I.items():
        assert abs(E0[R] - ref) < 2e-4, (R, E0[R], ref)


def test_pes_monotone_no_interior_minimum(physical_build):
    E0, _, _ = physical_build
    e = np.array([E0[R] for R in RS])
    assert np.all(np.diff(e) > 0), e      # E rises with R everywhere
    assert int(np.argmin(e)) == 0         # minimum at the innermost grid point


def test_large_R_below_separated_atoms(physical_build):
    E0, _, bsse = physical_build
    b = bsse[8.0]
    e_sep = float(b["E_A_own"] + b["E_B_own"])
    assert abs(e_sep - E_SEP_PAPER) < 1e-3, e_sep
    assert E0[8.0] < e_sep - 0.2, (E0[8.0], e_sep)


def test_bsse_value_and_R_independence(physical_build):
    _, _, bsse = physical_build
    vals = [float(bsse[R]["BSSE"]) for R in BSSE_RS]
    for R, v in zip(BSSE_RS, vals):
        assert abs(v - BSSE_PAPER) < 1e-3, (R, v)
    assert max(vals) - min(vals) < 1e-6, vals


def test_fock_weighted_lambda_01_interior_minimum(physical_build):
    _, E01, _ = physical_build
    e = np.array([E01[R] for R in LAM_RS])
    i = int(np.argmin(e))
    assert 0 < i < len(LAM_RS) - 1, (i, e)  # interior: repulsive inner wall
    x = np.array(LAM_RS[i - 1:i + 2])
    y = e[i - 1:i + 2]
    c2, c1, _ = np.polyfit(x, y, 2)
    r_eq = -c1 / (2.0 * c2)
    assert c2 > 0
    assert 3.5 <= r_eq <= 4.0, r_eq
