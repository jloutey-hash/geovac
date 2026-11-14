"""LiH STO-3G raw Jordan-Wigner baseline (Papers 14, 19, 20; group-4 synthesis).

Pins the Gaussian comparison row (631 terms, 16.48 Ha, 177 QWC groups at R = 3.015 bohr)
to a build that exists in the repository: geovac.gaussian_reference.lih_sto3g.

Wrong answers rejected:
  * integrals that disagree with an independent code -- the 1.45 A build must reproduce
    OpenFermion's shipped LiH STO-3G file (RHF energy, terms, 1-norm, groups);
  * the round-off insertion order -- it gives 179 groups, not 177;
Record: docs/predictions/2026-11-13_lih_sto3g_rebuild.md.
"""
import os

import numpy as np
import pytest

openfermion = pytest.importorskip("openfermion")

from geovac.gaussian_reference import lih_sto3g
from geovac.measurement_grouping import count_qwc_groups
from geovac.qubit_encoding import build_fermion_op_from_integrals, exact_jordan_wigner

BOHR_A = 0.529177210903


def _measure(h1, eri, v_nn):
    qop = exact_jordan_wigner(build_fermion_op_from_integrals(h1, eri, v_nn))
    lam = float(sum(abs(c) for c in qop.terms.values()))
    return len(qop.terms), lam, count_qwc_groups(qop)


@pytest.fixture(scope="module")
def baseline():
    return lih_sto3g(3.015)


@pytest.fixture(scope="module")
def control():
    from openfermion import MolecularData
    d = os.path.join(os.path.dirname(openfermion.__file__), "testing", "data")
    mol = MolecularData(filename=os.path.join(d, "H1-Li1_sto-3g_singlet_1.45"))
    mol.load()
    return lih_sto3g(1.45 / BOHR_A), mol


def test_control_matches_openfermion_file(control):
    ours, mol = control
    assert abs(ours["hf_energy"] - mol.hf_energy) < 1e-8
    ref_eri = np.asarray(mol.two_body_integrals).transpose(0, 3, 1, 2)
    t_ref, l_ref, g_ref = _measure(mol.one_body_integrals, ref_eri, float(mol.nuclear_repulsion))
    t, lam, g = _measure(ours["h1"], ours["eri"], ours["nuclear_repulsion"])
    assert (t, g) == (t_ref, g_ref) == (631, 178)
    assert abs(lam - l_ref) < 1e-6


def test_baseline_values(baseline):
    assert abs(baseline["hf_energy"] - (-7.862009)) < 1e-5
    t, lam, g = _measure(baseline["h1"], baseline["eri"], baseline["nuclear_repulsion"])
    assert t == 631
    assert abs(lam - 16.4768) < 1e-3
    assert g == 177
