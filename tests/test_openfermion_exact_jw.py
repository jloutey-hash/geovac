"""Guard: the Jordan-Wigner builders keep terms that OpenFermion's EQ_TOLERANCE (1e-8) would delete.

OpenFermion deletes a term whenever its running coefficient falls below 1e-8 while operators are added, and a
single added term below 1e-8 never enters. He cc-pVTZ lost 16 real terms (21,623 -> 21,607), He n_max = 4 lost
144 and n_max = 5 lost 11,688 (CHANGELOG v6.7.23; Paper 14 tab:pauli). The fix: dict accumulation in
geovac.qubit_encoding's fermion builders and exact_jordan_wigner (transform at 2**20, scale back).

Each test names the wrong answer it rejects.
"""
import importlib.util
import os

import numpy as np
import pytest
from openfermion import FermionOperator, jordan_wigner

from geovac.qubit_encoding import build_fermion_op_from_integrals, exact_jordan_wigner

_spec = importlib.util.spec_from_file_location(
    "water_ccsd_closed_form", os.path.join(os.path.dirname(__file__), "test_water_ccsd.py"))
_w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_w)


def test_term_assembled_from_sub_tolerance_pieces_survives():
    """Two 6e-9 pieces of one fermion term: the true coefficient 1.2e-8 must survive.
    Rejects the plain OpenFermion path (each partial sum is deleted, the term vanishes)."""
    f = FermionOperator()
    for _ in range(2):
        f.terms[((0, 1), (1, 0))] = f.terms.get(((0, 1), (1, 0)), 0.0) + 6e-9
        f.terms[((1, 1), (0, 0))] = f.terms.get(((1, 1), (0, 0)), 0.0) + 6e-9
    q = exact_jordan_wigner(f)
    # c (a+_0 a_1 + a+_1 a_0) -> (c/2)(X0 X1 + Y0 Y1), c = 1.2e-8
    assert abs(q.terms.get(((0, "X"), (1, "X")), 0.0) - 6e-9) < 1e-20
    assert abs(q.terms.get(((0, "Y"), (1, "Y")), 0.0) - 6e-9) < 1e-20


def test_single_sub_tolerance_integral_is_kept():
    """h1[0,1] = h1[1,0] = 5e-9 enters the fermion operator and the qubit operator (X0 Z1 X2 coefficient 2.5e-9).
    Rejects fermion-side deletion (OpenFermion's += drops a lone 5e-9 term before JW sees it)."""
    h1 = np.diag([-1.0, -0.5])
    h1[0, 1] = h1[1, 0] = 5e-9
    eri = np.zeros((2, 2, 2, 2))
    f = build_fermion_op_from_integrals(h1, eri, 0.0)
    assert abs(f.terms[((0, 1), (2, 0))] - 5e-9) < 1e-24
    q = exact_jordan_wigner(f)
    assert abs(q.terms[((0, "X"), (1, "Z"), (2, "X"))] - 2.5e-9) < 1e-20


def _random_small_hamiltonian(seed: int, M: int = 4):
    """Real symmetric h1 and 8-fold symmetric chemist ERIs with entries spread over 1e-10..1e-6
    (where the drop bites), plus O(1) diagonals so cancellations occur."""
    rng = np.random.default_rng(seed)
    h1 = rng.standard_normal((M, M)) * 10.0 ** rng.uniform(-10, -6, (M, M))
    h1 = h1 + h1.T + np.diag(rng.uniform(-2, -0.5, M))
    eri = rng.standard_normal((M,) * 4) * 10.0 ** rng.uniform(-10, -6, (M,) * 4)
    for perm in [(1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)]:
        eri = eri + eri.transpose(perm)
    eri = (eri + eri.transpose(1, 0, 3, 2)) / 2
    for p in range(M):
        for q in range(M):
            eri[p, p, q, q] += 0.3 + 0.1 * (p + q)
    return h1, eri


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_count_and_norm_match_closed_form(seed):
    """Term count and 1-norm equal the closed-form counter (no fermion-operator expansion) on Hamiltonians
    whose integrals sit in 1e-10..1e-6. Rejects any accumulation drop, and a cutoff set above the real terms
    (e.g. 1e-8 at the end, which deletes the 1e-10..1e-8 terms the closed form keeps)."""
    h1, eri = _random_small_hamiltonian(seed)
    q = exact_jordan_wigner(build_fermion_op_from_integrals(h1, eri, 0.7))
    st = _w.jw_pauli_stats(h1, eri, 0.7, tol=1e-11)
    lam = sum(abs(c) for k, c in q.terms.items() if k)
    assert len(q.terms) == st["n_pauli"]
    assert abs(lam - st["one_norm"]) < 1e-12
    # the plain OpenFermion path loses terms on these inputs (shows the test exercises the drop)
    assert len(jordan_wigner(build_fermion_op_from_integrals(h1, eri, 0.7)).terms) < st["n_pauli"]


def test_clean_hamiltonian_bit_identical_to_plain_jw():
    """With all coefficients O(1), exact_jordan_wigner returns exactly OpenFermion's operator, coefficient for
    coefficient. Rejects a scale that is not a power of two, or a scale-back that alters clean outputs."""
    rng = np.random.default_rng(7)
    M = 3
    h1 = rng.uniform(0.2, 1.0, (M, M)); h1 = h1 + h1.T
    eri = rng.uniform(0.2, 1.0, (M,) * 4)
    eri = eri + eri.transpose(1, 0, 2, 3) + eri.transpose(0, 1, 3, 2) + eri.transpose(2, 3, 0, 1)
    f = build_fermion_op_from_integrals(h1, eri, 1.3)
    a, b = exact_jordan_wigner(f).terms, jordan_wigner(f).terms
    assert a.keys() == b.keys()
    assert all(a[k] == b[k] for k in a)


@pytest.mark.slow
def test_he_ccpvtz_count():
    """He cc-pVTZ: 21,623 terms incl. identity, = the closed-form count. Rejects the retired 21,607."""
    from geovac.gaussian_reference import build_qubit_hamiltonian, he_cc_pvtz
    d = he_cc_pvtz()
    _, q, n = build_qubit_hamiltonian(d)
    assert n == 21623
    assert n == _w.jw_pauli_stats(d["h1"], d["eri"], d["nuclear_repulsion"])["n_pauli"]


class _StubIndex:
    """Minimal LatticeIndex stand-in for JordanWignerEncoder: spatial H1 (sparse) and a physicist-notation
    ERI dict <ab|cd> = (ac|bd)."""

    def __init__(self, h1, chem_eri):
        import scipy.sparse as sps
        M = h1.shape[0]
        self.n_sp = 2 * M
        self.V_NN = 0.0
        self._H1_spatial = sps.csr_matrix(h1)
        self._eri = {(a, b, c, d): float(chem_eri[a, c, b, d])
                     for a in range(M) for b in range(M) for c in range(M) for d in range(M)
                     if chem_eri[a, c, b, d] != 0.0}


@pytest.mark.parametrize("seed", [4, 5])
def test_encoder_path_matches_closed_form(seed):
    """JordanWignerEncoder (the He library path behind the n_max = 4, 5 counts) gives the closed-form count and
    1-norm on 1e-10..1e-6 integrals. Rejects the encoder reverting to OpenFermion's += or to plain jordan_wigner."""
    from geovac.qubit_encoding import JordanWignerEncoder
    h1, eri = _random_small_hamiltonian(seed)
    q = JordanWignerEncoder(_StubIndex(h1, eri)).build_qubit_operator()
    st = _w.jw_pauli_stats(h1, eri, 0.0, tol=1e-11)
    lam = sum(abs(c) for k, c in q.terms.items() if k)
    assert len(q.terms) == st["n_pauli"]
    assert abs(lam - st["one_norm"]) < 1e-12
