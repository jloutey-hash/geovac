"""Paper 14 (qubit encoding) scaling-exponent backing tests.

Closes the NO-TEST gap flagged in docs/claim_test_matrix.md (group4 pre-work,
2026-06-28): the four headline scaling exponents were never pinned from GeoVac data
--- the only prior test (test_qubit_encoding.test_fit_scaling) feeds fit_pauli_scaling
SYNTHETIC data (p = 2*q**3) and so validates the FITTER, not the framework's exponents.

These tests recompute each exponent from the live He JW encoding over n_max=2,3,4
(Q=10,28,60) and assert it matches Paper 14's published value within a band that
(a) brackets the published full-range (n_max=2..5) figure and the tractable 3-point
subset, and (b) is tight enough to catch a broken encoder (which would push N_Pauli
toward the Gaussian ~Q^4 or the naive ~Q^2). Paper 14 sec:scaling published values
(fit over n_max=2..5, Q=10..110):
    N_Pauli ~ Q^3.15 (R^2 0.9995),  lambda(1-norm) ~ Q^1.69 (R^2 0.997),
    QWC ~ Q^3.36 (R^2 1.000).
Recomputed 3-point (n_max=2..4) values, this pass: N_Pauli ~ Q^3.10, 1-norm ~ Q^1.67,
QWC ~ Q^3.356 (see debug/qa/group4_scaling_recompute.json). The n_max=5 point (Q=110)
is omitted: its build + QWC partition is too slow for a regression test.

All marked slow (the n_max=4 He build is ~20 s).
"""
from __future__ import annotations
import numpy as np
import pytest

from geovac.lattice_index import LatticeIndex
from geovac.qubit_encoding import JordanWignerEncoder, fit_pauli_scaling
from geovac.measurement_grouping import count_qwc_groups


NMAX = (2, 3, 4)  # Q = 10, 28, 60


def _l1_nonid(qop) -> float:
    return float(sum(abs(c) for t, c in qop.terms.items() if t))


def _he_encoders(nmax_values=NMAX):
    out = []
    for nmax in nmax_values:
        li = LatticeIndex(n_electrons=2, max_n=nmax, nuclear_charge=2,
                          vee_method='slater_full', h1_method='hybrid')
        out.append(JordanWignerEncoder(li))
    return out


@pytest.fixture(scope="module")
def he_sweep():
    """(Q, N_pauli, l1, qubit_op) for He at n_max=2,3,4."""
    rows = []
    for enc in _he_encoders():
        an = enc.analyze()
        qop = enc.build_qubit_operator()
        rows.append({'Q': an.n_qubits, 'N_pauli': an.n_pauli_terms,
                     'l1': _l1_nonid(qop), 'qop': qop})
    return rows


@pytest.mark.slow
def test_pauli_term_scaling_exponent(he_sweep):
    """N_Pauli ~ Q^alpha with alpha in [2.9, 3.3] (paper 3.15; 3-pt ~3.10)."""
    Q = np.array([r['Q'] for r in he_sweep], float)
    P = np.array([r['N_pauli'] for r in he_sweep], float)
    alpha, _ = fit_pauli_scaling(Q, P)
    assert 2.9 <= alpha <= 3.3, f"N_Pauli exponent {alpha:.3f} outside [2.9,3.3]"
    # the structural claim: well below the Gaussian molecular regime (3.9-4.3)
    assert alpha < 3.5, f"N_Pauli exponent {alpha:.3f} not below the Gaussian gap"


@pytest.mark.slow
def test_one_norm_scaling_exponent(he_sweep):
    """1-norm lambda ~ Q^alpha, sub-quadratic (paper 1.69; 3-pt ~1.67)."""
    Q = np.array([r['Q'] for r in he_sweep], float)
    L = np.array([r['l1'] for r in he_sweep], float)
    alpha, _ = fit_pauli_scaling(Q, L)
    assert 1.5 <= alpha <= 1.9, f"1-norm exponent {alpha:.3f} outside [1.5,1.9]"
    # the key fault-tolerant claim: sub-quadratic
    assert alpha < 2.0, f"1-norm exponent {alpha:.3f} not sub-quadratic"


@pytest.mark.slow
def test_qwc_group_scaling_exponent(he_sweep):
    """QWC measurement groups ~ Q^alpha with alpha in [3.1, 3.6] (paper 3.36)."""
    Q = np.array([r['Q'] for r in he_sweep], float)
    G = np.array([count_qwc_groups(r['qop']) for r in he_sweep], float)
    alpha, _ = fit_pauli_scaling(Q, G)
    assert 3.1 <= alpha <= 3.6, f"QWC exponent {alpha:.3f} outside [3.1,3.6]"


@pytest.mark.slow
def test_composed_pauli_scaling_exponent():
    """Composed (within-molecule, vary max_n) N_Pauli ~ Q^2.5 (paper)."""
    from geovac.composed_qubit import composed_lih_scaling_sweep
    sw = composed_lih_scaling_sweep(max_n_values=[1, 2, 3], verbose=False)
    Q = np.array([d['Q'] for d in sw['sweep_data']], float)
    P = np.array([d['N_pauli'] for d in sw['sweep_data']], float)
    alpha, _ = fit_pauli_scaling(Q, P)
    assert 2.2 <= alpha <= 2.8, f"composed exponent {alpha:.3f} outside [2.2,2.8] (paper 2.5)"
