"""Paper 14 (qubit encoding) scaling-exponent backing tests.

Closes the NO-TEST gap flagged in docs/claim_test_matrix.md (group4 pre-work,
2026-06-28): the four headline scaling exponents were never pinned from GeoVac data
--- the only prior test (test_qubit_encoding.test_fit_scaling) feeds fit_pauli_scaling
SYNTHETIC data (p = 2*q**3) and so validates the FITTER, not the framework's exponents.

These tests recompute each exponent from the live He JW encoding over n_max=2,3,4
(Q=10,28,60) and assert it matches Paper 14's published value within a band that
(a) brackets the published full-range (n_max=2..5) figure and the tractable 3-point
subset, and (b) is tight enough to catch a broken encoder (which would push N_Pauli
toward the Gaussian ~Q^4 or the naive ~Q^2). Paper 14 published values, CORRECTED 2026-08-30 under the exact
global-M_L ERI rule (fit over n_max=2..5, Q=10..110):
    N_Pauli ~ Q^3.773 (R^2 1.0000, 4 pts)
    lambda  ~ Q^1.774 (R^2 0.9997, 4 pts)
    QWC     ~ Q^4.013 (R^2 0.9976, 3 pts -- the Q=110 point carries 2.4M
              Pauli terms against an O(N^2) grouping and is infeasible)
Retired pair-diagonal values were 3.15 / 1.69 / 3.36.  Note the QWC exponent
now sits ABOVE the Pauli exponent and inside the Gaussian O(Q^4-5) band, so
no measurement-group scaling advantage is claimed.
Recomputed 3-point (n_max=2..4) subset, which is what these tests run:
N_Pauli ~ Q^3.779, 1-norm ~ Q^1.792, QWC ~ Q^4.013
(see debug/data/exact_rule_lambda_qwc.json).  The n_max=5 point (Q=110) is
omitted here: its build alone is ~28 min.

All marked slow.  MEASURED runtime of the four slow tests: ~69 min
(4175 s and 4144 s on two independent runs), dominated by the
O(N^2) QWC grouping at n_max=4 (250,403 terms, 86,224 groups).  The
earlier '~20 s' estimate predated the exact-rule correction, which
multiplied the term counts several-fold.  A fast default-run pin on
the two cheapest points is provided separately below.
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


def test_small_basis_points_are_pinned():
    """Fast default-run guard on the two cheapest points of the fit.

    The exponent tests above are slow-marked (the n_max=4 QWC grouping runs
    ~1 h), so without this the headline exponents would have no
    default-run protection at all -- backing that exists but never
    executes.  These two points are the inputs the fit is built from: if
    the encoder regresses, they move and this fires in seconds.

    MEASURED 2026-08-30 under the exact ERI rule.  Retired pair-diagonal
    values were N=120 / lambda=11.29 at Q=10 and N=2,659 / lambda=78.36 at
    Q=28, so the pins also lock out a silent revert to the retired rule.
    """
    import warnings
    expected = {2: (10, 287, 11.175), 3: (28, 14078, 74.207)}
    for nmax, (q, n_ref, lam_ref) in expected.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            li = LatticeIndex(n_electrons=2, max_n=nmax, nuclear_charge=2,
                              vee_method='slater_full', h1_method='hybrid')
        enc = JordanWignerEncoder(li)
        n_pauli = sum(1 for k in enc.build_qubit_operator().terms if k)
        lam = _l1_nonid(enc.build_qubit_operator())
        assert n_pauli == n_ref, (
            f"n_max={nmax} (Q={q}): N_Pauli={n_pauli}, expected {n_ref}"
        )
        assert lam == pytest.approx(lam_ref, rel=1e-4), (
            f"n_max={nmax} (Q={q}): lambda={lam:.3f}, expected {lam_ref}"
        )


@pytest.mark.slow
def test_pauli_term_scaling_exponent(he_sweep):
    """N_Pauli ~ Q^alpha, exact rule: 3.773 (4 pts) / 3.779 (3-pt subset).

    Retired pair-diagonal value 3.15; band re-centred 2026-08-30.
    """
    Q = np.array([r['Q'] for r in he_sweep], float)
    P = np.array([r['N_pauli'] for r in he_sweep], float)
    alpha, _ = fit_pauli_scaling(Q, P)
    assert alpha >= 3.6, f"N_Pauli exponent {alpha:.3f} below 3.6"
    # The structural claim is now NARROW, not comfortable: the exponent sits
    # just below the Gaussian molecular regime (3.9-4.3), where the retired
    # rule had put it 0.75 below.  Guard the side that matters.
    assert alpha < 3.9, f"N_Pauli exponent {alpha:.3f} reached the Gaussian band"


@pytest.mark.slow
def test_one_norm_scaling_exponent(he_sweep):
    """1-norm lambda ~ Q^alpha, exact rule: 1.774 (4 pts) / 1.792 (3-pt).

    Retired pair-diagonal value 1.69.  Sub-quadratic scaling is the one
    fault-tolerant-relevant advantage that survives the correction.
    """
    Q = np.array([r['Q'] for r in he_sweep], float)
    L = np.array([r['l1'] for r in he_sweep], float)
    alpha, _ = fit_pauli_scaling(Q, L)
    # Band EXCLUDES the retired 1.69: a regression restoring the q-sign
    # error would put lambda back there, and a band admitting it would give
    # this test no discrimination at all.  Live: 1.774 (4 pts) / 1.792 (3).
    assert 1.72 <= alpha <= 1.90, (
        f"1-norm exponent {alpha:.3f} outside [1.72,1.90] "
        f"(retired pair-diagonal value 1.69 is deliberately excluded)"
    )
    # The fault-tolerant claim is not merely "sub-quadratic" -- that is
    # implied by the band -- but that lambda grows strictly slower than the
    # term count, which is what makes Trotter depth favourable.
    P = np.array([r['N_pauli'] for r in he_sweep], float)
    alpha_terms, _ = fit_pauli_scaling(Q, P)
    assert alpha < alpha_terms - 1.5, (
        f"1-norm exponent {alpha:.3f} is not well below the term exponent "
        f"{alpha_terms:.3f} -- the coefficient-concentration argument fails"
    )


@pytest.mark.slow
def test_qwc_group_scaling_exponent(he_sweep):
    """QWC groups ~ Q^alpha, exact rule: 4.013 (3 pts).

    Retired pair-diagonal value 3.36.  The band moved by more than its own
    width, which is why this test had to be re-centred rather than widened:
    the old [3.1, 3.6] gate would FAIL on the live value.
    """
    Q = np.array([r['Q'] for r in he_sweep], float)
    G = np.array([count_qwc_groups(r['qop']) for r in he_sweep], float)
    alpha, _ = fit_pauli_scaling(Q, G)
    assert 3.8 <= alpha <= 4.2, f"QWC exponent {alpha:.3f} outside [3.8,4.2]"
    # The load-bearing structural fact: QWC is STEEPER than the term count,
    # so no measurement-group scaling advantage exists under the exact rule.
    P = np.array([r['N_pauli'] for r in he_sweep], float)
    alpha_p, _ = fit_pauli_scaling(Q, P)
    assert alpha > alpha_p, (
        f"QWC exponent {alpha:.3f} is not above the Pauli exponent "
        f"{alpha_p:.3f} -- the paper's 'no measurement-group advantage' "
        f"claim rests on this ordering"
    )


@pytest.mark.slow
def test_composed_pauli_scaling_exponent():
    """Composed within-molecule N_Pauli ~ Q^alpha, exact rule.

    This test fits THREE points (n_max=1,2,3) and measures 3.1685, matching
    the paper's small-basis 3.17.  Do not confuse it with the TWO-point
    (n_max=1,2) value 2.8163, which is exact and identical across molecules
    because it is forced by 1 + log_5(27.90/1.5) -- a different fit over a
    different range, and outside this band by construction.
    Retired pair-diagonal value: 2.5.
    """
    from geovac.composed_qubit import composed_lih_scaling_sweep
    sw = composed_lih_scaling_sweep(max_n_values=[1, 2, 3], verbose=False)
    Q = np.array([d['Q'] for d in sw['sweep_data']], float)
    P = np.array([d['N_pauli'] for d in sw['sweep_data']], float)
    alpha, _ = fit_pauli_scaling(Q, P)
    # MEASURED 3-point (n_max=1,2,3): 3.1685, matching the paper's
    # within-molecule 3.17.  Retired pair-diagonal band was [2.2,2.8].
    assert 3.0 <= alpha <= 3.35, (
        f"composed exponent {alpha:.3f} outside [3.0,3.35] "
        f"(paper 3.17; retired pair-diagonal value was 2.5)"
    )
