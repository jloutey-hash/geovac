"""Paper 62, matched accuracy (Sec. matched, Table VI): regression guard.

Record: docs/predictions/2026-10-31_matched_accuracy_gaussian.md (MA, GO).  Written as a
separate step after the result was homed (CLAUDE.md Sec. 9, guard-writing rule) and
fire-tested against the wrong answers named in each test's docstring.

What is recomputed here and what is read.  The aug-cc-pVDZ water run is RECOMPUTED end
to end from the basis file stored beside this test: RHF, fc-CCSD(T), and the
Jordan-Wigner cost (Pauli count, non-identity 1-norm, qubits) over the RHF canonical
orbitals with the oxygen 1s frozen.  The cc-pVDZ water run is RECOMPUTED the same way,
so the mean-field share of the augmentation step and the equal-size Pauli count the
ratio is measured against are own numbers, not copies.  The Slater DZP target, the
cc-pVDZ+TZpol rung and the cc-pVTZ anchor are READ as the registered literals (registry
keys p62_dzp_fc, p62_dzp_lambda_fc, p62_tzpol_fc, p62_ccpvtz_fc; C21 pins them to the
paper), because each is its own run of five to twenty minutes and the ladder guard's
rule is one recomputation per record.

Cost: the aug-cc-pVDZ integrals take about five minutes and the solve about one; the
cc-pVDZ run is about a minute.  Everything here is @pytest.mark.slow.
"""
import importlib.util
import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_ma', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

AUG_BASIS = os.path.join(HERE, 'data', 'paper62_aug_cc_pvdz_H_O_Ne_2026_10_31.json')
DZ_BASIS = os.path.join(HERE, 'data', 'paper62_cc_pvdz_H_O_Ne_2026_10_31.json')

# READ (registered literals, C21-pinned; measured in their own runs, not here)
DZP_FC = -76.281739          # p62_dzp_fc, the target (CHANGELOG v5.20.11)
DZP_LAMBDA = 706.3           # p62_dzp_lambda_fc, non-identity 1-norm at Q = 46
TZPOL_FC = -76.266498        # p62_tzpol_fc (record MA)
CCPVTZ_FC = -76.332191       # p62_ccpvtz_fc, the anchor below the target

# RECOMPUTED here and asserted
AUG_FC = -76.273743          # p62_aug_fc
AUG_RHF = -76.041422         # record MA table
AUG_PAULI_FC = 982785        # p62_aug_pauli_fc (identity included)
AUG_LAMBDA_FC = 1901.1       # p62_aug_lambda_fc (identity excluded)
AUG_Q_FC = 80
CCPVDZ_FC = -76.241025       # p62_ccpvdz_fc
CCPVDZ_PAULI_FC = 107382     # h2o_ccpvdz_pauli, Trenev et al. 2025 Table 5
CCPVDZ_LAMBDA_FC = 674.3     # p62_ccpvdz_lambda_fc
CCPVDZ_Q_FC = 46

_CACHE = {}


def water_run(path):
    """RHF, fc-CCSD(T) and the fc Jordan-Wigner cost for water in the basis at `path`,
    the route of debug/matched_accuracy_gaussian.py and of CHANGELOG v5.20.11."""
    if path in _CACHE:
        return _CACHE[path]
    nuclei = w.water_nuclei()
    fns, T = w.gaussian_basis(path, nuclei)
    S, h, g = w.ao_integrals(fns, nuclei, T)
    enuc = w.vnn(nuclei)
    r = w.run(S, h, g, 10, enuc, frozen_core=1)
    _, eps, C = w.rhf(S, h, g, 10, enuc)
    hm = C.T @ h @ C
    gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
    hh, gg, ecore = w.freeze_core(hm, gm, 1)
    st = w.jw_pauli_stats(hh, gg, enuc + ecore)
    out = {'nbf': T.shape[1], 'RHF': float(r['RHF']), 'CCSD(T)_fc': float(r['CCSD(T)']),
           'qubits': 2 * hh.shape[0], 'n_pauli': int(st['n_pauli']),
           'one_norm': float(st['one_norm']), 'identity': float(st['identity'])}
    _CACHE[path] = out
    return out


@pytest.mark.slow
def test_aug_cc_pvdz_water_energy_is_the_registered_value():
    """The aug-cc-pVDZ water fc-CCSD(T) energy recomputed from the stored basis file is
    the registered -76.273743 (RHF -76.041422) at 41 functions.

    Wrong answers rejected: a basis file whose exponents or contractions have moved, a
    solver that has drifted (a changed (T) weight shows at 5e-5 Ha), a registry literal
    mis-transcribed from the run, and the all-electron energy reported as frozen-core
    (2.3 mHa apart).
    """
    a = water_run(AUG_BASIS)
    assert a['nbf'] == 41, a['nbf']
    assert abs(a['RHF'] - AUG_RHF) < 5e-6, a['RHF']
    assert abs(a['CCSD(T)_fc'] - AUG_FC) < 5e-6, a['CCSD(T)_fc']


@pytest.mark.slow
def test_aug_cc_pvdz_is_short_of_the_target_and_the_crossing_sits_before_cc_pvtz():
    """MA1 and MA4 and the crossing clause: the recomputed aug-cc-pVDZ energy is 8.0 mHa
    ABOVE the Slater DZP target, below the cc-pVDZ+TZpol rung, and the crossing lies
    between aug-cc-pVDZ and cc-pVTZ at 0.137 of the energy interval.

    Wrong answers rejected: a gap of the wrong sign (aug-cc-pVDZ reported at or below the
    target, which would move the crossing to the smaller basis), a target literal that
    has drifted, the TZpol rung not above aug-cc-pVDZ, and a crossing fraction that puts
    the matched-accuracy basis near cc-pVTZ instead of near aug-cc-pVDZ.
    """
    a = water_run(AUG_BASIS)
    gap = 1e3 * (a['CCSD(T)_fc'] - DZP_FC)
    assert gap > 0, gap
    assert abs(gap - 8.00) < 0.1, gap
    assert a['CCSD(T)_fc'] < TZPOL_FC, (a['CCSD(T)_fc'], TZPOL_FC)
    below = 1e3 * (DZP_FC - CCPVTZ_FC)
    assert below > 0
    frac = gap / (gap + below)
    assert 0.13 < frac < 0.145, frac


@pytest.mark.slow
def test_augmentation_step_is_correlation_led_with_rhf_share_below_half():
    """MA3 (reported): of the gain from cc-pVDZ to aug-cc-pVDZ at fc-CCSD(T), the RHF
    part is 0.45, under the registered bound 0.5 (the Slater step is 0.71 mean-field).
    Both runs are recomputed here, and the cc-pVDZ run must reproduce v5.20.11.

    Wrong answers rejected: the share read from all-electron energies, a cc-pVDZ run that
    no longer gives -76.241025 (data drift under the baseline), and a share that would
    make augmentation a mean-field effect.
    """
    a, d = water_run(AUG_BASIS), water_run(DZ_BASIS)
    assert d['nbf'] == 24
    assert abs(d['CCSD(T)_fc'] - CCPVDZ_FC) < 5e-6, d['CCSD(T)_fc']
    share = (a['RHF'] - d['RHF']) / (a['CCSD(T)_fc'] - d['CCSD(T)_fc'])
    assert share < 0.5, share
    assert abs(share - 0.45) < 0.01, share


@pytest.mark.slow
def test_aug_cc_pvdz_jw_cost_is_nine_times_the_terms_and_2_7_times_the_1_norm():
    """The reported qubit cost: aug-cc-pVDZ fc on 80 qubits has 982,785 Pauli terms
    (identity included) and a non-identity 1-norm of 1,901.1 Ha; the cc-pVDZ fc count
    recomputed beside it is Trenev's 107,382 on 46 qubits with 1-norm 674.3; the ratios
    are 9.15 (terms, against the equal-size count) and 2.69 (1-norm, against the
    Slater DZP's 706.3).

    Wrong answers rejected: the identity counted into the 1-norm (adds 125 Ha), the
    frozen core not applied in the qubit step (1,091,076 terms on 82 qubits), a cc-pVDZ
    count that drifts from the Trenev value, and the term ratio measured against the
    wrong baseline.
    """
    a, d = water_run(AUG_BASIS), water_run(DZ_BASIS)
    assert a['qubits'] == AUG_Q_FC and d['qubits'] == CCPVDZ_Q_FC, (a['qubits'], d['qubits'])
    assert a['n_pauli'] == AUG_PAULI_FC, a['n_pauli']
    assert d['n_pauli'] == CCPVDZ_PAULI_FC, d['n_pauli']
    assert abs(a['one_norm'] - AUG_LAMBDA_FC) < 0.05, a['one_norm']
    assert abs(d['one_norm'] - CCPVDZ_LAMBDA_FC) < 0.05, d['one_norm']
    assert abs(a['identity']) > 100.0          # the identity is a separate, excluded number
    assert abs(a['n_pauli'] / d['n_pauli'] - 9.15) < 0.01
    assert abs(a['one_norm'] / DZP_LAMBDA - 2.69) < 0.01
