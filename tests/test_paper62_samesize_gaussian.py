"""Paper 62, same-size Gaussian control (Sec. dzp, Table samesize): regression guard.

Record: docs/predictions/2026-11-11_samesize_gaussian_control.md (SZ, all five held).
Written as a separate step after the result was homed (CLAUDE.md Sec. 9, guard-writing
rule) and fire-tested against the wrong answers named in each test's docstring.

What is recomputed here and what is read.  The ANO-R1, pcseg-1 and def2-SVP water runs
are RECOMPUTED end to end from the basis files stored beside this test (Basis Set
Exchange downloads of 2026-11-11, exponents and coefficients only): RHF, fc-CCSD(T),
and the Jordan-Wigner cost over the RHF canonical orbitals with the oxygen 1s frozen.
The Slater DZP, cc-pVDZ and the best atom-set ladder rung are READ as registered
literals (registry keys p62_dzp_fc, p62_ccpvdz_fc, p62_r5h3d_fc, p62_dzp_lambda_fc,
p62_ccpvdz_lambda_fc; C21 pins them to the paper); cc-pVDZ is recomputed by
tests/test_paper62_matched_accuracy.py.

Cost: each 24-function run is about half a minute.  Everything here is
@pytest.mark.slow.
"""
import importlib.util
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_sz', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

ANO_BASIS = os.path.join(HERE, 'data', 'paper62_ano_r1_H_O_Ne_2026_11_11.json')
PCSEG_BASIS = os.path.join(HERE, 'data', 'paper62_pcseg_1_H_O_Ne_2026_11_11.json')
DEF2_BASIS = os.path.join(HERE, 'data', 'paper62_def2_svp_H_O_Ne_2026_11_11.json')

# READ (registered literals, C21-pinned; measured in their own runs, not here)
DZP_FC = -76.281739          # p62_dzp_fc
DZP_LAMBDA = 706.3           # p62_dzp_lambda_fc
CCPVDZ_FC = -76.241025       # p62_ccpvdz_fc
CCPVDZ_LAMBDA = 674.3        # p62_ccpvdz_lambda_fc
R5H3D_FC = -76.260011        # p62_r5h3d_fc, best atom-set rung (3d set on the atom)

# RECOMPUTED here and asserted
ANO_FC = -76.264667          # p62_anor1_fc
ANO_RHF = -76.047820
ANO_LAMBDA = 625.36          # p62_anor1_lambda_fc
PCSEG_FC = -76.217263        # p62_pcseg1_fc
DEF2_FC = -76.174699         # p62_def2svp_fc
PAULI_24 = 107382            # every 24-function basis, fc, Q = 46

_CACHE = {}


def water_run(path):
    """RHF, fc-CCSD(T) and the fc Jordan-Wigner cost for water in the basis at `path`,
    the route of debug/samesize_gaussian_control.py (via debug/matched_accuracy_gaussian.py)."""
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
           'one_norm': float(st['one_norm'])}
    _CACHE[path] = out
    return out


@pytest.mark.slow
def test_ano_r1_water_energy_is_the_registered_value():
    """ANO-R1 water recomputed from the stored basis file: 24 functions, RHF -76.047820,
    fc-CCSD(T) -76.264667.

    Wrong answers rejected: a general contraction read as its first column only (fewer
    functions), a basis file whose exponents or coefficients have moved, a drifted
    solver, and the all-electron energy reported as frozen-core.
    """
    a = water_run(ANO_BASIS)
    assert a['nbf'] == 24, a['nbf']
    assert abs(a['RHF'] - ANO_RHF) < 5e-6, a['RHF']
    assert abs(a['CCSD(T)_fc'] - ANO_FC) < 5e-6, a['CCSD(T)_fc']


@pytest.mark.slow
def test_slater_dzp_margin_over_ano_r1_is_17_1_mha():
    """SZ2 and SZ3: ANO-R1 is 23.6 mHa below cc-pVDZ, and the Slater DZP is 17.1 mHa
    below ANO-R1, the margin the paper states against the best same-size Gaussian.

    Wrong answers rejected: the margin quoted against cc-pVDZ (40.7), a sign flip
    (ANO-R1 at or below the Slater DZP, which would erase the equal-size result), and
    ANO-R1 reported above cc-pVDZ.
    """
    a = water_run(ANO_BASIS)
    below_dz = 1e3 * (CCPVDZ_FC - a['CCSD(T)_fc'])
    margin = 1e3 * (a['CCSD(T)_fc'] - DZP_FC)
    assert abs(below_dz - 23.64) < 0.05, below_dz
    assert margin > 0, margin
    assert abs(margin - 17.07) < 0.05, margin


@pytest.mark.slow
def test_best_atom_set_rung_is_4_7_mha_above_ano_r1():
    """The best Slater rung set on free atoms (R5H, 3d set on the atom's correlated
    energy) is 4.7 mHa above ANO-R1.

    Wrong answers rejected: the comparison made against R5H without the atom-set 3d
    (16.6 mHa, the figure first reported to the PI and corrected), and a sign flip
    that would put the atom-set rung ahead of ANO-R1.
    """
    a = water_run(ANO_BASIS)
    gap = 1e3 * (R5H3D_FC - a['CCSD(T)_fc'])
    assert gap > 0, gap
    assert abs(gap - 4.66) < 0.05, gap


@pytest.mark.slow
def test_pcseg1_and_def2svp_land_above_cc_pvdz():
    """SZ1: pcseg-1 (-76.217263) and def2-SVP (-76.174699) are both above cc-pVDZ at
    fc-CCSD(T), def2-SVP the higher, each at 24 functions.

    Wrong answers rejected: either basis reported below cc-pVDZ, the order of the two
    reversed, a function count other than 24, and a drifted basis file or solver.
    """
    p, d = water_run(PCSEG_BASIS), water_run(DEF2_BASIS)
    assert p['nbf'] == 24 and d['nbf'] == 24, (p['nbf'], d['nbf'])
    assert abs(p['CCSD(T)_fc'] - PCSEG_FC) < 5e-6, p['CCSD(T)_fc']
    assert abs(d['CCSD(T)_fc'] - DEF2_FC) < 5e-6, d['CCSD(T)_fc']
    assert CCPVDZ_FC < p['CCSD(T)_fc'] < d['CCSD(T)_fc'], (p['CCSD(T)_fc'], d['CCSD(T)_fc'])


@pytest.mark.slow
def test_ano_r1_same_pauli_count_lower_1_norm_than_cc_pvdz():
    """SZ5 and the 1-norm observation: ANO-R1 fc has the 24-function count, 107,382
    Pauli terms on 46 qubits, and a non-identity 1-norm of 625.36 Ha, below cc-pVDZ's
    674.3 while more accurate; the Slater DZP's 1-norm is 12.9% above ANO-R1's.

    Wrong answers rejected: a Pauli count that depends on the radial basis, the identity
    counted into the 1-norm (moves it by over 100 Ha), the frozen core not applied in the
    qubit step (48 qubits), and a 1-norm that rises with accuracy across the two
    Gaussian families.
    """
    a = water_run(ANO_BASIS)
    assert a['qubits'] == 46, a['qubits']
    assert a['n_pauli'] == PAULI_24, a['n_pauli']
    assert abs(a['one_norm'] - ANO_LAMBDA) < 0.05, a['one_norm']
    assert a['one_norm'] < CCPVDZ_LAMBDA and a['CCSD(T)_fc'] < CCPVDZ_FC
    assert abs(100 * (DZP_LAMBDA / a['one_norm'] - 1) - 12.94) < 0.05
