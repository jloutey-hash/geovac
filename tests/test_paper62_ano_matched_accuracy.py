"""Paper 62, matched accuracy in the ANO family (Sec. matched, Table ano): regression guard.

Record: docs/predictions/2026-11-11_ano_matched_accuracy.md (AM; AM1 missed, AM2-AM4 held).
Written as a separate step after the paper was restated (CLAUDE.md Sec. 9, guard-writing
rule) and fire-tested against the wrong answers named in each test's docstring.

What is recomputed here and what is read.  The rungs are ANO-R2 kept to its first
contracted functions per shell.  The truncation is RECOMPUTED here from the ANO-R2 file
stored beside this test (Basis Set Exchange download of 2026-11-11, exponents and
coefficients only) and checked against the stored ANO-R1 file.  The 35-function rung's
water run is RECOMPUTED end to end: RHF, fc-CCSD(T), and the Jordan-Wigner cost over the
RHF canonical orbitals with the oxygen 1s frozen.  The Slater DZP target and the
24-function Pauli count are READ as registered literals (p62_dzp_fc, p62_dzp_lambda_fc,
h2o_ccpvdz_pauli; C21 pins them to the paper).  The 41-function rung is not recomputed
(about ten minutes); its row is pinned by C21 only.

Cost: the 35-function run is about six minutes; it is @pytest.mark.slow.  The
truncation checks are fast.
"""
import importlib.util
import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_am', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

R2_FILE = os.path.join(HERE, 'data', 'paper62_ano_r2_H_O_Ne_2026_11_11.json')
R1_FILE = os.path.join(HERE, 'data', 'paper62_ano_r1_H_O_Ne_2026_11_11.json')

# contracted functions kept per angular momentum (s, p, d, f), per element; the record's rungs
RUNGS = {
    'A24': {'1': (2, 1), '8': (3, 2, 1)},
    'A35': {'1': (2, 2), '8': (3, 2, 2)},
    'A41': {'1': (3, 2), '8': (4, 3, 2)},
}

# READ (registered literals, C21-pinned)
DZP_FC = -76.281739          # p62_dzp_fc
DZP_LAMBDA = 706.3           # p62_dzp_lambda_fc
PAULI_24 = 107382            # h2o_ccpvdz_pauli, every 24-function basis, fc, Q = 46
AUG_PAULI = 982785           # p62_aug_pauli_fc, the cc-pV ladder's 9.2x rung

# RECOMPUTED here and asserted
A35_FC = -76.287896          # p62_a35_fc
A35_RHF = -76.053557
A35_PAULI = 505467           # p62_a35_pauli_fc
A35_LAMBDA = 1632.93         # p62_a35_lambda_fc (registered at display precision 1632.9)


def truncate(data: dict, keep: dict) -> dict:
    """ANO-R2 kept to its first n contracted functions per shell (the driver's rule)."""
    out = {'elements': {}}
    for z, n_by_l in keep.items():
        shells = []
        for sh in data['elements'][z]['electron_shells']:
            l = sh['angular_momentum'][0]
            n = n_by_l[l] if l < len(n_by_l) else 0
            if n:
                assert n <= len(sh['coefficients']), (z, l, n)
                shells.append({**sh, 'coefficients': sh['coefficients'][:n]})
        out['elements'][z] = {'electron_shells': shells}
    return out


def _floats(xs):
    return np.array([float(x) for x in xs])


def _write(tmp_path, data, name):
    p = tmp_path / name
    p.write_text(json.dumps(data))
    return str(p)


def test_ano_r1_is_the_24_function_truncation_of_ano_r2():
    """The ladder's premise: ANO-R1 equals ANO-R2 kept to [3s2p1d]/[2s1p], the same
    primitives and the same leading contraction columns, shell by shell.

    Wrong answers rejected: a rung built from the wrong columns (the last ones, or a
    re-contraction), a primitive set that differs between R1 and R2, and a shell count
    other than ANO-R1's.
    """
    r2 = json.load(open(R2_FILE))
    r1 = json.load(open(R1_FILE))
    a24 = truncate(r2, RUNGS['A24'])
    for z in ('1', '8'):
        got = a24['elements'][z]['electron_shells']
        ref = r1['elements'][z]['electron_shells']
        assert len(got) == len(ref), (z, len(got), len(ref))
        for g, r in zip(got, ref):
            assert g['angular_momentum'] == r['angular_momentum']
            assert np.allclose(_floats(g['exponents']), _floats(r['exponents']), rtol=1e-12, atol=0)
            assert len(g['coefficients']) == len(r['coefficients'])
            for cg, cr in zip(g['coefficients'], r['coefficients']):
                assert np.allclose(_floats(cg), _floats(cr), rtol=1e-10, atol=1e-14)


def test_rung_function_counts_on_water(tmp_path):
    """A24, A35 and A41 have 24, 35 and 41 functions on water, the sizes of cc-pVDZ,
    cc-pVDZ+TZpol and aug-cc-pVDZ.

    Wrong answers rejected: Cartesian instead of spherical d functions, a truncation
    that keeps a shell the rung omits (an H d shell, an O f shell), and miscounted
    general contractions.
    """
    r2 = json.load(open(R2_FILE))
    nuclei = w.water_nuclei()
    for rung, n in (('A24', 24), ('A35', 35), ('A41', 41)):
        path = _write(tmp_path, truncate(r2, RUNGS[rung]), rung + '.json')
        _, T = w.gaussian_basis(path, nuclei)
        assert T.shape[1] == n, (rung, T.shape[1])


_CACHE = {}


def _a35_run(tmp_path_factory):
    if 'a35' in _CACHE:
        return _CACHE['a35']
    r2 = json.load(open(R2_FILE))
    path = _write(tmp_path_factory.mktemp('am'), truncate(r2, RUNGS['A35']), 'A35.json')
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
    _CACHE['a35'] = out
    return out


@pytest.mark.slow
def test_a35_reaches_the_slater_dzp(tmp_path_factory):
    """AM1 MISSED: the 35-function ANO rung's water fc-CCSD(T) is -76.287896, 6.2 mHa
    BELOW the Slater DZP (-76.281739); RHF -76.053557.

    Wrong answers rejected: the registered prediction (A35 short of the DZP, positive
    gap), a drifted basis file or solver, and the all-electron energy reported as
    frozen-core.
    """
    a = _a35_run(tmp_path_factory)
    assert a['nbf'] == 35, a['nbf']
    assert abs(a['RHF'] - A35_RHF) < 5e-6, a['RHF']
    assert abs(a['CCSD(T)_fc'] - A35_FC) < 5e-6, a['CCSD(T)_fc']
    below = 1e3 * (DZP_FC - a['CCSD(T)_fc'])
    assert below > 0, below
    assert abs(below - 6.16) < 0.05, below


@pytest.mark.slow
def test_matched_accuracy_ceiling_is_4_7x_terms_2_3x_1_norm(tmp_path_factory):
    """The paper's restated headline: the cheapest Gaussian basis measured that reaches
    the Slater DZP costs 505,467 Pauli terms on 68 qubits (4.71x the 24-function count)
    and a non-identity 1-norm of 1632.9 Ha (2.31x the DZP's).

    Wrong answers rejected: the retired "at least 9.2x" (aug-cc-pVDZ's ratio quoted as
    the ceiling), the identity counted into the 1-norm, the frozen core not applied in
    the qubit step (70 qubits, 570,822 terms), and a ratio taken against the wrong
    baseline count.
    """
    a = _a35_run(tmp_path_factory)
    assert a['qubits'] == 68, a['qubits']
    assert a['n_pauli'] == A35_PAULI, a['n_pauli']
    ratio = a['n_pauli'] / PAULI_24
    assert abs(ratio - 4.707) < 0.005, ratio
    assert ratio < AUG_PAULI / PAULI_24
    assert abs(a['one_norm'] - A35_LAMBDA) < 0.05, a['one_norm']
    assert abs(a['one_norm'] / DZP_LAMBDA - 2.312) < 0.005
