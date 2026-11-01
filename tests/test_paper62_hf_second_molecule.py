"""Paper 62, Sec. hf (Table tab:hf): regression guard for the second molecule, hydrogen fluoride.

Record: docs/predictions/2026-10-31_hf_slater_second_molecule.md (HS, GO).  Written as a separate step after the
result was homed (CLAUDE.md Sec. 9, guard-writing rule) and fire-tested against the wrong answers named below.

Recomputed here: the free F atom's radial average-of-configuration HF energy at the stored primitives, and that they
are a stationary minimum (fast); cc-pVDZ (as published and with its polarization re-optimized on HF) and
aug-cc-pVDZ on HF, end to end: RHF, fc-CCSD(T), and the fc Jordan-Wigner cost (slow, about 3 min).
Read, not recomputed: the Slater DZP energies (23 min at 12 Gaussians per function), cc-pVDZ+TZpol and cc-pVTZ,
as the registered literals (registry `p62_hf_*`, C21-pinned to the paper).  Basis files and fitted exponents are
stored beside this test (tests/data/paper62_hf_*), so nothing is read from debug/.
"""
import importlib.util
import json
import os

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location('water_ccsd_for_p62_hf', os.path.join(HERE, 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(w)

DATA = json.load(open(os.path.join(HERE, 'data', 'paper62_hf_2026_10_31.json')))
DZ = os.path.join(HERE, 'data', 'paper62_hf_cc_pvdz_H_F_2026_10_31.json')
AUG = os.path.join(HERE, 'data', 'paper62_hf_aug_cc_pvdz_H_F_2026_10_31.json')
S_N = (1, 1, 2, 2)
NUC = [(np.zeros(3), 9.0), (np.array([0.0, 0.0, DATA['R_HF']]), 1.0)]

# READ (registered literals; measured in their own runs)
DZP_FC, DZP_RHF, DZP_LAMBDA = -100.287713, -100.060707, 432.6
TZPOL_FC, CCPVTZ_FC = -100.255091, -100.338356
# RECOMPUTED here and asserted
CCPVDZ_FC, CCPVDZ_RHF, CCPVDZ_LAMBDA = -100.228146, -100.019419, 433.6
POLOPT_FC, POLOPT_RHF = -100.228217, -100.019359
AUG_FC, AUG_LAMBDA = -100.263583, 1169.9
PAULI_19, PAULI_AUG = 63535, 591674

_CACHE = {}


def basis_with_pol(path, d_exp, p_exp):
    data = json.load(open(path))
    for Z, val, am in (('9', d_exp, 2), ('1', p_exp, 1)):
        for sh in data['elements'][Z]['electron_shells']:
            if sh['angular_momentum'][0] == am and len(sh['exponents']) == 1:
                sh['exponents'] = [str(val)]
    tmp = os.path.join(HERE, 'data', '_paper62_hf_polopt_tmp.json')
    json.dump(data, open(tmp, 'w'))
    try:
        return w.gaussian_basis(tmp, NUC)
    finally:
        os.remove(tmp)


def run(name):
    if name in _CACHE:
        return _CACHE[name]
    if name == 'polopt':
        fns, T = basis_with_pol(DZ, *DATA['F_gto_pol'])
    else:
        fns, T = w.gaussian_basis({'ccpvdz': DZ, 'aug': AUG}[name], NUC)
    S, h, g = w.ao_integrals(fns, NUC, T)
    enuc = w.vnn(NUC)
    r = w.run(S, h, g, 10, enuc, frozen_core=1)
    _, eps, C = w.rhf(S, h, g, 10, enuc)
    hm = C.T @ h @ C
    gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
    hh, gg, ecore = w.freeze_core(hm, gm, 1)
    st = w.jw_pauli_stats(hh, gg, enuc + ecore)
    out = {'nbf': T.shape[1], 'RHF': float(r['RHF']), 'fc': float(r['CCSD(T)']), 'qubits': 2 * hh.shape[0],
           'n_pauli': int(st['n_pauli']), 'one_norm': float(st['one_norm'])}
    _CACHE[name] = out
    return out


def atom_energy(z):
    return w.radial_atom([(n, x) for n, x in zip(S_N, z[:4])], [(2, x) for x in z[4:]], Z=9.0, q_p=5)[0]


def test_fluorine_primitives_are_the_free_atom_optimum():
    """The stored F primitives give the registered free-atom energy -99.408295 and are a stationary minimum of the
    radial average-of-configuration HF objective (Z = 9, five p electrons): the design rule, not a molecular fit.

    Wrong answers rejected: primitives optimized for another occupation (oxygen's four p electrons), exponents
    tuned on the molecule, and a stored energy that does not match its exponents.
    """
    z = np.array(DATA['F_atom']['s'] + DATA['F_atom']['p'])
    E = atom_energy(z)
    assert abs(E - (-99.408295)) < 1e-6, E
    assert abs(E - DATA['F_atom']['E_radial_avg']) < 1e-9
    for k in range(len(z)):
        hstep = 1e-3 * z[k]
        up, dn = z.copy(), z.copy()
        up[k] += hstep
        dn[k] -= hstep
        Eu, Ed = atom_energy(up), atom_energy(dn)
        assert abs((Eu - Ed) / (2 * hstep)) < 5e-5, (k, (Eu - Ed) / (2 * hstep))
        assert Eu > E - 1e-9 and Ed > E - 1e-9, k


@pytest.mark.slow
def test_gaussian_rungs_recompute_to_the_registered_values():
    """cc-pVDZ (published and re-optimized polarization) and aug-cc-pVDZ on HF, recomputed end to end: energies,
    RHF, fc Pauli count, non-identity 1-norm and qubits match the registered values.

    Wrong answers rejected: a basis file that moved, a solver that drifted, the all-electron energy reported as
    frozen-core, the identity counted into the 1-norm, and the frozen core dropped from the qubit step.
    """
    d, p, a = run('ccpvdz'), run('polopt'), run('aug')
    assert d['nbf'] == 19 and a['nbf'] == 32
    assert abs(d['fc'] - CCPVDZ_FC) < 5e-6 and abs(d['RHF'] - CCPVDZ_RHF) < 5e-6, d
    assert abs(p['fc'] - POLOPT_FC) < 5e-6 and abs(p['RHF'] - POLOPT_RHF) < 5e-6, p
    assert abs(a['fc'] - AUG_FC) < 5e-6, a
    assert d['qubits'] == 36 and a['qubits'] == 62
    assert d['n_pauli'] == PAULI_19 and a['n_pauli'] == PAULI_AUG, (d['n_pauli'], a['n_pauli'])
    assert abs(d['one_norm'] - CCPVDZ_LAMBDA) < 0.05 and abs(a['one_norm'] - AUG_LAMBDA) < 0.05


@pytest.mark.slow
def test_hs_predictions_follow_from_the_recomputed_rungs():
    """HS1-HS4 derived from the recomputed Gaussian rungs and the registered Slater DZP: the equal-size gap is 59.50
    mHa (> 20, inside [39.6, 81.4]), its RHF share 0.695 (> 0.5); aug-cc-pVDZ (and the TZpol rung) lie above the Slater
    DZP and cc-pVTZ below; aug-cc-pVDZ costs 9.31x the terms and 2.70x the 1-norm.

    Wrong answers rejected: a gap of the wrong sign or below the bound, a share read from all-electron energies, the
    gap measured against the published rather than the re-optimized cc-pVDZ, and aug-cc-pVDZ at or below the DZP.
    """
    p, a = run('polopt'), run('aug')
    gap = 1e3 * (p['fc'] - DZP_FC)
    assert 20 < gap and 39.6 <= gap <= 81.4 and abs(gap - 59.50) < 0.01, gap
    share = (p['RHF'] - DZP_RHF) / (p['fc'] - DZP_FC)
    assert share > 0.5 and abs(share - 0.695) < 0.001, share
    assert a['fc'] > DZP_FC and TZPOL_FC > DZP_FC and CCPVTZ_FC < DZP_FC
    assert abs(a['n_pauli'] / PAULI_19 - 9.31) < 0.01
    assert abs(a['one_norm'] / DZP_LAMBDA - 2.70) < 0.01
