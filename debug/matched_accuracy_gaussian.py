"""Matched-accuracy Gaussian comparison for water (record MA, 2026-10-31).

Which Gaussian basis first reaches the Slater DZP's fc-CCSD(T) energy (-76.281739 Ha at
24 functions, CHANGELOG v5.20.11), and what does it cost in qubits, Pauli terms and 1-norm?
Ladder by size: cc-pVDZ (24) < cc-pVDZ+TZpol (35) < aug-cc-pVDZ (41) < cc-pVTZ (58).
cc-pVDZ+TZpol = cc-pVDZ's s and p shells with cc-pVTZ's polarization shells (O 2d, H 2p).

    python debug/matched_accuracy_gaussian.py prep          # build the two basis files
    python debug/matched_accuracy_gaussian.py ne <basis>    # control: Ne atom, non-scored
    python debug/matched_accuracy_gaussian.py water <basis> # scored run (after registration)

Basis names: ccpvdz, ccpvdz_tzpol, augccpvdz, ccpvtz.  Downloaded files (untrusted) live in
debug/lit_scan/dl_bse_augccpvdz/; `prep` re-serializes them into debug/data/basis/.
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, 'data')
BASIS = os.path.join(DATA, 'basis')
DL = os.path.join(HERE, 'lit_scan', 'dl_bse_augccpvdz')
spec = importlib.util.spec_from_file_location('w', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)

FILES = {
    'ccpvdz': 'cc-pvdz_H_O_Ne.json',
    'ccpvdz_tzpol': 'ccpvdz_tzpol_H_O_Ne.json',
    'augccpvdz': 'aug-cc-pvdz_H_O_Ne.json',
    'ccpvtz': 'cc-pvtz_H_O_Ne.json',
}
SLATER_DZP_FC = -76.281739


def _shells(data, z):
    return data['elements'][str(z)]['electron_shells']


def _clean(data):
    """Keep only what gaussian_basis reads: exponents and coefficients per shell."""
    out = {'elements': {}}
    for z, el in data['elements'].items():
        out['elements'][z] = {'electron_shells': [
            {'angular_momentum': [int(a) for a in sh['angular_momentum']],
             'exponents': [str(float(x)) for x in sh['exponents']],
             'coefficients': [[str(float(c)) for c in col] for col in sh['coefficients']]}
            for sh in el['electron_shells']]}
    return out


def prep():
    aug = _clean(json.load(open(os.path.join(DL, 'aug-cc-pvdz_H_O_Ne.json'))))
    dz_ho = json.load(open(os.path.join(BASIS, 'cc-pvdz_H_O.json')))
    dz_ne = json.load(open(os.path.join(DL, 'cc-pvdz_Ne.json')))
    tz_ho = json.load(open(os.path.join(BASIS, 'cc-pvtz_H_O.json')))
    tz_ne = json.load(open(os.path.join(DL, 'cc-pvtz_Ne.json')))
    dz = _clean({'elements': {**dz_ho['elements'], **dz_ne['elements']}})
    tz = _clean({'elements': {**tz_ho['elements'], **tz_ne['elements']}})
    merged = {'elements': {}}
    for z in ('1', '8', '10'):
        lmax_sp = 0 if z == '1' else 1
        sp = [sh for sh in _shells(dz, z) if sh['angular_momentum'][0] <= lmax_sp]
        pol = [sh for sh in _shells(tz, z) if sh['angular_momentum'][0] == lmax_sp + 1]
        merged['elements'][z] = {'electron_shells': sp + pol}
    for name, obj in (('aug-cc-pvdz_H_O_Ne.json', aug), ('cc-pvdz_H_O_Ne.json', dz),
                      ('cc-pvtz_H_O_Ne.json', tz), ('ccpvdz_tzpol_H_O_Ne.json', merged)):
        with open(os.path.join(BASIS, name), 'w') as f:
            json.dump(obj, f, indent=1)
        print('wrote', name)


def count(name, nuclei):
    fns, T = w.gaussian_basis(os.path.join(BASIS, FILES[name]), nuclei)
    return T.shape[1]


def ne(name):
    nuclei = [(np.zeros(3), 10.0)]
    t0 = time.time()
    fns, T = w.gaussian_basis(os.path.join(BASIS, FILES[name]), nuclei)
    S, h, g = w.ao_integrals(fns, nuclei, T)
    t1 = time.time()
    r = w.run(S, h, g, 10, 0.0)
    rf = w.run(S, h, g, 10, 0.0, frozen_core=1)
    row = {'case': 'Ne', 'basis': name, 'nbf': T.shape[1], **{k: float(v) for k, v in r.items()},
           **{k + '_fc': float(v) for k, v in rf.items() if k != 'RHF'},
           't_integrals': round(t1 - t0, 1), 't_solve': round(time.time() - t1, 1)}
    print(json.dumps(row), flush=True)
    with open(os.path.join(DATA, 'matched_accuracy_controls.jsonl'), 'a') as f:
        f.write(json.dumps(row) + '\n')
    return row


def water(name):
    nuclei = w.water_nuclei()
    t0 = time.time()
    fns, T = w.gaussian_basis(os.path.join(BASIS, FILES[name]), nuclei)
    S, h, g = w.ao_integrals(fns, nuclei, T)
    enuc = w.vnn(nuclei)
    t1 = time.time()
    r = w.run(S, h, g, 10, enuc)
    rf = w.run(S, h, g, 10, enuc, frozen_core=1)
    t2 = time.time()
    # qubit resources over the RHF canonical orbitals, as in v5.20.11
    # (debug/water_qubit_resources.py): chemist (pq|rs) in the MO basis, identity
    # counted in N_Pauli and excluded from the 1-norm by jw_pauli_stats itself.
    _, eps, C = w.rhf(S, h, g, 10, enuc)
    hm = C.T @ h @ C
    gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
    res = {}
    for fc in (0, 1):
        hh, gg, ecore = w.freeze_core(hm, gm, fc)
        st = w.jw_pauli_stats(hh, gg, enuc + ecore)
        res['fc%d' % fc] = {'qubits': 2 * hh.shape[0],
                            **{k: (float(v) if not isinstance(v, int) else v) for k, v in st.items()}}
    row = {'case': 'water', 'basis': name, 'nbf': T.shape[1],
           **{k: float(v) for k, v in r.items()},
           **{k + '_fc': float(v) for k, v in rf.items() if k != 'RHF'},
           'qubit': res, 't_integrals': round(t1 - t0, 1), 't_solve': round(t2 - t1, 1)}
    print(json.dumps(row), flush=True)
    with open(os.path.join(DATA, 'matched_accuracy.jsonl'), 'a') as f:
        f.write(json.dumps(row) + '\n')
    return row


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'prep':
        prep()
    elif cmd == 'ne':
        ne(sys.argv[2])
    elif cmd == 'water':
        water(sys.argv[2])
    elif cmd == 'count':
        print({n: count(n, w.water_nuclei()) for n in FILES})
