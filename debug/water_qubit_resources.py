"""Step 3 of debug/water_barrier_memo.md: qubit resources (Pauli count, 1-norm) of the
water Hamiltonian in three bases -- Slater DZP (24 fns), cc-pVDZ (24), cc-pVTZ (58).

  integrals <key>  : AO integrals -> RHF -> MO integrals saved to debug/data/water_mo_<key>.npz
                     keys: sto_dzp | ccpvdz | ccpvtz
  count            : JW Pauli count + 1-norm for every saved .npz, all-electron and frozen-core;
                     rows -> debug/data/water_qubit_resources.json

Counter: tests/test_water_ccsd.py::jw_pauli_stats (validated there against openfermion).
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('w', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_s)
_s.loader.exec_module(w)
_e = importlib.util.spec_from_file_location('e', os.path.join(ROOT, 'debug', 'water_sto_exponents.py'))
e = importlib.util.module_from_spec(_e)
_e.loader.exec_module(e)
DATA = os.path.join(ROOT, 'debug', 'data')


def build(key, nuc):
    if key == 'sto_dzp':
        d = e.load()
        return w.slater_dzp_contracted(nuc, e.contracted_ex(dict(d, O3d=d['O3d_fc'], H2p=d['H2p_fc'])), 12, ng_pol=8)
    return w.gaussian_basis(os.path.join(DATA, 'basis', key[:2] + '-' + key[2:] + '_H_O.json'), nuc)


if __name__ == '__main__':
    t0 = time.time()
    if sys.argv[1] == 'integrals':
        key = sys.argv[2]
        nuc = w.water_nuclei()
        fns, T = build(key, nuc)
        S, h, g = w.ao_integrals(fns, nuc, T)
        enuc = w.vnn(nuc)
        E, eps, C = w.rhf(S, h, g, 10, enuc)
        hmo = C.T @ h @ C
        gmo = np.einsum('pqrs,pi,qj,rk,sl->ijkl', g, C, C, C, C, optimize=True)
        np.savez(os.path.join(DATA, 'water_mo_%s.npz' % key), h=hmo, g=gmo, enuc=enuc, E_rhf=E, eps=eps)
        print(key, 'nbf', C.shape[1], 'RHF %.8f' % E, 'done in %.0fs' % (time.time() - t0))
    elif sys.argv[1] == 'count':
        rows = []
        for key in ('sto_dzp', 'ccpvdz', 'ccpvtz'):
            p = os.path.join(DATA, 'water_mo_%s.npz' % key)
            if not os.path.exists(p):
                print('missing', key); continue
            z = np.load(p)
            for fc in (0, 1):
                h, g, c0 = w.freeze_core(z['h'], z['g'], fc)
                st = w.jw_pauli_stats(h, g, c0 + float(z['enuc']))
                row = {'basis': key, 'frozen_core': fc, 'qubits': 2 * h.shape[0], **st}
                print(json.dumps(row)); sys.stdout.flush()
                rows.append(row)
        with open(os.path.join(DATA, 'water_qubit_resources.json'), 'w') as f:
            json.dump(rows, f, indent=1)
        print('count done in %.0fs' % (time.time() - t0))
