"""Exponent rule for the Slater DZP water basis (memo debug/water_barrier_memo.md):
  O s,p  : primitives 1s,1s',2s,2s' / 2p,2p',2p'' on the exact radial average-of-
           configuration HF energy of O 1s2 2s2 2p4; contracted like cc-pVDZ to
           [atomic 1s, atomic 2s, diffuse s | atomic 2p, diffuse p]
  H s,s' : RHF of H2 at R = 1.4 bohr
  pol    : O 3d + H 2p on water CCSD; cc-pVDZ's d/p exponents re-optimised the same way.
(Superseded: 'atoms' = uncontracted 3s2p on the atom, whose optimum spends the second
 s function on the core -- zetas 10.03/8.33/2.32 -- leaving single-zeta valence s.)
Engine: tests/test_water_ccsd.py.  Results -> debug/data/water_sto_exponents.json.

Usage: python debug/water_sto_exponents.py oatom | h2 | pol_sto | pol_gto | pol_gto_fc
       python debug/water_sto_exponents.py final [ng] [basis-key] [ng_pol]
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('w', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_s)
_s.loader.exec_module(w)
GE = w.GE
OUT = os.path.join(ROOT, 'debug', 'data', 'water_sto_exponents.json')
NG = 12


def load():
    return json.load(open(OUT)) if os.path.exists(OUT) else {}


def save(new):
    d = load()            # re-read: other modes may have saved meanwhile
    d.update(new)
    json.dump(d, open(OUT, 'w'), indent=1)


def oatom_energy(logz):
    z = np.exp(logz)
    sp = [(n, zz) for n, zz in zip(w.O_S_N, z[:4])]
    pp = [(2, zz) for zz in z[4:]]
    try:
        return w.radial_atom(sp, pp)[0]
    except (RuntimeError, np.linalg.LinAlgError):
        return 0.0


def contracted_ex(d):
    sp = [(n, zz) for n, zz in zip(w.O_S_N, d['Os'])]
    E, Cs, cp = w.radial_atom(sp, [(2, zz) for zz in d['Op']])
    return dict(d, cs=Cs, cp=cp)


def o_atom(logz):
    z = np.exp(logz)
    fns, blocks = [], []
    O = np.zeros(3)
    for (n, l), zz in zip(((1, 0), (2, 0), (2, 0), (2, 1), (2, 1)), z):
        w.sto_shell(O, n, l, zz, NG, fns, blocks)
    S, h, g = GE.integral_set_md(fns, [(O, 8.0)])
    s = [i for i, f in enumerate(fns) if f.lmn == (0, 0, 0)]
    ps = [[i for i, f in enumerate(fns) if f.lmn == c] for c in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    E = w.atom_average_energy(S, h, g, [(s, [2, 2])] + [(p, [4 / 3]) for p in ps], [1, 2, 3], 4)
    print('  O', np.round(z, 5), '%.8f' % E); sys.stdout.flush()
    return E


def h2(logz):
    z = np.exp(logz)
    nuc = [(np.zeros(3), 1.0), (np.array([0, 0, 1.4]), 1.0)]
    fns, blocks = [], []
    for C, _ in nuc:
        for zz in z:
            w.sto_shell(C, 1, 0, zz, NG, fns, blocks)
    S, h, g = GE.integral_set_md(fns, nuc)
    return w.rhf(S, h, g, 2, w.vnn(nuc))[0]


def water_ccsd_sto(logz, ex):
    e = dict(ex, O3d=float(np.exp(logz[0])), H2p=float(np.exp(logz[1])))
    nuc = w.water_nuclei()
    fns, T = w.slater_dzp_contracted(nuc, contracted_ex(e), NG)
    S, h, g = w.ao_integrals(fns, nuc, T)
    r = w.run(S, h, g, 10, w.vnn(nuc), triples=False)
    print('  sto pol', np.round(np.exp(logz), 5), '%.8f' % r['CCSD']); sys.stdout.flush()
    return r['CCSD']


def gto_basis_with_pol(nuc, d_exp, p_exp):
    data = json.load(open(os.path.join(ROOT, 'debug', 'data', 'basis', 'cc-pvdz_H_O.json')))
    for Z, val, am in (('8', d_exp, 2), ('1', p_exp, 1)):
        for sh in data['elements'][Z]['electron_shells']:
            if sh['angular_momentum'][0] == am and len(sh['exponents']) == 1:
                sh['exponents'] = [str(val)]
    tmp = os.path.join(ROOT, 'debug', 'data', 'basis', '_ccpvdz_polopt_tmp.json')
    json.dump(data, open(tmp, 'w'))
    return w.gaussian_basis(tmp, nuc)


def water_ccsd_gto(logz, fc=0):
    nuc = w.water_nuclei()
    fns, T = gto_basis_with_pol(nuc, float(np.exp(logz[0])), float(np.exp(logz[1])))
    S, h, g = w.ao_integrals(fns, nuc, T)
    r = w.run(S, h, g, 10, w.vnn(nuc), triples=False, frozen_core=fc)
    print('  gto pol', 'fc' if fc else 'ae', np.round(np.exp(logz), 5), '%.8f' % r['CCSD']); sys.stdout.flush()
    return r['CCSD']


if __name__ == '__main__':
    cmd = sys.argv[1]
    d = load()
    t0 = time.time()
    new = {}
    if cmd == 'oatom':
        x = np.log([10.0, 7.0, 2.6, 1.6, 5.0, 2.4, 1.2])
        for rnd in range(4):     # restarted Nelder-Mead
            res = minimize(oatom_energy, x, method='Nelder-Mead',
                           options={'xatol': 1e-6, 'fatol': 1e-11, 'maxiter': 20000, 'maxfev': 20000})
            x = res.x
            print('round', rnd, '%.10f' % res.fun, np.round(np.exp(x), 5)); sys.stdout.flush()
        z = np.exp(x)
        new = {'Os': list(z[:4]), 'Op': list(z[4:]), 'E_O_atom_avg_radial': float(res.fun)}
    elif cmd == 'h2':
        r2 = minimize(h2, np.log([1.4, 0.9]), method='Nelder-Mead', options={'xatol': 1e-4, 'fatol': 1e-10})
        zh = np.exp(r2.x)
        print('H2 RHF  E = %.8f  zetas %s' % (r2.fun, np.round(zh, 5)))
        new = {'H1s': float(zh[0]), "H1s'": float(zh[1]), 'E_H2_RHF': float(r2.fun)}
    elif cmd == 'pol_sto':
        res = minimize(water_ccsd_sto, np.log([1.8, 1.0]), args=(d,), method='Nelder-Mead',
                       options={'xatol': 2e-3, 'fatol': 2e-6, 'maxiter': 60})
        new = {'O3d': float(np.exp(res.x[0])), 'H2p': float(np.exp(res.x[1])), 'E_CCSD_sto_polopt': float(res.fun)}
    elif cmd == 'pol_gto':
        res = minimize(water_ccsd_gto, np.log([1.185, 0.727]), method='Nelder-Mead',
                       options={'xatol': 2e-3, 'fatol': 2e-6, 'maxiter': 60})
        new = {'gto_d': float(np.exp(res.x[0])), 'gto_p': float(np.exp(res.x[1])), 'E_CCSD_gto_polopt': float(res.fun)}
    elif cmd == 'pol_gto_fc':
        res = minimize(water_ccsd_gto, np.log([1.185, 0.727]), args=(1,), method='Nelder-Mead',
                       options={'xatol': 2e-3, 'fatol': 2e-6, 'maxiter': 60})
        new = {'gto_d_fc': float(np.exp(res.x[0])), 'gto_p_fc': float(np.exp(res.x[1])),
               'E_CCSD_fc_gto_polopt': float(res.fun)}
    elif cmd == 'final':
        nuc = w.water_nuclei()
        ng = int(sys.argv[2]) if len(sys.argv) > 2 else NG
        ngp = int(sys.argv[4]) if len(sys.argv) > 4 else ng
        # fc keys: polarization exponents optimized on frozen-core CCSD (the comparison criterion)
        rows = {'sto_dzp_contracted_ng%d_pol%d' % (ng, ngp):
                    lambda: w.slater_dzp_contracted(nuc, contracted_ex(dict(d, O3d=d['O3d_fc'], H2p=d['H2p_fc'])), ng, ng_pol=ngp),
                'ccpvdz': lambda: w.gaussian_basis(os.path.join(ROOT, 'debug', 'data', 'basis', 'cc-pvdz_H_O.json'), nuc),
                'ccpvdz_polopt_fc': lambda: gto_basis_with_pol(nuc, d['gto_d_fc'], d['gto_p_fc'])}
        if len(sys.argv) > 3:
            rows = {k: v for k, v in rows.items() if k == sys.argv[3] or k.startswith(sys.argv[3] + '_')}
        for name, build in rows.items():
            fns, T = build()
            S, h, g = w.ao_integrals(fns, nuc, T)
            r = w.run(S, h, g, 10, w.vnn(nuc))
            rf = w.run(S, h, g, 10, w.vnn(nuc), frozen_core=1)
            row = {'basis': name, 'nbf': T.shape[1], **{k: float(v) for k, v in r.items()},
                   **{k + '_fc': float(v) for k, v in rf.items() if k != 'RHF'}}
            print(json.dumps(row)); sys.stdout.flush()
            with open(os.path.join(ROOT, 'debug', 'data', 'water_ccsd.jsonl'), 'a') as f:
                f.write(json.dumps(row) + '\n')
    if new:
        save(new)
    print(cmd, 'done in %.0fs' % (time.time() - t0), json.dumps(new))
