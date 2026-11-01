"""Thread 3 (2026-10-31): does Paper 62's Slater-DZP result carry to a second molecule?

Hydrogen fluoride at r_e = 1.7325 bohr (0.9168 A), 19 functions as cc-pVDZ.  The Slater DZP is
built by the rule Paper 62 used for water (debug/water_barrier_memo.md, CHANGELOG v5.20.11):
  heavy atom s,p : primitives 1s,1s',2s,2s' / 2p,2p',2p'' optimized on the exact radial
                   average-of-configuration HF energy of the free atom; contracted like cc-pVDZ to
                   [atomic 1s, atomic 2s, most diffuse s | atomic 2p, most diffuse p]
  H s,s'         : H2 RHF at R = 1.4 (1.212, 0.627; unchanged from water, H-only rule)
  polarization   : heavy-atom 3d + H 2p optimized on the molecule's frozen-core CCSD, and
                   cc-pVDZ's d/p re-optimized the same way (both bases tuned equally).
Neon (isoelectronic closed shell) is the non-scored control: same rule, no hydrogens.

    python debug/hf_slater_matched.py prep                      # Gaussian basis files for H/F
    python debug/hf_slater_matched.py atom F|Ne                 # primitives on the radial atom
    python debug/hf_slater_matched.py pol_sto F|Ne              # Slater 3d (+ H 2p) on fc-CCSD
    python debug/hf_slater_matched.py pol_gto F|Ne              # cc-pVDZ d (+ p) on fc-CCSD
    python debug/hf_slater_matched.py final F|Ne <basis>        # RHF, CCSD(T) fc/ae, JW cost
basis: sto_dzp | ccpvdz | ccpvdz_polopt | ccpvdz_tzpol | augccpvdz | ccpvtz
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
_s = importlib.util.spec_from_file_location('w', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_s)
_s.loader.exec_module(w)
GE = w.GE
DATA = os.path.join(HERE, 'data')
BASIS = os.path.join(DATA, 'basis')
DL = os.path.join(HERE, 'lit_scan', 'dl_bse_fluorine')
EXP = os.path.join(DATA, 'hf_slater_exponents.json')
OUT = os.path.join(DATA, 'hf_slater_matched.jsonl')
R_HF = 1.7325
H_EX = {'H1s': 1.212, "H1s'": 0.627}          # H2 RHF, R = 1.4 (water memo)
SYSTEMS = {'F': dict(Z=9.0, q_p=5, nelec=10), 'Ne': dict(Z=10.0, q_p=6, nelec=10)}
NG, NG_OPT, NGP_OPT = 12, 8, 6
S_N = (1, 1, 2, 2)


def nuclei(sysname):
    X = (np.zeros(3), SYSTEMS[sysname]['Z'])
    return [X, (np.array([0.0, 0.0, R_HF]), 1.0)] if sysname == 'F' else [X]


def load():
    return json.load(open(EXP)) if os.path.exists(EXP) else {}


def save(new):
    d = load()
    d.update(new)
    json.dump(d, open(EXP, 'w'), indent=1)


def write(row):
    print(json.dumps(row), flush=True)
    with open(OUT, 'a') as f:
        f.write(json.dumps(row) + '\n')


# ------------------------------------------------------------------ Gaussian basis files
def _clean(data):
    return {'elements': {z: {'electron_shells': [
        {'angular_momentum': [int(a) for a in sh['angular_momentum']],
         'exponents': [str(float(x)) for x in sh['exponents']],
         'coefficients': [[str(float(c)) for c in col] for col in sh['coefficients']]}
        for sh in el['electron_shells']]} for z, el in data['elements'].items()}}


def prep():
    files = {b: _clean(json.load(open(os.path.join(DL, b + '_H_F.json'))))
             for b in ('cc-pvdz', 'aug-cc-pvdz', 'cc-pvtz')}
    dz, tz = files['cc-pvdz'], files['cc-pvtz']
    merged = {'elements': {}}
    for z in ('1', '9'):
        lsp = 0 if z == '1' else 1
        sp = [sh for sh in dz['elements'][z]['electron_shells'] if sh['angular_momentum'][0] <= lsp]
        pol = [sh for sh in tz['elements'][z]['electron_shells'] if sh['angular_momentum'][0] == lsp + 1]
        merged['elements'][z] = {'electron_shells': sp + pol}
    files['ccpvdz_tzpol'] = merged
    for name, obj in files.items():
        p = os.path.join(BASIS, name + '_H_F.json')
        json.dump(obj, open(p, 'w'), indent=1)
        print('wrote', p)


def ne_gaussian_file(name):
    """Neon files already on disk from record MA (H/O/Ne)."""
    return os.path.join(BASIS, {'ccpvdz': 'cc-pvdz_H_O_Ne.json', 'augccpvdz': 'aug-cc-pvdz_H_O_Ne.json',
                                'ccpvtz': 'cc-pvtz_H_O_Ne.json', 'ccpvdz_tzpol': 'ccpvdz_tzpol_H_O_Ne.json'}[name])


def f_gaussian_file(name):
    return os.path.join(BASIS, {'ccpvdz': 'cc-pvdz', 'augccpvdz': 'aug-cc-pvdz', 'ccpvtz': 'cc-pvtz',
                                'ccpvdz_tzpol': 'ccpvdz_tzpol'}[name] + '_H_F.json')


def gaussian(sysname, name, pol=None):
    path = f_gaussian_file(name) if sysname == 'F' else ne_gaussian_file(name)
    if pol is not None:
        data = json.load(open(path))
        heavy = '9' if sysname == 'F' else '10'
        for Z, val, am in ((heavy, pol[0], 2), ('1', pol[1] if len(pol) > 1 else None, 1)):
            if Z not in data['elements'] or val is None:
                continue
            for sh in data['elements'][Z]['electron_shells']:
                if sh['angular_momentum'][0] == am and len(sh['exponents']) == 1:
                    sh['exponents'] = [str(val)]
        path = os.path.join(BASIS, '_hf_polopt_tmp.json')
        json.dump(data, open(path, 'w'))
    return w.gaussian_basis(path, nuclei(sysname))


# ------------------------------------------------------------------ Slater DZP
def atom_energy(logz, Z, q_p):
    z = np.exp(logz)
    try:
        return w.radial_atom([(n, zz) for n, zz in zip(S_N, z[:4])], [(2, zz) for zz in z[4:]], Z=Z, q_p=q_p)[0]
    except (RuntimeError, np.linalg.LinAlgError):
        return 0.0


def slater_dzp(sysname, ex, n_gauss, ng_pol, d_exp, p_exp=None):
    sysd = SYSTEMS[sysname]
    nuc = nuclei(sysname)
    X = nuc[0][0]
    Es, Cs, cp = w.radial_atom([(n, z) for n, z in zip(S_N, ex['s'])], [(2, z) for z in ex['p']],
                               Z=sysd['Z'], q_p=sysd['q_p'])
    fns, blocks = [], []
    for n, z in zip(S_N, ex['s']):
        w.sto_shell(X, n, 0, z, n_gauss, fns, blocks)
    for z in ex['p']:
        w.sto_shell(X, 2, 1, z, n_gauss, fns, blocks)
    nX = len(fns)
    w.sto_shell(X, 3, 2, d_exp, ng_pol, fns, blocks)
    for C, Zc in nuc[1:]:
        w.sto_shell(C, 1, 0, H_EX['H1s'], n_gauss, fns, blocks)
        w.sto_shell(C, 1, 0, H_EX["H1s'"], n_gauss, fns, blocks)
        w.sto_shell(C, 2, 1, p_exp, ng_pol, fns, blocks)
    cols = []
    s = [i for i in range(nX) if fns[i].lmn == (0, 0, 0)]
    for k in range(2):
        v = np.zeros(len(fns)); v[s] = Cs[:, k]; cols.append(v)
    v = np.zeros(len(fns)); v[s[int(np.argmin(ex['s']))]] = 1.0; cols.append(v)
    jd = int(np.argmin(ex['p']))
    for c in w.CART[1]:
        p = [i for i in range(nX) if fns[i].lmn == c]
        v = np.zeros(len(fns)); v[p] = cp; cols.append(v)
        v = np.zeros(len(fns)); v[p[jd]] = 1.0; cols.append(v)
    rest = w.spherical_map(fns, [b for b in blocks if b[1] >= nX])
    return fns, np.hstack([np.array(cols).T, rest])


def fc_ccsd(sysname, fns, T):
    nuc = nuclei(sysname)
    S, h, g = w.ao_integrals(fns, nuc, T)
    return w.run(S, h, g, SYSTEMS[sysname]['nelec'], w.vnn(nuc), triples=False, frozen_core=1)['CCSD']


def pol_objective_sto(logz, sysname, ex):
    z = np.exp(logz)
    fns, T = slater_dzp(sysname, ex, NG_OPT, NGP_OPT, z[0], z[1] if len(z) > 1 else None)
    E = fc_ccsd(sysname, fns, T)
    print('  sto pol', np.round(z, 5), '%.8f' % E, flush=True)
    return E


def pol_objective_gto(logz, sysname):
    z = np.exp(logz)
    fns, T = gaussian(sysname, 'ccpvdz', pol=list(z))
    E = fc_ccsd(sysname, fns, T)
    print('  gto pol', np.round(z, 5), '%.8f' % E, flush=True)
    return E


def final(sysname, name):
    d = load()
    ex = d[sysname + '_atom']
    nuc = nuclei(sysname)
    t0 = time.time()
    if name == 'sto_dzp':
        pol = d[sysname + '_sto_pol']
        fns, T = slater_dzp(sysname, ex, NG, NG, pol[0], pol[1] if len(pol) > 1 else None)
    elif name == 'ccpvdz_polopt':
        fns, T = gaussian(sysname, 'ccpvdz', pol=d[sysname + '_gto_pol'])
    else:
        fns, T = gaussian(sysname, name)
    S, h, g = w.ao_integrals(fns, nuc, T)
    enuc = w.vnn(nuc)
    ne = SYSTEMS[sysname]['nelec']
    r = w.run(S, h, g, ne, enuc)
    rf = w.run(S, h, g, ne, enuc, frozen_core=1)
    row = {'system': 'HF' if sysname == 'F' else 'Ne', 'basis': name, 'nbf': T.shape[1],
           **{k: float(v) for k, v in r.items()}, **{k + '_fc': float(v) for k, v in rf.items() if k != 'RHF'}}
    if sysname == 'F':
        _, eps, C = w.rhf(S, h, g, ne, enuc)
        hm = C.T @ h @ C
        gm = np.einsum('pqrs,pa,qb,rc,sd->abcd', g, C, C, C, C, optimize=True)
        hh, gg, ecore = w.freeze_core(hm, gm, 1)
        st = w.jw_pauli_stats(hh, gg, enuc + ecore)
        row['qubit_fc'] = {'qubits': 2 * hh.shape[0], 'n_pauli': int(st['n_pauli']),
                           'one_norm': float(st['one_norm'])}
    row['t'] = round(time.time() - t0, 1)
    write(row)


if __name__ == '__main__':
    cmd = sys.argv[1]
    t0 = time.time()
    if cmd == 'prep':
        prep()
    elif cmd == 'atom':
        sysname = sys.argv[2]
        Z, q_p = SYSTEMS[sysname]['Z'], SYSTEMS[sysname]['q_p']
        x = np.log(np.array([10.0, 7.0, 2.6, 1.6, 5.0, 2.4, 1.2]) * Z / 8.0)
        for rnd in range(4):
            res = minimize(atom_energy, x, args=(Z, q_p), method='Nelder-Mead',
                           options={'xatol': 1e-6, 'fatol': 1e-11, 'maxiter': 20000, 'maxfev': 20000})
            x = res.x
            print('round', rnd, '%.10f' % res.fun, np.round(np.exp(x), 5), flush=True)
        z = np.exp(x)
        save({sysname + '_atom': {'s': list(z[:4]), 'p': list(z[4:]), 'E_radial_avg': float(res.fun)}})
    elif cmd == 'pol_sto':
        sysname = sys.argv[2]
        x0 = np.log([2.756 * SYSTEMS[sysname]['Z'] / 8.0] + ([1.579] if sysname == 'F' else []))
        res = minimize(pol_objective_sto, x0, args=(sysname, load()[sysname + '_atom']), method='Nelder-Mead',
                       options={'xatol': 2e-3, 'fatol': 2e-6, 'maxiter': 60})
        save({sysname + '_sto_pol': list(np.exp(res.x)), sysname + '_sto_pol_E_fc_ccsd': float(res.fun)})
    elif cmd == 'pol_gto':
        sysname = sys.argv[2]
        x0 = np.log([1.6 if sysname == 'F' else 2.2] + ([0.727] if sysname == 'F' else []))
        res = minimize(pol_objective_gto, x0, args=(sysname,), method='Nelder-Mead',
                       options={'xatol': 2e-3, 'fatol': 2e-6, 'maxiter': 60})
        save({sysname + '_gto_pol': list(np.exp(res.x)), sysname + '_gto_pol_E_fc_ccsd': float(res.fun)})
    elif cmd == 'final':
        final(sys.argv[2], sys.argv[3])
    print(cmd, 'done in %.0f s' % (time.time() - t0), flush=True)
