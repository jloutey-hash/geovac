"""Scale ladder (record SL): water at 24 functions with the two-scale Sturmian kit's primitive
types, releasing the number of distinct radial exponents on O (and, in one rung, on H).  Every
exponent is set on the O atom (radial average-of-configuration HF, tests/test_water_ccsd.radial_atom)
or on H2 (RHF); nothing is tuned on water.

    python debug/scale_ladder.py atom            O-atom ladder (control; prints exponents)
    python debug/scale_ladder.py h2              H2 ladder (control)
    python debug/scale_ladder.py water <rung>    one water run -> debug/data/scale_ladder.jsonl

O primitives (n, l): s (1,0) (1,0) (2,0); p (2,1) (3,1); d (3,2) at the 3p exponent.
H primitives: (1,0) (2,0) (2,1), one exponent kH unless the H rung releases them.
Rung layouts (which primitives share an exponent) are in RUNGS.
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

# numpy 2.4.0 on Python 3.14.0 leaks memory in `.flat` slice assignment, which np.diag uses to build a
# matrix from a vector (measured: ~2.8 KB per call, never freed; the atom optimizer reached 18 GB).
# Same result through fancy indexing, which does not leak.  Process-local patch, 1-D input only.
_np_diag = np.diag


def _diag(v, k=0):
    v = np.asarray(v)
    if v.ndim != 1 or k != 0:
        return _np_diag(v, k)
    out = np.zeros((len(v), len(v)), dtype=v.dtype)
    i = np.arange(len(v))
    out[i, i] = v
    return out


np.diag = _diag

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, 'data')
spec = importlib.util.spec_from_file_location('w', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)
GE = w.GE

KC, KV = 7.79324757, 2.28124155          # v5.20.12 two-scale values (atom-optimized)
O_PRIMS = [(1, 0), (1, 0), (2, 0), (2, 1), (3, 1)]
# rung -> group index of each O primitive (same index = shared exponent)
RUNGS = {
    'R2': [0, 1, 1, 1, 1],
    'R3': [0, 1, 1, 2, 2],
    'R4': [0, 1, 2, 3, 3],
    'R5': [0, 1, 2, 3, 4],
}
H_PRIMS = [(1, 0), (2, 0), (2, 1)]


def atom_energy(ks, groups):
    s = [(n, ks[g]) for (n, l), g in zip(O_PRIMS, groups) if l == 0]
    p = [(n, ks[g]) for (n, l), g in zip(O_PRIMS, groups) if l == 1]
    try:
        return w.radial_atom(s, p)[0]
    except RuntimeError:                      # a wild restart; the radial HF did not converge
        return 1e3


def optimize_atom(groups, starts=None):
    ng = max(groups) + 1
    base = [KC, KV] + [KV] * (ng - 2)
    cand = [base]
    rng = np.random.default_rng(0)
    for _ in range(40):
        cand.append(list(np.array(base) * np.exp(rng.normal(0, 0.5, ng))))
    best = None
    for x0 in cand:
        r = minimize(lambda k: atom_energy(np.abs(k), groups), x0, method='Nelder-Mead',
                     options=dict(xatol=1e-8, fatol=1e-11, maxiter=6000))
        if best is None or r.fun < best.fun:
            best = r
    return best.fun, np.abs(best.x)


def h2_energy(ks, groups, prims=H_PRIMS, R=1.4, ng=12):
    nuc = [(np.zeros(3), 1.0), (np.array([0, 0, R]), 1.0)]
    fns, blocks = [], []
    for C, _ in nuc:
        for (n, l), g in zip(prims, groups):
            w.sto_shell(C, n, l, ks[g], ng, fns, blocks)
    S, h, g = GE.integral_set_md(fns, nuc)
    return w.rhf(S, h, g, 2, w.vnn(nuc))[0]


def optimize_h2(groups, n_gauss=6):
    """Exponents set at n_gauss Gaussians per STO (12 costs ~40-80 s per call with 2s and 2p on H);
    the caller re-evaluates the optimum at 12."""
    ng = max(groups) + 1
    best = None
    for x0 in ([1.0] * ng, [1.3, 0.8, 1.2][:ng], [0.9, 1.4, 1.6][:ng]):
        r = minimize(lambda k: h2_energy(np.abs(k), groups, ng=n_gauss), x0, method='Nelder-Mead',
                     options=dict(xatol=1e-5, fatol=1e-10, maxiter=2000))
        if best is None or r.fun < best.fun:
            best = r
    return best.fun, np.abs(best.x)


def kit(nuclei, o_ks, o_groups, h_ks, h_groups, n_gauss=12, ng_pol=8):
    """O: the five primitives at their rung exponents plus 3d at the 3p exponent; H: three
    primitives.  Same builder as the v5.20.12 kits (_shell_kit with one primitive per function)."""
    fns, cols_of = [], []
    k3p = o_ks[o_groups[4]]
    for center, Z in nuclei:
        if Z == 8.0:
            for (n, l), g in zip(O_PRIMS, o_groups):
                w._shell_kit(center, [(l, {n: 1.0}, o_ks[g])], n_gauss, ng_pol, fns, [], cols_of)
            w._shell_kit(center, [(2, {3: 1.0}, k3p)], n_gauss, ng_pol, fns, [], cols_of)
        else:
            for (n, l), g in zip(H_PRIMS, h_groups):
                w._shell_kit(center, [(l, {n: 1.0}, h_ks[g])], n_gauss, ng_pol, fns, [], cols_of)
    return fns, w._kit_columns(fns, cols_of)


def water(name, o_ks, o_groups, h_ks, h_groups):
    t0 = time.time()
    nuc = w.water_nuclei()
    fns, T = kit(nuc, o_ks, o_groups, h_ks, h_groups)
    S, h, g = w.ao_integrals(fns, nuc, T)
    enuc = w.vnn(nuc)
    t1 = time.time()
    r = w.run(S, h, g, 10, enuc)
    rf = w.run(S, h, g, 10, enuc, frozen_core=1)
    row = {'rung': name, 'nbf': T.shape[1], 'o_ks': list(map(float, o_ks)), 'o_groups': o_groups,
           'h_ks': list(map(float, h_ks)), 'h_groups': h_groups,
           **{k: float(v) for k, v in r.items()}, **{k + '_fc': float(v) for k, v in rf.items() if k != 'RHF'},
           't_integrals': round(t1 - t0, 1), 't_solve': round(time.time() - t1, 1)}
    print(json.dumps(row), flush=True)
    with open(os.path.join(DATA, 'scale_ladder.jsonl'), 'a') as f:
        f.write(json.dumps(row) + '\n')
    return row


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'atom':                         # optional rung names; results merge into the json
        path = os.path.join(DATA, 'scale_ladder_atom.json')
        out = json.load(open(path)) if os.path.exists(path) else {}
        for name in (sys.argv[2:] or list(RUNGS)):
            t0 = time.time()
            E, ks = optimize_atom(RUNGS[name])
            out[name] = dict(E=float(E), ks=list(map(float, ks)), groups=RUNGS[name])
            print(name, '%.6f' % E, np.round(ks, 4), '%.0f s' % (time.time() - t0), flush=True)
            with open(path + '.tmp', 'w') as f:
                json.dump(out, f, indent=1)
            os.replace(path + '.tmp', path)
    elif cmd == 'h2':
        out = {}
        path = os.path.join(DATA, 'scale_ladder_h2.json')
        out = json.load(open(path)) if os.path.exists(path) else {}
        H_RUNGS = {'H1': [0, 0, 0], 'H2': [0, 1, 1], 'H3': [0, 1, 2]}
        for name in (sys.argv[2:] or list(H_RUNGS)):
            t0 = time.time()
            groups = H_RUNGS[name]
            E6, ks = optimize_h2(groups)
            E12 = h2_energy(ks, groups)
            out[name] = dict(E_ng6=float(E6), E=float(E12), ks=list(map(float, ks)), groups=groups)
            print(name, 'ng6 %.6f ng12 %.6f' % (E6, E12), np.round(ks, 4), '%.0f s' % (time.time() - t0), flush=True)
            with open(path + '.tmp', 'w') as f:
                json.dump(out, f, indent=1)
            os.replace(path + '.tmp', path)
        if 'k1' not in out:
            out['k1'] = dict(E=float(h2_energy([1.0], [0, 0, 0])))
            print('H1 at k = 1 (the kit): %.6f' % out['k1']['E'], flush=True)
            with open(path, 'w') as f:
                json.dump(out, f, indent=1)
    elif cmd == 'water':
        rung = sys.argv[2]
        if rung == 'R2rep':                       # control: reproduce the v5.20.12 two-scale run
            water('R2rep', [KC, KV], RUNGS['R2'], [1.0], [0, 0, 0])
            sys.exit(0)
        atom = json.load(open(os.path.join(DATA, 'scale_ladder_atom.json')))
        if rung in RUNGS:
            water(rung, atom[rung]['ks'], RUNGS[rung], [1.0], [0, 0, 0])
        elif rung == 'R5H':
            h2 = json.load(open(os.path.join(DATA, 'scale_ladder_h2.json')))['H3']
            water('R5H', atom['R5']['ks'], RUNGS['R5'], h2['ks'], h2['groups'])
