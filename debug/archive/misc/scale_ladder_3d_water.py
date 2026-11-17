"""SL-3d ladder, water side (record docs/predictions/2026-11-01_scale_ladder_3d.md).

Water at 24 functions, same geometry, builder and solver as the registered ladder (debug/scale_ladder.py).
The only change from a registered rung: the O 3d exponent is the atom-set value in
debug/data/scale_ladder_3d_atom_<rung>.json (fc-CCSD(T) optimum of the free O atom, s/p fixed), not the
3p exponent.  Nothing is tuned on water.

    python debug/scale_ladder_3d_water.py R2d     control (not scored; run in Phase 1)
    python debug/scale_ladder_3d_water.py R3d | R4d | R5d | R5Hd     scored (Phase 2 only)

Rows append to debug/data/scale_ladder_3d.jsonl.
"""
import importlib.util
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('sl', os.path.join(HERE, 'scale_ladder.py'))
sl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sl)          # patches np.diag, loads tests/test_water_ccsd as sl.w
w = sl.w
DATA = sl.DATA
NAMES = {'R2d': 'R2', 'R3d': 'R3', 'R4d': 'R4', 'R5d': 'R5', 'R5Hd': 'R5'}


def kit(nuclei, o5, k3d, h_ks, h_groups, n_gauss=12, ng_pol=8):
    fns, cols_of = [], []
    for center, Z in nuclei:
        if Z == 8.0:
            for (n, l), k in zip(sl.O_PRIMS, o5):
                w._shell_kit(center, [(l, {n: 1.0}, k)], n_gauss, ng_pol, fns, [], cols_of)
            w._shell_kit(center, [(2, {3: 1.0}, k3d)], n_gauss, ng_pol, fns, [], cols_of)
        else:
            for (n, l), g in zip(sl.H_PRIMS, h_groups):
                w._shell_kit(center, [(l, {n: 1.0}, h_ks[g])], n_gauss, ng_pol, fns, [], cols_of)
    return fns, w._kit_columns(fns, cols_of)


if __name__ == '__main__':
    name = sys.argv[1]
    rung = NAMES[name]
    at = json.load(open(os.path.join(DATA, 'scale_ladder_3d_atom_%s.json' % rung)))
    o5, k3d = at['s_p_exponents'], at['k3d_new']
    if name == 'R5Hd':
        h2 = json.load(open(os.path.join(DATA, 'scale_ladder_h2.json')))['H3']
        h_ks, h_groups = h2['ks'], h2['groups']
    else:
        h_ks, h_groups = [1.0], [0, 0, 0]
    t0 = time.time()
    nuc = w.water_nuclei()
    fns, T = kit(nuc, o5, k3d, h_ks, h_groups)
    S, h, g = w.ao_integrals(fns, nuc, T)
    enuc = w.vnn(nuc)
    t1 = time.time()
    r = w.run(S, h, g, 10, enuc)
    rf = w.run(S, h, g, 10, enuc, frozen_core=1)
    row = {'rung': name, 'nbf': T.shape[1], 'o_k_1s_1s_2s_2p_3p': o5, 'k3d': k3d,
           'h_ks': list(map(float, h_ks)), 'h_groups': h_groups,
           **{k: float(v) for k, v in r.items()}, **{k + '_fc': float(v) for k, v in rf.items() if k != 'RHF'},
           't_integrals': round(t1 - t0, 1), 't_solve': round(time.time() - t1, 1)}
    print(json.dumps(row), flush=True)
    with open(os.path.join(DATA, 'scale_ladder_3d.jsonl'), 'a') as f:
        f.write(json.dumps(row) + '\n')
