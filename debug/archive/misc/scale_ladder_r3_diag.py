"""Diagnostic for the R3 non-monotone step of the scale ladder (record SL).  See
debug/scale_ladder_r3_diag_memo.md.  Same builders as debug/scale_ladder.py; the only new freedom is
that the 3p and 3d exponents can be set separately from the 2p, and each value used is one the ladder
already holds (R2's 2.28124 or R3's own values).  Nothing is tuned on water.

    python debug/scale_ladder_r3_diag.py oion            O^2- closed-shell one-centre probe (fast)
    python debug/scale_ladder_r3_diag.py water <name>    one water run -> debug/data/scale_ladder_r3_diag.jsonl

Names: R3dB (R3 with 3d at 2.28124), R3vB (3p and 3d at 2.28124), R3sB (s-valence at 2.28124).
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('sl', os.path.join(HERE, 'scale_ladder.py'))
sl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sl)          # patches np.diag, loads tests/test_water_ccsd as sl.w
w = sl.w
DATA = sl.DATA
KV = sl.KV
ATOM = json.load(open(os.path.join(DATA, 'scale_ladder_atom.json')))
R3 = ATOM['R3']['ks']                 # [1s core, s-valence, p]
KC3, KS3, KP3 = R3

# name -> (k1s, k1s', k2s, k2p, k3p, k3d)
SETS = {
    'R2':   (sl.KC, KV, KV, KV, KV, KV),
    'R3':   (KC3, KS3, KS3, KP3, KP3, KP3),
    'R3dB': (KC3, KS3, KS3, KP3, KP3, KV),
    'R3vB': (KC3, KS3, KS3, KP3, KV, KV),
    'R3sB': (KC3, KV, KV, KP3, KP3, KP3),
}
O_PRIMS = [(1, 0), (1, 0), (2, 0), (2, 1), (3, 1)]


def kit(nuclei, oks, h_k=1.0, n_gauss=12, ng_pol=8):
    fns, cols_of = [], []
    for center, Z in nuclei:
        if Z == 8.0:
            for (n, l), k in zip(O_PRIMS, oks[:5]):
                w._shell_kit(center, [(l, {n: 1.0}, k)], n_gauss, ng_pol, fns, [], cols_of)
            w._shell_kit(center, [(2, {3: 1.0}, oks[5])], n_gauss, ng_pol, fns, [], cols_of)
        else:
            for (n, l) in sl.H_PRIMS:
                w._shell_kit(center, [(l, {n: 1.0}, h_k)], n_gauss, ng_pol, fns, [], cols_of)
    return fns, w._kit_columns(fns, cols_of)


def energies(nuclei, oks, nelec):
    fns, T = kit(nuclei, oks)
    S, h, g = w.ao_integrals(fns, nuclei, T)
    enuc = w.vnn(nuclei)
    r = w.run(S, h, g, nelec, enuc)
    rf = w.run(S, h, g, nelec, enuc, frozen_core=1)
    return {'nbf': T.shape[1], **{k: float(v) for k, v in r.items()},
            **{k + '_fc': float(v) for k, v in rf.items() if k != 'RHF'}}


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'oion':
        nuc = [(np.zeros(3), 8.0)]
        out = {}
        for name in ('R2', 'R3', 'R3dB', 'R3vB', 'R3sB'):
            t0 = time.time()
            out[name] = energies(nuc, SETS[name], 10)
            e = out[name]
            print('%-5s RHF %.6f  fcCCSD %.6f  fcCCSD(T) %.6f  corr_fc %.6f  %.0f s' % (
                name, e['RHF'], e['CCSD_fc'], e['CCSD(T)_fc'], e['CCSD(T)_fc'] - e['RHF'], time.time() - t0),
                flush=True)
        with open(os.path.join(DATA, 'scale_ladder_r3_diag_oion.json'), 'w') as f:
            json.dump({'sets': SETS, 'energies': out}, f, indent=1)
    elif cmd == 'water':
        name = sys.argv[2]
        t0 = time.time()
        nuc = w.water_nuclei()
        row = {'rung': name, 'o_k_1s_1s_2s_2p_3p_3d': list(SETS[name]), 'h_k': 1.0}
        row.update(energies(nuc, SETS[name], 10))
        row['t_total'] = round(time.time() - t0, 1)
        print(json.dumps(row), flush=True)
        with open(os.path.join(DATA, 'scale_ladder_r3_diag.jsonl'), 'a') as f:
            f.write(json.dumps(row) + '\n')
