"""Independent check against Stark & Zimmerman, arXiv:2510.15204 (2025), Table 3: HF heat-bath CI correlation
energy in their TZ Slater basis, -0.2558 Ha (frozen core, eps1 = 1e-4, eps2 = 1e-7, RI integrals, not extrapolated).

Their TZ basis (SI Table S1, read from the downloaded PDF, layout text): uncontracted, nodeless (every s is 1s-type,
p 2p-type, d 3d-type), spherical harmonics.
  F: s 15.520 8.466 4.619 2.520 1.375 0.750 | p 5.667 3.063 1.656 | d 4.000 2.000          (25 functions)
  H: s 2.280 1.520 1.013 | p 2.100 1.500                                                   (9 functions)
Geometry (SI Table S2): H at the origin, F at 0.9170 A on z: R = 1.73288 bohr.
Ours: the validated RHF + spin-orbital CCSD(T) of tests/test_water_ccsd.py, F 1s frozen, Slater functions as
Gaussian fits (NG per function).  Comparison is of the CORRELATION energy only (they report no SCF energies).

    python debug/stark_zimmerman_hf_tz.py control   # evaluator control: HF Slater DZP at 8 Gaussians vs 12 (-100.287713)
    python debug/stark_zimmerman_hf_tz.py tz 8|12   # their TZ set
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
_s = importlib.util.spec_from_file_location('w', os.path.join(ROOT, 'tests', 'test_water_ccsd.py'))
w = importlib.util.module_from_spec(_s)
_s.loader.exec_module(w)
OUT = os.path.join(HERE, 'data', 'stark_zimmerman_hf_tz.jsonl')
R_SZ = 0.9170 / 0.529177210903
F_S = [15.520, 8.466, 4.619, 2.520, 1.375, 0.750]
F_P = [5.667, 3.063, 1.656]
F_D = [4.000, 2.000]
H_S = [2.280, 1.520, 1.013]
H_P = [2.100, 1.500]
SZ_ECORR = -0.2558


def write(row):
    print(json.dumps(row), flush=True)
    with open(OUT, 'a') as f:
        f.write(json.dumps(row) + '\n')


def tz_basis(ng):
    F = np.zeros(3)
    H = np.array([0.0, 0.0, R_SZ])
    fns, blocks = [], []
    for z in F_S:
        w.sto_shell(F, 1, 0, z, ng, fns, blocks)
    for z in F_P:
        w.sto_shell(F, 2, 1, z, ng, fns, blocks)
    for z in F_D:
        w.sto_shell(F, 3, 2, z, ng, fns, blocks)
    for z in H_S:
        w.sto_shell(H, 1, 0, z, ng, fns, blocks)
    for z in H_P:
        w.sto_shell(H, 2, 1, z, ng, fns, blocks)
    nuc = [(F, 9.0), (H, 1.0)]
    return nuc, fns, w.spherical_map(fns, blocks)


def solve(nuc, fns, T, label, extra):
    t0 = time.time()
    S, h, g = w.ao_integrals(fns, nuc, T)
    enuc = w.vnn(nuc)
    r = w.run(S, h, g, 10, enuc, frozen_core=1)
    row = {'case': label, 'nbf': T.shape[1], 'RHF': float(r['RHF']), 'CCSD(T)_fc': float(r['CCSD(T)']),
           'Ecorr_fc': float(r['CCSD(T)'] - r['RHF']), 't': round(time.time() - t0, 1), **extra}
    write(row)
    return row


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'control':
        _h = importlib.util.spec_from_file_location('hf', os.path.join(HERE, 'hf_slater_matched.py'))
        hfm = importlib.util.module_from_spec(_h)
        _h.loader.exec_module(hfm)
        d = hfm.load()
        pol = d['F_sto_pol']
        fns, T = hfm.slater_dzp('F', d['F_atom'], 8, 8, pol[0], pol[1])
        solve(hfm.nuclei('F'), fns, T, 'control_hf_slater_dzp_ng8', {'ng': 8, 'ng12_reference_fc': -100.287713})
    elif cmd == 'tz':
        ng = int(sys.argv[2])
        nuc, fns, T = tz_basis(ng)
        assert T.shape[1] == 34, T.shape
        solve(nuc, fns, T, 'stark_zimmerman_tz', {'ng': ng, 'R_bohr': R_SZ, 'their_Ecorr_HBCI': SZ_ECORR})
