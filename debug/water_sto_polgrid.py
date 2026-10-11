"""Slater DZP polarization exponents (O 3d, H 2p) on water CCSD by a 3x3 grid +
quadratic fit; records all-electron and frozen-core (O 1s) CCSD (replaces pol_sto's Nelder-Mead, ~20 min per point at 12 Gaussians).
Evaluator: 8 Gaussians on s/p (a near-constant core-cusp offset), 6 on polarization.
Usage: python debug/water_sto_polgrid.py [d1 p1 d2 p2 ...]   -> appends debug/data/water_sto_polgrid.jsonl"""
import itertools, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_s = importlib.util.spec_from_file_location('x', os.path.join(ROOT, 'debug', 'water_sto_exponents.py'))
x = importlib.util.module_from_spec(_s); _s.loader.exec_module(x)
w = x.w
OUT = os.path.join(ROOT, 'debug', 'data', 'water_sto_polgrid.jsonl')
NG, NGP = 8, 6


def point(dz, pz):
    d = dict(x.load(), O3d=dz, H2p=pz)
    nuc = w.water_nuclei()
    fns, T = w.slater_dzp_contracted(nuc, x.contracted_ex(d), NG, ng_pol=NGP)
    S, h, g = w.ao_integrals(fns, nuc, T)
    E = w.run(S, h, g, 10, w.vnn(nuc), triples=False)['CCSD']
    Efc = w.run(S, h, g, 10, w.vnn(nuc), triples=False, frozen_core=1)['CCSD']
    row = {'O3d': dz, 'H2p': pz, 'ng': NG, 'ng_pol': NGP, 'CCSD': float(E), 'CCSD_fc': float(Efc)}
    print(json.dumps(row)); sys.stdout.flush()
    open(OUT, 'a').write(json.dumps(row) + '\n')
    return E


if __name__ == '__main__':
    t0 = time.time()
    if len(sys.argv) >= 3:      # pairs: d1 p1 d2 p2 ...
        a = [float(v) for v in sys.argv[1:]]
        for dz, pz in zip(a[::2], a[1::2]):
            point(dz, pz)
    else:
        for dz, pz in itertools.product((1.5, 1.8, 2.2), (0.8, 1.0, 1.25)):
            point(dz, pz)
    print('done %.0fs' % (time.time() - t0))
