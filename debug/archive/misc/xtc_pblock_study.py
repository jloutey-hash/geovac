"""
DECISIVE p-block xTC sparsity study.
====================================
Deliverable (parallels the Li p-inclusive study, now with a p-BLOCK reference so
the leg-3 L'=0 collapse is switched OFF):
  (A) full radial-weighted s+p engine on carbon (+ Be s-ref control): gates +
      angular density / fill-in / 1-norm / Pauli, PLAIN vs xTC.
  (B) radial-free EXACT angular fill-in across s+p / s+p+d / s+p+d+f bases -- the
      model-independent core (angular selection is radial-independent).
Contrast throughout with Li's / s-reference's 0 fill-in.
"""
import os, sys, json, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import xtc_poc_li_pinclusive as P
from xtc_pblock_engine import REFERENCES, ref_angular
from xtc_pblock_fast_angular import fast_angular, orbs, REF_ANG
from xtc_poc_li_pinclusive_study import (deliverable, threebody_density,
                                         gate_v3_mconservation, gate_contraction_fidelity)


def run_engine(name, gamma=1.0, k=2.0, Ng=800, nx=160):
    spec = REFERENCES[name]
    o = P.assemble(max_n=2, k=k, gamma=gamma, Z=spec['Z'], n_elec=spec['n_elec'],
                   Ng=Ng, nx=nx, want_exact3=True, ref_occ=spec['ref_occ'])
    badm, nzv3 = gate_v3_mconservation(o)
    fid = gate_contraction_fidelity(o)
    D = deliverable(o)
    tb = threebody_density(o)
    return o, dict(
        name=name, Z=spec['Z'], n_elec=spec['n_elec'], ref_occ=list(spec['ref_occ']),
        ndet=o['ndet'], gamma=gamma,
        E_plain=o['E_plain'], E_TC2=o['E_TC2'], E_exactTC=o['E_exactTC'],
        E_xTC=o['E_xTC'], imag_xTC=o['imag_xTC'],
        v3_m_violations=badm, v3_nonzero=nzv3, contraction_fidelity=fid,
        threebody_density=tb['density'],
        coulomb_spatial=D['coulomb_spatial'],
        v2_vs_coulomb=D['v2_vs_coulomb'],
        op_metrics=D['op_metrics'], pauli=D['pauli'])


def geminal_zero(name, k=2.0):
    spec = REFERENCES[name]; rows = []
    for g in (4.0, 10.0, 25.0):
        o = P.assemble(max_n=2, k=k, gamma=g, Z=spec['Z'], n_elec=spec['n_elec'],
                       Ng=700, nx=140, want_exact3=True, ref_occ=spec['ref_occ'])
        rows.append(dict(gamma=g, E_plain=o['E_plain'], E_xTC=o['E_xTC'],
                         imag=o['imag_xTC'], diff=abs(o['E_xTC'] - o['E_plain'])))
    return rows


def main():
    t0 = time.time()
    rep = {}

    print("=" * 72)
    print("(A) RADIAL-WEIGHTED s+p ENGINE  (carbon p-ref vs Be s-ref control)")
    print("=" * 72)
    rep['engine'] = {}
    engines = {}
    for name in ('Be_1S', 'C_3P', 'O_3P'):
        t = time.time()
        o, r = run_engine(name)
        engines[name] = o
        rep['engine'][name] = r
        om = r['op_metrics']
        vc = r['v2_vs_coulomb']
        print(f"\n[{name}] Z={r['Z']} n_elec={r['n_elec']} ndet={r['ndet']}  "
              f"({time.time()-t:.0f}s)")
        print(f"  energies: plain={r['E_plain']:.5f}  TC2={r['E_TC2']:.5f}  "
              f"exactTC={r['E_exactTC']:.5f}  xTC={r['E_xTC']:.5f} (imag {r['imag_xTC']:.1e})")
        print(f"  GATES: m-viol {r['v3_m_violations']}/{r['v3_nonzero']}   "
              f"contraction fidelity |xTC-exactTC| = {r['contraction_fidelity']:.2e}")
        print(f"  3-body L3 raw density = {r['threebody_density']*100:.2f}%")
        print(f"  Coulomb spatial 2-body: {r['coulomb_spatial']['nnz']}/"
              f"{r['coulomb_spatial']['total']} ({r['coulomb_spatial']['density']*100:.2f}%)")
        print(f"  >>> FILL-IN (v2 nonzero where Coulomb zero): {vc['n_fill_in']}  "
              f"(v2 support {vc['n_supp_v2']}, Coulomb {vc['n_supp_coul']}) <<<")
        print(f"  1-norm[2b]  plain={om['plain']['l1_2body']:.2f}  xTC={om['xTC']['l1_2body']:.2f}"
              f"  ratio {om['xTC']['l1_2body']/om['plain']['l1_2body']:.3f}")
        print(f"  1-norm[tot] plain={om['plain']['l1_total']:.2f}  xTC={om['xTC']['l1_total']:.2f}"
              f"  ratio {om['xTC']['l1_total']/om['plain']['l1_total']:.3f}")
        if 'error' not in r['pauli']:
            pp, px = r['pauli']['plain'], r['pauli']['xTC']
            print(f"  Pauli terms plain={pp['n_pauli']} xTC={px['n_pauli']} "
                  f"(ratio {px['n_pauli']/pp['n_pauli']:.3f});  "
                  f"Pauli L1 ratio {px['l1_pauli']/pp['l1_pauli']:.3f}")

    print("\n" + "=" * 72)
    print("[G5] geminal->0 sanity (E_xTC -> E_plain) for carbon C_3P")
    gz = geminal_zero('C_3P')
    for row in gz:
        print(f"   g={row['gamma']:5.1f}: plain={row['E_plain']:.5f} xTC={row['E_xTC']:.5f}"
              f"  |d|={row['diff']:.2e}  imag={row['imag']:.1e}")
    rep['geminal_zero_C'] = gz

    print("\n" + "=" * 72)
    print("(B) EXACT ANGULAR FILL-IN across bases (radial-free, model-independent)")
    print("    s-reference vs p-reference; this is where the L'=0 collapse is OFF")
    print("=" * 72)
    rep['fast_angular'] = {}
    for lmax in (1, 2, 3):
        tag = {1: 's+p', 2: 's+p+d', 3: 's+p+d+f'}[lmax]
        rep['fast_angular'][tag] = {}
        print(f"\n  basis {tag} ({len(orbs(lmax))} angular orbitals):")
        for name in ('Be_s2 (s-ref)', 'C_2p2_Hund (p-ref)', 'O_2p4 (p-ref)'):
            r = fast_angular(lmax, REF_ANG[name])
            rep['fast_angular'][tag][name] = dict(
                n_coulomb=r['n_coulomb'], n_l3=r['n_l3'], n_fill_in=r['n_fill_in'],
                total=r['total'], coul_density=r['coul_density'],
                l3_density=r['l3_density'], fill_frac=r['fill_frac'],
                refmult=[list(x) for x in r['refmult']],
                fill_examples=r['fill_examples'][:6])
            print(f"    {name:22s} refmult={r['refmult']}  "
                  f"Coul {r['coul_density']*100:5.2f}%  L3 {r['l3_density']*100:5.2f}%  "
                  f"FILL {r['n_fill_in']:3d} ({r['fill_frac']*100:.3f}%)")

    rep['walltime_s'] = time.time() - t0
    os.makedirs('debug/data', exist_ok=True)
    with open('debug/data/xtc_pblock_study.json', 'w') as f:
        json.dump(rep, f, indent=2,
                  default=lambda x: list(x) if isinstance(x, set) else float(x))
    print(f"\nwrote debug/data/xtc_pblock_study.json  ({rep['walltime_s']:.0f}s)")


if __name__ == '__main__':
    main()
