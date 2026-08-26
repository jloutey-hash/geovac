"""
xTC PoC study on Li (1s^2 2s doublet) in the Coulomb-Sturmian basis.

Produces:
  (1) gamma-scan at fixed (ns,k): TC2 / exact-3body-TC / xTC energies
  (2) accuracy-vs-basis: plain FCI vs xTC across ns (s-only), k optimized for plain
  (3) qubit / sparsity metrics of the effective 2-body operator (plain vs xTC)
  (4) plain full-l FCI (SturmianCI, includes p/d) for the qubit-count comparison
Reference: Li non-rel exact -7.478060 Ha ; HF limit -7.432727 Ha (corr -45.3 mHa).
"""
import json, time
import numpy as np
import xtc_poc_li as X

TRUE_LI = -7.478060
HF_LI = -7.432727
NG = 1000
NX = 160

def gamma_scan(ns=3, k=1.5, gammas=(0.6, 0.8, 1.0, 1.2, 1.5, 2.0)):
    rows = []
    for g in gammas:
        o = X.assemble(ns, k, g, Z=3, n_elec=3, Ng=NG, nx=NX, want_exact3=True)
        rows.append(dict(gamma=g, E_plain=o['E_plain'], E_TC2=o['E_TC2'],
                         E_exactTC=o['E_exactTC'], E_xTC=o['E_xTC'],
                         imag_xTC=o['imag_xTC'],
                         d3_exact_mHa=(o['E_exactTC']-o['E_TC2'])*1e3,
                         d3_xtc_mHa=(o['E_xTC']-o['E_TC2'])*1e3))
        print(f"  g={g:.2f}: plain={o['E_plain']:.5f} TC2={o['E_TC2']:.5f} "
              f"exactTC={o['E_exactTC']:.5f} xTC={o['E_xTC']:.5f}  "
              f"3b(exact/xtc)={rows[-1]['d3_exact_mHa']:+.2f}/{rows[-1]['d3_xtc_mHa']:+.2f} mHa")
    return rows

def optimize_k_plain(ns, ks):
    best = (1e9, None, None)
    for k in ks:
        E, m, nso, nd = X.plain_fci(ns, k, Ng=NG, nx=NX)
        if E < best[0]:
            best = (E, k, m)
    return best  # (E, k*, metrics)

def accuracy_vs_basis(gamma, ns_list=(2, 3, 4, 5), ks=(1.1, 1.3, 1.5, 1.7, 1.9)):
    rows = []
    for ns in ns_list:
        Ep, kstar, mp = optimize_k_plain(ns, ks)
        o = X.assemble(ns, kstar, gamma, Z=3, n_elec=3, Ng=NG, nx=NX, want_exact3=(ns <= 5))
        row = dict(ns=ns, qubits=2*ns, ndet=o['ndet'], k=kstar, gamma=gamma,
                   E_plain=Ep, E_TC2=o['E_TC2'], E_exactTC=o.get('E_exactTC'),
                   E_xTC=o['E_xTC'], imag_xTC=o['imag_xTC'],
                   err_plain_mHa=(Ep-TRUE_LI)*1e3, err_xTC_mHa=(o['E_xTC']-TRUE_LI)*1e3,
                   corr_plain_mHa=(Ep-HF_LI)*1e3, corr_xTC_mHa=(o['E_xTC']-HF_LI)*1e3,
                   metrics_plain=mp, metrics_xTC=o['metrics_xTC'])
        rows.append(row)
        print(f"  ns={ns} q={2*ns} k*={kstar}: plain={Ep:.5f} ({row['err_plain_mHa']:+.1f}mHa) "
              f"xTC={o['E_xTC']:.5f} ({row['err_xTC_mHa']:+.1f}mHa)  "
              f"gain={ (o['E_xTC']-Ep)*1e3:+.2f}mHa  "
              f"L1[2b] plain={mp['l1_2body']:.2f} xTC={o['metrics_xTC']['l1_2body']:.2f}")
    return rows

def s_limit(ns_list=(5, 6), ks=(1.3, 1.5, 1.7, 1.9)):
    out = []
    for ns in ns_list:
        E, k, m = optimize_k_plain(ns, ks)
        out.append(dict(ns=ns, k=k, E_plain=E, err_mHa=(E-TRUE_LI)*1e3))
        print(f"  plain s-only ns={ns}: E={E:.6f} (k={k}) err={(E-TRUE_LI)*1e3:+.1f} mHa")
    return out

def plain_full_l():
    """Plain full-l FCI via production SturmianCI (includes p/d) for qubit comparison."""
    import sys
    sys.path.insert(0, '..')
    from geovac.sturmian_solver import SturmianCI
    out = []
    for mx in [2]:                      # max_n>=3 FCI is too slow in pure-Python SC; q-count noted in memo
        try:
            t0 = time.time()
            ci = SturmianCI(Z=3, n_electrons=3, max_n=mx)
            if ci.n_sd > 5000:
                print(f"  full-l max_n={mx}: n_det={ci.n_sd} too big, skip"); continue
            res = ci.optimize_k(k_range=(1.0, 3.0), n_scan=10)
            out.append(dict(max_n=mx, n_spatial=ci.n_spatial, qubits=2*ci.n_spatial,
                            ndet=ci.n_sd, E=res['energy'], k=res['k'],
                            err_mHa=(res['energy']-TRUE_LI)*1e3))
            print(f"  full-l max_n={mx}: q={2*ci.n_spatial} n_det={ci.n_sd} E={res['energy']:.5f} "
                  f"err={(res['energy']-TRUE_LI)*1e3:+.1f} mHa  [{time.time()-t0:.0f}s]")
        except Exception as e:
            print(f"  full-l max_n={mx}: FAILED {e}")
    return out

if __name__ == '__main__':
    t0 = time.time()
    report = {'true_Li': TRUE_LI, 'HF_Li': HF_LI, 'grid': {'Ng': NG, 'nx': NX}}
    print("\n[1] gamma-scan (ns=3, k=1.5)")
    report['gamma_scan'] = gamma_scan(ns=3, k=1.5)
    # pick representative gamma = min exact-3body-TC energy (best correlation proxy)
    gstar = min(report['gamma_scan'], key=lambda r: r['E_exactTC'])['gamma']
    report['gamma_star'] = gstar
    print(f"  -> representative gamma* = {gstar}")

    print("\n[2] accuracy vs basis (s-only, k optimized per ns for plain)")
    report['accuracy'] = accuracy_vs_basis(gstar)

    print("\n[3] plain s-only saturation (s-limit)")
    report['s_limit'] = s_limit()

    print("\n[4] plain full-l FCI (production SturmianCI, includes p/d)")
    report['full_l'] = plain_full_l()

    report['walltime_s'] = time.time() - t0
    with open('data/xtc_poc_li_study.json', 'w') as f:
        json.dump(report, f, indent=2, default=float)
    print(f"\nwrote debug/data/xtc_poc_li_study.json  ({report['walltime_s']:.0f}s)")
