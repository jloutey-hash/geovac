"""
Accuracy validation of the "symmetrization escape hatch"
========================================================

Gate for the escape-hatch finding of sprint_tc_quantum_cost_measure_memo.md:
symmetrizing the full non-Hermitian TC operator  H~_sym = (H~ + H~^H)/2  drops the
anti-Hermitian part of the convective K (worth ~33 mHa at gamma=1 on Li).  Does that
move the energy TOWARD or AWAY from the exact non-relativistic energy?

TC IS NON-VARIATIONAL, so accuracy CANNOT be read off "lowest energy wins".  The
principled operating point is gamma-STATIONARITY (dE_TC/dgamma = 0), and a genuine
cusp fix must (i) move E toward exact WITHOUT overshooting past it and (ii) have its
correction PERSIST/GROW as the basis grows.

Systems: He (1s^2, s-only, 2-body-only: H~ = D + K) and Li (1s^2 2s, s-only,
genuine 3-body: H~ = D + K + xTC-L3).  Atomic only -- K is built/validated only in
the s-only engine.

Engines reused READ-ONLY:  xtc_poc_li.py (low-level builders + FCI + ground).
This driver only ORCHESTRATES them (caches gamma-independent one-body pieces, loops
gamma).  It does not modify any engine and re-implements no physics.

Outputs:  debug/data/xtc_sym_accuracy.json  + a table on stdout.
"""
import os, sys, json, time
import numpy as np
import scipy.linalg as sla

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
_ROOT = os.path.abspath(os.path.join(_HERE, '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import xtc_poc_li as X  # READ-ONLY reuse

# exact non-relativistic ground energies (Ha)
EXACT = {'He': -2.90372, 'Li': -7.47806}
SYS = {  # name: (Z, n_elec, with_L3, k-grid to optimise plain over)
    'He': dict(Z=2, n_elec=2, with_L3=False, kgrid=[1.6, 1.8, 2.0, 2.2, 2.4, 2.6]),
    'Li': dict(Z=3, n_elec=3, with_L3=True,  kgrid=[1.4, 1.5, 1.6, 1.7, 1.8]),
}
GAMMAS = [0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0]
NS_LIST = [3, 4, 5]
NG = 600
NX = 96


# ----------------------------------------------------------------------------
# gamma-independent one-body prep (cached per (system, ns, k))
# ----------------------------------------------------------------------------
def prep(ns, k, Z, ne, Ng=NG):
    r, wr = X.make_grid(k, Ng=Ng)
    S, h1s, Rtab, W = X.build_one_body(ns, r, wr, k, Z)
    Xm = X.lowdin(S)
    h1o = X.transform_1(h1s, Xm)
    nso = 2 * ns
    dets, didx = X.make_dets(nso, ne)
    hso = X.h_spin(h1o, nso)
    de = np.diag(h1o); order = np.argsort(de)
    ref = (2 * order[0], 2 * order[0] + 1) if ne == 2 else \
          (2 * order[0], 2 * order[0] + 1, 2 * order[1])
    return dict(r=r, wr=wr, S=S, Rtab=Rtab, W=W, Xm=Xm, h1o=h1o, nso=nso,
                dets=dets, didx=didx, hso=hso, ref=ref, ns=ns, ne=ne)


def plain_energy(P, nx=NX):
    """plain Coulomb FCI on the SAME grid as the TC scan (clean deltas)."""
    Km = X.build_kernels(P['r'], 1.0, nx=nx)   # coul is gamma-independent
    eri_coul, _, _ = X.two_body(P['ns'], P['Rtab'], P['W'], Km)
    asym = X.asym_from_phys(X.transform_2(eri_coul, P['Xm']), P['nso'])
    H = X.build_H(P['dets'], P['didx'], P['hso'], asym, P['nso'])
    E, _ = X.ground(H, hermitian=True)
    return E


def tc_energies(P, g, with_L3, nx=NX):
    """Full non-Hermitian TC ground E_TC (+imag) and symmetrized ground E_sym."""
    Km = X.build_kernels(P['r'], g, nx=nx)
    _, eri_w, eri_K = X.two_body(P['ns'], P['Rtab'], P['W'], Km)
    Xm = P['Xm']; nso = P['nso']
    asym_w = X.asym_from_phys(X.transform_2(eri_w, Xm), nso)
    asym_K = X.asym_from_phys(X.transform_2(eri_K, Xm), nso)
    if with_L3:
        V3o = X.transform_3(X.three_body(P['ns'], P['Rtab'], P['W'], Km), Xm)
        v2, v1, v0 = X.xtc_contract(V3o, nso, P['ref'])
        hso = P['hso'] + v1
        asym = asym_w + asym_K + v2
    else:
        hso = P['hso']; asym = asym_w + asym_K; v0 = 0.0
    H = X.build_H(P['dets'], P['didx'], hso, asym, nso, v0=v0)
    E_TC, im = X.ground(H, hermitian=False)
    Hs = 0.5 * (H + H.conj().T)
    E_sym = float(sla.eigh(Hs, eigvals_only=True)[0])
    return E_TC, im, E_sym


def optimise_k(ns, Z, ne, kgrid):
    best = (None, 1e9)
    for k in kgrid:
        E, _, _, _ = X.plain_fci(ns, k, Z=Z, n_elec=ne, Ng=700, nx=NX)
        if E < best[1]:
            best = (k, E)
    return best[0]


def find_stationary(gs, es):
    """Return (gamma*, dEdg at *, is_monotone).  Central differences; a stationary
    point is a sign change of dE/dgamma.  If none -> monotone, report the flattest."""
    gs = np.array(gs); es = np.array(es)
    dEdg = np.gradient(es, gs)
    sign = np.sign(dEdg)
    xings = np.where(np.diff(sign) != 0)[0]
    if xings.size:
        i = xings[0]
        # linear interp of the zero crossing
        g0, g1 = gs[i], gs[i + 1]
        d0, d1 = dEdg[i], dEdg[i + 1]
        gstar = float(g0 - d0 * (g1 - g0) / (d1 - d0)) if d1 != d0 else float(g0)
        return gstar, 0.0, False
    j = int(np.argmin(np.abs(dEdg)))
    return None, float(dEdg[j]), True, float(gs[j])


def interp_cross(gs, ys, target):
    """gamma at which y(gamma) == target (first crossing), or None."""
    gs = np.array(gs); ys = np.array(ys)
    d = ys - target
    for i in range(len(gs) - 1):
        if d[i] == 0:
            return float(gs[i])
        if d[i] * d[i + 1] < 0:
            g0, g1 = gs[i], gs[i + 1]
            return float(g0 - d[i] * (g1 - g0) / (d[i + 1] - d[i]))
    return None


def run_system(name):
    cfg = SYS[name]; Z, ne, wl3 = cfg['Z'], cfg['n_elec'], cfg['with_L3']
    exact = EXACT[name]
    out = dict(system=name, exact=exact, with_L3=wl3, per_ns={})
    for ns in NS_LIST:
        kstar = optimise_k(ns, Z, ne, cfg['kgrid'])
        P = prep(ns, kstar, Z, ne)
        E_plain = plain_energy(P)
        rows = []
        for g in GAMMAS:
            E_TC, im, E_sym = tc_energies(P, g, wl3)
            rows.append(dict(gamma=g, E_TC=E_TC, imag_TC=im, E_sym=E_sym,
                             d_plain_TC=E_plain - E_TC,
                             d_plain_sym=E_plain - E_sym,
                             sym_minus_TC=E_sym - E_TC,
                             dTC_exact=E_TC - exact,
                             dsym_exact=E_sym - exact))
            print("  %s ns=%d k=%.2f g=%4.2f  Ep=%.5f Etc=%.5f im=%.0e "
                  "Esym=%.5f  sym-tc=%+7.2f  Etc-ex=%+7.2f Esym-ex=%+7.2f mHa"
                  % (name, ns, kstar, g, E_plain, E_TC, im, E_sym,
                     (E_sym - E_TC) * 1e3, (E_TC - exact) * 1e3,
                     (E_sym - exact) * 1e3))
        gs = [r['gamma'] for r in rows]
        etc = [r['E_TC'] for r in rows]
        esym = [r['E_sym'] for r in rows]
        stat = find_stationary(gs, etc)
        gstar = stat[0]; monotone = stat[2]
        flat_g = stat[3] if len(stat) > 3 else None
        out['per_ns'][ns] = dict(
            k=kstar, E_plain=E_plain,
            gamma_star=gstar, monotone=monotone, flattest_gamma=flat_g,
            gamma_TC_crosses_exact=interp_cross(gs, etc, exact),
            gamma_sym_crosses_exact=interp_cross(gs, esym, exact),
            rows=rows)
    return out


def converged_plain(name):
    cfg = SYS[name]
    best = (None, 1e9)
    for k in [x for x in np.arange(1.3, 3.01, 0.1)]:
        E, _, _, _ = X.plain_fci(6, round(float(k), 2), Z=cfg['Z'],
                                 n_elec=cfg['n_elec'], Ng=800, nx=NX)
        if E < best[1]:
            best = (round(float(k), 2), E)
    return dict(ns=6, k=best[0], E_plain_converged=best[1])


def main():
    t0 = time.time()
    rep = dict(params=dict(Ng=NG, nx=NX, gammas=GAMMAS, ns_list=NS_LIST,
                           exact=EXACT),
               note=('E_TC = full non-Hermitian TC ground (scipy.linalg.eig, real GS). '
                     'E_sym = ground of (H~+H~^H)/2 (drops anti-Herm part of K). '
                     's-only plain FCI converges to a floor ~25-35 mHa ABOVE exact '
                     '(the residual is angular l>0 correlation, unreachable by a '
                     'radial cusp geminal). Any E below that floor / below exact is '
                     'non-variational drift, not a basis-honest cusp fix.'))
    for name in ('He', 'Li'):
        print('=== %s ===' % name)
        rep[name] = run_system(name)
        rep[name]['converged_plain_floor'] = converged_plain(name)
        print('  converged plain floor:', rep[name]['converged_plain_floor'])
    rep['walltime_s'] = time.time() - t0

    os.makedirs(os.path.join(_ROOT, 'debug', 'data'), exist_ok=True)
    outp = os.path.join(_ROOT, 'debug', 'data', 'xtc_sym_accuracy.json')
    with open(outp, 'w') as f:
        json.dump(rep, f, indent=2, default=float)
    print('\nwrote', outp, '  (%.0fs)' % rep['walltime_s'])

    # ---- summary table sliced at the memo's reference gamma=1.0 -----------
    print('\n' + '=' * 118)
    print('SUMMARY @ reference gamma=1.0  (deltas to exact in mHa; +above exact / -overshoot below)')
    print('%-4s %-4s %-6s %-10s %-10s %-10s %-10s | %-9s %-9s %-9s' % (
        'sys', 'ns', 'k', 'E_plain', 'E_TC', 'E_sym', 'exact',
        'd_plain', 'd_TC', 'd_sym'))
    print('-' * 118)
    for name in ('He', 'Li'):
        ex = EXACT[name]
        for ns in NS_LIST:
            d = rep[name]['per_ns'][ns]
            row = [r for r in d['rows'] if abs(r['gamma'] - 1.0) < 1e-9][0]
            print('%-4s %-4d %-6.2f %-10.5f %-10.5f %-10.5f %-10.5f | %+9.2f %+9.2f %+9.2f'
                  % (name, ns, d['k'], d['E_plain'], row['E_TC'], row['E_sym'], ex,
                     (d['E_plain'] - ex) * 1e3, row['dTC_exact'] * 1e3,
                     row['dsym_exact'] * 1e3))
    print('=' * 118)
    print('gamma-stationarity (dE_TC/dgamma=0) and exact-crossing gammas:')
    for name in ('He', 'Li'):
        for ns in NS_LIST:
            d = rep[name]['per_ns'][ns]
            print('  %s ns=%d: monotone=%s gamma*=%s  TC crosses exact @g=%s  sym crosses exact @g=%s'
                  % (name, ns, d['monotone'], d['gamma_star'],
                     d['gamma_TC_crosses_exact'], d['gamma_sym_crosses_exact']))
    print('basis trend of the TC/sym correction at gamma=1.0 and 1.5 (mHa):')
    for name in ('He', 'Li'):
        for g in (1.0, 1.5):
            vals = []
            for ns in NS_LIST:
                r = [x for x in rep[name]['per_ns'][ns]['rows']
                     if abs(x['gamma'] - g) < 1e-9][0]
                vals.append((ns, r['d_plain_TC'] * 1e3, r['d_plain_sym'] * 1e3))
            s = '  '.join('ns%d:TC%+.1f/sym%+.1f' % v for v in vals)
            print('  %s g=%.1f  (E_plain-E_TC / E_plain-E_sym):  %s' % (name, g, s))


if __name__ == '__main__':
    main()
