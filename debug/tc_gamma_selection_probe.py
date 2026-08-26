"""Stage-2: is there a PRINCIPLED, INTERNAL criterion that selects the accurate
geminal width gamma for the native non-Hermitian TC operator -- WITHOUT knowing
the exact energy (no external VMC/Jastrow optimization)?

Stage-1 (debug/tc_accuracy_stage1_he.py) established that the non-Hermitian TC
matches the R12-CI accuracy anchor (~2.6 mHa) only at a hand-picked gamma~0.9 on
He, is non-variationally fragile elsewhere (gamma=0.5 -> -18.9 mHa overshoot), and
E_TC(gamma) is monotone (no energy stationarity).  This probe asks: does ANY
internal gamma-picker land the accurate regime?

Criteria tested (each is a gamma-scan; report the gamma it SELECTS with NO knowledge
of exact, then the accuracy E_TC(gamma_selected) - exact):

  (1a) variance-min, aufbau ref :   argmin_g  sigma^2 = <Phi|(H~-E0)+(H~-E0)|Phi>,
       Phi = aufbau reference determinant, E0 = <Phi|H~|Phi> = H[ref,ref].
       (standard TC/VMC internal criterion: the gamma that makes the reference det
        closest to an eigenstate of the transcorrelated H~.)
  (1b) variance-min, dominant det : same, Phi = dominant determinant of the ground
        RIGHT eigenvector.
  (1c) projected-energy stationarity : gamma where dE_ref/dgamma = 0,
        E_ref(g) = H[ref,ref] (single-determinant / projected TC energy).
  (2a) non-normality-min : argmin_g  ||[H~+,H~]||_F / ||H~||_F^2.
  (2b) left-right consistency : argmin_g (1 - |<L0|R0>|) for the ground bi-orthogonal
        eigenpair  ( = argmin Bauer-Fike kappa_0 ).
  (3)  range/cusp shape condition : gamma fixed by geminal shape, not energy.
        The Slater geminal u'(0)=1/2 fixes the Kato cusp for ANY gamma, so the cusp
        determines nothing; the concrete representative shape-match is gamma = k
        (geminal decay = orbital decay).

CALIBRATION per system: ORACLE gamma = argmin |E_TC - exact|; for He also the
R12-CI-anchored "good gamma" ~0.9.

DECISION GATE: GO if >=1 internal criterion robustly selects a gamma giving E_TC
within ~2 mHa of exact on BOTH He and Li (and near the oracle gamma).  STOP if
every internal criterion misses.

Engines reused READ-ONLY (geovac/transcorrelated_sturmian.py, debug/ctf12_r12ci_he.py);
no engine/physics modified.
"""
import os, sys, json, time
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, '..'))
for p in (_HERE, _ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

from geovac import transcorrelated_sturmian as T   # tracked TC engine (read-only)
import ctf12_r12ci_he as R                          # He R12-CI anchor (read-only)

EXACT = {'He': -2.903724, 'Li': -7.47806}

# systems: name -> (ns, k, Z, n_elec, with_L3)
SYSTEMS = [
    ('He',     dict(ns=3, k=1.9, Z=2, ne=2, wl3=False)),
    ('Li_ns3', dict(ns=3, k=1.6, Z=3, ne=3, wl3=True)),   # as instructed (basis-limited)
    ('Li_ns5', dict(ns=5, k=1.6, Z=3, ne=3, wl3=True)),   # fair criterion test (accurate regime exists)
]
GAMMAS = [0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10,
          1.20, 1.40, 1.60, 2.00, 2.50, 3.00]
NG, NX = 600, 96


def _refidx(sysobj):
    return sysobj.didx[tuple(sorted(sysobj.ref_occ))]


def scan_system(ns, k, Z, ne, wl3):
    """gamma-scan; per gamma collect E_TC, projected E_ref, variance (aufbau + dominant),
    non-normality, ground left-right inconsistency."""
    rows = []
    for g in GAMMAS:
        s = T.build_atomic_system(ns, k, g, Z, ne, Ng=NG, nx=NX, with_L3=wl3)
        variant = s.xtc_full() if s.v2 is not None else s.tc2()
        H = T.build_fci_matrix(s, variant)
        E_TC, im = T.ground(H, hermitian=False)

        ref = _refidx(s)
        E_ref = float(H[ref, ref])                      # projected / single-det TC energy
        col = H[:, ref].copy()                          # (H~ |Phi>) column
        # variance around the projected mean E0 = H[ref,ref] (standard local-E variance)
        off = col.copy(); off[ref] = 0.0
        var_ref = float(np.sum(off * off))              # sum_{J!=ref} |H[J,ref]|^2

        # dominant determinant of the ground RIGHT eigenvector
        w, VL, VR = T._eig_lr(H)
        order = np.argsort(w.real)
        i0 = order[0]
        v0 = VR[:, i0].real
        dom = int(np.argmax(np.abs(v0)))
        cold = H[:, dom].copy()
        E0d = float(H[dom, dom])
        offd = cold.copy(); offd[dom] = 0.0
        var_dom = float(np.sum(offd * offd))

        # ground bi-orthogonal left-right consistency
        l0 = VL[:, i0]; r0 = VR[:, i0]
        ov = abs(np.vdot(l0, r0)) / (np.linalg.norm(l0) * np.linalg.norm(r0))
        lr_incons = float(1.0 - ov)                     # 0 = normal, ->1 = degenerate biorthogonal

        # non-normality  ||[H+,H]||_F / ||H||_F^2
        Hd = H.conj().T
        comm = Hd @ H - H @ Hd
        nonnorm = float(np.linalg.norm(comm, 'fro') / np.linalg.norm(H, 'fro') ** 2)

        rows.append(dict(gamma=g, E_TC=E_TC, im=im, E_ref=E_ref,
                         var_ref=var_ref, var_dom=var_dom,
                         nonnorm=nonnorm, lr_incons=lr_incons,
                         dom_is_ref=int(dom == ref)))
    return rows


def _parabolic_min(gs, ys):
    """discrete argmin + parabolic refinement (returns refined gamma, clamped to grid)."""
    gs = np.asarray(gs, float); ys = np.asarray(ys, float)
    i = int(np.argmin(ys))
    if 0 < i < len(gs) - 1:
        g0, g1, g2 = gs[i-1], gs[i], gs[i+1]
        y0, y1, y2 = ys[i-1], ys[i], ys[i+1]
        denom = (y0 - 2*y1 + y2)
        if abs(denom) > 1e-300:
            # vertex of parabola through 3 pts (non-uniform ok via Lagrange derivative)
            # use uniform-ish estimate on index then map — grid is non-uniform so
            # fall back to fit in g directly:
            A = np.array([[g0*g0, g0, 1], [g1*g1, g1, 1], [g2*g2, g2, 1]], float)
            a, b, c = np.linalg.solve(A, np.array([y0, y1, y2]))
            if a > 0:
                gv = -b / (2*a)
                if g0 <= gv <= g2:
                    return float(gv), 'interior'
    edge = 'left-edge' if i == 0 else ('right-edge' if i == len(gs)-1 else 'interior')
    return float(gs[i]), edge


def _stationary(gs, ys):
    """gamma where dy/dgamma = 0 (first sign change of gradient); else monotone."""
    gs = np.asarray(gs, float); ys = np.asarray(ys, float)
    d = np.gradient(ys, gs)
    sg = np.sign(d)
    x = np.where(np.diff(sg) != 0)[0]
    if x.size:
        i = x[0]
        g0, g1 = gs[i], gs[i+1]; d0, d1 = d[i], d[i+1]
        gstar = float(g0 - d0 * (g1 - g0) / (d1 - d0)) if d1 != d0 else float(g0)
        return gstar, 'stationary'
    return None, 'monotone(%s)' % ('increasing' if d.mean() > 0 else 'decreasing')


def _etc_at(ns, k, Z, ne, wl3, g):
    s = T.build_atomic_system(ns, k, g, Z, ne, Ng=NG, nx=NX, with_L3=wl3)
    variant = s.xtc_full() if s.v2 is not None else s.tc2()
    H = T.build_fci_matrix(s, variant)
    E, im = T.ground(H, hermitian=False)
    return E


def he_r12_anchor(k):
    """R12-CI-anchored 'good gamma' for He = gamma where E_TC best matches the R12-CI
    variational anchor (the accuracy the correlated method actually reaches)."""
    # R12-CI variational optimum gamma (its own lowest energy above exact)
    best = (1e9, None)
    for g in GAMMAS:
        Er = R.he_energy(3, k, g, n_gem=1, gem_refs=[1], Ng=NG, nx=NX)[0]
        if Er > EXACT['He'] - 5e-3 and Er < best[0]:
            best = (Er, g)
    return best  # (E_R12, gamma_R12opt)


def main():
    t0 = time.time()
    report = dict(params=dict(gammas=GAMMAS, Ng=NG, nx=NX, exact=EXACT),
                  systems={})
    table = []   # flat rows for the final {system x criterion} table

    for name, cfg in SYSTEMS:
        ns, k, Z, ne, wl3 = cfg['ns'], cfg['k'], cfg['Z'], cfg['ne'], cfg['wl3']
        ex = EXACT['He' if name == 'He' else 'Li']
        print('=== %s  (ns=%d k=%.2f Z=%d ne=%d L3=%s) exact=%.6f ===' %
              (name, ns, k, Z, ne, wl3, ex))
        rows = scan_system(ns, k, Z, ne, wl3)
        gs = [r['gamma'] for r in rows]
        etc = np.array([r['E_TC'] for r in rows])

        # oracle
        d = np.abs(etc - ex)
        io = int(np.argmin(d))
        g_oracle = gs[io]
        print('  ORACLE gamma=%.2f  |E_TC-exact|=%.2f mHa' % (g_oracle, d[io]*1e3))

        # criteria selections
        sel = {}
        sel['1a_var_ref'] = _parabolic_min(gs, [r['var_ref'] for r in rows])
        sel['1b_var_dom'] = _parabolic_min(gs, [r['var_dom'] for r in rows])
        sel['1c_Eref_stat'] = _stationary(gs, [r['E_ref'] for r in rows])
        sel['2a_nonnorm'] = _parabolic_min(gs, [r['nonnorm'] for r in rows])
        sel['2b_lr_incons'] = _parabolic_min(gs, [r['lr_incons'] for r in rows])
        sel['3_shape_g=k'] = (float(k), 'shape')   # gamma = k

        crit_out = {}
        for cname, (gsel, kind) in sel.items():
            if gsel is None:
                E_at = None; dmH = None
            else:
                gcl = min(max(gsel, GAMMAS[0]), GAMMAS[-1])
                E_at = _etc_at(ns, k, Z, ne, wl3, gcl)
                dmH = (E_at - ex) * 1e3
            crit_out[cname] = dict(gamma_sel=gsel, kind=kind,
                                   E_TC=E_at, dTC_mHa=dmH)
            print('  %-14s gamma_sel=%s (%s)  E_TC-exact=%s mHa' %
                  (cname, ('%.3f' % gsel) if gsel is not None else 'None',
                   kind, ('%+.2f' % dmH) if dmH is not None else 'n/a'))
            table.append(dict(system=name, criterion=cname,
                              gamma_sel=(None if gsel is None else round(gsel, 3)),
                              dTC_mHa=(None if dmH is None else round(dmH, 2)),
                              g_oracle=g_oracle,
                              oracle_dTC_mHa=round(float(d[io]*1e3), 2)))

        extra = {}
        if name == 'He':
            Er, g_r12 = he_r12_anchor(k)
            extra['R12_anchor'] = dict(gamma_R12opt=g_r12, E_R12=Er,
                                       dR12_mHa=(Er-ex)*1e3)
            print('  R12-CI anchor: variational-opt gamma=%.2f  E=%.6f (%+.2f mHa)'
                  % (g_r12, Er, (Er-ex)*1e3))

        report['systems'][name] = dict(cfg=cfg, exact=ex, rows=rows,
                                       g_oracle=g_oracle,
                                       oracle_dTC_mHa=float(d[io]*1e3),
                                       selections={k2: dict(gamma_sel=v[0], kind=v[1])
                                                   for k2, v in sel.items()},
                                       criteria=crit_out, extra=extra)
        print()

    report['flat_table'] = table
    report['walltime_s'] = time.time() - t0

    os.makedirs(os.path.join(_HERE, 'data'), exist_ok=True)
    outp = os.path.join(_HERE, 'data', 'tc_gamma_selection.json')
    with open(outp, 'w') as f:
        json.dump(report, f, indent=2, default=float)

    # ---- final compact table ----
    print('=' * 92)
    print('%-8s %-14s %8s %10s %8s %10s' %
          ('system', 'criterion', 'gamma', 'dTC(mHa)', 'oracle_g', 'oracle_dTC'))
    print('-' * 92)
    for row in table:
        print('%-8s %-14s %8s %10s %8.2f %10.2f' %
              (row['system'], row['criterion'],
               ('%.3f' % row['gamma_sel']) if row['gamma_sel'] is not None else 'None',
               ('%+.2f' % row['dTC_mHa']) if row['dTC_mHa'] is not None else 'n/a',
               row['g_oracle'], row['oracle_dTC_mHa']))
    print('=' * 92)
    print('wrote', outp, ' (%.0fs)' % report['walltime_s'])


if __name__ == '__main__':
    main()
