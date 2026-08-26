"""
xTC quantum-algorithm cost measurement on GeoVac's non-Hermitian effective
2-body operators.
==========================================================================

Task: the classical xTC PoC (sprint_xtc_{poc_li,pblock}_memo.md) diagonalized the
effective operator with scipy.linalg.eig but never measured the operator
properties that drive quantum-algorithm cost for a NON-HERMITIAN object.  Here we
measure them on the ACTUAL operators, reusing the validated engines READ-ONLY:

    xtc_poc_li.py            (s-only; builds the convective K term  -> He, Li)
    xtc_poc_li_pinclusive.py (s+p;   D + xTC-L3 only, NO K)         -> C (p-ref)

Full non-Hermitian effective 2-body H~ = D (Herm) + K (non-Herm convective)
                                          + xTC-contracted-L3.

Structural fact established first (probe): the xTC 3-body->2-body contraction v2
is HERMITIAN to machine precision; D (asym_w) is Hermitian; the ONLY non-Hermitian
term is the generic 2-body convective K (asym_K, present in ANY transcorrelated H).

Per system we measure:
  (1) spectrum of the FCI H~ (scipy.linalg.eig): all-real? Im spread? GS real?
  (2) departure from normality  ||H^H H - H H^H||_F / ||H||_F^2
  (3) eigenvector condition number kappa_V = cond(V), full matrix AND restricted
      to the low-lying (ground + few excited) invariant subspace  [LOAD-BEARING]
      + per-eigenvalue Bauer-Fike condition numbers  kappa_i = 1/|<u_i|v_i>|
  (4) spectral gap  Delta = Re(E1)-Re(E0)
  (5) non-Hermitian JW LCU 1-norm  lambda_nonHerm  vs plain Hermitian lambda_Herm
  (6) escape hatch: symmetrized  H~_sym = (H~+H~^H)/2  ground vs true TC ground:
      fraction of the cusp win retained
  (7) fold  H~^H H~ : condition number + gap vs H~'s
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

import xtc_poc_li as X
import xtc_poc_li_pinclusive as P
from xtc_pblock_engine import REFERENCES

TOL = 1e-9


# ---------------------------------------------------------------------------
# matrix-level cost metrics for a (possibly non-Hermitian) dense FCI matrix
# ---------------------------------------------------------------------------
def measure_matrix(H, m_low=4, deg_tol=2e-3):
    """deg_tol (absolute Ha): eigenvalues within deg_tol of the ground are treated
    as one (near-)degenerate manifold; the physical gap is measured to the first
    level beyond it.  (The xTC single-determinant reference is spin-broken, so it
    artificially splits the exactly-degenerate ground spin/term manifold by
    ~0.1 mHa; deg_tol merges that spurious split.)"""
    H = np.asarray(H)
    n = H.shape[0]
    m_low = min(m_low, n)
    Hd = H.conj().T
    # (2) departure from normality
    comm = Hd @ H - H @ Hd
    fH2 = np.linalg.norm(H, 'fro') ** 2
    nonnorm = float(np.linalg.norm(comm, 'fro') / fH2) if fH2 > 0 else 0.0

    # (1) spectrum + (3) eigenvectors (right + left, matched ordering)
    w, VL, VR = sla.eig(H, left=True, right=True)
    order = np.argsort(w.real)
    w = w[order]; VR = VR[:, order]; VL = VL[:, order]
    im_spread = float(np.max(np.abs(w.imag)))
    all_real = bool(im_spread < 1e-6 * (1 + np.max(np.abs(w.real))))
    E0 = float(w[0].real); im0 = float(abs(w[0].imag))

    # (4) gap: raw, physical (deg-merged), + ground-manifold splitting
    re = w.real
    gap_raw = float(re[1] - re[0]) if n > 1 else 0.0
    gap_merged = 0.0
    for i in range(1, n):
        if re[i] - re[0] > deg_tol:
            gap_merged = float(re[i] - re[0]); break
    in_mani = re[re - re[0] <= deg_tol]
    gs_split = float(in_mani.max() - in_mani.min())
    gs_degeneracy = int(in_mani.size)

    # (3a) full eigenvector condition number
    kV_full = float(np.linalg.cond(VR))
    # (3b) restricted to low-lying subspace: cond of the n x m_low right-evec block
    Vlow = VR[:, :m_low]
    kV_low = float(np.linalg.cond(Vlow))
    # (3c) per-eigenvalue Bauer-Fike condition numbers for the lowest m_low
    kappa_eig = []
    for i in range(m_low):
        r = VR[:, i]; l = VL[:, i]
        nr = np.linalg.norm(r); nl = np.linalg.norm(l)
        ov = abs(np.vdot(l, r))
        kappa_eig.append(float(nr * nl / ov) if ov > 0 else float('inf'))
    max_kappa_eig_low = float(max(kappa_eig))

    return dict(n=int(n), nonnormality=nonnorm, im_spread=im_spread,
                all_real=all_real, E0=E0, im0=im0,
                gap_raw=gap_raw, gap_merged=gap_merged,
                gs_split=gs_split, gs_degeneracy=gs_degeneracy,
                kV_full=kV_full, kV_low=kV_low,
                kappa_eig_low=kappa_eig, max_kappa_eig_low=max_kappa_eig_low)


def sym_ground(H):
    """Hermitian ground energy of (H + H^H)/2."""
    Hs = 0.5 * (H + H.conj().T)
    wr = sla.eigh(Hs, eigvals_only=True)
    return float(wr[0])


def fold_stats(H, tol_deg=1e-6):
    """Fold G = H^H H (Hermitian PSD): condition number + gap (two smallest evals)."""
    G = H.conj().T @ H
    ev = sla.eigh(G, eigvals_only=True)
    ev = np.sort(ev.real)
    lo = ev[ev > 1e-14]
    condG = float(ev[-1] / lo[0]) if lo.size else float('inf')
    gapG = 0.0
    for i in range(1, ev.size):
        if ev[i] - ev[0] > tol_deg * (1 + abs(ev[0])):
            gapG = float(ev[i] - ev[0]); break
    return dict(cond=condG, gap=gapG, eval_min=float(ev[0]), eval_max=float(ev[-1]))


# ---------------------------------------------------------------------------
# JW LCU 1-norm of a second-quantized operator (h1 + antisymmetrized 2-body)
#   H = sum h[p,q] a+p a_q + (1/4) sum asym[p,q,r,s] a+p a+q a_s a_r  (+ v0)
# complex coefficients allowed (non-Hermitian friendly)
# ---------------------------------------------------------------------------
def lcu_lambda(hso, asym, nso, tol=1e-9):
    import openfermion as of
    ferm = of.FermionOperator()
    for p in range(nso):
        for q in range(nso):
            c = hso[p, q]
            if abs(c) > tol:
                ferm += of.FermionOperator(((p, 1), (q, 0)), complex(c))
    for p in range(nso):
        for q in range(nso):
            for r in range(nso):
                for s in range(nso):
                    c = 0.25 * asym[p, q, r, s]
                    if abs(c) > tol:
                        ferm += of.FermionOperator(((p, 1), (q, 1), (s, 0), (r, 0)),
                                                   complex(c))
    qop = of.jordan_wigner(ferm)
    qop.compress(1e-9)
    terms = {t: c for t, c in qop.terms.items() if t != ()}
    lam = float(sum(abs(c) for c in terms.values()))
    max_imag = float(max((abs(c.imag) for c in terms.values()), default=0.0))
    return dict(n_pauli=len(terms), lam=lam, max_imag_coeff=max_imag)


# ---------------------------------------------------------------------------
# hermiticity deviation of an antisym tensor <pq||rs>: max|A - A^dag|
# ---------------------------------------------------------------------------
def herm_dev(A):
    return float(np.max(np.abs(A - np.conjugate(np.transpose(A, (2, 3, 0, 1))))))


def dagger_asym(A):
    """Hermitian conjugate of an antisym 2-body tensor <pq||rs> -> <rs||pq>*."""
    return np.conjugate(np.transpose(A, (2, 3, 0, 1)))


# ===========================================================================
# System assembly (s-only engine, exposes ALL pieces incl. convective K)
# ===========================================================================
def build_sonly(ns, k, gamma, Z, n_elec, Ng, nx, with_L3):
    r, wr = X.make_grid(k, Ng=Ng)
    S, h1s, Rtab, W = X.build_one_body(ns, r, wr, k, Z)
    Km = X.build_kernels(r, gamma, nx=nx)
    eri_coul, eri_w, eri_K = X.two_body(ns, Rtab, W, Km)
    V3 = X.three_body(ns, Rtab, W, Km) if with_L3 else None
    Xm = X.lowdin(S)
    h1o = X.transform_1(h1s, Xm)
    eri_coul_o = X.transform_2(eri_coul, Xm)
    eri_w_o = X.transform_2(eri_w, Xm)
    eri_K_o = X.transform_2(eri_K, Xm)
    V3o = X.transform_3(V3, Xm) if with_L3 else None
    nso = 2 * ns
    dets, didx = X.make_dets(nso, n_elec)
    hso = X.h_spin(h1o, nso)
    asym_coul = X.asym_from_phys(eri_coul_o, nso)
    asym_w = X.asym_from_phys(eri_w_o, nso)
    asym_K = X.asym_from_phys(eri_K_o, nso)
    # aufbau reference occupation
    diag_e = np.diag(h1o); order = np.argsort(diag_e)
    if n_elec == 2:
        o0 = order[0]; ref_occ = (2 * o0, 2 * o0 + 1)
    else:
        o0, o1 = order[0], order[1]; ref_occ = (2 * o0, 2 * o0 + 1, 2 * o1)
    v2 = v1 = None; v0 = 0.0
    if with_L3:
        v2, v1, v0 = X.xtc_contract(V3o, nso, ref_occ)
    return dict(nso=nso, dets=dets, didx=didx, hso=hso, asym_coul=asym_coul,
                asym_w=asym_w, asym_K=asym_K, v2=v2, v1=v1, v0=v0,
                ref_occ=ref_occ, ndet=len(dets))


def build_pblock(name, k, gamma, Ng, nx):
    spec = REFERENCES[name]
    o = P.assemble(max_n=2, k=k, gamma=gamma, Z=spec['Z'], n_elec=spec['n_elec'],
                   Ng=Ng, nx=nx, want_exact3=True, ref_occ=spec['ref_occ'])
    a = o['_arr']
    nso = a['nso']
    dets, didx = X.make_dets(nso, spec['n_elec'])
    return dict(nso=nso, dets=dets, didx=didx, hso=a['hso'],
                asym_coul=a['asym_coul'], asym_w=a['asym_w'],
                v2=a['v2'], v1=a['v1'], v0=a['v0'],
                asym_x=a['asym_x'], hso_x=a['hso_x'],
                E_plain=o['E_plain'], E_xTC=o['E_xTC'], E_TC2=o['E_TC2'],
                E_exactTC=o['E_exactTC'], ndet=o['ndet'])


def measure_variant(bits, hso, asym, label, is_xTC=False, m_low=4, plain_lambda=None,
                    E_plain=None, want_sym=False, want_fold=False):
    nso = bits['nso']; dets = bits['dets']; didx = bits['didx']
    v0 = bits.get('v0', 0.0) if is_xTC else 0.0
    H = X.build_H(dets, didx, hso, asym, nso, v0=v0)
    mm = measure_matrix(H, m_low=m_low)
    lam = lcu_lambda(hso, asym, nso)
    # tensor-level 1-norms (memo-style op_metrics: sum|h1| , sum|<pq||rs>|)
    l1_1b = float(np.sum(np.abs(hso)))
    l1_2b = float(np.sum(np.abs(asym)))
    rec = dict(label=label, is_xTC=is_xTC, hdev_asym=herm_dev(asym), **mm, lcu=lam,
               l1_1body=l1_1b, l1_2body=l1_2b)
    if plain_lambda is not None:
        rec['lam_ratio_vs_plain'] = float(lam['lam'] / plain_lambda)
    if want_sym and E_plain is not None:
        E_sym = sym_ground(H)
        win_tot = E_plain - mm['E0']
        win_sym = E_plain - E_sym
        rec['E_sym'] = E_sym
        rec['cusp_win_total_mHa'] = float(win_tot * 1e3)
        rec['cusp_win_sym_mHa'] = float(win_sym * 1e3)
        rec['sym_win_retained'] = float(win_sym / win_tot) if abs(win_tot) > 1e-12 else None
    if want_fold:
        rec['fold'] = fold_stats(H)
    return rec


def run_He(Ng=800, nx=128):
    b = build_sonly(ns=3, k=1.7, gamma=1.0, Z=2, n_elec=2, Ng=Ng, nx=nx, with_L3=False)
    plain = measure_variant(b, b['hso'], b['asym_coul'], 'plain(Coulomb)')
    lamp = plain['lcu']['lam']; plain['lam_ratio_vs_plain'] = 1.0
    E_plain = plain['E0']
    donly = measure_variant(b, b['hso'], b['asym_w'], 'D only (Herm)', plain_lambda=lamp)
    asym_TC2 = b['asym_w'] + b['asym_K']
    tc2 = measure_variant(b, b['hso'], asym_TC2, 'D+K (2body TC)', plain_lambda=lamp,
                          E_plain=E_plain, want_sym=True, want_fold=True)
    asym_sym = 0.5 * (asym_TC2 + dagger_asym(asym_TC2))
    symv = measure_variant(b, b['hso'], asym_sym, 'sym(D+K) [Herm]', plain_lambda=lamp)
    return dict(system='He', basis='s-only (1s,2s,3s)', ndet=b['ndet'], nso=b['nso'],
                E_plain=E_plain, primary='D+K (2body TC)',
                asymK_hdev=herm_dev(b['asym_K']), asymw_hdev=herm_dev(b['asym_w']),
                variants=[plain, donly, tc2, symv])


def run_Li(Ng=800, nx=128):
    b = build_sonly(ns=3, k=1.5, gamma=1.0, Z=3, n_elec=3, Ng=Ng, nx=nx, with_L3=True)
    plain = measure_variant(b, b['hso'], b['asym_coul'], 'plain(Coulomb)')
    lamp = plain['lcu']['lam']; plain['lam_ratio_vs_plain'] = 1.0
    E_plain = plain['E0']
    hso_x = b['hso'] + b['v1']
    asym_DL3 = b['asym_w'] + b['v2']
    dL3 = measure_variant(b, hso_x, asym_DL3, 'D+xTC-L3 (noK)', is_xTC=True, plain_lambda=lamp)
    asym_full = b['asym_w'] + b['asym_K'] + b['v2']
    full = measure_variant(b, hso_x, asym_full, 'D+K+xTC-L3 (FULL)', is_xTC=True,
                           plain_lambda=lamp, E_plain=E_plain, want_sym=True, want_fold=True)
    asym_TC2 = b['asym_w'] + b['asym_K']
    tc2 = measure_variant(b, b['hso'], asym_TC2, 'D+K (2body TC)', plain_lambda=lamp)
    asym_sym = 0.5 * (asym_full + dagger_asym(asym_full))
    symv = measure_variant(b, hso_x, asym_sym, 'sym(FULL) [Herm]', is_xTC=True,
                           plain_lambda=lamp)
    return dict(system='Li', basis='s-only (1s,2s,3s)', ndet=b['ndet'], nso=b['nso'],
                E_plain=E_plain, primary='D+K+xTC-L3 (FULL)',
                v2_hdev=herm_dev(b['v2']), asymK_hdev=herm_dev(b['asym_K']),
                asymw_hdev=herm_dev(b['asym_w']),
                variants=[plain, tc2, dL3, full, symv])


def run_C(Ng=800, nx=160):
    b = build_pblock('C_3P', k=2.0, gamma=1.0, Ng=Ng, nx=nx)
    plain = measure_variant(b, b['hso'], b['asym_coul'], 'plain(Coulomb)')
    lamp = plain['lcu']['lam']; plain['lam_ratio_vs_plain'] = 1.0
    E_plain = plain['E0']
    full = measure_variant(b, b['hso_x'], b['asym_x'], 'D+xTC-L3 (noK, s+p)', is_xTC=True,
                           plain_lambda=lamp, E_plain=E_plain, want_sym=True, want_fold=True)
    return dict(system='C', basis='s+p (1s,2s,2p) [K omitted: not in validated engine]',
                ndet=b['ndet'], nso=b['nso'], E_plain=E_plain,
                primary='D+xTC-L3 (noK, s+p)',
                v2_hdev=herm_dev(b['v2']), variants=[plain, full])


def main():
    t0 = time.time()
    rep = dict(params=dict(gamma=1.0,
                           He='s-only ns=3 k=1.7 Z=2 2e',
                           Li='s-only ns=3 k=1.5 Z=3 3e',
                           C='s+p max_n=2 k=2.0 Z=6 6e (C_3P Hund ref)'),
               note=('Non-Hermiticity source: only the 2-body convective K '
                     '(asym_K); D and the xTC-L3 contraction v2 are Hermitian. '
                     'The p-block engine does not build K, so C is measured '
                     'without K (its xTC operator is Hermitian).'))
    print('=== He ==='); rep['He'] = run_He(); print('  done %.0fs' % (time.time() - t0))
    print('=== Li ==='); rep['Li'] = run_Li(); print('  done %.0fs' % (time.time() - t0))
    print('=== C  ==='); rep['C'] = run_C();  print('  done %.0fs' % (time.time() - t0))
    rep['walltime_s'] = time.time() - t0

    os.makedirs(os.path.join(_ROOT, 'debug', 'data'), exist_ok=True)
    outp = os.path.join(_ROOT, 'debug', 'data', 'xtc_quantum_cost.json')
    with open(outp, 'w') as f:
        json.dump(rep, f, indent=2, default=float)
    print('wrote', outp)

    print('\n' + '=' * 126)
    print('%-4s %-22s %-6s %-9s %-8s %-7s %-8s %-9s %-8s %-6s %-6s' %
          ('sys', 'variant', 'real?', 'nonnorm', 'kV_low', 'kV_ful', 'gap',
           'lam/lamH', 'symret', 'nPaul', 'lam'))
    print('-' * 126)
    for s in ('He', 'Li', 'C'):
        for v in rep[s]['variants']:
            print('%-4s %-22s %-6s %-9.2e %-8.3g %-7.3g %-8.4f %-9s %-8s %-6d %-6.2f' % (
                s, v['label'], 'Y' if v['all_real'] else 'N',
                v['nonnormality'], v['kV_low'], v['kV_full'], v['gap_merged'],
                ('%.3f' % v['lam_ratio_vs_plain']) if 'lam_ratio_vs_plain' in v else '-',
                ('%.2f' % v['sym_win_retained']) if v.get('sym_win_retained') is not None else '-',
                v['lcu']['n_pauli'], v['lcu']['lam']))
    print('=' * 126)
    print('gap column = physical gap (deg-merged); ground-manifold spurious split (gs_split):')
    for s in ('He', 'Li', 'C'):
        pv = [v for v in rep[s]['variants'] if v['label'] == rep[s]['primary']][0]
        print('  %-4s primary %-22s  E0=%.5f  gap=%.4f  gs_split=%.2e (deg=%d)  max_kappa_eig_low=%.3g' % (
            s, pv['label'], pv['E0'], pv['gap_merged'], pv['gs_split'],
            pv['gs_degeneracy'], pv['max_kappa_eig_low']))


if __name__ == '__main__':
    main()
