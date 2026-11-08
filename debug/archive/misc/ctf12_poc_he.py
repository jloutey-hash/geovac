"""CT-F12 PoC on a Coulomb-Sturmian basis for He (singlet ground state).

Compares plain FCI (bare 1/r12) vs canonical/Hermitian transcorrelated F12
(bare H + two-body geminal dressing w(r12), Kato-cusp fixed amplitude) at
matched basis size. Reuses the validated SturmianCI angular/Loewdin/Slater-
Condon machinery; only the two-electron radial kernel is swapped.

  w(r) = (1-e^{-g r})/r + (g/2) e^{-g r} - (1/4) e^{-2 g r}   (see ctf12_validate.py)

Basis: Coulomb-Sturmians S_{nl}(r)=hydrogenic_radial(r,n,l,Z=n*k) (shared decay
e^{-kr}); size set by max_n (l runs 0..n-1). k = single scale parameter.
"""
import numpy as np
import json, time
from scipy.special import eval_legendre
from scipy.linalg import eigh
from geovac.sturmian_solver import (hydrogenic_radial, _ck_coefficient,
                                     SturmianCI)

EXACT_HE = -2.903724377   # Ha, non-relativistic infinite-mass He ground state


# ---------------------------------------------------------------- kernels
def w_kernel(r, g):
    r = np.asarray(r, float)
    small = r < 1e-8
    rr = np.where(small, 1.0, r)
    out = -np.expm1(-g * rr) / rr + (g / 2) * np.exp(-g * rr) - 0.25 * np.exp(-2 * g * rr)
    out = np.where(small, 1.5 * g - 0.25, out)
    return out


# ---------------------------------------------------------------- grid + radial
def make_grid(k, max_n, Ng=1400):
    """Quadratic radial grid r = r_max * t^2 (clusters near 0 for the cusp);
    nonuniform trapezoid weights. ~5-10 uHa on (1s1s|1s1s) at Ng>=1200."""
    r_max = 45.0 / k + 6.0 * max_n
    t = np.linspace(0.0, 1.0, Ng)
    r = r_max * t ** 2
    r[0] = 1e-9
    wr = np.zeros(Ng)
    wr[1:-1] = (r[2:] - r[:-2]) / 2.0
    wr[0] = (r[1] - r[0]) / 2.0
    wr[-1] = (r[-1] - r[-2]) / 2.0
    return r, wr


def radial_products(states, k, r):
    uniq = sorted(set((n, l) for n, l, m in states))
    return {(n, l): hydrogenic_radial(r, n, l, n * k) for (n, l) in uniq}


def multipole_kernels(r, kmax, kind, g=None, nx=96):
    """g_k[i,j] on the (r_i,r_j) grid for k=0..kmax.
       kind='coulomb' -> r<^k/r>^{k+1};  kind='ctf12' -> ((2k+1)/2) int P_k w dx."""
    R1, R2 = np.meshgrid(r, r, indexing='ij')
    ker = {}
    if kind == 'coulomb':
        rl = np.minimum(R1, R2)
        rg = np.maximum(R1, R2)
        for k in range(kmax + 1):
            ker[k] = rl ** k / rg ** (k + 1)
    else:
        xs, ws = np.polynomial.legendre.leggauss(nx)
        acc = {k: np.zeros_like(R1) for k in range(kmax + 1)}
        for x, wx in zip(xs, ws):
            r12 = np.sqrt(R1 * R1 + R2 * R2 - 2 * R1 * R2 * x)
            wv = w_kernel(r12, g)
            for k in range(kmax + 1):
                acc[k] += wx * eval_legendre(k, x) * wv
        for k in range(kmax + 1):
            ker[k] = (2 * k + 1) / 2.0 * acc[k]
    return ker


# ---------------------------------------------------------------- ERI builder
def build_eri_grid(states, k, r, wr, kind, g=None):
    n_sp = len(states)
    lmax = max(l for n, l, m in states)
    kmax = 2 * lmax
    Rtab = radial_products(states, k, r)
    ker = multipole_kernels(r, kmax, kind, g)
    r2w = r * r * wr

    ac_k_map = {}
    for a in range(n_sp):
        la, ma = states[a][1], states[a][2]
        for c in range(n_sp):
            lc, mc = states[c][1], states[c][2]
            lst = []
            for kk in range(0, la + lc + 1):
                if (la + lc + kk) % 2:
                    continue
                v = _ck_coefficient(la, ma, lc, mc, kk)
                if abs(v) > 1e-15:
                    lst.append((kk, v))
            if lst:
                ac_k_map[(a, c)] = lst

    rk_cache = {}

    def get_G(a, c, b, d, kk):
        na, la = states[a][0], states[a][1]
        nc, lc = states[c][0], states[c][1]
        nb, lb = states[b][0], states[b][1]
        nd, ld = states[d][0], states[d][1]
        key = (na, la, nc, lc, nb, lb, nd, ld, kk)
        if key in rk_cache:
            return rk_cache[key]
        f1 = Rtab[(na, la)] * Rtab[(nc, lc)] * r2w
        f2 = Rtab[(nb, lb)] * Rtab[(nd, ld)] * r2w
        G = float(f1 @ ker[kk] @ f2)
        rk_cache[key] = G
        return G

    eri = {}
    for (a, c), acl in ac_k_map.items():
        ma = states[a][2]
        mc = states[c][2]
        for (b, d), bdl in ac_k_map.items():
            mb = states[b][2]
            md = states[d][2]
            if ma + mb != mc + md:
                continue
            val = 0.0
            for kac, cac in acl:
                for kbd, cbd in bdl:
                    if kac != kbd:
                        continue
                    val += cac * cbd * get_G(a, c, b, d, kac)
            if abs(val) > 1e-14:
                eri[(a, b, c, d)] = val
    return eri


# ---------------------------------------------------------------- FCI driver
def solve(Z, max_n, k, kind, g=None, Ng=600):
    ci = SturmianCI(Z, 2, max_n)
    states = ci.states
    S = ci._build_overlap(k)
    h1 = ci._build_h1_sturmian(k, S)
    r, wr = make_grid(k, max_n, Ng)
    eri = build_eri_grid(states, k, r, wr, kind, g)
    h1o, erio = ci._lowdin_transform(S, h1, eri)
    Hfci = ci._build_fci_hamiltonian(h1o, erio)
    ev = eigh(Hfci, eigvals_only=True)
    return float(ev[0]), states


def n_orb_for(Z, max_n):
    return len(SturmianCI(Z, 2, max_n).states)


# ---------------------------------------------------------------- end-to-end geminal validation
def validate_geminal_1111(k, g):
    """(1s1s|w|1s1s) via grid assembly vs direct (r1,r2,x) quadrature (Laguerre)."""
    states = [(1, 0, 0)]
    r, wr = make_grid(k, 1, 900)
    eri = build_eri_grid(states, k, r, wr, 'ctf12', g)
    grid_val = eri[(0, 0, 0, 0)]

    def R1s(rr):
        return hydrogenic_radial(np.array([rr]), 1, 0, k)[0]

    def inner(r1, r2):
        xs, ws = np.polynomial.legendre.leggauss(200)
        r12 = np.sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * xs)
        return 0.5 * np.sum(ws * w_kernel(r12, g))

    from numpy.polynomial.laguerre import laggauss
    xa, wa = laggauss(80)
    r1g = xa / (2 * k)
    tot = 0.0
    for i, r1 in enumerate(r1g):
        f1 = R1s(r1) ** 2 * r1 * r1 * np.exp(2 * k * r1) * wa[i] / (2 * k)
        for j, r2 in enumerate(r1g):
            f2 = R1s(r2) ** 2 * r2 * r2 * np.exp(2 * k * r2) * wa[j] / (2 * k)
            tot += f1 * f2 * inner(r1, r2)
    return grid_val, tot


# ---------------------------------------------------------------- main sweep
if __name__ == '__main__':
    t0 = time.time()
    gv, dv = validate_geminal_1111(1.6875, 1.2)
    print(f"[validate] (1s1s|w|1s1s): grid={gv:.8f}  directquad={dv:.8f}  |d|={abs(gv - dv):.2e}\n")

    Z = 2
    results = {'system': 'He', 'exact': EXACT_HE,
               'validation_1111': {'grid': gv, 'direct': dv, 'absdiff': abs(gv - dv)},
               'runs': []}
    kscan = np.round(np.arange(1.4, 2.31, 0.1), 3)
    gammas = [1.0, 1.2, 1.5]
    print(f"{'max_n':>5} {'n_orb':>5} {'method':>8} {'best_par':>16} {'E':>12} {'err_mHa':>9}")
    for max_n in [1, 2, 3, 4, 5]:
        no = n_orb_for(Z, max_n)
        best = (1e9, None)
        for k in kscan:
            E, _ = solve(Z, max_n, float(k), 'coulomb')
            if E < best[0]:
                best = (E, {'k': float(k)})
        err = (best[0] - EXACT_HE) * 1000
        print(f"{max_n:5d} {no:5d} {'FCI':>8} {'k=%.2f' % best[1]['k']:>16} {best[0]:12.6f} {err:9.3f}")
        results['runs'].append({'max_n': max_n, 'n_orb': no, 'qubit_proxy': 2 * no,
                                'method': 'FCI', 'E': best[0], 'err_mHa': err, 'par': best[1]})
        bestct = (1e9, None)
        for g in gammas:
            for k in kscan:
                E, _ = solve(Z, max_n, float(k), 'ctf12', g)
                if E < bestct[0]:
                    bestct = (E, {'k': float(k), 'gamma': float(g)})
        errct = (bestct[0] - EXACT_HE) * 1000
        print(f"{max_n:5d} {no:5d} {'CT-F12':>8} {'k=%.2f g=%.1f' % (bestct[1]['k'], bestct[1]['gamma']):>16} {bestct[0]:12.6f} {errct:9.3f}")
        results['runs'].append({'max_n': max_n, 'n_orb': no, 'qubit_proxy': 2 * no,
                                'method': 'CT-F12', 'E': bestct[0], 'err_mHa': errct, 'par': bestct[1]})
    results['walltime_s'] = time.time() - t0
    with open('debug/data/ctf12_poc_he.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nwall {time.time() - t0:.1f}s -> debug/data/ctf12_poc_he.json")
