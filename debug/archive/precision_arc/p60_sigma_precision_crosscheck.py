"""Independent extended-precision cross-check of geovac/sturmian_sigma_basis.py (Paper 60 sigma basis).

Route (independent of the float64 Gauss quadrature):
  every sigma Coulomb-Sturmian on centre A (z=-R/2) or B (z=+R/2) is  poly(xi, eta) * exp(-k r_c),
  because r_c = h(xi +/- eta), z_c = h(xi*eta +/- 1) and r^l P_l(z/r) is a polynomial in (z, r^2).
  The volume element over 1/r_c cancels exactly, so every S/W entry is a finite sum of
      X_p(alpha) * Y_q(b),  X_p = int_1^inf xi^p e^{-alpha xi},  Y_q = int_{-1}^1 eta^q e^{-b eta},
  evaluated in closed form with mpmath at DPS digits (exact polynomial coefficients, no quadrature).
  Eigenvalues by mpmath eigsy at the same precision.  Third route spot check: mpmath 2-D quad.

Output: debug/data/p60_sigma_precision_crosscheck.json.  Diagnostic only (writes nothing outside debug/).
"""
import json
import sys
import time
from math import comb
from pathlib import Path

import mpmath as mp
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import geovac.sturmian_sigma_basis as fb  # noqa: E402

DPS = 60
mp.mp.dps = DPS


# ---------------------------------------------------------------- bivariate polynomials in (xi, eta)
def pmul(a, b):
    out = {}
    for (p1, q1), c1 in a.items():
        for (p2, q2), c2 in b.items():
            key = (p1 + p2, q1 + q2)
            out[key] = out.get(key, 0) + c1 * c2
    return out


def padd(a, b, s=1):
    out = dict(a)
    for kk, c in b.items():
        out[kk] = out.get(kk, 0) + s * c
    return out


def pscale(a, s):
    return {kk: c * s for kk, c in a.items()}


def ppow(a, e):
    out = {(0, 0): mp.mpf(1)}
    for _ in range(e):
        out = pmul(out, a)
    return out


def sigma_function(n, l, k, h, centre):
    """Polynomial part of chi_{n l 0} on centre 0 (A) or 1 (B); full function = poly * exp(-k r_c)."""
    sg = 1 if centre == 0 else -1
    r = {(1, 0): h, (0, 1): sg * h}               # r_c = h(xi +/- eta)
    z = {(1, 1): h, (0, 0): sg * h}               # z_c = h(xi*eta +/- 1)
    r2 = pmul(r, r)
    # solid harmonic r^l P_l(z/r) = 2^-l sum_j (-1)^j C(l,j) C(2l-2j,l) z^{l-2j} r^{2j}
    sol = {}
    for j in range(l // 2 + 1):
        c = mp.mpf((-1) ** j * comb(l, j) * comb(2 * l - 2 * j, l)) / mp.mpf(2) ** l
        sol = padd(sol, pscale(pmul(ppow(z, l - 2 * j), ppow(r2, j)), c))
    # Laguerre L_m^{(2l+1)}(2 k r)
    m, a = n - l - 1, 2 * l + 1
    lag = {}
    for i in range(m + 1):
        c = mp.mpf((-1) ** i * comb(m + a, m - i)) / mp.factorial(i) * (2 * k) ** i
        lag = padd(lag, pscale(ppow(r, i), c))
    nrm = mp.sqrt((2 * k) ** 3 * mp.factorial(n - l - 1) / (2 * n * mp.factorial(n + l)))
    nrm *= (2 * k) ** l * mp.sqrt(mp.mpf(2 * l + 1) / (4 * mp.pi))
    return pscale(pmul(sol, lag), nrm)


class Moments:
    def __init__(self, alpha):
        self.alpha = alpha
        self.X, self.Y = {}, {}

    def x(self, p):
        if p not in self.X:
            a = self.alpha
            self.X[p] = mp.factorial(p) * mp.e ** (-a) * mp.fsum(
                1 / (mp.factorial(j) * a ** (p - j + 1)) for j in range(p + 1))
        return self.X[p]

    def y(self, q, b):
        key = (q, b)
        if key not in self.Y:
            if b == 0:
                self.Y[key] = mp.mpf(2) / (q + 1) if q % 2 == 0 else mp.mpf(0)
            else:   # series in b: sum_m (-b)^m/m! * int eta^{q+m}
                tot, m, term = mp.mpf(0), 0, mp.mpf(1)
                while True:
                    if (q + m) % 2 == 0:
                        tot += term * 2 / (q + m + 1)
                    m += 1
                    term *= -b / m
                    if m > 20 and abs(term) < mp.mpf(10) ** (-DPS - 5):
                        break
                self.Y[key] = tot
        return self.Y[key]


def build_mp(N, k, R, Z=1):
    k, R = mp.mpf(k), mp.mpf(R)
    h, alpha = R / 2, k * R
    labels = fb.basis(N)
    funcs = [(sigma_function(n, l, k, h, 0), 0, n) for (n, l) in labels] + \
            [(sigma_function(n, l, k, h, 1), 1, n) for (n, l) in labels]
    M = Moments(alpha)
    wa = {(1, 0): 2 * mp.pi * h ** 2, (0, 1): -2 * mp.pi * h ** 2}   # dV / r_A
    wb = {(1, 0): 2 * mp.pi * h ** 2, (0, 1): 2 * mp.pi * h ** 2}    # dV / r_B

    def integ(poly, b):
        return mp.fsum(c * M.x(p) * M.y(q, b) for (p, q), c in poly.items())

    m = len(funcs)
    S, W = mp.matrix(m, m), mp.matrix(m, m)
    va_tab, vb_tab = {}, {}
    for i in range(m):
        fi, ci, _ = funcs[i]
        for j in range(i, m):
            fj, cj, _ = funcs[j]
            b = {(0, 0): alpha, (1, 1): -alpha}.get((ci, cj), mp.mpf(0))  # eta exponent coefficient
            prod = pmul(fi, fj)
            va = integ(pmul(prod, wa), b)
            vb = integ(pmul(prod, wb), b)
            va_tab[(i, j)] = va_tab[(j, i)] = va
            vb_tab[(i, j)] = vb_tab[(j, i)] = vb
    for i in range(m):
        for j in range(m):
            _, cj, nj = funcs[j]
            S[i, j] = (nj / k) * (va_tab[(i, j)] if cj == 0 else vb_tab[(i, j)])
            W[i, j] = Z * (va_tab[(i, j)] + vb_tab[(i, j)]) / k
    Ssym = (S + S.T) / 2
    return Ssym, W, labels, S


def eig_mp(A):
    ev = mp.eigsy(A, eigvals_only=True)
    return sorted(ev)


def gerade_block_mp(S, labels):
    m = len(labels)
    U = mp.matrix(2 * m, m)
    for i, (_, l) in enumerate(labels):
        U[i, i] = 1 / mp.sqrt(2)
        U[m + i, i] = (-1) ** l / mp.sqrt(2)
    return U.T * S * U


def top_generalized_mp(W, S):
    ev, V = mp.eigsy(S)
    m = S.rows
    X = mp.matrix(m, m)
    for j in range(m):
        for i in range(m):
            X[i, j] = V[i, j] / mp.sqrt(ev[j])
    Wsym = (W + W.T) / 2
    return max(mp.eigsy(X.T * Wsym * X, eigvals_only=True))


def galerkin_energy_mp(N, R, k0):
    def f(kk):
        S, W, _, _ = build_mp(N, kk, R)
        return top_generalized_mp(W, S) - kk
    kroot = mp.findroot(f, (mp.mpf(k0) - mp.mpf('1e-4'), mp.mpf(k0) + mp.mpf('1e-4')),
                        solver='secant', tol=mp.mpf(10) ** (-30))
    return -kroot ** 2 / 2, kroot


def to_np(A):
    return np.array([[float(A[i, j]) for j in range(A.cols)] for i in range(A.rows)])


def main():
    t0 = time.time()
    out = {"meta": {"dps": DPS, "route": "exact polynomial x closed-form prolate moments (mpmath)",
                    "float_route": f"geovac/sturmian_sigma_basis.py GL{fb.NX} x GLeg{fb.NE} float64",
                    "date": "2026-10-09"}}

    # --- 1. entry-level agreement + intra-block identity + third-route spot check
    Smp, Wmp, labels, _ = build_mp(6, 1, 3)
    Sf, Wf, _ = fb.build(6, 1.0, 3.0)
    m = len(labels)
    dS = np.max(np.abs(to_np(Smp) - Sf)); dW = np.max(np.abs(to_np(Wmp) - Wf))
    intra = max(abs(Smp[i, j] - (1 if i == j else 0)) for i in range(m) for j in range(m))
    out["entries_N6_kR3"] = {"max_abs_dS_float_vs_mp": float(dS), "max_abs_dW_float_vs_mp": float(dW),
                             "mp_intra_block_dev_from_identity": float(intra)}
    print("entries N=6 kR=3:", out["entries_N6_kR3"], flush=True)

    # third route: direct 2-D mpmath quad of one cross element S[A 4f, B 6h]-type pair
    mp.mp.dps = 25
    iA, jB = labels.index((4, 3)), m + labels.index((6, 5))
    k_, R_ = mp.mpf(1), mp.mpf(3); h_ = R_ / 2
    fA = sigma_function(4, 3, k_, h_, 0); fB = sigma_function(6, 5, k_, h_, 1)

    def ev_poly(P, x, e):
        return mp.fsum(c * x ** p * e ** q for (p, q), c in P.items())

    def integrand(x, e):   # <A|1/r_B|B> dV, radial/angle evaluated pointwise from scipy-free formula
        ra, rb = h_ * (x + e), h_ * (x - e)
        val = ev_poly(fA, x, e) * ev_poly(fB, x, e) * mp.e ** (-k_ * (ra + rb))
        return val * 2 * mp.pi * h_ ** 3 * (x ** 2 - e ** 2) / rb
    q3 = mp.quad(integrand, [1, 3, 10, mp.inf], [-1, 0, 1]) * 6 / k_
    mp.mp.dps = DPS
    out["third_route_quad"] = {"entry": "S[A(4,3), B(6,5)] kR=3", "mp_quad": float(q3),
                               "mp_closed": float(Smp[iA, jB]), "float64": float(Sf[iA, jB])}
    print("third route:", out["third_route_quad"], flush=True)

    # --- 2. metric statistics, every quoted number
    out["metric"] = []
    for s in (1, 3):
        for N in (2, 3, 4, 5, 6):
            Smp, _, labels, _ = build_mp(N, 1, s)
            ev = eig_mp(Smp)
            evg = eig_mp(gerade_block_mp(Smp, labels))
            mm = len(labels)
            deficit = 1 - mp.fsum(Smp[i, mm] ** 2 for i in range(mm))
            Sf, _, _ = fb.build(N, 1.0, float(s))
            evf = np.linalg.eigvalsh(Sf)
            Sg, _ = fb.gerade_ungerade(Sf, labels)
            rec = {"kR": s, "N": N, "fns_per_centre": mm,
                   "cond_mp": float(ev[-1] / ev[0]), "cond_float": float(evf[-1] / evf[0]),
                   "cond_gerade_mp": float(evg[-1] / evg[0]),
                   "cond_gerade_float": float(np.linalg.cond(Sg)),
                   "soft_1e-2_mp": int(sum(1 for e in ev if e < mp.mpf('0.01'))),
                   "soft_1e-2_float": int((evf < 1e-2).sum()),
                   "deficit_mp": float(deficit), "deficit_float": fb.bessel_deficit(N, float(s)),
                   "smallest8_mp": [mp.nstr(e, 12) for e in ev[:8]],
                   "smallest8_float": [float(e) for e in evf[:8]],
                   "lam_max_mp": float(ev[-1])}
            out["metric"].append(rec)
            print(f"kR={s} N={N}: cond mp {rec['cond_mp']:.4e} fl {rec['cond_float']:.4e} | "
                  f"ger mp {rec['cond_gerade_mp']:.4e} fl {rec['cond_gerade_float']:.4e} | "
                  f"soft {rec['soft_1e-2_mp']}/{rec['soft_1e-2_float']} | "
                  f"def {rec['deficit_mp']:.5f}/{rec['deficit_float']:.5f}", flush=True)

    # --- 3. H2+ R=2 Galerkin energies: mp (no tolerance needed) vs float at several CANON_TOL
    exact = -1.102634214
    out["h2plus_R2"] = []
    for N in (2, 3, 5):
        kf = np.sqrt(-2 * fb.galerkin_energy(N, 2.0))
        Emp, kmp = galerkin_energy_mp(N, 2, kf)
        rec = {"N": N, "fns_per_centre": len(fb.basis(N)), "E_mp": mp.nstr(Emp, 20),
               "err_mp": float(Emp - exact), "float_by_tol": {}}
        for tol in (0.0, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6):
            fb.CANON_TOL = tol
            try:
                rec["float_by_tol"][str(tol)] = fb.galerkin_energy(N, 2.0) - exact
            except Exception as exc:  # noqa: BLE001
                rec["float_by_tol"][str(tol)] = f"fail: {exc}"
        fb.CANON_TOL = 1e-12
        out["h2plus_R2"].append(rec)
        print(f"H2+ N={N}: mp err {rec['err_mp']:.6e}  float {rec['float_by_tol']}", flush=True)

    out["meta"]["wall_s"] = time.time() - t0
    p = Path(__file__).parent / "data" / "p60_sigma_precision_crosscheck.json"
    p.write_text(json.dumps(out, indent=1))
    print("wall", out["meta"]["wall_s"], "->", p)


if __name__ == "__main__":
    main()
