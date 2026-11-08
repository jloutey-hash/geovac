"""Q1 of the STO -> shared-k Coulomb-Sturmian evaluator scoping (diagnostic, debug only).

Avery & Avery (2017, STO appendix): with s = k r and rho = zeta / k,

    s^l e^{-rho s} = sum_{n>l} c_{n,l}(rho) S_{n,l}(s),
    c_{n,l}(rho) = 2^{l+1} n sqrt((n+l)!/(n (n-1-l)!)) (rho-1)^{n-1-l} / (rho+1)^{n+1+l},

S_{n,l} the L2-normalised Coulomb Sturmian (hydrogenic radial at Z = n k).  Higher-N STOs
r^{N-1} e^{-zeta r} (N > l+1) follow by (-d/drho)^{N-1-l}.  Everything depends on rho only.

Exact tools used (no quadrature):
  * Sturmian L2 Gram is tridiagonal: G_nn = 1, G_{n,n+1} = -1/2 sqrt((n-l)(n+l+1)/(n(n+1)));
  * Sturmian kinetic matrix: <S_m|T|S_n> = k^2 delta_mn - (k^2/2) G_mn  (from the Sturmian eq.);
  * so the L2 and kinetic (H1-seminorm) norms of the truncation residual are tail quadratic
    forms in the c's -- exact once the tail is summed far enough (n_max = 2000 here).

Checks: (a) Gram formula and (b) the derivative coefficients vs exact Fraction projection.
Outputs: debug/data/sto_sturmian_eval_convergence.json
"""
from __future__ import annotations

import json
import os
from fractions import Fraction as Fr
from math import comb, factorial

import mpmath as mp

mp.mp.dps = 40
NMAX = 2000
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                   "sto_sturmian_eval_convergence.json")


# ----------------------------------------------------------------- exact cross-checks

def _lag(m: int, a: int):
    return [Fr((-1) ** j * comb(m + a, m - j), factorial(j)) for j in range(m + 1)]


def _stpoly(n: int, l: int):
    """S_{n,l}(r) / (N e^{-r}) at k = 1 as {power: coeff}."""
    return {l + j: cj * 2 ** (l + j) for j, cj in enumerate(_lag(n - l - 1, 2 * l + 1))}


def _N2(n: int, l: int) -> Fr:
    return Fr(8 * factorial(n - l - 1), 2 * n * factorial(n + l))


def _inner(p1, p2, w: int, a: Fr) -> Fr:
    s = Fr(0)
    for i, ci in p1.items():
        for j, cj in p2.items():
            e = i + j + w
            s += ci * cj * Fr(factorial(e)) / a ** (e + 1)
    return s


def _mf(q: Fr):
    return mp.mpf(q.numerator) / q.denominator


def gram_off(n: int, l: int):
    return -mp.sqrt(mp.mpf((n - l) * (n + l + 1)) / (n * (n + 1))) / 2


def check_gram() -> float:
    worst = 0.0
    for l in range(3):
        for n in range(l + 1, l + 6):
            for m in range(l + 1, l + 6):
                g = _inner(_stpoly(n, l), _stpoly(m, l), 2, Fr(2))
                g = mp.mpf(g.numerator) / g.denominator * mp.sqrt(_mf(_N2(n, l) * _N2(m, l)))
                want = 1 if n == m else (gram_off(min(n, m), l) if abs(n - m) == 1 else 0)
                worst = max(worst, abs(g - want))
    return float(worst)


# ----------------------------------------------------------------- coefficients

def coeffs(Nr: int, l: int, rho, nmax: int = NMAX):
    """c_n for r^{Nr-1} e^{-rho r} (unnormalised, k = 1), n = l+1..nmax."""
    p = Nr - 1 - l
    rho = mp.mpf(rho)
    out = []
    for n in range(l + 1, nmax + 1):
        a, b = n - 1 - l, n + 1 + l
        A = 2 ** (l + 1) * n * mp.sqrt(mp.factorial(n + l) / (n * mp.factorial(n - 1 - l)))
        tot = mp.mpf(0)
        for j in range(p + 1):          # Leibniz on (rho-1)^a (rho+1)^{-b}
            if j > a:
                continue
            d1 = mp.factorial(a) / mp.factorial(a - j) * (rho - 1) ** (a - j)
            d2 = (-1) ** (p - j) * mp.rf(b, p - j) * (rho + 1) ** (-b - (p - j))
            tot += comb(p, j) * d1 * d2
        out.append((-1) ** p * A * tot)
    return out


def check_coeffs() -> float:
    worst = 0.0
    for (Nr, l) in ((2, 0), (3, 0), (3, 1), (4, 2)):
        for rho in (Fr(1, 3), Fr(5, 2)):
            c = coeffs(Nr, l, mp.mpf(rho.numerator) / rho.denominator, nmax=l + 6)
            f = {Nr - 1: Fr(1)}
            for i, n in enumerate(range(l + 1, l + 7)):
                pr = _inner(_stpoly(n, l), f, 1, 1 + rho) * n
                pr = mp.mpf(pr.numerator) / pr.denominator * mp.sqrt(_mf(_N2(n, l)))
                worst = max(worst, abs(pr - c[i]) / (abs(pr) + 1e-30))
    return float(worst)


# ----------------------------------------------------------------- residual norms

def _nmax(l: int, rho) -> int:
    t = abs(float(rho) - 1) / (float(rho) + 1)
    if t < 1e-14:
        return l + 8
    return min(NMAX, l + 80 + int(80.0 / -mp.log(t)))


def residual_profile(Nr: int, l: int, rho):
    """Relative L2 and kinetic residual norms after truncation at N terms (tail to ~1e-35)."""
    rho = mp.mpf(rho)
    nmax = _nmax(l, rho)
    c = coeffs(Nr, l, rho, nmax=nmax)
    m = len(c)
    ns = list(range(l + 1, nmax + 1))
    g = [gram_off(n, l) for n in ns[:-1]]
    # f norms: ||f||^2 = (2Nr)!/(2rho)^{2Nr+1};  <T> from Sturmian identity on full expansion
    f2 = mp.factorial(2 * Nr) / (2 * rho) ** (2 * Nr + 1)
    # tails: S2[i] = sum_{j>=i} c_j^2 ; X[i] = sum_{j>=i} g_j c_j c_{j+1}
    S2 = [mp.mpf(0)] * (m + 1)
    X = [mp.mpf(0)] * (m + 1)
    for i in range(m - 1, -1, -1):
        S2[i] = S2[i + 1] + c[i] ** 2
        X[i] = X[i + 1] + (g[i] * c[i] * c[i + 1] if i < m - 1 else 0)
    t_full = S2[0] - (S2[0] + 2 * X[0]) / 2      # <f|T|f> = sum c^2 - 1/2 c^T G c
    l2_full = S2[0] + 2 * X[0]
    rel_l2, rel_t = [], []
    for N in range(1, m):           # keep first N terms; residual = indices >= N
        l2 = S2[N] + 2 * X[N]
        tt = S2[N] - l2 / 2
        rel_l2.append(mp.sqrt(max(l2, 0) / f2))
        rel_t.append(mp.sqrt(max(tt, 0) / t_full))
    return rel_l2, rel_t, float(abs(l2_full / f2 - 1))


def terms_needed(prof, eps: float) -> int:
    for i, v in enumerate(prof):
        if v < eps:
            return i + 1
    return -1


TYPES = {"1s": (1, 0), "2s": (2, 0), "2p": (2, 1), "3s": (3, 0), "3p": (3, 1), "3d": (3, 2)}
RHOS = [1 / 16, 1 / 8, 1 / 4, 1 / 3, 1 / 2, 2 / 3, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 8.0, 16.0]

# Shell sets (exponents as used in noci_engine-based Paper 58 work; O = Clementi-Raimondi)
SHELLS = {
    "NaH": {"Na": [("1s", 10.63), ("2s", 3.3), ("2p", 3.44), ("3s", 0.836)], "H": [("1s", 1.0)]},
    "H2O": {"O": [("1s", 7.6579), ("2s", 2.2458), ("2p", 2.2266)], "H": [("1s", 1.0)]},
    "LiH": {"Li": [("1s", 2.69), ("2s", 0.65)], "H": [("1s", 1.0)]},
}


def main() -> None:
    res = {"checks": {"gram_max_abs_err": check_gram(), "deriv_coeff_max_rel_err": check_coeffs()}}
    print("checks", res["checks"])
    cache = {}

    def N_of(kind, rho, eps, which):
        key = (kind, round(float(rho), 12))
        if key not in cache:
            Nr, l = TYPES[kind]
            a, b, comp = residual_profile(Nr, l, rho)
            cache[key] = (a, b, comp)
        a, b, _ = cache[key]
        return terms_needed(a if which == "l2" else b, eps)

    tab = []
    for kind in TYPES:
        for rho in RHOS:
            row = {"type": kind, "rho": rho, "t": abs(rho - 1) / (rho + 1)}
            for eps in (1e-6, 1e-10):
                row[f"N_l2_{eps:g}"] = N_of(kind, rho, eps, "l2")
                row[f"N_kin_{eps:g}"] = N_of(kind, rho, eps, "kin")
            row["completeness_err"] = cache[(kind, round(float(rho), 12))][2]
            tab.append(row)
            print(row)
    res["table"] = tab

    # optimal k per centre and global, for the shell sets
    import numpy as np
    shells_out = {}
    for mol, cent in SHELLS.items():
        all_z = [z for c in cent.values() for (_, z) in c]
        out = {}

        def cost(shell_list, k, eps, which):
            ns = [N_of(kind, z / k, eps, which) for kind, z in shell_list]
            return ns

        for eps in (1e-6, 1e-10):
            for which in ("l2", "kin"):
                tag = f"{which}_{eps:g}"
                # per-centre optimum
                per = {}
                for cname, sl in cent.items():
                    zs = [z for _, z in sl]
                    grid = np.geomspace(min(zs) * 0.7, max(zs) * 1.3, 41)
                    best = None
                    for k in grid:
                        ns = cost(sl, k, eps, which)
                        if min(ns) < 0:
                            continue
                        score = (max(ns), sum(ns))
                        if best is None or score < best[0]:
                            best = (score, float(k), ns)
                    per[cname] = {"k_opt": best[1], "terms": dict(zip([s for s, _ in sl], best[2])),
                                  "max_terms": best[0][0]}
                # global shared k
                sl_all = [(f"{c}:{s}", z) for c, sl in cent.items() for s, z in sl]
                grid = np.geomspace(min(all_z) * 0.7, max(all_z) * 1.3, 41)
                best = None
                for k in grid:
                    ns = [N_of(s.split(":")[1], z / k, eps, which) for s, z in sl_all]
                    if min(ns) < 0:
                        continue
                    score = (max(ns), sum(ns))
                    if best is None or score < best[0]:
                        best = (score, float(k), ns)
                out[tag] = {"per_centre": per,
                            "global": {"k_opt": best[1], "terms": dict(zip([s for s, _ in sl_all], best[2])),
                                       "max_terms": best[0][0]},
                            "geo_mean_k": float(np.sqrt(min(all_z) * max(all_z)))}
                print(mol, tag, out[tag])
        shells_out[mol] = out
    res["shell_sets"] = shells_out
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
