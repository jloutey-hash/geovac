"""Atlas gap memo, door A: exact Rayleigh-Schroedinger series for the hydrogen
ground state in a uniform electric field (Stark) and in a uniform magnetic
field (diamagnetic term only; the linear Zeeman term vanishes for m = 0).

Method (exact, Fractions only). psi_k = e^{-r} sum_{i,j} c_{ij} r^i P_j(u),
u = cos(theta). With L = H0 + 1/2,
    L[e^{-r} r^i P_j] = e^{-r} [ i r^{i-1} + (j(j+1) - i(i+1))/2 r^{i-2} ] P_j,
which is triangular in i, so each order is solved top-down in exact
rationals. Intermediate normalization <psi0|psi_k> = 0 (k >= 1).
E_k = <psi0|V|psi_{k-1}> / <psi0|psi0>.

Stark:  V = F z = F r u          (series in F; odd orders vanish)
Diamag: V = (B^2/8)(x^2+y^2) = (beta/8) r^2 (1-u^2), beta = B^2 (series in beta)

Output: debug/data/atlas_gap_field_series.json
"""
from fractions import Fraction as Fr
from math import factorial, pi
import json, sys, time

def u_times(poly):
    """multiply {(i,j): c} by u in the Legendre index: u P_j = ((j+1)P_{j+1} + j P_{j-1})/(2j+1)"""
    out = {}
    for (i, j), c in poly.items():
        a = Fr(j + 1, 2 * j + 1)
        out[(i, j + 1)] = out.get((i, j + 1), 0) + c * a
        if j > 0:
            b = Fr(j, 2 * j + 1)
            out[(i, j - 1)] = out.get((i, j - 1), 0) + c * b
    return {k: v for k, v in out.items() if v != 0}

def r_times(poly, p):
    return {(i + p, j): c for (i, j), c in poly.items()}

def add(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v != 0}

def scale(a, s):
    return {k: v * s for k, v in a.items() if v * s != 0}

def V_stark(poly):
    return u_times(r_times(poly, 1))

def V_diamag(poly):
    # (1/8) r^2 (1 - u^2)
    r2 = r_times(poly, 2)
    return scale(add(r2, u_times(u_times(r2)), -1), Fr(1, 8))

def solve_L(rhs):
    """solve L g = rhs, g = sum a_{ij} r^i P_j, for each j separately (top-down)."""
    byj = {}
    for (i, j), c in rhs.items():
        byj.setdefault(j, {})[i] = c
    sol = {}
    for j, b in byj.items():
        top = max(b)
        a = {}
        cst = lambda i: Fr(j * (j + 1) - i * (i + 1), 2)
        # coefficient of r^m: (m+1) a_{m+1} + cst(m+2) a_{m+2} = b_m
        for m in range(top, j - 2, -1):
            if m + 1 < j:
                break
            val = b.get(m, Fr(0)) - cst(m + 2) * a.get(m + 2, Fr(0))
            if m + 1 == 0:
                continue
            a[m + 1] = val / (m + 1)
        # consistency: coefficient of r^{j-2} (j>=1) or r^{-1} (j=0)
        lowest = j - 2 if j >= 1 else -1
        resid = b.get(lowest, Fr(0)) - cst(lowest + 2) * a.get(lowest + 2, Fr(0)) - (lowest + 1) * a.get(lowest + 1, Fr(0))
        if j == 0:
            resid = b.get(-1, Fr(0)) - cst(1) * a.get(1, Fr(0))
        assert resid == 0, (j, resid)
        for i, v in a.items():
            if v != 0:
                sol[(i, j)] = v
    # check no negative powers in rhs other than handled
    return sol

def overlap0(poly):
    """<e^{-r}|e^{-r} poly> / <e^{-r}|e^{-r}>; only P_0 survives; int r^{i+2} e^{-2r} = (i+2)!/2^{i+3}"""
    tot = Fr(0)
    for (i, j), c in poly.items():
        if j == 0:
            tot += c * Fr(factorial(i + 2), 2 ** (i + 3))
    return tot / Fr(2, 8)  # norm: 2!/2^3 = 1/4

def series(V, K):
    psi = [{(0, 0): Fr(1)}]
    E = [Fr(-1, 2)]
    for k in range(1, K + 1):
        Ek = overlap0(V(psi[k - 1]))
        E.append(Ek)
        rhs = scale(V(psi[k - 1]), -1)
        for jj in range(1, k + 1):
            if E[jj] != 0:
                rhs = add(rhs, scale(psi[k - jj], E[jj]))
        # rhs must contain no r^{-1} etc. except as consistent
        g = solve_L(rhs)
        # intermediate normalization
        g[(0, 0)] = g.get((0, 0), Fr(0)) - overlap0(g)
        g = {kk: v for kk, v in g.items() if v != 0}
        psi.append(g)
    E.append(overlap0(V(psi[K])))
    return E

if __name__ == "__main__":
    t0 = time.time()
    KS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    KD = int(sys.argv[2]) if len(sys.argv) > 2 else 25
    Es = series(V_stark, KS)
    Ed = series(V_diamag, KD)
    out = {"stark_coeffs_F^k": {str(k): str(v) for k, v in enumerate(Es)},
           "diamag_coeffs_beta^k(beta=B^2)": {str(k): str(v) for k, v in enumerate(Ed)}}
    print("Stark:")
    for k in range(0, 9):
        print(k, Es[k])
    print("Diamag:")
    for k in range(0, 5):
        print(k, Ed[k])
    # large order: Stark a_k (coefficient of F^{2k}) vs -(6/pi)(3/2)^{2k}(2k)!
    rat = []
    for k in range(1, len(Es) // 2):
        a = Es[2 * k]
        pred = -(6 / pi) * (1.5) ** (2 * k) * factorial(2 * k)
        rat.append((k, float(a) / pred))
    print("Stark a_k / [-(6/pi)(3/2)^{2k}(2k)!]:", [(k, round(x, 6)) for k, x in rat[-6:]])
    # Richardson-like: ratio test for diamag: a_k/a_{k-1} / k^2 -> ?
    dr = []
    for k in range(3, len(Ed)):
        dr.append((k, float(Ed[k] / Ed[k - 1]) / (k * k)))
    print("Diamag (a_k/a_{k-1})/k^2:", [(k, round(x, 6)) for k, x in dr[-6:]], " -4/pi^2 =", -4 / pi ** 2)
    out["stark_largeorder_ratio_to_minus6overpi"] = [[k, x] for k, x in rat]
    out["diamag_ratio_over_k2"] = [[k, x] for k, x in dr]
    out["seconds"] = time.time() - t0
    json.dump(out, open("debug/data/atlas_gap_field_series.json", "w"), indent=1)
    print("time", time.time() - t0)
