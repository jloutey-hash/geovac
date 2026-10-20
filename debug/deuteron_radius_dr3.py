"""DR3: number class of the matrix elements (basis-potential pairing).

(a) Laguerre, registered orthonormal basis phi_n: are all entries of H (MeV) and <r^2> (fm^2)
    rational for rational lambda?  Exact check (entry = rational / sqrt(S_mm S_nn)).
(b) Laguerre, unnormalized f_n (S exactly diagonal rational): all entries rational (Fraction).
(c) HO <0s|e^{-mu r}/r|0s>: closed form with e^{y^2} erfc(y), y = mu b/2; nonzero coefficient;
    closed form vs quadrature to 1e-30 at three b.
(d) Minnesota Gaussian e^{-k r^2}: algebraic in HO (powers of (1+k b^2)^{-1/2}); carries
    e^{z^2} erfc(z) in Laguerre.
"""
from __future__ import annotations

import json
from fractions import Fraction as Fr
from math import isqrt

import mpmath as mp
import sympy as sp

import deuteron_radius_laguerre as L


def is_square_rational(q: Fr) -> bool:
    if q < 0:
        return False
    n, d = q.numerator, q.denominator
    return isqrt(n) ** 2 == n and isqrt(d) ** 2 == d


def laguerre_literal(N: int, lam: Fr):
    M = L.exact_mats(N, lam, L.MT3)
    assert L.check_orthogonality(M["S"], lam)
    res = {}
    for key in ("H", "R2", "T", "V"):
        irr = 0
        tot = 0
        first = None
        for i in range(N):
            for j in range(N):
                tot += 1
                q = M[key][i][j]
                if q == 0:
                    continue
                if not is_square_rational(M["S"][i][i] * M["S"][j][j]):
                    irr += 1
                    if first is None:
                        first = (i, j, str(q), str(M["S"][i][i] * M["S"][j][j]))
        res[key] = {"entries": tot, "irrational_in_orthonormal_basis": irr, "first_irrational": first}
    # unnormalized: rational by construction; confirm types
    res["unnormalized_all_Fraction"] = all(isinstance(M[k][i][j], Fr) for k in ("H", "R2", "S")
                                           for i in range(N) for j in range(N))
    res["unnormalized_S_diagonal_exact"] = True
    return res


def ho_yukawa_0s():
    r, b, mu = sp.symbols("r b mu", positive=True)
    N2 = 4 / (sp.sqrt(sp.pi) * b ** 3)
    expr = sp.simplify(sp.integrate(N2 * r ** 2 * sp.exp(-r ** 2 / b ** 2) * sp.exp(-mu * r) / r, (r, 0, sp.oo)))
    y = mu * b / 2
    closed = 2 / (sp.sqrt(sp.pi) * b) * (1 - sp.sqrt(sp.pi) * y * sp.exp(y ** 2) * sp.erfc(y))
    sym_diff = sp.simplify(sp.expand_func(expr - closed).rewrite(sp.erfc))
    out = {"sympy_integral": str(expr), "closed_form": str(closed), "sympy_difference": str(sym_diff),
           "erfc_coefficient": "-2*y/b (nonzero for mu, b > 0)"}
    mp.mp.dps = 45
    checks = []
    for hw in ("10", "20", "40"):
        bb = mp.sqrt(mp.mpf("41.47") * 2 / mp.mpf(hw))
        for m in ("1.55", "3.11"):
            mm = mp.mpf(m)
            yy = mm * bb / 2
            cf = 2 / (mp.sqrt(mp.pi) * bb) * (1 - mp.sqrt(mp.pi) * yy * mp.exp(yy ** 2) * mp.erfc(yy))
            q = mp.quad(lambda x: 4 / (mp.sqrt(mp.pi) * bb ** 3) * x * mp.exp(-x ** 2 / bb ** 2 - mm * x),
                        [0, bb, 2 * bb, 4 * bb, 8 * bb, mp.inf])
            no_erfc = 2 / (mp.sqrt(mp.pi) * bb)
            checks.append({"hw": hw, "mu": m, "b": mp.nstr(bb, 20), "closed": mp.nstr(cf, 35),
                           "quad": mp.nstr(q, 35), "abs_diff": mp.nstr(abs(cf - q), 5),
                           "erfc_term_value": mp.nstr(cf - no_erfc, 20)})
    out["numeric"] = checks
    out["PASS"] = all(mp.mpf(c["abs_diff"]) < mp.mpf("1e-30") for c in checks)
    return out


def minnesota_pairing():
    r, b, k, lam = sp.symbols("r b kappa lambda", positive=True)
    out = {}
    # HO s-wave functions n = 0..2 (unnormalized: r L_n^{1/2}(r^2/b^2) e^{-r^2/(2b^2)})
    ho = [r * sp.assoc_laguerre(n, sp.Rational(1, 2), r ** 2 / b ** 2) * sp.exp(-r ** 2 / (2 * b ** 2))
          for n in range(3)]
    norms = [sp.integrate(f ** 2, (r, 0, sp.oo)) for f in ho]
    tab = {}
    t = sp.symbols("t", positive=True)  # t = (1 + kappa b^2)^{-1/2}
    alg = True
    for m in range(3):
        for n in range(m, 3):
            e = sp.integrate(ho[m] * ho[n] * sp.exp(-k * r ** 2), (r, 0, sp.oo)) / sp.sqrt(norms[m] * norms[n])
            e = sp.simplify(e)
            et = sp.simplify(e.subs(k, (1 / t ** 2 - 1) / b ** 2))
            ok = et.is_polynomial(t) if et.free_symbols <= {t} else False
            alg = alg and bool(ok)
            tab[f"{m}{n}"] = {"element": str(e), "in_t": str(sp.expand(et)), "polynomial_in_t": bool(ok)}
    out["HO_minnesota"] = tab
    out["HO_algebraic"] = alg
    # Laguerre phi_0 = N0 x e^{-x/2}, x = 2 lam r, N0^2 = lam
    lag0 = sp.sqrt(lam) * 2 * lam * r * sp.exp(-lam * r)
    e = sp.simplify(sp.integrate(lag0 ** 2 * sp.exp(-k * r ** 2), (r, 0, sp.oo)))
    out["Laguerre_00_minnesota"] = str(e)
    out["Laguerre_contains_erfc"] = e.has(sp.erfc) or e.has(sp.erf)
    # numeric check of the erfc closed form at Minnesota triplet kappa_T = 0.639, lambda = 1.3949
    mp.mp.dps = 45
    kk = mp.mpf("0.639")
    ll = mp.mpf("1.3949")
    fn = sp.lambdify((lam, k), e, "mpmath")
    cf = fn(ll, kk)
    q = mp.quad(lambda x: 4 * ll ** 3 * x ** 2 * mp.exp(-2 * ll * x - kk * x ** 2), [0, 1, 2, 4, 8, mp.inf])
    z = ll / mp.sqrt(kk)
    out["Laguerre_numeric"] = {"closed": mp.nstr(cf, 35), "quad": mp.nstr(q, 35),
                               "abs_diff": mp.nstr(abs(cf - q), 5), "z": mp.nstr(z, 15)}
    return out


if __name__ == "__main__":
    out = {}
    for lam in (Fr(13949, 10000), Fr(1, 2)):
        out[f"laguerre_literal_N40_lam_{lam}"] = laguerre_literal(40, lam)
    out["ho_yukawa_0s"] = ho_yukawa_0s()
    out["minnesota"] = minnesota_pairing()
    json.dump(out, open("data/deuteron_radius_dr3.json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
