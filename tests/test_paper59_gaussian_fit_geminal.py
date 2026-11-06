"""Guard for Paper 59 sec:density_ring, the Gaussian-geminal fit priced on Y (record
docs/predictions/2026-10-25_gaussian_fit_geminal.md, GF1-GF3; GF1 missed, GF2 and GF3 held).

Y at density scales (1, 2) and apex t, with exp(-gamma r) geminals on the two edges (gamma = 3/2, 3), each replaced by
an n-term fit sum_i c_i exp(-alpha_i gamma^2 r^2) made once at gamma = 1:
Y_n = sum_ij c_i c_j Y_GG(alpha_i g1^2, alpha_j g2^2), E_n(t) = Y_n/Y - 1.  Two fits of exp(-r):
  W1: min int_0^oo r^2 [e^{-r} - fit]^2 dr (plain norm);  W2: the same with the extra weight e^{-r^2}.
Claims guarded, as the paper states them ([MEASURED]):
  * W1 at three and six terms gives the STO-3G and STO-6G exponents;
  * six-term W1 makes Y too small by 3.6e-5 at t = 1/2 down to 2.4e-6 at t = 4;
  * the error falls like e^{-b sqrt(n)}: ln|E_n| (W1, t = 1, n = 2..10) fits a + b sqrt(n) with RSS ratio 0.073 against
    a + b n (matrix row);
  * W2 loses its tail, and its six-term error is 5.4e-2 to 1.7e-2;  so the error is set by the tail, not the cusp.
Data: tests/data/paper59_gf_2026_11_01.json, a copy of debug/data/gf_fits.json and gf_results.json built by
debug/build_gf_guard_data_2026_11_01.py (no debug/ import here).  Y_exact there is the 60-digit grid's Y, an independent
route from the Gaussian engine below (a copy of debug/gf_engine.py, as in test_paper63_production_geminal_fit.py).
Independent recomputations: the published STO-3G/6G tables (Hehre-Stewart-Pople 1969), a from-scratch float refit at
n = 3 (erfcx moments, Nelder-Mead), re-solved coefficients and gradients at 40 digits, pointwise fit errors,
the RSS fits, and Y_1 (fast) and Y_6 (slow) through the engine against the exact Y.
Scope: one configuration, model densities; the production weight is not reproduced (the paper's cap).

Fire-tested 2026-11-01 (debug/qa/fire_test.py), every plant fired:
  P1 stored W1 n=6 exponent 23.1030314881 -> 23.1130314881 (a non-STO-6G fit)    -> test_w1_fit_is_sto_ng
  P2 refit in the 1-D norm (r^2 weight dropped)                                   -> test_w1_fit_is_sto_ng
  P3 W2 checked against weight exp(-2 r^2) (fit labelled with the wrong weight)   -> test_stored_fits_are_least_squares_optima
  P4 a stored W2 n=6 coefficient perturbed by 5e-6 relative (not re-solved)     -> test_stored_fits_are_least_squares_optima
     (a 4.5e-13 perturbation, below the 1e-12 tolerance, did not fire: the leg resolves c to 1e-12, not better)
  P5 stored E sign flipped at W1 n=6 t=1/2 (fit makes Y too large)               -> test_six_term_errors
  P6 W1 and W2 six-term rows swapped                                              -> test_six_term_errors
  P7 GF3 reversed: W1 n=6 t=1/2 and t=4 errors swapped                            -> test_six_term_errors
  P8 geometric rate: W1 t=1 errors replaced by exp(-1.3 n)                        -> test_rate_is_root_exponential
  P9 W2 fit given W1's exponents (W2 coefficients kept)                           -> test_w2_loses_its_tail
  P10 engine apex density scale doubled (rho(t) -> rho(2t))                       -> test_one_term_y_recomputes
  P11 the second geminal at gamma = 2 for 3                                       -> test_one_term_y_recomputes
"""
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest
from scipy.optimize import minimize
from scipy.special import dawsn, erfcx, gamma

DATA = json.loads((Path(__file__).parent / "data" / "paper59_gf_2026_11_01.json").read_text())
PI = np.pi
TS = DATA["t"]
G1, G2 = 1.5, 3.0
Y_EXACT = np.array(DATA["Y_exact"])
E = {w: {int(n): np.array(v) for n, v in byn.items()} for w, byn in DATA["E"].items()}
FIT = {w: {int(n): ([mp.mpf(a) for a in f["alpha"]], [mp.mpf(c) for c in f["c"]]) for n, f in byn.items()}
       for w, byn in DATA["fits"].items()}
SHIFT = {"W1": 0, "W2": 1}  # extra weight exp(-SHIFT r^2)

# Hehre, Stewart, Pople, J. Chem. Phys. 51, 2657 (1969), zeta = 1: (exponent, contraction coefficient of the
# normalized primitive (2a/pi)^(3/4) e^{-a r^2}); the contracted function is normalized.
STO3G = [(0.109818, 0.4446345422), (0.405771, 0.5353281423), (2.22766, 0.1543289673)]
STO6G = [(0.06510953954, 0.1303340841), (0.1580884151, 0.4164915298), (0.4070988982, 0.3705627997),
         (1.185056519, 0.1685383049), (4.235915534, 0.04936149294), (23.10303149, 0.009163596281)]


def f64(n, w="W1"):
    a, c = FIT[w][n]
    return np.array([float(x) for x in a]), np.array([float(x) for x in c])


def sto_compare(alpha, c, table):
    """Relative differences of exponents and of normalized coefficients against a published table."""
    pa = np.array([p[0] for p in table])
    pc = np.array([p[1] for p in table]) * (2 * pa / PI) ** 0.75 * np.sqrt(PI)  # in units of e^{-r}/sqrt(pi)
    norm2 = PI ** 1.5 * np.sum(np.outer(c, c) / (alpha[:, None] + alpha[None, :]) ** 1.5) / PI
    return alpha / pa - 1, c / np.sqrt(norm2) / pc - 1


def refit_float(n, r_power=2):
    """From-scratch variable-projection fit of e^{-r} in the norm int r^r_power (.)^2 dr, float64, closed-form
    moments (erfcx); it does not read the stored fits."""
    def parts(al):
        b = []
        for a in al:
            I = [np.sqrt(PI / (4 * a)) * erfcx(1 / (2 * np.sqrt(a)))]
            I.append((1 - I[0]) / (2 * a))
            for k in range(1, r_power):
                I.append((k * I[k - 1] - I[k]) / (2 * a))
            b.append(I[r_power])
        s = al[:, None] + al[None, :]
        return gamma((r_power + 1) / 2) / (2 * s ** ((r_power + 1) / 2)), np.array(b)

    def F(la):
        G, b = parts(np.exp(la))
        return -b @ np.linalg.solve(G, b)

    best = min((minimize(F, np.log(a0) + np.log(4.0) * np.arange(n), method="Nelder-Mead",
                         options={"xatol": 1e-10, "fatol": 1e-16, "maxiter": 6000}) for a0 in (0.08, 0.15)),
               key=lambda r: r.fun)
    al = np.sort(np.exp(best.x))
    G, b = parts(al)
    return al, np.linalg.solve(G, b)


def test_w1_fit_is_sto_ng():
    """W1 at n = 3 and n = 6 is the STO-3G / STO-6G fit.  Rejects: a stored fit that is not the plain-norm optimum
    (exponents off the published table), and a fit made in the 1-D norm without r^2 (the from-scratch refit, which
    does not read the stored fit, must land on STO-3G too)."""
    da, dc = sto_compare(*f64(3), STO3G)
    assert np.all(np.abs(da) < 1e-5) and np.all(np.abs(dc) < 1e-5), (da, dc)
    da, dc = sto_compare(*f64(6), STO6G)
    assert np.all(np.abs(da) < 1e-8) and np.all(np.abs(dc) < 1e-8), (da, dc)
    al, c = refit_float(3, r_power=2)
    da, dc = sto_compare(al, c, STO3G)
    assert np.all(np.abs(da) < 1e-5) and np.all(np.abs(dc) < 1e-5), (al, da, dc)


def Ik(a, kmax):
    a = mp.mpf(a)
    I = [mp.sqrt(mp.pi / (4 * a)) * mp.exp(1 / (4 * a)) * mp.erfc(1 / (2 * mp.sqrt(a)))]
    I.append((1 - I[0]) / (2 * a))
    for k in range(1, kmax):
        I.append((k * I[k - 1] - I[k]) / (2 * a))
    return I


def test_stored_fits_are_least_squares_optima():
    """Every stored fit, n = 2..10 under both weights, is a stationary point of its own stated norm: the linear
    coefficients re-solved at 40 digits reproduce the stored ones and the gradient in log alpha vanishes.  Rejects:
    a W2 fit made with a different Gaussian weight (e.g. exp(-2 r^2)) and labelled exp(-r^2), and coefficients that
    were not solved for (a perturbed c)."""
    with mp.workdps(40):
        for w, shift in SHIFT.items():
            f2 = mp.quad(lambda r: r * r * mp.exp(-2 * r - shift * r * r), [0, 1, 5, mp.inf])
            for n in range(2, 11):
                al, cs = FIT[w][n]
                G = mp.matrix(n, n); G4 = mp.matrix(n, n); b = mp.matrix(n, 1); b4 = mp.matrix(n, 1)
                for i in range(n):
                    I = Ik(al[i] + shift, 4)
                    b[i], b4[i] = I[2], I[4]
                    for j in range(n):
                        s = al[i] + al[j] + shift
                        G[i, j] = mp.sqrt(mp.pi) / (4 * s ** mp.mpf(1.5))
                        G4[i, j] = 3 * mp.sqrt(mp.pi) / (8 * s ** mp.mpf(2.5))
                c = mp.lu_solve(G, b)
                rel_c = max(abs(c[i] / cs[i] - 1) for i in range(n))
                assert rel_c < 1e-12, (w, n, rel_c)
                F = f2 - (b.T * c)[0]
                grad = max(abs(2 * c[k] * al[k] * (b4[k] - sum(c[j] * G4[j, k] for j in range(n))))
                           for k in range(n))
                assert grad < 1e-6 * F, (w, n, grad, F)


def test_six_term_errors():
    """The six-term errors the paper prints, their sign, and their order in t.  Rejects: a fit that makes Y too large
    (sign), W1 and W2 swapped, and GF3 reversed (|E| growing as the apex becomes compact)."""
    assert np.allclose(E["W1"][6], [-3.6e-5, -1.36e-5, -4.1e-6, -2.4e-6], rtol=0.03), E["W1"][6]
    assert np.allclose(E["W2"][6], [-5.4e-2, -4.1e-2, -2.8e-2, -1.7e-2], rtol=0.03), E["W2"][6]
    for w in E:
        for n in range(1, 11):
            assert np.all(E[w][n] < 0), (w, n)
            assert np.all(np.diff(np.abs(E[w][n])) < 0), (w, n)  # falls as t grows
        assert all(np.all(np.abs(E[w][n + 1]) < np.abs(E[w][n])) for n in range(1, 10)), w
    assert abs(E["W1"][6][0]) / abs(E["W1"][6][3]) > 10  # record: a factor of 15


def rss(y, x):
    X = np.c_[np.ones_like(x), x]
    return float(np.sum((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2))


def test_rate_is_root_exponential():
    """ln|E_n|, W1, t = 1, n = 2..10: RSS(a + b sqrt n) / RSS(a + b n) = 0.073 (registered bound 0.25), slope b in the
    record's -6.1..-7.8.  Rejects: a geometric rate e^{-c n}, for which the linear fit is exact and the ratio large."""
    n = np.arange(2, 11, dtype=float)
    y = np.log(np.abs([E["W1"][int(k)][1] for k in n]))
    ratio = rss(y, np.sqrt(n)) / rss(y, n)
    assert abs(ratio - 0.073) < 0.003, ratio
    b = np.polyfit(np.sqrt(n), y, 1)[0]
    assert -7.8 <= b <= -6.1, b
    for t in range(4):
        yt = np.log(np.abs([E["W1"][int(k)][t] for k in n]))
        assert rss(yt, np.sqrt(n)) / rss(yt, n) < 0.25, TS[t]


def test_w2_loses_its_tail():
    """Pointwise relative error of the six-term fits, from the stored coefficients: W1 is worst at the cusp (-6.0e-2
    at r = 0) and below 1e-2 for 1/2 <= r <= 6; W2 is better at the cusp (-2.4e-2) and loses the tail (-1.1e-2 at
    r = 3, -0.19 at 4, -0.80 at 6).  Rejects: a W2 fit with W1's tail, and the reading reversed (W2 worse at the
    cusp)."""
    def rel(w, r):
        a, c = f64(6, w)
        return np.sum(c[:, None] * np.exp(-a[:, None] * r * r), axis=0) * np.exp(r) - 1
    pts = np.array([0.0, 3.0, 4.0, 6.0])
    assert np.allclose(rel("W2", pts), [-2.4e-2, -1.1e-2, -0.19, -0.80], rtol=0.05), rel("W2", pts)
    assert abs(rel("W1", np.array([0.0]))[0] + 6.0e-2) < 1e-3
    r = np.linspace(0.5, 6.0, 5501)
    assert np.max(np.abs(rel("W1", r))) < 1e-2
    tail = np.linspace(4.0, 6.0, 2001)
    assert np.min(np.abs(rel("W2", tail))) > 10 * np.max(np.abs(rel("W1", tail)))


# ---- the Gaussian engine (copy of debug/gf_engine.py) ----

def rho(c, k):
    return c ** 3 * (k * k + c * c) ** -1.5


def G_loop(R, mu):
    x = R / (2 * np.sqrt(mu))
    small = x < 1e-6
    out = np.empty_like(R)
    out[~small] = 8 * mu * PI ** 1.5 * dawsn(x[~small]) / R[~small]
    out[small] = 8 * mu * PI ** 1.5 / (2 * np.sqrt(mu))
    return out


def gl(n, lo, hi):
    x, w = np.polynomial.legendre.leggauss(n)
    return (hi - lo) / 2 * x + (hi + lo) / 2, (hi - lo) / 2 * w


def half_line(n, lo, scale):
    u, wu = gl(n, 0.0, 1.0)
    return lo + scale * (u / (1 - u)) ** 2, wu * scale * 2 * u / (1 - u) ** 3


def geo(lo, hi, n, h0=0.5):
    br, h = [lo], h0
    while br[-1] + h < hi:
        br.append(lo + h)
        h *= 2
    br.append(hi)
    xs, ws = zip(*[gl(n, br[j], br[j + 1]) for j in range(len(br) - 1)])
    return np.concatenate(xs), np.concatenate(ws)


def two_sided(lo, hi, n, h0=0.5):
    m = (lo + hi) / 2
    x1, w1 = geo(lo, m, n, h0)
    x2, w2 = geo(-hi, -m, n, h0)
    return np.concatenate([x1, -x2]), np.concatenate([w1, w2])


def Y_GG(A, B, ts, a=1.0, b=2.0, n=24):
    cut = 4.0 + 12 * np.sqrt(A + B)
    x1, w1 = geo(0.0, cut, n)
    x2, w2 = half_line(2 * n, cut, cut)
    w, ww = np.concatenate([x1, x2]), np.concatenate([w1, w2])
    mu, c = A * B / (A + B), A / (A + B)
    pref_T = (2 * PI) ** -3 * 4 * PI * (PI * PI / (A * B)) ** 1.5
    ph = np.empty_like(w)
    for i, wi in enumerate(w):
        cu = wi + 8.0 + 12 * np.sqrt(mu)
        s1, ws1 = geo(wi, cu, n)
        s2, ws2 = half_line(2 * n, cu, cu)
        sg, wsg = np.concatenate([s1, s2]), np.concatenate([ws1, ws2])
        dl, wdl = two_sided(-wi, wi, n, min(0.5, wi / 2))
        S, D = np.meshgrid(sg, dl, indexing="ij")
        p, q = (S + D) / 2, (S - D) / 2
        R = np.sqrt(np.maximum((1 - c) * p * p + c * q * q - c * (1 - c) * wi * wi, 0.0))
        T = pref_T * np.exp(-wi * wi / (4 * (A + B))) * G_loop(R, mu)
        ph[i] = (2 * PI) ** -3 * 2 * PI / wi * np.sum(np.outer(wsg, wdl) * (S * S - D * D) / 8 * rho(a, p) * rho(b, q) * T)
    return np.array([1 / (2 * PI * PI) * np.sum(ww * w * w * rho(t, w) * ph) for t in ts])


def E_fit(w, n):
    al, c = f64(n, w)
    Yn = sum(c[i] * c[j] * Y_GG(al[i] * G1 ** 2, al[j] * G2 ** 2, TS) for i in range(n) for j in range(n))
    return Yn / Y_EXACT - 1


def test_one_term_y_recomputes():
    """E_1 for both weights recomputed through the engine against the exact Y of the 60-digit grid (an independent
    route) reproduces the stored -1.3e-1..-9.7e-2 (W1) and -5.2e-1..-3.6e-1 (W2).  Rejects: a wrong apex density
    scale and a wrong geminal exponent in the pipeline."""
    for w in ("W1", "W2"):
        got = E_fit(w, 1)
        assert np.allclose(got, E[w][1], rtol=1e-8, atol=0), (w, got, E[w][1])


@pytest.mark.slow
def test_six_term_y_recomputes():
    """E_6 for both weights recomputed (72 engine calls) against the exact Y: W1 reproduces -3.6e-5..-2.4e-6, which
    checks the engine against the independent 60-digit Y to about 1e-7; W2 reproduces -5.4e-2..-1.7e-2."""
    for w in ("W1", "W2"):
        got = E_fit(w, 6)
        assert np.allclose(got, E[w][6], rtol=1e-6, atol=1e-10), (w, got, E[w][6])
