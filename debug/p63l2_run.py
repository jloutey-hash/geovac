"""Scoring driver for record P63L2 (Paper 63, l > 0 attainment beyond the p_z-leaf chain): shapes B1, B2.

Objects and derivations: debug/p63odd_design_memo.md (B.1, B.2); integrands: debug/p63odd_engine.py.
  S1(eps) = int_1^oo I2(y) x^2 (2/3) Phi2x[f_p E](x) dx                           (B1, d root)
  S2(eps) = int_1^oo x^2 [(I0(y) + I1(y)/y) Phi0x - I2(y) Phi2x](x) dx             (B2, mixed root)
  y = (2 - eps) x.  Derived: eps^{3/2} S1 -> A1 = 1/(12 pi), kappa1 = 15/8;  eps^{1/2} S2 -> A2 = 3/(8 pi), kappa2 = -3/8.

Two routes per eps:
  float : p63odd_engine.branch{1,2}_f (log-grid trapezoid, 400001 nodes, Gauss-Laguerre potentials).
  referee: x in [1, 32]: 5 dyadic panels x 40 Gauss-Legendre nodes; e^{2x}Phi_Lx at those nodes from mp30 quadrature
           (p63odd_engine.phix_mp), computed once (mode `cache`, eps-independent) and frozen;
           x >= 32: the local-rule asymptotic series of e^{2x}Phi_Lx (p63odd_tail, NSER terms, the derivation's input),
           integrated with mp.quad at dps 30.
Modes:
  cache  : compute the referee nodes (parallel), write debug/data/p63l2_cache.json. No branch evaluated.
  check  : instrument checks only (series vs cached mp nodes on [16, 32]; GL vs trapezoid with the float integrand on
           [1, 32] at eps = 0.5). No branch at eps < 0.1 evaluated.
  run    : the scored run, eps = 1e-2 .. 1e-4 (9 points); writes debug/data/p63l2_run.json.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import p63odd_engine as EN  # noqa: E402
import p63odd_tail as TL  # noqa: E402

PANELS = [(1, 2), (2, 4), (4, 8), (8, 16), (16, 32)]
NGL = 40
XCUT = 32
NSER = 30
CACHE = HERE / "data" / "p63l2_cache.json"
CLOSED = {
    "B1": {"A": 1 / (12 * np.pi), "A_expr": "1/(12 pi)", "kappa": 15 / 8, "pw": 1.5, "powers": [0, 2, 3, 4]},
    "B2": {"A": 3 / (8 * np.pi), "A_expr": "3/(8 pi)", "kappa": -3 / 8, "pw": 0.5, "powers": [0, 1, 2, 3, 4]},
}


def gl_nodes():
    t, w = np.polynomial.legendre.leggauss(NGL)
    xs, ws = [], []
    for a, b in PANELS:
        xs += list((b - a) / 2 * t + (a + b) / 2)
        ws += list((b - a) / 2 * w)
    return xs, ws


def _node(args):
    x, L = args
    mp.mp.dps = 30
    t0 = time.time()
    v = EN.phix_mp(x, L)
    return x, L, mp.nstr(v, 28), time.time() - t0


def mode_cache(workers: int) -> None:
    xs, _ = gl_nodes()
    jobs = [(x, L) for L in (0, 2) for x in xs]
    done = {}
    if CACHE.exists():
        done = json.loads(CACHE.read_text())
    todo = [j for j in jobs if f"{j[1]}:{j[0]!r}" not in done]
    print(f"{len(todo)} nodes to compute ({len(jobs)} total)", flush=True)
    t0 = time.time()
    with Pool(workers) as pool:
        for i, (x, L, v, dt) in enumerate(pool.imap_unordered(_node, todo)):
            done[f"{L}:{x!r}"] = v
            if i % 10 == 0:
                CACHE.write_text(json.dumps(done, indent=0))
                print(f"{i+1}/{len(todo)}  L={L} x={x:.4f}  {dt:.0f}s  [{time.time()-t0:.0f}s]", flush=True)
    CACHE.write_text(json.dumps(done, indent=0))
    print("cache complete", len(done))


def series_coeffs():
    """e^{2x} Phi_Lx[f_p E](x) ~ x^a sum_j c_j x^-j, L = 0, 2 (NSER terms, mp floats)."""
    TL.N = NSER
    f_p = TL.Ser(sp.sqrt(sp.pi / 2), 1, sp.Rational(-1, 2), TL.a_coef(1, NSER)).scale(-1 / (2 * sp.pi**2))
    E_full = TL.local_potential(f_p, 1)
    w0 = f_p * E_full
    out = {}
    for L in (0, 2):
        Ph = TL.local_potential(w0, L)
        assert Ph.k == 2
        out[L] = (mp.mpf(sp.N(Ph.a, 40)), [mp.mpf(sp.N(Ph.pref * c, 40)) for c in Ph.c])
    return out


def phi_series(x, ser) -> mp.mpf:
    a, c = ser
    x = mp.mpf(x)
    return x**a * mp.fsum(cj * x ** (-j) for j, cj in enumerate(c))


def kern(y, which):
    """(scaled) radial kernels: I2 e^{-y};  (I0 + I1/y) e^{-y}."""
    if which == 2:
        return mp.besseli(2, y) * mp.exp(-y)
    return (mp.besseli(0, y) + mp.besseli(1, y) / y) * mp.exp(-y)


def integrand_mp(shape, x, eps, p0, p2):
    y = (2 - eps) * x
    damp = mp.exp(-eps * x) * x**2
    if shape == "B1":
        return kern(y, 2) * damp * mp.mpf(2) / 3 * p2
    return damp * (kern(y, 0) * p0 - kern(y, 2) * p2)


def referee(shape, eps, cache, ser) -> mp.mpf:
    xs, ws = gl_nodes()
    tot = mp.mpf(0)
    for x, w in zip(xs, ws):
        p0 = mp.mpf(cache[f"0:{x!r}"])
        p2 = mp.mpf(cache[f"2:{x!r}"])
        tot += w * integrand_mp(shape, mp.mpf(x), eps, p0, p2)
    f = lambda x: integrand_mp(shape, x, eps, phi_series(x, ser[0]), phi_series(x, ser[2]))  # noqa: E731
    pts = [XCUT, 64, 128, 256, 512, 1024, 4096, 16384, 65536, mp.inf]
    tot += mp.quad(f, pts)
    return tot


def mode_check() -> dict:
    mp.mp.dps = 30
    cache = json.loads(CACHE.read_text())
    t0 = time.time()
    ser = series_coeffs()
    out = {"series_build_s": time.time() - t0, "series_vs_cache": {}}
    xs, ws = gl_nodes()
    for L in (0, 2):
        a, c = ser[L]
        out["series_vs_cache"][L] = {
            "last_term_rel_at_32": float(abs(c[-1] * mp.mpf(32) ** (-(len(c) - 1)) / c[0])),
            "rel_at_nodes_16_32": [float(abs(phi_series(x, ser[L]) / mp.mpf(cache[f"{L}:{x!r}"]) - 1))
                                   for x in xs if x > 16][::8]}
    # GL panel adequacy with the float integrand, eps = 0.5 (not a scored eps)
    eps = 0.5
    xg = np.array(xs)
    yy = (2 - eps) * xg
    from scipy.special import i0e, i1e, ive
    h0, h2 = EN.h2_pair_f(xg)
    g1 = float(np.dot(ws, ive(2, yy) * EN.h1_f(xg) * np.exp(-eps * xg)))
    g2 = float(np.dot(ws, ((i0e(yy) + i1e(yy) / yy) * h0 - ive(2, yy) * h2) * np.exp(-eps * xg)))
    yl = np.linspace(0, np.log(32.0), 200001)
    xl = np.exp(yl)
    yy = (2 - eps) * xl
    h0, h2 = EN.h2_pair_f(xl)
    t1 = float(np.trapezoid(ive(2, yy) * EN.h1_f(xl) * np.exp(-eps * xl) * xl, yl))
    t2 = float(np.trapezoid(((i0e(yy) + i1e(yy) / yy) * h0 - ive(2, yy) * h2) * np.exp(-eps * xl) * xl, yl))
    out["gl_vs_trapezoid_1_32_eps0.5"] = {"B1": abs(g1 / t1 - 1), "B2": abs(g2 / t2 - 1)}
    # referee vs float at eps = 0.5 (not scored)
    for shape, bf in (("B1", EN.branch1_f), ("B2", EN.branch2_f)):
        r = referee(shape, mp.mpf(eps), cache, ser)
        f = bf(eps)
        out[f"referee_vs_float_eps0.5_{shape}"] = float(abs(f / r - 1))
    out["seconds"] = time.time() - t0
    print(json.dumps(out, indent=1))
    (HERE / "data" / "p63l2_check.json").write_text(json.dumps(out, indent=1, default=str))
    return out


def mode_run() -> None:
    mp.mp.dps = 30
    cache = json.loads(CACHE.read_text())
    ser = series_coeffs()
    eps_list = list(np.logspace(-2, -4, 9))
    res = {"rows": {"B1": [], "B2": []}}
    t0 = time.time()
    for shape, bf in (("B1", EN.branch1_f), ("B2", EN.branch2_f)):
        C = CLOSED[shape]
        for eps in eps_list:
            f = bf(eps)
            r = referee(shape, mp.mpf(eps), cache, ser)
            rf = float(r)
            res["rows"][shape].append({"eps": eps, "S_float64": f, "S_referee": mp.nstr(r, 22),
                                       "rel_float_vs_referee": (f - rf) / rf, "scaled_over_A": eps ** C["pw"] * rf / C["A"]})
            print(shape, f"eps={eps:.3e} S_ref={mp.nstr(r, 15)} rel {(f - rf) / rf:.1e} "
                  f"scaled/A {eps ** C['pw'] * rf / C['A']:.7f} [{time.time() - t0:.0f}s]", flush=True)
        e = np.array(eps_list)
        S = np.array([float(mp.mpf(row["S_referee"])) for row in res["rows"][shape]])
        le, lS = np.log(e), np.log(np.abs(S))
        u = np.sqrt(e)
        y = e ** C["pw"] * S

        def fit(sel):
            M = np.stack([u[sel] ** p for p in C["powers"]], 1)
            c, *_ = np.linalg.lstsq(M, y[sel], rcond=None)
            return dict(zip(C["powers"], c))
        fa, fs = fit(slice(None)), fit(slice(3, None))
        # reported, not scored: the same fit on the float route (full integrand from x = 1, independent of the series)
        Sf = np.array([row["S_float64"] for row in res["rows"][shape]])
        yf = e ** C["pw"] * Sf
        Mf = np.stack([u[3:] ** p for p in C["powers"]], 1)
        cf, *_ = np.linalg.lstsq(Mf, yf[3:], rcond=None)
        cf = dict(zip(C["powers"], cf))
        res[shape] = {
            "closed_form": {k: C[k] for k in ("A", "A_expr", "kappa")},
            "sign_of_S": [int(np.sign(s)) for s in S],
            "log_slope_last5": float(np.polyfit(le[-5:], lS[-5:], 1)[0]),
            "log_slope_all": float(np.polyfit(le, lS, 1)[0]),
            "local_log_slopes": [float((lS[i + 1] - lS[i]) / (le[i + 1] - le[i])) for i in range(len(e) - 1)],
            "fit_all9": {"A_ratio": float(fa[0] / C["A"]), "u2_over_Akappa": float(fa[2] / (C["A"] * C["kappa"]))},
            "fit_smallest6": {"A_ratio": float(fs[0] / C["A"]), "u2_over_Akappa": float(fs[2] / (C["A"] * C["kappa"]))},
            "reported_float_route_fit_smallest6": {"A_ratio": float(cf[0] / C["A"]),
                                                   "u2_over_Akappa": float(cf[2] / (C["A"] * C["kappa"]))},
            "max_rel_float_vs_referee": float(max(abs(r_["rel_float_vs_referee"]) for r_ in res["rows"][shape])),
        }
        print(shape, json.dumps(res[shape], indent=1), flush=True)
    res["runtime_s"] = time.time() - t0
    (HERE / "data" / "p63l2_run.json").write_text(json.dumps(res, indent=1))
    print("wrote debug/data/p63l2_run.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["cache", "check", "run"])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.mode == "cache":
        mode_cache(a.workers)
    elif a.mode == "check":
        mode_check()
    else:
        mode_run()
