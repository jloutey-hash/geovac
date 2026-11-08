"""Probe A -- analysis: tables, log-log fits with residuals, Bernstein floors."""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
from geovac.sturmian_sigma_law import sw_cross_block, sigma_spectrum  # noqa: E402

ROWS = json.load(open(os.path.join(HERE, "data", "probeA_scan.json")))


def spectrum(s, n):
    sig = sigma_spectrum(sw_cross_block(s, n))
    return np.sort(np.concatenate([1.0 - sig, 1.0 + sig]))


def bernstein_floor(lam, mode, conv, headroom=0.5):
    """Rigorous lower bound on deg p from |p'(x)| <= d/sqrt(1-x^2) (Bernstein).

    Between two accuracy points, |p(x_b)-p(x_a)| <= d |arcsin x_b - arcsin x_a|,
    so d >= |y_a - y_b| / |arcsin x_b - arcsin x_a|  (eps -> 0 limit).
    """
    lam = np.sort(lam); lo, hi = lam[0], lam[-1]
    if mode == "aware":
        lm = lam
    else:
        lm = lo * (hi / lo) ** np.linspace(0, 1, 20001)
    x = lm / hi if conv == "p60" else (2 * lm - (lo + hi)) / (hi - lo)
    y = headroom * np.sqrt(lo / lm)
    a = np.arcsin(np.clip(x, -1, 1))
    return float(np.max(np.abs(np.diff(y)) / np.maximum(np.abs(np.diff(a)), 1e-300)))


def fit(xs, ys):
    """log-log least squares; returns slope, intercept, R^2, max|resid| in dex."""
    lx, ly = np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float))
    A = np.vstack([lx, np.ones_like(lx)]).T
    sol, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ sol
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    return float(sol[0]), float(np.exp(sol[1])), 1 - ss_res / max(ss_tot, 1e-300), \
        float(np.max(np.abs(ly - pred)) / np.log(10)), (ly - pred) / np.log(10)


def sel(**kw):
    out = [r for r in ROWS if all(r[k] == v for k, v in kw.items()) and r["d"] > 0]
    return sorted(out, key=lambda r: r["n"])


ARMS = [("aware", "p60"), ("generic", "p60"), ("aware", "shift"), ("generic", "shift")]

print("=" * 100)
print("TABLE 1 -- minimal degree d (LP), per arm.  kappa = cond(S) = (1+sig_max)/(1-sig_max)")
print("=" * 100)
for eps in (1e-3, 1e-5):
    print(f"\n  eps = {eps:g}")
    print(f"  {'s':>4} {'n':>3} {'N':>3} {'kappa':>9} | " +
          " ".join(f"{m[:4]}/{c:<5}" for m, c in ARMS) + f" | {'d_model':>9}")
    for s in (1.4, 2.0, 3.0):
        for n in (4, 6, 8, 10, 12, 14, 16, 18):
            got = {}
            for m, c in ARMS:
                rr = [r for r in ROWS if r["s"] == s and r["n"] == n and r["mode"] == m
                      and r["conv"] == c and r["eps"] == eps]
                got[(m, c)] = rr[0] if rr else None
            if not any(got.values()):
                continue
            any_r = next(r for r in got.values() if r)
            cells = []
            for a in ARMS:
                r = got[a]
                cells.append("    --   " if r is None else
                             ("   >max  " if r["d"] < 0 else f"{r['d']:>8d} "))
            print(f"  {s:>4} {n:>3} {2*n:>3} {any_r['kappa']:>9.2f} | " + " ".join(cells) +
                  f" | {any_r['d_model']:>9.0f}")

print()
print("=" * 100)
print("TABLE 2 -- log-log fits.   d = A * kappa^p   (pooled over s; kappa = 4.86 n^2/s^2)")
print("=" * 100)
for eps in (1e-3, 1e-5):
    print(f"\n  eps = {eps:g}")
    print(f"  {'arm':<16} {'#pts':>5} {'kappa range':>16} {'exponent p':>11} {'A':>9} "
          f"{'R^2':>8} {'max|res| dex':>13}")
    for m, c in ARMS:
        rs = [r for r in ROWS if r["mode"] == m and r["conv"] == c and r["eps"] == eps
              and r["d"] > 0]
        if len(rs) < 3:
            continue
        p, A, r2, mx, _ = fit([r["kappa"] for r in rs], [r["d"] for r in rs])
        kr = f"{min(r['kappa'] for r in rs):.0f}-{max(r['kappa'] for r in rs):.0f}"
        print(f"  {m+'/'+c:<16} {len(rs):>5} {kr:>16} {p:>11.3f} {A:>9.3f} {r2:>8.4f} {mx:>13.3f}")

print()
print("=" * 100)
print("TABLE 3 -- log-log fits in N = 2n at FIXED s   (the pre-registered gate variable)")
print("=" * 100)
for eps in (1e-3, 1e-5):
    print(f"\n  eps = {eps:g}")
    print(f"  {'arm':<16} {'s':>5} {'#pts':>5} {'N range':>10} {'exponent':>10} {'R^2':>8} "
          f"{'max|res| dex':>13}")
    for m, c in ARMS:
        for s in (1.4, 2.0, 3.0):
            rs = sel(mode=m, conv=c, eps=eps, s=s)
            if len(rs) < 3:
                continue
            p, A, r2, mx, _ = fit([r["dim"] for r in rs], [r["d"] for r in rs])
            print(f"  {m+'/'+c:<16} {s:>5} {len(rs):>5} "
                  f"{min(r['dim'] for r in rs):>4}-{max(r['dim'] for r in rs):<4} "
                  f"{p:>10.3f} {r2:>8.4f} {mx:>13.3f}")

print()
print("=" * 100)
print("TABLE 4 -- ratios (the deliverable): what spectrum-awareness / the shift buy")
print("=" * 100)
print(f"  {'s':>4} {'n':>3} {'kappa':>9} {'eps':>7} | {'gen/aware (p60)':>16} "
      f"{'p60/shift (gen)':>16} {'model/aware(p60)':>17}")
for eps in (1e-3, 1e-5):
    for s in (1.4, 2.0, 3.0):
        for n in (4, 6, 8, 10, 12, 14, 16, 18):
            def g(m, c):
                rr = [r for r in ROWS if r["s"] == s and r["n"] == n and r["mode"] == m
                      and r["conv"] == c and r["eps"] == eps and r["d"] > 0]
                return rr[0]["d"] if rr else None
            ap, gp, gs = g("aware", "p60"), g("generic", "p60"), g("generic", "shift")
            if ap is None:
                continue
            kap = next(r["kappa"] for r in ROWS if r["s"] == s and r["n"] == n)
            dm = next(r["d_model"] for r in ROWS if r["s"] == s and r["n"] == n
                      and r["eps"] == eps)
            c1 = f"{gp/ap:>16.2f}" if gp else " " * 16
            c2 = f"{gp/gs:>16.2f}" if (gp and gs) else " " * 16
            print(f"  {s:>4} {n:>3} {kap:>9.2f} {eps:>7g} | {c1} {c2} {dm/ap:>17.1f}")

print()
print("=" * 100)
print("TABLE 5 -- rigorous Bernstein lower bounds on the degree (eps -> 0)")
print("=" * 100)
print(f"  {'s':>4} {'n':>3} {'kappa':>9} | " + " ".join(f"{m[:4]}/{c:<5}" for m, c in ARMS)
      + " |  floor ratio gen/aware (p60)")
for s in (1.4, 2.0, 3.0):
    for n in (4, 8, 12, 18):
        lam = spectrum(s, n)
        f = {a: bernstein_floor(lam, *a) for a in ARMS}
        print(f"  {s:>4} {n:>3} {lam[-1]/lam[0]:>9.2f} | " +
              " ".join(f"{f[a]:>9.1f}" for a in ARMS) +
              f" |  {f[('generic','p60')]/f[('aware','p60')]:>6.2f}")

print()
print("=" * 100)
print("TABLE 6 -- per-s fits in kappa, and the asymptotic ratio d_generic/d_aware")
print("=" * 100)
for eps in (1e-3, 1e-5):
    print(f"\n  eps = {eps:g}")
    for m, c in ARMS:
        for s in (1.4, 2.0, 3.0):
            rs = sel(mode=m, conv=c, eps=eps, s=s)
            if len(rs) < 3:
                continue
            p, A, r2, mx, _ = fit([r["kappa"] for r in rs], [r["d"] for r in rs])
            tail = rs[-1]
            print(f"  {m+'/'+c:<16} s={s:<4} pts={len(rs):<3} d = {A:6.3f} kappa^{p:.3f} "
                  f"(R^2={r2:.4f}, max|res|={mx:.3f} dex)   d/kappa at largest kappa "
                  f"({tail['kappa']:.0f}) = {tail['d']/tail['kappa']:.3f}")
print()
print("  Asymptotic large-kappa reading (s=1.4 branch, kappa 51-864):")
for eps in (1e-3,):
    aw = sel(mode="aware", conv="p60", eps=eps, s=1.4)
    ge = [r for r in ROWS if r["mode"] == "generic" and r["conv"] == "p60"
          and r["eps"] == eps and r["d"] > 0]
    pg, Ag, *_ = fit([r["kappa"] for r in ge], [r["d"] for r in ge])
    for r in aw:
        dg = Ag * r["kappa"] ** pg
        print(f"    kappa={r['kappa']:7.1f}  d_aware={r['d']:4d} ({r['d']/r['kappa']:.3f} kappa)"
              f"   d_generic(fit)={dg:7.0f} ({dg/r['kappa']:.2f} kappa)"
              f"   ratio={dg/r['d']:5.2f}   d_model/d_generic={r['d_model']/dg:5.2f}")
