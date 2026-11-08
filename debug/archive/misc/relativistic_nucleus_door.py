"""Relativistic finite nucleus: tests RD1-RD5 of docs/predictions/2026-10-18_relativistic_nucleus.md.

Registered in commit aef8ba6d before this script ran. Writes debug/data/relativistic_nucleus_door.json.

Dirac radial equations for G = r g, F = r f, total energy W (rest energy included), atomic units:
  dG/dr = -(kappa/r) G + [(W - V + c^2)/c] F,     dF/dr = (kappa/r) F - [(W - V - c^2)/c] G.
Inside a uniformly charged sphere V = -V0 + B r^2 (V0 = 3Z/(2R), B = Z/(2R^3)): exact power series.
Outside V = -Z/r: integrated inward in s = ln r from the decaying solution, matched at r = R.
"""
from __future__ import annotations

import json
from typing import Dict, List

import mpmath as mp

C = mp.mpf("137.035999084")
OUT: Dict[str, object] = {}


def interior_ratio(W, kappa: int, Z, R, tol_digits: int = 40):
    """F/G at r = R for the regular interior solution, by the power-series recurrence."""
    V0 = 3 * Z / (2 * R)
    B = Z / (2 * R**3)
    P = W + V0 + C**2
    M = W + V0 - C**2
    N = 400
    a = [mp.mpf(0)] * (N + 4)
    b = [mp.mpf(0)] * (N + 4)
    # leading index: kappa < 0 -> G ~ r^(-kappa)... use the k where (k + kappa) = 0 (G) or (k - kappa) = 0 (F)
    if kappa < 0:
        k0 = -kappa
        a[k0] = mp.mpf(1)
    else:
        k0 = kappa
        b[k0] = mp.mpf(1)
    for k in range(k0 + 1, N):
        bm1 = b[k - 1]
        bm3 = b[k - 3] if k >= 3 else 0
        am1 = a[k - 1]
        am3 = a[k - 3] if k >= 3 else 0
        if k + kappa != 0:
            a[k] = (P * bm1 - B * bm3) / C / (k + kappa)
        if k - kappa != 0:
            b[k] = -(M * am1 - B * am3) / C / (k - kappa)
    G = mp.fsum(a[k] * R**k for k in range(N))
    F = mp.fsum(b[k] * R**k for k in range(N))
    last = abs(a[N - 1] * R ** (N - 1)) + abs(b[N - 1] * R ** (N - 1))
    assert last < mp.mpf(10) ** (-tol_digits) * (abs(G) + abs(F)), "interior series not converged"
    return F / G


def exterior_start(W, kappa: int, Z):
    """Starting point (s_max, G, F) on the decaying exterior solution."""
    if W > -C**2:
        lam = mp.sqrt(C**4 - W**2) / C
        r_max = 45 / lam
        G0, F0 = mp.mpf(1), -mp.sqrt((C**2 - W) / (C**2 + W))
    else:
        r_max = (mp.mpf(60) ** 2) / (8 * Z)
        dlogG = -mp.sqrt(2 * Z / r_max) - mp.mpf(1) / (4 * r_max)
        G0 = mp.mpf(1)
        F0 = (C * r_max / Z) * (dlogG + kappa / r_max) * G0
    return mp.log(r_max), G0, F0


def exterior_rhs(W, kappa: int, Z):
    def rhs(t, Y):  # t = -s
        r = mp.e ** (-t)
        G, F = Y
        dG = -kappa * G + ((W + C**2) * r + Z) / C * F
        dF = kappa * F - ((W - C**2) * r + Z) / C * G
        return [-dG, -dF]
    return rhs


def exterior_segment(W, kappa: int, Z, s_from, Y_from, s_to):
    """Integrate from s_from down to s_to with a fresh solver, so memory stays bounded."""
    sol = mp.odefun(exterior_rhs(W, kappa, Z), -s_from, list(Y_from))
    return sol(-s_to)


def exterior_on_grid(W, kappa: int, Z, s_grid):
    """Exterior (G, F) at each s in s_grid (descending), integrated segment by segment and rescaled."""
    s0, G0, F0 = exterior_start(W, kappa, Z)
    out = []
    s_prev, Y = s0, [G0, F0]
    for sv in s_grid:
        Y = exterior_segment(W, kappa, Z, s_prev, Y, sv)
        nrm = mp.sqrt(Y[0] ** 2 + Y[1] ** 2)
        Y = [Y[0] / nrm, Y[1] / nrm]
        out.append(Y)
        s_prev = sv
    return out


def exterior_solution(W, kappa: int, Z, s_min):
    """Integrate the decaying exterior solution inward; returns a function s -> (G, F)."""
    if W > -C**2:
        lam = mp.sqrt(C**4 - W**2) / C
        r_max = 45 / lam
        ratio = -mp.sqrt((C**2 - W) / (C**2 + W))  # F/G at large r
        G0, F0 = mp.mpf(1), ratio
    else:  # W = -c^2: G ~ exp(-sqrt(8 Z r)); F from dG/dr = -kappa G/r + Z F/(c r)
        r_max = (mp.mpf(60) ** 2) / (8 * Z)
        dlogG = -mp.sqrt(2 * Z / r_max) - mp.mpf(1) / (4 * r_max)
        G0 = mp.mpf(1)
        F0 = (C * r_max / Z) * (dlogG + kappa / r_max) * G0
    s_max = mp.log(r_max)

    def rhs(t, Y):  # t = -s, integrate forward in t from -s_max
        s = -t
        r = mp.e**s
        G, F = Y
        dG = -kappa * G + ((W + C**2) * r + Z) / C * F
        dF = kappa * F - ((W - C**2) * r + Z) / C * G
        return [-dG, -dF]

    sol = mp.odefun(rhs, -s_max, [G0, F0])
    return lambda s: sol(-s)


def match(W, kappa, Z, R, ext=None):
    if ext is None:
        ext = exterior_solution(W, kappa, Z, mp.log(R))
    Gs, Fs = ext(mp.log(R))
    # Wronskian form (smooth in W, R and Z; no poles): interior normalized to G = 1
    return (interior_ratio(W, kappa, Z, R) * Gs - Fs) / mp.sqrt(Gs**2 + Fs**2)


def W_point(n: int, kappa: int, Z):
    g = mp.sqrt(kappa**2 - (Z / C) ** 2)
    nr = n - abs(kappa)
    return C**2 / mp.sqrt(1 + (Z / C) ** 2 / (nr + g) ** 2)


def W_bound(kappa, Z, R, guess):
    f = lambda W: match(W, kappa, Z, R)
    eps = (C**2 - guess) * mp.mpf(10) ** -14 + mp.mpf(10) ** -10
    f0, f1 = f(guess), f(guess + eps)
    W1 = guess - f0 * eps / (f1 - f0)  # D is smooth and nearly linear: one secant step, then refine
    return mp.findroot(f, (W1, W1 + eps))


def W_bound_bracket(kappa, Z, R, hi, lo_frac=mp.mpf("0.05")):
    """Bracketing root search below hi (used where the finite-size shift is not small, Z alpha = 1)."""
    f = lambda W: match(W, kappa, Z, R)
    fhi = f(hi)
    W = hi
    step = (hi - lo_frac * hi) / 40
    while W > lo_frac * hi:
        Wn = W - step
        fn = f(Wn)
        if mp.sign(fn) != mp.sign(fhi):
            return mp.findroot(f, (Wn, W), solver="illinois")
        W, fhi = Wn, fn
    raise RuntimeError("no sign change")


def bracket_root(f, a, fa, b, fb, steps: int = 40, xtol=mp.mpf(10) ** -12):
    """Illinois regula falsi with a fixed step cap and a progress print (mpmath's own stalled here)."""
    side = 0
    for _ in range(steps):
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)
        print("  root step", mp.nstr(c, 15), mp.nstr(fc, 3), flush=True)
        if fc == 0 or abs(b - a) < xtol:
            return c
        if mp.sign(fc) == mp.sign(fb):
            b, fb = c, fc
            if side == -1:
                fa /= 2
            side = -1
        else:
            a, fa = c, fc
            if side == 1:
                fb /= 2
            side = 1
    return c


def run(stages) -> None:
    mp.mp.dps = 25
    if "rd2" in stages:
        stage_rd2()
    if "rd3" in stages:
        stage_rd3()
    if "rd4" in stages:
        stage_rd4()
    if "rd5" in stages:
        stage_rd5()
    path = "debug/data/relativistic_nucleus_door.json"
    try:
        with open(path, encoding="utf-8") as fh:
            old = json.load(fh)
    except (OSError, ValueError):
        old = {}
    old.update(OUT)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(old, fh, indent=1, default=str)


def stage_rd2() -> None:
    # RD check: Sommerfeld limit and RD2 exponent
    rd2: Dict[str, object] = {}
    for za in (mp.mpf("0.3"), mp.mpf("0.6"), mp.mpf("0.9")):
        Z = za * C
        g = mp.sqrt(1 - za**2)
        W0 = W_point(1, -1, Z)
        Rs = [mp.mpf(10) ** (-e) for e in (8, 7, 6, 5)]
        dW = []
        for R in Rs:
            W = W_bound(-1, Z, R, W0)
            dW.append(W - W0)
        slopes = [mp.log(dW[i + 1] / dW[i]) / mp.log(Rs[i + 1] / Rs[i]) for i in range(len(Rs) - 1)]
        rd2[f"Zalpha={za}"] = {"2gamma": mp.nstr(2 * g, 10), "slopes_small_to_large_R": [mp.nstr(sl, 10) for sl in slopes],
                               "dW_at_R=1e-8": mp.nstr(dW[0], 6)}
        print("RD2", za, rd2[f"Zalpha={za}"], flush=True)
    # p 1/2 (kappa = +1): RD1, the level exists and shifts with R
    Z = mp.mpf("0.6") * C
    W0p = W_point(2, 1, Z)
    Wp = W_bound(1, Z, mp.mpf(10) ** -5, W0p)
    rd2["2p1/2_Zalpha=0.6_R=1e-5_shift"] = mp.nstr(Wp - W0p, 6)
    OUT["RD2"] = rd2


def stage_rd3() -> None:
    # RD3: Z alpha = 1
    Z = C
    rd3 = {}
    Rs = [mp.mpf(10) ** (-e) for e in (3, 4, 5, 6, 7, 8, 9)]
    Ws: List = []
    guess = C**2 * mp.mpf("0.6")  # the bracket search walks down from here to the first sign change
    for R in Rs:
        W = W_bound_bracket(-1, Z, R, guess)
        Ws.append(W)
        guess = W
        print("RD3", mp.nstr(R, 3), mp.nstr(W / C**2, 12), flush=True)
    rd3["W_over_c2"] = {mp.nstr(R, 3): mp.nstr(W / C**2, 12) for R, W in zip(Rs, Ws)}
    rd3["local_slopes_dlnW_dlnR"] = [mp.nstr(mp.log(Ws[i + 1] / Ws[i]) / mp.log(Rs[i + 1] / Rs[i]), 6) for i in range(len(Rs) - 1)]
    rd3["W*ln(1/R)/c2 (not registered)"] = [mp.nstr(W * mp.log(1 / R) / C**2, 6) for R, W in zip(Rs, Ws)]
    OUT["RD3"] = rd3


def stage_rd4() -> None:
    # RD4: critical charge with R = 1.2 (2.5 Z)^(1/3) fm
    bohr_fm = mp.mpf("52917.72109")
    Rz = lambda Z: mp.mpf("1.2") * (mp.mpf("2.5") * Z) ** (mp.mpf(1) / 3) / bohr_fm
    Wc = -C**2
    Dz = lambda Z: match(Wc, -1, Z, Rz(Z))
    scan = {}
    prev = None
    zc = None
    for Zv in range(160, 186, 2):
        d = Dz(mp.mpf(Zv))
        scan[Zv] = mp.nstr(d, 5)
        print("RD4 scan", Zv, scan[Zv], flush=True)
        if prev is not None and mp.sign(d) != mp.sign(prev[1]) and zc is None:
            zc = bracket_root(Dz, mp.mpf(prev[0]), prev[1], mp.mpf(Zv), d)
        prev = (Zv, d)
    OUT["RD4"] = {"scan_D(Z)": scan, "Z_cr": mp.nstr(zc, 8) if zc is not None else None, "D(Z_cr)": mp.nstr(Dz(zc), 3) if zc is not None else None,
                  "R(Z_cr)_bohr": mp.nstr(Rz(zc), 6) if zc is not None else None}
    print("RD4", OUT["RD4"], flush=True)


def stage_rd5() -> None:
    # RD5: Z alpha = 1.2, W = -c^2: zeros of the matching function in R.
    # The exterior is integrated one grid cell at a time (bounded memory); a root inside a cell is
    # refined by integrating only across that cell.
    mp.mp.dps = 30
    Wc = -C**2
    za = mp.mpf("1.2")
    Z = za * C
    sigma = mp.sqrt(za**2 - 1)
    grid = [mp.mpf(10) ** (-mp.mpf(e) / 8) for e in range(8, 8 * 14 + 1)]  # R = 1e-1 .. 1e-14
    sg = [mp.log(R) for R in grid]
    ext = exterior_on_grid(Wc, -1, Z, sg)
    D = lambda R, Y: (interior_ratio(Wc, -1, Z, R) * Y[0] - Y[1]) / mp.sqrt(Y[0] ** 2 + Y[1] ** 2)
    vals = [D(R, Y) for R, Y in zip(grid, ext)]
    zeros = []
    for i in range(len(grid) - 1):
        if mp.sign(vals[i]) != mp.sign(vals[i + 1]):
            def f(sv, i=i):
                Y = exterior_segment(Wc, -1, Z, sg[i], ext[i], sv)
                return D(mp.e**sv, Y)
            sv = mp.findroot(f, (sg[i], sg[i + 1]), solver="illinois", tol=mp.mpf(10) ** -16, verify=False)
            zeros.append(mp.e**sv)
            print("RD5 zero", mp.nstr(zeros[-1], 8), flush=True)
    ratios = [zeros[k] / zeros[k + 1] for k in range(len(zeros) - 1)]
    OUT["RD5"] = {"sigma": mp.nstr(sigma, 8), "e^(pi/sigma)": mp.nstr(mp.e ** (mp.pi / sigma), 8),
                  "zeros_R": [mp.nstr(z_, 6) for z_ in zeros], "ratios": [mp.nstr(q, 6) for q in ratios]}
    print("RD5", OUT["RD5"], flush=True)


if __name__ == "__main__":
    import sys
    run(sys.argv[1:] or ["rd2", "rd3", "rd4", "rd5"])
