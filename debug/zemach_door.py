"""Pre-registered ZM1-ZM4: the hyperfine door, a finite magnetization distribution.

Record: docs/predictions/2026-10-19_zemach_door.md (committed be16fccd before this script ran).

Non-relativistic s states, Z = 1, x = R. Hyperfine factor H_n = int rho_M |psi_n|^2 / |psi_n^point(0)|^2.
  zm2:   point charge, uniform-ball magnetization, exact sympy, n = 1..6.
  zmnum: uniform-sphere charge; magnetization point / ball / shell; 40 digits; fits in x.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp
import sympy.physics.hydrogen  # noqa: F401

OUT = Path(__file__).resolve().parent / "data" / "zemach_door.json"


# ---------------- Zemach moments (exact) ----------------
def zemach_moments():
    r1, r2, R = sp.symbols("r1 r2 R", positive=True)
    # angular average of |r1 - r2| for fixed radii: r_> + r_<^2 / (3 r_>)
    def mean_dist_to_ball(r):  # point at radius r <= R to a uniform ball of radius R
        inner = sp.integrate(3 * r2**2 / R**3 * (r + r2**2 / (3 * r)), (r2, 0, r))
        outer = sp.integrate(3 * r2**2 / R**3 * (r2 + r**2 / (3 * r2)), (r2, r, R))
        return sp.simplify(inner + outer)
    d1 = mean_dist_to_ball(r1)
    ball_ball = sp.simplify(sp.integrate(3 * r1**2 / R**3 * d1, (r1, 0, R)) / R)
    ball_point = sp.simplify(sp.integrate(3 * r1**2 / R**3 * r1, (r1, 0, R)) / R)
    ball_shell = sp.simplify(sp.integrate(3 * r1**2 / R**3 * (R + r1**2 / (3 * R)), (r1, 0, R)) / R)
    return {"ball_point": ball_point, "ball_ball": ball_ball, "ball_shell": ball_shell}


# ---------------- ZM2: point charge, exact ----------------
def stage_zm2(nmax=6, order=6):
    r, x = sp.symbols("r x", positive=True)
    out = {}
    for n in range(1, nmax + 1):
        Rn = sp.simplify(sp.physics.hydrogen.R_nl(n, 0, r, 1))
        dens = sp.simplify(Rn**2 / Rn.subs(r, 0)**2)
        H = 3 / x**3 * sp.integrate(r**2 * dens, (r, 0, x))
        ser = sp.series(sp.simplify(H), x, 0, order + 1).removeO()
        coeffs = [sp.nsimplify(ser.coeff(x, k)) for k in range(0, order + 1)]
        has_log = ser.has(sp.log)
        out[n] = {"coeffs": [str(c) for c in coeffs], "has_log": bool(has_log),
                  "x1_ok": coeffs[1] == sp.Rational(-3, 2),
                  "x2_ok": sp.simplify(coeffs[2] - (1 + sp.Rational(1, 5) / n**2)) == 0}
        print("ZM2 n", n, out[n], flush=True)
    # polynomial in 1/n^2: c_k(n) fitted by degree floor(k/2) in t = 1/n^2 must reproduce every n
    t = sp.symbols("t")
    poly_ok = {}
    for k in range(1, order + 1):
        deg = k // 2
        pts = [(sp.Rational(1, n**2), sp.Rational(out[n]["coeffs"][k])) for n in range(1, nmax + 1)]
        P = sp.interpolate(pts[:deg + 1], t)
        poly_ok[k] = {"degree": deg, "polynomial": str(sp.expand(P)),
                      "reproduces_all_n": all(sp.simplify(P.subs(t, tt) - v) == 0 for tt, v in pts)}
    return {"by_n": out, "polynomial_in_1_over_n2": poly_ok}


# ---------------- numerics: uniform-sphere charge ----------------
def interior(R, E, tol):
    """u = sum a_k r^k, odd k, for u'' = (-3/R - 2E + r^2/R^3) u, a_1 = 1. Returns u(R), u'(R), int_0^R u^2."""
    c0 = -3 / R - 2 * E
    a = {0: mp.mpf(0), 1: mp.mpf(1), 2: mp.mpf(0), 3: c0 / 6}
    k = 3
    while True:
        k += 1
        a[k] = (c0 * a[k - 2] + a[k - 4] / R**3) / (k * (k - 1))
        if k > 12 and k % 2 == 1 and abs(a[k] * R**k) < tol * R and abs(a[k - 2] * R**(k - 2)) < tol * R:
            break
    ks = sorted(a)
    uR = mp.fsum(a[j] * R**j for j in ks)
    dR = mp.fsum(j * a[j] * R**(j - 1) for j in ks if j)
    odd = [j for j in ks if j % 2 == 1]
    sq = mp.fsum(a[i] * a[j] * R**(i + j + 1) / (i + j + 1) for i in odd for j in odd)
    return uR, dR, sq


def W(nu, z):
    return mp.whitw(nu, mp.mpf(1) / 2, z)


def ext_logderiv(nu, R):
    z = 2 * R / nu
    w0, w1 = W(nu, z), W(nu + 1, z)
    dz = ((z / 2 - nu) * w0 - w1) / z
    return (2 / nu) * dz / w0, w0


def match(nu, R, tol):
    E = -1 / (2 * nu**2)
    uR, dR, _ = interior(R, E, tol)
    ld, _ = ext_logderiv(nu, R)
    return dR / uR - ld


def solve_nu(n, R, tol):
    a = mp.mpf(n) + mp.mpf(2) / 5 * R**2 / 1   # quantum defect guess, nu = n - mu
    b = a + mp.mpf(10)**(-6) * R**2 + mp.mpf(10)**(-25)
    fa, fb = match(a, R, tol), match(b, R, tol)
    for it in range(60):
        c = b - fb * (b - a) / (fb - fa)
        fc = match(c, R, tol)
        a, fa, b, fb = b, fb, c, fc
        if abs(b - a) < mp.mpf(10)**(-(mp.mp.dps - 6)):
            break
    return b


def bound_state(n, R, tol):
    nu = solve_nu(n, R, tol)
    E = -1 / (2 * nu**2)
    uR, dR, sq_in = interior(R, E, tol)
    _, wR = ext_logderiv(nu, R)
    c = uR / wR
    pts = [R, 2 * R, 10 * R, mp.mpf(1) / 10, 1, 3 * nu, 10 * nu, 30 * nu, 80 * nu]
    pts = sorted(set(p for p in pts if p >= R))
    sq_out = c**2 * mp.quad(lambda r: W(nu, 2 * r / nu)**2, pts + [mp.inf])
    Nrm = sq_in + sq_out
    return {"nu": nu, "uR": uR, "sq_in": sq_in, "N": Nrm}


def hyperfine(n, R, bs):
    n3 = mp.mpf(n)**3
    return {"point": n3 / (4 * bs["N"]) - 1,
            "ball": 3 * n3 * bs["sq_in"] / (4 * R**3 * bs["N"]) - 1,
            "shell": n3 * bs["uR"]**2 / (4 * R**2 * bs["N"]) - 1}


def lsq(xs, ys, cols):
    A = mp.matrix([[f(x) for f in cols] for x in xs])
    b = mp.matrix(ys)
    sol, res = mp.qr_solve(A, b)
    return [sol[i] for i in range(len(cols))], res


COLS = {
    "x": lambda x: x,
    "x2lnx": lambda x: x**2 * mp.log(x),
    "x2": lambda x: x**2,
    "x3lnx": lambda x: x**3 * mp.log(x),
    "x3": lambda x: x**3,
    "x4": lambda x: x**4,
    "x4lnx": lambda x: x**4 * mp.log(x),
}


def stage_num(ns=(1, 2), npts=12, dps=40):
    mp.mp.dps = dps
    tol = mp.mpf(10)**(-(dps + 5))
    xs = [mp.mpf(10)**(-5) * (mp.mpf(300))**(mp.mpf(i) / (npts - 1)) for i in range(npts)]
    raw, fits = {}, {}
    for n in ns:
        rows = []
        for x in xs:
            bs = bound_state(n, x, tol)
            h = hyperfine(n, x, bs)
            rows.append((x, h, bs["nu"]))
            print("ZM n", n, "x", mp.nstr(x, 4), "mu", mp.nstr(n - bs["nu"], 8),
                  {k: mp.nstr(v, 12) for k, v in h.items()}, flush=True)
        raw[n] = [{"x": mp.nstr(x, 20), "mu": mp.nstr(n - nu, 25),
                   **{k: mp.nstr(v, 30) for k, v in h.items()}} for x, h, nu in rows]
        for prof in ("point", "ball", "shell"):
            ys = [h[prof] for _, h, _ in rows]
            out = {}
            for label, names in {"registered": ["x", "x2lnx", "x2", "x3lnx", "x3", "x4"],
                                 "drop_last": ["x", "x2lnx", "x2", "x3lnx", "x3"],
                                 "plus_x4lnx": ["x", "x2lnx", "x2", "x3lnx", "x3", "x4lnx", "x4"],
                                 "no_logs": ["x", "x2", "x3", "x4"]}.items():
                c, res = lsq(xs, ys, [COLS[k] for k in names])
                out[label] = {"coeffs": {k: mp.nstr(v, 12) for k, v in zip(names, c)}, "residual": mp.nstr(res, 3)}
            fits[f"n{n}_{prof}"] = out
            print(prof, out["registered"], flush=True)
    return {"raw": raw, "fits": fits}


def run(stages):
    data = json.loads(OUT.read_text()) if OUT.exists() else {}
    if "moments" in stages:
        data["zemach_moments_over_R"] = {k: str(v) for k, v in zemach_moments().items()}
        print(data["zemach_moments_over_R"])
    if "zm2" in stages:
        data["ZM2"] = stage_zm2()
    if "num" in stages:
        data["ZM1_ZM3"] = stage_num()
    OUT.write_text(json.dumps(data, indent=1, default=str))


if __name__ == "__main__":
    run(sys.argv[1:] or ["moments", "zm2", "num"])
