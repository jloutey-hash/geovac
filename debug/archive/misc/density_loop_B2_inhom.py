"""Loop test PT2, step B2: the Mellin-route operators of step B annihilate G and C only up to
terms the regularized Mellin transform cannot see (pure powers t^lam log^k t). Measure L G and L C
at several t with a cancellation-free theta^j psi (asymptotic series at large argument) and fit
the remainder to powers.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import density_loop_B_ode as B  # noqa: E402

OUT = Path(__file__).parent / "data" / "density_loop_B2.json"
u, y = B.u, B.y


def theta_j_psi_factory(J):
    al, be, ga = B.theta_powers(J, 1 + y**2, -2 / sp.pi)
    lam = [tuple(sp.lambdify(y, e, "mpmath") for e in (al[j], be[j], ga[j])) for j in range(J + 1)]

    def th(yy, j):
        yy = mp.mpf(yy)
        if yy < B.X0:
            with mp.workdps(mp.mp.dps + int(yy / 2.2) + 15 + 6 * j):
                a, b, c = lam[j]
                return +(a(yy) * B.psi(yy) + b(yy) * B.theta_psi(yy) + c(yy))
        tot, m = mp.mpf(0), 0
        while True:
            cc = (-1) ** m * mp.factorial(2 * m) * mp.binomial(mp.mpf(-1.5), m) * yy ** (-2 * m - 2) * (-2 * m - 2) ** j
            tot += cc
            if abs(cc) < mp.eps * abs(tot) or m > yy / 2:
                return 2 / mp.pi * tot
            m += 1
    return th


def LG(Q, t, th):
    """sum_i t^{2i} Q_i(theta_t) G, theta_t^j [t^2 psi(tx)] = t^2 [(theta+2)^j psi](tx)."""
    t = mp.mpf(t)
    J = max(sp.Poly(q, u).degree() for q in Q)
    # operator on the t-factor: sum_i t^{2i} Q_i(theta+2) applied to psi at y = t x
    coef = {}
    for i, q in enumerate(Q):
        pc = sp.Poly(sp.expand(q.subs(u, u + 2)), u).all_coeffs()[::-1]
        coef[i] = [mp.mpf(sp.Rational(c).p) / sp.Rational(c).q for c in pc]

    def h(x):
        yy = t * x
        val = sum(t ** (2 * i) * sum(c * th(yy, j) for j, c in enumerate(cs)) for i, cs in coef.items())
        return x * x * B.psi(x) ** 2 * val
    pts = [0, mp.mpf(1) / 4, 1, 4, 16, 64, B.X0, 4 * B.X0, mp.inf]
    return mp.pi / 2 * t * t * mp.quad(h, pts)


def LC(Q, t):
    t = mp.mpf(t)
    J = max(sp.Poly(q, u).degree() for q in Q)
    al, be, _ = B.theta_powers(J, y**2, 0)
    coefs = {}
    for i, q in enumerate(Q):
        pc = sp.Poly(sp.expand(q.subs(u, u + 3)), u).all_coeffs()[::-1]
        A = sum(sp.Rational(c) * al[j] for j, c in enumerate(pc))
        Bb = sum(sp.Rational(c) * be[j] for j, c in enumerate(pc))
        coefs[i] = (sp.lambdify(y, sp.expand(A), "mpmath"), sp.lambdify(y, sp.expand(Bb), "mpmath"))

    def h(x):
        yy = t * x
        k0, k1 = mp.besselk(0, yy), mp.besselk(1, yy)
        val = sum(t ** (2 * i) * (a(yy) * k0 - b(yy) * yy * k1) for i, (a, b) in coefs.items())
        return B.chi(x) ** 2 * val
    return 8 * t**3 / mp.pi**3 * mp.quad(h, [0, mp.mpf(1) / 4, 1, 4, 16, 64, mp.inf])


def main(dps: int) -> None:
    mp.mp.dps = dps
    QG = [sp.sympify(q) for q in json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())["QG"]]
    QC = [sp.sympify(q) for q in json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())["QC"]]
    th = theta_j_psi_factory(12)
    ts = ["0.5", "0.7", "1.3", "1.6", "2.5", "3.5"]
    rec = {"dps": dps, "LG": {}, "LC": {}}
    for t in ts:
        g, c = LG(QG, t, th), LC(QC, t)
        rec["LG"][t], rec["LC"][t] = mp.nstr(g, dps - 3), mp.nstr(c, dps - 3)
        print("t", t, " L_G G =", mp.nstr(g, 18), "  L_C C =", mp.nstr(c, 18), flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
