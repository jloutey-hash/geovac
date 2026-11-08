"""Loop test, PT1 and PT3 (record docs/predictions/2026-10-17_density_loop.md, committed ade983ab).

PT1a  the 3-D massless unique triangle  int d^3l / prod |l - x_i|^2 = pi^3 / (|x12| |x23| |x31|),
      checked by the independent Feynman-parameter route
      I = (pi^2 / 2) int_simplex (a1 a2 x12^2 + a1 a3 x13^2 + a2 a3 x23^2)^(-3/2).
PT1b  the reduced formula T = (1/pi) int_{|p-q|<=w<=p+q} rho~_a(p) rho~_b(q) rho~_c(w)
      for one-centre normalised 1s densities rho = zeta^3/pi exp(-2 zeta r)
      (rho~ = 16 zeta^4 / (k^2 + 4 zeta^2)^2, the P3c convention), against
      Remiddi's closed form for f(0,0,0;0,0,0) (Pachucki-Puchalski-Remiddi 2004, eq. 7,
      d^3r/(4 pi) normalisation):  T = 64 prod zeta^3 * (-d_w1)(-d_w2)(-d_w3) f000,  w_i = 2 zeta_i.
PT3   the one-centre values at (1,1,1), (1,2,3), (2,1,3) and a PSLQ rationality test.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import mpmath as mp

OUT = Path(__file__).parent / "data" / "density_loop_A.json"


# ------------------------------------------------------------------ PT1a
def triangle_feynman(x1, x2, x3):
    d2 = lambda u, v: sum((ui - vi) ** 2 for ui, vi in zip(u, v))
    s12, s13, s23 = d2(x1, x2), d2(x1, x3), d2(x2, x3)

    def f(a1, a2):
        a3 = 1 - a1 - a2
        if a3 <= 0:
            return mp.mpf(0)
        return (a1 * a2 * s12 + a1 * a3 * s13 + a2 * a3 * s23) ** mp.mpf(-1.5)

    # simplex as a1 in [0,1], a2 in [0, 1-a1]
    val = mp.quad(lambda a1: mp.quad(lambda a2: f(a1, a2), [0, (1 - a1) / 2, 1 - a1]),
                  [0, mp.mpf(1) / 2, 1])
    return mp.pi ** 2 / 2 * val


def triangle_unique(x1, x2, x3):
    d = lambda u, v: mp.sqrt(sum((ui - vi) ** 2 for ui, vi in zip(u, v)))
    return mp.pi ** 3 / (d(x1, x2) * d(x2, x3) * d(x3, x1))


# ------------------------------------------------------------------ PT1b / PT3
def t_reduced(za, zb, zc):
    """(1/pi) int_0^inf int_0^inf rho_a(p) rho_b(q) [Phi_c(p+q) - Phi_c(|p-q|)],
    Phi_c(w) = int_0^w rho~_c, closed form."""
    za, zb, zc = map(mp.mpf, (za, zb, zc))
    rt = lambda k, z: 16 * z ** 4 / (k * k + 4 * z * z) ** 2
    g = 2 * zc

    def Phi(w):
        return 16 * zc ** 4 * (w / (2 * g ** 2 * (w * w + g * g)) + mp.atan(w / g) / (2 * g ** 3))

    def inner(p):
        h = lambda q: rt(q, zb) * (Phi(p + q) - Phi(abs(p - q)))
        return mp.quad(h, [0, p / 2, p, 2 * p + 1, mp.inf]) if p > 0 else mp.mpf(0)

    val = mp.quad(lambda p: rt(p, za) * inner(p), [0, mp.mpf(1) / 2, 2, 8, mp.inf])
    return val / mp.pi


def f000(w1, w2, w3):
    """Remiddi's master, as read from PPR 2004 eq. (7)."""
    def term(a, b, c):           # a = the singled-out w, b + c the pair
        x = a / (b + c)
        return mp.log(x) * mp.log(1 + x) + mp.polylog(2, -x) + mp.polylog(2, 1 - x)
    return -(term(w3, w1, w2) + term(w2, w3, w1) + term(w1, w2, w3)) / (2 * w1 * w2 * w3)


def t_remiddi(za, zb, zc):
    w = [2 * mp.mpf(z) for z in (za, zb, zc)]
    d3 = mp.diff(lambda a, b, c: f000(a, b, c), tuple(w), (1, 1, 1))
    return 64 * (mp.mpf(za) * zb * zc) ** 3 * (-d3)


def main(dps: int) -> None:
    mp.mp.dps = dps
    rec = {"dps": dps, "pt1a": [], "pt1b_pt3": {}}
    random.seed(17)
    for _ in range(5):
        xs = [[mp.mpf(random.uniform(-2, 2)) for _ in range(3)] for _ in range(3)]
        a, b = triangle_feynman(*xs), triangle_unique(*xs)
        rec["pt1a"].append({"feynman": mp.nstr(a, 20), "unique": mp.nstr(b, 20),
                            "rel": mp.nstr(abs(a / b - 1), 3)})
        print("PT1a", mp.nstr(a, 15), mp.nstr(b, 15), "rel", mp.nstr(abs(a / b - 1), 3), flush=True)
    for zs in [(1, 1, 1), (1, 2, 3), (2, 1, 3)]:
        tr = t_reduced(*zs)
        trm = t_remiddi(*zs)
        rel = abs(tr / trm - 1)
        rat = mp.identify(tr) if dps <= 40 else None
        pslq = mp.pslq([tr, mp.mpf(1)], maxcoeff=10 ** 12, maxsteps=10 ** 6)
        rec["pt1b_pt3"][str(zs)] = {"reduced": mp.nstr(tr, dps - 5), "remiddi": mp.nstr(trm, dps - 5),
                                    "rel": mp.nstr(rel, 3), "pslq_Q": str(pslq)}
        print("PT1b/PT3", zs, mp.nstr(tr, 25), mp.nstr(trm, 25), "rel", mp.nstr(rel, 3),
              "pslq[T,1]", pslq, flush=True)
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30)
