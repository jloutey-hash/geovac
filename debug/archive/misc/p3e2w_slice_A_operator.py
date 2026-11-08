"""P3e-3 step A (pre-registered: docs/predictions/2026-10-19_p3e2_widened.md, commit bb3ffe11).

Slice G(t;1,1) = int_0^inf x K1(t x) J(x)^2 dx, J = 3L - 3 x K0 - x^2 K1 (Slater weight), t on the
central density.  Mellin route of density_loop_B_ode.py (imported, not modified):
  theta J = x^3 K2, ((theta-3)^2 - 4 - x^2) theta J = 0  =>  theta^3 J = (z - 5) theta J + 6 theta^2 J,
  z = x^2;  J^2 obeys the symmetric square;  Phi(s) = M[G](s) = A(2 - s) E(s),  A = M[J^2],
  E = M[K1], E(s+2) = (s^2 - 1) E(s).
Self-check: the same code reproduces F_3's operator 4 th^3 - t^2 (th+1)^3 for int K0(x)^2 K0(tx) dx.
Numerical check: L G at several t, theta_t^j applied under the integral (theta^2 K1 = (1+y^2) K1).
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_B_ode as B  # noqa: E402
from three_electron_p3e_check import L as Lfun  # noqa: E402

u, y, z = B.u, B.y, B.z
OUT = Path(__file__).parent / "data" / "p3e2w_slice_A.json"


def Jf(x):
    return 3 * Lfun(x) - 3 * x * mp.besselk(0, x) - x * x * mp.besselk(1, x)


def LGnum(Q, t, Jcache):
    """sum_i t^{2i} Q_i(theta_t) G(t), theta_t^j K1(t x) = (theta_y^j K1)(t x)."""
    t = mp.mpf(t)
    J = max(sp.Poly(q, u).degree() for q in Q)
    al, be, _ = B.theta_powers(J, 1 + y**2, 0)
    coefs = []
    for i, q in enumerate(Q):
        pc = sp.Poly(sp.expand(q), u).all_coeffs()[::-1]
        A = sum(sp.Rational(c) * al[j] for j, c in enumerate(pc))
        Bb = sum(sp.Rational(c) * be[j] for j, c in enumerate(pc))
        coefs.append((i, sp.lambdify(y, sp.expand(A), "mpmath"), sp.lambdify(y, sp.expand(Bb), "mpmath")))

    def h(x):
        yy = t * x
        k0, k1 = mp.besselk(0, yy), mp.besselk(1, yy)
        th1 = -yy * k0 - k1                       # theta K1 = y K1' = -y K0 - K1
        val = sum(t ** (2 * i) * (a(yy) * k1 + b(yy) * th1) for i, a, b in coefs)
        return x * Jf(x) ** 2 * val
    return mp.quad(h, [0, mp.mpf(1) / 8, 1, 4, 16, 64, mp.inf])


def main(dps: int) -> None:
    rec = {}
    vK0 = B.sym_power([z, 0], 2, 2)
    QF3 = B.mellin_t_ode(vK0, 1, 0, lambda sg: sg**2)
    print("self-check F3:", QF3, flush=True)
    assert [sp.expand(q) for q in QF3] in ([sp.expand(4 * u**3), sp.expand(-(u + 1) ** 3)],
                                           [sp.expand(-4 * u**3), sp.expand((u + 1) ** 3)])
    vJ = B.sym_power([0, z - 5, 6], 3, 2)
    print("J^2 relation:", vJ, flush=True)
    QG = B.mellin_t_ode(vJ, 2, 0, lambda sg: sg**2 - 1)
    print("L_G order", max(sp.Poly(q, u).degree() for q in QG), " t-degree", 2 * (len(QG) - 1), flush=True)
    for i, q in enumerate(QG):
        print(f"  t^{2*i}:", q, flush=True)
    lead = sum(sp.Poly(q, u).coeff_monomial(u ** max(sp.Poly(qq, u).degree() for qq in QG)) * sp.Symbol("t") ** (2 * i)
               for i, q in enumerate(QG))
    print("  leading coefficient:", sp.factor(lead), flush=True)
    rec["QG"] = [str(q) for q in QG]
    rec["lead"] = str(sp.factor(lead))
    mp.mp.dps = dps
    rec["check"] = {}
    for ts in ["0.5", "0.9", "1.3", "2.5", "3.5"]:
        t0 = time.time()
        r = LGnum(QG, ts, None)
        g = LGnum([sp.Integer(1)], ts, None)
        rec["check"][ts] = {"LG": mp.nstr(r, 15), "G": mp.nstr(g, 25)}
        print(f"  t={ts}: G = {mp.nstr(g, 20)}   L_G G = {mp.nstr(r, 8)}  ({time.time()-t0:.0f}s)", flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20)
