"""FP step P: numerical check of B P + q(t^2) P_h4 = 0 (step O2), P_h4 = (16/pi^4) t^3 int K0(t x) x^6 K0(x)^3 dx
(by quadrature), theta^j P from the step-E cache; and A4 P_h4 = 0 (Mellin operator of P_h4)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_M_eval as ME  # noqa: E402
import four_density_star_ops as O  # noqa: E402
u, t = O.u, sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_P.json"


def Ph4_thetas(tv, J):
    tv = mp.mpf(tv)
    al, be = [sp.Integer(1)], [sp.Integer(0)]
    y = sp.Symbol("y")
    for _ in range(J):
        a, b = al[-1], be[-1]
        al.append(sp.expand(y * sp.diff(a, y) + b * y**2)); be.append(sp.expand(a + y * sp.diff(b, y)))
    I = []
    for j in range(J + 1):
        fa, fb = sp.lambdify(y, al[j], "mpmath"), sp.lambdify(y, be[j], "mpmath")
        f = lambda x: x**6 * mp.besselk(0, x)**3 * (fa(tv*x) * mp.besselk(0, tv*x) - fb(tv*x) * tv*x * mp.besselk(1, tv*x))
        I.append(mp.quad(f, [0, mp.mpf(1)/4, 1, 4, 16, 64]))
    N = 16 / mp.pi**4
    return [N * tv**3 * sum(mp.binomial(k, m) * 3**(k-m) * I[m] for m in range(k+1)) for k in range(J+1)]


def main():
    dps = 60
    d = json.loads((HERE / "data" / "four_density_path_O2.json").read_text())
    BQ = [sp.sympify(q).subs(sp.Symbol("u"), u) for q in d["B_Q"]]
    q = sp.sympify(d["q"]).subs(sp.Symbol("t"), t)
    A4 = O.from_Q([sp.sympify(x).subs(sp.Symbol("u"), u) for x in d["A4_Q"]])
    B = O.from_Q(BQ)
    rec = {}
    for tv in ["0.5", "1.3", "2.5", "3.5"]:
        mp.mp.dps = dps
        bp, big, _ = ME.evaluate(B, tv, dps)
        th4 = Ph4_thetas(tv, 4)
        qv = mp.mpf(sp.Rational(q.subs(t, sp.Rational(tv))).p) / sp.Rational(q.subs(t, sp.Rational(tv))).q
        res = bp + qv * th4[0]
        wv = sp.Rational(tv) ** 2
        a4 = [mp.mpf(sp.Rational(c.as_expr().subs(O.w_sym, wv)).p) / sp.Rational(c.as_expr().subs(O.w_sym, wv)).q * th4[k] / 2**k for k, c in enumerate(A4)]
        rec[tv] = {"BP": mp.nstr(bp, 20), "qPh4": mp.nstr(qv * th4[0], 20), "rel": mp.nstr(abs(res) / max(big, abs(qv*th4[0])), 4),
                   "A4_Ph4_rel": mp.nstr(abs(sum(a4)) / max(abs(x) for x in a4), 4)}
        print(tv, rec[tv], flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
