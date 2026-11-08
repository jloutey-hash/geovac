"""QG step 5c: recompute the null vector at a given (r, d) from the long series and recognise it as an
integer operator (normalised on the leading theta-coefficient of z^d or t^(2d)), then re-check the
recurrence exactly on all N coefficients and report the leading polynomial (singular points).
usage: python ring_ext_QG_5c_rational.py <gamma> <phi|ylog> <tag> <N> <dps> <r> <d> <digits> <maxden-exp>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import mpmath as mp
import sympy as sp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_guess as Gs
import ring_ext_QG_5b_guess_big as B

DATA = Path(__file__).parent / "data"
u, z = sp.symbols("u z")


def main():
    g, which, tag, N, dps, r, d, digits, mde = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), \
        int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9])
    tag = "" if tag == "-" else tag
    c, lam, step = B.series(g, which, tag, N, dps)
    p, worst = Gs.try_operator(c, lam, step, r, d, n_to=N, scale=4, tol_digits=digits)
    print("residual", worst.mid().str(3, radius=False), flush=True)
    mp.mp.dps = dps
    keys = sorted(p)
    vals = {k: mp.mpf(p[k].mid().str(dps, radius=False)) for k in keys}
    # normalise on the largest-|.| entry among the top theta-power row
    ref = max([k for k in keys if k[1] == r], key=lambda k: abs(vals[k]))
    fr = {}
    for k in keys:
        x = vals[k] / vals[ref]
        q = Fraction(mp.nstr(x, digits + 5, min_fixed=-10**9, max_fixed=10**9)).limit_denominator(10 ** mde)
        if abs(mp.mpf(q.numerator) / q.denominator - x) > mp.mpf(10) ** (-digits) * max(1, abs(x)):
            print("rationalisation FAILED at", k, flush=True)
            return
        fr[k] = q
    L = 1
    for q in fr.values():
        L = lcm(L, q.denominator)
    ints = {k: int(q * L) for k, q in fr.items()}
    G = 0
    for v in ints.values():
        G = gcd(G, v)
    P = {k: v // G for k, v in ints.items()}
    ex = Gs.exact_residuals(P, c, lam, step, N)
    print("integer operator: max |coeff| ~ 10^%d, exact recheck %s" % (len(str(max(abs(v) for v in P.values()))),
                                                                       ex.mid().str(3, radius=False)), flush=True)
    Q = [sp.expand(sum(sp.Integer(P[(i, j)]) * u**j for j in range(r + 1) if (i, j) in P)) for i in range(d + 1)]
    lead = sp.factor(sum(Q[i].coeff(u, r) * z ** i for i in range(d + 1)))
    print("leading theta^r polynomial in z (or t^2):", lead, flush=True)
    print("Q_0 =", sp.factor(Q[0]), flush=True)
    print("Q_d =", sp.factor(Q[d]), flush=True)
    out = DATA / f"ring_ext_QG_5c_{which}_g{g}_r{r}_d{d}.json"
    out.write_text(json.dumps({"r": r, "d": d, "P": {f"{i},{j}": str(v) for (i, j), v in P.items()},
                               "exact_recheck": ex.mid().str(3, radius=False), "lead": str(lead)}, indent=1))


if __name__ == "__main__":
    main()
