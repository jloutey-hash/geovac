"""FD1-FD3 step F: algebraic (hyperexponential) solutions = order-1 right factors, mod p, of an
operator in Q-form sum_i t^(2i) Q_i(theta_t), following debug/p3e2w_slice_E_hyperexp.py but with the
third singular point t = 3 of the star.  An algebraic solution is log-free at t = 0, so lies in the
span of the log-free Frobenius solutions (free slots 0).  Ansatz:
  y (1 - t^2)^(-al) (1 - t^2/4)^(-be) (1 - t^2/9)^(-ga) = polynomial of degree <= mu_max - 2(al+be+ga).
Control: the chain's M6 = L_C/(theta-2) must return t^3/(1-t^2)^2 (al = -2)."""
from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import density_loop_C_factor as C  # noqa: E402

u = sp.Symbol("u")
P = C.P


def fr(x):
    x = Fraction(x)
    return x.numerator % P * pow(x.denominator % P, P - 2, P) % P


def frob_free(Q, lam, NT=400):
    """log-free solution t^lam sum c_n t^(2n), mod p, free slots set to 0 (D6.frob_free, copied)."""
    Qp = [sp.Poly(q, u) for q in Q]

    def ev(i, x):
        val = Fraction(0)
        for (m,), c in zip(Qp[i].monoms(), Qp[i].coeffs()):
            val += Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) * Fraction(x) ** m
        return fr(val)
    c = [1]
    free = []
    for n in range(1, NT // 2):
        s_ = sum(ev(i, lam + 2 * (n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0) % P
        dd = ev(0, lam + 2 * n)
        if dd == 0:
            if s_ != 0:
                return None, None
            c.append(0)
            free.append(lam + 2 * n)
        else:
            c.append((-s_) * pow(dd, P - 2, P) % P)
    return {lam + 2 * n: cn for n, cn in enumerate(c) if cn}, free


def binom_series(alpha, scale, N):
    out, c = {}, Fraction(1)
    for n in range(N // 2 + 1):
        if c:
            out[2 * n] = fr(c)
        c = c * (Fraction(alpha) - n) / (n + 1) * (-Fraction(scale))
    return out


def mul(a, b, top):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            if e1 + e2 < top:
                out[e1 + e2] = (out.get(e1 + e2, 0) + c1 * c2) % P
    return out


def algebraic_solutions(Q, e1, e2, e3, mu_max, N=400, verbose=True):
    roots = sorted(r for r in set(sp.roots(sp.Poly(Q[0], u).as_expr(), u)) if r.is_rational)
    sols = []
    for lam in roots:
        ser, free = frob_free(Q, int(lam) if lam.is_integer else lam, NT=N)
        if ser is not None:
            sols.append((lam, ser))
    lo = min(l for l, _ in sols)
    found = []
    for al, be, ga in itertools.product(e1, e2, e3):
        w = mul(mul(binom_series(-Fraction(al), 1, N), binom_series(-Fraction(be), Fraction(1, 4), N), N + 4),
                binom_series(-Fraction(ga), Fraction(1, 9), N), N + 4)
        dmax = mu_max - 2 * (al + be + ga)
        if dmax < lo:
            continue
        top = N + int(lo) - 4
        cols = [mul(s_, w, top) for _, s_ in sols]
        rows = [[col.get(e, 0) for col in cols] for e in range(int(sp.floor(dmax)) + 1, top)]
        if not rows:
            continue
        nul = len(cols) - C.rank_mod(rows)
        if nul:
            found.append((str(al), str(be), str(ga), nul, str(dmax)))
            if verbose:
                print(f"   al={al} be={be} ga={ga}: {nul} algebraic solution(s), degree <= {dmax}", flush=True)
    return [str(l) for l, _ in sols], found


if __name__ == "__main__":
    sys.path.insert(0, str(HERE.parent / "tests"))
    import test_paper59_density_ring as DR
    import four_density_star_ops as O
    QC = [sp.sympify(str(q)).subs(sp.Symbol("u"), u) for q in DR.QC]
    M6, rem = O.rdiv(O.from_Q(QC), O.from_Q([u - 2]))
    QM6, _ = O.to_Q(M6)
    print("control: chain M6", flush=True)
    print(algebraic_solutions(QM6, [0, -1, -2, -3], [0, -1, -2], [0], 3))
