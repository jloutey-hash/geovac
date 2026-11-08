"""Loop test PT2, step C: structure of the triangle operator L_G (step B, corrected) and of the
chain control's operator.

L_G = sum_i t^{2i} Q_i(theta_t), annihilates G up to the elementary remainder c t^2 (measured here).
(1) remainder: L_G G / t^2 at three t, 18 digits;
(2) mod-p minimal-annihilator search (orders 1..r-1, degree <= 8, 20, 40; control: order r finds
    the operator) for the Frobenius solutions at t = 0, in w = t^2;
(3) quadric count of the solution space at a regular rational point (exact Fractions).
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
u = sp.Symbol("u")
P = (1 << 61) - 1
OUT = Path(__file__).parent / "data" / "density_loop_C.json"


def frob_series(Q, lam, NT):
    """Solution w^lam sum c_n w^n of sum_i w^i Q_i(2 theta_w) (theta_t = 2 theta_w), mod p.
    Recurrence: sum_i Q_i(2(lam + n - i)) c_{n-i} = 0."""
    Qp = [sp.Poly(q, u) for q in Q]

    def ev(i, x):  # Q_i(x) mod p, x a Fraction
        val = Fraction(0)
        for (m,), c in zip(Qp[i].monoms(), Qp[i].coeffs()):
            val += Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) * x**m
        return val.numerator % P * pow(val.denominator % P, P - 2, P) % P
    c = [1]
    for n in range(1, NT):
        s_ = sum(ev(i, 2 * (lam + n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0) % P
        d = ev(0, 2 * (lam + n))
        assert d != 0, ("resonance", n)
        c.append((-s_) * pow(d, P - 2, P) % P)
    return c


def rank_mod(rows):
    rows = [x[:] for x in rows]
    rk = 0
    for col in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col]), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = pow(rows[rk][col], P - 2, P)
        rows[rk] = [x * iv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col]:
                f = rows[i][col]
                rows[i] = [(x - f * yv) % P for x, yv in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def annihilators(c, lam, order, D, NT):
    """Count operators sum_{j<=order, d<=D} a_{jd} w^d theta_w^j killing w^lam sum c_n w^n.
    theta_w^j w^{lam+n} = (lam+n)^j w^{lam+n}: coefficient of w^{lam+m}: sum a_{jd} (lam+m-d)^j c_{m-d}."""
    lamp = lam.numerator % P * pow(lam.denominator % P, P - 2, P) % P
    rows = []
    for m in range(NT):
        row = []
        for j in range(order + 1):
            for d in range(D + 1):
                row.append(pow((lamp + m - d) % P, j, P) * c[m - d] % P if m - d >= 0 else 0)
        rows.append(row)
    return (order + 1) * (D + 1) - rank_mod(rows)


def taylor_basis(Q, t0, N):
    """Exact Taylor solutions at a regular point t0 of sum_i t^{2i} Q_i(theta)."""
    t = sp.Symbol("t")
    f = sp.Function("f")(t)
    e = 0
    for i, q in enumerate(Q):
        cs = sp.Poly(q, u).all_coeffs()[::-1]
        cur, acc = f, 0
        for cc in cs:
            acc += cc * cur
            cur = t * sp.diff(cur, t)
        e += t ** (2 * i) * acc
    e = sp.expand(e)
    r = max(k for k in range(10) if e.coeff(sp.diff(f, t, k)) != 0) if True else 0
    a = [sp.expand(e.coeff(sp.diff(f, t, k))) if k else sp.expand(e.subs({sp.diff(f, t, j): 0 for j in range(1, r + 1)}).coeff(f))
         for k in range(r + 1)]
    x = sp.Symbol("x")
    A = [sp.Poly(sp.expand(ak.subs(t, t0 + x)), x) for ak in a]
    Ac = [[Fraction(int(sp.Rational(cc).p), int(sp.Rational(cc).q)) for cc in Ak.all_coeffs()[::-1]] for Ak in A]
    sols = []
    for k0 in range(r):
        y = [Fraction(0)] * N
        y[k0] = Fraction(1)
        # sum_k a_k(x) y^{(k)} = 0, y = sum y_n x^n; coefficient of x^m
        for m in range(N - r):
            # find y_{m+r}
            tot = Fraction(0)
            lead = None
            for k in range(r + 1):
                for d, cf in enumerate(Ac[k]):
                    n = m - d + k  # y^{(k)} coefficient x^{m-d} uses y_{m-d+k}
                    if m - d < 0 or n >= N:
                        continue
                    fac = Fraction(1)
                    for q_ in range(k):
                        fac *= (m - d + k - q_)
                    if n == m + r and k == r and d == 0:
                        lead = cf * fac
                        continue
                    tot += cf * fac * y[n]
            y[m + r] = -tot / lead
        sols.append(y)
    return sols, r


def quadric_count(sols, N):
    r = len(sols)
    prods = []
    for i in range(r):
        for j in range(i, r):
            pr = [sum(sols[i][k] * sols[j][n - k] for k in range(n + 1)) for n in range(N)]
            prods.append(pr)
    M = sp.Matrix([[int(v.numerator) % P * pow(int(v.denominator) % P, P - 2, P) % P for v in row] for row in prods])
    rows = [list(map(int, M.row(i))) for i in range(M.rows)]
    return len(prods) - rank_mod(rows)


def main() -> None:
    QG = [sp.sympify(q) for q in json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())["QG"]]
    rec = {"QG": [str(q) for q in QG]}
    # (1) remainder
    import density_loop_B2_inhom as B2
    mp.mp.dps = 22
    th = B2.theta_j_psi_factory(6)
    rem = {}
    for t in ["0.6", "1.3", "2.7"]:
        r_ = B2.LG(QG, t, th) / mp.mpf(t) ** 2
        rem[t] = mp.nstr(r_, 18)
        print("L_G G / t^2 at t =", t, ":", mp.nstr(r_, 18), flush=True)
    rec["remainder_over_t2"] = rem
    # (2) annihilator search on Frobenius solutions at t = 0 (exponents in w: 0 and 3/2)
    NT = 340
    res = {}
    for lam in [Fraction(0), Fraction(3, 2)]:
        try:
            c = frob_series(QG, lam, NT)
        except AssertionError as e:
            print("lam", lam, "resonant:", e, flush=True)
            continue
        for order in (1, 2, 3):
            for D in ((8, 20, 40) if order < 3 else (2, 4)):
                a = annihilators(c, lam, order, D, NT)
                res[f"lam={lam},order={order},deg={D}"] = a
                print(f"lam {lam}: order {order}, degree {D}: annihilators {a}", flush=True)
    rec["annihilators"] = res
    # (3) quadric count at t0 = 1/2
    sols, r = taylor_basis(QG, sp.Rational(1, 2), 70)
    qc = quadric_count(sols, 70)
    print("order", r, " quadric count at t0 = 1/2:", qc, flush=True)
    rec["quadric_count"] = qc
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
