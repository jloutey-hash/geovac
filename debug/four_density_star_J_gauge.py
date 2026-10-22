"""FD1-FD3 step J: mod-p rational-gauge tests (density_loop_D5_rational.test generalized to the third
singular point t = 3).  Sub test (R): exists A = sum_{j<r} a_j(t) theta^j, a_j = p_j(t)/den,
den = (t^2-1)^b (t^2-4)^c (t^2-9)^e, p_j Laurent in t with exponents in [-a, Dg], such that
M (A f) = 0 for a solution f of a reference operator of order r: then A maps Sol(ref) into Sol(M)
(ref irreducible).  Quotient test (L): the same with the adjoint of M and the adjoint reference.
References: f3 (F_3's operator, 3F2 at t^2/4), f4 (F_4's operator Q4, exponent-2 Frobenius solution),
Legendre at t^2, t^2/4, t^2/9."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import density_loop_C_factor as C  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402
import four_density_star_F_hyperexp as F  # noqa: E402

u = sp.Symbol("u")
P = C.P
t = sp.Symbol("t")
Q4 = [9 * u**2 * (u - 2) ** 2, -2 * u**2 * (5 * u**2 + 1), (u + 1) ** 4]
L3 = [4 * u**3, -(u + 1) ** 3]
QSUN = [9 * u**2, -2 * (5 * u**2 + 10 * u + 6), (u + 2) ** 2]


def adjoint(Q):
    return [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]


def apply_L(Q, ser):
    Qp = [sp.Poly(q, u) for q in Q]
    out = {}
    for i, qp in enumerate(Qp):
        cs = [(m[0], F.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)))) for m, c in zip(qp.monoms(), qp.coeffs())]
        for e, c in ser.items():
            val = sum(cc * pow(e % P, m, P) for m, cc in cs) % P
            if val:
                out[e + 2 * i] = (out.get(e + 2 * i, 0) + val * c) % P
    return {e: c for e, c in out.items() if c}


def inv_series(den, N):
    cs = [F.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))) for c in sp.Poly(den, t).all_coeffs()[::-1]]
    out = [pow(cs[0], P - 2, P)]
    for n in range(1, N):
        s_ = sum(cs[k] * out[n - k] for k in range(1, min(n, len(cs) - 1) + 1)) % P
        out.append((-s_) * out[0] % P)
    return {n: c for n, c in enumerate(out) if c}


def refs(NT):
    f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4))
    f4, _ = F.frob_free(Q4, 2, NT=NT)
    out = {"F3": (f3, 3, L3), "F4": (f4, 4, Q4)}
    # odd-weight four-Bessel moment int x K0^3 K0(t x) dx: Mellin operator of order 2 (sunrise type)
    fs, _ = F.frob_free(QSUN, 0, NT=NT)
    out["Sun"] = (fs, 2, QSUN)
    for nm, sc in (("Leg(t^2)", 1), ("Leg(t^2/4)", Fraction(1, 4)), ("Leg(t^2/9)", Fraction(1, 9))):
        out[nm] = (D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT, scale=Fraction(sc)), 2,
                   [u**2, -sp.Rational(sc) * (u + 1) ** 2])
    return out


def adjoint_ref(name, NT):
    """a solution series of the adjoint of the reference operator (log-free Frobenius, lowest
    exponent available)."""
    _, r, Qr = refs(NT)[name]
    Qa = adjoint(Qr)
    for lam in sorted(set(sp.roots(sp.Poly(Qa[0], u).as_expr(), u))):
        ser, _ = F.frob_free(Qa, int(lam), NT=NT)
        if ser is not None:
            return ser
    raise RuntimeError("no log-free adjoint solution")


def test(Q, f, r, Dg, b, c, e, a=6, NT=700, par=0, extra=1):
    w = inv_series(sp.expand((t**2 - 1) ** b * (t**2 - 4) ** c * (t**2 - 9) ** e * extra), NT)
    top = min(f) + NT - 2 * len(Q) - Dg - 8
    cols = []
    for j in range(r):
        g = F.mul(w, D.theta(f, j), NT + min(f))
        for d in range(-a + ((a + par) % 2), Dg + 1, 2):
            cols.append(apply_L(Q, D.shift(g, d)))
    exps = sorted(x for x in set().union(*cols) if x < top)
    rows = [[col.get(x, 0) for col in cols] for x in exps]
    return len(cols) - C.rank_mod(rows), len(cols), len(rows)


def scan(Q, f, r, Dg, maxpow=3, a=6, NT=700, verbose=False):
    for tot in range(0, 3 * maxpow + 1):
        for b in range(maxpow + 1):
            for c in range(maxpow + 1):
                e = tot - b - c
                if e < 0 or e > maxpow:
                    continue
                for par in (0, 1):
                    n, nc, nr = test(Q, f, r, Dg, b, c, e, a=a, NT=NT, par=par)
                    if verbose:
                        print(f"      b={b} c={c} e={e} parity {par}: nullity {n} ({nc} unknowns, {nr} rows)", flush=True)
                    if n:
                        return n, (b, c, e, par)
    return 0, None


if __name__ == "__main__":
    # controls: F_3's operator contains f3 (R and L), F_4's operator contains f4 and not f3
    NT = 500
    R = refs(NT)
    print("control F3 op vs F3:", scan(L3, R["F3"][0], 3, 8, maxpow=1, NT=NT), flush=True)
    print("control F3 op vs F3 (L):", scan(adjoint(L3), adjoint_ref("F3", NT), 3, 8, maxpow=1, NT=NT), flush=True)
    print("control F4 op vs F4:", scan(Q4, R["F4"][0], 4, 8, maxpow=1, NT=NT), flush=True)
    print("control F4 op vs F3:", scan(Q4, R["F3"][0], 3, 16, maxpow=2, NT=NT), flush=True)
    print("control F3 op vs F4:", scan(L3, R["F4"][0], 4, 16, maxpow=2, NT=NT), flush=True)
