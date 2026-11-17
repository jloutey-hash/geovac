"""Paper 63 owed item (B): tail derivations of the branch constants at t = -(n-1), exact (sympy).

LOCAL RULE.  For a source s(r) P_L(cos) with s = e^{-k r} r^a sum_j s_j r^-j, the e^{-k r} part of its channel-L
Coulomb potential, (4pi/(2L+1)) [r^L int_r^oo u^{1-L} s - r^{-L-1} int_r^oo u^{L+2} s], is the unique formal
solution e^{-k r} r^a sum_j phi_j r^-j of the radial Poisson equation
      (theta - L)(theta + L + 1) Phi = -4 pi r^2 s,
      (theta^2 + theta - L(L+1)) e^{-k r} r^b = [k^2 r^2 - 2k(b+1) r + b(b+1) - L(L+1)] e^{-k r} r^b.
Leading: phi_0 = -4 pi s_0 / k^2 (L-independent); L enters at relative order r^-2.

CONTROLS (analytic; must reproduce committed numbers):
  T1  P63L1 leaf tail E ~ sqrt(2/pi) r^-1/2 e^-r [1 + 11/8 r^-1 + 97/128 r^-2 - 463/1024 r^-3] (record C3).
  T2  P63L1 A = 1/(3 pi^2), kappa = 47/8;   T3  record SA / P63L1 toy (s chain of three): A = 1/pi^2, kappa = 31/8.
  T4  star of five (prop:star_attain): eps^(1/2) int I0 R^4 x^-2 -> c_5 Gamma(1/2) = pi^2/(8 sqrt 2).
DERIVED (new shapes; prediction constants only, from tails; no evaluation of any branch integral):
  B1  chain p_z leaf - p_z middle - d root (h_2(d) t^3 K0(t r)/(2 pi^2), radial t^5 K2(t r) P2), t on the root:
      S1(eps) = int_1^oo I2((2-eps)x) x^2 Phi2[(2/3) f_p E](x) dx (e^-2x part);  eps^(3/2) S1 -> A1 (1 + kappa1 eps).
  B2  same chain, root (lap - (5/2) h_2(d)) t^3 K0/(2 pi^2): leading angular weight cancels by construction;
      S2(eps) = int_1^oo x^2 [ (I0(y) + I1(y)/y) Phi0[f_p E] - I2(y) Phi2[f_p E] ] dx,  y = (2-eps)x;
      eps^(1/2) S2 -> A2 (1 + kappa2 eps) + odd powers of sqrt(eps).
Output: debug/data/p63odd_tail.json
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
N = 8  # series terms kept
pi = sp.pi
eps = sp.symbols("epsilon", positive=True)


def a_coef(nu, N=N):
    out = [sp.Integer(1)]
    for k in range(1, N):
        out.append(out[-1] * sp.Rational(4 * nu**2 - (2 * k - 1) ** 2, 8 * k))
    return out


class Ser:
    """pref * e^{-k x} x^a * sum_j c_j x^-j."""

    def __init__(self, pref, k, a, c):
        self.pref, self.k, self.a, self.c = sp.nsimplify(pref) if not isinstance(pref, sp.Expr) else pref, k, sp.nsimplify(a), list(c)[:N]

    def __mul__(self, o):
        c = [sum(self.c[i] * o.c[j - i] for i in range(j + 1)) for j in range(N)]
        return Ser(self.pref * o.pref, self.k + o.k, self.a + o.a, c)

    def scale(self, s):
        return Ser(self.pref * s, self.k, self.a, self.c)

    def xpow(self, p):
        return Ser(self.pref, self.k, self.a + p, self.c)

    def normalized(self):
        c0 = self.c[0]
        return Ser(self.pref * c0, self.k, self.a, [x / c0 for x in self.c])


def K_ser(nu):
    return Ser(sp.sqrt(pi / 2), 1, sp.Rational(-1, 2), a_coef(nu))


def local_potential(s: Ser, L: int) -> Ser:
    """e^{-k r} part of the channel-L potential of the source s (formal series)."""
    k, a = s.k, s.a
    rhs = [-4 * pi * s.pref * c for c in s.c]  # coefficient of e^{-kr} r^{a+2-j}
    phi = []
    for j in range(N):
        # r^{a+2-j}: k^2 phi_j - 2k(b+1) phi_{j-1} + (b(b+1) - L(L+1)) phi_{j-2}, b = a - (j-1), a - (j-2)
        acc = rhs[j]
        if j >= 1:
            b = a - (j - 1)
            acc += 2 * k * (b + 1) * phi[j - 1]
        if j >= 2:
            b = a - (j - 2)
            acc -= (b * (b + 1) - L * (L + 1)) * phi[j - 2]
        phi.append(sp.simplify(acc / k**2))
    return Ser(sp.Integer(1), k, a, phi)


def kernel_I(nu, extra_I1_over_y=False):
    """I_nu(y) ~ (2 pi)^-1/2 y^-1/2 e^y sum (-1)^k a_k y^-k (+ I1/y if requested): returns list c_j in y^-j."""
    c = [(-1) ** k * x for k, x in enumerate(a_coef(nu))]
    if extra_I1_over_y:
        c1 = [(-1) ** k * x for k, x in enumerate(a_coef(1))]
        c = [c[j] + (c1[j - 1] if j >= 1 else 0) for j in range(N)]
    return c


def branch_expansion(terms, order=2):
    """terms: list of (sign, kernel coeffs c_j in y^-j with y = (2-eps)x, H Ser with k = 2).
    Returns the singular part of S(eps) = sum sign int_1^oo k(y) H dx as a series in eps (half-integer powers),
    using int_1^oo x^q e^{-eps x} = Gamma(q+1) eps^{-q-1} + analytic (q half-integer)."""
    tot = 0
    for sgn, c, H in terms:
        assert H.k == 2
        for j in range(order + 2):
            for kk in range(order + 2):
                q = H.a - sp.Rational(1, 2) - j - kk
                tot += (sgn * H.pref * H.c[kk] * c[j] * (2 * pi) ** sp.Rational(-1, 2)
                        * sp.series((2 - eps) ** (sp.Rational(-1, 2) - j), eps, 0, order + 2).removeO() * sp.gamma(q + 1) * eps ** (-q - 1))
    return sp.expand(tot)


def leading(expr, nterms=2):
    s = sp.series(expr, eps, 0, None)
    out = []
    for term in sp.Add.make_args(sp.expand(expr)):
        out.append(term)
    # collect by power of eps
    d = sp.collect(sp.expand(expr), eps, evaluate=False)
    items = sorted(d.items(), key=lambda kv: float(sp.Poly(kv[0], eps).degree() if kv[0].is_polynomial(eps) else sp.limit(sp.log(kv[0]) / sp.log(eps), eps, 0)))
    return items[:nterms]


def powers(expr):
    d = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        c, p = term.as_coeff_exponent(eps)
        d[p] = sp.simplify(d.get(p, 0) + c)
    return dict(sorted(d.items(), key=lambda kv: kv[0]))


def main() -> None:
    out = {}
    # densities (fingerprint weights, scale 1)
    f_p = K_ser(1).scale(-1 / (2 * pi**2))           # p_z radial factor -K1/(2 pi^2)
    rho_s = K_ser(0).scale(1 / (2 * pi**2))          # s density K0/(2 pi^2)
    # T1: E = e^{-r} part of the L = 1 potential of f_p
    E = local_potential(f_p, 1).normalized()
    out["T1_E_pref"] = str(sp.simplify(E.pref))
    out["T1_E_coeffs"] = [str(x) for x in E.c[:4]]
    out["T1_pass"] = (sp.simplify(E.pref - sp.sqrt(2 / pi)) == 0
                      and E.c[:4] == [1, sp.Rational(11, 8), sp.Rational(97, 128), sp.Rational(-463, 1024)])
    # s leaf potential e^{-r} part = -(2/pi) R/r: check against R's series d_m (p63l1_engine)
    Vs = local_potential(rho_s, 0).normalized()
    import sys
    sys.path.insert(0, str(HERE))
    from p63l1_engine import _R_series_coeffs  # exact Fractions
    d = _R_series_coeffs(N)
    Rser = Ser(sp.sqrt(pi / 2), 1, sp.Rational(1, 2), [sp.Rational(x.numerator, x.denominator) for x in d])
    Rr = Rser.xpow(-1).scale(-2 / pi).normalized()
    out["T1b_s_leaf_tail_pass"] = bool(sp.simplify(Vs.pref - Rr.pref) == 0 and all(sp.simplify(x - y) == 0 for x, y in zip(Vs.c, Rr.c)))

    def report(name, S_sing, lead_pow):
        P = powers(S_sing)
        lead = sp.nsimplify(sp.simplify(P.get(lead_pow, 0)))
        nxt = sp.simplify(P.get(lead_pow + 1, 0) / lead) if lead != 0 else None
        out[name] = {"powers": {str(k): str(v) for k, v in list(P.items())[:4]},
                     "A": str(lead), "A_float": float(lead), "kappa": str(nxt), "kappa_float": float(nxt) if nxt is not None else None}
        print(name, out[name]["A"], out[name]["kappa"])

    # T2: P63L1, H = (2/(3pi)) x^2 E^2, kernel I0
    H = (E * E).xpow(2).scale(2 / (3 * pi))
    report("T2_P63L1", branch_expansion([(1, kernel_I(0), H)]), sp.Rational(-3, 2))
    out["T2_pass"] = sp.simplify(sp.sympify(out["T2_P63L1"]["A"]) - 1 / (3 * pi**2)) == 0 and out["T2_P63L1"]["kappa"] == "47/8"
    # T3: s chain, H_sing = (2/pi)^3 R^2, kernel I0
    H = (Rser * Rser).scale((2 / pi) ** 3)
    report("T3_s_chain", branch_expansion([(1, kernel_I(0), H)]), sp.Rational(-3, 2))
    out["T3_pass"] = sp.simplify(sp.sympify(out["T3_s_chain"]["A"]) - 1 / pi**2) == 0 and out["T3_s_chain"]["kappa"] == "31/8"
    # T4: star of five, leading constant c5 Gamma(1/2) of int I0((4-eps)x) R^4 x^-2 (k = 4: leading only, by hand)
    R4 = Rser * Rser * Rser * Rser
    c5G = R4.pref * (2 * pi * 4) ** sp.Rational(-1, 2) * sp.gamma(sp.Rational(1, 2))
    out["T4_star5_const"] = str(sp.simplify(c5G))
    out["T4_pass"] = sp.simplify(c5G - pi**2 / (8 * sp.sqrt(2))) == 0

    # B1: d root; w = (2/3) f_p E (e^{-2x} part of the middle source in channel 2), Phi2, kernel I2
    E_full = local_potential(f_p, 1)            # un-normalised e^{-r} part of g
    w = (f_p * E_full).scale(sp.Rational(2, 3))
    Phi2 = local_potential(w, 2)
    H1 = Phi2.xpow(2)
    report("B1_d_root", branch_expansion([(1, kernel_I(2), H1)]), sp.Rational(-3, 2))
    # B2: mixed root; channel weights +1 (L = 0, kernel I0 + I1/y) and -1 (L = 2, kernel I2)
    w0 = f_p * E_full
    H0 = local_potential(w0, 0).xpow(2)
    H2 = local_potential(w0, 2).xpow(2)
    S2 = branch_expansion([(1, kernel_I(0, True), H0), (-1, kernel_I(2), H2)], order=3)
    P = powers(S2)
    out["B2_powers"] = {str(k): str(sp.simplify(v)) for k, v in list(P.items())[:4]}
    report("B2_mixed_root", S2, sp.Rational(-1, 2))
    out["B2_leading_eps^-3/2_coeff"] = str(sp.simplify(P.get(sp.Rational(-3, 2), 0)))
    out["angular_check_B2"] = str(sp.integrate(sp.Symbol("c") ** 2 * (1 - sp.Rational(5, 2) * (3 * sp.Symbol("c") ** 2 - 1) / 2), (sp.Symbol("c"), -1, 1)))
    json.dump(out, open(HERE / "data" / "p63odd_tail.json", "w"), indent=1, default=str)
    for k in ("T1_pass", "T1b_s_leaf_tail_pass", "T2_pass", "T3_pass", "T4_pass"):
        print(k, out[k])


if __name__ == "__main__":
    main()
