"""PM independent recheck of GY1-GY3 (docs/predictions/2026-10-23_gaussian_geminal_triangle.md).

Own code path; nothing imported from the agent's gy_* scripts.
(1) Exact: F1 right-divides L9G with zero remainder, and the quotient is F2 up to a rational left factor.
(2) Exact: Newton polygon of L9G at x = infinity and its characteristic polynomial.
(3) Exact: indicial exponents of L9G at x = 1, 4, the roots of the degree-7 factor (numerically located) and
    whether each is a regular singular / apparent point.
(4) Numerical, independent of the agent's routes: Y_G(t) by a position-space reduction.
    The geminal pair: g(x13) g(x23) = exp(-2 gamma |x3 - m|^2 - gamma |x12|^2 / 2), m = (x1 + x2)/2,
    so the apex integral is h(|m|) = int rho_3(x) exp(-2 gamma |x - m|^2) d^3x (one radial quadrature);
    then with w = |x12| as the angular variable, d mu / w = dw / (r1 r2), and
    Y = 8 pi^2 int int r1 r2 rho(r1) rho(r2) int_{|r1-r2|}^{r1+r2} exp(-gamma w^2/2) h(m(w)) dw dr1 dr2,
    m^2 = (2 r1^2 + 2 r2^2 - w^2)/4, rho_a(r) = a^3 K0(a r)/(2 pi^2) (transform (k^2/a^2 + 1)^(-3/2)).
    Compared with the agent's grid values of Y_G at three t."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.special import k0

ROOT = Path(__file__).parent
D = json.loads((ROOT / "data" / "gy_pieces.json").read_text())
x = sp.Symbol("x")


def op(lst):
    return [sp.Poly(sp.sympify(s.replace("^", "**")), x) for s in lst]


def th(p, m=1):
    """(2x d/dx)^m p"""
    for _ in range(m):
        p = sp.Poly(2 * x * p.diff(x).as_expr(), x)
    return p


def compose(A, B):
    """A o B with theta = t d/dt = 2x d/dx: theta^i b = sum_m C(i,m) (theta^m b) theta^(i-m)."""
    out = [sp.Poly(0, x) for _ in range(len(A) + len(B) - 1)]
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            for m in range(i + 1):
                out[i - m + j] += a * sp.binomial(i, m) * th(b, m)
    return out


def rdiv(L, R):
    """pseudo right division over Q(x): L = Q o R + rem; returns (Q as rational functions, rem)."""
    L = [c.as_expr() for c in L]
    Q = [0] * (len(L) - len(R) + 1)
    Rx = [c.as_expr() for c in R]
    rem = list(L)
    while len(rem) >= len(R) and any(sp.simplify(c) != 0 for c in rem):
        k = len(rem) - len(R)
        q = sp.cancel(rem[-1] / Rx[-1])
        Q[k] = q
        mon = [sp.Poly(0, x)] * k + [sp.Poly(1, x)]
        # theta^k o R, then multiply by q on the left
        tk = compose([sp.Poly(0, x)] * k + [sp.Poly(1, x)], R)
        for i, c in enumerate(tk):
            rem[i] = sp.cancel(rem[i] - q * c.as_expr())
        while rem and sp.simplify(rem[-1]) == 0:
            rem.pop()
        if len(rem) < len(R):
            break
    return Q, rem


def main():
    out = {}
    L = op(D["operator_L9G"]["theta_form"])
    F1 = op(D["factor_1"])
    F2 = op(D["factor_2"])
    print("orders", len(L) - 1, len(F1) - 1, len(F2) - 1, "x-degrees", max(c.degree() for c in L))
    Q, rem = rdiv(L, F1)
    out["F1_right_divides_L9G"] = (len(rem) == 0)
    ratios = {sp.factor(sp.cancel(q / f.as_expr())) for q, f in zip(Q, F2) if f.as_expr() != 0}
    out["quotient_over_F2_ratios"] = [str(r) for r in ratios]
    out["quotient_is_F2_up_to_left_factor"] = (len(ratios) == 1)
    print("F1 | L9G exactly:", out["F1_right_divides_L9G"], " quotient/F2:", out["quotient_over_F2_ratios"])

    # Newton polygon at infinity
    pts = [(k, c.degree()) for k, c in enumerate(L)]
    out["newton_points"] = pts
    on_edge = [k for k, dg in pts if dg + k == 18]
    C = sp.Symbol("C")
    cp = sp.factor(sum(sp.LC(L[k]) * C**k for k in on_edge))
    out["edge_k_plus_deg_18"] = on_edge
    out["char_poly"] = str(cp)
    print("points", pts, " edge", on_edge, " char", cp)

    # local exponents: convert to D = d/dx form
    Dsym = sp.Symbol("D")
    # theta^k = sum_j S(k,j) 2^k ... build via operator algebra with D: theta = 2 x D
    Pd = [0] * len(L)
    thk = {0: {0: sp.Integer(1)}}       # theta^k as {j: coefficient poly in x of D^j}
    for k in range(1, len(L)):
        prev = thk[k - 1]
        new = {}
        for j, c in prev.items():       # 2x D (c D^j) = 2x c' D^j + 2x c D^(j+1)
            new[j] = new.get(j, 0) + 2 * x * sp.diff(c, x)
            new[j + 1] = new.get(j + 1, 0) + 2 * x * c
        thk[k] = {j: sp.expand(v) for j, v in new.items()}
    for k, c in enumerate(L):
        for j, v in thk[k].items():
            Pd[j] = sp.expand(Pd[j] + c.as_expr() * v)
    lead = sp.factor(Pd[-1])
    out["lead_D_form"] = str(lead)
    r = sp.Symbol("r")

    def indicial(x0):
        s = sp.Symbol("s")
        P = [sp.Poly(sp.expand(p.subs(x, x0 + s)), s) for p in Pd]
        v = [min((m[0] for m in p.monoms()), default=10**9) for p in P]
        delta = min(v[k] - k for k in range(len(P)))
        poly = 0
        for k in range(len(P)):
            if v[k] - k == delta:
                ff = sp.Integer(1)
                for i in range(k):
                    ff *= (r - i)
                poly += P[k].coeff_monomial(s ** v[k]) * ff
        poly = sp.Poly(sp.expand(poly), r)
        return delta, v[-1] - (len(P) - 1), poly

    for x0 in (sp.Integer(1), sp.Integer(4), sp.Integer(0)):
        delta, dlead, poly = indicial(x0)
        regular = (delta == dlead)
        roots = sp.roots(poly)
        out[f"x0={x0}"] = {"fuchsian_here": regular, "indicial_degree": poly.degree(),
                           "roots": {str(k): m for k, m in roots.items()}}
        print("x0 =", x0, "regular-singular" if regular else "IRREGULAR", "deg", poly.degree(), "roots", roots)
    # the degree-7 factor: all its roots real? apparent test numerically = exponents are distinct nonneg integers
    f7 = [f for f, m in sp.factor_list(Pd[-1])[1] if sp.degree(f, x) == 7][0]
    rr = [complex(z) for z in sp.Poly(f7, x).nroots(n=30)]
    out["deg7_roots"] = [str(z) for z in rr]
    print("degree-7 factor roots:", [f"{z.real:.4f}{z.imag:+.4f}i" for z in rr])

    # (4) independent Y_G values
    gam = 1.0

    def h(R, t, n=400):
        # radial quadrature on [0, Rmax] with r = s^2 substitution near 0
        s, w = np.polynomial.legendre.leggauss(n)
        smax = np.sqrt(float(np.max(R)) + 12.0)
        sv = (s + 1) / 2 * smax
        ws = w / 2 * smax
        rv = sv**2
        jac = 2 * sv
        rho3 = t**3 * k0(t * rv) / (2 * np.pi**2)
        RR = np.asarray(R)[..., None]
        a = 4 * gam * rv * RR
        # exp(-2g(r^2+R^2)) sinh(4grR)/(4grR) * 2, stably
        e1 = np.exp(-2 * gam * (rv - RR) ** 2)
        e2 = np.exp(-2 * gam * (rv + RR) ** 2)
        small = a < 1e-8
        kern = np.where(small, 2 * np.exp(-2 * gam * (rv**2 + RR**2)), (e1 - e2) / np.where(small, 1, a))
        return 2 * np.pi * np.sum(ws * jac * rho3 * rv**2 * kern, axis=-1)

    def Y(t, n=160, nw=48):
        """outer r1 on r = s^2; for each r1 the r2 integral is split at r2 = r1 (kink of |r1 - r2|):
        [0, r1] with r2 = r1 u^2, [r1, r1 + 30] with r2 = r1 + v^2."""
        g, gw_ = np.polynomial.legendre.leggauss(n)
        smax = np.sqrt(30.0)
        sv = (g + 1) / 2 * smax
        rv = sv**2
        wr = gw_ / 2 * smax * 2 * sv
        u01 = (g + 1) / 2
        w01 = gw_ / 2
        gw, gwt = np.polynomial.legendre.leggauss(nw)
        tot = 0.0
        for i in range(n):
            r1 = rv[i]
            pieces = [(r1 * u01**2, w01 * 2 * u01 * r1),
                      (r1 + (u01 * np.sqrt(30.0))**2, w01 * np.sqrt(30.0) * 2 * u01 * np.sqrt(30.0))]
            acc = 0.0
            for r2, w2 in pieces:
                lo = np.abs(r1 - r2)
                hi = r1 + r2
                wv = lo[:, None] + (gw[None, :] + 1) / 2 * (hi - lo)[:, None]
                ww = gwt[None, :] / 2 * (hi - lo)[:, None]
                m = np.sqrt(np.maximum(2 * r1**2 + 2 * r2[:, None] ** 2 - wv**2, 0)) / 2
                inner = np.sum(ww * np.exp(-gam * wv**2 / 2) * h(m, t), axis=1)
                acc += np.sum(w2 * r2 * k0(r2) / (2 * np.pi**2) * inner)
            tot += wr[i] * r1 * k0(r1) / (2 * np.pi**2) * acc
        return 8 * np.pi**2 * tot

    agent = D["on_Y"]["grid_h_1_128"]
    out["Y_check"] = {}
    for ts in ("0.8", "1.3", "2.5"):
        t = float(ts)
        y1 = Y(t)
        y2 = Y(t, n=200, nw=64)
        ya = float(agent[ts]["Y"])
        out["Y_check"][ts] = {"pm": y2, "pm_coarse_vs_fine": abs(y1 - y2) / abs(y2), "agent": ya,
                              "rel_diff": abs(y2 - ya) / abs(ya)}
        print(f"t = {ts}: PM {y2:.15e} (coarse/fine {abs(y1-y2)/abs(y2):.1e})  agent {ya:.15e}  rel {abs(y2-ya)/abs(ya):.2e}")
    (ROOT / "data" / "gy_pm_check.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
