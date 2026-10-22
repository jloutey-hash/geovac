"""WG1-WG3 (docs/predictions/2026-10-21_wg_closure.md, committed 61d5b6c9 + fd4221cb before computing).

Source of a Mellin-route operator, exactly.  For F(t) = pref * int_0^inf x^a F_x(x) K_nu(t x) dx the
Mellin transform is Phi(s) = A(alpha - s) E(s + beta), A = M[f^m] (f = J or chi), E = M[K_nu].  The
term t^(2i) Q_i(theta) F is an integral on the line c - 2i; moving it to c crosses the poles of
Phi(w), w = s + 2i, at w = p with c < p < c + 2i.  With the recurrence sum_i Q_i(-s-2i) Phi(s+2i) = 0
(exact, see the record), L F = S(t) = -sum_i sum_p Res_{s = p - 2i} [Q_i(-s-2i) Phi(s+2i) t^(-s)].
Near w = p = alpha + n the singular part of A(alpha - w) is -sum_m f_{n,m} m!/eps^(m+1), eps = w - p,
where f^m = sum f_{n,m} x^n log^m x at x = 0.

Oscillation: N2 y = a2 theta^2 y + a1 theta y + a0 y; y'' + f y' + g y = 0 with f = (a2 + a1)/(a2 t),
g = a0/(a2 t^2); q = g - f'/2 - f^2/4 (Liouville normal form v'' + q v = 0).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(ROOT / "tests"))
import density_loop_B_ode as B  # noqa: E402

z, s, u, y = B.z, B.s, B.u, B.y
x = sp.Symbol("x", positive=True)
t = sp.Symbol("t", positive=True)
Lx, Lt, eps = sp.symbols("Lx Lt eps")
OUT = Path(__file__).parent / "data" / "wg_closure_A.json"


# ---------------------------------------------------------------- series at x = 0
def _int0x(expr: sp.Expr) -> sp.Expr:
    """int_0^x of a polynomial in x and Lx = log x (all powers of x >= 0)."""
    out = sp.Integer(0)
    P = sp.Poly(sp.expand(expr), x, Lx)
    for (n, m), c in P.terms():
        X = sp.Symbol("X", positive=True)
        out += c * sp.integrate(X**n * sp.log(X) ** m, (X, 0, x))
    return sp.expand(sp.expand_log(out, force=True).subs(sp.log(x), Lx))


def _ser(f: sp.Expr, N: int) -> sp.Expr:
    e = sp.series(f, x, 0, N).removeO()
    return sp.expand(sp.expand_log(e, force=True).subs(sp.log(x), Lx))


def J_series(N: int) -> sp.Expr:
    """J = 3 int_0^x K0 - 3 x K0 - x^2 K1, through x^(N-1)."""
    k0, k1 = _ser(sp.besselk(0, x), N + 2), _ser(sp.besselk(1, x), N + 2)
    return _trunc(3 * _int0x(k0) - 3 * x * k0 - x**2 * k1, N)


def chi_series(N: int) -> sp.Expr:
    """chi = int_0^x s K1(s) ds."""
    k1 = _ser(sp.besselk(1, x), N + 2)
    return _trunc(_int0x(sp.expand(x * k1)), N)


def _trunc(e: sp.Expr, N: int) -> sp.Expr:
    P = sp.Poly(sp.expand(e), x, Lx)
    return sp.expand(sum(c * x**n * Lx**m for (n, m), c in P.terms() if n < N))


def coeffs(e: sp.Expr, N: int) -> Dict[Tuple[int, int], sp.Expr]:
    P = sp.Poly(_trunc(e, N), x, Lx)
    return {(n, m): sp.nsimplify(c) if c.is_number and c.is_rational else c for (n, m), c in P.terms()}


# ---------------------------------------------------------------- Mellin data
def E_K(nu: int) -> Callable[[sp.Expr], sp.Expr]:
    return lambda w: 2 ** (w - 2) * sp.gamma((w - nu) / 2) * sp.gamma((w + nu) / 2)


def source(Q: List[sp.Expr], fc: Dict[Tuple[int, int], sp.Expr], alpha: int, beta: int,
           Ef: Callable[[sp.Expr], sp.Expr], c: sp.Rational, pref: sp.Expr = 1,
           drop: Optional[int] = None) -> Tuple[sp.Expr, List]:
    """S(t) = L F by the crossed residues.  drop: leave out the pole w = drop (negative control)."""
    ns = sorted({n for (n, m) in fc})
    S = sp.Integer(0)
    crossed = []
    for i, Qi in enumerate(Q):
        if i == 0:
            continue
        for n in ns:
            p = alpha + n
            if not (c < p < c + 2 * i) or p == drop:
                continue
            s0 = p - 2 * i
            mmax = max(m for (nn, m) in fc if nn == n)
            Asing = -sum(fc[(n, m)] * sp.factorial(m) / eps ** (m + 1) for m in range(mmax + 1) if (n, m) in fc)
            Qs = sp.expand(Qi.subs(u, -(s0 + eps) - 2 * i))
            Et = sp.series(Ef(p + beta + eps), eps, 0, mmax + 1).removeO()
            tt = sp.series(sp.exp(-eps * Lt), eps, 0, mmax + 1).removeO()
            expr = sp.expand(Asing * Qs * Et * tt)
            res = expr.coeff(eps, -1)
            S += -pref * res * t ** (-s0)
            crossed.append((i, p, s0, mmax))
    return sp.expand(S), crossed


def simplify_coeffs(S: sp.Expr) -> Dict[str, Tuple[str, str]]:
    """Coefficients of t^k Lt^j, simplified exactly, with a 40-digit value."""
    out = {}
    S = sp.expand(sp.expand_func(S))
    if S == 0:
        return out
    P = sp.Poly(sp.expand(S * t**6), t, Lt)
    for (k, j), cf in P.terms():
        cs = sp.simplify(sp.expand_func(cf))
        out[f"t^{k - 6} Lt^{j}"] = (str(cs), mp.nstr(sp.N(cf, 40), 12))
    return out


# ---------------------------------------------------------------- quadrature of L F
def LF_num(Q: List[sp.Expr], Fx: Callable, nu: int, a: int, tshift: int, tv, pref=1):
    """(sum_i t^(2i) Q_i(theta)) [pref t^tshift int_0^inf x^a Fx(x) K_nu(t x) dx] at t = tv."""
    rel2 = 1 + y**2 if nu == 1 else y**2
    J = max(sp.Poly(q, u).degree() for q in Q)
    al, be, _ = B.theta_powers(J, rel2, 0)
    parts = []
    for i, q in enumerate(Q):
        pc = sp.Poly(sp.expand(q.subs(u, u + tshift)), u).all_coeffs()[::-1]
        A_ = sp.expand(sum(cf * al[j] for j, cf in enumerate(pc)))
        B_ = sp.expand(sum(cf * be[j] for j, cf in enumerate(pc)))
        parts.append((i, sp.lambdify(y, A_, "mpmath"), sp.lambdify(y, B_, "mpmath")))
    tv = mp.mpf(tv)
    big = [mp.mpf(0)]

    def h(xx):
        yy = tv * xx
        k0, k1 = mp.besselk(0, yy), mp.besselk(1, yy)
        g, thg = (k1, -yy * k0 - k1) if nu == 1 else (k0, -yy * k1)
        terms = [tv ** (2 * i) * (A_(yy) * g + B_(yy) * thg) for i, A_, B_ in parts]
        w = xx**a * Fx(xx)
        big[0] = max(big[0], max(abs(w * term) for term in terms))
        return w * sum(terms)
    val = pref * tv**tshift * mp.quad(h, [0, mp.mpf(1) / 8, 1, 4, 16, 64, mp.inf])
    return val


def Jnum(xx):
    from three_electron_p3e_check import L as Lfun  # noqa: E402
    return 3 * Lfun(xx) - 3 * xx * mp.besselk(0, xx) - xx * xx * mp.besselk(1, xx)


def chinum(xx):
    xx = mp.mpf(xx)
    if xx == 0:
        return mp.mpf(0)
    return mp.pi * xx / 2 * (mp.besselk(1, xx) * mp.struvel(0, xx) + mp.besselk(0, xx) * mp.struvel(1, xx))


def rebind(e):
    """Replace foreign symbols named t, u (the test files build them with other assumptions)."""
    e = sp.sympify(e)
    m = {q: {"t": t, "u": u}[q.name] for q in e.free_symbols if q.name in ("t", "u")}
    return e.xreplace(m)


# ---------------------------------------------------------------- oscillation
def oscillation(a: List[sp.Expr], q0: sp.Rational, Tmax: int = 10) -> Dict:
    a0, a1, a2 = [sp.together(rebind(e)) for e in a]
    f = sp.cancel((a2 + a1) / (a2 * t))
    g = sp.cancel(a0 / (a2 * t**2))
    q = sp.cancel(g - sp.diff(f, t) / 2 - f**2 / 4)
    qinf = sp.limit(q, t, sp.oo)
    num, den = sp.fraction(sp.cancel(q - q0))
    roots = [r for r in sp.Poly(num, t).real_roots()] + [r for r in sp.Poly(den, t).real_roots()]
    sing = [r for r in sp.Poly(sp.fraction(sp.cancel(a2))[0], t).real_roots()]
    T = max([sp.Integer(0)] + [r for r in roots + sing if r > 0])
    Tf = float(T)
    sign_after = sp.sign(sp.N((q - q0).subs(t, T + 1), 30)) if Tf < 1e6 else None
    big = sp.series(q.subs(t, 1 / sp.Symbol("e", positive=True)), sp.Symbol("e", positive=True), 0, 4).removeO()
    qp = sp.cancel(sp.diff(q, t))
    qn, qd = sp.fraction(qp)
    qp_order = (sp.Poly(qn, t).degree() - sp.Poly(qd, t).degree()) if qp != 0 else -sp.oo
    ok = bool(qinf > q0) and Tf <= Tmax and sign_after == 1 and qp_order <= -3
    return {"q_inf": str(qinf), "T": str(sp.nsimplify(T)) if T.is_rational else f"{Tf:.6f}", "T_float": Tf,
            "q_minus_q0_sign_past_T": str(sign_after), "q_series_1_over_t": str(big),
            "q_prime_order": str(qp_order), "held": ok}


# ---------------------------------------------------------------- main
def main() -> None:
    import test_paper59_chain_slater_slice as SL
    import test_paper59_density_ring as DR
    mp.mp.dps = 30
    rec: Dict = {}
    N = 14

    # series, with the odd-power check
    Js, chis = J_series(N), chi_series(N)
    fJ = coeffs(Js, N)
    print("J exponents:", sorted(fJ), flush=True)

    # C3: representations
    vJ1 = B.sym_power([0, z - 5, 6], 3, 1)
    vJ2 = B.sym_power([0, z - 5, 6], 3, 2)
    vJ3 = B.sym_power([0, z - 5, 6], 3, 3)
    vchi = B.sym_power([0, z - 3, 4], 3, 2)
    QG = B.mellin_t_ode(vJ2, 2, 0, lambda sg: sg**2 - 1)
    QCc = B.mellin_t_ode(vchi, -2, 3, lambda sg: sg**2)

    def same(A_, B_):
        r = sp.cancel(A_[0] / B_[0])
        B_ = [rebind(b_) for b_ in B_]
        r = sp.cancel(A_[0] / B_[0])
        return len(A_) == len(B_) and all(sp.expand(a_ - r * b_) == 0 for a_, b_ in zip(A_, B_))
    rec["C3"] = {"slice_QG": same(QG, SL.QG), "chain_QC": same(QCc, DR.QC),
                 "v0_J2": str(vJ2[0]), "v0_chi2": str(vchi[0]), "v0_J1": str(vJ1[0])}
    print("C3:", rec["C3"], flush=True)

    # C1': G1 = int x K1(t x) J dx; strip 2 < Re s < 3
    Q1 = B.mellin_t_ode(vJ1, 2, 0, lambda sg: sg**2 - 1)
    S1, cr1 = source(Q1, fJ, 2, 0, E_K(1), sp.Rational(5, 2))
    S1n = sp.lambdify(t, sp.expand_func(S1).subs(Lt, sp.log(t)), "mpmath")
    rows = {}
    for tv in ("0.7", "1.9"):
        q = LF_num(Q1, Jnum, 1, 1, 0, tv)
        rows[tv] = {"quad": mp.nstr(q, 25), "engine": mp.nstr(S1n(mp.mpf(tv)), 25)}
    rec["C1p"] = {"Q1_terms": len(Q1), "crossed": cr1, "S1": simplify_coeffs(S1), "values": rows}
    print("C1':", json.dumps(rec["C1p"], indent=1), flush=True)
    # C1'': drop the highest crossed pole
    pdrop = max(p for (_, p, _, _) in cr1)
    S1b, _ = source(Q1, fJ, 2, 0, E_K(1), sp.Rational(5, 2), drop=pdrop)
    S1bn = sp.lambdify(t, sp.expand_func(S1b).subs(Lt, sp.log(t)), "mpmath")
    rec["C1pp"] = {"dropped_pole": pdrop,
                   "values": {tv: mp.nstr(S1bn(mp.mpf(tv)), 20) for tv in ("0.7", "1.9")}}
    print("C1'':", rec["C1pp"], flush=True)
    OUT.write_text(json.dumps(rec, indent=1, default=str))

    # WG1: slice
    fJ2 = coeffs(sp.expand(Js**2), N)
    S, cr = source(QG, fJ2, 2, 0, E_K(1), sp.Integer(3))
    rec["WG1"] = {"crossed": cr, "S": simplify_coeffs(S), "S_is_zero": sp.simplify(sp.expand_func(S)) == 0}
    print("WG1:", json.dumps(rec["WG1"], indent=1), flush=True)
    OUT.write_text(json.dumps(rec, indent=1, default=str))

    # extra check, not registered: the double-pole branch, G0 = int x K0(t x) J^2 dx, strip 2 < Re s < 4
    Q0 = B.mellin_t_ode(vJ2, 2, 0, lambda sg: sg**2)
    S0, cr0 = source(Q0, fJ2, 2, 0, E_K(0), sp.Integer(3))
    S0n = sp.lambdify(t, sp.expand_func(S0).subs(Lt, sp.log(t)), "mpmath")
    rec["X_logpole_G0"] = {"crossed": cr0, "S0": simplify_coeffs(S0),
                           "values": {tv: {"quad": mp.nstr(LF_num(Q0, lambda xx: Jnum(xx) ** 2, 0, 1, 0, tv), 25),
                                           "engine": mp.nstr(S0n(mp.mpf(tv)), 25)} for tv in ("0.7", "1.9")}}
    print("X G0:", json.dumps(rec["X_logpole_G0"], indent=1), flush=True)
    OUT.write_text(json.dumps(rec, indent=1, default=str))

    # WG3a: chain, strip: A(-2 - s), chi^2 ~ x^2: -2 < -2 - s < 0 -> -2 < s < 0; E0(s + 3) needs s > -3
    fc2 = coeffs(sp.expand(chis**2), N)
    Sc, crc = source(QCc, fc2, -2, 3, E_K(0), sp.Integer(-1), pref=8 / sp.pi**3)
    rec["WG3_source"] = {"crossed": crc, "S": simplify_coeffs(Sc), "S_is_zero": sp.simplify(sp.expand_func(Sc)) == 0}
    print("WG3 source:", json.dumps(rec["WG3_source"], indent=1), flush=True)
    OUT.write_text(json.dumps(rec, indent=1, default=str))

    # C2 + WG2 + WG3b: oscillation
    rec["C2"] = {"bessel": oscillation([t**2 - 1, 0, 1], sp.Rational(1, 2)),
                 "modified": oscillation([-(t**2 + 1), 0, 1], sp.Rational(1, 2))}
    print("C2:", rec["C2"], flush=True)
    rec["WG2"] = oscillation(SL.N2, sp.Integer(5))
    print("WG2:", rec["WG2"], flush=True)
    N2c = DR._chain_pieces()[3]
    rec["WG3_N2"] = [str(e) for e in N2c]
    rec["WG3_osc"] = oscillation(N2c, sp.Integer(2))
    print("WG3 osc:", rec["WG3_osc"], flush=True)
    OUT.write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
