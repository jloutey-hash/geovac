"""nosplit step 1: exact log structure of local solutions at x = 0 (t = 0) and x = oo.
Generalised Frobenius: y = sum_n x^(rho0+n) q_n(log x), q_n polynomials; P(delta) acts on
x^rho q(L) as x^rho P(rho + d/dL) q.  Exact linear algebra on layers 0..n_max (beyond the
last exponent of the group the layers are forced and cannot raise the log degree).
Reports, per operator and point, the dimension of the exponent group and the maximal log degree.
A split W = Sol(E) (+) V with V ~ Sol(R3) forces maxlog(W) <= max(maxlog(E), maxlog(R3))."""
from __future__ import annotations
import json, sys
import sympy as sp
import nosplit_core as N

x, r = N.x, N.r
Lg = sp.Symbol("Lg")


def layers(Lop, p):
    """Return (P_i polynomials in r) with Lop = sum_i z^i P_i(delta_z), z the local variable
    (x at 0, w = 1/x at oo; delta_w = -delta_x), after clearing a common factor."""
    w = sp.Symbol("w")
    if p == 0:
        cs = [sp.cancel(c) for c in Lop]; z = x; sgn = 1
    else:
        cs = [sp.cancel(c.subs(x, 1 / w)) for c in Lop]; z = w; sgn = -1
    den = sp.lcm([sp.denom(sp.together(c)) for c in cs])
    cs = [sp.expand(sp.cancel(c * den)) for c in cs]
    lo = min(min(m[0] for m in sp.Poly(c, z).monoms()) for c in cs if c != 0)
    P = {}
    for k, c in enumerate(cs):
        pc = sp.Poly(c, z)
        for (e,), v in zip(pc.monoms(), pc.coeffs()):
            P[e - lo] = P.get(e - lo, 0) + v * (sgn * r) ** k
    return {i: sp.expand(q) for i, q in P.items()}


def act(P, rho, q):
    """P(rho + D) q, D = d/dLg."""
    Pp = sp.Poly(P, r)
    out = 0
    shifted = sp.Poly(Pp.as_expr().subs(r, r + rho), r)
    for (k,), c in zip(shifted.monoms(), shifted.coeffs()):
        out += c * sp.diff(q, Lg, k)
    return sp.expand(out)


def log_structure(Lop, p, rho0, nmax, K=5):
    P = layers(Lop, p)
    unk, ys = [], []
    for n in range(nmax + 1):
        cs = sp.symbols(f"a{n}_0:{K + 1}")
        unk += cs
        ys.append(sum(c * Lg**j for j, c in enumerate(cs)))
    eqs = []
    for n in range(nmax + 1):
        tot = 0
        for i, Pi in P.items():
            if n - i >= 0:
                tot += act(Pi, rho0 + n - i, ys[n - i])
        tot = sp.expand(tot)
        eqs += [tot.coeff(Lg, j) for j in range(K + 1)]
    M = sp.Matrix([[sp.diff(e, u) for u in unk] for e in eqs])
    ns = M.nullspace()
    dim = len(ns)
    maxlog = -1
    for v in ns:
        sub = dict(zip(unk, v))
        for n in range(nmax + 1):
            q = sp.expand(ys[n].subs(sub))
            if q != 0:
                maxlog = max(maxlog, sp.Poly(q, Lg).degree())
    return dim, maxlog


def main():
    ops = N.chain_ops()
    E = N.compose(ops["R1"], ops["T2"])
    W = N.compose(ops["R3"], E)
    W4 = N.compose(ops["R3"], ops["R1"])
    rec = {}
    cases = [("W", W, 0, sp.Rational(1, 2), 2), ("W", W, 0, sp.Integer(1), 0),
             ("R3", ops["R3"], 0, sp.Rational(1, 2), 2), ("E", E, 0, sp.Integer(1), 1),
             ("L3", ops["L3"], 0, sp.Integer(0), 0),
             ("W", W, "oo", sp.Integer(-1), 3), ("W", W, "oo", sp.Rational(1, 2), 0),
             ("R3", ops["R3"], "oo", sp.Integer(-1), 2), ("E", E, "oo", sp.Integer(-1), 0),
             ("W4=R3oR1", W4, 0, sp.Rational(1, 2), 2), ("W4=R3oR1", W4, "oo", sp.Integer(-1), 2)]
    for name, Lop, p, rho0, nmax in cases:
        dim, ml = log_structure(Lop, p, rho0, nmax)
        print(f"{name:10s} at {str(p):3s} group rho0={rho0}, layers 0..{nmax}: dim {dim}, max log degree {ml}", flush=True)
        rec[f"{name}@{p}@{rho0}"] = {"dim": dim, "maxlog": ml}
    (N.HERE / "data" / "nosplit_local_logs.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
