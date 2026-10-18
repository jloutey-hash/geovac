"""nosplit step 4: exact argument that the chain C lies on the order-4 factor U = ker P.

Facts used:  Sol(W) = U (+) <g2>, g2 the artanh solution of E (log at x = 1);  T_1 - 1 (local
monodromy at t = 1) has rank 1 on Sol(W) (step-1 reasoning; checked here by log degrees).
If every local solution of P at x = 1 is log-free (integer exponents, no logs), then U is
contained in ker(T_1 - 1), which has dimension 4, so U = ker(T_1 - 1) = {solutions of W that
are single-valued around t = 1}.  The chain C(t) = (8t^3/pi^3) int chi^2 K0(t x) dx is analytic
on Re t > 0, hence single-valued around t = 1, hence C in U, hence P C = 0 exactly.
Log structure at finite p via the local Euler form  z^n L = sum_k a_k(p+z) z^{n-k} [theta_z]_k."""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import nosplit_core as N
import nosplit_search as NS
import nosplit_local_logs as LL
import nosplit_order4 as O4

x, r = N.x, N.r
z = sp.Symbol("z")


def local_layers(L, p):
    a = NS.to_d_form(L)
    n = len(a) - 1
    terms = {}
    for k, c in enumerate(a):
        ff = sp.prod([r - j for j in range(k)])
        terms[k] = (sp.cancel(c.subs(x, p + z) * z ** (n - k)), ff)
    den = sp.lcm([sp.denom(sp.together(c)) for c, _ in terms.values()])
    P = {}
    for k, (c, ff) in terms.items():
        pc = sp.Poly(sp.expand(sp.cancel(c * den)), z)
        for (e,), v in zip(pc.monoms(), pc.coeffs()):
            P[e] = P.get(e, 0) + v * ff
    lo = min(P)
    return {i - lo: sp.expand(q) for i, q in P.items()}


def log_structure_at(L, p, rho0, nmax, K=5):
    LL_layers = LL.layers
    LL.layers = lambda Lop, pp: local_layers(Lop, p)
    try:
        return LL.log_structure(L, "finite", rho0, nmax, K)
    finally:
        LL.layers = LL_layers


def main():
    ops = N.chain_ops()
    R3, R1, T2 = ops["R3"], ops["R1"], ops["T2"]
    E = N.compose(R1, T2)
    W = N.compose(R3, E)
    S1 = O4.S1_exact()
    P = [sp.cancel(a - b) for a, b in zip(T2 + [0] * 3, N.compose(S1, E))]
    rec = {}
    for name, L in [("E", E), ("W", W), ("P", P), ("R3", R3)]:
        I, reg, ex = N.exponents(L, sp.Integer(1))
        lo, hi = min(ex), max(ex)
        dim, ml = log_structure_at(L, sp.Integer(1), lo, int(hi - lo))
        print(f"{name} at x=1: exponents {sorted(ex)}, group dim {dim}, max log degree {ml}", flush=True)
        rec[name] = {"exps": [str(e) for e in sorted(ex)], "dim": dim, "maxlog": ml}
    # also P's exponents everywhere (for step 5)
    print("P singular factors:", N.singular_points(P), flush=True)
    (N.HERE / "data" / "nosplit_membership.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
