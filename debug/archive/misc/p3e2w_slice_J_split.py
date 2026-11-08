"""P3e-3 step J: does the K3 piece split off?  W = R3 o R1 o (theta+2) (step I) has K3 (R3) as its
top factor.  A split (K3 also a sub) means a rational gauge B = b0 + b1 th + b2 th^2 with W(B f) = 0
for f in Sol(L3), L3 = 4 th^3 - t^2 (th+1)^3.  Search, mod p, with ONE maximal denominator
  t^a (t^2-1)^b (t^2-4)^c q(t)^e   (q = apparent-singularity factor(s) of R3)
and numerators of degree <= Dg (a maximal box contains every smaller one).  The box is compared with
pole/degree bounds from the local exponents of W and L3 (printed).
Positive control: E2 o R3 (elementary piece below K3 reversed: the split exists by construction)
must be found by the same code at the same box.  Second control: the chain's W (density loop D12,
no split there) must give none at the same box."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402
import density_loop_D10_k3_exact_check as D10  # noqa: E402
import density_loop_D11_chain_minimal as D11  # noqa: E402
import p3e2w_slice_D_exponents as DX  # noqa: E402

t, u = D10.t, D10.u
P = C.P
HERE = Path(__file__).parent
NT = 700
F4.NT = NT


def inv_series(den, N):
    cs = [D.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))) for c in sp.Poly(den, t).all_coeffs()[::-1]]
    assert cs[0] != 0
    out = [pow(cs[0], P - 2, P)]
    for n in range(1, N):
        s_ = sum(cs[k] * out[n - k] for k in range(1, min(n, len(cs) - 1) + 1)) % P
        out.append((-s_) * out[0] % P)
    return {n: c for n, c in enumerate(out) if c}


def mul(a, b, top):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            if e1 + e2 < top:
                out[e1 + e2] = (out.get(e1 + e2, 0) + c1 * c2) % P
    return {e: c for e, c in out.items() if c}


def test(Q, den, a, Dg):
    f = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4))
    w = inv_series(den, NT)
    top = NT - 2 * len(Q) - Dg - a - 8
    cols = []
    for j in range(3):
        g = mul(w, D.theta(f, j), NT)
        for d in range(-a, Dg + 1):
            cols.append(F4.apply_L(Q, D.shift(g, d)))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    return len(cols) - C.rank_mod(rows), len(cols), len(rows)


def main():
    H = json.loads((HERE / "data" / "p3e2w_slice_H.json").read_text())
    R3 = [sp.sympify(x) for x in H["R3"]]
    R1 = [-(1 + 7 * t**2), 1 - t**2]
    E2 = D10.compose(R1, [sp.Integer(2), sp.Integer(1)])
    W = D10.compose(R3, E2)
    QW, _ = D11.to_Q(W)
    Qc, _ = D11.to_Q(D10.compose(E2, R3))
    QchainW, _ = D11.to_Q(D11.W)
    # apparent-singularity factors of R3 (denominator factors other than t, t^2-1, t^2-4)
    dens = sp.factor_list(sp.lcm([sp.denom(sp.together(c)) for c in R3]))[1]
    q = sp.Integer(1)
    for fac, m in dens:
        if sp.Poly(fac, t).degree() and sp.rem(sp.expand((t**2 - 1) * (t**2 - 4) * t), fac, t) != 0:
            q *= fac
    print("R3 denominator factors:", dens, " apparent part q =", q, flush=True)
    # local exponents (printed for the bound)
    for name, Q in [("W", QW), ("L3", [4 * u**3, -(u + 1) ** 3])]:
        Lc = DX.theta_coeffs(Q)
        print(f"{name}: t=0 {sp.factor(Q[0])}, inf {sp.factor(Q[-1])}", flush=True)
        for t0 in (1, 2):
            print(f"    t={t0}: {DX.indicial_at(Lc, t0)[0]}", flush=True)
        if q != 1:
            for r_ in sp.Poly(q, t).all_roots()[:1]:
                pass
    a, b, c, e, Dg = (int(x) for x in sys.argv[1:6]) if len(sys.argv) > 5 else (12, 6, 6, 3, 40)
    den = sp.expand((t**2 - 1) ** b * (t**2 - 4) ** c * q**e)
    rec = {"box": {"a": a, "b": b, "c": c, "e": e, "Dg": Dg, "q": str(q)}}
    for name, Q in [("positive control E2 o R3", Qc), ("W = R3 R1 (th+2)", QW), ("chain W (D12)", QchainW)]:
        t0 = time.time()
        n, nc, nr = test(Q, den, a, Dg)
        rec[name] = n
        print(f"{name:28s}: nullity {n}   [{nc} unknowns, {nr} rows, {time.time()-t0:.0f}s]", flush=True)
    (HERE / "data" / "p3e2w_slice_J.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
