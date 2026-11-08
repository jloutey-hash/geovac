"""GS step C (record docs/predictions/2026-10-21_generic_slice.md): E2 = LCLM(R1a, R1b) (order 2),
E3 = E2 (theta + 2), M14 = L / E3 exactly; then does V's operator R_V (step A, order 4) sit in M14 as a
sub (R side: a rational gauge maps Sol(R_V) into Sol(M14)) or a quotient (L side, on adjoints)?
Mod-p series test as density_loop_D5_rational, one maximal denominator
t^-a (t^2-1)^B (t^2-4)^B (t^2-9)^B, numerators to degree Dg (a maximal box contains every smaller one).
Reference: the exponent-2 Frobenius solution of R_V at t = 0 (rational coefficients).
Controls: R_V itself (both sides must hit); F_4's operator (Paper 59, neither side)."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D4_factors as F4  # noqa: E402
from generic_slice_B import coeff_list, Q, rdiv, compose, cth, t, u  # noqa: E402

HERE = Path(__file__).parent
P = C.P
NT = 760
F4.NT = NT
A = json.loads((HERE / "data" / "generic_slice_A.json").read_text())
QV = [sp.sympify(q, locals={"u": u}) for q in A["QV"]]


def lclm2(ha, hb):
    """Monic-free order-2 operator killing ha, hb: det [[f, th f, th^2 f], [ha ...], [hb ...]]."""
    rows = [[ha, cth(ha), cth(cth(ha))], [hb, cth(hb), cth(cth(hb))]]
    c0 = sp.cancel(rows[0][1] * rows[1][2] - rows[0][2] * rows[1][1])
    c1 = sp.cancel(-(rows[0][0] * rows[1][2] - rows[0][2] * rows[1][0]))
    c2 = sp.cancel(rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0])
    den = sp.lcm([sp.fraction(c)[1] for c in (c0, c1, c2)])
    out = [sp.factor(sp.cancel(c * den)) for c in (c0, c1, c2)]
    g = sp.gcd_list([sp.fraction(sp.cancel(c))[0] for c in out])
    return [sp.expand(sp.cancel(c / g)) for c in out]


def to_Q(Lc):
    den = sp.lcm([sp.fraction(sp.cancel(c))[1] for c in Lc])
    pol = [sp.expand(sp.cancel(c * den)) for c in Lc]
    deg = max(sp.Poly(p, t).degree() for p in pol if p != 0)
    assert all(sp.Poly(p, t).as_expr().subs(t, -t).expand() == p for p in pol), "not even in t"
    Qo = []
    for i in range(deg // 2 + 1):
        Qo.append(sp.expand(sum(sp.Poly(p, t).coeff_monomial(t ** (2 * i)) * u**k for k, p in enumerate(pol))))
    while Qo and Qo[-1] == 0:
        Qo.pop()
    return Qo


def inv_series(den, N):
    cs = [D.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))) for c in sp.Poly(den, t).all_coeffs()[::-1]]
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


def test(Qo, f, r, Bp, a, Dg):
    den = sp.expand((t**2 - 1) ** Bp * (t**2 - 4) ** Bp * (t**2 - 9) ** Bp)
    w = inv_series(den, NT)
    top = min(f) + NT - 2 * len(Qo) - Dg - a - 8
    cols = []
    for j in range(r):
        g = mul(w, D.theta(f, j), NT)
        for d in range(-a, Dg + 1):
            cols.append(F4.apply_L(Qo, D.shift(g, d)))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    return len(cols) - C.rank_mod(rows), len(cols), len(rows)


def main(Bp=3, a=6, Dg=30):
    rec = {}
    t0 = time.time()
    L = coeff_list(Q)
    M16, rem = rdiv(L, [sp.Integer(2), sp.Integer(1)])
    assert all(r == 0 for r in rem)
    E2 = lclm2(t / (1 - t**2) ** 4, t / (4 - t**2) ** 4)
    print("E2:", [sp.factor(c) for c in E2], flush=True)
    M14, rem = rdiv(M16, E2)
    ok = all(sp.cancel(r) == 0 for r in rem)
    print(f"L = M14 o E2 o (theta+2) exactly: {ok}; M14 order {len(M14)-1}  ({time.time()-t0:.0f}s)", flush=True)
    rec["E2"] = [str(c) for c in E2]
    rec["E3_exact"] = ok
    QM = to_Q(M14)
    rec["M14_Q"] = [str(q) for q in QM]
    print("M14 in Q-form: t-degree", 2 * (len(QM) - 1), flush=True)
    (HERE / "data" / "generic_slice_C.json").write_text(json.dumps(rec, indent=1))
    fV = D.frob_t(QV, 2, NT=NT)
    LCd = json.loads((HERE / "data" / "linking_ceiling_B_ode.json").read_text())
    Q4 = [sp.sympify(q, locals={"u": u}) for q in LCd["4"]["Q"]]
    for name, Qo in [("control R_V", QV), ("control F_4", Q4), ("M14", QM), ("L/(th+2)", to_Q(M16)), ("L", Q)]:
        for side, QQ in [("R", Qo), ("L", F4.adjoint(Qo))]:
            t1 = time.time()
            n, nc, nr = test(QQ, fV, 4, Bp, a, Dg)
            rec[f"{name}|{side}"] = n
            print(f"{name:10s} ({side}): nullity {n}  [{nc} cols, {nr} rows, {time.time()-t1:.0f}s]", flush=True)
            (HERE / "data" / "generic_slice_C.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
