"""QG step 6: rung tests on the Mellin-route operator L_Y of the geminal triangle (apex), gamma given.
Reference families as t-series mod p (density_loop_D_gauge): f2 = 2F1(1/2,1/2;1;t^2) (Legendre,
the two-density elliptic rung), f3 = 3F2(1/2,1/2,1/2;1,1;t^2/4) (F_3's operator 4 th^3 - t^2 (th+1)^3,
the K3 rung).  Side (R): exists A = sum_{j<r} a_j theta^j, a_j = p_j(t)/den, with L (A f) = 0 (the
reference is a sub-module = right factor up to gauge); side (L): the same on the adjoint (quotient).
Denominators den = t^a (t^2-4)^c (t^2-12)^e (t^2+20)^h.  Positive controls: F_3's operator and Legendre
must find themselves; the chain operator W = R3 R1 (theta-2) (Paper 59) must contain f3 on side (R);
F_4's operator must find neither.
usage: python ring_ext_QG_6_factor.py <path-to-operator-json-or-'mellin'> <gamma>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C
import density_loop_D_gauge as D
import density_loop_D4_factors as F4
import ring_ext_QG_ops as O

u = sp.Symbol("u")
t = sp.Symbol("t")
P = C.P
NT = F4.NT


def inv_series(den_poly_t, N):
    cs = [D.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))) for c in sp.Poly(den_poly_t, t).all_coeffs()[::-1]]
    out = [pow(cs[0], P - 2, P)]
    for n in range(1, N):
        s_ = sum(cs[k] * out[n - k] for k in range(1, min(n, len(cs) - 1) + 1)) % P
        out.append((-s_) * out[0] % P)
    return {n: c for n, c in enumerate(out) if c}


def mul(a, b, top):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            e = e1 + e2
            if e < top:
                out[e] = (out.get(e, 0) + c1 * c2) % P
    return {e: c for e, c in out.items() if c}


def test(Q, f, r, Dg, den, a=6):
    w = inv_series(sp.expand(den), NT)
    top = NT - 2 * len(Q) - Dg - 6
    cols = []
    for j in range(r):
        g = mul(w, D.theta(f, j), NT)
        for d in range(-a, Dg + 1):
            cols.append(F4.apply_L(Q, D.shift(g, d)))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    return len(cols) - C.rank_mod(rows)


DENS = [sp.Integer(1), (t**2 - 4), (t**2 - 4) ** 2, (t**2 - 4) ** 3, (t**2 - 4) ** 2 * (t**2 - 12) * (t**2 + 20),
        (t**2 - 4) ** 3 * (t**2 - 12) ** 2 * (t**2 + 20) ** 2, (t**2 - 4) ** 4 * (t**2 - 12) ** 3 * (t**2 + 20) ** 3,
        (t**2 - 4) ** 2 * (4 * t**2 + 11), (t**2 - 4) ** 3 * (4 * t**2 + 11) ** 2]


def scan(name, Q, refs, Dg=34):
    res = {}
    for ref, f, r in refs:
        for side, QQ in [("R", Q), ("L", F4.adjoint(Q))]:
            best = 0
            for den in DENS:
                n = test(QQ, f, r, Dg, den)
                if n:
                    best = n
                    break
            res[f"{ref}|{side}"] = best
            print(f"{name:28s} {ref} ({side}): {'FOUND (nullity %d)' % best if best else 'none'}", flush=True)
    return res


def main():
    f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT)
    f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4))
    refs = [("f3", f3, 3), ("f2", f2, 2)]
    rec = {}
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("controls", "all"):
        Q3 = [4 * u**3, -(u + 1) ** 3]
        Q2 = [u**2, -(u + 1) ** 2]
        LCd = json.loads((Path(__file__).parent / "data" / "linking_ceiling_B_ode.json").read_text())
        Q4 = [sp.sympify(q) for q in LCd["4"]["Q"]]
        R3 = [sp.sympify(x) for x in json.loads((Path(__file__).parent / "data" / "density_loop_D10.json").read_text())["R3"]]
        import density_loop_D10_k3_exact_check as X
        import density_loop_D11_chain_minimal as D11
        Wc = X.compose(X.compose(R3, [-(3 + X.t**2), 1 - X.t**2]), [sp.Integer(-2), sp.Integer(1)])
        QW, _ = D11.to_Q(Wc)
        QW = [sp.expand(q.subs(D11.u, u)) for q in QW]
        for nm, Q in [("control F_3 op", Q3), ("control Legendre", Q2), ("control F_4 op", Q4), ("control chain W", QW)]:
            rec[nm] = scan(nm, Q, refs)
    if which in ("op", "all"):
        g = sys.argv[2] if len(sys.argv) > 2 else "1"
        for r in (2, 4):
            path = Path(__file__).parent / "data" / f"ring_ext_QG_5_guess_phi_g{g}.json"
            try:
                Qphi = O.load_Q(str(path), r)
            except Exception:
                continue
            QY = O.mellin_Y_from_phi(Qphi)
            rec[f"L_Y(from L_phi order {r})"] = scan(f"L_Y g={g} (L_phi ord {r})", QY, refs)
    out = Path(__file__).parent / "data" / f"ring_ext_QG_6_factor_{which}.json"
    out.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
