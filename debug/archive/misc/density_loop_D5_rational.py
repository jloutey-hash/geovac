"""Loop test PT2, step D5: the factor test of D4 with rational change-of-basis coefficients,
a_j = p_j(t) / (t^a (t^2-1)^b (t^2-4)^c), p_j polynomials (Laurent shifts d >= -a), b, c <= 3.
Controls as in D4, plus: the triangle must show f2 as a quotient (known from step D)."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C
import density_loop_D_gauge as D
import density_loop_D4_factors as F4
u = sp.Symbol("u"); P = C.P; NT = F4.NT
OUT = Path(__file__).parent / "data" / "density_loop_D5.json"

def inv_series(den_poly_t, N):
    """1/den as t-series mod p, den a polynomial in t with den(0) != 0."""
    cs = [D.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))) for c in sp.Poly(den_poly_t, sp.Symbol("t")).all_coeffs()[::-1]]
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

def test(Q, f, r, Dg, b, c, a=4):
    t = sp.Symbol("t")
    w = inv_series(sp.expand((t**2 - 1) ** b * (t**2 - 4) ** c), NT)
    top = NT - 2 * len(Q) - Dg - 6
    cols = []
    for j in range(r):
        g = mul(w, D.theta(f, j), NT)
        for d in range(-a, Dg + 1):
            cols.append(F4.apply_L(Q, D.shift(g, d)))
    exps = sorted(e for e in set().union(*cols) if e < top)
    rows = [[col.get(e, 0) for col in cols] for e in exps]
    return len(cols) - C.rank_mod(rows)

def main():
    B = json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())
    QG = [sp.sympify(q) for q in B["QG"]]; QC = [sp.sympify(q) for q in B["QC"]]
    LCd = json.loads((Path(__file__).parent / "data" / "linking_ceiling_B_ode.json").read_text())
    Q4 = [sp.sympify(q) for q in LCd["4"]["Q"]]
    f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT)
    f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4))
    rec = {}
    for name, Q in [("control F_4 op", Q4), ("triangle L_G", QG), ("chain L_C", QC)]:
        for ref, f, r in [("f3", f3, 3), ("f2", f2, 2)]:
            for side, QQ in [("R", Q), ("L", F4.adjoint(Q))]:
                best = 0
                for b in range(4):
                    for c in range(4):
                        n = test(QQ, f, r, 10, b, c)
                        best = max(best, n)
                        if n:
                            print(f"{name:16s} {ref} ({side}) denominators (t^2-1)^{b} (t^2-4)^{c}: {n}", flush=True)
                            break
                    if best:
                        break
                rec[f"{name}|{ref}|{side}"] = best
                print(f"{name:16s} {ref} ({side}): {'FOUND' if best else 'none'}", flush=True)
    OUT.write_text(json.dumps(rec, indent=1))

if __name__ == "__main__":
    main()
