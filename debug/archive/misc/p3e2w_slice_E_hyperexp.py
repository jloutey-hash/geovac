"""P3e-3 step E: algebraic (hyperexponential) solutions of an operator sum_i t^{2i} Q_i(theta), i.e.
its order-1 right factors with algebraic solutions.  Any such y is log-free at t = 0, so y lies in
the span of the log-free Frobenius solutions y_lam (free slots set to 0, density_loop_D6.frob_free).
Ansatz y (1-t^2)^{-al} (1-t^2/4)^{-be} = polynomial, al in exps(1), be in exps(2) (step D), degree
<= mu_max - 2 al - 2 be with mu_max the largest exponent at infinity.  Linear in the span
coefficients; mod p (p = 2^61-1), series length N.
Control: the chain operator's M6 (density loop) must return R1's solution t^3/(1-t^2)^2."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C  # noqa: E402
import density_loop_D_gauge as D  # noqa: E402
import density_loop_D6_peel as F6  # noqa: E402

u = sp.Symbol("u")
P = C.P
HERE = Path(__file__).parent


def binom_series(alpha, scale, N):
    """(1 - scale t^2)^alpha as dict t-exp -> coeff mod p (alpha a Fraction)."""
    out, c = {}, Fraction(1)
    for n in range(N // 2 + 1):
        if c:
            out[2 * n] = D.fr(c)
        c = c * (Fraction(alpha) - n) / (n + 1) * (-Fraction(scale))
    return out


def mul(a, b, top):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            if e1 + e2 < top:
                out[e1 + e2] = (out.get(e1 + e2, 0) + c1 * c2) % P
    return out


def algebraic_solutions(Q, e1, e2, einf, N=360):
    roots = sorted(set(sp.roots(sp.Poly(Q[0], u).as_expr(), u)))
    sols = []
    for lam in roots:
        ser, free = F6.frob_free(Q, lam, NT=N)
        if ser is not None:
            sols.append((lam, ser))
    lo = min(l for l, _ in sols)
    found = []
    for al in e1:
        for be in e2:
            w = mul(binom_series(-Fraction(al), 1, N), binom_series(-Fraction(be), Fraction(1, 4), N), N + lo + 2)
            dmax = max(einf) - 2 * al - 2 * be
            cols = [mul(s_, w, N + lo - 4) for _, s_ in sols]
            rows = [[col.get(e, 0) for col in cols] for e in range(int(sp.floor(dmax)) + 1, N + int(lo) - 4)]
            if not rows:
                continue
            nul = len(cols) - C.rank_mod(rows)
            if nul:
                found.append((str(al), str(be), nul))
                print(f"   alpha={al} beta={be}: {nul} algebraic solution(s), degree <= {dmax}", flush=True)
    return [l for l, _ in sols], found


def main():
    rec = {}
    Bc = json.loads((HERE / "data" / "density_loop_B.json").read_text())
    import density_loop_D7_m5 as F7  # noqa: E402
    M6c = F7.rdiv([sp.sympify(q) for q in Bc["QC"]], [u - 2])
    print("control: chain M6 (expects alpha=-2 i.e. (1-t^2)^-2, beta=0)", flush=True)
    einf_c = [int(r) for r in sp.roots(sp.Poly(M6c[-1], u).as_expr(), u)]
    print("   chain M6 exponents at infinity:", einf_c, flush=True)
    rec["control"] = algebraic_solutions(M6c, [0, 1, 2, -1, -2, -3], [0, 1, 2, -1, -2], einf_c)
    M6 = [sp.sympify(q) for q in json.loads((HERE / "data" / "p3e2w_slice_B.json").read_text())["M6"]]
    print("M6 = L_G/(theta+2):", flush=True)
    rec["M6"] = algebraic_solutions(M6, [0, 1, 2, 3, 4, -4], [0, 1, 2, 3, 4, Fraction(-11, 2)], [-4, -7, -8, -12])
    (HERE / "data" / "p3e2w_slice_E.json").write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
