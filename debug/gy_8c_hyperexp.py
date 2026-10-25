"""GY stage 8c: order-one right factors of an IRREGULAR piece F, allowing the exponential factor its Newton polygon permits.
A hyperexponential solution s of F has the form exp(c x) t^f0 (1 - x/4)^e4 P(t), with c a root of the slope -1 edge's
characteristic polynomial (exp(c t^2)), f0 a t = 0 exponent, e4 an exponent at x = 4 (t = +-2), P a polynomial (the finite
singular points are regular, so no other exponential part).  Conjugating, F_c = exp(-c x) o F o exp(c x) = F(theta + 2 c x)
(theta_t (c x) = 2 c x), and s' = s exp(-c x) solves F_c; the exponent-bounded test of lib4 (test_candidate) is then run on
F_c with a STATED degree bound Nmax for P (the infinity bound of lib4 is not valid for an irregular piece).
Applied to F and to F* (order-one left factors).  With F = F1, F2 of Y_G's chain this closes the escape of an F_3-type piece
sitting in the middle of F2 (which would force an order-one sub- or quotient factor).
Control: the Slater G4 at gamma = 1 (order 4 = R3 o R1, R1 with the algebraic solution t^3 (...)/(t^2-4)^(9/2)) must give a hit
at c = 0.
usage: python gy_8c_hyperexp.py <pieces json> <key> <Nmax> <c list> <out tag> | python gy_8c_hyperexp.py control"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib4 as S
import gy_lib as A


def twist(L, c):
    """L(theta + 2 c x) as an fmpz_poly operator (primitive)"""
    c = Fraction(c)
    if c == 0:
        return L
    Lq = [flint.fmpq_poly([flint.fmpq(int(v)) for v in ck.coeffs()]) if not ck.is_zero() else flint.fmpq_poly([]) for ck in L]
    T = [flint.fmpq_poly([0, flint.fmpq(2 * c.numerator, c.denominator)]), flint.fmpq_poly([1])]   # theta + 2 c x
    res = [flint.fmpq_poly([])]
    Pk = [flint.fmpq_poly([1])]
    for k, ck in enumerate(Lq):
        if not ck.is_zero():
            res = G.add_op(res, G.scale_op(ck, Pk))
        if k < len(Lq) - 1:
            Pk = G.compose(T, Pk)
    return G.op_to_fmpz(G.trim(res))


def times_inverse(L, den, m):
    """polynomial operator K with L o den^(-m) = den^(-(m + ord L)) K  (theta^i (1/D) = N_i / D^(i+1), D = den^m; as in
    correlated_triangle_g32_lib2.verify_map)"""
    from math import comb
    D = den ** m
    K_ = len(L) - 1
    Nl = [flint.fmpz_poly([1])]
    tD = G.Dth(D)
    for i in range(K_):
        Nl.append(G.Dth(Nl[-1]) * D - (i + 1) * Nl[-1] * tD)
    out = [flint.fmpz_poly([]) for _ in range(K_ + 1)]
    for k, c in enumerate(L):
        if c.is_zero():
            continue
        for i in range(k + 1):
            out[k - i] = out[k - i] + c * comb(k, i) * Nl[i] * D ** (K_ - i)
    return G.content_free(G.trim(out))


def class_mins(vals):
    return S.class_mins(vals)


def search(L, cs, Nmax, label):
    out = []
    for c in cs:
        Lc = twist(L, c)
        e0 = class_mins(A.Lc.exps_t0(Lc)[0])
        e4 = class_mins(A.Lc.exps_finite_t(Lc, Fraction(4))[0])
        hits = []
        ncand = 0
        for f0 in e0:
            for e in e4:
                combo = {a: (e if abs(a) == 2 else Fraction(0)) for a in S.TPTS}
                ncand += 1
                nul = S.test_candidate(Lc, f0, combo, Nmax, "sol")
                if nul:
                    hits.append({"f0": str(f0), "e4": str(e), "nullity": nul})
        print(f"   {label}: c = {c}: {ncand} exponent choices, P degree <= {Nmax}: hits {hits}", flush=True)
        out.append({"c": str(c), "candidates": ncand, "Nmax": Nmax, "hits": hits})
    return out


if __name__ == "__main__":
    A.set_points(["1", "2"])
    if sys.argv[1] == "control":
        G4 = G.load_G4()
        r = search(G4, ["0", "1/8"], 60, "Slater G4 (gamma = 1)")
        ok = any(h["hits"] for h in r if h["c"] == "0")
        print("control: order-one factor of G4 found at c = 0:", ok, flush=True)
        (G.DATA / "gy_8c_hyperexp_control.json").write_text(json.dumps({"runs": r, "found_at_c0": ok}, indent=1))
    else:
        d = json.loads(Path(sys.argv[1]).read_text())
        key, Nmax, cs, tag = sys.argv[2], int(sys.argv[3]), sys.argv[4].split(","), sys.argv[5]
        F = G.parse_op(d[key])
        t0 = time.time()
        Fa = G.content_free(G.adjoint(F))
        rec = {"key": key, "F": search(F, cs, Nmax, key), "F_adjoint": search(Fa, [str(-Fraction(c)) for c in cs], Nmax, key + "*")}
        # adjoint solutions may have poles at F's apparent points (exponents down to -mpole there): allow app^(-m)
        import cf_g2_lib as C
        app = flint.fmpz_poly([1])
        for f, e in C.lead_factors(Fa):
            if f.degree() >= 2:
                app = app * f ** e
        mpole = int(sys.argv[6]) if len(sys.argv) > 6 else 2
        rec["F_adjoint_with_app_poles"] = {}
        for m in range(1, mpole + 1):
            Km = times_inverse(Fa, app, m)
            rec["F_adjoint_with_app_poles"][str(m)] = search(Km, [str(-Fraction(c)) for c in cs], Nmax, key + f"* o app^-{m}")
        rec["apparent_factor_used"] = str(app)
        rec["seconds"] = time.time() - t0
        (G.DATA / f"gy_8c_hyperexp_{tag}.json").write_text(json.dumps(rec, indent=1))
        print(f"({time.time()-t0:.0f}s)", flush=True)
