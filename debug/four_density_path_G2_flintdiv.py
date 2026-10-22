"""FP step G2: exact right division L = M o W12 over Q(w) with python-flint (fmpq_poly), w = t^2, delta = w d/dw.
Operators in delta form with polynomial coefficients (from the Q-form: sum_i w^i Q_i(2 delta)).
R^(m) = delta^m o R kept with polynomial coefficients; the running dividend is N_k / lead^e (common exponent e).
Remainder zero <=> all remainder numerators are zero (exact)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import flint
import sympy as sp

HERE = Path(__file__).parent
u = sp.Symbol("u")


def to_delta_poly(Q):
    """Q-form -> list over delta-powers of fmpq_poly in w."""
    n = max(sp.Poly(sp.sympify(q).subs(sp.Symbol("u"), u), u).degree() for q in Q) + 1
    coeffs = [dict() for _ in range(n)]
    for i, q in enumerate(Q):
        p = sp.Poly(sp.expand(sp.sympify(q).subs(sp.Symbol("u"), u).subs(u, 2 * u)), u)
        for (k,), c in p.terms():
            c = sp.Rational(c)
            coeffs[k][i] = coeffs[k].get(i, 0) + c
    out = []
    for d in coeffs:
        m = max(d) if d else 0
        lst = [flint.fmpq(0)] * (m + 1)
        for i, c in d.items():
            lst[i] = flint.fmpq(int(c.p), int(c.q))
        out.append(flint.fmpq_poly(lst))
    return out


def dpoly(a):
    """w a'(w)"""
    return a.derivative() * flint.fmpq_poly([0, 1])


def delta_compose_left(R):
    """delta o R (R list of polys): coeff_k -> d(a_k) + a_{k-1}."""
    out = [flint.fmpq_poly([0])] * (len(R) + 1)
    for k, a in enumerate(R):
        out[k] = out[k] + dpoly(a)
        out[k + 1] = out[k + 1] + a
    return out


def rdiv_exact(L, R, verbose=True):
    t0 = time.time()
    r = len(R) - 1
    lead = R[-1]
    Rm = [R]
    for m in range(len(L) - len(R)):
        Rm.append(delta_compose_left(Rm[-1]))
    N = [a for a in L]       # numerators over lead^e
    e = 0
    for k in range(len(L) - 1, r - 1, -1):
        m = k - r
        ck = N[k]            # quotient coeff = ck / lead^(e+1)
        if ck == 0:
            continue
        # N_i <- N_i * lead - ck * Rm[m][i]
        Rmm = Rm[m]
        newN = []
        for i in range(len(N)):
            v = N[i] * lead
            if i < len(Rmm):
                v = v - ck * Rmm[i]
            newN.append(v)
        N = newN
        e += 1
        assert N[k] == 0
        # remove common content periodically
        if e % 4 == 0:
            g = None
            for a in N[:k]:
                if a != 0:
                    g = a if g is None else g.gcd(a)
            if g is not None and g.degree() > 0:
                # only divide by the part of g dividing lead^e: keep it simple, divide numerators by gcd(g, lead^e)
                pass
        if verbose and (k % 8 == 0):
            print(f"   k={k}: e={e}, max deg {max(a.degree() for a in N[:k] if a != 0) if any(a != 0 for a in N[:k]) else -1}"
                  f" ({time.time()-t0:.0f}s)", flush=True)
    rem = N[:r]
    return all(a == 0 for a in rem), e


def main() -> None:
    QL = json.loads((HERE / "data" / "four_density_path_A.json").read_text())["L"]["Q"]
    QW = json.loads((HERE / "data" / "four_density_path_D.json").read_text())["W12"]["Q"]
    L, W = to_delta_poly(QL), to_delta_poly(QW)
    t0 = time.time()
    ok, e = rdiv_exact(L, W)
    print("L = M o W12 exact (flint):", ok, f"steps {e} ({time.time()-t0:.0f}s)", flush=True)
    rec = {"L_div_W12_zero_flint": ok}
    # control: the known exact division W12 = Y4 o R8 must also come out zero, and a perturbed W12 must not
    V = json.loads((HERE / "data" / "four_density_path_V.json").read_text())
    R8 = to_delta_poly(V["R8_Q"])
    ok2, _ = rdiv_exact(W, R8, verbose=False)
    Wp = list(W)
    Wp[3] = Wp[3] * flint.fmpq(201, 200)
    ok3, _ = rdiv_exact(L, Wp, verbose=False)
    rec["control_W12_div_R8_zero"] = ok2
    rec["control_L_div_perturbed_W12_zero"] = ok3
    print("control W12 / R8 zero:", ok2, "  L / perturbed W12 zero (must be False):", ok3, flush=True)
    (HERE / "data" / "four_density_path_G2.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
