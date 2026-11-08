"""GS step B (record docs/predictions/2026-10-21_generic_slice.md, commit a606d28f): exact peeling of
the elementary right factors of L (step A), b = 2.  theta-form coefficient lists, low -> high.
Candidates: theta + 2 (solution t^-2); R1a with solution t/(1-t^2)^4 (the b = 1 slice's R1);
R1b with solution t/(4-t^2)^4 (its rescaling to the leaf at scale 2)."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp

t, u = sp.symbols("t u")
HERE = Path(__file__).parent
A = json.loads((HERE / "data" / "generic_slice_A.json").read_text())
Q = [sp.sympify(q, locals={"u": u}) for q in A["Q"]]


def coeff_list(Q):
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    out = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            out[k] += c * t ** (2 * i)
    return [sp.expand(c) for c in out]


def cth(c):
    return sp.cancel(t * sp.diff(c, t))


def compose(A, B):
    out = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        if a == 0:
            continue
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = cth(bi)
    return [sp.cancel(c) for c in out]


def rdiv(L, R):
    L = [sp.cancel(c) for c in L]
    r = len(R) - 1
    Qo = [sp.Integer(0)] * (len(L) - r)
    for k in range(len(L) - 1, r - 1, -1):
        c = sp.cancel(L[k] / R[-1])
        Qo[k - r] = c
        if c == 0:
            continue
        sub = compose([0] * (k - r) + [c], R)
        L = [sp.cancel(L[i] - (sub[i] if i < len(sub) else 0)) for i in range(len(L))]
    return Qo, L[:r]


def first_order(h):
    """theta - theta(h)/h, scaled to polynomial coefficients."""
    lam = sp.cancel(cth(h) / h)
    n, d = sp.fraction(lam)
    return [sp.expand(-n), sp.expand(d)]


def main():
    rec = {}
    L = coeff_list(Q)
    print("order", len(L) - 1, flush=True)
    cands = {"th+2": [sp.Integer(2), sp.Integer(1)],
             "R1a": first_order(t / (1 - t**2) ** 4),
             "R1b": first_order(t / (4 - t**2) ** 4)}
    print({k: v for k, v in cands.items()}, flush=True)
    cur, peeled = L, []
    for name in ["th+2", "R1a", "R1b"]:
        t0 = time.time()
        qo, rem = rdiv(cur, cands[name])
        ok = all(sp.cancel(r) == 0 for r in rem)
        print(f"peel {name}: remainder zero = {ok}  ({time.time()-t0:.0f}s)", flush=True)
        rec[name] = ok
        if ok:
            cur, peeled = qo, peeled + [name]
        else:
            # try the other order later
            pass
    if not rec.get("R1b", False) or not rec.get("R1a", False):
        # second attempt: R1b before R1a
        L1, _ = rdiv(L, cands["th+2"])
        for order in (["R1b", "R1a"],):
            c2 = L1
            res = []
            for name in order:
                qo, rem = rdiv(c2, cands[name])
                ok = all(sp.cancel(r) == 0 for r in rem)
                res.append((name, ok))
                if ok:
                    c2 = qo
            print("alternative order", res, flush=True)
            rec["alt"] = res
    rec["peeled"] = peeled
    rec["quotient_order"] = len(cur) - 1
    rec["quotient"] = [str(sp.factor(c)) for c in cur]
    (HERE / "data" / "generic_slice_B.json").write_text(json.dumps(rec, indent=1))
    print("quotient order", len(cur) - 1, flush=True)
    print("quotient leading coefficient:", sp.factor(cur[-1]), flush=True)


if __name__ == "__main__":
    main()
