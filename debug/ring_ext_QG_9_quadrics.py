"""QG step 9: quadric counter (method of debug/linking_ceiling_C_quadrics.py) on the exact factors:
constant-coefficient quadratic relations among a Taylor basis of solutions at rational regular points.
Reading (Zariski-dense monodromy): SL2 order 2 -> 0; Sym^2 (K3 of F_3) -> 1; Sp4 / SL4 -> 0.
Controls: F_3's operator (expect 1), Legendre (0).
usage: python ring_ext_QG_9_quadrics.py <gamma>"""
from __future__ import annotations
import json, sys
from fractions import Fraction as Fr
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_ops as O
from ring_ext_QG_ops import u, tt as t

DATA = Path(__file__).parent / "data"
h = sp.Symbol("h")


def d_form_from_theta(Lc):
    """function-left theta coefficients -> D-form coefficients a_k(t) (polynomial after clearing)."""
    from sympy.functions.combinatorial.numbers import stirling
    r = len(Lc) - 1
    a = [sp.Integer(0)] * (r + 1)
    for j, c in enumerate(Lc):
        for k in range(j + 1):
            a[k] += c * stirling(j, k, kind=2) * t**k
    den = sp.lcm([sp.denom(sp.together(x)) for x in a])
    return [sp.expand(sp.cancel(x * den)) for x in a]


def taylor_basis(a, t0, N):
    r = len(a) - 1
    A = [[Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in
          reversed(sp.Poly(sp.expand(x.subs(t, t0 + h)), h).all_coeffs())] for x in a]
    assert A[r][0] != 0, "t0 singular"
    ff = lambda n, k: 1 if k == 0 else n * ff(n - 1, k - 1)
    basis = []
    for e in range(r):
        c = [Fr(0)] * (N + r)
        for k in range(r):
            c[k] = Fr(1 if k == e else 0) / ff(k, k)
        for M in range(N):
            tot = Fr(0)
            for k in range(r + 1):
                for l, alk in enumerate(A[k]):
                    idx = M - l + k
                    if alk == 0 or idx < 0 or (k == r and l == 0):
                        continue
                    tot += alk * c[idx] * ff(idx, k)
            c[M + r] = -tot / (A[r][0] * ff(M + r, r))
        basis.append(c[:N])
    return basis


def rank_Q(rows):
    rows = [r[:] for r in rows]
    rk = 0
    for col in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col] != 0), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][col] != 0:
                f = rows[i][col] / rows[rk][col]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def quadrics(Lc, t0, N=40):
    a = d_form_from_theta(Lc)
    r = len(a) - 1
    B = taylor_basis(a, t0, N)
    prods = [[sum(B[i][k] * B[j][M - k] for k in range(M + 1)) for M in range(N)] for i in range(r) for j in range(i, r)]
    return r, len(prods), len(prods) - rank_Q(prods)


def main():
    g = sys.argv[1]
    rec = {}
    ops = {"control F_3 op": O.Q_to_coeffs([4 * u**3, -(u + 1) ** 3]),
           "control Legendre": O.Q_to_coeffs([u**2, -(u + 1) ** 2])}
    k3 = json.loads((DATA / f"ring_ext_QG_7_k3map_g{g}.json").read_text())
    ops["R3 (K3 factor of G4)"] = [sp.sympify(x, locals={"t": t}) for x in k3["R3_monic"]]
    d = json.loads((DATA / f"ring_ext_QG_6_gcrd_g{g}.json").read_text())
    ops["G4 (minimal op of Y)"] = [sp.sympify(x, locals={"t": t}) for x in d["G4"]]
    for nm, Lc in ops.items():
        for t0 in (sp.Rational(1, 2), sp.Rational(7, 5)):
            r, pairs, q = quadrics(Lc, t0)
            print(f"{nm:24s} order {r} t0={t0}: pairs {pairs}, quadratic relations {q}", flush=True)
            rec[f"{nm}@{t0}"] = {"order": r, "pairs": pairs, "relations": q}
    (DATA / f"ring_ext_QG_9_quadrics_g{g}.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
