"""EC1-EC3 positive control (docs/predictions/2026-10-22_c2b_identification.md): Legendre's family.
(1) Stage C with three poles (0, 1, infinity): every splitting of 6, must recover
    j = 256 (l^2 - l + 1)^3 / (l^2 (l - 1)^2).
(2) Stage E: from the j-line hypergeometric operator pulled back by that j and gauged by P^(-1/4),
    an exact order-one map onto the gauged Legendre operator of c2b_ident_ops (A = 1 + (l^2 + 3) d/dl
    applied to K(l)), found by the same linear algebra as for C2b."""
import itertools
import json
import sys
from pathlib import Path

import flint
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import c2b_ident_C_jsolve as C  # noqa: E402
import c2b_ident_ops as O  # noqa: E402
from c2b_ident_D_compare import pullback_op  # noqa: E402
from c2b_ident_E_map import fq, gauge, apply_c2b_to_image, int_nullspace  # noqa: E402

x = sp.Symbol("x")
DATA = Path(__file__).parent / "data"


def stage_c():
    found = []
    for n, a in itertools.product(range(1, 5), repeat=2):
        c = 6 - n - a
        if c < 1:
            continue
        D = sp.expand(x**n * (x - 1) ** a)
        sols, note, degs = C.solve_D(D)
        ex = [C.exact_check(D, s["identified"]) for s in sols if s["real"]]
        ex = [e for e in ex if e]
        print(f"Legendre widths ({n},{a},{c}): {len(sols)} solutions, exact: {[e['J'] for e in ex]}", flush=True)
        found += [(n, a, c, e) for e in ex]
    target = sp.factor(256 * (x**2 - x + 1) ** 3 / (x**2 * (x - 1) ** 2))
    ok = any(sp.simplify(sp.sympify(e["J"]) - target) == 0 for *_, e in found)
    print("recovered Legendre j:", ok)
    return ok


def stage_e(k=0, dN=4):
    P = x**2 - x + 1
    J = 256 * P**3 / (x**2 * (x - 1) ** 2)
    ops = pullback_op(J)
    src = gauge(ops, sp.together(-sp.Rational(1, 4) * sp.diff(P, x) / P))
    g = sp.gcd_list([sp.expand(v) for v in src])
    src = [sp.expand(sp.cancel(v / g)) for v in src]
    b = [fq(v) for v in src]
    at, _ = O.legendre_gauged()
    a = [flint.fmpq_poly([flint.fmpq(int(c)) for c in aj]) for aj in at]
    E = fq(sp.Integer(1) if k == 0 else (x * (x - 1)) ** k)
    X = flint.fmpq_poly([0, 1])
    zero = flint.fmpq_poly([])
    cols = [apply_c2b_to_image(a, b, X**j, zero, E) for j in range(dN + 1)]
    cols += [apply_c2b_to_image(a, b, zero, X**j, E) for j in range(dN + 2)]
    deg = max(max(u.degree(), v.degree()) for u, v in cols) + 1
    rows = []
    for r in range(deg):
        rows.append([u.coeffs()[r] if r < len(u.coeffs()) else 0 for u, v in cols])
        rows.append([v.coeffs()[r] if r < len(v.coeffs()) else 0 for u, v in cols])
    NS, nul = int_nullspace(rows, len(cols))
    print(f"Legendre stage E (k = {k}, dN = {dN}): nullity {nul}")
    return nul


if __name__ == "__main__":
    okc = stage_c()
    nul = stage_e()
    json.dump({"stage_c_recovered": okc, "stage_e_nullity": nul}, open(DATA / "c2b_ident_F_controls.json", "w"))
