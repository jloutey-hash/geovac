"""From the derived closed form of K (pachucki_bound_1s.py): the two coordinate relations,
the minimal operator over Q(R), its singular data, and comparison with the committed
mod-p operator (debug/data/exch_run_1_2_<p>_EXS.json; read-only).

Run: python debug/pachucki_bound_operator.py [A B [prime:suffix ...]]   (default 1 2)
"""
from __future__ import annotations

import json
import sys

import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, "debug")
import pachucki_bound_core as C  # noqa: E402
from pachucki_bound_core import r, LC  # noqa: E402

# rebuild K exactly (same code path as pachucki_bound_1s, imported as a module run)
import runpy  # noqa: E402
PAIR = sys.argv[1:3] if len(sys.argv) > 2 else ["1", "2"]
PRIMES = [(int(x.split(":")[0]), x.split(":")[1]) for x in sys.argv[3:]] or [(67108859, "EXS"), (33554393, "EXS")]
sys.argv = ["x"] + PAIR
ns = runpy.run_path("debug/pachucki_bound_1s.py")
Kfun, closure = ns["Kfun"], ns["closure"]

rows, cur = [], Kfun
for i in range(10):
    rows.append([sum(c * r**p for p, c in cur.get(k, {}).items()) for k in closure])
    cur = C.fdiff(cur)
M = sp.Matrix(rows)                          # 10 x 11
dom = sp.QQ.frac_field(r, LC)

# (a) right null space of the first 9 rows' span: the 2 coordinate relations
Md = DomainMatrix.from_Matrix(M).convert_to(dom)
ns_right = Md.nullspace().to_Matrix()
print("coordinate relations (rows = vectors over keys", [str(k) for k in closure], "):")
for v in ns_right.tolist():
    print("  ", [sp.factor(x) for x in v])

# (b) minimal operator: left null vector of the 10 x 11 matrix
lv = DomainMatrix.from_Matrix(M.T).convert_to(dom).nullspace().to_Matrix()
assert lv.shape[0] == 1, lv.shape
c = [sp.factor(x) for x in lv.row(0)]
den = sp.lcm([sp.fraction(sp.together(x))[1] for x in c])
c = [sp.factor(sp.together(x * den)) for x in c]
import functools
g = functools.reduce(sp.gcd, [x for x in c if x != 0])
c = [sp.factor(sp.cancel(x / g)) for x in c]
print("operator depends on LC:", any(x.has(LC) for x in c))
polys = [sp.Poly(sp.expand(x), r) for x in c]
print("order", len(c) - 1, " coefficient degrees", [p.degree() for p in polys])

# leading coefficient factorisation, R = 0 indicial exponents
lc = sp.factor(polys[-1].as_expr())
print("leading coefficient:", lc)
k = sp.Symbol("k")
# indicial polynomial at R=0: sum_i c_i(R) d^i, lowest weight val(c_i) - i
vals = [(p.as_expr().as_poly(r).monoms()[-1][0] if not p.is_zero else None) for p in polys]
wmin = min(v - i for i, v in enumerate(vals) if v is not None)
ind = sum(polys[i].as_expr().coeff(r, wmin + i) * sp.ff(k, i) for i in range(len(polys)))
print("indicial roots at R=0:", sp.roots(sp.Poly(sp.expand(ind), k)))
# rates at infinity: top-degree part
dtop = max(p.degree() for p in polys)
char = sum(polys[i].as_expr().coeff(r, dtop) * sp.Symbol("x")**i for i in range(len(polys)))
print("characteristic poly at infinity:", sp.factor(char))

# (c) compare with the committed mod-p operator
out = {"coordinate_relations": [[str(sp.factor(x)) for x in v] for v in ns_right.tolist()],
       "order": len(c) - 1, "degrees": [p.degree() for p in polys],
       "indicial_roots": {str(kk): vv for kk, vv in sp.roots(sp.Poly(sp.expand(ind), k)).items()},
       "char_inf": str(sp.factor(char))}
TAG = "_".join(x.replace("/", "_") for x in PAIR)
for prime, sfx in PRIMES:
    J = json.load(open(f"debug/data/exch_run_{TAG}_{prime}_{sfx}.json"))
    cm = J["minimal_operator"]["coeffs_mod_p"]
    ours = [[int(sp.Rational(co).p * pow(sp.Rational(co).q, -1, prime)) % prime
             for co in reversed(p.all_coeffs())] for p in polys]
    # scale: match leading coefficient of order-9 term's top power
    s = cm[-1][-1] * pow(ours[-1][-1], -1, prime) % prime
    ours_s = [[x * s % prime for x in row] for row in ours]
    same = all(len(ra) == len(rb) and ra == rb for ra, rb in zip(ours_s, cm))
    print("prime", prime, ": derived operator == committed operator (up to scalar):", same)
    out[f"equal_mod_{prime}"] = same
json.dump(out, open(f"debug/data/pachucki_bound_operator_{TAG}.json", "w"), indent=1)
