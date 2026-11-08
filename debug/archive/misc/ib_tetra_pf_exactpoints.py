"""PF2 Y clause on the exact off-slice tetrahedron operator (data/ib_tetra_pf_op_Q.json, form_D on S = Y/t^3 in
x = t^2): factor the leading coefficient over Z, exact local exponents over Q at each rational root, and write the
operator reduced mod a fresh prime (not one it was built from) in the format of ib_tetra_pf_apparent.py.
usage: python ib_tetra_pf_exactpoints.py"""
import json
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
from flint import fmpz_poly, fmpq_poly, fmpq

DATA = Path(__file__).parent / "data"
d = json.loads((DATA / "ib_tetra_pf_op_Q.json").read_text())
ops = [fmpz_poly([int(v) for v in p]) for p in d["polys"]]
order = len(ops) - 1
lead = ops[-1]
c, facs = lead.factor()
print("order", order, "leading degree", lead.degree())
lin = []
for f, e in facs:
    print(" factor degree", f.degree(), "mult", e, "coeffs" if f.degree() <= 2 else "", f.coeffs() if f.degree() <= 2 else "")
    if f.degree() == 1:
        a0, a1 = [int(v) for v in f.coeffs()]
        lin.append(F(-a0, a1))
print("rational roots:", sorted(lin))
x = sp.Symbol("x"); th = sp.Symbol("th")


def exponents(rho):
    """indicial polynomial at x = rho (rho = 0 handled the same way), exact over Q."""
    sh = fmpq_poly([fmpq(rho.numerator, rho.denominator), 1])
    Ps = [fmpq_poly(q)(sh) for q in ops]
    v = [min((i for i, cc in enumerate(q.coeffs()) if cc != 0), default=10**9) for q in Ps]
    delta = min(v[k] - k for k in range(order + 1))
    ind = 0
    for k in range(order + 1):
        if v[k] - k == delta:
            ff = 1
            for j in range(k):
                ff *= (th - j)
            cc = Ps[k].coeffs()[v[k]]
            ind += sp.Rational(int(cc.p), int(cc.q)) * ff
    return sp.roots(sp.Poly(sp.expand(ind), th))


for rho in [F(0)] + sorted(lin) + [F(25, 4), F(25), F(3)]:
    r = exponents(rho)
    print(f"x = {rho}: exponents", sorted(((sp.nsimplify(k), m) for k, m in r.items()), key=lambda z: float(z[0])),
          "total", sum(r.values()), flush=True)

# reduce mod a fresh prime for the extension-field apparent test
used = set(int(json.loads(Path(f).read_text())["P"]) for f in
           [DATA / "ib_tetra_pf_op.json", DATA / "ib_tetra_pf_p1_op.json", DATA / "ib_tetra_pf_p2_op.json"])
q = int(sp.prevprime(2**61 - 12345))
while q in used:
    q = int(sp.prevprime(q))
out = {"P": q, "polys": [[int(v) % q for v in p.coeffs()] for p in ops],
       "note": "exact operator ib_tetra_pf_op_Q.json reduced mod a fresh prime"}
(DATA / "ib_tetra_pf_exact_modq.json").write_text(json.dumps(out))
print("fresh prime", q, "wrote ib_tetra_pf_exact_modq.json")
