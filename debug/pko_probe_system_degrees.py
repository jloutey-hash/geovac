"""PK record, phase-1 feasibility probe (disclosed in docs/predictions/2026-10-26_parametric_kite_operator.md).
Reads the saved mod-P differential system of the PF kite at (1, 2; 3/2, 3) (debug/data/ib_pf_p5_op.json) and
reports, per entry, the numerator and denominator degrees in z and the roots of the common denominator.
Scores nothing; reads existing data only."""
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from flint import nmod_poly

DATA = Path(__file__).parent / "data"
d = json.loads((DATA / "ib_pf_p5_op.json").read_text())
P = d["P"]


def fr(x):
    x = F(x)
    return x.numerator % P * pow(x.denominator, P - 2, P) % P


cand = {}
for s in ["0", "-1/4", "-1", "-9/4", "-4", "-9", "-49/4", "-16", "-81/4", "-25/4", "-25"]:
    cand[fr(s)] = s
dn, dd = Counter(), Counter()
den_lcm = nmod_poly([1], P)
tot = []
for row in d["A"] + [d["c"]]:
    ents = row if isinstance(row[0][0], list) else [row]
    for nu, de in ents:
        N, D = nmod_poly(nu, P), nmod_poly(de, P)
        if N.is_zero():
            continue
        dn[N.degree()] += 1
        dd[D.degree()] += 1
        tot.append(N.degree() + D.degree())
        den_lcm = den_lcm * D // den_lcm.gcd(D)
print("entries:", len(tot), "max num deg", max(dn), "max den deg", max(dd), "max total", max(tot))
print("num degree histogram", sorted(dn.items()))
print("den degree histogram", sorted(dd.items()))
c, fac = den_lcm.factor()
print("lcm of denominators: degree", den_lcm.degree())
for f, e in fac:
    if f.degree() == 1:
        r = (-int(f.coeffs()[0])) % P
        print("  linear factor, root", cand.get(r, "OTHER"), "multiplicity", e)
    else:
        print("  factor of degree", f.degree(), "multiplicity", e)
