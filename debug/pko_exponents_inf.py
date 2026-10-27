"""PK record: exponents at z = infinity of an integer d/dz-form operator sum_k p_k(z) D^k.
Solutions behave like z^rho at infinity, rho a root of sum_k lc_k ff(rho, k) over the k with
deg p_k - k maximal.  usage: python pko_exponents_inf.py <op json>.  On the base point this only reads
existing data (disclosure in docs/predictions/2026-10-26_parametric_kite_operator.md)."""
import json
import sys
import sympy as sp

d = json.load(open(sys.argv[1]))
polys = d["polys"]
r = sp.Symbol("r")
top = max(len(p) - 1 - k for k, p in enumerate(polys))
ind = 0
for k, p in enumerate(polys):
    if len(p) - 1 - k == top:
        ind += sp.Integer(int(p[-1])) * sp.ff(r, k)
roots = sp.roots(sp.Poly(sp.expand(ind), r))
print("rho at infinity (solutions ~ z^rho):", dict(sorted(roots.items())), "found", sum(roots.values()), "of", len(polys) - 1)
