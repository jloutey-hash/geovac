"""PF1, PF3, GR1, TG1: the fixed map search (pf_12_hom.py with source apparent points in the denominator) at twelve
primes instead of two.  A factor equivalent to a source only over a quadratic field Q(sqrt d) gives maps mod p exactly
at the primes where d is a square, so two primes can both miss it.  Box a 8, b 6, numerator degree 300 (R_V 300).
usage: python pf_12e_hom_primes.py <out tag>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_12_hom as P

DATA = Path(__file__).parent / "data"


def main():
    tag = sys.argv[1]
    ps = []
    q = 2 ** 30
    while len(ps) < 12:
        q = int(sp.prevprime(q))
        ps.append(q)
        q -= 10007
    P.PRIMES = ps
    pieces = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F8 = G.parse_op(pieces["factor_8"])
    src = {"R_V": __import__("pf_10_rvmap").RV}
    for nm in ("z1", "z2", "d1", "d2"):
        src[nm] = P.legendre_source(nm)
    out = {"primes": ps, "box": [8, 6, 300]}
    t0 = time.time()
    for nm, S in src.items():
        r = P.hom(S, F8, 8, 6, 300)
        l = P.hom(G.adjoint(S), G.adjoint(F8), 8, 6, 300)
        out[f"Hom({nm}, F_8)"] = r
        out[f"Hom({nm}*, F_8*)"] = l
        print(nm, "right", {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")},
              "left", {k: v["nullity_two_primes"] for k, v in l.items() if k.startswith("nu")}, f"({time.time()-t0:.0f}s)", flush=True)
        (DATA / f"pf_12e_hom_primes_{tag}.json").write_text(json.dumps(out, indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
