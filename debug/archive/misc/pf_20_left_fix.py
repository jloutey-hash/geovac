"""Second gap in the map search (PF addenda 3-5), left end: the adjoint of an operator has poles at the operator's
apparent points (L with exponents 0, 2 at a point: L* has -1, 1), so a map into Sol(F*) may need poles at F's apparent
points.  pf_12_hom allowed poles at the cuts and the SOURCE's own extra points only.  Right-end searches (into F) are
not affected: F's solutions are analytic at its apparent points.
Fix: the denominator also carries the target's extra factors (its leading coefficient minus x and the cuts), to power k.
Discriminating control: the sector filtration (pf_14, pf_15) proves the sunrise block, equivalent to R_V (pf_17), is a
submodule of Y's cyclic module, so R_V is a left factor of L: Hom(R_V*, L*) must be nonzero.  Old denominator vs new.
usage: python pf_20_left_fix.py <out tag> <target: L19|L15> [k M amax]
Third gap found here: the adjoint has t-exponents down to about -25 at x = 0 (L has 3..25), beyond amax = 8."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_lib  # noqa: F401
import pf_12_hom as P

DATA = Path(__file__).parent / "data"
_orig_den = P.den_coeffs
TEXTRA = {"poly": None, "k": 0}


def den_with_target(b, extra=None):
    d = flint.fmpz_poly(_orig_den(b, extra))
    if TEXTRA["poly"] is not None:
        d = d * TEXTRA["poly"] ** TEXTRA["k"]
    return [int(v) for v in d.coeffs()]


def target_extra(F):
    return P.source_extra(F)


def nul(r):
    return {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")}


def load_target(name):
    if name == "L19":
        return G.parse_op(json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())["chain_R8"])
    return G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])


def main():
    tag, tname = sys.argv[1], sys.argv[2]
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    M = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    amax = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    t0 = time.time()
    L = load_target(tname)
    Ls = G.adjoint(L)
    te = target_extra(L)
    print(f"target {tname}: order {len(L)-1}, x-degree {G.op_degree(L)}; extra factor degree {te.degree()}", flush=True)
    RV = __import__("pf_10_rvmap").RV
    out = {"target": tname, "k": k, "M": M, "amax": amax, "target_extra_degree": te.degree()}

    def rec(name, r):
        out[name] = r
        print(f"{name}: {nul(r)} ({time.time()-t0:.0f}s)", flush=True)
        (DATA / f"pf_20_left_fix_{tag}.json").write_text(json.dumps(out, indent=1, default=str))

    P.den_coeffs = _orig_den
    rec("old denominator: Hom(R_V*, L*)", P.hom(G.adjoint(RV), Ls, amax, 6, M))
    P.den_coeffs = den_with_target
    TEXTRA.update(poly=te, k=k)
    rec("new denominator: Hom(R_V*, L*)", P.hom(G.adjoint(RV), Ls, amax, 6, M))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
