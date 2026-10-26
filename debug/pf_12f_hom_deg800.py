"""PF1, PF3, GR1, TG1 and the induced hypothesis (pf_13_induced.py): the fixed map search (pf_12_hom.py, source
apparent points in the denominator) at numerator degree 800 instead of 300 (R_V 500), for every source into F_8 and
F_8*.  F_8 has x-degree 282, mostly apparent points, and a map onto a factor of it could need a large numerator.
Control at this degree: the R_V and z1 planted right factor F_8 o S(theta + 1) and planted left factor S o F_8.
usage: python pf_12f_hom_deg800.py <out tag> [M]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_12_hom as P

DATA = Path(__file__).parent / "data"


def main():
    tag = sys.argv[1]
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 800
    pieces = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F8 = G.parse_op(pieces["factor_8"])
    ind = json.loads((DATA / "pf_13_induced_v1.json").read_text())
    src = {"R_V": __import__("pf_10_rvmap").RV}
    for nm in ("z1", "z2", "d1", "d2"):
        src[nm] = P.legendre_source(nm)
    for tw, rec in ind.items():
        src[f"induced E+ twist {tw}"] = G.parse_op(rec["theta_form"])
    out = {"M": M}
    t0 = time.time()
    for nm in ("R_V", "z1"):
        S = src[nm]
        r = P.hom(S, G.compose(F8, G.shift_theta(S, 1)), 8, 6, M)
        l = P.hom(G.adjoint(S), G.adjoint(G.compose(S, F8)), 8, 6, M)
        out[f"planted {nm}"] = {"right": r, "left": l}
        print("planted", nm, {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")},
              {k: v["nullity_two_primes"] for k, v in l.items() if k.startswith("nu")}, f"({time.time()-t0:.0f}s)", flush=True)
    for nm, S in src.items():
        r = P.hom(S, F8, 8, 6, M)
        l = P.hom(G.adjoint(S), G.adjoint(F8), 8, 6, M)
        out[f"Hom({nm}, F_8)"] = r
        out[f"Hom({nm}*, F_8*)"] = l
        print(nm, {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")},
              {k: v["nullity_two_primes"] for k, v in l.items() if k.startswith("nu")}, f"({time.time()-t0:.0f}s)", flush=True)
        (DATA / f"pf_12f_hom_deg800_{tag}.json").write_text(json.dumps(out, indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
