"""PF1, PF3, GR1: the planted controls of pf_12_hom.py at the large box and with targets as large as F_8: the source S
as a right factor of F_8 o S(theta + 1) and as a left factor of S o F_8 (adjoint search).  Both must be found.
usage: python pf_12b_hom_planted_big.py <out tag> [amax bmax M]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_12_hom as P

DATA = Path(__file__).parent / "data"


def main():
    tag = sys.argv[1]
    amax, bmax, M = (int(v) for v in sys.argv[2:5]) if len(sys.argv) > 4 else (8, 6, 300)
    t0 = time.time()
    pieces = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F8 = G.parse_op(pieces["factor_8"])
    src = {"R_V": __import__("pf_10_rvmap").RV, "z1": P.legendre_source("z1"), "z2": P.legendre_source("z2")}
    out = {"box": {"amax": amax, "bmax": bmax, "M": M}}
    for nm, S in src.items():
        r = P.hom(S, G.compose(F8, G.shift_theta(S, 1)), amax, bmax, M)
        l = P.hom(G.adjoint(S), G.adjoint(G.compose(S, F8)), amax, bmax, M)
        out[f"planted right: Hom({nm}, F_8 o {nm}(theta+1))"] = r
        out[f"planted left: Hom({nm}*, ({nm} o F_8)*)"] = l
        print(nm, "right", {k: v["nullity_two_primes"] for k, v in r.items()}, "left", {k: v["nullity_two_primes"] for k, v in l.items()},
              f"({time.time()-t0:.0f}s)", flush=True)
    (DATA / f"pf_12b_hom_planted_big_{tag}.json").write_text(json.dumps(out, indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
