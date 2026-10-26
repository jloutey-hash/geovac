"""Which sector blocks share composition factors with the ends of F_8?  And is the density sunrise R_V?

Sources: the six sector blocks of debug/pf_16_blocks_Q.py (theta_t form; their own apparent points go into the map
denominator, pf_12_hom.source_extra).  Hom(S, F) != 0: a left factor of S is equivalent to a right factor of F.
Hom(S*, F*) != 0: a right factor of S is equivalent to a left factor of F.  Box a 8, b 6, numerator degree 300,
two primes, both t-parities (pf_12_hom.hom).
Tests: sunrise block {g1} against R_V both ways; {g1} = {g2} as operators; each block into F_8 and F_8*.
Controls: planted F_8 o S(theta + 1) and S o F_8 for the sunrise and for the {g2, C} block.
usage: python pf_17_block_hom.py <out tag>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_lib  # noqa: F401  (sets the cut points used by the libraries)
import pf_12_hom as P

DATA = Path(__file__).parent / "data"


def nul(r):
    return {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")}


def main():
    tag = sys.argv[1]
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    t0 = time.time()
    pieces = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F8 = G.parse_op(pieces["factor_8"])
    blocks = json.loads((DATA / "pf_16_blocks_Q.json").read_text())
    B = {k: G.parse_op(v["theta_form"]) for k, v in blocks.items()}
    RV = __import__("pf_10_rvmap").RV
    out = {"box": [8, 6, M]}

    def rec(name, r):
        out[name] = r
        print(f"{name}: {nul(r)} ({time.time()-t0:.0f}s)", flush=True)
        (DATA / f"pf_17_block_hom_{tag}.json").write_text(json.dumps(out, indent=1, default=str))

    out["sunrise blocks equal"] = G.op_to_strings(G.content_free(B["(3,)"])) == G.op_to_strings(G.content_free(B["(4,)"]))
    print("sunrise blocks {g1} and {g2} equal as operators:", out["sunrise blocks equal"], flush=True)
    sun = B["(3,)"]
    rec("Hom(sunrise, R_V)", P.hom(sun, RV, 8, 6, M))
    rec("Hom(R_V, sunrise)", P.hom(RV, sun, 8, 6, M))
    for nm in ("(3,)", "(4, 5)"):
        S = B[nm]
        rec(f"planted right {nm}", P.hom(S, G.compose(F8, G.shift_theta(S, 1)), 8, 6, M))
        rec(f"planted left {nm}", P.hom(G.adjoint(S), G.adjoint(G.compose(S, F8)), 8, 6, M))
    for nm, S in B.items():
        if nm == "(4,)":
            continue
        rec(f"Hom({nm}, F_8)", P.hom(S, F8, 8, 6, M))
        rec(f"Hom({nm}*, F_8*)", P.hom(G.adjoint(S), G.adjoint(F8), 8, 6, M))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
