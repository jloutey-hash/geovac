"""Is R_V the only irreducible left factor of L15?  If so the socle of D/D L15 is R_V, every proper quotient kills R_V,
and since R_S (the right factor without R_V) leaves Y nonzero (pf_23_needs), L15 is Y's minimal operator.
Candidates for another simple submodule are the other composition factors: z1, z2 (exact maps into the middle
piece, pf_24) and F_1..F_7 (pf_pieces_tetra_exact).  Test Hom(S*, L15*) for each, box a 20 (L15* has t-exponents
down to -31 at 0), b 12, the source's extra points to b, L15's extra factor (degree 77) to power 2, M 400.
Positive controls in the same box: R_V (expected 1) and the planted F_3 o L15 for an order-one source.
usage: python pf_26_socle.py <tag>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_lib  # noqa: F401
import pf_12_hom as P
import pf_20_left_fix as LF

DATA = Path(__file__).parent / "data"
AMAX, BMAX, M = 20, 12, 400


def left(S, F):
    P.den_coeffs = LF.den_with_target
    LF.TEXTRA.update(poly=P.source_extra(F), k=2)
    r = P.hom(G.adjoint(S), G.adjoint(F), AMAX, BMAX, M)
    P.den_coeffs = LF._orig_den
    return r


def nul(r):
    return {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")}


def main():
    tag = sys.argv[1]
    t0 = time.time()
    L15 = G.parse_op(json.loads((DATA / "pf_19_quotient_Q.json").read_text())["L15_lam_1/2"])
    pc = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    Fs = {f"F_{i}": G.parse_op(pc[f"factor_{i}"]) for i in range(1, 8)}
    src = {"R_V": __import__("pf_10_rvmap").RV, "z1": P.legendre_source("z1"), "z2": P.legendre_source("z2"), **Fs}
    out = {"box": [AMAX, BMAX, M]}

    def rec(name, r):
        out[name] = r
        print(f"{name}: {nul(r)} ({time.time()-t0:.0f}s)", flush=True)
        (DATA / f"pf_26_socle_{tag}.json").write_text(json.dumps(out, indent=1, default=str))

    rec("planted left: Hom(F_3*, (F_3 o L15)*)", left(Fs["F_3"], G.compose(Fs["F_3"], L15)))
    for nm, S in src.items():
        rec(f"Hom({nm}*, L15*)", left(S, L15))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
