"""GS2's cap: is R_V (the K3 operator of the chain (t; 1, 2), equivalent to the sunrise sector block) split by an
order-two factor?  The Frobenius-solution probe of pf_lib (minimal annihilators of every log-free Frobenius solution
of R_V and of R_V*, at 0 and the cut points, order < 4, coefficient degree <= dmax).  Control: the same probe on the
planted product z1 o R_V-free... here on z1 o z2 (pf_12_hom Legendre sources), which must split.
usage: python pf_25_rv_split.py <dmax>"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_lib as PL
import pf_12_hom as P

DATA = Path(__file__).parent / "data"


def main():
    dmax = int(sys.argv[1])
    t0 = time.time()
    RV = __import__("pf_10_rvmap").RV
    ctrl = G.compose(P.legendre_source("z1"), P.legendre_source("z2"))
    out = {"dmax": dmax}
    for nm, L in (("control z1 o z2", ctrl), ("R_V", RV)):
        c = PL.probe(L, dmax, verbose=False)
        out[nm] = [{k: str(v) for k, v in d.items()} for d in c if d["r"] < PL.order(L)]
        print(f"{nm}: proper factor candidates {out[nm]} ({time.time()-t0:.0f}s)", flush=True)
    (DATA / f"pf_25_rv_split_d{dmax}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
