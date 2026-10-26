"""PF1, PF3, GR1, TG1: re-run of pf_12_hom.py with the denominator fixed to allow poles at the source's own singular and
apparent points (the first run allowed the eight cuts only; reading TG1 showed the gap).  New discriminating control:
the source planted with a pole at its apparent point, F_1 o (S o (x - x_a)), whose map y -> y/(x - x_a) the old
denominator cannot express (must give 0) and the new one must find.  R_V (no extra points) is re-run at a larger
numerator degree.
usage: python pf_12d_hom_rerun.py <out tag> [amax bmax M MRV]"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import pf_12_hom as P

DATA = Path(__file__).parent / "data"


def main():
    tag = sys.argv[1]
    amax, bmax, M, MRV = (int(v) for v in sys.argv[2:6]) if len(sys.argv) > 5 else (8, 6, 300, 500)
    t0 = time.time()
    pieces = json.loads((DATA / "pf_pieces_tetra_exact.json").read_text())
    F8 = G.parse_op(pieces["factor_8"])
    F1 = G.parse_op(pieces["factor_1"])
    L = G.parse_op(pieces["chain_R8"])
    out = {"box": {"amax": amax, "bmax": bmax, "M": M, "M_RV": MRV, "den": "eight cuts times the source's other lead factors"}}
    ctrl, tests = {}, {}
    for nm in ("z1", "z2", "d1", "d2"):
        S = P.legendre_source(nm)
        ex = P.source_extra(S)
        lin = [f for f, _ in ex.factor()[1] if f.degree() == 1]
        if lin:
            planted = G.compose(F1, G.compose(S, [lin[0]]))
            old = P.hom(S, planted, 2, 1, 20, extra=None)
            new = P.hom(S, planted, 2, 1, 20)
            ctrl[f"apparent-pole planted, Hom({nm}, F_1 o {nm} o ({lin[0]})): old den"] = old
            ctrl[f"apparent-pole planted, Hom({nm}, F_1 o {nm} o ({lin[0]})): new den"] = new
            print(nm, "apparent-pole planted control: old den", {k: v["nullity_two_primes"] for k, v in old.items() if k.startswith("nu")},
                  "new den", {k: v["nullity_two_primes"] for k, v in new.items() if k.startswith("nu")}, flush=True)
        for tname, T in (("F_8", F8), ("L", L)):
            r = P.hom(S, T, amax, bmax, M)
            l = P.hom(G.adjoint(S), G.adjoint(T), amax, bmax, M)
            tests[f"Hom({nm}, {tname})"] = r
            tests[f"Hom({nm}*, {tname}*)"] = l
            print(nm, tname, "right", {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")},
                  "left", {k: v["nullity_two_primes"] for k, v in l.items() if k.startswith("nu")}, f"({time.time()-t0:.0f}s)", flush=True)
        (DATA / f"pf_12d_hom_rerun_{tag}.json").write_text(json.dumps({**out, "controls": ctrl, "tests": tests}, indent=1))
    RV = __import__("pf_10_rvmap").RV
    for tname, T in (("F_8", F8), ("L", L)):
        r = P.hom(RV, T, amax, bmax, MRV)
        l = P.hom(G.adjoint(RV), G.adjoint(T), amax, bmax, MRV)
        tests[f"Hom(R_V, {tname}) M={MRV}"] = r
        tests[f"Hom(R_V*, {tname}*) M={MRV}"] = l
        print("R_V", tname, f"M={MRV}", "right", {k: v["nullity_two_primes"] for k, v in r.items() if k.startswith("nu")},
              "left", {k: v["nullity_two_primes"] for k, v in l.items() if k.startswith("nu")}, f"({time.time()-t0:.0f}s)", flush=True)
    out.update({"controls": ctrl, "tests": tests, "seconds": time.time() - t0})
    (DATA / f"pf_12d_hom_rerun_{tag}.json").write_text(json.dumps(out, indent=1))
    print("saved", flush=True)


if __name__ == "__main__":
    main()
