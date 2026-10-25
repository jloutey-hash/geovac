"""PF1-PF4 stage 6 (copy of debug/ga3_6_survey.py with pf_lib): survey of the minimal annihilator L: local data at t = 0, infinity and every rational root of the leading coefficient; apparent-point check; order-1 right factors of L and L*; minimal annihilators of the log-free Frobenius solutions.
Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).
usage: python pf_6_survey.py <op json> <tag> [rmax] [dmax]"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as Lc
import correlated_triangle_g32_lib4 as S
import cf_g2_lib as C
import pf_lib as A


def main():
    L = C.load_Q(sys.argv[1])
    tag = sys.argv[2]
    rmax = int(sys.argv[3]) if len(sys.argv) > 3 else len(L) - 2
    dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 80
    rec = {}
    t0 = time.time()
    rr = C.rational_roots(L)
    print("rational roots of the leading coefficient (x, multiplicity):", [(str(a), e) for a, e in rr], flush=True)
    xpts = sorted(set([a for a, e in rr if a != 0] + A.XPTS))
    rec["L"] = C.local_report(L, xpts, name="L")
    La = G.adjoint(L)
    rec["Ladj"] = C.local_report(La, xpts, name="L*")
    rec["apparent"] = {}
    for f, e in C.lead_factors(L):
        if f.degree() >= 2:
            ac = C.apparent_check(L, f, nprimes=2)
            rec["apparent"][str(f)[:40]] = ac
            print(f"   apparent check, factor of degree {f.degree()} (mult {e}):", ac, flush=True)
    print(f"({time.time()-t0:.0f}s)", flush=True)
    rec["hyperexp_L"] = S.split_search(L, "sol")
    rec["hyperexp_Ladj"] = S.split_search(La, "sol")
    print(f"({time.time()-t0:.0f}s)", flush=True)
    rec["minann"] = {}
    for name, M in (("L", L), ("L*", La)):
        for x0 in [Fraction(0)] + xpts:
            for cl in Lc.monodromy_type(M, x0):
                for e, m, ok in cl["exponents"]:
                    if not ok:
                        continue
                    res = C.fast_min_annihilator(M, x0, Fraction(e), rmax, dmax)
                    key = f"{name} x0={x0} e={e}"
                    rec["minann"][key] = res
                    print(f"   min annihilator of {key}: {res}", flush=True)
    print(f"({time.time()-t0:.0f}s)", flush=True)
    (C.DATA / f"pf_6_survey_{tag}.json").write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
