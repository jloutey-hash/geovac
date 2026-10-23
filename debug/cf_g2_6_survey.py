"""CF5-CF7 stage 6: survey of the minimal annihilator L of Y_log at gamma = 2.
 (a) local data at t = 0, infinity and every rational root of the leading coefficient; apparent-point
     check (mod p) of the irreducible factors;
 (b) order-1 right factors of L and of L* (hyperexponential solutions; exponent-bounded search, lib4);
 (c) minimal annihilators (mod p) of the log-free Frobenius solutions of L and of L* at the true
     singular points (right factors of L, and of L*, i.e. left factors of L).
usage: python cf_g2_6_survey.py <op json> <tag> [rmax] [dmax] [t points, default 1,2,3,4]"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as Lc
import correlated_triangle_g32_lib4 as S
import cf_g2_lib as C

TPOS = [Fraction(1), Fraction(2), Fraction(3), Fraction(4)]


def main():
    L = C.load_Q(sys.argv[1])
    tag = sys.argv[2]
    rmax = int(sys.argv[3]) if len(sys.argv) > 3 else len(L) - 2
    dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 70
    tpos = [Fraction(v) for v in sys.argv[5].split(",")] if len(sys.argv) > 5 else TPOS
    C.set_points(tpos)
    rec = {}
    t0 = time.time()
    rr = C.rational_roots(L)
    print("rational roots of the leading coefficient (x, multiplicity):", [(str(a), e) for a, e in rr], flush=True)
    xpts = sorted(set([a for a, e in rr if a != 0] + [v * v for v in tpos]))
    rec["L"] = C.local_report(L, xpts, name="L")
    La = G.adjoint(L)
    rec["Ladj"] = C.local_report(La, xpts, name="L*")
    rec["apparent"] = {}
    for f, e in C.lead_factors(L):
        if f.degree() >= 2:
            ac = C.apparent_check(L, f, nprimes=2)
            rec["apparent"][str(f)[:40]] = ac
            print(f"   apparent check, factor of degree {f.degree()}:", ac, flush=True)
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
    (C.DATA / f"cf_g2_6_survey_{tag}.json").write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
