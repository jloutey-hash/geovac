"""Mechanism probe for the JOB-A floor: partial-wave (l_max) increments at fixed n_max.

If the K->inf floor at lmax=3 is the partial-wave truncation error, the l-increments
DeltaE_l = E(lmax=l) - E(lmax=l-1) should follow the Schwartz law ~ (l+1/2)^-4 and the
extrapolated tail sum_{l>3} should reproduce the fitted floor.
"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import debug.p60_engine as E
import geovac.sturmian_secular as S

GROUND, EXC1, CHEM = -2.903724377, -2.145974046, 1.5936014616
out = {}
for nmax in [int(x) for x in (sys.argv[1:] or ["10", "12"])]:
    E.set_grid(max(80.0, 5.0 * nmax * nmax), 24000, "grade", 2.0)
    prev0 = prev1 = None
    rows = []
    for lmax in range(0, 4):
        t0 = time.time()
        cfgs = S.build_configs(E.family(nmax, lmax))
        M = S.build_M(cfgs, Z=2.0)
        p = np.sort(np.linalg.eigvalsh(M))[::-1]
        e0, e1 = -p[0] ** 2 / 2, -p[1] ** 2 / 2
        d0 = None if prev0 is None else (e0 - prev0) * 1000
        d1 = None if prev1 is None else (e1 - prev1) * 1000
        rows.append(dict(lmax=lmax, K=len(cfgs), e0=e0, e1=e1,
                         gap0=(e0 - GROUND) * 1000, gap1=(e1 - EXC1) * 1000,
                         d0_mHa=d0, d1_mHa=d1, wall=time.time() - t0))
        print("nmax=%2d lmax=%d K=%4d | E0=%.9f gap0=%8.4f mHa dE0=%s | E1=%.9f gap1=%7.4f mHa dE1=%s [%.0fs]"
              % (nmax, lmax, len(cfgs), e0, (e0 - GROUND) * 1000,
                 "   --  " if d0 is None else "%+8.4f" % d0, e1, (e1 - EXC1) * 1000,
                 "   --  " if d1 is None else "%+7.4f" % d1, time.time() - t0))
        sys.stdout.flush()
        prev0, prev1 = e0, e1
    # Schwartz (l+1/2)^-4 tail from the last increment
    for tag, key in (("gnd", "d0_mHa"), ("2^1S", "d1_mHa")):
        d3 = rows[3][key]
        tail = sum(d3 * ((3 + 0.5) / (l + 0.5)) ** 4 for l in range(4, 4000))
        print("   %s: DeltaE_3 = %+.4f mHa -> Schwartz (l+1/2)^-4 tail beyond lmax=3 = %+.4f mHa (%.2fx chem)"
              % (tag, d3, tail, abs(tail) / CHEM))
        rows[3][tag + "_schwartz_tail_mHa"] = tail
    out[str(nmax)] = rows
json.dump(out, open("debug/data/p60_lmax_increments.json", "w"), indent=1)
print("wrote debug/data/p60_lmax_increments.json")
