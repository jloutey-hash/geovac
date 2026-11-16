"""Building-block check (zeta-free): generator X_{j1j2}(p) series == sympy
ordered_xi_general series, tau <= TMAX, exact over Q."""
import sys, time
sys.path.insert(0, "debug"); sys.path.insert(0, "."); sys.path.insert(0, "tests")
from exch_gen import *
from exch_scope_series import series_of, R, G
from geovac.two_center_eri import ordered_xi_general

TMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 4
TMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 0
NS = 25
r = RingQ()
name = {sp.Integer(1): "1", G: "g", sp.log(2): "l2"}
allok = True
for tau in range(TMIN, TMAX + 1):
    for (j1, j2) in KERNELS:
        t0 = time.time()
        Xs = ordered_xi_general(tau, 0, 0, 0, j1, j2, R, R)
        ser = series_of(Xs, NS + 2 * tau + 12)  # atom series carry poles to p^-(2tau+6)
        t1 = time.time()
        mine = generate_X(tau, j1, j2, NS - 8, r)
        ref = {}
        for (ce, lp), d in ser.items():
            key = "lp" if lp == 1 else name.get(ce, str(ce))
            assert not (lp == 1 and ce != 1), (ce, lp)
            for e, v in d.items():
                if e <= NS - 8:
                    ref.setdefault(key, {})[e] = v
        ok = True
        for key in set(ref) | set(mine):
            a, b = ref.get(key, {}), mine.get(key, {})
            for e in set(a) | set(b):
                if a.get(e, 0) != b.get(e, 0):
                    ok = False
        allok &= ok
        print(f"tau={tau} ({j1},{j2}) match={ok}  comps={sorted(mine)} lead={min(min(d) for d in mine.values() if d)} sympy {t1-t0:.1f}s", flush=True)
print("ALL", allok)
