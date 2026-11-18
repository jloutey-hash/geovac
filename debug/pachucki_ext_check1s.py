"""Check the jet/direct-solver engine against the reviewed helper tests/_pachucki_exchange.py at (1,2), (1,3/2)."""
import sys, time
from fractions import Fraction
sys.path.insert(0, "debug"); sys.path.insert(0, "tests")
import pachucki_ext_core as E
import _pachucki_exchange as PX

for a, b in ((1, 2), (1, Fraction(3, 2))):
    t0 = time.time()
    ctx = E.Ctx("Q", a, b)
    D = E.downset([(1, 1, 1, 1)])
    calc, fS, Fs, P2, nfree = E.solve_all(ctx, D)
    t1 = time.time()
    alg = PX.Algebra(a, b)
    fS2, free2, Fs2, P22 = alg.solve_all()
    ok_src = all(Fs[al] == Fs2[al] and P2[al] == P22[al] for al in D)
    ok_f = all(fS[al] == fS2[al] for al in D)
    print((a, b), "sources equal:", ok_src, " solutions equal:", ok_f, " free:", nfree, "%.1fs (helper %.1fs)" % (t1 - t0, time.time() - t1))
