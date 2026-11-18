"""Reviewer check: the jet engine's sources for multi-indices with REPEATED exponents (the part of
debug/pachucki_ext_core.py that the reviewed helper never exercised) against plain sympy differentiation of the
helper's transcription of labels 20, 30-34, exactly, at two pairs.  Wrong answer excluded: a jet rule (Ei, exp,
log, inverse, factorial normalisation) that is right on subsets (all orders <= 1) but wrong at order >= 2.
"""
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "debug")
sys.path.insert(0, "tests")
import pachucki_ext_core as E  # noqa: E402
import _pachucki_exchange as PX  # noqa: E402

ALS = [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2), (3, 0, 0, 0), (0, 3, 0, 0), (2, 1, 0, 0), (1, 2, 1, 0),
       (2, 0, 2, 0), (0, 2, 0, 2), (2, 1, 1, 1), (1, 1, 1, 2), (3, 1, 1, 1), (1, 3, 1, 1), (2, 2, 1, 1), (3, 1, 3, 1)]
for a, b in ((Fraction(1), Fraction(2)), (Fraction(4, 5), Fraction(3, 2))):
    t0 = time.time()
    ctx = E.Ctx("Q", a, b)
    D = E.downset(ALS)
    Fs, P2 = E.sources(ctx, D)
    alg = PX.Algebra(a, b)
    sub = PX.exchange_point(sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator))
    okF = okP = True
    for al in ALS:
        e, q = PX._dsym(al)
        F_ref = alg.parse(e.subs(sub))
        p_ref = PX._frac(sp.nsimplify(q.subs(sub)))
        f1 = {k: {p: Fraction(c) for p, c in d.items()} for k, d in Fs[al].items()}
        sF, sP = (PX.fadd(f1, F_ref, -1) == {}), (P2[al] == p_ref)
        okF, okP = okF and sF, okP and sP
        print("  ", (str(a), str(b)), al, "F jet == sympy:", sF, " p^2 jet == sympy:", sP, flush=True)
    print((str(a), str(b)), "all F:", okF, "all p^2:", okP, "%.0fs" % (time.time() - t0), flush=True)
