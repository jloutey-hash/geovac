"""Reviewer check (c)/(g): the author's engine (debug/pachucki_ext_core.py) evaluated through the reviewer's own
transcription of the l > 0 map, at a pair and R the author did not use, against the reviewer's independent
Neumann integrals (debug/pachucki_ext_review_neumann.py).  Also: f0 (label 59/60 quadrature) against the direct
Neumann value of the master integral, and the same with the sign of the e^{2bR}-e^{2aR} pair of FM1 flipped.
Usage: python debug/pachucki_ext_review_closedform.py a b R case [case ...]
"""
import sys
import time
from fractions import Fraction

import mpmath as mp
import sympy as sp

sys.path.insert(0, "debug")
sys.path.insert(0, "tests")
import pachucki_ext_core as E  # noqa: E402
import _pachucki_exchange as PX  # noqa: E402
from pachucki_ext_review_neumann import neumann  # noqa: E402

a, b, R = Fraction(sys.argv[1]), Fraction(sys.argv[2]), Fraction(sys.argv[3])
cases = sys.argv[4:]
DPS = 40
mp.mp.dps = DPS + 10
aS, bS, RS = [sp.Rational(x.numerator, x.denominator) for x in (a, b, R)]
X, Y = sp.symbols("X Y")              # X = rA^2, Y = rB^2
zA = (X + RS ** 2 - Y) / (2 * RS)
zB = zA - RS
SPEC = {
    "1s": (sp.Integer(1), aS ** 3 * bS ** 3),
    "2psig": (zA, aS ** 5 * bS ** 3),
    "3dsig": (3 * zA ** 2 - X, aS ** 7 * bS ** 3 / 18),
    "ppsig": (zA * zB, aS ** 5 * bS ** 5),
    "pppi0": (X - zA ** 2, aS ** 5 * bS ** 5 / 4),
}
alg = PX.Algebra(a, b)
x = mp.mpf(R.numerator) / R.denominator
f0 = PX.f0_numeric(alg, x)
f0_direct, _ = neumann("f0", a, b, R, DPS)
print("f0 via label 59/60 quadrature:", mp.nstr(f0, 35))
print("f0 direct (Neumann)          :", mp.nstr(f0_direct, 35), " rel diff", mp.nstr(abs(f0 / f0_direct - 1), 3))
# FM1 with the e^{2bR}E1 - e^{2aR}E1 pair's sign flipped (the bug the author reports having had in the limit script)
bad = {k: ({p: (-c if (k[0] == "Ei" and k[1] > 0) else c) for p, c in d.items()}) for k, d in alg.FM1.items()}
f0_bad = mp.quad(lambda s: PX.evaluate(alg, bad, s, 0) / s, [x, x + 1, x + 2, x + 4, x + 10, x + 30, mp.inf])
print("f0 with that sign flipped    :", mp.nstr(f0_bad, 35), " rel diff", mp.nstr(abs(f0_bad / f0_direct - 1), 3))

for case in cases:
    t0 = time.time()
    P, N2pi2 = SPEC[case]
    Pp = sp.Poly(sp.expand(P), X, Y)
    mono = list(zip(Pp.monoms(), Pp.coeffs()))
    # r^n e^{-u r} = (-1)^{n+1} d_u^{n+1} [e^{-u r}/r]; n = 2i (even) -> order 2i+1, total sign (+) over 4 factors
    terms = [((2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1), c1 * c2) for (i1, j1), c1 in mono for (i2, j2), c2 in mono]
    D = E.downset([al for al, _ in terms])
    ctx = E.Ctx("Q", a, b)
    calc, fS, Fs, P2, nfree = E.solve_all(ctx, D, check=True)
    tot = mp.mpf(0)
    for al, c in terms:
        tot += mp.mpf(c.p) / c.q * PX.evaluate(alg, fS[al], x, f0)
    Kcf = mp.mpf(sp.Rational(16) * N2pi2 / RS) * tot      # N^2 (4 pi)^2 / R
    Knum, per = neumann(case, a, b, R, DPS)
    print("%-6s (a,b,R)=(%s,%s,%s) closed form %s" % (case, a, b, R, mp.nstr(Kcf, 35)))
    print("       Neumann (independent)      %s   rel diff %s  [%d idx, free %d, %.0fs]"
          % (mp.nstr(Knum, 35), mp.nstr(abs(Kcf / Knum - 1), 3), len(D), nfree, time.time() - t0), flush=True)
