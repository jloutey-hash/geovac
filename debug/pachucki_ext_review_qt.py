"""Reviewer check (a): what the Q(t) run needs in order to be a proof at EVERY real t > 1.

For each multi-index alpha (not only for K) collect
  * the irreducible factors of every coefficient denominator of f_alpha, of the source F_alpha and of d_alpha p^2;
  * every function key of f_alpha with its decay rate (lam for E keys, lam + mu for Ei keys) as a linear form in t,
    and whether that rate is < 0 for all t > 1 (uniqueness-by-decay is used at every alpha, pointwise in t);
  * zero Ei arguments mu (Ei(0) is undefined) and E keys with rate 0 (non-decaying).
Also: rank of K..K^(8) in the nine coordinates at several t (including the key-collision values t = 2, 3).
Usage: python debug/pachucki_ext_review_qt.py CASE     (CASE = 1s | 2psig | ...)
"""
import sys
import time

import sympy as sp

sys.path.insert(0, "debug")
import pachucki_ext_core as E  # noqa: E402

case = sys.argv[1] if len(sys.argv) > 1 else "1s"
t0 = time.time()
ctx = E.Ctx("Qt")
KK, tS = ctx.K, ctx.t
GEN = {"1s": [(1, 1, 1, 1)], "2psig": [(3, 1, 3, 1), (1, 3, 1, 3), (3, 1, 1, 3), (1, 3, 3, 1), (3, 1, 1, 1),
                                        (1, 1, 3, 1), (1, 3, 1, 1), (1, 1, 1, 3), (1, 1, 1, 1)]}
D = E.downset(GEN[case])
calc, fS, Fs, P2, nfree = E.solve_all(ctx, D, check=True)
print(case, "multi-indices", len(D), "free", nfree, "solved %.0fs" % (time.time() - t0), flush=True)


def den_factors(x, acc):
    d = KK.to_sympy(x)
    for f, _ in sp.factor_list(sp.denom(sp.together(d)))[1]:
        acc.add(str(f))


dens_f, dens_F, dens_p = set(), set(), set()
keys = set()
for al in D:
    for k, d in fS[al].items():
        keys.add(k)
        for c in d.values():
            den_factors(c, dens_f)
    for k, d in Fs[al].items():
        keys.add(k)
        for c in d.values():
            den_factors(c, dens_F)
    if P2[al] != 0:
        den_factors(P2[al], dens_p)
print("denominator factors: f_alpha", sorted(dens_f), "| F_alpha", sorted(dens_F), "| d p^2", sorted(dens_p))
bad = []
for k in sorted(keys, key=str):
    if k[0] == "f":
        continue
    lam = KK.to_sympy(k[1])
    rate = lam if k[0] == "E" else sp.expand(lam + KK.to_sympy(k[2]))
    r1, slope = rate.subs(tS, 1), sp.diff(rate, tS)
    assert sp.diff(rate, tS, 2) == 0
    decays = (r1 <= 0 and slope <= 0 and not (r1 == 0 and slope == 0))
    extra = ""
    if k[0] == "Ei":
        mu = KK.to_sympy(k[2])
        roots = [r for r in sp.solve(mu, tS) if r.is_real and r > 1]
        if roots or mu == 0:
            extra = " Ei ARGUMENT VANISHES at t=%s" % roots
    if not decays or extra:
        bad.append((k, rate, extra))
    print("  key", (k[0], KK.to_sympy(k[1]), KK.to_sympy(k[2]) if k[0] == "Ei" else k[2]), " rate", rate,
          " decays for all t>1:", decays, extra)
print("keys failing decay / with vanishing Ei argument for some t > 1:", bad)
print("done %.0fs" % (time.time() - t0))
