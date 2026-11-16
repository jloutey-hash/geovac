import sys; sys.path.insert(0,"debug"); sys.path.insert(0,"."); sys.path.insert(0,"tests")
from exch_scope_series import *
for p in (1000003, 998244353):
    f = Q.exchange_closed_form(Fraction(1), S1, S1, Fraction(1), S1, S1, R, tau_max=6)
    ser = series_of(f, 70)
    for r, d in [(3, 20), (4, 7), (4, 8), (4, 10), (5, 8)]:
        ns = guess_mod_p(ser, r, d, p, mmax=60)
        msg = f"p={p} r={r} d={d}: nullity {len(ns)}"
        if (r, d) == (4, 8) and len(ns) == 1:
            msg += f"; prop. OP_EXCH_1: {proportional(ns[0], op_mod_p(OP_EXCH_1, p, 8), p)}"
        print(msg)
    # how few terms suffice for (4,8)?
    for mm in (20, 25, 30, 40):
        print(" mmax", mm, "nullity", len(guess_mod_p(ser, 4, 8, p, mmax=mm)))
# planted non-holonomic-looking control: tan(R) + lnR*R^3/(e^R-1)? use series of tan
