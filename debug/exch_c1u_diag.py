"""Locate the C1u disagreement at unequal zeta WITHOUT printing any coefficient:
evaluate (a) generator series, (b) scope sympy series, (c) the sympy closed form
itself, numerically at R = 0.4 (tau <= 6 throughout); print only agreement flags and
log10 of relative differences."""
import sys
sys.path.insert(0, "debug"); sys.path.insert(0, "."); sys.path.insert(0, "tests")
import mpmath as mp
from fractions import Fraction
from exch_gen import RingQ, generate
from exch_scope_series import series_of, R, G
from exch_controls import map_sym, to_frac_dict, S1
from geovac import qfd_core as Q
import sympy as sp

mp.mp.dps = 60
consts = {"1": mp.mpf(1), "g": +mp.euler, "l2": mp.log(2), "l3": mp.log(3), "l5": mp.log(5)}


def ev(ser, x):
    t = mp.mpf(0)
    for (c, m), d in ser.items():
        t += consts[c] * mp.log(x) ** m * sum(mp.mpf(v.numerator) / v.denominator * x ** k
                                               for k, v in d.items())
    return t


for (za, zb) in ((1, Fraction(3, 2)), (1, 2)):
    x = mp.mpf("0.4")
    gen = to_frac_dict(generate(za, zb, 60, RingQ(), tau_cap=6)[0])
    f = Q.exchange_closed_form(Fraction(za), S1, S1, Fraction(zb), S1, S1, R, tau_max=6)
    fx = sp.N(f.subs(R, sp.Rational(2, 5)), 50)
    ref = map_sym(series_of(f, 80), 60)
    a, b, c = ev(gen, x), ev(ref, x), mp.mpf(str(fx))
    rel = lambda u, v: mp.nstr(mp.log10(abs(u - v) / abs(v) + mp.mpf(10) ** -60), 3)
    print(za, zb, "log10 rel: gen-vs-closed", rel(a, c), " sympyseries-vs-closed", rel(b, c), flush=True)
