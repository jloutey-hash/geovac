"""BBBG (arXiv:0801.0891) conjecture 11, eq. (c80): c_{8,0} =? pi^2 int I0^2 (6 K0^2 - pi^2 I0^2) K0^4 dt.
Low-precision check (algebraic t^-4 tail; tanh-sinh on [0, oo))."""
import json
import mpmath as mp
mp.mp.dps = 40
MOM = json.load(open("debug/data/ising_class_moments_arb.json"))["moments"]
c80 = mp.mpf(MOM["8,0"]["value"])
def f(t):
    i0 = mp.besseli(0, t); k0 = mp.besselk(0, t)
    return i0**2 * (6 * k0**2 - mp.pi**2 * i0**2) * k0**4
pts = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
v = mp.quad(f, pts)
# tail beyond 512: integrand ~ sum a_j t^-(4+2j); use the large-t expansion I0 K0 ~ 1/(2t)(1 + 1/(8t^2)+...) numerically fitted is avoided:
# instead integrate the tail with the substitution t = 512/u on (0,1]
tail = mp.quad(lambda u: f(512 / u) * 512 / u**2, [0, 0.5, 1])
val = mp.pi**2 * (v + tail)
rel = abs(val / c80 - 1)
print("c80", mp.nstr(c80, 30)); print("rhs", mp.nstr(val, 30)); print("rel", mp.nstr(rel, 3))
json.dump({"c80": mp.nstr(c80, 35), "rhs": mp.nstr(val, 35), "rel": mp.nstr(rel, 3), "dps": 40},
          open("debug/data/ising_class_c80_sumrule.json", "w"), indent=1)
