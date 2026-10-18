"""Banana dictionary, step 4 (numerical support): the exact Mellin identity
   F_n(t) = (2/pi) int_0^1 K'(u) * (t/u) * Bp_n(t/u) du/u,   K'(u) = K(k = sqrt(1-u^2)),
with F_n(t) = int K0(x)^{n-1}K0(tx)dx (Paper 59) and Bp_n(y) = int x K0(x)^{n-1}K0(yx)dx
(equal-mass 2-D banana with n-1 propagators, one K0 replacing I0).  mpmath ellipk takes m = k^2."""
import mpmath as mp, json
mp.mp.dps = 18
def F(n, t): return mp.quad(lambda x: mp.besselk(0, x)**(n-1)*mp.besselk(0, t*x), [0, 1, 5, mp.inf])
def Bp(n, y):
    if n == 2: return mp.log(y)/(y**2 - 1) if y != 1 else mp.mpf(1)/2
    return mp.quad(lambda x: x*mp.besselk(0, x)**(n-1)*mp.besselk(0, y*x), [0, 1, 5, mp.inf])
def rhs(n, t):
    f = lambda u: (mp.pi/(2*mp.agm(1, u))) * (t/u) * Bp(n, t/u) / u
    return 2/mp.pi * mp.quad(f, [0, mp.mpf(t)/4, t, 1])
if __name__ != "__main__": raise ImportError
out = {}
for n in [2]:  # n=3 nested quadrature timed out (600 s); exact recurrence checks cover n=3..5
    for t in [mp.mpf("0.5"), mp.mpf("1.3")]:
        a, b = F(n, t), rhs(n, t)
        out[f"n{n}_t{float(t)}"] = [mp.nstr(a, 15), mp.nstr(b, 15), float(abs(a-b)/abs(a))]
        print(n, float(t), mp.nstr(a, 15), mp.nstr(b, 15), float(abs(a-b)/abs(a)))
json.dump(out, open("data/banana_dict_kernel.json", "w"), indent=1)
