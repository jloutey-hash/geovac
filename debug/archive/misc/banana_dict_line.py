"""Banana dictionary, step 7 (numerical support): F_n(t) = int_0^inf dp B_n(p; t,1,..,1), where
B_n(p) = int_0^inf x J0(p x) K0(t x) K0(x)^{n-1} dx is the n-propagator 2-D banana at Euclidean
momentum p (masses t,1,...,1).  Formal proof: int_0^inf J0(p x) dp = 1/x.  Checked at n = 2 with the
closed-form bubble  B_2(p) = ln((a^2+b^2+p^2+sqrt(D))/(2ab))/sqrt(D),  D = (a^2+b^2+p^2)^2 - 4a^2b^2."""
import mpmath as mp, json
mp.mp.dps = 25
out = {}
for t in [mp.mpf("0.5"), mp.mpf("1.3")]:
    a, b = mp.mpf(1), t
    def B2(p):
        D = (a*a + b*b + p*p)**2 - 4*a*a*b*b
        return mp.log((a*a + b*b + p*p + mp.sqrt(D))/(2*a*b))/mp.sqrt(D)
    lhs = mp.quad(lambda x: mp.besselk(0, x)*mp.besselk(0, t*x), [0, 1, mp.inf])
    rhs = mp.quad(B2, [0, 1, 10, mp.inf])
    out[str(t)] = [mp.nstr(lhs, 20), mp.nstr(rhs, 20), float(abs(lhs - rhs)/lhs)]
    print(float(t), out[str(t)])
json.dump(out, open("data/banana_dict_line.json", "w"), indent=1)
