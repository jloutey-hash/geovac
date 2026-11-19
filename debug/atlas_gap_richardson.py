"""Richardson extrapolation of the large-order ratios from atlas_gap_field_series (exact coefficients)."""
import json, sys
from fractions import Fraction as Fr
from math import factorial
import mpmath as mp
sys.path.insert(0, "debug")
from atlas_gap_field_series import series, V_stark, V_diamag
mp.mp.dps = 50
KS, KD = 120, 60
Es = series(V_stark, KS); Ed = series(V_diamag, KD)
def rich(seq, N=6):
    # seq: list of (k, x_k) with x_k = L + c1/k + ...; N-th order Richardson at the last point
    ks = [k for k, _ in seq][-N-1:]; xs = [x for _, x in seq][-N-1:]
    k0 = ks[0]; tot = mp.mpf(0)
    for i in range(N+1):
        k = ks[i]
        tot += xs[i] * mp.mpf(k)**N * (-1)**(i+N) / (mp.factorial(i) * mp.factorial(N-i))
    return tot
st = [(k, mp.mpf(Es[2*k].numerator)/Es[2*k].denominator / (-(6/mp.pi)*mp.mpf(1.5)**(2*k)*mp.factorial(2*k))) for k in range(1, KS//2+1)]
# diamag: test a_k / [(-1)^{k+1} C pi^{-2k} Gamma(2k+3/2)] -> constant C; print the bare ratio with C=1
dm = [(k, mp.mpf(Ed[k].numerator)/Ed[k].denominator / ((-1)**(k+1) * mp.pi**(-2*k) * mp.gamma(2*k+mp.mpf(3)/2))) for k in range(2, KD+1)]
res = {"stark_ratio_last": str(st[-1][1]), "stark_richardson": [str(rich(st, N)) for N in (2,4,6,8)],
       "diamag_C_last": str(dm[-1][1]), "diamag_C_richardson": [str(rich(dm, N)) for N in (2,4,6,8)],
       "candidates": {"(4/pi)^(5/2)": str((4/mp.pi)**2.5), "4/pi^(5/2)*?": str(4/mp.pi**2.5)}}
print(json.dumps(res, indent=1))
json.dump(res, open("debug/data/atlas_gap_richardson.json", "w"), indent=1)
