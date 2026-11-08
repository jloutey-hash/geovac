"""NP record: read the Stage-C formula, homogenize with a^2, factor, restrict to F1 and compare with Stage A.
Reads data only."""
import json
import sympy as sp
from fractions import Fraction as F
d = json.load(open("data/np_stageC_dense.json"))
A, B, G1, G2, s = sp.symbols("A B G1 G2 s")  # squared rates a^2, b^2, g1^2, g2^2
X = (B, G1, G2)
def poly(dct, deg):
    out = 0
    for k, v in dct.items():
        m = json.loads(k.replace("(", "[").replace(")", "]"))
        out += sp.Rational(v) * X[0]**m[0] * X[1]**m[1] * X[2]**m[2] * A**(deg - sum(m))
    return sp.expand(out)
N = poly(d["N"], 3); D = poly(d["D"], 2)
print("N =", sp.factor(N)); print("D =", sp.factor(D))
z9 = sp.cancel(N / D)
print("z9 (homogeneous) =", sp.factor(z9))
f1 = sp.factor(z9.subs({A: 1, B: 4, G1: sp.Rational(9, 4), G2: s}))
print("restricted to F1:", f1, " Stage A:", sp.factor((s - 9) * (4 * s - 21) / (5 * (s - 4))),
      " identical:", sp.simplify(f1 - (s - 9) * (4 * s - 21) / (5 * (s - 4))) == 0)
print("at P2:", z9.subs({A: 1, B: sp.Rational(9, 4), G1: sp.Rational(49, 9), G2: sp.Rational(64, 9)}))
print("at 2*P2:", z9.subs({A: 4, B: 9, G1: sp.Rational(196, 9), G2: sp.Rational(256, 9)}))
print("symmetry (a<->b, g1<->g2):", sp.simplify(z9 - z9.subs({A: B, B: A, G1: G2, G2: G1}, simultaneous=True)) == 0)
print("symmetry (a<->g1? b<->g2):", sp.simplify(z9 - z9.subs({A: G1, G1: A, B: G2, G2: B}, simultaneous=True)) == 0)
