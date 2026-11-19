"""R2 (top constants) in float64 with exponentially scaled Bessel functions and Gauss-Laguerre on [x, oo).
B1: h_2 root - p_z - p_z (predicted C = pi/15); D: s root with two h_2 leaves (predicted C = pi/10).
G^top / (x e^{-2x}) -> C, exponent (5-n)/2 = 1 at n = 3. Reviewer's own code."""
import json
import numpy as np
from scipy.special import kve, roots_laguerre

vn, vw = roots_laguerre(80)


def lag(f, x, k):
    """int_x^oo f(s) e^{-k (s-x)} ds, f smooth (already scaled)."""
    v = vn / k
    return np.sum(vw / k * f(x + v))


def E1s(s):  # e^{s} * exterior tail of leaf p_z potential (channel 1), W1 = -K1
    s = np.atleast_1d(s)
    out = []
    for si in s:
        A = lag(lambda u: u ** 3 * kve(1, u), si, 1.0)       # e^{s} int_s^oo u^3 K1
        out.append((si ** -2 * A - si * kve(0, si)) / 3)     # (1/3)[-s^-2 int u^3 W1 + s int W1], int_s K1 = K0
    return np.array(out)


def W2s(s):  # e^{2s} * W2^top
    return (2.0 / 3.0) * (-kve(1, s)) * E1s(s)


def Gtop_B1(x):  # e^{2x} * G^top
    A = lag(lambda s: s ** 4 * W2s(s), x, 2.0)
    B = lag(lambda s: s ** -1 * W2s(s), x, 2.0)
    E = (-x ** -3 * A + x ** 2 * B) / 5
    dE = (3 * x ** -4 * A + 2 * x * B) / 5
    return (-x ** 2 * W2s(np.array([x]))[0] + 3 * x * dE + 9 * E) / 5


def E2leaf(x):  # e^{x} * exterior tail of h_2 leaf (channel 2, W = K2)
    A = lag(lambda s: s ** 4 * kve(2, s), x, 1.0)
    B = lag(lambda s: s ** -1 * kve(2, s), x, 1.0)
    return (-x ** -3 * A + x ** 2 * B) / 5


def Gtop_D(x):
    return x ** 2 * E2leaf(x) ** 2 / 5


out = {}
for name, fn, pred in (('B1', Gtop_B1, np.pi / 15), ('D', Gtop_D, np.pi / 10)):
    xs = [50.0, 100.0, 200.0, 400.0]
    vals = [fn(x) / x for x in xs]
    r1 = [2 * vals[i + 1] - vals[i] for i in range(3)]
    r2 = [(4 * r1[i + 1] - r1[i]) / 3 for i in range(2)]
    r3 = (8 * r2[1] - r2[0]) / 7
    slope = np.log(fn(800.0) / fn(400.0)) / np.log(2)
    out[name] = dict(raw=vals, richardson3=r3, predicted=pred, rel=abs(r3 - pred) / pred,
                     local_exponent_400_800=slope,
                     wrong_principal_no_n1_factor=pred / 4 if name == 'B1' else None,
                     wrong_exponent_shift_by_degree=1 + (2 if name == 'B1' else 4))
    print(name, out[name], flush=True)
json.dump(out, open('debug/data/tree_attain_l_review_tailfast.json', 'w'), indent=1)
