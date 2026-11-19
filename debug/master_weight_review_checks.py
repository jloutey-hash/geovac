"""Cold review of debug/master_weight_memo.md: independent re-implementation and checks.

Independent pieces (none imported from debug/master_weight_closed_form.py):
  Ein(x) = Ei(x) - gamma - ln|x| (real x), via mpmath ei
  E2(x)  = int_0^1 Ein(x u) du / u                  (quadrature, not the series)
  Phi(x,y) = int_0^1 e^{x u} Ein(y u) du / u        (quadrature, not the double series)
  K(lam, mu): my own derivation of the Mellin finite part (same formula shape) AND a numeric
     inversion K_num = T_quad(r0) + g ln r0 + ln^2 r0 / 2 + g Ein + ln r0 Ein - E2 + Phi at small r0
  FM1: my own transcription of label 60 of arXiv:0905.1635 (source fetched by the reviewer).
Writes debug/data/master_weight_review_checks.json.
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tests"))


def Ein(x):
    x = mp.mpf(x)
    if x == 0:
        return mp.mpf(0)
    if abs(x) < mp.mpf("1e-3"):
        return mp.nsum(lambda k: x ** k / (k * mp.factorial(k)), [1, mp.inf])
    return mp.ei(x) - mp.euler - mp.log(abs(x))


def E2(x):
    return mp.quad(lambda u: Ein(x * u) / u, [0, 0.5, 1])


def Phi(x, y):
    return mp.quad(lambda u: mp.exp(x * u) * Ein(y * u) / u, [0, 0.25, 0.5, 0.75, 1])


def Kform(lam, mu):
    lam, mu = mp.mpf(lam), mp.mpf(mu)
    p = -lam
    if mu < 0:
        q = -mu
        return -((mp.euler + mp.log(q)) ** 2 / 2 + mp.pi ** 2 / 12 + mp.polylog(2, -p / q))
    return -((mp.euler + mp.log(mu)) ** 2 / 2 - 5 * mp.pi ** 2 / 12 + mp.re(mp.polylog(2, p / mu)))


def T_rest(lam, mu, R):
    """T - K."""
    g = mp.euler + mp.log(abs(mu))
    L = mp.log(R)
    e = Ein(lam * R)
    return -g * L - L ** 2 / 2 - g * e - L * e + E2(lam * R) - Phi(lam * R, mu * R)


def T_quad(lam, mu, R):
    lam, mu, R = mp.mpf(lam), mp.mpf(mu), mp.mpf(R)
    pts = [R, R + mp.mpf(1) / 2, R + 2, R + 6, R + 20, R + 60, mp.inf]
    return mp.quad(lambda s: mp.exp(lam * s) * mp.ei(mu * s) / s, pts)


def FM1_mine(a, b, s):
    """Label 60, u = a, w = b, transcribed from the fetched source."""
    u, w = a, b
    lc = mp.log(abs((w - u) / (w + u)))
    pre = 1 / (2 * (u * u - w * w))
    t1 = (mp.exp(2 * s * w) - mp.exp(2 * s * u)) * mp.ei(-2 * s * (u + w))
    t2 = -mp.exp(-2 * s * w) * (mp.ei(2 * s * (w - u)) - 2 * mp.ei(s * (w - u)) + lc)
    t3 = mp.exp(-2 * s * u) * (mp.ei(2 * s * (u - w)) - 2 * mp.ei(s * (u - w)) + lc)
    return pre * (t1 + t2 + t3)


def f0_quad_mine(a, b, R):
    pts = [R, R + mp.mpf(1) / 2, R + 2, R + 6, R + 20, R + 60, mp.inf]
    return mp.quad(lambda s: FM1_mine(a, b, s) / s, pts)


def terms(a, b):
    return [(1, 2 * b, -2 * (a + b)), (-1, 2 * a, -2 * (a + b)), (-1, -2 * b, 2 * (b - a)),
            (2, -2 * b, b - a), (1, -2 * a, 2 * (a - b)), (-2, -2 * a, a - b)]


def Cab(a, b):
    c = 1 / (2 * (a * a - b * b))
    LC = mp.log((a + b) / (b - a))
    return c * (sum(k * Kform(l, m) for k, l, m in terms(a, b)) + LC * mp.log(a / b))


def f0_closed_mine(a, b, R):
    c = 1 / (2 * (a * a - b * b))
    LC = mp.log((a + b) / (b - a))
    s = sum(k * (Kform(l, m) + T_rest(l, m, R)) for k, l, m in terms(a, b))
    s += LC * (mp.e1(2 * b * R) - mp.e1(2 * a * R))
    return c * s


def W2(a, b, R):
    c = 1 / (2 * (a * a - b * b))
    return -c * sum(k * Phi(l * R, m * R) for k, l, m in terms(a, b))


def W2_limit(a, R):
    """Memo's a -> b limit of W2: R/(2a) * d_x Phi(x, y) at x = 2aR, y = -4aR, with
    d_x Phi = int_0^1 e^{xu} Ein(yu) du = e^x Ein(y)/x - (Ein(x+y) - Ein(x))/x (integration by parts)."""
    x, y = 2 * a * R, -4 * a * R
    dxPhi = mp.exp(x) * Ein(y) / x - (Ein(x + y) - Ein(x)) / x
    dxPhi_q = mp.quad(lambda u: mp.exp(x * u) * Ein(y * u), [0, 1])
    return R / (2 * a) * dxPhi, abs(dxPhi - dxPhi_q)


def main():
    mp.mp.dps = 50
    t0 = time.time()
    out = {"dps": 50, "T": [], "Knum": [], "f0": [], "C": [], "limit": []}
    # 1. per-term T at fresh (lam, mu), both signs of mu
    for lam, mu in ((mp.mpf(-5) / 3, mp.mpf(-1) / 2), (mp.mpf(3), mp.mpf(-9) / 2),
                    (mp.mpf(-7) / 2, mp.mpf(3) / 2), (mp.mpf(-5) / 2, mp.mpf(5) / 4)):
        for R in (mp.mpf("0.25"), mp.mpf("2.3")):
            q = T_quad(lam, mu, R)
            cf = Kform(lam, mu) + T_rest(lam, mu, R)
            out["T"].append({"lam": mp.nstr(lam, 8), "mu": mp.nstr(mu, 8), "R": mp.nstr(R, 5),
                             "quad": mp.nstr(q, 40), "closed": mp.nstr(cf, 40), "absdiff": mp.nstr(abs(q - cf), 5)})
        r0 = mp.mpf("0.02")
        Kn = T_quad(lam, mu, r0) - T_rest(lam, mu, r0)
        out["Knum"].append({"lam": mp.nstr(lam, 8), "mu": mp.nstr(mu, 8), "K_formula": mp.nstr(Kform(lam, mu), 40),
                            "K_numeric_r0_0.02": mp.nstr(Kn, 40), "absdiff": mp.nstr(abs(Kn - Kform(lam, mu)), 5)})
    print("T done", time.time() - t0, flush=True)
    # 2. assembled f0 at new pairs vs my own quadrature of label 60 and vs the helper
    import _pachucki_extension as PE
    for (fa, fb) in ((Fraction(3, 4), Fraction(2)), (Fraction(1), Fraction(5))):
        a, b = mp.mpf(fa.numerator) / fa.denominator, mp.mpf(fb.numerator) / fb.denominator
        for R in (mp.mpf("0.2"), mp.mpf("0.9"), mp.mpf("2.2")):
            mine_q = f0_quad_mine(a, b, R)
            helper = PE.f0_quadrature(fa, fb, R)
            cf = f0_closed_mine(a, b, R)
            out["f0"].append({"a": str(fa), "b": str(fb), "R": mp.nstr(R, 4), "closed": mp.nstr(cf, 40),
                              "my_quad": mp.nstr(mine_q, 40), "helper_quad": mp.nstr(helper, 40),
                              "rel_closed_vs_myquad": mp.nstr(abs(cf - mine_q) / abs(mine_q), 5),
                              "rel_myquad_vs_helper": mp.nstr(abs(helper - mine_q) / abs(mine_q), 5)})
    print("f0 done", time.time() - t0, flush=True)
    # 3. C(a, b) = 0 at new pairs (incl. irrational ratio) + f0 at tiny R by my quadrature
    for a, b in ((mp.mpf(3) / 4, mp.mpf(2)), (mp.mpf(1), mp.mpf(5)), (mp.mpf(2), mp.mpf(3)),
                 (mp.mpf(1), mp.sqrt(2)), (mp.mpf(1), mp.mpf(3))):
        out["C"].append({"a": mp.nstr(a, 10), "b": mp.nstr(b, 10), "C": mp.nstr(Cab(a, b), 5)})
    out["f0_tiny_R"] = mp.nstr(f0_quad_mine(mp.mpf(3) / 4, mp.mpf(2), mp.mpf("1e-10")), 10)
    # control: C with one sign flipped in K must NOT vanish (guard sensitivity)
    a, b = mp.mpf(3) / 4, mp.mpf(2)
    c = 1 / (2 * (a * a - b * b))
    bad = c * (sum(k * Kform(l, m) for k, l, m in terms(a, b)) - mp.log((a + b) / (b - a)) * mp.log(a / b))
    out["C_control_LC_sign_flipped"] = mp.nstr(bad, 10)
    print("C done", time.time() - t0, flush=True)
    # 4. a -> b: W2 tends to the weight-one limit; (W2 - lim)/eps should be bounded
    a = mp.mpf(1)
    for R in (mp.mpf("0.7"), mp.mpf("1.9")):
        lim, ibp_err = W2_limit(a, R)
        for eps in (mp.mpf("1e-2"), mp.mpf("1e-3"), mp.mpf("1e-4"), mp.mpf("1e-5")):
            w = W2(a, a + eps, R)
            out["limit"].append({"R": mp.nstr(R, 3), "eps": mp.nstr(eps, 2), "W2": mp.nstr(w, 25),
                                 "limit": mp.nstr(lim, 25), "diff_over_eps": mp.nstr((w - lim) / eps, 10),
                                 "ibp_check": mp.nstr(ibp_err, 5)})
    out["seconds"] = round(time.time() - t0, 1)
    with open(os.path.join(HERE, "data", "master_weight_review_checks.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
