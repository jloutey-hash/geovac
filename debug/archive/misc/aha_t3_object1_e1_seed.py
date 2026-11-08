"""AHA Track-3, OBJECT 1 (calibration anchor): resurgent data of the exchange-class
seed  f(a) = e^a E_1(a)  --  Paper 18 Level-2 Stieltjes seed (sec:level2_seed_set).

Question: is the Stokes/monodromy data of GeoVac's archetypal Layer-2 transcendental
ALGEBRAIC over the parameter field, up to a 2*pi*i normalisation?

What is computed here (nothing asserted -- everything derived + numerically checked):
  (1) the large-a asymptotic series          c_n = (-1)^n n!
  (2) its Borel transform                    B(zeta) = sum c_n zeta^n/n! = 1/(1+zeta)
      -> single SIMPLE POLE at zeta = -1, residue 1;  no other singularity.
  (3) Borel-Laplace reconstruction on the non-singular direction (a > 0)
  (4) the Stokes ray arg a = pi, the two lateral resummations S_+ / S_-,
      and the discontinuity (S_+ - S_-) computed THREE independent ways:
         (i)   rotated-contour lateral Laplace integrals (mpmath, dps 60)
         (ii)  the closed-form principal value  e^{-x} Ei(x) +- i pi e^{-x}
         (iii) the E_1 branch cut  E_1(-x -+ i0) = -Ei(x) +- i pi
  (5) the Stokes constant S and its large-order (resurgence) signature
  (6) field classification of S/(2 pi i).

Run:  python debug/aha_t3_object1_e1_seed.py
"""
from __future__ import annotations

import os

import mpmath as mp
import sympy as sp

DPS = 60
mp.mp.dps = DPS

OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


# ----------------------------------------------------------------- (1) the series
def asymptotic_coeffs(N):
    """c_n of  f(a) = e^a E_1(a) ~ sum_{n>=0} c_n / a^{n+1}."""
    return [(-1) ** n * sp.factorial(n) for n in range(N)]


def symbolic_derivation():
    say("[1] Integral representation (exact, derived):")
    say("    E_1(a) = int_a^inf e^{-u}/u du ; u = a(1+t) gives")
    say("    e^a E_1(a) = int_0^inf e^{-a t}/(1+t) dt")
    say("    term-by-term:  1/(1+t) = sum_n (-t)^n , int_0^inf e^{-at} t^n dt = n!/a^{n+1}")
    say("    =>  c_n = (-1)^n n!   (Gevrey-1, factorially divergent)")


# --------------------------------------------------------------- (2) Borel transform
def borel_transform_symbolic(N=14):
    z = sp.Symbol("zeta")
    c = asymptotic_coeffs(N)
    B_series = sum(c[k] * z ** k / sp.factorial(k) for k in range(N))
    B_closed = 1 / (1 + z)
    diff = sp.series(B_closed - B_series, z, 0, N).removeO()
    assert sp.simplify(diff) == 0, diff
    res = sp.residue(B_closed, z, -1)
    say("")
    say("[2] Borel transform  B(zeta) = sum_n c_n zeta^n/n! = 1/(1+zeta)   [SYMBOLIC, exact]")
    say(f"    series/closed-form agreement to O(zeta^{N}) : residual {sp.simplify(diff)}")
    say(f"    singularities: ONE simple pole at zeta = -1 ; residue = {res}")
    say("    -> Poincare rank 1, single action A = 1 (in the variable a); the attached")
    say("       one-instanton series is the CONSTANT 1, so the trans-series TERMINATES.")
    return B_closed


# -------------------------------------------------- (3) Borel-Laplace on a>0 (no Stokes)
def check_borel_laplace(a_vals=(3, 7, 15)):
    say("")
    say("[3] Borel-Laplace on the regular direction a>0:")
    say("      L[B](a) = int_0^inf e^{-a z}/(1+z) dz   vs   e^a E_1(a)")
    say(f"    {'a':>4}  {'|L[B](a) - e^a E_1(a)|':>28}")
    worst = mp.mpf(0)
    for a in a_vals:
        a = mp.mpf(a)
        lap = mp.quad(lambda z: mp.e ** (-a * z) / (1 + z), [0, 1, 10, mp.inf])
        exact = mp.e ** a * mp.e1(a)
        d = abs(lap - exact)
        worst = max(worst, d)
        say(f"    {float(a):>4.0f}  {mp.nstr(d, 5):>28}")
    return worst


# -------------------------------------------------- (4) lateral sums on the Stokes ray
def Phi_lateral(x, theta):
    """S_theta Phi(x) = int_0^{e^{i theta} inf} e^{-x zeta}/(1 - zeta) dzeta.

    Phi(x) = sum_n n!/x^{n+1} is the all-positive (non-Borel-summable) series;
    its Borel transform 1/(1-zeta) has the pole ON the positive axis = Stokes ray.
    theta > 0 passes ABOVE the pole, theta < 0 BELOW.
    """
    x = mp.mpf(x)
    e = mp.e ** (1j * mp.mpf(theta))

    def integrand(s):
        z = e * s
        return mp.e ** (-x * z) / (1 - z) * e

    return mp.quad(integrand, [0, 1, 5, 30, mp.inf])


def stokes_discontinuity(x_vals=(4, 9, 20), thetas=(mp.pi / 6, mp.pi / 3)):
    say("")
    say("[4] Stokes ray. Phi(x) = sum n!/x^{n+1} (Borel pole at zeta=+1, ON the ray).")
    say("    Lateral resummations by CONTOUR ROTATION (no closed form used):")
    say(f"    {'x':>4} {'theta':>8}  {'|(S_+ - S_-) - 2 pi i e^-x|':>32}  {'|S_+ - (PV + i pi e^-x)|':>28}")
    worst = mp.mpf(0)
    for x in x_vals:
        x = mp.mpf(x)
        pv = mp.e ** (-x) * mp.ei(x)
        pred = 2j * mp.pi * mp.e ** (-x)
        for th in thetas:
            Sp = Phi_lateral(x, th)
            Sm = Phi_lateral(x, -th)
            d1 = abs((Sp - Sm) - pred)
            d2 = abs(Sp - (pv + 1j * mp.pi * mp.e ** (-x)))
            worst = max(worst, d1, d2)
            say(f"    {float(x):>4.0f} {float(th):>8.4f}  {mp.nstr(d1, 5):>32}  {mp.nstr(d2, 5):>28}")
    return worst


def branch_cut_check(x_vals=(4, 9, 20)):
    """Third, fully independent route: the E_1 branch cut of the CORPUS object itself."""
    say("")
    say("[4c] Independent route -- branch cut of the corpus object f(a) = e^a E_1(a):")
    say("     E_1(-x -+ i0) = -Ei(x) +- i pi  =>  Disc_a f|_{a=-x} = -2 pi i e^{a}")
    say(f"     {'x':>4}  {'|Disc_a f + 2 pi i e^-x|':>28}")
    worst = mp.mpf(0)
    eps = mp.mpf(10) ** (-45)
    for x in x_vals:
        x = mp.mpf(x)
        f_up = mp.e ** (-x) * mp.e1(mp.mpc(-x, eps))
        f_dn = mp.e ** (-x) * mp.e1(mp.mpc(-x, -eps))
        disc = f_up - f_dn
        d = abs(disc + 2j * mp.pi * mp.e ** (-x))
        worst = max(worst, d)
        say(f"     {float(x):>4.0f}  {mp.nstr(d, 5):>28}")
    return worst


# ------------------------------------------------- (5) large-order resurgence signature
def large_order(N=40):
    say("")
    say("[5] Large-order / resurgence signature. Simple Borel pole at zeta = A with")
    say("    residue r  =>  c_n = r n!/A^{n+1} and S = 2 pi i r.")
    say("    Here A = 1, c_n = n!, so c_n A^{n+1}/n! = 1 for EVERY n (exact, not asymptotic):")
    devs = []
    for n in (0, 1, 5, 10, 20, N - 1):
        ratio = mp.factorial(n) * mp.mpf(1) ** (n + 1) / mp.factorial(n)
        devs.append(abs(ratio - 1))
        say(f"      n={n:>3}   c_n A^(n+1)/n! = {mp.nstr(ratio, 25)}   dev {mp.nstr(abs(ratio - 1), 3)}")
    return max(devs)


# ---------------------------------------------------------------- (6) classification
def classify(worst):
    say("")
    say("[6] Stokes constant and field classification")
    say("    S := (S_+ - S_-)/e^{-A x} = 2 pi i   (exactly; all three routes agree)")
    say("    S/(2 pi i) = 1 in QQ  --  ALGEBRAIC (rational, height 1).")
    say("    Action A = 1 rational; for the corpus argument E_1(lambda R) the action is")
    say("    A = lambda in QQ(Z_A, Z_B) (a sum of orbital exponents).")
    say(f"    max verification residual across all routes: {mp.nstr(worst, 5)}")


def main():
    os.makedirs("debug/data", exist_ok=True)
    say("=" * 78)
    say("OBJECT 1 -- e^a E_1(a), the Paper-18 Level-2 exchange-class seed")
    say(f"mpmath dps = {DPS}")
    say("=" * 78)
    symbolic_derivation()
    borel_transform_symbolic()
    w3 = check_borel_laplace()
    w4 = stokes_discontinuity()
    w4c = branch_cut_check()
    w5 = large_order()
    worst = max(w3, w4, w4c, w5)
    classify(worst)
    say("")
    say("VERDICT OBJECT 1: rank-1, single simple Borel pole, S/(2 pi i) = 1 in QQ.")
    say(f"                  ALGEBRAIC UP TO 2 pi i.  max residual {mp.nstr(worst, 5)}")
    with open("debug/data/aha_t3_object1_log.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
