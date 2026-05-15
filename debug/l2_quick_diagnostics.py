"""Quick diagnostic checks against the c value at 14-15 digit precision.

Uses MR-C 80-dps c = 4.1093214674877940927579607260741005838057691088362615503253972964276017819113301...

Computes c against a few targeted candidates that aren't in the basis but
might be obvious closed forms:
  - 4 * (some natural constant)
  - small rational
  - involves Glaisher A, Khintchine, etc.
"""
from __future__ import annotations

import sys
from pathlib import Path
PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ))

import mpmath
from mpmath import mp, mpf, pi, log, sqrt

mp.dps = 100
c = mpf("4.1093214674877940927579607260741005838057691088362615503253972964276017819113301")

print(f"c = {mpmath.nstr(c, 40)}")
print()

# Decompose
print(f"c - 4 = {mpmath.nstr(c - 4, 40)}")
print(f"c - 4*log(2) = {mpmath.nstr(c - 4*log(mpf(2)), 40)}")
print(f"c / pi = {mpmath.nstr(c / pi, 40)}")
print(f"c * pi = {mpmath.nstr(c * pi, 40)}")
print(f"c * pi^2 = {mpmath.nstr(c * pi**2, 40)}")
print(f"c / log(2) = {mpmath.nstr(c / log(mpf(2)), 40)}")
print(f"c - pi = {mpmath.nstr(c - pi, 40)}")
print(f"c - 4/pi - pi = {mpmath.nstr(c - 4/pi - pi, 40)}")
print(f"c - 4*log(2) - pi/2 = {mpmath.nstr(c - 4*log(mpf(2)) - pi/2, 40)}")
print()

# Some L-values
print("Try simple closed forms:")
candidates = [
    ("(4/pi) * (gamma_E + log(8*pi))", (4/pi) * (mpmath.euler + log(8*pi))),
    ("(4/pi) * (1 + log(8*pi)) - some", (4/pi) * (1 + log(8*pi))),
    ("4*log(pi*e^gamma_E/2)/pi", 4*log(pi*mpmath.exp(mpmath.euler)/2)/pi),
    ("4/pi + 4/pi * log(2)", 4/pi + 4/pi * log(mpf(2))),
    ("4 - pi/exp(1)", mpf(4) - pi/mpmath.exp(1)),
    ("4*log(2)/pi + 4*gamma_E/pi", 4*log(mpf(2))/pi + 4*mpmath.euler/pi),
    ("(4/pi)*(1 + 2*log(2))", (4/pi)*(1 + 2*log(mpf(2)))),
    ("4*log(2) + 4*gamma_E/pi - 1/pi^2", 4*log(mpf(2)) + 4*mpmath.euler/pi - 1/pi**2),
    # Some involving Catalan:
    ("4/pi * (gamma_E + 2*log(2))", (4/pi)*(mpmath.euler + 2*log(mpf(2)))),
    ("4*Catalan/pi + 4/pi * log(2)", 4*mpmath.catalan/pi + 4/pi * log(mpf(2))),
    ("4 - 4*Catalan/pi^2", 4 - 4*mpmath.catalan/pi**2),
    # 4 * ln(asymptotic constant in something)
    ("4 + (gamma_E + log(pi))/pi", 4 + (mpmath.euler + log(pi))/pi),
    # SU(2) Haar normalization
    ("4*log(pi)/pi", 4*log(pi)/pi),
    ("4/pi*log(pi)", 4/pi*log(pi)),
    ("4/pi*log(2*pi)", 4/pi*log(2*pi)),
    ("4/pi*(2*log(2)+gamma_E)", 4/pi*(2*log(mpf(2)) + mpmath.euler)),
    # Glaisher-Kinkelin
    ("4*log(A)/pi", 4*mpmath.log(mpmath.glaisher)/pi),
    ("12*log(A)/pi", 12*mpmath.log(mpmath.glaisher)/pi),
    ("c -log(A)*16/pi", c - 16*mpmath.log(mpmath.glaisher)/pi),
]
for label, v in candidates:
    delta = c - v
    rel = abs(delta) / abs(c)
    if abs(delta) < mpf("1e-10"):
        flag = " <-- HIT"
    elif abs(delta) < mpf("1e-5"):
        flag = " <-- close"
    else:
        flag = ""
    print(f"  c - ({label}) = {mpmath.nstr(delta, 12)}{flag}")

# Check c - 4*gamma_E/pi - 4*log(2)/pi - 4/pi vs. something
print()
combo = (4 + 4*mpmath.euler + 4*log(mpf(2))) / pi
print(f"(4 + 4*gamma_E + 4*log(2))/pi = {mpmath.nstr(combo, 40)}")
print(f"vs c                          = {mpmath.nstr(c, 40)}")
print(f"diff                          = {mpmath.nstr(c - combo, 40)}")
