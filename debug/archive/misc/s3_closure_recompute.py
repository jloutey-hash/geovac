"""S^(3) closure sprint — corrected canonical value of the decomposition.

The stage-1 figure S3 = 30.61538052881950... (s3_decomp_numerics.json) was
computed with the broken trailing-1 evaluator: the four table entries
t3(2,1,1), t3(2,3,1), t3(4,1,1), t3(4,3,1) were good to only ~12 digits,
and with table coefficients up to 3072 the assembled value was good to
~1e-9 at best.

This driver re-assembles the decomposition entirely from the 220-dps
cache (debug/data/s3_pslq_cache.json):
  - 26 double t-values     (Levin on pure algebraic decay — trusted),
  - 8 triple t-values b3>=2 (same),
  - 4 triple t-values b3=1  (corrected Abel-summation values written by
                             s3_closure_trailing1.py),
  - lambda / lambda-product / Hurwitz closed forms,
  - S_min via its identified closed form (Paper 28 eq:smin_closed) AND
    the Levin Hurwitz-tail anchor, cross-checked.

Effective precision ~200 digits (cache strings are dps-5 = 215 digits).

Output: debug/data/s3_closure_recompute.json  (canonical S^(3) value)
"""
from __future__ import annotations

import ast
import json
from pathlib import Path

import mpmath

HERE = Path(__file__).parent
OUT = {}

DPS = 200
mpmath.mp.dps = DPS
pi = mpmath.pi

CACHE = json.loads((HERE / "data" / "s3_pslq_cache.json").read_text())


def cval(key):
    k = "%s@220" % key
    if k not in CACHE:
        raise KeyError("cache miss: %s" % k)
    return mpmath.mpf(CACHE[k])


def lam(s):
    return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)


tab = json.loads((HERE / "data" / "s3_decomp_table.json").read_text())
lin = {ast.literal_eval(k): mpmath.mpf(int(v[0])) / int(v[1])
       for k, v in tab["S3_linear_table"].items()}


def sym_val(k):
    if k == ('one',):
        return mpmath.mpf(1)
    if k[0] == 'lam':
        return lam(k[1])
    if k[0] == 'lam2':
        return lam(k[1]) * lam(k[2])
    if k[0] == 't2':
        return cval("t2(%d,%d)" % k[1:])
    if k[0] == 't3':
        return cval("t3(%d,%d,%d)" % k[1:])
    raise ValueError(k)


S3_linear = mpmath.fsum(c * sym_val(k) for k, c in lin.items())

hz = lambda s: mpmath.zeta(s, mpmath.mpf(5) / 2)
P = pi ** 2 - pi ** 4 / 12 - mpmath.mpf(64) / 81
Q = 4 * hz(4) - 2 * hz(6) + mpmath.mpf(1) / 4 * hz(8)
R3 = 8 * hz(6) - 6 * hz(8) + mpmath.mpf(3) / 2 * hz(10) \
    - mpmath.mpf(1) / 8 * hz(12)

# S_min two ways: identified closed form (Paper 28 eq:smin_closed) + Levin
ln2, z3, z5 = mpmath.log(2), mpmath.zeta(3), mpmath.zeta(5)
S_min_closed = (8 * pi ** 2 * ln2 - mpmath.mpf(2) / 3 * pi ** 4 * ln2
                - 3 * pi ** 2 * z3 + mpmath.mpf(1) / 2 * pi ** 4 * z3
                - mpmath.mpf(5) / 2 * pi ** 2 * z5
                + pi ** 6 / 4 - mpmath.mpf(3) / 2 * pi ** 4 - pi ** 8 / 96)
S_min_levin = mpmath.nsum(
    lambda k: (2 * mpmath.zeta(2, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)
               - mpmath.mpf(1) / 2
               * mpmath.zeta(4, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)) ** 2,
    [1, mpmath.inf], method='levin')
d_smin = abs(S_min_closed - S_min_levin)
OUT["smin_closed_vs_levin"] = mpmath.nstr(d_smin, 4)
assert d_smin < mpmath.mpf(10) ** -(DPS - 30), "S_min cross-check FAILED"
print("S_min closed-form vs Levin: %s  PASS" % mpmath.nstr(d_smin, 4))

S3 = S3_linear - 16 * S_min_closed * P - 8 * P ** 3 + 8 * P * Q + R3

STALE = mpmath.mpf("30.61538052881950")
print("\nS^(3) corrected  = %s" % mpmath.nstr(S3, 60))
print("S^(3) stale      = %s" % mpmath.nstr(STALE, 20))
print("|corrected - stale| = %s" % mpmath.nstr(abs(S3 - STALE), 4))

OUT["S3_canonical_200dps"] = mpmath.nstr(S3, 200)
OUT["S3_minus_stale"] = mpmath.nstr(abs(S3 - STALE), 4)
OUT["S3_linear_part"] = mpmath.nstr(S3_linear, 60)
(HERE / "data" / "s3_closure_recompute.json").write_text(
    json.dumps(OUT, indent=2))
print("Saved:", HERE / "data" / "s3_closure_recompute.json")
