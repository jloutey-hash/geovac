"""S^(3) round-2 decomposition — Step 3: high-precision numerical verification.

Two INDEPENDENT evaluations of S^(3) at physical exponents (a=4, p=1):

  (A) Decomposition route: the exact rational table from s3_decomp_engine
      (12 Hoffman triple t-values + double t-values + lambda terms) with
      every t-value evaluated by iterated Hurwitz tails + mpmath nsum
      (Levin u), plus the closed-form product terms
      -16 S_min P - 8 P^3 + 8 P Q + R3.

  (B) Factorized closed-form anchor: S3 = sum_m phi(m) L(m)^2 with
      L(m) = 2 A(m) - P - phi(m) and A(m) in closed digamma/Hurwitz form
      (the channel-count closed form turns Paper 28's O(N^3) factorized
      algorithm into an O(1)-per-term single sum). Summed by Euler-Maclaurin
      (mpmath sumem) with Levin cross-check.

NO brute-force triple sum is used anywhere (charter: the raw sum converges
as N^-1.31; the decomposition IS the method).

Numerical formulas (o2 = 2j+1):
  t3(b1,b2,b3) = sum_{j>=1} (2j+1)^-b2 * tau_{b1}(j) * pL_{b3}(j)
     tau_b(j)  = 2^-b zeta(b, j+3/2)              [tail over largest var]
     pL_c(j)   = lam(c) - 2^-c zeta(c, j+1/2)     [partial over smallest]
     pL_1(j)   = (psi0(j+1/2) - psi0(1/2))/2
  t2(b,c)      = sum_{j>=0} (2j+1)^-c * 2^-b zeta(b, j+3/2)

Output: debug/data/s3_decomp_numerics.json
"""
from __future__ import annotations

import ast
import json
import time
from pathlib import Path

import mpmath

DPS = 130
mpmath.mp.dps = DPS
OUT = {"dps": DPS}
t_start = time.time()

pi = mpmath.pi
HERE = Path(__file__).parent


def lam(s):
    return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)


PSI0_HALF = mpmath.digamma(mpmath.mpf(1) / 2)
LAM = {s: lam(s) for s in range(1, 13) if s != 1}


def t2val(b, c):
    def term(j):
        j = int(j)
        return (mpmath.mpf(2 * j + 1) ** -c * mpmath.power(2, -b)
                * mpmath.zeta(b, mpmath.mpf(j) + mpmath.mpf(3) / 2))
    return mpmath.nsum(term, [0, mpmath.inf], method='levin')


def t3val(b1, b2, b3):
    lam3 = LAM.get(b3)

    def term(j):
        j = int(j)
        tau = mpmath.power(2, -b1) * mpmath.zeta(
            b1, mpmath.mpf(j) + mpmath.mpf(3) / 2)
        if b3 == 1:
            pL = (mpmath.digamma(mpmath.mpf(j) + mpmath.mpf(1) / 2)
                  - PSI0_HALF) / 2
        else:
            pL = lam3 - mpmath.power(2, -b3) * mpmath.zeta(
                b3, mpmath.mpf(j) + mpmath.mpf(1) / 2)
        return mpmath.mpf(2 * j + 1) ** -b2 * tau * pL
    return mpmath.nsum(term, [1, mpmath.inf], method='levin')


# ---------------------------------------------------------------- load table
tab = json.loads((HERE / "data" / "s3_decomp_table.json").read_text())
lin = {ast.literal_eval(k): mpmath.mpf(int(v[0])) / int(v[1])
       for k, v in tab["S3_linear_table"].items()}

t3_syms = sorted(k for k in lin if k[0] == 't3')
t2_syms = sorted(k for k in lin if k[0] == 't2')

print("Evaluating %d double t-values (Levin, %d dps) ..." %
      (len(t2_syms), DPS))
T2 = {}
for k in t2_syms:
    t0 = time.time()
    T2[k] = t2val(k[1], k[2])
    print("  t2%s = %s  (%.0fs)" % (k[1:], mpmath.nstr(T2[k], 25),
                                    time.time() - t0))

print("Evaluating %d triple t-values (Levin, %d dps) ..." %
      (len(t3_syms), DPS))
T3 = {}
for k in t3_syms:
    t0 = time.time()
    T3[k] = t3val(k[1], k[2], k[3])
    print("  t3%s = %s  (%.0fs)" % (k[1:], mpmath.nstr(T3[k], 25),
                                    time.time() - t0))

OUT["t2_values"] = {repr(k): mpmath.nstr(v, 110) for k, v in T2.items()}
OUT["t3_values"] = {repr(k): mpmath.nstr(v, 110) for k, v in T3.items()}


def sym_val(k):
    if k == ('one',):
        return mpmath.mpf(1)
    if k[0] == 'lam':
        return LAM[k[1]]
    if k[0] == 'lam2':
        return LAM[k[1]] * LAM[k[2]]
    if k[0] == 't2':
        return T2[k]
    if k[0] == 't3':
        return T3[k]
    raise ValueError(k)


S3_linear = sum(c * sym_val(k) for k, c in lin.items())

# product terms
hz = lambda s: mpmath.zeta(s, mpmath.mpf(5) / 2)
P = pi ** 2 - pi ** 4 / 12 - mpmath.mpf(64) / 81
Q = 4 * hz(4) - 2 * hz(6) + mpmath.mpf(1) / 4 * hz(8)
R3 = 8 * hz(6) - 6 * hz(8) + mpmath.mpf(3) / 2 * hz(10) \
    - mpmath.mpf(1) / 8 * hz(12)

print("S_min anchor (Levin) ...")
S_min = mpmath.nsum(
    lambda k: (2 * mpmath.zeta(2, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)
               - mpmath.mpf(1) / 2
               * mpmath.zeta(4, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)) ** 2,
    [1, mpmath.inf], method='levin')
FROZEN80 = ("2.4799369380342225544135795008293821446879257866172884583787"
            "98726559552777818374")
assert mpmath.nstr(S_min, 80).startswith(FROZEN80[:75]), "S_min guard FAILED"
OUT["S_min"] = mpmath.nstr(S_min, 110)

S3_dec = S3_linear - 16 * S_min * P - 8 * P ** 3 + 8 * P * Q + R3
print("\n(A) S3 decomposition route = %s" % mpmath.nstr(S3_dec, 70))
OUT["S3_decomposition"] = mpmath.nstr(S3_dec, 110)


# ---------------------------------------------------------------- anchor (B)
def A_closed(m):
    m = mpmath.mpf(int(m))
    z2m = mpmath.zeta(2, m + mpmath.mpf(5) / 2)
    z3m = mpmath.zeta(3, m + mpmath.mpf(5) / 2)
    z4m = mpmath.zeta(4, m + mpmath.mpf(5) / 2)
    return (2 * (mpmath.digamma(m + mpmath.mpf(5) / 2)
                 - mpmath.digamma(mpmath.mpf(5) / 2))
            - 3 * (hz(2) - z2m) - mpmath.mpf(1) / 2 * (hz(3) - z3m)
            + mpmath.mpf(3) / 4 * (hz(4) - z4m)
            + m * (2 * z2m - mpmath.mpf(1) / 2 * z4m))


def phi_f(m):
    x = mpmath.mpf(int(m)) + mpmath.mpf(3) / 2
    return 2 / x ** 2 - mpmath.mpf(1) / 2 / x ** 4


# closed-form sanity: A(1) must equal P + ... check directly vs partial sums
A_dir = lambda m, N: sum(min(n, m) * phi_f(n) for n in range(1, N + 1))
for m in (1, 2, 5):
    d = abs(A_closed(m) - (A_dir(m, 100000) + m * mpmath.zeta(
        2, mpmath.mpf(100000) + mpmath.mpf(5) / 2) * 2
        - m * mpmath.mpf(1) / 2 * mpmath.zeta(
            4, mpmath.mpf(100000) + mpmath.mpf(5) / 2)))
    assert d < mpmath.mpf(10) ** -25, ("A_closed check m=%d: %s" % (m, d))
print("A(m) closed form checked against direct partial sums: PASS")


def anchor_term(m):
    L = 2 * A_closed(m) - P - phi_f(m)
    return phi_f(m) * L * L


print("(B) factorized anchor: Euler-Maclaurin tail + direct head ...")
t0 = time.time()
M0 = 200
head = mpmath.fsum(anchor_term(m) for m in range(1, M0))
tail = mpmath.sumem(anchor_term, [M0, mpmath.inf])
S3_anchor_em = head + tail
print("  sumem:  S3 = %s  (%.0fs)" % (mpmath.nstr(S3_anchor_em, 70),
                                      time.time() - t0))
t0 = time.time()
S3_anchor_lv = mpmath.nsum(anchor_term, [1, mpmath.inf], method='levin')
print("  levin:  S3 = %s  (%.0fs)" % (mpmath.nstr(S3_anchor_lv, 70),
                                      time.time() - t0))
agree = abs(S3_anchor_em - S3_anchor_lv)
print("  |sumem - levin| = %s" % mpmath.nstr(agree, 5))
OUT["S3_anchor_sumem"] = mpmath.nstr(S3_anchor_em, 110)
OUT["S3_anchor_levin"] = mpmath.nstr(S3_anchor_lv, 110)
OUT["anchor_method_agreement"] = mpmath.nstr(agree, 5)

res = abs(S3_dec - S3_anchor_em)
print("\nDECOMPOSITION RESIDUAL |A - B| = %s" % mpmath.nstr(res, 5))
OUT["decomposition_residual"] = mpmath.nstr(res, 5)
OUT["elapsed_sec"] = time.time() - t_start

outp = HERE / "data" / "s3_decomp_numerics.json"
outp.write_text(json.dumps(OUT, indent=2))
print("Saved:", outp)
