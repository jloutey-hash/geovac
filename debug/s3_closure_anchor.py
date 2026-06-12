"""S^(3) closure sprint — follow-on 2: global anchor via rigorous bracket.

The stage-1 anchor attempts failed for the SAME reason as the trailing-1
t3 evaluator: the factorized summand phi(m) L(m)^2 carries (ln m)^2 growth
through the digamma in A(m), and both mpmath.sumem (diverged, -1.8e+403)
and nsum-levin (30.2197, wrong) silently mishandle log-modulated algebraic
decay.  This driver uses NO series acceleration anywhere:

  S3 = sum_{m>=1} phi(m) L(m)^2,         phi(m) = 8((2m+3)^2-1)/(2m+3)^4,
  L(m) = sum_{n>=1} I(n,m) phi(n) = 2 A(m) - P - phi(m),
  A(m) = sum_{n<m} n phi(n) + m T(m),    T(m) = P - sum_{n<m} phi(n),
  P    = sum_{n>=1} phi(n) = pi^2 - pi^4/12 - 64/81.

All terms positive, so the partial sum S_head(N) is a rigorous LOWER
bound, and an explicit integral bound gives a rigorous UPPER bound:

  L(m) <= 2 A(m) <= 4 (ln m + 1)        [n phi(n) < 2/(n+1), H-bound;
                                          m T(m) < 2 by integral comparison]
  phi(m) < 2/(m+1)^2
  tail(N) < 32 * [(ln N + 1)^2 + 2(ln N + 1) + 2] / N
            [closed-form integral of (ln x + 1)^2 / x^2, integrand
             decreasing for x > 1].

Gates:
  GA (closed form): incremental A(m) vs the digamma/Hurwitz closed form
     A_closed(m) from s3_decomp_numerics.py at m in {1,2,5,100,1000}.
  GB (exact rational): the cube-truncated triple sum at N=300 computed in
     exact Fractions equals the mpf cube sum to 1e-45 — this IS the
     memo's "exact-rational partial triple sum at moderate N".
  GC (bound validity): L(m) <= 4(ln m + 1) and phi(m) <= 2/(m+1)^2
     verified numerically on sampled m.
  VERDICT: bracket [S_head(N), S_head(N) + tail_bound(N)] must contain
     the decomposition value 30.61538... and exclude the bad Levin
     anchor 30.2197...

Output: debug/data/s3_closure_anchor.json
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

import mpmath

HERE = Path(__file__).parent
OUT = {}
t_start = time.time()

S3_DEC_STALE = mpmath.mpf("30.61538052881950")   # 130-dps run, stale t3(*,*,1)
S3_LEVIN_BAD = mpmath.mpf("30.2197")

mpmath.mp.dps = 60
pi = mpmath.pi
P_const = pi ** 2 - pi ** 4 / 12 - mpmath.mpf(64) / 81


def phi_mpf(m):
    x = mpmath.mpf(2 * m + 3)
    x2 = x * x
    return 8 * (x2 - 1) / (x2 * x2)


def phi_frac(n):
    o = 2 * n + 3
    return Fraction(8 * (o * o - 1), o ** 4)


# ------------------------------------------------------------------ GA
def A_closed(m):
    hz = lambda s: mpmath.zeta(s, mpmath.mpf(5) / 2)
    m_ = mpmath.mpf(m)
    z2m = mpmath.zeta(2, m_ + mpmath.mpf(5) / 2)
    z3m = mpmath.zeta(3, m_ + mpmath.mpf(5) / 2)
    z4m = mpmath.zeta(4, m_ + mpmath.mpf(5) / 2)
    return (2 * (mpmath.digamma(m_ + mpmath.mpf(5) / 2)
                 - mpmath.digamma(mpmath.mpf(5) / 2))
            - 3 * (hz(2) - z2m) - mpmath.mpf(1) / 2 * (hz(3) - z3m)
            + mpmath.mpf(3) / 4 * (hz(4) - z4m)
            + m_ * (2 * z2m - mpmath.mpf(1) / 2 * z4m))


def A_incremental_single(m):
    cum_phi = mpmath.mpf(0)
    cum_nphi = mpmath.mpf(0)
    for n in range(1, m):
        p = phi_mpf(n)
        cum_phi += p
        cum_nphi += n * p
    return cum_nphi + m * (P_const - cum_phi)


print("=== GA: incremental A(m) vs digamma/Hurwitz closed form ===")
ga = {}
for m in (1, 2, 5, 100, 1000):
    d = abs(A_incremental_single(m) - A_closed(m))
    ga["m=%d" % m] = mpmath.nstr(d, 4)
    print("  m=%d: |incr - closed| = %s" % (m, mpmath.nstr(d, 4)))
GA_PASS = all(mpmath.mpf(v) < mpmath.mpf(10) ** -45 for v in ga.values())
OUT["GA"] = ga
OUT["GA_pass"] = GA_PASS
print("  GA:", "PASS" if GA_PASS else "FAIL")

# ------------------------------------------------------------------ GB
print("=== GB: exact-Fraction cube sum at N=300 vs mpf ===")
N_CUBE = 300
t0 = time.time()
P_N = sum(phi_frac(n) for n in range(1, N_CUBE + 1))
cube_frac = Fraction(0)
cum_phi_f = Fraction(0)
cum_nphi_f = Fraction(0)
for m in range(1, N_CUBE + 1):
    pm = phi_frac(m)
    # A_N(m) = sum_{n<=N} min(n,m) phi(n) = cum_nphi(m) + m*(P_N - cum_phi(m))
    A_N = cum_nphi_f + m * (P_N - cum_phi_f)
    L_N = 2 * A_N - P_N - pm
    cube_frac += pm * L_N * L_N
    cum_phi_f += pm
    cum_nphi_f += m * pm

cube_mpf = mpmath.mpf(0)
P_N_m = mpmath.fsum(phi_mpf(n) for n in range(1, N_CUBE + 1))
cum_phi_m = mpmath.mpf(0)
cum_nphi_m = mpmath.mpf(0)
for m in range(1, N_CUBE + 1):
    pm = phi_mpf(m)
    A_N = cum_nphi_m + m * (P_N_m - cum_phi_m)
    L_N = 2 * A_N - P_N_m - pm
    cube_mpf += pm * L_N * L_N
    cum_phi_m += pm
    cum_nphi_m += m * pm

d_gb = abs(mpmath.mpf(cube_frac.numerator) / cube_frac.denominator - cube_mpf)
GB_PASS = d_gb < mpmath.mpf(10) ** -45
OUT["GB_cube_N"] = N_CUBE
OUT["GB_frac_vs_mpf"] = mpmath.nstr(d_gb, 4)
OUT["GB_cube_value"] = mpmath.nstr(cube_mpf, 30)
OUT["GB_pass"] = GB_PASS
print("  |Fraction - mpf| = %s (%.0fs)  GB: %s"
      % (mpmath.nstr(d_gb, 4), time.time() - t0,
         "PASS" if GB_PASS else "FAIL"))

# ------------------------------------------------------------------ head loop
print("=== Head partial sums (no acceleration) ===")
CHECKPOINTS = [10 ** 4, 10 ** 5, 10 ** 6, 4 * 10 ** 6]
mpmath.mp.dps = 30
pi30 = mpmath.pi
P30 = pi30 ** 2 - pi30 ** 4 / 12 - mpmath.mpf(64) / 81

head = mpmath.mpf(0)
cum_phi = mpmath.mpf(0)
cum_nphi = mpmath.mpf(0)
results = {}
gc_samples = {}
t0 = time.time()
m = 0
for cp in CHECKPOINTS:
    while m < cp:
        m += 1
        x = mpmath.mpf(2 * m + 3)
        x2 = x * x
        pm = 8 * (x2 - 1) / (x2 * x2)
        A = cum_nphi + m * (P30 - cum_phi)
        L = 2 * A - P30 - pm
        head += pm * L * L
        cum_phi += pm
        cum_nphi += m * pm
        if m in (10, 100, 10 ** 4, 10 ** 6):
            gc_samples[m] = (L, pm)
    # rigorous tail bound: 32*[(ln N + 1)^2 + 2(ln N + 1) + 2]/N
    u = mpmath.log(cp) + 1
    tb = 32 * (u * u + 2 * u + 2) / cp
    results[cp] = (head, tb)
    print("  N=%-9d head=%s  tail_bound=%s  (%.0fs)"
          % (cp, mpmath.nstr(head, 20), mpmath.nstr(tb, 4),
             time.time() - t0))

# ------------------------------------------------------------------ GC
print("=== GC: bound validity on sampled m ===")
gc = {}
for ms, (L, pm) in sorted(gc_samples.items()):
    okL = L <= 4 * (mpmath.log(ms) + 1)
    okP = pm <= mpmath.mpf(2) / (ms + 1) ** 2
    gc["m=%d" % ms] = bool(okL and okP)
    print("  m=%d: L<=4(ln m+1): %s   phi<=2/(m+1)^2: %s" % (ms, okL, okP))
GC_PASS = all(gc.values())
OUT["GC"] = gc
OUT["GC_pass"] = GC_PASS

# ------------------------------------------------------------------ verdict
print("=== VERDICT ===")
ver = {}
for cp, (h, tb) in results.items():
    lo, hi = h, h + tb
    in_dec = lo <= S3_DEC_STALE <= hi
    out_bad = not (lo <= S3_LEVIN_BAD <= hi)
    ver[str(cp)] = {
        "lower": mpmath.nstr(lo, 22), "upper": mpmath.nstr(hi, 22),
        "width": mpmath.nstr(tb, 4),
        "contains_decomposition_30.615": bool(in_dec),
        "excludes_levin_30.2197": bool(out_bad)}
    print("  N=%-9d  S3 in [%s, %s]   dec-in: %s   levin-out: %s"
          % (cp, mpmath.nstr(lo, 14), mpmath.nstr(hi, 14), in_dec, out_bad))
OUT["bracket"] = ver
best = ver[str(CHECKPOINTS[-1])]
OUT["anchor_closed"] = (GA_PASS and GB_PASS and GC_PASS
                        and best["contains_decomposition_30.615"]
                        and best["excludes_levin_30.2197"])
OUT["elapsed_sec"] = time.time() - t_start
(HERE / "data" / "s3_closure_anchor.json").write_text(
    json.dumps(OUT, indent=2))
print("\nANCHOR CLOSED:", OUT["anchor_closed"])
print("Saved:", HERE / "data" / "s3_closure_anchor.json")
