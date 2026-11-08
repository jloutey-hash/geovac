"""S^(3) anchor discriminator.

The two candidate values from s3_decomp_numerics.py disagree:
    (A) decomposition route : 30.61538052881950...
    (B) nsum-levin anchor   : 30.21972625247710...   [suspect]
All anchor terms phi(m) L(m)^2 are POSITIVE, so the partial sums S(M)
increase monotonically to S3. If S(M) ever exceeds (B), (B) is refuted.

Method: exact recurrence  T(m) = T(m-1) - phi(m-1),  A(m) = A(m-1) + T(m)
(A(1) = T(1) = P), term = phi(m) (2A(m) - P - phi(m))^2, summed in double
with Kahan compensation and re-anchored from closed forms (mpmath, 60 dps)
every 200k steps. Then a log-polynomial Richardson fit
    S(M) = S_inf - (a ln^2 M + b ln M + c)/M - (d ln^3 M)/M^2 ...
across checkpoints estimates S_inf to ~8-12 digits.

Also: positive-term truncation windows for the 12 triple t-values
(partial sum at O=2*10^5 plus crude tail bound) to catch any gross
Levin error in the per-symbol evaluations of route (A).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import mpmath
import numpy as np

mpmath.mp.dps = 60
HERE = Path(__file__).parent
OUT = {}

pi = mpmath.pi
P_mp = pi ** 2 - pi ** 4 / 12 - mpmath.mpf(64) / 81
P = float(P_mp)


def phi(m: float) -> float:
    x = m + 1.5
    return 2.0 / (x * x) - 0.5 / (x ** 4)


def T_closed(m: int) -> mpmath.mpf:
    a = mpmath.mpf(m) + mpmath.mpf(3) / 2
    return 2 * mpmath.zeta(2, a) - mpmath.mpf(1) / 2 * mpmath.zeta(4, a)


def A_closed(m: int) -> mpmath.mpf:
    mm = mpmath.mpf(m)
    z2m = mpmath.zeta(2, mm + mpmath.mpf(5) / 2)
    z3m = mpmath.zeta(3, mm + mpmath.mpf(5) / 2)
    z4m = mpmath.zeta(4, mm + mpmath.mpf(5) / 2)
    hz = lambda s: mpmath.zeta(s, mpmath.mpf(5) / 2)
    return (2 * (mpmath.digamma(mm + mpmath.mpf(5) / 2)
                 - mpmath.digamma(mpmath.mpf(5) / 2))
            - 3 * (hz(2) - z2m) - mpmath.mpf(1) / 2 * (hz(3) - z3m)
            + mpmath.mpf(3) / 4 * (hz(4) - z4m)
            + mm * (2 * z2m - mpmath.mpf(1) / 2 * z4m))


print("Monotone partial sums of the anchor (positive terms) ...")
t0 = time.time()
M_END = 6_400_000
CHECK = [100_000, 200_000, 400_000, 800_000, 1_600_000, 3_200_000, 6_400_000]
REANCHOR = 200_000

S = 0.0
comp = 0.0  # Kahan
m = 1
T = float(T_closed(1))
A = float(A_closed(1))
checkpoints = {}
while m <= M_END:
    L = 2.0 * A - P - phi(m)
    y = phi(m) * L * L - comp
    t = S + y
    comp = (t - S) - y
    S = t
    if m in CHECK or m == M_END:
        checkpoints[m] = S
        print("  S(%9d) = %.14f   (%.0fs)" % (m, S, time.time() - t0))
    m += 1
    if m % REANCHOR == 0:
        T = float(T_closed(m - 1))
        A = float(A_closed(m - 1))
    # recurrence step to m
    T = T - phi(m - 1)
    A = A + T

B_LEVIN = 30.21972625247710
refuted = any(v > B_LEVIN for v in checkpoints.values())
print("  levin anchor 30.2197... refuted by monotonicity: %s" % refuted)
OUT["partial_sums"] = {str(k): repr(v) for k, v in checkpoints.items()}
OUT["levin_anchor_refuted"] = bool(refuted)

# Richardson-style fit:  S(M) = S_inf - (a ln^2 M + b ln M + c)/M
#                                 - (d ln^3 M + e ln^2 M)/M^2
Ms = np.array(sorted(checkpoints), dtype=float)
Sv = np.array([checkpoints[int(k)] for k in sorted(checkpoints)])
lo = np.log(Ms)
X = np.column_stack([
    np.ones_like(Ms), -lo ** 2 / Ms, -lo / Ms, -1.0 / Ms,
    -lo ** 3 / Ms ** 2, -lo ** 2 / Ms ** 2,
])
coef, res_, rank_, sv_ = np.linalg.lstsq(X, Sv, rcond=None)
S_inf = coef[0]
fit_resid = float(np.max(np.abs(X @ coef - Sv)))
print("  log-Richardson fit: S_inf = %.12f  (max fit residual %.2e)"
      % (S_inf, fit_resid))
OUT["S_inf_fit"] = repr(S_inf)
OUT["fit_residual"] = fit_resid

A_DEC = 30.61538052881950170777397460662109
print("  decomposition route 30.615380528819... consistent: |diff| = %.2e"
      % abs(S_inf - A_DEC))
OUT["fit_vs_decomposition"] = repr(abs(S_inf - A_DEC))

# ------------------------------------------------------------------
# truncation windows for the 12 t3 values (positive terms)
# ------------------------------------------------------------------
print("\nTruncation windows for t3 (partial @ O=2e5 + tail bound) ...")
mpmath.mp.dps = 30
LAM = {s: (1 - mpmath.power(2, -s)) * mpmath.zeta(s) for s in range(2, 13)}
PSI0_HALF = mpmath.digamma(mpmath.mpf(1) / 2)

num = json.loads((HERE / "data" / "s3_decomp_numerics.json").read_text())
import ast
t3_levin = {ast.literal_eval(k): mpmath.mpf(v)
            for k, v in num["t3_values"].items()}

JMAX = 100_000
windows = {}
for (tag, b1, b2, b3), lev in sorted(t3_levin.items()):
    s = mpmath.mpf(0)
    # iterate over o2 = 2j+1
    for j in range(1, JMAX + 1):
        tau = mpmath.power(2, -b1) * mpmath.zeta(
            b1, mpmath.mpf(j) + mpmath.mpf(3) / 2)
        if b3 == 1:
            pL = (mpmath.digamma(mpmath.mpf(j) + mpmath.mpf(1) / 2)
                  - PSI0_HALF) / 2
        else:
            pL = LAM[b3] - mpmath.power(2, -b3) * mpmath.zeta(
                b3, mpmath.mpf(j) + mpmath.mpf(1) / 2)
        s += mpmath.mpf(2 * j + 1) ** -b2 * tau * pL
    # tail bound: term(j) <= (2j+1)^-b2 * tau_max * pL_max with
    # tau(j) <= 2^-b1 * ((b1-1)^-1 (j+1/2)^(1-b1));  pL <= lam(b3) (b3>=2)
    # or pL <= ln(2*JMAX+3) (b3=1, crude).
    o = mpmath.mpf(2 * JMAX + 3)
    pLmax = (mpmath.log(o) if b3 == 1 else LAM[b3])
    # sum_{o2 > O} o2^-b2 * tau_{b1}(o2) <= pLmax * integral bound
    expo = b1 + b2 - 1
    tail = pLmax * mpmath.mpf(2) / (expo - 1) / o ** (expo - 1)
    lo_, hi_ = s, s + tail
    ok = (lo_ <= lev <= hi_)
    windows[str((b1, b2, b3))] = {
        "partial": mpmath.nstr(lo_, 20), "tail_bound": mpmath.nstr(tail, 5),
        "levin": mpmath.nstr(lev, 20), "inside": bool(ok)}
    print("  t3(%d,%d,%d): partial %s  tail<%s  levin %s  %s"
          % (b1, b2, b3, mpmath.nstr(lo_, 12), mpmath.nstr(tail, 3),
             mpmath.nstr(lev, 12), "OK" if ok else "OUTSIDE"))

OUT["t3_windows"] = windows
OUT["all_t3_inside"] = all(w["inside"] for w in windows.values())
print("all t3 inside windows:", OUT["all_t3_inside"])

outp = HERE / "data" / "s3_anchor_check.json"
outp.write_text(json.dumps(OUT, indent=2))
print("Saved:", outp)
