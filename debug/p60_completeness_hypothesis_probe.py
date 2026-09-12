"""Is overcompleteness ELIMINABLE?  The decisive test is hypothesis V1.

Proposition A (frames scan, 2026-09-11) says: if g != 0 lies in the closed span
of {f_i}, then lam_min(G_N) <= dist(g, V_N)^2 -> 0.  Completeness of the
ONE-CENTRE set alone forces the Gram to degenerate; the second centre is
incidental.  If that hypothesis holds here, "overcompleteness cannot be
eliminated" is a THEOREM for this basis, not a measurement.

The scan flagged V1 as the one load-bearing thing it could not verify: Coulomb
Sturmians are complete in the 1/r_A-weighted space, and L^2(1/r_A) and L^2(1/r_B)
are NOT the same space, so completeness in the MOLECULAR (Shibuya-Wulfman)
metric does not follow formally.

It is directly measurable, and cheaply, because of a fact Paper 60 already
proves: in the SW metric the intra-centre block is EXACTLY the identity.  So the
A-set is orthonormal, the projection of a displaced Sturmian chi^B_1 onto
span{chi^A_1..chi^A_N} has coefficients C_{i1}, and by Bessel

    eps_N^2  =  1 - sum_{i=1..N} C_{i1}^2   >=  0,

with eps_N -> 0 IFF chi^B_1 lies in the closed span of the A-set.  Prop A then
predicts lam_min <= eps_N^2, and the measured lam_min ~ N^-2 predicts
eps_N ~ N^-1.  Three things settle at once: the hypothesis, Prop A's
applicability, and whether the N^-2 law is explained by it.
"""
import numpy as np
from geovac.sturmian_sigma_law import sw_cross_block

M_QUAD = 400_001


def bessel_deficit(s, N):
    """eps_N^2 = 1 - sum_i <chi^A_i, chi^B_1>^2, in the SW metric."""
    C = sw_cross_block(s, N, M=M_QUAD)
    return float(1.0 - np.sum(C[:, 0] ** 2))


print("eps_N^2 = 1 - sum_i C_{i1}^2   (0 <=> the displaced Sturmian IS in the A-span)")
print(f"  {'N':>5}", "".join(f"{f'kR={s}':>16}" for s in (1.0, 2.0, 4.0)))
Ns = (4, 8, 16, 32, 64, 128, 256)
data = {s: [] for s in (1.0, 2.0, 4.0)}
for N in Ns:
    row = ""
    for s in (1.0, 2.0, 4.0):
        d = bessel_deficit(s, N)
        data[s].append(d)
        row += f"{d:16.3e}"
    print(f"  {N:5d}{row}")

print("\nrates (eps_N^2 ~ N^p, and eps_N ~ N^(p/2)):")
for s in (1.0, 2.0, 4.0):
    y = np.array(data[s])
    m = y > 1e-13
    p = np.polyfit(np.log(np.array(Ns)[m]), np.log(y[m]), 1)[0]
    loc = np.log(y[-1] / y[-2]) / np.log(Ns[-1] / Ns[-2])
    print(f"  kR={s}:  global p = {p:+.3f},  local p (last pair) = {loc:+.3f}"
          f"   => eps_N ~ N^{loc/2:+.3f}")

print("\nProp A's bound, against the measured lam_min:")
print(f"  {'N':>5} {'eps_N^2 (bound)':>17} {'lam_min = 1-sigma_max':>23} {'bound holds?':>13}")
for N in (8, 32, 128):
    C = sw_cross_block(2.0, N, M=M_QUAD)
    lam = float(1.0 - np.linalg.svd(C, compute_uv=False)[0])
    eps2 = bessel_deficit(2.0, N)
    print(f"  {N:5d} {eps2:17.3e} {lam:23.3e} {'YES' if lam <= eps2 + 1e-12 else 'NO':>13}")
