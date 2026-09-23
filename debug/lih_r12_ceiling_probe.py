"""LiH R12 F12-SYNTHESIS CEILING PROBE (diagnostic-before-engineering, 2026-09-22).

Question the diagnostic answers: how much of the Route-C residual can an r12 geminal
recover before the expensive two-center 4e mpf build?  Route C plateaus at E=-8.012
(58 mHa above exact) at the FCI determinant wall; the dominant residual is the Li 1s^2
core-core correlation, which the single analytic core orbital leaves at ZERO (so the
geminal's core contribution is DOUBLE-COUNTING-FREE).  The isolated Li+ (He-like Z=3)
core is the clean, transferable proxy for that pair.

Model: {Phi0, b_k Phi0} generalized eigenproblem, Phi0 = 1s(zeta)^2 with zeta=2.6875
(= Route C's frozen analytic core exponent -- NOT re-optimized, the realistic scenario),
in Hylleraas coordinates (r1, r2, r12); volume element 8 pi^2 r1 r2 r12; kinetic energy
via the symmetric gradient form (first derivatives only; the corpus 'bounded gradient').
Hylleraas primitives b(r1,r2,r12) in {1, u=r12, u2, t2=(r1-r2)^2, s=r1+r2, ut2}.

Finding (see debug/sprint_lih_r12_ceiling_diagnostic_memo.md):
  - r12-only {u,u2}      -> ~52% of the core deficit (PLATEAUS; more u-terms do not help)
  - radial-only {t2,s}   -> ~45-51%  (an ORBITAL degree of freedom, NOT a geminal one)
  - both {u,t2}          -> ~74-79%  (near-additive: the two halves are complementary)
  - rich {u,t2,s,u2,ut2} -> ~95-98%  (full core correlation with a small correlated basis)
So a geminal ALONE caps at ~half the core correlation; the other half needs a 2nd core
ORBITAL (radial flexibility), which Route C's ladder can supply with existing code.

Validated against known He (Z=2) Hylleraas: {u,t2}=-2.892 (classic 3-term -2.9024),
5-term -> 98% of the 57 mHa He correlation.  Float64, ~seconds.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

N1 = N2 = 72; NU = 56; RMAX = 7.0
_x1, _w1 = leggauss(N1)
r1 = 0.5 * RMAX * (_x1 + 1.0); wr1 = 0.5 * RMAX * _w1
r2 = r1.copy(); wr2 = wr1.copy()
_xu, _wu = leggauss(NU)
R1 = r1[:, None, None]; R2 = r2[None, :, None]
ulo = np.abs(R1 - R2); uhi = R1 + R2
R12 = ulo + (uhi - ulo) * 0.5 * (_xu[None, None, :] + 1.0)
WVOL = (8 * np.pi**2 * (wr1[:, None, None] * R1) * (wr2[None, :, None] * R2)
        * (0.5 * (uhi - ulo) * _wu[None, None, :]) * R12)
C1 = (R1**2 - R2**2 + R12**2) / (2 * R1 * R12)     # rhat1 . (r1-r2)/r12
C2 = (R2**2 - R1**2 + R12**2) / (2 * R2 * R12)     # rhat2 . (r2-r1)/r12


def integ(g):
    return np.sum(WVOL * g)


def basisfn(kind):
    """return (b, db/dr1, db/dr2, db/dr12) arrays for a Hylleraas primitive."""
    O = np.ones_like(R12); Z = np.zeros_like(R12)
    if kind == '1':   return O, Z, Z, Z
    if kind == 'u':   return R12, Z, Z, O
    if kind == 'u2':  return R12**2, Z, Z, 2 * R12
    if kind == 't2':  return (R1 - R2)**2, 2 * (R1 - R2) * O, -2 * (R1 - R2) * O, Z
    if kind == 's':   return (R1 + R2) * O, O, O, Z
    if kind == 'ut2': return R12 * (R1 - R2)**2, 2 * (R1 - R2) * R12, -2 * (R1 - R2) * R12, (R1 - R2)**2
    raise ValueError(kind)


def energy(Z, zeta, kinds):
    N = zeta**3 / np.pi
    Phi2 = (N * np.exp(-zeta * (R1 + R2)))**2
    V = -Z / R1 - Z / R2 + 1.0 / R12
    B = [basisfn(k) for k in (['1'] + kinds)]
    n = len(B)
    S = np.zeros((n, n)); T = np.zeros((n, n)); Vm = np.zeros((n, n))
    for i in range(n):
        bi, di1, di2, diu = B[i]
        for j in range(i, n):
            bj, dj1, dj2, dju = B[j]
            S[i, j] = S[j, i] = integ(Phi2 * bi * bj)
            Vm[i, j] = Vm[j, i] = integ(Phi2 * bi * bj * V)
            # electron 1: (di1 - z bi)(dj1 - z bj) + diu dju + [(di1-z bi)dju+(dj1-z bj)diu] C1
            a_i1 = di1 - zeta * bi; a_j1 = dj1 - zeta * bj
            a_i2 = di2 - zeta * bi; a_j2 = dj2 - zeta * bj
            e1 = a_i1 * a_j1 + diu * dju + (a_i1 * dju + a_j1 * diu) * C1
            e2 = a_i2 * a_j2 + diu * dju + (a_i2 * dju + a_j2 * diu) * C2
            T[i, j] = T[j, i] = 0.5 * integ(Phi2 * (e1 + e2))
    H = T + Vm
    sval, svec = np.linalg.eigh(S)
    Xs = svec @ np.diag(1 / np.sqrt(np.clip(sval, 1e-13, None))) @ svec.T
    w = np.linalg.eigvalsh(Xs @ H @ Xs)
    return H[0, 0] / S[0, 0], w[0], np.linalg.cond(S)


# ---- quadrature self-check (norm=1, <1/r1>=zeta, T00=zeta^2) --------------- #
def _validate(zeta=2.6875, Z=3.0):
    N = zeta**3 / np.pi
    Phi2 = (N * np.exp(-zeta * (R1 + R2)))**2
    norm = integ(Phi2)
    inv_r1 = integ(Phi2 / R1) / norm
    _, _, _ = energy(Z, zeta, ['u'])
    print("quadrature self-check (zeta=2.6875):")
    print(f"  norm={norm:.8f} (1)   <1/r1>={inv_r1:.6f} (zeta={zeta})   "
          f"E0={energy(Z, zeta, ['u'])[0]:.6f} (zeta^2-2Z zeta+5zeta/8={zeta**2-2*Z*zeta+5*zeta/8:.6f})\n")


CASES = {'He (Z=2)': (2.0, 27/16, -2.90372), 'Li+ (Z=3)': (3.0, 2.6875, -7.27991)}
SETS = [
    ("r12 only        {u}",          ['u']),
    ("r12 only        {u,u2}",       ['u', 'u2']),
    ("radial only     {t2}",         ['t2']),
    ("radial only     {t2,s}",       ['t2', 's']),
    ("BOTH  {u,t2}",                 ['u', 't2']),
    ("BOTH  {u,t2,s,u2,ut2}",        ['u', 't2', 's', 'u2', 'ut2']),
]
if __name__ == "__main__":
    _validate()
    for name, (Z, zeta, exact) in CASES.items():
        E0, _, _ = energy(Z, zeta, ['u'])
        corr = E0 - exact
        print("=" * 76)
        print(f"{name}: ref E0={E0:.5f} exact={exact} corr={corr*1e3:.1f} mHa")
        print("=" * 76)
        for label, ks in SETS:
            _, ER, cond = energy(Z, zeta, ks)
            print(f"  {label:26s} E={ER:.5f}  rec={ (E0-ER)*1e3:6.1f} mHa  {100*(E0-ER)/corr:5.0f}%  cond={cond:.0e}")
        print()
