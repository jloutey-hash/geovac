"""Corrected locality measurement for the preconditioner probe.

The first pass fitted exp(-d/L) to an off-diagonal profile that the symbol
analysis says is ALGEBRAIC (|c_j| ~ j^-5/4 from the UV chirp).  Fitting an
exponential to a power law returns an L that scales with the fit window, hence
with n -- so that measurement could not distinguish the two cases and its
apparent "decay length grows like 0.14 n" was an artifact, not a result.

The operationally meaningful question is: to reach a fixed RELATIVE accuracy,
how does the required bandwidth scale with n?

  S^{-1/2} (ungerade):  symbol (1-a)^{-1/2} blows up at chi = pi and is not in
                        L^1.  The singularity SHARPENS as n grows.
                        Predict: bandwidth GROWS with n.
  G^{-1/2}:             symbol is bounded and smooth at chi = pi (the ratio
                        (1-a)/g tends to (kR)^2/24 != 0) and carries only the
                        chi -> 0 chirp, which does not depend on n.
                        Predict: bandwidth n-INDEPENDENT.
"""
import numpy as np
from p60_preconditioner_probe import cross_block, toeplitz_minus_hankel, inv_sqrt

kR = 2.0
PCOEF = {0: 2.0, 1: 1.0}


def bandwidth_for(M, tol):
    """Smallest b with ||M - band_b(M)||_F <= tol * ||M||_F."""
    n = M.shape[0]
    total = np.linalg.norm(M)
    idx = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    for b in range(n):
        if np.linalg.norm(M * (idx > b)) <= tol * total:
            return b
    return n


def profile_exponent(M):
    """Power-law exponent of the band-averaged off-diagonal profile."""
    n = M.shape[0]
    d = np.arange(1, n // 2)
    prof = np.array([np.abs(np.diag(M, k)).mean() for k in d])
    m = (d >= 3) & (prof > 0)
    return np.polyfit(np.log(d[m]), np.log(prof[m]), 1)[0]


print("Bandwidth needed for a fixed RELATIVE Frobenius accuracy")
for tol in (1e-2, 1e-3):
    print(f"\n  tol = {tol:.0e}")
    print(f"  {'n':>5} {'b(S^-1/2 unger)':>17} {'b/n':>7} {'b(G^-1/2)':>11} {'b/n':>7}"
          f" {'b(P^-1/2)':>11}")
    for n in (32, 64, 128, 256):
        C = cross_block(n, kR)
        A = np.eye(n) - C
        P = toeplitz_minus_hankel(PCOEF, n)
        Pis = inv_sqrt(P)
        G = Pis @ A @ Pis
        bS, bG, bP = (bandwidth_for(inv_sqrt(A), tol),
                      bandwidth_for(inv_sqrt(G), tol),
                      bandwidth_for(Pis, tol))
        print(f"  {n:5d} {bS:17d} {bS/n:7.3f} {bG:11d} {bG/n:7.3f} {bP:11d}")

print("\nPower-law exponent of the off-diagonal profile (chirp predicts ~ -5/4 where"
      "\nthe chi=pi singularity is absent):")
print(f"  {'n':>5} {'p(S^-1/2 unger)':>17} {'p(G^-1/2)':>12}")
for n in (64, 128, 256):
    C = cross_block(n, kR)
    A = np.eye(n) - C
    P = toeplitz_minus_hankel(PCOEF, n)
    Pis = inv_sqrt(P)
    G = Pis @ A @ Pis
    print(f"  {n:5d} {profile_exponent(inv_sqrt(A)):17.3f} {profile_exponent(inv_sqrt(G)):12.3f}")
