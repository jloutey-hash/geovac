"""End-to-end pricing of the preconditioner lever: does the QSVT DEGREE win
survive once subnormalization is counted?

Paper 60's model: whiten to X^dag W X, apply the metric factor by QSVT with a
polynomial approximating x^{-1/2}; the metric multiplies the query count by
(1 + d_inv), d_inv ~ kappa ln(kappa/eps).  The lever drops kappa from ~1e4 to
2.23.  But a block-encoding also carries a SUBNORMALIZATION alpha, and composing
block-encodings MULTIPLIES alphas, which can be far larger than the norm of the
product.  That is the failure mode that killed an earlier exponent-changing
variant in this same paper, so it must be checked and not assumed.

THE INVARIANCE THAT DECIDES IT.  Any X with X^dag A X = I satisfies
X = A^{-1/2} U for some unitary U, hence ||X|| = ||A^{-1/2}|| EXACTLY --
independent of how X is factored.  So no factorization can lower the amplitude
floor; the preconditioner can only buy DEPTH.  Checked numerically below rather
than asserted.
"""
import numpy as np
from geovac.sturmian_sigma_law import sw_cross_block

KR = 2.0
EPS = 1.6e-3


def tridiag(n):
    return (np.diag(2.0 * np.ones(n)) + np.diag(np.ones(n - 1), 1)
            + np.diag(np.ones(n - 1), -1))


def inv_sqrt(M):
    ev, U = np.linalg.eigh(M)
    return U @ np.diag(ev ** -0.5) @ U.T


def d_inv(kappa):
    """QSVT degree for x^{-1/2} on [1/kappa, 1] (Paper 60's model)."""
    return kappa * np.log(kappa / EPS)


print("A. the amplitude floor is factorization-invariant")
print(f"  {'n':>5} {'||A^-1/2||':>12} {'||P^-1/2 G^-1/2||':>19} {'rel diff':>10}")
for n in (20, 40, 80, 160):
    A = np.eye(n) - sw_cross_block(KR, n, M=100_001)
    P_is = inv_sqrt(tridiag(n))
    G = P_is @ A @ P_is
    X = P_is @ inv_sqrt(G)
    a, b = np.linalg.norm(inv_sqrt(A), 2), np.linalg.norm(X, 2)
    print(f"  {n:5d} {a:12.4f} {b:19.4f} {abs(a-b)/a:10.2e}")

print("\nB. what each route actually costs")
print(f"  {'n':>5} {'kappa(A)':>10} {'d_inv naive':>12} {'alpha naive':>12}"
      f" {'kappa(G)':>9} {'d_inv pre':>10} {'alpha pre':>11} {'net':>8}")
for n in (20, 40, 80, 160):
    A = np.eye(n) - sw_cross_block(KR, n, M=100_001)
    P = tridiag(n)
    P_is = inv_sqrt(P)
    G = P_is @ A @ P_is

    kA, kG = np.linalg.cond(A), np.linalg.cond(G)
    # naive: block-encode A (alpha ~ ||A||), QSVT to A^{-1/2}; the achieved
    # amplitude is the invariant floor ||A^{-1/2}||
    alpha_naive = np.linalg.norm(inv_sqrt(A), 2)
    # preconditioned by COMPOSITION: alpha_G = alpha_P^2 * alpha_A, then G^{-1/2}
    alpha_P = np.linalg.norm(P_is, 2)                    # = lam_min(P)^{-1/2}
    alpha_A = np.linalg.norm(A, 2)
    alpha_pre = alpha_P * np.linalg.norm(inv_sqrt(G), 2) * 1.0   # ||X|| achieved
    alpha_G_composed = alpha_P ** 2 * alpha_A
    dn, dp = d_inv(kA), d_inv(kG)
    net = (alpha_naive * dn) / (alpha_G_composed * dp)
    print(f"  {n:5d} {kA:10.1f} {dn:12.0f} {alpha_naive:12.1f}"
          f" {kG:9.3f} {dp:10.1f} {alpha_G_composed:11.1f} {net:8.0f}x")

print("\nC. the composition penalty, stated plainly")
n = 160
A = np.eye(n) - sw_cross_block(KR, n, M=100_001)
P_is = inv_sqrt(tridiag(n))
G = P_is @ A @ P_is
alpha_P = np.linalg.norm(P_is, 2)
print(f"  n = {n}")
print(f"    ||P^-1/2||            = {alpha_P:10.2f}   (analytic (n+1)/pi = {(n+1)/np.pi:.2f})")
print(f"    ||G||                 = {np.linalg.norm(G, 2):10.4f}")
print(f"    alpha_G if COMPOSED   = {alpha_P**2 * np.linalg.norm(A, 2):10.1f}"
      f"   <- {alpha_P**2 * np.linalg.norm(A,2)/np.linalg.norm(G,2):.0f}x larger than ||G||")
print(f"    amplitude floor ||X|| = {np.linalg.norm(P_is @ inv_sqrt(G), 2):10.2f}"
      f"   (invariant; the naive route already achieves it)")
print("\n  => the lever buys DEPTH, not amplitude; and a naive COMPOSED encoding of")
print("     G overpays amplitude by ~alpha_P^2.  A direct block-encoding of G is")
print("     the open item gating the end-to-end claim.")
