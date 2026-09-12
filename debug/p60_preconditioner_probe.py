"""Does the Toeplitz-preconditioner escape survive independent verification?

The 2026-09-11 Toeplitz literature scan reported that Serra's theorem
(Math. Comp. 66 (1997) 651) applies to Paper 60's two-centre metric: a
trigonometric polynomial matching the symbol's zero gives a preconditioned
spectrum inside a fixed interval for EVERY n.  It also reported a constructive
factorization in which all the non-locality sits in P^{-1/2} for P tridiagonal
(an exact DST), while G = P^{-1/2}(I-C)P^{-1/2} is well conditioned with
n-independent decay.

The agent reproducing its own number is not an independent route.  This is the
independent route: everything below is built from the verified symbol routine
of debug/p60_symbol_pole_decay.py, with no number taken on report.

Structure of the check:
  1. cond(I -+ C): ill-conditioning confined to the UNGERADE block (paper's claim).
  2. P = Toeplitz(2 + 2cos chi) has a quadratic zero at chi = pi, exactly where
     1 - a(chi) does, and in this basis the Hankel part vanishes identically,
     so P is EXACTLY tridiagonal(1, 2, 1).
  3. P is diagonalized exactly by the DST-I -- checked against the closed-form
     eigenpairs, not against a numerical eigensolver.
  4. cond(G) flat in n?  This is the load-bearing claim.
  5. Locality: S^{-1/2} decay length vs n (should DIVERGE, per the holomorphic-
     calculus obstruction) against G's (should be n-independent).
  6. What it would buy on Paper 60's own resource model, d_inv ~ kappa ln(kappa/eps).
"""
import numpy as np

GL_N = 48
_xg, _wg = np.polynomial.legendre.leggauss(GL_N)
_cache = {}


def c_j(j, kR, Smax=2000.0):
    """Cosine coefficient of a(chi) = j0(kR cot(chi/2)); see p60_symbol_pole_decay."""
    key = (j, kR)
    if key in _cache:
        return _cache[key]
    edges, s = [0.0], 0.0
    while s < Smax:
        s += min(np.pi / 2, np.pi / max(2 * j * kR / (kR * kR + s * s), 1e-300))
        edges.append(min(s, Smax))
    e = np.asarray(edges)
    mid, half = 0.5 * (e[:-1] + e[1:]), 0.5 * (e[1:] - e[:-1])
    sv = (mid[:, None] + half[:, None] * _xg[None, :]).ravel()
    wv = (half[:, None] * _wg[None, :]).ravel()
    g = np.cos(2 * j * np.arctan2(kR, sv)) * np.sinc(sv / np.pi) \
        * 2 * kR / (kR * kR + sv * sv)
    _cache[key] = float(np.dot(wv, g) / np.pi)
    return _cache[key]


def cross_block(n, kR):
    """C_{nm} = c_{n-m} - c_{n+m}: Toeplitz MINUS Hankel, not Toeplitz."""
    c = {j: c_j(j, kR) for j in range(0, 2 * n + 2)}
    return np.array([[c[abs(a - b)] - c[a + b] for b in range(1, n + 1)]
                     for a in range(1, n + 1)])


def toeplitz_minus_hankel(coeffs, n):
    """Same basis convention, for an arbitrary cosine-coefficient dict."""
    g = lambda j: coeffs.get(j, 0.0)
    return np.array([[g(abs(a - b)) - g(a + b) for b in range(1, n + 1)]
                     for a in range(1, n + 1)])


def inv_sqrt(M):
    ev, U = np.linalg.eigh(M)
    return U @ np.diag(ev ** -0.5) @ U.T


def decay_length(M):
    """Fit |M_{i,i+d}| ~ exp(-d/L) over the band-averaged profile; large L = delocal."""
    n = M.shape[0]
    prof = np.array([np.abs(np.diag(M, d)).mean() for d in range(1, n // 2)])
    prof = np.maximum(prof, 1e-300)
    d = np.arange(1, len(prof) + 1)
    m = prof > prof[0] * 1e-12
    sl = np.polyfit(d[m], np.log(prof[m]), 1)[0]
    return -1.0 / sl if sl < 0 else np.inf


kR = 2.0
# Preconditioner symbol g(chi) = 2 + 2 cos chi: zero at chi = pi, quadratic.
# Cosine coefficients: g_0 = 2, g_1 = 1.  n+m >= 2 always, so the Hankel part
# hits only g_{>=2} = 0 -> P is exactly tridiagonal.
PCOEF = {0: 2.0, 1: 1.0}

print("1-2. conditioning by parity, and the preconditioned ungerade block")
print(f"  {'n':>5} {'cond(I+C)':>12} {'cond(I-C)':>14} {'cond(G)':>10}  {'G symbol ratio':>14}")
rows = []
for n in (10, 20, 40, 80, 160):
    C = cross_block(n, kR)
    I = np.eye(n)
    P = toeplitz_minus_hankel(PCOEF, n)
    assert np.allclose(P, np.triu(np.tril(P, 1), -1)), "P is not tridiagonal"
    Pinv_sqrt = inv_sqrt(P)
    G = Pinv_sqrt @ (I - C) @ Pinv_sqrt
    rows.append((n, np.linalg.cond(I + C), np.linalg.cond(I - C), np.linalg.cond(G)))
    print(f"  {n:5d} {rows[-1][1]:12.4f} {rows[-1][2]:14.1f} {rows[-1][3]:10.4f}"
          f"  {(kR**2/24):14.5f}")

print("\n3. is P exactly DST-I diagonalizable?  (closed form, not an eigensolver)")
n = 64
P = toeplitz_minus_hankel(PCOEF, n)
k = np.arange(1, n + 1)
V = np.sqrt(2 / (n + 1)) * np.sin(np.outer(k, k) * np.pi / (n + 1))   # DST-I
lam_closed = 2 + 2 * np.cos(k * np.pi / (n + 1))
resid = np.abs(P @ V - V * lam_closed).max()
print(f"   max |P V - V diag(lambda)| = {resid:.3e}   (DST-I basis, closed-form lambda)")
print(f"   V orthogonal to {np.abs(V @ V.T - np.eye(n)).max():.3e}")

print("\n5. locality: decay length vs n  (the holomorphic obstruction vs the repair)")
print(f"  {'n':>5} {'L(S^-1/2 ungerade)':>20} {'L(G^-1/2)':>12} {'L(P^-1/2)':>12}")
for n in (32, 64, 128):
    C = cross_block(n, kR)
    A = np.eye(n) - C
    P = toeplitz_minus_hankel(PCOEF, n)
    Pis = inv_sqrt(P)
    G = Pis @ A @ Pis
    print(f"  {n:5d} {decay_length(inv_sqrt(A)):20.2f} {decay_length(inv_sqrt(G)):12.2f}"
          f" {decay_length(Pis):12.2f}")

print("\n6. what it would buy on Paper 60's own resource model, d_inv ~ kappa ln(kappa/eps)")
eps = 1.6e-3
for n, cg, cu, cG in rows:
    d_naive = cu * np.log(cu / eps)
    d_pre = cG * np.log(cG / eps)
    print(f"  n={n:4d}  kappa(I-C)={cu:12.1f} -> d_inv={d_naive:12.1f}    "
          f"kappa(G)={cG:6.3f} -> d_inv={d_pre:7.1f}   ratio {d_naive/d_pre:9.1f}x")
