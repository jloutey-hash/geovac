"""H2 in prolate-spheroidal, re-based monomial -> Laguerre(xi) x Legendre(eta).

Strategy that avoids reimplementing the physics:
  * Work in the UNSYMMETRIZED product basis g(j,l,k,m; mu) -- one .terms entry
    each -- so geovac's own (validated) one_body + neumann_vee_general_m.vee_matrix
    build S, T+V_ne, V_ee with ACCURATE entries (V_ee is mpmath-seeded, dps=30).
    The physical H is 1<->2 symmetric, so its lowest eigenvalue over the full
    product space is the spatially-symmetric (singlet) ground state -- no
    hand symmetrization needed.
  * Re-base each electron's radial monomial xi^j -> L_j(2a(xi-1)) and angular
    eta^l -> P_l(eta). Same SPAN => identical exact-arithmetic energy; only the
    conditioning changes. C = T_xi (x) T_eta (x) T_xi (x) T_eta, block-diagonal
    in mu.
  * Solve H_orth c = E S_orth c in the well-conditioned basis.

Validate against the known monomial result at (3,3), mu<=1 (99.09% of D_e),
then push (j_max,l_max) and watch: climb toward 99.97% (geometry reaches
spectroscopic accuracy) or plateau (the e-e cusp / Layer-3 wall).
"""
import time
import numpy as np
import mpmath as mp

from geovac import prolate_general_m as pg
from geovac import neumann_vee_general_m as ngm

mp.mp.dps = 40
DE_EXACT = 0.174475
E_EXACT = -1.174475
R = pg.R_DEFAULT


# ---- one-electron re-basing transforms (exact coeffs, downcast to float64) ----
def _pmul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out

def laguerre_row(n, alpha, width):
    """coeff of xi^j in L_n(2a(xi-1)), j=0..width-1 (float64)."""
    s = 2 * alpha
    cx = [(mp.mpf(-1) ** k) * mp.binomial(n, k) / mp.factorial(k) for k in range(n + 1)]
    lin = [mp.mpf(-s), mp.mpf(s)]           # s*xi - s
    out = [mp.mpf(0)]
    xpow = [mp.mpf(1)]
    for k, ck in enumerate(cx):
        if k > 0:
            xpow = _pmul(xpow, lin)
        term = [ck * t for t in xpow]
        if len(term) > len(out):
            out += [mp.mpf(0)] * (len(term) - len(out))
        for i, t in enumerate(term):
            out[i] += t
    row = [0.0] * width
    for i in range(min(len(out), width)):
        row[i] = float(out[i])
    return row

def legendre_row(l, width):
    """coeff of eta^l' in P_l(eta) (float64)."""
    c = list(ngm._leg_coeffs(l))
    row = [0.0] * width
    for i in range(min(len(c), width)):
        row[i] = float(c[i])
    return row


class UnsymP:
    """single-term (unsymmetrized) product basis function."""
    __slots__ = ("j", "k", "l", "m", "mu", "alpha")
    def __init__(self, j, l, k, m, mu, alpha):
        self.j, self.l, self.k, self.m, self.mu, self.alpha = j, l, k, m, mu, alpha
    @property
    def terms(self):
        return [(self.j, self.l, self.k, self.m)]


def run(j_max, l_max, mu_max, alpha=1.0, l_neumann=16, verbose=False):
    t0 = time.time()
    # unsym product index set (gerade: l+m even), same set for mono & orth
    idx = [(j, l, k, m, mu)
           for mu in range(mu_max + 1)
           for j in range(j_max + 1)
           for l in range(l_max + 1)
           for k in range(j_max + 1)
           for m in range(l_max + 1)
           if (l + m) % 2 == 0]
    N = len(idx)

    # monomial matrices from geovac (accurate entries)
    mono = [UnsymP(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
    n_mom = 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 30
    momA = pg.Moments(2.0 * alpha, n_mom)
    S, H1 = pg.one_body(mono, R, 1.0, momA)
    V = ngm.vee_matrix(mono, R, l_neumann=l_neumann, verbose=verbose)
    H = H1 + V + (1.0 / R) * S
    t_build = time.time() - t0

    # re-basing transforms
    Tx = np.array([laguerre_row(n, alpha, j_max + 1) for n in range(j_max + 1)])
    Te = np.array([legendre_row(n, l_max + 1) for n in range(l_max + 1)])

    # C[I, u] = d(muI==muu) Tx[A,j] Te[B,l] Tx[Cc,k] Te[D,m]
    pos = {t: i for i, t in enumerate(idx)}
    C = np.zeros((N, N))
    for I, (A, B, Cc, D, MU) in enumerate(idx):
        for u, (j, l, k, m, mu) in enumerate(idx):
            if mu != MU:
                continue
            C[I, u] = Tx[A, j] * Te[B, l] * Tx[Cc, k] * Te[D, m]

    S_o = C @ S @ C.T
    H_o = C @ H @ C.T

    # robust solve: symmetric orthogonalization on the (well-conditioned) S_o
    w, Uv = np.linalg.eigh(0.5 * (S_o + S_o.T))
    keep = w > 1e-12 * w[-1]
    X = Uv[:, keep] / np.sqrt(w[keep])
    Hp = X.T @ (0.5 * (H_o + H_o.T)) @ X
    E = float(np.linalg.eigvalsh(Hp)[0])

    condS_mono = np.linalg.cond(S)
    condS_orth = np.linalg.cond(S_o)
    de = 100.0 * (-1.0 - E) / DE_EXACT
    var = "" if E > E_EXACT - 5e-6 else "  <-NON-VARIATIONAL"
    print(f"  ({j_max},{l_max}) mu<={mu_max}  N={N:5d}  kept={int(keep.sum()):5d}  "
          f"E={E:12.7f}  D_e%={de:7.3f}  "
          f"cond(S):mono={condS_mono:.1e} orth={condS_orth:.1e}  "
          f"[{t_build:.0f}s]{var}")
    return de, E, N


if __name__ == "__main__":
    import sys
    print("field reference: TMR grid l_max=6 -> 99.97% ; Wolniewicz exact.")
    print("Paper 12 monomial ceiling: (3,3) mu<=1 = 99.09% (conditioning-capped)")
    print()
    cases = {
        "val": [(2, 2, 0), (3, 3, 0), (2, 2, 1), (3, 3, 1)],
        "push1": [(4, 4, 1), (5, 5, 1)],
        "delta": [(3, 3, 2), (4, 4, 2)],
        "push2": [(6, 6, 1)],
    }
    if len(sys.argv) == 4:            # single case: j l mu
        j, l, mu = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
        run(j, l, mu, verbose=True)
    else:
        which = sys.argv[1] if len(sys.argv) > 1 else "val"
        for (j, l, mu) in cases[which]:
            run(j, l, mu, verbose=True)
