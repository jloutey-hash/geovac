"""#3 follow-on: DIRECT float64 one-body build in the orthogonal basis, no monomials.

Builds the prolate two-electron one-body matrices S and H1 = T + V_ne DIRECTLY in
the Laguerre(xi) x Legendre(eta) basis via banded recurrence operators, and validates
against the mpf re-basing pipeline (geovac.prolate_recondition).

This pass: the SIGMA sector (mu = 0), which exercises every mechanism --
  * radial overlaps r0/r1/r2 via the Laguerre multiply-by-z operator Z + orthogonality,
  * radial kinetic  K_rad = alpha^2 (2D-I)^T (r2 - r0) (2D-I),  D = d/dz operator,
  * angular overlap a0 (Legendre norm), a2 (eta^2, tridiagonal Y), and the DIAGONAL
    angular kinetic K_ang = l(l+1) a0 (Legendre ODE eigenvalue).
mu > 0 adds the azimuthal mu^2 term (shifted weight) -- the next block.

Ground truth: geovac.prolate_recondition.{one_body_mp, _factored_cob} (mpf, exact).
"""
from __future__ import annotations

import time
import numpy as np
import mpmath as mp

from geovac import prolate_recondition as pr

R = pr.R_DEFAULT


# ==================== radial operators (Laguerre L^{(0)}) ====================
def _Z(n: int) -> np.ndarray:
    """Symmetric 'multiply by z' on L^{(0)}: <L_i|z|L_j>_{e^{-z}} (h=1)."""
    Z = np.zeros((n, n))
    for j in range(n):
        Z[j, j] = 2 * j + 1
        if j + 1 < n:
            Z[j, j + 1] = Z[j + 1, j] = -(j + 1)
    return Z


def _D(n: int) -> np.ndarray:
    """d/dz on L^{(0)}: L_n' = -sum_{k<n} L_k  ->  D[k,n] = -1 for k<n."""
    D = np.zeros((n, n))
    for col in range(n):
        for row in range(col):
            D[row, col] = -1.0
    return D


def radial_blocks_mu0(n_r: int, alpha: float, pad: int = 4):
    """r0, r1, r2 (xi^0,1,2) and K_rad for mu=0, float64, banded.

    Operators are built PADDED (size n+pad) and the polynomial products (xi^2 = Z^2,
    the kinetic) formed there, then truncated to the basis size -- otherwise the
    z-coupling to index n_r+1 is dropped and r2[n_r,n_r] is undercounted.
    """
    n = n_r + 1
    m = n + pad
    C = float(mp.e ** (-2 * alpha) / (2 * alpha))     # int L_iL_j e^{-2a xi}dxi = C delta
    Z = _Z(m)
    I = np.eye(m)
    xi = I + Z / (2.0 * alpha)                          # xi = 1 + z/2a  (operator)
    r0 = C * I
    r1 = C * xi
    r2 = C * (xi @ xi)
    V = 2.0 * _D(m) - I                                 # (2D - I)
    K_rad = alpha * alpha * (V.T @ (r2 - r0) @ V)
    r0, r1, r2 = r0[:n, :n], r1[:n, :n], r2[:n, :n]
    K_rad = K_rad[:n, :n]
    return r0, r1, r2, 0.5 * (K_rad + K_rad.T)


# ==================== angular operators (Legendre P_l) =======================
def angular_blocks_mu0(l_max: int, pad: int = 4):
    """a0 (norm), a2 (eta^2), K_ang (l(l+1) a0) for mu=0, float64.

    Y (multiply-by-eta) built PADDED so Y^2 keeps the coupling to index l_max+1,
    then truncated -- else a2[l_max,l_max] is undercounted.
    """
    n = l_max + 1
    m = n + pad
    h = np.array([2.0 / (2 * l + 1) for l in range(m)])
    Y = np.zeros((m, m))                                # eta P_l = [(l+1)P_{l+1}+l P_{l-1}]/(2l+1)
    for l in range(m):
        if l + 1 < m:
            Y[l + 1, l] = (l + 1) / (2 * l + 1)
        if l - 1 >= 0:
            Y[l - 1, l] = l / (2 * l + 1)
    Y2 = Y @ Y
    a2 = (h[:, None] * Y2)                              # <P_i|eta^2 P_j> = h_i (Y^2)[i,j]
    a2 = 0.5 * (a2 + a2.T)[:n, :n]
    a0 = np.diag(h[:n])
    K_ang = np.diag([l * (l + 1) for l in range(n)]) @ a0
    return a0, a2, K_ang


# ==================== assemble two-electron S, H1 (mu=0) =====================
def build_direct_mu0(j_max: int, l_max: int, alpha: float):
    r0, r1, r2, K_rad = radial_blocks_mu0(j_max, alpha)
    a0, a2, K_ang = angular_blocks_mu0(l_max)

    idx = pr._product_index(j_max, l_max, 0)           # (j,l,k,m,mu=0)
    N = len(idx)
    h6 = float((mp.mpf(R) / 2) ** 6)
    cc = 4.0 * np.pi ** 2                                # phi_cc(0)
    pref_T = 0.5 * (4.0 / R ** 2) * h6
    pref_V = -(4.0 / R) * h6

    # one-electron blocks over (radial a, angular b):
    #   ov = r2[a,c] a0[b,d] - r0[a,c] a2[b,d]
    #   vne = r1[a,c] a0[b,d]
    #   kin(grad) = K_rad[a,c] a0[b,d] + r0[a,c] K_ang[b,d]
    def ov(ra, la, rc, lc):
        return r2[ra, rc] * a0[la, lc] - r0[ra, rc] * a2[la, lc]

    def vne(ra, la, rc, lc):
        return r1[ra, rc] * a0[la, lc]

    def kin(ra, la, rc, lc):
        return K_rad[ra, rc] * a0[la, lc] + r0[ra, rc] * K_ang[la, lc]

    S = np.zeros((N, N))
    H = np.zeros((N, N))
    for i, (ji, li, ki, mi, _) in enumerate(idx):
        for j, (jj, lj, kj, mj, _) in enumerate(idx):
            o1 = ov(ji, li, jj, lj)                     # electron 1: radial j, angular l
            o2 = ov(ki, mi, kj, mj)                     # electron 2: radial k, angular m
            S[i, j] = h6 * cc * o1 * o2
            k1 = kin(ji, li, jj, lj)
            k2 = kin(ki, mi, kj, mj)
            v1 = vne(ji, li, jj, lj)
            v2 = vne(ki, mi, kj, mj)
            H[i, j] = (pref_T * cc * (k1 * o2 + o1 * k2)
                       + pref_V * cc * (v1 * o2 + o1 * v2))
    return S, H


# ============================================================================
# GENERAL-mu one-body engine: exact mpf 1D blocks (tiny -> fast) + float64 assembly.
# Each one-electron block factors radial x angular; blocks reuse the validated
# _ov/_vne/_kin polynomial logic (prolate_recondition mpf helpers), summed over
# the orthogonal-basis monomial expansions.  The O(N^2) two-electron assembly is
# float64 -- so the engine is fast and free of the monomial dynamic-range blowup.
# ============================================================================
def _mx_poly(jx, alpha):
    """mu=0 radial derivative poly: d/dxi[xi^jx e^{-a xi}] = mx * e^{-a xi}."""
    a = mp.mpf(alpha)
    t = pr._shift([-a], jx)
    if jx > 0:
        t = pr._pa(t, pr._shift([mp.mpf(jx)], jx - 1))
    return t


def _nx_poly(jx, mu, alpha):
    """mu>0 radial derivative poly (mirrors prolate_general_m._kin nx)."""
    a, m = mp.mpf(alpha), mp.mpf(mu)
    t = pr._shift([m], jx + 1)
    t = pr._ps(t, pr._pm(pr._shift([a], jx), pr._xi2m1(1)))
    if jx > 0:
        t = pr._pa(t, pr._pm(pr._shift([mp.mpf(jx)], jx - 1), pr._xi2m1(1)))
    return t


def _my_poly(lx):
    """mu=0 angular derivative poly: d/deta[eta^lx]."""
    return pr._shift([mp.mpf(lx)], lx - 1) if lx > 0 else [mp.mpf(0)]


def _ny_poly(lx, mu):
    """mu>0 angular derivative poly (mirrors prolate_general_m._kin ny)."""
    m = mp.mpf(mu)
    t = pr._shift([-m], lx + 1)
    if lx > 0:
        t = pr._pa(t, pr._pm(pr._shift([mp.mpf(lx)], lx - 1), pr._meta2(1)))
    return t


def _combine(coeffs, poly_fn):
    """sum_p coeffs[p] * poly_fn(p)  (linearity of the derivative over the basis)."""
    out = [mp.mpf(0)]
    for p, cp in enumerate(coeffs):
        if cp == 0:
            continue
        out = pr._pa(out, [cp * t for t in poly_fn(p)])
    return out


def radial_blocks(n_r, mu, alpha, A):
    """r0,r1,r2 (weight (xi^2-1)^mu), K_rad, and r0',r2' (weight (xi^2-1)^{mu-1})."""
    n = n_r + 1
    Lc = [pr.laguerre_coeffs(a, alpha, n) if mu == 0
          else pr.assoc_laguerre_coeffs(a, mu, alpha, n) for a in range(n)]
    xw = pr._xi2m1(mu)
    xw1 = pr._xi2m1(mu - 1) if mu >= 1 else pr._xi2m1(1)   # shifted weight for K/azi
    wk = 1 if mu == 0 else mu - 1                          # (xi^2-1) power in K_rad

    def blk(xi_pow, weight_poly):
        M = np.empty((n, n), object)
        for a in range(n):
            for c in range(a, n):
                prod = pr._pm(pr._pm(Lc[a], Lc[c]), weight_poly)
                M[a, c] = M[c, a] = pr._mom_xi(pr._shift(prod, xi_pow), A)
        return M

    r0, r1, r2 = blk(0, xw), blk(1, xw), blk(2, xw)
    NX = [_combine(Lc[a], (lambda p, mu=mu, alpha=alpha:
                           _mx_poly(p, alpha) if mu == 0 else _nx_poly(p, mu, alpha)))
          for a in range(n)]
    K_rad = np.empty((n, n), object)
    xwk = pr._xi2m1(wk)
    for a in range(n):
        for c in range(a, n):
            prod = pr._pm(pr._pm(NX[a], NX[c]), xwk)
            K_rad[a, c] = K_rad[c, a] = pr._mom_xi(prod, A)
    # azi radial (mu>0): weight (xi^2-1)^{mu-1}
    if mu >= 1:
        r0p = blk(0, xw1)
        r2p = blk(2, xw1)
    else:
        r0p = r2p = None
    return r0, r1, r2, K_rad, r0p, r2p


def angular_blocks(l_max, mu):
    """a0,a2 (weight (1-eta^2)^mu), K_ang, and a0',a2' (weight (1-eta^2)^{mu-1})."""
    n = l_max + 1
    Gc = [pr.legendre_coeffs(b, n) if mu == 0
          else pr.gegenbauer_coeffs(b, mp.mpf(mu) + mp.mpf('0.5'), n) for b in range(n)]
    yw = pr._meta2(mu)
    yw1 = pr._meta2(mu - 1) if mu >= 1 else pr._meta2(1)
    wk = 1 if mu == 0 else mu - 1

    def blk(eta_pow, weight_poly):
        M = np.empty((n, n), object)
        for b in range(n):
            for d in range(b, n):
                prod = pr._pm(pr._pm(Gc[b], Gc[d]), weight_poly)
                M[b, d] = M[d, b] = pr._mom_eta(pr._shift(prod, eta_pow))
        return M

    a0, a2 = blk(0, yw), blk(2, yw)
    NY = [_combine(Gc[b], (lambda r, mu=mu: _my_poly(r) if mu == 0 else _ny_poly(r, mu)))
          for b in range(n)]
    K_ang = np.empty((n, n), object)
    ywk = pr._meta2(wk)
    for b in range(n):
        for d in range(b, n):
            prod = pr._pm(pr._pm(NY[b], NY[d]), ywk)
            K_ang[b, d] = K_ang[d, b] = pr._mom_eta(prod)
    if mu >= 1:
        a0p = blk(0, yw1)
        a2p = blk(2, yw1)
    else:
        a0p = a2p = None
    return a0, a2, K_ang, a0p, a2p


def build_direct_full(j_max, l_max, mu_max, alpha):
    """Complete one-body S, H1 = T + V_ne for all mu, float64, via exact mpf blocks."""
    with mp.workdps(pr.DEFAULT_DPS):
        A = pr.ngm._mono_moments(2.0 * alpha, 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 20)
        h6 = float((mp.mpf(R) / 2) ** 6)
        pref_T = 0.5 * (4.0 / R ** 2) * h6
        pref_V = -(4.0 / R) * h6
        # per-mu float64 block sets
        RB, AB, ccv, ssv = {}, {}, {}, {}
        for mu in range(mu_max + 1):
            r0, r1, r2, Kr, r0p, r2p = radial_blocks(j_max, mu, alpha, A)
            a0, a2, Ka, a0p, a2p = angular_blocks(l_max, mu)
            f = lambda M: None if M is None else np.array([[float(M[i, j])
                        for j in range(M.shape[1])] for i in range(M.shape[0])])
            RB[mu] = tuple(f(M) for M in (r0, r1, r2, Kr, r0p, r2p))
            AB[mu] = tuple(f(M) for M in (a0, a2, Ka, a0p, a2p))
            ccv[mu] = float(pr._phi_cc(mu))
            ssv[mu] = float(pr._phi_ss(mu))

    idx = pr._product_index(j_max, l_max, mu_max)
    N = len(idx)
    S = np.zeros((N, N))
    H = np.zeros((N, N))
    L1 = l_max + 1
    idx_arr = np.array(idx)
    for mu in range(mu_max + 1):
        r0, r1, r2, Kr, r0p, r2p = RB[mu]
        a0, a2, Ka, a0p, a2p = AB[mu]
        cc, ss = ccv[mu], ssv[mu]
        # single-electron blocks over se = radial*(l_max+1) + angular  (Kronecker)
        OV = np.kron(r2, a0) - np.kron(r0, a2)
        GR = np.kron(Kr, a0) + np.kron(r0, Ka)
        VN = np.kron(r1, a0)
        AZ = (mu * mu) * (np.kron(r2p, a0p) - np.kron(r0p, a2p)) if r0p is not None \
            else np.zeros_like(OV)
        # this-mu basis functions and their two single-electron indices
        rows = np.where(idx_arr[:, 4] == mu)[0]
        se1 = idx_arr[rows, 0] * L1 + idx_arr[rows, 1]      # (j, l)
        se2 = idx_arr[rows, 2] * L1 + idx_arr[rows, 3]      # (k, m)
        OV1, OV2 = OV[np.ix_(se1, se1)], OV[np.ix_(se2, se2)]
        GR1, GR2 = GR[np.ix_(se1, se1)], GR[np.ix_(se2, se2)]
        VN1, VN2 = VN[np.ix_(se1, se1)], VN[np.ix_(se2, se2)]
        AZ1, AZ2 = AZ[np.ix_(se1, se1)], AZ[np.ix_(se2, se2)]
        Sb = h6 * cc * OV1 * OV2
        Hb = (pref_T * (cc * (GR1 * OV2 + OV1 * GR2) + ss * (AZ1 * OV2 + OV1 * AZ2))
              + pref_V * cc * (VN1 * OV2 + OV1 * VN2))
        S[np.ix_(rows, rows)] = Sb
        H[np.ix_(rows, rows)] = Hb
    return S, H


def ground_truth_full(j_max: int, l_max: int, mu_max: int, alpha: float):
    """mpf S_o, H1_o (one-body, all mu) from prolate_recondition, downcast."""
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        n_mom = 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 20
        A = pr.ngm._mono_moments(2.0 * alpha, n_mom)
        S, H1 = pr.one_body_mp(fns, alpha, R, A)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1) if (l + m) % 2 == 0])
        # same (gegenbauer / mu-adapted) basis as build_direct_full
        Tr, Ta = pr._transforms_per_mu("gegenbauer", j_max, l_max, mu_max, alpha)
        S_o = pr._factored_cob(S, mu_max + 1, Nr, Na, Tr, Ta)
        H_o = pr._factored_cob(H1, mu_max + 1, Nr, Na, Tr, Ta)
        N = len(fns)
        Sf = np.array([[float(S_o[i, j]) for j in range(N)] for i in range(N)])
        Hf = np.array([[float(H_o[i, j]) for j in range(N)] for i in range(N)])
    return Sf, Hf


def ground_truth_mu0(j_max: int, l_max: int, alpha: float):
    """mpf S_o, H1_o (one-body only) from prolate_recondition, downcast to float64."""
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, 0)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        n_mom = 4 * j_max + 4 + 4 * l_max + 16
        A = pr.ngm._mono_moments(2.0 * alpha, n_mom)
        S, H1 = pr.one_body_mp(fns, alpha, R, A)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1) if (l + m) % 2 == 0])
        Tr, Ta = pr._transforms_per_mu("laguerre_legendre", j_max, l_max, 0, alpha)
        S_o = pr._factored_cob(S, 1, Nr, Na, Tr, Ta)
        H_o = pr._factored_cob(H1, 1, Nr, Na, Tr, Ta)
        N = len(fns)
        Sf = np.array([[float(S_o[i, j]) for j in range(N)] for i in range(N)])
        Hf = np.array([[float(H_o[i, j]) for j in range(N)] for i in range(N)])
    return Sf, Hf


def maxrel(A, B):
    d = 0.0
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            b = B[i, j]
            if abs(b) > 1e-12:
                d = max(d, abs(A[i, j] - b) / abs(b))
    return d


if __name__ == "__main__":
    alpha = 1.0
    print("A. CORRECTNESS: complete direct float64 one-body (all mu) vs mpf\n", flush=True)
    for (j_max, l_max, mu_max) in [(2, 2, 1), (2, 2, 2), (3, 3, 2)]:
        t0 = time.time()
        Sd, Hd = build_direct_full(j_max, l_max, mu_max, alpha)
        t_direct = time.time() - t0
        t1 = time.time()
        Sg, Hg = ground_truth_full(j_max, l_max, mu_max, alpha)
        t_mpf = time.time() - t1
        print(f"  ({j_max},{l_max}) mu<={mu_max} N={Sd.shape[0]:4d}  "
              f"relS={maxrel(Sd, Sg):.2e} relH={maxrel(Hd, Hg):.2e}  "
              f"[direct {t_direct:.3f}s vs mpf {t_mpf:.1f}s]", flush=True)

    print("\nB. SPEED (vectorized build) at production sizes\n", flush=True)
    for (j_max, l_max, mu_max) in [(5, 5, 2), (7, 7, 2)]:
        t0 = time.time()
        Sd, Hd = build_direct_full(j_max, l_max, mu_max, alpha)
        t_direct = time.time() - t0
        w = np.linalg.eigvalsh(0.5 * (Sd + Sd.T))
        condS = w[-1] / w[0] if w[0] > 0 else float('inf')
        print(f"  ({j_max},{l_max}) mu<={mu_max} N={Sd.shape[0]:5d}  "
              f"raw cond(S)={condS:.2e} (norm-spread; ~9e4 normalized)  "
              f"[build {t_direct:.2f}s]", flush=True)
