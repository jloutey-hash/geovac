"""QR step 2: partial-wave (position-space) evaluator for Coulomb rings of n densities.
Each 1/|x - y| = sum_l g_l(r, r') P_l(cos), g_l = r<^l / r>^(l+1); the angular integrals around
a ring give (4 pi)^n / (2l+1)^(n-1), so
    Ring_n = (4 pi)^n sum_l (2l+1)^(1-n) R_l,   R_l = int prod_v r_v^2 rho_v(r_v) dr_v prod_edges g_l.
R_l is a sum over the n! orderings of nested one-dimensional integrals (the kernels are
semi-separable).  In x = log r each nesting level is an exponential smoothing
G_k(x) = int_{-inf}^x e^{-p_k (x-y)} w_k(y) G_{k-1}(y) dy, done panel-wise with exact
exponential weights on Gauss-Legendre panels (float64).  Open chains are the l = 0 term alone
with the ring edge removed (no angular coupling survives).
Densities: rho_v(r) = a_v^3 K0(a_v r)/(2 pi^2) (fingerprint, unit charge); vertex weights can be
replaced by theta_t-derivatives for the ODE checks."""
from __future__ import annotations
import itertools
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.signal import lfilter
from scipy.special import k0, k1

NG = 12
XG, WG = leggauss(NG)
S01 = (XG + 1) / 2          # nodes in [0, 1]
W01 = WG / 2


class Grid:
    def __init__(self, xmin=-36.0, xmax=np.log(90.0), H=0.04):
        self.H = H
        self.np_ = int(np.ceil((xmax - xmin) / H))
        self.x0 = xmin + H * np.arange(self.np_)
        self.x = (self.x0[:, None] + H * S01[None, :])     # (panels, NG)
        self.r = np.exp(self.x)
        self.cache = {}

    def mats(self, p):
        """M[i, j] = int_0^{s_i} e^{-pH(s_i - s)} L_j(s) ds (Lagrange basis on S01), times H;
        m[j] = same to s = 1; d = e^{-pH s_i}; D = e^{-pH}."""
        if p in self.cache:
            return self.cache[p]
        H = self.H
        xq, wq = leggauss(60)
        # Lagrange basis values at arbitrary points
        def lag(s):
            V = np.ones((len(s), NG))
            for j in range(NG):
                for m in range(NG):
                    if m != j:
                        V[:, j] *= (s - S01[m]) / (S01[j] - S01[m])
            return V
        M = np.zeros((NG, NG))
        for i in range(NG):
            s = S01[i] * (xq + 1) / 2
            w = S01[i] * wq / 2 * np.exp(-p * H * (S01[i] - s))
            M[i] = H * (w @ lag(s))
        s = (xq + 1) / 2
        m = H * ((wq / 2 * np.exp(-p * H * (1 - s))) @ lag(s))
        d = np.exp(-p * H * S01)
        out = (M, m, d, np.exp(-p * H))
        self.cache[p] = out
        return out

    def smooth(self, F, p):
        """G(x) = int_{-inf}^x e^{-p(x-y)} F(y) dy on the panel nodes (F zero below xmin)."""
        M, m, d, D = self.mats(p)
        c = F @ m                               # contribution of each panel to its end
        ends = lfilter([1.0], [1.0, -D], c)      # G at panel ends
        starts = np.concatenate([[0.0], ends[:-1]])
        return starts[:, None] * d[None, :] + F @ M.T

    def integrate(self, F):
        return self.H * np.sum(F @ W01)


def rho(a, r):
    return a**3 * k0(a * r) / (2 * np.pi**2)


def ring_l(grid, w, l, edges, n):
    """sum over orderings of nested integrals of prod_v r_v^{e_v} w_v; w[v] arrays r^2 rho on the grid
    (in x-variable the measure dr = r dx is included here)."""
    nbr = {v: [] for v in range(n)}
    for (u, v) in edges:
        nbr[u].append(v); nbr[v].append(u)
    tot = 0.0
    r = grid.r
    for perm in itertools.permutations(range(n)):
        rank = {v: i for i, v in enumerate(perm)}
        e = {v: sum((l if rank[v] < rank[u] else -l - 1) for u in nbr[v]) for v in range(n)}
        G = np.ones_like(r); P = 0
        for k, v in enumerate(perm):
            pk = P + e[v] + 3
            # H_k = r^{pk} G_k ; integrand s^{e} w H_{k-1} ds = s^{e + P + 1} w G ds/s... in x:
            F = w[v] * r**(e[v] + P + 1 - pk) * G      # = w r^{-2} G : powers absorbed
            if k < n - 1:
                G = grid.smooth(F, pk)
                P = pk
            else:
                tot += grid.integrate(F * r**pk)
    return tot


def ring(grid, avals, n, lmax, closed=True, weights=None):
    edges = [(i, (i + 1) % n) for i in range(n)] if closed else [(i, i + 1) for i in range(n - 1)]
    w = [grid.r**2 * rho(a, grid.r) if weights is None or weights[v] is None else weights[v]
         for v, a in enumerate(avals)]
    if not closed:
        return (4 * np.pi)**n * ring_l(grid, w, 0, edges, n), None
    terms = []
    for l in range(lmax + 1):
        terms.append((2 * l + 1)**(1 - n) * ring_l(grid, w, l, edges, n))
    return (4 * np.pi)**n * np.array(terms), terms


def tail_sum(terms, n, kfit=6, lfit=40):
    """Sum of the l-series with the tail fitted by sum_k c_k (2l+1)^{-k}, k >= 2n-2."""
    T = np.array(terms); L = len(T) - 1
    ls = np.arange(L - lfit + 1, L + 1)
    ks = np.arange(2 * n - 2, 2 * n - 2 + kfit)
    A = np.array([[(2 * l + 1.0)**(-k) for k in ks] for l in ls])
    c = np.linalg.lstsq(A, T[ls], rcond=None)[0]
    tail = 0.0
    for k, ck in zip(ks, c):
        # sum_{l > L} (2l+1)^{-k} via Hurwitz: 2^{-k} zeta(k, L + 3/2)
        from mpmath import zeta
        tail += ck * float(2.0**(-k) * zeta(k, L + 1.5))
    return T.sum() + tail, T.sum()


# ---------------------------------------------------------------- batched version (theta-derivatives)
def theta_bessel(J):
    """theta^j [y^3 K0(y)] = sum c[(m, nu)] y^m K_nu(y), theta = y d/dy."""
    out = [{(3, 0): 1.0}]
    for _ in range(J):
        new = {}
        for (m, nu), c in out[-1].items():
            if nu == 0:     # theta(y^m K0) = m y^m K0 - y^{m+1} K1
                new[(m, 0)] = new.get((m, 0), 0) + m * c
                new[(m + 1, 1)] = new.get((m + 1, 1), 0) - c
            else:           # theta(y^m K1) = (m-1) y^m K1 - y^{m+1} K0
                new[(m, 1)] = new.get((m, 1), 0) + (m - 1) * c
                new[(m + 1, 0)] = new.get((m + 1, 0), 0) - c
        out.append(new)
    return out


def d_weights(grid, ts, J):
    """array (len(ts)*(J+1), panels, NG): r^2 theta_t^j [t^3 K0(t r)] / (2 pi^2)."""
    tb = theta_bessel(J)
    r = grid.r
    out = []
    for t in ts:
        y = t * r
        K = {0: k0(y), 1: k1(y)}
        for j in range(J + 1):
            val = sum(c * y**m * K[nu] for (m, nu), c in tb[j].items())
            out.append(r**2 * val / r**3 / (2 * np.pi**2))
    return np.array(out)


def smooth_b(grid, F, p):
    M, m, d, D = grid.mats(p)
    c = F @ m
    ends = lfilter([1.0], [1.0, -D], c, axis=-1)
    starts = np.concatenate([np.zeros(ends.shape[:-1] + (1,)), ends[..., :-1]], axis=-1)
    return starts[..., None] * d + F @ M.T


def ring_l_b(grid, w, l, edges, n):
    nbr = {v: [] for v in range(n)}
    for (u, v) in edges:
        nbr[u].append(v); nbr[v].append(u)
    tot = 0.0
    r = grid.r
    for perm in itertools.permutations(range(n)):
        rank = {v: i for i, v in enumerate(perm)}
        e = {v: sum((l if rank[v] < rank[u] else -l - 1) for u in nbr[v]) for v in range(n)}
        G = 1.0; P = 0
        for k, v in enumerate(perm):
            pk = P + e[v] + 3
            F = w[v] * r**(e[v] + P + 1 - pk) * G
            if k < n - 1:
                G = smooth_b(grid, F, pk)
                P = pk
            else:
                tot = tot + grid.H * np.sum((F * r**pk) @ W01, axis=-1)
    return tot


def batch(grid, n, tpos, ts, J, lmax, closed=True):
    """values[(t, j)] of theta_t^j Ring (closed) or Chain (open), t-density at vertex tpos,
    other densities scale 1.  Returns array (len(ts), J+1) and the per-l terms."""
    wd = d_weights(grid, ts, J)
    w1 = grid.r**2 * rho(1.0, grid.r)
    w = [wd if v == tpos else w1 for v in range(n)]
    if closed:
        edges = [(i, (i + 1) % n) for i in range(n)]
        terms = [(2 * l + 1)**(1 - n) * ring_l_b(grid, w, l, edges, n) for l in range(lmax + 1)]
        terms = (4 * np.pi)**n * np.array(terms)            # (L+1, nt*(J+1))
        vals = np.array([tail_sum(terms[:, i], n)[0] for i in range(terms.shape[1])])
    else:
        edges = [(i, i + 1) for i in range(n - 1)]
        vals = (4 * np.pi)**n * ring_l_b(grid, w, 0, edges, n)
        terms = None
    return vals.reshape(len(ts), J + 1), terms
