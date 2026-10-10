"""m-seam rung 4 (tests + engine): the converged referee for rungs 2/2S/3.

Same model (one electron, central Z_c, two unit charges at distance R in
directions (s, 0, +-c)), solved with a complete radial representation
(finite-element DVR, element boundary at r = R where the multipole kernel
kinks) and a one-centre partial-wave expansion to high l.  At linear geometry
m is exact, so each |m| block is small; the s^2 coefficient a2 of
E_in - E_out is exact second-order perturbation theory:

    a2(x) = <pi|V2|pi> + <pi|V1|psi1>,   (H0 - E_pi) psi1 = -V1 pi  (Sternheimer)

with V1 coupling the in-plane pi (m = 1c) to m = 0 and m = 2c, and the
out-of-plane pi (1s) to 2s.  No sum over states, so nothing is truncated but
l_max and the radial grid.  Expansion of the two-proton multipole sum in s
(theta_n = arcsin s, so d/ds = d/dtheta_n through second order):

  P_k(n+.r) + P_k(n-.r) = 2 P_k(u)                       (k even)
                        + 2 s P_k'(u) sin(t) cos(phi)     (k odd)
                        + s^2 [P_k''(u) sin^2 t cos^2 phi - u P_k'(u)]  (k even)

Driver: debug/renner_teller_rung4_grid.py (sweep data in debug/data/).
"""
import numpy as np
import pytest
from numpy.polynomial import legendre as npleg
from scipy.special import lpmv, gammaln

R = 2.0
ZC = 1.0


# ---------------------------------------------------------------- radial FE-DVR
def lobatto(n):
    """n-point Gauss-Lobatto nodes/weights on [-1, 1]."""
    c = np.zeros(n); c[-1] = 1
    inner = npleg.legroots(npleg.legder(c))
    x = np.concatenate(([-1.0], np.sort(inner), [1.0]))
    Pn1 = npleg.legval(x, c)
    w = 2.0 / (n * (n - 1) * Pn1 ** 2)
    return x, w


def lagrange_deriv(x):
    n = len(x)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                p = 1.0 / (x[j] - x[i])
                for k in range(n):
                    if k != i and k != j:
                        p *= (x[i] - x[k]) / (x[j] - x[k])
                D[i, j] = p
    # D[i, j] = L_j'(x_i)
    for i in range(n):
        D[i, i] = sum(1.0 / (x[i] - x[k]) for k in range(n) if k != i)
    return D


def fedvr(bounds, order):
    x, w = lobatto(order)
    D = lagrange_deriv(x)
    nodes, weights, blocks = [], [], []
    for e in range(len(bounds) - 1):
        a, b = bounds[e], bounds[e + 1]
        h = (b - a) / 2
        r = a + h * (x + 1)
        we = h * w
        Te = 0.5 * (D.T * we) @ D / h ** 2      # (1/2) int L_i' L_j'
        start = len(nodes) - 1 if e > 0 else 0
        if e == 0:
            nodes.extend(r); weights.extend(we)
        else:
            weights[-1] += we[0]
            nodes.extend(r[1:]); weights.extend(we[1:])
        blocks.append((start, Te))
    N = len(nodes)
    T = np.zeros((N, N))
    for start, Te in blocks:
        T[start:start + order, start:start + order] += Te
    W = np.array(weights)
    T = T / np.sqrt(np.outer(W, W))
    keep = slice(1, N - 1)                    # u(0) = u(rmax) = 0
    return np.array(nodes)[keep], T[keep, keep]


BOUNDS = [0, 0.15, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.9, 3.5, 4.3, 5.5, 7, 9, 12, 16, 21, 28, 37, 48, 62]


# ---------------------------------------------------------------- angular
NU = 200          # Gauss-Legendre points in u (exact to degree 399)
NPHI = 24         # trapezoid in phi (exact to trig degree 23)
U, WU = npleg.leggauss(NU)
PHI = 2 * np.pi * np.arange(NPHI) / NPHI
WPHI = np.full(NPHI, 2 * np.pi / NPHI)


def real_Y(l, m, kind):
    """Real orthonormal harmonic on the (U, PHI) grid, shape (NU, NPHI)."""
    norm = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - m + 1) - gammaln(l + m + 1)))
    P = norm * lpmv(m, l, U) * (-1) ** m
    if m == 0:
        return np.outer(P, np.ones(NPHI))
    f = np.cos(m * PHI) if kind == 'c' else np.sin(m * PHI)
    return np.sqrt(2) * np.outer(P, f)


def harmonics(m, kind, lmax, zpar=None):
    """Labels in an |m| block, optionally restricted to z-parity (-1)^(l+m)."""
    out = []
    for l in range(m, lmax + 1):
        if zpar is None or (-1) ** (l + m) == zpar:
            out.append((l, m, kind))
    return out


def ang_matrix(bra, ket, func):
    """<Y_a|func(u, phi)|Y_b> by exact quadrature; func returns (NU, NPHI)."""
    Ya = np.array([real_Y(*a) for a in bra])
    Yb = np.array([real_Y(*b) for b in ket])
    Wgt = np.outer(WU, WPHI) * func
    return np.einsum('aij,ij,bij->ab', Ya, Wgt, Yb)


SIN = np.sqrt(1 - U ** 2)
COSPHI = np.cos(PHI)


def Pk(k, der=0):
    c = np.zeros(k + 1); c[k] = 1
    for _ in range(der):
        c = npleg.legder(c)
    return npleg.legval(U, c) if len(c) else np.zeros_like(U)


def fk(r, k):
    rl, rg = np.minimum(r, R), np.maximum(r, R)
    return rl ** k / rg ** (k + 1)


# ---------------------------------------------------------------- operators
def h0_block(r, T, labs, kmax):
    """Linear-geometry H0 on one m block: dense (len(labs)*Nr)^2."""
    Nr, L = len(r), len(labs)
    H = np.zeros((L, Nr, L, Nr))
    lvals = np.array([a[0] for a in labs])
    for a in range(L):
        H[a, :, a, :] += T
        H[a, :, a, :] += np.diag(lvals[a] * (lvals[a] + 1) / (2 * r ** 2) - ZC / r)
    for k in range(0, kmax + 1, 2):
        A = ang_matrix(labs, labs, np.outer(2 * Pk(k), np.ones(NPHI)))
        if np.max(np.abs(A)) < 1e-14:
            continue
        f = fk(r, k)
        for a in range(L):
            for b in range(L):
                if A[a, b] != 0.0:
                    H[a, :, b, :] -= np.diag(A[a, b] * f)
    return H.reshape(L * Nr, L * Nr)


def radial_op(r, bra, ket, kset, angfun):
    """-sum_k f_k(r) <bra|angfun(k)|ket>, diagonal in r."""
    Nr = len(r)
    M = np.zeros((len(bra), Nr, len(ket), Nr))
    for k in kset:
        A = ang_matrix(bra, ket, angfun(k))
        if np.max(np.abs(A)) < 1e-14:
            continue
        f = fk(r, k)
        for a in range(len(bra)):
            for b in range(len(ket)):
                if abs(A[a, b]) > 0:
                    M[a, :, b, :] -= np.diag(A[a, b] * f)
    return M.reshape(len(bra) * Nr, len(ket) * Nr)


def v1_fun(k):
    return np.outer(2 * Pk(k, 1) * SIN, COSPHI)


def v2_fun(k):
    return (np.outer(Pk(k, 2) * SIN ** 2, COSPHI ** 2)
            - np.outer(U * Pk(k, 1), np.ones(NPHI)))


def a2_linear(lmax, order=10, bounds=BOUNDS, nsigma=8):
    """Returns dict with E_pi, a2 (in - out), its pieces, and the sigma/delta
    neighbours of pi with their V1 couplings (the 2x2 data)."""
    r, T = fedvr(bounds, order)
    kmax = 2 * lmax
    # pi is z-even: (-1)^(l+m) = +1 -> l odd for m = 1.  V0 and V1, V2 preserve z-parity.
    pi_c = harmonics(1, 'c', lmax, +1)
    pi_s = harmonics(1, 's', lmax, +1)
    sig = harmonics(0, 'c', lmax, +1)
    del_c = harmonics(2, 'c', lmax, +1)
    del_s = harmonics(2, 's', lmax, +1)
    H1 = h0_block(r, T, pi_c, kmax)
    ev, V = np.linalg.eigh(H1)
    Epi, pi = ev[0], V[:, 0]
    odd_k = range(1, kmax + 1, 2)
    even_k = range(0, kmax + 1, 2)
    out = {'lmax': lmax, 'Nr': len(r), 'E_pi': Epi}
    pieces = {}
    for tag, pilab, targets in (('in', pi_c, [('sigma', sig), ('delta', del_c)]),
                                ('out', pi_s, [('delta', del_s)])):
        direct = pi @ radial_op(r, pilab, pilab, even_k, v2_fun) @ pi
        relay = 0.0
        for name, lab in targets:
            H0t = h0_block(r, T, lab, kmax)
            rhs = radial_op(r, lab, pilab, odd_k, v1_fun) @ pi
            psi1 = np.linalg.solve(H0t - Epi * np.eye(len(H0t)), -rhs)
            relay += rhs @ psi1
            pieces[tag + '_' + name] = float(rhs @ psi1)
            if tag == 'in' and name == 'sigma':
                es, Vs = np.linalg.eigh(H0t)
                cpl = Vs.T @ rhs
                near = np.argsort(np.abs(es - Epi))[:nsigma]
                out['sigma_neighbours'] = [{'E': float(es[t]), 'gap': float(es[t] - Epi),
                                            'coupling': float(cpl[t]),
                                            'contrib': float(-cpl[t] ** 2 / (es[t] - Epi))}
                                           for t in near]
                out['sigma_sum_check'] = float(np.sum(-cpl ** 2 / (es - Epi)))
        pieces[tag + '_direct'] = float(direct)
        pieces[tag] = float(direct + relay)
    out['pieces'] = pieces
    out['a2'] = pieces['in'] - pieces['out']
    return out


def full_bent(lmax, svals=(0.04, 0.02, 0.01), order=10, bounds=BOUNDS):
    """Referee for the perturbative formulas: all m, finite bend, Richardson."""
    r, T = fedvr(bounds, order)
    Nr = len(r)
    labs_c = [(l, m, 'c') for l in range(lmax + 1) for m in range(0, l + 1)]
    labs_s = [(l, m, 's') for l in range(lmax + 1) for m in range(1, l + 1)]

    def lowest_pi(labs, sv, Eref):
        c = np.sqrt(1 - sv ** 2)
        L = len(labs)
        H = np.zeros((L, Nr, L, Nr))
        for a, lab in enumerate(labs):
            H[a, :, a, :] += T + np.diag(lab[0] * (lab[0] + 1) / (2 * r ** 2) - ZC / r)
        for k in range(0, 2 * lmax + 1):
            cosg_p = c * U[:, None] + sv * SIN[:, None] * COSPHI[None, :]
            cosg_m = -c * U[:, None] + sv * SIN[:, None] * COSPHI[None, :]
            ck = np.zeros(k + 1); ck[k] = 1
            func = npleg.legval(cosg_p, ck) + npleg.legval(cosg_m, ck)
            A = ang_matrix(labs, labs, func)
            f = fk(r, k)
            for a in range(L):
                for b in range(L):
                    if abs(A[a, b]) > 1e-14:
                        H[a, :, b, :] -= np.diag(A[a, b] * f)
        ev = np.linalg.eigvalsh(H.reshape(L * Nr, L * Nr))
        return ev[np.argmin(np.abs(ev - Eref))]

    pert = a2_linear(lmax, order, bounds)
    Epi = pert['E_pi']
    q = []
    for sv in svals:
        d = lowest_pi(labs_c, sv, Epi) - lowest_pi(labs_s, sv, Epi)
        q.append(d / sv ** 2)
    r1, r2 = (4 * q[1] - q[0]) / 3, (4 * q[2] - q[1]) / 3
    return r2, abs(r1 - r2), pert['a2']




# ================================================================ tests
COARSE = [0, 0.5, 1.2, 2.0, 3.0, 4.5, 7, 11, 17, 26, 40]


def test_radial_grid_is_exact_for_hydrogen():
    r, T = fedvr(BOUNDS, 10)
    for l, exact in ((0, -0.5), (1, -0.125)):
        E = np.linalg.eigvalsh(T + np.diag(l * (l + 1) / (2 * r ** 2) - 1 / r))[0]
        assert abs(E - exact) < 1e-10


def test_perturbative_a2_matches_a_real_bend():
    """Exact second-order a2 at linear geometry == finite-bend Richardson with
    all m mixed.  Excludes a wrong V1/V2 expansion or a wrong relay channel."""
    rich, spread, pert = full_bent(5, order=7, bounds=COARSE)
    assert spread < 2e-4
    assert abs(rich - pert) < 5 * spread


def test_delta_relay_cancels_in_the_splitting():
    """The |m| = 2 relay shifts both pi components equally, so the in/out
    splitting at s^2 is direct + sigma relay only."""
    p = a2_linear(9, order=7, bounds=COARSE)['pieces']
    assert abs(p['in_delta']) > 0.05
    assert abs(p['in_delta'] - p['out_delta']) < 1e-12
    assert abs(p['in'] - p['out'] - (p['in_direct'] - p['out_direct'] + p['in_sigma'])) < 1e-12


def test_dominant_sigma_partner_sits_below_pi():
    """One sigma state (the first excited z-even sigma) carries most of the
    sigma relay; it lies ~0.09 Ha BELOW pi at R = 2, not near-degenerate."""
    o = a2_linear(15, order=8, bounds=BOUNDS)
    nb = max(o['sigma_neighbours'], key=lambda n: abs(n['contrib']))
    assert -0.12 < nb['gap'] < -0.07
    assert nb['contrib'] > 0.8 * o['pieces']['in_sigma'] > 0
    assert abs(o['sigma_sum_check'] - o['pieces']['in_sigma']) < 1e-8


def test_a2_converges_in_l():
    """a2 and E_pi converge monotonically from above as l_max grows;
    the l -> infinity values (fit, l <= 55) are 0.44941 and -0.7385753."""
    vals = [a2_linear(lm, order=8, bounds=BOUNDS) for lm in (11, 15, 19)]
    a = [v['a2'] for v in vals]
    e = [v['E_pi'] for v in vals]
    assert a[0] > a[1] > a[2] > 0.4494
    assert (a[1] - a[2]) < (a[0] - a[1])
    assert e[0] > e[1] > e[2] > -0.738576
    assert abs(a[2] - 0.44941) < 0.006
