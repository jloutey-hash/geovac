r"""Backing test for Paper 61 sec:modular: the second-cusp Stokes amplitude of the
scale-integrated two-scale Bessel twist is the closed form

    2|A| = (c1 + b^2)^{3/2} / (4 sqrt(pi) sqrt(c1) b).

Added v5.12.5 to close a coverage gap: this [MEASURED] claim previously had NO
permanent backing (only transient debug/ drivers _stokes_amp*.py). The claim was
validated there only to ~1% (fit-to-closed-form ratio 0.988-0.989), which this
test HARDENS to the convergence statement it actually is.

Object (Paper 59/61): the fibre J(c_s,c_t,b)=sum_n F_n c_t^{n+1} m_n(c_s,b), with
F_n the Taylor coeffs of e^{-sqrt(1+e)}[(1+e)^{-3/2}+3(1+e)^{-2}+3(1+e)^{-5/2}]
and m_n(c_s,b)=int_0^inf k^{2n} j0(kb) P(c_s,k) dk the one-mass moments. The
one-mass (c_t->0) expansion's coefficients a_n = F_n m_n are Gevrey-(2n)!; their
dominant Borel singularity is at z* = -(sqrt(c1) - i b)^2 with a (z*-w)^{3/2}
branch (exponent p=-5/2), so

    a_n |z*|^n / (2n)! * n^{5/2}  ->  P cos(n theta) + Q sin(n theta) + O(1/n),
    theta = arg(z*),   2|A|_fit = sqrt(P^2 + Q^2).

VERIFIED (this test, at all four prose points (c1,b) = (0.2,0.5), (0.1,0.5),
(0.3,0.5), (0.2,0.7); rewritten 2026-11-03 from one point at 3.5%): with the
leading + 1/n fit, 2|A|_fit rises MONOTONICALLY toward the closed form from below
across windows ending at n = 40..60; with the 1/n^2 deficit included in the fit,
the top window and its 1/n^2 Richardson limit are within 5e-4 of the closed form
(prose: ~4e-4); and no algebraic factor (c1+b^2, b^2/(c1+b^2), c1/(c1+b^2),
sqrt(2/pi)) is consistent with the measured ratio -- only the factor 1 is.

Self-contained (no debug/ import), per the permanent-record policy. Slow (~7 min per
point; the four points are independent and parallelize under -n).
"""
import mpmath as mp
import numpy as np
import pytest


def _fast_gl(N):
    """Gauss-Legendre nodes/weights on [-1,1] at mpmath precision (numpy seed)."""
    x0, _ = np.polynomial.legendre.leggauss(N)
    xs, ws = [], []
    for xd in x0:
        x = mp.mpf(float(xd))
        for _ in range(6):
            p0, p1 = mp.mpf(1), x
            for k in range(2, N + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dP = N * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dP
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-(mp.mp.dps + 6)):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, N + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dP = N * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dP * dP))
    return xs, ws


def _Pcs(cs, k):
    d = mp.sqrt(cs * k * k + 1)
    return cs * mp.e ** (-d) * (d ** -3 + 3 * d ** -4 + 3 * d ** -5)


def _F_taylor(N):
    f = lambda e: mp.e ** (-mp.sqrt(1 + e)) * (
        (1 + e) ** (mp.mpf(-3) / 2) + 3 * (1 + e) ** -2 + 3 * (1 + e) ** (mp.mpf(-5) / 2))
    return mp.taylor(f, 0, N)


def _dcoeffs(cs, A, M):
    """Laurent tail coeffs of P(cs,k) e^{Ak} k^3 (A = sqrt(cs) = tail decay rate)."""
    k0 = max(mp.mpf(60), 35 / mp.sqrt(cs))
    ks = [k0 * (i + 1) for i in range(M + 1)]
    gs = [_Pcs(cs, k) * mp.e ** (A * k) * k ** 3 for k in ks]
    V = mp.matrix(M + 1, M + 1)
    for i, k in enumerate(ks):
        for m in range(M + 1):
            V[i, m] = 1 / k ** m
    return list(mp.lu_solve(V, mp.matrix(gs)))


def _m_hi(cs, b, n, d, A, xs, ws, K, M):
    """m_n = int_0^inf k^{2n} j0(kb) P(cs,k) dk: bounded GL on [0,K] + gammainc tail."""
    bnd = mp.mpf(0)
    for x, w in zip(xs, ws):
        k = K * (x + 1) / 2
        wk = K * w / 2
        j0 = mp.sin(k * b) / (k * b) if k * b > 1e-40 else mp.mpf(1)
        bnd += wk * (k ** (2 * n)) * j0 * _Pcs(cs, k)
    z = A - 1j * b
    t = mp.mpc(0)
    for j in range(M + 1):
        o = 2 * n - 3 - j
        t += d[j] * z ** (-o) * mp.gammainc(o, z * K)
    return bnd + (t / b).imag




# ----------------------------------------------------------------------------------------------
# Four-point guard (2026-11-03, group3 guard pass).  The prose claims the closed form at FOUR
# points (c1, b) = (0.2,0.5), (0.1,0.5), (0.3,0.5), (0.2,0.7) with a Richardson limit of 1 to
# ~4e-4; the earlier test covered one point at 3.5%.  Reaching 4e-4 needs coefficients to n = 60
# and the 1/n^2 term in the fit, so the bounded quadrature is done ONCE per node for all n
# (identical object to _m_hi, cross-checked below) instead of once per n.
# ----------------------------------------------------------------------------------------------

POINTS = [('0.2', '0.5'), ('0.1', '0.5'), ('0.3', '0.5'), ('0.2', '0.7')]
N_MAX = 60


def _all_moments(cs, b, N, xs, ws, K, M, d, A):
    """m_n for n = 0..N in one pass over the GL nodes (same bounded + gammainc tail as _m_hi)."""
    bnd = [mp.mpf(0)] * (N + 1)
    for x, w in zip(xs, ws):
        k = K * (x + 1) / 2
        wk = K * w / 2
        j0 = mp.sin(k * b) / (k * b) if k * b > 1e-40 else mp.mpf(1)
        base = wk * j0 * _Pcs(cs, k)
        k2 = k * k
        p = mp.mpf(1)
        for n in range(N + 1):
            bnd[n] += base * p
            p *= k2
    z = A - 1j * b
    out = []
    for n in range(N + 1):
        t = mp.mpc(0)
        for j in range(M + 1):
            o = 2 * n - 3 - j
            t += d[j] * z ** (-o) * mp.gammainc(o, z * K)
        out.append(bnd[n] + (t / b).imag)
    return out


def _coeffs(cs, b, N=N_MAX, Nq=2400, K=mp.mpf(110), M=12):
    cs, b = mp.mpf(cs), mp.mpf(b)
    A = mp.sqrt(cs)
    d = _dcoeffs(cs, A, M)
    Fn = _F_taylor(N)
    xs, ws = _fast_gl(Nq)
    m = _all_moments(cs, b, N, xs, ws, K, M, d, A)
    return [Fn[n] * m[n] for n in range(N + 1)], (d, xs, ws, K, M, A)


def _fit(a, cs, b, nlo, nhi, orders):
    """2|A| from a_n |z*|^n n^{5/2}/(2n)! = sum_{j<orders} (P_j cos + Q_j sin)(n theta)/n^j."""
    cs, b = mp.mpf(cs), mp.mpf(b)
    zst = -(mp.sqrt(cs) - 1j * b) ** 2
    absz, th = abs(zst), mp.arg(zst)
    rows, rhs = [], []
    for n in range(nlo, nhi + 1):
        y = a[n] * absz ** n / mp.factorial(2 * n) * mp.mpf(n) ** (mp.mpf(5) / 2)
        c, s = mp.cos(n * th), mp.sin(n * th)
        r = []
        for j in range(orders):
            r += [c / mp.mpf(n) ** j, s / mp.mpf(n) ** j]
        rows.append(r)
        rhs.append(y)
    Mx = mp.matrix(rows)
    sol = mp.lu_solve(Mx.T * Mx, Mx.T * mp.matrix(rhs))
    return mp.sqrt(sol[0] ** 2 + sol[1] ** 2)


def _closed_form(cs, b):
    cs, b = mp.mpf(cs), mp.mpf(b)
    return (cs + b * b) ** (mp.mpf(3) / 2) / (4 * mp.sqrt(mp.pi) * mp.sqrt(cs) * b)


@pytest.mark.slow
@pytest.mark.parametrize("cs,b", POINTS)
def test_paper61_stokes_amplitude_four_points(cs, b):
    mp.mp.dps = 40
    a, (d, xs, ws, K, M, A) = _coeffs(cs, b)
    # the one-pass moments are the same object as the per-n _m_hi evaluator
    for n in (5, 40):
        ref = mp.mpf(_m_hi(mp.mpf(cs), mp.mpf(b), n, d, A, xs, ws, K, M))
        assert abs(a[n] / _F_taylor(n)[n] - ref) <= mp.mpf(10) ** -25 * abs(ref), n
    pred = _closed_form(cs, b)

    # (1) the prose's fit (leading + 1/n): monotone from below toward the closed form
    lo = [_fit(a, cs, b, nhi - 12, nhi, 2) / pred for nhi in (40, 44, 48, 52, 56, 60)]
    assert all(x < y for x, y in zip(lo, lo[1:])), [mp.nstr(x, 8) for x in lo]
    assert all(x < 1 for x in lo), [mp.nstr(x, 8) for x in lo]

    # (2) with the 1/n^2 finite-order deficit resolved, the amplitude is 1 to ~4e-4:
    #     the top-window fit and its 1/n^2 Richardson limit both sit within 5e-4 of 1
    w = [(nhi - 12, nhi) for nhi in (54, 60)]
    amps = [_fit(a, cs, b, lo_, hi_, 3) for lo_, hi_ in w]
    n1, n2 = [mp.mpf(lo_ + hi_) / 2 for lo_, hi_ in w]
    C = (amps[1] - amps[0]) / (1 / n1 ** 2 - 1 / n2 ** 2)
    rich = (amps[1] + C / n2 ** 2) / pred
    assert abs(amps[1] / pred - 1) < mp.mpf('5e-4'), mp.nstr(amps[1] / pred, 10)
    assert abs(rich - 1) < mp.mpf('5e-4'), mp.nstr(rich, 10)

    # (3) only the algebraic factor 1 is consistent with the measured ratio
    csm, bm = mp.mpf(cs), mp.mpf(b)
    fac = amps[1] / pred
    for name, val in (('c1+b2', csm + bm * bm), ('b2/(c1+b2)', bm * bm / (csm + bm * bm)),
                      ('c1/(c1+b2)', csm / (csm + bm * bm)), ('sqrt2/sqrtpi', mp.sqrt(2 / mp.pi))):
        assert abs(fac - val) > mp.mpf('0.05'), (name, mp.nstr(fac, 8), mp.nstr(val, 8))
