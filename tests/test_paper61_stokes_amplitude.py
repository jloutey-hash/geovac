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

VERIFIED (this test): 2|A|_fit rises MONOTONICALLY toward the closed form from
below across increasing-n windows, its 1/n^2 Richardson limit equals the closed
form, and no algebraic factor (c1+b^2, b^2/(c1+b^2), c1/(c1+b^2)) is consistent
with the measured ratio -- only the factor 1 is. A wrong closed form (any of the
algebraic factors, or a sqrt(2)/pi slip) fails the window ratio and the limit.

Self-contained (no debug/ import), per the permanent-record policy. Slow (~2 min).
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


def _amp_windows(cs, b, N, windows, Nq=2400, K=mp.mpf(110), M=12):
    cs, b = mp.mpf(cs), mp.mpf(b)
    A = mp.sqrt(cs)
    d = _dcoeffs(cs, A, M)
    Fn = _F_taylor(N)
    xs, ws = _fast_gl(Nq)
    a = [Fn[n] * _m_hi(cs, b, n, d, A, xs, ws, K, M) for n in range(N + 1)]
    zst = -(A - 1j * b) ** 2
    absz, th = abs(zst), mp.arg(zst)
    amps = []
    for nlo, nhi in windows:
        rows, rhs = [], []
        for n in range(nlo, nhi + 1):
            y = a[n] * absz ** n / mp.factorial(2 * n) * mp.mpf(n) ** (mp.mpf(5) / 2)
            c, s = mp.cos(n * th), mp.sin(n * th)
            rows.append([c, s, c / n, s / n])
            rhs.append(y)
        Mx = mp.matrix(rows)
        sol = mp.lu_solve(Mx.T * Mx, Mx.T * mp.matrix(rhs))
        amps.append(mp.sqrt(sol[0] ** 2 + sol[1] ** 2))
    pred = (cs + b * b) ** (mp.mpf(3) / 2) / (4 * mp.sqrt(mp.pi) * mp.sqrt(cs) * b)
    return amps, pred


@pytest.mark.slow
def test_paper61_stokes_amplitude_closed_form():
    mp.mp.dps = 40
    cs, b = mp.mpf('0.2'), mp.mpf('0.5')
    windows = [(26, 38), (31, 43), (34, 46)]
    centers = [mp.mpf(lo + hi) / 2 for lo, hi in windows]
    amps, pred = _amp_windows(cs, b, 46, windows)

    ratios = [a / pred for a in amps]
    # (1) monotone increase toward the closed form, from below
    assert amps[0] < amps[1] < amps[2], [mp.nstr(a, 8) for a in amps]
    assert all(r < 1 for r in ratios), [mp.nstr(r, 8) for r in ratios]
    # (2) the top window is already within ~4% and climbing
    assert mp.mpf('0.955') < ratios[-1] < mp.mpf('1.0'), mp.nstr(ratios[-1], 8)

    # (3) 1/n^2 Richardson limit equals the closed form (the "converges to 1" claim,
    #     quantified): fit amp(n_c) = amp_inf - C/n_c^2 on the two widest-separated
    #     centers and check amp_inf/pred -> 1.
    n1, n2 = centers[0], centers[-1]
    a1, a2 = amps[0], amps[-1]
    # amp(n) = amp_inf - C/n^2 with C>0 (amp rises to the limit), so the limit is
    # amp_inf = amp(n2) + C/n2^2 -- extrapolate UP from the top window, not down.
    C = (a2 - a1) / (1 / n1 ** 2 - 1 / n2 ** 2)
    amp_inf = a2 + C / n2 ** 2
    assert amp_inf > amps[-1], mp.nstr(amp_inf, 10)          # extrapolated toward 1
    assert abs(amp_inf / pred - 1) < mp.mpf('0.035'), mp.nstr(amp_inf / pred, 10)

    # (4) NO algebraic factor is consistent: only the factor 1 fits. Each of the
    #     candidate factors would drive amp/pred far from the measured ~0.99.
    fac = ratios[-1]
    for name, val in (('c1+b2', cs + b * b),
                      ('b2/(c1+b2)', b * b / (cs + b * b)),
                      ('c1/(c1+b2)', cs / (cs + b * b))):
        # if the true amplitude carried this extra factor, fac would equal `val`
        assert abs(fac - val) > mp.mpf('0.3'), (name, mp.nstr(fac, 8), mp.nstr(val, 8))
    # ...and fac itself is uniquely near 1
    assert abs(fac - 1) < mp.mpf('0.05'), mp.nstr(fac, 8)
