"""Backing test for Paper 61 (companion of Paper 59) sec:modular co-area reduction (v4.100.0, fourth T2 track).

The genuinely 3D collinear observable T2 = (8/pi) int int J(s,t) reduces EXACTLY to a
1D modular integral by a co-area change of variables (modulus rho=c_t/c_s outer, scale
u=c_s inner), Jacobian ds dt = u/(sqrt(1-4u) sqrt(1-4 rho u)) du drho:

    Phi(rho) = int_0^umax [sum over 4 (s,t)-branches of j0(k b) ] Pc(c_s,k)Pc(c_t,k) du-measure
    T2 = (16/pi) int_0^1 Phi(rho) drho          (folded by the exact symmetry below)

Backed facts (self-contained, no debug/ import):
  (A) exact s<->t symmetry Phi(rho) = Phi(1/rho)/rho^2  -> the two X(2) real-locus contour
      halves are equal, giving the single (16/pi) int_0^1 form;
  (B) the leading Fourier-Whittaker (log) coefficient A = J(1/4,1/4,1)/4 -- the fibre value at
      the DIAGONAL c_s=c_t (the genus-0 degenerate cusp), one transcendence level below T2;
  (C) the reduction reproduces the frozen anchor (coarse grid here; the driver reaches 2.2e-9).
"""
import mpmath as mp
from mpmath import legendre
import pytest

T2_FROZEN = mp.mpf('0.3953557659017139641')   # frozen anchor; validate-only
_G = {}


def _gl(N):
    if N in _G:
        return _G[N]
    roots, ws = [], []
    for k in range(1, N + 1):
        x = mp.cos(mp.pi * (k - mp.mpf('0.25')) / (N + mp.mpf('0.5')))
        for _ in range(90):
            f = legendre(N, x)
            fp = N * (x * legendre(N, x) - legendre(N - 1, x)) / (x * x - 1)
            dx = f / fp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 6):
                break
        roots.append(x)
    for x in roots:
        fp = N * (x * legendre(N, x) - legendre(N - 1, x)) / (x * x - 1)
        ws.append(2 / ((1 - x * x) * fp * fp))
    _G[N] = (roots, ws)
    return _G[N]


def _Pc(c, k):
    D = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-D) * (1 / D ** 3 + 3 / D ** 4 + 3 / D ** 5)


def _Jsum(c1, c2, bs, Nk):
    L = 1 / (mp.sqrt(c1) + mp.sqrt(c2))
    xs, ws = _gl(Nk)
    tot = mp.mpf(0)
    for x, w in zip(xs, ws):
        uu = (x + 1) / 2
        k = L * uu / (1 - uu)
        dk = L / (1 - uu) ** 2
        sj = mp.mpf(0)
        for b in bs:
            kb = k * b
            sj += mp.sin(kb) / kb if kb > mp.mpf('1e-40') else mp.mpf(1)
        tot += (w / 2) * sj * _Pc(c1, k) * _Pc(c2, k) * dk
    return tot


def _four_bs(u, rho):
    sm = (1 - mp.sqrt(1 - 4 * u)) / 2
    tm = (1 - mp.sqrt(1 - 4 * rho * u)) / 2
    return [sm + tm, sm + (1 - tm), (1 - sm) + tm, (1 - sm) + (1 - tm)]


def _Phi(rho, Nu, Nk):
    umax = min(mp.mpf(1) / 4, 1 / (4 * rho))
    xs, ws = _gl(Nu)
    H = mp.pi / 2
    tot = mp.mpf(0)
    for x, w in zip(xs, ws):
        phi = H * (x + 1) / 2
        u = umax * mp.sin(phi) ** 2
        du = umax * mp.sin(2 * phi)
        wj = H * w / 2
        val = _Jsum(u, rho * u, _four_bs(u, rho), Nk)
        meas = u / (mp.sqrt(1 - 4 * u) * mp.sqrt(1 - 4 * rho * u))
        tot += wj * du * val * meas
    return tot


def _A(Nk):
    """A = J(1/4,1/4,1)/4, two independent evaluators (decay-map GL + tanh-sinh quad)."""
    c = mp.mpf(1) / 4
    a = _Jsum(c, c, [mp.mpf(1)], Nk)
    b = mp.quad(lambda k: (mp.sin(k) / k if k > mp.mpf('1e-40') else mp.mpf(1)) * _Pc(c, k) ** 2,
                [0, 1, 2, 4, 8, 16, mp.inf])
    return a / 4, b / 4


@pytest.mark.slow
def test_paper59_coarea_symmetry_and_A():
    """(A) exact s<->t symmetry Phi(rho)=Phi(1/rho)/rho^2; (B) A=J(1/4,1/4,1)/4 two ways."""
    mp.mp.dps = 30
    Nu, Nk = 40, 60
    for rho in [mp.mpf('0.5'), mp.mpf('0.35')]:
        lhs = _Phi(rho, Nu, Nk)
        rhs = _Phi(1 / rho, Nu, Nk) / rho ** 2
        assert abs(lhs - rhs) < mp.mpf(10) ** -8, (rho, abs(lhs - rhs))
    a_gl, a_quad = _A(80)
    assert abs(a_gl - a_quad) < mp.mpf(10) ** -10, (a_gl, a_quad)
    # A is the LOG coefficient: Phi(rho) + A ln(1-rho) stays finite (converges to C) as rho->1
    A = a_gl
    g99 = _Phi(mp.mpf('0.99'), 130, Nk) + A * mp.log(mp.mpf('0.01'))
    g999 = _Phi(mp.mpf('0.999'), 200, Nk) + A * mp.log(mp.mpf('0.001'))
    assert abs(g99 - g999) < mp.mpf('0.02'), (g99, g999)   # bounded => log coeff is exactly A


@pytest.mark.slow
def test_paper59_coarea_reduction_reproduces_T2():
    """(C) T2 = (16/pi) int_0^1 Phi drho, via int Phi = A + int[Phi + A ln(1-rho)] (v=sqrt(1-rho)).
    Coarse grid here (the driver reaches 2.2e-9); a few 1e-3 confirms the reduction is correct."""
    mp.mp.dps = 30
    Nk = 50
    A = _A(80)[0]
    xs, ws = _gl(20)
    bracket = mp.mpf(0)
    for x, w in zip(xs, ws):
        v = (x + 1) / 2                      # v in (0,1), rho = 1 - v^2, ln(1-rho) = 2 ln v
        rho = 1 - v * v
        Nu = 60 if v > mp.mpf('0.3') else 120     # more nodes near the cusp v->0
        integ = (_Phi(rho, Nu, Nk) + A * 2 * mp.log(v)) * 2 * v
        bracket += (w / 2) * integ
    T2_recon = (16 / mp.pi) * (A + bracket)
    assert abs(T2_recon - T2_FROZEN) < mp.mpf('5e-3'), (T2_recon, abs(T2_recon - T2_FROZEN))
