"""Shared helpers for test QG of docs/predictions/2026-10-17_ring_extensions.md.

The 3-D one-loop triangle with masses: closed form by the solid-angle identity
  int_simplex (alpha^T M alpha)^(-3/2) d alpha = Omega / sqrt(det M),
  tan(Omega/2) = sqrt(det M) / (m1 m2 m3 + M12 m3 + M13 m2 + M23 m1),
M_ij = (m_i^2 + m_j^2 + s_ij)/2 (M_ii = m_i^2), s_ij = |x_i - x_j|^2, so
  I_111 = int d^3 l / prod_i (|l - x_i|^2 + m_i^2) = pi^2 arctan(sqrt(det M)/D) / sqrt(det M).
I_221 (propagators 1 and 2 squared) = d/d(m1^2) d/d(m2^2) I_111, by a bi-dual number.
All functions take mpmath numbers (or flint arb via duck typing where noted).
"""
from __future__ import annotations

import mpmath as mp


class BD:
    """Bi-dual number v + a e1 + b e2 + c e1 e2 (e1^2 = e2^2 = 0)."""
    __slots__ = ("v", "a", "b", "c")

    def __init__(self, v, a=0, b=0, c=0):
        self.v, self.a, self.b, self.c = v, a, b, c

    @staticmethod
    def lift(x):
        return x if isinstance(x, BD) else BD(x)

    def __add__(self, o):
        o = BD.lift(o)
        return BD(self.v + o.v, self.a + o.a, self.b + o.b, self.c + o.c)
    __radd__ = __add__

    def __neg__(self):
        return BD(-self.v, -self.a, -self.b, -self.c)

    def __sub__(self, o):
        return self + (-BD.lift(o))

    def __rsub__(self, o):
        return BD.lift(o) + (-self)

    def __mul__(self, o):
        o = BD.lift(o)
        return BD(self.v * o.v, self.v * o.a + self.a * o.v, self.v * o.b + self.b * o.v,
                  self.v * o.c + self.a * o.b + self.b * o.a + self.c * o.v)
    __rmul__ = __mul__

    def fn(self, f0, f1, f2):
        """f(self) given f, f', f'' at self.v."""
        return BD(f0, f1 * self.a, f1 * self.b, f1 * self.c + f2 * self.a * self.b)

    def inv(self):
        iv = 1 / self.v
        return self.fn(iv, -iv * iv, 2 * iv * iv * iv)

    def __truediv__(self, o):
        return self * BD.lift(o).inv()

    def __rtruediv__(self, o):
        return BD.lift(o) * self.inv()

    def sqrt(self):
        s = mp.sqrt(self.v)
        return self.fn(s, 1 / (2 * s), -1 / (4 * s * self.v))


def A_fun(y):
    """A(y) = arctan(sqrt y)/sqrt y and its first two derivatives (analytic in y > -1)."""
    if abs(y) < mp.mpf("0.05"):
        a0 = a1 = a2 = mp.mpf(0)
        k = 0
        while True:
            term = (-y) ** k / (2 * k + 1)
            a0 += term
            if k >= 1:
                a1 += k * (-1) ** k * y ** (k - 1) / (2 * k + 1)
            if k >= 2:
                a2 += k * (k - 1) * (-1) ** k * y ** (k - 2) / (2 * k + 1)
            if k > 4 and abs(term) < mp.eps * mp.mpf(10) ** -5:
                break
            k += 1
        return a0, a1, a2
    if y > 0:
        s = mp.sqrt(y)
        a0 = mp.atan(s) / s
    else:
        s = mp.sqrt(-y)
        a0 = mp.atanh(s) / s
    a1 = (1 / (1 + y) - a0) / (2 * y)
    a2 = -1 / (2 * y * (1 + y) ** 2) - 3 * a1 / (2 * y)
    return a0, a1, a2


def tri_from_sq(a, b, c, s12, s13, s23):
    """I_111 as a function of squared masses a = m1^2, b = m2^2, c = m3^2 (any may be BD)."""
    A_, B_, C_ = BD.lift(a), BD.lift(b), BD.lift(c)
    m1, m2 = A_.sqrt(), B_.sqrt()
    m3 = BD.lift(mp.sqrt(c)) if not isinstance(c, BD) else C_.sqrt()
    M12 = (A_ + B_ + s12) * mp.mpf(0.5)
    M13 = (A_ + C_ + s13) * mp.mpf(0.5)
    M23 = (B_ + C_ + s23) * mp.mpf(0.5)
    det = A_ * (B_ * C_ - M23 * M23) - M12 * (M12 * C_ - M23 * M13) + M13 * (M12 * M23 - B_ * M13)
    D = m1 * m2 * m3 + M12 * m3 + M13 * m2 + M23 * m1
    y = det / (D * D)
    Ay = y.fn(*A_fun(y.v))
    return mp.pi ** 2 * Ay / D


def I111(m1, m2, m3, s12, s13, s23):
    return tri_from_sq(mp.mpf(m1) ** 2, mp.mpf(m2) ** 2, mp.mpf(m3) ** 2, s12, s13, s23).v


def I221(m1, m2, m3, s12, s13, s23):
    """int d^3l / [(|l-x1|^2+m1^2)^2 (|l-x2|^2+m2^2)^2 (|l-x3|^2+m3^2)] = d_a d_b I_111."""
    a = BD(mp.mpf(m1) ** 2, 1, 0, 0)
    b = BD(mp.mpf(m2) ** 2, 0, 1, 0)
    return tri_from_sq(a, b, mp.mpf(m3) ** 2, s12, s13, s23).c


# ---------------------------------------------------------------------------------------------
# Double-bubble form (route B).  Fourier convention F~(k) = int e^{-ik.x} F(x) d^3x.
# rho~ = (k^2+1)^(-3/2) (normalized density, scale 1), f = exp(-gamma r), f~ = 8 pi gamma/(k^2+gamma^2)^2.
# Apex t: Y(t) = (2 pi)^-3 int d^3L rho~_k(L;t) Phi~(L) = (1/(2 pi^2)) int_0^inf w^2 t^3 (w^2+t^2)^(-3/2) Phi~(w) dw,
# Phi~(L) = (2pi)^-6 int int [f~(a) f~(a+L)] [rho~(k) rho~(k+L)] (4 pi/|k-a|^2)
#        = 2 gamma^2 int_0^1 int_0^1 du dv sqrt(u(1-u)) v(1-v) K / (M^3 N^5),
# M = sqrt(1 + u(1-u) w^2), N = sqrt(gamma^2 + v(1-v) w^2), s = M + N, k2 = (u-v)^2 w^2,
# K = (1/kappa) int_0^inf e^{-s r} (1+Mr)(3+3Nr+N^2 r^2) sin(kappa r) dr
#   = 3/D + 6 s^2/D^2 + 2(N^2+3MN)(3s^2-k2)/D^3 + 24 M N^2 s (s^2-k2)/D^4,  D = s^2 + k2.
# All-Yukawa control (e^{-gamma r}/r on all three edges):
# Phi~_Y = 2 int int du dv sqrt(u(1-u)) /(M^3 N) [1/D + 2 M s/D^2], N = sqrt(gamma^2 + v(1-v) w^2),
#          s = M + N + gamma.
# ---------------------------------------------------------------------------------------------
from flint import arb, acb, ctx as fctx  # noqa: E402

_NODES = {}


def cheb2_nodes(n):
    """int_0^1 g(u) sqrt(u(1-u)) du ~ sum W_k g(u_k)  (Gauss-Chebyshev 2nd kind)."""
    key = ("c2", n, fctx.prec)
    if key not in _NODES:
        pi = arb.pi()
        out = []
        for k in range(1, n + 1):
            th = pi * k / (n + 1)
            x = th.cos()
            out.append(((1 - x) / 2, pi / (n + 1) * th.sin() ** 2 / 4))
        _NODES[key] = out
    return _NODES[key]


def gl_nodes01(n):
    key = ("gl", n, fctx.prec)
    if key not in _NODES:
        out = []
        for k in range(n):
            x, wt = arb.legendre_p_root(n, k, weight=True)
            out.append(((1 - x) / 2, wt / 2))
        _NODES[key] = out
    return _NODES[key]


def phi_gem(w, g, nu=120, nv=120):
    """Phi~(w) of the exponential-geminal kite (apex), w real or complex (acb), gamma = g (arb)."""
    z = acb(w) * acb(w)
    g = arb(g)
    us = cheb2_nodes(nu)
    vs = gl_nodes01(nv)
    Mu = []
    for u, Wu in us:
        M = (1 + u * (1 - u) * z).sqrt()
        Mu.append((u, Wu / (M * M * M), M))
    Nv = []
    for v, Wv in vs:
        V = v * (1 - v)
        N = (g * g + V * z).sqrt()
        N2 = N * N
        Nv.append((v, Wv * V / (N2 * N2 * N), N, N2))
    tot = acb(0)
    for u, cu, M in Mu:
        for v, cv, N, N2 in Nv:
            s = M + N
            k2 = (u - v) * (u - v) * z
            s2 = s * s
            iD = 1 / (s2 + k2)
            iD2 = iD * iD
            K = 3 * iD + 6 * s2 * iD2 + 2 * (N2 + 3 * M * N) * (3 * s2 - k2) * iD2 * iD \
                + 24 * M * N2 * s * (s2 - k2) * iD2 * iD2
            tot += cu * cv * K
    return 2 * g * g * tot


def phi_yuk3(w, g, nu=120, nv=120):
    """Phi~_Y(w): Yukawa e^{-g r}/r on all three edges (apex), GL in v (use small g with care)."""
    z = acb(w) * acb(w)
    g = arb(g)
    tot = acb(0)
    Mu = [(u, Wu, (1 + u * (1 - u) * z).sqrt()) for u, Wu in cheb2_nodes(nu)]
    for v, Wv in gl_nodes01(nv):
        N = (g * g + v * (1 - v) * z).sqrt()
        for u, Wu, M in Mu:
            s = M + N + g
            iD = 1 / (s * s + (u - v) * (u - v) * z)
            tot += Wu * Wv / (M * M * M * N) * (iD + 2 * M * s * iD * iD)
    return 2 * tot


def to_mp(x):
    """flint arb/acb (midpoint) -> mpmath mpf/mpc at the current mp precision."""
    if isinstance(x, acb):
        re = mp.mpf(x.real.mid().str(fctx.dps + 5, radius=False))
        im = mp.mpf(x.imag.mid().str(fctx.dps + 5, radius=False))
        return mp.mpc(re, im) if im != 0 else re
    return mp.mpf(x.mid().str(fctx.dps + 5, radius=False))


def ts_nodes01(h, eps_digits):
    """Tanh-sinh nodes on (0,1): (u, 1-u, weight), truncated when the weight < 10^-eps_digits."""
    key = ("ts", str(h), eps_digits, fctx.prec)
    if key in _NODES:
        return _NODES[key]
    h = arb(h)
    pi = arb.pi()
    out = []
    tiny = arb(10) ** (-eps_digits)
    k = 0
    while True:
        tau = h * k
        sh = tau.sinh()
        e = (pi * sh).exp()
        um = 1 / (1 + e)              # 1 - u for tau >= 0
        up = e / (1 + e)
        wt = pi * h * tau.cosh() * e / (1 + e) ** 2
        if k == 0:
            out.append((up, um, wt))
        else:
            out.append((up, um, wt))
            out.append((um, up, wt))
        if wt < tiny or um < tiny ** 2:
            break
        k += 1
    _NODES[key] = out
    return out


def phi_gem_ts(w, g, h="1/16", eps=40):
    """phi_gem with tanh-sinh in u and v (robust for large real w)."""
    z = acb(w) * acb(w)
    g = arb(g)
    nodes = ts_nodes01(arb(h) if isinstance(h, str) is False else arb(h), eps)
    Mu = []
    for u, um, Wu in nodes:
        U = u * um
        M = (1 + U * z).sqrt()
        Mu.append((u, Wu * U.sqrt() / (M * M * M), M))
    Nv = []
    for v, vm, Wv in nodes:
        V = v * vm
        N = (g * g + V * z).sqrt()
        N2 = N * N
        Nv.append((v, Wv * V / (N2 * N2 * N), N, N2))
    tot = acb(0)
    for u, cu, M in Mu:
        for v, cv, N, N2 in Nv:
            s = M + N
            k2 = (u - v) * (u - v) * z
            s2 = s * s
            iD = 1 / (s2 + k2)
            iD2 = iD * iD
            K = 3 * iD + 6 * s2 * iD2 + 2 * (N2 + 3 * M * N) * (3 * s2 - k2) * iD2 * iD \
                + 24 * M * N2 * s * (s2 - k2) * iD2 * iD2
            tot += cu * cv * K
    return 2 * g * g * tot


def phi_yuk3_ts(w, g, h="1/16", eps=40):
    z = acb(w) * acb(w)
    g = arb(g)
    nodes = ts_nodes01(arb(h), eps)
    Mu = []
    for u, um, Wu in nodes:
        U = u * um
        M = (1 + U * z).sqrt()
        Mu.append((u, Wu * U.sqrt() / (M * M * M), M))
    tot = acb(0)
    for v, vm, Wv in nodes:
        N = (g * g + v * vm * z).sqrt()
        for u, cu, M in Mu:
            s = M + N + g
            iD = 1 / (s * s + (u - v) * (u - v) * z)
            tot += cu * Wv / N * (iD + 2 * M * s * iD * iD)
    return 2 * tot


def phi_gem_z(z, g, nu, nv):
    """phi(z) = Phi~(sqrt z) for the geminal kite (apex), z complex (acb)."""
    z = acb(z)
    g = arb(g)
    Mu = []
    for u, Wu in cheb2_nodes(nu):
        M = (1 + u * (1 - u) * z).sqrt()
        Mu.append((u, Wu / (M * M * M), M))
    Nv = []
    for v, Wv in gl_nodes01(nv):
        V = v * (1 - v)
        N = (g * g + V * z).sqrt()
        N2 = N * N
        Nv.append((v, Wv * V / (N2 * N2 * N), N, N2))
    tot = acb(0)
    for u, cu, M in Mu:
        for v, cv, N, N2 in Nv:
            s = M + N
            k2 = (u - v) * (u - v) * z
            s2 = s * s
            iD = 1 / (s2 + k2)
            iD2 = iD * iD
            K = 3 * iD + 6 * s2 * iD2 + 2 * (N2 + 3 * M * N) * (3 * s2 - k2) * iD2 * iD \
                + 24 * M * N2 * s * (s2 - k2) * iD2 * iD2
            tot += cu * cv * K
    return 2 * g * g * tot


def phi_coul_bubble(w, nu, nphi):
    """gamma = 0 member of the all-Yukawa family (Coulomb on all three edges), v = (1-cos phi)/2:
    Phi~_C(w) = (2/w) int_0^pi dphi int_0^1 du sqrt(u(1-u)) M^-3 [1/D + 2 M s/D^2], N = w sin(phi)/2."""
    w = arb(w)
    z = w * w
    pi = arb.pi()
    tot = arb(0)
    Mu = [(u, Wu, (1 + u * (1 - u) * z).sqrt()) for u, Wu in cheb2_nodes(nu)]
    for x, wx in gl_nodes01(nphi):
        ph = pi * x
        v = (1 - ph.cos()) / 2
        N = w * ph.sin() / 2
        for u, Wu, M in Mu:
            s = M + N
            iD = 1 / (s * s + (u - v) * (u - v) * z)
            tot += Wu * wx * pi / (M * M * M) * (iD + 2 * M * s * iD * iD)
    return 2 * tot / w


def phi1_gem_z(z, g, nu, ns):
    """Secondary placement (t on density i: the vertex with f13 and the Coulomb kernel 12).
    Kite with bubbles {rho_j, rho_k} (scales 1, 1) and {v12, f13}, middle line f23:
    Phi~_1 = 2 gamma int du sqrt(u(1-u)) M^-3 int_0^1 2 dsigma Nt^-3 K1,   v = sigma^2,
    Nt = sqrt(gamma^2 + (1 - v) z), N = sigma Nt, s = M + N + gamma, D = s^2 + (u-v)^2 z,
    K1 = 2s/D^2 + 2(M+N)(3s^2 - k2)/D^3 + 24 M N s (s^2 - k2)/D^4."""
    z = acb(z)
    g = arb(g)
    Mu = [(u, Wu / ((1 + u * (1 - u) * z).sqrt()) ** 3, (1 + u * (1 - u) * z).sqrt()) for u, Wu in cheb2_nodes(nu)]
    Sv = []
    for sg, Ws in gl_nodes01(ns):
        v = sg * sg
        Nt = (g * g + (1 - v) * z).sqrt()
        Sv.append((v, 2 * Ws / (Nt * Nt * Nt), sg * Nt))
    tot = acb(0)
    for u, cu, M in Mu:
        for v, cv, N in Sv:
            s = M + N + g
            k2 = (u - v) * (u - v) * z
            s2 = s * s
            iD = 1 / (s2 + k2)
            iD2 = iD * iD
            K1 = 2 * s * iD2 + 2 * (M + N) * (3 * s2 - k2) * iD2 * iD + 24 * M * N * s * (s2 - k2) * iD2 * iD2
            tot += cu * cv * K1
    return 2 * g * tot
