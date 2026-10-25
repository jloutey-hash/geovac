"""PF1-PF4 (docs/predictions/2026-10-23_price_formula_slots.md, committed b84e1231 before computing): shared helpers.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12.
Conventions of record QG (debug/sprint_ring_ext_QG_memo.md): rho~ = (c k^2 + 1)^(-3/2), scale = inverse length,
f~ = 8 pi gamma/(k^2+gamma^2)^2, Y(t) = (1/(4 pi^2)) t^3 int z^(1/2) (z+t^2)^(-3/2) phi(z) dz, x = t^2, operators in x
with theta = x d/dx in the ga3 Q-form conventions.

Kite (route B), derived from the definition with the lines labelled as in the record (O1 = a, O2 = b, 13 = gamma_1,
23 = gamma_2, 12 = Coulomb).  Momentum space: Y = (2 pi)^-3 int d^3L rho~_t(L) Phi~(L),
  Phi~(L) = (2 pi)^-6 int d^3k d^3q [rho~_a(k) rho~_b(k - L)] [f~_1(q) f~_2(q - L)] 4 pi/|k - q|^2.
Feynman parameter u on the density pair, weighting the mass-a line (rho~_a rho~_b = a^3 b^3 (8/pi) int du
sqrt(u(1-u)) [(k - (1-u)L)^2 + M^2]^-3, M^2 = u a^2 + (1-u) b^2 + u(1-u) z); v on the geminal pair, weighting the
gamma_2 line (f~_1 f~_2 = 384 pi^2 gamma_1 gamma_2 int dv v(1-v) [(q - vL)^2 + N^2]^-4,
N^2 = (1-v) gamma_1^2 + v gamma_2^2 + v(1-v) z).  After the shifts the Coulomb line carries k - q + (1-u-v) L, so the
position-space convolution gives kappa^2 = (1-u-v)^2 z, and
  Phi~ = 2 gamma_1 gamma_2 a^3 b^3 int_0^1 int_0^1 du dv sqrt(u(1-u)) v(1-v) K / (M^3 N^5),
  K = (1/kappa) int_0^oo e^{-s r} (1 + M r)(3 + 3 N r + N^2 r^2) sin(kappa r) dr
    = 3/D + 6 s^2/D^2 + 2 (N^2 + 3 M N)(3 s^2 - k2)/D^3 + 24 M N^2 s (s^2 - k2)/D^4,  s = M + N, D = s^2 + k2.
On the slice a = b the weight sqrt(u(1-u)) is symmetric under u -> 1-u and (1-u-v)^2 <-> (u-v)^2, which is why
record QG's k2 = (u-v)^2 z was correct there; off the slice the assignment matters (assign='wrong' keeps the slice
form with the same M and is the discriminating control).  Route A (the momentum triangle with the closed-form loop
I_221(gamma_1, gamma_2, 0; w^2, p^2, q^2) of ring_ext_QG_common, built directly from the definition) is the arbiter.
Nothing in ring_ext_QG_common / ga1_common is modified; their node tables and loop are imported."""
from __future__ import annotations
import sys
from pathlib import Path
import mpmath as mp
import numpy as np
from flint import arb, acb, ctx as fctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm  # noqa: E402

DATA = Path(__file__).parent / "data"
# the registered configuration
A_SC, B_SC, G1, G2 = "1", "2", "3/2", "3"


def params(a=A_SC, b=B_SC, g1=G1, g2=G2):
    return arb(a), arb(b), arb(g1), arb(g2)


def kernel(M, N, k2):
    """K(M, N, s, k2) of the exp-geminal kite (identical to ring_ext_QG_common.phi_gem_z's K)."""
    N2 = N * N
    s = M + N
    s2 = s * s
    iD = 1 / (s2 + k2)
    iD2 = iD * iD
    return 3 * iD + 6 * s2 * iD2 + 2 * (N2 + 3 * M * N) * (3 * s2 - k2) * iD2 * iD \
        + 24 * M * N2 * s * (s2 - k2) * iD2 * iD2


def phi_z(z, nu, nv, a=A_SC, b=B_SC, g1=G1, g2=G2, assign="correct"):
    """phi(z) = Phi~(sqrt z) of the kite, z complex (acb); Gauss-Chebyshev (2nd kind) in u, Gauss-Legendre in v.
    assign: 'correct' k2 = (1-u-v)^2 z (derived above); 'wrong' k2 = (u-v)^2 z with the same M (slice form)."""
    z = acb(z)
    a, b, g1, g2 = params(a, b, g1, g2)
    a2, b2, g12, g22 = a * a, b * b, g1 * g1, g2 * g2
    Mu = []
    for u, Wu in Cm.cheb2_nodes(nu):
        M = (a2 * u + b2 * (1 - u) + u * (1 - u) * z).sqrt()
        Mu.append((u, Wu / (M * M * M), M))
    Nv = []
    for v, Wv in Cm.gl_nodes01(nv):
        V = v * (1 - v)
        N = (g12 * (1 - v) + g22 * v + V * z).sqrt()
        N2 = N * N
        Nv.append((v, Wv * V / (N2 * N2 * N), N))
    tot = acb(0)
    for v, cv, N in Nv:
        acc = acb(0)
        for u, cu, M in Mu:
            d = (1 - u - v) if assign == "correct" else (u - v)
            acc += cu * kernel(M, N, d * d * z)
        tot += cv * acc
    return 2 * g1 * g2 * a2 * a * b2 * b * tot


def phi_ts(w, h, eps, a=A_SC, b=B_SC, g1=G1, g2=G2):
    """phi at real w with tanh-sinh in u and v (for the Y grid at large real w; w real arb)."""
    z = acb(w) * acb(w)
    a, b, g1, g2 = params(a, b, g1, g2)
    a2, b2, g12, g22 = a * a, b * b, g1 * g1, g2 * g2
    nodes = Cm.ts_nodes01(arb(h), eps)
    Mu = []
    for u, um, Wu in nodes:
        U = u * um
        M = (a2 * u + b2 * um + U * z).sqrt()
        Mu.append((u, Wu * U.sqrt() / (M * M * M), M))
    Nv = []
    for v, vm, Wv in nodes:
        V = v * vm
        N = (g12 * vm + g22 * v + V * z).sqrt()
        N2 = N * N
        Nv.append((v, Wv * V / (N2 * N2 * N), N))
    tot = acb(0)
    for v, cv, N in Nv:
        acc = acb(0)
        for u, cu, M in Mu:
            d = 1 - u - v
            acc += cu * kernel(M, N, d * d * z)
        tot += cv * acc
    return 2 * g1 * g2 * a2 * a * b2 * b * tot


# ---------------------------------------------------------------- route A: the momentum triangle
def phi_routeA(w, a=A_SC, b=B_SC, g1=G1, g2=G2):
    """Phi~(w) = (2 pi)^-3 int d^3k rho~_a(k) rho~_b(k - L) T(k, L),
    T = (2 pi)^-3 (4 pi) (8 pi g1)(8 pi g2) I_221(g1, g2, 0; s12 = w^2, s13 = p^2, s23 = q^2), p = |k|, q = |k - L|;
    the angular integral gives (2 pi / w) int int_{|p-q| <= w <= p+q} p q dp dq.  mpmath quadrature."""
    w = mp.mpf(w)
    a, b, g1, g2 = (mpq(s) for s in (a, b, g1, g2))
    ra = lambda k: a ** 3 * (k * k + a * a) ** mp.mpf(-1.5)
    rb = lambda k: b ** 3 * (k * k + b * b) ** mp.mpf(-1.5)
    pref = (2 * mp.pi) ** -3 * 4 * mp.pi * (8 * mp.pi * g1) * (8 * mp.pi * g2)

    def inner(p):
        return mp.quad(lambda q: q * rb(q) * Cm.I221(g1, g2, 0, w * w, p * p, q * q), [abs(p - w), p + w])
    val = mp.quad(lambda p: p * ra(p) * inner(p), [0, w, 2 * w, 4 * w + 2, mp.inf])
    return (2 * mp.pi) ** -3 * 2 * mp.pi / w * pref * val


# ---------------------------------------------------------------- the loop against direct 3-D quadrature
def direct3d(x1, x2, g1, g2, Nang=160, Nr=240):
    """int d^3 l g_1(l - x1) g_2(l - x2) / |l|^2 (x3 = 0), g_i(q) = (q^2 + g_i^2)^-2, float64 Gauss rules,
    independent of the closed form (copy of ga3_1_checks.direct with two masses)."""
    x1, x2 = np.array(x1), np.array(x2)
    uu, wu = np.polynomial.legendre.leggauss(Nang)
    ph = np.linspace(0, 2 * np.pi, Nang, endpoint=False)
    wph = 2 * np.pi / Nang
    st = np.sqrt(1 - uu**2)
    n = np.stack([np.outer(st, np.cos(ph)), np.outer(st, np.sin(ph)), np.outer(uu, np.ones_like(ph))], -1)
    d1, d2 = n @ x1, n @ x2
    R = max(np.linalg.norm(x1), np.linalg.norm(x2)) + 2
    sr, wr = np.polynomial.legendre.leggauss(Nr)

    def f_at(ri):
        q1 = ri**2 - 2 * ri * d1 + x1 @ x1
        q2 = ri**2 - 2 * ri * d2 + x2 @ x2
        return (q1 + g1 * g1) ** -2 * (q2 + g2 * g2) ** -2
    tot = 0.0
    for lo, hi in [(0, R / 3), (R / 3, 2 * R / 3), (2 * R / 3, R)]:
        r = (hi - lo) / 2 * sr + (hi + lo) / 2
        wgt = (hi - lo) / 2 * wr
        for ri, wi in zip(r, wgt):
            tot += wi * np.sum(wu[:, None] * f_at(ri)) * wph
    s = (sr + 1) / 2
    ws = wr / 2
    for si, wi in zip(s, ws):
        ri = R / si
        tot += wi * R / si**2 * np.sum(wu[:, None] * f_at(ri)) * wph
    return tot


def mpq(s):
    """'3/2' -> mpf"""
    s = str(s)
    if "/" in s:
        n, d = s.split("/")
        return mp.mpf(n) / mp.mpf(d)
    return mp.mpf(s)
