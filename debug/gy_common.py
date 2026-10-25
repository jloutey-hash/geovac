"""GY1-GY3 (docs/predictions/2026-10-23_gaussian_geminal_triangle.md): shared helpers for the Gaussian geminal
g(r) = exp(-gamma r^2) on both edges of the correlated triangle (apex placement t on k, scales (1, 1, t),
fingerprint fibre rho~ = (k^2+1)^(-3/2)).  Conventions of record QG: F~(k) = int exp(-ik.x) F d^3x,
g~(k) = (pi/gamma)^(3/2) exp(-k^2/(4 gamma)).  Nothing in the lineage (ring_ext_QG_common, ga1_common) is modified.

Derivation (checked symbolically and numerically in debug/gy_0_kernel.py):
 * Gaussian pair:  g~(a) g~(a+L) = (pi/gamma)^3 exp(-(a + L/2)^2/(2 gamma)) exp(-L^2/(8 gamma)), no Feynman parameter.
 * Loop (route A): int d^3l g(l-x1) g(l-x2) / |l-x3|^2 with g(q) = exp(-q^2/(4 gamma)) equals
       exp(-s12/(8 gamma)) * 4 gamma pi^(3/2) Daw(R/sqrt(2 gamma)) / R,   R^2 = (s13 + s23)/2 - s12/4,
   Daw the Dawson function (from 1/q^2 = int_0^oo dlam e^{-lam q^2} and a Gaussian in l).
 * Kite (route B, convolution theorem of ga1_0_kernel.py with the Coulomb line 1/r in position space):
       Phi~(w) = sqrt(pi/8) gamma^(-3/2) exp(-z/(8 gamma)) int_0^1 du sqrt(u(1-u)) M^-3 K_G(M, kappa),
       M = sqrt(1 + u(1-u) z), kappa^2 = (u - 1/2)^2 z, z = w^2,
       K_G = (1/kappa) int_0^oo e^{-M r - gamma r^2/2} (1 + M r) sin(kappa r) dr
           = [ H(b-) - H(b+) ] / (2 i kappa),   b+- = M +- i kappa,   H(b) = I(b) (1 - M b/gamma),
       I(b) = int_0^oo e^{-gamma r^2/2 - b r} dr = sqrt(pi/(2 gamma)) e^{b^2/(2 gamma)} erfc(b/sqrt(2 gamma)),
   using int r e^{-gamma r^2/2 - b r} dr = (1 - b I(b))/gamma.
 * Split of the transcendental content: with E = I(1), P(b) = e^{(b^2-1)/(2 gamma)}, Qh(b) = int_1^b e^{-(y^2-1)/(2 gamma)} dy,
       I(b) = E P(b) - (1/gamma) P(b) Qh(b),
   so Phi~ = E Phi_a + Phi_b with Phi_a, Phi_b in pi^(3/2) Q[[z]] (E = sqrt(pi/2) e^{1/2} erfc(1/sqrt 2) at gamma = 1).
   An operator over Q(z) annihilating Phi~ annihilates Phi_a and Phi_b separately (E irrational).
 * Y(t) = (1/(2 pi^2)) int_0^oo w^2 Phi~(w) t^3 (w^2 + t^2)^(-3/2) dw   (ring_ext_QG_3_values.Y_theta, kernel shape-independent)."""
from __future__ import annotations
import sys
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx as fctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm


# ---------------------------------------------------------------- route A: the Gaussian loop (mpmath)
def dawson(x):
    """Dawson function D(x) = e^{-x^2} int_0^x e^{y^2} dy = x 1F1(1; 3/2; -x^2) (mpmath, any precision)."""
    x = mp.mpf(x)
    return x * mp.hyp1f1(1, mp.mpf(3) / 2, -x * x)


def loop_gauss(g, s12, s13, s23):
    """int d^3l exp(-|l-x1|^2/(4g)) exp(-|l-x2|^2/(4g)) / |l-x3|^2, squared distances s_ij."""
    g = mp.mpf(g)
    R2 = (mp.mpf(s13) + mp.mpf(s23)) / 2 - mp.mpf(s12) / 4
    pre = mp.exp(-mp.mpf(s12) / (8 * g)) * 4 * g * mp.pi ** mp.mpf(1.5)
    if R2 <= 0:
        # R = 0: Daw(x)/x -> 1 / sqrt(2g) ... D(R/sqrt(2g))/R -> 1/sqrt(2g)
        return pre / mp.sqrt(2 * g)
    R = mp.sqrt(R2)
    return pre * dawson(R / mp.sqrt(2 * g)) / R


def phi_routeA(w, g):
    """Phi~(w) = (2pi)^-3 (2pi/w) int int_{|p-q|<=w<=p+q} p q rho~(p) rho~(q) T(p,q,w),
    T = (2pi)^-3 (4 pi) (pi/g)^3 loop_gauss(w^2, p^2, q^2)   (same skeleton as ga1_common.phi_routeA)."""
    w = mp.mpf(w)
    g = mp.mpf(g)
    rt = lambda k: (k * k + 1) ** mp.mpf(-1.5)
    pref = (2 * mp.pi) ** -3 * (mp.pi / g) ** 3 * 4 * mp.pi

    def inner(p):
        return mp.quad(lambda q: q * rt(q) * loop_gauss(g, w * w, p * p, q * q), [abs(p - w), p + w])
    val = mp.quad(lambda p: p * rt(p) * inner(p), [0, w / 2, w, 2 * w, 4 * w + 2, mp.inf])
    return (2 * mp.pi) ** -3 * 2 * mp.pi / w * pref * val


# ---------------------------------------------------------------- route B: the kite (flint acb)
def I_fun(b, g):
    """I(b) = sqrt(pi/(2g)) e^{b^2/(2g)} erfc(b/sqrt(2g)), b acb, g arb"""
    return (arb.pi() / (2 * g)).sqrt() * (b * b / (2 * g)).exp() * (b / (2 * g).sqrt()).erfc()


def P_fun(b, g):
    return ((b * b - 1) / (2 * g)).exp()


def Qh_fun(b, g):
    """Qh(b) = int_1^b e^{-(y^2-1)/(2g)} dy = e^{1/(2g)} sqrt(g pi/2) [erf(b/sqrt(2g)) - erf(1/sqrt(2g))]"""
    s2g = (2 * g).sqrt()
    return (1 / (2 * g)).exp() * (g * arb.pi() / 2).sqrt() * ((b / s2g).erf() - (1 / s2g).erf())


def E_const(g):
    return I_fun(acb(1), g).real


def KG_parts(M, kap2, Uz, g, which, kap=None):
    """K_G built from a given 'G-function': which = 'full' (I), 'a' (P), 'b' (-(1/g) P Qh).
    K_G = [H(b-) - H(b+)]/(2 i kappa), H(b) = Gf(b) (1 - M b/g); even in kappa, kappa^2 = kap2 (acb).
    At kappa = 0 the limit is used: K_G(M, 0) = int (1 + M r) r e^{-M r - g r^2/2} dr."""
    def Gf(b):
        if which == "full":
            return I_fun(b, g)
        if which == "a":
            return P_fun(b, g)
        return -(1 / g) * P_fun(b, g) * Qh_fun(b, g)
    if kap2 == 0:
        # kappa -> 0: K_G(M,0) = J1(M) + M J2(M), J1 = -d/db Gf, J2 = d^2/db^2 Gf  (Gf' = (b/g) Gf - c, c = 1/g for I, else 0)
        # use a symmetric finite step in the analytic variable? no: exact via the ODE each piece satisfies:
        # Gf' = (b/g) Gf - c with c = 1/g ('full', 'b'), 0 ('a').
        c = 0 if which == "a" else 1 / g
        G0 = Gf(acb(M))
        G1 = M / g * G0 - c
        G2 = G0 / g + M / g * G1
        return -G1 + M * G2
    if kap is None:
        kap = acb(kap2).sqrt()
    bp = M + acb(0, 1) * kap
    bm = M - acb(0, 1) * kap
    Hp = Gf(bp) * (1 - M * bp / g)
    Hm = Gf(bm) * (1 - M * bm / g)
    return (Hm - Hp) / (2 * acb(0, 1) * kap)


def _phi_nodes(nodes, z, g, which):
    """sum over (u, 1-u, weight) nodes of weight * M^-3 K_G; 'nodes' yields (u, um, W) with um = 1 - u."""
    half = arb(1) / 2
    tot = acb(0)
    zs = z.sqrt()
    for u, um, Wt in nodes:
        U = u * um
        Uz = U * z
        M = (1 + Uz).sqrt()
        W = (u - half) * (u - half)
        if W == 0:
            tot += Wt / (M * M * M) * KG_parts(M, 0, Uz, g, which)
        else:
            tot += Wt / (M * M * M) * KG_parts(M, W * z, Uz, g, which, kap=W.sqrt() * zs)
    return tot


def prefactor(z, g):
    """sqrt(pi/8) g^(-3/2) exp(-z/(8g))"""
    return (arb.pi() / 8).sqrt() * g ** (arb(-3) / 2) * (-z / (8 * g)).exp()


def phi_z(z, g, nu, which="full"):
    """phi(z) = Phi~(sqrt z), z complex (acb), gamma = g (arb); Gauss-Chebyshev (2nd kind) in u, nu even
    (no node at u = 1/2).  which = 'full' | 'a' | 'b' (the E-split pieces, phi = E phi_a + phi_b)."""
    z = acb(z)
    g = arb(g)
    assert nu % 2 == 0
    nodes = [(u, 1 - u, Wu) for u, Wu in Cm.cheb2_nodes(nu)]
    return prefactor(z, g) * _phi_nodes(nodes, z, g, which)


def phi_ts(w, g, h="1/16", eps=40, which="full"):
    """phi at real w with tanh-sinh in u (for the Y grid, large real w)."""
    z = acb(w) * acb(w)
    g = arb(g)
    nodes = [(u, um, Wt * (u * um).sqrt()) for u, um, Wt in Cm.ts_nodes01(arb(h), eps)]
    return prefactor(z, g) * _phi_nodes(nodes, z, g, which)


def phi_all_z(z, g, nu):
    """(phi, phi_a, phi_b) at complex z in one pass over the Gauss-Chebyshev nodes (nu even)."""
    z = acb(z)
    g = arb(g)
    assert nu % 2 == 0
    half = arb(1) / 2
    tots = [acb(0), acb(0), acb(0)]
    zs = z.sqrt()
    for u, Wu in Cm.cheb2_nodes(nu):
        U = u * (1 - u)
        Uz = U * z
        M = (1 + Uz).sqrt()
        W = (u - half) * (u - half)
        cu = Wu / (M * M * M)
        kap = W.sqrt() * zs
        for i, which in enumerate(("full", "a", "b")):
            tots[i] += cu * KG_parts(M, W * z, Uz, g, which, kap=kap)
    pre = prefactor(z, g)
    return [pre * t for t in tots]


def KG_direct(M, kappa, g, dps=40):
    """K_G by direct r-quadrature (mpmath), real M, kappa (kappa may be imaginary: pass mpc)."""
    mp.mp.dps = dps
    M, g = mp.mpf(M), mp.mpf(g)
    f = lambda r: mp.exp(-M * r - g * r * r / 2) * (1 + M * r) * (mp.sin(kappa * r) / kappa if kappa != 0 else r)
    return mp.quad(f, [0, 1, 2, 4, 8, 16, mp.inf])
