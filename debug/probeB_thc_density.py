"""Probe B -- ANALYTIC momentum transition densities for exponential-type s-orbitals.

The object under test is the momentum factorization of the ERI (Paper 59
`sec:f12` kernel-swap identity; the same object as the three-centre momentum
reduction in ``debug/routeC_momentum_poc.py``):

    (pq|rs) = (1/(2 pi)^3) int d^3k (4 pi / k^2) rho~_pq(k)^* rho~_rs(k)
    rho~_pq(k) = int phi_p(r) phi_q(r) e^{-i k.r} d^3r

For hydrogenic s-orbitals the transition density has an EXACT, fit-free analytic
representation obtained by the Yukawa/Feynman reduction (no Gaussian fit, no
Barnett--Coulson L-truncation).  Writing

    e^{-a r_A} = -d/da ( e^{-a r_A} / r_A ),   FT[e^{-a r_A}/r_A] = 4 pi e^{-i k.A}/(k^2+a^2)

the convolution theorem plus one Feynman parameter gives, for the UNNORMALISED
1s x 1s product e^{-a r_A} e^{-b r_B},

    rho~^{(00)}(k) = 2 pi  d^2/da db  int_0^1 dt  e^{-i k.P(t)}  e^{-D Delta}/Delta
    Delta(t,k) = sqrt( t(1-t) k^2 + t a^2 + (1-t) b^2 ),   D = |A-B|
    P(t)       = (1-t) A + t B                    (Feynman-interpolated centre)

and the a,b derivatives are done in CLOSED FORM (sympy, once).  Higher radial
powers (the hydrogenic 2s node) come from extra a/b derivatives:

    int r_A^{m1} r_B^{m2} e^{-a r_A - b r_B} e^{-i k.r}
        = (-1)^{m1+m2} d^{m1}_a d^{m2}_b rho~^{(00)}

AZIMUTHAL REDUCTION.  All centres sit on the z axis, so Delta depends on |k| only
and the only angular dependence is the phase e^{-i k mu P_z(t)}, mu = cos theta_k.
Hence rho~ is a function of (k, mu) alone and the k-grid is 2-D.

The Coulomb weight also collapses: with d^3k = 2 pi k^2 dk dmu,

    (pq|rs) = (1/pi) int_0^inf dk int_{-1}^{1} dmu  rho~_pq^*(k,mu) rho~_rs(k,mu)

-- the 4pi/k^2 kernel is EXACTLY cancelled by the k^2 Jacobian, so the quadrature
weight is flat and there is no k->0 singularity at all (the protocol's "small-k
audit" is therefore about where the MASS sits, not about a divergence).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss


# ---------------------------------------------------------------- orbital spec
@dataclass(frozen=True)
class SOrb:
    """Hydrogenic s-orbital on the z axis:  phi(r) = N sum_m c_m r^m e^{-a r}."""
    z: float          # centre (z coordinate)
    a: float          # decay rate  (hydrogenic: a = Z/n)
    kind: str         # '1s', '2s' or '3s'
    label: str = ""

    @property
    def norm(self) -> float:
        """N with int |N sum_m c_m r^m e^{-a r}|^2 d^3r = 1 (exact moments)."""
        c = 2.0 * self.a
        s = 0.0
        for m1, c1 in self.terms:
            for m2, c2 in self.terms:
                s += c1 * c2 * math.factorial(m1 + m2 + 2) / c ** (m1 + m2 + 3)
        return 1.0 / math.sqrt(4.0 * math.pi * s)

    @property
    def terms(self) -> Tuple[Tuple[int, float], ...]:
        """Hydrogenic ns radial polynomial (up to normalisation), a = Z/n."""
        if self.kind == "1s":
            return ((0, 1.0),)
        if self.kind == "2s":
            return ((0, 1.0), (1, -self.a))
        if self.kind == "3s":
            return ((0, 1.0), (1, -2.0 * self.a), (2, 2.0 * self.a ** 2 / 3.0))
        raise ValueError(self.kind)

    def radial(self, r) -> np.ndarray:
        r = np.asarray(r, float)
        out = np.zeros_like(r)
        for m, c in self.terms:
            out = out + c * r ** m
        return self.norm * out * np.exp(-self.a * r)


# ------------------------------------------- closed-form a/b derivatives of F
_a, _b, _t, _k, _D = sp.symbols("a b t k D", positive=True)
_Delta = sp.sqrt(_t * (1 - _t) * _k ** 2 + _t * _a ** 2 + (1 - _t) * _b ** 2)
_F = sp.exp(-_D * _Delta) / _Delta

_DERIV_CACHE: Dict[Tuple[int, int], Callable] = {}


def _dF(p: int, q: int) -> Callable:
    """Lambdified  d^p/da^p d^q/db^q [ e^{-D Delta}/Delta ]  (numpy, broadcasting)."""
    key = (p, q)
    if key not in _DERIV_CACHE:
        expr = sp.diff(_F, _a, p, _b, q)
        _DERIV_CACHE[key] = sp.lambdify((_a, _b, _t, _k, _D), sp.powsimp(expr, force=True), "numpy")
    return _DERIV_CACHE[key]


# ------------------------------------------------------------- t (Feynman) grid
def feynman_grid(y_max: float = 24.0, panel: float = 0.4, n_g: int = 10):
    """Composite Gauss-Legendre in the logistic variable  t = 1/(1+e^{-y}).

    The t-integrand concentrates near t=0 and t=1 (where Delta is smallest) with
    width ~ rate^2/k^2, so a plain [0,1] grid fails at large k.  In y the two
    endpoint structures become O(1)-wide bumps at y ~ -+ 2 ln(k/rate), which a
    uniform-panel composite rule resolves for every k in range.
    Returns (t, w) with  int_0^1 f(t) dt = sum w_i f(t_i)  (dt = t(1-t) dy).
    """
    x, wx = leggauss(n_g)
    n_pan = int(round(2 * y_max / panel))
    edges = np.linspace(-y_max, y_max, n_pan + 1)
    ys, ws = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
        ys.append(mid + half * x)
        ws.append(half * wx)
    y = np.concatenate(ys)
    w = np.concatenate(ws)
    t = 1.0 / (1.0 + np.exp(-y))
    return t, w * t * (1.0 - t)


# ------------------------------------------------------- the density evaluator
def rho_tilde(p: SOrb, q: SOrb, k, mu, tgrid=None, chunk: int = 48) -> np.ndarray:
    """rho~_pq(k, mu) on the outer product grid k x mu.  Shape (nk, nmu), complex.

    Exact analytic evaluation (Feynman reduction + closed-form a,b derivatives);
    the only numerics is the 1-D t quadrature, converged to ~1e-12.
    """
    if tgrid is None:
        tgrid = feynman_grid()
    t, wt = tgrid
    k = np.atleast_1d(np.asarray(k, float))
    mu = np.atleast_1d(np.asarray(mu, float))
    a, b = p.a, q.a
    D = abs(p.z - q.z)
    Pz = (1.0 - t) * p.z + t * q.z

    if D < 1e-14:                   # EXACT closed form; phase factors out
        rad = rho_tilde_onecentre(p, q, k)
        return rad[:, None] * np.exp(-1j * np.outer(k, mu) * p.z)
    out = np.zeros((k.size, mu.size), dtype=complex)
    onecentre = False
    for i0 in range(0, k.size, chunk):
        ks = np.maximum(k[i0:i0 + chunk], 1e-12)
        G = np.zeros((t.size, ks.size))
        for (m1, c1) in p.terms:
            for (m2, c2) in q.terms:
                f = _dF(m1 + 1, m2 + 1)
                val = f(a, b, t[:, None], ks[None, :], D)
                G = G + ((-1) ** (m1 + m2)) * c1 * c2 * np.asarray(val, float)
        G = G * wt[:, None]                                   # (nt, nk)
        if onecentre:
            rad = G.sum(axis=0)                               # (nk,) real
            out[i0:i0 + chunk, :] = (rad[:, None]
                                     * np.exp(-1j * np.outer(ks, mu) * p.z))
        else:
            ph = np.exp(-1j * (ks[None, :, None] * mu[None, None, :])
                        * Pz[:, None, None])
            out[i0:i0 + chunk, :] = np.einsum("tk,tkm->km", G, ph)
    return 2.0 * math.pi * p.norm * q.norm * out


# ------------------------------------------ CLOSED-FORM one-centre density FT
def rho_tilde_onecentre(p: SOrb, q: SOrb, k) -> np.ndarray:
    """EXACT rho~ for a ONE-CENTRE product (both orbitals on the same centre).

    rho_pq(r) = N_p N_q sum_m c_m r^m e^{-c r},  c = a_p + a_q, and

        rho~(k) = (4 pi / k) N_p N_q sum_m c_m Im[ (m+1)! / (c - i k)^{m+2} ]

    (elementary; no quadrature).  The full rho~(k,mu) is this times e^{-i k mu z0}.
    Sanity: 1s x 1s (rate a) gives the textbook 16 a^4/(k^2+4a^2)^2.
    """
    assert abs(p.z - q.z) < 1e-14
    k = np.atleast_1d(np.asarray(k, float))
    c = p.a + q.a
    prod: Dict[int, float] = {}
    for m1, c1 in p.terms:
        for m2, c2 in q.terms:
            prod[m1 + m2] = prod.get(m1 + m2, 0.0) + c1 * c2
    z = c - 1j * k
    out = np.zeros(k.size)
    for m, cm in prod.items():
        out = out + cm * math.factorial(m + 1) * np.imag(z ** (-(m + 2)))
    return (4.0 * math.pi / np.maximum(k, 1e-300)) * p.norm * q.norm * out


# --------------------------------------------------- independent direct-FT check
def rho_tilde_direct(p: SOrb, q: SOrb, kval: float, muval: float,
                     n_g: int = 40, n_pan: int = 30, n_eta: int = 260) -> complex:
    """Brute-force FT of phi_p phi_q by 2-D prolate-spheroidal quadrature.

    Foci at the two centres; the azimuthal integral is done analytically (J0).
    Same-centre pairs fall back to the 1-D spherical (j0) form.  Completely
    independent of the Feynman representation above.
    """
    from scipy.special import j0 as bessel_j0
    kz = kval * muval
    kperp = kval * math.sqrt(max(0.0, 1.0 - muval ** 2))
    if abs(p.z - q.z) < 1e-14:                                  # one centre
        rate = p.a + q.a
        x, wx = leggauss(n_eta)
        edges = np.linspace(0.0, 70.0 / rate, n_pan + 1)
        tot = 0.0
        for lo, hi in zip(edges[:-1], edges[1:]):
            m, h = 0.5 * (lo + hi), 0.5 * (hi - lo)
            r = m + h * x
            w = h * wx
            kr = kval * r
            jj = np.where(kr > 1e-12, np.sin(kr) / np.where(kr > 1e-12, kr, 1.0), 1.0)
            tot += 4 * math.pi * float(np.sum(w * p.radial(r) * q.radial(r) * jj * r ** 2))
        return tot * np.exp(-1j * kz * p.z)
    c = 0.5 * abs(p.z - q.z)
    hi_z, lo_z = (p, q) if p.z > q.z else (q, p)                # focus F1 = higher z
    M = 0.5 * (p.z + q.z)
    rate = hi_z.a + lo_z.a
    x, wx = leggauss(n_g)
    xi_max = 1.0 + 70.0 / (c * rate)
    edges = 1.0 + (xi_max - 1.0) * np.linspace(0, 1, n_pan + 1) ** 1.5
    xe, we = leggauss(n_eta)
    tot = 0.0 + 0j
    for lo, hi in zip(edges[:-1], edges[1:]):
        m, h = 0.5 * (lo + hi), 0.5 * (hi - lo)
        xi = (m + h * x)[:, None]
        w_xi = (h * wx)[:, None]
        eta = xe[None, :]
        w_eta = we[None, :]
        rA = c * (xi - eta)                                     # dist to F1 (higher z)
        rB = c * (xi + eta)                                     # dist to F2 (lower z)
        f = hi_z.radial(rA) * lo_z.radial(rB)
        rho_perp = c * np.sqrt(np.maximum((xi ** 2 - 1) * (1 - eta ** 2), 0.0))
        integ = ((xi ** 2 - eta ** 2) * f * np.exp(-1j * kz * c * xi * eta)
                 * bessel_j0(kperp * rho_perp))
        tot += np.sum(w_xi * w_eta * integ)
    return 2 * math.pi * c ** 3 * tot * np.exp(-1j * kz * M)


# ------------------------------------------------------------ systems + Loewdin
def system(name: str):
    """(orbitals, nuclei) for the pre-registered probe systems."""
    if name == "H2":
        R = 1.4
        return ([SOrb(+R / 2, 1.0, "1s", "H_A 1s"), SOrb(-R / 2, 1.0, "1s", "H_B 1s")],
                [(+R / 2, 1.0), (-R / 2, 1.0)])
    if name == "LiH":
        R = 3.015
        return ([SOrb(0.0, 3.0, "1s", "Li 1s"), SOrb(0.0, 1.5, "2s", "Li 2s"),
                 SOrb(R, 1.0, "1s", "H 1s")],
                [(0.0, 3.0), (R, 1.0)])
    if name == "H2_4o":
        R = 1.4
        return ([SOrb(+R / 2, 1.0, "1s", "H_A 1s"), SOrb(+R / 2, 0.5, "2s", "H_A 2s"),
                 SOrb(-R / 2, 1.0, "1s", "H_B 1s"), SOrb(-R / 2, 0.5, "2s", "H_B 2s")],
                [(+R / 2, 1.0), (-R / 2, 1.0)])
    if name == "LiH_5o":
        R = 3.015
        return ([SOrb(0.0, 3.0, "1s", "Li 1s"), SOrb(0.0, 1.5, "2s", "Li 2s"),
                 SOrb(0.0, 1.0, "3s", "Li 3s"), SOrb(R, 1.0, "1s", "H 1s"),
                 SOrb(R, 0.5, "2s", "H 2s")],
                [(0.0, 3.0), (R, 1.0)])
    raise ValueError(name)


def overlap_matrix(orbs: Sequence[SOrb], tgrid=None) -> np.ndarray:
    """S_pq = rho~_pq(0) -- the k->0 value of the transition density."""
    n = len(orbs)
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            v = rho_tilde(orbs[i], orbs[j], np.array([1e-9]), np.array([0.0]), tgrid)[0, 0]
            S[i, j] = S[j, i] = v.real
    return S


def lowdin_X(S: np.ndarray) -> np.ndarray:
    w, U = np.linalg.eigh(S)
    keep = w > 1e-10
    return U[:, keep] @ np.diag(1.0 / np.sqrt(w[keep])) @ U[:, keep].T


if __name__ == "__main__":
    print("Probe B -- analytic momentum transition density: VALIDATION\n")
    tg = feynman_grid()
    print(f"  Feynman t-grid: {tg[0].size} nodes (logistic composite GL)\n")

    for sysname in ("H2", "LiH", "H2_4o"):
        orbs, _ = system(sysname)
        print(f"--- {sysname}: {[o.label for o in orbs]}")
        S = overlap_matrix(orbs, tg)
        print("  S = rho~(0):")
        for row in S:
            print("     " + "  ".join(f"{v:+.10f}" for v in row))
        if sysname == "LiH":
            print(f"  [check] <Li1s|Li2s> hydrogenic orthogonality: {S[0,1]:+.2e} (exact 0)")
        worst = 0.0
        for i in range(len(orbs)):
            for j in range(i, len(orbs)):
                for (kv, mv) in ((0.7, 0.3), (2.5, -0.8), (6.0, 1.0), (13.0, 0.45)):
                    fast = rho_tilde(orbs[i], orbs[j], np.array([kv]),
                                     np.array([mv]), tg)[0, 0]
                    ref = rho_tilde_direct(orbs[i], orbs[j], kv, mv)
                    err = abs(fast - ref)
                    worst = max(worst, err)
                    if err > 1e-8:
                        print(f"     !! ({i},{j}) k={kv} mu={mv}: "
                              f"feynman {fast:.12g}  direct {ref:.12g}  |d|={err:.2e}")
        print(f"  max |rho~_feynman - rho~_direct| over pairs x 4 k-points: {worst:.3e}\n")
