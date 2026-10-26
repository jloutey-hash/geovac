"""GF1-GF3 (docs/predictions/2026-10-25_gaussian_fit_geminal.md, committed 9cac600d before computing): the engine.

Y_GG(A, B; t) = <rho_i rho_j rho_k exp(-A r_ik^2) exp(-B r_jk^2) / r_ij>, density scales a (i), b (j), t (apex k),
fingerprint fibre rho~_c(k) = c^3 (k^2 + c^2)^(-3/2), conventions of record QG as in debug/pf_common.py:
  Y = (1/(2 pi^2)) int_0^oo dw w^2 rho~_t(w) Phi~(w),
  Phi~(w) = (2 pi)^-3 (2 pi / w) int int_{|p-q| <= w <= p+q} p q rho~_a(p) rho~_b(q) T(p, q, w) dp dq,
  T = (2 pi)^-3 int d^3q' g~_A(q') g~_B(q' - L) 4 pi / |k - q'|^2,   g~_A(q) = (pi/A)^(3/2) exp(-q^2/(4A)).
Gaussian pair:  g~_A(q) g~_B(q - L) = (pi^2/(A B))^(3/2) exp(-(q - cL)^2/(4 mu)) exp(-L^2/(4(A+B))),
  mu = A B/(A+B), c = A/(A+B).
Coulomb against one Gaussian:  G(R; mu) = int d^3s exp(-s^2/(4 mu)) / |s - R|^2 = 8 mu pi^(3/2) Daw(R/(2 sqrt mu)) / R
  (record GY's loop at mu = g/2; checked here against direct quadrature by gf_checks.py),
  R^2 = |k - cL|^2 = (1-c) p^2 + c q^2 - c(1-c) w^2.
Inner variables sigma = p + q in (w, oo), delta = p - q in (-w, w): p q dp dq = (sigma^2 - delta^2)/4 * dsigma ddelta/2.
float64, Gauss-Legendre with algebraic maps; Phi~ on the w nodes is shared by every t."""
from __future__ import annotations
import numpy as np
from scipy.special import dawsn

PI = np.pi


def rho(c, k):
    return c ** 3 * (k * k + c * c) ** -1.5


def G_loop(R, mu):
    x = R / (2 * np.sqrt(mu))
    small = x < 1e-6
    out = np.empty_like(R)
    out[~small] = 8 * mu * PI ** 1.5 * dawsn(x[~small]) / R[~small]
    out[small] = 8 * mu * PI ** 1.5 / (2 * np.sqrt(mu))  # Daw(x)/R -> 1/(2 sqrt mu)
    return out


def gl(n, lo, hi):
    x, w = np.polynomial.legendre.leggauss(n)
    return (hi - lo) / 2 * x + (hi + lo) / 2, (hi - lo) / 2 * w


def half_line(n, lo, scale):
    """nodes on (lo, oo): s = lo + scale * u^2/(1-u)^2 ... use u in (0,1), s = lo + scale * u/(1-u), cubed map
    for the algebraic tail."""
    u, wu = gl(n, 0.0, 1.0)
    s = lo + scale * (u / (1 - u)) ** 2
    ds = scale * 2 * u / (1 - u) ** 3
    return s, wu * ds


def geo(lo, hi, n, h0=0.5):
    """Gauss-Legendre on [lo, hi] split at lo + h0 * 2^k."""
    br = [lo]
    h = h0
    while br[-1] + h < hi:
        br.append(lo + h)
        h *= 2
    br.append(hi)
    xs, ws = zip(*[gl(n, br[j], br[j + 1]) for j in range(len(br) - 1)])
    return np.concatenate(xs), np.concatenate(ws)


def two_sided(lo, hi, n, h0=0.5):
    """[lo, hi] refined geometrically towards both ends."""
    m = (lo + hi) / 2
    x1, w1 = geo(lo, m, n, h0)
    x2, w2 = geo(-hi, -m, n, h0)
    return np.concatenate([x1, -x2]), np.concatenate([w1, w2])


def w_nodes(nw, A, B):
    """w in (0, oo): geometric pieces to 4 + 12 sqrt(A+B), then an algebraic tail."""
    cut = 4.0 + 12 * np.sqrt(A + B)
    x1, w1 = geo(0.0, cut, nw)
    x2, w2 = half_line(2 * nw, cut, cut)
    return np.concatenate([x1, x2]), np.concatenate([w1, w2])


def phi_w(w, A, B, a=1.0, b=2.0, ns=24, nd=24):
    """Phi~ at each w (array)."""
    mu = A * B / (A + B)
    c = A / (A + B)
    pref_T = (2 * PI) ** -3 * 4 * PI * (PI * PI / (A * B)) ** 1.5
    out = np.empty_like(w)
    for i, wi in enumerate(w):
        cut = wi + 8.0 + 12 * np.sqrt(mu)
        s1, ws1 = geo(wi, cut, ns)
        s2, ws2 = half_line(2 * ns, cut, cut)
        sg = np.concatenate([s1, s2]); wsg = np.concatenate([ws1, ws2])
        dl, wdl = two_sided(-wi, wi, nd, min(0.5, wi / 2))
        S, D = np.meshgrid(sg, dl, indexing="ij")
        WS = np.outer(wsg, wdl)
        p = (S + D) / 2
        q = (S - D) / 2
        R2 = (1 - c) * p * p + c * q * q - c * (1 - c) * wi * wi
        R = np.sqrt(np.maximum(R2, 0.0))
        T = pref_T * np.exp(-wi * wi / (4 * (A + B))) * G_loop(R, mu)
        f = (S * S - D * D) / 8 * rho(a, p) * rho(b, q) * T
        out[i] = (2 * PI) ** -3 * 2 * PI / wi * np.sum(WS * f)
    return out


def Y_GG(A, B, ts, a=1.0, b=2.0, nw=24, ns=24, nd=24):
    w, ww = w_nodes(nw, A, B)
    ph = phi_w(w, A, B, a, b, ns, nd)
    return np.array([1 / (2 * PI * PI) * np.sum(ww * w * w * rho(t, w) * ph) for t in ts])
