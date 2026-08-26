"""Core high-precision evaluator for the Paper 59 collinear observable T2 in the (KW) frame:

    T2 = (8/pi) int_0^inf dk int_0^1 dw cos(k w) R(k,w)^2,
    R(k,w) = int_0^1 cos(k w (s-1/2)) P(s,k) ds = 2 int_0^{1/2} cos(k w (1/2 - s)) P(s,k) ds,
    P(x,k) = c e^{-D}(D^-3+3D^-4+3D^-5), c=x(1-x), D=sqrt(c k^2+1).

Pieces
  * R via s = (1/2) y^p graded GL (p resolves the 1/k^2 endpoint layer analytically).
  * F(k) = int_0^1 cos(kw) R^2 dw by GL (integrand entire in w, bandwidth 2k -> N_w ~ 0.6k+40).
  * int_0^K F dk by panelled GL.
  * int_K^inf F dk ANALYTICALLY: Watson expansion of R in x=1/k (Sec.2 of the findings memo).
"""
from __future__ import annotations
import sys
sys.path.insert(0, 'debug')
import mpmath as mp
from _fastgl import fast_gl


# ---------------------------------------------------------------- basic pieces
def P(x, k):
    c = x * (1 - x)
    D = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-D) * (1 / D ** 3 + 3 / D ** 4 + 3 / D ** 5)


_SNODES = {}
def s_nodes(Ns, p):
    """GL nodes/weights for  2 int_0^{1/2} f(s) ds  under s = (1/2) y^p, y in [0,1]."""
    key = (Ns, p, mp.mp.dps)
    if key in _SNODES:
        return _SNODES[key]
    xs, ws = fast_gl(Ns)
    S, W = [], []
    for xg, wg in zip(xs, ws):
        y = (xg + 1) / 2
        s = y ** p / 2
        jac = p * y ** (p - 1) / 2          # ds/dy
        S.append(s)
        W.append(2 * (wg / 2) * jac)        # the outer factor 2 of R folded in
    _SNODES[key] = (S, W)
    return S, W


def R_setup(k, Ns, p):
    """Precompute (1/2 - s_j) and W_j*P(s_j,k) for fixed k."""
    S, W = s_nodes(Ns, p)
    a = []; c = []
    for s, w in zip(S, W):
        a.append((mp.mpf(1) / 2 - s) * k)
        c.append(w * P(s, k))
    return a, c


def R_of_w(a, c, w):
    tot = mp.mpf(0)
    for aj, cj in zip(a, c):
        tot += cj * mp.cos(aj * w)
    return tot


def F_of_k(k, Ns, Nw, p):
    """F(k) = int_0^1 cos(kw) R(k,w)^2 dw."""
    a, c = R_setup(k, Ns, p)
    xs, ws = fast_gl(Nw)
    tot = mp.mpf(0)
    for xg, wg in zip(xs, ws):
        w = (xg + 1) / 2
        Rv = R_of_w(a, c, w)
        tot += (wg / 2) * mp.cos(k * w) * Rv * Rv
    return tot


def nw_for(k, base=40, slope=mp.mpf('0.62')):
    return int(slope * float(k)) + base


def int_F(a, b, npan, nnode, Ns, p=6, nw=None, prog=None):
    """Panelled Gauss-Legendre for int_a^b F(k) dk."""
    a = mp.mpf(a); b = mp.mpf(b)
    xs, ws = fast_gl(nnode)
    h = (b - a) / npan
    tot = mp.mpf(0)
    for ip in range(npan):
        lo = a + ip * h
        for xg, wg in zip(xs, ws):
            k = lo + h * (xg + 1) / 2
            NW = nw(k) if nw else nw_for(k)
            NS = Ns(k) if callable(Ns) else Ns
            tot += (h * wg / 2) * F_of_k(k, NS, NW, p)
        if prog: prog(ip, npan, tot)
    return tot
