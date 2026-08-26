"""STRUCTURALLY INDEPENDENT route: the ORIGINAL (s,t)-outer / k-inner frame,
    T2 = (8/pi) int_0^1 int_0^1 J(s,t) ds dt,   J = int_0^inf j0(kb) P(s,k)P(t,k) dk,
with an ANALYTIC large-k fibre tail (no Vandermonde fit) and the sigma^2-Duffy corner.

Fibre tail: P(x,k) = k^-3 g(1/k) e^{-sqrt(c) k}  with
    g(u) = c exp(-u/(sqrt c + sqrt(c+u^2))) [ (c+u^2)^{-3/2} + 3u(c+u^2)^{-2} + 3u^2(c+u^2)^{-5/2} ]
analytic at u=0 (the exponent's removable singularity is resolved in closed form).  Hence
    P(s,k)P(t,k) = e^{-Ak} sum_m d_m k^{-6-m},  A = sqrt(c_s)+sqrt(c_t), d = taylor(g_s g_t),
    int_Kf^inf j0(kb) PP dk = (1/b) Im[ sum_m d_m z^{6+m} Gamma(-6-m, z Kf) ],  z = A - i b.
Shares NOTHING with the (KW) frame except the definition of P.
"""
from __future__ import annotations
import sys
sys.path.insert(0, 'debug')
import mpmath as mp
from _fastgl import fast_gl
from beta2_t2_kw_core import P
from _ser import smul, sinv, sexp

_GT = {}
def gtaylor(c, M):
    """Taylor coefficients of  g(u) = c exp(-u/(sqrt c + sqrt(c+u^2)))
    [ (c+u^2)^{-3/2} + 3u (c+u^2)^{-2} + 3u^2 (c+u^2)^{-5/2} ]  at u=0,
    by explicit power-series arithmetic (mp.taylor's numerical differentiation is far too slow)."""
    key = (mp.nstr(c, 32), M, mp.mp.dps)
    if key in _GT: return _GT[key]
    sc = mp.sqrt(c)
    def Apow(p):                      # (1 + u^2/c)^p  as a series in u
        out = [mp.mpf(0)]*(M+1); out[0] = mp.mpf(1)
        for m in range(1, M//2 + 1):
            out[2*m] = mp.binomial(p, m)/c**m
        return out
    sqrtA = Apow(mp.mpf(1)/2)
    den = [sc*v for v in sqrtA]; den[0] += sc          # sqrt c + sqrt(c+u^2)
    E = smul([mp.mpf(0), mp.mpf(-1)] + [mp.mpf(0)]*(M-1), sinv(den, M), M)
    eE = sexp(E, M)
    T = [mp.mpf(0)]*(M+1)
    a32 = Apow(mp.mpf(-3)/2); a2 = Apow(-2); a52 = Apow(mp.mpf(-5)/2)
    f32 = c**(mp.mpf(-3)/2); f2 = c**mp.mpf(-2); f52 = c**(mp.mpf(-5)/2)
    for i in range(M+1):
        T[i] += f32*a32[i]
        if i >= 1: T[i] += 3*f2*a2[i-1]
        if i >= 2: T[i] += 3*f52*a52[i-2]
    g = [c*v for v in smul(eE, T, M)]
    _GT[key] = g
    return g


def fibre_tail(cs, ct, b, Kf, M):
    """(1/b) Im[ sum_m d_m T_m ],  T_m = int_Kf^inf k^{-(7+m)} e^{-zk} dk,  z = A - i b.
    T_M seeded from Gamma(-6-M, zKf); the recurrence is run DOWNWARD in m
    (T_m = Kf^{-(7+m)}e^{-zKf}/z - ((7+m)/z) T_{m+1}), which is the CANCELLATION-FREE
    direction (upward loses ~log10(|z|Kf/(6+m)) digits per step)."""
    ds = gtaylor(cs, M); dt = gtaylor(ct, M)
    d = [sum(ds[i] * dt[m - i] for i in range(m + 1)) for m in range(M + 1)]
    A = mp.sqrt(cs) + mp.sqrt(ct)
    z = mp.mpc(A, -b)
    x = z * Kf
    ex = mp.e ** (-x)
    T = [None] * (M + 1)
    T[M] = z ** (6 + M) * mp.gammainc(-6 - M, x)
    for m in range(M - 1, -1, -1):
        T[m] = Kf ** (-(7 + m)) * ex / z - ((7 + m) / z) * T[m + 1]
    tot = mp.mpc(0)
    for m in range(M + 1):
        tot += d[m] * T[m]
    return (tot / b).imag


_NQ = (128, 160, 200, 256, 320, 400, 512, 640, 800, 1024, 1300, 1700, 2200)


def J_acc(s, t, M=26, kfac=mp.mpf(15), acut=mp.mpf(140)):
    """int_0^inf j0(kb) P(s,k)P(t,k) dk.

    Two regimes, both exponentially safe:
      * tail regime  Kf = kfac/sqrt(cmin) <= acut/A : bounded GL on [0,Kf] + ANALYTIC tail
        (series radius sqrt(cmin) => truncation ~ kfac^-M).
      * cut regime   acut/A < kfac/sqrt(cmin) : the integrand is already < e^{-acut} there,
        so integrate to Kf = acut/A and drop the tail (neglected < e^{-125} ~ 5e-55).
    Node count resolves BOTH the exponential decay (A*Kf <= acut) and the j0 oscillation (Kf*b).
    """
    cs = s * (1 - s); ct = t * (1 - t)
    if cs == 0 or ct == 0: return mp.mpf(0)
    b = s + t
    cmin = min(cs, ct)
    A = mp.sqrt(cs) + mp.sqrt(ct)
    Kt = kfac / mp.sqrt(cmin)
    Kc = acut / A
    use_tail = (Kt <= Kc)
    Kf = Kt if use_tail else Kc
    span = max(float(A * Kf), float(Kf * b))
    npan = max(1, int(span / 3.5) + 1)
    nn = 44
    xs, ws = fast_gl(nn)
    h = Kf / npan
    tot = mp.mpf(0)
    for ip in range(npan):
        lo = ip * h
        for xg, wg in zip(xs, ws):
            k = lo + h * (xg + 1) / 2
            kb = k * b
            j0 = mp.sin(kb) / kb if kb > mp.mpf('1e-40') else mp.mpf(1)
            tot += (h * wg / 2) * j0 * P(s, k) * P(t, k)
    if use_tail:
        tot += fibre_tail(cs, ct, b, Kf, M)
    return tot


def T2_st(delta, Nc, Nt, M=8, Jf=None):
    """sigma^2-Duffy corner on {s+t<delta} + two kink-free blocks (same decomposition as
    routeC_T2_highprec, different fibre)."""
    J = Jf or J_acc
    delta = mp.mpf(delta)
    H = mp.pi / 2
    xs, ws = fast_gl(Nc); xa, wa = fast_gl(Nc)
    smax = mp.sqrt(delta); tot = mp.mpf(0)
    for xi, wi in zip(xs, ws):                       # corner
        sig = smax * (xi + 1) / 2; wsig = smax * wi / 2
        for xj, wj in zip(xa, wa):
            psi = H * (xj + 1) / 2; wpsi = H * wj / 2
            al = mp.sin(psi) ** 2; jal = mp.sin(2 * psi)
            s = sig * sig * al; t = sig * sig * (1 - al)
            tot += wsig * wpsi * jal * 2 * sig ** 3 * J(s, t, M)
    xs, ws = fast_gl(Nt); xa, wa = fast_gl(Nt)
    for xi, wi in zip(xs, ws):                       # block 1: s in [0,delta]
        phis = H * (xi + 1) / 2; s = delta * mp.sin(phis) ** 2
        js = delta * mp.sin(2 * phis); wphis = H * wi / 2
        lo = delta - s
        for xj, wj in zip(xa, wa):
            phit = H * (xj + 1) / 2
            t = lo + (1 - lo) * mp.sin(phit) ** 2
            jt = (1 - lo) * mp.sin(2 * phit); wphit = H * wj / 2
            tot += wphis * js * wphit * jt * J(s, t, M)
    for xi, wi in zip(xs, ws):                       # block 2: s in [delta,1]
        phis = H * (xi + 1) / 2; s = delta + (1 - delta) * mp.sin(phis) ** 2
        js = (1 - delta) * mp.sin(2 * phis); wphis = H * wi / 2
        for xj, wj in zip(xa, wa):
            phit = H * (xj + 1) / 2; t = mp.sin(phit) ** 2
            jt = mp.sin(2 * phit); wphit = H * wj / 2
            tot += wphis * js * wphit * jt * J(s, t, M)
    return (8 / mp.pi) * tot
