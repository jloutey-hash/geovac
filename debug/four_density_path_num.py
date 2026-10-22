"""FP1-FP3 numerics (record docs/predictions/2026-10-22_four_density_path.md, commit 66af9822).

Position space (derivation; conventions of four_density_star_B_numerics.py):
  rho_a(r) = a^3 K0(a r)/(2 pi^2),  phi_a(r) = (2/pi) chi(a r)/r  (unit leaf).
  Path a - b - c - d, scales (1, t, 1, 1):
  E = int int [rho_b phi_a](r) [rho_c phi_d](r') / |r - r'| d^3r d^3r'
    = (4 pi)^2 int int r^2 r'^2 [t^3 K0(t r)/(2 pi^2)] (2/pi) chi(r)/r [K0(r')/(2 pi^2)] (2/pi) chi(r')/r' / max(r, r')
    = (16/pi^4) t^3 int int min(r, r') K0(t r) chi(r) K0(r') chi(r') dr dr'
    = N t^3 int_0^oo K0(t r) chi(r) H(r) dr,   N = 16/pi^4,
  H(r) = int_0^r r' F(r') dr' + r int_r^oo F(r') dr',  F = K0 chi.
H is computed by its two one-sided integrals on the double-exponential grid x = exp(s - e^-s):
Gauss-Legendre (m points) on each grid panel in the variable s, cumulative sums.
theta_t^k P = N t^3 sum_m C(k,m) 3^(k-m) I_m,  I_m = int g(x) theta^m[K0(t x)] dx,
theta^m K0(y) = a_m(y) K0 + b_m(y) theta K0 (integer polynomials, theta^2 K0 = y^2 K0).

Momentum control (independent of chi, H, K0): with f(k) = (1+k^2)^(-3/2)/k^2,
  E = (4 pi)^3 (2 pi)^-9 I = I/(8 pi^6),  I = int d^3q q^-2 M_t(q) M_1(q) = 4 pi int_0^oo M_t M_1 dq,
  M_t = f * g_t (3-D convolution), g_t(k) = (1 + k^2/t^2)^(-3/2)  (four_density_star_B_numerics._conv).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))


def chi(x):
    x = mp.mpf(x)
    if x == 0:
        return mp.mpf(0)
    with mp.extradps(20):
        return +(mp.pi * x / 2 * (mp.besselk(1, x) * mp.struvel(0, x) + mp.besselk(0, x) * mp.struvel(1, x)))


def theta_k0_polys(n):
    """integer coefficient lists a_m, b_m (index = power of y), m = 0..n."""
    al, be = [[1]], [[0]]
    for _ in range(n):
        a, b = al[-1], be[-1]
        na = [0] * (max(len(a), len(b) + 2) + 1)
        nb = [0] * (max(len(a), len(b)) + 1)
        for i, c in enumerate(a):        # y a'
            na[i] += i * c
            nb[i] += c                   # a -> b_{m+1}
        for i, c in enumerate(b):
            na[i + 2] += c               # y^2 b
            nb[i] += i * c               # y b'
        while len(na) > 1 and na[-1] == 0:
            na.pop()
        while len(nb) > 1 and nb[-1] == 0:
            nb.pop()
        al.append(na)
        be.append(nb)
    return al, be


def horner(cs, y):
    acc = mp.mpf(0)
    for c in reversed(cs):
        acc = acc * y + c
    return acc


class Grid:
    def __init__(self, dps, hinv, smin=-5.5, smax=None, m=10, verbose=True):
        mp.mp.dps = dps
        self.dps, self.hinv = dps, hinv
        if smax is None:
            smax = float(mp.log(2.4 * dps / 0.45 + 60))   # K0(t x) at t >= 0.45 negligible beyond
        h = mp.mpf(1) / hinv
        ks = list(range(int(smin * hinv), int(smax * hinv) + 2))
        sv = [k * h for k in ks]
        xs = [mp.exp(s_ - mp.exp(-s_)) for s_ in sv]
        ws = [h * x * (1 + mp.exp(-s_)) for s_, x in zip(sv, xs)]
        gl_x, gl_w = zip(*[(mp.mpf(a), mp.mpf(b)) for a, b in mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(
            self._gl_degree(m), mp.mp.prec)])
        t0 = time.time()
        cum1, cum0 = [mp.mpf(0)], [mp.mpf(0)]   # int_0^{x_k} r F, int_0^{x_k} F (from x_0, below negligible)
        Fx = []
        for i in range(len(sv) - 1):
            a, b = sv[i], sv[i + 1]
            mid, half = (a + b) / 2, (b - a) / 2
            s1 = s0 = mp.mpf(0)
            for gx, gw in zip(gl_x, gl_w):
                s_ = mid + half * gx
                x = mp.exp(s_ - mp.exp(-s_))
                jac = x * (1 + mp.exp(-s_))
                F = mp.besselk(0, x) * chi(x)
                s0 += gw * half * jac * F
                s1 += gw * half * jac * F * x
            cum1.append(cum1[-1] + s1)
            cum0.append(cum0[-1] + s0)
        self.Hinf = cum1[-1]
        Ftot = cum0[-1]
        self.nodes = []
        for x, w, c1, c0 in zip(xs, ws, cum1, cum0):
            H = c1 + x * (Ftot - c0)
            cx = chi(x)
            self.nodes.append((x, w, cx * H, cx, H))
        self.ginf = mp.pi / 2 * self.Hinf
        if verbose:
            print(f"grid dps={dps} hinv={hinv} s=[{smin},{float(smax):.2f}] nodes={len(self.nodes)} "
                  f"Hinf={mp.nstr(self.Hinf, 25)} ({time.time()-t0:.0f}s)", flush=True)

    @staticmethod
    def _gl_degree(m):
        # mpmath GaussLegendre degree d gives 3*2^(d-1) nodes
        d = 1
        while 3 * 2 ** (d - 1) < m:
            d += 1
        return d

    def thetas(self, tv, J, which="g"):
        """[theta_t^k P](tv), k = 0..J, P = N t^3 int K0(t x) G(x) dx with G = g (default)."""
        tv = mp.mpf(tv)
        al, be = theta_k0_polys(J)
        I = [mp.mpf(0)] * (J + 1)
        for x, w, gx, cx, H in self.nodes:
            G = {"g": gx, "chi": cx, "H": H}[which] if isinstance(which, str) else which(x, gx, cx, H)
            if G == 0:
                continue
            y = tv * x
            k0 = mp.besselk(0, y)
            if abs(k0) * abs(G) * w < mp.mpf(10) ** (-self.dps - 30):
                continue
            thk0 = -y * mp.besselk(1, y)
            wg = w * G
            for j in range(J + 1):
                I[j] += wg * (horner(al[j], y) * k0 + horner(be[j], y) * thk0)
        N = 16 / mp.pi**4
        out = []
        for k in range(J + 1):
            out.append(N * tv**3 * sum(mp.binomial(k, m_) * 3 ** (k - m_) * I[m_] for m_ in range(k + 1)))
        return out


def P_mom(tv):
    import four_density_star_B_numerics as SB
    tv = mp.mpf(tv)
    kf = lambda k: (1 + k * k) ** mp.mpf(-1.5) / k

    def Bf(q):
        if q == 0:
            return -mp.inf
        w = mp.sqrt(1 + q * q)
        return 1 / w + mp.log(q / (w + 1))
    Bt = lambda q: -tv**2 / mp.sqrt(1 + q * q / tv**2)
    B1 = lambda q: -1 / mp.sqrt(1 + q * q)
    integrand = lambda P: SB._conv(kf, Bt, P) * SB._conv(kf, B1, P) if P != 0 else mp.mpf(0)
    return 4 * mp.pi * mp.quad(integrand, [0, mp.mpf(1) / 4, 1, 3, 10, mp.inf]) / (8 * mp.pi**6)


def evaluate(L, thv, tv):
    """L function-left in delta = theta/2, w = t^2 (four_density_star_ops); returns (value, largest, rel)."""
    import sympy as sp
    import four_density_star_ops as O
    wv = sp.Rational(str(tv)) ** 2
    terms = []
    for k, c in enumerate(L):
        if c == 0:
            continue
        num = sp.Poly(c.numer.as_expr(), O.w_sym).eval(wv) if c.numer.as_expr().free_symbols else sp.Rational(c.numer.as_expr())
        den = sp.Poly(c.denom.as_expr(), O.w_sym).eval(wv) if c.denom.as_expr().free_symbols else sp.Rational(c.denom.as_expr())
        cv = mp.mpf(sp.Rational(num).p) / sp.Rational(num).q / (mp.mpf(sp.Rational(den).p) / sp.Rational(den).q)
        terms.append(cv * thv[k] / mp.mpf(2) ** k)
    tot = sum(terms)
    big = max(abs(x) for x in terms)
    return tot, big, abs(tot) / big
