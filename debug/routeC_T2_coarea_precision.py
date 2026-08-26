"""G1: high-precision co-area twist Phi(rho) via the analytic-tail fibre, and the 1D integral
T2 = (16/pi) int_0^1 Phi(rho) drho pushed past the ~13-16 digit walls of the (s,t)-quadrature and
the fixed-Nk fibre.  Phi(rho) = int_0^umax [sum_4 branch fibre] * u/(sqrt(1-4u) sqrt(1-4 rho u)) du,
umax=min(1/4,1/(4rho)).  The 4 branches share c_s=u, c_t=rho*u (P depends only on c), so ONE
bounded int_0^K pass (4 j0's summed) + 4 analytic tails per u-node.  Fibre accuracy ~37-40 digits
(validated in _jtail_check).  Reuses routeC_T2_highprec.P and _fastgl.fast_gl.
"""
from __future__ import annotations
import sys
sys.path.insert(0, 'debug')
import mpmath as mp
import routeC_T2_highprec as H
from _fastgl import fast_gl


def _dcoeffs(sm, tm, A, M):
    """1/k-series coeffs of g(k)=P(sm,k)P(tm,k) e^{Ak} k^6 (depends on c_s,c_t only)."""
    cmin = min(sm * (1 - sm), tm * (1 - tm))
    k0 = max(mp.mpf(60), 35 / mp.sqrt(cmin))
    ks = [mp.mpf(k0) * (i + 1) for i in range(M + 1)]
    gs = [H.P(sm, k) * H.P(tm, k) * mp.e ** (A * k) * k ** 6 for k in ks]
    V = mp.matrix(M + 1, M + 1)
    for i, k in enumerate(ks):
        for m in range(M + 1):
            V[i, m] = 1 / k ** m
    d = mp.lu_solve(V, mp.matrix(gs))
    return [d[m] for m in range(M + 1)]


def _tail_b(d, A, b, K, M):
    """int_K^inf j0(kb) [PP] dk via int_K^inf sin(kb) e^{-Ak} k^{-(6+m+1)} dk = Im[z^{...} Gamma],
    z=A-ib, using the shared d coeffs (g=PP e^{Ak} k^6 => PP = sum d_m e^{-Ak} k^{-(6+m)},
    and j0=sin(kb)/(kb) adds one more 1/k)."""
    z = A - 1j * b
    if abs(z) * K < mp.mpf('0.5'):
        return mp.mpf(0)
    x = z * K; ex = mp.e ** (-x)
    Gs = []; G = mp.gammainc(-6, x)
    for m in range(M + 1):
        Gs.append(G); a = -6 - m
        G = (G - x ** (a - 1) * ex) / (a - 1)
    tot = mp.mpc(0)
    for m in range(M + 1):
        tot += d[m] * z ** (6 + m) * Gs[m]
    val = (tot / b).imag
    return val if abs(val) < 1 else mp.mpf(0)


_NQ_LADDER = (300, 400, 500, 650, 850, 1100, 1500, 2000, 2800, 4000)


def branch_fibre(u, rho, M=6, guard=mp.mpf('1e-7')):
    """sum over the 4 co-area branches at c_s=u, c_t=rho*u, with ADAPTIVE K and Nq (drawn from a
    cached ladder so fast_gl is not recomputed).  K must exceed ~1/sqrt(cmin) for tail validity;
    Nq must resolve the K*bmax oscillation.  Guard where the fibre*weight is negligible."""
    sm = (1 - mp.sqrt(1 - 4 * u)) / 2
    tm = (1 - mp.sqrt(1 - 4 * rho * u)) / 2
    cs = u; ct = rho * u
    cmin = min(cs, ct)
    if cmin < guard:
        return mp.mpf(0)
    bs = [sm + tm, sm + (1 - tm), (1 - sm) + tm, (1 - sm) + (1 - tm)]
    A = mp.sqrt(cs) + mp.sqrt(ct)
    bmax = max(bs)
    K = max(mp.mpf(95), mp.mpf('3.2') / mp.sqrt(cmin))
    want = max(300, 2.2 * float(K * bmax) + 200)
    Nq = next((g for g in _NQ_LADDER if g >= want), _NQ_LADDER[-1])
    xs, ws = fast_gl(Nq)                       # cached by (Nq, dps)
    bnd = mp.mpf(0)
    for xg, wg in zip(xs, ws):
        k = K * (xg + 1) / 2; wk = K * wg / 2
        sj = mp.mpf(0)
        for b in bs:
            kb = k * b
            sj += mp.sin(kb) / kb if kb > mp.mpf('1e-40') else mp.mpf(1)
        bnd += wk * sj * H.P(sm, k) * H.P(tm, k)
    d = _dcoeffs(sm, tm, A, M)
    tl = sum(_tail_b(d, A, b, K, M) for b in bs)
    return bnd + tl


def Phi_acc(rho, Nu, M=6):
    """co-area twist via sin^2-map u-integration + adaptive analytic-tail branch fibre."""
    rho = mp.mpf(rho)
    umax = min(mp.mpf(1) / 4, 1 / (4 * rho))
    xs, ws = fast_gl(Nu); Hh = mp.pi / 2; tot = mp.mpf(0)
    for xg, wg in zip(xs, ws):
        phi = Hh * (xg + 1) / 2
        u = umax * mp.sin(phi) ** 2
        du = umax * mp.sin(2 * phi)
        wj = Hh * wg / 2
        bf = branch_fibre(u, rho, M=M)
        integ = bf * u / (mp.sqrt(1 - 4 * u) * mp.sqrt(1 - 4 * rho * u))
        tot += wj * du * integ
    return tot


if __name__ == '__main__':
    import time
    mp.mp.dps = 44
    for rho in ['0.5', '0.9']:
        prev = None
        print(f"rho={rho}:", flush=True)
        for Nu in (40, 80, 160):
            t0 = time.time(); v = Phi_acc(rho, Nu); dt = time.time() - t0
            d = '' if prev is None else f"  |dNu|={mp.nstr(abs(v - prev), 3)}"
            print(f"   Nu={Nu:4d}: {mp.nstr(v, 34)}{d}  ({dt:.1f}s)", flush=True)
            prev = v
