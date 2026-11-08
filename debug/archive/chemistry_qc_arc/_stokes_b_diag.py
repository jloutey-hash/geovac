r"""Diagnose the (c1,b)=(0.2,1.0) blowup: precision, or a competing Borel singularity?

Recompute the amplitude extraction at (0.2,1.0) at HIGH precision and across b, and
also report the sub-dominant singularity distances. If amp/pred converges to ~0.99
at high dps -> it was precision (b-scaling OK). If it stays garbage -> a second
singularity competes at larger b (scope limit on the single-singularity formula).
"""
import sys; sys.path.insert(0, 'debug')
import mpmath as mp
from _fastgl import fast_gl
import _watson_fibre as W

def Pcs(cs, k):
    d = mp.sqrt(cs * k * k + 1); return cs * mp.e ** (-d) * (d ** -3 + 3 * d ** -4 + 3 * d ** -5)

def dcoeffs(cs, A, M):
    k0 = max(mp.mpf(60), 35 / mp.sqrt(cs)); ks = [k0 * (i + 1) for i in range(M + 1)]
    gs = [Pcs(cs, k) * mp.e ** (A * k) * k ** 3 for k in ks]; V = mp.matrix(M + 1, M + 1)
    for i, k in enumerate(ks):
        for m in range(M + 1): V[i, m] = 1 / k ** m
    return list(mp.lu_solve(V, mp.matrix(gs)))

def m_hi(cs, b, n, d, A, xs, ws, K, M=12):
    bnd = mp.mpf(0)
    for x, w in zip(xs, ws):
        k = K * (x + 1) / 2; wk = K * w / 2
        j0 = mp.sin(k * b) / (k * b) if k * b > 1e-40 else mp.mpf(1)
        bnd += wk * (k ** (2 * n)) * j0 * Pcs(cs, k)
    z = A - 1j * b; t = mp.mpc(0)
    for j in range(M + 1):
        o = 2 * n - 3 - j; t += d[j] * z ** (-o) * mp.gammainc(o, z * K)
    return bnd + (t / b).imag

def run(cs, b, dps, K, Nq, N, windows):
    mp.mp.dps = dps
    cs = mp.mpf(cs); b = mp.mpf(b); A = mp.sqrt(cs)
    d = dcoeffs(cs, A, 12); Fn = W.F_taylor(N); xs, ws = fast_gl(Nq)
    a = [Fn[n] * m_hi(cs, b, n, d, A, xs, ws, mp.mpf(K)) for n in range(N + 1)]
    zst = -(A - 1j * b) ** 2; absz = abs(zst); th = mp.arg(zst)
    pred = (cs + b * b) ** (mp.mpf(3) / 2) / (4 * mp.sqrt(mp.pi) * mp.sqrt(cs) * b)
    # dominant vs sub-dominant Borel singularities (branch points of the fibre):
    # z*_1 = -(A-ib)^2 (source+scale), and the pure-scale one at |z|=? report |z*| and 4c1 (the x=-1 image ~ 2A scale)
    print(f"(c1,b)=({mp.nstr(cs,4)},{mp.nstr(b,4)}) dps={dps} K={K} Nq={Nq}: |z*|={mp.nstr(absz,8)}  pred={mp.nstr(pred,10)}", flush=True)
    for nlo, nhi in windows:
        rows = []; rhs = []
        for n in range(nlo, nhi + 1):
            y = a[n] * absz ** n / mp.factorial(2 * n) * mp.mpf(n) ** (mp.mpf(5) / 2)
            c = mp.cos(n * th); s = mp.sin(n * th); rows.append([c, s, c / n, s / n]); rhs.append(y)
        Mx = mp.matrix(rows); sol = mp.lu_solve(Mx.T * Mx, Mx.T * mp.matrix(rhs))
        amp = mp.sqrt(sol[0] ** 2 + sol[1] ** 2)
        print(f"    window {(nlo,nhi)}: amp={mp.nstr(amp,10)}  ratio={mp.nstr(amp/pred,8)}", flush=True)

# (A) the blowup point at escalating precision -- precision or model?
run('0.2', '1.0', 60, 150, 3200, 54, [(30,42),(36,48),(40,54)])
# (B) a robustness check: intermediate b at good precision (does b-scaling hold between 0.5 and 1.0?)
run('0.2', '0.7', 60, 140, 3000, 54, [(30,42),(36,48),(40,54)])
