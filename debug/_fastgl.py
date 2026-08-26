"""Fast Gauss-Legendre nodes/weights at mpmath precision: numpy double-precision seed +
2-4 high-precision Newton steps via the O(N) Legendre recurrence. ~100x faster than
Newton-on-mpmath.legendre for large N. Drop-in replacement for routeC_T2_highprec.gl."""
import mpmath as mp
import numpy as np

def _legendre_and_deriv(N, x):
    # P_N(x), P_N'(x) via recurrence (all mpf)
    p0 = mp.mpf(1); p1 = x
    if N == 0: return p0, mp.mpf(0)
    for k in range(2, N+1):
        p0, p1 = p1, ((2*k-1)*x*p1 - (k-1)*p0)/k
    pN = p1
    dP = N*(x*pN - p0)/(x*x - 1)
    return pN, dP

_CACHE = {}
def fast_gl(N):
    key = (N, mp.mp.dps)
    if key in _CACHE: return _CACHE[key]
    x0, w0 = np.polynomial.legendre.leggauss(N)   # double precision seed
    xs = []; ws = []
    for xd in x0:
        x = mp.mpf(float(xd))
        for _ in range(5):
            pN, dP = _legendre_and_deriv(N, x)
            dx = pN/dP; x -= dx
            if abs(dx) < mp.mpf(10)**(-(mp.mp.dps+6)): break
        pN, dP = _legendre_and_deriv(N, x)
        xs.append(x); ws.append(2/((1-x*x)*dP*dP))
    _CACHE[key] = (xs, ws)
    return xs, ws
