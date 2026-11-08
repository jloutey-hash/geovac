"""Minimal truncated power-series arithmetic over mpf (lists of coefficients, index=power)."""
import mpmath as mp

def smul(a, b, M):
    out = [mp.mpf(0)]*(M+1)
    for i, ai in enumerate(a):
        if i > M or ai == 0: continue
        for j, bj in enumerate(b):
            if i+j > M: break
            out[i+j] += ai*bj
    return out

def sinv(a, M):
    """1/a for a[0] != 0."""
    out = [mp.mpf(0)]*(M+1)
    out[0] = 1/a[0]
    for n in range(1, M+1):
        s = mp.mpf(0)
        for j in range(1, n+1):
            if j < len(a): s += a[j]*out[n-j]
        out[n] = -s/a[0]
    return out

def sexp(a, M):
    """exp(a) for a[0] = 0."""
    out = [mp.mpf(0)]*(M+1); out[0] = mp.mpf(1)
    for n in range(1, M+1):
        s = mp.mpf(0)
        for j in range(1, n+1):
            if j < len(a): s += j*a[j]*out[n-j]
        out[n] = s/n
    return out
