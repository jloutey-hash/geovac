"""Guess the linear ODE in the bond angle x = cos(gamma) satisfied by f(x) = sum_l F_l P_l(x), directly from
high-precision Legendre coefficients (record AO, docs/predictions/2026-10-28_angle_operator.md).

Ansatz  sum_{i<=R, j<=Q} a_ij x^j T^i f = 0,  T = (1 - x^2) d/dx.  Both x and T act locally on Legendre coefficients:
    (x f)_l = l/(2l-1) F_{l-1} + (l+1)/(2l+3) F_{l+1}
    (T f)_l = (l+1)(l+2)/(2l+3) F_{l+1} - (l-1) l/(2l-1) F_{l-1}
so every Legendre coefficient of the left side is a finite combination of the F_l (no truncated tail).  One
equation per l; the a_ij are a null vector.  A true ODE shows as a singular value at the working-precision floor,
separated from the next by many orders, and validates on rows held out of the fit.
"""
import mpmath as mp


def X(s):
    n = len(s) - 1
    out = []
    for l in range(n):
        a = mp.mpf(l) / (2 * l - 1) * s[l - 1] if l >= 1 else 0
        out.append(a + mp.mpf(l + 1) / (2 * l + 3) * s[l + 1])
    return out


def T(s):
    n = len(s) - 1
    out = []
    for l in range(n):
        a = mp.mpf((l - 1) * l) / (2 * l - 1) * s[l - 1] if l >= 1 else 0
        out.append(mp.mpf((l + 1) * (l + 2)) / (2 * l + 3) * s[l + 1] - a)
    return out


def basis(F, R, Q):
    cols = {}
    ti = list(F)
    for i in range(R + 1):
        xj = ti
        for j in range(Q + 1):
            cols[(i, j)] = xj
            xj = X(xj)
        ti = T(ti)
    return cols


def guess(F, R, Q, nfit):
    cols = basis(F, R, Q)
    keys = sorted(cols)
    nrow = min(len(c) for c in cols.values())
    rows = []
    for l in range(nrow):
        row = [cols[k][l] for k in keys]
        m = max(abs(v) for v in row)
        rows.append([v / m for v in row])
    fit, held = rows[:nfit], rows[nfit:]
    A = mp.matrix(fit)
    cs = [max(abs(A[i, j]) for i in range(A.rows)) or 1 for j in range(A.cols)]
    for j in range(A.cols):
        for i in range(A.rows):
            A[i, j] /= cs[j]
    U, S, V = mp.svd_r(A)
    k = min(range(len(S)), key=lambda i: S[i])
    s = sorted(S[i] for i in range(len(S)))
    v = [V[k, j] / cs[j] for j in range(A.cols)]
    res = [abs(sum(a * b for a, b in zip(r, v))) / sum(abs(a * b) for a, b in zip(r, v)) for r in held]
    return dict(R=R, Q=Q, keys=keys, v=v, smin=s[0], s2=s[1], n_fit=len(fit), n_held=len(held),
                held_max=max(res) if res else None)
