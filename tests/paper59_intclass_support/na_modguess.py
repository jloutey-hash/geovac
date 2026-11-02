"""Exact (mod p) search for the bond-angle operator of a Legendre-coefficient sequence (record NA).
Ansatz sum_{i<=R, j<=Q} a_ij x^j T^i f = 0, T = (1 - x^2) d/dx, with x and T acting locally on the coefficients
(as debug/ao_guess.py), every entry mod p.  The null space is exact: no noise floor, so no masking.
nullity(R, Q) = number of independent operators of order <= R and degree <= Q; for a true operator of order R and
degree Q0, nullity(R, Q) >= Q - Q0 + 1 for Q >= Q0 (x^m L), so the first Q with nullity >= 1 is the minimal degree.
"""
import flint


def X(s, p):
    out = []
    for l in range(len(s) - 1):
        v = flint.nmod(l, p) / flint.nmod(2 * l - 1, p) * s[l - 1] if l >= 1 else flint.nmod(0, p)
        out.append(v + flint.nmod(l + 1, p) / flint.nmod(2 * l + 3, p) * s[l + 1])
    return out


def T(s, p):
    out = []
    for l in range(len(s) - 1):
        v = flint.nmod((l - 1) * l, p) / flint.nmod(2 * l - 1, p) * s[l - 1] if l >= 1 else flint.nmod(0, p)
        out.append(flint.nmod((l + 1) * (l + 2), p) / flint.nmod(2 * l + 3, p) * s[l + 1] - v)
    return out


class Basis:
    """Cache of x^j T^i s, extended on demand."""

    def __init__(self, s, p):
        self.p = p
        self.t = [list(s)]
        self.cols = {}

    def col(self, i, j):
        if (i, j) in self.cols:
            return self.cols[(i, j)]
        while len(self.t) <= i:
            self.t.append(T(self.t[-1], self.p))
        c = self.t[i] if j == 0 else X(self.col(i, j - 1), self.p)
        self.cols[(i, j)] = c
        return c


def nullity(B, R, Q, margin):
    keys = [(i, j) for i in range(R + 1) for j in range(Q + 1)]
    cols = [B.col(i, j) for i, j in keys]
    nrow = min(len(c) for c in cols)
    if nrow < len(keys) + margin:
        return None, keys, None
    M = flint.nmod_mat(nrow, len(keys), [cols[k][l] for l in range(nrow) for k in range(len(keys))], B.p)
    N, nul = M.nullspace()
    return nul, keys, (N, nrow)


def min_degree(B, R, Qmax, margin):
    """Smallest Q <= Qmax with nullity >= 1 at order R (linear scan from Q = 0), or None."""
    for Q in range(Qmax + 1):
        nul, keys, info = nullity(B, R, Q, margin)
        if nul is None:
            return None, Q
        if nul >= 1:
            return Q, nul
    return None, Qmax


def min_degree_bisect(B, R, margin, nmax):
    """Smallest Q with nullity >= 1 at order R by bisection (nullity is non-decreasing in Q: the ansatz space grows),
    Q up to the largest the rows allow (nmax = number of coefficients); returns (Q0, nullity at Q0) or (None, Qtop)."""
    Qtop = 0
    while (R + 1) * (Qtop + 2) + margin <= nmax - R - (Qtop + 1):
        Qtop += 1
    nul, _, _ = nullity(B, R, Qtop, margin)
    if not nul:
        return None, Qtop
    lo, hi = -1, Qtop                     # nullity(lo) = 0 (or lo = -1), nullity(hi) >= 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        n, _, _ = nullity(B, R, mid, margin)
        if n:
            hi = mid
        else:
            lo = mid
    n, _, _ = nullity(B, R, hi, margin)
    return hi, n
