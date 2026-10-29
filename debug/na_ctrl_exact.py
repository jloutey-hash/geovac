"""Exact positive control for debug/na_modguess.py (record NA): f(x) = sum_{k=1..K} x^k e^{c_k x}, c_k = k/3 + 1/2,
a known D-finite function of minimal order K (LCLM of K first-order operators x f' = (k + c_k x) f).  Its Legendre
coefficients are exact combinations of e^{+-c_k} with rational coefficients, the shape of the integral's; the
sequence is mod p with the constants e^{+-c_k} replaced by random residues (a generic combination), and also returned
per component.
    F_l = (2l+1)/2 sum_k int_{-1}^{1} x^k e^{c x} P_l(x) dx,  int_{-1}^{1} x^n e^{c x} dx = e^c U_n - e^{-c} V_n.
"""
from fractions import Fraction

import flint


def legendre_coeffs(lmax):
    """P_l as lists of exact coefficients (fmpq), l = 0..lmax."""
    P = [[flint.fmpq(1)], [flint.fmpq(0), flint.fmpq(1)]]
    for l in range(1, lmax):
        a = [flint.fmpq(0)] + [flint.fmpq(2 * l + 1, l + 1) * c for c in P[l]]
        b = P[l - 1] + [flint.fmpq(0)] * (len(a) - len(P[l - 1]))
        P.append([u - flint.fmpq(l, l + 1) * v for u, v in zip(a, b)])
    return P


def moments(c, nmax):
    """int_{-1}^1 x^n e^{cx} dx = e^c U_n + e^{-c} V_n (exact U, V): I_n = ([x^n e^{cx}]_{-1}^{1} - n I_{n-1}) / c."""
    c = flint.fmpq(Fraction(c).numerator, Fraction(c).denominator)
    U, V = [1 / c], [-1 / c]
    for n in range(1, nmax + 1):
        U.append((1 - n * U[-1]) / c)
        V.append((-(-1) ** n - n * V[-1]) / c)
    return U, V


def control(K, lmax, p, weights):
    """Generic mod-p sequence and per-component sequences {('+'|'-', k): [...]}."""
    P = legendre_coeffs(lmax)
    comps = {}
    for k in range(1, K + 1):
        ck = Fraction(k, 3) + Fraction(1, 2)
        U, V = moments(ck, lmax + k)
        for sgn, W in (('+', U), ('-', V)):
            seq = []
            for l in range(lmax + 1):
                s = sum((a * W[n + k] for n, a in enumerate(P[l]) if a != 0), flint.fmpq(0))
                s = flint.fmpq(2 * l + 1, 2) * s
                seq.append(flint.nmod(int(s.p), p) / flint.nmod(int(s.q), p))
            comps[(sgn, k)] = seq
    gen = [sum((flint.nmod(weights[i], p) * comps[key][l] for i, key in enumerate(sorted(comps))), flint.nmod(0, p))
           for l in range(lmax + 1)]
    return gen, comps
