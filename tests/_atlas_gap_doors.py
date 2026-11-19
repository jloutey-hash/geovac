"""Helper for tests/test_atlas_gap_doors.py (proposed atlas doors 28-30, memo debug/atlas_gap_memo.md).

Copied (not imported) from the memo's scripts debug/atlas_gap_field_series.py, debug/atlas_gap_richardson.py,
debug/atlas_gap_linear_stark.py, debug/atlas_gap_lens.py and debug/atlas_gap_lens_F.py, so the guard does not
depend on the transient debug/ tree.  Everything exact unless the name says otherwise.

Field series (doors 28, 29).  psi_k = e^{-r} sum c_ij r^i P_j(cos theta).  With L = H0 + 1/2,
    L[e^{-r} r^i P_j] = e^{-r} [ i r^{i-1} + (j(j+1) - i(i+1))/2 r^{i-2} ] P_j,
triangular in i, so each Rayleigh-Schroedinger order is solved top-down in exact rationals
(intermediate normalization).  Stark V = r u (series in F); diamagnetic V = (1/8) r^2 (1 - u^2)
(series in beta = B^2).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import comb, factorial
from typing import Callable, Dict, List, Sequence, Tuple

import mpmath as mp
import sympy as sp

Poly = Dict[Tuple[int, int], Fr]


# ----------------------------------------------------------------------------------------------
# Rayleigh-Schroedinger series for the hydrogen ground state in a uniform field
# ----------------------------------------------------------------------------------------------

def u_times(poly: Poly) -> Poly:
    """Multiply by u = cos(theta) in the Legendre index: u P_j = ((j+1) P_{j+1} + j P_{j-1})/(2j+1)."""
    out: Poly = {}
    for (i, j), c in poly.items():
        out[(i, j + 1)] = out.get((i, j + 1), Fr(0)) + c * Fr(j + 1, 2 * j + 1)
        if j > 0:
            out[(i, j - 1)] = out.get((i, j - 1), Fr(0)) + c * Fr(j, 2 * j + 1)
    return {k: v for k, v in out.items() if v != 0}


def r_times(poly: Poly, p: int) -> Poly:
    return {(i + p, j): c for (i, j), c in poly.items()}


def add(a: Poly, b: Poly, s: int = 1) -> Poly:
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, Fr(0)) + s * v
    return {k: v for k, v in out.items() if v != 0}


def scale(a: Poly, s: Fr) -> Poly:
    return {k: v * s for k, v in a.items() if v * s != 0}


def V_stark(poly: Poly) -> Poly:
    """z = r cos(theta)."""
    return u_times(r_times(poly, 1))


def V_diamag(poly: Poly) -> Poly:
    """(1/8) r^2 sin^2(theta) = (1/8) r^2 (1 - u^2); multiplies beta = B^2."""
    r2 = r_times(poly, 2)
    return scale(add(r2, u_times(u_times(r2)), -1), Fr(1, 8))


def solve_L(rhs: Poly) -> Poly:
    """Solve L g = rhs, g = sum a_ij r^i P_j, each j separately (top-down); assert solvability."""
    byj: Dict[int, Dict[int, Fr]] = {}
    for (i, j), c in rhs.items():
        byj.setdefault(j, {})[i] = c
    sol: Poly = {}
    for j, b in byj.items():
        top = max(b)
        a: Dict[int, Fr] = {}

        def cst(i: int, j: int = j) -> Fr:
            return Fr(j * (j + 1) - i * (i + 1), 2)

        # coefficient of r^m: (m+1) a_{m+1} + cst(m+2) a_{m+2} = b_m
        for m in range(top, j - 2, -1):
            if m + 1 < j:
                break
            val = b.get(m, Fr(0)) - cst(m + 2) * a.get(m + 2, Fr(0))
            if m + 1 == 0:
                continue
            a[m + 1] = val / (m + 1)
        lowest = j - 2 if j >= 1 else -1
        resid = (b.get(lowest, Fr(0)) - cst(lowest + 2) * a.get(lowest + 2, Fr(0))
                 - (lowest + 1) * a.get(lowest + 1, Fr(0)))
        if j == 0:
            resid = b.get(-1, Fr(0)) - cst(1) * a.get(1, Fr(0))
        assert resid == 0, (j, resid)
        for i, v in a.items():
            if v != 0:
                sol[(i, j)] = v
    return sol


def overlap0(poly: Poly) -> Fr:
    """<e^{-r}| e^{-r} poly> / <e^{-r}|e^{-r}>; only P_0 survives; int r^{i+2} e^{-2r} = (i+2)!/2^{i+3}."""
    tot = Fr(0)
    for (i, j), c in poly.items():
        if j == 0:
            tot += c * Fr(factorial(i + 2), 2 ** (i + 3))
    return tot / Fr(1, 4)


def series(V: Callable[[Poly], Poly], K: int) -> List[Fr]:
    """Energy coefficients E_0..E_{K+1} (E_k multiplies lambda^k, V scaled by lambda)."""
    psi: List[Poly] = [{(0, 0): Fr(1)}]
    E: List[Fr] = [Fr(-1, 2)]
    for k in range(1, K + 1):
        Ek = overlap0(V(psi[k - 1]))
        E.append(Ek)
        rhs = scale(V(psi[k - 1]), Fr(-1))
        for jj in range(1, k + 1):
            if E[jj] != 0:
                rhs = add(rhs, scale(psi[k - jj], E[jj]))
        g = solve_L(rhs)
        g[(0, 0)] = g.get((0, 0), Fr(0)) - overlap0(g)
        psi.append({kk: v for kk, v in g.items() if v != 0})
    E.append(overlap0(V(psi[K])))
    return E


def richardson(xs: Sequence[mp.mpf], ks: Sequence[int], N: int) -> mp.mpf:
    """N-th order Richardson at the last point, x_k = L + c1/k + ... (consecutive k)."""
    ks = list(ks)[-N - 1:]
    xs = list(xs)[-N - 1:]
    tot = mp.mpf(0)
    for i in range(N + 1):
        tot += xs[i] * mp.mpf(ks[i]) ** N * (-1) ** (i + N) / (mp.factorial(i) * mp.factorial(N - i))
    return tot


def stark_ratios(E: List[Fr], kmax: int) -> Tuple[List[int], List[mp.mpf]]:
    """a_k / [-(6/pi)(3/2)^{2k}(2k)!], a_k the coefficient of F^{2k}."""
    ks = list(range(1, kmax + 1))
    xs = [mp.mpf(E[2 * k].numerator) / E[2 * k].denominator
          / (-(6 / mp.pi) * mp.mpf(3) ** (2 * k) / mp.mpf(2) ** (2 * k) * mp.factorial(2 * k)) for k in ks]
    return ks, xs


def diamag_ratios(E: List[Fr], kmax: int) -> Tuple[List[int], List[mp.mpf]]:
    """a_k / [(-1)^{k+1} pi^{-2k} Gamma(2k + 3/2)], a_k the coefficient of beta^k; tends to (4/pi)^{5/2}."""
    ks = list(range(2, kmax + 1))
    xs = [mp.mpf(E[k].numerator) / E[k].denominator
          / ((-1) ** (k + 1) * mp.pi ** (-2 * k) * mp.gamma(2 * k + mp.mpf(3) / 2)) for k in ks]
    return ks, xs


# ----------------------------------------------------------------------------------------------
# First-order (linear) Stark matrices in a shell
# ----------------------------------------------------------------------------------------------

_r = sp.symbols("r", positive=True)


def z_block(n: int, m: int) -> sp.Matrix:
    """<n l1 m| r cos(theta) |n l2 m>, l1, l2 = |m|..n-1, exact (radial integral x Gaunt)."""
    from sympy.physics.hydrogen import R_nl
    from sympy.physics.wigner import gaunt

    ls = list(range(abs(m), n))

    def elem(l1: int, l2: int) -> sp.Expr:
        if abs(l1 - l2) != 1:
            return sp.Integer(0)
        rad = sp.integrate(R_nl(n, l1, _r, 1) * R_nl(n, l2, _r, 1) * _r ** 3, (_r, 0, sp.oo))
        ang = (-1) ** m * gaunt(l1, 1, l2, -m, 0, m)
        return sp.nsimplify(sp.simplify(rad * sp.sqrt(4 * sp.pi / 3) * ang))

    return sp.Matrix(len(ls), len(ls), lambda i, j: elem(ls[i], ls[j]))


def parabolic_shifts(n: int, m: int) -> List[sp.Rational]:
    """(3/2) n (n1 - n2) with n1 + n2 = n - |m| - 1."""
    s = n - abs(m) - 1
    return sorted(sp.Rational(3, 2) * n * (n1 - (s - n1)) for n1 in range(s + 1))


# ----------------------------------------------------------------------------------------------
# Lens spaces L(p,1)
# ----------------------------------------------------------------------------------------------

def mult_weights(p: int, n: int) -> int:
    """Route 1: count the n^2 parabolic states (n1, n2, m) of shell n with L_z + A_z = m + (n1 - n2) = 0 mod p."""
    cnt = 0
    for m in range(-(n - 1), n):
        s = n - abs(m) - 1
        for n1 in range(s + 1):
            n2 = s - n1
            if (m + n1 - n2) % p == 0:
                cnt += 1
    return cnt


def mult_character(p: int, n: int) -> int:
    """Route 2: (1/p) sum_t n chi_j(2 pi t/p), chi_j(theta) = sin(n theta)/sin(theta), exact cyclotomics."""
    w = sp.exp(2 * sp.pi * sp.I / p)
    tot = sp.Integer(0)
    for t in range(p):
        # chi_j(theta) = sum_{k = -(n-1), -(n-3), ..., n-1} e^{i k theta}; at theta = 2 pi t/p, e^{i theta} = w^t
        chi = sum(w ** (t * k) for k in range(-(n - 1), n, 2))
        tot += n * chi
    val = sp.nsimplify(sp.simplify(sp.expand_complex(tot / p)))
    assert val.is_Integer, val
    return int(val)


def molien_coeffs(p: int, order: int) -> List[int]:
    """Route 3: (1 - t^2) x group average of 1/det(1 - t g) on R^4 = C^2, as a power series (exact integers).
    det(1 - t g)^{-1} = (1 - t w)^{-2} (1 - t w^{-1})^{-2} = sum (a+1)(b+1) t^{a+b} w^{a-b}; averaging keeps a = b mod p."""
    c = [0] * (order + 3)
    for a in range(order + 3):
        for b in range(order + 3 - a):
            if (a - b) % p == 0:
                c[a + b] += (a + 1) * (b + 1)
    return [c[k] - (c[k - 2] if k >= 2 else 0) for k in range(order + 1)]


T = sp.symbols("t")

#: generating functions sum_n m_n t^{n-1} as stated in the memo's table
GEN_FUNCS: Dict[int, sp.Expr] = {
    2: (1 + 6 * T ** 2 + T ** 4) / (1 - T ** 2) ** 3,
    3: (1 + T) * (1 - 2 * T + 5 * T ** 2 - 2 * T ** 3 + T ** 4) / ((1 - T) ** 3 * (1 + T + T ** 2) ** 2),
    5: (1 + T) * (1 - 2 * T + 5 * T ** 2 - 8 * T ** 3 + 13 * T ** 4 - 8 * T ** 5 + 5 * T ** 6 - 2 * T ** 7 + T ** 8)
       / ((1 - T) ** 3 * (1 + T + T ** 2 + T ** 3 + T ** 4) ** 2),
}


def series_coeffs(G: sp.Expr, order: int) -> List[int]:
    s = sp.series(G, T, 0, order + 1).removeO()
    return [int(s.coeff(T, k)) for k in range(order + 1)]


def c_count(p: int, n: int) -> int:
    """c(n) = m_n / n: number of k in {-(n-1), ..., n-1} (step 2) divisible by p."""
    return sum(1 for k in range(-(n - 1), n, 2) if k % p == 0)


def F_lens_numeric(p: int, dps: int = 70) -> mp.mpf:
    """F_p = sum_n m_n / n^4 = sum c(n)/n^3 via Hurwitz zeta on residue classes mod 2p.
    Requires c(n) = alpha_a n + beta_a on each class n = a mod 2p (checked by the caller)."""
    q = 2 * p
    with mp.workdps(dps):
        tot = mp.mpf(0)
        for a in range(1, q + 1):
            c1, c2 = c_count(p, a), c_count(p, a + q)
            alpha = mp.mpf(c2 - c1) / q
            beta = c1 - alpha * a
            tot += alpha * mp.zeta(2, mp.mpf(a) / q) / q ** 2 + beta * mp.zeta(3, mp.mpf(a) / q) / q ** 3
        return +tot


def cot_poly_power_sum(p: int) -> Fr:
    """Exact sum of squares of the roots of P(x) = Im((x + i)^p) / 1, whose roots are cot(pi k/p), k = 1..p-1.
    P(x) = sum_k (-1)^k C(p, 2k+1) x^{p-1-2k}; e1 = 0, e2 = -C(p,3)/C(p,1); sum x^2 = e1^2 - 2 e2."""
    e1 = Fr(0)
    e2 = Fr(-comb(p, 3), comb(p, 1))
    return e1 * e1 - 2 * e2


def cot_poly_value(p: int, x: mp.mpf) -> mp.mpf:
    return sum((-1) ** k * comb(p, 2 * k + 1) * x ** (p - 1 - 2 * k) for k in range((p - 1) // 2 + 1))
