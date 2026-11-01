"""Paper 59 sec:density_ring: the order-19 tetrahedron operator for the off-slice correlated triangle Y, and its
singular points, at (a, b, g1, g2) = (1, 2, 3/2, 3).

Record docs/predictions/2026-10-23_price_formula_slots.md (status note 2, disposition addendum 2: PF2 Y clause).
Written as a separate step from the results; no debug/ import.  Data tests/data/paper59_offslice_y_2026_11_01.json
(built by debug/build_offslice_y_guard_data_2026_11_01.py): the exact operator L for S = Y/t^3 in x = t^2 as
ascending integer coefficients of polys[k](x), L = sum_k polys[k](x) (d/dx)^k; theta_t^k Y, k = 0..19, at
t = 0.8, 1.7, 2.6, 3.8 from the 60-digit grid.  Y's minimal operator M13 is read from the neighbouring guard's data,
tests/data/paper59_y_2026_10_26.json (theta_t form on Y).

Six guards, each with the wrong answer it rejects.
(1) L has order 19, x-degree 102 + k in the coefficient of (d/dx)^k, integers up to 934 bits.  Rejected: another
    order or degree range.
(2) The leading coefficient factors over Z as x^7 (4x-1)(x-4)(4x-49)(x-16) [(x-1)(x-9)(4x-9)(4x-81)]^3 f, f
    irreducible of degree 98.  Rejected: another multiplicity, a ninth cut factor (e.g. 4x - 25), f reducible.
(3) L annihilates Y at t = 0.8, 1.7, 2.6, 3.8, below 1e-40 of the largest t^(2i) group (the theta form is rebuilt
    here from the D form); the dominant monomial times 1 + 1/200 leaves above 1e-2.  Rejected: an operator that
    does not annihilate Y (the perturbed one).
(4) The roots of f are apparent: at a root of f in F_p (two primes at which f has a linear factor) the 19 local
    exponents are distinct non-negative integers and no logarithm is forced, so all 19 solutions are analytic.
    Rejected: "a root of f is a genuine point of L" (the same routine finds the cut x = 9 non-trivial).
(5) Each of the eight cuts is genuine for Y's minimal operator M13: its exact local exponents there are
    non-integer, or integer with a forced logarithm (mod two primes).  Independent of the record's loop test
    (numerical continuation of Y, not re-run here).  Rejected: "x0 is apparent for M" at any cut; the routine
    returns trivial monodromy at the ordinary control x = 3.
(6) M13 divides L on the right (the step behind the paper's first arm: every genuine point of M is one of L): 13
    independent power-series solutions of M13 at x = 3, mod p, are annihilated by L to 130 orders.  Rejected: a
    non-factor (M13 with one coefficient changed mod p leaves non-zero terms).
Not re-run here: the reduction (34 masters), the lift from two primes and the check on a third, the loop test.
Fire-tested 2026-11-01 (debug/fire_offslice_y_guard_2026_11_01.sh), 11/11 FIRED, each plant against its wrong
answer: P1 order 18 claimed; P2 multiplicity 2 on (x - 9); P3 a ninth cut (4x - 25) claimed; P4 f claimed
reducible (49 + 49); P5 a non-annihilating operator (leading coefficient times 201); P6 the 1/200 perturbation
removed from the control; P7 apparentness tested at the cut x = 9 instead of a root of f; P8 the ordinary
x = 25/4 claimed genuine in place of 49/4; P9 the log obstruction ignored (log-type cuts 1/4, 4, 49/4, 16 read
as apparent); P10 M13 + 1 in its lowest coefficient as the right factor; P11 L19 + 1 in its leading coefficient
in the division test.
Measured here (2026-11-01): at the cuts M13's exponents are half-integral at 1, 9/4, 9, 81/4 and integral with a
forced log at 1/4, 4, 49/4, 16; at roots of f the exponents of L are 0..17, 19.  Run time about 2 s.
"""
from __future__ import annotations

import json
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import pytest
from flint import fmpq_poly, fmpz_poly, nmod_poly

HERE = Path(__file__).parent / "data"
DATA = json.loads((HERE / "paper59_offslice_y_2026_11_01.json").read_text())
YDATA = json.loads((HERE / "paper59_y_2026_10_26.json").read_text())
PRIMES = [1073741789, 1073740783, 1073741827, 1073741831, 1073741833, 1073741839, 1073741843, 1073741857,
          1073741891, 1073741909, 1073741939, 1073741953, 1073741969, 1073741971, 1073741987, 1073741993]
CUTS = [Fr(s) for s in DATA["cuts_x"]]
ORDINARY = [Fr(s) for s in DATA["ordinary_x"]]

L19 = [fmpz_poly([int(v) for v in p]) for p in DATA["polys"]]


# ---------------------------------------------------------------- operator conversions
def stirling2(k: int, j: int) -> int:
    if k == j:
        return 1
    if j == 0 or j > k:
        return 0
    return j * stirling2(k - 1, j) + stirling2(k - 1, j - 1)


def theta_t_on_y_to_d_on_s(c):
    """M = sum_k c_k(x) theta_t^k acting on Y; with Y = x^(3/2) S and theta_t = 2 theta_x,
    x^(-3/2) M x^(3/2) = sum_k c_k(x) (2 theta_x + 3)^k; theta_x^m = sum_j S(m, j) x^j D^j.  Returns D-form."""
    R = len(c) - 1
    out = [fmpz_poly([]) for _ in range(R + 1)]
    for k, ck in enumerate(c):
        # (2 th + 3)^k = sum_m binom(k, m) 2^m 3^(k-m) th^m
        from math import comb
        for m in range(k + 1):
            w = comb(k, m) * 2 ** m * 3 ** (k - m)
            for j in range(m + 1):
                s = stirling2(m, j)
                if s:
                    out[j] = out[j] + ck * fmpz_poly([0] * j + [w * s])
    return out


def d_on_s_to_theta_t_on_y(polys):
    """x^19 x^(3/2) sum_k p_k(x) D^k x^(-3/2) = sum_k p_k(x) x^(19-k) prod_{j<k} (u/2 - 3/2 - j), u = theta_t.
    Returns {(theta power, x power): Fraction}."""
    R = len(polys) - 1
    out: dict = {}
    for k, p in enumerate(polys):
        P = [Fr(1)]                                  # coefficients in u, ascending
        for j in range(k):
            a = Fr(-3, 2) - j                        # multiply by (u/2 + a)
            Q = [Fr(0)] * (len(P) + 1)
            for i, v in enumerate(P):
                Q[i] += v * a
                Q[i + 1] += v / 2
            P = Q
        for i, cv in enumerate(p.coeffs()):
            cv = int(cv)
            if cv:
                for m, v in enumerate(P):
                    if v:
                        key = (m, i + R - k)
                        out[key] = out.get(key, Fr(0)) + v * cv
    return out


M13 = theta_t_on_y_to_d_on_s([fmpz_poly([int(v) for v in c]) for c in YDATA["M13"]])


# ---------------------------------------------------------------- (1), (2): shape and leading coefficient
def test_operator_shape():
    """Rejects another order, another degree range, or another coefficient size than the record's."""
    assert len(L19) - 1 == 19
    assert [p.degree() for p in L19] == [102 + k for k in range(20)]
    bits = max(abs(int(v)).bit_length() for p in L19 for v in p.coeffs())
    assert bits == 934


def test_leading_coefficient_factorization():
    """Rejects another multiplicity at a cut, a ninth cut factor, a missing x^7, and a reducible f."""
    expect = {(0, 1): 7, (-1, 4): 1, (-4, 1): 1, (-49, 4): 1, (-16, 1): 1,
              (-1, 1): 3, (-9, 1): 3, (-9, 4): 3, (-81, 4): 3}
    _, facs = L19[-1].factor()
    lin, big = {}, []
    for g, m in facs:
        c = [int(v) for v in g.coeffs()]
        if g.degree() == 1:
            if c[1] < 0:
                c = [-c[0], -c[1]]
            lin[(c[0], c[1])] = m
        else:
            big.append((g.degree(), m))
    assert lin == expect
    assert big == [(98, 1)]


# ---------------------------------------------------------------- (3): L on Y
THETA = d_on_s_to_theta_t_on_y(L19)


def residual(t: str, perturb: bool = False):
    mp.mp.dps = 80
    th = [mp.mpf(v) for v in DATA["thetaY"][t]]
    t2 = mp.mpf(t) ** 2
    mono = {}
    for (m, i), v in THETA.items():
        mono[(m, i)] = mp.mpf(v.numerator) / v.denominator * t2 ** i * th[m]
    if perturb:
        kmax = max(mono, key=lambda key: abs(mono[key]))
        mono[kmax] *= mp.mpf(201) / 200
    D = max(i for _, i in mono)
    groups = [mp.fsum(v for (m, i), v in mono.items() if i == j) for j in range(D + 1)]
    return abs(mp.fsum(groups)) / max(abs(g) for g in groups)


@pytest.mark.parametrize("t", ["0.8", "1.7", "2.6", "3.8"])
def test_l19_annihilates_y(t):
    """Rejects an operator that does not annihilate Y: the dominant monomial times 1 + 1/200 must leave > 1e-2."""
    assert residual(t) < mp.mpf(10) ** -40
    assert residual(t, perturb=True) > mp.mpf(10) ** -2


# ---------------------------------------------------------------- local analysis mod p
def shifted(polys, x0: int, p: int):
    """q_k(s) = polys[k](x0 + s) mod p, as lists of ints."""
    s = nmod_poly([x0 % p, 1], p)
    out = []
    for q in polys:
        qp = nmod_poly([int(v) % p for v in q.coeffs()], p)
        out.append([int(v) for v in qp(s).coeffs()])
    return out


def ff(n: int, k: int, p: int) -> int:
    r = 1
    for i in range(k):
        r = r * (n - i) % p
    return r


def local_trivial_modp(polys, x0: int, p: int, exps=None):
    """Local monodromy at x0 (an element of F_p) of a Fuchsian D-form operator, from the reduction mod p.
    exps: the local exponents as Fractions (exact), or None to find them as integers in [-60, 400] mod p.
    Returns 'nonint', 'repeated', 'log' (non-trivial monodromy) or 'trivial'."""
    q = shifted(polys, x0, p)
    R = len(q) - 1
    val = [next((j for j, v in enumerate(c) if v), None) for c in q]
    delta = min(v - k for k, v in enumerate(val) if v is not None)
    assert val[R] - R == delta, "not Fuchsian at this point"

    def ind(n):
        return sum(q[k][val[k]] * ff(n, k, p) for k in range(R + 1)
                   if val[k] is not None and val[k] - k == delta) % p

    if exps is None:
        ints = [n for n in range(-60, 400) if ind(n) == 0]
        # the indicial polynomial has degree R; confirm these are all its roots (with multiplicity one each)
        if len(ints) != R:
            return "repeated" if 0 < len(ints) < R else "nonint"
        e = sorted(ints)
    else:
        if any(x.denominator != 1 for x in exps):
            return "nonint"
        e = sorted(int(x) for x in exps)
        if len(set(e)) < len(e):
            return "repeated"
    rho = e[0]
    res = [x - rho for x in e]
    N = res[-1]
    a = []
    for n in range(N + 1):
        rhs = [0] * R
        for k in range(R + 1):
            ck = q[k]
            for d in range(1, n + 1):
                j = d + k + delta
                if 0 <= j < len(ck) and ck[j]:
                    w = ck[j] * ff(rho + n - d, k, p) % p
                    prev = a[n - d]
                    for i in range(R):
                        if prev[i]:
                            rhs[i] = (rhs[i] + w * prev[i]) % p
        if n in res:
            assert ind(rho + n) == 0
            if any(rhs):
                return "log"
            v = [0] * R
            v[res.index(n)] = 1
            a.append(v)
        else:
            I = ind(rho + n)
            assert I, "accidental zero of the indicial polynomial mod p"
            inv = pow(I, -1, p)
            a.append([(-r * inv) % p for r in rhs])
    return "trivial"


def exact_exponents(polys, x0: Fr):
    """Local exponents at the rational point x0, exactly: roots of the indicial polynomial over Q."""
    lin = fmpz_poly([-x0.numerator, x0.denominator])
    vals, lead = [], []
    for c in polys:
        if c.is_zero():
            vals.append(None)
            lead.append(None)
            continue
        v, r = 0, c
        while True:
            qq, rr = divmod(r, lin)
            if not rr.is_zero():
                break
            r, v = qq, v + 1
        vals.append(v)
        lead.append(sum(Fr(int(cv)) * x0 ** i for i, cv in enumerate(r.coeffs())) * x0.denominator ** v)
    R = len(polys) - 1
    delta = min(v - k for k, v in enumerate(vals) if v is not None)
    assert vals[R] - R == delta, "not Fuchsian"
    I = fmpq_poly([0])
    for k in range(R + 1):
        if vals[k] is not None and vals[k] - k == delta:
            f = fmpq_poly([1])
            for i in range(k):
                f = f * fmpq_poly([-i, 1])
            I = I + f * fmpq_poly([lead[k].numerator]) / lead[k].denominator
    assert I.degree() == R
    _, facs = I.factor()
    exps = []
    for g, m in facs:
        if g.degree() == 1:
            c = g.coeffs()
            exps += [Fr(-int(c[0].p) * int(c[1].q), int(c[0].q) * int(c[1].p))] * m
        else:
            exps += [Fr(1, 10 ** 9)] * (g.degree() * m)       # irrational: marked non-integer
    return sorted(exps)


def red(x: Fr, p: int) -> int:
    return x.numerator * pow(x.denominator, -1, p) % p


# ---------------------------------------------------------------- (4): f apparent
def int_exponents_modp(polys, x0: int, p: int):
    """Integer roots in [-60, 400] of the indicial polynomial at x0, mod p."""
    q = shifted(polys, x0, p)
    R = len(q) - 1
    val = [next((j for j, v in enumerate(c) if v), None) for c in q]
    delta = min(v - k for k, v in enumerate(val) if v is not None)
    return [n for n in range(-60, 400)
            if sum(q[k][val[k]] * ff(n, k, p) for k in range(R + 1)
                   if val[k] is not None and val[k] - k == delta) % p == 0]


def f_factor():
    for g, m in L19[-1].factor()[1]:
        if g.degree() == 98:
            return g
    raise AssertionError("no degree-98 factor")


def test_roots_of_f_apparent():
    """Rejects "a root of f is a genuine point of L": at roots of f in F_p (two primes with a linear factor of f)
    the 19 exponents must be distinct integers >= 0 with no forced log; the cut x = 9 must come out non-trivial."""
    f = f_factor()
    found = 0
    for p in PRIMES:
        roots = [r for r, _ in nmod_poly([int(v) % p for v in f.coeffs()], p).roots()]
        if not roots:
            continue
        x0 = int(roots[0])
        assert local_trivial_modp(L19, x0, p) == "trivial"
        assert int_exponents_modp(L19, x0, p) == list(range(18)) + [19]
        found += 1
        if found == 2:
            break
    assert found == 2
    p = PRIMES[0]
    assert local_trivial_modp(L19, red(Fr(9), p), p, exact_exponents(L19, Fr(9))) != "trivial"


# ---------------------------------------------------------------- (5): the eight cuts genuine for M13
GENUINE = CUTS


@pytest.mark.parametrize("x0", [str(c) for c in GENUINE])
def test_cut_genuine_for_minimal_operator(x0):
    """Rejects "x0 is apparent (trivial local monodromy) for Y's minimal operator" at any of the eight cuts."""
    x0 = Fr(x0)
    e = exact_exponents(M13, x0)
    for p in PRIMES[:2]:
        assert local_trivial_modp(M13, red(x0, p), p, e) != "trivial"


def test_ordinary_control_for_minimal_operator():
    """The routine of (5) must return trivial monodromy where M13 and L are ordinary (x = 3, 25/4, 25)."""
    for x0 in ORDINARY:
        assert sum(Fr(int(c)) * x0 ** i for i, c in enumerate(L19[-1].coeffs())) != 0
        e = exact_exponents(M13, x0)
        assert e == [Fr(i) for i in range(13)]
        assert local_trivial_modp(M13, red(x0, PRIMES[0]), PRIMES[0], e) == "trivial"


# ---------------------------------------------------------------- (6): M13 divides L19 on the right
def series_solutions(polys, x0: int, p: int, N: int):
    """The R power-series solutions of an operator ordinary at x0, mod p, to N terms (initial data a_i = delta)."""
    q = shifted(polys, x0, p)
    R = len(q) - 1
    assert q[R][0] % p
    sols = [[0] * N for _ in range(R)]
    for i in range(R):
        sols[i][i] = 1
    for n in range(N - R):
        # coefficient of s^n: sum_k sum_j q_{k,j} ff(n - j + k, k) a_{n-j+k} = 0; unknown a_{n+R} from k=R, j=0
        lead = q[R][0] * ff(n + R, R, p) % p
        inv = pow(lead, -1, p)
        for sidx in range(R):
            a = sols[sidx]
            tot = 0
            for k in range(R + 1):
                ck = q[k]
                for j in range(len(ck)):
                    m = n - j + k
                    if m < 0:
                        break
                    if (k == R and j == 0) or m >= N or not ck[j]:
                        continue
                    if a[m]:
                        tot += ck[j] * ff(m, k, p) * a[m]
            a[n + R] = (-tot * inv) % p
    return sols


def apply_op(polys, x0: int, p: int, series, n_check: int):
    q = shifted(polys, x0, p)
    S = nmod_poly(series, p)
    acc = nmod_poly([], p)
    d = S
    for k, c in enumerate(q):
        if k:
            d = d.derivative()
        acc = acc + nmod_poly(c, p) * d
    co = [int(v) for v in acc.coeffs()] + [0] * (n_check + 1)
    return co[:n_check]


def test_m13_right_divides_l19():
    """Rejects "M13 is not a right factor of L": every solution of M13 at the ordinary point x = 3 must be
    annihilated by L to 130 orders (a non-factor leaves non-zero terms)."""
    p = PRIMES[0]
    M = list(M13)
    N = 150
    sols = series_solutions(M, 3, p, N)
    for s in sols:
        out = apply_op(L19, 3, p, s, N - 20)
        assert not any(out)
