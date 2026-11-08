"""GT1-GT3 (docs/predictions/2026-10-22_correlated_triangle_g32.md): operator library.

Operators L = sum_k c_k(x) theta^k, theta = theta_t = t d/dt, x = t^2, coefficients function-left,
stored as a list of flint polynomials in x (fmpz_poly / fmpq_poly exact, nmod_poly mod p).
On functions of x, theta_t = 2 x d/dx.  The QG "Q-form" sum_i t^(2i) Q_i(theta) has c_k = sum_i [u^k]Q_i x^i.

Exact work over Q: composition, adjoint (theta* = -theta), pseudo right division (exact remainder test),
local indicial polynomials.  Mod p work: Frobenius series at x0, minimal-annihilator searches,
multi-prime CRT + rational reconstruction.
"""
from __future__ import annotations

import json
from fractions import Fraction
from math import comb
from pathlib import Path

import flint
import sympy as sp

HERE = Path(__file__).parent
DATA = HERE / "data"
u_sym = sp.Symbol("u")


# ----------------------------------------------------------------------------------------------
# polynomial helpers (generic over flint poly types)

def X_of(p):
    """the generator x in the ring of p"""
    if isinstance(p, flint.nmod_poly):
        return flint.nmod_poly([0, 1], p.modulus())
    return type(p)([0, 1])


def zero_like(p):
    if isinstance(p, flint.nmod_poly):
        return flint.nmod_poly([], p.modulus())
    return type(p)([])


def const_like(p, c):
    if isinstance(p, flint.nmod_poly):
        return flint.nmod_poly([c], p.modulus())
    return type(p)([c])


def Dth(p):
    """theta_t applied to a function of x: 2 x p'(x)"""
    if p.is_zero():
        return p
    return 2 * X_of(p) * p.derivative()


def trim(L):
    L = list(L)
    while len(L) > 1 and L[-1].is_zero():
        L.pop()
    return L


def is_zero_op(L):
    return all(c.is_zero() for c in L)


def theta_left(B):
    """theta o B"""
    out = [zero_like(B[0]) for _ in range(len(B) + 1)]
    for j, b in enumerate(B):
        if b.is_zero():
            continue
        out[j] = out[j] + Dth(b)
        out[j + 1] = out[j + 1] + b
    return out


def add_op(A, B):
    n = max(len(A), len(B))
    z = zero_like(A[0])
    return [(A[i] if i < len(A) else z) + (B[i] if i < len(B) else z) for i in range(n)]


def scale_op(f, B):
    """f(x) * B (left multiplication by a function)"""
    return [f * b for b in B]


def compose(A, B):
    """A o B"""
    res = [zero_like(B[0])]
    cur = list(B)
    for k, a in enumerate(A):
        if not a.is_zero():
            res = add_op(res, scale_op(a, cur))
        if k < len(A) - 1:
            cur = theta_left(cur)
    return trim(res)


def adjoint(L):
    """L* = sum_k (-theta)^k o c_k"""
    res = [zero_like(L[0])]
    for k, c in enumerate(L):
        cur = [c]
        for _ in range(k):
            cur = theta_left(cur)
        sgn = -1 if k % 2 else 1
        res = add_op(res, [sgn * x for x in cur])
    return trim(res)


def shift_theta(L, s):
    """L o t^s = t^s * L(theta + s): returns the operator L(theta + s) (function part t^s factored out)"""
    r = len(L) - 1
    out = [zero_like(L[0]) for _ in range(r + 1)]
    for k, c in enumerate(L):
        for i in range(k + 1):
            out[i] = out[i] + c * (comb(k, i) * (s ** (k - i)))
    return trim(out)


def pseudo_rdiv(L, R):
    """r^m * L = Qt o R + rem with r = lead(R); returns (m, Qt, rem)."""
    L = trim(L)
    R = trim(R)
    r = len(R) - 1
    lead = R[-1]
    Qt = [zero_like(L[0])]
    m = 0
    while len(L) - 1 >= r and not is_zero_op(L):
        k = len(L) - 1
        lk = L[-1]
        L = scale_op(lead, L)
        Qt = scale_op(lead, Qt)
        m += 1
        term = [zero_like(lk)] * (k - r) + [lk]
        sub = compose(term, R)
        L = trim(add_op(L, [-x for x in sub]))
        Qt = add_op(Qt, term)
        if len(L) - 1 == k and not L[-1].is_zero():
            raise RuntimeError("leading term did not cancel")
    return m, trim(Qt), trim(L)


def content_free(L):
    """divide an fmpz_poly operator by the gcd of all its coefficient polynomials (incl. integer content)"""
    g = None
    for c in L:
        if c.is_zero():
            continue
        g = c if g is None else g.gcd(c)
    if g is None:
        return L
    if g.leading_coefficient() < 0:
        g = -g
    return [flint.fmpz_poly(c // g) if not c.is_zero() else c for c in L]


def op_to_fmpz(L):
    """fmpq_poly operator -> primitive fmpz_poly operator (scalar multiple)"""
    den = 1
    for c in L:
        den = sp.ilcm(den, int(c.denom()))
    out = [flint.fmpz_poly([int(v) for v in (c * den).numer().coeffs()]) if not c.is_zero() else flint.fmpz_poly([]) for c in L]
    return content_free(out)


def op_mod(L, p):
    return [flint.nmod_poly([int(v) % p for v in c.coeffs()], p) for c in L]


def op_degree(L):
    return max(c.degree() for c in L)


# ----------------------------------------------------------------------------------------------
# conversions

def from_Q(Q):
    """QG Q-form (list of sympy polys in u, Q_i <-> t^(2i)) -> fmpz_poly operator"""
    r = max(sp.Poly(q, u_sym).degree() for q in Q)
    c = [[0] * len(Q) for _ in range(r + 1)]
    for i, q in enumerate(Q):
        for (j,), v in sp.Poly(q, u_sym).terms():
            c[j][i] = int(v)
    return [flint.fmpz_poly(row) for row in c]


def to_Q_lists(L):
    """operator -> {i: [coeff of u^0..u^r]} as python ints/str, for json"""
    r = len(L) - 1
    d = op_degree(L)
    out = []
    for i in range(d + 1):
        out.append([str(int(L[k].coeffs()[i])) if i < len(L[k].coeffs()) else "0" for k in range(r + 1)])
    return out


def op_to_strings(L):
    return [str(c) for c in L]


def load_L9():
    rec = json.loads((DATA / "ring_ext_QG_5c_ylog_g1.5_r9_d25.json").read_text())
    Q = [sp.sympify(q, locals={"u": u_sym}) for q in rec["Q"]]
    return from_Q(Q)


def load_G4():
    d = json.loads((DATA / "ring_ext_QG_6_gcrd_g1.json").read_text())
    Q = [sp.sympify(q, locals={"u": u_sym}) for q in d["QG4"]]
    return from_Q(Q)


def L3_op():
    """F_3's operator 4 th^3 - t^2 (th+1)^3"""
    return from_Q([4 * u_sym**3, -(u_sym + 1) ** 3])


def save_op(path, ops: dict, extra=None):
    rec = {k: [str(c) for c in v] for k, v in ops.items()}
    if extra:
        rec.update(extra)
    Path(path).write_text(json.dumps(rec, indent=1))


def parse_op(strs):
    xs = sp.Symbol("x")
    out = []
    for s in strs:
        e = sp.sympify(s.replace("^", "**"), locals={"x": xs})
        if e == 0:
            out.append(flint.fmpz_poly([]))
            continue
        pc = sp.Poly(e, xs).all_coeffs()[::-1]
        out.append(flint.fmpz_poly([int(v) for v in pc]))
    return out


# ----------------------------------------------------------------------------------------------
# local forms

def stirling2(n, k):
    return int(sp.functions.combinatorial.numbers.stirling(n, k, kind=2))


def dform(L):
    """a_j(x) with L = sum_j a_j(x) d^j/dx^j  (theta_t^k = 2^k sum_j S(k,j) x^j d^j)"""
    r = len(L) - 1
    z = zero_like(L[0])
    X = X_of(L[0])
    a = []
    for j in range(r + 1):
        acc = z
        for k in range(j, r + 1):
            if not L[k].is_zero():
                acc = acc + L[k] * ((2 ** k) * stirling2(k, j))
        a.append(acc * X ** j if not acc.is_zero() else acc)
    return a


def valuation(p):
    cs = p.coeffs()
    for i, c in enumerate(cs):
        if c != 0:
            return i
    return 10 ** 9


def hform(L, x0):
    """operator at x0 in h = x - x0, normalised: h^m L = sum_s h^s B_s(theta_h).
    Returns list of coefficient lists of B_s (in u = theta_h), entries in the base ring (int mod p or fmpq).
    x0: Fraction (exact, uses fmpq) or int mod p (for nmod)."""
    a = dform(L)
    mod = isinstance(L[0], flint.nmod_poly)
    if mod:
        p = L[0].modulus()
        shift = flint.nmod_poly([x0 % p, 1], p)
        ah = [ak(shift) if not ak.is_zero() else ak for ak in a]
    else:
        x0 = Fraction(x0)
        shift = flint.fmpq_poly([flint.fmpq(x0.numerator, x0.denominator), 1])
        ah = [flint.fmpq_poly(ak)(shift) if not ak.is_zero() else flint.fmpq_poly([]) for ak in a]
    r = len(L) - 1
    ords = [valuation(ah[j]) for j in range(r + 1)]
    m = max(j - ords[j] for j in range(r + 1) if ords[j] < 10 ** 9)
    S = max(ah[j].degree() - j + m for j in range(r + 1) if not ah[j].is_zero())
    # falling factorial coefficients ff_j(u) = u(u-1)...(u-j+1)
    ff = [[1]]
    for j in range(1, r + 1):
        prev = ff[-1]
        nxt = [0] * (j + 1)
        for i, cf in enumerate(prev):
            nxt[i + 1] += cf
            nxt[i] -= (j - 1) * cf
        ff.append(nxt)
    B = []
    for s in range(S + 1):
        row = [0] * (r + 1)
        for j in range(r + 1):
            idx = s - m + j
            cs = ah[j].coeffs()
            if 0 <= idx < len(cs) and cs[idx] != 0:
                for i, f in enumerate(ff[j]):
                    row[i] = row[i] + cs[idx] * f
        B.append(row)
    return B, m


def indicial_at(L, x0):
    """exact indicial polynomial at finite x0 in theta_h (h = x - x0), as a sympy poly in u.
    At x0 = 0 the local exponent in t is 2 * (root)."""
    B, _ = hform(L, Fraction(x0))
    return sp.Poly(sum(sp.Rational(str(c)) * u_sym ** i for i, c in enumerate(B[0])), u_sym)


def indicial_t0(L):
    """indicial polynomial at t = 0 in u = theta_t (roots = t-exponents)"""
    acc = 0
    for k, c in enumerate(L):
        cs = c.coeffs()
        v = valuation(c)
        acc = acc  # placeholder
    vmin = min(valuation(c) for c in L if not c.is_zero())
    acc = 0
    for k, c in enumerate(L):
        cs = c.coeffs()
        if vmin < len(cs):
            acc += sp.Rational(str(cs[vmin])) * u_sym ** k
    return sp.Poly(acc, u_sym)


def indicial_inf(L):
    """growth exponents mu (y ~ t^mu at t = infinity): top-x-degree part, in u = theta_t"""
    D = op_degree(L)
    acc = 0
    for k, c in enumerate(L):
        if c.degree() == D:
            acc += sp.Rational(str(c.coeffs()[D])) * u_sym ** k
    return sp.Poly(acc, u_sym)


def exps_of(poly):
    """roots with multiplicity of a sympy poly (rational roots + irreducible remainder note)"""
    f = sp.factor_list(poly.as_expr())
    out = []
    rest = []
    for fac, m in f[1]:
        pf = sp.Poly(fac, u_sym)
        if pf.degree() == 1:
            out += [sp.Rational(-pf.all_coeffs()[1], pf.all_coeffs()[0])] * m
        else:
            rest.append((str(fac), m))
    return sorted(out), rest


# ----------------------------------------------------------------------------------------------
# mod p series

def inv(a, p):
    return pow(a % p, p - 2, p)


def frac_mod(e, p):
    e = Fraction(e)
    return e.numerator % p * inv(e.denominator, p) % p


def frob_mod(L, x0, e, N, p, allow_resonance_nolog=True):
    """Frobenius series h^e sum_{n<N} c_n h^n (c_0 = 1) of L at x = x0 (mod p), e rational."""
    Lp = op_mod(L, p) if not isinstance(L[0], flint.nmod_poly) else L
    B, _ = hform(Lp, frac_mod(x0, p))
    ep = frac_mod(e, p)
    Bc = [[int(c) % p for c in row] for row in B]

    def ev(s, x):
        v = 0
        for c in reversed(Bc[s]):
            v = (v * x + c) % p
        return v
    assert ev(0, ep) == 0, "e is not an exponent"
    c = [1]
    S = len(Bc) - 1
    for n in range(1, N):
        acc = 0
        for s in range(1, min(n, S) + 1):
            if c[n - s]:
                acc = (acc + ev(s, (ep + n - s) % p) * c[n - s]) % p
        d = ev(0, (ep + n) % p)
        if d == 0:
            if acc == 0 and allow_resonance_nolog:
                c.append(0)
                continue
            raise ZeroDivisionError(f"log resonance at n = {n}")
        c.append((-acc) * inv(d, p) % p)
    return c


def taylor_basis_mod(L, x0, N, p):
    """basis of solutions at a regular point x0 (exponents 0..r-1), mod p"""
    r = len(L) - 1
    return [frob_mod(L, x0, j, N, p) for j in range(r)]


def theta_series(s, e, x0, p):
    """theta_t on h^e sum s_n h^n at x0: returns series at base e - 1 (same length)"""
    N = len(s)
    out = [0] * N
    x0p = frac_mod(x0, p)
    for n in range(N):
        v = 2 * x0p * ((e + n) % p) % p * s[n] % p
        if n >= 1:
            v = (v + 2 * ((e + n - 1) % p) * s[n - 1]) % p
        out[n] = v
    return out


def columns(series_list, x0, p, r, d):
    """for each series (e, s): columns of x^i theta^k y, i <= d, k <= r, aligned at base e - r.
    Returns rows (list of row lists) for the stacked homogeneous system, and the column labels."""
    x0p = frac_mod(x0, p)
    allrows = []
    labels = [(i, k) for i in range(d + 1) for k in range(r + 1)]
    for e, s in series_list:
        ep = frac_mod(e, p)
        N = len(s)
        ths = [s]
        for k in range(1, r + 1):
            ths.append(theta_series(ths[-1], (ep - (k - 1)) % p, x0, p))
        # (x0 + h)^i as coefficient lists
        xpow = [[1]]
        for i in range(1, d + 1):
            prev = xpow[-1]
            nxt = [0] * (len(prev) + 1)
            for j, cf in enumerate(prev):
                nxt[j] = (nxt[j] + cf * x0p) % p
                nxt[j + 1] = (nxt[j + 1] + cf) % p
            xpow.append(nxt)
        cols = []
        for (i, k) in labels:
            base = [0] * (r - k) + ths[k][: N - (r - k)]
            pol = flint.nmod_poly(base, p) * flint.nmod_poly(xpow[i], p)
            cs = [int(v) for v in pol.coeffs()][:N]
            cs += [0] * (N - len(cs))
            cols.append(cs)
        valid = N - r - 2
        for n in range(valid):
            allrows.append([col[n] for col in cols])
    return allrows, labels


def nullspace_mod(rows, ncols, p):
    M = flint.nmod_mat(len(rows), ncols, [v for row in rows for v in row], p)
    X, nul = M.nullspace()
    vecs = []
    for j in range(nul):
        vecs.append([int(X[i, j]) for i in range(ncols)])
    return vecs


def min_annihilator_mod(series_list, x0, p, r, dmax, dmin=0, margin=30):
    """smallest x-degree d <= dmax with an order-<=r annihilator (exact order r unless lower exists)
    of all series; returns (d, vecs, labels) or None"""
    for d in range(dmin, dmax + 1):
        U = (r + 1) * (d + 1)
        totrows = sum(len(s) - r - 2 for _, s in series_list)
        if totrows < U + margin:
            return ("rows", d)
        rows, labels = columns(series_list, x0, p, r, d)
        vecs = nullspace_mod(rows, U, p)
        if vecs:
            return d, vecs, labels
    return None


def vec_to_op_mod(v, labels, r, d, p):
    c = [[0] * (d + 1) for _ in range(r + 1)]
    for val, (i, k) in zip(v, labels):
        c[k][i] = val
    return [flint.nmod_poly(row, p) for row in c]


# ----------------------------------------------------------------------------------------------
# CRT + rational reconstruction

def ratrec(a, m):
    """rational reconstruction of a mod m: r/s with |r|, s < sqrt(m/2); None if fails"""
    a %= m
    bound = int(flint.fmpz(m // 2).isqrt())
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    if s1 < 0:
        r1, s1 = -r1, -s1
    return Fraction(r1, s1)


def primes_from(start=(1 << 62), count=10):
    out = []
    q = flint.fmpz(start)
    while len(out) < count:
        q = q - 1
        while not q.is_prime():
            q -= 1
        out.append(int(q))
    return out


def crt_vectors(vecs_by_prime):
    """vecs_by_prime: list of (p, vec normalised). Returns (modulus, residues)."""
    m = 1
    res = None
    for p, v in vecs_by_prime:
        if res is None:
            res = [x % p for x in v]
            m = p
            continue
        mi = inv(m % p, p)
        res = [r + m * (((x - r) % p) * mi % p) for r, x in zip(res, v)]
        m *= p
    return m, res


def reconstruct(m, res):
    out = []
    for a in res:
        q = ratrec(a, m)
        if q is None:
            return None
        out.append(q)
    return out


def fracs_to_op(fr, labels, r, d):
    den = 1
    for q in fr:
        den = sp.ilcm(den, q.denominator)
    c = [[0] * (d + 1) for _ in range(r + 1)]
    for q, (i, k) in zip(fr, labels):
        c[k][i] = int(q * den)
    return content_free([flint.fmpz_poly(row) for row in c])


def normalise(v, idx, p):
    iv = inv(v[idx], p)
    return [x * iv % p for x in v]
