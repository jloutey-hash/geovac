"""Paper 33 Sec. 2, route 2: coordinates on S^3 in R^4 (independent of route 1).

Durable port (2026-11-07) of the VV stage-A route-2 driver, plus the vertex-
tensor assembly (orthogonal bases, sum_m |V|^2, diagonal vertex) that the
stage-A driver did inline.  Record: docs/predictions/2026-11-07_physical_vector_vertex.md.
Imports nothing from route 1 or from geovac.

* Spinors: positive-chirality half S^+ = C^2 of R^4 restricted to S^3 (ambient
  trivialization by constant spinors).  gamma_a = [[0, tau_a],[tau_a^dag, 0]],
  tau = (1, i s1, i s2, i s3).  On S^+, gamma_a gamma_b = tau_a tau_b^dag.
  Clifford multiplication on S^3: c(v) = xslash vslash = sum x_a v_b tau_a tau_b^dag.
  Dirac: D = sum_{a<b} tau_a tau_b^dag L_ab + c0,  L_ab = x_a d_b - x_b d_a, with the
  constant c0 = +-3/2 (Gauss formula); the sign is fixed by requiring the
  spectrum to be symmetric (control C1 tests both).
* 1-forms: tangent polynomial 4-vectors A_a(x), homogeneous degree q.  (*dA)_a =
  eps_abcd x_b d_c A_d (pullback-independent); delta A = -div A.  Eigenforms of *d
  with eigenvalue +-(q+1), checked tangent and coclosed.
* SU(2) x SU(2) content: total so(4) generators G_ab = L_ab + (spin or vector
  part); Casimirs of the two su(2) factors measured on each eigenspace.
* Integrals: exact monomial moments on S^3, normalized volume 1:
  <x^(2b)> = prod (2b_i)!/b_i! / (4^|b| (|b|+1)!).
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement, permutations
from math import factorial
from typing import Dict, List, Tuple

from sympy import QQ, QQ_I, I as sI, Rational
from sympy.polys.matrices import DomainMatrix

NV = 4

# ---------------------------------------------------------------------------
# monomials and polynomials (dict exp-tuple -> Gaussian rational (re, im))
# ---------------------------------------------------------------------------
_mon_cache: Dict[int, List[Tuple[int, ...]]] = {}


def monomials(d: int) -> List[Tuple[int, ...]]:
    if d not in _mon_cache:
        out = []
        for c in combinations_with_replacement(range(NV), d):
            e = [0] * NV
            for i in c:
                e[i] += 1
            out.append(tuple(e))
        _mon_cache[d] = sorted(set(out))
    return _mon_cache[d]


def mon_index(d: int) -> Dict[Tuple[int, ...], int]:
    return {m: i for i, m in enumerate(monomials(d))}


# Gaussian rational as complex-like pair of Fractions
class G:
    __slots__ = ("re", "im")

    def __init__(self, re=Fraction(0), im=Fraction(0)):
        self.re = Fraction(re)
        self.im = Fraction(im)

    def __add__(self, o):
        return G(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return G(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        if isinstance(o, G):
            return G(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
        return G(self.re * o, self.im * o)

    def conj(self):
        return G(self.re, -self.im)

    def iszero(self):
        return self.re == 0 and self.im == 0

    def abs2(self):
        return self.re * self.re + self.im * self.im

    def __repr__(self):
        return f"({self.re}+{self.im}i)"


ONE, ZERO, IU = G(1), G(0), G(0, 1)

Poly = Dict[Tuple[int, ...], G]


def p_add(a: Poly, b: Poly, s: G = ONE) -> Poly:
    out = dict(a)
    for k, v in b.items():
        w = out.get(k, ZERO) + v * s
        if w.iszero():
            out.pop(k, None)
        else:
            out[k] = w
    return out


def p_scale(a: Poly, s) -> Poly:
    out = {}
    for k, v in a.items():
        w = v * s
        if not w.iszero():
            out[k] = w
    return out


def p_mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            w = out.get(k, ZERO) + va * vb
            if w.iszero():
                out.pop(k, None)
            else:
                out[k] = w
    return out


def p_conj(a: Poly) -> Poly:
    return {k: v.conj() for k, v in a.items()}


def p_x(a: Poly, i: int) -> Poly:
    out = {}
    for k, v in a.items():
        e = list(k)
        e[i] += 1
        out[tuple(e)] = v
    return out


def p_d(a: Poly, i: int) -> Poly:
    out: Poly = {}
    for k, v in a.items():
        if k[i] == 0:
            continue
        e = list(k)
        c = e[i]
        e[i] -= 1
        out[tuple(e)] = v * c
    return out


def p_L(a: Poly, i: int, j: int) -> Poly:
    """L_ij = x_i d_j - x_j d_i."""
    return p_add(p_x(p_d(a, j), i), p_x(p_d(a, i), j), G(-1))


# ---------------------------------------------------------------------------
# sphere moments (normalized volume 1)
# ---------------------------------------------------------------------------
def moment(e: Tuple[int, ...]) -> Fraction:
    if any(x % 2 for x in e):
        return Fraction(0)
    b = [x // 2 for x in e]
    num = 1
    for bi in b:
        num *= factorial(2 * bi) // factorial(bi)
    B = sum(b)
    return Fraction(num, 4 ** B * factorial(B + 1))


def integrate(p: Poly) -> G:
    acc = G()
    for k, v in p.items():
        m = moment(k)
        if m:
            acc = acc + v * m
    return acc


# ---------------------------------------------------------------------------
# Pauli / tau matrices (2x2 over Z[i])
# ---------------------------------------------------------------------------
def mat(a, b, c, d):
    return [[a, b], [c, d]]


S1 = mat(G(0), G(1), G(1), G(0))
S2 = mat(G(0), G(0, -1), G(0, 1), G(0))
S3 = mat(G(1), G(0), G(0), G(-1))
ID = mat(G(1), G(0), G(0), G(1))


def mmul(A, B):
    return [[A[i][0] * B[0][j] + A[i][1] * B[1][j] for j in range(2)] for i in range(2)]


def mdag(A):
    return [[A[j][i].conj() for j in range(2)] for i in range(2)]


def mscale(A, s):
    return [[A[i][j] * s for j in range(2)] for i in range(2)]


TAU = [ID, mscale(S1, IU), mscale(S2, IU), mscale(S3, IU)]
TT = [[mmul(TAU[a], mdag(TAU[b])) for b in range(4)] for a in range(4)]  # tau_a tau_b^dag

# Levi-Civita on 4 indices
EPS4: Dict[Tuple[int, int, int, int], int] = {}
for perm in permutations(range(4)):
    sgn = 1
    p = list(perm)
    for i in range(4):
        for j in range(i + 1, 4):
            if p[i] > p[j]:
                sgn = -sgn
    EPS4[perm] = sgn

# ---------------------------------------------------------------------------
# spinor fields: list of 2 Polys ; vector fields: list of 4 Polys
# ---------------------------------------------------------------------------
Spinor = List[Poly]
Vector = List[Poly]


def mat_apply(M, psi: Spinor) -> Spinor:
    return [p_add(p_scale(psi[0], M[i][0]), p_scale(psi[1], M[i][1])) for i in range(2)]


def D0(psi: Spinor) -> Spinor:
    """sum_{a<b} tau_a tau_b^dag L_ab psi (no constant)."""
    out: Spinor = [{}, {}]
    for a in range(4):
        for b in range(a + 1, 4):
            Lpsi = [p_L(psi[0], a, b), p_L(psi[1], a, b)]
            t = mat_apply(TT[a][b], Lpsi)
            out = [p_add(out[0], t[0]), p_add(out[1], t[1])]
    return out


def G_spinor(psi: Spinor, a: int, b: int) -> Spinor:
    """so(4) generator on spinor fields: L_ab + (1/2) tau_a tau_b^dag."""
    Lpsi = [p_L(psi[0], a, b), p_L(psi[1], a, b)]
    t = mat_apply(TT[a][b], psi)
    return [p_add(Lpsi[i], t[i], G(Fraction(1, 2))) for i in range(2)]


def star_d(A: Vector) -> Vector:
    out: Vector = [{}, {}, {}, {}]
    for (a, b, c, d), s in EPS4.items():
        t = p_x(p_d(A[d], c), b)
        out[a] = p_add(out[a], t, G(s))
    return out


def div(A: Vector) -> Poly:
    out: Poly = {}
    for a in range(4):
        out = p_add(out, p_d(A[a], a))
    return out


def x_dot(A: Vector) -> Poly:
    out: Poly = {}
    for a in range(4):
        out = p_add(out, p_x(A[a], a))
    return out


def G_vector(A: Vector, a: int, b: int) -> Vector:
    """so(4) generator on vector fields: L_ab + E_ab, (E_ab A)_c = d_ca A_b - d_cb A_a."""
    out = [p_L(A[c], a, b) for c in range(4)]
    out[a] = p_add(out[a], A[b])
    out[b] = p_add(out[b], A[a], G(-1))
    return out


def G_scalar(f: Poly, a: int, b: int) -> Poly:
    return p_L(f, a, b)


# ---------------------------------------------------------------------------
# exact linear algebra helpers
# ---------------------------------------------------------------------------
def g_to_sym(g: G):
    return Rational(g.re.numerator, g.re.denominator) + sI * Rational(g.im.numerator, g.im.denominator)


def sym_to_g(x) -> G:
    from sympy import re, im, nsimplify
    r = Rational(re(x))
    i = Rational(im(x))
    return G(Fraction(int(r.p), int(r.q)), Fraction(int(i.p), int(i.q)))


def flatten(fields: List[List[Poly]], d: int) -> List[List[G]]:
    """columns: each field -> vector over (component, monomial of degree d)."""
    idx = mon_index(d)
    nm = len(idx)
    cols = []
    for f in fields:
        v = [ZERO] * (len(f) * nm)
        for c, p in enumerate(f):
            for k, val in p.items():
                assert sum(k) == d, (k, d)
                v[c * nm + idx[k]] = val
        cols.append(v)
    return cols


def nullspace_QQI(rows: List[List[G]]) -> List[List[G]]:
    dom = QQ_I
    data = [[dom.from_sympy(g_to_sym(x)) for x in r] for r in rows]
    M = DomainMatrix(data, (len(rows), len(rows[0])), dom)
    N = M.nullspace()  # rows are basis vectors
    out = []
    for r in N.to_Matrix().tolist():
        out.append([sym_to_g(x) for x in r])
    return out


def columns_to_rows(cols: List[List[G]]) -> List[List[G]]:
    return [list(r) for r in zip(*cols)]


# ---------------------------------------------------------------------------
# harmonic polynomials
# ---------------------------------------------------------------------------
def harmonic_basis(k: int) -> List[Poly]:
    mons = monomials(k)
    if k < 2:
        return [{m: ONE} for m in mons]
    lap_cols = []
    for m in mons:
        p = {m: ONE}
        lp: Poly = {}
        for a in range(4):
            lp = p_add(lp, p_d(p_d(p, a), a))
        lap_cols.append(lp)
    cols = flatten([[c] for c in lap_cols], k - 2)
    rows = columns_to_rows(cols)
    ns = nullspace_QQI(rows)
    basis = []
    for vec in ns:
        basis.append({mons[i]: vec[i] for i in range(len(mons)) if not vec[i].iszero()})
    return basis


# ---------------------------------------------------------------------------
# inner products, Gram-Schmidt
# ---------------------------------------------------------------------------
def inner(f: List[Poly], g: List[Poly]) -> G:
    acc = G()
    for a, b in zip(f, g):
        acc = acc + integrate(p_mul(p_conj(a), b))
    return acc


def field_add(f, g, s: G):
    return [p_add(a, b, s) for a, b in zip(f, g)]


def gram_schmidt(fields: List[List[Poly]]) -> Tuple[List[List[Poly]], List[Fraction]]:
    """Orthogonal (not normalized) basis and squared norms (exact)."""
    out, norms = [], []
    for f in fields:
        g = f
        for h, nh in zip(out, norms):
            c = inner(h, g)
            if not c.iszero():
                g = field_add(g, h, c * Fraction(-1) * (Fraction(1) / nh))
        n2 = inner(g, g)
        assert n2.im == 0 and n2.re > 0
        out.append(g)
        norms.append(n2.re)
    return out, norms


# ---------------------------------------------------------------------------
# eigenspaces
# ---------------------------------------------------------------------------
def spinor_eigenspaces(k: int, c0: Fraction):
    """Eigenspaces of D = D0 + c0 on S^+ (x) H_k. Returns dict eigenvalue -> list of spinors,
    plus total dimension of S^+ (x) H_k."""
    H = harmonic_basis(k)
    basis: List[Spinor] = []
    for h in H:
        basis.append([h, {}])
        basis.append([{}, h])
    Dcols = flatten([D0(b) for b in basis], k)
    Ecols = flatten(basis, k)
    out = {}
    # candidate eigenvalues of D0 on S^+ (x) H_k: -k and k+2 (found, then the
    # dimension count in C1 certifies there are no others); D = D0 + c0.
    for lam0 in (Fraction(-k), Fraction(k + 2)):
        lam = lam0 + c0
        shift = c0 - lam
        cols = [[dc[i] + ec[i] * shift for i in range(len(dc))] for dc, ec in zip(Dcols, Ecols)]
        ns = nullspace_QQI(columns_to_rows(cols))
        fields = []
        for vec in ns:
            psi: Spinor = [{}, {}]
            for coef, b in zip(vec, basis):
                if not coef.iszero():
                    psi = [p_add(psi[0], b[0], coef), p_add(psi[1], b[1], coef)]
            fields.append(psi)
        out[lam] = fields
    return out, len(basis)


def check_spinor_eigen(psi: Spinor, c0: Fraction, lam: Fraction) -> bool:
    Dpsi = D0(psi)
    for i in range(2):
        r = p_add(p_add(Dpsi[i], psi[i], G(c0)), psi[i], G(-lam))
        if r:
            return False
    return True


def vector_eigenspace(q: int, lam: int) -> List[Vector]:
    mons = monomials(q)
    basis: List[Vector] = []
    for a in range(4):
        for m in mons:
            v: Vector = [{}, {}, {}, {}]
            v[a] = {m: ONE}
            basis.append(v)
    Scols = flatten([star_d(b) for b in basis], q)
    Ecols = flatten(basis, q)
    cols = [[sc[i] - ec[i] * lam for i in range(len(sc))] for sc, ec in zip(Scols, Ecols)]
    ns = nullspace_QQI(columns_to_rows(cols))
    fields = []
    for vec in ns:
        A: Vector = [{}, {}, {}, {}]
        for coef, b in zip(vec, basis):
            if not coef.iszero():
                A = [p_add(A[i], b[i], coef) for i in range(4)]
        fields.append(A)
    return fields


# ---------------------------------------------------------------------------
# su(2) x su(2) Casimirs
# ---------------------------------------------------------------------------
# indices 0..3 = x1..x4 ; J^{+-}_k = 1/2 (G_{0k} +- 1/2 eps_kij G_ij), k = 1..3
EPS3 = {(1, 2, 3): 1, (2, 3, 1): 1, (3, 1, 2): 1, (1, 3, 2): -1, (3, 2, 1): -1, (2, 1, 3): -1}


def casimir(field, Gop, sign: int):
    """Sum_k (J^sign_k)^2 field."""
    def J(f, k):
        t = Gop(f, 0, k)
        out = [dict(p) for p in t] if isinstance(t, list) else dict(t)
        for (kk, i, j), e in EPS3.items():
            if kk != k:
                continue
            g = Gop(f, i, j)
            if isinstance(out, list):
                out = [p_add(o, gg, G(Fraction(sign * e, 2))) for o, gg in zip(out, g)]
            else:
                out = p_add(out, g, G(Fraction(sign * e, 2)))
        if isinstance(out, list):
            return [p_scale(o, Fraction(1, 2)) for o in out]
        return p_scale(out, Fraction(1, 2))

    tot = None
    for k in (1, 2, 3):
        t = J(J(field, k), k)
        if tot is None:
            tot = t
        elif isinstance(tot, list):
            tot = [p_add(a, b) for a, b in zip(tot, t)]
        else:
            tot = p_add(tot, t)
    return tot


def eigen_ratio(f, g):
    """If g = r f (fields, lists of polys or a poly), return r, else None."""
    if not isinstance(f, list):
        f, g = [f], [g]
    r = None
    for a, b in zip(f, g):
        keys = set(a) | set(b)
        for k in keys:
            va = a.get(k, ZERO)
            vb = b.get(k, ZERO)
            if va.iszero():
                if not vb.iszero():
                    return None
                continue
            # r = vb / va
            den = va.abs2()
            rr = vb * va.conj() * (Fraction(1) / den)
            if r is None:
                r = rr
            elif not (r - rr).iszero():
                return None
    return r


# ---------------------------------------------------------------------------
# integer-scaled fields and the vertex tensor
# ---------------------------------------------------------------------------
from math import lcm as _lcm
import numpy as _np


def to_int_field(f: List[Poly]) -> Tuple[List[Dict[Tuple[int, ...], Tuple[int, int]]], int]:
    """Scale a field to Gaussian-integer coefficients; returns (field, scale L)."""
    L = 1
    for p in f:
        for v in p.values():
            L = _lcm(L, v.re.denominator, v.im.denominator)
    out = []
    for p in f:
        out.append({k: (int(v.re * L), int(v.im * L)) for k, v in p.items()})
    return out, L


def ipoly_mul_conj(a, b):
    """conj(a) * b for Gaussian-integer polys {exp: (re, im)}."""
    out: Dict[Tuple[int, ...], List[int]] = {}
    for ka, (ar, ai) in a.items():
        for kb, (br, bi) in b.items():
            k = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2], ka[3] + kb[3])
            # (ar - i ai)(br + i bi)
            re = ar * br + ai * bi
            im = ar * bi - ai * br
            if k in out:
                out[k][0] += re
                out[k][1] += im
            else:
                out[k] = [re, im]
    return out


def _tt_int():
    T = [[None] * 4 for _ in range(4)]
    for a in range(4):
        for b in range(4):
            T[a][b] = [[(int(TT[a][b][s][t].re), int(TT[a][b][s][t].im)) for t in range(2)] for s in range(2)]
    return T


TTI = _tt_int()


def current(u, v, dK: int):
    """K_b = sum_a x_a sum_st conj(u_s) (tau_a tau_b^dag)_st v_t, b = 0..3, as flat
    re/im integer arrays over (b, monomial of degree dK)."""
    idx = mon_index(dK)
    nm = len(idx)
    re = _np.zeros(4 * nm, dtype=object)
    im = _np.zeros(4 * nm, dtype=object)
    re[:] = 0
    im[:] = 0
    P = [[ipoly_mul_conj(u[s], v[t]) for t in range(2)] for s in range(2)]
    for a in range(4):
        for b in range(4):
            for s in range(2):
                for t in range(2):
                    cr, ci = TTI[a][b][s][t]
                    if cr == 0 and ci == 0:
                        continue
                    for k, (pr, pi) in P[s][t].items():
                        e = list(k)
                        e[a] += 1
                        j = b * nm + idx[tuple(e)]
                        re[j] += cr * pr - ci * pi
                        im[j] += cr * pi + ci * pr
    return re, im


def moment_tilde(e) -> int:
    if any(x % 2 for x in e):
        return 0
    num = 1
    for x in e:
        num *= factorial(x) // factorial(x // 2)
    return num


def moment_den(total_deg: int) -> int:
    B = total_deg // 2
    return 4 ** B * factorial(B + 1)


def wm_matrix(wfields_int, q: int, dK: int):
    """Columns j: (Mtilde w_j) over (b, monomial deg dK); Mtilde[m, m'] = moment_tilde(m+m')."""
    mK = monomials(dK)
    nm = len(mK)
    out = _np.zeros((4 * nm, len(wfields_int)), dtype=object)
    out[:] = 0
    for j, w in enumerate(wfields_int):
        for b in range(4):
            for kw, (wr, wi) in w[b].items():
                assert wi == 0
                for i, m in enumerate(mK):
                    mt = moment_tilde((m[0] + kw[0], m[1] + kw[1], m[2] + kw[2], m[3] + kw[3]))
                    if mt:
                        out[b * nm + i, j] += mt * wr
    return out


# ---------------------------------------------------------------------------
# assembly (ported from the stage-A driver; exact throughout)
# ---------------------------------------------------------------------------
C0 = Fraction(-3, 2)   # constant term of D; C1 shows it is the one symmetric choice


def spinor_spaces(nmax: int, c0: Fraction = C0) -> Dict[Tuple[int, int], List[Spinor]]:
    """(sign of D, level n) -> eigenspinors with D = sign (n + 3/2), every n <= nmax.

    Raises AssertionError if the polynomial degrees k <= nmax + 1 do not
    decompose into eigenspaces of the stated form."""
    spaces: Dict[Tuple[int, int], List[Spinor]] = {}
    for k in range(nmax + 2):
        es, dim = spinor_eigenspaces(k, c0)
        assert sum(len(v) for v in es.values()) == dim, ("eigenspaces do not exhaust", k)
        for lam, fields in es.items():
            if not fields:
                continue
            s = 1 if lam > 0 else -1
            n = abs(lam) - Fraction(3, 2)
            assert n.denominator == 1 and n >= 0, ("eigenvalue not +-(n+3/2)", lam)
            assert (s, int(n)) not in spaces
            spaces[(s, int(n))] = fields
    return {k_: v for k_, v in spaces.items() if k_[1] <= nmax}


def vector_spaces(qmax: int) -> Dict[Tuple[int, int], List[Vector]]:
    """(sign of *d, q) -> polynomial eigenforms with *d = sign (q + 1)."""
    return {(sd, q): vector_eigenspace(q, sd * (q + 1)) for q in range(1, qmax + 1) for sd in (1, -1)}


def _calibration() -> Dict[int, Fraction]:
    x1 = {(1, 0, 0, 0): ONE}
    return {sg: eigen_ratio(x1, casimir(x1, G_scalar, sg)).re for sg in (1, -1)}


def content(fields, Gop, cal) -> Tuple[int, int]:
    """(2 j_+, 2 j_-) measured by the two su(2) Casimirs (scalar on every field)."""
    res = []
    for sg in (1, -1):
        vals = set()
        for f in fields:
            r = eigen_ratio(f, casimir(f, Gop, sg))
            assert r is not None and r.im == 0, "Casimir not scalar"
            vals.add(r.re)
        assert len(vals) == 1
        jj = vals.pop() / cal[sg] * Fraction(3, 4)
        t = 0
        while Fraction(t * (t + 2), 4) < jj:
            t += 1
        assert Fraction(t * (t + 2), 4) == jj
        res.append(t)
    return tuple(res)


class Route2:
    """Orthogonal integer-scaled bases on route 2's range and the vertex sums."""

    def __init__(self, nmax: int, qmax: int):
        self.nmax, self.qmax = nmax, qmax
        self.SP = spinor_spaces(nmax)
        self.VEC = vector_spaces(qmax)
        cal = _calibration()
        cont = {("s",) + k_: content(v, G_spinor, cal) for k_, v in self.SP.items()}
        c_plus0 = cont[("s", 1, 0)]
        first = 0 if c_plus0[0] > c_plus0[1] else 1
        fs = lambda c: (c[first], c[1 - first])  # noqa: E731
        self.spinor_content = {k_[1:]: fs(v) for k_, v in cont.items()}
        self.comp_sign: Dict[Tuple[int, str], int] = {}
        for (sd, q), flds in self.VEC.items():
            c = fs(content(flds, G_vector, cal))
            if c == (q + 1, q - 1):
                self.comp_sign[(q, "A")] = sd
            elif c == (q - 1, q + 1):
                self.comp_sign[(q, "B")] = sd
            else:
                raise AssertionError(("vector content", sd, q, c))
        self.SPB = {k_: self._ortho(v) for k_, v in self.SP.items()}
        self.VCB = {qc: self._ortho(self.VEC[(sd, qc[0])]) for qc, sd in self.comp_sign.items()}
        self._wm: Dict = {}

    @staticmethod
    def _ortho(fields):
        ob, nrm = gram_schmidt(fields)
        out = []
        for f, n2 in zip(ob, nrm):
            fi, L = to_int_field(f)
            out.append((fi, n2 * L * L))
        return out

    def deg(self, s: int, n: int) -> int:
        for p in self.SP[(s, n)][0]:
            for k_ in p:
                return sum(k_)
        raise AssertionError

    def wm(self, q: int, comp: str, dK: int):
        ck = (q, comp, dK)
        if ck not in self._wm:
            self._wm[ck] = wm_matrix([w for (w, _) in self.VCB[(q, comp)]], q, dK)
        return self._wm[ck]

    def cell(self, s1: int, n1: int, comp: str, q: int, s2: int, n2: int) -> Tuple[bool, Fraction]:
        """(open, sum over all m of |<psi1| c(A) |psi2>|^2), Haar volume one, orthonormal bases."""
        U, Vv = self.SPB[(s1, n1)], self.SPB[(s2, n2)]
        dK = self.deg(s1, n1) + self.deg(s2, n2) + 1
        WM = self.wm(q, comp, dK)
        W = self.VCB[(q, comp)]
        C = moment_den(dK + q)
        S, nz = Fraction(0), False
        for (u, nu) in U:
            for (v, nv) in Vv:
                re, im = current(u, v, dK)
                vr, vi = re.dot(WM), im.dot(WM)
                part = Fraction(0)
                for jw, (w, nw) in enumerate(W):
                    a2 = vr[jw] * vr[jw] + vi[jw] * vi[jw]
                    if a2:
                        nz = True
                        part += Fraction(a2, nw)
                if part:
                    S += part / (C * C * nu * nv)
        return nz, S

    def diagonal_ground(self, s: int) -> Dict[str, List[Fraction]]:
        """Per n = 0 basis state a: sum over q = 1 photon m of |<a| c(A) |a>|^2, per component."""
        U = self.SPB[(s, 0)]
        d = self.deg(s, 0)
        out: Dict[str, List[Fraction]] = {"A": [], "B": []}
        for (u, nu) in U:
            re, im = current(u, u, 2 * d + 1)
            for comp in ("A", "B"):
                WM = self.wm(1, comp, 2 * d + 1)
                C = moment_den(2 * d + 2)
                vr, vi = re.dot(WM), im.dot(WM)
                tot = Fraction(0)
                for jw, (w, nw) in enumerate(self.VCB[(1, comp)]):
                    tot += Fraction(vr[jw] ** 2 + vi[jw] ** 2, C * C * nu * nu * nw)
                out[comp].append(tot)
        return out


def tangential_gradient(f: Poly, p: int) -> Vector:
    """Tangential gradient of a degree-p harmonic f on S^3, homogenized: r^2 grad f - p x f."""
    out = []
    for a in range(4):
        g = p_d(f, a)
        t: Poly = {}
        for b in range(4):
            t = p_add(t, p_x(p_x(g, b), b))
        out.append(p_add(t, p_x(f, a), G(-p)))
    return out


def clifford(Avec: Vector, psi: Spinor) -> Spinor:
    """c(A) psi = sum_ab x_a A_b tau_a tau_b^dag psi."""
    out: Spinor = [{}, {}]
    for a in range(4):
        for b in range(4):
            xa = [p_x(p_mul(Avec[b], psi[0]), a), p_x(p_mul(Avec[b], psi[1]), a)]
            t = mat_apply(TT[a][b], xa)
            out = [p_add(out[0], t[0]), p_add(out[1], t[1])]
    return out
