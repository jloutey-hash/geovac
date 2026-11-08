"""IB record (docs/predictions/2026-10-24_ibp_kite.md): a small Laporta IBP reducer, exact modulo a prime.

A family: loop momenta L, external momenta E (here one, w, with w.w = z), propagators
D_e = p_e^2 + m_e^2 with p_e an integer combination of L and E, and base indices nu_e (Fraction; the
density lines carry 3/2 + Z).  Integrals I(n) = int prod_e D_e^-(nu_e + n_e), n an integer shift vector.
"Integer" lines (base index an integer) can drop out (index <= 0: a numerator); half-integer lines never do.

Everything is evaluated at one numeric point (z, D, masses) modulo P, so a reduction is a sparse linear
solve over GF(P).  Rational functions of z are rebuilt from many such points elsewhere.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product

import os
P = int(os.environ.get("IB_PRIME", str((1 << 61) - 1)))


def inv(a: int) -> int:
    return pow(a % P, P - 2, P)


def fr(x) -> int:
    x = Fraction(x)
    return x.numerator % P * inv(x.denominator) % P


class Family:
    """loops, externals: names; props: list of (coeff dict over names, mass^2 Fraction, base index Fraction);
    ext_sq: {(e1, e2): value} symbolic tag 'z' handled by the caller giving z numerically."""

    def __init__(self, loops, externals, props, zero_sector=None):
        self.L = list(loops)
        self.E = list(externals)
        self.props = props
        self.ne = len(props)
        names = self.L + self.E
        # scalar products involving at least one loop momentum
        self.sp = [(a, b) for i, a in enumerate(names) for b in names[i:] if a in self.L or b in self.L]
        assert len(self.sp) == self.ne, "propagators must span the scalar products"
        self.integer = [Fraction(nu).denominator == 1 for _, _, nu in props]
        self.zero_sector = zero_sector or (lambda pos: False)
        self.nu2 = [int(2 * Fraction(nu)) for _, _, nu in props]
        self._zc = {}
        # A: D_e = sum_s A[e][s] sp_s + (external-squared part) + m_e^2
        self.A = []
        self.ext_part = []   # list of dicts {(e1,e2): coeff}
        for c, m2, nu in props:
            row = [Fraction(0)] * self.ne
            ext = {}
            for i, a in enumerate(names):
                for j, b in enumerate(names):
                    ca, cb = c.get(a, 0), c.get(b, 0)
                    if ca * cb == 0:
                        continue
                    key = (a, b) if i <= j else (b, a)
                    if key in self.sp:
                        row[self.sp.index(key)] += ca * cb
                    else:
                        ext[key] = ext.get(key, 0) + ca * cb
            self.A.append(row)
            self.ext_part.append(ext)
        self.Ainv = _inverse(self.A)

    def setup(self, ext_vals: dict, D, masses2=None):
        """numeric point: ext_vals {(e1,e2): Fraction or int mod-P value}, D, optional override of m^2."""
        self.Dm = fr(D)
        self.ev = {k: (v % P if isinstance(v, int) else fr(v)) for k, v in ext_vals.items()}
        m2 = masses2 or [m for _, m, _ in self.props]
        self.c = []
        for e in range(self.ne):
            val = fr(m2[e]) if not isinstance(m2[e], int) else m2[e] % P
            for k, v in self.ext_part[e].items():
                val = (val + fr(v) * self.ev[k]) % P
            self.c.append(val)
        self.AinvP = [[fr(x) for x in row] for row in self.Ainv]
        self.nu = [fr(nu) for _, _, nu in self.props]
        self.nuF = [Fraction(nu) for _, _, nu in self.props]
        self._vp = {}

    # sp_s = sum_e Ainv[s][e] (D_e - c_e)  ->  {e: coeff} plus a constant (key None)
    def sp_in_D(self, s):
        out = {}
        const = 0
        for e in range(self.ne):
            a = self.AinvP[s][e]
            if a:
                out[e] = a
                const = (const - a * self.c[e]) % P
        out[None] = const
        return out

    def dot(self, x: dict, y: dict):
        """x.y for momenta given as coefficient dicts: {e or None: coeff} in D's plus constant."""
        out = {}
        names = self.L + self.E
        for a, ca in x.items():
            for b, cb in y.items():
                if not ca or not cb:
                    continue
                i, j = names.index(a), names.index(b)
                key = (a, b) if i <= j else (b, a)
                f = ca * cb % P
                if key in self.sp:
                    for e, v in self.sp_in_D(self.sp.index(key)).items():
                        out[e] = (out.get(e, 0) + f * v) % P
                else:
                    out[None] = (out.get(None, 0) + f * self.ev[key]) % P
        return out

    def index(self, n, e):
        return self.nuF[e] + n[e]

    def is_zero(self, n):
        z = self._zc.get(n)
        if z is None:
            pos = tuple((not self.integer[e]) or self.nu2[e] + 2 * n[e] > 0 for e in range(self.ne))
            z = self._zc[n] = self.zero_sector(pos)
        return z

    def deriv_terms(self, n, a, v):
        """int d/d(l_a) . (v * prod D^-nu) as {shift tuple: coeff}.  a: loop name; v: momentum dict."""
        out = {}
        if v == {a: 1}:
            out[n] = self.Dm
        for e, (c, _, _) in enumerate(self.props):
            ca = c.get(a, 0)
            if not ca:
                continue
            nue = (self.nu[e] + n[e]) % P
            if not nue:
                continue
            pref = (-2 * ca * nue) % P
            ck = (tuple(sorted(v.items())), e)
            vp = self._vp.get(ck)
            if vp is None:
                vp = self._vp[ck] = self.dot(v, c)
            raised = list(n)
            raised[e] += 1
            for f, val in vp.items():
                key = list(raised)
                if f is not None:
                    key[f] -= 1
                key = tuple(key)
                out[key] = (out.get(key, 0) + pref * val) % P
        return {k: v for k, v in out.items() if v}

    def ibp(self, n):
        vecs = [{x: 1} for x in self.L + self.E]
        return [self.deriv_terms(n, a, v) for a in self.L for v in vecs]

    def complexity(self, n):
        # doubled indices, integers only
        t = r = s = 0
        for e in range(self.ne):
            nu2 = self.nu2[e] + 2 * n[e]
            base2 = 2 if self.integer[e] else 1
            if nu2 >= base2:
                if self.integer[e]:
                    t += 1
                r += nu2 - base2
            else:
                s += (base2 - nu2) if not self.integer[e] else -nu2
        return (t, r + s, s, r, tuple(n))


def _inverse(M):
    n = len(M)
    A = [list(map(Fraction, row)) + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        piv = A[c][c]
        A[c] = [x / piv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


class Reducer:
    """Laporta: pivot on the most complex integral of each equation."""

    def __init__(self, fam: Family):
        self.f = fam
        self.rows = {}   # pivot key -> {key: coeff} with pivot coeff 1 (pivot not stored)
        self._key = {}

    def key(self, n):
        k = self._key.get(n)
        if k is None:
            k = self.f.complexity(n)
            self._key[n] = k
        return k

    def add(self, eq: dict):
        eq = {n: c for n, c in eq.items() if c and not self.f.is_zero(n)}
        while eq:
            top = max(eq, key=self.key)
            row = self.rows.get(top)
            if row is None:
                c = inv(eq.pop(top))
                self.rows[top] = {n: v * c % P for n, v in eq.items()}
                return True
            f = eq.pop(top)
            for n, v in row.items():
                nv = (eq.get(n, 0) - f * v) % P
                if nv:
                    eq[n] = nv
                else:
                    eq.pop(n, None)
        return False

    def seed(self, shifts, sym=None):
        for n in sorted(shifts, key=self.key):
            for eq in self.f.ibp(n):
                self.add(dict(eq))
            if sym is not None:
                m = sym(n)
                if m != n:
                    self.add({n: 1, m: P - 1})

    def reduce(self, eq: dict, memo=None):
        """express sum c_n I(n) in non-pivot integrals (masters)."""
        memo = {} if memo is None else memo
        out = {}
        stack = [(n, c) for n, c in eq.items() if c and not self.f.is_zero(n)]
        work = {}
        for n, c in stack:
            work[n] = (work.get(n, 0) + c) % P
        while work:
            top = max(work, key=self.key)
            c = work.pop(top)
            if not c:
                continue
            row = self.rows.get(top)
            if row is None:
                out[top] = (out.get(top, 0) + c) % P
                continue
            for n, v in row.items():
                if self.f.is_zero(n):
                    continue
                work[n] = (work.get(n, 0) - c * v) % P
        return {n: c for n, c in out.items() if c}


def shifts_box(ne, dots, nums, integer, sectors=None):
    """shift vectors relative to the base indices (1 on integer lines, 1/2 on half-integer lines) with
    total dots <= dots and total numerators <= nums, counted per sector: an integer line at index 0 is
    absent, not a numerator; index -k is k numerators; a half-integer line below 1/2 counts its depth."""
    out = []
    rng = range(-nums - 1, dots + 1)
    for n in product(rng, repeat=ne):
        r = s = 0
        for e in range(ne):
            if integer[e]:
                idx = 1 + n[e]
                if idx >= 1:
                    r += idx - 1
                else:
                    s += -idx
            else:
                if n[e] >= 0:
                    r += n[e]
                else:
                    s += -n[e]
        if r <= dots and s <= nums:
            out.append(n)
    return out


def ratrec(a: int, bound: int = 1 << 30):
    """Wang rational reconstruction of a mod P."""
    r0, r1, s0, s1 = P, a % P, 0, 1
    while r1 >= bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) >= bound:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def dz_terms(fam: Family, n, ext: str, z: int):
    """d/dz I(n), z = ext.ext: (1/(2z)) ext.d/d(ext) on the integrand."""
    out = {}
    hz = inv(2 * z)
    for e, (c, _, _) in enumerate(fam.props):
        de = c.get(ext, 0)
        if not de:
            continue
        nue = (fam.nu[e] + n[e]) % P
        pref = (-2 * de * nue) % P * hz % P
        vp = fam.dot({ext: 1}, c)
        raised = list(n)
        raised[e] += 1
        for f, val in vp.items():
            key = list(raised)
            if f is not None:
                key[f] -= 1
            key = tuple(key)
            out[key] = (out.get(key, 0) + pref * val) % P
    return {k: v for k, v in out.items() if v}
