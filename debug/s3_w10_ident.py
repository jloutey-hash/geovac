"""Sprint S^(3) W10 identification (2026-06-12).

Closes the v4.7.0 open item: identify the two depth-3 generators
t3_GV(2,4,4), t3_GV(3,3,4) of the W10 block at derivable grade, using

  (1) Charlton-Hoffman Symmetry Theorem 2.21 (arXiv:2204.14183), depth
      m=3 at phi=0: generating-series identity whose coefficient at
      y1^(a-1) y2^(b-1) y3^(c-1) gives the explicit-product form of
      t_CH(a,b,c) + t_CH(c,b,a) = products  (weight 10 even).
  (2) Four clean lam*t2 descending stuffles whose pairwise differences
      eliminate the auxiliary triples t3(4,2,4) / t3(3,4,3) and give the
      pair DIFFERENCES over the catalogued ring.
  (3) Level-1 depth-2 double shuffle at w in {4,6,8} (word engine,
      conventions numerically calibrated) + classical Euler trailing-1
      rows (numerically gated), reducing every zeta double to the ring
      Q[zeta(2),zeta(3),...] + one w8 generator zeta53 := zeta_GV(5,3).

CONVENTIONS (two in play, never mixed silently):
  GV/Hoffman descending: t_GV(s1,...,sk) = sum_{n1>...>nk odd} prod ni^-si
                         zeta_GV(s1,..)  = sum_{n1>...>nk}     prod ni^-si
  CH ascending (the 2204.14183 paper): t_CH(n1,...,nm), zeta_CH likewise,
                         k1 < k2 < ... < km, FIRST index on SMALLEST.
  Conversion: t_CH(a,b,c) = t_GV(c,b,a); zeta_CH(a,b) = zeta_GV(b,a).

No PSLQ anywhere. All coefficients exact Fractions; numerics are gates.

Subcommands: extract | verify | assemble | all
"""
from __future__ import annotations
import json
import sys
import itertools
from fractions import Fraction as F
from collections import defaultdict

from mpmath import mp, mpf, zeta as mpzeta, euler as mpeuler, bernoulli as mpbern

HERE = r"C:\Users\jlout\Desktop\Project_Geometric\debug"
DATA = HERE + r"\data"

TOTDEG = 7          # we extract coefficient of monomials of total degree 7
WEIGHT = 10

# ----------------------------------------------------------------------
# SECTION 1: coefficient algebra.
# A "term" is a tuple of symbol-atoms, sorted; an atom is a tuple like
#   ('t1', n) ('t2', n1, n2) ('t3', n1, n2, n3)   [CH ascending t-values]
#   ('z1', n) ('z2', n1, n2)                      [CH ascending MZVs, T=0 reg]
# A "coef" is dict[term] -> Fraction.  A "poly" is dict[(i,j,k)] -> coef
# over monomials y1^i y2^j y3^k with i+j+k <= TOTDEG.
# ----------------------------------------------------------------------

def cmul(c1: dict, c2: dict) -> dict:
    out = defaultdict(F)
    for t1, f1 in c1.items():
        for t2, f2 in c2.items():
            out[tuple(sorted(t1 + t2))] += f1 * f2
    return {t: f for t, f in out.items() if f}


def cadd_into(acc: dict, c: dict, scale: F = F(1)) -> None:
    for t, f in c.items():
        acc[t] = acc.get(t, F(0)) + f * scale
        if not acc[t]:
            del acc[t]


def pmul(p1: dict, p2: dict) -> dict:
    out = {}
    for m1, c1 in p1.items():
        for m2, c2 in p2.items():
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])
            if sum(m) > TOTDEG:
                continue
            c = cmul(c1, c2)
            if m in out:
                cadd_into(out[m], c)
            else:
                out[m] = dict(c)
    return {m: c for m, c in out.items() if c}


def padd(*ps, scales=None) -> dict:
    out = {}
    scales = scales or [F(1)] * len(ps)
    for p, s in zip(ps, scales):
        for m, c in p.items():
            if m not in out:
                out[m] = {}
            cadd_into(out[m], c, s)
    return {m: c for m, c in out.items() if c}


def const_poly(term, coef=F(1)) -> dict:
    return {(0, 0, 0): {term: coef}}


def var_poly(idx: int, sign: int = 1) -> dict:
    """The polynomial  sign * y_idx  (idx in 0,1,2)."""
    m = [0, 0, 0]
    m[idx] = 1
    return {tuple(m): {(): F(sign)}}


def ppow(p: dict, k: int) -> dict:
    out = {(0, 0, 0): {(): F(1)}}
    for _ in range(k):
        out = pmul(out, p)
    return out


# ----------------------------------------------------------------------
# SECTION 2: generating series of Theorem 2.21 (m = 3, phi = 0), CH conv.
# ----------------------------------------------------------------------

def t1_series(arg: dict) -> dict:
    """Li^t depth 1 in polynomial argument `arg`: sum_n t1[n] arg^(n-1).
    t1[1] = T = log 2 (stuffle reg) -- kept as symbol ('T',)."""
    out = {}
    for n in range(1, WEIGHT + 1):
        term = ('T',) if n == 1 else ('t1', n)
        out = padd(out, pmul(const_poly((term,)), ppow(arg, n - 1)))
    return out


def t2_series(arg1: dict, arg2: dict) -> dict:
    out = {}
    for n1 in range(1, WEIGHT):
        a1 = ppow(arg1, n1 - 1)
        for n2 in range(1, WEIGHT + 1 - n1):
            # CH ascending: divergent iff LAST index 1 -> stuffle reg symbol
            atom = ('t2r', n1, n2) if n2 == 1 else ('t2', n1, n2)
            out = padd(out, pmul(const_poly((atom,)),
                                 pmul(a1, ppow(arg2, n2 - 1))))
    return out


def t3_series(arg1: dict, arg2: dict, arg3: dict) -> dict:
    out = {}
    for n1 in range(1, WEIGHT - 1):
        a1 = ppow(arg1, n1 - 1)
        for n2 in range(1, WEIGHT - n1):
            a12 = pmul(a1, ppow(arg2, n2 - 1))
            for n3 in range(1, WEIGHT + 1 - n1 - n2):
                atom = (('t3r', n1, n2, n3) if n3 == 1
                        else ('t3', n1, n2, n3))
                out = padd(out, pmul(const_poly((atom,)),
                                     pmul(a12, ppow(arg3, n3 - 1))))
    return out


def z1_series(arg: dict) -> dict:
    """MZV depth-1 generating series, T=0 stuffle reg: zeta(1) -> 0."""
    out = {}
    for n in range(2, WEIGHT + 1):
        out = padd(out, pmul(const_poly((('z1', n),)), ppow(arg, n - 1)))
    return out


def z2_series(arg1: dict, arg2: dict) -> dict:
    out = {}
    for n1 in range(1, WEIGHT):
        a1 = ppow(arg1, n1 - 1)
        for n2 in range(1, WEIGHT + 1 - n1):
            atom = ('z2r', n1, n2) if n2 == 1 else ('z2', n1, n2)
            out = padd(out, pmul(const_poly((atom,)),
                                 pmul(a1, ppow(arg2, n2 - 1))))
    return out


def bt_series(idx: int) -> dict:
    """B^t(0|y_idx) = 2 sum_{k>=1} t(2k) y^(2k-1)  (no T dependence)."""
    out = {}
    y = var_poly(idx)
    for k in range(1, WEIGHT // 2 + 1):
        out = padd(out, pmul(const_poly((('t1', 2 * k),), F(2)),
                             ppow(y, 2 * k - 1)))
    return out


def half_diff(i: int, j: int) -> dict:
    """(1/2)(y_i - y_j) as a poly (i, j are 0-based variable indices)."""
    p = padd(var_poly(i), var_poly(j, -1))
    return {m: {t: f / 2 for t, f in c.items()} for m, c in p.items()}


def build_theorem_m3() -> dict:
    """Full LHS of Thm 2.21 at m=3, phi=0 (equals 0 identically)."""
    y1, y2, y3 = (var_poly(i) for i in range(3))
    n1, n2, n3 = (var_poly(i, -1) for i in range(3))
    # line 1: sum_{j=0..3} (-1)^j Li^t(y_{j+1..m}) Li^t(-y_j..-y_1)
    line1 = padd(
        t3_series(y1, y2, y3),                              # j=0
        pmul(t2_series(y2, y3), t1_series(n1)),             # j=1 (sign -)
        pmul(t1_series(y3), t2_series(n2, n1)),             # j=2 (+)
        t3_series(n3, n2, n1),                              # j=3 (-)
        scales=[F(1), F(-1), F(1), F(-1)],
    )
    # line 2: -(1/4) sum_{j=1..3} (-1)^(j-1) Li(..)B(y_j)Li(..)
    l2_1 = pmul(bt_series(0), z2_series(half_diff(1, 0), half_diff(2, 0)))
    l2_2 = pmul(z1_series(half_diff(1, 0)),
                pmul(bt_series(1), z1_series(half_diff(2, 1))))
    l2_3 = pmul(z2_series(half_diff(2, 1), half_diff(2, 0)), bt_series(2))
    line2 = padd(l2_1, l2_2, l2_3, scales=[F(1), F(-1), F(1)])
    return padd(line1, line2, scales=[F(1), F(-1, 4)])


def build_theorem_m2() -> dict:
    """Thm 2.21 at m=2, phi=0; equals (1/2!)(i pi/2)^2 = -pi^2/8 times
    the constant monomial.  We only use coefficients of nonconstant
    monomials, where the RHS is 0."""
    global TOTDEG
    y1, y2 = var_poly(0), var_poly(1)
    m1, m2 = var_poly(0, -1), var_poly(1, -1)
    line1 = padd(
        t2_series(y1, y2),                                  # j=0
        pmul(t1_series(y2), t1_series(m1)),                 # j=1 (-)
        t2_series(m2, m1),                                  # j=2 (+)
        scales=[F(1), F(-1), F(1)],
    )
    l2_1 = pmul(bt_series(0), z1_series(half_diff(1, 0)))
    l2_2 = pmul(z1_series(half_diff(1, 0)), bt_series(1))
    line2 = padd(l2_1, l2_2, scales=[F(1), F(-1)])
    return padd(line1, line2, scales=[F(1), F(-1, 2)])


def extract(p: dict, a: int, b: int, c: int = None) -> dict:
    """Coefficient (as coef-dict) of y1^(a-1) y2^(b-1) [y3^(c-1)]."""
    if c is None:
        return p.get((a - 1, b - 1, 0), {})
    return p.get((a - 1, b - 1, c - 1), {})


# ----------------------------------------------------------------------
# SECTION 3: descending stuffle for the pair-difference relations.
# ----------------------------------------------------------------------

def stuffle_1_into_2(x: int, y: int, z: int) -> dict:
    """t_GV(x) * t_GV(y,z) expansion as dict[GV-atom] -> Fraction.
    Atoms: ('T3', a,b,c), ('T2', a,b)  in GV descending convention."""
    out = defaultdict(F)
    out[('T3', x, y, z)] += 1
    out[('T3', y, x, z)] += 1
    out[('T3', y, z, x)] += 1
    out[('T2', x + y, z)] += 1
    out[('T2', y, x + z)] += 1
    return dict(out)


def pair_difference_relations() -> dict:
    """Return the two clean difference relations (GV convention):
       t3(4,4,2) - t3(2,4,4) = lam2*t2(4,4) - lam4*t2(2,4) - t2(4,6) + t2(2,8)
       t3(3,3,4) - t3(4,3,3) = (1/2)[lam3*(t2(3,4)-t2(4,3))
                                     - t2(6,4) - t2(3,7) + t2(7,3) + t2(4,6)]
    derived by eliminating t3(4,2,4) (resp. t3(3,4,3)) from two stuffles."""
    RA = stuffle_1_into_2(2, 4, 4)   # lam2*t2(4,4)
    RB = stuffle_1_into_2(4, 2, 4)   # lam4*t2(2,4)
    RC = stuffle_1_into_2(3, 3, 4)   # lam3*t2(3,4)
    RD = stuffle_1_into_2(3, 4, 3)   # lam3*t2(4,3)
    return {'RA': RA, 'RB': RB, 'RC': RC, 'RD': RD}


# ----------------------------------------------------------------------
# SECTION 4: level-1 word-shuffle engine, weight-w depth-2 reduction.
# ----------------------------------------------------------------------

def shuffle_words(w1: tuple, w2: tuple) -> dict:
    """Shuffle product of two words; returns dict[word]->int multiplicity."""
    out = defaultdict(int)
    if not w1:
        out[w2] += 1
        return out
    if not w2:
        out[w1] += 1
        return out
    for w, mult in shuffle_words(w1[1:], w2).items():
        out[(w1[0],) + w] += mult
    for w, mult in shuffle_words(w1, w2[1:]).items():
        out[(w2[0],) + w] += mult
    return dict(out)


def word_of_zeta_desc(comp: tuple) -> tuple:
    """Word for descending zeta_GV(s1,...,sk): innermost-first reading,
    word = concat over i = k..1 of ( [1] + [0]*(s_i - 1) )."""
    w = []
    for s in reversed(comp):
        w.append(1)
        w.extend([0] * (s - 1))
    return tuple(w)


def zeta_of_word_desc(w: tuple) -> tuple:
    """Inverse of word_of_zeta_desc; returns None if word is not
    admissible (must start with 1 and end with 0 for convergence...
    here we only require parseability; convergence checked by caller)."""
    if not w or w[0] != 1:
        return None
    comps = []
    i = 0
    while i < len(w):
        assert w[i] == 1
        j = i + 1
        while j < len(w) and w[j] == 0:
            j += 1
        comps.append(j - i)   # s = 1 + number of zeros
        i = j
    return tuple(reversed(comps))


def level1_relations(w: int) -> list:
    """Conv shuffle + conv stuffle relations among zeta_GV doubles of
    weight w, as dicts mapping atoms to Fractions, atoms:
      ('Z2', a, b) [descending], ('P', p, q) [zeta(p)*zeta(q)],
      ('Z1', w) [zeta(w)].  Each relation: sum coef*atom = 0."""
    rels = []
    for p in range(2, w - 1):
        q = w - p
        if q < 2 or q < p:
            continue
        # stuffle: zeta(p)zeta(q) = Z2(p,q)+Z2(q,p)+zeta(w)
        rel = defaultdict(F)
        rel[('P', p, q)] += 1
        rel[('Z2', p, q)] -= 1
        rel[('Z2', q, p)] -= 1
        rel[('Z1', w)] -= 1
        rels.append(dict(rel))
        # shuffle via word engine
        rel = defaultdict(F)
        rel[('P', p, q)] += 1
        sh = shuffle_words(word_of_zeta_desc((p,)), word_of_zeta_desc((q,)))
        for word, mult in sh.items():
            comp = zeta_of_word_desc(word)
            assert comp is not None and len(comp) == 2
            rel[('Z2', comp[0], comp[1])] -= mult
        rels.append(dict(rel))
    return rels


def euler_trailing1(w: int) -> dict:
    """Euler 1775 (classical, numerically gated in verify):
    zeta_GV(w-1, 1) = (w-1)/2 zeta(w) - 1/2 sum_{k=2}^{w-2} zeta(k)zeta(w-k)
    as a relation dict."""
    rel = defaultdict(F)
    rel[('Z2', w - 1, 1)] += 1
    rel[('Z1', w)] -= F(w - 1, 2)
    for k in range(2, w - 1):
        p, q = min(k, w - k), max(k, w - k)
        rel[('P', p, q)] += F(1, 2)
    return dict(rel)


def reduce_level1(w: int, basis_keep=()) -> dict:
    """Gaussian-eliminate the level-1 weight-w depth-2 system; return
    dict ('Z2',a,b) -> dict of known-atoms (P/Z1/kept-Z2) -> Fraction."""
    unknowns = [('Z2', a, w - a) for a in range(2, w)]   # b = w-a >= 1
    unknowns = [u for u in unknowns if u not in basis_keep]
    rels = level1_relations(w) + [euler_trailing1(w)]
    # eliminate: treat kept basis + P + Z1 as knowns (move right)
    rows = []
    for rel in rels:
        lhs = {u: rel.get(u, F(0)) for u in unknowns if rel.get(u)}
        rhs = {k: -v for k, v in rel.items() if k not in lhs}
        rows.append((dict(lhs), rhs))
    solved = {}
    for u in list(unknowns):
        piv = None
        for i, (lhs, rhs) in enumerate(rows):
            if lhs.get(u):
                piv = i
                break
        if piv is None:
            continue
        lhs, rhs = rows.pop(piv)
        cu = lhs.pop(u)
        # u = (rhs - lhs_rest)/cu
        expr = defaultdict(F)
        for k, v in rhs.items():
            expr[k] += v / cu
        for k, v in lhs.items():
            expr[(('UNK',) + k)] -= v / cu   # placeholder, resolved below
        solved[u] = dict(expr)
        # substitute into remaining rows
        for lhs2, rhs2 in rows:
            if lhs2.get(u):
                c = lhs2.pop(u)
                for k, v in solved[u].items():
                    if k[0] == 'UNK':
                        lhs2[k[1:]] = lhs2.get(k[1:], F(0)) + c * v
                    else:
                        rhs2[k] = rhs2.get(k, F(0)) - c * v
    # back-substitute UNK placeholders
    changed = True
    while changed:
        changed = False
        for u, expr in solved.items():
            new = defaultdict(F)
            for k, v in expr.items():
                if k[0] == 'UNK' and k[1:] in solved:
                    for k2, v2 in solved[k[1:]].items():
                        new[k2] += v * v2
                    changed = True
                else:
                    new[k] += v
            solved[u] = {k: v for k, v in new.items() if v}
    # final sanity: no UNK left except kept basis members
    for u, expr in solved.items():
        for k in expr:
            if k[0] == 'UNK':
                raise RuntimeError(f"unresolved unknown {k} in {u}")
    return solved


# ----------------------------------------------------------------------
# SECTION 5: numerics -- evaluators and caches.
# ----------------------------------------------------------------------

_cache = {}

def load_caches():
    for fn in ('s3_pslq_cache.json', 's3_w10_value_cache.json'):
        try:
            d = json.load(open(DATA + "\\" + fn))
            for k, v in d.items():
                _cache[k] = v if isinstance(v, str) else v.get('value', v)
        except FileNotFoundError:
            pass


def cached_t(name: str, dps: int):
    """Look up t-value from cache at >= dps digits, e.g. 't3(2,4,4)'."""
    best = None
    bestp = -1
    for k, v in _cache.items():
        if k.startswith(name + '@'):
            p = int(k.split('@')[1])
            if p > bestp:
                best, bestp = v, p
    if best is None or bestp < dps:
        return None
    return mpf(best) if not isinstance(best, str) else mp.mpf(best)


def lam_num(k: int) -> mpf:
    return (1 - mpf(2) ** (-k)) * mpzeta(k)


def t2_gv_num(a: int, b: int, N: int = 400) -> mpf:
    """t_GV(a,b) = sum_{n1>n2, odd} n1^-a n2^-b, Levin-safe via Hurwitz
    tails (log-free; requires b >= 2 handled, b=1 via cache)."""
    from mpmath import nsum
    if b == 1:
        v = cached_t(f"t2({a},1)", 50)
        if v is not None:
            return v
        raise ValueError("t2(a,1) needs cache")
    # inner prefix S_b(m) over odd k < (2m+1): use exact partial sums via
    # zeta_h: sum_{odd k >= K} k^-b = 2^-b zeta_h(b, K/2 ... ) -- do directly:
    # sum over odd n2 < n1:  S(n1) = lam(b) - tail_b(n1) with
    # tail_b(n1) = sum_{odd k >= n1} k^-b = 2^-b zeta(b, n1/2) ... careful:
    # odd k = 2j+1 >= n1  (n1 odd)  -> j >= (n1-1)/2:
    # tail = sum_{j >= (n1-1)/2} (2j+1)^-b = 2^-b zeta_h(b, (n1+1)/2 - 1/2)
    #      = 2^-b zeta_h(b, n1/2)
    lamb = lam_num(b)
    def term(i):
        n1 = 2 * i + 1            # odd, i >= 1 means n1 >= 3
        tail = mpf(2) ** (-b) * mpzeta(b, n1 / mpf(2))
        return n1 ** mpf(-a) * (lamb - tail)
    return nsum(term, [1, mp.inf], method='levin')


def zeta2_gv_num(a: int, b: int) -> mpf:
    """zeta_GV(a,b) descending, a >= 2, b >= 2: Levin-safe."""
    from mpmath import nsum
    zb = mpzeta(b)
    def term(n):
        return n ** mpf(-a) * (zb - mpzeta(b, n))
    return nsum(term, [2, mp.inf], method='levin')


def zeta2_gv_trailing1_num(s: int, N: int = 4000, kmax: int = 60) -> mpf:
    """zeta_GV(s,1) = sum_{n>=2} n^-s H_{n-1}; partial sum + asymptotic
    tail via d/ds Hurwitz (NO acceleration on the log-modulated tail)."""
    g = mpeuler
    head = mpf(0)
    H = mpf(0)
    for n in range(2, N + 1):
        H += mpf(1) / (n - 1)
        head += H * mpf(n) ** (-s)
    # tail: H_{n-1} = ln n + g - 1/(2n) - sum_{k>=1} B_2k/(2k n^2k)
    a0 = N + 1
    tail = -mpzeta(s, a0, 1) + g * mpzeta(s, a0) - mpzeta(s + 1, a0) / 2
    for k in range(1, kmax + 1):
        tail -= mpbern(2 * k) / (2 * k) * mpzeta(s + 2 * k, a0)
    return head + tail


def alt1_num(n: int) -> mpf:
    """zeta_GV(n; -) = sum (-1)^m m^-n = -eta(n) = (2^(1-n)-1) zeta(n);
    n=1 -> -ln 2."""
    if n == 1:
        return -mp.log(2)
    return (mpf(2) ** (1 - n) - 1) * mpzeta(n)


# ----------------------------------------------------------------------
# SECTION 6: symbolic catalogue (sympy) and assembly.
# ----------------------------------------------------------------------

def sympy_catalogue():
    import sympy as sp
    pi = sp.pi
    z3, z5, z7, z9, z53 = sp.symbols('z3 z5 z7 z9 z53', positive=True)
    ln2 = sp.symbols('ln2', positive=True)
    a82, a73, a64 = sp.symbols('a82 a73 a64')
    Z = {2: sp.zeta(2), 4: sp.zeta(4), 6: sp.zeta(6), 8: sp.zeta(8),
         10: sp.zeta(10), 3: z3, 5: z5, 7: z7, 9: z9}
    lam = {k: (1 - sp.Rational(1, 2 ** k)) * Z[k] for k in range(2, 11)}
    cat = {}
    # --- depth-1
    for k in range(2, 11):
        cat[('t1', k)] = lam[k]
        cat[('z1', k)] = Z[k]
    cat[('T',)] = ln2
    # --- GV t-doubles, ground-truthed Hoffman App-A rows (descending):
    T = {}
    T[(2, 2)] = (lam[2] ** 2 - lam[4]) / 2
    T[(3, 3)] = (lam[3] ** 2 - lam[6]) / 2
    T[(4, 4)] = (lam[4] ** 2 - lam[8]) / 2
    T[(5, 5)] = (lam[5] ** 2 - lam[10]) / 2
    T[(4, 2)] = -sp.Rational(1, 7) * lam[6] + sp.Rational(1, 7) * lam[3] ** 2
    T[(2, 4)] = sp.Rational(11, 28) * lam[6] - sp.Rational(1, 7) * lam[3] ** 2
    T[(4, 3)] = (-sp.Rational(1, 2) * lam[7]
                 + sp.Rational(6, 7) * lam[3] * lam[4]
                 - sp.Rational(10, 31) * lam[2] * lam[5])
    T[(3, 4)] = lam[3] * lam[4] - lam[7] - T[(4, 3)]
    T[(2, 3)] = (-sp.Rational(1, 2) * lam[5]
                 + sp.Rational(4, 7) * lam[2] * lam[3])
    T[(3, 2)] = lam[2] * lam[3] - lam[5] - T[(2, 3)]
    # w10 doubles via antisymmetric generators:
    A = {(8, 2): a82, (7, 3): a73, (6, 4): a64}
    for (x, y), sym in A.items():
        s = lam[x] * lam[y] - lam[10]
        T[(x, y)] = (s + sym) / 2
        T[(y, x)] = (s - sym) / 2
    for k, v in T.items():
        cat[('T2gv',) + k] = v
    # --- level-1 zeta doubles, weights 4, 6, 8 (basis zeta53 at w8):
    red = {}
    for w in (4, 6):
        red.update(reduce_level1(w))
    red.update(reduce_level1(8, basis_keep=(('Z2', 5, 3),)))
    def resolve(atomdict):
        e = sp.Integer(0)
        for k, v in atomdict.items():
            fr = sp.Rational(v.numerator, v.denominator)
            if k[0] == 'P':
                e += fr * Z[k[1]] * Z[k[2]]
            elif k[0] == 'Z1':
                e += fr * Z[k[1]]
            elif k == ('Z2', 5, 3):
                e += fr * z53
            else:
                raise RuntimeError(f"unexpected atom {k}")
        return e
    for u, expr in red.items():
        cat[('Z2gv', u[1], u[2])] = resolve(expr)
    cat[('Z2gv', 5, 3)] = z53
    # w2 reg + small cases used by L2 substitutions:
    cat[('Z2gv_r11',)] = -Z[2] / 2          # zeta_CH-reg(1,1) at T=0
    syms = dict(pi=pi, z3=z3, z5=z5, z7=z7, z9=z9, z53=z53, ln2=ln2,
                a82=a82, a73=a73, a64=a64, lam=lam, Z=Z)
    return cat, syms


def coef_to_sympy(coef: dict, cat: dict):
    """Convert an extracted coefficient (CH-convention atoms) to sympy,
    using the conversion t_CH(a,b)=t_GV(b,a), z_CH(a,b)=z_GV(b,a)."""
    import sympy as sp
    expr = sp.Integer(0)
    unresolved = []
    for term, frac in coef.items():
        fr = sp.Rational(frac.numerator, frac.denominator)
        prod = sp.Integer(1)
        ok = True
        for atom in term:
            if atom[0] == 't1' or atom == ('T',) or atom[0] == 'z1':
                prod *= cat[atom]
            elif atom[0] == 't2':
                prod *= cat[('T2gv', atom[2], atom[1])]
            elif atom[0] == 'z2':
                key = ('Z2gv', atom[2], atom[1])
                if key in cat:
                    prod *= cat[key]
                else:
                    ok = False
            elif atom[0] == 'z2r':
                # zeta_CH-reg(n,1)|_{T=0} = -zeta_CH(1,n) - zeta(n+1)
                #                        = -zeta_GV(n,1) - zeta(n+1)
                n = atom[1]
                if n == 1:
                    prod *= cat[('Z2gv_r11',)]
                else:
                    prod *= (-cat[('Z2gv', n, 1)] - cat[('z1', n + 1)])
            elif atom[0] in ('t3', 't3r', 't2r'):
                ok = False
            else:
                raise RuntimeError(f"unknown atom {atom}")
        if ok:
            expr += fr * prod
        else:
            unresolved.append((term, frac))
    return expr, unresolved


# ----------------------------------------------------------------------
# SECTION 7: subcommands.
# ----------------------------------------------------------------------

def run_extract():
    print("building Theorem 2.21 m=3 phi=0 series (CH conv)...")
    P3 = build_theorem_m3()
    out = {}
    for (a, b, c) in [(2, 4, 4), (3, 3, 4), (6, 2, 2)]:
        coef = extract(P3, a, b, c)
        out[f"({a},{b},{c})"] = {repr(t): str(f) for t, f in coef.items()}
        print(f"\n=== coefficient at (a,b,c)=({a},{b},{c}) "
              f"[identity: == 0] ===")
        for t, f in sorted(coef.items()):
            print(f"   {f}  *  {t}")
    json.dump(out, open(DATA + r"\s3_w10_ident_extraction.json", 'w'),
              indent=1)
    print("\nwrote s3_w10_ident_extraction.json")
    return P3


# --- numeric resolution of atoms (CH convention) ----------------------

_numcache = {}

def atom_num(atom) -> mpf:
    key = (atom, mp.dps)
    if key in _numcache:
        return _numcache[key]
    if atom == ('T',):
        v = mp.log(2)
    elif atom[0] == 't1':
        v = lam_num(atom[1])
    elif atom[0] == 'z1':
        v = mpzeta(atom[1])
    elif atom[0] == 't2':
        # CH(n1,n2) = GV(n2,n1); prefer cache, else evaluator
        a, b = atom[2], atom[1]          # GV indices
        cv = cached_t(f"t2({a},{b})", 60)
        v = cv if cv is not None else t2_gv_num(a, b)
    elif atom[0] == 't3':
        a, b, c = atom[3], atom[2], atom[1]   # GV indices (reversed)
        cv = cached_t(f"t3({a},{b},{c})", 60)
        if cv is None:
            raise RuntimeError(f"no cache for t3({a},{b},{c})")
        v = cv
    elif atom[0] == 'z2':
        a, b = atom[2], atom[1]          # GV indices
        v = zeta2_gv_trailing1_num(a) if b == 1 else zeta2_gv_num(a, b)
    elif atom[0] == 'z2r':
        n = atom[1]
        if n == 1:
            v = -mpzeta(2) / 2
        else:
            v = -zeta2_gv_trailing1_num(n) - mpzeta(n + 1)
    else:
        raise RuntimeError(f"no numeric for {atom}")
    _numcache[key] = v
    return v


def coef_num(coef: dict) -> mpf:
    tot = mpf(0)
    for term, frac in coef.items():
        prod = mpf(frac.numerator) / mpf(frac.denominator)
        for atom in term:
            prod *= atom_num(atom)
        tot += prod
    return tot


def gate(name: str, residual, tol) -> bool:
    ok = abs(residual) < tol
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: residual = "
          f"{mp.nstr(abs(residual), 5)}  (gate {tol})")
    return ok


def run_verify():
    results = []
    print("\n--- gate W: word-engine convention (zeta(2)zeta(3)) ---")
    sh = shuffle_words(word_of_zeta_desc((2,)), word_of_zeta_desc((3,)))
    tot = mpf(0)
    for word, mult in sh.items():
        comp = zeta_of_word_desc(word)
        a, b = comp
        tot += mult * (zeta2_gv_num(a, b) if b >= 2
                       else zeta2_gv_trailing1_num(a))
    results.append(gate("shuffle zeta(2)*zeta(3)",
                        tot - mpzeta(2) * mpzeta(3), mpf(10) ** (-mp.dps + 15)))

    print("\n--- gate E: Euler trailing-1 rows w=4,6,8 ---")
    for w in (4, 6, 8):
        rel = euler_trailing1(w)
        tot = mpf(0)
        for k, v in rel.items():
            fv = mpf(v.numerator) / mpf(v.denominator)
            if k[0] == 'Z2':
                tot += fv * (zeta2_gv_trailing1_num(k[1]) if k[2] == 1
                             else zeta2_gv_num(k[1], k[2]))
            elif k[0] == 'Z1':
                tot += fv * mpzeta(k[1])
            else:
                tot += fv * mpzeta(k[1]) * mpzeta(k[2])
        results.append(gate(f"Euler zeta_GV({w-1},1)", tot,
                            mpf(10) ** (-mp.dps + 15)))

    print("\n--- gate L: level-1 reductions w=4,6,8 (vs evaluators) ---")
    red = {}
    for w in (4, 6):
        red.update(reduce_level1(w))
    red.update(reduce_level1(8, basis_keep=(('Z2', 5, 3),)))
    z53num = zeta2_gv_num(5, 3)
    for u, expr in sorted(red.items()):
        a, b = u[1], u[2]
        lhs = zeta2_gv_trailing1_num(a) if b == 1 else zeta2_gv_num(a, b)
        rhs = mpf(0)
        for k, v in expr.items():
            fv = mpf(v.numerator) / mpf(v.denominator)
            if k[0] == 'P':
                rhs += fv * mpzeta(k[1]) * mpzeta(k[2])
            elif k[0] == 'Z1':
                rhs += fv * mpzeta(k[1])
            elif k == ('Z2', 5, 3):
                rhs += fv * z53num
            else:
                raise RuntimeError(f"atom {k}")
        results.append(gate(f"reduce zeta_GV({a},{b})", lhs - rhs,
                            mpf(10) ** (-mp.dps + 15)))

    print("\n--- gate D: pair-difference stuffles (vs caches) ---")
    st = pair_difference_relations()
    # RA: lam2*t2GV(4,4) = t3(2,4,4)+t3(4,2,4)+t3(4,4,2)+t2(6,4)+t2(4,6)
    def stuffle_num(x, y, z):
        tot = -lam_num(x) * (cached_t(f"t2({y},{z})", 60) or t2_gv_num(y, z))
        for atom, fr in stuffle_1_into_2(x, y, z).items():
            fv = mpf(fr.numerator) / mpf(fr.denominator)
            if atom[0] == 'T3':
                cv = cached_t(f"t3({atom[1]},{atom[2]},{atom[3]})", 60)
                if cv is None:
                    raise RuntimeError(f"no cache t3{atom[1:]}")
                tot += fv * cv
            else:
                tot += fv * (cached_t(f"t2({atom[1]},{atom[2]})", 60)
                             or t2_gv_num(atom[1], atom[2]))
        return tot
    for (x, y, z), nm in [((2, 4, 4), 'RA'), ((4, 2, 4), 'RB'),
                          ((3, 3, 4), 'RC'), ((3, 4, 3), 'RD')]:
        results.append(gate(f"stuffle {nm} t({x})*t2({y},{z})",
                            stuffle_num(x, y, z), mpf(10) ** (-78)))

    print("\n--- gate M2: Theorem 2.21 m=2 extraction at (4,4), (2,6) ---")
    P2 = build_theorem_m2()
    for (a, b) in [(4, 4), (2, 6)]:
        coef = extract(P2, a, b)
        results.append(gate(f"Thm2.21 m=2 at ({a},{b})", coef_num(coef),
                            mpf(10) ** (-78)))

    print("\n--- gate M3: Theorem 2.21 m=3 identities ---")
    P3 = build_theorem_m3()
    for (a, b, c) in [(6, 2, 2), (2, 4, 4), (3, 3, 4)]:
        coef = extract(P3, a, b, c)
        results.append(gate(f"Thm2.21 m=3 at ({a},{b},{c})",
                            coef_num(coef), mpf(10) ** (-78)))

    nfail = sum(1 for r in results if not r)
    print(f"\n==== verify: {len(results) - nfail}/{len(results)} PASS ====")
    return nfail == 0


def run_assemble():
    import sympy as sp
    cat, syms = sympy_catalogue()
    P3 = build_theorem_m3()
    pi, z3, z5, z7, z9, z53, ln2 = (syms[k] for k in
                                    ('pi', 'z3', 'z5', 'z7', 'z9',
                                     'z53', 'ln2'))
    a82, a73, a64 = syms['a82'], syms['a73'], syms['a64']
    lam = syms['lam']

    def T2(a, b):
        return cat[('T2gv', a, b)]

    # --- the two symmetry sums P1, P2 (solve extraction == 0 for t3 sum)
    sums = {}
    for (a, b, c), key in [((2, 4, 4), 'P1'), ((3, 3, 4), 'P2')]:
        coef = extract(P3, a, b, c)
        expr, unresolved = coef_to_sympy(coef, cat)
        # unresolved must be exactly the two t3 atoms with coeff +1
        atoms = sorted(t for t, f in unresolved)
        assert all(f == 1 for _, f in unresolved), unresolved
        assert len(unresolved) == 2, unresolved
        sums[key] = sp.expand(-expr)   # t3+t3rev = -expr
        print(f"{key}: t3 pair sum [GV {tuple(reversed((a,b,c)))} + "
              f"{(a,b,c)}] =")
        sp.pprint(sp.nsimplify(sums[key]))
    P1, P2 = sums['P1'], sums['P2']

    # --- pair differences from RA-RB and RC-RD (exact stuffles)
    D1 = sp.expand(lam[2] * T2(4, 4) - lam[4] * T2(2, 4)
                   - T2(4, 6) + T2(2, 8))
    D2 = sp.expand(sp.Rational(1, 2) * (lam[3] * (T2(3, 4) - T2(4, 3))
                   - T2(6, 4) - T2(3, 7) + T2(7, 3) + T2(4, 6)))
    # NOTE: T2(3,7)/T2(7,3) enter via a73; T2(6,4)/T2(4,6) via a64; etc.

    t3_442 = sp.expand((P1 + D1) / 2)
    t3_244 = sp.expand((P1 - D1) / 2)
    t3_334 = sp.expand((P2 + D2) / 2)
    t3_433 = sp.expand((P2 - D2) / 2)
    print("\nt3_GV(2,4,4) ="); sp.pprint(t3_244)
    print("\nt3_GV(3,3,4) ="); sp.pprint(t3_334)

    # --- assembly-normalized corollary (v4.7.0 boxed identity, GV conv)
    a43 = sp.expand(T2(4, 3) - T2(3, 4))
    a42 = sp.expand(T2(4, 2) - T2(2, 4))
    W10a_over_256 = sp.expand(
        6 * lam[10] - lam[2] * lam[8] - 2 * lam[3] * lam[7]
        - 3 * lam[4] * lam[6] - 4 * lam[3] * a43 - 2 * lam[4] * a42
        + a82 + 2 * a73 - 3 * a64 - 4 * t3_244 - 8 * t3_334)
    print("\nW10_assembly/256 (fully substituted) =")
    sp.pprint(sp.collect(W10a_over_256, [a82, a73, a64, z53]))

    # --- Ident from the assembly sprint (26 terms)
    Z = syms['Z']
    Ident = (
        128 * pi**2 * ln2**2 - sp.Rational(32, 3) * pi**4 * ln2**2
        + sp.Rational(7232, 81) * pi**2 * ln2
        - sp.Rational(35504, 243) * pi**4 * ln2
        + sp.Rational(352, 15) * pi**6 * ln2
        - sp.Rational(314, 315) * pi**8 * ln2
        - 96 * pi**2 * ln2 * z3 + 16 * pi**4 * ln2 * z3
        + 20 * pi**2 * z3**2 - sp.Rational(2200, 27) * pi**2 * z3
        + sp.Rational(4988, 81) * pi**4 * z3
        - sp.Rational(66, 5) * pi**6 * z3
        + sp.Rational(157, 210) * pi**8 * z3
        - 80 * pi**2 * ln2 * z5
        + sp.Rational(7460, 81) * pi**2 * z5
        + sp.Rational(85, 3) * pi**4 * z5
        - sp.Rational(11, 3) * pi**6 * z5
        - 7 * pi**2 * z7 + sp.Rational(91, 12) * pi**4 * z7
        - 30 * pi**2 * z9
        + sp.Rational(4096, 6561) * pi**2
        - sp.Rational(109240, 19683) * pi**4
        + sp.Rational(27427, 1215) * pi**6
        - sp.Rational(288859, 51030) * pi**8
        + sp.Rational(1381, 2835) * pi**10
        - sp.Rational(17851, 1247400) * pi**12)

    S3_closed = sp.expand(Ident + 256 * W10a_over_256)
    S3_collected = sp.collect(S3_closed, [a82, a73, a64, z53])
    print("\n===== S^(3) FULL CLOSED FORM =====")
    sp.pprint(S3_collected)

    # --- numeric gate at 200 dps
    mp.dps = 210
    _numcache.clear()
    a82n = (cached_t("t2(8,2)", 200) - cached_t("t2(2,8)", 200))
    a73n = (cached_t("t2(7,3)", 200) - cached_t("t2(3,7)", 200))
    a64n = (cached_t("t2(6,4)", 200) - cached_t("t2(4,6)", 200))
    z53n = zeta2_gv_num(5, 3)
    subs = {pi: mp.pi, z3: mpzeta(3), z5: mpzeta(5), z7: mpzeta(7),
            z9: mpzeta(9), z53: z53n, ln2: mp.log(2),
            a82: a82n, a73: a73n, a64: a64n}
    # evaluate via lambdify-free walk (sympy evalf with mp-backed floats):
    expr_num = S3_closed
    val = mpf(0)
    poly = sp.Poly(S3_closed, pi, z3, z5, z7, z9, z53, ln2, a82, a73, a64)
    for monom, coeff in poly.terms():
        term = mpf(str(sp.Rational(coeff)))
        for base, e in zip((mp.pi, mpzeta(3), mpzeta(5), mpzeta(7),
                            mpzeta(9), z53n, mp.log(2), a82n, a73n, a64n),
                           monom):
            if e:
                term *= base ** e
        val += term
    rec = json.load(open(DATA + r"\s3_closure_recompute.json"))
    canon = mp.mpf(rec['S3_canonical_200dps'])
    print(f"\nclosed form numeric: {mp.nstr(val, 40)}")
    print(f"canonical:           {mp.nstr(canon, 40)}")
    resid = val - canon
    ok = gate("FINAL S^(3) closed-form vs canonical", resid, mpf(10)**(-180))

    # --- individual triple gates (validates the sum/difference split
    #     member-by-member, not just in the assembled combination)
    print("\n--- individual t3 closed-form gates (200 dps) ---")
    gens = (mp.pi, mpzeta(3), mpzeta(5), mpzeta(7), mpzeta(9), z53n,
            mp.log(2), a82n, a73n, a64n)
    def eval_expr(e):
        p = sp.Poly(e, pi, z3, z5, z7, z9, z53, ln2, a82, a73, a64)
        v = mpf(0)
        for monom, coeff in p.terms():
            t = mpf(str(sp.Rational(coeff)))
            for base, ee in zip(gens, monom):
                if ee:
                    t *= base ** ee
            v += t
        return v
    indiv = {}
    for nm, e in [('t3(2,4,4)', t3_244), ('t3(3,3,4)', t3_334),
                  ('t3(4,4,2)', t3_442), ('t3(4,3,3)', t3_433)]:
        cv = cached_t(nm.replace('t3', 't3'), 190)
        r = eval_expr(e) - cv
        indiv[nm] = mp.nstr(abs(r), 5)
        ok = gate(f"closed form {nm} vs cache", r, mpf(10) ** (-180)) and ok
    out = {
        'P1': str(P1), 'P2': str(P2), 'D1': str(D1), 'D2': str(D2),
        't3_244': str(t3_244), 't3_334': str(t3_334),
        't3_442': str(t3_442), 't3_433': str(t3_433),
        'S3_closed_form': str(S3_collected),
        'final_residual': mp.nstr(abs(resid), 8),
        'final_gate': 'PASS' if ok else 'FAIL',
    }
    json.dump(out, open(DATA + r"\s3_w10_ident_result.json", 'w'), indent=1)
    print("wrote s3_w10_ident_result.json")
    return ok


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all'
    mp.dps = 100
    load_caches()
    if cmd in ('extract',):
        run_extract()
    elif cmd == 'verify':
        run_verify()
    elif cmd == 'assemble':
        run_assemble()
    elif cmd == 'all':
        run_extract()
        ok = run_verify()
        if ok:
            run_assemble()


if __name__ == '__main__':
    main()
