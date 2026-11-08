"""VV stage A, route 1: group theory on S^3 = SU(2), left-trivialized frame.

Record: docs/predictions/2026-11-07_physical_vector_vertex.md.

Construction (all exact):
  * Functions on SU(2): Wigner functions D^j_{m m'}(g); right translation
    acts on the second index m' as the standard spin-j ket action, left
    translation on the first index m.  X_i = -2i R_i (left-invariant fields,
    unit S^3, [X_i, X_j] = 2 eps_ijk X_k).
  * Spinors: f(g) (x) chi, chi in C^2 (frame index).  Dirac operator in the
    left-invariant frame:  D = sum_i (i sigma_i) X_i + 3/2 = 2 sigma.R + 3/2.
    Eigenstates are standard CG couplings of (j, 1/2) -> J; J = j+1/2 gives
    D = 2j + 3/2, J = j-1/2 gives D = -(2j + 1/2).
  * 1-forms: f(g) (x) e_mu (spin-1 frame index, spherical basis).  *d =
    -2 S.R - 2; coupled J = j+1 gives -(q+1) (q = 2j+1), J = j-1 gives +(q+1)
    (q = 2j-1), J = j gives 0 (exact forms).
  * Content labels: record's "first" factor = J (right SU(2), frame
    included); "second" = j (left SU(2)).  With this, D = +(n+3/2) is
    ((n+1)/2, n/2) and component A = ((q+1)/2, (q-1)/2): exactly the
    record's convention.
  * Vertex V = int <psi1, c(A) psi2> dg (Haar normalized), c(A) = sum_i a_i sigma_i
    (the factor i of gamma^i = i sigma_i is dropped; irrelevant to zero/|.|^2).
    The integral factorizes:
       V = sqrt((2j1+1)(2jA+1)(2j2+1))/(2j1+1) * <jA mA j2 m2|j1 m1>  *  Rt(M1, MA, M2)
    with Rt the frame-side contraction (computed by brute-force CG sums).

Numbers: elements of the multi-quadratic field as {squarefree s: Fraction c}
meaning sum c sqrt(s).  Zero means exactly zero.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import factorial, gcd
from typing import Dict, Iterable, List, Tuple

# ---------------------------------------------------------------------------
# exact surd arithmetic
# ---------------------------------------------------------------------------
_PRIMES = [p for p in range(2, 200) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


@lru_cache(maxsize=None)
def squarefree_split(n: int) -> Tuple[int, int]:
    """n = s * t^2 with s squarefree; returns (s, t). n > 0, smooth."""
    s, t = 1, 1
    for p in _PRIMES:
        if p * p > n and n > 1 and n < 200 * 200:
            pass
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            t *= p ** (e // 2)
            if e % 2:
                s *= p
        if n == 1:
            break
    if n != 1:
        raise ValueError("non-smooth integer in squarefree_split")
    return s, t


Surd = Dict[int, Fraction]


def surd_sqrt(F: Fraction, sign: Fraction = Fraction(1)) -> Surd:
    """sign * sqrt(F) as a Surd (F >= 0 rational)."""
    if F == 0 or sign == 0:
        return {}
    n, d = F.numerator, F.denominator
    s, t = squarefree_split(n * d)
    return {s: sign * Fraction(t, d)}


def s_add(a: Surd, b: Surd) -> Surd:
    out = dict(a)
    for k, v in b.items():
        w = out.get(k, Fraction(0)) + v
        if w == 0:
            out.pop(k, None)
        else:
            out[k] = w
    return out


def s_mul(a: Surd, b: Surd) -> Surd:
    out: Surd = {}
    for s1, c1 in a.items():
        for s2, c2 in b.items():
            g = gcd(s1, s2)
            s = (s1 // g) * (s2 // g)
            c = c1 * c2 * g
            w = out.get(s, Fraction(0)) + c
            if w == 0:
                out.pop(s, None)
            else:
                out[s] = w
    return out


def s_scale(a: Surd, r: Fraction) -> Surd:
    if r == 0:
        return {}
    return {k: v * r for k, v in a.items()}


def s_rational(a: Surd) -> Fraction:
    """Return a as a rational; raise if it has irrational part."""
    for k in a:
        if k != 1:
            raise ValueError(f"not rational: {a}")
    return a.get(1, Fraction(0))


# ---------------------------------------------------------------------------
# Clebsch-Gordan (Racah formula), doubled-integer arguments
# ---------------------------------------------------------------------------
def _half(x2: int) -> int:
    assert x2 % 2 == 0
    return x2 // 2


@lru_cache(maxsize=None)
def cg2(tj1: int, tm1: int, tj2: int, tm2: int, tJ: int, tM: int) -> Tuple[Fraction, Fraction]:
    """<j1 m1 j2 m2 | J M> = r * sqrt(F); returns (r, F). Doubled args."""
    zero = (Fraction(0), Fraction(0))
    if tm1 + tm2 != tM:
        return zero
    if abs(tm1) > tj1 or abs(tm2) > tj2 or abs(tM) > tJ:
        return zero
    if (tj1 + tm1) % 2 or (tj2 + tm2) % 2 or (tJ + tM) % 2:
        return zero
    if (tj1 + tj2 + tJ) % 2:
        return zero
    if tJ < abs(tj1 - tj2) or tJ > tj1 + tj2:
        return zero
    a = _half(tJ + tj1 - tj2)
    b = _half(tJ - tj1 + tj2)
    c = _half(tj1 + tj2 - tJ)
    d = _half(tj1 + tj2 + tJ) + 1
    F = Fraction((tJ + 1) * factorial(a) * factorial(b) * factorial(c), factorial(d))
    F *= (factorial(_half(tJ + tM)) * factorial(_half(tJ - tM))
          * factorial(_half(tj1 - tm1)) * factorial(_half(tj1 + tm1))
          * factorial(_half(tj2 - tm2)) * factorial(_half(tj2 + tm2)))
    r = Fraction(0)
    k = 0
    while True:
        args = [k, _half(tj1 + tj2 - tJ) - k, _half(tj1 - tm1) - k, _half(tj2 + tm2) - k,
                _half(tJ - tj2 + tm1) + k, _half(tJ - tj1 - tm2) + k]
        if args[1] < 0 or args[2] < 0 or args[3] < 0:
            break
        if min(args) >= 0:
            den = 1
            for x in args:
                den *= factorial(x)
            r += Fraction((-1) ** k, den)
        k += 1
    if r == 0:
        return zero
    # absorb |r| into F to keep a single surd: CG = sign(r) sqrt(F r^2)
    return (Fraction(1 if r > 0 else -1), F * r * r)


def cg_surd(tj1: int, tm1: int, tj2: int, tm2: int, tJ: int, tM: int) -> Surd:
    r, F = cg2(tj1, tm1, tj2, tm2, tJ, tM)
    return surd_sqrt(F, r)


def mrange(tj: int) -> List[int]:
    return list(range(-tj, tj + 1, 2))


# sigma_mu (spherical) matrix elements <s1| sigma_mu |s2>, doubled s, mu in {-1,0,1}
# sigma_{+1} = -(sx + i sy)/sqrt2 = -sqrt2 |up><dn| ; sigma_0 = sz ; sigma_{-1} = +sqrt2 |dn><up|
def sigma_el(ts1: int, mu: int, ts2: int) -> Surd:
    if mu == 0:
        if ts1 == ts2:
            return {1: Fraction(1 if ts1 == 1 else -1)}
        return {}
    if mu == 1 and ts1 == 1 and ts2 == -1:
        return {2: Fraction(-1)}
    if mu == -1 and ts1 == -1 and ts2 == 1:
        return {2: Fraction(1)}
    return {}


# ---------------------------------------------------------------------------
# content
# ---------------------------------------------------------------------------
def spinor_content(sign: int, n: int) -> Tuple[int, int]:
    """(tJ, tj) doubled: record-first = J (right+frame), second = j (left)."""
    if sign > 0:
        return (n + 1, n)       # J = (n+1)/2, j = n/2  (J = j + 1/2) -> D = n + 3/2
    return (n, n + 1)           # J = n/2, j = (n+1)/2 (J = j - 1/2) -> D = -(n + 3/2)


def vector_content(comp: str, q: int) -> Tuple[int, int]:
    if comp == "A":
        return (q + 1, q - 1)   # J = j + 1, *d = -(q+1)
    return (q - 1, q + 1)       # J = j - 1, *d = +(q+1)


# ---------------------------------------------------------------------------
# frame-side contraction Rt and the vertex summaries
# ---------------------------------------------------------------------------
def rt_tensor(tJ1: int, tj1: int, tJA: int, tjA: int, tJ2: int, tj2: int) -> Dict[Tuple[int, int, int], Surd]:
    """Rt(M1, MA, M2) = sum <jA a j2 b|j1 m1'> <j1 m1' 1/2 s1|J1 M1>
       <jA a 1 mu|JA MA> <j2 b 1/2 s2|J2 M2> <s1|sigma_mu|s2>  (nonzero entries only)."""
    out: Dict[Tuple[int, int, int], Surd] = {}
    for tM1 in mrange(tJ1):
        for tMA in mrange(tJA):
            for tM2 in mrange(tJ2):
                acc: Surd = {}
                for ta in mrange(tjA):
                    tmu = tMA - ta
                    if tmu not in (-2, 0, 2):
                        continue
                    mu = tmu // 2
                    cA = cg_surd(tjA, ta, 2, tmu, tJA, tMA)
                    if not cA:
                        continue
                    for ts2 in (-1, 1):
                        tb = tM2 - ts2
                        if abs(tb) > tj2:
                            continue
                        c2 = cg_surd(tj2, tb, 1, ts2, tJ2, tM2)
                        if not c2:
                            continue
                        tm1p = ta + tb
                        cL = cg_surd(tjA, ta, tj2, tb, tj1, tm1p)
                        if not cL:
                            continue
                        ts1 = tM1 - tm1p
                        if ts1 not in (-1, 1):
                            continue
                        c1 = cg_surd(tj1, tm1p, 1, ts1, tJ1, tM1)
                        if not c1:
                            continue
                        sg = sigma_el(ts1, mu, ts2)
                        if not sg:
                            continue
                        term = s_mul(s_mul(s_mul(cA, c2), s_mul(cL, c1)), sg)
                        acc = s_add(acc, term)
                if acc:
                    out[(tM1, tMA, tM2)] = acc
    return out


def left_open(tj1: int, tjA: int, tj2: int) -> bool:
    """Is <jA mA j2 m2|j1 m1> nonzero for some m (computed, not assumed)."""
    for tmA in mrange(tjA):
        for tm2 in mrange(tj2):
            tm1 = tmA + tm2
            if abs(tm1) <= tj1 and cg2(tjA, tmA, tj2, tm2, tj1, tm1)[1] != 0:
                return True
    return False


def vertex_route1(s1: int, n1: int, comp: str, q: int, s2: int, n2: int) -> Dict:
    """Open/closed and Sum_m |V|^2 for one (sign pair, component, n1, n2, q)."""
    tJ1, tj1 = spinor_content(s1, n1)
    tJA, tjA = vector_content(comp, q)
    tJ2, tj2 = spinor_content(s2, n2)
    lo = left_open(tj1, tjA, tj2)
    rt = rt_tensor(tJ1, tj1, tJA, tjA, tJ2, tj2) if lo else {}
    open_ = lo and bool(rt)
    if not open_:
        return {"open": False, "S": Fraction(0)}
    sum_rt2: Surd = {}
    for v in rt.values():
        sum_rt2 = s_add(sum_rt2, s_mul(v, v))
    S = Fraction((tjA + 1) * (tj2 + 1)) * s_rational(sum_rt2)
    # m-conservation check in both factors (rule 4): Rt nonzero only at M1 = MA + M2
    mcons = all(k[0] == k[1] + k[2] for k in rt)
    return {"open": True, "S": S, "m_conserved": mcons}


def diagonal_ground(sign: int) -> Dict:
    """SA4: <a| c(A) |a> for the n = 0 states a, q = 1, both components, all photon m.
    Full V including the left CG (n = 0 has j or J = 0 on one side)."""
    tJ, tj = spinor_content(sign, 0)
    res = {}
    for comp in ("A", "B"):
        tJA, tjA = vector_content(comp, 1)
        if not left_open(tj, tjA, tj):
            res[comp] = {"left_open": False}
            continue
        rt = rt_tensor(tJ, tj, tJA, tjA, tJ, tj)
        per_state = {}
        for tmL in mrange(tj):
            for tM in mrange(tJ):
                tot: Surd = {}
                for tmAL in mrange(tjA):
                    L = cg_surd(tjA, tmAL, tj, tmL, tj, tmL + tmAL) if abs(tmL + tmAL) <= tj else {}
                    if tmAL != 0:
                        L = {}  # diagonal: same left label on both sides => photon mL = 0
                    if not L:
                        continue
                    for tMA in mrange(tJA):
                        v = rt.get((tM, tMA, tM))
                        if v is None:
                            continue
                        amp = s_mul(L, v)
                        tot = s_add(tot, s_mul(amp, amp))
                norm = Fraction((tj + 1) * (tjA + 1) * (tj + 1), (tj + 1) ** 2)
                per_state[f"mL={tmL}/2,M={tM}/2"] = str(norm * s_rational(tot)) if tot else "0"
        res[comp] = {"left_open": True, "sum_mA_|V(a,A,a)|^2": per_state}
    return res


# ---------------------------------------------------------------------------
# C1 / C2 explicit operator checks (route 1)
# ---------------------------------------------------------------------------
def _spin_matrices(tj: int):
    """Standard J_z, J_+, J_- as dicts {(tm_row, tm_col): Surd} for spin tj/2."""
    Jz, Jp, Jm = {}, {}, {}
    for tm in mrange(tj):
        Jz[(tm, tm)] = {1: Fraction(tm, 2)} if tm else {}
        if tm + 2 <= tj:
            # <m+1|J+|m> = sqrt(j(j+1) - m(m+1))
            F = Fraction(tj * (tj + 2) - tm * (tm + 2), 4)
            Jp[(tm + 2, tm)] = surd_sqrt(F)
            Jm[(tm, tm + 2)] = surd_sqrt(F)
    return Jz, Jp, Jm


def _apply(op: Dict, vec: Dict) -> Dict:
    out: Dict = {}
    for (r, c), val in op.items():
        if c in vec and val:
            out[r] = s_add(out.get(r, {}), s_mul(val, vec[c]))
    return {k: v for k, v in out.items() if v}


def dot_SR_on_coupled(tj: int, ts: int, tJ: int, tM: int) -> Tuple[Dict, Dict]:
    """Apply S.R (= Sz Rz + (S+R- + S-R+)/2) to the coupled state |(j s) J M>.
    Returns (coupled vector, S.R vector) as dicts over (tm_orb, tm_spin)."""
    vec = {}
    for tm in mrange(tj):
        tms = tM - tm
        if abs(tms) <= ts:
            c = cg_surd(tj, tm, ts, tms, tJ, tM)
            if c:
                vec[(tm, tms)] = c
    Rz, Rp, Rm = _spin_matrices(tj)
    Sz, Sp, Sm = _spin_matrices(ts)

    def kron_apply(A, B, v):
        out = {}
        for (r1, c1), a in A.items():
            if not a:
                continue
            for (r2, c2), b in B.items():
                if not b:
                    continue
                if (c1, c2) in v:
                    key = (r1, r2)
                    out[key] = s_add(out.get(key, {}), s_mul(s_mul(a, b), v[(c1, c2)]))
        return out

    t1 = kron_apply(Rz, Sz, vec)
    t2 = kron_apply(Rm, Sp, vec)
    t3 = kron_apply(Rp, Sm, vec)
    res = dict(t1)
    for t in (t2, t3):
        for k, v in t.items():
            res[k] = s_add(res.get(k, {}), s_scale(v, Fraction(1, 2)))
    res = {k: v for k, v in res.items() if v}
    return vec, res


def check_eigen(tj: int, ts: int, tJ: int, coef: Fraction, const: Fraction, expect: Fraction) -> bool:
    """(coef * S.R + const) |(j s)J M> == expect |(j s)J M> for all M (exact)."""
    for tM in mrange(tJ):
        vec, sr = dot_SR_on_coupled(tj, ts, tJ, tM)
        keys = set(vec) | set(sr)
        for k in keys:
            lhs = s_add(s_scale(sr.get(k, {}), coef), s_scale(vec.get(k, {}), const))
            rhs = s_scale(vec.get(k, {}), expect)
            if s_add(lhs, s_scale(rhs, Fraction(-1))):
                return False
    return True
