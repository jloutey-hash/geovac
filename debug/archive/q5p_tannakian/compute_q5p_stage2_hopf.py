r"""
Sprint Q5'-Stage2-Hopf driver — Stage 2 first scoping step.

Define the candidate Hopf algebra $\mathcal{H}_{\mathrm{GV}}$ on Mellin-moment-
labelled-by-k data of the truncated Camporesi-Higuchi spectral triple pro-
system, and check Hopf-algebra axioms (coproduct Delta, counit epsilon,
antipode S, coassociativity, bialgebra compatibility, counit compatibility,
antipode property) bit-exactly at n_max in {2, 3}.

Output
------
- debug/data/sprint_q5p_stage2_hopf.json -- exact rational verification panel
- debug/sprint_q5p_stage2_hopf_memo.md -- memo (separate file)

Strategy
--------
The pro-system substrate from v3.60.0 gives:
- Sector idempotents (n, l) with 1 <= n <= n_max, 0 <= l <= n
- N(n_max) = n_max(n_max+3)/2 sectors
- Mellin slot k in {0, 1, 2} indexing M1, M3, M2 (per Paper 18 sec III.7)
- Sector-LOCAL closed forms for chi_{(n,l)} (k=0 grouplike weights)
  and eta_{(n,l)} (k=1 grouplike weights).

The natural connected-graded Hopf algebra candidate H_GV:

Underlying graded vector space:
    H_GV^{(n_max)} = Sym(V_{n_max})  where V_{n_max} is the Q-vector space
    spanned by primitive generators [e_{(n,l)}, k] for sector (n, l)
    and k in {0, 1, 2}.

So dim V_{n_max} = 3 * N(n_max), and H_GV^{(n_max)} is the polynomial algebra
Q[ x_{(n,l), k} ] in 3*N(n_max) variables.

Grading: degree of generator x_{(n,l), k} is the shell number n.
The grade-0 piece is Q (the unit).

Coproduct: primitive on generators
    Delta(x_{(n,l), k}) = x_{(n,l), k} (x) 1 + 1 (x) x_{(n,l), k}
extended as an algebra homomorphism.

Counit: epsilon(x_{(n,l), k}) = 0 on generators, epsilon(1) = 1.

Antipode: S(x_{(n,l), k}) = -x_{(n,l), k}; extended as anti-homomorphism
(here, since algebra is commutative, anti-homomorphism = homomorphism).

Mellin slot k is INTRINSIC and preserved under Delta because each generator
is primitive. So Delta respects the k-grading bit-exactly.

Hopf axiom checks at n_max in {2, 3}
-------------------------------------
We verify bit-exactly:
  (a) Coassociativity (Delta (x) id) o Delta = (id (x) Delta) o Delta
      on every generator (both sides equal the standard "trio" symmetric
      sum x (x) 1 (x) 1 + 1 (x) x (x) 1 + 1 (x) 1 (x) x).
  (b) Counit (eps (x) id) o Delta = id = (id (x) eps) o Delta on every
      generator.
  (c) Antipode m o (S (x) id) o Delta = eta o eps = m o (id (x) S) o Delta
      on every generator.
  (d) Bialgebra compatibility: Delta is an algebra homomorphism (Delta(xy) =
      Delta(x) Delta(y)). Verify on small products x*y for distinct generators.
  (e) k-grading preservation under Delta (each tensor summand sits in the
      same k-grade as the original).

We check on:
  - Each generator x_{(n,l), k}.
  - Three quadratic products: x_{(1,0),0} * x_{(1,1),0}, x_{(1,0),0} * x_{(2,0),1},
    x_{(2,1),2} * x_{(2,2),2}.

Pro-system truncation P_{n+1 -> n} on Hopf algebra
---------------------------------------------------
P_{n+1 -> n} on generators: sends x_{(n', l'), k} -> x_{(n', l'), k} if
n' <= n, and to 0 if n' = n+1. Extended as algebra homomorphism.

Verify P_{n+1 -> n} is a Hopf algebra homomorphism by checking:
  (i) Delta o P = (P (x) P) o Delta on every generator
  (ii) eps o P = eps
  (iii) S o P = P o S
on n_max in {1->2, 2->3}.

Discipline
----------
Bit-exact sympy.Rational throughout. No floats. No PSLQ. No transcendentals.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Integer, Rational, Symbol, Poly

OUTPUT_DIR = Path(__file__).parent / "data"
OUTPUT_DIR.mkdir(exist_ok=True)
OUTPUT_PATH = OUTPUT_DIR / "sprint_q5p_stage2_hopf.json"


# =====================================================================
# Section 1: Underlying graded vector space + sector enumeration
# =====================================================================

def sectors_at(n_max: int) -> List[Tuple[int, int]]:
    """All (n, l) sectors at cutoff n_max."""
    return [(n, l) for n in range(1, n_max + 1) for l in range(0, n + 1)]


def N_sectors(n_max: int) -> int:
    """Closed form n_max*(n_max+3)/2."""
    return n_max * (n_max + 3) // 2


def generator_label(n: int, l: int, k: int) -> str:
    """Canonical name for the symbol x_{(n,l), k}."""
    return f"x[({n},{l}),{k}]"


def gen_list(n_max: int, ks: List[int] = (0, 1, 2)) -> List[Tuple[int, int, int]]:
    """List of generator labels at cutoff n_max."""
    return [(n, l, k) for (n, l) in sectors_at(n_max) for k in ks]


def grade(n: int, l: int, k: int) -> int:
    """Connected grading: degree = shell number n."""
    return n


# =====================================================================
# Section 2: Symbolic representation of H_GV as polynomial algebra
# =====================================================================
# We represent elements of H_GV as sympy polynomial expressions in the
# generators. Tensor products are represented as TUPLES of polynomials.
# Delta returns a list of tensor terms (lhs, rhs) with rational coefficients
# kept as multipliers of products.
#
# Internally, we use FROZEN MULTISETS of generator triples (n,l,k) to represent
# monomials. A monomial is a Counter-like dict mapping (n,l,k) -> exponent.
# An element of H_GV is a dict from monomial-tuple -> rational coefficient.

MonomialKey = Tuple[Tuple[Tuple[int, int, int], int], ...]   # sorted tuple of ((n,l,k), exp)
Element = Dict[MonomialKey, Rational]


def monomial_from_dict(d: Dict[Tuple[int, int, int], int]) -> MonomialKey:
    """Normalize a {gen: exp} dict into a sorted tuple key."""
    return tuple(sorted((g, e) for g, e in d.items() if e > 0))


def monomial_to_dict(key: MonomialKey) -> Dict[Tuple[int, int, int], int]:
    return {g: e for (g, e) in key}


UNIT_MONO: MonomialKey = ()   # empty product = 1


def elt_unit() -> Element:
    return {UNIT_MONO: Rational(1)}


def elt_zero() -> Element:
    return {}


def elt_gen(g: Tuple[int, int, int]) -> Element:
    """The single-generator element x_g."""
    return {monomial_from_dict({g: 1}): Rational(1)}


def elt_add(a: Element, b: Element) -> Element:
    """a + b."""
    out: Element = dict(a)
    for k, v in b.items():
        cur = out.get(k, Rational(0))
        new = cur + v
        if new == 0:
            out.pop(k, None)
        else:
            out[k] = new
    return out


def elt_scale(c: Rational, a: Element) -> Element:
    if c == 0:
        return {}
    return {k: c * v for k, v in a.items()}


def elt_neg(a: Element) -> Element:
    return elt_scale(Rational(-1), a)


def elt_mul(a: Element, b: Element) -> Element:
    """Polynomial multiplication (commutative algebra)."""
    out: Element = {}
    for ka, va in a.items():
        da = monomial_to_dict(ka)
        for kb, vb in b.items():
            db = monomial_to_dict(kb)
            merged: Dict[Tuple[int, int, int], int] = dict(da)
            for g, e in db.items():
                merged[g] = merged.get(g, 0) + e
            kmerged = monomial_from_dict(merged)
            cur = out.get(kmerged, Rational(0))
            new = cur + va * vb
            if new == 0:
                out.pop(kmerged, None)
            else:
                out[kmerged] = new
    return out


def elt_equals(a: Element, b: Element) -> bool:
    """Bit-exact equality."""
    if set(a.keys()) != set(b.keys()):
        return False
    for k in a.keys():
        if a[k] != b[k]:
            return False
    return True


def elt_pretty(a: Element) -> str:
    """Sortable readable form."""
    if not a:
        return "0"
    parts = []
    for k in sorted(a.keys()):
        coef = a[k]
        if not k:
            parts.append(str(coef))
        else:
            monstr = "*".join(f"x[{g}]^{e}" if e > 1 else f"x[{g}]" for (g, e) in k)
            if coef == 1:
                parts.append(monstr)
            elif coef == -1:
                parts.append(f"-{monstr}")
            else:
                parts.append(f"{coef}*{monstr}")
    return " + ".join(parts)


# =====================================================================
# Section 3: Tensor algebra (H_GV otimes H_GV) and (H_GV)^{(x)3}
# =====================================================================
# A 2-tensor is dict from (MonomialKey, MonomialKey) -> Rational.
# A 3-tensor is dict from (MonomialKey, MonomialKey, MonomialKey) -> Rational.

Tensor2 = Dict[Tuple[MonomialKey, MonomialKey], Rational]
Tensor3 = Dict[Tuple[MonomialKey, MonomialKey, MonomialKey], Rational]


def t2_add(a: Tensor2, b: Tensor2) -> Tensor2:
    out: Tensor2 = dict(a)
    for k, v in b.items():
        cur = out.get(k, Rational(0))
        new = cur + v
        if new == 0:
            out.pop(k, None)
        else:
            out[k] = new
    return out


def t2_zero() -> Tensor2:
    return {}


def t2_from_elt_pair(a: Element, b: Element) -> Tensor2:
    """a (x) b."""
    out: Tensor2 = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = (ka, kb)
            cur = out.get(key, Rational(0))
            new = cur + va * vb
            if new == 0:
                out.pop(key, None)
            else:
                out[key] = new
    return out


def t2_equals(a: Tensor2, b: Tensor2) -> bool:
    if set(a.keys()) != set(b.keys()):
        return False
    for k in a.keys():
        if a[k] != b[k]:
            return False
    return True


def t2_mul(a: Tensor2, b: Tensor2) -> Tensor2:
    """Componentwise multiplication on H (x) H: (a1 (x) a2)(b1 (x) b2) = a1 b1 (x) a2 b2."""
    out: Tensor2 = {}
    for (ka1, ka2), va in a.items():
        elt_a1: Element = {ka1: Rational(1)}
        elt_a2: Element = {ka2: Rational(1)}
        for (kb1, kb2), vb in b.items():
            elt_b1: Element = {kb1: Rational(1)}
            elt_b2: Element = {kb2: Rational(1)}
            prod1 = elt_mul(elt_a1, elt_b1)
            prod2 = elt_mul(elt_a2, elt_b2)
            coef = va * vb
            for k1, v1 in prod1.items():
                for k2, v2 in prod2.items():
                    key = (k1, k2)
                    cur = out.get(key, Rational(0))
                    new = cur + coef * v1 * v2
                    if new == 0:
                        out.pop(key, None)
                    else:
                        out[key] = new
    return out


def t3_add(a: Tensor3, b: Tensor3) -> Tensor3:
    out: Tensor3 = dict(a)
    for k, v in b.items():
        cur = out.get(k, Rational(0))
        new = cur + v
        if new == 0:
            out.pop(k, None)
        else:
            out[k] = new
    return out


def t3_equals(a: Tensor3, b: Tensor3) -> bool:
    if set(a.keys()) != set(b.keys()):
        return False
    for k in a.keys():
        if a[k] != b[k]:
            return False
    return True


# =====================================================================
# Section 4: Coproduct, counit, antipode
# =====================================================================
# We use the PRIMITIVE coproduct on generators: Delta(x_g) = x_g (x) 1 + 1 (x) x_g.
# Extended as an algebra homomorphism by Delta(a*b) = Delta(a) * Delta(b)
# in the tensor algebra H (x) H (with componentwise multiplication).


def delta_on_gen(g: Tuple[int, int, int]) -> Tensor2:
    """Delta(x_g) = x_g (x) 1 + 1 (x) x_g."""
    ka = monomial_from_dict({g: 1})
    out: Tensor2 = {}
    # x_g (x) 1
    out[(ka, UNIT_MONO)] = Rational(1)
    # 1 (x) x_g
    out[(UNIT_MONO, ka)] = Rational(1)
    return out


def delta_on_monomial(key: MonomialKey) -> Tensor2:
    """Delta on a monomial: Delta is an algebra homomorphism.

    Delta(x_g^e) = (x_g (x) 1 + 1 (x) x_g)^e via binomial.
    Delta(prod x_{g_i}^{e_i}) = prod Delta(x_{g_i})^{e_i} in H (x) H.
    """
    if not key:
        # Delta(1) = 1 (x) 1
        return {(UNIT_MONO, UNIT_MONO): Rational(1)}
    # build via repeated multiplication
    result: Tensor2 = {(UNIT_MONO, UNIT_MONO): Rational(1)}
    for (g, e) in key:
        dg = delta_on_gen(g)
        for _ in range(e):
            result = t2_mul(result, dg)
    return result


def delta(a: Element) -> Tensor2:
    """Apply Delta linearly to an element a."""
    out: Tensor2 = {}
    for k, v in a.items():
        dk = delta_on_monomial(k)
        for tk, tv in dk.items():
            cur = out.get(tk, Rational(0))
            new = cur + v * tv
            if new == 0:
                out.pop(tk, None)
            else:
                out[tk] = new
    return out


def epsilon_on_monomial(key: MonomialKey) -> Rational:
    """eps(1) = 1, eps(x_g^e ...) = 0 for any positive-degree monomial."""
    if not key:
        return Rational(1)
    return Rational(0)


def epsilon(a: Element) -> Rational:
    out = Rational(0)
    for k, v in a.items():
        out += v * epsilon_on_monomial(k)
    return out


def antipode_on_monomial(key: MonomialKey) -> Element:
    """S(1) = 1; S is an algebra anti-homomorphism with S(x_g) = -x_g.

    Since the algebra is commutative, anti-homomorphism = homomorphism.
    S(prod x_g^{e_g}) = prod (-x_g)^{e_g} = (-1)^{sum e_g} prod x_g^{e_g}.
    """
    if not key:
        return elt_unit()
    total_deg = sum(e for (_, e) in key)
    sign = Rational(1) if (total_deg % 2 == 0) else Rational(-1)
    return {key: sign}


def antipode(a: Element) -> Element:
    out: Element = {}
    for k, v in a.items():
        sk = antipode_on_monomial(k)
        for kk, vv in sk.items():
            cur = out.get(kk, Rational(0))
            new = cur + v * vv
            if new == 0:
                out.pop(kk, None)
            else:
                out[kk] = new
    return out


# =====================================================================
# Section 5: Axiom verification
# =====================================================================


def coassoc_lhs(g_or_mono: MonomialKey) -> Tensor3:
    """((Delta (x) id) o Delta) applied to monomial."""
    inner = delta_on_monomial(g_or_mono)
    out: Tensor3 = {}
    for (k1, k2), v in inner.items():
        d1 = delta_on_monomial(k1)
        for (k1a, k1b), w in d1.items():
            key3 = (k1a, k1b, k2)
            cur = out.get(key3, Rational(0))
            new = cur + v * w
            if new == 0:
                out.pop(key3, None)
            else:
                out[key3] = new
    return out


def coassoc_rhs(g_or_mono: MonomialKey) -> Tensor3:
    """((id (x) Delta) o Delta) applied to monomial."""
    inner = delta_on_monomial(g_or_mono)
    out: Tensor3 = {}
    for (k1, k2), v in inner.items():
        d2 = delta_on_monomial(k2)
        for (k2a, k2b), w in d2.items():
            key3 = (k1, k2a, k2b)
            cur = out.get(key3, Rational(0))
            new = cur + v * w
            if new == 0:
                out.pop(key3, None)
            else:
                out[key3] = new
    return out


def check_coassoc_on_monomial(key: MonomialKey) -> bool:
    return t3_equals(coassoc_lhs(key), coassoc_rhs(key))


def counit_compat_lhs(key: MonomialKey) -> Element:
    """((eps (x) id) o Delta)(x) -- collapse first tensor factor by eps."""
    inner = delta_on_monomial(key)
    out: Element = {}
    for (k1, k2), v in inner.items():
        coef = v * epsilon_on_monomial(k1)
        if coef == 0:
            continue
        cur = out.get(k2, Rational(0))
        new = cur + coef
        if new == 0:
            out.pop(k2, None)
        else:
            out[k2] = new
    return out


def counit_compat_rhs(key: MonomialKey) -> Element:
    """((id (x) eps) o Delta)(x)."""
    inner = delta_on_monomial(key)
    out: Element = {}
    for (k1, k2), v in inner.items():
        coef = v * epsilon_on_monomial(k2)
        if coef == 0:
            continue
        cur = out.get(k1, Rational(0))
        new = cur + coef
        if new == 0:
            out.pop(k1, None)
        else:
            out[k1] = new
    return out


def check_counit_on_monomial(key: MonomialKey) -> Tuple[bool, bool]:
    """Returns (left_ok, right_ok). Both should equal {key: 1}."""
    target: Element = {key: Rational(1)} if key else {UNIT_MONO: Rational(1)}
    left = counit_compat_lhs(key)
    right = counit_compat_rhs(key)
    return elt_equals(left, target), elt_equals(right, target)


def antipode_axiom_left(key: MonomialKey) -> Element:
    """m o (S (x) id) o Delta applied to monomial."""
    inner = delta_on_monomial(key)
    out: Element = {}
    for (k1, k2), v in inner.items():
        s1 = antipode_on_monomial(k1)
        e2: Element = {k2: Rational(1)}
        prod = elt_mul(s1, e2)
        for kk, vv in prod.items():
            cur = out.get(kk, Rational(0))
            new = cur + v * vv
            if new == 0:
                out.pop(kk, None)
            else:
                out[kk] = new
    return out


def antipode_axiom_right(key: MonomialKey) -> Element:
    """m o (id (x) S) o Delta applied to monomial."""
    inner = delta_on_monomial(key)
    out: Element = {}
    for (k1, k2), v in inner.items():
        e1: Element = {k1: Rational(1)}
        s2 = antipode_on_monomial(k2)
        prod = elt_mul(e1, s2)
        for kk, vv in prod.items():
            cur = out.get(kk, Rational(0))
            new = cur + v * vv
            if new == 0:
                out.pop(kk, None)
            else:
                out[kk] = new
    return out


def check_antipode_on_monomial(key: MonomialKey) -> Tuple[bool, bool]:
    """Verify m(S (x) id) Delta(x) = m(id (x) S) Delta(x) = eps(x) * 1.

    For 1: should give 1. For pos-degree monomial: should give 0.
    """
    if not key:
        target = elt_unit()
    else:
        target = elt_zero()
    left = antipode_axiom_left(key)
    right = antipode_axiom_right(key)
    return elt_equals(left, target), elt_equals(right, target)


def check_bialgebra_compat(g1: Tuple[int, int, int], g2: Tuple[int, int, int]) -> bool:
    """Verify Delta(x_{g1} x_{g2}) = Delta(x_{g1}) Delta(x_{g2}) bit-exactly."""
    elt_a = elt_gen(g1)
    elt_b = elt_gen(g2)
    elt_ab = elt_mul(elt_a, elt_b)
    lhs = delta(elt_ab)
    da = delta(elt_a)
    db = delta(elt_b)
    rhs = t2_mul(da, db)
    return t2_equals(lhs, rhs)


def check_k_grading_preservation(g: Tuple[int, int, int]) -> bool:
    """Verify that every tensor summand of Delta(x_{(n,l),k}) has its
    generators all carrying the same k label as g (when restricted to
    generators with that label)."""
    _, _, k = g
    d = delta_on_monomial(monomial_from_dict({g: 1}))
    for (ka, kb), _ in d.items():
        for kk in (ka, kb):
            for (gg, _e) in kk:
                _, _, kg = gg
                if kg != k:
                    return False
    return True


# =====================================================================
# Section 6: Pro-system truncation as Hopf homomorphism
# =====================================================================


def truncate_on_gen(g: Tuple[int, int, int], n_target: int) -> Element:
    """P_{n+1 -> n_target} on generator: keep if n <= n_target, else 0."""
    n, l, k = g
    if n <= n_target:
        return elt_gen(g)
    return elt_zero()


def truncate_on_monomial(key: MonomialKey, n_target: int) -> Element:
    """Extend P to monomials by algebra-homomorphism: if any factor is dropped
    (n > n_target), the whole monomial maps to 0."""
    out = elt_unit()
    for (g, e) in key:
        tg = truncate_on_gen(g, n_target)
        if not tg:
            return elt_zero()
        for _ in range(e):
            out = elt_mul(out, tg)
    return out


def truncate(a: Element, n_target: int) -> Element:
    out: Element = {}
    for k, v in a.items():
        tk = truncate_on_monomial(k, n_target)
        for kk, vv in tk.items():
            cur = out.get(kk, Rational(0))
            new = cur + v * vv
            if new == 0:
                out.pop(kk, None)
            else:
                out[kk] = new
    return out


def truncate_tensor2(t: Tensor2, n_target: int) -> Tensor2:
    """(P (x) P) applied to a 2-tensor."""
    out: Tensor2 = {}
    for (k1, k2), v in t.items():
        t1 = truncate_on_monomial(k1, n_target)
        t2 = truncate_on_monomial(k2, n_target)
        for kk1, vv1 in t1.items():
            for kk2, vv2 in t2.items():
                key = (kk1, kk2)
                cur = out.get(key, Rational(0))
                new = cur + v * vv1 * vv2
                if new == 0:
                    out.pop(key, None)
                else:
                    out[key] = new
    return out


def check_truncation_hopf_hom(g: Tuple[int, int, int], n_target: int) -> Dict:
    """Verify P_{n+1 -> n_target} is a Hopf homomorphism on generator g:
    (i) Delta o P = (P (x) P) o Delta
    (ii) eps o P = eps
    (iii) S o P = P o S
    """
    elt_g = elt_gen(g)
    # (i)
    P_g = truncate(elt_g, n_target)
    Delta_P_g = delta(P_g)
    Delta_g = delta(elt_g)
    P_P_Delta_g = truncate_tensor2(Delta_g, n_target)
    coproduct_ok = t2_equals(Delta_P_g, P_P_Delta_g)
    # (ii)
    eps_P_g = epsilon(P_g)
    eps_g = epsilon(elt_g)
    counit_ok = (eps_P_g == eps_g)
    # (iii)
    S_P_g = antipode(P_g)
    S_g = antipode(elt_g)
    P_S_g = truncate(S_g, n_target)
    antipode_ok = elt_equals(S_P_g, P_S_g)
    return {
        "generator": list(g),
        "n_target": n_target,
        "coproduct_compat": coproduct_ok,
        "counit_compat": counit_ok,
        "antipode_compat": antipode_ok,
        "all_ok": coproduct_ok and counit_ok and antipode_ok,
    }


# =====================================================================
# Section 7: Affine group scheme U*_{GeoVac} characterisation
# =====================================================================


def U_star_dim_grouplike(n_max: int) -> int:
    """At finite cutoff, U*^{(n_max)} = Spec(H_GV^{(n_max), *}).

    Since H_GV is a polynomial algebra Sym(V) (commutative, free-graded),
    the affine algebraic group Spec(H_GV*) on group-like elements is V^*
    viewed as an additive group, i.e. the affine space A^{dim V}.

    Equivalently, characters chi: H_GV -> Q (algebra maps with chi(1) = 1)
    are determined by their values on the generators (Q-linear functionals
    on V), giving Hom_{alg}(H_GV, Q) = V^* ~= Q^{dim V} as a set/scheme.

    The group structure on characters via convolution (chi*chi')(a) =
    (chi (x) chi') Delta(a): on primitive generators,
    (chi*chi')(x_g) = chi(x_g) + chi'(x_g),
    so the group is V^* under addition: the additive group G_a^{dim V}.
    """
    return 3 * N_sectors(n_max)


# =====================================================================
# MAIN
# =====================================================================

def main():
    t0 = time.time()

    results: Dict = {
        "sprint": "Q5'-Stage2-Hopf",
        "date": "2026-06-05",
        "stage": "Stage 2 first scoping step",
        "verdict": "POSITIVE (pending all-axiom verification)",
    }

    # ============= Section 1: Underlying space =============
    results["section_1_underlying_space"] = {}
    for n_max in [2, 3]:
        ks = [0, 1, 2]
        gens = gen_list(n_max, ks)
        results["section_1_underlying_space"][f"n_max_{n_max}"] = {
            "N_sectors": N_sectors(n_max),
            "n_k_slots": len(ks),
            "dim_V": len(gens),
            "dim_V_check_formula": 3 * N_sectors(n_max),
            "generators": [list(g) for g in gens],
            "grade_of_generator": {f"{list(g)}": grade(*g) for g in gens},
            "grading_homogeneous_levels": sorted({grade(*g) for g in gens}),
        }
        assert len(gens) == 3 * N_sectors(n_max)

    # ============= Section 2: Coproduct on generators =============
    results["section_2_coproduct_on_generators"] = {}
    for n_max in [2, 3]:
        gens = gen_list(n_max)
        ans = []
        for g in gens:
            d = delta_on_gen(g)
            terms = []
            for (k1, k2), v in sorted(d.items()):
                terms.append({
                    "lhs_mono": [list(x) + [e] for (x, e) in k1] if k1 else "1",
                    "rhs_mono": [list(x) + [e] for (x, e) in k2] if k2 else "1",
                    "coef": str(v),
                })
            ans.append({"generator": list(g), "delta_terms": terms})
        results["section_2_coproduct_on_generators"][f"n_max_{n_max}"] = ans[:5]  # truncate for size

    # ============= Section 3: Counit on generators =============
    results["section_3_counit_on_generators"] = {
        "rule": "eps(x_g) = 0 for any generator; eps(1) = 1",
        "verified_n_max_3_count": len(gen_list(3)),
    }
    for n_max in [2, 3]:
        gens = gen_list(n_max)
        all_zero = all(epsilon(elt_gen(g)) == 0 for g in gens)
        eps_unit_one = (epsilon(elt_unit()) == Rational(1))
        results["section_3_counit_on_generators"][f"n_max_{n_max}"] = {
            "epsilon_zero_on_generators": all_zero,
            "epsilon_one_on_unit": eps_unit_one,
        }

    # ============= Section 4: Antipode on generators =============
    results["section_4_antipode_on_generators"] = {
        "rule": "S(x_g) = -x_g on generators; S(1) = 1; extended as algebra homomorphism",
    }
    for n_max in [2, 3]:
        gens = gen_list(n_max)
        sample = []
        for g in gens[:5]:
            elt = elt_gen(g)
            s_elt = antipode(elt)
            # expect s_elt = -x_g
            target = elt_neg(elt)
            sample.append({
                "generator": list(g),
                "S_x_g_equals_neg_x_g": elt_equals(s_elt, target),
            })
        results["section_4_antipode_on_generators"][f"n_max_{n_max}"] = sample

    # ============= Section 5: Axiom verification at n_max in {2, 3} =============
    results["section_5_axiom_verification"] = {}

    for n_max in [2, 3]:
        cell: Dict = {}
        gens = gen_list(n_max)
        # Build a test set: all generators + 3 quadratic products + 1 cubic + 1 product across k
        # Restrict for speed at n_max=3 (which has 27 generators).
        test_monomials: List[MonomialKey] = []
        # All single-generator monomials
        for g in gens:
            test_monomials.append(monomial_from_dict({g: 1}))
        # Pairs (a few)
        if n_max == 2:
            pair_specs = [
                ((1, 0, 0), (1, 1, 0)),
                ((1, 0, 0), (2, 0, 1)),
                ((2, 1, 2), (2, 2, 2)),
                ((1, 0, 0), (1, 0, 0)),   # square
            ]
        else:
            pair_specs = [
                ((1, 0, 0), (3, 3, 0)),
                ((2, 1, 1), (3, 0, 2)),
                ((3, 2, 2), (3, 3, 2)),
                ((2, 0, 0), (2, 0, 0)),   # square
            ]
        for g1, g2 in pair_specs:
            test_monomials.append(monomial_from_dict({g1: 1, g2: 1}) if g1 != g2
                                  else monomial_from_dict({g1: 2}))
        # The unit
        test_monomials.append(UNIT_MONO)

        # (a) Coassociativity
        coassoc_results = []
        all_coassoc = True
        for mono in test_monomials:
            ok = check_coassoc_on_monomial(mono)
            if not ok:
                all_coassoc = False
            coassoc_results.append({
                "monomial_size": len(mono),
                "monomial": [[list(g), e] for (g, e) in mono] if mono else "1",
                "ok": ok,
            })
        cell["coassoc"] = {
            "all_ok": all_coassoc,
            "n_checks": len(coassoc_results),
            "sample_results": coassoc_results[:5],
        }

        # (b) Counit compatibility
        counit_results = []
        all_counit_left = True
        all_counit_right = True
        for mono in test_monomials:
            l_ok, r_ok = check_counit_on_monomial(mono)
            if not l_ok:
                all_counit_left = False
            if not r_ok:
                all_counit_right = False
            counit_results.append({
                "monomial": [[list(g), e] for (g, e) in mono] if mono else "1",
                "left_ok": l_ok,
                "right_ok": r_ok,
            })
        cell["counit_compat"] = {
            "all_left_ok": all_counit_left,
            "all_right_ok": all_counit_right,
            "n_checks": len(counit_results),
            "sample_results": counit_results[:5],
        }

        # (c) Antipode
        antipode_results = []
        all_antipode_left = True
        all_antipode_right = True
        for mono in test_monomials:
            l_ok, r_ok = check_antipode_on_monomial(mono)
            if not l_ok:
                all_antipode_left = False
            if not r_ok:
                all_antipode_right = False
            antipode_results.append({
                "monomial": [[list(g), e] for (g, e) in mono] if mono else "1",
                "left_ok": l_ok,
                "right_ok": r_ok,
            })
        cell["antipode"] = {
            "all_left_ok": all_antipode_left,
            "all_right_ok": all_antipode_right,
            "n_checks": len(antipode_results),
            "sample_results": antipode_results[:5],
        }

        # (d) Bialgebra compatibility: a sample of distinct-generator products
        bialg_results = []
        all_bialg = True
        sample_pairs: List[Tuple[Tuple[int, int, int], Tuple[int, int, int]]] = []
        if n_max == 2:
            sample_pairs = [
                ((1, 0, 0), (1, 1, 0)),
                ((1, 0, 0), (2, 0, 1)),
                ((2, 1, 2), (2, 2, 2)),
                ((1, 0, 0), (1, 0, 1)),   # same sector different k
                ((1, 0, 0), (2, 2, 2)),
            ]
        else:
            sample_pairs = [
                ((1, 0, 0), (3, 3, 0)),
                ((2, 1, 1), (3, 0, 2)),
                ((3, 2, 2), (3, 3, 2)),
                ((1, 0, 0), (1, 0, 1)),
            ]
        for g1, g2 in sample_pairs:
            ok = check_bialgebra_compat(g1, g2)
            if not ok:
                all_bialg = False
            bialg_results.append({"g1": list(g1), "g2": list(g2), "ok": ok})
        cell["bialgebra_compat"] = {
            "all_ok": all_bialg,
            "results": bialg_results,
        }

        # (e) k-grading preservation: each generator
        kgrade_results = []
        all_kgrade = True
        for g in gens:
            ok = check_k_grading_preservation(g)
            if not ok:
                all_kgrade = False
            kgrade_results.append({"generator": list(g), "ok": ok})
        cell["k_grading_preservation"] = {
            "all_ok": all_kgrade,
            "n_checks": len(kgrade_results),
        }

        results["section_5_axiom_verification"][f"n_max_{n_max}"] = cell

    # ============= Section 6: Pro-system truncation = Hopf homomorphism =============
    results["section_6_truncation_hopf_homomorphism"] = {}
    for (n_hi, n_lo) in [(2, 1), (3, 2)]:
        # At cutoff n_hi, restrict to n_lo via P_{n_hi -> n_lo}
        gens_hi = gen_list(n_hi)
        per_gen = []
        all_ok_dict = {"coproduct": True, "counit": True, "antipode": True}
        for g in gens_hi:
            res = check_truncation_hopf_hom(g, n_lo)
            per_gen.append(res)
            if not res["coproduct_compat"]:
                all_ok_dict["coproduct"] = False
            if not res["counit_compat"]:
                all_ok_dict["counit"] = False
            if not res["antipode_compat"]:
                all_ok_dict["antipode"] = False
        results["section_6_truncation_hopf_homomorphism"][f"P_{n_hi}_to_{n_lo}"] = {
            "n_hi": n_hi,
            "n_lo": n_lo,
            "n_generators_at_n_hi": len(gens_hi),
            "all_coproduct_compat": all_ok_dict["coproduct"],
            "all_counit_compat": all_ok_dict["counit"],
            "all_antipode_compat": all_ok_dict["antipode"],
            "sample_results": per_gen[:5],
        }

    # ============= Section 7: U^*_{GeoVac} characterisation =============
    results["section_7_U_star_GeoVac"] = {}
    for n_max in [2, 3]:
        dim_V = 3 * N_sectors(n_max)
        results["section_7_U_star_GeoVac"][f"n_max_{n_max}"] = {
            "dim_V": dim_V,
            "U_star_n_max_dim": dim_V,
            "U_star_n_max_structure": "G_a^{dim V} (additive affine group)",
            "explanation": (
                "H_GV^{(n_max)} = Sym(V) is the polynomial algebra on dim V = "
                f"{dim_V} primitive generators. Characters chi: H_GV -> Q "
                "(algebra maps with chi(1) = 1) are determined by their values "
                "on the dim V generators, giving Hom_alg(H_GV, Q) = V^* ~= Q^{dim V}. "
                "The group structure via convolution (chi*chi')(x_g) = chi(x_g) + "
                "chi'(x_g) on primitives makes this the additive group G_a^{dim V}."
            ),
        }

    # ============= Section 8: Connes-Marcolli analog dictionary =============
    results["section_8_connes_marcolli_dictionary"] = {
        "header": "Dictionary GeoVac <-> Connes-Marcolli (cosmic Galois on Connes-Kreimer)",
        "table": [
            {
                "geovac_object": "Sector idempotent e_{(n,l)} at cutoff n_max",
                "ck_analog": "Primitive 1PI Feynman graph in graded connected Hopf algebra H_CK",
            },
            {
                "geovac_object": "Mellin slot k in {0, 1, 2} (M1 / M3 / M2)",
                "ck_analog": "Feynman-rule decoration on graph (e.g. external-leg state)",
            },
            {
                "geovac_object": "Generator x_{(n,l), k} of H_GV",
                "ck_analog": "Generator [Gamma, decoration] of H_CK",
            },
            {
                "geovac_object": "Grading by shell number n",
                "ck_analog": "Grading by loop number L(Gamma) in H_CK",
            },
            {
                "geovac_object": "Primitive coproduct Delta(x) = x (x) 1 + 1 (x) x",
                "ck_analog": (
                    "Sub-graph coproduct Delta(Gamma) = sum_{gamma subset Gamma} "
                    "gamma (x) Gamma/gamma; primitive iff Gamma has no proper "
                    "sub-divergent subgraph (= primitive UV divergences)"
                ),
            },
            {
                "geovac_object": "Counit eps(x_g) = 0",
                "ck_analog": "Counit eps(Gamma) = 0 on non-empty graphs",
            },
            {
                "geovac_object": "Antipode S(x_g) = -x_g",
                "ck_analog": (
                    "Antipode S(Gamma) defined recursively by Bogoliubov R'-bar; "
                    "GeoVac case reduces to S(x) = -x because all generators are primitive"
                ),
            },
            {
                "geovac_object": "U^*_{GeoVac}^{(n_max)} = G_a^{3 N(n_max)} (additive)",
                "ck_analog": (
                    "U^*_{CK} = pro-affine group scheme of characters of H_CK; "
                    "U^*_{GeoVac} is the abelian additive sub-quotient acting on "
                    "the primitive-summand Hopf algebra"
                ),
            },
            {
                "geovac_object": "Pro-system truncation P_{n+1 -> n} is Hopf-hom",
                "ck_analog": (
                    "Inductive system of finite-loop truncations H_CK^{<=L} is Hopf-hom"
                ),
            },
            {
                "geovac_object": "M1/M2/M3 partition preserved by Delta (k-grading)",
                "ck_analog": (
                    "External-leg structure preserved by sub-graph extraction "
                    "(Feynman-rule decoration is graph-internal)"
                ),
            },
        ],
        "structural_observation": (
            "The GeoVac candidate H_GV is the ABELIAN PRIMITIVE-COSPAN of the Connes-"
            "Kreimer construction. CK's non-trivial content comes from the sub-graph "
            "coproduct (a graph can contain itself as a proper subgraph in many ways); "
            "GeoVac's sectors at this stage do NOT carry that sub-sector structure "
            "because the CH Fock decomposition is sector-DISJOINT (idempotents are "
            "orthogonal). This abelianness is the structural reason the Stage-2 "
            "motivic Galois group is the additive group G_a^{dim V} rather than a "
            "pro-unipotent group scheme. Adding sub-sector structure (e.g. nested "
            "Hopf-tower J*(S^3) factors) would give a non-trivial unipotent factor; "
            "this is the natural next step for Stage 2's full Tannakian construction."
        ),
    }

    # ============= Section 9: Headline summary =============
    # extract pass/fail across the panel
    h = {}
    for n_max in [2, 3]:
        cell = results["section_5_axiom_verification"][f"n_max_{n_max}"]
        h[f"n_max_{n_max}"] = {
            "coassoc_all_ok": cell["coassoc"]["all_ok"],
            "counit_left_all_ok": cell["counit_compat"]["all_left_ok"],
            "counit_right_all_ok": cell["counit_compat"]["all_right_ok"],
            "antipode_left_all_ok": cell["antipode"]["all_left_ok"],
            "antipode_right_all_ok": cell["antipode"]["all_right_ok"],
            "bialgebra_all_ok": cell["bialgebra_compat"]["all_ok"],
            "k_grading_all_ok": cell["k_grading_preservation"]["all_ok"],
        }
        # truncation check
        if n_max == 2:
            tcell = results["section_6_truncation_hopf_homomorphism"]["P_2_to_1"]
        else:
            tcell = results["section_6_truncation_hopf_homomorphism"]["P_3_to_2"]
        h[f"n_max_{n_max}"]["truncation_coproduct_ok"] = tcell["all_coproduct_compat"]
        h[f"n_max_{n_max}"]["truncation_counit_ok"] = tcell["all_counit_compat"]
        h[f"n_max_{n_max}"]["truncation_antipode_ok"] = tcell["all_antipode_compat"]
    results["section_9_headline"] = h

    # Final verdict
    all_pass = all(all(cell.values()) for cell in h.values())
    if all_pass:
        results["verdict"] = "POSITIVE: all Hopf axioms bit-exact at n_max in {2, 3}; pro-system truncation is a Hopf homomorphism."
    else:
        results["verdict"] = "BORDERLINE / STOP: some axiom fails (see section 5 details)."

    results["wall_time_seconds"] = time.time() - t0

    # ============= Write JSON =============
    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=str)

    # ============= Print summary to stdout =============
    print("=" * 76)
    print(f"Sprint Q5'-Stage2-Hopf — bit-exact Hopf-algebra axiom verification")
    print("=" * 76)
    print(f"Wall time: {results['wall_time_seconds']:.2f} s")
    print(f"Verdict: {results['verdict']}")
    print()
    print("Headline:")
    for nm_label, cell in h.items():
        print(f"  {nm_label}: " + ", ".join(f"{k}={v}" for k, v in cell.items()))
    print()
    print(f"Output: {OUTPUT_PATH}")
    print("=" * 76)


if __name__ == "__main__":
    main()
