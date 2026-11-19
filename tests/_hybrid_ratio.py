"""Helper for tests/test_paper58_hybrid_every_ratio.py -- the quotient count order(f) = d_E + k on PARAMETRIC
closed forms of four (AA|AB) hybrids, one positive parameter u per exponent region.

Own code path.  Nothing is imported from debug/: the parametric atom forms (built once, 2026-11-19, by running
geovac's `hybrid_closed_form_shell` algebra with symbolic exponents; normalisation dropped) are frozen in
tests/data/paper58_hybrid_ratio_2026_11_19.json, and every count below is RECOMPUTED here from those forms.  The
frozen forms are tied back to the live evaluator by the test file (exact radial_norm product at fresh u).

Atoms (Lemma GI labels): e^{cR} (kind '1'), e^{cR} ln R ('log'), e^{cR} Ei(kappa R) ('Ei'; E1(kR) = -Ei(-kR)),
c and kappa linear expressions in u.

Proposition (quotient count; hand proof, cold-reviewed 2026-11-19).  With W the free C(R)-module on the
d/dR-closed atom set and W_exp its exponential atoms,
    order(f) = d_E + k,
  d_E = sum over E-shifts c of dim_C span{q_j : q_j the coefficient of an E-atom (Ei or ln R) with shift c},
  k   = number of exponential rates with a NONZERO coefficient in L f, L the order-d_E operator (Wronskian
        determinant) whose kernel is span{q_j e^{c_j R}} over a basis of those spans.
`quotient_count` implements exactly this, for any field of coefficients carrying R as a generator; the same code
is run on the review's eight toys (against a direct rank over Q(R)) and on the frozen hybrid forms.
"""
from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Dict, List, Tuple

import mpmath as mp
import sympy as sp
from sympy.polys.fields import field

DATA_PATH = Path(__file__).parent / "data" / "paper58_hybrid_ratio_2026_11_19.json"
U = sp.Symbol("u", positive=True)
RS = sp.Symbol("R", positive=True)

ORB = {"1s": (1, 0, 0), "2s": (2, 0, 0), "2p0": (2, 1, 0), "2p1": (2, 1, 1), "3d0": (3, 2, 0)}
CASES = {"fam1": ("2p0", "2p0", "1s", "1s"), "H2": ("2p1", "2p1", "2p1", "2p1"),
         "H4": ("2p0", "2p0", "2p0", "3d0"), "H7": ("2p0", "3d0", "2s", "1s")}
REGIONS = ("lt", "mid", "gt", "eq1", "eq2")
OPEN = ("lt", "mid", "gt")

Label = Tuple[sp.Expr, str, sp.Expr]          # (c, kind, kappa), kind in {'1', 'log', 'Ei'}


def load_data() -> dict:
    return json.loads(DATA_PATH.read_text())


def _sym(s: str) -> sp.Expr:
    return sp.expand(sp.sympify(s, locals={"u": U, "R": RS}))


# =============================================================== parse a frozen form into a field
class Form:
    """f = sum over labels of coeff * atom, coeff in K = Q(u, R, C_0, ...) (C_i = integer-power constants:
    Euler gamma, pi, logs).  A fractional-power constant (sqrt 3, ...) must multiply EVERY term and is dropped as
    a common nonzero factor (asserted)."""

    def __init__(self, entry: dict):
        self.info = {k: _sym(v) for k, v in entry["info"].items()}
        raw = []
        for ks, cs in entry["form"].items():
            c, kind, kap, consts = ast.literal_eval(ks)
            raw.append((_sym(c), kind, _sym(kap), [(n, sp.Rational(e)) for n, e in consts], cs))
        frac_sets = [frozenset((n, e) for n, e in consts if e.q != 1) for *_x, consts, _cs in raw]
        assert len(set(frac_sets)) == 1, "a fractional constant that is not a common factor"
        self.common_fractional = sorted(frac_sets[0])
        names = sorted({n for *_x, consts, _cs in raw for n, e in consts if e.q == 1})
        self.const_names = names
        self.csyms = [sp.Symbol(f"C{i}") for i in range(len(names))]
        self.K, *gens = field([U, RS] + self.csyms, sp.QQ)
        self.UG, self.RG = gens[0], gens[1]
        cg = dict(zip(names, gens[2:]))
        self.coords: Dict[Label, object] = {}
        for c, kind, kap, consts, cs in raw:
            m = self.K.one
            for n, e in consts:
                if e.q == 1:
                    m *= cg[n] ** int(e)
            val = self.K.from_expr(sp.sympify(cs, locals={"u": U, "R": RS})) * m
            lab = (c, kind, kap)
            self.coords[lab] = self.coords.get(lab, self.K.zero) + val
        self.coords = {k: v for k, v in self.coords.items() if v != 0}

    def const_free(self, v) -> bool:
        """The field element involves no constant generator (numerator and denominator)."""
        nvars = 2
        for poly in (v.numer, v.denom):
            for mon in poly.monoms():
                if any(mon[nvars:]):
                    return False
        return True


# =============================================================== generic calculus on {label: field element}
def derivative(coords: Dict, K, RG) -> Dict:
    out: Dict = {}

    def add(key, v):
        out[key] = out.get(key, K.zero) + v

    for (c, kind, kap), p in coords.items():
        add((c, kind, kap), p.diff(RG) + K.from_expr(c) * p)
        if kind == "log":
            add((c, "1", sp.Integer(0)), p / RG)
        elif kind == "Ei":
            add((sp.expand(c + kap), "1", sp.Integer(0)), p / RG)
    return {k: v for k, v in out.items() if v != 0}


def twisted(q, c, j: int, K, RG):
    """e^{-cR} d^j/dR^j (q e^{cR}) = (d/dR + c)^j q."""
    cc = K.from_expr(c)
    for _ in range(j):
        q = q.diff(RG) + cc * q
    return q


def det(rows: List[List[object]], K):
    n = len(rows)
    if n == 0:
        return K.one
    if n == 1:
        return rows[0][0]
    tot = K.zero
    for j in range(n):
        if rows[0][j] == 0:
            continue
        minor = [r[:j] + r[j + 1:] for r in rows[1:]]
        tot += (-1) ** j * rows[0][j] * det(minor, K)
    return tot


def span_basis(qs: List[object], K, RG) -> List[object]:
    """Greedy basis of span_C{q_j}: a q joins iff the Wronskian of the chosen set with it is nonzero."""
    chosen: List[object] = []
    for q in qs:
        trial = chosen + [q]
        m = len(trial)
        W = [[twisted(t, 0, j, K, RG) for j in range(m)] for t in trial]
        if det(W, K) != 0:
            chosen = trial
    return chosen


def quotient_count(coords: Dict, K, RG) -> dict:
    """d_E, L's coefficients, L f and k, by the Proposition.  Raises if L does not kill the E-part."""
    E = {lb: v for lb, v in coords.items() if lb[1] != "1"}
    shifts = sorted({lb[0] for lb in E}, key=str)
    basis: List[Tuple[sp.Expr, object]] = []
    per_shift = {}
    for c in shifts:
        qs = [v for lb, v in sorted(E.items(), key=lambda t: str(t[0])) if lb[0] == c]
        b = span_basis(qs, K, RG)
        per_shift[c] = len(b)
        basis += [(c, q) for q in b]
    d = len(basis)
    ds = [coords]
    for _ in range(d):
        ds.append(derivative(ds[-1], K, RG))
    A = []
    for j in range(d + 1):
        cols = [jj for jj in range(d + 1) if jj != j]
        A.append((-1) ** (d + j) * det([[twisted(q, c, jj, K, RG) for jj in cols] for (c, q) in basis], K))
    Lf: Dict = {}
    for j in range(d + 1):
        for key, v in ds[j].items():
            Lf[key] = Lf.get(key, K.zero) + A[j] * v
    Lf = {k: v for k, v in Lf.items() if v != 0}
    leftover = [k for k in Lf if k[1] != "1"]
    assert not leftover, f"L did not kill the E-part: {leftover}"
    rates = sorted({lb[0] for lb in Lf}, key=str)
    return {"d_E": d, "per_shift": per_shift, "A": A, "Lf": Lf, "k": len(rates), "Lf_rates": rates,
            "order": d + len(rates)}


# =============================================================== zero tests in u
def u_gcd(v, K) -> sp.Expr:
    """gcd over Q[u] of v's numerator coefficients w.r.t. R and the constants (v in K, lowest terms)."""
    groups: Dict[tuple, sp.Expr] = {}
    for mon, cf in v.numer.terms():
        groups[mon[1:]] = groups.get(mon[1:], sp.Integer(0)) + sp.Rational(cf) * U ** mon[0]
    g = sp.Integer(0)
    for p in groups.values():
        p = sp.expand(p)
        g = p if g == 0 else sp.gcd(g, p)
    return sp.factor(g)


def u_denominator_factors(v) -> List[sp.Expr]:
    den = v.denom.as_expr()
    return [fa for fa, _e in sp.factor_list(den)[1] if fa.has(U)]


def positive_roots(p: sp.Expr) -> List[sp.Expr]:
    """Real roots u > 0 of a nonzero polynomial in u (exact).  A nonzero constant has none."""
    p = sp.expand(p)
    assert p != 0, "zero polynomial"
    if not p.has(U):
        return []
    P = sp.Poly(p, U)
    return [r for r in P.real_roots() if r > 0]


def lin_vanishes_for_some_u(e: sp.Expr) -> bool:
    """True iff the (polynomial in u) expression e is zero at some u > 0, or identically."""
    e = sp.expand(e)
    if e == 0:
        return True
    return bool(positive_roots(e)) if e.has(U) else False


# =============================================================== toys (the cold review's eight)
def toy_coords(terms, K, RG) -> Dict:
    """terms: (coeff in R, (c, kind, kappa)), kind '1' | 'y' (ln R if kappa = 0, else Ei(kappa R))."""
    out: Dict = {}
    for cf, (c, kind, k) in terms:
        k = sp.Rational(k)
        kk = "1" if kind == "1" else ("log" if k == 0 else "Ei")
        lab = (sp.Rational(c), kk, k)
        out[lab] = out.get(lab, K.zero) + K.from_expr(sp.sympify(cf))
    return {k: v for k, v in out.items() if v != 0}


def toys():
    R = RS
    q = R ** 2 + 1
    qm = q.subs(R, -R)
    return {
        "T1 generic": [(q, (-1, "y", 1)), (-q, (-1, "y", 2)), (R, (-1, "1", 0)), (1, (-3, "1", 0))],
        "T2 r in kernel": [(q, (-1, "y", 1)), (-q, (-1, "y", 2)), (5 * q, (-1, "1", 0)), (1, (-3, "1", 0))],
        "T3 closure-rate exp": [(q, (-1, "y", 1)), (-q, (-1, "y", 2)), (R ** 3, (0, "1", 0))],
        "T4 independent at one shift": [(q, (-1, "y", 1)), (R, (-1, "y", 2)), (1, (-3, "1", 0))],
        "T5 dependent triple": [(q, (-1, "y", 1)), (R, (-1, "y", 2)), (q + R, (-1, "y", 3)), (R, (-1, "1", 0))],
        "T6 mirror pair": [(q, (-1, "y", sp.Rational(-1, 2))), (-q, (-1, "y", sp.Rational(1, 3))),
                           (qm, (1, "y", -3)), (-qm, (1, "y", -sp.Rational(5, 2))),
                           (q, (-1, "1", 0)), (R + 2, (-2, "1", 0)), (1 / R, (-sp.Rational(3, 2), "1", 0))],
        "T7 log + Ei": [(q, (-1, "y", 0)), (-q, (-1, "y", 2)), (R, (-1, "1", 0)), (1, (1, "1", 0))],
        "T8 near-miss": [(q, (-1, "y", 1)), (-q / R, (0, "1", 0)), (-q, (-1, "y", 2))],
    }


# =============================================================== numerical value of a frozen form at u0
def const_value(name: str, e: sp.Rational, u0) -> mp.mpf:
    if name == "PI":
        return mp.pi ** (mp.mpf(e.p) / e.q)
    if name == "GAM":
        return mp.euler ** int(e)
    if name.startswith("P:"):
        return mp.mpf(int(name[2:])) ** (mp.mpf(e.p) / e.q)
    if name.startswith("LOG:"):
        arg = sp.Rational(_sym(name[4:]).subs(U, u0))
        assert arg > 0, (name, arg)
        return mp.log(mp.mpf(arg.p) / arg.q) ** int(e)
    raise ValueError(name)


def form_value(entry: dict, u0, Rv, dps: int = 80) -> mp.mpf:
    """The frozen (normalisation-free) form at u = u0, R = Rv, from its JSON entry (all constants kept)."""
    with mp.workdps(dps):
        x = mp.mpf(sp.Rational(Rv).p) / sp.Rational(Rv).q
        tot = mp.mpf(0)
        for ks, cs in entry["form"].items():
            c, kind, kap, consts = ast.literal_eval(ks)
            cv = sp.Rational(_sym(c).subs(U, u0))
            kv = sp.Rational(_sym(kap).subs(U, u0))
            coef = sp.sympify(cs, locals={"u": U, "R": RS}).subs({U: u0, RS: sp.Rational(Rv)})
            assert coef.is_Rational, coef
            val = mp.mpf(int(coef.p)) / int(coef.q)
            for n, e in consts:
                val *= const_value(n, sp.Rational(e), u0)
            atom = mp.exp(mp.mpf(cv.p) / cv.q * x)
            if kind == "log":
                atom *= mp.log(x)
            elif kind == "Ei":
                assert kv != 0
                atom *= mp.ei(mp.mpf(kv.p) / kv.q * x)
            tot += val * atom
        return tot


__all__ = ["load_data", "Form", "derivative", "twisted", "quotient_count", "span_basis", "u_gcd",
           "u_denominator_factors", "positive_roots", "lin_vanishes_for_some_u", "toy_coords", "toys",
           "form_value", "CASES", "ORB", "REGIONS", "OPEN", "U", "RS"]
