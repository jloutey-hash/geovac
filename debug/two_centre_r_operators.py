"""Two-centre Slater-1s integrals as functions of R: minimal linear ODEs and singular points.

Chemistry period inventory (docs/chemistry_period_inventory.md), rows 2, 3a, 3b, 3c,
criteria 2 (governing equation) and 3 (singular points).

Method (exact, no floating-point data enters the operator):
  1. Closed form f(R) from geovac (qfd_core one-electron S/T/V and exchange;
     two_center_eri aabb / hybrid).
  2. f is decomposed as  sum_atoms  a(R; consts) * e^{cR} * S,  S in {1, ln R, E1(kR), Ei(kR)},
     a rational in R and polynomial in the constants {gamma_E, ln 2, ln 3, ...}.
     The derivative closes on this atom set (holonomic closure), so f^(j) are vectors
     over Q(R) indexed by (atom, constant-monomial).
  3. The first r with a Q(R)-linear dependency among f, ..., f^(r) gives the minimal-order
     operator with coefficients in Q[R] (constants split, i.e. gamma_E, ln 2 treated as
     independent of Q(R); minimality also relies on the standard independence of the atoms
     over Q(R) -- distinct exponentials, E1's at distinct rates, ln R).  The order over
     K(R) with the constants kept as field elements is reported alongside.
  4. Verification: sympy applies the operator to the closed form itself (sympy's own
     derivatives of expint/Ei/log) and simplifies to zero; plus a 60-digit mpmath check.
  5. Singular points: R = 0 (Fuchs test + indicial roots + whether f's own local
     expansion carries ln R), finite nonzero roots of the leading coefficient (exponents;
     apparent iff they are distinct non-negative integers -- guaranteed here because every
     Galois conjugate of f is a combination of the atoms, analytic on C minus {0}),
     R = infinity (irregular; exponential rates from the top-degree characteristic
     polynomial, compared with the atom rates c and c - k).

Writes debug/data/two_centre_r_operators.json.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from geovac import qfd_core as Q  # noqa: E402
from geovac.two_center_eri import R_s, aabb_closed_form, hybrid_closed_form  # noqa: E402

R = R_s
QR = sp.QQ.frac_field(R)
S1 = (1, 0, 0)

# ----------------------------------------------------------------- decomposition


def _const_symbols(expr):
    """Map non-rational constants (gamma_E, log of numbers, surds) to symbols."""
    reps = {}
    k = 0
    for a in sorted(expr.atoms(sp.log, sp.Pow, sp.EulerGamma.__class__), key=str):
        if a.free_symbols:
            continue
        if isinstance(a, sp.log) or a == sp.EulerGamma or (
                isinstance(a, sp.Pow) and not a.exp.is_integer):
            reps[a] = sp.Symbol(f"k{k}")
            k += 1
    if expr.has(sp.EulerGamma):
        reps[sp.EulerGamma] = sp.Symbol("gE")
    return reps


def decompose(expr):
    """f -> ({(c, kind, k): coeff}, const_symbol_map).  kind in '1','log','E1','Ei'."""
    expr = sp.expand(sp.expand_log(sp.expand(expr), force=True))
    # pull a global surd / constant factor shared by every term (normalization)
    reps = _const_symbols(expr)
    expr = sp.expand(expr.xreplace(reps))
    atoms: dict = {}
    for term in sp.Add.make_args(expr):
        c, kind, k, coeff = sp.Integer(0), "1", sp.Integer(0), sp.Integer(1)
        for fct in sp.Mul.make_args(term):
            if isinstance(fct, sp.exp):
                c += sp.diff(fct.args[0], R)
                assert sp.diff(fct.args[0], R, 2) == 0
                assert fct.args[0].subs(R, 0) == 0, fct
            elif isinstance(fct, sp.Pow) and isinstance(fct.base, sp.exp):
                arg = fct.base.args[0] * fct.exp
                c += sp.diff(arg, R)
            elif isinstance(fct, sp.expint):
                assert fct.args[0] == 1 and kind == "1", term
                kind, k = "E1", sp.diff(fct.args[1], R)
            elif isinstance(fct, sp.Ei):
                assert kind == "1", term
                kind, k = "Ei", sp.diff(fct.args[0], R)
            elif fct == sp.log(R):
                assert kind == "1"
                kind = "log"
            elif isinstance(fct, sp.Pow) and fct.base == sp.log(R):
                raise ValueError("ln^2 R")
            else:
                assert not fct.has(sp.exp, sp.expint, sp.Ei, sp.log), fct
                coeff *= fct
        key = (sp.nsimplify(c), kind, sp.nsimplify(k))
        atoms[key] = atoms.get(key, 0) + coeff
    atoms = {a: sp.together(v) for a, v in atoms.items()}
    atoms = {a: v for a, v in atoms.items() if sp.simplify(v) != 0}
    return atoms, reps


FR = sp.QQ.frac_field(R)
RG = FR.gens[0]


def to_field(atoms, reps, mode="split"):
    """Atoms -> {(atom, const-monomial): Q(R) field element}.

    mode='split': constants split off as independent (coefficients over Q(R));
    mode='spec':  constants specialised to fixed rationals (generic K(R) test).
    Splitting commutes with d/dR because the constants are constants.
    """
    csyms = sorted(set(reps.values()), key=str)
    cvals = {s: sp.Rational(7 + 3 * i, 13 + 5 * i) for i, s in enumerate(csyms)}
    out = {}
    for a, v in atoms.items():
        if mode == "spec" or not csyms:
            w = sp.together(v.subs(cvals)) if csyms else sp.together(v)
            parts = {(): w}
        else:
            num, den = sp.fraction(sp.together(v))
            assert not den.has(*csyms)
            P = sp.Poly(sp.expand(num), *csyms)
            parts = {m: cf / den for m, cf in P.terms()}
        for m, w in parts.items():
            el = FR.from_sympy(sp.together(w))
            if el != 0:
                out[(a, m)] = out.get((a, m), FR.zero) + el
    return out


def d_atoms(atoms):
    """d/dR on {((c, kind, k), mono): Q(R)} -- holonomic closure of the atom set."""
    out: dict = {}

    def add(key, v):
        out[key] = out.get(key, FR.zero) + v

    for ((c, kind, k), m), a in atoms.items():
        cq = FR.from_sympy(c)
        add(((c, kind, k), m), a.diff(RG) + cq * a)
        if kind == "log":
            add(((c, "1", 0), m), a / RG)
        elif kind == "E1":
            add(((c - k, "1", 0), m), -a / RG)
        elif kind == "Ei":
            add(((c + k, "1", 0), m), a / RG)
    return {key: v for key, v in out.items() if v != 0}


def _rank_at(derivs, x0):
    keys = sorted(set().union(*[set(D) for D in derivs]), key=str)
    rows = []
    for kk in keys:
        row = []
        for D in derivs:
            v = D.get(kk, FR.zero)
            row.append(sp.QQ.from_sympy(sp.Rational(FR.to_sympy(v).subs(R, x0))))
        rows.append(row)
    return DomainMatrix(rows, (len(keys), len(derivs)), sp.QQ).rank()


def min_order(atoms, reps, split_consts=True, rmax=16):
    """Smallest r with f..f^(r) dependent over Q(R) (generic rank at two rational points)."""
    derivs = [to_field(atoms, reps, "split" if split_consts else "spec")]
    for r in range(1, rmax + 1):
        derivs.append(d_atoms(derivs[-1]))
        if all(_rank_at(derivs, x0) < r + 1
               for x0 in (sp.Rational(37, 59), sp.Rational(101, 43))):
            return r, derivs
    raise RuntimeError("no operator up to rmax")


def minimal_operator(atoms, reps, dmax=40):
    """Minimal-order operator over Q(R), minimal coefficient degree, exact ansatz."""
    r, derivs = min_order(atoms, reps, split_consts=True)
    keys = sorted(set().union(*[set(D) for D in derivs]), key=str)
    eqrows = []                        # per key: list of Poly P_j (denominators cleared)
    for kk in keys:
        fr = [(sp.Poly(FR.numer(D.get(kk, FR.zero)).as_expr(), R),
               sp.Poly(FR.denom(D.get(kk, FR.zero)).as_expr(), R)) for D in derivs]
        L = fr[0][1]
        for _n, dd in fr[1:]:
            L = L.lcm(dd)
        eqrows.append([n * L.exquo(dd) for n, dd in fr])
    for d in range(0, dmax + 1):
        nunk = (r + 1) * (d + 1)
        rows = []
        for Ps in eqrows:
            emax = max(P.degree() for P in Ps) + d
            for e in range(emax + 1):
                row = [sp.QQ(0)] * nunk
                nz = False
                for j, P in enumerate(Ps):
                    for m in range(d + 1):
                        c = P.coeff_monomial(R ** (e - m)) if e - m >= 0 else 0
                        if c != 0:
                            row[j * (d + 1) + m] = sp.QQ.from_sympy(c)
                            nz = True
                if nz:
                    rows.append(row)
        M = DomainMatrix(rows, (len(rows), nunk), sp.QQ)
        ns = M.nullspace()
        if ns.shape[0] > 0:
            v = ns.to_Matrix().row(0)
            polys = [sp.expand(sum(v[j * (d + 1) + m] * R ** m for m in range(d + 1)))
                     for j in range(r + 1)]
            g = sp.gcd_list([p for p in polys if p != 0])
            polys = [sp.expand(sp.cancel(p / g)) for p in polys]
            lc = sp.Poly(polys[-1], R).LC()
            polys = [sp.expand(p / lc) for p in polys]
            m = sp.lcm([sp.fraction(c)[1] for p in polys for c in sp.Poly(p, R).coeffs()])
            polys = [sp.expand(p * m) for p in polys]
            return r, polys, int(ns.shape[0])
    raise RuntimeError("degree cap")


# ----------------------------------------------------------------- verification


def verify(expr, polys):
    res = sum(p * sp.diff(expr, R, j) for j, p in enumerate(polys))
    res = sp.expand(sp.expand_func(res))
    # symbolic: sympy's own derivatives of exp/expint/Ei/log, collected by atom; every
    # rational coefficient must cancel identically (const symbols kept symbolic)
    ra, _ = decompose(res)
    sym = len(ra) == 0
    mp.mp.dps = 60
    f = sp.lambdify(R, res, "mpmath")
    vals = [abs(f(mp.mpf(x))) for x in ("0.7", "1.4", "3.3")]
    scale = [abs(sp.lambdify(R, expr, "mpmath")(mp.mpf(x))) for x in ("0.7", "1.4", "3.3")]
    return sym, float(max(v / (s if s else 1) for v, s in zip(vals, scale)))


# ----------------------------------------------------------------- singular points


def _ord0(p):
    P = sp.Poly(p, R)
    return min(m[0] for m in P.monoms()) if not P.is_zero else None


def indicial_at_zero(polys):
    r = len(polys) - 1
    th = sp.Symbol("theta")
    vals = [(_ord0(p) - k) if p != 0 else None for k, p in enumerate(polys)]
    v = min(x for x in vals if x is not None)
    fuchs = vals[r] == v
    ind = 0
    for k, p in enumerate(polys):
        if vals[k] == v:
            low = sp.Poly(p, R).coeff_monomial(R ** _ord0(p))
            ind += low * sp.ff(th, k)
    ind = sp.expand(ind)
    roots = sp.roots(sp.Poly(ind, th)) if ind != 0 else {}
    return fuchs, {str(k): int(m) for k, m in roots.items()}, sp.factor(ind)


def finite_singular(polys):
    r = len(polys) - 1
    lead = sp.Poly(polys[-1], R)
    out = []
    for fac, mult in sp.factor_list(lead.as_expr())[1]:
        if sp.Poly(fac, R).degree() == 0:
            continue
        if sp.Poly(fac, R).as_expr() == R:
            continue
        rec = {"factor": str(fac), "multiplicity": int(mult),
               "roots": [str(sp.N(x, 12)) for x in sp.Poly(fac, R).nroots(n=15)]}
        if mult == 1:
            g = sp.Poly(fac, R)
            dlead = sp.Poly(sp.diff(lead.as_expr(), R), R)
            inv = sp.invert(dlead.as_expr(), fac, R)
            val = sp.rem(sp.expand(-sp.Poly(polys[r - 1], R).as_expr() * inv), fac, R)
            e_last = sp.simplify(r - 1 + val)
            rec["exponents"] = list(range(r - 1)) + [str(e_last)]
            rec["apparent"] = bool(e_last.is_Integer and e_last >= r - 1)
        else:
            # numeric indicial at each root, Frobenius in t = R - rho
            t = sp.Symbol("t")
            exps = []
            for rho in sp.Poly(fac, R).nroots(n=30):
                sh = [sp.expand(p.subs(R, rho + t)) for p in polys]
                vals = []
                for k, p in enumerate(sh):
                    P = sp.Poly(sp.N(p, 30), t)
                    terms = [(m[0], c) for m, c in P.terms() if abs(c) > 1e-18]
                    o = min(m for m, _ in terms) if terms else None
                    vals.append((o - k) if o is not None else None)
                v = min(x for x in vals if x is not None)
                th = sp.Symbol("theta")
                ind = 0
                for k, p in enumerate(sh):
                    if vals[k] == v:
                        P = sp.Poly(sp.N(p, 30), t)
                        ind += P.coeff_monomial(t ** (v + k)) * sp.ff(th, k)
                exps.append([str(sp.N(x, 10)) for x in sp.Poly(sp.expand(ind), th).nroots()])
            rec["exponents_numeric"] = exps
        out.append(rec)
    return out


def rates_at_infinity(polys):
    degs = [sp.Poly(p, R).degree() if p != 0 else -10 ** 9 for p in polys]
    D = max(degs)
    lam = sp.Symbol("lam")
    chi = sum(sp.Poly(p, R).LC() * lam ** k for k, p in enumerate(polys) if degs[k] == D)
    roots = sp.roots(sp.Poly(chi, lam))
    return {str(k): int(m) for k, m in roots.items()}, sp.factor(chi)


def atom_rates(atoms):
    rs = set()
    for (c, kind, k) in atoms:
        rs.add(c)
        if kind == "E1":
            rs.add(c - k)
        if kind == "Ei":
            rs.add(c + k)
    return sorted(rs)


def local_log_part(atoms, reps):
    """Coefficient B(R) of ln R in f's expansion at R = 0 (elementary)."""
    B = 0
    for (c, kind, k), a in atoms.items():
        if kind == "log":
            B += a * sp.exp(c * R)
        elif kind == "E1":
            B -= a * sp.exp(c * R)
        elif kind == "Ei":
            B += a * sp.exp(c * R)
    return sp.simplify(B)


def local_series(expr, reps, n=10):
    e = expr
    try:
        s = sp.series(e, R, 0, n).removeO()
    except Exception as exc:  # pragma: no cover
        return str(exc)
    return s


# ----------------------------------------------------------------- panel


def A(Z):
    return ("A", Fraction(Z), 1)


def B(Z):
    return ("B", Fraction(Z), 1)


PAIRS = [(1, 1), (1, 2), (2, 1), (Fraction(3, 2), 1)]


def build_panel():
    cases = []
    for za, zb in PAIRS:
        tag = f"zA={za},zB={zb}"
        cases.append(("S", "2", tag, lambda za=za, zb=zb: Q.overlap(A(za), B(zb), R)))
        cases.append(("T", "2", tag, lambda za=za, zb=zb: Q.kinetic(A(za), B(zb), R)))
        cases.append(("V <a|1/r_A|b>", "2", tag,
                      lambda za=za, zb=zb: Q._inv_r(A(za), B(zb), "A", R)))
        cases.append(("V <a|1/r_B|a>", "2", tag,
                      lambda za=za, zb=zb: Q._inv_r(A(za), A(za), "B", R)))
        cases.append(("(AA|BB) 1s", "3a", tag,
                      lambda za=za, zb=zb: aabb_closed_form(Fraction(za), S1, S1,
                                                            Fraction(zb), S1, S1)))
        rel = "eq" if za == zb else ("a_d>a_c" if zb > za else "a_d<a_c")
        cases.append((f"(AA|AB) 1s [{rel}]", "3b", tag,
                      lambda za=za, zb=zb: hybrid_closed_form(Fraction(za), S1, S1, S1,
                                                              Fraction(zb), S1)))
    # hybrid with a p-type one-centre pair: the class that carries E1 / Ei / gamma+ln R
    p0 = (2, 1, 0)
    for za, zb, rel in ((3, 3, "eq"), (3, 4, "a_d>a_c"), (3, 2, "a_d<a_c")):
        cases.append((f"(2p0 2p0|1s 1s_B) [{rel}]", "3b", f"ZA={za},ZB={zb}",
                      lambda za=za, zb=zb: hybrid_closed_form(Fraction(za), p0, p0, S1,
                                                              Fraction(zb), S1)))
    for z in (1, Fraction(3, 2), 2):
        cases.append(("(AB|AB) 1s", "3c", f"zA=zB={z}",
                      lambda z=z: Q.exchange_closed_form(Fraction(z), S1, S1, Fraction(z),
                                                         S1, S1, R, tau_max=6)))
    return cases


def analyse(name, row, tag, fn):
    t0 = time.time()
    expr = fn()
    atoms, reps = decompose(expr)
    r, polys, nsdim = minimal_operator(atoms, reps)
    rK, _ = min_order(atoms, reps, split_consts=False)
    sym, num = verify(expr, polys)
    fuchs, ind_roots, ind = indicial_at_zero(polys)
    finite = finite_singular(polys)
    inf_roots, chi = rates_at_infinity(polys)
    B = local_log_part(atoms, reps)
    Blead = None
    if B != 0:
        ser = sp.series(B.xreplace({v: k for k, v in reps.items()}), R, 0, 16).removeO()
        Blead = str(sp.Poly(ser, R).terms()[-1]) if ser != 0 else "0 to O(R^16)"
    r0_singular = sp.Poly(polys[-1], R).eval(0) == 0
    rec = {
        "class": name, "row": row, "params": tag,
        "atoms": [f"{c}|{kind}|{k}" for (c, kind, k) in atoms],
        "constants": {str(v): str(k) for k, v in reps.items()},
        "order_QR": r, "nullspace_dim_at_min_degree": nsdim, "order_KR_consts_in_field": rK,
        "coeff_degrees": [sp.Poly(p, R).degree() if p != 0 else None for p in polys],
        "operator": [str(sp.factor(p)) for p in polys],
        "verified_symbolic": bool(sym), "verified_numeric_rel": num,
        "R0_singular": bool(r0_singular), "R0_fuchsian": bool(fuchs), "lnR_coeff_leading_(power,coeff)": Blead, "R0_indicial": str(ind), "R0_exponents": ind_roots,
        "f_has_lnR_at_0": bool(B != 0), "lnR_coefficient": str(B),
        "finite_singular": finite,
        "inf_char_poly": str(chi), "inf_rates": inf_roots,
        "atom_rates": [str(x) for x in atom_rates(atoms)],
        "seconds": round(time.time() - t0, 2),
    }
    return rec, expr


def main() -> None:
    out = {"cases": [], "cross_checks": {}}
    for name, row, tag, fn in build_panel():
        rec, expr = analyse(name, row, tag, fn)
        out["cases"].append(rec)
        print(f"[{row}] {name:28s} {tag:16s} order {rec['order_QR']} (K:{rec['order_KR_consts_in_field']})"
              f" deg {max(d for d in rec['coeff_degrees'] if d is not None)}"
              f" sym={rec['verified_symbolic']} num={rec['verified_numeric_rel']:.1e}"
              f" R0 {rec['R0_exponents']} fuchs={rec['R0_fuchsian']} ln={rec['f_has_lnR_at_0']}"
              f" inf {rec['inf_rates']} fin {[(f['factor'], f.get('exponents', f.get('exponents_numeric'))) for f in rec['finite_singular']]}"
              f" ({rec['seconds']}s)", flush=True)
        if name.startswith("(AA|BB)") and tag == "zA=1,zB=1":
            ser = sp.series(expr, R, 0, 11).removeO()
            out["cross_checks"]["silverstone_odd"] = {
                "odd_coeffs_R1_R3_R5": [str(ser.coeff(R, k)) for k in (1, 3, 5)],
                "R7": str(ser.coeff(R, 7)), "R9": str(ser.coeff(R, 9)),
                "match": [ser.coeff(R, k) for k in (1, 3, 5)] == [0, 0, 0]
                and ser.coeff(R, 7) == sp.Rational(2, 315)
                and ser.coeff(R, 9) == sp.Rational(16, 14175)}
            print("  Silverstone odd family:", out["cross_checks"]["silverstone_odd"])
            # the same coefficients from the OPERATOR alone: the analytic solution at R=0
            # (exponents -1, 0) normalised by J(0) = 5/8
            P = [sp.sympify(p, locals={"R": R}) for p in rec["operator"]]
            N = 14
            a = sp.symbols(f"a0:{N}")
            y = sum(a[i] * R ** i for i in range(N))
            res = sp.expand(sum(p * sp.diff(y, R, j) for j, p in enumerate(P)))
            sol = list(sp.linsolve([res.coeff(R, k) for k in range(N - 4)]
                                   + [a[0] - sp.Rational(5, 8)], a))[0]
            out["cross_checks"]["silverstone_from_operator"] = {
                "taylor_R0_to_R10": [str(x) for x in sol[:11]],
                "match": sol[7] == sp.Rational(2, 315) and sol[9] == sp.Rational(16, 14175)
                and sol[1] == sol[3] == sol[5] == 0}
            print("  from operator:", out["cross_checks"]["silverstone_from_operator"])
    # Pachucki / Paper 59 Bessel limit: r f'' + f' - p^2 r f = 0
    p = sp.Symbol("p", positive=True)
    bess = [-p ** 2 * R, sp.Integer(1), R]
    fu, ir, ind = indicial_at_zero([sp.Poly(x, R).as_expr() for x in bess]) if False else (None, None, None)
    th = sp.Symbol("theta")
    out["cross_checks"]["pachucki_bessel"] = {
        "operator": "R D^2 + D - p^2 R", "R0_indicial": str(sp.factor(th * (th - 1) + th)),
        "R0_exponents": "{0: 2} (log solution K_0)", "inf_rates": "+p, -p",
        "check_K0": bool(sp.simplify(R * sp.diff(sp.besselk(0, p * R), R, 2)
                                     + sp.diff(sp.besselk(0, p * R), R)
                                     - p ** 2 * R * sp.besselk(0, p * R)) == 0)}
    path = ROOT / "debug" / "data" / "two_centre_r_operators.json"
    path.write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
    print("wrote", path)


if __name__ == "__main__":
    main()
