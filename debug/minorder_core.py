"""Minimal order in R of two-centre integrals (Paper 58 sec:census): exact lower bounds with TRUE constants.

Shared machinery for debug/minorder_run.py.  Nothing here edits a paper, test or doc.

Given a closed form f(R) = sum over atoms  coeff(R) * e^{cR} * {1, ln R, E1(kR), Ei(kR)}  with coeff in
Q(constants)(R) (constants: Euler gamma, logs of rationals, normalisation surds), this module

  * reads the atom decomposition with the guarded routine of tests/test_paper58_r_operators.py (loaded, not
    copied; that file's `decompose`, `to_field`, `d_atoms` are cross-checked against sympy's derivatives there);
  * closes the atom set under d/dR (structurally) and reports its size n: order <= n;
  * builds the TRUE coordinate matrix of f, f', ..., f^(r-1) on the atoms at a rational R0, entries polynomials
    in the constant symbols over Q, picks an r x r minor, computes it exactly as a polynomial in the constants,
    and evaluates that polynomial at the true constants in interval arithmetic (mpmath.iv).  An interval that
    excludes 0 proves rank r over C(R), with no independence assumption on gamma, ln 2, ... .

The independence of the atoms over C(R) is the hand lemma GI of debug/minorder_two_centre_memo.md.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import mpmath as mp
import sympy as sp
from sympy.polys.matrices import DomainMatrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

_spec = importlib.util.spec_from_file_location("_r_ops_minorder", ROOT / "tests" / "test_paper58_r_operators.py")
RO = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(RO)

R = RO.R
FR = RO.FR


def atom_label(a) -> Tuple[sp.Rational, str, sp.Rational]:
    """(c, kind, kappa): the atom is e^{cR} * y,  y = 1 | ln R (kappa = 0) | Ei(kappa R) (E1(kR) = -Ei(-kR))."""
    c, kind, k = a
    c, k = sp.Rational(c), sp.Rational(k)
    if kind == "1":
        return (c, "1", sp.Integer(0))
    if kind == "log":
        return (c, "y", sp.Integer(0))
    assert k > 0, a
    return (c, "y", -k if kind == "E1" else k)


def closure(atom_keys) -> set:
    """Smallest set of atoms containing atom_keys and closed under d/dR:
    (e^{cR} ln R)' brings e^{cR};  (e^{cR}E1(kR))' brings e^{(c-k)R};  (e^{cR}Ei(kR))' brings e^{(c+k)R}."""
    out = set(atom_keys)
    for (c, kind, k) in list(atom_keys):
        if kind == "log":
            out.add((c, "1", sp.Integer(0)))
        elif kind == "E1":
            out.add((sp.nsimplify(c - k), "1", sp.Integer(0)))
        elif kind == "Ei":
            out.add((sp.nsimplify(c + k), "1", sp.Integer(0)))
    return out


def true_constant_intervals(reps: Dict[sp.Expr, sp.Symbol], dps: int = 60):
    """{symbol: mpmath interval containing the true constant}; width 1e-(dps-10) relative."""
    out, mids = {}, {}
    for expr, sym in reps.items():
        v = sp.N(expr, dps + 10)
        assert v.is_real, (expr, v)
        x = mp.mpf(str(v))
        w = abs(x) * mp.mpf(10) ** (-(dps - 10))
        out[sym] = mp.iv.mpf([x - w, x + w])
        mids[sym] = x
    return out, mids


def iv_eval(poly: sp.Poly, ivals) -> "mp.iv.mpf":
    tot = mp.iv.mpf(0)
    gens = poly.gens
    for mono, cf in poly.terms():
        cf = sp.Rational(cf)
        term = mp.iv.mpf(int(cf.p)) / mp.iv.mpf(int(cf.q))
        for g, e in zip(gens, mono):
            if e:
                term = term * ivals[g] ** int(e)
        tot = tot + term
    return tot


def analyse(f: sp.Expr, r: int, R0=sp.Rational(37, 59), dps: int = 60, upper: bool = True) -> dict:
    """Exact data for 'f satisfies no ODE of order < r over C(R)' (given lemma GI)."""
    atoms, reps = RO.decompose(f)
    csyms = sorted(set(reps.values()), key=str)
    ds = [RO.to_field(atoms, reps, "split")]
    for _ in range(max(r, 1)):
        ds.append(RO.d_atoms(ds[-1]))
    f_atoms = sorted({a for (a, _m) in ds[0]}, key=str)
    clos = sorted(closure(f_atoms), key=str)
    seen = sorted({a for D in ds for (a, _m) in D}, key=str)
    assert set(seen) <= set(clos), (seen, clos)
    labels = [atom_label(a) for a in clos]
    assert len(set(labels)) == len(labels), "atom labels not distinct"
    rates = sorted({lab[0] for lab in labels}, key=lambda x: float(x))
    n_y = sum(1 for lab in labels if lab[1] == "y")

    # true coordinate matrix at R0: rows = atoms of the closure, columns = derivatives 0..r-1
    mp.iv.dps = dps
    with mp.workdps(dps):
        ivals, mids = true_constant_intervals(reps, dps)

        def entry(a, j):
            tot = sp.Integer(0)
            for (aa, m), w in ds[j].items():
                if aa != a:
                    continue
                val = sp.Rational(FR.to_sympy(w).subs(R, R0))
                mono = sp.Integer(1)
                for s, e in zip(csyms, m):
                    mono *= s ** e
                tot += val * mono
            return sp.expand(tot)

        T = [[entry(a, j) for j in range(r)] for a in clos]

        def num(e):
            return mp.mpf(str(sp.N(e.subs({s: sp.Float(mp.nstr(mids[s], dps), dps) for s in csyms}), dps))) if csyms \
                else mp.mpf(sp.Rational(e).p) / mp.mpf(sp.Rational(e).q)

        Tn = [[num(e) for e in row] for row in T]
        # greedy row choice by numerical elimination with pivoting
        rows_left = list(range(len(clos)))
        work = {i: Tn[i][:] for i in rows_left}
        chosen = []
        for col in range(r):
            if not rows_left:
                break
            piv = max(rows_left, key=lambda i: abs(work[i][col]))
            if abs(work[piv][col]) < mp.mpf(10) ** (-(dps // 2)):
                break
            chosen.append(piv)
            rows_left.remove(piv)
            for i in rows_left:
                fct = work[i][col] / work[piv][col]
                work[i] = [x - fct * y for x, y in zip(work[i], work[piv])]
        result = {
            "f_atoms": [str(a) for a in f_atoms],
            "closure_atoms": [str(a) for a in clos],
            "n_closure": len(clos),
            "n_rates_closure": len(rates),
            "rates": [str(x) for x in rates],
            "n_E_type_atoms": n_y,
            "constants": {str(v): str(k) for k, v in reps.items()},
            "claimed_order": r,
            "R0": str(R0),
        }
        if len(chosen) < r:
            result.update(rank_true=len(chosen), minor_interval=None, proved_lower=False)
            return result
        if csyms:
            dom = sp.QQ[tuple(csyms)]
            M = DomainMatrix([[dom.from_sympy(T[i][j]) for j in range(r)] for i in chosen], (r, r), dom)
            det = sp.Poly(dom.to_sympy(M.det()), *csyms)
            ivd = iv_eval(det, ivals)
            lo, hi = ivd.a, ivd.b
            nterms = len(det.terms())
            result["minor_monomials_in_constants"] = [
                " * ".join(f"{g}^{e}" for g, e in zip(det.gens, mono) if e) or "1" for mono, _ in det.terms()][:6]
            result["minor_free_of_gamma_and_logs"] = all(
                all(e == 0 for g, e in zip(det.gens, mono)
                    if str(g) == "gE" or str(dict((v, k) for k, v in reps.items())[g]).startswith("log"))
                for mono, _ in det.terms())
        else:
            M = DomainMatrix([[sp.QQ.from_sympy(T[i][j]) for j in range(r)] for i in chosen], (r, r), sp.QQ)
            d = sp.Rational(sp.QQ.to_sympy(M.det()))
            lo = hi = mp.mpf(d.p) / mp.mpf(d.q)
            nterms = 1
            det = d
        excl = (lo > 0) or (hi < 0)
        result.update(
            rank_true=r if excl else None,
            minor_rows=[str(clos[i]) for i in chosen],
            minor_interval=[mp.nstr(lo, 12), mp.nstr(hi, 12)],
            minor_is_exact_rational=not csyms,
            minor_terms_in_constants=nterms,
            proved_lower=bool(excl),
            upper_by_closure=len(clos),
            order_exact_by_closure=bool(excl) and len(clos) == r,
        )
        if len(clos) > r and upper:
            # order <= r: generic rank over Q(R) of the split (atom, constant-monomial) x (r+1) matrix equals r.
            # A Q(R)-relation among the split vectors of f..f^(r) is a relation among the true functions.
            keys = sorted(set().union(*[set(D) for D in ds[:r + 1]]), key=str)
            Mfull = DomainMatrix([[D.get(kk, FR.zero) for D in ds[:r + 1]] for kk in keys],
                                 (len(keys), r + 1), FR)
            result["generic_rank_r_plus_1_columns"] = int(Mfull.rank())
            result["upper_proved"] = result["generic_rank_r_plus_1_columns"] <= r
        else:
            result["upper_proved"] = len(clos) <= r
        return result
