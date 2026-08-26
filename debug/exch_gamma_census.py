"""Track B / exchange class -- STEP 1: exact term census of the EXCHANGE two-centre
ERI kernel's closed form, and the structural laws that set up the resurgence.

Object: geovac.two_center_eri.ordered_xi_closed(p1, p2)  (tau=0, sigma=0, j1=j2=0)
        with p_i = rate_i * R, rate_i rational -- Increment 3c, seeds {E_1, ln, gamma}.
Superset check: ordered_xi_general(tau, sigma, H1, H2, j1, j2, p1, p2) -- Increment 3e.

The T3 parser (aha_t3_hybrid_lib) RAISES on an R-dependent log; that raise is
exactly the boundary this track crosses, so the parser is re-implemented here with
the log split into  {R-free log}  and  {log R}  and with EulerGamma pulled out.

Run:  python debug/exch_gamma_census.py
"""
from __future__ import annotations

import os
from collections import defaultdict

import sympy as sp

from geovac.two_center_eri import R_s, ordered_xi_closed, ordered_xi_general

OUT = []
G = sp.Symbol("Gam")          # stand-in for EulerGamma (sympy will not .coeff() on it)


def say(s=""):
    print(s)
    OUT.append(s)


# rate pairs (alpha, beta) = (p1/R, p2/R); includes the equal-rate degenerate case
CASES = [
    (sp.Rational(3, 2), sp.Integer(1)),
    (sp.Integer(2), sp.Rational(5, 2)),
    (sp.Rational(1, 2), sp.Rational(7, 2)),
    (sp.Integer(1), sp.Integer(1)),
    (sp.Rational(5, 2), sp.Rational(3, 2)),
    (sp.Rational(7, 3), sp.Rational(4, 5)),
]


class T:
    """coeff * R^k * exp(-c R) * X,  X in {1, E1(lam R), lnR, ln(q)}; gam = gamma flag."""

    __slots__ = ("coeff", "k", "c", "lam", "logq", "lnR", "gam")

    def __init__(self, coeff, k, c, lam, logq, lnR, gam):
        self.coeff, self.k, self.c = sp.sympify(coeff), sp.Integer(k), sp.sympify(c)
        self.lam = None if lam is None else sp.sympify(lam)
        self.logq = None if logq is None else sp.sympify(logq)
        self.lnR, self.gam = bool(lnR), bool(gam)

    @property
    def kind(self):
        if self.lam is not None:
            return "E1"
        if self.lnR:
            return "lnR"
        if self.gam:
            return "gamma"
        if self.logq is not None:
            return "logq"
        return "elem"

    @property
    def action(self):
        return self.c + (self.lam if self.lam is not None else 0)

    def __repr__(self):
        x = ("E1(%sR)" % self.lam if self.lam is not None else
             "lnR" if self.lnR else "gamma" if self.gam else
             "log(%s)" % self.logq if self.logq is not None else "1")
        return "T(%s, R^%s, e^(-%sR), %s)" % (self.coeff, self.k, self.c, x)


def parse(expr, R=R_s):
    """Exact parse. Raises on ANY shape outside {exp, E_1, log, gamma} * R^k."""
    expr = sp.expand(sp.expand(expr).subs(sp.EulerGamma, G))
    out = []
    for t in sp.Add.make_args(expr):
        rest, lam, logq, lnR, gam, c = t, None, None, False, False, sp.Integer(0)

        e1s = [f for f in t.atoms(sp.Function) if isinstance(f, sp.expint)]
        if len(e1s) > 1:
            raise ValueError("two E_1 in one term: %s" % t)
        if e1s:
            f = e1s[0]
            assert f.args[0] == 1, "expint order != 1: %s" % f
            lam = sp.simplify(f.args[1] / R)
            if lam.has(R):
                raise ValueError("E_1 arg not linear in R: %s" % f)
            rest = rest / f

        lgs = list(t.atoms(sp.log))
        if len(lgs) > 1:
            raise ValueError("two logs in one term: %s" % t)
        if lgs:
            f = lgs[0]
            if f.args[0] == R:
                lnR = True
            elif f.args[0].has(R):
                raise ValueError("log arg neither R nor R-free: %s" % f)
            else:
                logq = f.args[0]
            rest = rest / f

        if rest.has(G):
            gam = True
            rest = rest / G
            assert not rest.has(G)

        exs = [f for f in t.atoms(sp.exp)]
        if len(exs) > 1:
            raise ValueError("two exps in one term: %s" % t)
        if exs:
            f = exs[0]
            c = sp.simplify(-f.args[0] / R)
            if c.has(R):
                raise ValueError("exp arg not linear in R: %s" % f)
            rest = rest / f

        rest = sp.powsimp(sp.cancel(sp.simplify(rest)))
        coeff, k = rest.as_coeff_exponent(R)
        if coeff.has(R):
            raise ValueError("unparsed R dependence: %s -> %s" % (t, rest))
        out.append(T(coeff, k, c, lam, logq, lnR, gam))
    return out


def summarize(terms):
    """(sector actions, {action: kinds}, gamma-coeffs, lnR-coeffs)."""
    acts = sorted({t.action for t in terms}, key=lambda z: float(z))
    kinds = defaultdict(set)
    for t in terms:
        kinds[t.action].add(t.kind)
    cg, cl = defaultdict(lambda: sp.Integer(0)), defaultdict(lambda: sp.Integer(0))
    for t in terms:
        if t.gam:
            cg[(t.c, t.k)] += t.coeff
        if t.lnR:
            cl[(t.c, t.k)] += t.coeff
    cg = {k: v for k, v in cg.items() if v != 0}
    cl = {k: v for k, v in cl.items() if v != 0}
    return acts, dict(kinds), cg, cl


def main():
    os.makedirs("debug/data", exist_ok=True)
    say("=" * 96)
    say("EXCHANGE CLASS -- STEP 1: exact term census (ordered_xi_closed, tau=0)")
    say("=" * 96)

    say("  %12s %7s %24s %30s %9s" % ("(a,b)", "#terms", "sectors", "kinds@action", "gam==lnR"))
    rows = []
    all_single, all_bundle = True, True
    for al, be in CASES:
        terms = parse(ordered_xi_closed(al * R_s, be * R_s))
        acts, kinds, cg, cl = summarize(terms)
        single = (len(acts) == 1)
        bundle = (cg == cl) and len(cg) > 0
        all_single &= single
        all_bundle &= bundle
        kd = ",".join(sorted(kinds[acts[0]])) if single else "MULTI"
        say("  %12s %7d %24s %30s %9s"
            % ("(%s,%s)" % (al, be), len(terms), [str(a) for a in acts], kd, bundle))
        rows.append((al, be, terms))

    n = len(CASES)
    say("")
    say("  L1  every term parses as coeff*R^k*e^(-cR)*X, X in {1,E1,lnR,log q,gamma}"
        "   : PASS %d/%d" % (n, n))
    say("  L2  the WHOLE object sits in ONE trans-series sector (single action A=a+b)"
        "  : %s %d/%d" % ("PASS" if all_single else "FAIL", n, n))
    say("  L3  coeff(gamma) == coeff(ln R) identically (the Ein bundle)"
        "                : %s %d/%d" % ("PASS" if all_bundle else "FAIL", n, n))

    say("")
    say("  Reduced normal form (checked as an EXACT symbolic identity, per case):")
    say("      F = 1/(2ab R^2) * [ e^{-(a+b)R}(gamma + ln(kappa R))")
    say("                          + e^{(a-b)R}E1(2aR) + e^{(b-a)R}E1(2bR)")
    say("                          - e^{(a+b)R}E1(2(a+b)R) ],   kappa = 2ab/(a+b)")
    say("")
    say("  %12s %10s %14s %14s %26s"
        % ("(a,b)", "kappa", "(2a)(2b)/(2A)", "NF residual", "E1 rates"))
    ok_nf, ok_kap = True, True
    R = R_s
    for al, be, terms in rows:
        A = al + be
        kap = 2 * al * be / A
        kap2 = (2 * al) * (2 * be) / (2 * A)
        nf = (1 / (2 * al * be * R ** 2)) * (
            sp.exp(-A * R) * (sp.EulerGamma + sp.log(kap * R))
            + sp.exp((al - be) * R) * sp.E1(2 * al * R)
            + sp.exp((be - al) * R) * sp.E1(2 * be * R)
            - sp.exp(A * R) * sp.E1(2 * A * R))
        res = sp.simplify(sp.expand(sp.expand(ordered_xi_closed(al * R, be * R))
                                    - sp.expand(nf)))
        rates = sorted({t.lam for t in terms if t.lam is not None}, key=lambda z: float(z))
        ok_nf &= (res == 0)
        ok_kap &= (sp.simplify(kap - kap2) == 0)
        say("  %12s %10s %14s %14s %26s"
            % ("(%s,%s)" % (al, be), kap, kap2, res, [str(r) for r in rates]))
    say("")
    say("  L4  reduced normal form EXACT (residual identically 0)"
        "                       : %s %d/%d" % ("PASS" if ok_nf else "FAIL", n, n))
    say("  L5  kappa = (2a)(2b)/(2A): charge-weighted product of the other three"
        "        : %s %d/%d" % ("PASS" if ok_kap else "FAIL", n, n))

    say("")
    say("  DECOY CONTROLS on L5 (a wrong charge assignment MUST fail):")
    decoys = [
        ("D1  (2a)(2b)(2A)          [charge +1 on 2A]",
         lambda a, b: (2 * a) * (2 * b) * (2 * (a + b))),
        ("D2  (2a)(2A)/(2b)         [charges permuted]",
         lambda a, b: (2 * a) * (2 * (a + b)) / (2 * b)),
        ("D3  (2a+2b)/2             [additive, not multiplicative]",
         lambda a, b: (2 * a + 2 * b) / 2),
        ("D4  (2a)(2b)/(2A)^2       [wrong charge weight]",
         lambda a, b: (2 * a) * (2 * b) / (2 * (a + b)) ** 2),
        ("D5  sqrt((2a)(2b))        [geometric mean]",
         lambda a, b: sp.sqrt((2 * a) * (2 * b))),
    ]
    for name, f in decoys:
        hits = sum(1 for al, be, _ in rows
                   if sp.simplify(2 * al * be / (al + be) - f(al, be)) == 0)
        verdict = "FAIL(as required)" if hits < len(rows) else "PASSED -- NOT A CONTROL"
        say("    %-46s matches %d/%d  -> %s" % (name, hits, len(rows), verdict))

    say("")
    say("=" * 96)
    say("STEP 1b -- ordered_xi_general: does the general kernel stay inside the census?")
    say("=" * 96)
    GEN = [(0, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0), (1, 1, 1, 1, 0, 0),
           (2, 0, 1, 0, 1, 0), (2, 2, 2, 2, 0, 1), (3, 1, 1, 2, 1, 1)]
    say("  %22s %7s %26s %6s %9s %8s"
        % ("(tau,sig,H1,H2,j1,j2)", "#terms", "sectors", "lnR?", "gam==lnR", "logq?"))
    gen_ok, gen_bundle = True, True
    al, be = sp.Rational(3, 2), sp.Integer(1)
    for prm in GEN:
        e = ordered_xi_general(*prm, al * R_s, be * R_s)
        try:
            terms = parse(e)
        except ValueError as ex:
            say("  %22s  PARSE RAISED: %s" % (str(prm), ex))
            gen_ok = False
            continue
        acts, kinds, cg, cl = summarize(terms)
        haslnR = any(t.lnR for t in terms)
        haslq = any(t.logq is not None for t in terms)
        bundle = (cg == cl)
        gen_bundle &= bundle
        say("  %22s %7d %26s %6s %9s %8s"
            % (str(prm), len(terms), [str(a) for a in acts], haslnR, bundle, haslq))
    say("")
    say("  L6  general kernel parses inside the same census"
        "                              : %s" % ("PASS" if gen_ok else "FAIL"))
    say("  L7  gamma/lnR bundle survives general (tau,sigma,H,j)"
        "                        : %s" % ("PASS" if gen_bundle else "FAIL"))

    with open("debug/data/exch_gamma_census.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
