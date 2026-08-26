"""AHA Track-3, OBJECT 2: resurgent data of the HYBRID two-centre ERI class
{E_1, ln} (Paper 58 / docs/neumann_general_m_build_plan.md sections 8.4-8.4.2).

The test: does GeoVac's Layer-2 chemistry transcendental carry ALGEBRAIC Stokes /
monodromy data over the parameter field Q(Z_A, Z_B) (up to a 2*pi*i normalisation),
with the genuinely transcendental content (the logs) confined to a boundary term?

Pipeline (exact-first, PSLQ only on the residual numeric constants):
  [A] exact term census of hybrid_closed_form for a sweep of quartets
  [B] trans-series regrouping F(R) = sum_A e^{-A R} phi_A(1/R); action lattice
  [C] Borel singularity table: every sector's local Borel poles sit at
      xi = -(A -+ a_d), i.e. global zeta = -+ a_d  (the far-centre orbital exponent)
  [D] Stokes constants, TWO independent numeric routes to >= 25 digits:
        (D1) blind large-order extrapolation of phi_A's coefficients
        (D2) exact lateral discontinuity of the assembled sector function
      then PSLQ against the algebraic candidates over Q(Z_A, Z_B)
  [E] the discriminating audit: do the LOGS ever enter the Stokes data?
  [F] closed-form validation against the independent quadrature route.

Run:  python debug/aha_t3_object2_hybrid.py
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction

import mpmath as mp
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from aha_t3_hybrid_lib import (  # noqa: E402
    group_elementary, group_logs, group_prefactors, parse_closed_form,
    phi_coefficients, sector_actions,
)
from geovac.two_center_eri import (  # noqa: E402
    R_s, hybrid_closed_form, hybrid_quadrature, radial_poly,
)

DPS = 80
mp.mp.dps = DPS
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


# quartets: (label, ZA, oa, ob, oc, ZB, od).  oa/ob = the one-centre multipole pair
# (needs l > 0 for the {E_1, ln} route), oc = third A orbital, od = the lone B orbital.
QUARTETS = [
    ("(2p0 2p0|1s   1s_B)  ZA=3 ZB=1", Fraction(3), (2, 1, 0), (2, 1, 0), (1, 0, 0), Fraction(1), (1, 0, 0)),
    ("(2p0 2p0|2p0  1s_B)  ZA=3 ZB=1", Fraction(3), (2, 1, 0), (2, 1, 0), (2, 1, 0), Fraction(1), (1, 0, 0)),
    ("(2p1 2p1|1s   1s_B)  ZA=3 ZB=1", Fraction(3), (2, 1, 1), (2, 1, 1), (1, 0, 0), Fraction(1), (1, 0, 0)),
    ("(2p1 2p1|2p1 2p1_B)  ZA=3 ZB=2", Fraction(3), (2, 1, 1), (2, 1, 1), (2, 1, 1), Fraction(2), (2, 1, 1)),
    ("(3d0 3d0|1s   1s_B)  ZA=3 ZB=1", Fraction(3), (3, 2, 0), (3, 2, 0), (1, 0, 0), Fraction(1), (1, 0, 0)),
    ("(2p0 2p0|1s   2s_B)  ZA=4 ZB=3", Fraction(4), (2, 1, 0), (2, 1, 0), (1, 0, 0), Fraction(3), (2, 0, 0)),
    ("(2p0 2p0|1s   1s_B)  ZA=2 ZB=1", Fraction(2), (2, 1, 0), (2, 1, 0), (1, 0, 0), Fraction(1), (1, 0, 0)),
    ("(2p0 2p0|2s   1s_B)  ZA=5 ZB=2", Fraction(5), (2, 1, 0), (2, 1, 0), (2, 0, 0), Fraction(2), (1, 0, 0)),
    ("(3d1 3d1|2p0  1s_B)  ZA=3 ZB=1", Fraction(3), (3, 2, 1), (3, 2, 1), (2, 1, 0), Fraction(1), (1, 0, 0)),
]


def rate(Z, o):
    return sp.Rational(radial_poly(Z, o[0], o[1])[1])


# --------------------------------------------------------------------- helpers
def poly_dict_str(d, g):
    return {int(k): sp.radsimp(sp.cancel(v / g)) for k, v in sorted(d.items(), key=lambda x: -int(x[0]))}


def common_algebraic_factor(terms):
    """The single algebraic generator g = sqrt(d) (d squarefree) with every
    coefficient in Q + Q*g.  Returns None if the coefficients need more."""
    g = sp.Integer(1)
    for t in terms:
        c = sp.sympify(t.coeff)
        if c.is_Rational:
            continue
        sq = sp.cancel(sp.expand(c ** 2))
        if not sq.is_Rational:
            return None
        num, den = sp.Rational(sq).p, sp.Rational(sq).q
        d = sp.Integer(1)
        for pr, ex in sp.factorint(num * den).items():
            if ex % 2:
                d *= pr
        cand = sp.sqrt(d)
        if g == 1:
            g = cand
        elif sp.simplify(g - cand) != 0:
            return None
    for t in terms:
        if not sp.radsimp(sp.cancel(t.coeff / g)).is_Rational:
            return None
    return g


def minpoly_degrees(terms):
    degs = set()
    for t in terms:
        c = sp.sympify(t.coeff)
        try:
            degs.add(sp.degree(sp.minimal_polynomial(c, sp.Symbol("x")), sp.Symbol("x")))
        except Exception:
            degs.add(-1)
    return sorted(degs)


# ------------------------------------------------------------- [D1] large order
def large_order_constant(terms, A, lam_min, N0=300, M=6):
    """Blind extraction of  p0 = lim_N u_N ,  u_N := a_N lam_min^N /((-1)^{N-1}(N-1)!).

    Uses ONLY the coefficient sequence a_N of phi_A (what an experimentalist has),
    with the standard large-order ansatz
        u_N = p0 + sum_{m>=1} pi_m / ((N-1)(N-2)...(N-m))  + O((lam_min/lam_max)^N)
    solved EXACTLY over Q on M+1 nodes (exact rational linear algebra -- no
    conditioning loss).  Also returns the raw un-extrapolated u_Nmax and a
    Neville-extrapolated cross-check.
    """
    Nmax = N0 + M
    a = phi_coefficients(terms, A, Nmax)
    us = {}
    for N in range(N0, Nmax + 1):
        us[N] = a[N] * lam_min ** N / ((-1) ** (N - 1) * sp.factorial(N - 1))
    rows, rhs = [], []
    for j in range(M + 1):
        N = N0 + j
        basis = [sp.Integer(1)]
        prod = sp.Integer(1)
        for m in range(1, M + 1):
            prod *= sp.Integer(N - m)
            basis.append(sp.Rational(1, 1) / prod)
        rows.append(basis)
        rhs.append(us[N])
    sol = sp.Matrix(rows).solve(sp.Matrix(rhs))
    p0_fit = sp.simplify(sol[0])
    return (mp.mpf(sp.N(p0_fit, DPS)),
            mp.mpf(sp.N(us[Nmax], DPS)))


# ------------------------------------------------------- [D2] exact discontinuity
def sector_function_numeric(terms, A, R):
    """f_A(R) = sum over terms of action A, evaluated in mpmath at complex R."""
    tot = mp.mpc(0)
    for t in terms:
        if t.action != A:
            continue
        co = mp.mpf(sp.N(t.coeff, DPS))
        val = co * R ** int(t.k) * mp.e ** (-mp.mpf(sp.N(t.c, DPS)) * R)
        if t.lam is not None:
            val *= mp.e1(mp.mpf(sp.N(t.lam, DPS)) * R)
        elif t.logarg is not None:
            val *= mp.log(mp.mpf(sp.N(t.logarg, DPS)))
        tot += val
    return tot


def full_function_numeric(terms, R):
    tot = mp.mpc(0)
    for t in terms:
        co = mp.mpf(sp.N(t.coeff, DPS))
        val = co * R ** int(t.k) * mp.e ** (-mp.mpf(sp.N(t.c, DPS)) * R)
        if t.lam is not None:
            val *= mp.e1(mp.mpf(sp.N(t.lam, DPS)) * R)
        elif t.logarg is not None:
            val *= mp.log(mp.mpf(sp.N(t.logarg, DPS)))
        tot += val
    return tot


def predicted_disc_numeric(terms, A, R):
    """-2 pi i * sum_{(c,lam): c+lam=A} P_{c,lam}(R) e^{-c R}   (E_1 cut = -2 pi i)."""
    tot = mp.mpc(0)
    for t in terms:
        if t.action != A or t.lam is None:
            continue
        co = mp.mpf(sp.N(t.coeff, DPS))
        tot += co * R ** int(t.k) * mp.e ** (-mp.mpf(sp.N(t.c, DPS)) * R)
    return -2j * mp.pi * tot



# --------------------------------------------- [E2] the log cross-ratio law
def log_crossratio_law(res):
    """Is the {ln} half of the seed set determined by the Borel actions?

    The four E_1 rates are lam in {a_c -+ a_d, mu -+ a_d} (the Borel singularity
    positions measured from their sectors).  Test whether the WHOLE logarithmic
    content of the closed form is

        L(R) = s * P_(a_d, mu-a_d)(R) * log(Lambda) ,
        Lambda = (a_c - a_d)(mu + a_d) / [(a_c + a_d)(mu - a_d)]

    i.e. the log of the multiplicative cross-ratio of the four Borel actions,
    carried by an E_1 prefactor polynomial -- no independent transcendental.
    """
    a_c, a_d, mu = res["a_c"], res["a_d"], res["mu"]
    Lam = sp.Rational(1) * ((a_c - a_d) * (mu + a_d)) / ((a_c + a_d) * (mu - a_d))
    all_at_ad = all(c == a_d for (c, _arg) in res["lgs"])
    Ld = {}
    for (c, arg), d in res["lgs"].items():
        for k, v in d.items():
            Ld[k] = Ld.get(k, sp.Integer(0)) + v * sp.log(arg)
    P = res["pre"].get((a_d, mu - a_d))
    if P is None:
        return dict(ok=False, why="no (a_d, mu-a_d) E_1 group", Lam=Lam)
    ratios = []
    for k, v in Ld.items():
        if k in P and sp.simplify(P[k]) != 0:
            ratios.append(sp.simplify(sp.expand_log(v / P[k], force=True)))
    if not ratios:
        return dict(ok=False, why="no overlapping R-powers", Lam=Lam)
    uniform = all(sp.simplify(r - ratios[0]) == 0 for r in ratios)
    s_sign = None
    for sgn in (1, -1):
        if sp.simplify(ratios[0] - sgn * sp.log(Lam)) == 0:
            s_sign = sgn
            break
    return dict(ok=(uniform and s_sign is not None), uniform=uniform, sign=s_sign,
                Lam=Lam, ratio=ratios[0], logs_all_in_leading_sector=all_at_ad)


# ------------------------------------------------------------------- the analysis
def analyse(label, ZA, oa, ob, oc, ZB, od, verbose=True):
    t0 = time.time()
    expr = hybrid_closed_form(ZA, oa, ob, oc, ZB, od, R_s)
    terms = parse_closed_form(expr)

    b_A = rate(ZA, oa) + rate(ZA, ob)
    a_c = rate(ZA, oc)
    a_d = rate(ZB, od)
    mu = b_A + a_c

    acts = sector_actions(terms)
    pre = group_prefactors(terms)
    lgs = group_logs(terms)
    elem = group_elementary(terms)
    g = common_algebraic_factor(terms) or sp.Integer(1)

    res = dict(label=label, terms=terms, acts=acts, pre=pre, lgs=lgs, elem=elem,
               b_A=b_A, a_c=a_c, a_d=a_d, mu=mu, g=g, secs=time.time() - t0)

    if not verbose:
        return res

    say("")
    say("-" * 78)
    say(f"QUARTET {label}")
    say("-" * 78)
    say(f"  orbital rates:  b_A(pair) = {b_A}   a_c = {a_c}   mu = b_A + a_c = {mu}   a_d(far) = {a_d}")
    say(f"  {len(terms)} closed-form terms; common algebraic factor g = {g}; "
        f"minpoly degrees over Q = {minpoly_degrees(terms)}")
    say(f"  exponential actions present: {acts}")
    say("")
    say("  [A] E_1 groups   (c, lam) -> prefactor P_(c,lam)(R) coefficients in units of g")
    for key in sorted(pre, key=lambda kv: (float(kv[0]), float(kv[1]))):
        c, lam = key
        say(f"      c={str(c):>6}  lam={str(lam):>6}  action A=c+lam={str(c + lam):>6}   P = {poly_dict_str(pre[key], g)}")
    say("  [A] log groups   (c, arg) -> coefficients in units of g")
    for key in sorted(lgs, key=lambda kv: (float(kv[0]), float(kv[1]))):
        say(f"      c={str(key[0]):>6}  arg={str(key[1])}   {poly_dict_str(lgs[key], g)}")
    say("  [A] elementary groups  c -> coefficients in units of g")
    for c in sorted(elem, key=float):
        say(f"      c={str(c):>6}   {poly_dict_str(elem[c], g)}")
    return res


def structural_checks(res):
    """The exact structural statements, checked symbolically."""
    terms, pre, lgs = res["terms"], res["pre"], res["lgs"]
    a_d, b_A, mu, g = res["a_d"], res["b_A"], res["mu"], res["g"]
    a_c = res["a_c"]
    ok = {}

    # (S1) no term carries BOTH a log and an E_1  -> logs cannot enter Stokes data
    ok["S1 no term has log*E1"] = all(not (t.lam is not None and t.logarg is not None) for t in terms)
    # (S2) every log argument is R-INDEPENDENT and rational (a ratio of rates)
    ok["S2 log args rational, R-free"] = all(sp.sympify(t.logarg).is_Rational
                                             for t in terms if t.logarg is not None)
    # (S3) every E_1 rate lam and every exp rate c is rational in the orbital exponents
    ok["S3 rates rational"] = all(sp.sympify(t.c).is_Rational and
                                  (t.lam is None or sp.sympify(t.lam).is_Rational) for t in terms)
    # (S4) the exp rate accompanying every E_1 is exactly -+ a_d
    ok["S4 exp rate of E1 terms = +- a_d"] = all(abs(sp.sympify(t.c)) == a_d
                                                 for t in terms if t.lam is not None)
    # (S5) actions of the E_1 sectors are exactly the A-side aggregate rates {b_A, mu}
    e1_actions = {sp.sympify(t.action) for t in terms if t.lam is not None}
    ok["S5 E1 sector actions = {a_c, mu}"] = e1_actions <= {a_c, mu}
    # (S6) cut cancellation: for each c, sum over lam of P_(c,lam) = 0 identically in R
    cancel = True
    per_c = {}
    for (c, lam), d in pre.items():
        per_c.setdefault(c, {})
        for k, v in d.items():
            per_c[c][k] = per_c[c].get(k, sp.Integer(0)) + v
    for c, d in per_c.items():
        if any(sp.simplify(v) != 0 for v in d.values()):
            cancel = False
    ok["S6 sum_lam P_(c,lam) = 0 (cut cancellation)"] = cancel
    # (S7) all coefficients lie in Q(g), g algebraic of degree <= 2
    ok["S7 coeffs in Q(g)"] = all(sp.radsimp(sp.cancel(t.coeff / g)).is_Rational for t in terms)
    # (S8) the four E_1 rates are exactly {a_c -+ a_d, mu -+ a_d}
    lam_set = {sp.sympify(t.lam) for t in terms if t.lam is not None}
    ok["S8 E1 rates = {a_c-+a_d, mu-+a_d}"] = lam_set == {a_c - a_d, a_c + a_d, mu - a_d, mu + a_d}
    # (S9) leading sector carries NO E_1 -> its series TERMINATES (exactly summable)
    A_lead = min((sp.sympify(t.action) for t in terms), key=float)
    ok["S9 leading sector E_1-free (terminates)"] = all(
        t.lam is None for t in terms if sp.sympify(t.action) == A_lead)
    # (S10) every log sits in the leading sector
    ok["S10 all logs in leading sector"] = all(
        sp.sympify(t.c) == A_lead for t in terms if t.logarg is not None)
    return ok, per_c


def main():
    os.makedirs("debug/data", exist_ok=True)
    say("=" * 78)
    say("OBJECT 2 -- the HYBRID two-centre ERI class {E_1, ln} (Paper 58)")
    say(f"mpmath dps = {DPS}")
    say("=" * 78)

    results = []
    for q in QUARTETS:
        results.append(analyse(*q))

    # ---------------------------------------------------------------- structural
    say("")
    say("=" * 78)
    say("[B/E] STRUCTURAL CHECKS (exact, symbolic) across the quartet sweep")
    say("=" * 78)
    keys = None
    for r in results:
        ok, per_c = structural_checks(r)
        r["per_c"] = per_c
        if keys is None:
            keys = list(ok)
            say(f"{'quartet':<34} " + "  ".join(f"S{i + 1}" for i in range(len(keys))))
        say(f"{r['label']:<34} " + "  ".join((" OK" if ok[k] else "FAIL") for k in keys))
        r["ok"] = ok
    say("")
    for i, k in enumerate(keys):
        say(f"   S{i + 1}: {k}")

    # ------------------------------------------------------------------- [C] Borel
    say("")
    say("=" * 78)
    say("[C] BOREL SINGULARITY TABLE")
    say("=" * 78)
    say("  For a sector at action A, the divergent piece is P(R) e^{-cR} E_1(lam R) with")
    say("  c + lam = A, so phi_A ~ P(R) sum_n (-1)^n n!/(lam R)^{n+1}: local Borel transform")
    say("  singular at xi = -lam, i.e. GLOBAL position zeta = A - lam = c.")
    say("")
    say(f"  {'quartet':<34} {'A':>8} {'lam':>8} {'xi_sing':>9} {'zeta_glob':>10} {'= a_d?':>8} {'type':>18}")
    for r in results:
        for (c, lam), d in sorted(r["pre"].items(), key=lambda kv: (float(kv[0][0] + kv[0][1]), float(kv[0][1]))):
            A = c + lam
            kmin = min(int(k) for k in d)          # most negative power of R present
            typ = "simple pole" if max(int(k) for k in d) == 0 and len(d) == 1 else f"pole/log, R-powers {kmin}..{max(int(k) for k in d)}"
            say(f"  {r['label']:<34} {str(A):>8} {str(lam):>8} {str(-lam):>9} {str(c):>10} "
                f"{('yes' if abs(c) == r['a_d'] else 'NO'):>8} {typ:>18}")

    # ------------------------------------------------------- [D] Stokes constants
    say("")
    say("=" * 78)
    say("[D] STOKES CONSTANTS -- two independent numeric routes + PSLQ")
    say("=" * 78)
    say("  Convention: S_A := 2 pi i * p0 , p0 = lim_N a_N lam_min^N /((-1)^{N-1}(N-1)!)")
    say("  where a_N are the coefficients of phi_A and lam_min = A - a_d is the")
    say("  smallest Borel action in that sector.")
    say("")
    stokes_rows = []
    for r in results[:4]:                       # the four cheapest quartets
        terms, a_d = r["terms"], r["a_d"]
        for A in sorted({sp.sympify(t.action) for t in terms if t.lam is not None}, key=float):
            lams = sorted({sp.sympify(t.lam) for t in terms if t.lam is not None and t.action == A}, key=float)
            lam_min = lams[0]
            # exact p0 (route D0, from the closed form)
            p0_exact = sum(t.coeff for t in terms
                           if t.lam == lam_min and t.action == A and int(t.k) == 0)
            p0_exact = sp.radsimp(sp.expand(p0_exact))
            # blind large-order extraction (route D1)
            p0_num, p0_raw = large_order_constant(terms, A, lam_min)
            err = abs(p0_num - mp.mpf(sp.N(p0_exact, DPS)))
            rel = err / max(abs(p0_num), mp.mpf(1))
            stokes_rows.append((r["label"], A, lam_min, p0_exact, p0_num, rel, r["g"]))
            say(f"  {r['label']}   sector A={A}, lam_min={lam_min}")
            say(f"     D1 blind large-order p0 = {mp.nstr(p0_num, 32)}")
            say(f"     D0 exact from closed form = {p0_exact} = {mp.nstr(mp.mpf(sp.N(p0_exact, DPS)), 32)}")
            say(f"     agreement (relative)      = {mp.nstr(rel, 5)}   [un-extrapolated u_N = {mp.nstr(p0_raw, 12)}]")
            # PSLQ against the algebraic candidates over Q(Z_A, Z_B)
            g = mp.mpf(sp.N(r["g"], DPS))
            rel_pslq = mp.pslq([p0_num, mp.mpf(1), g], tol=mp.mpf(10) ** (-30),
                               maxcoeff=10 ** 10, maxsteps=10 ** 6)
            say(f"     PSLQ [p0, 1, g={r['g']}] -> {rel_pslq}")
            say(f"     => S_A = 2 pi i * ({sp.radsimp(sp.cancel(p0_exact / r['g']))}) * {r['g']}   ALGEBRAIC")
            say("")

    # ------------------------------- [D2] exact lateral discontinuity, numerically
    say("=" * 78)
    say("[D2] LATERAL DISCONTINUITY of each sector function across arg R = pi")
    say("=" * 78)
    say("  Disc f_A = -2 pi i * sum_{c+lam=A} P_(c,lam)(R) e^{-cR}   (E_1 cut = -2 pi i, exact)")
    say(f"  {'quartet':<34} {'A':>6} {'x':>5} {'rel. residual':>16}")
    worst = mp.mpf(0)
    eps = mp.mpf(10) ** (-50)
    for r in results[:4]:
        terms = r["terms"]
        for A in sorted({sp.sympify(t.action) for t in terms if t.lam is not None}, key=float):
            for x in (mp.mpf("1.3"), mp.mpf("2.7")):
                Rp, Rm = mp.mpc(-x, eps), mp.mpc(-x, -eps)
                disc = sector_function_numeric(terms, A, Rp) - sector_function_numeric(terms, A, Rm)
                pred = predicted_disc_numeric(terms, A, mp.mpc(-x, 0))
                rel = abs(disc - pred) / abs(pred)
                worst = max(worst, rel)
                say(f"  {r['label']:<34} {str(A):>6} {float(x):>5.1f} {mp.nstr(rel, 6):>16}")
    say(f"  worst relative residual: {mp.nstr(worst, 5)}")

    # --------------------------------------------- the physical function is cut-free
    say("")
    say("=" * 78)
    say("[E] IS THE PHYSICAL ERI ITSELF CUT-FREE?  (do the algebraic Stokes constants cancel)")
    say("=" * 78)
    say(f"  {'quartet':<34} {'x':>5} {'|Disc F| / |F|':>18}")
    for r in results[:4]:
        for x in (mp.mpf("1.3"), mp.mpf("2.7")):
            Rp, Rm = mp.mpc(-x, eps), mp.mpc(-x, -eps)
            F = full_function_numeric(r["terms"], mp.mpc(-x, eps))
            d = full_function_numeric(r["terms"], Rp) - full_function_numeric(r["terms"], Rm)
            say(f"  {r['label']:<34} {float(x):>5.1f} {mp.nstr(abs(d) / abs(F), 6):>18}")

    # --------------------------------------------- [E2] the log cross-ratio law
    say("")
    say("=" * 78)
    say("[E2] WHERE DOES THE ln HALF OF THE SEED SET LIVE?")
    say("=" * 78)
    say("  Claim under test: the whole logarithmic content is")
    say("     L(R) = s * P_(a_d, mu-a_d)(R) * log(Lambda),")
    say("     Lambda = (a_c-a_d)(mu+a_d) / [(a_c+a_d)(mu-a_d)]")
    say("  = log of the multiplicative cross-ratio of the FOUR Borel actions,")
    say("    carried by an E_1 prefactor polynomial. If true, the ln is not an")
    say("    independent transcendental: it is a function of the Stokes data.")
    say("")
    say(f"  {'quartet':<34} {'Lambda':>12} {'sign':>5} {'uniform':>8} {'logs in leading sector':>24} {'VERDICT':>9}")
    for r in results:
        law = log_crossratio_law(r)
        r["law"] = law
        say(f"  {r['label']:<34} {str(law['Lam']):>12} {str(law.get('sign')):>5} "
            f"{str(law.get('uniform')):>8} {str(law.get('logs_all_in_leading_sector')):>24} "
            f"{('PASS' if law['ok'] else 'FAIL'):>9}")
    say("")
    say("  DECOY CONTROL -- the same test with the two 'minus' rates swapped for the")
    say("  'plus' ones, Lambda_decoy = (a_c+a_d)(mu+a_d)/[(a_c-a_d)(mu-a_d)].  A test that")
    say("  passes for a wrong Lambda would be vacuous.")
    say(f"  {'quartet':<34} {'Lambda_decoy':>14} {'VERDICT':>9}")
    for r in results:
        a_c, a_d, mu = r["a_c"], r["a_d"], r["mu"]
        dec = sp.Rational(1) * ((a_c + a_d) * (mu + a_d)) / ((a_c - a_d) * (mu - a_d))
        rat = r["law"].get("ratio")
        hit = rat is not None and any(
            sp.simplify(rat - sgn * sp.log(dec)) == 0 for sgn in (1, -1))
        say(f"  {r['label']:<34} {str(dec):>14} {('PASS (BAD)' if hit else 'FAIL (good)'):>9}")
    say("")
    say("  (a_c = third A-orbital rate, mu = a_c + b_A, a_d = far-centre rate;")
    say("   the four E_1 rates are lam = a_c -+ a_d and mu -+ a_d.)")

    # -------------------------------------------------------------- [F] validation
    say("")
    say("=" * 78)
    say("[F] VALIDATION -- parsed closed form vs the independent quadrature route")
    say("=" * 78)
    say(f"  {'quartet':<34} {'R':>5} {'closed form':>22} {'quadrature':>22} {'|diff|':>12}")
    for (label, ZA, oa, ob, oc, ZB, od), r in zip(QUARTETS, results):
        Rv = 2.5
        cf = float(full_function_numeric(r["terms"], mp.mpf(Rv)).real)
        qd = hybrid_quadrature(ZA, oa, ob, oc, ZB, od, Rv)
        say(f"  {label:<34} {Rv:>5.1f} {cf:>22.14e} {qd:>22.14e} {abs(cf - qd):>12.2e}")

    with open("debug/data/aha_t3_object2_log.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
