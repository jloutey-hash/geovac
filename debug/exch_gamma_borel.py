"""Track B / exchange class -- STEPS 2-5: Borel transform, singularity table,
generalized Stokes data by THREE independent routes, classification, decoys.

Object (validated by tests/test_two_center_eri_aabb.py::test_ordered_xi_*):

    F(R) = ordered_xi_closed(a R, b R)
         = 1/(2ab R^2) [ e^{-(a+b)R}(gamma + ln(kappa R))
                         + e^{(a-b)R}E1(2aR) + e^{(b-a)R}E1(2bR)
                         - e^{(a+b)R}E1(2(a+b)R) ],   kappa = 2ab/(a+b)

(normal form established EXACTLY in debug/exch_gamma_census.py, L4 6/6).

The whole object is ONE trans-series sector, action A = a+b.  Write

    F(R) = e^{-A R} phi(R),    phi(R) = int_0^oo e^{-xi R} Bor(xi) d xi.

CLAIM (route 1, closed form) -- the LOCAL Borel transform is elementary:

    Bor(xi) = (1/(2ab)) [ -xi ln(xi/kappa)
                          + (xi+2a) ln((xi+2a)/(2a))
                          + (xi+2b) ln((xi+2b)/(2b))
                          - (xi+2A) ln((xi+2A)/(2A)) ]
            = (1/(2ab)) sum_j q_j [ (xi - xi_j) ln(xi - xi_j) - (-xi_j)ln(-xi_j) ]

  xi_j in {0, -2a, -2b, -2A},   q_j in {-1, +1, +1, -1},   sum_j q_j = 0.

So the R-DEPENDENT log is the q = -1 charge sitting at xi = 0 -- the DEGENERATE
member (lambda -> 0) of the same four-point family.  Its normalisation is forced:
kappa = prod_{j != 0} lambda_j^{q_j}.

Run:  python debug/exch_gamma_borel.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction

import mpmath as mp
import sympy as sp

from geovac.two_center_eri import R_s, ordered_xi_closed, ordered_xi_general

from exch_gamma_census import parse            # same directory

OUT = []
DPS = 90


def say(s=""):
    print(s)
    OUT.append(s)


CASES = [
    (sp.Rational(3, 2), sp.Integer(1)),
    (sp.Integer(2), sp.Rational(5, 2)),
    (sp.Rational(1, 2), sp.Rational(7, 2)),
    (sp.Integer(1), sp.Integer(1)),          # degenerate: 2a = 2b
    (sp.Rational(5, 2), sp.Rational(3, 2)),
    (sp.Rational(7, 3), sp.Rational(4, 5)),
]


# ----------------------------------------------------------------- generic tools

def sector_and_coeffs(expr, M):
    """From the PARSED closed form: (action A, {m: (rational, lnR-rational)}) for
    phi(R) = e^{A R} F(R) ~ sum_m (a_m + b_m ln R) R^{-m}, exact over Q.

    Deliberately generic: it reads the term census, not the hand normal form.
    gamma is kept symbolic (sp.EulerGamma).
    """
    terms = parse(expr)
    acts = {t.action for t in terms}
    assert len(acts) == 1, "not a single-sector object: %s" % acts
    A = acts.pop()
    a = {}
    b = {}

    def add(d, m, v):
        if 1 <= m <= M:
            d[m] = d.get(m, sp.Integer(0)) + v

    for t in terms:
        if t.kind == "E1":
            lam = t.lam
            for n in range(0, M + 3):
                m = n + 1 - int(t.k)
                if m > M:
                    break
                add(a, m, t.coeff * (-1) ** n * sp.factorial(n) / lam ** (n + 1))
        else:
            m = -int(t.k)
            sym = (sp.EulerGamma if t.gam else
                   sp.log(t.logq) if t.logq is not None else sp.Integer(1))
            if t.lnR:
                add(b, m, t.coeff)
            else:
                add(a, m, t.coeff * sym)
    return A, a, b


def bor_closed(al, be):
    """(list of (xi_j, q_j), Bor(xi) as an mpmath callable, prefactor 1/(2ab))."""
    A = al + be
    sings = [(sp.Integer(0), sp.Integer(-1)),
             (-2 * al, sp.Integer(1)),
             (-2 * be, sp.Integer(1)),
             (-2 * A, sp.Integer(-1))]
    kap = 2 * al * be / A
    pref = mp.mpf(1) / mp.mpf(sp.Rational(2 * al * be))

    def B(xi):
        xi = mp.mpf(xi)
        tot = -xi * mp.log(xi / mp.mpf(sp.Rational(kap)))
        for lam in (2 * al, 2 * be, 2 * A):
            s = 1 if lam != 2 * A else -1
            L = mp.mpf(sp.Rational(lam))
            tot += s * (xi + L) * mp.log((xi + L) / L)
        return pref * tot
    return sings, B


def phi_numeric(al, be, R):
    """phi(R) = e^{AR} F(R), evaluated from the closed form at mp.dps precision."""
    A = mp.mpf(sp.Rational(al + be))
    a, b = mp.mpf(sp.Rational(al)), mp.mpf(sp.Rational(be))
    R = mp.mpf(R)
    kap = 2 * a * b / A
    pre = 1 / (2 * a * b * R ** 2)
    return pre * (mp.euler + mp.log(kap * R)
                  + mp.e ** (2 * a * R) * mp.e1(2 * a * R)
                  + mp.e ** (2 * b * R) * mp.e1(2 * b * R)
                  - mp.e ** (2 * A * R) * mp.e1(2 * A * R))


# ------------------------------------------------------- blind large-order (Prony)

def _mpf(expr):
    """sympy exact -> mpf at the ambient precision (handles radicals)."""
    return mp.mpf(str(sp.N(expr, mp.mp.dps + 5)))



def blind_type_and_scale(a_exact, m0):
    """Blind (s, lambda_min) from three consecutive coefficients.

    Model  a_m ~ C (-1)^{m-1} (m-s)! lambda^{...}  =>  a_{m+1}/a_m = -(m+1-s)/lambda.
    Two consecutive ratios determine s and lambda with no prior assumption.
    """
    r1 = _mpf(a_exact[m0 + 1] / a_exact[m0])
    r2 = _mpf(a_exact[m0 + 2] / a_exact[m0 + 1])
    t = r1 / r2
    s = (m0 + 1) - t / (1 - t)
    lam = -(m0 + 1 - s) / r1
    return s, lam


def prony_exact(w, dmax=4):
    """Exact Prony over Q: given w_k = sum_j q_j x_j^k (k = 1..len), return the
    minimal d, the monic characteristic polynomial coefficients, and the residual
    check on the UNUSED tail. Returns (d, poly_coeffs_low_to_high, tail_ok)."""
    n = len(w)
    for d in range(1, dmax + 1):
        if n < 2 * d + 2:
            break
        Mt = sp.Matrix([[w[k + i] for i in range(d)] for k in range(d)])
        rhs = sp.Matrix([-w[k + d] for k in range(d)])
        if Mt.det() == 0:
            continue
        c = Mt.solve(rhs)                     # w_{k+d} + sum c_i w_{k+i} = 0
        ok = True
        for k in range(d, n - d):
            r = w[k + d] + sum(c[i] * w[k + i] for i in range(d))
            if sp.simplify(r) != 0:
                ok = False
                break
        if ok:
            return d, [c[i] for i in range(d)], True
    return None, None, False


def prony_extract(a_exact, s, mlo, mhi):
    """Blind extraction of (lambda_j, q_j) from the exact coefficient sequence,
    assuming only the integer type shift s found by blind_type_and_scale.

        v_m := a_m (-1)^{m-1} / (m-s)!  = sum_j q_j x_j^{m-s+1},  x_j = 1/lambda_j
    """
    w = []
    for m in range(mlo, mhi + 1):
        w.append(sp.together(a_exact[m] * (-1) ** (m - 1) / sp.factorial(m - s)))
    d, c, ok = prony_exact(w)
    if d is None:
        return None
    x = sp.Symbol("x")
    P = x ** d + sum(c[i] * x ** i for i in range(d))
    roots = sp.roots(sp.Poly(P, x))
    xs = []
    for r, mult in roots.items():
        xs += [sp.radsimp(r)] * mult
    # Vandermonde solve for the charges, on the FIRST d rows
    k0 = mlo - s + 1
    V = sp.Matrix([[xs[j] ** (k0 + i) for j in range(d)] for i in range(d)])
    q = V.solve(sp.Matrix(w[:d]))
    lam = [sp.radsimp(1 / xi) for xi in xs]
    # residual on the whole window
    resid = sp.Integer(0)
    for i, m in enumerate(range(mlo, mhi + 1)):
        pred = sum(q[j] * xs[j] ** (m - s + 1) for j in range(d))
        resid += sp.Abs(sp.simplify(pred - w[i]))
    return lam, [q[j] for j in range(d)], sp.simplify(resid), d, ok


# ----------------------------------------------------------------------- main

def main():
    mp.mp.dps = DPS
    os.makedirs("debug/data", exist_ok=True)

    say("=" * 100)
    say("EXCHANGE CLASS -- STEPS 2-5: Borel / Stokes / classification   (mp.dps = %d)" % DPS)
    say("=" * 100)

    # ---------------------------------------------------------------- STEP 2
    say("")
    say("STEP 2 -- the LOCAL Borel transform in closed form, verified by Laplace")
    say("  phi(R) = e^{AR} F(R) = int_0^oo e^{-xi R} Bor(xi) d xi  ;  digits agreed:")
    say("")
    say("  %12s %8s %26s %26s %8s"
        % ("(a,b)", "R", "phi (closed form)", "Laplace of Bor", "digits"))
    lap_min = 10 ** 9
    for al, be in CASES:
        sings, B = bor_closed(al, be)
        brk = sorted({float(sp.Rational(2 * al)), float(sp.Rational(2 * be)),
                      float(sp.Rational(2 * (al + be)))})
        for R in (mp.mpf(1), mp.mpf("2.5"), mp.mpf(7)):
            lhs = phi_numeric(al, be, R)
            pts = [mp.mpf(0)] + [mp.mpf(p) for p in brk] + [mp.inf]
            rhs = mp.quad(lambda xi: mp.e ** (-xi * R) * B(xi), pts)
            dig = -mp.log10(abs(lhs - rhs) / abs(lhs))
            lap_min = min(lap_min, float(dig))
            say("  %12s %8s %26s %26s %8.1f"
                % ("(%s,%s)" % (al, be), mp.nstr(R, 4), mp.nstr(lhs, 20),
                   mp.nstr(rhs, 20), float(dig)))
    say("")
    say("  L8  closed-form Borel transform verified by Laplace, worst agreement "
        "%.1f digits : %s" % (lap_min, "PASS" if lap_min > 25 else "FAIL"))

    # ---------------------------------------------------------------- STEP 3
    say("")
    say("=" * 100)
    say("STEP 3 -- Borel singularity table + generalized Stokes data")
    say("=" * 100)
    say("")
    say("  local xi_j = global zeta_j - A.  Local form q_j (xi-xi_j) ln(xi-xi_j):")
    say("  monodromy xi-xi_j -> e^{2 pi i}(xi-xi_j) gives Disc Bor = 2 pi i q_j (xi-xi_j)/(2ab),")
    say("  i.e. the alien derivative attaches the SINGLE term 2 pi i q_j /(2ab R^2).")
    say("")
    say("  %12s %10s %12s %8s %22s %10s"
        % ("(a,b)", "xi_j", "zeta_j", "q_j", "type", "S/(2 pi i)"))
    for al, be in CASES:
        sings, _ = bor_closed(al, be)
        A = al + be
        for xj, qj in sings:
            typ = "log branch (u ln u)" + (" [AT ORIGIN]" if xj == 0 else "")
            say("  %12s %10s %12s %8s %22s %10s"
                % ("(%s,%s)" % (al, be), xj, A + xj, qj, typ,
                   sp.nsimplify(qj / (2 * al * be))))

    # ---------------- route (ii): BLIND large order, with synthetic validation
    say("")
    say("  ROUTE (ii) BLIND large-order extraction from the coefficient sequence only.")
    say("  SYNTHETIC VALIDATION FIRST (known log-Borel data, incl. irrational charges):")
    say("")
    say("  %-46s %10s %10s %s" % ("synthetic model", "s (blind)", "lam_min", "recovered (lam_j, q_j)"))
    syn = [
        ("q(xi+L)ln(xi+L): L={2,5,9}, q={1,-3/7,sqrt(2)}", 3,
         [(sp.Integer(2), sp.Integer(1)), (sp.Integer(5), sp.Rational(-3, 7)),
          (sp.Integer(9), sp.sqrt(2))]),
        ("q ln(xi+L): L={3,8}, q={2,-5}  [pure log, s=2]", 2,
         [(sp.Integer(3), sp.Integer(2)), (sp.Integer(8), sp.Integer(-5))]),
        ("q/(xi+L): L={4,11}, q={1,7/5}  [simple pole, s=1]", 1,
         [(sp.Integer(4), sp.Integer(1)), (sp.Integer(11), sp.Rational(7, 5))]),
    ]
    syn_ok = True
    for name, s_true, data in syn:
        M = 46
        aa = {}
        for m in range(1, M + 1):
            tot = sp.Integer(0)
            for L, q in data:
                if s_true == 3 and m >= 3:
                    tot += q * (-1) ** (m - 1) * sp.factorial(m - 3) * L ** (2 - m)
                elif s_true == 2 and m >= 2:
                    tot += q * (-1) ** m * sp.factorial(m - 2) * L ** (1 - m)
                elif s_true == 1 and m >= 1:
                    tot += q * (-1) ** (m - 1) * sp.factorial(m - 1) * L ** (-m)
            if tot != 0:
                aa[m] = tot
        s_est, lam_est = blind_type_and_scale(aa, 34)
        s_int = int(mp.nint(s_est))
        got = prony_extract(aa, s_int, max(s_int + 1, 4), max(s_int + 1, 4) + 13)
        pairs = sorted(zip([sp.nsimplify(l) for l in got[0]], got[1]),
                       key=lambda p: float(p[0]))
        exp_pairs = sorted(data, key=lambda p: float(p[0]))
        sgn = (-1) ** (s_true + 1)      # v_m convention: qhat = (-1)^{s+1} q
        match = all(sp.simplify(g[0] - e[0]) == 0
                    and sp.simplify(g[1] - sgn * e[1]) == 0
                    for g, e in zip(pairs, exp_pairs)) and s_int == s_true
        syn_ok &= match
        say("  %-46s %10s %10s %s  -> %s"
            % (name, mp.nstr(s_est, 6), mp.nstr(lam_est, 6),
               [(str(l), str(q)) for l, q in pairs], "OK" if match else "MISMATCH"))
    say("")
    say("  L9  blind extractor validated on synthetic log-Borel data (type s AND")
    say("      charges recovered exactly, including irrational q)               : %s"
        % ("PASS" if syn_ok else "FAIL"))

    say("")
    say("  Now BLIND on the real object (coefficients taken from the PARSED closed form):")
    say("")
    say("  %12s %10s %10s %28s %16s"
        % ("(a,b)", "s (blind)", "lam_min", "blind (lambda_j, q_j)", "Prony residual"))
    blind_ok = True
    blind_store = {}
    for al, be in CASES:
        A, aa, bb = sector_and_coeffs(ordered_xi_closed(al * R_s, be * R_s), 46)
        aa = {m: v for m, v in aa.items() if not v.has(sp.EulerGamma) and not v.has(sp.log)}
        s_est, lam_est = blind_type_and_scale(aa, 34)
        s_int = int(mp.nint(s_est))
        lam, q, resid, d, ok = prony_extract(aa, s_int, 4, 24)
        pairs = sorted(zip(lam, q), key=lambda p: float(p[0]))
        exp = sorted([(2 * al, sp.Rational(1, 2 * al * be)),
                      (2 * be, sp.Rational(1, 2 * al * be)),
                      (2 * (al + be), sp.Rational(-1, 2 * al * be))],
                     key=lambda p: float(p[0]))
        # merge duplicates (degenerate 2a == 2b)
        merged = {}
        for L, qq in exp:
            merged[L] = merged.get(L, 0) + qq
        exp = sorted(merged.items(), key=lambda p: float(p[0]))
        match = (s_int == 3 and resid == 0 and len(pairs) == len(exp)
                 and all(sp.simplify(g[0] - e[0]) == 0 and sp.simplify(g[1] - e[1]) == 0
                         for g, e in zip(pairs, exp)))
        blind_ok &= match
        blind_store[(al, be)] = pairs
        say("  %12s %10s %10s %28s %16s  %s"
            % ("(%s,%s)" % (al, be), mp.nstr(s_est, 8), mp.nstr(lam_est, 8),
               [(str(l), str(qq)) for l, qq in pairs], resid,
               "OK" if match else "MISMATCH"))
    say("")
    say("  L10 blind large-order route reproduces route (i) EXACTLY (residual 0)  : %s %d/%d"
        % ("PASS" if blind_ok else "FAIL", len(CASES), len(CASES)))
    say("      type shift s = 3 blind  =>  singularity is (xi-xi_j) ln(xi-xi_j),")
    say("      NOT a pole (s=1) and NOT a bare log (s=2).")

    # ---------------- route (iii): the R-plane branch cut
    say("")
    say("  ROUTE (iii) exact discontinuity across arg R = pi (the E_1 and ln R cuts).")
    say("  Prediction:  Disc F = -(2 pi i)/(2ab R^2) * sum_j q_j e^{-zeta_j R}.")
    say("")
    say("  %12s %8s %30s %30s %10s"
        % ("(a,b)", "x", "Disc F  (numeric)", "prediction", "rel resid"))
    disc_worst = 0.0
    off = mp.mpf(10) ** (-60)
    for al, be in CASES:
        a, b = mp.mpf(sp.Rational(al)), mp.mpf(sp.Rational(be))
        A = a + b
        kap = 2 * a * b / A

        def F(R):
            return (1 / (2 * a * b * R ** 2)) * (
                mp.e ** (-A * R) * (mp.euler + mp.log(kap * R))
                + mp.e ** ((a - b) * R) * mp.e1(2 * a * R)
                + mp.e ** ((b - a) * R) * mp.e1(2 * b * R)
                - mp.e ** (A * R) * mp.e1(2 * A * R))

        for x in (mp.mpf("1.3"), mp.mpf("2.7")):
            Rp, Rm = -x + 1j * off, -x - 1j * off
            disc = F(mp.mpc(Rp)) - F(mp.mpc(Rm))
            Rc = mp.mpc(-x)
            pred = -(2j * mp.pi) / (2 * a * b * Rc ** 2) * (
                (-1) * mp.e ** (-A * Rc) + mp.e ** (-(b - a) * Rc)
                + mp.e ** (-(a - b) * Rc) + (-1) * mp.e ** (A * Rc))
            rel = abs(disc - pred) / abs(pred)
            disc_worst = max(disc_worst, float(rel))
            say("  %12s %8s %30s %30s %10s"
                % ("(%s,%s)" % (al, be), mp.nstr(x, 3), mp.nstr(disc, 14),
                   mp.nstr(pred, 14), mp.nstr(rel, 4)))
    say("")
    say("  L11 branch-cut route matches the charge/position table, worst rel resid "
        "%.2e : %s" % (disc_worst, "PASS" if disc_worst < 1e-40 else "FAIL"))

    # ---------------------------------------------------------------- STEP 4
    say("")
    say("=" * 100)
    say("STEP 4 -- classification")
    say("=" * 100)
    say("")
    say("  (a) FIELD OF THE STOKES DATA")
    say("      S_j = 2 pi i * q_j / (2ab),  q_j in {-1,+1,+1,-1} (INTEGERS).")
    say("      => S_j in 2 pi i * Q(rates).  ALGEBRAIC up to a single power of pi.")
    say("      No radical is needed at all here (the hybrid class needed Q(sqrt d)).")
    say("")
    say("  (b) IS THE BOUNDARY DATUM FORCED BY THE SKELETON?")
    say("      Boundary datum = the single number (gamma + ln kappa) multiplying")
    say("      e^{-AR}/(2ab R^2) alongside ln R.  Forcing law:")
    say("           ln kappa = sum_{j != 0} q_j ln lambda_j      (charge-weighted)")
    say("      and charge neutrality sum_j q_j = 0 is what removes the regular")
    say("      xi-linear remainder from Bor.  Test as an exact identity + decoys:")
    say("")
    say("  %12s %12s %14s %10s %30s"
        % ("(a,b)", "kappa", "prod lam^q", "equal?", "Bor xi-linear remainder"))
    xi = sp.Symbol("xi", positive=True)
    force_ok = True
    for al, be in CASES:
        A = al + be
        kap = 2 * al * be / A
        prod = (2 * al) ** 1 * (2 * be) ** 1 * (2 * A) ** (-1)
        # exact symbolic Bor, then subtract the pure charge sum -> must be a CONSTANT
        Bsym = (-xi * sp.log(xi / kap)
                + (xi + 2 * al) * sp.log((xi + 2 * al) / (2 * al))
                + (xi + 2 * be) * sp.log((xi + 2 * be) / (2 * be))
                - (xi + 2 * A) * sp.log((xi + 2 * A) / (2 * A)))
        pure = (-xi * sp.log(xi)
                + (xi + 2 * al) * sp.log(xi + 2 * al)
                + (xi + 2 * be) * sp.log(xi + 2 * be)
                - (xi + 2 * A) * sp.log(xi + 2 * A))
        rem = sp.expand(Bsym - pure)
        lin = sp.simplify(sp.expand_log(sp.together(sp.diff(rem, xi)), force=True))
        # independent check: rem must take the SAME value at two rational points
        r1 = sp.simplify(sp.expand_log(rem.subs(xi, sp.Rational(1, 3)), force=True))
        r2 = sp.simplify(sp.expand_log(rem.subs(xi, sp.Rational(11, 7)), force=True))
        flat = sp.simplify(sp.expand_log(r1 - r2, force=True)) == 0
        eq = sp.simplify(kap - prod) == 0
        force_ok &= eq and flat
        say("  %12s %12s %14s %10s %30s"
            % ("(%s,%s)" % (al, be), kap, sp.cancel(prod), eq,
               "d/dxi|_{xi=1/3} = %s ; xi-flat(exact) = %s"
               % (mp.nstr(_mpf(lin.subs(xi, sp.Rational(1, 3))), 3), flat)))
    say("")
    say("  L12 kappa = prod lambda_j^{q_j} AND Bor - (pure charge sum) is xi-CONSTANT")
    say("      (=> no free regular linear term; gamma+ln kappa is forced)        : %s %d/%d"
        % ("PASS" if force_ok else "FAIL", len(CASES), len(CASES)))

    say("")
    say("  (c) CUT-FREENESS / BOREL SUMMABILITY")
    say("      Bor(xi) is analytic on the whole Laplace ray (0, oo): the three")
    say("      non-trivial singularities sit at xi = -2a, -2b, -2A < 0 and the")
    say("      fourth is the ENDPOINT xi = 0.  So the object is Borel summable")
    say("      along R > 0 with NO lateral ambiguity -- unlike the hybrid class,")
    say("      whose sector constants had to cancel pairwise to achieve this.")
    say("      BUT: Disc_{arg R = pi} F != 0 (L11), so F is genuinely MULTIVALUED")
    say("      in complex R -- the hybrid ERI was single-valued.  That is exactly")
    say("      the price of the R-dependent log, and the discontinuity it creates")
    say("      is 2 pi i x (rational) x elementary -- no new transcendental.")

    say("")
    say("  DECOY CONTROLS on the charge table (each MUST fail):")
    bad = [("all +1", (1, 1, 1)), ("swap sign on 2a", (-1, 1, -1)),
           ("swap sign on 2b", (1, -1, -1)), ("charge 2 on 2A", (1, 1, -2)),
           ("charges (2,1,-1)", (2, 1, -1)), ("charges (0,1,-1)", (0, 1, -1))]
    say("    %-22s %s" % ("decoy charge vector", "cases matching kappa (must be < 6/6)"))
    dec_ok = True
    for nm, ch in bad:
        hits, where = 0, []
        for al, be in CASES:
            A = al + be
            k = (2 * al) ** ch[0] * (2 * be) ** ch[1] * (2 * A) ** ch[2]
            if sp.simplify(k - 2 * al * be / A) == 0:
                hits += 1
                where.append("(%s,%s)" % (al, be))
        dec_ok &= (hits < len(CASES))
        say("    %-22s %d/%d %s" % (nm, hits, len(CASES),
                                    ("  accidental at " + ",".join(where)) if where else ""))
    say("")
    say("    L12c every decoy charge vector FAILS somewhere                       : %s"
        % ("PASS" if dec_ok else "FAIL"))
    say("    NOTE (honest): the 'swap sign on 2a' decoy accidentally matches at")
    say("    (a,b) = (1/2,7/2) because there 2a = 1 and 1^{+-1} = 1 -- a degenerate")
    say("    rate, not a law.  It fails on the other five, which is why the sweep")
    say("    is run over six rate pairs rather than one.")

    say("")
    say("  Lateral-Laplace check (numerical proof of no ambiguity on the ray R > 0):")
    say("  %12s %8s %14s %14s %12s"
        % ("(a,b)", "R", "theta", "|L_+ - L_-|", "rel"))
    lat_worst = 0.0
    for al, be in CASES[:3]:
        sings, B = bor_closed(al, be)

        def Bc(z, al=al, be=be):
            A = sp.Rational(al + be)
            kap = 2 * al * be / (al + be)
            z = mp.mpc(z)
            tot = -z * mp.log(z / mp.mpf(sp.Rational(kap)))
            for lam, sg in ((2 * al, 1), (2 * be, 1), (2 * (al + be), -1)):
                L = mp.mpf(sp.Rational(lam))
                tot += sg * (z + L) * mp.log((z + L) / L)
            return tot / mp.mpf(sp.Rational(2 * al * be))

        for R in (mp.mpf(2), mp.mpf(5)):
            for th in (mp.pi / 6, mp.pi / 4):
                Lp = mp.quad(lambda u: mp.e ** (-u * mp.e ** (1j * th) * R)
                             * Bc(u * mp.e ** (1j * th)) * mp.e ** (1j * th),
                             [0, 1, 10, mp.inf])
                Lm = mp.quad(lambda u: mp.e ** (-u * mp.e ** (-1j * th) * R)
                             * Bc(u * mp.e ** (-1j * th)) * mp.e ** (-1j * th),
                             [0, 1, 10, mp.inf])
                rel = abs(Lp - Lm) / abs(Lp)
                lat_worst = max(lat_worst, float(rel))
                say("  %12s %8s %14s %14s %12s"
                    % ("(%s,%s)" % (al, be), mp.nstr(R, 3), mp.nstr(th, 6),
                       mp.nstr(abs(Lp - Lm), 4), mp.nstr(rel, 4)))
    say("")
    say("  L12b lateral resummations agree (no Stokes ambiguity on R > 0), worst")
    say("       relative gap %.2e                                             : %s"
        % (lat_worst, "PASS" if lat_worst < 1e-30 else "FAIL"))

    say("")
    say("  (d) WHERE THE ln R SITS IN THE EXPANSION (which orders R^-m carry it)")
    say("  %12s %26s %26s"
        % ("(a,b)", "orders m with ln R", "orders m with a_m (power)"))
    for al, be in CASES:
        A, aa, bb = sector_and_coeffs(ordered_xi_closed(al * R_s, be * R_s), 10)
        say("  %12s %26s %26s"
            % ("(%s,%s)" % (al, be), sorted(bb), sorted(aa)))
    say("      => the log-monomial appears at EXACTLY ONE order (m = 2, the leading")
    say("         one): a single weight-1 log trans-monomial, no (ln R)^2 anywhere.")

    # ---------------------------------------------------------------- STEP 4b
    say("")
    say("=" * 100)
    say("STEP 4b -- genericity: the general (tau, sigma, H, j) exchange kernel")
    say("=" * 100)
    GEN = [(0, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0), (1, 1, 1, 1, 0, 0),
           (2, 0, 1, 0, 1, 0), (2, 2, 2, 2, 0, 1), (3, 1, 1, 2, 1, 1)]
    al, be = sp.Rational(3, 2), sp.Integer(1)
    A = al + be
    allowed = {2 * al, 2 * be, 2 * A}
    say("")
    say("  %22s %26s %18s %18s"
        % ("(tau,sig,H1,H2,j1,j2)", "E1 rates lambda", "R powers on E1", "R powers on lnR"))
    gen_pos_ok = True
    for prm in GEN:
        terms = parse(ordered_xi_general(*prm, al * R_s, be * R_s))
        lam = sorted({t.lam for t in terms if t.lam is not None}, key=lambda z: float(z))
        ke = sorted({int(t.k) for t in terms if t.kind == "E1"})
        kl = sorted({int(t.k) for t in terms if t.lnR})
        gen_pos_ok &= set(lam) <= allowed
        say("  %22s %26s %18s %18s"
            % (str(prm), [str(x) for x in lam], ke, kl))
    say("")
    say("  L13 general kernel's Borel singularity POSITIONS stay in {2a, 2b, 2A}")
    say("      (i.e. global zeta in {+-A, +-(a-b)}) for every tested (tau,sigma,H,j) : %s"
        % ("PASS" if gen_pos_ok else "FAIL"))
    say("      R powers k <= -2 on the E_1 terms => (xi+lam)^{|k|-1} ln(xi+lam)")
    say("      higher log-branch orders; k = -1 would be a simple pole.  Both are")
    say("      2 pi i x rational monodromy, so the classification is unchanged.")

    with open("debug/data/exch_gamma_borel.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
