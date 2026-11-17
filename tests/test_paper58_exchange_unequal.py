"""Paper 58 -- the 1s two-centre exchange integral at UNEQUAL exponents (records EX + EXb).

Guard for the "[MEASURED] At unequal exponents the 1s exchange integral has no closed form" paragraph of
Paper 58 (papers/group2_quantum_chemistry/paper_58_abelian_residue.tex), records
docs/predictions/2026-11-16_exchange_unequal_zeta.md (EX, VOID on CRT headroom at (1, 3/2), every clause held) and
docs/predictions/2026-11-16_exchange_unequal_zeta_b.md (EXb, GO; joint verdict EX + EXb = GO for inventory row 3c).
Written as a SEPARATE pass from the derivation and scoring (CLAUDE.md Sec. 9).

f(R) = (1s_A 1s_B | 1s_A 1s_B), Slater 1s exponent a on A, b on B, (a, b) = (1, 3/2) and (1, 2).  Claim: an operator
L = sum_{i<=9} p_i(R) D^i, deg p_i <= 16, over Q, with

  * L f = 0 (exact, on the Neumann log-Taylor series at R = 0);
  * its Frobenius solution from the data at the R = 0 exponents is f (checked against the production route);
  * every finite nonzero singular point apparent;
  * R = 0 Fuchsian, exponents {-1 (x2), 0, 1, 2, 3, 4, 6, 9}; f's ln R part starts at R^9 ln R, coefficient
    (8/3375) a^5 b^5 = 9/500 and 256/3375;
  * rates at infinity {-(a+b), 2a, 2b, 0 (x2), -2a (x2), -2b (x2)}.

PROVENANCE OF THE FROZEN OPERATORS.  `OP_EXCH_1_3HALF` is `operator_Q_crt` of debug/data/exch_heldout_1_3_2_EXB.json
(record EXb: CRT + rational reconstruction over 32 primes, 831 bits; sha256 of the json d0bf8b0b...8d22) and
`OP_EXCH_1_2` is `operator_Q_crt` of debug/data/exch_heldout_1_2_EXS.json (record EX: four primes, 103 bits; sha256
05988bd9...82f5), copied verbatim (integer polynomials, content removed by the driver).  Nothing here reads debug/:
the series is rebuilt in-test by a Q-only port of the frozen generator debug/exch_gen.py (`generate`, record EX sha256
5a17efc6...96d5), and the numeric check uses geovac.qfd_core.exchange_hp, an independent production route.

What each test rejects (the wrong answer is named in each docstring; plants at the end of the file):

* annihilation: a stored operator that does not kill f's exact series (any coefficient off, the wrong pair's operator,
  the equal-exponent order-4 operator), or a series generator with a or b wrong;
* R = 0: a wrong exponent set (e.g. 5 for 6, a single -1), a non-Fuchsian point, a log-free local solution space;
* infinity: a wrong rate multiset (e.g. an a - b rate, unsplit rates), a rank other than one;
* ln R: an onset other than R^9, a coefficient other than (8/3375) a^5 b^5 -- read both from the series and as the
  value FORCED by the operator from the non-log data below R^9;
* apparent points: a genuine singular point (non-integer or colliding exponent, or a log solution) at a root of the
  degree-11 factor of p_9, read mod primes where that factor has a linear factor (the record's own reading; the
  factor is irreducible over Q);
* numeric: a Frobenius solution that is not the exchange integral.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from typing import Dict, List, Tuple

import numpy as np
import pytest
import sympy as sp

R = sp.Symbol("R")

# ---------------------------------------------------------------- frozen operators
# Coefficient lists [p_0, p_1, ..., p_9] of  L = sum_i p_i(R) D^i.

OP_EXCH_1_3HALF = [
    "72*(200475000000*R**14 - 6875204062500*R**13 - 180159946593750*R**12 - 1229894981681250*R**11 - 4941251945656875*R**10 - 13167040488458750*R**9 - 24934634749690500*R**8 - 35125795088488000*R**7 - 36778160139208200*R**6 - 29039006176905720*R**5 - 17401845162786696*R**4 - 7907628289573440*R**3 - 2629512797006592*R**2 - 572694479884800*R - 53320646615040)",
    "-24*(601425000000*R**15 - 20625612187500*R**14 - 542217289781250*R**13 - 3627956988168750*R**12 - 13329788433470625*R**11 - 30983628923538750*R**10 - 49074520366783375*R**9 - 54015152733520875*R**8 - 35432732764931600*R**7 - 4773000549615660*R**6 + 14189065299239412*R**5 + 13385949966465300*R**4 + 6007737229345440*R**3 + 2310800758680960*R**2 + 984312426601728*R + 153090362880000)",
    "6*(801900000000*R**16 - 28489826250000*R**15 - 691594992000000*R**14 - 3864392338612500*R**13 - 9871695584707500*R**12 - 8280982721098750*R**11 + 24717017688967125*R**10 + 105363904676014500*R**9 + 216563000634717250*R**8 + 294813711073636560*R**7 + 280436006512528608*R**6 + 185155514172876960*R**5 + 78684081683455728*R**4 + 14819482371450240*R**3 - 3620654698848768*R**2 - 2528096925312000*R - 533206466150400)",
    "6*(989010000000*R**16 - 32519694375000*R**15 - 931731674362500*R**14 - 6907515790920000*R**13 - 27276221279061250*R**12 - 71158059811152875*R**11 - 132659552222280250*R**10 - 182205211877194050*R**9 - 183066644395126440*R**8 - 125067419246700392*R**7 - 41221873991751480*R**6 + 20306850107832336*R**5 + 36619920210525120*R**4 + 22799029301116416*R**3 + 7792904356715520*R**2 + 1941544147931136*R + 188418908160000)",
    "R*(668250000000*R**15 - 19170691875000*R**14 - 694686866250000*R**13 - 6041289724387500*R**12 - 25106895294161875*R**11 - 74841835546147500*R**10 - 176278157991715250*R**9 - 321904322906193000*R**8 - 461962006249297800*R**7 - 522808559631745920*R**6 - 456566447177294256*R**5 - 291726221592902400*R**4 - 127421810735454720*R**3 - 38560773912837120*R**2 - 10049645489135616*R - 1130513448960000)",
    "-R*(1822095000000*R**15 - 61068602812500*R**14 - 1676550771912500*R**13 - 11982651297505625*R**12 - 47913058039447500*R**11 - 126857811875525750*R**10 - 236033276468598000*R**9 - 323237281565705400*R**8 - 330662557655556480*R**7 - 245941535626125264*R**6 - 119805661924431840*R**5 - 23731439888080512*R**4 + 12795115534479360*R**3 + 12114156566734848*R**2 + 4738925584834560*R + 744908089540608)",
    "-R**2*(735075000000*R**14 - 24465096562500*R**13 - 673550136937500*R**12 - 4673788295734375*R**11 - 17801235220633125*R**10 - 46175290131157500*R**9 - 86955286089320000*R**8 - 123622710688075500*R**7 - 136478451267393900*R**6 - 118962226507726080*R**5 - 82136954585055744*R**4 - 44413952170141440*R**3 - 18155179952750592*R**2 - 5115763401154560*R - 744908089540608)",
    "R**3*(49005000000*R**13 - 1638728437500*R**12 - 46737432962500*R**11 - 366582506194375*R**10 - 1487156521940000*R**9 - 3651174464565500*R**8 - 6078385078122000*R**7 - 7308600999024300*R**6 - 6537832385325120*R**5 - 4451038340538816*R**4 - 2446091479975680*R**3 - 1105706037325824*R**2 - 310378370641920*R - 41383782752256)",
    "R**4*(66825000000*R**12 - 2247184687500*R**11 - 60790383187500*R**10 - 405112436396875*R**9 - 1400122490731250*R**8 - 3141615710512500*R**7 - 5016404375967000*R**6 - 5916170695639500*R**5 - 5186909374021800*R**4 - 3333207175424640*R**3 - 1521333629194752*R**2 - 457260125253120*R - 67248646972416)",
    "2*R**5*(4455000000*R**11 - 151000312500*R**10 - 3992292087500*R**9 - 24878273313125*R**8 - 76755983840000*R**7 - 148036260295500*R**6 - 196259782122000*R**5 - 185067612112500*R**4 - 123712823733120*R**3 - 57263380050816*R**2 - 17435951205120*R - 2586486422016)",
]

OP_EXCH_1_2 = [
    "384*(69672960*R**14 - 2694083328*R**13 - 29956029120*R**12 - 148887629280*R**11 - 466729434648*R**10 - 1003963993560*R**9 - 1565271506208*R**8 - 1832839402176*R**7 - 1618827097900*R**6 - 1088301620502*R**5 - 558794688435*R**4 - 215596001268*R**3 - 58966431048*R**2 - 10330662888*R - 746746560)",
    "-96*(278691840*R**15 - 10776333312*R**14 - 120520846080*R**13 - 567854893440*R**12 - 1594375746912*R**11 - 2935822330080*R**10 - 3720936308232*R**9 - 3210554266524*R**8 - 1507051780960*R**7 + 195415986852*R**6 + 853018254470*R**5 + 617739228873*R**4 + 251193073788*R**3 + 82440679608*R**2 + 25932788280*R + 3225096000)",
    "16*(557383680*R**16 - 22156498944*R**15 - 219183197184*R**14 - 821838911616*R**13 - 1454341260096*R**12 - 298591472592*R**11 + 4751245206516*R**10 + 13160560943880*R**9 + 20987098707962*R**8 + 23120073164808*R**7 + 18149890922592*R**6 + 9981063923724*R**5 + 3555402388839*R**4 + 573693629076*R**3 - 103389745032*R**2 - 63233804088*R - 11201198400)",
    "16*(603832320*R**16 - 22555224576*R**15 - 286175251584*R**14 - 1462704317568*R**13 - 4459784794992*R**12 - 9325917716364*R**11 - 14164825761756*R**10 - 15984765999322*R**9 - 13235654711628*R**8 - 7393749472674*R**7 - 1807115269401*R**6 + 1139878966569*R**5 + 1491098942034*R**4 + 768645653688*R**3 + 221425089480*R**2 + 45904927824*R + 3870115200)",
    "8*R*(69672960*R**15 - 2084444928*R**14 - 46711264896*R**13 - 289438695360*R**12 - 1014567754644*R**11 - 2634015119976*R**10 - 5300302924830*R**9 - 8198987234508*R**8 - 9926222573690*R**7 - 9479532896352*R**6 - 7003658424369*R**5 - 3811571354232*R**4 - 1434154020984*R**3 - 373166182872*R**2 - 80608657248*R - 7740230400)",
    "-24*R*(110315520*R**15 - 4196926656*R**14 - 49352654592*R**13 - 243371242548*R**12 - 733573115928*R**11 - 1519318534014*R**10 - 2275036995456*R**9 - 2550212393786*R**8 - 2160947131920*R**7 - 1347081899149*R**6 - 560941193601*R**5 - 104907546090*R**4 + 33751543968*R**3 + 30308930856*R**2 + 10127536272*R + 1340288208)",
    "-2*R**2*(452874240*R**14 - 17145758592*R**13 - 200800847232*R**12 - 957652529400*R**11 - 2777904247788*R**10 - 5656766413896*R**9 - 8539786060008*R**8 - 9873098300430*R**7 - 8942359814161*R**6 - 6430354470144*R**5 - 3676626348768*R**4 - 1650624834456*R**3 - 561271821696*R**2 - 131850742464*R - 16083458496)",
    "6*R**3*(5806080*R**13 - 219668544*R**12 - 2879983296*R**11 - 15828684876*R**10 - 49872521808*R**9 - 100193004036*R**8 - 139704592536*R**7 - 142619149931*R**6 - 109520387928*R**5 - 64248505532*R**4 - 30035999880*R**3 - 11384698944*R**2 - 2680576416*R - 297841824)",
    "9*R**4*(5806080*R**12 - 221281344*R**11 - 2541783168*R**10 - 11535214236*R**9 - 30512541372*R**8 - 54696822852*R**7 - 71219148404*R**6 - 69334783167*R**5 - 50599461551*R**4 - 27206250120*R**3 - 10404758512*R**2 - 2617282584*R - 322661976)",
    "R**5*(5806080*R**11 - 222571584*R**10 - 2467592640*R**9 - 10438506396*R**8 - 24713371152*R**7 - 38221242084*R**6 - 41491515672*R**5 - 32453435903*R**4 - 18146025648*R**3 - 7043999400*R**2 - 1795425912*R - 223381368)",
]
# The EQUAL-exponent operator (zeta = 1), order 4, degree 8: OP_EXCH_1 of tests/test_paper58_r_operators.py,
# verbatim.  Used only as a negative control (and as a fire-test substitute).
OP_EXCH_1_EQUAL = [
    "24*R*(32*R**6 + 192*R**5 + 576*R**4 + 1380*R**3 + 2070*R**2 + 1890*R + 945)",
    "-16*(16*R**8 + 80*R**7 + 72*R**6 - 246*R**5 - 1035*R**4 - 2430*R**3 - 3915*R**2 - 4050*R - 2025)",
    "-2*R*(64*R**7 + 768*R**6 + 2976*R**5 + 5400*R**4 + 6480*R**3 + 5940*R**2 + 4050*R + 2025)",
    "2*R**2*(32*R**6 + 64*R**5 - 480*R**4 - 1980*R**3 - 3690*R**2 - 4050*R - 2025)",
    "R**3*(32*R**5 + 240*R**4 + 720*R**3 + 1260*R**2 + 1350*R + 675)",
]

PAIRS = {
    "1_3half": (Fraction(1), Fraction(3, 2), OP_EXCH_1_3HALF),
    "1_2": (Fraction(1), Fraction(2), OP_EXCH_1_2),
}
LN_COEF = Fraction(8, 3375)                      # (8/3375) a^5 b^5 at the R^9 ln R onset
EXPONENTS_0 = {-1: 2, 0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 6: 1, 9: 1}


def rates_inf(a: Fraction, b: Fraction) -> Dict[sp.Rational, int]:
    q = lambda x: sp.Rational(x.numerator, x.denominator)
    return {q(-(a + b)): 1, q(2 * a): 1, q(2 * b): 1, sp.Integer(0): 2, q(-2 * a): 2, q(-2 * b): 2}


@lru_cache(maxsize=None)
def op_coeffs(key: str) -> Tuple[Tuple[Fraction, ...], ...]:
    """p_i as ascending integer coefficient tuples."""
    out = []
    for s in PAIRS[key][2]:
        P = sp.Poly(sp.sympify(s, locals={"R": R}), R)
        cs = [Fraction(0)] * (P.degree() + 1)
        for (j,), c in P.terms():
            cs[j] = Fraction(int(c.p), int(c.q))
        out.append(tuple(cs))
    return tuple(out)


def op_polys(key: str) -> List[sp.Expr]:
    return [sp.expand(sp.sympify(s, locals={"R": R})) for s in PAIRS[key][2]]


# ---------------------------------------------------------------- series generator (Q-only port of debug/exch_gen.py)
#
#   f(R) = (a^3 b^3 / 8) R^5 sum_tau (2 tau + 1) [beta_0^2 X22 - 2 beta_0 beta_2 X20 + beta_2^2 X00],
#   p = (a+b) R/2, q = (a-b) R/2; beta_{tau,k}(q) from Legendre moments; X(p) from tau-independent Laplace tables of
#   poly x {1, ln v, ln(v+2), ln(v+4), ln 2}.  Components (const, m): const in {'1', 'g' (gamma_E), 'l<prime>'},
#   m = power of ln R.  Truncation tau <= N/2 rests on val_R(tau-term) >= 2 tau, checked per tau (val_ok).


def _vec(n: int) -> np.ndarray:
    a = np.empty(n, dtype=object)
    a[:] = [Fraction(0)] * n
    return a


def _mat(n: int, m: int) -> np.ndarray:
    a = np.empty((n, m), dtype=object)
    for i in range(n):
        a[i, :] = [Fraction(0)] * m
    return a


def _stack(rows) -> np.ndarray:
    M = np.empty((len(rows), len(rows[0])), dtype=object)
    for i, v in enumerate(rows):
        M[i, :] = v
    return M


def legendre_PV(tmax: int):
    """P_n and V_n (Q_n = P_n Q_0 - V_n) coefficient lists, n = 0..tmax."""
    P = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    V = [[Fraction(0)], [Fraction(1)]]
    for n in range(1, tmax):
        for L in (P, V):
            a, b = L[n], L[n - 1]
            new = [Fraction(0)] * (len(a) + 1)
            for i, c in enumerate(a):
                new[i + 1] += Fraction(2 * n + 1, n + 1) * c
            for i, c in enumerate(b):
                new[i] -= Fraction(n, n + 1) * c
            L.append(new)
    return P[: tmax + 1], V[: tmax + 1]


def leg_moment(n: int, tau: int) -> Fraction:
    """int_{-1}^1 x^n P_tau(x) dx."""
    if n < tau or (n - tau) % 2:
        return Fraction(0)
    return Fraction(2 ** (tau + 1) * factorial(n) * factorial((n + tau) // 2),
                    factorial((n - tau) // 2) * factorial(n + tau + 1))


class Tables:
    """int_0^oo e^{-pv} v^n B(v) dv as Laurent arrays in p, n = 0..nmax, exponents emin..emax."""

    def __init__(self, nmax: int, emax: int):
        self.nmax, self.emax = nmax, emax
        self.emin = -(nmax + 1)
        W = emax - self.emin + 1
        self.W = W
        fact = [factorial(n) for n in range(nmax + 1)]
        H = [Fraction(0)]
        for n in range(1, nmax + 1):
            H.append(H[-1] + Fraction(1, n))
        powT, lnv1 = _mat(nmax + 1, W), _mat(nmax + 1, W)
        for n in range(nmax + 1):
            i = -n - 1 - self.emin
            powT[n, i] = Fraction(fact[n])
            lnv1[n, i] = Fraction(fact[n]) * H[n]
        neg = -powT
        self.T = {"pow": {"1": powT}, "lnv": {"1": lnv1, "g": neg, "lp": neg.copy()},
                  "l2": {"l2": powT.copy()}, "a2": self._tn(2, 1), "a4": self._tn(4, 2)}

    def _tn(self, a: int, lmult: int):
        """T_n^a = int e^{-pv} v^n ln(v+a); T_0 = (1/p)[ln a + e^{ap} E1(ap)]; T_n = (-d/dp)^n T_0."""
        nmax, emin = self.nmax, self.emin
        top = self.emax + nmax + 2
        L = top + 2
        ea = [Fraction(a) ** k / factorial(k) for k in range(L + 1)]
        S = [Fraction(0)] + [Fraction(-a) ** k / (k * factorial(k)) for k in range(1, L + 1)]
        eS = [sum(ea[i] * S[k - i] for i in range(k + 1)) for k in range(L + 1)]
        Wx = top - emin + 1

        def arr(series_nonneg, shift):
            v = _vec(Wx)
            for k, c in enumerate(series_nonneg):
                e = k + shift
                if emin <= e <= top and c != 0:
                    v[e - emin] = Fraction(c)
            return v
        one = arr([-c for c in eS], -1)
        g = arr([-c for c in ea], -1)
        la = arr([Fraction(1) - ea[0]] + [-c for c in ea[1:]], -1)
        lp = arr([-c for c in ea], -1)
        em = np.array([Fraction(int(e)) for e in range(emin, top + 1)], dtype=object)

        def d(v):
            w = v * em
            out = _vec(Wx)
            out[:-1] = w[1:]
            return out

        rows = {k: [] for k in ("1", "g", "l2", "lp")}
        cur = {"1": one, "g": g, "la": la, "lp": lp}
        for _n in range(nmax + 1):
            for k in ("1", "g", "lp"):
                rows[k].append(cur[k][: self.W])
            rows["l2"].append(cur["la"][: self.W] * lmult)
            lp_over_p = _vec(Wx)
            lp_over_p[:-1] = cur["lp"][1:]
            cur = {"1": -(d(cur["1"]) + lp_over_p), "g": -d(cur["g"]), "la": -d(cur["la"]),
                   "lp": -d(cur["lp"])}
        return {k: _stack(v) for k, v in rows.items()}


def _biv_eval_uminus1(Tm: np.ndarray) -> np.ndarray:
    """sum_{i,k} T[i,k] u^i (u-1)^k -> coefficients in u."""
    Du, Dy = Tm.shape
    res = _vec(Du + Dy)
    for k in range(Dy - 1, -1, -1):
        sh = _vec(Du + Dy)
        sh[1:] = res[:-1]
        res = sh - res
        res[:Du] = res[:Du] + Tm[:, k]
    return res


def _biv_eval_uhalf(Tm: np.ndarray) -> np.ndarray:
    """sum_{i,k} T[i,k] u^i (u/2)^k -> coefficients in u."""
    Du, Dy = Tm.shape
    res = _vec(Du + Dy)
    w = Fraction(1)
    for k in range(Dy):
        res[k:k + Du] = res[k:k + Du] + Tm[:, k] * w
        w = w * Fraction(1, 2)
    return res


def _shift2(c: np.ndarray) -> np.ndarray:
    """q(u) -> q(2+v) coefficients in v."""
    D = len(c)
    res = _vec(D)
    for i in range(D - 1, -1, -1):
        sh = _vec(D)
        sh[1:] = res[:-1]
        res = sh + res * 2
        res[0] = res[0] + c[i]
    return res


def g_pieces(tau: int, j1: int, j2: int, Pc, Vc):
    """g(2+v) = pow + lnv ln v + a2 ln(v+2) + a4 ln(v+4) + l2 ln 2 (polynomials in v)."""
    piv = np.array([Fraction(x) for x in Pc[tau]], dtype=object)
    viv = np.array([Fraction(x) for x in Vc[tau]], dtype=object)
    deg = tau + j1 + j2
    Du = Dy = deg + 1
    Pxy = _mat(Du, Dy)                                  # P_tau(u - y)
    for m, cm in enumerate(Pc[tau]):
        if cm == 0:
            continue
        for i in range(m + 1):
            Pxy[i, m - i] = Fraction(cm * comb(m, i) * (-1) ** (m - i))

    def mul_umy(A):
        B = _mat(*A.shape)
        B[1:, :] = A[:-1, :]
        B[:, 1:] = B[:, 1:] - A[:, :-1]
        return B

    def mul_y(A, j):
        B = _mat(*A.shape)
        if j:
            B[:, j:] = A[:, :-j]
        else:
            B[:, :] = A
        return B

    def xpow(A, j):
        for _ in range(j):
            A = mul_umy(A)
        return A
    K = mul_y(xpow(Pxy, j1), j2) + mul_y(xpow(Pxy, j2), j1)
    Dy2 = Dy + tau
    M, N = _mat(Du, Dy2), _mat(Du, Dy2)
    for i in range(Du):
        M[i, :] = np.convolve(K[i, :], piv)[:Dy2]
        cv = np.convolve(K[i, :], viv)
        N[i, :len(cv)] = cv[:Dy2]
    Dy3 = Dy2 + 1
    ks = np.arange(Dy2)
    invk1 = np.array([Fraction(1, k + 1) for k in range(Dy2)], dtype=object)
    s = M * invk1[None, :]
    Nint = _mat(Du, Dy3)
    Nint[:, 1:] = N * invk1[None, :]
    out_lc, out_rc = {}, {}
    for c in (+1, -1):
        mc = -c
        L = _mat(Du, Dy3)
        L[:, 1:] = s
        sign_k1 = np.array([Fraction(mc ** (int(k) + 1)) for k in ks], dtype=object)
        L[:, 0] = L[:, 0] - (s * sign_k1[None, :]).sum(axis=1)
        w = s * sign_k1[None, :]
        suf = np.flip(np.cumsum(np.flip(w, axis=1), axis=1), axis=1)
        Rm = _mat(Du, Dy3)
        for j in range(1, Dy3):
            Rm[:, j] = suf[:, j - 1] * Fraction(-1, j) * (mc ** j)
        out_lc[c], out_rc[c] = L, Rm
    srcs = (("L+", out_lc[1]), ("L-", out_lc[-1]), ("R+", out_rc[1]), ("R-", out_rc[-1]), ("N", Nint))
    up = {k: _biv_eval_uminus1(v) for k, v in srcs}
    lo = {k: _biv_eval_uhalf(v) for k, v in srcs}
    h = lambda a: a * Fraction(1, 2)
    pieces = {
        "a2": h(up["L+"]),                              # ln u    = ln(v+2)
        "a4": -h(lo["L+"]),                             # ln(u+2) = ln(v+4)
        "lnv": h(lo["L-"]) - h(up["L-"]),               # ln(u-2) = ln v
        "l2": h(lo["L+"]) - h(lo["L-"]),
        "pow": h(up["R+"] - lo["R+"]) - h(up["R-"] - lo["R-"]) - (up["N"] - lo["N"]),
    }
    out = {}
    for k, v in pieces.items():
        nz = np.nonzero(v != 0)[0]
        top = int(nz[-1]) + 1 if len(nz) else 1
        assert top <= 2 * tau + 6, (k, top, tau)
        out[k] = _shift2(v[:top])
    return out


KERNELS = ((2, 2), (2, 0), (0, 0))


def _prime_logs(x: Fraction) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for pr, e in sp.factorint(x.numerator).items():
        out[int(pr)] = out.get(int(pr), 0) + int(e)
    for pr, e in sp.factorint(x.denominator).items():
        out[int(pr)] = out.get(int(pr), 0) - int(e)
    return {k: v for k, v in out.items() if v}


def generate(alpha: Fraction, beta: Fraction, N: int):
    """Exact series of f to R^N: ({(const, m): object array over exponents 0..N}, val_ok)."""
    a, b = Fraction(alpha), Fraction(beta)
    kap, lam = (a + b) / 2, (a - b) / 2
    taumax = N // 2 if lam != 0 else 2
    Pc, Vc = legendre_PV(max(taumax, 2))
    nmax = 2 * taumax + 5
    emax = N + 2
    TB = Tables(nmax, emax)
    E0 = TB.emin
    RLO, RHI = -5, N - 5
    WR = RHI - RLO + 1
    acc = {k: _vec(WR) for k in ("1", "g", "l2", "lp")}
    val_ok = True
    kscale = _vec(TB.W)
    for i in range(TB.W):
        kscale[i] = kap ** (E0 + i)
    for tau in range(0, taumax + 1):
        bser = {}
        for k in (0, 2):
            v = _vec(N + 1)
            for m in range(N + 1):
                if lam == 0 and m > 0:
                    break
                mo = leg_moment(k + m, tau)
                if mo:
                    v[m] = Fraction(-lam) ** m / factorial(m) * mo
            bser[k] = v
        w = {(2, 2): np.convolve(bser[0], bser[0])[:N + 1],
             (2, 0): -2 * np.convolve(bser[0], bser[2])[:N + 1],
             (0, 0): np.convolve(bser[2], bser[2])[:N + 1]}
        if all(not np.any(w[kk] != 0) for kk in KERNELS):
            continue
        for kk in KERNELS:
            lo = max(0, 2 * tau - 4)
            if np.any(w[kk][:lo] != 0):
                val_ok = False
        Stau = {k: _vec(WR) for k in acc}
        for (j1, j2) in KERNELS:
            pcs = g_pieces(tau, j1, j2, Pc, Vc)
            Xh: Dict[str, np.ndarray] = {}
            for basis, poly in pcs.items():
                n = len(poly)
                assert n <= nmax + 1, (n, nmax)
                for comp, tab in TB.T[basis].items():
                    contrib = poly.dot(tab[:n])
                    Xh[comp] = contrib if comp not in Xh else Xh[comp] + contrib
            vmin = -(j1 + j2 + 1)
            for comp, v in Xh.items():
                if np.any(v[: vmin - E0] != 0):
                    val_ok = False
                seg = (v * kscale)[vmin - E0:]
                prod = np.convolve(w[(j1, j2)], seg)
                off = vmin - RLO
                nn = min(len(prod), WR - off)
                Stau[comp][off:off + nn] = Stau[comp][off:off + nn] + prod[:nn]
        cut = 2 * tau - 5 - RLO
        for comp, v in Stau.items():
            if cut > 0 and np.any(v[:cut] != 0):
                val_ok = False
        for comp in acc:
            acc[comp] = acc[comp] + Stau[comp] * (2 * tau + 1)
    em = _vec(N + 1)
    for n in range(N + 1):
        em[n] = Fraction(-2 * kap) ** n / factorial(n)
    pref = a ** 3 * b ** 3 / 8
    out: Dict[Tuple[str, int], np.ndarray] = {}

    def put(key, v):
        out[key] = out[key] + v if key in out else v
    pl = _prime_logs(kap)
    for comp, v in acc.items():
        full = np.convolve(v, em)[:WR] * pref
        if comp == "lp":                                # ln p = ln kappa + ln R
            put(("1", 1), full)
            for pr, e in pl.items():
                put((f"l{pr}", 0), full * e)
        else:
            put((comp, 0), full)
    return out, val_ok


N_FAST = 40


@lru_cache(maxsize=None)
def series(key: str, N: int = N_FAST):
    a, b, _ = PAIRS[key]
    ser, ok = generate(a, b, N)
    comps = {k: [Fraction(x) for x in v] for k, v in ser.items() if any(x != 0 for x in v)}
    return comps, ok


# ---------------------------------------------------------------- operator on log-Taylor series


def _ff_poly(i: int) -> List[int]:
    """Ascending integer coefficients of k(k-1)...(k-i+1)."""
    c = [1]
    for s in range(i):
        new = [0] * (len(c) + 1)
        for e, v in enumerate(c):
            new[e + 1] += v
            new[e] -= s * v
        c = new
    return c


@lru_cache(maxsize=None)
def _ffd_coeffs(i: int, t: int) -> Tuple[int, ...]:
    c = _ff_poly(i)
    for _ in range(t):
        c = [e * c[e] for e in range(1, len(c))] or [0]
    return tuple(c)


def ffd(i: int, t: int, k: int) -> int:
    """t-th k-derivative of the falling factorial k^(i), at integer k."""
    v = 0
    for c in reversed(_ffd_coeffs(i, t)):
        v = v * k + c
    return v


def apply_operator(op, comps, tmax: int) -> Dict[Tuple[str, int, int], Fraction]:
    """Coefficients of L f at R^t (ln R)^m', t = -r..tmax, per constant: {(const, m', t): value}.

    R^j D^i [R^k ln^m R] = sum_t C(m, t) (d/dk)^t k^(i) R^{k-i+j} ln^{m-t} R.
    Exact for t <= N - r when the series is exact to R^N.
    """
    out: Dict[Tuple[str, int, int], Fraction] = {}
    for (c, m), vals in comps.items():
        for i, p in enumerate(op):
            for j, cij in enumerate(p):
                if not cij:
                    continue
                for k, fk in enumerate(vals):
                    if not fk:
                        continue
                    t_ord = k - i + j
                    if t_ord > tmax:
                        break
                    for t in range(m + 1):
                        w = comb(m, t) * ffd(i, t, k)
                        if w:
                            key = (c, m - t, t_ord)
                            out[key] = out.get(key, Fraction(0)) + w * cij * fk
    return out


def frobenius(op, data: Dict[Tuple[str, int], Dict[int, Fraction]], kmax: int):
    """Frobenius continuation of the series from its values at the R = 0 exponents.

    Returns (F, roots, root_residuals): F = {(const, m): {k: coeff}}, and the left side of every equation at an
    exponent (where the indicial coefficient vanishes), which must be 0 for a consistent solution.
    """
    cij = {(i, j): v for i, p in enumerate(op) for j, v in enumerate(p) if v}
    delta = min(j - i for (i, j) in cij)
    k = sp.Symbol("k")
    I = sum(cv * sp.ff(k, i) for (i, j), cv in cij.items() if j - i == delta)
    roots = sorted(int(rt) for rt in sp.roots(sp.Poly(I, k)) if rt.is_integer)
    consts = sorted({c for c, m in data})
    F, resid = {}, {}
    for c in consts:
        mmax = max(m for cc, m in data if cc == c)
        Fc: Dict[int, Dict[int, Fraction]] = {m: {} for m in range(mmax + 1)}
        for kk in range(min(roots + [0]), kmax + 1):
            for m in range(mmax, -1, -1):
                lhs, top = Fraction(0), Fraction(0)
                for (i, j), cv in cij.items():
                    k2 = kk + delta + i - j
                    for t in range(mmax - m + 1):
                        w = comb(m + t, t) * ffd(i, t, k2)
                        if not w:
                            continue
                        if k2 == kk and t == 0:
                            top += w * cv
                        else:
                            lhs += w * cv * Fc[m + t].get(k2, Fraction(0))
                if top == 0:
                    Fc[m][kk] = data.get((c, m), {}).get(kk, Fraction(0))
                    resid[(c, m, kk)] = lhs
                else:
                    Fc[m][kk] = -lhs / top
        for m in Fc:
            F[(c, m)] = Fc[m]
    return F, roots, resid


def indicial(op) -> Tuple[bool, sp.Poly]:
    """(Fuchsian at R = 0, indicial polynomial in theta)."""
    th = sp.Symbol("theta")
    r = len(op) - 1
    ords = [next(j for j, v in enumerate(p) if v) for p in op]
    vals = [o - i for i, o in enumerate(ords)]
    v = min(vals)
    ind = sum(sp.Rational(op[i][ords[i]].numerator, op[i][ords[i]].denominator) * sp.ff(th, i)
              for i in range(r + 1) if vals[i] == v)
    return vals[r] == v, sp.Poly(sp.expand(ind), th)


def infinity_rank_and_rates(op) -> Tuple[sp.Rational, Dict[sp.Expr, int], int]:
    """Newton-polygon rank at infinity and the characteristic roots of the top-degree edge (with multiplicity)."""
    r = len(op) - 1
    degs = [len(p) - 1 for p in op]
    w = [degs[j] - j for j in range(r + 1)]
    rank = max(sp.Rational(w[j] - w[r], r - j) for j in range(r))
    D = max(degs)
    lam = sp.Symbol("lam")
    chi = sum(sp.Rational(op[k][D].numerator, op[k][D].denominator) * lam ** k
              for k in range(r + 1) if degs[k] == D)
    P = sp.Poly(chi, lam)
    return rank, dict(sp.roots(P)), P.degree()


# ---------------------------------------------------------------- apparent points, mod p


def _primes_from(lo: int):
    n = lo
    while True:
        if sp.isprime(n):
            yield n
        n += 1


def _roots_mod_p(coeffs: List[int], p: int) -> List[int]:
    """All roots in F_p of the integer polynomial (ascending coefficients), by vectorised Horner."""
    x = np.arange(p, dtype=np.int64)
    v = np.zeros(p, dtype=np.int64)
    for c in reversed(coeffs):
        v = (v * x + (c % p)) % p
    return [int(t) for t in np.nonzero(v == 0)[0]]


def _rank_mod_p(rows: List[List[int]], p: int) -> int:
    M = [[x % p for x in row] for row in rows]
    rank, ncol = 0, len(M[0])
    for col in range(ncol):
        piv = next((i for i in range(rank, len(M)) if M[i][col]), None)
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = pow(M[rank][col], p - 2, p)
        M[rank] = [x * inv % p for x in M[rank]]
        for i in range(len(M)):
            if i != rank and M[i][col]:
                f = M[i][col]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rank])]
        rank += 1
    return rank


def apparent_at_root_mod_p(op_int: List[List[int]], rho: int, p: int, emax: int = 400):
    """At t = R - rho (rho a simple root of p_r mod p): exponents {0..r-2, e} and the analytic-solution count.

    Returns (e or None, dim).  Apparent  <=>  e is an integer >= r-1 (all exponents distinct non-negative
    integers) and dim == r (no log solution).  Over F_p, "e is an integer" is read as e mod p in [r-1, emax].
    """
    r = len(op_int) - 1
    tay = []                                            # tay[i][l] = coefficient of t^l in p_i(rho + t)
    for pc in op_int:
        d = len(pc) - 1
        tay.append([sum(pc[j] * comb(j, l) * pow(rho, j - l, p) for j in range(l, d + 1)) % p
                    for l in range(d + 1)])
    assert tay[r][0] == 0 and tay[r][1] != 0, "rho not a simple root of p_r mod p"
    e = (r - 1 - tay[r - 1][0] * pow(tay[r][1], p - 2, p)) % p
    if not (r - 1 <= e <= emax):
        return None, None
    rows = []
    for s in range(1 - r, e + 2 - r):                   # coefficient of t^s; unknowns a_0..a_e
        row = []
        for k in range(e + 1):
            acc = 0
            for i in range(r + 1):
                l = s - k + i
                if 0 <= l < len(tay[i]) and k >= i:
                    acc += tay[i][l] * ffd(i, 0, k)
            row.append(acc % p)
        rows.append(row)
    return e, (e + 1) - _rank_mod_p(rows, p)


@lru_cache(maxsize=None)
def apparent_readings(key: str, n_primes: int = 2, lo: int = 100003):
    """The record's reading: on primes where the degree-11 factor of p_9 has a root, the root is apparent."""
    op = op_coeffs(key)
    op_int = [[int(c) for c in p] for p in op]          # integer polynomials (content removed by the driver)
    assert all(c.denominator == 1 for p in op for c in p)
    p9 = sp.Poly(sum(int(c) * R ** j for j, c in enumerate(op_int[-1])), R)
    _, facs = sp.factor_list(p9.as_expr())
    big = [f for f, m in facs if sp.Poly(f, R).degree() > 1]
    assert len(big) == 1, facs
    fac = sp.Poly(big[0], R)
    fac_c = [int(c) for c in reversed(fac.all_coeffs())]
    out = []
    for p in _primes_from(lo):
        if fac_c[-1] % p == 0:
            continue
        for rho in _roots_mod_p(fac_c, p):
            if rho == 0:
                continue
            e, dim = apparent_at_root_mod_p([[c % p for c in pc] for pc in op_int], rho, p)
            out.append((p, rho, e, dim))
        if len({q for q, *_ in out}) >= n_primes:
            break
    return fac.degree(), [m for f, m in facs if f == R], out


# ================================================================ tests


KEYS = list(PAIRS)


@pytest.mark.parametrize("key", KEYS)
def test_operator_shape(key):
    """Rejects: an operator of the wrong order or degree (the equal-exponent OP_EXCH_1 is order 4, degree 8)."""
    op = op_coeffs(key)
    assert len(op) - 1 == 9
    assert max(len(p) - 1 for p in op) == 16
    assert all(any(p) for p in op)


@pytest.mark.parametrize("key", KEYS)
def test_operator_annihilates_exact_series(key):
    """(a) L f = 0 exactly over Q, every (constant, ln-power) component, orders -9 .. N-9 (N = 40).

    Rejects: any perturbed operator coefficient (every c_ij enters the order t = j at k = i, below N - 9), the other
    pair's or the equal-exponent operator, and a generator with a, b or the Neumann weights wrong.  Non-vacuous: the
    series has the five components the record names and the residual system has >= 4 (N - 9) checked entries.
    """
    comps, ok = series(key)
    assert ok, "per-tau valuation check failed"
    a, b, _ = PAIRS[key]
    kap = (a + b) / 2
    primes = sorted(f"l{q}" for q in _prime_logs(kap) if q != 2)
    assert set(comps) == {("1", 0), ("1", 1), ("g", 0), ("l2", 0)} | {(c, 0) for c in primes}
    res = apply_operator(op_coeffs(key), comps, N_FAST - 9)
    assert len(res) >= 4 * (N_FAST - 9)
    bad = {k: v for k, v in res.items() if v != 0}
    assert not bad, sorted(bad)[:5]


def test_operators_are_pair_specific():
    """Negative control: each pair's operator FAILS on the other pair's series, and so does the equal-exponent
    OP_EXCH_1 (order 4) of tests/test_paper58_r_operators.py.  Rejects an annihilation check that passes anything."""
    s32, _ = series("1_3half")
    s2, _ = series("1_2")
    assert any(v != 0 for v in apply_operator(op_coeffs("1_2"), s32, 20).values())
    assert any(v != 0 for v in apply_operator(op_coeffs("1_3half"), s2, 20).values())
    op1 = [tuple(Fraction(int(c)) for c in reversed(sp.Poly(sp.sympify(t, locals={"R": R}), R).all_coeffs()))
           for t in OP_EXCH_1_EQUAL]
    assert len(op1) == 5
    assert any(v != 0 for v in apply_operator(op1, s32, 20).values())


@pytest.mark.parametrize("key", KEYS)
def test_R0_fuchsian_exponents(key):
    """(b) R = 0 regular singular; exponents {-1 (x2), 0, 1, 2, 3, 4, 6, 9} with multiplicity.

    Rejects: 5 in place of 6, a simple -1, any other exponent set, an irregular (non-Fuchsian) point.
    """
    fuchs, ind = indicial(op_coeffs(key))
    assert fuchs
    assert ind.degree() == 9
    assert {int(k): int(m) for k, m in sp.roots(ind).items()} == EXPONENTS_0


@pytest.mark.parametrize("key", KEYS)
def test_R0_log_solution_and_data_consistency(key):
    """The local solution space at R = 0 is NOT log-free, and f's own data are consistent with it.

    Frobenius from f's values at the exponents reproduces the generator's series to N = 40 (orders between
    exponents are forced), and every equation at an exponent is satisfied.  Rejects an operator whose Frobenius
    recursion disagrees with the series between the data orders, and one for which f's data are inconsistent.
    """
    comps, _ = series(key)
    op = op_coeffs(key)
    data = {kk: {k: v for k, v in enumerate(vals) if v} for kk, vals in comps.items()}
    roots_data = {kk: {k: v for k, v in d.items() if k in EXPONENTS_0} for kk, d in data.items()}
    F, roots, resid = frobenius(op, roots_data, N_FAST)
    assert roots == sorted(EXPONENTS_0)
    assert all(v == 0 for v in resid.values()), {k: v for k, v in resid.items() if v}
    for kk, d in data.items():
        for k in range(N_FAST + 1):
            assert F[kk].get(k, 0) == d.get(k, 0), (kk, k)


@pytest.mark.parametrize("key", KEYS)
def test_lnR_onset_and_coefficient_from_series(key):
    """(d) f's ln R part starts at R^9 ln R with coefficient (8/3375) a^5 b^5 = 9/500 at (1, 3/2), 256/3375 at (1, 2).

    Rejects: 4/3375 (or any other prefactor), an onset other than R^9, the equal-exponent value (8/3375) zeta^10 read
    at the wrong zeta.  The literal values are checked as well as the formula.
    """
    a, b, _ = PAIRS[key]
    comps, _ = series(key)
    ln = comps[("1", 1)]
    k0 = next(k for k, v in enumerate(ln) if v)
    assert k0 == 9
    assert ln[9] == LN_COEF * a ** 5 * b ** 5
    assert ln[9] == {"1_3half": Fraction(9, 500), "1_2": Fraction(256, 3375)}[key]


@pytest.mark.parametrize("key", KEYS)
def test_lnR_coefficient_forced_by_operator(key):
    """(d') The R^9 ln R coefficient is FORCED by the operator from f's non-log data below R^9.

    At the exponent 9 the ln-free equation is a consistency condition that is affine in the R^9 ln R coefficient x;
    solving it gives x = (8/3375) a^5 b^5.  Independent of the generator's own ln R output.  Rejects an operator that
    would carry a different log coefficient, and a wrong expected constant.
    """
    a, b, _ = PAIRS[key]
    comps, _ = series(key)
    base = {kk: {k: v for k, v in enumerate(vals) if v and k in EXPONENTS_0 and k < 9}
            for kk, vals in comps.items() if kk[1] == 0}
    vals = []
    for x in (Fraction(0), Fraction(1)):
        d = dict(base)
        d[("1", 1)] = {9: x}
        _, _, resid = frobenius(op_coeffs(key), d, 9)
        vals.append(resid[("1", 0, 9)])
    assert vals[1] != vals[0]
    x = -vals[0] / (vals[1] - vals[0])
    assert x == LN_COEF * a ** 5 * b ** 5


@pytest.mark.parametrize("key", KEYS)
def test_infinity_rank_one_rates(key):
    """(c) Irregular rank exactly 1 at infinity; the top-degree edge carries all nine behaviours with rates
    {-(a+b), 2a, 2b, 0 (x2), -2a (x2), -2b (x2)} -- six distinct, no a - b rate.

    Rejects: an inserted a - b rate, unsplit rates (the equal-exponent {0, 2z, -2z x2}), a rate with the wrong
    multiplicity, fewer than nine behaviours on the edge, a rank other than 1.
    """
    a, b, _ = PAIRS[key]
    rank, roots, deg = infinity_rank_and_rates(op_coeffs(key))
    assert rank == 1
    assert deg == 9
    assert roots == rates_inf(a, b)
    assert len(roots) == 6
    assert sp.Rational((a - b).numerator, (a - b).denominator) not in roots


@pytest.mark.parametrize("key", KEYS)
def test_finite_points_apparent_mod_p(key):
    """(e) p_9 = c R^5 F(R), F of degree 11 irreducible over Q; at every root of F in F_p (two primes with roots,
    from 100003 up), exponents {0..7, e} with e a non-negative integer >= 8 (the same e on every prime and root:
    Galois-conjugate points share it), and nine analytic solutions (no log solution).

    Rejects: a genuine singular point (non-integer or colliding exponent, a log solution), which any perturbation of
    p_9 or p_8 produces at a root of F mod p with overwhelming probability.  Mod-p reading, as in the record.
    """
    degF, rpow, out = apparent_readings(key)
    assert degF == 11
    assert rpow == [5]
    assert len({p for p, *_ in out}) >= 2
    es = {e for *_, e, _d in out}
    assert None not in es, out
    assert len(es) == 1, out
    assert all(dim == 9 for *_, dim in out), out


# ---------------------------------------------------------------- numeric

FROB_K = 90
# exchange_hp at dps 80 (record numeric step, `exchange_hp_dps80`, tau_max = 32), for the cheap record check below.
RECORD_DPS80 = {
    ("1_3half", "1.4"): "0.272600564598402583271880658867844616019656149590167032857415",
    ("1_3half", "3.015"): "0.0283454185592348529587015928859179181830153685961323648344686",
    ("1_2", "1.4"): "0.20632030202807079004021890670600820962904556376580789324225",
    ("1_2", "3.015"): "0.0152431491774871955879056964903847142468863512279614803377193",
}


def frobenius_value(key: str, Rv: str, K: int = FROB_K, dps: int = 60):
    """Value at R = Rv of the Frobenius solution of the frozen operator, data = f's values at the exponents only."""
    from mpmath import mp
    comps, _ = series(key)
    data = {kk: {k: v for k, v in enumerate(vals) if v and k in EXPONENTS_0} for kk, vals in comps.items()}
    F, _, resid = frobenius(op_coeffs(key), data, K)
    assert all(v == 0 for v in resid.values())
    with mp.workdps(dps):
        x = mp.mpf(Rv)

        def cval(c):
            return {"1": mp.mpf(1), "g": +mp.euler}.get(c) or mp.log(int(c[1:]))
        tot = mp.mpf(0)
        for (c, m), d in F.items():
            tot += cval(c) * mp.log(x) ** m * mp.fsum(mp.mpf(v.numerator) / v.denominator * x ** k
                                                    for k, v in d.items() if v)
        last = max(abs(mp.mpf(v.numerator) / v.denominator * x ** k)
                   for (c, m), d in F.items() for k, v in d.items() if k >= K - 3 and v)
        return +tot, +last


def _exchange_hp(key: str, Rv: str, tau_max: int, dps: int):
    from geovac import qfd_core as Q
    a, b, _ = PAIRS[key]
    S1 = (1, 0, 0)
    return Q.exchange_hp(a, S1, S1, b, S1, S1, sp.Rational(Rv), tau_max=tau_max, dps=dps)


def test_frobenius_matches_production_exchange_hp_fast():
    """(f) Frobenius solution of OP_EXCH_1_3HALF from 8 data orders vs the live production route
    geovac.qfd_core.exchange_hp (independent: per-tau closed forms of the xi and eta integrals, mpmath sum) at
    R = 1.4, dps 30, tau <= 10 (the omitted tail is below 1e-28 there), to 1e-20.

    Rejects: an operator or data whose Frobenius continuation is not the exchange integral -- including a correct
    operator for the WRONG pair, which moves the value at the 1e-1 level.
    """
    from mpmath import mp
    vf, last = frobenius_value("1_3half", "1.4", dps=40)
    assert last < mp.mpf(10) ** -30
    vh, per = _exchange_hp("1_3half", "1.4", tau_max=10, dps=30)
    assert abs(per[-1]) < mp.mpf(10) ** -24
    assert abs(vf - vh) < mp.mpf(10) ** -20, (vf, vh)


@pytest.mark.parametrize("key,Rv", list(RECORD_DPS80), ids=[f"{k}-{r}" for k, r in RECORD_DPS80])
def test_frobenius_matches_record_values(key, Rv):
    """Frobenius at both pairs and both record radii vs the record's 80-digit exchange_hp values, to 1e-40
    (the record measured 1e-82 or better).  Cheap; the live production check is the test above and the slow tier."""
    from mpmath import mp
    vf, last = frobenius_value(key, Rv, K=160, dps=70)
    assert last < mp.mpf(10) ** -50
    with mp.workdps(70):
        assert abs(vf - mp.mpf(RECORD_DPS80[(key, Rv)])) < mp.mpf(10) ** -40


# ---------------------------------------------------------------- slow tier


@pytest.mark.slow
@pytest.mark.parametrize("key", KEYS)
def test_operator_annihilates_exact_series_N120(key):
    """The heldout step's exact check: L f = 0 over Q on orders -9..111 (N = 120), every one unused by the fit
    (the run fits mod p; this is exact over Q).  ~4-5 min per pair."""
    a, b, _ = PAIRS[key]
    ser, ok = generate(a, b, 120)
    assert ok
    comps = {k: [Fraction(x) for x in v] for k, v in ser.items() if any(x != 0 for x in v)}
    res = apply_operator(op_coeffs(key), comps, 111)
    assert len(res) >= 4 * 111
    assert all(v == 0 for v in res.values())


@pytest.mark.slow
@pytest.mark.parametrize("key,Rv", [("1_2", "1.4"), ("1_3half", "3.015"), ("1_2", "3.015")])
def test_frobenius_matches_production_exchange_hp_slow(key, Rv):
    """(f) at the other pair / radius, live: exchange_hp at dps 40 with tau_max = 20 vs Frobenius, to 1e-25."""
    from mpmath import mp
    vf, last = frobenius_value(key, Rv, K=160, dps=60)
    assert last < mp.mpf(10) ** -40
    vh, per = _exchange_hp(key, Rv, tau_max=20, dps=40)
    assert abs(per[-1]) < mp.mpf(10) ** -28
    assert abs(vf - vh) < mp.mpf(10) ** -25, (vf, vh)


# TIMINGS (2026-11-16, this machine): fast tier 22 tests in 46 s (live exchange_hp at R = 1.4, tau <= 10, dps 30:
# 25 s; exact series to N = 40: 6-9 s per pair; everything else < 1.2 s).  Slow tier 5 tests in 1110 s (--slow):
# exact annihilation to N = 120, ~240 s per pair; live exchange_hp tau <= 20, dps 40: 180-220 s per point.
# Measured here: apparent points mod p (p >= 100003) all have exponents {0..7, 9} and nine analytic solutions; the
# rates at (1, 3/2) are {-5/2, 3, 2, 0 x2, -2 x2, -3 x2}.
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper58_exchange_unequal.py --plant "OLD=>NEW" -k <sel>),
# 2026-11-16, written and fired as a separate step after the guard; 28 runs, 26 fired; the 2 that did not are
# explained below.  ANN = "annihilates_exact_series and not N120".
#     p_0 coefficient +1                  "72*(200475000000*R**14" -> "...001*R**14"          -k "ANN and 1_3half"        FIRED
#     p_9 coefficient +1                  "151000312500*R**10" -> "...501*R**10" (in p_9)      -k "apparent and 1_3half"   FIRED
#     p_8 coefficient +1                  "2247184687500*R**11" -> "...501*R**11" (in p_8)     -k "apparent and 1_3half"   FIRED
#     p_8 coefficient +1                  same                                                 -k "ANN and 1_3half"        FIRED
#     kappa, lambda swapped in generator  "(a + b) / 2, (a - b) / 2" -> "(a - b) / 2, (a + b) / 2"   -k ANN             FIRED
#     pair (1, 3/2) -> (1, 5/4)           "Fraction(3, 2)" in PAIRS -> "Fraction(5, 4)"        -k "ANN and 1_3half"        FIRED
#     ln R coefficient 8/3375 -> 4/3375   "LN_COEF = Fraction(8, 3375)" -> "(4, 3375)"         -k lnR_onset                FIRED
#     same, operator-forced route         same                                                 -k lnR_coefficient_forced   FIRED
#     exponent 6 -> 5                     "6: 1, 9: 1}" -> "5: 1, 9: 1}"                       -k fuchsian                 FIRED
#     -1 simple instead of double         "{-1: 2," -> "{-1: 1,"                               -k fuchsian                 FIRED
#     rate a - b inserted                 "q(-(a + b)): 1" -> "q(a - b): 1"                    -k infinity                 FIRED
#     rates unsplit (2b -> 2a)            "q(2 * b): 1" -> "q(2 * a): 1"                       -k infinity                 FIRED
#     equal-exponent OP_EXCH_1 for (1,3/2)  OP_EXCH_1_3HALF -> OP_EXCH_1_EQUAL in PAIRS        -k "1_3half and (ANN)"      FIRED
#     same                                same                                                 -k "operator_shape and 1_3half"   FIRED
#     (1, 2) operator for (1, 3/2)        OP_EXCH_1_3HALF -> OP_EXCH_1_2 in PAIRS              -k "ANN and 1_3half"        FIRED
#     prefactor a^3 b^3/8 -> /4           "pref = a ** 3 * b ** 3 / 8" -> "/ 4"                 -k lnR_onset                FIRED
#     same                                same                                                 -k exchange_hp_fast         FIRED
#     Neumann weight 2 tau + 1 -> + 3     "Stau[comp] * (2 * tau + 1)" -> "(2 * tau + 3)"      -k ANN                      FIRED
#     e^{-2 kappa R} -> e^{-2 a R}        "Fraction(-2 * kap) ** n" -> "Fraction(-2 * a) ** n" -k ANN                      FIRED
#     ln-R derivative of k^(i) dropped    "[e * c[e] for e ..." -> "[c[e] for e ..."           -k ANN                      FIRED
#     exchange_hp at (a, a)               "Q.exchange_hp(a, S1, S1, b," -> "..., a,"            -k exchange_hp_fast         FIRED
#     exchange_hp tau_max 10 -> 3         "tau_max=10, dps=30)" -> "tau_max=3, dps=30)"        -k exchange_hp_fast         FIRED
#     Frobenius K 90 -> 20                "FROB_K = 90" -> "FROB_K = 20"                       -k exchange_hp_fast         FIRED
#     record value off at 1e-34           "...8678446" -> "...8678447" (1_3half, R = 1.4)      -k record_values            FIRED
#     apparent exponent sign              "e = (r - 1 - tay[r - 1][0]" -> "+ tay"              -k apparent                 FIRED
#     p_0 coefficient +1, slow tier       as the first plant                                   -k "N120 and 1_3half"       FIRED
#   Did not fire, and why (not a guard defect):
#     generate(a, b) -> generate(b, a)    -k "ANN or lnR or record_values or exchange_hp_fast": f is symmetric under
#       a <-> b (swapping the exponents swaps the centres; q -> -q and beta_tau(-q) = (-1)^tau beta_tau(q) enters
#       squared), so the swap is not a wrong answer.  The swap that IS wrong (kappa <-> lambda) fires above.
#     prefactor /8 -> /4 with -k ANN alone: L is linear, so L(2f) = 0; annihilation is blind to an overall scale by
#       construction.  The ln R coefficient test and the live numeric check fire on the same plant.
