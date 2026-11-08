"""Round-1 diagnostic dossier driver: S_min structure, decomposition, reduction.

Six-question diagnostic (no paper edits, no production changes, no bulk PSLQ
sweeps).  Pipeline:

  A. Recompute S_min anchor at 120 dps (mpmath.nsum Levin, per v2.23.1 fix).
  B. Exact (Fraction) finite check of the min-weight interchange identity:
       sum_{k=1}^N (sum_{n=k}^N phi(n))^2
         == sum_{n1,n2<=N} min(n1,n2) phi(n1) phi(n2),
     phi(n) = g_n/|lambda_n|^4 = 8(o^2-1)/o^4 at o = 2n+3.
  C. Decomposition of S_min into odd-integer double sums (Hoffman multiple
     t-values) + boundary + diagonal depth-1 terms; 100-dps verification.
  D. Stuffle reduction of the even-weight t-values (exact identities).
  E. Odd-weight t-value closed forms via parity-theorem ansatz; coefficients
     pinned by minimal (<=7-element, weight-homogeneous) PSLQ as a DERIVATION
     AID, each verified to >=90 dps.  (Not an irreducibility sweep.)
  F. Symbolic assembly of the full closed form for S_min; final >=100-dps
     verification against the anchor.
  G. CG channel-count audit: I(n1,n2) = sum_q W(n1,n2,q) vs 2*min(n1,n2)
     (how exact is the min-weight model of the true CG-weighted sunset?).

Output: debug/data/smin_dossier_round1.json
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

import mpmath
import sympy as sp

mpmath.mp.dps = 130
OUT = {}

FROZEN_80 = ("2.4799369380342225544135795008293821446879257866172884583787"
             "98726559552777818374912328")


# =====================================================================
# A. Anchor
# =====================================================================
def T_tail(k):
    a = mpmath.mpf(k) + mpmath.mpf(3) / 2
    return 2 * mpmath.hurwitz(2, a) - mpmath.mpf(1) / 2 * mpmath.hurwitz(4, a)


t0 = time.time()
S_anchor = mpmath.nsum(lambda k: T_tail(int(k)) ** 2, [1, mpmath.inf],
                       method='levin')
dt = time.time() - t0
match_frozen = abs(S_anchor - mpmath.mpf(FROZEN_80)) < mpmath.mpf(10) ** -78
print("A. S_min anchor (120 dps, Levin, %.1fs):" % dt)
print("   %s" % mpmath.nstr(S_anchor, 100))
print("   matches frozen v2.23.1 value to >=78 digits: %s" % match_frozen)
OUT["anchor"] = {"S_min_120dps": mpmath.nstr(S_anchor, 120),
                 "matches_frozen_80": bool(match_frozen),
                 "time_sec": dt}
assert match_frozen


# =====================================================================
# B. Exact finite interchange check (pure Fraction arithmetic)
# =====================================================================
def phi_frac(n: int) -> Fraction:
    """phi(n) = g_n/|lambda_n|^4 = 2(n+1)(n+2)/(n+3/2)^4 = 8(o^2-1)/o^4."""
    o = 2 * n + 3
    return Fraction(8 * (o * o - 1), o ** 4)


N = 40
lhs = sum((sum(phi_frac(n) for n in range(k, N + 1))) ** 2
          for k in range(1, N + 1))
rhs = sum(min(n1, n2) * phi_frac(n1) * phi_frac(n2)
          for n1 in range(1, N + 1) for n2 in range(1, N + 1))
print("\nB. Exact finite interchange identity at N=%d (Fraction): %s"
      % (N, lhs == rhs))
OUT["interchange_exact_N40"] = (lhs == rhs)
assert lhs == rhs

# odd-integer form: same sum in o = 2n+3 variables, weight (min(o1,o2)-3)/2
rhs_odd = sum(Fraction(min(o1, o2) - 3, 2)
              * Fraction(64 * (o1 * o1 - 1) * (o2 * o2 - 1),
                         o1 ** 4 * o2 ** 4)
              for o1 in range(5, 2 * N + 4, 2)
              for o2 in range(5, 2 * N + 4, 2))
print("   odd-integer reindexing exact: %s" % (rhs == rhs_odd))
OUT["odd_reindex_exact_N40"] = (rhs == rhs_odd)
assert rhs == rhs_odd


# =====================================================================
# C. Decomposition into multiple t-values + boundary + diagonal
# =====================================================================
# Off-diagonal (o1 < o2): 64*(o1-3)(o1^2-1)(o2^2-1)/(o1^4 o2^4)
# (o-3)(o^2-1)/o^4 = o^-1 - 3 o^-2 - o^-3 + 3 o^-4   (sympy check below)
# (o^2-1)/o^4      = o^-2 - o^-4
o = sp.symbols('o', positive=True)
p1 = sp.expand((o - 3) * (o ** 2 - 1))
c_coeffs = {c: int(p1.coeff(o, 3 - (c - 1) - 0))
            for c in range(1, 5)}  # exponent o^(4-c) in numerator -> o^-c
# robust extraction:
c_coeffs = {}
for c in range(1, 5):
    c_coeffs[c] = int(p1.coeff(o, 4 - c))
b_coeffs = {2: 1, 4: -1}
assert c_coeffs == {1: 1, 2: -3, 3: -1, 4: 3}, c_coeffs
print("\nC. inner-index coefficients (c -> eps_c):", c_coeffs,
      " outer (b -> eps_b):", b_coeffs)

# diagonal: 32*(o-3)(o^2-1)^2/o^8 = 32 * sum_s d_s o^-s
pd = sp.expand((o - 3) * (o ** 2 - 1) ** 2)
d_coeffs = {s: int(pd.coeff(o, 8 - s)) for s in range(3, 9)}
assert d_coeffs == {3: 1, 4: -3, 5: -2, 6: 6, 7: 1, 8: -3}, d_coeffs
print("   diagonal coefficients (s -> d_s):", d_coeffs)


def lam(s):
    """lambda(s) = sum over odd o>=1 of o^-s = (1 - 2^-s) zeta(s)."""
    return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)


def tval(b, c):
    """Hoffman double t-value t(b,c) = sum_{o2 > o1 >= 1, both odd}
    o2^-b o1^-c.  Inner index o1 = 2j+1, outer tail via Hurwitz."""
    def term(j):
        j = int(j)
        return (mpmath.mpf(2 * j + 1) ** -c
                * mpmath.power(2, -b) * mpmath.hurwitz(b, mpmath.mpf(j) + mpmath.mpf(3) / 2))
    return mpmath.nsum(term, [0, mpmath.inf], method='levin')


print("   computing the 8 double t-values t(b,c), b in {2,4}, c in {1..4} ...")
t0 = time.time()
TV = {}
for b in (2, 4):
    for c in (1, 2, 3, 4):
        TV[(b, c)] = tval(b, c)
print("   done in %.1fs" % (time.time() - t0))
for key in sorted(TV):
    print("   t%s = %s" % (key, mpmath.nstr(TV[key], 40)))
OUT["t_values_60dps"] = {str(k): mpmath.nstr(v, 60) for k, v in TV.items()}

# boundary: t~(c,b) [o1>=5] = t(b,c) - (lam(b)-1) - 3^-c (lam(b)-1-3^-b)
def t_restricted(c, b):
    return (TV[(b, c)] - (lam(b) - 1)
            - mpmath.power(3, -c) * (lam(b) - 1 - mpmath.power(3, -b)))


S_offdiag = mpmath.mpf(0)
for c, ec in c_coeffs.items():
    for b, eb in b_coeffs.items():
        S_offdiag += 64 * ec * eb * t_restricted(c, b)

S_diag = 32 * sum(ds * (lam(s) - 1 - mpmath.power(3, -s))
                  for s, ds in d_coeffs.items())

S_rebuild = S_offdiag + S_diag
resC = abs(S_rebuild - S_anchor)
print("   S_min rebuilt from t-values + diagonal: %s"
      % mpmath.nstr(S_rebuild, 60))
print("   |rebuild - anchor| = %s  (target < 1e-100)" % mpmath.nstr(resC, 5))
OUT["decomposition_residual"] = mpmath.nstr(resC, 5)
assert resC < mpmath.mpf(10) ** -100


# =====================================================================
# D. Stuffle reduction of even-weight t-values (exact identities)
# =====================================================================
stuffle_checks = {
    "t(2,2) = (lam(2)^2 - lam(4))/2":
        abs(TV[(2, 2)] - (lam(2) ** 2 - lam(4)) / 2),
    "t(4,4) = (lam(4)^2 - lam(8))/2":
        abs(TV[(4, 4)] - (lam(4) ** 2 - lam(8)) / 2),
    "t(2,4) + t(4,2) = lam(2)lam(4) - lam(6)":
        abs(TV[(2, 4)] + TV[(4, 2)] - (lam(2) * lam(4) - lam(6))),
}
print("\nD. Stuffle identities (residuals):")
for k, v in stuffle_checks.items():
    print("   %s : %s" % (k, mpmath.nstr(v, 5)))
    assert v < mpmath.mpf(10) ** -100
OUT["stuffle_residuals"] = {k: mpmath.nstr(v, 5)
                            for k, v in stuffle_checks.items()}
# crucial coefficient symmetry: eps(2,4) == eps(4,2) so the symmetric stuffle
# combination is exactly what appears
assert c_coeffs[2] * b_coeffs[4] == c_coeffs[4] * b_coeffs[2] == 3


# =====================================================================
# E. Odd-weight t-value closed forms (parity-theorem ansatz + tiny PSLQ
#    as derivation aid; each verified to >= 90 dps)
# =====================================================================
pi = mpmath.pi
ln2 = mpmath.log(2)
z3, z5, z7 = mpmath.zeta(3), mpmath.zeta(5), mpmath.zeta(7)

ansatz = {
    (2, 1): [("zeta3", z3), ("pi2*ln2", pi ** 2 * ln2),
             ("ln2^3", ln2 ** 3)],
    (4, 1): [("zeta5", z5), ("pi2*zeta3", pi ** 2 * z3),
             ("pi4*ln2", pi ** 4 * ln2), ("zeta3*ln2^2", z3 * ln2 ** 2),
             ("pi2*ln2^3", pi ** 2 * ln2 ** 3), ("ln2^5", ln2 ** 5)],
    (2, 3): [("zeta5", z5), ("pi2*zeta3", pi ** 2 * z3),
             ("pi4*ln2", pi ** 4 * ln2), ("zeta3*ln2^2", z3 * ln2 ** 2),
             ("pi2*ln2^3", pi ** 2 * ln2 ** 3), ("ln2^5", ln2 ** 5)],
    (4, 3): [("zeta7", z7), ("pi2*zeta5", pi ** 2 * z5),
             ("pi4*zeta3", pi ** 4 * z3), ("pi6*ln2", pi ** 6 * ln2),
             ("zeta5*ln2^2", z5 * ln2 ** 2), ("zeta3^2*ln2", z3 ** 2 * ln2)],
}

closed = {}
print("\nE. Odd-weight t-value identifications (derivation aid):")
for key, basis in ansatz.items():
    vals = [TV[key]] + [v for _, v in basis]
    rel = mpmath.pslq(vals, tol=mpmath.mpf(10) ** -90, maxcoeff=10 ** 10)
    if rel is None or rel[0] == 0:
        print("   t%s : NO identification in ansatz -> OBSTRUCTION" % (key,))
        closed[key] = None
        continue
    coeffs = {basis[i - 1][0]: Fraction(-rel[i], rel[0])
              for i in range(1, len(rel)) if rel[i] != 0}
    recon = sum(mpmath.mpf(fr.numerator) / fr.denominator * dict(basis)[name]
                for name, fr in coeffs.items())
    res = abs(TV[key] - recon)
    print("   t%s = %s   (residual %s)"
          % (key, " + ".join("%s*%s" % (fr, name)
                             for name, fr in coeffs.items()),
             mpmath.nstr(res, 5)))
    assert res < mpmath.mpf(10) ** -90
    closed[key] = {name: str(fr) for name, fr in coeffs.items()}
OUT["odd_weight_closed_forms"] = {str(k): v for k, v in closed.items()}


# =====================================================================
# F. Symbolic assembly of the closed form (only if all four identified)
# =====================================================================
if all(closed.get(k) is not None for k in ansatz):
    spi = sp.pi
    sln2 = sp.log(2)
    sz3, sz5, sz7 = sp.zeta(3), sp.zeta(5), sp.zeta(7)
    sym_map = {
        "zeta3": sz3, "zeta5": sz5, "zeta7": sz7,
        "pi2*ln2": spi ** 2 * sln2, "pi4*ln2": spi ** 4 * sln2,
        "pi6*ln2": spi ** 6 * sln2, "pi2*zeta3": spi ** 2 * sz3,
        "pi2*zeta5": spi ** 2 * sz5, "pi4*zeta3": spi ** 4 * sz3,
        "ln2^3": sln2 ** 3, "zeta3*ln2^2": sz3 * sln2 ** 2,
        "pi2*ln2^3": spi ** 2 * sln2 ** 3, "ln2^5": sln2 ** 5,
        "zeta5*ln2^2": sz5 * sln2 ** 2, "zeta3^2*ln2": sz3 ** 2 * sln2,
    }

    def slam(s):
        # lambda(s) symbolic; even s in terms of pi
        zeta_even = {2: spi ** 2 / 6, 4: spi ** 4 / 90, 6: spi ** 6 / 945,
                     8: spi ** 8 / 9450}
        zs = zeta_even.get(s, sp.zeta(s) if s % 2 == 0 else
                           {3: sz3, 5: sz5, 7: sz7}[s])
        return (1 - sp.Rational(1, 2 ** s)) * zs

    tsym = {
        (2, 2): (slam(2) ** 2 - slam(4)) / 2,
        (4, 4): (slam(4) ** 2 - slam(8)) / 2,
    }
    for key in ansatz:
        tsym[key] = sum(sp.Rational(fr) * sym_map[name]
                        for name, fr in closed[key].items())
    # the (2,4)+(4,2) pair enters with equal coefficient +3:
    pair_24 = slam(2) * slam(4) - slam(6)

    Ssym = sp.Integer(0)
    for c, ec in c_coeffs.items():
        for b, eb in b_coeffs.items():
            if (b, c) in ((2, 4), (4, 2)):
                continue  # handled as a pair
            boundary = ((slam(b) - 1)
                        + sp.Rational(1, 3 ** c) * (slam(b) - 1
                                                    - sp.Rational(1, 3 ** b)))
            Ssym += 64 * ec * eb * (tsym[(b, c)] - boundary)
    # symmetric pair (coefficient +3 each):
    bnd_24 = ((slam(2) - 1) + sp.Rational(1, 3 ** 4) * (slam(2) - 1 - sp.Rational(1, 9)))
    bnd_42 = ((slam(4) - 1) + sp.Rational(1, 3 ** 2) * (slam(4) - 1 - sp.Rational(1, 81)))
    Ssym += 64 * 3 * (pair_24 - bnd_24 - bnd_42)
    # diagonal:
    Ssym += 32 * sum(ds * (slam(s) - 1 - sp.Rational(1, 3 ** s))
                     for s, ds in d_coeffs.items())

    Ssym = sp.expand(sp.simplify(sp.expand(Ssym)))
    print("\nF. Closed form for S_min:")
    print("   S_min =", Ssym)
    Snum = mpmath.mpf(sp.N(Ssym, 120).__str__())
    resF = abs(Snum - S_anchor)
    print("   closed form evaluated: %s" % mpmath.nstr(Snum, 60))
    print("   |closed - anchor| = %s" % mpmath.nstr(resF, 5))
    OUT["closed_form"] = {"sympy": str(Ssym),
                          "latex": sp.latex(Ssym),
                          "residual_vs_anchor": mpmath.nstr(resF, 5),
                          "verified_100dps": bool(resF < mpmath.mpf(10) ** -100)}
else:
    print("\nF. SKIPPED: at least one odd-weight t-value not identified.")
    OUT["closed_form"] = None


# =====================================================================
# G. CG channel-count audit: I(n1,n2) vs 2*min(n1,n2)
# =====================================================================
def so4_count(n1: int, n2: int, q: int) -> int:
    j1_L = Fraction(n1 + 1, 2); j1_R = Fraction(n1, 2)
    j2_L = Fraction(n2, 2); j2_R = Fraction(n2 + 1, 2)
    cnt = 0
    for jg_L, jg_R in ((Fraction(q + 1, 2), Fraction(q - 1, 2)),
                       (Fraction(q - 1, 2), Fraction(q + 1, 2))):
        if jg_L < 0 or jg_R < 0:
            continue
        if (abs(j1_L - jg_L) <= j2_L <= j1_L + jg_L
                and abs(j1_R - jg_R) <= j2_R <= j1_R + jg_R):
            cnt += 1
    return cnt


def I_channel(n1: int, n2: int) -> int:
    tot = 0
    for q in range(max(1, abs(n1 - n2)), n1 + n2 + 1):
        if (n1 + n2 + q) % 2 == 0:
            continue
        tot += so4_count(n1, n2, q)
    return tot


print("\nG. CG channel count I(n1,n2) vs 2*min(n1,n2)  (delta = I - 2min):")
delta_table = {}
for n1 in range(0, 13):
    row = []
    for n2 in range(0, 13):
        d = I_channel(n1, n2) - 2 * min(n1, n2)
        row.append(d)
        delta_table[(n1, n2)] = d
    print("   n1=%2d: %s" % (n1, row))
nonzero = {k: v for k, v in delta_table.items() if v != 0}
print("   nonzero deltas:", nonzero if len(nonzero) < 40 else
      "%d entries" % len(nonzero))
OUT["cg_audit"] = {"delta_nonzero": {str(k): v for k, v in nonzero.items()}}

# quantitative: true CG-weighted sunset partial vs 2*S_min partial (N=40)
def phi_mp(n):
    a = mpmath.mpf(n) + mpmath.mpf(3) / 2
    return 2 / a ** 2 - mpmath.mpf(1) / 2 / a ** 4


sunset_partial = sum(I_channel(n1, n2) * phi_mp(n1) * phi_mp(n2)
                     for n1 in range(0, 41) for n2 in range(0, 41))
smin_partial = sum(2 * min(n1, n2) * phi_mp(n1) * phi_mp(n2)
                   for n1 in range(0, 41) for n2 in range(0, 41))
print("   partial (N=40): sunset = %s, 2*S_min-model = %s, diff = %s"
      % (mpmath.nstr(sunset_partial, 12), mpmath.nstr(smin_partial, 12),
         mpmath.nstr(sunset_partial - smin_partial, 8)))
OUT["cg_audit"]["sunset_partial_N40"] = mpmath.nstr(sunset_partial, 30)
OUT["cg_audit"]["two_smin_partial_N40"] = mpmath.nstr(smin_partial, 30)
OUT["cg_audit"]["difference_N40"] = mpmath.nstr(
    sunset_partial - smin_partial, 30)

outp = Path(__file__).parent / "data" / "smin_dossier_round1.json"
outp.parent.mkdir(parents=True, exist_ok=True)
outp.write_text(json.dumps(OUT, indent=2))
print("\nSaved:", outp)
