"""Round-1b: full level-2 (alternating-MZV) basis identification of the two
remaining odd-weight double t-values t(4,1) [w5] and t(4,3) [w7], and final
closed-form assembly of S_min, all at 260 dps.

Weight-5 level-2 conjectural basis (dim 8, Deligne/BBV):
    zeta(5), pi^2 zeta(3), pi^4 ln2, zeta(3) ln2^2, pi^2 ln2^3, ln2^5,
    Li4(1/2) ln2, Li5(1/2)
Weight-7 level-2 basis (dim 21 = 17 product monomials + zeta(7) + Li7(1/2)
    + 2 depth-2 irreducibles; we include 2 named depth-2 alternating
    candidates A1 = zeta(-6,1), A2 = zeta(-5,2) in a second pass).

Output: debug/data/smin_dossier_round1b.json
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

import mpmath

mpmath.mp.dps = 260
OUT = {}

pi = mpmath.pi
ln2 = mpmath.log(2)
z3, z5, z7 = mpmath.zeta(3), mpmath.zeta(5), mpmath.zeta(7)
Li4 = mpmath.polylog(4, mpmath.mpf(1) / 2)
Li5 = mpmath.polylog(5, mpmath.mpf(1) / 2)
Li6 = mpmath.polylog(6, mpmath.mpf(1) / 2)
Li7 = mpmath.polylog(7, mpmath.mpf(1) / 2)


def lam(s):
    return (1 - mpmath.power(2, -s)) * mpmath.zeta(s)


def tval(b, c):
    def term(j):
        j = int(j)
        return (mpmath.mpf(2 * j + 1) ** -c
                * mpmath.power(2, -b)
                * mpmath.hurwitz(b, mpmath.mpf(j) + mpmath.mpf(3) / 2))
    return mpmath.nsum(term, [0, mpmath.inf], method='levin')


print("Computing the 8 t-values at 260 dps ...")
t0 = time.time()
TV = {}
for b in (2, 4):
    for c in (1, 2, 3, 4):
        TV[(b, c)] = tval(b, c)
        print("  t(%d,%d) done (%.0fs)" % (b, c, time.time() - t0))
OUT["t_values_200dps"] = {str(k): mpmath.nstr(v, 200) for k, v in TV.items()}

print("Anchor at 260 dps ...")
S_anchor = mpmath.nsum(
    lambda k: (2 * mpmath.hurwitz(2, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)
               - mpmath.mpf(1) / 2
               * mpmath.hurwitz(4, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)) ** 2,
    [1, mpmath.inf], method='levin')
print("  S_min = %s..." % mpmath.nstr(S_anchor, 60))


def identify(name, value, basis, maxcoeff=10 ** 8, tol_exp=-220):
    vals = [value] + [v for _, v in basis]
    rel = mpmath.pslq(vals, tol=mpmath.mpf(10) ** tol_exp, maxcoeff=maxcoeff)
    if rel is None or rel[0] == 0:
        print("  %s: NO relation (basis size %d)" % (name, len(basis)))
        return None
    coeffs = {basis[i - 1][0]: Fraction(-rel[i], rel[0])
              for i in range(1, len(rel)) if rel[i] != 0}
    recon = sum(mpmath.mpf(fr.numerator) / fr.denominator * dict(basis)[nm]
                for nm, fr in coeffs.items())
    res = abs(value - recon)
    print("  %s = %s" % (name, " + ".join(
        "(%s)*%s" % (fr, nm) for nm, fr in coeffs.items())))
    print("    residual %s" % mpmath.nstr(res, 5))
    if res > mpmath.mpf(10) ** (tol_exp + 20):
        print("    REJECTED (residual too large)")
        return None
    return coeffs


# ---------------------------------------------------------------------
# t(4,1): full weight-5 level-2 basis (dim 8)
# ---------------------------------------------------------------------
w5_basis = [
    ("zeta5", z5), ("pi2*zeta3", pi ** 2 * z3), ("pi4*ln2", pi ** 4 * ln2),
    ("zeta3*ln2^2", z3 * ln2 ** 2), ("pi2*ln2^3", pi ** 2 * ln2 ** 3),
    ("ln2^5", ln2 ** 5), ("Li4*ln2", Li4 * ln2), ("Li5", Li5),
]
print("\nt(4,1) against full w5 level-2 basis:")
c41 = identify("t(4,1)", TV[(4, 1)], w5_basis)
OUT["t41"] = {k: str(v) for k, v in c41.items()} if c41 else None

# cross-check t(2,3) against the same full basis (should reproduce round-1a)
print("t(2,3) against full w5 level-2 basis (cross-check):")
c23 = identify("t(2,3)", TV[(2, 3)], w5_basis)
OUT["t23"] = {k: str(v) for k, v in c23.items()} if c23 else None

# ---------------------------------------------------------------------
# t(4,3): weight-7 level-2 product monomials + zeta7 + Li7  (19 elements)
# ---------------------------------------------------------------------
w7_products = [
    ("zeta7", z7), ("pi2*zeta5", pi ** 2 * z5), ("pi4*zeta3", pi ** 4 * z3),
    ("pi6*ln2", pi ** 6 * ln2), ("zeta5*ln2^2", z5 * ln2 ** 2),
    ("zeta3^2*ln2", z3 ** 2 * ln2), ("zeta3*pi2*ln2^2", z3 * pi ** 2 * ln2 ** 2),
    ("zeta3*ln2^4", z3 * ln2 ** 4), ("pi4*ln2^3", pi ** 4 * ln2 ** 3),
    ("pi2*ln2^5", pi ** 2 * ln2 ** 5), ("ln2^7", ln2 ** 7),
    ("Li4*zeta3", Li4 * z3), ("Li4*pi2*ln2", Li4 * pi ** 2 * ln2),
    ("Li4*ln2^3", Li4 * ln2 ** 3), ("Li5*pi2", Li5 * pi ** 2),
    ("Li5*ln2^2", Li5 * ln2 ** 2), ("Li6*ln2", Li6 * ln2), ("Li7", Li7),
]
print("\nt(4,3) against w7 level-2 product basis (19 elements):")
c43 = identify("t(4,3)", TV[(4, 3)], w7_products)
OUT["t43_products_only"] = ({k: str(v) for k, v in c43.items()}
                            if c43 else None)

if c43 is None:
    # second pass: add two named depth-2 alternating irreducible candidates
    print("  adding depth-2 alternating candidates zeta(-6,1), zeta(-5,2):")

    def alt_double(b, c):
        """zeta(bbar, c) = sum_{m>n>=1} (-1)^m / (m^b n^c)."""
        def term(n):
            n = int(n)
            # inner: sum over m > n of (-1)^m/m^b  (alternating Hurwitz tail)
            # (-1)^m/m^b for m>=n+1:  use Lerch transcendent:
            # sum_{m>=n+1} (-1)^m/m^b = (-1)^(n+1) * Phi(-1, b, n+1)
            tail = mpmath.lerchphi(-1, b, n + 1) * (-1) ** (n + 1)
            return tail / mpmath.mpf(n) ** c
        return mpmath.nsum(term, [1, mpmath.inf], method='levin')

    t0 = time.time()
    A1 = alt_double(6, 1)
    A2 = alt_double(5, 2)
    print("  A1 = zeta(-6,1) = %s (%.0fs)" % (mpmath.nstr(A1, 30),
                                              time.time() - t0))
    print("  A2 = zeta(-5,2) = %s" % mpmath.nstr(A2, 30))
    w7_ext = w7_products + [("zeta(-6,1)", A1), ("zeta(-5,2)", A2)]
    c43 = identify("t(4,3)", TV[(4, 3)], w7_ext, maxcoeff=10 ** 8)
    OUT["t43_with_depth2"] = ({k: str(v) for k, v in c43.items()}
                              if c43 else None)
    w7_basis_used = w7_ext
else:
    w7_basis_used = w7_products

# ---------------------------------------------------------------------
# Final assembly (numeric, exact Fraction coefficients)
# ---------------------------------------------------------------------
print("\nFinal assembly:")
c_coeffs = {1: 1, 2: -3, 3: -1, 4: 3}
b_coeffs = {2: 1, 4: -1}
d_coeffs = {3: 1, 4: -3, 5: -2, 6: 6, 7: 1, 8: -3}


def t_restricted(c, b):
    return (TV[(b, c)] - (lam(b) - 1)
            - mpmath.power(3, -c) * (lam(b) - 1 - mpmath.power(3, -b)))


S_offdiag = sum(64 * ec * eb * t_restricted(c, b)
                for c, ec in c_coeffs.items()
                for b, eb in b_coeffs.items())
S_diag = 32 * sum(ds * (lam(s) - 1 - mpmath.power(3, -s))
                  for s, ds in d_coeffs.items())
res_decomp = abs(S_offdiag + S_diag - S_anchor)
print("  decomposition residual at 260 dps: %s" % mpmath.nstr(res_decomp, 5))
OUT["decomposition_residual_260"] = mpmath.nstr(res_decomp, 5)

# closed-form substitution: replace each t-value by its reduction
t_closed = {
    (2, 2): (lam(2) ** 2 - lam(4)) / 2,
    (4, 4): (lam(4) ** 2 - lam(8)) / 2,
}
basis_lookup = dict(w5_basis + w7_products
                    + [("zeta3", z3), ("pi2*ln2", pi ** 2 * ln2),
                       ("ln2^3", ln2 ** 3)])
if 'A1' in dir():
    basis_lookup["zeta(-6,1)"] = A1
    basis_lookup["zeta(-5,2)"] = A2

w3_basis = [("zeta3", z3), ("pi2*ln2", pi ** 2 * ln2), ("ln2^3", ln2 ** 3)]
print("t(2,1) re-identification at 260 dps:")
c21 = identify("t(2,1)", TV[(2, 1)], w3_basis)
OUT["t21"] = {k: str(v) for k, v in c21.items()} if c21 else None

all_identified = all(x is not None for x in (c21, c23, c41, c43))
OUT["fully_reduced"] = bool(all_identified)
if all_identified:
    def closed_val(coeffs):
        return sum(mpmath.mpf(fr.numerator) / fr.denominator
                   * basis_lookup[nm] for nm, fr in coeffs.items())

    t_closed[(2, 1)] = closed_val(c21)
    t_closed[(2, 3)] = closed_val(c23)
    t_closed[(4, 1)] = closed_val(c41)
    t_closed[(4, 3)] = closed_val(c43)
    # symmetric pair via stuffle:
    pair24 = lam(2) * lam(4) - lam(6)

    S_closed = mpmath.mpf(0)
    for c, ec in c_coeffs.items():
        for b, eb in b_coeffs.items():
            bnd = ((lam(b) - 1)
                   + mpmath.power(3, -c) * (lam(b) - 1 - mpmath.power(3, -b)))
            if (b, c) in ((2, 4), (4, 2)):
                S_closed += 64 * ec * eb * (-bnd)  # t-part via pair below
            else:
                S_closed += 64 * ec * eb * (t_closed[(b, c)] - bnd)
    S_closed += 64 * 3 * pair24  # +3*t(2,4) +3*t(4,2) combined
    S_closed += S_diag

    res_final = abs(S_closed - S_anchor)
    print("\n  S_min closed form evaluated: %s" % mpmath.nstr(S_closed, 60))
    print("  |closed - anchor| = %s" % mpmath.nstr(res_final, 5))
    OUT["final_residual"] = mpmath.nstr(res_final, 5)
    OUT["final_verified_200dps"] = bool(res_final < mpmath.mpf(10) ** -200)
else:
    print("\n  Closed form NOT assembled (obstruction stands).")

outp = Path(__file__).parent / "data" / "smin_dossier_round1b.json"
outp.write_text(json.dumps(OUT, indent=2))
print("\nSaved:", outp)
