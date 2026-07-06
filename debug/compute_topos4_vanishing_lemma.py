"""Sprint Topos-4 — the general-parameter vanishing lemma for inter-frame
hydrogenic overlaps.

Topos-2 left the overlap-vanishing question as a biconditional whose (=>)
direction was only panel-scoped: rate coincidence Z1/n1 = Z2/n2 with
|n1-n2| >= 2 forces a zero (exhaustively verified), but sporadic
off-coincidence zeros exist beyond the panel (boundary case Z'=2, l=3,
n=n'=5; rates 1/5 vs 2/5).

This driver characterizes ALL zeros exactly.  Claim:

  For same-l inter-frame overlaps, in the rate variable
      t = (p - q) / (p + q),   p = Z1/n1,  q = Z2/n2,   t in (-1, 1),
  the exact unnormalized overlap factors as
      S(p, q) = (positive monomial prefactor) * P_{n1,n2,l}(t),
  where P is an explicit polynomial of degree D = (n1-l-1) + (n2-l-1).
  Hence for physical rates (p, q > 0):
      S = 0  <=>  t is a real root of P in (-1, 1).
  The rate-coincidence zeros are the single root t = 0 (present <=>
  |n1 - n2| >= 2); every "sporadic" zero is a nonzero root of P.

We further identify P with a classical Jacobi polynomial and verify the
boundary case lands exactly on one of its roots.

Output: debug/data/sprint_topos4_vanishing_lemma.json
Backing test: tests/test_topos4_vanishing_lemma.py (self-contained, exact).
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb, factorial

import sympy as sp


# ----------------------------------------------------------------------
# exact overlap, symbolic in the two decay rates p, q
# (same convention as debug/compute_topos2_family1_meetjoin.py; positive
#  monomial rescalings are support-neutral and do not move physical zeros)
# ----------------------------------------------------------------------

def genlag_coeffs(k: int, alpha: int):
    """c_j of L_k^{alpha}(x) = sum_j c_j x^j, exact rationals."""
    return [sp.Rational((-1) ** j * comb(k + alpha, k - j), factorial(j))
            for j in range(k + 1)]


def overlap_symbolic(n1: int, n2: int, l: int, p, q):
    """Exact symbolic  int_0^inf R_{n1 l}(p) R_{n2 l}(q) r^2 dr,
    R_{nl}(rate) ~ sum_j lag[j] (2*rate)^j r^{l+j} exp(-rate r)."""
    lag1 = genlag_coeffs(n1 - l - 1, 2 * l + 1)
    lag2 = genlag_coeffs(n2 - l - 1, 2 * l + 1)
    total = sp.Integer(0)
    for j1, c1 in enumerate(lag1):
        a1 = c1 * (2 * p) ** j1
        for j2, c2 in enumerate(lag2):
            a2 = c2 * (2 * q) ** j2
            k = (l + j1) + (l + j2) + 2
            total += a1 * a2 * sp.factorial(k) / (p + q) ** (k + 1)
    return sp.simplify(total)


def overlap_exact(Z1: int, n1: int, Z2: int, n2: int, l: int) -> Fraction:
    """Exact numeric overlap for integer Z (matches the Topos-2 driver)."""
    lag1 = genlag_coeffs(n1 - l - 1, 2 * l + 1)
    lag2 = genlag_coeffs(n2 - l - 1, 2 * l + 1)
    p = Fraction(Z1, n1)
    q = Fraction(Z2, n2)
    rate = p + q
    total = Fraction(0)
    for j1, c1 in enumerate(lag1):
        a1 = Fraction(int(c1.p), int(c1.q)) * (2 * p) ** j1
        for j2, c2 in enumerate(lag2):
            a2 = Fraction(int(c2.p), int(c2.q)) * (2 * q) ** j2
            k = (l + j1) + (l + j2) + 2
            total += a1 * a2 * Fraction(factorial(k)) / rate ** (k + 1)
    return total


# ----------------------------------------------------------------------
# factor S into (prefactor) * P(t),  t = (p - q)/(p + q)
# ----------------------------------------------------------------------

def overlap_numerator(n1: int, n2: int, l: int, p, q):
    """Homogeneous-of-degree-D numerator N(p,q) of the overlap, obtained by
    clearing the common denominator (p+q)^{n1+n2+1}.  S = N(p,q)/(p+q)^{K+1}
    with K = n1+n2, so S=0 (for p,q>0) <=> N(p,q)=0.  Pure polynomial
    arithmetic (no simplify): each term is coeff * p^j1 q^j2 (p+q)^{D-j1-j2},
    total degree D = n1+n2-2l-2 for every term."""
    lag1 = genlag_coeffs(n1 - l - 1, 2 * l + 1)
    lag2 = genlag_coeffs(n2 - l - 1, 2 * l + 1)
    D = (n1 - l - 1) + (n2 - l - 1)          # = n1+n2-2l-2
    N = sp.Integer(0)
    for j1, c1 in enumerate(lag1):
        for j2, c2 in enumerate(lag2):
            k = 2 * l + j1 + j2 + 2
            coeff = c1 * c2 * 2 ** (j1 + j2) * sp.factorial(k)
            N += coeff * p ** j1 * q ** j2 * (p + q) ** (D - j1 - j2)
    return sp.expand(N), D


def poly_in_t(n1: int, n2: int, l: int):
    """Return (P(t) as sympy Poly, degree D, symbol t).  Cheap: no roots.

    N(p,q) is homogeneous of degree D, so substituting the linear forms
    p = 1+t, q = 1-t (=> t = (p-q)/(p+q)) gives P(t) = N(1+t, 1-t), a
    polynomial in t whose roots in (-1,1) are EXACTLY the physical
    (p,q>0) vanishing loci.  No rational substitution, no simplify."""
    t = sp.symbols('t', real=True)
    N, D = overlap_numerator(n1, n2, l, 1 + t, 1 - t)
    P = sp.Poly(sp.expand(N), t)
    if P.degree() >= 1:
        P = sp.Poly(sp.primitive(P.as_expr())[1], t)
    return P, D, t


def interior_roots(P, t):
    """Exact roots of P in the open interval (-1, 1) with multiplicity."""
    if P.degree() < 1:
        return {}
    out = {}
    for r, mult in sp.roots(P).items():
        if r.is_real and sp.Abs(r) < 1:
            out[r] = mult
    return out


def jacobi_match(n1: int, n2: int, l: int, P, t):
    """Try to identify P(t) with a classical Jacobi polynomial P_d^{(a,b)}(t)
    up to a rational scalar, where d = deg P.  Returns (a, b, scalar) or
    None."""
    d = P.degree()
    if d < 1:
        return None
    for a in range(0, 2 * l + 6):
        for b in range(0, 2 * l + 6):
            Jp = sp.Poly(sp.expand(sp.jacobi(d, a, b, t)), t)
            if Jp.degree() != d or Jp.LC() == 0:
                continue
            scal = sp.nsimplify(P.LC() / Jp.LC())
            if sp.expand(P.as_expr() - scal * Jp.as_expr()) == 0:
                return (a, b, scal)
    return None


def main():
    out = {}

    # --- 0. reproduce the three exact pins (convention control) -----------
    pins = {
        "S(1s1|1s2)": (overlap_exact(1, 1, 2, 1, 0), Fraction(2, 27)),
        "S(1s1|2s2)": (overlap_exact(1, 1, 2, 2, 0), Fraction(-1, 4)),
        "S(2p1|2p2)": (overlap_exact(1, 2, 2, 2, 1), Fraction(256, 81)),
    }
    out["pins"] = {k: {"computed": str(v), "expected": str(e),
                       "match": v == e} for k, (v, e) in pins.items()}

    # --- 1. the boundary case: n=n'=5, l=3 --------------------------------
    P, deg, t = poly_in_t(5, 5, 3)
    interior = interior_roots(P, t)
    match = jacobi_match(5, 5, 3, P, t)
    # physical check: rates 1/5 (Z=1,n=5) and 2/5 (Z=2,n=5) => t = -1/3
    t_boundary = sp.Rational(sp.Rational(1, 5) - sp.Rational(2, 5),
                             1) / (sp.Rational(1, 5) + sp.Rational(2, 5))
    S_boundary = overlap_exact(1, 5, 2, 5, 3)
    out["boundary_n5n5l3"] = {
        "P(t)": str(P.as_expr()),
        "P(t)_factored": str(sp.factor(P.as_expr())),
        "degree": deg,
        "interior_roots": {str(r): int(m) for r, m in interior.items()},
        "jacobi_match (a,b,scalar)": (None if match is None else
                                      [int(match[0]), int(match[1]),
                                       str(match[2])]),
        "t_at_Z1_Zp2": str(t_boundary),
        "P(t_boundary)": str(sp.simplify(P.as_expr().subs(t, t_boundary))),
        "S_exact_at_boundary": str(S_boundary),
        "vanishes": S_boundary == 0,
    }

    # --- 2. the t=0 (rate-coincidence) criterion across a grid ------------
    #  claim: t=0 is a root of P  <=>  |n1 - n2| >= 2
    coincidence = []
    for n1 in range(1, 7):
        for n2 in range(1, 7):
            for l in range(0, min(n1, n2)):
                P, deg, t = poly_in_t(n1, n2, l)
                t0_root = (P.eval(0) == 0)
                coincidence.append({
                    "n1": n1, "n2": n2, "l": l,
                    "degree": deg,
                    "t0_is_root": bool(t0_root),
                    "abs_dn_ge_2": abs(n1 - n2) >= 2,
                    "consistent": bool(t0_root) == (abs(n1 - n2) >= 2),
                })
    out["rate_coincidence_criterion"] = {
        "rows": coincidence,
        "all_consistent": all(r["consistent"] for r in coincidence),
    }

    # --- 3. Jacobi identification across a family --------------------------
    jac = []
    for (n1, n2, l) in [(1, 1, 0), (2, 1, 0), (3, 1, 0), (2, 2, 0),
                        (3, 2, 0), (3, 3, 0), (2, 2, 1), (3, 3, 1),
                        (5, 5, 3), (4, 2, 1), (5, 3, 2)]:
        P, deg, t = poly_in_t(n1, n2, l)
        interior = interior_roots(P, t)
        m = jacobi_match(n1, n2, l, P, t)
        jac.append({
            "n1": n1, "n2": n2, "l": l, "degree": deg,
            "P(t)": str(P.as_expr()),
            "P(t)_factored": str(sp.factor(P.as_expr())),
            "jacobi (a,b)": (None if m is None else [int(m[0]), int(m[1])]),
            "interior_roots": [str(r) for r in interior],
        })
    out["jacobi_family"] = jac

    # --- 4. exhaustive enumeration of INTEGER-hit sporadic zeros ----------
    #  scan small integer (Z1, Z2, n1, n2, l), record every vanishing overlap
    #  and whether it is rate-coincidence or sporadic.
    #  every integer zero is uniquely one of the three root families:
    #    rate_coincidence  (Z1/n1 == Z2/n2, t = 0),
    #    orthogonality     (Z1 == Z2, n1 != n2, standard <n1 l|n2 l> = 0),
    #    residual          (neither; a nonzero non-orthogonality root of P).
    zeros = {"rate_coincidence": 0, "orthogonality": 0, "residual": []}
    N = 7
    for n1 in range(1, N):
        for n2 in range(n1, N):          # n1 <= n2 wlog
            for l in range(0, min(n1, n2)):
                for Z1 in range(1, N):
                    for Z2 in range(1, N):
                        if overlap_exact(Z1, n1, Z2, n2, l) != 0:
                            continue
                        p, q = Fraction(Z1, n1), Fraction(Z2, n2)
                        if p == q:
                            zeros["rate_coincidence"] += 1
                        elif Z1 == Z2:
                            zeros["orthogonality"] += 1
                        else:
                            zeros["residual"].append({
                                "Z1": Z1, "n1": n1, "Z2": Z2, "n2": n2,
                                "l": l, "t": str((p - q) / (p + q))})
    total = (zeros["rate_coincidence"] + zeros["orthogonality"]
             + len(zeros["residual"]))
    out["integer_zero_census"] = {
        "scan_range": f"Z,n < {N}",
        "total_zeros": total,
        "rate_coincidence": zeros["rate_coincidence"],
        "orthogonality": zeros["orthogonality"],
        "residual_count": len(zeros["residual"]),
        "residual_zeros": zeros["residual"],
        "fully_accounted": True,   # by construction the three are exhaustive
    }

    # --- 5. the three structural root families ----------------------------
    #  (i)  t = 0                       rate coincidence  <=> |dn| >= 2
    #  (ii) t = (n2-n1)/(n2+n1)         same-Z orthogonality  <n1 l|n2 l>=0
    #                                   (always a root for n1 != n2)
    #  (iii) residual roots R           genuinely extra; generically
    #                                   irrational; the "sporadic" zeros are
    #                                   exactly the rational residual roots.
    root_families = []
    for n1 in range(1, 7):
        for n2 in range(n1, 7):
            for l in range(0, min(n1, n2)):
                P, D, t = poly_in_t(n1, n2, l)
                info = {"n1": n1, "n2": n2, "l": l, "degree": D}
                # (i) rate coincidence
                info["t0_root"] = bool(P.eval(0) == 0)
                info["abs_dn_ge_2"] = abs(n1 - n2) >= 2
                info["coincidence_ok"] = info["t0_root"] == info["abs_dn_ge_2"]
                # (ii) orthogonality root (only meaningful n1 != n2)
                if n1 != n2:
                    t_orth = sp.Rational(n2 - n1, n2 + n1)
                    info["t_orth"] = str(t_orth)
                    info["orth_is_root"] = bool(P.eval(t_orth) == 0)
                else:
                    info["t_orth"] = None
                    info["orth_is_root"] = None
                # (iii) residual: strip t^a and (t - t_orth), report roots
                Pt = sp.Poly(P.as_expr(), t)
                # strip rate-coincidence factor(s)
                while Pt.degree() >= 1 and Pt.eval(0) == 0:
                    Pt = sp.div(Pt, sp.Poly(t, t))[0]
                if n1 != n2:
                    t_orth = sp.Rational(n2 - n1, n2 + n1)
                    while Pt.degree() >= 1 and Pt.eval(t_orth) == 0:
                        Pt = sp.div(Pt, sp.Poly(t - t_orth, t))[0]
                # rational residual roots only (ground_roots is fast; it does
                # NOT attempt to solve the generically-irrational residual).
                grat = Pt.ground_roots() if Pt.degree() >= 1 else {}
                rat_interior = [r for r in grat if sp.Abs(r) < 1]
                info["residual_degree"] = Pt.degree()
                info["residual_rational_roots"] = [str(r) for r in rat_interior]
                info["residual_has_rational_root"] = len(rat_interior) > 0
                root_families.append(info)
    out["root_families"] = {
        "rows": root_families,
        "all_coincidence_ok": all(r["coincidence_ok"] for r in root_families),
        "all_orth_roots": all(r["orth_is_root"] for r in root_families
                              if r["orth_is_root"] is not None),
        "rational_residual_cases": [
            {"n1": r["n1"], "n2": r["n2"], "l": r["l"],
             "residual_rational_roots": r["residual_rational_roots"]}
            for r in root_families if r["residual_has_rational_root"]],
    }

    with open("debug/data/sprint_topos4_vanishing_lemma.json", "w") as fh:
        json.dump(out, fh, indent=1)

    # ---- console summary ----
    print("PINS:", {k: v["match"] for k, v in out["pins"].items()})
    b = out["boundary_n5n5l3"]
    print(f"\nBOUNDARY n=n'=5,l=3:  P(t) = {b['P(t)']}")
    print(f"  degree {b['degree']}  interior roots {b['interior_roots']}")
    print(f"  Jacobi (a,b,scalar) = {b['jacobi_match (a,b,scalar)']}")
    print("  t at (Z=1,Z'=2) = {};  P(t_boundary) = {};  S vanishes = {}".format(
        b['t_at_Z1_Zp2'], b['P(t_boundary)'], b['vanishes']))
    print(f"\nRATE-COINCIDENCE CRITERION (t=0 root <=> |dn|>=2):  "
          f"all consistent = {out['rate_coincidence_criterion']['all_consistent']}")
    print("\nJACOBI FAMILY:")
    for j in out["jacobi_family"]:
        print(f"  ({j['n1']},{j['n2']},l={j['l']}) deg {j['degree']}  "
              f"Jacobi(a,b)={j['jacobi (a,b)']}  roots={j['interior_roots']}")
    rf = out["root_families"]
    print(f"\nROOT FAMILIES (n,l grid):")
    print(f"  (i)  t=0 root  <=> |dn|>=2 :  all consistent = {rf['all_coincidence_ok']}")
    print(f"  (ii) t=(n2-n1)/(n2+n1) is a root (n1!=n2, orthogonality):  "
          f"all hold = {rf['all_orth_roots']}")
    print(f"  (iii) residual roots after stripping (i)+(ii):")
    print(f"     rational-residual cases (the genuine 'sporadic' zeros):")
    for c in rf["rational_residual_cases"]:
        print(f"       (n1={c['n1']},n2={c['n2']},l={c['l']}): "
              f"residual rational roots {c['residual_rational_roots']}")
    izc = out["integer_zero_census"]
    print(f"\nINTEGER ZERO CENSUS (Z,n<{N}):  {izc['total_zeros']} total = "
          f"{izc['rate_coincidence']} rate-coincidence + "
          f"{izc['orthogonality']} orthogonality + "
          f"{izc['residual_count']} residual")
    for r in izc["residual_zeros"]:
        print(f"    RESIDUAL (genuine sporadic): Z1={r['Z1']},n1={r['n1']} | "
              f"Z2={r['Z2']},n2={r['n2']}  l={r['l']}  t={r['t']}")


if __name__ == "__main__":
    main()
