"""Track RH-F: Arithmetic / Galois structure of the Ihara zeros.

For every non-trivial polynomial factor appearing in the closed-form Ihara
zetas of the GeoVac Hopf graphs (scalar S^3 Coulomb at max_n=2,3; scalar
S^5 Bargmann-Segal at N_max=2,3; Dirac-S^3 Rule A / Rule B at n_max=1,2,3),
this script computes:

  (a) irreducible factorization over Q, Q(i), Q(sqrt 2), Q(sqrt 5),
      Q(omega)=Q(sqrt -3), and Q(sqrt disc) where disc is the polynomial's
      discriminant;
  (b) Galois group of the polynomial (if irreducible over Q);
  (c) minimal splitting field degree over Q;
  (d) ring-classification of the zero set of each closed form.

Run: python debug/compute_galois_ihara.py
Output: debug/data/galois_ihara.json
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
from sympy.polys.numberfields import galois_group

s = sp.symbols('s')

# ---------------------------------------------------------------------------
# Closed-form factorizations: polynomial -> symbolic expression.
# Extracted from debug/data/ihara_zeta_geovac_hopf.json and
# debug/data/ihara_zeta_dirac_s3.json. Trivial (s±1) factors are dropped:
# those contribute s = ±1 in Q by inspection.
# ---------------------------------------------------------------------------

CLOSED_FORMS = {
    # Scalar S^3 Coulomb
    "S3_max_n_3_cyc_minus": s**2 - s + 1,
    "S3_max_n_3_cyc_plus": s**2 + s + 1,
    "S3_max_n_3_cubic_minus": 2*s**3 - s**2 + s - 1,
    "S3_max_n_3_cubic_plus": 2*s**3 + s**2 + s + 1,
    # Scalar S^5 Bargmann-Segal N_max=2
    "S5_N2_quad": 2*s**2 + 1,
    "S5_N2_biquad": 3*s**4 + 3*s**2 + 1,
    "S5_N2_sextic": 24*s**6 + 21*s**4 + s**2 - 1,
    # Scalar S^5 Bargmann-Segal N_max=3
    "S5_N3_P12": (432*s**12 + 666*s**10 + 374*s**8 + 135*s**6
                  + 47*s**4 + 11*s**2 + 1),
    "S5_N3_P22": (829440*s**22 + 2453184*s**20 + 3308104*s**18
                  + 2696682*s**16 + 1470640*s**14 + 557227*s**12
                  + 146654*s**10 + 25709*s**8 + 2630*s**6
                  + 79*s**4 - 12*s**2 - 1),
    # Dirac-S^3 Rule A n_max=2
    "DiracA_n2_quad": s**2 + 1,
    # Dirac-S^3 Rule B n_max=2
    "DiracB_n2_3s2p1": 3*s**2 + 1,
    "DiracB_n2_4s2p1": 4*s**2 + 1,
    "DiracB_n2_lin_minus": 2*s**2 - s + 1,
    "DiracB_n2_lin_plus": 2*s**2 + s + 1,
    "DiracB_n2_biquad": 12*s**4 + 9*s**2 - 1,
    # Dirac-S^3 Rule A n_max=3
    "DiracA_n3_quad": s**2 + 1,  # reappears
    "DiracA_n3_cyc_minus": s**2 - s + 1,
    "DiracA_n3_cyc_plus": s**2 + s + 1,
    "DiracA_n3_cubic_a": 2*s**3 - 2*s**2 + 2*s - 1,
    "DiracA_n3_cubic_b": 2*s**3 - s**2 + s - 1,
    "DiracA_n3_cubic_c": 2*s**3 + s**2 + s + 1,
    "DiracA_n3_cubic_d": 2*s**3 + 2*s**2 + 2*s + 1,
    "DiracA_n3_quartic_a": 2*s**4 - 2*s**3 + 2*s**2 - s + 1,
    "DiracA_n3_quartic_b": 2*s**4 + 2*s**3 + 2*s**2 + s + 1,
    # Dirac-S^3 Rule B n_max=3
    "DiracB_n3_9s2p1": 9*s**2 + 1,
    "DiracB_n3_P22": (538876800*s**22 + 1088750160*s**20 + 925413984*s**18
                      + 456734400*s**16 + 148182464*s**14 + 33353711*s**12
                      + 5288330*s**10 + 580341*s**8 + 41332*s**6
                      + 1589*s**4 + 10*s**2 - 1),
    "DiracB_n3_P24": (538876800*s**24 + 1108028160*s**22 + 934113744*s**20
                      + 450036324*s**18 + 145121264*s**16 + 35278007*s**14
                      + 7106943*s**12 + 1222735*s**10 + 170871*s**8
                      + 17793*s**6 + 1257*s**4 + 53*s**2 + 1),
}

# Candidate extension fields: Q, Q(i), Q(sqrt 2), Q(sqrt 5), Q(omega)=Q(sqrt -3),
# plus Q(sqrt 3) and Q(sqrt -1,-2)-composite for diagnostic completeness.
EXTENSIONS = {
    "Q": [],
    "Q(i)": [sp.I],
    "Q(sqrt 2)": [sp.sqrt(2)],
    "Q(sqrt 3)": [sp.sqrt(3)],
    "Q(sqrt 5)": [sp.sqrt(5)],
    "Q(sqrt -3) = Q(omega)": [sp.sqrt(-3)],
    "Q(sqrt -7)": [sp.sqrt(-7)],
    "Q(zeta_8)": [sp.sqrt(2), sp.I],
    "Q(zeta_12)": [sp.sqrt(3), sp.I],
}


def factor_over(extension, poly_expr):
    """Factor poly over Q adjoined the given extension generators.

    Returns (expr, n_poly_factors) where n_poly_factors counts the number of
    distinct polynomial-of-s factors (ignoring scalar leading coefficients).
    """
    try:
        if not extension:
            fexpr = sp.factor(poly_expr)
        else:
            fexpr = sp.factor(poly_expr, extension=extension)
    except Exception as e:
        return (f"ERROR: {e!s}", None)
    # Count polynomial-of-s factors, with multiplicity
    n_poly = 0
    for a in sp.Mul.make_args(fexpr):
        # a may be Pow(base, e)
        base, exp = a.as_base_exp()
        try:
            deg = sp.Poly(base, s).degree()
        except (sp.PolynomialError, sp.GeneratorsNeeded):
            deg = 0
        if deg > 0:
            # Count with multiplicity
            try:
                mult = int(exp)
            except TypeError:
                mult = 1
            n_poly += mult
    return (str(fexpr), n_poly)


def even_in_s_reduce(poly_expr):
    """If poly(s) is even in s, return poly(u) with u = s^2, else None."""
    poly = sp.Poly(poly_expr, s)
    coeffs = poly.all_coeffs()  # highest-to-lowest
    deg = poly.degree()
    if deg % 2 != 0:
        return None
    # Odd powers must all be zero
    for i, c in enumerate(coeffs):
        power = deg - i
        if power % 2 == 1 and c != 0:
            return None
    u = sp.symbols('u')
    even_coeffs = coeffs[::2]
    expr = sum(c * u**(len(even_coeffs) - 1 - i) for i, c in enumerate(even_coeffs))
    return expr


def compute_galois_data(name, poly_expr):
    """Compute Galois group, discriminant, degree, and factorizations."""
    poly = sp.Poly(poly_expr, s)
    degree = poly.degree()
    disc = sp.discriminant(poly_expr, s)

    # Square-free part of disc gives the quadratic resolvent field
    disc_rational = sp.Rational(disc)
    sf_disc = sp.sqrtdenest(sp.sqrt(disc_rational))
    # Extract square-free part explicitly
    try:
        sqfree = sp.core.numbers.Integer(disc).p
        # Factor into primes
        # Prefer sp.integer_nthroot
    except Exception:
        sqfree = int(disc)

    result = {
        "name": name,
        "polynomial": str(poly_expr),
        "degree": int(degree),
        "discriminant": int(disc) if disc.is_Integer else str(disc),
        "factorizations": {},
        "galois_group": None,
        "galois_order": None,
        "solvable": None,
        "is_abelian": None,
        "splitting_field_degree": None,
        "irreducible_over_Q": None,
    }

    # Irreducibility over Q: count s-polynomial factors (with multiplicity)
    _, n_poly_q = factor_over([], poly_expr)
    # Irreducible iff exactly one s-polynomial factor with multiplicity 1
    # Check via degree equality: compare single-factor degree to poly degree
    factored_q_expr = sp.factor(poly_expr)
    single_factor_deg = 0
    for a in sp.Mul.make_args(factored_q_expr):
        base, exp = a.as_base_exp()
        try:
            d = sp.Poly(base, s).degree()
        except (sp.PolynomialError, sp.GeneratorsNeeded):
            d = 0
        if d > single_factor_deg:
            single_factor_deg = d
    irreducible_over_Q = (single_factor_deg == degree) and (n_poly_q == 1)
    result["irreducible_over_Q"] = bool(irreducible_over_Q)

    # Factorization over each extension (skip expensive extensions for
    # degree > 12 polynomials; the u-reduced form is more informative).
    result["n_factors_over"] = {}
    heavy_thresh = 12
    for ext_name, ext in EXTENSIONS.items():
        if degree > heavy_thresh and ext_name in {
            "Q(zeta_8)", "Q(zeta_12)", "Q(sqrt -7)", "Q(sqrt 3)",
            "Q(sqrt 5)",
        }:
            result["factorizations"][ext_name] = "SKIPPED (deg > 12)"
            result["n_factors_over"][ext_name] = None
            continue
        factored, n_poly = factor_over(ext, poly_expr)
        result["factorizations"][ext_name] = str(factored)
        result["n_factors_over"][ext_name] = n_poly

    # Factor over Q(sqrt disc) — the quadratic resolvent, expected to
    # partially split a cubic with S_3 Galois group.
    # Skip if disc is huge (> 10^15) — factoring over Q(sqrt D) for giant D
    # is prohibitively slow in sympy and the result tells us little beyond
    # the fact that the cubic's Galois sub-action is the sign of D.
    if disc.is_Integer:
        try:
            disc_int = int(disc)
            if abs(disc_int) < 10**15:
                sign = -1 if disc_int < 0 else 1
                a = abs(disc_int)
                # extract square part
                sq_free = a
                k = 2
                while k * k <= sq_free:
                    while sq_free % (k * k) == 0:
                        sq_free //= k * k
                    k += 1
                sq_free = sign * sq_free
                if sq_free not in {1, -1}:
                    ext_gen = sp.sqrt(sp.Integer(sq_free))
                    factored, n_poly = factor_over([ext_gen], poly_expr)
                    result["factorizations"][f"Q(sqrt disc = sqrt {sq_free})"] = str(factored)
                    result["n_factors_over"][f"Q(sqrt disc = sqrt {sq_free})"] = n_poly
                    result["sq_free_disc"] = sq_free
                else:
                    result["sq_free_disc"] = sq_free
            else:
                result["sq_free_disc"] = "disc too large (skipped)"
        except Exception as e:
            result["sq_free_disc_error"] = str(e)

    # Galois group (only meaningful if irreducible over Q; deg <= 6)
    if irreducible_over_Q and degree <= 6:
        try:
            gg, is_alt = galois_group(poly)
            result["galois_group"] = str(gg)
            result["galois_order"] = int(gg.order())
            result["solvable"] = bool(gg.is_solvable)
            result["is_abelian"] = bool(gg.is_abelian)
            # Splitting field degree = |Gal| for irreducible polynomials
            result["splitting_field_degree"] = int(gg.order())
        except Exception as e:
            result["galois_group"] = f"ERROR: {e!s}"
    elif irreducible_over_Q and degree > 6:
        result["galois_group"] = f"N/A (sympy limit, deg > 6)"

    # Even-in-s reduction: if P(s) = Q(s^2), compute Galois on Q(u)
    u_reduced = even_in_s_reduce(poly_expr)
    if u_reduced is not None:
        u_deg = sp.Poly(u_reduced, sp.symbols('u')).degree()
        result["even_in_s"] = True
        result["u_reduced_polynomial"] = str(u_reduced)
        result["u_reduced_degree"] = int(u_deg)
        if u_deg <= 6:
            try:
                u_sym = sp.symbols('u')
                u_poly = sp.Poly(u_reduced, u_sym)
                # Check irreducibility of u-polynomial
                u_factored = sp.factor(u_reduced)
                u_single_deg = 0
                for a in sp.Mul.make_args(u_factored):
                    base, exp = a.as_base_exp()
                    try:
                        d = sp.Poly(base, u_sym).degree()
                    except (sp.PolynomialError, sp.GeneratorsNeeded):
                        d = 0
                    if d > u_single_deg:
                        u_single_deg = d
                u_irreducible = (u_single_deg == u_deg)
                result["u_irreducible_over_Q"] = bool(u_irreducible)
                if u_irreducible:
                    gg_u, _ = galois_group(u_poly)
                    result["u_galois_order"] = int(gg_u.order())
                    result["u_galois_group"] = str(gg_u)
                    result["u_solvable"] = bool(gg_u.is_solvable)
                    result["u_abelian"] = bool(gg_u.is_abelian)
                    # Splitting field of P(s) = P(u)(u=s^2) over splitting field
                    # of P(u) is degree 2^k where k is #distinct roots needing
                    # a sqrt. Upper bound: 2 * |Gal(P_u)|.
                    result["s_splitting_field_upper_bound"] = 2 * int(gg_u.order())
            except Exception as e:
                result["u_galois_error"] = str(e)
    else:
        result["even_in_s"] = False

    return result


def cyclotomic_check(poly_expr):
    """Detect if a polynomial is a cyclotomic polynomial (and which one)."""
    poly = sp.Poly(poly_expr, s)
    if poly.LC() != 1:
        return None
    # Compare against phi_n for small n
    for n_cyc in range(1, 60):
        phi = sp.cyclotomic_poly(n_cyc, s)
        if sp.simplify(poly.as_expr() - phi) == 0:
            return n_cyc
    return None


def main():
    import sys
    # Force line-buffered stdout
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass
    out = {
        "track": "RH-F",
        "description": ("Galois / number-field classification of Ihara-zero "
                        "polynomials from Paper 29 closed forms"),
        "candidate_fields": list(EXTENSIONS.keys()),
        "polynomials": {},
    }

    print("=" * 72)
    print("TRACK RH-F: GALOIS STRUCTURE OF IHARA ZEROS")
    print("=" * 72)

    for name, poly_expr in CLOSED_FORMS.items():
        print(f"\n--- {name} ---")
        print(f"  polynomial: {poly_expr}")
        data = compute_galois_data(name, poly_expr)
        cyc_n = cyclotomic_check(poly_expr)
        if cyc_n is not None:
            data["cyclotomic_n"] = cyc_n
            print(f"  CYCLOTOMIC: Phi_{cyc_n}(s)")
        else:
            data["cyclotomic_n"] = None
        print(f"  disc = {data['discriminant']}")
        print(f"  deg  = {data['degree']}")
        print(f"  irreducible over Q: {data['irreducible_over_Q']}")
        if data["galois_order"] is not None:
            print(f"  Galois order = {data['galois_order']} "
                  f"(solvable: {data['solvable']}, abelian: {data['is_abelian']})")
        # Print where the polynomial first splits (by n_poly_factors)
        n_q = data["n_factors_over"]["Q"]
        first_split_found = False
        for ext_name in ["Q(i)", "Q(sqrt 2)", "Q(sqrt 3)", "Q(sqrt 5)",
                         "Q(sqrt -3) = Q(omega)", "Q(sqrt -7)",
                         "Q(zeta_8)", "Q(zeta_12)"]:
            n_ext = data["n_factors_over"].get(ext_name)
            if n_ext is not None and n_ext > n_q:
                print(f"  SPLITS over {ext_name}: {n_q} -> {n_ext} factors")
                print(f"    {data['factorizations'][ext_name]}")
                first_split_found = True
                break
        if not first_split_found and data["irreducible_over_Q"]:
            print("  (remains irreducible over all tested quadratic extensions)")
        out["polynomials"][name] = data

    # Minimal splitting field per zero set
    # Zero set of a polynomial P with Gal(P/Q) = G: splitting field has
    # degree |G| over Q (if P irreducible). Composite splitting fields for
    # products use the compositum.
    out["summary"] = {
        "cyclotomic_factors": [
            name for name, data in out["polynomials"].items()
            if data.get("cyclotomic_n") is not None
        ],
        "abelian_factors": [
            name for name, data in out["polynomials"].items()
            if data.get("is_abelian") is True
        ],
        "S3_galois_factors": [
            name for name, data in out["polynomials"].items()
            if data.get("galois_order") == 6 and data.get("solvable")
        ],
        "factors_splitting_over_Qi": [],
        "factors_splitting_over_Qsqrt2": [],
    }

    for name, data in out["polynomials"].items():
        n_q = data["n_factors_over"]["Q"]
        for ext_name in EXTENSIONS.keys():
            if ext_name == "Q":
                continue
            n_ext = data["n_factors_over"].get(ext_name)
            if n_ext is not None and n_ext > n_q:
                key = f"factors_splitting_over_{ext_name}"
                if key not in out["summary"]:
                    out["summary"][key] = []
                out["summary"][key].append(name)

    # Write JSON
    out_path = Path("debug/data/galois_ihara.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)

    print(f"\n\nWrote {out_path}")


if __name__ == "__main__":
    main()
