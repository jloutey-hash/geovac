"""
Cross-check alpha-X's sympy-computed S^5 Seeley-DeWitt coefficients against
Kluth-Litim 2020 (EPJ C 80:269, arXiv:1910.00543) Table 1.

Kluth-Litim conventions (KL):
    KL Eq. 1:    dU/dt = (nabla^2 + E) U  (sign on E: POSITIVE)
    KL Eq. 11:   Tr_s U_E(t,sigma) = Vol/(4pi t)^{d/2} * sum_n [b^{(s)}_{2n}(E) t^n + c^{(s)}_{d+2n}(E) t^{d/2+n}]
    KL Eq. 9:    b^{(s)}_n(E) = (1/Vol) * Tr_s[tilde_b_n(E)]
    -- so KL's b coefficients are VOLUME-NORMALIZED, coincident-limit, trace-reduced.

Thus, alpha-X's a_k coefficients (which INCLUDE a Vol factor) correspond to
    alpha_X_b_k := alpha_X_a_k / Vol = dim_S * (scalar invariant)

For comparison, this script:
  1. Reads KL Table 1 entries for d=5, scalar Laplacian, E=0.
  2. Transforms scalar values to E != 0 via KL Eq. 13.
  3. For the Dirac operator D^2, identifies E_KL via the Lichnerowicz formula
     D^2 = -nabla^2 + R/4, noting sign convention: exp(-t D^2) = e^{t(nabla^2 - R/4)}
     so E_KL_Dirac = -R/4.
  4. Multiplies by dim_spin=4 (5D spinor dimension).
  5. Reports the match or any mismatch at b_0, b_2, b_4.

Output: comparison table + verdict.
"""
from sympy import Rational, symbols, simplify, Integer


def main():
    R_sym = symbols('R', positive=True)

    # ---------------------------------------------------
    # 1) KL Table 1 scalar values at d=5, E=0
    # ---------------------------------------------------
    # Directly read from KL Table 1 at d=5:
    b0_E0 = Integer(1)
    b2_E0 = R_sym / Rational(6)
    b4_E0 = R_sym**2 / Rational(75)

    print("=" * 68)
    print("KL Table 1 (scalar Laplacian, E=0) at d=5:")
    print("=" * 68)
    print(f"  b_0^(0) = {b0_E0}")
    print(f"  b_2^(0) = {b2_E0}")
    print(f"  b_4^(0) = {b4_E0}")
    print()

    # ---------------------------------------------------
    # 2) alpha-X scalar values, Vol-normalized
    # ---------------------------------------------------
    # alpha-X computed:
    #   a_0^scalar = pi^3 = Vol(S^5_{r=1})
    #   a_2^scalar = 10 pi^3 / 3
    #   a_4^scalar = 16 pi^3 / 3
    # Dividing by Vol = pi^3, these become:
    #   b_0^scalar_aX = 1
    #   b_2^scalar_aX = 10/3
    #   b_4^scalar_aX = 16/3
    alphaX_b0_scalar = Integer(1)
    alphaX_b2_scalar = Rational(10, 3)
    alphaX_b4_scalar = Rational(16, 3)

    print("=" * 68)
    print("alpha-X scalar values on unit S^5 (Vol-divided, R_scalar=20):")
    print("=" * 68)
    print(f"  b_0^scalar_aX = {alphaX_b0_scalar}")
    print(f"  b_2^scalar_aX = {alphaX_b2_scalar}")
    print(f"  b_4^scalar_aX = {alphaX_b4_scalar}")
    print()

    # KL values evaluated on unit S^5 (R=20, since R_scalar = d(d-1) = 20 at d=5, r=1)
    kl_b0_scalar_unit = b0_E0
    kl_b2_scalar_unit = b2_E0.subs(R_sym, 20)
    kl_b4_scalar_unit = b4_E0.subs(R_sym, 20)

    print("=" * 68)
    print("KL values substituted to unit S^5 (R=20):")
    print("=" * 68)
    print(f"  b_0^(0)_KL = {kl_b0_scalar_unit}")
    print(f"  b_2^(0)_KL = {kl_b2_scalar_unit}")
    print(f"  b_4^(0)_KL = {kl_b4_scalar_unit}")
    print()

    # Check scalar match
    scalar_match = (
        simplify(alphaX_b0_scalar - kl_b0_scalar_unit) == 0
        and simplify(alphaX_b2_scalar - kl_b2_scalar_unit) == 0
        and simplify(alphaX_b4_scalar - kl_b4_scalar_unit) == 0
    )
    print(f"SCALAR LAPLACIAN MATCH: {scalar_match}")
    print()

    # ---------------------------------------------------
    # 3) Dirac D^2 via KL Eq. 13: transform E=0 -> E_KL = -R/4
    # ---------------------------------------------------
    # KL's E convention: dU/dt = (nabla^2 + E) U, i.e., operator = -(nabla^2 + E) = -nabla^2 - E.
    # Dirac squared via Lichnerowicz: D^2 = -nabla^2 + R/4
    # So D^2 = -nabla^2 - (-R/4), i.e., P = -nabla^2 - E with E = -R/4.
    # Therefore E_KL_Dirac = -R/4.
    E_dirac_KL = -R_sym / Rational(4)

    # KL Eq. 13:
    #   b_2n^(s)(E) = sum_{k=0}^n (E - E_bar)^k / k! * b_{2(n-k)}^(s)(E_bar)
    # With E_bar = 0 (Table 1 reference), E = E_KL_Dirac = -R/4:
    b0_dirac_scalar = b0_E0
    b2_dirac_scalar = E_dirac_KL * b0_E0 + b2_E0
    b4_dirac_scalar = (E_dirac_KL**2 / 2) * b0_E0 + E_dirac_KL * b2_E0 + b4_E0

    b2_dirac_scalar_s = simplify(b2_dirac_scalar)
    b4_dirac_scalar_s = simplify(b4_dirac_scalar)

    print("=" * 68)
    print("KL-derived Dirac (scalar part) on S^5 via E_KL = -R/4:")
    print("=" * 68)
    print(f"  b_0^Dirac_scalar = {b0_dirac_scalar}")
    print(f"  b_2^Dirac_scalar = {b2_dirac_scalar_s}")
    print(f"  b_4^Dirac_scalar = {b4_dirac_scalar_s}")
    print()

    # Multiply by dim_spin = 4 (5D Dirac spinor dimension)
    dim_spin_5D = 4
    kl_b0_dirac = dim_spin_5D * b0_dirac_scalar
    kl_b2_dirac = simplify(dim_spin_5D * b2_dirac_scalar)
    kl_b4_dirac = simplify(dim_spin_5D * b4_dirac_scalar)

    print(f"Multiplied by dim_spin = 4 (general R_sym):")
    print(f"  b_0^Dirac_KL = {kl_b0_dirac}")
    print(f"  b_2^Dirac_KL = {kl_b2_dirac}")
    print(f"  b_4^Dirac_KL = {kl_b4_dirac}")
    print()

    # Substitute R=20 for unit S^5
    kl_b0_dirac_unit = kl_b0_dirac
    kl_b2_dirac_unit = kl_b2_dirac.subs(R_sym, 20)
    kl_b4_dirac_unit = kl_b4_dirac.subs(R_sym, 20)

    print(f"Substituted to unit S^5 (R=20):")
    print(f"  b_0^Dirac_KL = {kl_b0_dirac_unit}")
    print(f"  b_2^Dirac_KL = {simplify(kl_b2_dirac_unit)}")
    print(f"  b_4^Dirac_KL = {simplify(kl_b4_dirac_unit)}")
    print()

    # ---------------------------------------------------
    # 4) alpha-X Dirac values, Vol-normalized
    # ---------------------------------------------------
    alphaX_b0_dirac = Integer(4)
    alphaX_b2_dirac = Rational(-20, 3)
    alphaX_b4_dirac = Rational(14, 3)

    print("=" * 68)
    print("alpha-X Dirac values on unit S^5 (Vol-divided):")
    print("=" * 68)
    print(f"  b_0^Dirac_aX = {alphaX_b0_dirac}")
    print(f"  b_2^Dirac_aX = {alphaX_b2_dirac}")
    print(f"  b_4^Dirac_aX = {alphaX_b4_dirac}")
    print()

    # Check Dirac match
    dirac_b0_match = simplify(alphaX_b0_dirac - kl_b0_dirac_unit) == 0
    dirac_b2_match = simplify(alphaX_b2_dirac - kl_b2_dirac_unit) == 0
    dirac_b4_match = simplify(alphaX_b4_dirac - kl_b4_dirac_unit) == 0

    print(f"DIRAC OPERATOR MATCH:")
    print(f"  b_0: {dirac_b0_match}")
    print(f"  b_2: {dirac_b2_match}")
    print(f"  b_4: {dirac_b4_match}")
    print()

    # ---------------------------------------------------
    # 5) Final verdict
    # ---------------------------------------------------
    print("=" * 68)
    print("FINAL VERDICT")
    print("=" * 68)
    if scalar_match and dirac_b0_match and dirac_b2_match and dirac_b4_match:
        print("MATCH: All six values (b_0, b_2, b_4 for both scalar and Dirac)")
        print("match Kluth-Litim Table 1 (and its KL Eq.13 transformation to")
        print("E_KL = -R/4 for the Dirac Lichnerowicz endomorphism).")
        print()
        print("alpha-X's sympy computation INDEPENDENTLY VERIFIED.")
    elif scalar_match:
        print("PARTIAL: Scalar matches, Dirac differs -- check sign convention")
        print("on Lichnerowicz endomorphism or spinor dimension.")
    else:
        print("MISMATCH: further investigation needed.")


if __name__ == "__main__":
    main()
