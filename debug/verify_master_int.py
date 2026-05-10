"""Verify the Hylleraas master integral closed form symbolically."""
import sympy as sp


def verify_master(l, m, n, alpha_val=2):
    """Compare symbolic vs claimed closed form."""
    s, t, u, alpha = sp.symbols('s t u alpha', positive=True)

    integrand = sp.exp(-2 * alpha * s) * s**l * t**(2*m) * u**n * u * (s**2 - t**2)

    I_t = sp.integrate(integrand, (t, -u, u))
    I_u = sp.integrate(I_t, (u, 0, s))
    I_s = sp.integrate(I_u, (s, 0, sp.oo))

    I_sym = sp.simplify(I_s)

    # Claimed: 4 (n + 4m + 6) (l + n + 2m + 5)! / [(2m+1)(2m+3)(n+2m+3)(n+2m+5) (2 alpha)^(l+n+2m+6)]
    N_total = l + n + 2*m + 5
    claim = (4 * (n + 4*m + 6) * sp.factorial(N_total)
             / ((2*m + 1) * (2*m + 3) * (n + 2*m + 3) * (n + 2*m + 5)
                * (2*alpha)**(N_total + 1)))
    claim = sp.simplify(claim)

    diff = sp.simplify(I_sym - claim)
    return I_sym, claim, diff


for (l, m, n) in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                  (2, 0, 0), (1, 1, 0), (2, 1, 1)]:
    I_sym, claim, diff = verify_master(l, m, n)
    print(f"({l},{m},{n}): I = {I_sym}")
    print(f"          claim = {claim}")
    print(f"          diff  = {diff}")
    if diff != 0:
        print(f"          *** MISMATCH ***")
    print()
