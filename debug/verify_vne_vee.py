"""Verify Vne and Vee closed forms in (s,t,u) Hylleraas integrals."""
import sympy as sp


def verify_vne(L, M, N):
    """Vne integral check: integrand has s^L t^(2M) u^N exp(-2 alpha s) with
    1/r_1 + 1/r_2 = 4 s/(s^2 - t^2). Volume = 8 pi^2 u (s^2 - t^2). Combined:
    Integrand_total = -Z * 8 pi^2 * 4 * s^(L+1) t^(2M) u^(N+1) e^(-2 alpha s)
                       (no s^2 - t^2 factor).
    Integral over volume:
      8 pi^2 * (-4 Z) * integral exp s^(L+1) t^(2M) u^(N+1) ds dt du over s>=u>=|t|
    """
    s, t, u, alpha = sp.symbols('s t u alpha', positive=True)
    integrand = sp.exp(-2*alpha*s) * s**(L+1) * t**(2*M) * u**(N+1)
    I_t = sp.integrate(integrand, (t, -u, u))
    I_u = sp.integrate(I_t, (u, 0, s))
    I_s = sp.integrate(I_u, (s, 0, sp.oo))
    sym = sp.simplify(I_s)

    # Claim: 2 (L+N+2M+4)! / [(2M+1)(N+2M+3)(2 alpha)^(L+N+2M+5)]
    claim = (2 * sp.factorial(L + N + 2*M + 4)
             / ((2*M + 1) * (N + 2*M + 3) * (2*alpha)**(L + N + 2*M + 5)))
    claim = sp.simplify(claim)
    return sym, claim, sp.simplify(sym - claim)


def verify_vee(L, M, N):
    """Vee integral check: integrand s^L t^(2M) u^(N-1) e^(-2 alpha s) with
    volume u (s^2 - t^2) -> Integrand = s^L t^(2M) u^N (s^2 - t^2) e^(-2 alpha s).
    Note the lack of u beyond u^N: 1/u from V_ee cancels one u from volume.
    Claim: 2 (L+N+2M+4)! [1/((2M+1)(N+2M+2)) - 1/((2M+3)(N+2M+4))] / (2 alpha)^(L+N+2M+5)
    """
    s, t, u, alpha = sp.symbols('s t u alpha', positive=True)
    integrand = sp.exp(-2*alpha*s) * s**L * t**(2*M) * u**N * (s**2 - t**2)
    I_t = sp.integrate(integrand, (t, -u, u))
    I_u = sp.integrate(I_t, (u, 0, s))
    I_s = sp.integrate(I_u, (s, 0, sp.oo))
    sym = sp.simplify(I_s)

    bracket = sp.Rational(1, (2*M + 1)*(N + 2*M + 2)) - sp.Rational(1, (2*M + 3)*(N + 2*M + 4))
    claim = (2 * bracket * sp.factorial(L + N + 2*M + 4) / (2*alpha)**(L + N + 2*M + 5))
    claim = sp.simplify(claim)
    return sym, claim, sp.simplify(sym - claim)


print("=== Vne integrals (without 8 pi^2 * (-4Z) prefactor) ===")
for (L, M, N) in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 1, 1)]:
    sym, claim, diff = verify_vne(L, M, N)
    print(f"L={L},M={M},N={N}: sym = {sym}, claim = {claim}, diff = {diff}")

print("\n=== Vee integrals (without 8 pi^2 prefactor) ===")
for (L, M, N) in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 1, 1)]:
    sym, claim, diff = verify_vee(L, M, N)
    print(f"L={L},M={M},N={N}: sym = {sym}, claim = {claim}, diff = {diff}")
