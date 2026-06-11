"""Re-verify the triplet T matrix element analytically."""
import sympy as sp

s, t, u, a = sp.symbols('s t u a', positive=True)

# Full integrand for <phi_t|T|phi_t> with phi_t = exp(-alpha s) t.
# T phi = e^(-alpha s) [-alpha^2 t + 4 t/(s^2-t^2) + 4 alpha s t/(s^2-t^2)]
# phi * T phi = e^(-2 alpha s) [-alpha^2 t^2 + 4 t^2/(s^2-t^2) + 4 alpha s t^2/(s^2-t^2)]
# Multiply by volume pi^2 u(s^2-t^2):
# integrand = pi^2 e^(-2 alpha s) u [-alpha^2 t^2 (s^2-t^2) + 4 t^2 + 4 alpha s t^2]

# Without pi^2 (we add it at the end):
integrand = sp.exp(-2*a*s) * u * (
    -a**2 * t**2 * (s**2 - t**2)
    + 4 * t**2
    + 4 * a * s * t**2
)
I_t = sp.integrate(integrand, (t, -u, u))
I_u = sp.integrate(I_t, (u, 0, s))
I_s = sp.integrate(I_u, (s, 0, sp.oo))
T_int = sp.simplify(I_s)
T_full = sp.pi**2 * T_int
print(f"Symbolic T = {T_full}")
print(f"T at alpha=2 = {float(T_full.subs(a, 2))}")
print(f"Expected: 5*pi^2/(2*alpha^6) at alpha=2 = {float(5*sp.pi**2/(2*2**6))}")
