"""Wider PSLQ search for S^7 scalar zeta'(0) closed form."""

import mpmath as mp
import sympy as sp

mp.mp.dps = 200

# Recompute the value (same as in main script)
def zeta_S7_conf_prime(k_max: int = 100):
    leading = 2 * (mp.zeta(-6, derivative=1) - 5 * mp.zeta(-4, derivative=1) + 4 * mp.zeta(-2, derivative=1))
    series_sum = mp.mpf(0)
    for k in range(1, k_max + 1):
        term = (mp.zeta(2*k - 6) - 5 * mp.zeta(2*k - 4) + 4 * mp.zeta(2*k - 2)) / (k * mp.mpf(4)**k)
        series_sum += term
    return (leading + series_sum) / 360

target = zeta_S7_conf_prime()
print(f"Target value: {mp.nstr(target, 50)}")
print()

# Try several increasingly broad bases
bases_to_try = [
    # Try with rational scaling - maybe coefficients are over /360 or /5760 denominators
    (["log(2)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6"],
     [mp.log(2), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6]),
    # With 1 to catch rational constants
    (["log(2)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6", "1"],
     [mp.log(2), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6, mp.mpf(1)]),
    # With pi^2, pi^4, pi^6 in case those don't combine cleanly
    (["log(2)", "zeta(3)", "zeta(5)", "zeta(7)", "pi^2", "pi^4", "pi^6"],
     [mp.log(2), mp.zeta(3), mp.zeta(5), mp.zeta(7), mp.pi**2, mp.pi**4, mp.pi**6]),
    # With even higher zeta values
    (["log(2)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6", "zeta(9)/pi^8"],
     [mp.log(2), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6, mp.zeta(9)/mp.pi**8]),
]

for names, basis in bases_to_try:
    print(f"Testing basis: {names}")
    vec = [target] + basis
    for tol_exp in (30, 50, 80):
        for maxc in (10**8, 10**12, 10**16):
            r = mp.pslq(vec, tol=mp.mpf(f'1e-{tol_exp}'), maxcoeff=maxc)
            if r is not None and r[0] != 0:
                print(f"  tol=1e-{tol_exp}, maxcoeff=10^{int(mp.log10(maxc))}: {r}")
                # Decode
                a0 = r[0]
                print(f"    framework = ", end="")
                terms = []
                for i, n in enumerate(names):
                    rat = sp.Rational(-r[i+1], a0)
                    if rat != 0:
                        terms.append(f"({rat}) * {n}")
                print(" + ".join(terms))
                break
        if r is not None and r[0] != 0:
            break
    else:
        print(f"  no relation found")
    print()
