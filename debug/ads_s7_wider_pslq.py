"""Extended PSLQ search for S^7 scalar zeta'(0) closed form.

Round 2 follow-up to v3.2.2 negative (which tested only
{log 2, zeta(3)/pi^2, zeta(5)/pi^4, zeta(7)/pi^6}). The v3.2.3 finding
that the minimal scalar on S^5 introduces log 3 and log 5 from m=2
zero-mode boundary suggests the S^d ring may include log p for primes
p relevant to the multiplicity structure.

Try bases including:
- {log 2, log 3, log 5, zeta(3)/pi^2, zeta(5)/pi^4, zeta(7)/pi^6}
- + log 7
- + log p products
- + Dirichlet beta values (Catalan G, beta(4), beta(6))
- + zeta(9)/pi^8

At least one of these should give a clean integer relation for the S^7
value -0.001594661553469573861754924... computed via Hurwitz at 200 dps.
"""

import mpmath as mp
import sympy as sp

mp.mp.dps = 200


def s7_value(k_max=100):
    """Recompute S^7 scalar zeta'(0) value."""
    leading = 2 * (mp.zeta(-6, derivative=1) - 5 * mp.zeta(-4, derivative=1) + 4 * mp.zeta(-2, derivative=1))
    series_sum = mp.mpf(0)
    for k in range(1, k_max + 1):
        term = (mp.zeta(2*k - 6) - 5 * mp.zeta(2*k - 4) + 4 * mp.zeta(2*k - 2)) / (k * mp.mpf(4)**k)
        series_sum += term
    return (leading + series_sum) / 360


target = s7_value()
print(f"S^7 scalar zeta'(0) target value: {mp.nstr(target, 50)}")
print()


bases = [
    # Round 2a: add log 3
    (["log(2)", "log(3)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6"],
     [mp.log(2), mp.log(3), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6]),
    # Round 2b: log 3 + log 5
    (["log(2)", "log(3)", "log(5)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6"],
     [mp.log(2), mp.log(3), mp.log(5), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6]),
    # Round 2c: + log 7
    (["log(2)", "log(3)", "log(5)", "log(7)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6"],
     [mp.log(2), mp.log(3), mp.log(5), mp.log(7), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6]),
    # Round 2d: + Catalan G
    (["log(2)", "log(3)", "Catalan_G", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6"],
     [mp.log(2), mp.log(3), mp.catalan, mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6]),
    # Round 2e: + zeta(9)/pi^8 and log(15/4) = log 3 + log 5 - 2 log 2 (single combo)
    (["log(2)", "log(15/4)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6", "zeta(9)/pi^8"],
     [mp.log(2), mp.log(mp.mpf(15)/4), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6, mp.zeta(9)/mp.pi**8]),
    # Round 2f: minimal logarithmic basis + zetas
    (["log(2)", "log(3)", "log(5)", "log(7)", "zeta(3)", "zeta(5)", "zeta(7)", "pi^2", "pi^4", "pi^6"],
     [mp.log(2), mp.log(3), mp.log(5), mp.log(7), mp.zeta(3), mp.zeta(5), mp.zeta(7), mp.pi**2, mp.pi**4, mp.pi**6]),
]

for names, basis in bases:
    print(f"Testing basis: {names}")
    vec = [target] + basis
    for tol_exp in (30, 50, 80):
        for maxc in (10**8, 10**12, 10**15):
            try:
                r = mp.pslq(vec, tol=mp.mpf(f'1e-{tol_exp}'), maxcoeff=maxc)
            except (ValueError, TypeError):
                r = None
            if r is not None and r[0] != 0:
                print(f"  tol=1e-{tol_exp}, maxcoeff=10^{int(mp.log10(maxc))}: relation {r}")
                a0 = r[0]
                terms = []
                for i, n in enumerate(names):
                    rat = sp.Rational(-r[i+1], a0)
                    if rat != 0:
                        terms.append(f"({rat}) * {n}")
                print(f"    framework target = " + " + ".join(terms))
                # Verify numerically
                reconstr = sum(sp.Rational(-r[i+1], a0) * basis[i] for i in range(len(basis)))
                diff = abs(target - reconstr)
                print(f"    numerical residual: {mp.nstr(diff, 10)}")
                break
        if r is not None and r[0] != 0:
            break
    else:
        print(f"  NO relation in this basis")
    print()
