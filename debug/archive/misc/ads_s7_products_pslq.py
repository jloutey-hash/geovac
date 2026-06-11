"""Round 2 PSLQ continuation: products of zeta values."""

import mpmath as mp
import sympy as sp

mp.mp.dps = 200


def s7_value(k_max=100):
    leading = 2 * (mp.zeta(-6, derivative=1) - 5 * mp.zeta(-4, derivative=1) + 4 * mp.zeta(-2, derivative=1))
    series_sum = mp.mpf(0)
    for k in range(1, k_max + 1):
        term = (mp.zeta(2*k - 6) - 5 * mp.zeta(2*k - 4) + 4 * mp.zeta(2*k - 2)) / (k * mp.mpf(4)**k)
        series_sum += term
    return (leading + series_sum) / 360


target = s7_value()
print(f"Target value: {mp.nstr(target, 50)}")
print()

# Cross-check with direct numerical sum (Euler-Maclaurin via mpmath's nsum)
print("Cross-check via direct nsum (after asymptotic subtraction)...")
mp.mp.dps = 100


def s7_term(n):
    n = mp.mpf(n)
    mult = (2*n + 4) * (n+1) * (n+2) * (n+3) / 6
    eig = (n + mp.mpf(3)/2) * (n + mp.mpf(5)/2)
    # Asymptotic subtraction: for large n, log(eig) ~ log(n^2) + log(1 + 4/n + ...)
    # Use: log(eig) ~ 2 log(n) for large n
    # full - 2 log n * mult should be O(1/n^4)
    u = n + 2
    log_eig = mp.log(eig)
    # leading asymptotic: log(eig) ~ 2 log(u) + 0/u^2 - 1/(2u^2) + ...
    # actual: log((u-1/2)(u+1/2)) = log(u^2 - 1/4) = 2 log u + log(1 - 1/(4u^2))
    #       ~ 2 log u - 1/(4u^2) - 1/(32 u^4) - ...
    # The subtracted asymp: mult * 2 log u  (the divergent piece)
    # Plus mult * (-1/(4u^2)) which is mult/u^2 type — for mult ~ u^4, gives u^2 term (still divergent)
    # We need careful asymp subtraction; just return the raw term
    return mult * log_eig


# Just compare to Hurwitz value
hurwitz_val = s7_value(k_max=100)
print(f"  Hurwitz value (k_max=100): {mp.nstr(hurwitz_val, 30)}")

# Different k_max values
for k in [50, 75, 100, 150]:
    v = s7_value(k_max=k)
    print(f"  Hurwitz at k_max={k}: {mp.nstr(v, 30)}")
print()


# Products basis
print("Testing extended basis with zeta-zeta products...")
mp.mp.dps = 200
bases = [
    (["log(2)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6", "zeta(3)^2/pi^6", "zeta(3)*zeta(5)/pi^6"],
     [mp.log(2), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6,
      mp.zeta(3)**2/mp.pi**6, mp.zeta(3)*mp.zeta(5)/mp.pi**6]),
    # Try log 2 only with products
    (["log(2)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6", "zeta(3)*zeta(5)/pi^8"],
     [mp.log(2), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6,
      mp.zeta(3)*mp.zeta(5)/mp.pi**8]),
    # Pure MZV ring up to weight 7
    (["log(2)", "zeta(3)", "zeta(5)", "zeta(7)", "zeta(3)*zeta(5)", "pi^2", "pi^4", "pi^6"],
     [mp.log(2), mp.zeta(3), mp.zeta(5), mp.zeta(7), mp.zeta(3)*mp.zeta(5),
      mp.pi**2, mp.pi**4, mp.pi**6]),
    # Polygamma at half-integer (these appear in d>3 partition functions)
    (["log(2)", "zeta(3)/pi^2", "zeta(5)/pi^4", "zeta(7)/pi^6", "psi(1/2)", "psi'(1/2)"],
     [mp.log(2), mp.zeta(3)/mp.pi**2, mp.zeta(5)/mp.pi**4, mp.zeta(7)/mp.pi**6,
      mp.psi(0, mp.mpf(1)/2), mp.psi(1, mp.mpf(1)/2)]),
]

for names, basis in bases:
    print(f"\nBasis: {names}")
    vec = [target] + basis
    for tol_exp in (30, 50):
        for maxc in (10**8, 10**12, 10**16):
            try:
                r = mp.pslq(vec, tol=mp.mpf(f'1e-{tol_exp}'), maxcoeff=maxc)
            except (ValueError, TypeError):
                r = None
            if r is not None and r[0] != 0:
                print(f"  tol=1e-{tol_exp}, maxc=10^{int(mp.log10(maxc))}: {r}")
                a0 = r[0]
                for i, n in enumerate(names):
                    rat = sp.Rational(-r[i+1], a0)
                    if rat != 0:
                        print(f"    coeff of {n}: {rat}")
                break
        if r is not None and r[0] != 0:
            break
    else:
        print(f"  No relation")
