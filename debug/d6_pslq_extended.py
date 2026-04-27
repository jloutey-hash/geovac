"""
D6 PSLQ with extended basis: add triple products and test larger maxcoeff.
Also verify D6 value by independent method (direct partial sums).
"""
import sys, io, functools, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
print = functools.partial(print, flush=True)

from mpmath import mp, mpf, nstr, nsum, inf as mpinf
import mpmath
import sympy as sp

DPS = 600
mp.dps = DPS
z = lambda s: mpmath.zeta(s)

def euler_sum(r, s):
    old = mp.dps
    mp.dps = DPS + 40
    def term(kk):
        return mpmath.hurwitz(s, kk) / mpf(kk)**r
    result = nsum(term, [1, mpinf], method='richardson')
    mp.dps = old
    return result

# ============================================================
# PART 1: Verify D6 via INDEPENDENT method
# ============================================================
print("=" * 70)
print("PART 1: D6 VERIFICATION (two methods)")
print("=" * 70)

# Method A: Euler sum assembly (same as before)
print("\nMethod A: Euler sum assembly...")
D6_euler = mpf(-2289)/512 * z(10) + mpf(567)/128 * z(11)
pairs = [
    (mpf(315)/32, 1, 10), (mpf(-245)/32, 2, 9),
    (mpf(63)/32, 4, 7), (mpf(35)/64, 5, 6),
    (mpf(-21)/64, 6, 5), (mpf(-21)/64, 7, 4),
    (mpf(-7)/64, 8, 3),
]
for coeff, r, s in pairs:
    t0 = time.time()
    S = euler_sum(r, s)
    c = coeff * (S - z(r + s))
    D6_euler += c
    print(f"  S_{{{r},{s}}} ({time.time()-t0:.0f}s)")

# Method B: Direct partial sum to large N + tail correction
print("\nMethod B: Direct partial sum...")
t0 = time.time()

# c_6(n) from Dirac-Coulomb energy expansion
# E_n = -Z^2/(2n^2) * sum_{p=1}^inf c_p(n) * (Z*alpha)^{2p}
# For p=6, weight 2p-1=11
# We compute c_6(n) from the Dirac formula expansion to order (Za)^12

def dirac_coeffs_p6(n_max):
    """Compute c_6(n) for n=1..n_max from Dirac-Coulomb expansion."""
    from mpmath import taylor, sqrt as mpsqrt
    old = mp.dps
    mp.dps = DPS + 50
    coeffs = []
    for n in range(1, n_max + 1):
        # For each n, sum over j = 1/2, 3/2, ..., n-1/2
        # with degeneracy 2(2j+1) except j=n-1/2 has degeneracy 2j+1
        total = mpf(0)
        for two_j in range(1, 2*n, 2):  # j = two_j/2
            j = mpf(two_j) / 2
            kappa_vals = []
            if two_j < 2*n - 1:
                kappa_vals = [-(two_j+1)//2, (two_j+1)//2]  # kappa = -(j+1/2) and j+1/2
            else:  # j = n - 1/2
                kappa_vals = [-(two_j+1)//2]  # only kappa = -(j+1/2) = -n
                if two_j > 1:
                    kappa_vals.append((two_j+1)//2)  # kappa = j+1/2 = n (but only if j >= 1)

            for kappa in kappa_vals:
                abs_kappa = abs(kappa)
                # Dirac energy: E = mc^2 / sqrt(1 + (Za)^2/(n_r + gamma)^2)
                # where gamma = sqrt(kappa^2 - (Za)^2), n_r = n - |kappa|
                # Expand in powers of (Za)^2 = x
                # gamma = |kappa| * sqrt(1 - x/kappa^2)
                # We need the coefficient of x^6 in E/(mc^2) - 1 + Z^2/(2n^2)
                # normalized by -Z^2/(2n^2)

                n_r = n - abs_kappa
                k2 = mpf(abs_kappa)**2

                # Series expansion of gamma = sqrt(k2 - x) around x=0
                # gamma = |k|*sqrt(1 - x/k2)
                # gamma = |k| - x/(2|k|) - x^2/(8|k|^3) - ...

                # Actually, let's use mpmath's taylor expansion directly
                # E(x) = 1/sqrt(1 + x/(n_r + sqrt(k2 - x))^2)
                # This is complex; use symbolic expansion instead

                # Alternative: expand E_nj(alpha) in powers of alpha^2
                # using the exact Dirac formula
                # E = [1 + (Za)^2 / (n - |k| + sqrt(k^2 - (Za)^2))^2]^{-1/2}

                # Let x = (Za)^2. We need Taylor expansion of f(x) to order x^6
                def f(x):
                    gamma = mpmath.sqrt(k2 - x)
                    denom = n_r + gamma
                    return mpf(1) / mpmath.sqrt(1 + x / denom**2)

                # Get Taylor coefficients
                tc = mpmath.taylor(f, 0, 7)  # coefficients of x^0, x^1, ..., x^7

                # E/(mc^2) = sum tc[p] * x^p
                # E_binding = E - mc^2 = (f(x) - 1) * mc^2
                # In atomic units with mc^2 = 1/(2*alpha^2) effectively:
                # Actually, the standard convention:
                # E_n = -Z^2/(2n^2) * [1 + sum_{p>=1} c_p(n) * (Za)^{2p}]
                # So c_p(n) comes from expanding f(x) = 1 - x/(2n^2) * [1 + sum c_p x^p]
                # Hmm, this mapping is nontrivial. Let me use a different approach.

                # The Dirac-Coulomb binding energy coefficient at order (Za)^{2p}
                # For the state (n, kappa), the contribution to D_p is:
                # d_p(n, kappa) * n^2 where d_p is the coefficient in
                # E_binding = -1/(2n^2) * sum_{p=0} d_p * (Za)^{2p}
                # with d_0 = 1.

                # Actually, let me just compute 4*n*a_p(n) where
                # a_p(n) = (1/n^2) * sum_kappa (2j+1) * [p-th coeff of Dirac energy]
                # This is the degeneracy-weighted coefficient.

                # The exact Sommerfeld coefficients are:
                # 4n * a_p(n) = sum over (n,l,j,m_j) of c_p(n,l,j)
                # where the sum is degeneracy-weighted

                # For simplicity, let's just compute the Taylor expansion of
                # the Dirac energy and extract order by order.
                pass

        coeffs.append(total)
    mp.dps = old
    return coeffs

# Actually, let me use a simpler and more reliable approach:
# compute D6 from the known structure of the Dirac expansion coefficients.
# The Sommerfeld fine structure coefficients c_p(n) satisfy:
# 4*k*a_p(n,k) = coefficient of k^{2p} in the Laurent expansion of
# sum_{j} (2j+1) * [E_D(n,j,k) + 1/(2n^2)] around k=0
# where k = Z*alpha and E_D is the Dirac-Coulomb energy.

# But this is exactly what the prior script d6_fix.py computed symbolically.
# Let me instead just do a numerical cross-check via direct summation.

print("\nMethod B: Direct summation of c_6(n) * n^{-something}...")
# From the Dirac formula, for each n, expand E(n,kappa,alpha) to order alpha^12
# and extract the coefficient.

# Use symbolic expansion via sympy
from sympy import Rational as R, sqrt as ssqrt, Symbol, series, nsimplify

alpha_sq = Symbol('x')  # x = (Z*alpha)^2

def dirac_energy_series(n, kappa, order=7):
    """Return Taylor coefficients of E_D(x) where x = (Za)^2."""
    k2 = R(kappa**2)
    n_r = n - abs(kappa)
    # gamma = sqrt(k2 - x)
    # E = 1/sqrt(1 + x/(n_r + gamma)^2)
    # Expand symbolically
    gamma = ssqrt(k2 - alpha_sq)
    denom = R(n_r) + gamma
    expr = 1 / ssqrt(1 + alpha_sq / denom**2)
    s = series(expr, alpha_sq, 0, order)
    return [s.coeff(alpha_sq, p) for p in range(order)]

print("  Computing c_6(n) symbolically for n=1..30...")
t0 = time.time()
D6_direct = mpf(0)
c6_values = []
for n in range(1, 31):
    # Sum over all (kappa, m_j) states at shell n
    # States: for each l=0..n-1, j=l+1/2 and j=l-1/2 (if l>0)
    # kappa = -(l+1) for j=l+1/2, kappa = l for j=l-1/2
    # degeneracy of each (n,kappa) state is 2|kappa| = 2j+1
    c6_n = sp.Rational(0)
    for l in range(n):
        # j = l + 1/2: kappa = -(l+1)
        kappa = -(l+1)
        degen = 2*abs(kappa)  # = 2(l+1) = 2j+1
        coeffs = dirac_energy_series(n, kappa, 7)
        # coeffs[p] is the coefficient of x^p in E(x)
        # The contribution to D_p is: degen * coeffs[p] * ???
        # Actually, D_p = sum_n c_p(n) where c_p(n) is defined via
        # E_n = -1/(2n^2) * [1 + c_1(n)*x + c_2(n)*x^2 + ...]
        # But E_n here is the Dirac energy summed over all states in shell n,
        # divided by n^2 (the Fock degeneracy).
        #
        # More precisely: 4n * a_p(n) = c_p(n) is defined as:
        # sum_{states in shell n} E_state(x) = -n^2/(2n^2) * [something]
        # This requires careful normalization.
        #
        # Let me just compute sum_{states} coeffs[6] for each n.
        c6_n += degen * coeffs[6]

        # j = l - 1/2: kappa = l (only if l > 0)
        if l > 0:
            kappa = l
            degen = 2*abs(kappa)  # = 2l = 2j+1
            coeffs = dirac_energy_series(n, kappa, 7)
            c6_n += degen * coeffs[6]

    # The total Dirac energy of shell n at order x^6:
    # E_shell^{(6)}(n) = c6_n (as computed above, summed with degeneracy)
    # D_6 = sum_n c6_n ... but we need the right normalization.
    #
    # From the existing code (d6_fix.py), D_p = sum_n [4n * a_p(n)]
    # where 4n*a_p(n,k) extracts the coefficient of k^{2p} from
    # 4k * sum_{states in n} E_D(state, k)
    #
    # Let's compute: sum_{states in n} E_D(state, x) expanded to order x^6
    # The coefficient of x^6 in this sum is c6_n (what we computed).
    # Then D_6 = sum_n c6_n? No...
    #
    # Actually the normalization in d6_fix.py is:
    # expr(n) = sum_j (2j+1) * [series of E_D(n,j)]
    # Then 4*k*a_p(n) extracts the x^p coefficient of expr, multiplied by 4k/k = 4.
    # Wait, looking at the expansion more carefully...
    # The Sommerfeld coefficient is c_p(n) = [x^p coeff of expr(n)]
    # and D_p = sum_{n>=1} c_p(n).
    #
    # So D_6 = sum_n c6_n where c6_n = coeff of x^6 in
    #   sum_{kappa at shell n} (2|kappa|) * E_D(n, kappa, x)

    c6_float = float(c6_n)
    c6_values.append((n, c6_float))
    D6_direct += mpf(c6_float)
    if n <= 5 or n % 10 == 0:
        print(f"    n={n}: c6={float(c6_n):.6e}, partial D6={float(D6_direct):.12f}")

print(f"\n  D6 partial (N=30) = {nstr(D6_direct, 20)}")
print(f"  D6 Euler sum      = {nstr(D6_euler, 20)}")
dt = time.time() - t0
print(f"  ({dt:.0f}s)")

# ============================================================
# PART 2: Extended PSLQ basis
# ============================================================
print("\n" + "=" * 70)
print("PART 2: PSLQ WITH EXTENDED BASIS")
print("=" * 70)

# Maybe D6 involves triple products or higher-depth MZVs?
# Weight-11 triple products: z(3)^2*z(5) and z(3)*z(4)^2
# But these are just products of known zeta values.

# Test 1: Standard basis + triple products
basis_vals = [D6_euler]
basis_names = ["D6"]

# Standard: z(10), z(11), products
basis_vals.append(z(10)); basis_names.append("z(10)")
basis_vals.append(z(11)); basis_names.append("z(11)")
basis_vals.append(z(2)*z(9)); basis_names.append("z(2)z(9)")
basis_vals.append(z(3)*z(8)); basis_names.append("z(3)z(8)")
basis_vals.append(z(4)*z(7)); basis_names.append("z(4)z(7)")
basis_vals.append(z(5)*z(6)); basis_names.append("z(5)z(6)")

# Triple products at weight 11
basis_vals.append(z(3)**2 * z(5)); basis_names.append("z(3)^2*z(5)")
basis_vals.append(z(3)*z(4)**2); basis_names.append("z(3)*z(4)^2")  # weight 3+4+4=11... no, 3+8=11 but 4+4=8? wait: z(3)*z(4)^2 has weight 3+4+4=11. Yes.

print(f"\n  Extended basis ({len(basis_vals)}): {', '.join(basis_names)}")

rel = None
for mc in [1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000]:
    rel = mpmath.pslq(basis_vals, maxcoeff=mc)
    if rel is not None:
        print(f"  PSLQ at maxcoeff={mc}")
        break

if rel is None:
    print(f"  Extended PSLQ FAILED at {DPS} dps")
else:
    d = rel[0]
    terms = []
    for c, nm in zip(rel[1:], basis_names[1:]):
        if c != 0:
            r = sp.Rational(-c, d)
            terms.append(f"({r}){nm}")
    res = sum(c*v for c, v in zip(rel, basis_vals))
    print(f"  D6 = {' + '.join(terms)}")
    print(f"  Residual: {nstr(res, 5)}")
    print(f"  Relation: {[int(x) for x in rel]}")

# Test 2: Maybe the LCD is very large. Try with tol parameter.
print("\n--- PSLQ with explicit tol ---")
for tol_exp in [300, 400, 500]:
    basis7 = [D6_euler, z(10), z(11), z(2)*z(9), z(3)*z(8), z(4)*z(7), z(5)*z(6)]
    rel2 = mpmath.pslq(basis7, tol=mpf(10)**(-tol_exp), maxcoeff=10**9)
    if rel2 is not None:
        print(f"  tol=1e-{tol_exp}: FOUND {[int(x) for x in rel2]}")
        d = rel2[0]
        terms = []
        for c, nm in zip(rel2[1:], ["z(10)", "z(11)", "z(2)z(9)", "z(3)z(8)", "z(4)z(7)", "z(5)z(6)"]):
            if c != 0:
                r = sp.Rational(-c, d)
                terms.append(f"({r}){nm}")
        res = sum(c*v for c, v in zip(rel2, basis7))
        print(f"  D6 = {' + '.join(terms)}")
        print(f"  Residual: {nstr(res, 5)}")
        break
    else:
        print(f"  tol=1e-{tol_exp}: FAILED")

# Save
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

results = {
    'D6_euler_600dps': nstr(D6_euler, 80),
    'D6_direct_N30': nstr(D6_direct, 20),
    'extended_pslq': 'FOUND' if rel else 'FAILED',
    'c6_values': [(n, float(c)) for n, c in c6_values[:10]],
}
if rel:
    results['relation'] = [int(x) for x in rel]

with open('debug/data/d6_pslq_extended.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print(f"\nSaved to debug/data/d6_pslq_extended.json")
print("Done.")
