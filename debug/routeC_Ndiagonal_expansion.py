"""Diagnostic: expand N(D) (Paper 59 eq:laplace) around the Bessel point rho=0 in the
scale-mismatch rho=c2/c1; measure coefficient structure + convergence (do not assert).

L(D,rho) = int_0^inf e^{-D cosh t}/sqrt(rho cosh^2 t + 1 - rho) dt   (L = sqrt(c1) N; = K0(D) at rho=0)
Term-by-term: L = sum_n a_n rho^n K_n(D)/D^n,  a_n = (-1)^n ((2n-1)!!)^2 / (2^n n!).
"""
import mpmath as mp
mp.mp.dps = 30

def a_coeff(n):
    dfac = mp.mpf(1)
    for k in range(1, 2*n, 2):
        dfac *= k
    return (-1)**n * dfac**2 / (mp.mpf(2)**n * mp.factorial(n))

def L_direct(D, rho):
    # e^{-D cosh t} kills the tail; cosh(6)~202, e^{-D*202} negligible for D>=1.5
    f = lambda t: mp.e**(-D*mp.cosh(t)) / mp.sqrt(rho*mp.cosh(t)**2 + 1 - rho)
    return mp.quad(f, [0, 6])

def report(D, rho, Nmax=22):
    exact = L_direct(D, rho)
    print(f"\n=== D={D}, rho={rho}   exact L = {mp.nstr(exact,18)} ===")
    S = mp.mpf(0); best_err = None; best_N = None; errs = []
    for n in range(0, Nmax+1):
        term = a_coeff(n) * rho**n * mp.besselk(n, D) / D**n
        S += term
        err = abs(S - exact); errs.append(err)
        if best_err is None or err < best_err:
            best_err, best_N = err, n
        if n <= 10 or n % 3 == 0:
            print(f"  n={n:>2}  term={mp.nstr(term,7):>15}  |S_N-exact|={mp.nstr(err,4)}")
    grew = errs[-1] > 30*best_err
    print(f"  --> optimal truncation N={best_N}, min err {mp.nstr(best_err,4)};  "
          f"err at N={Nmax} = {mp.nstr(errs[-1],4)}  =>  "
          f"{'ASYMPTOTIC (diverges past optimal N)' if grew else 'convergent'}")

for D, rho in [(2.0, 0.3), (3.0, 0.2), (2.0, 0.1)]:
    report(D, rho)
print("\n[sanity] L(2,0)-K0(2) =", mp.nstr(L_direct(2.0,0.0)-mp.besselk(0,2.0),4))
print("[coeffs] a_0..a_5 =", [mp.nstr(a_coeff(n),6) for n in range(6)])
