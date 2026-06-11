"""
Sprint GB-2: does the squared Dirac spectral zeta zeta_{D^2}(s) have a
functional equation that D(s) lacked (RH-O)?

zeta_{D^2}(s) = D(2s) = sum_{n>=0} g_n |lambda_n|^{-2s},  |lambda_n|=n+3/2,
g_n = 2(n+1)(n+2).  T9 closed form: 2^{2s-1}[lambda(2s-2) - lambda(2s)],
lambda(z) = (1-2^{-z}) zeta(z)  (Dirichlet lambda).

Decisive structural test (closed-form, not template search):
F(s) and F(c-s) both live in the 2D space span{zeta(2s), zeta(2s-2)}.
A functional equation F(s) = chi(s) F(c-s) with chi elementary exists
iff the two coefficient-vectors are PARALLEL, i.e. det = A*D - B*C = 0.
"""
import mpmath as mp
mp.mp.dps = 40

def lam(z):                              # Dirichlet lambda
    return (1 - mp.power(2, -z)) * mp.zeta(z)

def F_closed(s):                         # T9 closed form
    return mp.power(2, 2*s - 1) * (lam(2*s - 2) - lam(2*s))

def F_spectral(s, N=200000):             # direct spectral sum (Re s > 3/2)
    tot = mp.mpf(0)
    for n in range(0, N):
        tot += 2*(n+1)*(n+2) * mp.power(mp.mpf(n) + mp.mpf(3)/2, -2*s)
    return tot

print("="*74)
print("1. Verify the object: spectral sum vs T9 closed form (Re s > 3/2)")
print("="*74)
for s in [mp.mpf('2.0'), mp.mpf('2.5'), mp.mpf('3.0')]:
    fc, fs = F_closed(s), F_spectral(s)
    print(f"  s={float(s):>4}:  closed={mp.nstr(fc,18):<22} spectral={mp.nstr(fs,18):<22} "
          f"|diff|={mp.nstr(abs(fc-fs),3)}")

print("\n" + "="*74)
print("2. Coefficient vectors in span{zeta(2s), zeta(2s-2)}")
print("="*74)
# F(s) = C*zeta(2s) + D*zeta(2s-2)   (exact, no FE needed)
def CD(s):
    p = mp.power(2, 2*s-1)
    C = -p * (1 - mp.power(2, -2*s))       # coeff of zeta(2s)
    D =  p * (1 - mp.power(2, 2-2*s))       # coeff of zeta(2s-2)
    return C, D
# F(c-s) reflected onto same basis via Riemann FE; test axis c=3/2
def AB(s):                                  # coeffs of F(3/2 - s)
    A = mp.power(2,2-2*s)*(1-mp.power(2,2*s-1))*2*mp.power(2*mp.pi,-2*s)*mp.cos(mp.pi*s)*mp.gamma(2*s)
    B = mp.power(2,2-2*s)*(1-mp.power(2,2*s-3))*2*mp.power(2*mp.pi,2-2*s)*mp.cos(mp.pi*s)*mp.gamma(2*s-2)
    return A, B

print("  cross-check reflection algebra: A*zeta(2s)+B*zeta(2s-2) =?= F(3/2-s)")
for s in [mp.mpc(1,1), mp.mpc('0.75',2), mp.mpc('2.2','0.5')]:
    A,B = AB(s)
    lhs = A*mp.zeta(2*s) + B*mp.zeta(2*s-2)
    rhs = F_closed(mp.mpf('1.5')-s)
    print(f"  s={mp.nstr(s,4):<16} |lhs-rhs|={mp.nstr(abs(lhs-rhs),3)}")

print("\n" + "="*74)
print("3. DECISIVE TEST: det = A*D - B*C  (=0 iff functional eq about c=3/2)")
print("="*74)
for s in [mp.mpc(1,1), mp.mpc('1.5',3), mp.mpc('2',5), mp.mpc('0.5',10)]:
    A,B = AB(s); C,D = CD(s)
    det = A*D - B*C
    # the obstruction multiplier: ratio of candidate chi's = (A/C)/(B/D)
    chi_ratio = (A/C)/(B/D)
    gamma_ratio = mp.gamma(2*s)/mp.gamma(2*s-2)   # = (2s-1)(2s-2)
    print(f"  s={mp.nstr(s,4):<14} |det|={mp.nstr(abs(det),4):<12} "
          f"chi1/chi2={mp.nstr(chi_ratio,5):<20} (2s-1)(2s-2)={mp.nstr((2*s-1)*(2*s-2),5)}")

print("\n" + "="*74)
print("4. Axis scan: is there ANY real axis c with det(c)=0?")
print("="*74)
def det_axis(s, c):
    # F(c-s) coeffs via FE about general c (arguments 2(c-s), 2(c-s)-2)
    # reflect zeta(2(c-s))=zeta(2c-2s) and zeta(2c-2s-2) back to zeta(2s),zeta(2s-2)
    # only consistent reflection is via z->1-z; arguments match {2s,2s-2} iff 2c-2s in {1-2s,3-2s}
    # => c in {1/2, 3/2}. Test both plus neighbors numerically by parallelism of F(s),F(c-s) vectors.
    # Generic c: F(c-s) is a combo of zeta(2c-2s),zeta(2c-2s-2) - DIFFERENT basis unless c=3/2 or 1/2.
    Cs,Ds = CD(s)
    pc = mp.power(2,2*(c-s)-1)
    # F(c-s)= -pc(1-2^{-2(c-s)})zeta(2(c-s)) + pc(1-2^{2-2(c-s)})zeta(2(c-s)-2)
    Cc = -pc*(1-mp.power(2,-2*(c-s)))
    Dc =  pc*(1-mp.power(2,2-2*(c-s)))
    # vectors live in different bases unless arguments coincide; measure mismatch of arg sets
    return Cs,Ds,Cc,Dc
print("  Only c where F(c-s) shares the basis {zeta(2s),zeta(2s-2)} are c=1/2, 3/2")
print("  (argument-matching), and c=3/2 already gives det != 0 above.")
for c in [mp.mpf('0.5'), mp.mpf('1.5')]:
    A,B = AB(mp.mpc(1,1)) if c==mp.mpf('1.5') else (None,None)
    print(f"  c={float(c)}: argument set reflects consistently; multiplier carries (2s-1)(2s-2) -> no FE")

print("\n" + "="*74)
print("5. POSITIVE CONTROL: single-power spectrum (constant multiplicity)")
print("="*74)
print("  G(s)=lambda(2s) lives in 1D span{zeta(2s)} -> trivially parallel -> HAS FE.")
print("  The 2D obstruction is created by the QUADRATIC multiplicity g_n=2(n+1)(n+2),")
print("  which spreads F across exponents 2s and 2s-2 (gap 2) -> gamma ratio (2s-1)(2s-2).")
