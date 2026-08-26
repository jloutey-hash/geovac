"""Build Phase 1 (diagnostic gate): extract N(D)'s Stokes constants and test whether they are
ALGEBRAIC, and verify the D=0 reduction to the elliptic period.

Borel transform of N(D)'s asymptotic series is psi(zeta)=[(zeta+2)(rho(1+zeta)^2+1-rho)]^{-1/2}
(psi_k = Taylor coeffs).  psi has SQUARE-ROOT singularities at zeta*=-2 and -1 +- i*omega,
omega=sqrt((1-rho)/rho).  The Stokes constant at a sqrt singularity is the local amplitude
   a_* = lim_{zeta->zeta*} psi(zeta)*sqrt(zeta-zeta*)     (up to universal factors).
Prediction: a_* is ALGEBRAIC in rho (psi is 1/sqrt(quartic)).  We extract a_* two ways:
(A) closed-form limit;  (B) from the large-order growth of psi_k (Richardson), and check they agree.
Then verify N(0)*sqrt(c1)=L(0,rho) equals the complete elliptic integral K(1-rho).
"""
import mpmath as mp
mp.mp.dps = 40

def run(rho):
    rho = mp.mpf(rho); omega = mp.sqrt((1-rho)/rho)
    psi = lambda z: 1/mp.sqrt((z+2)*(rho*(1+z)**2 + (1-rho)))
    print(f"\n===== rho = {mp.nstr(rho,6)}  (omega = {mp.nstr(omega,6)}) =====")

    # ---- (A) closed-form Stokes amplitudes at the three non-dominant branch points
    aA = {}
    for name, zc in [("zeta=-2", mp.mpf(-2)),
                     ("zeta=-1+iw", -1 + 1j*omega),
                     ("zeta=-1-iw", -1 - 1j*omega)]:
        amp = mp.limit(lambda h: psi(zc + h)*mp.sqrt(h), 0)
        aA[name] = amp
        print(f"  (A) amplitude a_* at {name:>10}: {mp.nstr(amp,12)}")

    # ---- (B) leading amplitude from large-order psi_k (nearest singularity zeta=-2, dist 2)
    M = 60
    psic = mp.taylor(psi, 0, M)
    # (zeta+2)^{-1/2} Taylor coeff = 2^{-1/2} binom(-1/2,k) 2^{-k};  psi_k / that -> a_*(-2)
    def ref_k(k):
        return mp.mpf(2)**(mp.mpf(-1)/2) * mp.binomial(mp.mpf(-1)/2, k) * mp.mpf(2)**(-k)
    seq = [psic[k]/ref_k(k) for k in range(20, M+1)]     # -> a_*(-2) with O(1/k) tail
    # Richardson on the tail
    def richardson(seq):
        s = list(seq)
        for _ in range(6):
            s = [ (2*s[i+1]-s[i]) for i in range(len(s)-1) ] if False else \
                [ ((i+1)*s[i+1]-s[i]) for i in range(len(s)-1) ]  # generic 1/k Richardson
        return s[-1]
    aB = richardson(seq)
    print(f"  (B) large-order extraction of a_*(-2):    {mp.nstr(aB,12)}   (closed-form {mp.nstr(aA['zeta=-2'],12)})")
    print(f"      agree? {mp.nstr(abs(aB-aA['zeta=-2']),3)}")

    # ---- algebraicity check: are the amplitudes algebraic in rho?  test squares.
    print(f"  a_*(-2)^2               = {mp.nstr(aA['zeta=-2']**2,10)}   (predict 1)")
    print(f"  a_*(-1+iw)^2            = {mp.nstr(aA['zeta=-1+iw']**2,10)}  (algebraic in rho?)")
    #   analytic predicted: near -1+iw, (z+2)->1+iw, other factor ~ 2 i rho omega (z-(-1+iw))
    pred = 1/(2j*rho*omega*(1+1j*omega))
    print(f"  predicted a_*(-1+iw)^2  = {mp.nstr(pred,10)}   diff {mp.nstr(abs(aA['zeta=-1+iw']**2-pred),3)}")

    # ---- (C) D=0 reduction to the elliptic period
    L0 = mp.quad(lambda x: 1/mp.sqrt((x**2-1)*(rho*x**2 + (1-rho))), [1, mp.inf])
    # mpmath ellipk(m) is K with parameter m=k^2
    print(f"  L(0,rho)                = {mp.nstr(L0,14)}")
    for label, val in [("K(1-rho) [m=1-rho]", mp.ellipk(1-rho)),
                       ("K(rho)   [m=rho]",   mp.ellipk(rho)),
                       ("(1/sqrt?)*K...",     mp.ellipk(1-rho))]:
        print(f"     vs {label:>20} = {mp.nstr(val,14)}   ratio {mp.nstr(L0/val,10)}")

for rho in ["1/5", "1/3", "2/7"]:
    run(rho)
