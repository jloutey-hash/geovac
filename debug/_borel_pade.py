"""Target B step 2: Borel-Pade the fibre c_t-series a_n=F_n m_n (n=0..24) to LOCATE the Borel
singularities.  a_n ~ (2n)!-type => Borel transform B(z)=sum a_n z^n/(2n)! is convergent; its poles
(Pade) approximate the singularities.  Predicted (curve branch points, c_s=0.2): a real one near
z~c_s and/or a complex-conjugate pair (the erratic sign-oscillating ratios)."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import _watson_fibre as W
mp.mp.dps=60
cs=mp.mpf('0.2'); b=mp.mpf('0.5'); N=24
Fn=W.F_taylor(N); ms=[W.m_n(cs,b,n) for n in range(N+1)]
a=[Fn[n]*ms[n] for n in range(N+1)]
# Borel transform coeffs b_n = a_n/(2n)!
bc=[a[n]/mp.factorial(2*n) for n in range(N+1)]
print("Borel coeffs b_n=a_n/(2n)! (first few):",[mp.nstr(x,6) for x in bc[:6]],flush=True)
# Pade [L/M]
for L,M in [(12,12),(11,11),(10,10),(13,11)]:
    try:
        p,q=mp.pade(bc, L, M)
        roots=mp.polyroots(q[::-1], maxsteps=200, extraprec=200)
        roots=sorted(roots, key=lambda z: abs(z))
        print(f"\nPade[{L}/{M}] nearest poles of B(z) (=Borel singularities z*):",flush=True)
        for z in roots[:5]:
            print(f"   z*={mp.nstr(z,10)}   |z*|={mp.nstr(abs(z),8)}",flush=True)
    except Exception as e:
        print(f"Pade[{L}/{M}] failed: {e}",flush=True)
# predictions from the curve y^2=(c_s k^2+1)(c_t k^2+1): branch k=+-i/sqrt(c_s) => c_t k^2=-c_t/c_s,
# F branch at c_t k^2=-1 => c_t=c_s (diagonal). Report c_s for comparison.
print(f"\npredicted diagonal-degeneration scale c_t=c_s={mp.nstr(cs,6)} (real);",flush=True)
print("complex pair expected from the j0(kb)/other-branch interference (erratic ratios).",flush=True)
