"""Target B step 3: Borel-Laplace RESUM the fibre c_t-series and test vs exact (beat optimal trunc).
a_n ~ (2n)! => J/c_t = sum a_n c_t^n = int_0^inf e^{-t} B(c_t t^2) dt, B(z)=sum a_n z^n/(2n)! (Pade'd).
Singularities of B are COMPLEX (0.054+-0.48i), none on z>=0 => Laplace contour clean => unambiguous."""
import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import _watson_fibre as W
mp.mp.dps=50
cs=mp.mpf('0.2'); b=mp.mpf('0.5'); N=24
Fn=W.F_taylor(N); ms=[W.m_n(cs,b,n) for n in range(N+1)]
a=[Fn[n]*ms[n] for n in range(N+1)]
bc=[a[n]/mp.factorial(2*n) for n in range(N+1)]
def J_exact(ct):
    f=lambda k:(mp.sin(k*b)/(k*b) if k>1e-30 else mp.mpf(1))*W.Pcs(cs,k)*W.Pcs(ct,k)
    return mp.quadosc(f,[0,mp.inf],period=2*mp.pi/b)
def J_borel(ct,L,M):
    p,q=mp.pade(bc,L,M)
    def B(z):
        zz=mp.mpf(z) if not isinstance(z,mp.mpc) else z
        num=sum(p[i]*zz**i for i in range(len(p))); den=sum(q[i]*zz**i for i in range(len(q)))
        return num/den
    # J/c_t = int_0^inf e^{-t} B(c_t t^2) dt
    val=mp.quad(lambda t: mp.e**(-t)*B(ct*t*t),[0,mp.inf])
    return ct*val
# also optimal-truncation for comparison
def J_opt(ct):
    S=mp.mpf(0); prev=None; best=None
    for n in range(N+1):
        term=a[n]*ct**(n+1); S+=term
        if prev is not None and abs(term)>abs(prev): break
        best=S; prev=term
    return best
print("c_t     |J_opt-Jexact|     |J_borel[12/12]-Jexact|   |J_borel[11/13]-Jexact|",flush=True)
for ct in ['0.05','0.02','0.01']:
    ct=mp.mpf(ct); Je=J_exact(ct)
    jo=J_opt(ct)
    try: jb1=J_borel(ct,12,12)
    except Exception as e: jb1=None
    try: jb2=J_borel(ct,11,13)
    except Exception as e: jb2=None
    s1=mp.nstr(abs(jb1-Je),3) if jb1 is not None else 'fail'
    s2=mp.nstr(abs(jb2-Je),3) if jb2 is not None else 'fail'
    print(f"  {mp.nstr(ct,4):>5}: {mp.nstr(abs(jo-Je),3):>12}      {s1:>16}   {s2:>16}",flush=True)
