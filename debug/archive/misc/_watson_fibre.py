"""The Watson small-c_t fibre  J(c_s,c_t,b) = sum_n F_n c_t^{n+1} m_n(c_s,b).
F(eps)=e^{-sqrt(1+eps)}[(1+eps)^{-3/2}+3(1+eps)^{-2}+3(1+eps)^{-5/2}], eps=c_t k^2, F_n=Taylor coeffs;
m_n(c_s,b)=int_0^inf k^{2n} j0(kb) P(c_s,k) dk (one-mass moments).  Functions only."""
import mpmath as mp

def F_taylor(N):
    f=lambda e: mp.e**(-mp.sqrt(1+e))*((1+e)**(mp.mpf(-3)/2)+3*(1+e)**-2+3*(1+e)**(mp.mpf(-5)/2))
    return mp.taylor(f,0,N)

def Pcs(cs,k):
    d=mp.sqrt(cs*k*k+1); return cs*mp.e**(-d)*(d**-3+3*d**-4+3*d**-5)

def m_n(cs,b,n):
    f=lambda k:(k**(2*n))*(mp.sin(k*b)/(k*b) if k>1e-30 else mp.mpf(1))*Pcs(cs,k)
    return mp.quadosc(f,[0,mp.inf],period=2*mp.pi/b)

def J_watson(cs,ct,b,Fn,ms,N=None):
    if N is None: N=len(Fn)-1
    terms=[Fn[n]*ct**(n+1)*ms[n] for n in range(N+1)]
    S=mp.mpf(0); best=None; bn=0
    for n,t in enumerate(terms):
        S+=t
        if n>0 and abs(t)>abs(terms[n-1]):
            break
        best=S; bn=n
    return best,bn
