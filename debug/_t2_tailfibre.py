"""Tail-analytic fibre: J(s,t)=int_0^K j0(kb)PP dk (bounded quadrature) + analytic tail.
Tail via incomplete gammas: int_K^inf sin(kb)e^{-Ak}/k^n dk = Im[z^{n-1} Gamma(1-n, zK)], z=A-ib.
g(k)=PP e^{Ak} k^6 is stable (e^{Ak} cancels e^{-ak}); its 1/k-series coeffs d_m weight the tail."""
import mpmath as mp, sys
sys.path.insert(0,'debug'); import routeC_T2_highprec as H
from _fastgl import fast_gl

def A_of(s,t):
    return mp.sqrt(s*(1-s))+mp.sqrt(t*(1-t))

def gfun(s,t,k):
    A=A_of(s,t)
    return H.P(s,k)*H.P(t,k)*mp.e**(A*k)*k**6   # -> 1/(as*at), stable (no cancellation)

def dcoeffs(s,t,M,k0=None):
    # STABLE: fit g(k)=sum_{m=0..M} d_m/k^m at k=k0,2k0,...,(M+1)k0 (Vandermonde in 1/k).
    # adaptive k0 >= ~35/sqrt(c_min) keeps the fit inside the k>1/sqrt(c) asymptotic regime.
    if k0 is None:
        cmin=min(s*(1-s),t*(1-t)); k0=max(mp.mpf(60), 35/mp.sqrt(cmin))
    ks=[mp.mpf(k0)*(i+1) for i in range(M+1)]
    gs=[gfun(s,t,k) for k in ks]
    V=mp.matrix(M+1,M+1)
    for i,k in enumerate(ks):
        for m in range(M+1): V[i,m]=1/k**m
    d=mp.lu_solve(V, mp.matrix(gs))
    return [d[m] for m in range(M+1)]

def tail(s,t,b,K,M):
    A=A_of(s,t); z=A-1j*b
    cmin=min(s*(1-s),t*(1-t))
    # deep-corner / tiny-scale guard: tail is negligible AND numerically unstable there (z->0)
    if cmin < mp.mpf('1e-8') or abs(z)*K < mp.mpf('0.5'):
        return mp.mpf(0)
    try:
        d=dcoeffs(s,t,M)
        x=z*K; ex=mp.e**(-x)
        Gs=[]; G=mp.gammainc(-6, x)
        for m in range(M+1):
            Gs.append(G); a=-6-m
            G=(G - x**(a-1)*ex)/(a-1)
        tot=mp.mpc(0)
        for m in range(M+1):
            tot+= d[m]*z**(6+m)*Gs[m]
        val=(tot/b).imag
        return val if abs(val) < 1 else mp.mpf(0)   # clamp spurious blowups (|J|<~0.5 physically)
    except Exception:
        return mp.mpf(0)

def bounded(s,t,b,K,Nq):
    # int_0^K j0(kb) P(s,k)P(t,k) dk via GL on [0,K]
    xs,ws=fast_gl(Nq); tot=mp.mpf(0)
    for x,w in zip(xs,ws):
        k=K*(x+1)/2; wk=K*w/2
        j0=mp.sin(k*b)/(k*b) if k*b>mp.mpf('1e-40') else mp.mpf(1)
        tot+=wk*j0*H.P(s,k)*H.P(t,k)
    return tot

def Jtail(s,t,K=None,Nq=None,M=6):
    b=s+t
    cs=s*(1-s); ct=t*(1-t)
    if cs==0 or ct==0: return mp.mpf(0)
    cmin=min(cs,ct)
    if cmin < mp.mpf('1e-9'): return mp.mpf(0)   # deep tip: J~0
    # ADAPTIVE K: must exceed 1/sqrt(cmin) so the tail asymptotic series is valid
    if K is None: K=max(mp.mpf(95), 3.2/mp.sqrt(cmin))
    # Nq resolves the K*b oscillations (small at the corner since b->0 there)
    if Nq is None:
        want=max(300, 2.2*float(K*b)+200)
        for g in (300,400,500,650,850,1100,1500,2000):
            if g>=want: Nq=g; break
        else: Nq=2000
    return bounded(s,t,b,K,Nq)+tail(s,t,b,K,M)

if __name__=='__main__':
    mp.mp.dps=40
    print('Validate tail-analytic fibre vs direct decay-map fibre (H.J, high Nk):',flush=True)
    for (s,t,tag) in [(mp.mpf('0.3'),mp.mpf('0.17'),'interior'),
                      (mp.mpf('0.05'),mp.mpf('0.28'),'edge s small'),
                      (mp.mpf('0.5'),mp.mpf('0.5'),'diagonal')]:
        ref=H.J(s,t,400)
        for (K,Nq,M) in [(30,400,12),(40,500,16)]:
            v=Jtail(s,t,K,Nq,M)
            print(f'  {tag} (s={mp.nstr(s,3)},t={mp.nstr(t,3)}) K={K} Nq={Nq} M={M}: |Jtail-ref|={mp.nstr(abs(v-ref),3)}',flush=True)
    print('DONE',flush=True)
