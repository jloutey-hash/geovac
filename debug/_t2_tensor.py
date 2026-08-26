"""FAST tensor T2: global k-grid + sin^2 outer maps => s,t are independent 1D node sets.
Precompute P(s_i,k_m), sin(k_m s_i), cos(k_m s_i) ONCE; inner double integral is pure
multiply-add (sin(k(s+t))=sin*cos+cos*sin). ~50x faster than per-point quadrature.
Validated against the 19-digit anchor.  (No Duffy corner: tests whether the sin^2 tensor
grid alone handles the (s+t)->0 rho^{3/2} corner; if not, add a corner patch.)"""
import mpmath as mp, sys, time
sys.path.insert(0,'debug')
from _fastgl import fast_gl

def Pcol(sval, kk):
    # P(s,k)=c e^{-D}(1/D^3+3/D^4+3/D^5), c=s(1-s), D=sqrt(c k^2+1); vector over kk for one s
    c=sval*(1-sval)
    out=[]
    for k in kk:
        D=mp.sqrt(c*k*k+1)
        out.append(c*mp.e**(-D)*(1/D**3+3/D**4+3/D**5))
    return out

def T2_tensor(Nc, Nk, Kmax):
    H=mp.pi/2
    # global k-grid on decay map k=Kmax u/(1-u)
    xu,wu=fast_gl(Nk)
    kk=[]; wk=[]
    for x,w in zip(xu,wu):
        u=(x+1)/2; k=Kmax*u/(1-u); kk.append(k); wk.append((w/2)*Kmax/(1-u)**2)
    # outer sin^2 nodes on [0,1]
    xp,wp=fast_gl(Nc)
    sv=[]; sj=[]
    for x,w in zip(xp,wp):
        phi=H*(x+1)/2; sv.append(mp.sin(phi)**2); sj.append(mp.sin(2*phi)*H*w/2)
    # precompute per-node transcendentals
    P=[Pcol(s,kk) for s in sv]
    S=[[mp.sin(k*s) for k in kk] for s in sv]
    C=[[mp.cos(k*s) for k in kk] for s in sv]
    invk=[1/k for k in kk]
    # double sum (inner loop: pure mult-add)
    tot=mp.mpf(0)
    for i in range(Nc):
        Pi=P[i]; Si=S[i]; Ci=C[i]; si=sv[i]; wi=sj[i]
        for j in range(Nc):
            Pj=P[j]; Sj=S[j]; Cj=C[j]; b=si+sv[j]
            acc=mp.mpf(0)
            for m in range(Nk):
                # j0(k b)=sin(kb)/(kb); sin(kb)=Si*Cj+Ci*Sj
                sinkb=Si[m]*Cj[m]+Ci[m]*Sj[m]
                acc+=wk[m]*invk[m]*sinkb*Pi[m]*Pj[m]
            tot+=wi*sj[j]*(acc/b)
    return (8/mp.pi)*tot

if __name__=='__main__':
    mp.mp.dps=40; ref=mp.mpf('0.3953557659017139641')
    for Nc,Nk,Kmax in [(60,400,40),(80,600,50),(100,800,60)]:
        t0=time.time(); v=T2_tensor(Nc,Nk,Kmax); dt=time.time()-t0
        print(f'Nc={Nc} Nk={Nk} Kmax={Kmax}: {mp.nstr(v,30)}  ({dt:.0f}s)  |anch|={mp.nstr(abs(v-ref),3)}',flush=True)
    print('DONE',flush=True)
