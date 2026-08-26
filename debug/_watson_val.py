import sys,time; sys.path.insert(0,'debug')
import mpmath as mp
import _watson_fibre as W
mp.mp.dps=40
cs=mp.mpf('0.2'); b=mp.mpf('0.5'); Nmax=16
Fn=W.F_taylor(Nmax); ms=[W.m_n(cs,b,n) for n in range(Nmax+1)]
def J_exact(ct):
    f=lambda k:(mp.sin(k*b)/(k*b) if k>1e-30 else mp.mpf(1))*W.Pcs(cs,k)*W.Pcs(ct,k)
    return mp.quadosc(f,[0,mp.inf],period=2*mp.pi/b)
print("c_t     J_watson(opt-trunc)         trunc-n   |Jw-Jexact|",flush=True)
for ct in ['0.05','0.02','0.01','0.005','0.002']:
    ct=mp.mpf(ct)
    Jw,bn=W.J_watson(cs,ct,b,Fn,ms)
    Je=J_exact(ct)
    print(f"  {mp.nstr(ct,4):>6}: {mp.nstr(Jw,24)}  n={bn:2d}   {mp.nstr(abs(Jw-Je),3)}",flush=True)
