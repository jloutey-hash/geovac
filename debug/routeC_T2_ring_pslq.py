"""General CM/Eisenstein ring PSLQ for T2 (and its twist-free sibling).

Generators carry (name, value, weight, signed?).  'signed' generators (periods varpi,P8)
admit NEGATIVE powers => quasiperiod (1/period) directions (Legendre E=pi/4varpi+varpi/2).
Weight-graded monomials up to total weight wmax; guarded by decoy + cross-precision.

Principled ring for a length-2, weight<=3 iterated-Eisenstein X(2) integral:
  pi(w1), varpi=K(1/2)(w1,signed), G=Catalan=beta(2)(w2), ln2(w1), P8 disc-8 period(w1,signed),
  zeta3(w3).  Justified: int_0^1 K = 2G, int_0^1 E = G+1/2, E(1/2)=pi/4varpi+varpi/2.
"""
from __future__ import annotations
import sys, itertools
import mpmath as mp
sys.path.insert(0,'debug')
from routeC_T2_pslq_decisive import guarded_pslq

def gens(dps, names):
    mp.mp.dps=dps+30
    pi=mp.pi
    varpi=mp.gamma(mp.mpf(1)/4)**2/(4*mp.sqrt(pi))
    P8=(mp.sqrt(1+mp.sqrt(2))*mp.gamma(mp.mpf(1)/8)*mp.gamma(mp.mpf(3)/8)
        /(mp.mpf(2)**(mp.mpf(13)/4)*mp.sqrt(pi)))
    table={'pi':(pi,1,False),'vp':(varpi,1,True),'G':(mp.catalan,2,False),
           'ln2':(mp.log(2),1,False),'P8':(P8,1,True),'z3':(mp.zeta(3),3,False)}
    out=[(n,)+table[n] for n in names]
    mp.mp.dps=dps
    return out  # list of (name,value,weight,signed)

def graded_ring(dps, wmax, names):
    G=gens(dps,names); mp.mp.dps=dps
    ring={'1':mp.mpf(1)}
    # exponent ranges: signed gens in [-wmax,wmax], others in [0, wmax//weight]
    ranges=[]
    for (n,v,w,signed) in G:
        hi=wmax  # cap by weight below
        ranges.append(range(-wmax,wmax+1) if signed else range(0,wmax+1))
    for exps in itertools.product(*ranges):
        wt=sum(abs(e)*G[i][2] for i,e in enumerate(exps))
        if 1<=wt<=wmax:
            key='*'.join(f'{G[i][0]}^{e}' for i,e in enumerate(exps) if e)
            val=mp.mpf(1)
            for i,e in enumerate(exps): val*=G[i][1]**e
            ring[key]=val
    return ring

def run(W_target_str, dps, wmax, names, maxcoeff=10**9, is_W=True):
    """W_target = the natural period (V*pi/8) if is_W, else raw target."""
    mp.mp.dps=dps+20
    W=mp.mpf(W_target_str)
    decoy=W*(1+mp.mpf(10)**(-11))+mp.euler/mp.mpf(10)**5
    ring=graded_ring(dps,wmax,names); mp.mp.dps=dps
    n=len(ring)
    print(f"\n=== ring={{{','.join(names)}}} wt<=({wmax})  dim={n}  dps={dps}  fp10^{mp.nstr(mp.mpf(dps)/max(1,n-1),3)} ===")
    real=guarded_pslq(W,ring,dps,maxcoeff,"REAL")
    dec =guarded_pslq(decoy,ring,dps,maxcoeff,"DECOY")
    def h(r): return None if r is None else max(abs(x) for x in r)
    hr,hd=h(real),h(dec)
    if hr is not None and hr<=80 and (hd is None or hd>8*hr):
        print(f"    >>> CANDIDATE (REAL h={hr}, decoy h={hd})")
        return ('CANDIDATE',real,ring)
    print(f"    >>> no clean hit (REAL h={hr}, decoy h={hd})")
    return ('NEG',real,ring)

if __name__=='__main__':
    W=sys.argv[1]; dps=int(sys.argv[2])
    # incremental principled rings
    for wmax,names in [(3,['pi','vp','G']),(3,['pi','vp','G','ln2']),
                       (3,['pi','vp','G','ln2','z3']),(3,['pi','vp','P8','G'])]:
        run(W,dps,wmax,names)
