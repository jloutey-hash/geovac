"""Hybrid: true = tensor(inf) + delta_inf, delta_inf = lim[tail_tensor(Nc) - tensor(Nc)] on a
MATCHED sin^2 grid (slow outer error cancels in the difference; only the smooth fibre correction
survives). Measures delta(Nc) at several Nc to confirm fast convergence."""
import mpmath as mp, sys, time
sys.path.insert(0,'debug')
import _t2_tensor as T           # fixed-grid fibre, sin^2 outer
import _t2_tailfibre as TF       # accurate tail fibre
from _fastgl import fast_gl

def tail_tensor(Nc):
    # SAME sin^2 outer as T.T2_tensor, but the accurate tail fibre
    H=mp.pi/2; xp,wp=fast_gl(Nc)
    sv=[]; sj=[]
    for x,w in zip(xp,wp):
        phi=H*(x+1)/2; sv.append(mp.sin(phi)**2); sj.append(mp.sin(2*phi)*H*w/2)
    tot=mp.mpf(0)
    for i in range(Nc):
        si=sv[i]; wi=sj[i]
        for j in range(Nc):
            tot+=wi*sj[j]*TF.Jtail(si,sv[j])
    return (8/mp.pi)*tot

if __name__=='__main__':
    mp.mp.dps=int(sys.argv[1]) if len(sys.argv)>1 else 50
    Nk,Km=1200,80
    ncs=[int(x) for x in (sys.argv[2].split(',') if len(sys.argv)>2 else ['50','70','90'])]
    print(f'delta(Nc)=tail_tensor - tensor(Nk={Nk},Km={Km}), matched sin^2 grid, dps={mp.mp.dps}',flush=True)
    prev=None
    for Nc in ncs:
        t0=time.time(); a=tail_tensor(Nc); tb=time.time(); b=T.T2_tensor(Nc,Nk,Km); dt=time.time()-t0
        d=a-b
        dd='' if prev is None else f'  |d delta|={mp.nstr(abs(d-prev),3)}'
        print(f'  Nc={Nc}: delta={mp.nstr(d,mp.mp.dps-6)}  (tail {tb-t0:.0f}s + tensor {time.time()-tb:.0f}s){dd}',flush=True)
        prev=d
    print('DONE',flush=True)
