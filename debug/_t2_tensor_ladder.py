"""Clean tensor convergence study: fix Nk,Kmax high; vary ONLY Nc (outer). Then a Nk/Kmax
saturation check. Enables clean Richardson/Shanks extrapolation to ~32+ digits."""
import mpmath as mp, sys, time
sys.path.insert(0,'debug')
import _t2_tensor as T
mp.mp.dps=48
ref=mp.mpf('0.3953557659017139641')
print('== Nc-ladder at FIXED Nk=1200, Kmax=80 ==',flush=True)
vals=[]; ncs=[80,120,160,200,240]
for Nc in ncs:
    t0=time.time(); v=T.T2_tensor(Nc,1200,80); dt=time.time()-t0
    d='' if not vals else f'  |dprev|={mp.nstr(abs(v-vals[-1]),3)}'
    print(f'  Nc={Nc}: {mp.nstr(v,38)}  ({dt:.0f}s){d}  |anch|={mp.nstr(abs(v-ref),3)}',flush=True)
    vals.append(v)
print('== Nk saturation at Nc=160, Kmax=80 ==',flush=True)
prev=None
for Nk in [800,1200,1800]:
    t0=time.time(); v=T.T2_tensor(160,Nk,80); dt=time.time()-t0
    d='' if prev is None else f'  |dNk|={mp.nstr(abs(v-prev),3)}'
    print(f'  Nk={Nk}: {mp.nstr(v,38)}  ({dt:.0f}s){d}',flush=True); prev=v
print('== Kmax saturation at Nc=160, Nk=1200 ==',flush=True)
prev=None
for Km in [60,80,110]:
    t0=time.time(); v=T.T2_tensor(160,1200,Km); dt=time.time()-t0
    d='' if prev is None else f'  |dKm|={mp.nstr(abs(v-prev),3)}'
    print(f'  Kmax={Km}: {mp.nstr(v,38)}  ({dt:.0f}s){d}',flush=True); prev=v
print('DONE',flush=True)
