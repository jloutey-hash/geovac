import mpmath as mp, time, sys
sys.path.insert(0,'debug'); import routeC_T2_highprec as H
mp.mp.dps=50
ref=mp.mpf('0.3953557659017139641')
print('Fix Nc=44, sweep Nk (is Nk the ~1e-15 floor?):',flush=True)
prev=None
for Nk in [120,220,320,440]:
    t0=time.time(); v=H.T2(mp.mpf('0.08'),44,44,Nk); dt=time.time()-t0
    d='' if prev is None else f'  |dNk|={mp.nstr(abs(v-prev),3)}'
    print(f'  Nk={Nk}: {mp.nstr(v,26)}  ({dt:.0f}s){d}  |v-ref19|={mp.nstr(abs(v-ref),3)}',flush=True)
    prev=v
print('DONE',flush=True)
