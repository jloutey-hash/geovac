import mpmath as mp, time, sys
sys.path.insert(0,'debug')
import routeC_T2_highprec as H
mp.mp.dps=35
ref=mp.mpf('0.3953557659017139641')
print('FIBRE convergence at (s,t)=(0.3,0.17), saturate Nk:',flush=True)
prevJ=None
for Nk in [60,100,140,200]:
    J=H.J(mp.mpf('0.3'),mp.mpf('0.17'),Nk)
    d='' if prevJ is None else f'  |dJ|={mp.nstr(abs(J-prevJ),3)}'
    print(f'  Nk={Nk}: {mp.nstr(J,28)}{d}',flush=True); prevJ=J
print('\nOUTER ladder, SATURATED Nk=120, refine only Nc=Nt (dps=35):',flush=True)
prev=None
for Nc in [20,28,36,44,52]:
    t0=time.time(); v=H.T2(mp.mpf('0.08'),Nc,Nc,120); dt=time.time()-t0
    d='' if prev is None else f'  |dprev|={mp.nstr(abs(v-prev),3)}'
    print(f'  Nc={Nc} Nk=120: {mp.nstr(v,22)}  ({dt:.0f}s){d}   |v-ref19|={mp.nstr(abs(v-ref),3)}',flush=True)
    prev=v
print('DONE',flush=True)
