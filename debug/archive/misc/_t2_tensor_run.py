"""Parametrized tensor Nc-ladder for the high-precision push.
Usage: python _t2_tensor_run.py <Nk> <Kmax> <dps> <Nc_csv>"""
import mpmath as mp, sys, time
sys.path.insert(0,'debug'); import _t2_tensor as T
Nk=int(sys.argv[1]); Kmax=int(sys.argv[2]); dps=int(sys.argv[3])
ncs=[int(x) for x in sys.argv[4].split(',')]
mp.mp.dps=dps
ref=mp.mpf('0.3953557659017139641')
print(f'== tensor Nc-ladder Nk={Nk} Kmax={Kmax} dps={dps} ==',flush=True)
vals=[]
for Nc in ncs:
    t0=time.time(); v=T.T2_tensor(Nc,Nk,Kmax); dt=time.time()-t0
    d='' if not vals else f'  |dprev|={mp.nstr(abs(v-vals[-1]),3)}'
    print(f'  Nc={Nc}: {mp.nstr(v,dps-4)}  ({dt:.0f}s){d}',flush=True)
    vals.append(v)
print('VALUES:',','.join(mp.nstr(v,dps-4) for v in vals),flush=True)
print('DONE',flush=True)
