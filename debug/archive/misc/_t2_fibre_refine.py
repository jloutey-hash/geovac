"""Push the fixed-grid fibre toward its limit: at fixed high Nc, sweep MATCHED (Nk,Kmax)
configs (density ~ Nk/Kmax roughly constant) and Shanks-extrapolate the fibre limit.
Nc chosen so Nc-error << fibre-error (Nc=200 gives Nc-error ~6e-27)."""
import mpmath as mp, sys, time
sys.path.insert(0,'debug'); import _t2_tensor as T
mp.mp.dps=int(sys.argv[1]) if len(sys.argv)>1 else 44
Nc=int(sys.argv[2]) if len(sys.argv)>2 else 200
# matched configs: raise Nk and Kmax together
configs=[(1600,90),(2200,105),(2800,120),(3400,135),(4000,150)]
if len(sys.argv)>3:
    configs=[tuple(int(y) for y in c.split(':')) for c in sys.argv[3].split(',')]
print(f'== fibre refine at Nc={Nc}, dps={mp.mp.dps} ==',flush=True)
vals=[]
for Nk,Km in configs:
    t0=time.time(); v=T.T2_tensor(Nc,Nk,Km); dt=time.time()-t0
    d='' if not vals else f'  |dprev|={mp.nstr(abs(v-vals[-1]),3)}'
    print(f'  Nk={Nk} Kmax={Km}: {mp.nstr(v,mp.mp.dps-4)}  ({dt:.0f}s){d}',flush=True)
    vals.append(v)
# Shanks tower on the (Nk,Kmax) sequence
def sh(s):
    o=[]
    for i in range(1,len(s)-1):
        d1=s[i+1]-s[i]; d0=s[i]-s[i-1]; den=d1-d0
        o.append(s[i+1]-d1*d1/den if den!=0 else s[i+1])
    return o
s=list(vals); lvl=0
print('-- Shanks tower --',flush=True)
while len(s)>=3:
    s=sh(s); lvl+=1
    print(f'  L{lvl}: {mp.nstr(s[-1],mp.mp.dps-4)}'+('' if len(s)<2 else f'   selfconsist={mp.nstr(abs(s[-1]-s[-2]),3)}'),flush=True)
print('DONE',flush=True)
