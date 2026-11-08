import sys; sys.path.insert(0,'debug'); import mpmath as mp
import _t2_tailfibre as TF; import routeC_T2_highprec as H
mp.mp.dps=40
# reference: adaptive high Nk decay-map fibre
def ref(s,t):
    cmin=min(s*(1-s),t*(1-t)); Nk=int(max(500, 60/mp.sqrt(cmin)))
    return H.J(s,t,min(Nk,3000))
worst=(mp.mpf(0),None)
pts=[]
import itertools
for s in ['0.5','0.3','0.1','0.03','0.01','0.003','0.001']:
    for t in ['0.5','0.2','0.05','0.01','0.002']:
        pts.append((mp.mpf(s),mp.mpf(t)))
print(f'grid validation, {len(pts)} points:',flush=True)
for s,t in pts:
    v=TF.Jtail(s,t); r=ref(s,t); e=abs(v-r)
    rel=e/abs(r) if abs(r)>0 else e
    if e>worst[0]: worst=(e,(mp.nstr(s,3),mp.nstr(t,3)),mp.nstr(rel,3))
print('  MAX abs error:', mp.nstr(worst[0],3), 'at', worst[1], ' rel', worst[2],flush=True)
print('DONE',flush=True)
