import sys; sys.path.insert(0,'debug')
import mpmath as mp
mp.mp.dps=45
from _t2_decayfibre import Jdecay
from _t2_fibre_stress import J_quadosc
pts=[(mp.mpf('0.5'),mp.mpf('0.999'),'edge cmin=1e-3 osc=2.8'),
     (mp.mpf('0.5'),mp.mpf('0.5'),'interior cmin=.25 osc=1'),
     (mp.mpf('0.3'),mp.mpf('0.97'),'osc=2.0 cmin=.029'),
     (mp.mpf('0.999'),mp.mpf('0.89'),'osc=5.5 cmin=1e-3')]
for (s,t,tag) in pts:
    ref=J_quadosc(s,t)
    cs=s*(1-s); ct=t*(1-t); osc=float((s+t)/(mp.sqrt(cs)+mp.sqrt(ct)))
    row=f'{tag:26} osc={osc:.2f}: '
    for Nk in (500,800,1100,1500):
        v=Jdecay(s,t,Nk); ag=float(-mp.log10(abs(v-ref)/abs(ref))) if v!=ref else 99
        row+=f'Nk{Nk}={ag:.0f} '
    print(row,flush=True)
