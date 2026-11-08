import sys; sys.path.insert(0,'debug')
import mpmath as mp
mp.mp.dps=50
from _t2_decayfibre import Jdecay
from _t2_fibre_stress import J_quadosc
from _t2_tailfibre import Jtail
for (s,t) in [(mp.mpf('0.5'),mp.mpf('0.999')), (mp.mpf('0.11'),mp.mpf('0.999'))]:
    jq=J_quadosc(s,t)
    jt=Jtail(s,t)
    jd2=Jdecay(s,t,1600); jd3=Jdecay(s,t,2600)
    print(f's={float(s)} t={float(t)}:',flush=True)
    print(f'   quadosc      {mp.nstr(jq,34)}',flush=True)
    print(f'   Jtail        {mp.nstr(jt,34)}  d(quad)={mp.nstr(abs(jt-jq),3)}',flush=True)
    print(f'   Jdecay 1600  {mp.nstr(jd2,34)}  d(quad)={mp.nstr(abs(jd2-jq),3)}',flush=True)
    print(f'   Jdecay 2600  {mp.nstr(jd3,34)}  d(quad)={mp.nstr(abs(jd3-jq),3)}  d(1600-2600)={mp.nstr(abs(jd2-jd3),3)}',flush=True)
