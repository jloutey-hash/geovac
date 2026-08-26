"""Is quadosc accurate at the DEEPEST oscillatory-corner nodes (cmin~1e-4..1e-5) that Nrect=72
samples? Cross-check quadosc vs high-Nk Jdecay vs Jtail. If quadosc caps ~14 dig here, the T2
plateau is a FIBRE issue at deep nodes, not outer-GL convergence."""
import sys, time; sys.path.insert(0,'debug')
import mpmath as mp
mp.mp.dps=44
from _t2_fibre_stress import J_quadosc
from _t2_decayfibre import Jdecay
from _t2_tailfibre import Jtail
for (s,t,tag) in [(mp.mpf('0.999'),mp.mpf('0.999'),'(1,1) cmin~1e-3'),
                  (mp.mpf('0.9997'),mp.mpf('0.9997'),'(1,1) cmin~3e-4'),
                  (mp.mpf('0.9999'),mp.mpf('0.9999'),'(1,1) cmin~1e-4')]:
    jq=J_quadosc(s,t)
    jt=Jtail(s,t)
    jd1=Jdecay(s,t,2000); jd2=Jdecay(s,t,3200)
    dconv=float(-mp.log10(abs(jd1-jd2)/abs(jd2))) if jd1!=jd2 else 99
    print(f'{tag}:',flush=True)
    print(f'   quadosc     {mp.nstr(jq,30)}',flush=True)
    print(f'   Jtail       {mp.nstr(jt,30)}  d(quad)={mp.nstr(abs(jt-jq),3)}',flush=True)
    print(f'   Jdecay3200  {mp.nstr(jd2,30)}  d(quad)={mp.nstr(abs(jd2-jq),3)}  Jdecay2000-3200conv={dconv:.1f}dig',flush=True)
