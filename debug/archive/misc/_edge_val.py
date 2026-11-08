import sys; sys.path.insert(0,'debug'); import mpmath as mp
import _t2_tailfibre as TF; import routeC_T2_highprec as H
mp.mp.dps=44
print('Edge/corner validation (tail-analytic vs decay-map ref):',flush=True)
for (s,t,tag) in [(mp.mpf('0.05'),mp.mpf('0.28'),'edge s->0'),
                  (mp.mpf('0.02'),mp.mpf('0.3'),'deeper edge'),
                  (mp.mpf('0.008'),mp.mpf('0.012'),'CORNER both->0')]:
    ref=H.J(s,t,700); b=s+t
    best=None
    for K,M in [(70,4),(110,5),(160,6)]:
        v=TF.bounded(s,t,b,K,700)+TF.tail(s,t,b,K,M)
        e=abs(v-ref)
        if best is None or e<best[0]: best=(e,K,M)
    print(f'  {tag}: best |Jtail-ref|={mp.nstr(best[0],3)} at K={best[1]},M={best[2]}   (ref=H.J Nk700)',flush=True)
print('DONE',flush=True)
