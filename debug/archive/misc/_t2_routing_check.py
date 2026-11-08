"""Verify Jfast (fast decay fibre, adaptive Nk) is accurate wherever the split routes to it:
the (0,0) Duffy region and the bulk (R1,R2,R3 incl. edge strips). Compare vs quadosc.
Any point below ~34 dig => that region must be re-routed to Josc."""
import sys
import mpmath as mp
sys.path.insert(0, 'debug')
from routeC_T2_corner_subtraction import Jfast
from _t2_fibre_stress import J_quadosc

mp.mp.dps = 45

def chk(s, t, tag):
    s = mp.mpf(s); t = mp.mpf(t)
    jf = Jfast(s, t); jq = J_quadosc(s, t)
    ag = float(-mp.log10(abs(jf-jq)/abs(jq))) if jf != jq else 99.0
    cmin = float(min(s*(1-s), t*(1-t))); b = float(s+t)
    flag = '  <-- LOW' if ag < 34 else ''
    print(f"  {tag:22} s={float(s):.4f} t={float(t):.4f} cmin={cmin:.1e} b={b:.3f} : agree={ag:5.1f}dig{flag}", flush=True)

d = 0.1
print("== (0,0) Duffy region (b small, non-osc) ==", flush=True)
chk(0.05, 0.05, '(0,0) interiorish')
chk(0.005, 0.09, '(0,0) edge')
chk(0.001, 0.02, '(0,0) deep')
chk(0.09, 0.005, '(0,0) edge2')
print("== bulk R1 [d,1-d]x[0,1] incl top/bottom edges ==", flush=True)
chk(0.5, 0.5, 'R1 center')
chk(0.5, 0.001, 'R1 bottom edge')
chk(0.5, 0.999, 'R1 top edge')
chk(0.11, 0.001, 'R1 near-corner bottom')
chk(0.11, 0.999, 'R1 near-corner top')
chk(0.89, 0.001, 'R1 right-bottom')
chk(0.3, 0.97, 'R1 upper')
print("== bulk R2 [0,d]x[d,1-d] left edge strip ==", flush=True)
chk(0.001, 0.5, 'R2 mid')
chk(0.001, 0.11, 'R2 low')
chk(0.001, 0.89, 'R2 high')
chk(0.05, 0.11, 'R2 inner')
print("== bulk R3 [1-d,1]x[d,1-d] right edge strip ==", flush=True)
chk(0.999, 0.5, 'R3 mid')
chk(0.999, 0.11, 'R3 low')
chk(0.999, 0.89, 'R3 high')
print("DONE", flush=True)
