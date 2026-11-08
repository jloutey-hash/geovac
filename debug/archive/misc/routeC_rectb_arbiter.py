"""Fast arbiter for RectB truth value (M=192 k-grid, already shown >=24-digit
accurate in RectB by routeC_rectb_resolve.py's k-grid check). GL refinement vs
tanh-sinh (independent outer family). Whichever they co-converge to is truth."""
from __future__ import annotations
import sys, time
import mpmath as mp
sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed
from routeC_probe8 import k_grid_sinh_paneled
from routeC_hp_evaluator import std_gl_nodes
from routeC_rectb_resolve import rectB_gl, rectB_tanhsinh

mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 35
delta = mp.mpf('0.05')
kn = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), 4)  # M=192
print(f"dps={mp.mp.dps} M={len(kn)}", flush=True)

print("tanh-sinh (independent family):", flush=True)
t0 = time.time(); ts = rectB_tanhsinh(delta, kn)
print(f"  RectB(TS) = {mp.nstr(ts, mp.mp.dps-4)}  ({time.time()-t0:.1f}s)", flush=True)

prev = None
for d in [4, 5, 6]:
    t0 = time.time(); v = rectB_gl(d, d, delta, kn)
    diff = "" if prev is None else mp.nstr(abs(v-prev), 4)
    print(f"GL deg={d}: {mp.nstr(v, mp.mp.dps-4)}  diff={diff}  |GL-TS|={mp.nstr(abs(v-ts),4)}  ({time.time()-t0:.1f}s)", flush=True)
    prev = v
