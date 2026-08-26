"""DECISIVE test for the RectB=[delta,1]x[0,1] discrepancy (Paper 59 T2).

Two candidate values disagree at digit 7-9:
  GL-family (single & dyadic panel):  RectB ~ 0.1536499096283...
  Mobius-family (deg 4):              RectB ~ 0.1536499726615...
and the GL-composite total (0.3953557662787...) disagrees with the trusted
~16-digit anchor (0.3953557659017139) at digit 9.  The whole discrepancy is
in RectB (T1 corner, T2far, RectA all agree across methods to >=10 digits).

This script settles it two ways:
  (1) k-GRID CHECK: J_fixed (fixed paneled-sinh k-grid) vs J_adaptive
      (adaptive mp.quad k-integral, extended breakpoints) at RectB edge +
      interior (s,t) points -- confirms the SHARED k-grid is not the bias
      (near s->1 slow k-decay, and t->0).
  (2) OUTER CONVERGENCE: RectB via GL at increasing degree, AND via tanh-sinh
      nested mp.quad (a genuinely different outer family from both GL and
      Mobius).  Whichever value GL and tanh-sinh CO-converge to is the truth.

Usage: python routeC_rectb_resolve.py [dps]
"""
from __future__ import annotations
import sys
import time
import mpmath as mp

sys.path.insert(0, r'C:\Users\jlout\Desktop\Project_Geometric\debug')
from routeC_probe9 import J_fixed, J_adaptive
from routeC_probe8 import k_grid_sinh_paneled
from routeC_hp_evaluator import std_gl_nodes


def rectB_gl(deg_s, deg_t, delta, kn):
    std_s = std_gl_nodes(deg_s, mp.mp.prec)
    std_t = std_gl_nodes(deg_t, mp.mp.prec)
    slo, shi, tlo, thi = delta, mp.mpf(1), mp.mpf(0), mp.mpf(1)
    hs, ms = (shi - slo) / 2, (shi + slo) / 2
    ht, mt = (thi - tlo) / 2, (thi + tlo) / 2
    tot = mp.mpf(0)
    for xs, ws in std_s:
        s = ms + hs * xs
        acc = mp.mpf(0)
        for xt, wt in std_t:
            t = mt + ht * xt
            acc += wt * J_fixed(s, t, kn)
        tot += ws * acc
    return (hs * ht) * tot


def rectB_tanhsinh(delta, kn):
    """Nested tanh-sinh. Inner t in [0,1] with endpoint breakpoints (mild
    zeros at t=0,1); outer s in [delta,1] with a breakpoint near s=1 (slow
    k-decay edge)."""
    def inner(s):
        return mp.quad(lambda t: J_fixed(s, t, kn), [0, mp.mpf('0.5'), 1])
    return mp.quad(inner, [delta, mp.mpf('0.5'), mp.mpf('0.9'), 1])


def main():
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    mp.mp.dps = dps
    delta = mp.mpf('0.05')
    kn = k_grid_sinh_paneled(mp.mpf(16), mp.mpf(2), 6)  # M=? (pdeg6 -> 96/panel)
    print(f"dps={dps} delta={delta} k_M={len(kn)}", flush=True)

    # (1) k-grid check: J_fixed vs J_adaptive at RectB probe points
    print("\n(1) k-GRID CHECK  J_fixed vs J_adaptive (rel.diff)", flush=True)
    probes = [
        ('interior', mp.mpf('0.5'), mp.mpf('0.5')),
        ('t->0    ', mp.mpf('0.5'), mp.mpf('1e-3')),
        ('t->0,s=d', delta,         mp.mpf('1e-3')),
        ('s->1    ', mp.mpf('0.98'), mp.mpf('0.3')),
        ('s->1,t->0', mp.mpf('0.98'), mp.mpf('1e-3')),
        ('s=d,t=.3', delta,         mp.mpf('0.3')),
    ]
    for name, s, t in probes:
        jf = J_fixed(s, t, kn)
        ja = J_adaptive(s, t)
        rd = abs(jf - ja) / abs(ja) if ja != 0 else abs(jf - ja)
        print(f"    {name}  s={mp.nstr(s,4)} t={mp.nstr(t,4)}  "
              f"Jfix={mp.nstr(jf,18)}  reldiff={mp.nstr(rd,4)}", flush=True)

    # (2) outer convergence: GL refinement
    print("\n(2a) RectB via GL, degree refinement", flush=True)
    prev = None
    for d in [4, 5, 6, 7]:
        t0 = time.time()
        v = rectB_gl(d, d, delta, kn)
        diff = "" if prev is None else mp.nstr(abs(v - prev), 4)
        print(f"    deg={d}  RectB={mp.nstr(v, dps-4)}  diff={diff}  ({time.time()-t0:.1f}s)", flush=True)
        prev = v
    gl_val = prev

    # (2b) tanh-sinh independent family
    print("\n(2b) RectB via nested tanh-sinh (mp.quad)", flush=True)
    t0 = time.time()
    ts_val = rectB_tanhsinh(delta, kn)
    print(f"    RectB(tanh-sinh)={mp.nstr(ts_val, dps-4)}  ({time.time()-t0:.1f}s)", flush=True)

    print(f"\n|GL - tanhsinh| = {mp.nstr(abs(gl_val - ts_val), 6)}", flush=True)
    print(f"GL      = {mp.nstr(gl_val, dps-4)}")
    print(f"tanhsinh= {mp.nstr(ts_val, dps-4)}")


if __name__ == '__main__':
    main()
