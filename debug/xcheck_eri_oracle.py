"""Cross-check TwoCenterLM cross-center l>0 ERIs vs the independent cross-center oracle."""
import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from two_center_grid_lm import TwoCenterLM
from eri_crosscenter_oracle import eri_ref

R = 1.5
for Lmax in (14, 22):
    eng = TwoCenterLM(R, nr=1400, nu=32, nphi=32, rmax=50.0, Lmax=Lmax)
    cases = [
        ((1.0,1,0,"A"),(1.0,0,0,"B"),(1.0,1,0,"A"),(1.0,0,0,"B")),   # (pzA sB|pzA sB) cross, l>0
        ((1.0,1,0,"A"),(1.2,0,0,"B"),(0.9,0,0,"A"),(1.1,1,0,"B")),   # mixed exp, cross, l>0
        ((1.0,1,1,"A"),(1.0,1,1,"B"),(1.0,1,1,"A"),(1.0,1,1,"B")),   # ppi cross-center, m=1
        ((1.1,1,0,"A"),(0.8,0,0,"B"),(0.8,0,0,"B"),(1.1,1,0,"A")),   # exchange-type cross
    ]
    worst = 0.0; print(f"Lmax={Lmax}:")
    for a, b, c, d in cases:
        mine = eng.eri(a, b, c, d)
        ref = complex(eri_ref(a, b, c, d, R)).real
        worst = max(worst, abs(mine - ref))
        print(f"  {a[1:3]}{a[3]} {b[1:3]}{b[3]}|{c[1:3]}{c[3]} {d[1:3]}{d[3]}: mine={mine:+.7f} oracle={ref:+.7f} diff={mine-ref:+.1e}")
    print(f"  worst={worst:.1e}\n")
