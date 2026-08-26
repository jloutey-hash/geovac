"""Cross-check TwoCenterLM one-electron l>0 integrals vs the independent oracle (machine-prec)."""
import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from two_center_grid_lm import TwoCenterLM
from one_electron_lm_oracle import Orbital, s_ref, t_ref, v_ref

R = 1.5
eng = TwoCenterLM(R, nr=1400, nu=32, nphi=32, rmax=50.0, Lmax=14)

# (zeta,l,m,center) for my engine  <->  Orbital(center,zeta,l,m) for oracle
cases = [
    ((1.0, 1, 0, "A"), (1.2, 1, 0, "B")),    # pz-pz cross-center, mixed exp
    ((1.0, 1, 0, "A"), (0.9, 0, 0, "B")),    # pz-s cross-center
    ((1.1, 1, 1, "A"), (0.8, 1, 1, "B")),    # ppi-ppi cross-center (m=+1)
    ((1.3, 1, -1, "A"), (1.3, 1, -1, "A")),  # ppi same-center diagonal
    ((1.0, 1, 1, "A"), (1.0, 1, 0, "B")),    # ppi-pz cross (should vanish by m)
    ((0.9, 1, 0, "B"), (1.4, 1, 0, "B")),    # pz-pz same-center B, mixed exp
]
print(f"one-electron l>0 cross-check vs oracle (R={R}):  |diff|")
print(f"  {'case':30s} {'S':>10} {'T':>10} {'V_ne':>10}")
worst = 0.0
for gi, gj in cases:
    oi = Orbital(gi[3], gi[0], gi[1], gi[2]); oj = Orbital(gj[3], gj[0], gj[1], gj[2])
    dS = abs(eng.overlap(gi, gj) - complex(s_ref(oi, oj, R)).real)
    dT = abs(eng.kinetic(gi, gj) - complex(t_ref(oi, oj, R)).real)
    dV = abs(eng.nuclear(gi, gj) - complex(v_ref(oi, oj, R)).real)
    worst = max(worst, dS, dT, dV)
    lab = f"{gi[:3]}{gi[3]}-{gj[:3]}{gj[3]}"
    print(f"  {lab:30s} {dS:>10.1e} {dT:>10.1e} {dV:>10.1e}")
print(f"\nworst residual: {worst:.1e}  -> {'PASS (<2e-4 grid)' if worst < 2e-4 else 'CHECK'}")
