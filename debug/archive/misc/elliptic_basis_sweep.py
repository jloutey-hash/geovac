"""beta-sweep: is a CM singular-modulus exponent ratio variationally special for H2?
Fast grid for the sweep; high grid re-check at beta* and the CM markers."""
import os, sys, time, json
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from elliptic_basis_pilot import h2_energy
from geovac.sturmian_integrals import GoscinskianIntegrals

R = 1.4
Z_LO = 1.0
FAST = dict(Lmax=14, nr=1600, nth=120, rmax=50.0)
HIGH = dict(Lmax=24, nr=3000, nth=200, rmax=60.0)

# CM ratio markers beta = z_max/z_min
m4, m8 = 0.5, float(3 - 2 * 2 ** 0.5)
CM = {
    "disc4_p2 (m=1-1/b^2)": (1 - m4) ** -0.5,   # sqrt2
    "disc8_p2":             (1 - m8) ** -0.5,
    "disc4_p1 (m=1-1/b)":   1 / (1 - m4),        # 2.0
    "disc8_p1":             1 / (1 - m8),
}

def run():
    gi_f = GoscinskianIntegrals(R, **FAST)
    # grid sanity vs the known high-grid point E([1.2,1.2*sqrt2]) = -1.152351
    chk = h2_energy([1.2, 1.2 * 2 ** 0.5], R, gi_f)
    print(f"grid check (fast) E[1.2,1.2sqrt2] = {chk:.6f}  (high-grid ref -1.152351, "
          f"diff {chk-(-1.152351):+.2e})")

    betas = np.round(np.concatenate([np.arange(1.05, 2.25, 0.075),
                                     np.arange(2.3, 3.01, 0.2)]), 4)
    rows = []
    t0 = time.time()
    for b in betas:
        e = h2_energy([Z_LO, Z_LO * b], R, gi_f)
        rows.append((float(b), e))
        print(f"  beta={b:5.3f}  E_tot={e:.6f}")
    print(f"[sweep {time.time()-t0:.0f}s]")

    bstar, estar = min(rows, key=lambda t: t[1])
    print(f"\nvariational optimum on grid: beta* = {bstar:.3f}  E = {estar:.6f}")
    print("CM markers:", {k: round(v, 4) for k, v in CM.items()})

    # high-grid re-check at beta* and each CM marker
    gi_h = GoscinskianIntegrals(R, **HIGH)
    print("\nhigh-grid re-check:")
    dec = {}
    for label, b in [("beta*", bstar)] + [(k, v) for k, v in CM.items()]:
        if 1.02 <= b <= 3.1:
            e = h2_energy([Z_LO, Z_LO * b], R, gi_h)
            dec[label] = (round(float(b), 4), e)
            print(f"  {label:22s} beta={b:5.3f}  E_tot={e:.6f}")
    json.dump({"R": R, "z_lo": Z_LO, "sweep": rows, "cm": CM, "decision": dec},
              open(os.path.join(REPO, "debug/data/elliptic_basis_sweep.json"), "w"), indent=1)
    print("\nsaved debug/data/elliptic_basis_sweep.json")

if __name__ == "__main__":
    run()
