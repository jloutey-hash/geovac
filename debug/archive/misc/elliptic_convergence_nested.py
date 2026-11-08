"""Clean NESTED radial convergence (add one function, keep the rest) + explicit channel split.
Removes the even-tempered shift confound. Quantifies: radial tail (where the elliptic period
lives) vs the angular/cusp gap (where the accuracy actually is)."""
import os, sys, time, json
import numpy as np
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "debug"))
from elliptic_basis_pilot import h2_energy
from geovac.sturmian_integrals import GoscinskianIntegrals

R = 1.4
gi = GoscinskianIntegrals(R, Lmax=14, nr=1600, nth=120, rmax=50.0)

# nested geometric sequence anchored at 1.15, ratio 1.5 (each K adds ONE new exponent)
base, ratio = 1.15, 1.5
print(f"NESTED 1s ladder (anchor {base}, ratio {ratio}); each K adds one exponent:")
rows = []
for K in range(1, 6):
    exps = [base * ratio ** i for i in range(K)]
    t0 = time.time()
    E = h2_energy(exps, R, gi)
    rows.append((K, E)); print(f"  K={K}  E_tot={E:.6f}  [{time.time()-t0:.0f}s]  exps={[round(e,2) for e in exps]}")

print("\nnested increments / ratios (clean = value of the K-th function):")
Es = [e for _, e in rows]; prevd = None
for K in range(2, len(Es) + 1):
    d = Es[K - 1] - Es[K - 2]
    rat = d / prevd if prevd else float('nan')
    print(f"  dE({K}) = {d:+.2e}   ratio {rat:.3f}" if not np.isnan(rat) else f"  dE({K}) = {d:+.2e}")
    prevd = d

# geometric extrap of the radial s-limit
r = (Es[-1] - Es[-2]) / (Es[-2] - Es[-3])
Einf = Es[-1] + (Es[-1] - Es[-2]) * r / (1 - r)
EXACT = -1.17447  # H2 R=1.4 BO (Kolos-Wolniewicz)
print(f"\n  radial s-limit  E_inf(s)  ~ {Einf:.6f}  (ratio r={r:.3f} -> {'smooth/geometric' if 0<r<1 else 'STRUCTURED'})")
print(f"  radial tail remaining at K=5: {abs(Es[-1]-Einf):.2e} Ha  <-- MOST a modular resummation could buy")
print(f"  ANGULAR/cusp gap  E_inf(s) - E_exact = {Einf - EXACT:+.4f} Ha  <-- where the accuracy actually is (period cannot reach it)")
print(f"  channel ratio: angular gap / radial tail ~ {abs(Einf-EXACT)/max(abs(Es[-1]-Einf),1e-9):.0f}x")
json.dump({"rows": rows, "Einf_s": Einf, "exact": EXACT}, open(os.path.join(REPO, "debug/data/elliptic_convergence_nested.json"), "w"), indent=1)
