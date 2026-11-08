"""
Fine-grid confirmation of the LiH balanced tilt/curvature autopsy (n_max=2, live).
Validates the cached analytical curve == current solver, and gives a robust
eps'(R_true) and E_c''(R_true) from a dense grid rather than a 6-point fit.
"""
from __future__ import annotations
import json, time
import numpy as np
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import lih_spec

R_TRUE, E_EXACT = 3.015, -8.0706
CM1_TO_HA = 1.0 / 219474.6313702
MU_ME = (7.016004*1.007825)/(7.016004+1.007825) * 1822.888486
K_TRUE = MU_ME * (1405.65*CM1_TO_HA)**2

def E_of_R(R, max_n=2, n_grid=8000):
    spec = lih_spec(R=R, max_n=max_n)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(spec, R=R, n_grid_vne=n_grid, L_max=4,
                                     screened_cross_center=False, verbose=False)
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    return float(fci['E_coupled'])

grid = np.round(np.arange(2.80, 3.66, 0.05), 3)
print(f"Running balanced LiH n_max=2 on {len(grid)} points (n_grid_vne=8000)...")
Es = []
t0 = time.perf_counter()
for i, R in enumerate(grid):
    ta = time.perf_counter()
    E = E_of_R(R)
    Es.append(E)
    print(f"  [{i+1:2d}/{len(grid)}] R={R:.3f}  E={E:.6f}  ({time.perf_counter()-ta:.1f}s)", flush=True)
Es = np.array(Es)
print(f"total {time.perf_counter()-t0:.1f}s")

# robust local fit about R_TRUE (quartic over the whole dense window)
x = grid - R_TRUE
c = np.polyfit(x, Es, 4)
p = np.poly1d(c); dp = p.deriv(1); ddp = p.deriv(2)
# also a low-order (cubic) cross-check on the near-min window
mnear = np.abs(x) <= 0.35
c3 = np.polyfit(x[mnear], Es[mnear], 3); p3 = np.poly1d(c3)
dp3 = p3.deriv(1); ddp3 = p3.deriv(2)

roots = dp.r; roots = roots[np.isreal(roots)].real
mins = [r for r in roots if ddp(r) > 0 and grid.min()-R_TRUE < r < grid.max()-R_TRUE]
r_eq_x = min(mins, key=lambda r: p(r)) if mins else x[np.argmin(Es)]
R_eq = r_eq_x + R_TRUE

res = {
    'grid': grid.tolist(), 'E': Es.tolist(),
    'R_eq': float(R_eq), 'R_eq_err_pct': abs(R_eq-R_TRUE)/R_TRUE*100,
    'eps_true_Ha': float(p(0.0) - E_EXACT),
    'tilt_quartic': float(dp(0.0)), 'curv_quartic': float(ddp(0.0)),
    'tilt_cubic_near': float(dp3(0.0)), 'curv_cubic_near': float(ddp3(0.0)),
    'k_true': float(K_TRUE),
}
with open('debug/data/pk_amplification_finegrid.json', 'w') as f:
    json.dump(res, f, indent=2)

print("\n--- fine-grid autopsy (balanced n_max=2, live solver) ---")
print(f"  R_eq            = {R_eq:.3f} bohr  ({res['R_eq_err_pct']:.1f}% err)  "
      f"[cached analytical: 3.227]")
print(f"  eps(R_true)     = {res['eps_true_Ha']:+.4f} Ha  ({abs(res['eps_true_Ha'])/abs(E_EXACT)*100:.2f}%)")
print(f"  tilt  eps'(R_t) = {res['tilt_quartic']:+.4f} (quartic) / "
      f"{res['tilt_cubic_near']:+.4f} (cubic-near) Ha/bohr   [6pt fit gave -0.030]")
print(f"  curv  E_c''(R_t)= {res['curv_quartic']:+.4f} (quartic) / "
      f"{res['curv_cubic_near']:+.4f} (cubic-near)   k_true={K_TRUE:.4f} "
      f"(ratio {res['curv_quartic']/K_TRUE:.2f}x)")
print("[saved] debug/data/pk_amplification_finegrid.json")
