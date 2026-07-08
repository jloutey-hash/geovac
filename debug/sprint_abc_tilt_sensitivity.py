"""
A/B/C tilt sensitivity battery -- balanced LiH n_max=2, live solver.
===================================================================
Tests the three ahha connections at their decisive point: WHICH truncation
knob carries the geometry-defect tilt eps'(R_true) = -0.030 Ha/bohr?

Reference-free anchor: at the TRUE minimum R_true=3.015, dE_true/dR=0, and the
nuclear repulsion V_NN(R)=Z_Li*Z_H/R = 3/R is EXACT, so
    d<electronic>_true/dR = -dV_NN/dR = +3/R^2 = +0.330 Ha/bohr   (required).
The computed tilt is  dE_c/dR = dV_NN/dR + d<el>_c/dR, and since V_NN is exact,
    tilt = d<el>_c/dR - d<el>_true/dR    (the pure electronic-gradient error).

Knobs (each isolates one connection):
  n_grid_vne : radial quadrature of cross-center V_ne  -> is the tilt just
               unconverged quadrature? (mundane, neither Layer-1 nor Layer-2)
  L_max      : multipole/angular order of cross-center V_ne -> ANGULAR/discrete
               (skeleton) truncation. Sensitivity => Layer-1 origin (B), and A's
               "irreducible free-side wall" is wrong (it's angular truncation).
  max_n      : one-particle BASIS completeness -> Pulay/basis-response (C).
               Frozen tilt+curv while energy heals => gradient-lags-energy (C).

Output: debug/data/abc_tilt_sensitivity.json  (+ console table).
"""
from __future__ import annotations
import json, time
import numpy as np
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import lih_spec

R_TRUE = 3.015
CM1_TO_HA = 1.0 / 219474.6313702
MU_ME = (7.016004*1.007825)/(7.016004+1.007825) * 1822.888486
K_TRUE = MU_ME * (1405.65*CM1_TO_HA)**2          # 0.0659 Ha/bohr^2
Z_LI, Z_H = 3.0, 1.0
DVNN_DR = -Z_LI*Z_H / R_TRUE**2                    # exact dV_NN/dR at R_true = -0.330
REQ_ELEC_GRAD = -DVNN_DR                           # required d<el>/dR = +0.330

H = 0.05
STENCIL = np.array([-2, -1, 0, 1, 2]) * H          # 5-pt about R_true


def energy(R: float, max_n: int, n_grid_vne: int, L_max: int):
    spec = lih_spec(R=R, max_n=max_n)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(spec, R=R, n_grid_vne=n_grid_vne, L_max=L_max,
                                     screened_cross_center=False, verbose=False)
    fci = coupled_fci_energy(ham, n_electrons=n_e, verbose=False)
    vnn = ham.get('nuclear_repulsion', None)
    return float(fci['E_coupled']), (float(vnn) if vnn is not None else None)


def tilt_curv(max_n: int, n_grid_vne: int, L_max: int, label: str):
    t0 = time.perf_counter()
    Es, vnns = [], []
    for dx in STENCIL:
        E, vnn = energy(R_TRUE + dx, max_n, n_grid_vne, L_max)
        Es.append(E); vnns.append(vnn)
    Es = np.array(Es)
    c = np.polyfit(STENCIL, Es, 3)                  # cubic over 5 pts
    p = np.poly1d(c); dp = p.deriv(1); ddp = p.deriv(2)
    tilt = float(dp(0.0)); curv = float(ddp(0.0))
    elec_grad = tilt - DVNN_DR                       # = tilt + 0.330
    # verify the nuclear_repulsion FIELD carries exact d(3/R)/dR = -3/R^2 across
    # the stencil (field = R-independent core offset + 3/R; derivative isolates 3/R)
    vnn_ok = None
    if all(v is not None for v in vnns):
        cv = np.polyfit(STENCIL, np.array(vnns), 3)
        vnn_dr = float(np.poly1d(cv).deriv(1)(0.0))
        vnn_ok = vnn_dr - DVNN_DR                    # should be ~0 (both = -0.330)
    dt = time.perf_counter() - t0
    row = {
        'label': label, 'max_n': max_n, 'n_grid_vne': n_grid_vne, 'L_max': L_max,
        'tilt': tilt, 'curv': curv, 'curv_over_ktrue': curv/K_TRUE,
        'elec_grad': elec_grad, 'elec_grad_required': REQ_ELEC_GRAD,
        'elec_grad_deficit_pct': (REQ_ELEC_GRAD - elec_grad)/REQ_ELEC_GRAD*100,
        'E_at_Rtrue': float(p(0.0)), 'vnn_max_err_vs_3overR': vnn_ok,
        'wall_s': dt,
    }
    print(f"  {label:24s}  tilt={tilt:+.4f}  curv={curv:+.4f} ({curv/K_TRUE:.2f}x k)"
          f"  elec_grad={elec_grad:+.4f} (need {REQ_ELEC_GRAD:+.3f}, "
          f"deficit {row['elec_grad_deficit_pct']:+.1f}%)  [{dt:.0f}s]", flush=True)
    return row


rows = []
print("="*88)
print(f"BASE reproduce (n_max=2, n_grid_vne=8000, L_max=4)  [finegrid gave tilt -0.0298, curv 0.148]")
print("="*88)
rows.append(tilt_curv(2, 8000, 4, "base n2/g8000/L4"))

print("\n" + "="*88)
print("KNOB 1: n_grid_vne (radial quadrature of cross V_ne)  -- is tilt converged?")
print("="*88)
for g in (2000, 4000, 16000):
    rows.append(tilt_curv(2, g, 4, f"n2/g{g}/L4"))

print("\n" + "="*88)
print("KNOB 2: L_max (angular/multipole order of cross V_ne)  -- Layer-1/skeleton knob")
print("="*88)
for L in (2, 3, 6, 8):
    rows.append(tilt_curv(2, 8000, L, f"n2/g8000/L{L}"))

out = {
    'meta': {'R_true': R_TRUE, 'k_true': K_TRUE, 'dVNN_dR': DVNN_DR,
             'required_elec_grad': REQ_ELEC_GRAD, 'stencil_h': H},
    'rows': rows,
}
json.dump(out, open('debug/data/abc_tilt_sensitivity.json', 'w'), indent=2)
print("\n[saved] debug/data/abc_tilt_sensitivity.json")

# ---- summary read-out ----
base = rows[0]
print("\n" + "="*88); print("VERDICT READ-OUT"); print("="*88)
gvals = [r for r in rows if r['L_max']==4 and r['max_n']==2]
Lvals = [r for r in rows if r['n_grid_vne']==8000 and r['max_n']==2]
print(f"  tilt vs n_grid_vne: " + ", ".join(f"g{r['n_grid_vne']}={r['tilt']:+.4f}" for r in sorted(gvals,key=lambda r:r['n_grid_vne'])))
print(f"  tilt vs L_max     : " + ", ".join(f"L{r['L_max']}={r['tilt']:+.4f}" for r in sorted(Lvals,key=lambda r:r['L_max'])))
print(f"  curv vs L_max     : " + ", ".join(f"L{r['L_max']}={r['curv']:+.4f}" for r in sorted(Lvals,key=lambda r:r['L_max'])))
print(f"  (V_NN=3/R exact check, max abs err: {base['vnn_max_err_vs_3overR']})")
