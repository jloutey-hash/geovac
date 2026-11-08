"""
Sprint: balanced-coupled LiH R_eq-drift MECHANISM PROOF (diagnostic only).
=========================================================================
Paper 19 (Sec. "Fixed-geometry energy versus well shape", v5.14.8) NAMES the
R_eq-drift mechanism -- "the fixed, unpolarized hydrogenic bond orbital
(Z_orb=1) carries no R-adaptive relaxation" -- but does not DEMONSTRATE it.
This driver demonstrates (or refutes) it.

The balanced LiH well is a strict two-term V_NN <-> cross-V_ne balance; the 8.8%
outward drift is the cross-V_ne R-slope being 8.8% too weak. Paper 19 attributes
that to the bond orbital not contracting toward the bond as R changes. Both bond
orbitals default to Z_orb = 1.0 (Li-side Z_center AND H-side Z_partner). We give
them ONE R-adaptive shared exponent -- "the bonding density contracts toward the
bond" -- and re-measure:

    E_relaxed(R) = min_z E_balanced(R, Z_center=Z_partner=z)
    R_eq_relaxed = argmin_R E_relaxed(R)

PREDICTION if the named mechanism is right: the outward drift collapses
(R_eq -> ~3.015) and the tilt -> ~0, with z*(R) CONTRACTING (increasing) as R
shrinks. By the envelope theorem, d/dR E_relaxed = <psi| dH/dR |psi> at z*, and
only V_NN and cross-V_ne are R-dependent, so tilt_relaxed -> 0 IS the cross-V_ne
slope correcting to the required +Z_Li Z_H/R^2. If radial contraction does NOT
heal it, the missing ingredient is sharper (angular/p-polarization).

DIAGNOSTIC, not a fix: an R-dependent exponent breaks the zero-parameter
construction and re-enters the PK/Loewdin/non-orthogonal walls (Sec. 3). The
deliverable is proving WHAT the missing physics is.

Cost note: FCI ~12s/point dominates (build ~1.6s, grid-independent). Physical
energy = E_raw - E_CORE (coupled_fci_energy carries a constant +E_CORE offset;
the original termdecomp driver subtracts it too). Offset is R- and z-independent,
so R_eq/tilt are unaffected.
"""
from __future__ import annotations
import json
import time
import numpy as np
from pathlib import Path

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.molecular_spec import lih_spec, _FIRST_ROW_CORE_ENERGY

R_TRUE = 3.015
E_EXACT = -8.071
Z_LI, Z_H = 3.0, 1.0
E_CORE = _FIRST_ROW_CORE_ENERGY[3]     # -7.2799, R- and z-independent
MAX_N = 2
N_GRID_VNE = 400                       # tilt bit-identical to ng=2000 (calibrated)
L_MAX = 4

R_GRID = np.array([2.70, 2.85, 3.015, 3.15, 3.30, 3.45])
Z_GRID = np.array([0.80, 1.00, 1.20, 1.40, 1.60, 1.90, 2.30])   # 1.00 = baseline


def E_phys(R: float, z: float | None = None) -> float:
    """Physical balanced-coupled total energy at bond length R, shared bond
    exponent z (None -> spec default = 1.0)."""
    spec = lih_spec(R=R, max_n=MAX_N)
    bond = next(b for b in spec.blocks if b.block_type == 'bond')
    if z is not None:
        bond.Z_center = bond.Z_partner = float(z)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(spec, R=R, n_grid_vne=N_GRID_VNE, L_max=L_MAX)
    res = {'M': ham['M'], 'h1': ham['h1_no_pk'] + ham['h1_cross_vne'],
           'eri': ham['eri'], 'nuclear_repulsion': ham['nuclear_repulsion']}
    return float(coupled_fci_energy(res, n_electrons=n_e, verbose=False)['E_coupled']) - E_CORE


def _min_over_z(Es_z):
    """Interior minimum of E(z) by quartic fit over the z-grid; returns (E*, z*)."""
    p = np.poly1d(np.polyfit(Z_GRID, Es_z, 4))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r for r in roots if ddp(r) > 0 and Z_GRID.min() <= r <= Z_GRID.max()]
    if not cand:                       # fall back to grid min
        i = int(np.argmin(Es_z)); return float(Es_z[i]), float(Z_GRID[i])
    zstar = float(min(cand, key=lambda z: p(z)))
    return float(p(zstar)), zstar


def _quartic_req(R, E, lo=2.5, hi=3.55):
    p = np.poly1d(np.polyfit(R - R_TRUE, E, 4))
    dp, ddp = p.deriv(1), p.deriv(2)
    roots = dp.r[np.isreal(dp.r)].real
    cand = [r + R_TRUE for r in roots if ddp(r) > 0 and lo < r + R_TRUE < hi]
    return (float(min(cand, key=lambda rr: p(rr - R_TRUE))) if cand else float('nan'))


def _tilt(R, E):
    i0 = int(np.argmin(np.abs(R - R_TRUE)))
    p = np.poly1d(np.polyfit(R[i0 - 1:i0 + 2], E[i0 - 1:i0 + 2], 2))
    return float(p.deriv(1)(R_TRUE)), float(p.deriv(2)(R_TRUE))


def main():
    t0 = time.perf_counter()
    print(f"required electronic force at R_true = +Z_Li Z_H/R^2 = "
          f"{Z_LI*Z_H/R_TRUE**2:+.4f} Ha/bohr")
    print(f"z-grid {Z_GRID.tolist()}  (1.00 = fixed-basis baseline)\n")

    Etab = np.zeros((len(R_GRID), len(Z_GRID)))
    for i, R in enumerate(R_GRID):
        for j, z in enumerate(Z_GRID):
            Etab[i, j] = E_phys(R, z)
        print(f"  R={R:.3f}  E(z=1)={Etab[i, list(Z_GRID).index(1.00)]:.5f}  "
              f"E(min over z)={Etab[i].min():.5f}  [{time.perf_counter()-t0:.0f}s]")

    j1 = list(Z_GRID).index(1.00)
    E_fixed = Etab[:, j1].copy()
    E_relax = np.empty(len(R_GRID)); z_star = np.empty(len(R_GRID))
    for i in range(len(R_GRID)):
        E_relax[i], z_star[i] = _min_over_z(Etab[i])

    print("\n" + "=" * 78)
    print("z*(R) trajectory  (prediction: z* INCREASES as R shrinks = contraction)")
    print("=" * 78)
    for i, R in enumerate(R_GRID):
        print(f"  R={R:.3f}  z*={z_star[i]:.3f}  E_fixed={E_fixed[i]:.5f}  "
              f"E_relax={E_relax[i]:.5f}  gain={1e3*(E_fixed[i]-E_relax[i]):.2f} mHa")

    print("\n" + "=" * 78)
    print("R_eq / tilt by relaxation level (n_max=2)")
    print("=" * 78)
    summary = []
    for name, E in [('fixed Z_orb=1 (baseline)', E_fixed),
                    ('R-adaptive shared exponent', E_relax)]:
        req = _quartic_req(R_GRID, E)
        tilt, curv = _tilt(R_GRID, E)
        err = abs(req - R_TRUE) / R_TRUE * 100 if np.isfinite(req) else float('nan')
        print(f"  {name:28s}  R_eq={req:.4f} ({err:5.1f}%)  tilt={tilt:+.4f}  curv={curv:+.4f}")
        summary.append({'variant': name, 'R_eq': req, 'R_eq_err_pct': err,
                        'tilt': tilt, 'curvature': curv})

    # verdict
    base_err = summary[0]['R_eq_err_pct']; relax_err = summary[1]['R_eq_err_pct']
    print("\nVERDICT:")
    if np.isfinite(relax_err) and relax_err < 0.4 * base_err:
        print(f"  MECHANISM CONFIRMED (radial): R-adaptive contraction collapses the "
              f"drift {base_err:.1f}% -> {relax_err:.1f}%.")
    elif np.isfinite(relax_err) and relax_err < 0.8 * base_err:
        print(f"  PARTIAL: contraction helps ({base_err:.1f}% -> {relax_err:.1f}%) but "
              f"does not close it -> residual is likely angular/p-polarization.")
    else:
        print(f"  NOT radial contraction: drift {base_err:.1f}% -> {relax_err:.1f}% "
              f"(little change) -> mechanism is NOT a shared radial exponent.")

    out = {
        'R_grid': R_GRID.tolist(), 'z_grid': Z_GRID.tolist(),
        'E_table': Etab.tolist(),
        'E_fixed': E_fixed.tolist(), 'E_relax': E_relax.tolist(),
        'z_star': z_star.tolist(),
        'summary': summary, 'required_force': Z_LI * Z_H / R_TRUE**2,
        'params': {'MAX_N': MAX_N, 'N_GRID_VNE': N_GRID_VNE, 'L_MAX': L_MAX,
                   'E_CORE': E_CORE},
    }
    Path('debug/data').mkdir(parents=True, exist_ok=True)
    with open('debug/data/balanced_reqdrift_relaxation.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] debug/data/balanced_reqdrift_relaxation.json "
          f"({time.perf_counter()-t0:.0f}s)")


if __name__ == '__main__':
    main()
