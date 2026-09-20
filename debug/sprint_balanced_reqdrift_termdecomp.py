"""
Sprint: balanced-coupled LiH R_eq-drift term localization (DIAGNOSTIC ONLY).
============================================================================
Question: the balanced-coupled LiH PES reaches ~0.20% ENERGY accuracy (n_max=3)
but ~9% R_eq error -- energy converges, geometry drifts. WHICH Hamiltonian term's
R-dependence drives the drift?

Prior work (debug/sprint_chem_error_projection_memo.md, 2026-07-05): established
that R_eq error = residual TILT eps'(R_true) / computed curvature, the tilt is
frozen ~-0.030 Ha/bohr while energy converges, and it is OUTWARD (favours long
bonds). The well-shape test (tests/test_paper19_well_shape.py) further split the
tilt into dV_NN/dR (exact, = -3/R^2) + electronic gradient, and found the
electronic gradient ~8.8% too weak. NEITHER decomposed the electronic gradient
BY TERM. That is this script's job.

Method (rigorous force decomposition via term-freezing):
  E(R) = min_psi <psi| H(R) |psi>. For any term T, define H_frozenT(R) by
  replacing T(R) with T(R_ref). Then by Hellmann-Feynman at R_ref (where
  H_frozenT(R_ref)=H(R_ref) so psi is identical):
        dE_frozenT/dR |_{R_ref} = <psi| d(H-T)/dR |psi> = tilt_full - F_T,
  where F_T = <psi| dT/dR |psi> is term T's exact contribution to the force.
  => F_T = tilt_full - tilt_frozenT,  and  sum_T F_T = tilt_full (freeze all -> 0).

  We also locate R_eq of each frozen variant to see which term's R-dependence,
  when removed, restores the true geometry.

R-dependent pieces of the balanced Hamiltonian:
  - V_NN            : nuclear_repulsion = E_core(const) + Z_Li*Z_H/R = E_core + 3/R
  - cross-center V_ne: h1_cross_vne(R)  (Shibuya-Wulfman, hydrogenic Z_orb basis)
  - cross-block ERI : eri_balanced(R) - eri_within(R)
  - within-block h1 : h1_no_pk          (expected R-INDEPENDENT -> verified)
  - within-block ERI: eri_within(R)     (expected R-INDEPENDENT -> verified)
"""
from __future__ import annotations
import json
import time
import numpy as np
from pathlib import Path

from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.composed_qubit import build_composed_hamiltonian
from geovac.molecular_spec import lih_spec, _FIRST_ROW_CORE_ENERGY

R_TRUE = 3.015
E_EXACT = -8.071          # Ha exact LiH total (Paper 19)
Z_LI, Z_H = 3.0, 1.0
E_CORE = _FIRST_ROW_CORE_ENERGY[3]   # -7.2799, He-like Li(2+) core, R-independent
MAX_N = 2

# Grid: wide enough to bracket the (outward-drifted) minima of the frozen variants
R_GRID = np.array([2.70, 2.85, 3.015, 3.15, 3.30, 3.45, 3.60, 3.80])
H_DERIV = 3.015  # central point for force decomposition uses the grid neighbours


def build_components(R: float):
    """Return the R-dependent component matrices at bond length R."""
    spec = lih_spec(R=R, max_n=MAX_N)
    n_e = sum(b.n_electrons for b in spec.blocks)
    ham = build_balanced_hamiltonian(spec, R=R, n_grid_vne=8000, L_max=4)
    comp = build_composed_hamiltonian(spec, pk_in_hamiltonian=False, verbose=False)
    eri_within = comp['eri']
    M = ham['M']
    vnn = ham['nuclear_repulsion']          # = E_core + 3/R
    return {
        'R': R, 'M': M, 'n_e': n_e,
        'h1_no_pk': ham['h1_no_pk'].copy(),
        'h1_cross_vne': ham['h1_cross_vne'].copy(),
        'eri_balanced': ham['eri'].copy(),
        'eri_within': eri_within.copy(),
        'nuclear_repulsion': vnn,
        'vnn_only': vnn - E_CORE,           # 3/R part
    }


def fci_E(M, n_e, h1, eri, vnn):
    res = {'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': vnn}
    return float(coupled_fci_energy(res, n_electrons=n_e, verbose=False)['E_coupled'])


def main():
    t0 = time.perf_counter()
    print("Building component matrices on the R grid ...")
    C = {}
    for R in R_GRID:
        C[R] = build_components(R)
        print(f"  R={R:.3f} built ({time.perf_counter()-t0:.0f}s)")

    # --- verify which within-block pieces are R-independent -------------------
    R0 = R_TRUE
    diffs = {}
    for R in R_GRID:
        if abs(R - R0) < 1e-9:
            continue
        diffs[R] = {
            'd_h1_no_pk': float(np.max(np.abs(C[R]['h1_no_pk'] - C[R0]['h1_no_pk']))),
            'd_eri_within': float(np.max(np.abs(C[R]['eri_within'] - C[R0]['eri_within']))),
            'd_h1_cross_vne': float(np.max(np.abs(C[R]['h1_cross_vne'] - C[R0]['h1_cross_vne']))),
            'd_eri_cross_block': float(np.max(np.abs(
                (C[R]['eri_balanced'] - C[R]['eri_within'])
                - (C[R0]['eri_balanced'] - C[R0]['eri_within'])))),
            'd_vnn': float(abs(C[R]['vnn_only'] - C[R0]['vnn_only'])),
        }
    print("\nR-dependence probe (max|Delta| vs R_true):")
    for R, d in diffs.items():
        print(f"  R={R:.3f}: h1_no_pk={d['d_h1_no_pk']:.2e} eri_within={d['d_eri_within']:.2e} "
              f"cross_vne={d['d_h1_cross_vne']:.2e} cross_eri={d['d_eri_cross_block']:.2e} "
              f"vnn={d['d_vnn']:.2e}")

    # --- variant Hamiltonian assemblers (freeze one term's R-dep at R_ref) ----
    Rref = R_TRUE

    def assemble(R, freeze):
        c = C[R]; cr = C[Rref]
        M, n_e = c['M'], c['n_e']
        h1 = c['h1_no_pk'] + c['h1_cross_vne']
        eri = c['eri_balanced'].copy()
        vnn = c['nuclear_repulsion']
        cross_block_R = c['eri_balanced'] - c['eri_within']
        cross_block_ref = cr['eri_balanced'] - cr['eri_within']
        if 'vnn' in freeze:
            vnn = E_CORE + cr['vnn_only']
        if 'cross_vne' in freeze:
            h1 = c['h1_no_pk'] + cr['h1_cross_vne']
        if 'cross_eri' in freeze:
            eri = c['eri_within'] + cross_block_ref
        if 'eri_within' in freeze:
            eri = cr['eri_within'] + cross_block_R
        if 'h1_no_pk' in freeze:
            h1 = cr['h1_no_pk'] + c['h1_cross_vne']
        return M, n_e, h1, eri, vnn

    variants = {
        'FULL': [],
        'freeze_vnn': ['vnn'],
        'freeze_cross_vne': ['cross_vne'],
        'freeze_cross_eri': ['cross_eri'],
        'freeze_eri_within': ['eri_within'],
        'freeze_h1_no_pk': ['h1_no_pk'],
        'freeze_ALL_electronic': ['cross_vne', 'cross_eri', 'eri_within', 'h1_no_pk'],
        'freeze_ALL': ['vnn', 'cross_vne', 'cross_eri', 'eri_within', 'h1_no_pk'],
    }

    print("\nComputing PES for each variant ...")
    PES = {}
    for name, fr in variants.items():
        Es = []
        for R in R_GRID:
            Es.append(fci_E(*assemble(R, fr)))
        PES[name] = np.array(Es)
        print(f"  {name:24s} done ({time.perf_counter()-t0:.0f}s)")

    # --- analyse: R_eq and tilt at R_true for each variant -------------------
    def local_fit(R, E, center, window=0.5):
        m = np.abs(R - center) <= window
        x = R[m] - center
        order = min(4, m.sum() - 1)
        p = np.poly1d(np.polyfit(x, E[m], order))
        return p

    def r_eq_of(E, window_lo=2.7, window_hi=3.8):
        # global quartic over the whole grid, pick interior minimum
        p = np.poly1d(np.polyfit(R_GRID - R_TRUE, E, 4))
        dp = p.deriv(1); ddp = p.deriv(2)
        roots = dp.r[np.isreal(dp.r)].real
        cand = [r + R_TRUE for r in roots
                if ddp(r) > 0 and window_lo < r + R_TRUE < window_hi]
        if not cand:
            return float('nan')
        return float(min(cand, key=lambda rr: p(rr - R_TRUE)))

    print("\n" + "=" * 78)
    print("RESULT 1 -- R_eq of each variant (which frozen term restores geometry?)")
    print("=" * 78)
    req = {}
    for name in variants:
        req[name] = r_eq_of(PES[name])
        err = abs(req[name] - R_TRUE) / R_TRUE * 100 if np.isfinite(req[name]) else float('nan')
        print(f"  {name:24s} R_eq = {req[name]:.4f} bohr  (err {err:5.1f}%)")

    # tilt at R_true (central difference from the two grid neighbours of R_true)
    i0 = int(np.argmin(np.abs(R_GRID - R_TRUE)))
    hL = R_TRUE - R_GRID[i0 - 1]
    hR = R_GRID[i0 + 1] - R_TRUE
    def tilt_at_true(E):
        # asymmetric 3-point derivative at R_true using neighbours i0-1,i0,i0+1
        xL, x0, xR = R_GRID[i0 - 1], R_GRID[i0], R_GRID[i0 + 1]
        # fit a parabola through the 3 nearest points and evaluate slope at R_true
        p = np.poly1d(np.polyfit([xL, x0, xR], [E[i0 - 1], E[i0], E[i0 + 1]], 2))
        return float(p.deriv(1)(R_TRUE))

    print("\n" + "=" * 78)
    print("RESULT 2 -- force decomposition at R_true: F_T = tilt_full - tilt_frozenT")
    print("=" * 78)
    tilt_full = tilt_at_true(PES['FULL'])
    print(f"  tilt_full = dE/dR|R_true = {tilt_full:+.4f} Ha/bohr  "
          f"(<0 => minimum pushed OUTWARD)")
    print(f"  required electronic force to cancel dV_NN/dR: +Z_LiZ_H/R^2 = "
          f"{Z_LI*Z_H/R_TRUE**2:+.4f} Ha/bohr")
    F = {}
    for name, fr in variants.items():
        if name == 'FULL':
            continue
        tf = tilt_at_true(PES[name])
        F[name] = tilt_full - tf
        print(f"  F[{name:22s}] = {F[name]:+.4f} Ha/bohr   (tilt_frozen={tf:+.4f})")
    # consistency: sum of single-term freezes vs freeze_ALL
    single = ['freeze_vnn', 'freeze_cross_vne', 'freeze_cross_eri',
              'freeze_eri_within', 'freeze_h1_no_pk']
    print(f"\n  sum_T F_T (single freezes) = {sum(F[s] for s in single):+.4f}")
    print(f"  F[freeze_ALL]              = {F['freeze_ALL']:+.4f}  (should ~= tilt_full)")
    print(f"  tilt_full                  = {tilt_full:+.4f}")

    # --- energy accuracy at R_true and at full R_eq --------------------------
    E_full_true = PES['FULL'][i0] - E_CORE
    print("\n" + "=" * 78)
    print("RESULT 3 -- energy vs geometry split (n_max=2)")
    print("=" * 78)
    print(f"  E(R_true) [TC] = {E_full_true:.4f} Ha  "
          f"(exact {E_EXACT}, err {abs(E_full_true-E_EXACT)/abs(E_EXACT)*100:.2f}%)")
    print(f"  R_eq(FULL)     = {req['FULL']:.4f} bohr "
          f"(exact {R_TRUE}, err {abs(req['FULL']-R_TRUE)/R_TRUE*100:.1f}%)")

    out = {
        'R_grid': R_GRID.tolist(),
        'PES': {k: v.tolist() for k, v in PES.items()},
        'R_eq': req,
        'tilt_full': tilt_full,
        'forces_F_T': F,
        'required_elec_force': Z_LI * Z_H / R_TRUE**2,
        'E_true_TC': E_full_true,
        'E_exact': E_EXACT,
        'rdep_probe': diffs,
    }
    Path('debug/data').mkdir(parents=True, exist_ok=True)
    with open('debug/data/balanced_reqdrift_termdecomp.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n[saved] debug/data/balanced_reqdrift_termdecomp.json "
          f"(total {time.perf_counter()-t0:.0f}s)")


if __name__ == '__main__':
    main()
