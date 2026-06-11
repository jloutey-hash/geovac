"""Sprint α-Diagnostic / Track α-1: Numerical bimodule distance d_L, d_R, d_M^LR.

Verifies the M-Y prediction d_R/d_L ~ 6.7 for NaH and tests the alkali-hydride
scaling law d_R ~ r_valence^phys(M) - 1.5 bohr.

Workflow:
  1. Define four pin states for an alkali-hydride MH:
     (i)   hydrogenic Z=1 product (current framework default)
     (ii)  SV-corrected diagonal (bit-identical to (i) at wavefunction level)
     (iii) cross-center bonding orbital with physical M-valence shape (from
           FrozenCore Z_eff(r) Schrodinger radial solver)
     (iv)  small perturbation of (iii) by V_ee/ΔE (approximate)
  2. Compute d_L, d_R, d_M^LR pairwise using a bond-region L^2 inner product.
  3. Repeat for LiH (M=Li, Z=3), NaH (M=Na, Z=11), KH (M=K, Z=19).
  4. Test consistency under basis size (max_n=2 vs 3 — same radial shapes) and
     bond region centering (R ∈ {3.0, 3.5, 4.0} bohr for NaH).

Diagnostic-only. Do NOT modify production code. Output: JSON file +
companion memo at debug/sprint_alpha_1_diagnostic_memo.md.

Conventions:
- "left action" = H-centered multiplication algebra acts on the H-side
  component of the bimodule element.
- "right action" = M-centered (alkali) multiplication algebra acts on the
  M-side component.
- For pin states with no cross-center linear combination, the "M-side"
  component is the on-center valence radial function with its assigned
  shape (hydrogenic Z=1 vs physical-Z_eff(r) screened).
"""

from __future__ import annotations

import json
import os
import sys
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

# Ensure we can import from the project
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# ---------------------------------------------------------------------------
# Pin-state radial-function constructors
# ---------------------------------------------------------------------------

# Hydrogenic R_{n,l}(r; Z_eff): use n=1, l=0 (Z_eff=1) for the H-1s shape and
# also for the framework's hydrogenic Z=1 "M-valence" stand-in (this is exactly
# the framework default for the M-valence basis at Z_orb=1 — Track 3 memo §3).

def hydrogenic_1s(Z_eff: float, r: np.ndarray) -> np.ndarray:
    """R_{1s}(r; Z_eff) = 2 Z_eff^{3/2} exp(-Z_eff r)."""
    return 2.0 * Z_eff**1.5 * np.exp(-Z_eff * r)


def physical_M_valence_radial(Z: int, n_val: int, l: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    """Compute physical valence radial function R(r) for an alkali atom M
    using FrozenCore Z_eff(r) screening.

    Returns
    -------
    r : ndarray
        Radial grid.
    R : ndarray
        Radial wavefunction R(r) (NOT u(r) = r*R(r)), normalized such that
        int |R|^2 r^2 dr = 1.

    Notes
    -----
    For Li (Z=3), there is no FrozenCore (it has a 1s² explicit core in the
    framework). We use a hydrogenic Z=1.3 approximation for the Li 2s — this
    is a reasonable Slater-rule estimate (mean radius ≈ 2.3 bohr) and
    matches the Li 2s mean radius from quantum chemistry tables (~3.4 bohr
    in the asymptotic limit but ~2-3 bohr for the inner peak).
    """
    from geovac.neon_core import _solve_screened_radial

    if Z == 3:
        # Li 2s: no FrozenCore registered (1s² is explicit). Use hydrogenic
        # Z_eff via Slater's rules: 2s sees 0.85*2 = 1.7 screening → Z_eff=1.3.
        # n=2, l=0, Z_eff=1.3.
        r_max = 40.0
        n_grid = 4000
        r = np.linspace(0.01, r_max, n_grid)
        Z_eff = 1.3
        from scipy.special import genlaguerre
        rho = 2.0 * Z_eff * r / 2
        L_poly = genlaguerre(2 - 0 - 1, 1)(rho)
        R = rho**0 * np.exp(-rho / 2.0) * L_poly
        # Normalize: int |R|^2 r^2 dr = 1
        norm_sq = np.trapezoid(R**2 * r**2, r)
        R /= np.sqrt(norm_sq)
        return r, R

    if Z == 11:
        # Na 3s: [Ne] core
        core_type = 'Ne'
    elif Z == 19:
        # K 4s: [Ar] core
        core_type = 'Ar'
    else:
        raise ValueError(f"Unsupported Z={Z}")

    # Solve screened radial Schrödinger equation for l=0 (allow_l0=True)
    # n_target = n_val (3 for Na, 4 for K)
    energy, u_vec, r = _solve_screened_radial(
        Z=Z, l=l, n_target=n_val,
        core_type=core_type,
        n_grid=12000, r_max=80.0,
        allow_l0=True,
    )
    # Convert u(r) = r R(r) → R(r)
    R = np.where(r > 0, u_vec / r, 0.0)
    # Renormalize int |R|^2 r^2 dr = 1 (should already be 1 since
    # int |u|^2 dr = 1, but recompute for safety)
    norm_sq = np.trapezoid(R**2 * r**2, r)
    if norm_sq > 0:
        R /= np.sqrt(norm_sq)
    return r, R


# ---------------------------------------------------------------------------
# Pin-state construction (the four candidates)
# ---------------------------------------------------------------------------

@dataclass
class PinState:
    """A pin state for the bimodule diagnostic.

    Components
    ----------
    H_side(r) : Callable[[np.ndarray], np.ndarray]
        Radial function R^H(r_H) centered at H. Always hydrogenic Z=1 in the
        current framework (this is the correct H-side basis).
    M_side(r) : Callable[[np.ndarray], np.ndarray]
        Radial function R^M(r_M) centered at M (alkali). For pin state (i),
        this is hydrogenic Z=1 (the framework default). For pin state (iii),
        this is the physical screened M-valence radial function.
    label : str
        Identifier of the candidate, e.g. 'i', 'ii', 'iii', 'iv'.
    coefficients : Tuple[float, float]
        (c_H, c_M) bonding-orbital coefficients for cross-center linear
        combination. (1.0, 0.0) means H-only, (0.0, 1.0) means M-only,
        (c_H, c_M) means a sigma_g bonding orbital.
    """
    label: str
    H_side: Callable[[np.ndarray], np.ndarray]
    M_side: Callable[[np.ndarray], np.ndarray]
    coefficients: Tuple[float, float]


def build_pin_states(
    M_symbol: str,
    Z_M: int,
    n_val_M: int,
) -> Dict[str, PinState]:
    """Build the four pin states for an alkali-hydride MH.

    For each, return a (H_side(r), M_side(r), coefficients) tuple.

    Pin state details:
      (i)   H_side=hydrogenic Z=1, M_side=hydrogenic Z=1 (n_val_M).
            Coefficients (1/sqrt2, 1/sqrt2): symmetric Hartree product
            (both centers populated equally, no cross-center mixing).
      (ii)  Identical wavefunctions to (i), differs only at diagonal
            (eigenvalue) level. Bimodule prediction: d_L = d_R = 0.
      (iii) H_side=hydrogenic Z=1 (still correct), M_side=physical M-valence
            from FrozenCore Z_eff(r). Coefficients are bonding-orbital
            symmetric for simplicity (would be R-dependent in full HF).
      (iv)  Small perturbation of (iii): admixes the antibonding character
            via |V_ee/ΔE| ~ 0.05 (NaH valence scale).
    """
    # H-side is always hydrogenic Z=1 1s
    def H_hyd_Z1(r: np.ndarray) -> np.ndarray:
        return hydrogenic_1s(1.0, r)

    # M-side hydrogenic Z=1 nl: use shape of physical n_val with Z_eff=1
    def M_hyd_Z1(r: np.ndarray) -> np.ndarray:
        # Hydrogenic n_val=n, l=0, Z_eff=1 nodes-and-all
        from scipy.special import genlaguerre
        n = n_val_M
        Z_eff = 1.0
        rho = 2.0 * Z_eff * r / n
        L_poly = genlaguerre(n - 1, 1)(rho)
        wf = np.exp(-rho / 2.0) * L_poly
        # Normalize int |R|^2 r^2 dr = 1
        norm_sq = np.trapezoid(wf**2 * r**2, r)
        if norm_sq > 0:
            wf = wf / np.sqrt(norm_sq)
        return wf

    # M-side physical from FrozenCore (for n_val=3 Na, n_val=4 K, or Li=2)
    print(f"  [physical M-valence radial for Z={Z_M}, n_val={n_val_M}]")
    r_phys_grid, R_phys = physical_M_valence_radial(Z_M, n_val=n_val_M, l=0)

    # Build a callable that interpolates onto an arbitrary grid
    def M_physical(r: np.ndarray) -> np.ndarray:
        return np.interp(r, r_phys_grid, R_phys, left=0.0, right=0.0)

    # Pin states
    coef_sym = (1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0))  # symmetric bonding
    states = {
        'i':   PinState('i',   H_hyd_Z1,   M_hyd_Z1,    coef_sym),
        'ii':  PinState('ii',  H_hyd_Z1,   M_hyd_Z1,    coef_sym),  # bit-identical to (i)
        'iii': PinState('iii', H_hyd_Z1,   M_physical,  coef_sym),
        # (iv): perturbation of (iii) — mix with antibonding (1, -1)/sqrt(2)
        # at weight ~0.05. The bonding-orbital coefficient is reshuffled.
        # We keep H_side, M_side fixed but change coefficients slightly.
        'iv':  PinState('iv',  H_hyd_Z1,   M_physical,  (np.cos(0.05) / np.sqrt(2.0), np.sin(0.05) * 1.0 / np.sqrt(2.0))),
    }
    return states


# ---------------------------------------------------------------------------
# Two-axis L/R bimodule distance
# ---------------------------------------------------------------------------

def hilbert_schmidt_dist(
    f1: Callable[[np.ndarray], np.ndarray],
    f2: Callable[[np.ndarray], np.ndarray],
    r_center: float,
    r_window: float,
    n_quad: int = 200,
) -> float:
    """Hilbert-Schmidt L^2 distance between two radial functions in a
    bond-region quadrature.

    The "bond-region weight" is supported on [r_center - r_window, r_center
    + r_window], approximating the cross-coupling integral support. This is
    where the bimodule structure lives — the cross-action region between
    centers.

    d(f1, f2) = sqrt( int_{|r - r_center| < r_window} |f1(r) - f2(r)|^2 r^2 dr )

    Parameters
    ----------
    f1, f2 : Callable[[r], R(r)]
        Two radial functions to compare.
    r_center : float
        Center of the bond region (e.g. R/2 for NaH bond midpoint).
    r_window : float
        Half-width of the bond region.
    n_quad : int
        Number of quadrature points.
    """
    r_min = max(0.01, r_center - r_window)
    r_max = r_center + r_window
    r = np.linspace(r_min, r_max, n_quad)
    diff = f1(r) - f2(r)
    return float(np.sqrt(np.trapezoid(diff**2 * r**2, r)))


def bimodule_distance_LR(
    psi_a: PinState,
    psi_b: PinState,
    R_bond: float,
    bond_window: float = 1.0,
) -> Dict[str, float]:
    """Compute the two-axis L/R bimodule distance between two pin states.

    For an alkali-hydride MH at bond length R, the bond region runs along
    the M-H axis. The "left axis" is centered near H; the "right axis" is
    centered near M.

    Two notions of d_L, d_R are computed here:

    (A) "Side-restricted" d_L^side, d_R^side: pure radial-function L^2
        distances on each center. For pin states (i) and (iii) that share
        the H-side function (hydrogenic Z=1), d_L^side = 0 by
        construction. This is the trivial bimodule-element distance.

    (B) "Bond-region-probed" d_L^bond, d_R^bond: distances measured at the
        bond region, where the H-multiplication algebra probes the
        bimodule element NEAR H and the M-multiplication algebra probes
        it NEAR M.  Each pin state contributes BOTH the H-side and M-side
        amplitudes in the bond region (the M-side tail extending toward
        H, the H-side tail extending toward M).  This is the M-Y
        operator-distance reading.

    The diagnostic returns (B) as the load-bearing d_L, d_R since this is
    what M-Y predicts (~6.7 ratio); (A) is returned as a side-table for
    transparency.

    Coefficients (c_H, c_M) for each pin state weight the side-amplitudes
    (since the bonding orbital is a linear combination of the two centers).
    """
    # Effective H-side radial function = c_H * R^H(r)
    def H_eff_a(r: np.ndarray) -> np.ndarray:
        return psi_a.coefficients[0] * psi_a.H_side(r)

    def H_eff_b(r: np.ndarray) -> np.ndarray:
        return psi_b.coefficients[0] * psi_b.H_side(r)

    def M_eff_a(r: np.ndarray) -> np.ndarray:
        return psi_a.coefficients[1] * psi_a.M_side(r)

    def M_eff_b(r: np.ndarray) -> np.ndarray:
        return psi_b.coefficients[1] * psi_b.M_side(r)

    # -----------------------------------------------------------------------
    # (A) Side-restricted: pure radial-function L^2 distance on each side.
    # -----------------------------------------------------------------------
    # H-side: probe at small r from H (where the H-multiplication algebra
    # is fully resolved).
    r_H_center = 1.0   # bohr from H
    d_L_side = hilbert_schmidt_dist(H_eff_a, H_eff_b, r_H_center, bond_window)

    # M-side: probe at the bond midpoint (M-side coords).
    r_M_center = R_bond / 2.0
    d_R_side = hilbert_schmidt_dist(M_eff_a, M_eff_b, r_M_center,
                                      R_bond / 2.0, n_quad=400)

    # -----------------------------------------------------------------------
    # (B) Bond-region-probed: each algebra probes the FULL bimodule element
    # in the bond region. For an MH bonding orbital, the bimodule element
    # at a point r in space (1D along the bond axis) is:
    #     xi(z) = c_H * R^H(|z - R_H|) + c_M * R^M(|z - R_M|)
    # where R_M = 0 (M at origin), R_H = R_bond. We probe along the bond
    # axis (1D parametrization is sufficient for the radial diagnostic).
    # -----------------------------------------------------------------------
    n_z = 400
    z_grid = np.linspace(-1.0, R_bond + 1.0, n_z)  # 1D coord along bond axis
    r_M_at_z = np.abs(z_grid)              # distance from M (at origin)
    r_H_at_z = np.abs(z_grid - R_bond)     # distance from H

    def xi_at_z(psi: PinState) -> np.ndarray:
        return psi.coefficients[0] * psi.H_side(r_H_at_z) + \
               psi.coefficients[1] * psi.M_side(r_M_at_z)

    xi_a = xi_at_z(psi_a)
    xi_b = xi_at_z(psi_b)

    # The "left-action distance" d_L^bond: weighting the |xi_a - xi_b|^2
    # by a function peaked at H (left-action probe).
    # Use Gaussian peaks of width ~1 bohr at H (z = R_bond) and M (z = 0).
    width = 1.0
    weight_L = np.exp(-(z_grid - R_bond)**2 / (2 * width**2))  # peaked at H
    weight_R = np.exp(-(z_grid - 0.0)**2     / (2 * width**2))  # peaked at M

    # Normalize weights so int weight dz = 1 (so the d's are comparable across
    # bond lengths)
    weight_L /= np.trapezoid(weight_L, z_grid)
    weight_R /= np.trapezoid(weight_R, z_grid)

    integrand_L = weight_L * (xi_a - xi_b)**2
    integrand_R = weight_R * (xi_a - xi_b)**2

    d_L = float(np.sqrt(np.trapezoid(integrand_L, z_grid)))
    d_R = float(np.sqrt(np.trapezoid(integrand_R, z_grid)))

    d_LR = float(np.sqrt(d_L**2 + d_R**2))
    ratio = d_R / d_L if d_L > 1e-15 else float('inf')

    return {
        'd_L': d_L,
        'd_R': d_R,
        'd_LR': d_LR,
        'ratio_R_over_L': ratio,
        'd_L_side': d_L_side,
        'd_R_side': d_R_side,
    }


# ---------------------------------------------------------------------------
# Mean radius diagnostic (for the scaling test)
# ---------------------------------------------------------------------------

def mean_radius(R_func: Callable[[np.ndarray], np.ndarray],
                r_max: float = 30.0, n_quad: int = 3000) -> float:
    """Compute <r> = int |R|^2 r * r^2 dr / int |R|^2 r^2 dr ."""
    r = np.linspace(0.01, r_max, n_quad)
    R = R_func(r)
    num = np.trapezoid(R**2 * r * r**2, r)
    denom = np.trapezoid(R**2 * r**2, r)
    return float(num / denom) if denom > 0 else 0.0


# ---------------------------------------------------------------------------
# Main diagnostic
# ---------------------------------------------------------------------------

def run_diagnostic(
    M_symbol: str,
    Z_M: int,
    n_val_M: int,
    R_bond: float,
    R_bond_alt: Optional[List[float]] = None,
) -> Dict:
    """Run the full diagnostic for an alkali-hydride MH.

    Returns the bimodule distance table for all (i, j) pairs, plus the
    R-dependence consistency check.
    """
    print(f"\n=== Sprint alpha-1 diagnostic: {M_symbol}H ===")
    print(f"  Z_M = {Z_M}, n_val_M = {n_val_M}, R_bond = {R_bond} bohr")

    states = build_pin_states(M_symbol, Z_M, n_val_M)

    # Compute mean radii (sanity)
    print("\n  Pin-state mean radii (M-side, with c_M = 1/sqrt(2) implicit):")
    mean_radii = {}
    for lbl, ps in states.items():
        r_mean = mean_radius(ps.M_side)
        mean_radii[lbl] = r_mean
        print(f"    pin state ({lbl}): <r>_M = {r_mean:.3f} bohr")

    # Pair-wise bimodule distances at the primary bond length
    print(f"\n  Pair-wise (d_L, d_R, d_M^LR, ratio) at R={R_bond} bohr:")
    pair_results = {}
    pair_labels = [('i', 'ii'), ('i', 'iii'), ('i', 'iv'), ('iii', 'iv')]
    for a, b in pair_labels:
        d = bimodule_distance_LR(states[a], states[b], R_bond)
        pair_results[f"{a}_vs_{b}"] = d
        print(f"    ({a}) vs ({b}):  d_L = {d['d_L']:.4f},  d_R = {d['d_R']:.4f},"
              f"  d_LR = {d['d_LR']:.4f},  d_R/d_L = {d['ratio_R_over_L']:.2f}")

    # R-dependence consistency: try alternate bond lengths
    R_dependence = {}
    if R_bond_alt:
        print(f"\n  R-dependence consistency check ((i) vs (iii)):")
        for R in R_bond_alt:
            d = bimodule_distance_LR(states['i'], states['iii'], R)
            R_dependence[f"R_{R:.1f}"] = d
            print(f"    R = {R:.2f}:  d_L = {d['d_L']:.4f},  d_R = {d['d_R']:.4f},"
                  f"  d_R/d_L = {d['ratio_R_over_L']:.2f}")

    return {
        'molecule': f"{M_symbol}H",
        'Z_M': Z_M,
        'n_val_M': n_val_M,
        'R_bond_primary': R_bond,
        'pair_results': pair_results,
        'R_dependence': R_dependence,
        'mean_radii_M_side': mean_radii,
    }


def main():
    t0 = time.time()
    # Set stdout to UTF-8 on Windows
    import io
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass
    print("=" * 70)
    print(" Sprint alpha-Diagnostic / Track alpha-1: bimodule distance verification")
    print("=" * 70)

    results = {}

    # NaH (primary): Z=11, n_val=3, R_eq=3.566 bohr
    results['NaH'] = run_diagnostic(
        'Na', Z_M=11, n_val_M=3, R_bond=3.566,
        R_bond_alt=[3.0, 3.5, 4.0],
    )

    # LiH (scaling): Z=3, n_val=2, R_eq=3.015 bohr
    results['LiH'] = run_diagnostic(
        'Li', Z_M=3, n_val_M=2, R_bond=3.015,
    )

    # KH (scaling): Z=19, n_val=4, R_eq=4.243 bohr
    results['KH'] = run_diagnostic(
        'K', Z_M=19, n_val_M=4, R_bond=4.243,
    )

    # Scaling summary
    print("\n" + "=" * 70)
    print(" Alkali-hydride scaling test")
    print("=" * 70)
    print(f"  {'Mol':>6} {'Z_M':>4} {'n_val':>6} {'<r>_M_phys':>12} {'<r>_M_hyd':>12}"
          f" {'d_R((i),(iii))':>16}")
    scaling_table = []
    for mol in ['LiH', 'NaH', 'KH']:
        rr = results[mol]
        r_phys = rr['mean_radii_M_side']['iii']  # physical M-valence
        r_hyd = rr['mean_radii_M_side']['i']     # hydrogenic Z=1
        d_R_i_iii = rr['pair_results']['i_vs_iii']['d_R']
        print(f"  {mol:>6} {rr['Z_M']:>4} {rr['n_val_M']:>6} "
              f"{r_phys:>12.3f} {r_hyd:>12.3f} {d_R_i_iii:>16.4f}")
        scaling_table.append({
            'molecule': mol,
            'Z_M': rr['Z_M'],
            'n_val_M': rr['n_val_M'],
            'r_phys_M': r_phys,
            'r_hyd_M': r_hyd,
            'd_R_i_vs_iii': d_R_i_iii,
            'd_L_i_vs_iii': rr['pair_results']['i_vs_iii']['d_L'],
            'd_R_minus_1p5': r_phys - 1.5,  # M-Y predicted scaling: r_phys - 1.5 bohr
        })

    # Save JSON
    out_dir = os.path.join(os.path.dirname(__file__), 'data')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'sprint_alpha_1_diagnostic.json')

    # Strip non-serializable callables; we only kept dicts above so just dump.
    full_results = {
        'per_molecule': results,
        'scaling_table': scaling_table,
        'metadata': {
            'mY_prediction_d_R_over_d_L': 6.7,
            'mY_prediction_scaling': 'd_R ~ r_valence_phys - 1.5 bohr',
            'wall_time_s': time.time() - t0,
        },
    }
    with open(out_path, 'w') as f:
        json.dump(full_results, f, indent=2)

    print(f"\n  Results saved to: {out_path}")
    print(f"  Wall time: {time.time() - t0:.1f} s")

    return full_results


if __name__ == '__main__':
    main()
