"""Sprint F4 Step 1 — algebraic bonding-vs-core overlap diagnostic.

Algebraic-first gate before any production wiring:

At NaH R = R_eq = 3.566 bohr with the full F3 stack (W1c screening +
multi-zeta Na 3s + cross-block h1), compute the bonding orbital
from cross-block h1 diagonalization in the {|Na 3s⟩, |H 1s⟩} subspace.
Compute its overlap with each Na [Ne] core orbital (1s, 2s, 2p).
Estimate the PK barrier magnitude:

    ΔH^PK_{bond,bond} = sum_c (E_v - E_c) |S_{bond, c}|^2

Gate:
- |overlap with at least one core| > 5% AND predicted barrier > 1 Ha
  → proceed to Step 2 (W1e is real, PK should close it)
- overlap < 1% → STOP, W1e isn't a Pauli problem (deeper mechanism)
- intermediate → proceed but flag weak expected impact
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from typing import Dict, Tuple

import numpy as np
from scipy.linalg import eigh

# Project root on sys.path
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from geovac.cross_block_h1 import (  # noqa: E402
    hydrogenic_R_nl_analytical,
    overlap_ss_axial,
)
from geovac.multi_zeta_orbitals import (  # noqa: E402
    MultiZetaOrbital,
    STO,
    get_physical_valence_orbitals,
)


# ----------------------------------------------------------------------
# Na [Ne] core orbital constructors
# Sourced consistently with geovac/phillips_kleinman_cross_center.py:
# Clementi-Raimondi single-zeta exponents from neon_core._CLEMENTI_ZETA_NE
# and Clementi-Roetti HF eigenvalues from _CORE_ORBITAL_ENERGIES_NE.
# ----------------------------------------------------------------------

# Na (Z=11) Clementi-Raimondi exponents (from neon_core._CLEMENTI_ZETA_NE[11])
NA_CR_EXPONENTS = {
    (1, 0): 10.62608,  # 1s (zeta_1s)
    (2, 0): 3.21766,   # 2s
    (2, 1): 3.45200,   # 2p
}

# Na [Ne]-core HF eigenvalues (Hartree) from Clementi-Roetti 1974
# (matches phillips_kleinman_cross_center._CORE_ORBITAL_ENERGIES_NE[11])
NA_CORE_ENERGIES = {
    (1, 0): -40.4787,  # 1s
    (2, 0): -2.7967,   # 2s
    (2, 1): -1.5181,   # 2p
}


def build_na_core_orbital(n_core: int, l_core: int):
    """Return callable r -> R(r) for a hydrogenic Na core orbital.

    Uses Clementi-Raimondi zeta to define the effective Z_eff = n * zeta
    consistent with the convention in phillips_kleinman_cross_center.
    """
    zeta = NA_CR_EXPONENTS[(n_core, l_core)]
    Z_eff = n_core * zeta  # Bethe-Salpeter Z convention for R_nl with screening

    def _R(r):
        return hydrogenic_R_nl_analytical(Z_eff, n_core, l_core, np.asarray(r))

    return _R


def main():
    R_EQ = 3.566  # NaH equilibrium bond length (bohr)
    Z_NA = 11.0
    Z_H = 1.0

    # Quadrature parameters (same as F3 Step 1 converged)
    QUAD = dict(n_rho=80, n_z=100, rho_max=20.0, z_max=20.0)

    print("=" * 78)
    print("Sprint F4 Step 1 — bonding-vs-core overlap diagnostic")
    print(f"NaH R = {R_EQ} bohr")
    print(f"Quadrature: {QUAD}")
    print("=" * 78)

    t0 = time.perf_counter()

    # ------------------------------------------------------------------
    # Reproduce the F3 step-1 h1 reconstruction in the {|Na 3s⟩, |H 1s⟩}
    # subspace to obtain the bonding orbital coefficients.
    # ------------------------------------------------------------------
    # Multi-zeta Na 3s (physical-fit screened-Schrödinger eigenstate)
    na_valence = get_physical_valence_orbitals(11)
    na_3s_mz = None
    for orb in na_valence:
        if orb.n_orbital == 3 and orb.l_orbital == 0:
            na_3s_mz = orb
            break
    assert na_3s_mz is not None, "Na 3s multi-zeta orbital not found"

    def na_3s_call(r):
        return na_3s_mz.evaluate(np.asarray(r))

    # H 1s hydrogenic
    def h_1s_call(r):
        return hydrogenic_R_nl_analytical(1.0, 1, 0, np.asarray(r))

    # On-site eigenvalues used in F3 step 1 multi-zeta branch:
    # Na 3s multi-zeta on-site = -0.170 Ha (Clementi-Roetti HF value)
    # H 1s = -0.5 Ha (hydrogenic exact)
    E_Na3s = -0.170
    E_H1s = -0.5

    # 2x2 h1 matrix in (Na 3s_mz, H 1s) basis from F3 step 1 multi-zeta data:
    # F3 step1: cross_block_h1_NaH_total_Ha = -1.3703 Ha
    # Re-compute here to be self-contained.
    from geovac.cross_block_h1 import (
        inv_r_ss_axial,
        vne_ss_axial,
    )

    # In z-frame, A=Na at z=0, B=H at z=R_EQ
    z_A, z_B = 0.0, R_EQ

    S_AB = overlap_ss_axial(
        na_3s_call, z_A, h_1s_call, z_B, **QUAD
    )
    # V_ne contributions: sum over BOTH nuclei
    Vne_Na = vne_ss_axial(
        na_3s_call, z_A, h_1s_call, z_B,
        z_C=z_A, Z_C=Z_NA, **QUAD
    )
    Vne_H = vne_ss_axial(
        na_3s_call, z_A, h_1s_call, z_B,
        z_C=z_B, Z_C=Z_H, **QUAD
    )
    # Kinetic via T|psi_B> = (E_B + Z_B/r_B)|psi_B>, with E_B the screened
    # eigenvalue (multi-zeta path)
    inv_rB = inv_r_ss_axial(
        na_3s_call, z_A, h_1s_call, z_B, z_C=z_B, **QUAD
    )
    E_B_kin = E_H1s  # H 1s is the B center here
    T_AB = E_B_kin * S_AB + Z_H * inv_rB
    h12 = T_AB + Vne_Na + Vne_H

    # Diagonal h11 / h22 = on-site eigenvalues (multi-zeta convention)
    h11 = E_Na3s
    h22 = E_H1s

    print(f"\nh1 matrix in (Na 3s_mz, H 1s) basis:")
    print(f"  h11 (Na 3s diag) = {h11:+.4f} Ha")
    print(f"  h22 (H 1s diag)  = {h22:+.4f} Ha")
    print(f"  h12 (cross)      = {h12:+.4f} Ha   "
          f"[S={S_AB:+.4f}, Vne_Na={Vne_Na:+.4f}, "
          f"Vne_H={Vne_H:+.4f}, T={T_AB:+.4f}]")

    # Generalized eigenvalue problem (with overlap matrix)
    h_2x2 = np.array([[h11, h12], [h12, h22]])
    s_2x2 = np.array([[1.0, S_AB], [S_AB, 1.0]])
    w, v = eigh(h_2x2, s_2x2)
    # Identify bonding eigenvector: same-sign coefficients
    if v[0, 0] * v[1, 0] > 0:
        bond_idx = 0
    else:
        bond_idx = 1
    bond_coeff = v[:, bond_idx]
    E_bond = w[bond_idx]

    print(f"\nGeneralized eigenvalues:  E_lower = {w[0]:+.4f} Ha, "
          f"E_upper = {w[1]:+.4f} Ha")
    print(f"Bonding eigvec (bonding orbital coeffs in {{Na 3s, H 1s}}):")
    print(f"  c_Na = {bond_coeff[0]:+.4f},  c_H = {bond_coeff[1]:+.4f}")
    print(f"  E_bond = {E_bond:+.4f} Ha")

    # ------------------------------------------------------------------
    # Compute overlap of bonding orbital with each Na [Ne] core orbital.
    # bonding(r) = c_Na * psi_Na3s(r - R_A) + c_H * psi_H1s(r - R_B)
    # Na core orbital lives at A (Na center)
    # ------------------------------------------------------------------
    print(f"\n{'-' * 78}")
    print("Bonding-vs-Na-core overlaps (s-orbital cores; 2p has m=0 component):")
    print(f"{'-' * 78}")

    core_overlaps = {}
    pk_barrier_per_core = {}
    pk_barrier_total = 0.0

    for (n_c, l_c) in [(1, 0), (2, 0), (2, 1)]:
        if l_c != 0:
            # 2p — only the m=0 component participates with the s-bonding orbital
            # by axial symmetry. The s-orbital is the m=0 spherical harmonic;
            # 2p_m=0 = pz which on the bond axis is z/r, axially symmetric.
            # Use the axial 2p_z overlap.
            # We need axial 2p_z evaluation; for simplicity restrict to s-cores
            # in this first-pass diagnostic and report 2p separately via
            # an axial-aware quadrature.
            S_bond_core = compute_bonding_2pz_overlap(
                bond_coeff, na_3s_call, z_A, h_1s_call, z_B,
                NA_CR_EXPONENTS[(2, 1)], **QUAD,
            )
        else:
            R_core = build_na_core_orbital(n_c, l_c)
            # Overlap of each component with the s-core (axial)
            S_Na3s_core = overlap_ss_axial(
                na_3s_call, z_A, R_core, z_A, **QUAD
            )
            S_H1s_core = overlap_ss_axial(
                h_1s_call, z_B, R_core, z_A, **QUAD
            )
            S_bond_core = bond_coeff[0] * S_Na3s_core + bond_coeff[1] * S_H1s_core

        E_c = NA_CORE_ENERGIES[(n_c, l_c)]
        # PK barrier contribution (E_v - E_c) |S|^2 with E_v = E_bond
        # (purely-repulsive PK uses E_v = 0 reference; we report both)
        dH_pk_Ev = (E_bond - E_c) * S_bond_core ** 2
        dH_pk_abs = (0.0 - E_c) * S_bond_core ** 2  # absolute PK, E_v=0
        pk_barrier_per_core[f"({n_c},{l_c})"] = {
            "S_bond_core": S_bond_core,
            "abs_S_pct": 100.0 * abs(S_bond_core),
            "E_c_Ha": E_c,
            "dH_pk_Ev_bond": dH_pk_Ev,
            "dH_pk_abs_Ev0": dH_pk_abs,
        }
        pk_barrier_total += dH_pk_Ev
        print(
            f"  ({n_c},{l_c}) core (E_c = {E_c:+.4f} Ha): "
            f"S_bond_core = {S_bond_core:+.4f}  (|S| = {100*abs(S_bond_core):.2f}%)  "
            f"dH^PK_(Ev=E_bond) = {dH_pk_Ev:+.4f} Ha  "
            f"dH^PK_(Ev=0)    = {dH_pk_abs:+.4f} Ha"
        )
        core_overlaps[f"({n_c},{l_c})"] = S_bond_core

    print(f"\nTotal PK barrier (Ev=E_bond): {pk_barrier_total:+.4f} Ha")
    abs_pk_total = sum(
        v["dH_pk_abs_Ev0"] for v in pk_barrier_per_core.values()
    )
    print(f"Total PK barrier (Ev=0):       {abs_pk_total:+.4f} Ha")

    # ------------------------------------------------------------------
    # Gate logic
    # ------------------------------------------------------------------
    max_abs_overlap_pct = max(
        v["abs_S_pct"] for v in pk_barrier_per_core.values()
    )
    barrier_threshold = 1.0  # Ha
    overlap_threshold_pct = 5.0
    weak_overlap_threshold_pct = 1.0

    print(f"\n{'=' * 78}")
    print("GATE LOGIC")
    print(f"{'=' * 78}")
    print(f"  Max |overlap with any core| = {max_abs_overlap_pct:.2f}%")
    print(f"  Predicted total PK barrier (Ev=0) = {abs_pk_total:+.4f} Ha")
    print()

    if (
        max_abs_overlap_pct > overlap_threshold_pct
        and abs_pk_total > barrier_threshold
    ):
        verdict = "PROCEED_TO_STEP_2"
        verdict_reason = (
            f"max|S|={max_abs_overlap_pct:.2f}% > {overlap_threshold_pct}% "
            f"AND barrier={abs_pk_total:.2f} Ha > {barrier_threshold} Ha. "
            "W1e is real; PK should close it."
        )
    elif max_abs_overlap_pct < weak_overlap_threshold_pct:
        verdict = "STOP_BONDING_ALREADY_ORTHOGONAL"
        verdict_reason = (
            f"max|S|={max_abs_overlap_pct:.2f}% < {weak_overlap_threshold_pct}%. "
            "Bonding orbital already orthogonal to core; W1e isn't a Pauli problem."
        )
    else:
        verdict = "PROCEED_WITH_WEAK_EXPECTED_IMPACT"
        verdict_reason = (
            f"max|S|={max_abs_overlap_pct:.2f}% in intermediate range "
            f"[{weak_overlap_threshold_pct}, {overlap_threshold_pct}]% "
            f"or barrier={abs_pk_total:.2f} Ha in intermediate range; "
            "proceed but flag weak expected impact."
        )

    print(f"VERDICT: {verdict}")
    print(f"REASON:  {verdict_reason}")

    t1 = time.perf_counter()
    print(f"\nWall time: {t1 - t0:.2f} s")

    # Save results
    out = {
        "sprint": "F4 Step 1 — bonding-vs-core overlap diagnostic",
        "date": "2026-05-23",
        "R_AB_bohr": R_EQ,
        "quadrature": QUAD,
        "h1_2x2": {
            "h11_Na3s_mz_diag_Ha": h11,
            "h22_H1s_diag_Ha": h22,
            "h12_cross_Ha": h12,
            "S_AB": S_AB,
        },
        "bonding_orbital": {
            "E_bond_Ha": E_bond,
            "c_Na_3s": float(bond_coeff[0]),
            "c_H_1s": float(bond_coeff[1]),
            "E_upper_Ha": float(w[1 - bond_idx]),
        },
        "core_overlaps": pk_barrier_per_core,
        "summary": {
            "max_abs_overlap_pct": max_abs_overlap_pct,
            "total_pk_barrier_Ev_eq_Ebond_Ha": pk_barrier_total,
            "total_pk_barrier_Ev_eq_0_Ha": abs_pk_total,
            "barrier_threshold_Ha": barrier_threshold,
            "overlap_threshold_pct": overlap_threshold_pct,
            "weak_overlap_threshold_pct": weak_overlap_threshold_pct,
        },
        "verdict": verdict,
        "verdict_reason": verdict_reason,
        "wall_time_s": t1 - t0,
    }

    out_path = os.path.join(_HERE, "data", "sprint_f4_step1_diagnostic.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_path}")

    return out


def compute_bonding_2pz_overlap(
    bond_coeff, na_3s_call, z_A, h_1s_call, z_B, zeta_2p,
    n_rho, n_z, rho_max, z_max,
):
    """Compute <bonding | Na 2p_z> overlap, axial geometry.

    Na 2p_z(r - R_A) = N_2p * (zeta_2p)^... * (z - z_A)/r' * exp(-Z r' / 2)
    Equivalently, the m=0 component of 2p is the pz orbital.

    For an s-orbital bonding state on the bond axis, the overlap with 2p_z
    can be computed by 2D axial Gauss-Legendre quadrature: the angular
    z/r factor on 2p_z survives without phi-integration vanishing.
    """
    from scipy.special import roots_legendre

    # Build axial grid centered at midpoint between A and B
    z_center = 0.5 * (z_A + z_B)
    u, wu = roots_legendre(n_rho)
    rho_pts = rho_max * (u + 1.0) / 2.0
    rho_w = wu * rho_max / 2.0

    v, wv = roots_legendre(n_z)
    z_pts = z_center + z_max * v
    z_w = wv * z_max

    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')

    # Distances
    r_A = np.sqrt(rho_grid ** 2 + (z_grid - z_A) ** 2)
    r_B = np.sqrt(rho_grid ** 2 + (z_grid - z_B) ** 2)

    # Bonding component values: c_Na * psi_Na3s(r_A) + c_H * psi_H1s(r_B)
    R_Na = na_3s_call(r_A.flatten()).reshape(r_A.shape)
    R_H = h_1s_call(r_B.flatten()).reshape(r_B.shape)
    bond_val = bond_coeff[0] * R_Na + bond_coeff[1] * R_H

    # Na 2p_z = R_2p(r_A) * cos(theta_A) = R_2p(r_A) * (z - z_A) / r_A
    # R_2p radial for hydrogenic with Z_eff = 2 * zeta_2p (Bethe-Salpeter conv.)
    Z_eff = 2 * zeta_2p
    R_2p_A = hydrogenic_R_nl_analytical(
        Z_eff, 2, 1, r_A.flatten()
    ).reshape(r_A.shape)
    cos_theta_A = np.divide(
        z_grid - z_A, np.where(r_A > 1e-12, r_A, 1e-12)
    )
    psi_2pz = R_2p_A * cos_theta_A

    # Angular normalization: Y_{1,0} = sqrt(3 / (4 pi)) cos(theta)
    # So 2pz = R_2p * cos(theta) -> need normalization sqrt(3/(4 pi))
    # For an s-component, Y_{0,0} = 1/sqrt(4 pi)
    # Phi-integration: 2*pi
    # Result: 2*pi * (1/sqrt(4 pi)) * sqrt(3/(4 pi)) = (1/2) * sqrt(3/(4 pi)) * 2
    # Cleaner: use full prefactor 2*pi * Y_0 * Y_{1,0} cos(theta) part
    # The s-orbital prefactor is 1/sqrt(4 pi); for 2pz it is sqrt(3/(4 pi)) * cos(theta)
    # Cross-overlap radial integral has dV = rho drho dz dphi factor
    # angular factor = (1/sqrt(4 pi)) * sqrt(3/(4 pi)) * cos(theta)
    # which we already encoded as cos_theta_A above.

    # Prefactor: 2*pi (phi) * 1/sqrt(4 pi) (s) * sqrt(3/(4 pi)) (2pz norm)
    prefactor = 2.0 * np.pi * (1.0 / np.sqrt(4.0 * np.pi)) * np.sqrt(
        3.0 / (4.0 * np.pi)
    )

    integrand = bond_val * psi_2pz * rho_grid
    weight_2d = np.outer(rho_w, z_w)
    raw = np.sum(weight_2d * integrand)
    return float(prefactor * raw)


if __name__ == "__main__":
    main()
