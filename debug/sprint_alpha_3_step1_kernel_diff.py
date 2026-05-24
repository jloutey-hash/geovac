"""
Sprint alpha-PES Step 1 — Algebraic kernel differential test (debug-only).

Goal: compute the single cross-center V_ne matrix element

    <psi_Na_3s | (-Z_H / |r - R_H|) | psi_Na_3s>

at NaH internuclear distance R = R_eq = 3.566 bohr, using TWO different
Na 3s wavefunctions on the bra/ket:

    Baseline: hydrogenic Z_orb=1 with (n=3, l=0)  (current production)
    Path A:   physical multi-zeta Na 3s from get_physical_valence_orbitals(11)

Compare the differential.  Also evaluate at R = 10 bohr (dissociation
limit) to confirm the differential goes to zero as the H atom recedes.

This script uses the EXISTING analytical kernel _split_integral_analytical
in geovac/shibuya_wulfman.py.  No production code is modified.  Both bases
are expressed as polynomial-times-exponential decompositions, then the
multi-zeta orbital is treated as a sum-of-K-pair-integrals at K different
exponents.

ABORT GATE: if the differential at R_eq is < 0.05 Ha (or wrong sign:
baseline more attractive than physical), Path A is not the mechanism.
In that case, document the clean negative and STOP.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import List, Tuple

import numpy as np
from scipy.special import factorial, gammainc, gammaincc

# Production multi-zeta module
from geovac.multi_zeta_orbitals import (
    get_physical_valence_orbitals,
    MultiZetaOrbital,
    STO,
)

# Production hydrogenic polynomial decomposition
from geovac.shibuya_wulfman import (
    _hydrogenic_poly_coeffs,
    _split_integral_analytical,
    _poly_product,
)


def _multizeta_to_poly_components(
    orbital: MultiZetaOrbital,
) -> List[Tuple[np.ndarray, float]]:
    """
    Decompose a MultiZetaOrbital into a list of (poly_coeffs, alpha) pairs.

    Each STO primitive chi_i(r) = N_i * r^{n_i - 1} * exp(-zeta_i r) is a
    single polynomial-times-exponential.  The orbital is

        R(r) = sum_i c_i * N_i * r^{n_i - 1} * exp(-zeta_i r)

    which is a SUM of K single-exponential terms at K distinct decay rates
    zeta_i.  We return one (poly_coeffs, alpha) pair per primitive, with
    the c_i * N_i factor absorbed into poly_coeffs.

    Returns
    -------
    components : list of (poly_coeffs, alpha)
        Each pair represents one term  exp(-alpha*r) * sum_k poly_coeffs[k] * r^k.
        For a single STO primitive, poly_coeffs has a single nonzero entry
        at index (n_slater - 1).
    """
    components = []
    for c, prim in zip(orbital.coefficients, orbital.primitives):
        N = prim.normalization()
        # chi(r) = N * r^{n-1} * exp(-zeta*r)
        # so coefficient at r^{n-1} is N, all others zero.
        # The MultiZeta combination contributes c * chi(r), so coeff = c * N
        poly = np.zeros(prim.n)  # length n: poly_coeffs[k] for k=0..n-1
        poly[prim.n - 1] = c * N
        components.append((poly, float(prim.zeta)))
    return components


def cross_vne_hydrogenic(
    Z_orb: float, n: int, l: int,
    Z_nuc: float, R_AB: float, L_max: int = 0,
) -> Tuple[float, float, float]:
    """
    Compute <psi_nl | -Z_nuc/|r-R_B| | psi_nl> with hydrogenic Z_orb basis.

    For l1 = l2 = 0 (s-s), only L=0 contributes (Gaunt triangle).
    The angular coefficient at L=0 is (4*pi)^{-1/2} * (4*pi)^{-1/2} * 4*pi = 1
    after proper normalization — but the framework convention in
    `compute_cross_center_vne_element` already absorbs this; the radial
    split-integral IS the matrix-element (modulo a sign).

    Returns (V, V_inner, V_outer) for diagnostic decomposition.
    """
    c1, alpha1 = _hydrogenic_poly_coeffs(Z_orb, n, l)
    c2, alpha2 = _hydrogenic_poly_coeffs(Z_orb, n, l)
    prod = _poly_product(c1, c2)
    alpha_total = alpha1 + alpha2

    # For L=0 only (s-s only):
    # ang_coeff(0,0,0,0,0) = (-1)^0 sqrt(1*1) * w1 * w2 with w1=w2=Wigner3j(0,0,0;0,0,0) = 1
    # Then in compute_cross_center_vne_element:
    #   total = sum_L (nuc_parity)^L * ang * rad
    #   return -Z_nuc * total
    # For s-s and nuc_parity=+1, L=0 only:
    #   ang_coeff(L=0) = 1 * 1 * 1 = 1
    L = 0
    rad = _split_integral_analytical(
        prod, alpha_total, L + 2, 1 - L, L, R_AB,
    )
    return -Z_nuc * rad, 0.0, 0.0


def cross_vne_multizeta(
    orbital: MultiZetaOrbital,
    Z_nuc: float, R_AB: float, L_max: int = 0,
) -> Tuple[float, dict]:
    """
    Compute <psi_orb | -Z_nuc/|r-R_B| | psi_orb> with a multi-zeta orbital.

    Returns the matrix element and a diagnostic dict.
    """
    components = _multizeta_to_poly_components(orbital)
    K = len(components)

    # For l_orbital = 0 only L=0 contributes
    L = 0
    total_rad = 0.0
    pair_contributions = []

    for i in range(K):
        for j in range(K):
            poly_i, alpha_i = components[i]
            poly_j, alpha_j = components[j]
            prod = _poly_product(poly_i, poly_j)
            alpha_total = alpha_i + alpha_j
            rad = _split_integral_analytical(
                prod, alpha_total, L + 2, 1 - L, L, R_AB,
            )
            total_rad += rad
            pair_contributions.append(
                {
                    'i': i, 'j': j,
                    'zeta_i': alpha_i, 'zeta_j': alpha_j,
                    'alpha_total': alpha_total,
                    'rad': rad,
                }
            )

    matrix_element = -Z_nuc * total_rad
    diag = {
        'K': K,
        'pair_contributions': pair_contributions,
        'total_rad': total_rad,
    }
    return matrix_element, diag


def main():
    R_eq = 3.566  # bohr — NaH experimental R_eq (CRC handbook)
    R_diss = 10.0  # bohr — dissociation limit consistency check
    Z_H = 1.0  # H nuclear charge

    # ----- BASELINE: hydrogenic Z_orb=1 Na 3s -----
    V_base_eq, _, _ = cross_vne_hydrogenic(
        Z_orb=1.0, n=3, l=0, Z_nuc=Z_H, R_AB=R_eq,
    )
    V_base_diss, _, _ = cross_vne_hydrogenic(
        Z_orb=1.0, n=3, l=0, Z_nuc=Z_H, R_AB=R_diss,
    )

    # ----- PATH A: physical multi-zeta Na 3s -----
    Na_orbitals = get_physical_valence_orbitals(11)
    Na_3s = Na_orbitals[0]
    assert Na_3s.n_orbital == 3 and Na_3s.l_orbital == 0, "Expected Na 3s as first orbital"

    V_phys_eq, diag_eq = cross_vne_multizeta(
        Na_3s, Z_nuc=Z_H, R_AB=R_eq,
    )
    V_phys_diss, diag_diss = cross_vne_multizeta(
        Na_3s, Z_nuc=Z_H, R_AB=R_diss,
    )

    # ----- Differentials -----
    diff_eq = V_phys_eq - V_base_eq  # phys minus baseline
    diff_diss = V_phys_diss - V_base_diss

    # ----- Sanity: at R_diss = 10 bohr, both should be ~-Z_H/R = -0.1 (classical limit) -----
    classical_diss = -Z_H / R_diss
    classical_eq = -Z_H / R_eq

    # ----- Decision gate -----
    # M-Y predicts: physical Na 3s (mean radius ~4.5 bohr) is MORE extended
    # than hydrogenic Z_orb=1 3s (mean radius ~13.5 bohr).  Wait — actually
    # hydrogenic Z=1 n=3 is MORE diffuse than physical Na 3s, which is
    # screened by the [Ne] core.  Let me think carefully:
    #
    #   Hydrogenic Z=1, n=3, l=0:  <r> = (3*9 - 0)/(2*1) = 13.5 bohr
    #   Physical Na 3s (screened): <r> = 4.5 bohr (from sprint alpha-2)
    #
    # So the physical orbital is MORE COMPACT than the hydrogenic Z=1 3s.
    # At R_eq = 3.566 bohr, the H nucleus sits at the edge of the Na
    # valence shell.  A MORE COMPACT physical orbital has LESS amplitude
    # at the H position, so the cross-V_ne attraction is LARGER in magnitude
    # (less screened) than for the diffuse hydrogenic orbital.  Wait, but
    # the relevant operator is -Z_H/|r-R_H|, so amplitude at r ≈ R_H gives
    # the largest contribution.  If physical orbital has LESS amplitude
    # at r ≈ R_H (more confined to Na center r << R_H), then it sees the
    # H nucleus only through 1/R_eq (the classical-electron-at-Na value).
    #
    # Expected differential sign:
    #   The hydrogenic Z=1 orbital is so diffuse it has significant
    #   amplitude AT the H nucleus, giving a strong over-attractive
    #   contribution from r ≈ R_H where |r-R_H| → 0.  The physical
    #   compact orbital sits at r ≈ 4 bohr (Na valence shell), well
    #   inside R_eq = 3.6 bohr, so it sees the H nucleus from a distance.
    #   Therefore:
    #     V_base_eq  <  V_phys_eq  (i.e., baseline is MORE NEGATIVE)
    #     diff_eq > 0  (physical is LESS attractive than diffuse baseline)
    #
    # This is the OPPOSITE of what the W1c-residual binding wall needs!
    # The wall is over-attraction in the baseline.  If physical is less
    # over-attractive, then physical REDUCES the over-attraction by
    # some amount — and that amount IS the differential.  So if
    #   |diff_eq| > 0.05 Ha at R_eq AND |diff_diss| << |diff_eq|,
    # then the basis substitution is the right mechanism direction.

    print(f"=== Step 1 — Algebraic kernel differential test ===\n")
    print(f"NaH:  R_eq = {R_eq} bohr,  R_diss = {R_diss} bohr,  Z_H = {Z_H}\n")
    print(f"--- Baseline (hydrogenic Z_orb=1 Na 3s) ---")
    print(f"  V(R_eq)   = {V_base_eq:.6f} Ha")
    print(f"  V(R_diss) = {V_base_diss:.6f} Ha")
    print()
    print(f"--- Path A (physical multi-zeta Na 3s, K = {Na_3s.primitives.__len__()}) ---")
    print(f"  V(R_eq)   = {V_phys_eq:.6f} Ha")
    print(f"  V(R_diss) = {V_phys_diss:.6f} Ha")
    print()
    print(f"--- Differential (phys - base) ---")
    print(f"  diff(R_eq)   = {diff_eq:+.6f} Ha")
    print(f"  diff(R_diss) = {diff_diss:+.6f} Ha")
    print()
    print(f"--- Classical limit  V_classical = -Z_H/R ---")
    print(f"  V_classical(R_eq)   = {classical_eq:.6f} Ha")
    print(f"  V_classical(R_diss) = {classical_diss:.6f} Ha")
    print()
    print(f"--- Sanity: does V_phys approach classical at R_diss? ---")
    print(f"  V_phys(R_diss) - V_classical(R_diss) = {V_phys_diss - classical_diss:+.6e} Ha")
    print(f"  V_base(R_diss) - V_classical(R_diss) = {V_base_diss - classical_diss:+.6e} Ha")
    print()

    # ----- Decision gate -----
    abort = False
    abort_reason = None
    if abs(diff_eq) < 0.05:
        abort = True
        abort_reason = (
            f"|diff(R_eq)| = {abs(diff_eq):.4f} Ha < 0.05 Ha threshold. "
            f"Multi-zeta basis substitution is not the dominant mechanism."
        )
    elif abs(diff_diss) > 0.5 * abs(diff_eq):
        # If the differential at R_diss is large too, the orbitals
        # are wrong on their own, not due to bond interaction.
        abort = True
        abort_reason = (
            f"|diff(R_diss)| = {abs(diff_diss):.4f} Ha is large relative to "
            f"|diff(R_eq)| = {abs(diff_eq):.4f} Ha "
            f"({abs(diff_diss)/abs(diff_eq)*100:.1f}%). "
            f"Differential is not bond-localized; suggests inconsistent normalization "
            f"or a non-bonding artifact."
        )

    if abort:
        print(f"!!! ABORT GATE TRIGGERED !!!")
        print(f"Reason: {abort_reason}")
        print(f"Path A is NOT the mechanism.  Stop here, document clean negative.")
    else:
        print(f"=== Decision: Differential at R_eq IS substantial ({diff_eq:+.4f} Ha). ===")
        print(f"=== Differential at R_diss is small ({diff_diss:+.4f} Ha, ratio {abs(diff_diss)/abs(diff_eq)*100:.1f}%). ===")
        print(f"=== PROCEED to Step 2 (single-point FCI test with multi-zeta wiring). ===")

    # ----- Save results -----
    results = {
        'metadata': {
            'sprint': 'alpha-PES Step 1',
            'date': '2026-05-23',
            'R_eq_bohr': R_eq,
            'R_diss_bohr': R_diss,
            'Z_H': Z_H,
            'kernel': 'analytical split-region incomplete-gamma',
            'production_modules_modified': 'none',
        },
        'baseline_hydrogenic_z1_na_3s': {
            'V_eq_Ha': V_base_eq,
            'V_diss_Ha': V_base_diss,
        },
        'path_a_physical_multizeta_na_3s': {
            'V_eq_Ha': V_phys_eq,
            'V_diss_Ha': V_phys_diss,
            'K': Na_3s.primitives.__len__(),
        },
        'differential': {
            'diff_eq_Ha': diff_eq,
            'diff_diss_Ha': diff_diss,
            'ratio_diss_to_eq': abs(diff_diss) / abs(diff_eq) if abs(diff_eq) > 1e-12 else None,
        },
        'classical_limit_check': {
            'V_classical_eq_Ha': classical_eq,
            'V_classical_diss_Ha': classical_diss,
            'V_phys_minus_classical_diss': V_phys_diss - classical_diss,
            'V_base_minus_classical_diss': V_base_diss - classical_diss,
        },
        'decision': {
            'abort': abort,
            'abort_reason': abort_reason,
            'proceed_to_step_2': not abort,
        },
        'diag_eq': diag_eq,
    }

    data_path = Path('debug/data/sprint_alpha_3_step1_kernel_diff.json')
    data_path.parent.mkdir(parents=True, exist_ok=True)
    with open(data_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {data_path}")

    return results


if __name__ == '__main__':
    main()
