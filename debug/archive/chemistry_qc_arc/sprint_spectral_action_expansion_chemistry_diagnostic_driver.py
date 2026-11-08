"""Sprint Spectral-Action-Expansion-Chemistry diagnostic (2026-06-07).

Question
========
At fixed R = R_eq, expand the Marcolli-vS spectral action

    S(D, Lambda) = Tr exp(-D^2 / Lambda^2)

of the GeoVac chemistry Dirac D = h1 (the M-vS-2 confirmed assembled
gauge-network Hamiltonian, bit-exact = Track CD balanced_coupled h1) as a
function of Lambda.  Does the expansion factorize into chemistry-meaningful
coefficients, analogous to the Chamseddine-Connes SM expansion?

Setup
=====
- LiH default lih_spec() at R = 3.015 bohr.  M = 15 spatial orbitals.
- H2 two-center pilot (from bratteli_h2_pilot_driver) at R = 1.4 bohr.
  M = 10.
- NaH default nah_spec() at R = 3.566 bohr (NIST CCCBDB).  M = 10.

For each system:
  1. Build h1 via build_balanced_hamiltonian(cross_block_h1=True).
  2. Compute Lambda-sweep S(D, Lambda) on log-spaced Lambda in [0.1, 100].
  3. Compute trace invariants Tr(D^{2k}) for k = 0..6.
  4. Substitution panel:
       D_full = h1
       D_diag = diag(h1)
       D_off  = h1 - diag(h1)
     compute S(Lambda) for each in the same Lambda grid.
  5. Verify the trivial large-Lambda Taylor expansion:
       S(Lambda) = M - Tr(D^2)/Lambda^2 + Tr(D^4)/(2*Lambda^4) - ...
  6. Compare leading non-trivial coefficient Tr(D^2) to chemistry
     trace invariants.

Verdict
=======
POSITIVE  : Fit coefficients have chemistry-specific structure beyond the
            trivial Taylor expansion.
NEGATIVE  : Fit coefficients ARE the trivial Tr(D^{2k}) Taylor terms; no
            chemistry-specific structural content.
PARTIAL   : Leading coefficient meaningful but sub-leading not.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import scipy.linalg as la

try:
    sys.stdout.reconfigure(encoding='utf-8')
except (AttributeError, Exception):
    pass

from geovac.molecular_spec import lih_spec, nah_spec, MolecularSpec, OrbitalBlock
from geovac.balanced_coupled import build_balanced_hamiltonian


# ----------------------------------------------------------------------
# H2 two-center spec (lifted from bratteli_h2_pilot_driver.py)
# ----------------------------------------------------------------------


def make_h2_two_center_spec(R: float, max_n: int = 2) -> MolecularSpec:
    """Two-center H_2 spec with one lone_pair block per H atom."""
    nuclei = [
        {'Z': 1.0, 'position': (0.0, 0.0, -R / 2.0), 'label': 'H_a'},
        {'Z': 1.0, 'position': (0.0, 0.0,  R / 2.0), 'label': 'H_b'},
    ]
    blocks = [
        OrbitalBlock(
            label='H_a_atomic',
            block_type='lone_pair',
            Z_center=1.0,
            n_electrons=1,
            max_n=max_n,
            center_nucleus_idx=0,
        ),
        OrbitalBlock(
            label='H_b_atomic',
            block_type='lone_pair',
            Z_center=1.0,
            n_electrons=1,
            max_n=max_n,
            center_nucleus_idx=1,
        ),
    ]
    return MolecularSpec(
        name='H2',
        blocks=blocks,
        nuclear_repulsion_constant=1.0 / R,
        description=f'H2 two-center pilot at R={R} bohr',
        nuclei=nuclei,
        R=R,
    )


# ----------------------------------------------------------------------
# Core observables
# ----------------------------------------------------------------------


def spectral_action(D: np.ndarray, Lambda: float) -> float:
    """S(D, Lambda) = Tr exp(-D^2/Lambda^2)."""
    A = -(D @ D) / (Lambda * Lambda)
    return float(np.real(np.trace(la.expm(A))))


def trace_powers(D: np.ndarray, k_max: int = 6) -> Dict[str, float]:
    """Trace(D^{2k}) for k = 0..k_max."""
    out = {'Tr_I': int(D.shape[0])}
    D2 = D @ D
    Dpow = D2.copy()
    out['Tr_D2'] = float(np.real(np.trace(Dpow)))
    for k in range(2, k_max + 1):
        Dpow = Dpow @ D2
        out[f'Tr_D{2*k}'] = float(np.real(np.trace(Dpow)))
    return out


def chemistry_invariants(h1: np.ndarray) -> Dict[str, float]:
    """Chemistry-side scalar invariants of h1."""
    h_diag = np.diag(np.diag(h1))
    h_off = h1 - h_diag
    cross = float(np.real(np.trace(h_diag @ h_off)))
    return {
        'sum_diag': float(np.sum(np.diag(h1))),  # = Tr(h1)
        'sum_sq_diag': float(np.sum(np.diag(h1) ** 2)),
        'sum_sq_off': float(np.sum(h_off ** 2)),
        'frob_norm_sq_total': float(np.sum(h1 ** 2)),
        'Tr_diag_off': cross,
    }


def taylor_predict_S(M: int, traces: Dict[str, float], Lambda: float,
                     k_max: int = 6) -> float:
    """Predicted S(Lambda) via the trivial Taylor expansion of exp(-D^2/L^2).

    S = sum_{k=0}^{k_max} (-1)^k Tr(D^{2k}) / (k! Lambda^{2k})
    """
    val = float(M)  # k=0 term
    for k in range(1, k_max + 1):
        val += ((-1) ** k) * traces[f'Tr_D{2*k}'] / (
            float(math.factorial(k)) * (Lambda ** (2 * k))
        )
    return val


# ----------------------------------------------------------------------
# Per-system computation
# ----------------------------------------------------------------------


def compute_system_diagnostic(
    name: str, spec, R: float, max_n: int = 2,
    Lambda_grid: np.ndarray = None,
    k_max_trace: int = 6,
) -> Dict[str, Any]:
    """Build h1 for `spec` at R, compute Lambda-sweep + invariants."""
    print(f"\n  Building {name} at R = {R} bohr, max_n = {max_n}...")
    t0 = time.time()
    result = build_balanced_hamiltonian(
        spec, R=R, cross_block_h1=True, verbose=False,
    )
    h1 = np.asarray(result['h1'])
    M = h1.shape[0]
    print(f"    M = {M}, build time {time.time() - t0:.1f}s")

    # Hermitize defensively (h1 should already be Hermitian)
    h1 = 0.5 * (h1 + h1.T)

    h_diag = np.diag(np.diag(h1))
    h_off = h1 - h_diag

    # ---- Trace invariants on the three substitutions
    print(f"    computing trace powers up to D^{2*k_max_trace}...")
    tr_full = trace_powers(h1, k_max=k_max_trace)
    tr_diag = trace_powers(h_diag, k_max=k_max_trace)
    tr_off = trace_powers(h_off, k_max=k_max_trace)

    chem_inv = chemistry_invariants(h1)

    # Verify the identity Tr(h_full^2) = Tr(h_diag^2) + Tr(h_off^2)
    tr_d2_check = tr_diag['Tr_D2'] + tr_off['Tr_D2']
    tr_d2_resid = abs(tr_full['Tr_D2'] - tr_d2_check)

    # ---- Spectral action on Lambda grid for all three substitutions
    print(f"    sweeping Lambda grid ({len(Lambda_grid)} points)...")
    t_sw = time.time()
    S_full = np.array([spectral_action(h1, L) for L in Lambda_grid])
    S_diag = np.array([spectral_action(h_diag, L) for L in Lambda_grid])
    S_off = np.array([spectral_action(h_off, L) for L in Lambda_grid])
    print(f"    Lambda sweep {time.time() - t_sw:.1f}s")

    # ---- Compare large-Lambda S(Lambda) to Taylor prediction
    L_big = Lambda_grid[Lambda_grid >= 5.0]
    S_pred_full = np.array([
        taylor_predict_S(M, tr_full, L, k_max=k_max_trace) for L in L_big
    ])
    S_actual_full_big = S_full[Lambda_grid >= 5.0]
    taylor_resid_max = float(np.max(np.abs(S_pred_full - S_actual_full_big)))
    taylor_resid_rel_max = float(np.max(
        np.abs(S_pred_full - S_actual_full_big) /
        np.maximum(np.abs(S_actual_full_big), 1e-12)
    ))

    # Also: Taylor at moderate Lambda (e.g. Lambda = max|eigenvalue|)
    eigs = la.eigvalsh(h1)
    lam_max = float(np.max(np.abs(eigs)))
    lam_min_nz = float(np.min(np.abs(eigs[np.abs(eigs) > 1e-12])))

    # Where does the Taylor series start failing? Find smallest Lambda
    # where |S_taylor - S_actual| < 1e-6
    taylor_converged_Lambda = None
    for i, L in enumerate(Lambda_grid):
        try:
            S_pred = taylor_predict_S(M, tr_full, L, k_max=k_max_trace)
            if abs(S_pred - S_full[i]) < 1e-6:
                taylor_converged_Lambda = float(L)
                break
        except (OverflowError, FloatingPointError):
            continue

    return {
        'name': name,
        'R_bohr': R,
        'max_n': max_n,
        'M': M,
        'h1_trace_full': float(np.trace(h1)),
        'h1_eigs_min_nz': lam_min_nz,
        'h1_eigs_max': lam_max,
        'h1_eigenvalues': sorted(eigs.tolist()),
        'tr_full': tr_full,
        'tr_diag': tr_diag,
        'tr_off': tr_off,
        'chem_invariants': chem_inv,
        'tr_D2_decomposition_check_resid': tr_d2_resid,
        'Lambda_grid': Lambda_grid.tolist(),
        'S_Lambda_full': S_full.tolist(),
        'S_Lambda_diag': S_diag.tolist(),
        'S_Lambda_off': S_off.tolist(),
        'S_Lambda_taylor_predict_full_for_large_L': S_pred_full.tolist(),
        'taylor_resid_max_for_Lambda_ge_5': taylor_resid_max,
        'taylor_resid_rel_max_for_Lambda_ge_5': taylor_resid_rel_max,
        'taylor_converged_smallest_Lambda': taylor_converged_Lambda,
    }


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------


def print_system_table(name: str, d: Dict[str, Any]) -> None:
    print(f"\n--- {name} (M = {d['M']}, R = {d['R_bohr']} bohr) ---")
    print(f"  h1 trace             = {d['h1_trace_full']:.6f}")
    print(f"  |lambda|_max         = {d['h1_eigs_max']:.4f}")
    print(f"  |lambda|_min (nz)    = {d['h1_eigs_min_nz']:.4e}")
    tr = d['tr_full']
    tr_d = d['tr_diag']
    tr_o = d['tr_off']
    chem = d['chem_invariants']
    print(f"  Tr(h1)               = {chem['sum_diag']:12.6f}")
    print(f"  Tr(h1^2)             = {tr['Tr_D2']:12.6f}")
    print(f"  Tr(h_diag^2)         = {tr_d['Tr_D2']:12.6f}  "
          f"(= sum h_pp^2)")
    print(f"  Tr(h_off^2)          = {tr_o['Tr_D2']:12.6f}  "
          f"(= sum_{{p!=q}} h_pq^2)")
    print(f"  Tr(h_diag * h_off)   = {chem['Tr_diag_off']:12.4e}  "
          f"(should be ~0 by orthogonality of diag/off)")
    print(f"  identity check       = {d['tr_D2_decomposition_check_resid']:.4e}")
    print(f"  Tr(h1^4)             = {tr['Tr_D4']:12.6f}")
    print(f"  Tr(h1^6)             = {tr['Tr_D6']:12.6f}")
    print(f"  Tr(h1^8)             = {tr['Tr_D8']:12.6f}")
    print(f"  Tr(h1^10)            = {tr['Tr_D10']:12.6f}")
    print(f"  Tr(h1^12)            = {tr['Tr_D12']:12.6f}")
    print(f"  Taylor expansion check at Lambda >= 5:")
    print(f"    max |S_actual - S_taylor|     = "
          f"{d['taylor_resid_max_for_Lambda_ge_5']:.4e}")
    print(f"    max rel |S_actual - S_taylor| = "
          f"{d['taylor_resid_rel_max_for_Lambda_ge_5']:.4e}")
    cvg = d['taylor_converged_smallest_Lambda']
    if cvg is not None:
        print(f"    smallest Lambda where Taylor agrees to 1e-6: {cvg:.3f}")
    print(f"    (|lambda|_max = {d['h1_eigs_max']:.3f}, so Taylor must")
    print(f"    converge for Lambda >> sqrt(spectral_radius(D^2)))")


def compare_systems_invariants(systems: List[Dict[str, Any]]) -> None:
    print("\n" + "=" * 78)
    print("Cross-system comparison of trace invariants")
    print("=" * 78)
    hdr = f"  {'invariant':>22} | " + " | ".join(
        f"{d['name']:>12}" for d in systems
    )
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for key, label in [
        ('M', 'M (dim)'),
        ('Tr_h1', 'Tr(h1)'),
        ('Tr_D2_full', 'Tr(h1^2)'),
        ('Tr_D2_diag', 'Tr(h_diag^2)'),
        ('Tr_D2_off', 'Tr(h_off^2)'),
        ('off_frac', 'Tr(h_off^2)/Tr(h^2)'),
        ('Tr_D4_full', 'Tr(h1^4)'),
        ('Tr_D6_full', 'Tr(h1^6)'),
    ]:
        row = f"  {label:>22} | "
        for d in systems:
            if key == 'M':
                v = d['M']
            elif key == 'Tr_h1':
                v = d['h1_trace_full']
            elif key == 'Tr_D2_full':
                v = d['tr_full']['Tr_D2']
            elif key == 'Tr_D2_diag':
                v = d['tr_diag']['Tr_D2']
            elif key == 'Tr_D2_off':
                v = d['tr_off']['Tr_D2']
            elif key == 'off_frac':
                v = d['tr_off']['Tr_D2'] / d['tr_full']['Tr_D2']
            elif key == 'Tr_D4_full':
                v = d['tr_full']['Tr_D4']
            elif key == 'Tr_D6_full':
                v = d['tr_full']['Tr_D6']
            row += f"{v:>12.4g} | "
        print(row.rstrip(' |'))


def print_S_vs_Taylor_table(d: Dict[str, Any], n_pts: int = 8) -> None:
    """Show S(Lambda) actual vs Taylor prediction at a few representative Lambda."""
    L = np.array(d['Lambda_grid'])
    S_full = np.array(d['S_Lambda_full'])
    S_pred = np.array(d['S_Lambda_taylor_predict_full_for_large_L'])
    L_big_mask = L >= 5.0
    L_big = L[L_big_mask]
    S_big = S_full[L_big_mask]
    print(f"\n  {d['name']} S(Lambda) actual vs k_max=6 Taylor (only Lambda >= 5):")
    print(f"    {'Lambda':>10} | {'S_actual':>14} | {'S_taylor':>14} | "
          f"{'|diff|':>12}")
    n_show = min(n_pts, len(L_big))
    idxs = np.linspace(0, len(L_big) - 1, n_show).astype(int)
    for i in idxs:
        diff = abs(S_big[i] - S_pred[i])
        print(f"    {L_big[i]:>10.3f} | {S_big[i]:>14.8f} | "
              f"{S_pred[i]:>14.8f} | {diff:>12.4e}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------


def main():
    print("=" * 78)
    print("Spectral-Action-Expansion-Chemistry diagnostic")
    print("Sprint date: 2026-06-07")
    print("=" * 78)
    print("\nQuestion: at fixed R = R_eq, does S(D, Lambda) expansion factorize")
    print("into chemistry-meaningful coefficients (cf. Chamseddine-Connes SM)?")
    print("Decision gate:")
    print("  POSITIVE if fit coefficients beyond trivial Tr(D^{2k}) Taylor")
    print("           series show chemistry-specific structure.")
    print("  NEGATIVE if S(Lambda) is exactly the matrix-Taylor expansion.")
    print("  PARTIAL  if leading coefficient meaningful, sub-leading not.")

    # Lambda grid: log-spaced over 3 orders of magnitude
    Lambda_grid = np.logspace(np.log10(0.1), np.log10(100.0), 50)

    print(f"\nLambda grid: {len(Lambda_grid)} points, log-spaced in "
          f"[{Lambda_grid[0]:.3f}, {Lambda_grid[-1]:.2f}]")

    systems = []

    # ---- 1. LiH default at R_eq
    print("\n[1/3] LiH default lih_spec() at R = 3.015 bohr")
    spec_lih = lih_spec(R=3.015, max_n=2)
    d_lih = compute_system_diagnostic(
        'LiH', spec_lih, R=3.015, max_n=2, Lambda_grid=Lambda_grid,
    )
    print_system_table('LiH', d_lih)
    systems.append(d_lih)

    # ---- 2. H2 two-center
    print("\n[2/3] H2 two-center pilot at R = 1.4 bohr")
    spec_h2 = make_h2_two_center_spec(R=1.4, max_n=2)
    d_h2 = compute_system_diagnostic(
        'H2', spec_h2, R=1.4, max_n=2, Lambda_grid=Lambda_grid,
    )
    print_system_table('H2', d_h2)
    systems.append(d_h2)

    # ---- 3. NaH default at R_eq
    print("\n[3/3] NaH default nah_spec() at R = 3.566 bohr")
    spec_nah = nah_spec(R=3.566, max_n=2)
    d_nah = compute_system_diagnostic(
        'NaH', spec_nah, R=3.566, max_n=2, Lambda_grid=Lambda_grid,
    )
    print_system_table('NaH', d_nah)
    systems.append(d_nah)

    # ---- Cross-system comparison
    compare_systems_invariants(systems)

    # ---- Show S vs Taylor table for each
    for d in systems:
        print_S_vs_Taylor_table(d)

    # ---- Substitution panel
    print("\n" + "=" * 78)
    print("Substitution panel: S(Lambda) under D = h1, diag(h1), off(h1)")
    print("=" * 78)
    L_test = [0.5, 1.0, 2.0, 4.0, 10.0]
    for d in systems:
        L_arr = np.array(d['Lambda_grid'])
        S_full = np.array(d['S_Lambda_full'])
        S_d = np.array(d['S_Lambda_diag'])
        S_o = np.array(d['S_Lambda_off'])
        print(f"\n  {d['name']} (M = {d['M']}):")
        print(f"    {'Lambda':>8} | {'S(full)':>12} | {'S(diag)':>12} | "
              f"{'S(off)':>12} | {'S(d)+S(o)-M':>12}")
        for Ltest in L_test:
            i = int(np.argmin(np.abs(L_arr - Ltest)))
            additive = S_d[i] + S_o[i] - d['M']
            print(f"    {L_arr[i]:>8.3f} | {S_full[i]:>12.6f} | "
                  f"{S_d[i]:>12.6f} | {S_o[i]:>12.6f} | "
                  f"{additive:>12.6f}")

    # ---- Verdict synthesis
    print("\n" + "=" * 78)
    print("VERDICT SYNTHESIS")
    print("=" * 78)

    taylor_resid_global = max(
        d['taylor_resid_rel_max_for_Lambda_ge_5'] for d in systems
    )
    print(f"\nTest 1: large-Lambda Taylor expansion captures S(Lambda)?")
    print(f"  Max rel residual at Lambda >= 5 across {len(systems)} systems: "
          f"{taylor_resid_global:.4e}")
    if taylor_resid_global < 1e-8:
        finding_1 = "trivial_taylor_at_large_lambda"
    else:
        finding_1 = "non_taylor_residual"

    print(f"\nTest 2: Tr(D^2) splits into diag + off?")
    for d in systems:
        full = d['tr_full']['Tr_D2']
        diag = d['tr_diag']['Tr_D2']
        off = d['tr_off']['Tr_D2']
        ratio_diag = diag / full if full != 0 else float('nan')
        ratio_off = off / full if full != 0 else float('nan')
        print(f"  {d['name']:>4}: Tr(h_diag^2)/Tr(h^2) = {ratio_diag:6.3f}, "
              f"Tr(h_off^2)/Tr(h^2) = {ratio_off:6.3f}, "
              f"residual {d['tr_D2_decomposition_check_resid']:.4e}")

    print(f"\nTest 3: is the diag/off split chemistry-meaningful?")
    for d in systems:
        off_frac = d['tr_off']['Tr_D2'] / d['tr_full']['Tr_D2']
        print(f"    {d['name']:>4}: bond-coupling fraction Tr(h_off^2)/Tr(h^2) "
              f"= {off_frac:6.3f}")

    out_path = (
        Path(__file__).parent / 'data' /
        'spectral_action_expansion_chemistry.json'
    )
    out_data = {
        'sprint': 'spectral_action_expansion_chemistry_diagnostic',
        'date': '2026-06-07',
        'Lambda_grid': Lambda_grid.tolist(),
        'systems': systems,
        'verdict': {
            'finding_1_taylor_at_large_lambda': finding_1,
            'taylor_resid_global_max_rel': taylor_resid_global,
            'final_verdict': 'NEGATIVE_with_partial_caveat',
            'rationale': (
                'S(D, Lambda) for finite-dim D is exactly the matrix-exp '
                'Taylor series sum (-1)^k Tr(D^{2k})/(k! Lambda^{2k}) when '
                'Lambda >> ||D||_op. Coefficients Tr(D_diag^2) and '
                'Tr(D_off^2) ARE chemistry-meaningful as (kinetic + same-'
                'center V_ne) and cross-center V_ne respectively, but this '
                'reflects the block structure of h1 itself, not an emergent '
                'property of the spectral action.'
            ),
        },
    }

    def _jsonable(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_jsonable(v) for v in obj]
        return obj

    out_data = _jsonable(out_data)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out_data, f, indent=2)
    print(f"\nData written to {out_path}")
    print("=" * 78)


if __name__ == '__main__':
    main()
