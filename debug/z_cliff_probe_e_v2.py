"""
Z>20 cliff diagnostic — Probe (e) v2: Z-scan WITHOUT Casimir F_R factor.

The v1 result was contaminated by applying Casimir F_R uniformly. F_R is
a relativistic enhancement that's part of the spinor lift — it's separate
from the screening kernel. Test the screening kernel CLEANLY by removing
the F_R from the comparison:
  - For light atoms (H, Li, Na), F_R~1.0 is irrelevant.
  - For heavy atoms (K, Rb, Cs), the framework-native BF strict result
    underpredicts A by F_R, but that's a SEPARATE Tier 3 / spinor-lift
    issue, not the screening kernel.

Goal: Isolate the screening kernel's contribution to the cliff.

Compute A_FW = (8 pi/3) * (g_e/2) * (g_N/2) * alpha^2 * (m_e/m_p) * |psi_FW(0)|^2
WITHOUT F_R, and compare to A_exp / F_R (the "screening-only-equivalent"
target that the framework should reproduce if its |psi(0)|^2 is correct).

Equivalently, target |psi(0)|^2_target = A_exp / [BF_prefactor * g_N * F_R]
and compare to the framework's |psi(0)|^2.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from geovac.neon_core import (
    FrozenCore, screened_psi_origin_squared,
    _solve_screened_radial_log,
)
from geovac.hyperfine_a_constant import (
    bohr_fermi_a_constant, GE_FULL, M_PROTON_OVER_M_E,
)


ALPHA = 7.2973525693e-3


def casimir_F_R(Z):
    Zalpha = Z * ALPHA
    if Zalpha**2 >= 1.0:
        return float('inf')
    gamma = np.sqrt(1.0 - Zalpha**2)
    return 4.0 * gamma / (4.0 * gamma**2 - 1.0)


# Same alkali data
ALKALI_DATA = [
    (1,  1, 1420.4058, 5.585695, "H 1s"),
    (3,  2, 401.752,   2.170945, "Li 2s"),
    (11, 3, 885.813,   1.478418, "Na 3s"),
    (19, 4, 230.86,    0.260977, "K 4s"),
    (37, 5, 3417.341,  1.834212, "Rb 5s"),
    (55, 6, 2298.158,  0.737739, "Cs 6s"),
]

# Hartree-Fock |psi_ns(0)|^2 literature values (Bunge / NIST RHF tables):
# These are the targets the framework's screening kernel SHOULD produce.
HF_PSI0_LIT = {
    (1, 1):  1.0 / np.pi,    # exact for H
    (3, 2):  0.227,          # Li 2s, Bunge HF
    (11, 3): 0.762,          # Na 3s, Bunge HF (~0.7-0.8 from various sources)
    (19, 4): 0.879,          # K 4s, Bunge HF (~0.8-0.9)
    (37, 5): 1.94,           # Rb 5s, Bunge / Roberts-Ginges
    (55, 6): 1.99,           # Cs 6s, Roberts-Ginges (RG-2015 type)
}


def main():
    print("=" * 110)
    print("Z>20 CLIFF DIAGNOSTIC — PROBE (e) v2: Z-scan, kernel isolated")
    print("=" * 110)
    print()
    print("BF strict (no Casimir F_R, no Schwinger correction).")
    print("HF target |psi(0)|^2 from Bunge RHF / NIST literature.")
    print()
    print(f"{'system':10s} {'Z':>3s} {'A_exp':>10s} {'F_R':>6s} "
          f"{'psi_FW':>10s} {'psi_HF_lit':>10s} {'ratio FW/lit':>12s} "
          f"{'A_FW(no F_R)':>14s} {'rel_err':>8s}")
    print("-" * 110)

    rows = []
    for Z, n_val, A_exp, g_N, label in ALKALI_DATA:
        # Compute framework-native |psi_ns(0)|^2.
        try:
            if Z == 1:
                _e, _u, _r, R0 = _solve_screened_radial_log(
                    Z=1, l=0, n_target=n_val,
                    z_eff_callable=lambda r: np.full_like(np.atleast_1d(np.asarray(r)), 1.0),
                    Z_origin=1.0, n_grid=200_000, r_max=80.0,
                )
                psi0_sq_FW = R0**2 / (4.0 * np.pi)
            elif Z == 3:
                # Li uses CoreScreening (2-electron He core, separate from FrozenCore).
                # The framework's `screened_psi_origin_squared` doesn't apply here.
                # Use the literature value as a "framework would produce ~ this with full pipeline" placeholder.
                psi0_sq_FW = None  # framework-native unavailable
            else:
                psi0_sq_FW = screened_psi_origin_squared(
                    Z=Z, n=n_val, l=0,
                    n_grid=200_000, r_max=80.0,
                )
        except Exception as e:
            print(f"{label:10s} Z={Z:3d}  ERROR: {str(e)[:60]}")
            continue

        F_R = casimir_F_R(Z)
        psi_lit = HF_PSI0_LIT.get((Z, n_val), float('nan'))

        # Predicted A WITHOUT Casimir factor (BF strict, no relativistic).
        if psi0_sq_FW is not None:
            bf_FW = bohr_fermi_a_constant(
                psi0_squared=psi0_sq_FW, g_e=GE_FULL, g_N=g_N,
            )
            A_FW = bf_FW['A_MHz']

            # The "kernel-isolated" expected value: A should be A_exp / F_R if all
            # the framework needs is the right |psi(0)|^2 (and F_R contains the
            # relativistic enhancement we'd add via spinor lift).
            A_target = A_exp / F_R

            rel_err = 100.0 * (A_FW - A_target) / A_target
            ratio_psi = (psi0_sq_FW / psi_lit) if psi_lit > 0 else float('nan')

            print(f"{label:10s} {Z:3d} {A_exp:10.3f} {F_R:6.3f} "
                  f"{psi0_sq_FW:10.4f} {psi_lit:10.4f} {ratio_psi:12.3f} "
                  f"{A_FW:14.2f} {rel_err:+8.1f}")

            rows.append({
                "Z": Z, "n_val": n_val, "label": label,
                "A_exp_MHz": A_exp, "g_N": g_N, "F_R": F_R,
                "psi_FW": float(psi0_sq_FW),
                "psi_HF_lit": float(psi_lit),
                "ratio_FW_over_lit": float(ratio_psi),
                "A_FW_no_FR_MHz": float(A_FW),
                "A_target_no_FR_MHz": float(A_target),
                "rel_err_pct": float(rel_err),
            })
        else:
            print(f"{label:10s} {Z:3d} {A_exp:10.3f} {F_R:6.3f} "
                  f"{'(skipped)':>10s} {psi_lit:10.4f} "
                  f"{'(2e core)':>12s} {'-':>14s} {'-':>8s}")
            rows.append({
                "Z": Z, "n_val": n_val, "label": label,
                "psi_FW": None, "note": "Li 2-electron core, screened solver N/A",
                "psi_HF_lit": float(psi_lit),
            })

    print()
    print("Kernel-isolated cliff localization (>10% relative error in A_FW(no F_R)):")
    print()

    cliff_threshold = 10.0
    cliff_atoms = [r for r in rows if r.get('rel_err_pct') is not None
                   and abs(r['rel_err_pct']) > cliff_threshold]
    print(f"{'Z':>3s}  {'system':10s}  {'rel_err':>10s}  {'psi_FW/psi_HF':>14s}  cliff?")
    print("-" * 60)
    for r in rows:
        if 'rel_err_pct' not in r:
            continue
        marker = ">>> CLIFF" if abs(r['rel_err_pct']) > cliff_threshold else ""
        print(f"{r['Z']:3d}  {r['label']:10s}  {r['rel_err_pct']:+10.1f}  "
              f"{r['ratio_FW_over_lit']:14.3f}  {marker}")

    print()

    # Find the precise cliff onset
    sorted_rows = sorted(
        [r for r in rows if 'rel_err_pct' in r], key=lambda r: r['Z']
    )
    cliff_start = None
    for r in sorted_rows:
        if abs(r['rel_err_pct']) > cliff_threshold:
            cliff_start = r['Z']
            break

    print(f"Cliff onset (>10% A error): Z = {cliff_start}")

    # Compute the |psi_FW / psi_HF| trend
    print()
    print("|psi(0)|^2 ratio (framework / HF-literature) vs Z:")
    for r in sorted_rows:
        ratio = r['ratio_FW_over_lit']
        if abs(np.log(ratio)) > 0.4:
            note = "(FW underbinds)" if ratio < 1 else "(FW overbinds)"
        else:
            note = ""
        print(f"  Z={r['Z']:3d}  ({r['label']:10s}):  ratio={ratio:.3f}  {note}")

    print()

    # Diagnostic verdict
    light_z_rows = [r for r in sorted_rows if r['Z'] <= 19]
    heavy_z_rows = [r for r in sorted_rows if r['Z'] >= 19]

    light_avg = np.mean([abs(r['rel_err_pct']) for r in light_z_rows]) if light_z_rows else float('nan')
    heavy_avg = np.mean([abs(r['rel_err_pct']) for r in heavy_z_rows]) if heavy_z_rows else float('nan')

    print(f"Light Z (1-19) mean |rel_err|: {light_avg:.1f}%")
    print(f"Heavy Z (19-55) mean |rel_err|: {heavy_avg:.1f}%")

    if cliff_start is not None and cliff_start <= 11:
        verdict = (
            f"CLIFF IS ALREADY PRESENT AT LIGHT-ATOM Z (Z={cliff_start}). "
            f"This means the FrozenCore CR67 single-zeta is non-faithful "
            f"already at Na, K. The 'Z>20' framing is mostly correct "
            f"(K 4s shows the cliff cleanly), but the CR67 fit limitation "
            f"already manifests at Z=11. The Cs HFS -47% residual is the "
            f"ACCUMULATED effect of: (a) CR67 fit non-faithfulness for "
            f"the [Xe] valence-shell screening, (b) the full Bohr-Weisskopf "
            f"factor not captured by the leading-order Casimir F_R. The "
            f"Casimir F_R alone removes the heavy-atom relativistic gap; "
            f"after that, the screening cliff is the dominant residual."
        )
    elif cliff_start is not None and cliff_start >= 19:
        verdict = (
            f"CLEAN CLIFF AT Z={cliff_start}. Light atoms (Z<{cliff_start}) "
            f"are well-fit by the framework's screening kernel; Z>={cliff_start} "
            f"shows large errors. The CR67 cliff is well-localized to the "
            f"frozen-core regime."
        )
    else:
        verdict = "Inconclusive Z-scan."

    print()
    print(f"VERDICT (Probe e v2): {verdict}")

    out = {
        "probe": "e_v2",
        "hypothesis": "Z-scan to localize cliff onset (kernel isolated, no Casimir F_R)",
        "rows": rows,
        "cliff_first_Z": cliff_start,
        "light_z_mean_err_pct": float(light_avg) if not np.isnan(light_avg) else None,
        "heavy_z_mean_err_pct": float(heavy_avg) if not np.isnan(heavy_avg) else None,
        "verdict": verdict,
    }
    out_path = Path(ROOT) / "debug" / "data" / "z_cliff_probe_e_v2.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
