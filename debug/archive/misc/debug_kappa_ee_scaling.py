#!/usr/bin/env python3
"""
debug_kappa_ee_scaling.py
=========================
Verify kappa_ee = 5Z/2 calibration for the chordal V_ee model.

Steps executed (matching the mission specification):
  Step 1 — Algebraic verification (printed inline below)
  Step 4 — S^3 coordinates for five canonical states + degeneracy check
  Step 3 — He (Z=2, n_el=2) and Li (Z=3, n_el=3, max_n=4) energy audits

Author: GeoVac Development Team, February 2026
"""

import os
import sys
import time
import warnings

import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geovac.lattice_index import LatticeIndex


# -----------------------------------------------------------------------
# Step 4 — S^3 coordinate check
# -----------------------------------------------------------------------

def _s3_coord(n: int, l: int, m: int) -> np.ndarray:
    """Fock S^3 embedding for state (n, l, m), p0=1 normalised."""
    if l > 0:
        cos_th = np.clip(float(m) / np.sqrt(float(l * (l + 1))), -1.0, 1.0)
        sin_th = np.sqrt(max(0.0, 1.0 - cos_th * cos_th))
    else:
        cos_th = 1.0
        sin_th = 0.0
    fac = 2.0 * n / (n * n + 1.0)
    return np.array([
        fac * sin_th,
        0.0,
        fac * cos_th,
        (1.0 - n * n) / (1.0 + n * n),
    ])


def print_step4() -> None:
    print("\n" + "=" * 65)
    print("STEP 4 -- S^3 coordinates for canonical states (p0=1 normalised)")
    print("=" * 65)

    states_5 = [(1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 1, 1), (2, 1, -1)]
    coords_5 = [_s3_coord(*s) for s in states_5]

    print(f"\n  {'(n,l,m)':>14}   {'xi_1':>9}  {'xi_2':>9}  {'xi_3':>9}  {'xi_4':>9}  {'|xi|':>7}")
    print("  " + "-" * 65)
    for s, xi in zip(states_5, coords_5):
        print(f"  ({s[0]},{s[1]},{s[2]:+d}):  "
              f"{xi[0]:>9.5f}  {xi[1]:>9.5f}  {xi[2]:>9.5f}  {xi[3]:>9.5f}  "
              f"{np.linalg.norm(xi):>7.6f}")

    # Pairwise chordal distances among the 5 states
    print("\n  Pairwise d^2 = 2(1 - xi_i.xi_j) for the five states:")
    labels = [f"({s[0]},{s[1]},{s[2]:+d})" for s in states_5]
    col_w = 12
    hdr = "  " + " " * 14 + "".join(f"{lb:>{col_w}}" for lb in labels)
    print(hdr)
    print("  " + "-" * (14 + col_w * len(states_5) + 2))
    for i, (si, xi) in enumerate(zip(states_5, coords_5)):
        row = f"  {labels[i]:>14}"
        for j, (sj, xj) in enumerate(zip(states_5, coords_5)):
            if i == j:
                row += f"{'[d^2=4]':>{col_w}}"
            else:
                d_sq = 2.0 * (1.0 - np.dot(xi, xj))
                row += f"{d_sq:>{col_w}.5f}"
        print(row)

    # Degeneracy flag for the 5 states
    degenerate = []
    for i in range(len(states_5)):
        for j in range(i + 1, len(states_5)):
            d_sq = 2.0 * (1.0 - np.dot(coords_5[i], coords_5[j]))
            if d_sq < 1e-6:
                degenerate.append((states_5[i], states_5[j], d_sq))

    if degenerate:
        print("\n  [WARNING] Degenerate pairs among the 5 states (d^2 < 1e-6):")
        for (si, sj, d_sq) in degenerate:
            print(f"    {si} -- {sj}: d^2 = {d_sq:.2e}")
    else:
        print("\n  All 5 states are distinct on S^3 (no degeneracy).")

    # --- Full max_n=4 basis degeneracy audit ---
    print("\n  Full max_n=4 basis: scanning for degenerate S^3 pairs...")
    all_states = [
        (n, l, m)
        for n in range(1, 5)
        for l in range(n)
        for m in range(-l, l + 1)
    ]
    all_coords = np.array([_s3_coord(*s) for s in all_states])
    gram = all_coords @ all_coords.T
    d_sq_full = 2.0 * (1.0 - gram)
    np.fill_diagonal(d_sq_full, 4.0)

    degen_pairs = []
    for i in range(len(all_states)):
        for j in range(i + 1, len(all_states)):
            if d_sq_full[i, j] < 1e-6:
                degen_pairs.append((all_states[i], all_states[j], d_sq_full[i, j]))

    if degen_pairs:
        print(f"  Found {len(degen_pairs)} degenerate pair(s) in max_n=4 basis:")
        for (si, sj, d_sq) in degen_pairs:
            print(f"    {si} == {sj}   d^2 = {d_sq:.2e}")
        print("\n  [NOTE] These map to the same S^3 point (same-n, m=0, different l>0).")
        print("  The d_sq clip at 1e-14 prevents div/zero but V_ee -> 1/1e-14 ~ 1e14 Ha.")
        print("  Affected V_ee matrix entries are physically uninterpretable.")
        print("  For n <= 2 (He, H), no degeneracy exists. For Li (max_n=4), impacts")
        print("  V_ee for excited-state pairs — ground state dominated by n=1,2.")
    else:
        print("  No degenerate pairs found in max_n=4 basis.")


# -----------------------------------------------------------------------
# Step 3 — He and Li energy audits
# -----------------------------------------------------------------------

HE_EXACT_HA: float = -2.9037   # Pekeris 1959
LI_EXACT_HA: float = -7.4781   # Pekeris 1959


def _run_fci_audit(
    label: str,
    n_electrons: int,
    nuclear_charge: int,
    max_n: int,
    exact_energy: float,
    vee_method: str = 'chordal',
) -> dict:
    """Run FCI ground-state audit for one species. Returns result dict."""
    kappa_ee = 5.0 * nuclear_charge / 2.0
    print(f"\n  {label}: Z={nuclear_charge}, n_el={n_electrons}, "
          f"max_n={max_n}, vee={vee_method}, kappa_ee={kappa_ee:.3f} Ha")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t0 = time.perf_counter()
        idx = LatticeIndex(
            n_electrons=n_electrons,
            max_n=max_n,
            nuclear_charge=nuclear_charge,
            vee_method=vee_method,
        )
        H = idx.assemble_hamiltonian()
        t_asm = time.perf_counter() - t0

    k = min(1, idx.n_sd - 2)
    eigvals_raw, _ = eigsh(H, k=k, which="SA")
    E = float(eigvals_raw[np.argsort(eigvals_raw)][0])

    abs_err = abs(E - exact_energy)
    pct_err = abs_err / abs(exact_energy) * 100.0

    print(f"    n_sd = {idx.n_sd:,}  |  assembly = {t_asm:.2f}s")
    print(f"    E_GeoVac = {E:.6f} Ha")
    print(f"    E_exact  = {exact_energy:.6f} Ha")
    print(f"    |dE|     = {abs_err:.4f} Ha  ({pct_err:.2f}%)")

    return {
        "label": label,
        "Z": nuclear_charge,
        "n_el": n_electrons,
        "max_n": max_n,
        "kappa_ee": kappa_ee,
        "n_sd": idx.n_sd,
        "E_geovac": E,
        "E_exact": exact_energy,
        "abs_err": abs_err,
        "pct_err": pct_err,
    }


def print_step3() -> None:
    print("\n" + "=" * 65)
    print("STEP 3 -- Energy audit: He and Li with kappa_ee = 5Z/2")
    print("=" * 65)

    results = []

    # Helium: n_electrons=2, Z=2, max_n=4 -> C(60,2) = 1770 SDs (fast)
    r_he = _run_fci_audit(
        "Helium", n_electrons=2, nuclear_charge=2, max_n=4,
        exact_energy=HE_EXACT_HA,
    )
    results.append(r_he)

    # Lithium: n_electrons=3, Z=3, max_n=4 -> 34,220 SDs (benchmarked)
    r_li = _run_fci_audit(
        "Lithium", n_electrons=3, nuclear_charge=3, max_n=4,
        exact_energy=LI_EXACT_HA,
    )
    results.append(r_li)

    # Summary table
    print("\n" + "=" * 65)
    print("SUMMARY -- kappa_ee = 5Z/2 audit results")
    print("=" * 65)
    print(f"  {'Species':>8}  {'Z':>3}  {'kappa_ee':>9}  "
          f"{'E_GeoVac (Ha)':>14}  {'E_exact (Ha)':>13}  {'err %':>8}")
    print("  " + "-" * 62)
    for r in results:
        print(f"  {r['label']:>8}  {r['Z']:>3}  {r['kappa_ee']:>9.3f}  "
              f"{r['E_geovac']:>14.6f}  {r['E_exact']:>13.4f}  "
              f"{r['pct_err']:>7.2f}%")


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main() -> None:
    print("\n" + "=" * 65)
    print("kappa_ee = 5Z/2 calibration audit")
    print("=" * 65)

    print("\nSTEP 1 -- Algebraic derivation:")
    print("  Exact 1s-1s Hartree integral:  F^0(1s,1s) = 5Z/8 Ha")
    print("  Chordal model, same-orbital:   V_ee = kappa_ee / d^2 = kappa_ee / 4")
    print("  Setting equal:  kappa_ee / 4 = 5Z/8  =>  kappa_ee = 5Z/2")
    print()
    for Z in [2, 3]:
        kappa = 5.0 * Z / 2.0
        v_same = kappa / 4.0
        f0 = 5.0 * Z / 8.0
        ok = "OK" if abs(v_same - f0) < 1e-12 else "FAIL"
        print(f"  Z={Z}: kappa={kappa:.4f}, V_ee(1s,1s)={v_same:.4f} Ha, "
              f"F^0={f0:.4f} Ha  [{ok}]")

    print_step4()
    print_step3()
    print("\n[DONE]")


if __name__ == "__main__":
    main()
