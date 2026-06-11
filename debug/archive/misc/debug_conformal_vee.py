#!/usr/bin/env python3
"""
debug_conformal_vee.py
=======================
Sanity checks for conformal-weighted V_ee: V_ee = kappa_ee / d^2_phys
where d^2_phys = Omega_i * Omega_j * |p_i - p_j|^2  (Paper 7 Eq. chordal).

Checks:
  1. V_ee(1s, 2s) for Li must be in range 1-5 Ha.
  2. Same-orbital repulsion (antipodal d^2=4) > cross-shell repulsion.
  3. V_ee decreases monotonically with shell separation.

If check 1 passes, runs He and Li energy audits.

Author: GeoVac Development Team, February 2026
"""

import os
import sys
import time
import warnings

import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# -----------------------------------------------------------------------
# Sanity checks (standalone, not importing LatticeIndex yet)
# -----------------------------------------------------------------------

def _momentum_vec(n: int, l: int, m: int, Z: float) -> np.ndarray:
    """Representative momentum vector for state (n,l,m) at nuclear charge Z."""
    p_mag = Z / float(n)
    if l > 0:
        cos_th = np.clip(float(m) / np.sqrt(float(l * (l + 1))), -1.0, 1.0)
        sin_th = np.sqrt(max(0.0, 1.0 - cos_th ** 2))
    else:
        cos_th = 1.0
        sin_th = 0.0
    return np.array([p_mag * sin_th, 0.0, p_mag * cos_th])


def _omega(n: int, Z: float) -> float:
    """Conformal factor Omega_k = n/Z (on-shell: |p|=p0=Z/n => Omega=1/p0=n/Z)."""
    return float(n) / Z


def _d2_phys(n1: int, l1: int, m1: int,
             n2: int, l2: int, m2: int,
             Z: float) -> float:
    """Conformal-weighted flat-space distance squared."""
    p1 = _momentum_vec(n1, l1, m1, Z)
    p2 = _momentum_vec(n2, l2, m2, Z)
    o1 = _omega(n1, Z)
    o2 = _omega(n2, Z)
    dp = p1 - p2
    return o1 * o2 * float(np.dot(dp, dp))


def _vee(n1: int, l1: int, m1: int,
         n2: int, l2: int, m2: int,
         Z: float, same_orbital: bool = False) -> float:
    """V_ee = kappa_ee / d^2, kappa_ee = 5Z/2."""
    kappa = 5.0 * Z / 2.0
    if same_orbital:
        return kappa / 4.0   # antipodal override
    d2 = _d2_phys(n1, l1, m1, n2, l2, m2, Z)
    if d2 < 1e-14:
        return float('inf')
    return kappa / d2


def run_sanity_checks(Z: float = 3.0) -> bool:
    """Run all three sanity checks for nuclear charge Z. Returns True if check 1 passes."""
    kappa = 5.0 * Z / 2.0
    print(f"\n{'=' * 65}")
    print(f"SANITY CHECKS — Conformal-weighted V_ee for Z={Z:.0f}")
    print(f"{'=' * 65}")
    print(f"  kappa_ee = 5Z/2 = {kappa:.3f} Ha")

    # --- Check 1: V_ee(1s, 2s) for Li ---
    print(f"\n  CHECK 1: V_ee(1s, 2s) must be in [1, 5] Ha")
    p_1s = _momentum_vec(1, 0, 0, Z)
    p_2s = _momentum_vec(2, 0, 0, Z)
    o_1s = _omega(1, Z)
    o_2s = _omega(2, Z)
    dp = p_1s - p_2s
    flat_d2 = float(np.dot(dp, dp))
    d2_phys = o_1s * o_2s * flat_d2
    v_1s_2s = kappa / d2_phys if d2_phys > 1e-14 else float('inf')

    print(f"    p_1s = ({p_1s[0]:.4f}, {p_1s[1]:.4f}, {p_1s[2]:.4f}),  |p_1s| = {np.linalg.norm(p_1s):.4f}")
    print(f"    p_2s = ({p_2s[0]:.4f}, {p_2s[1]:.4f}, {p_2s[2]:.4f}),  |p_2s| = {np.linalg.norm(p_2s):.4f}")
    print(f"    Omega_1s = {o_1s:.6f}")
    print(f"    Omega_2s = {o_2s:.6f}")
    print(f"    |p_1s - p_2s|^2 = {flat_d2:.6f}")
    print(f"    d^2_phys = Omega_1s * Omega_2s * |dp|^2 = {d2_phys:.6f}")
    print(f"    V_ee(1s,2s) = kappa / d^2_phys = {v_1s_2s:.4f} Ha")

    check1_pass = 1.0 <= v_1s_2s <= 5.0
    print(f"    {'PASS' if check1_pass else 'FAIL'}: V_ee(1s,2s) = {v_1s_2s:.4f} Ha "
          f"({'in' if check1_pass else 'NOT in'} [1, 5] Ha)")

    # --- Check 2: Same-orbital > cross-shell ---
    print(f"\n  CHECK 2: Same-orbital repulsion > cross-shell repulsion")
    v_same_1s = _vee(1, 0, 0, 1, 0, 0, Z, same_orbital=True)
    print(f"    V_ee(1s, 1s) [antipodal d^2=4] = {v_same_1s:.4f} Ha")
    print(f"    V_ee(1s, 2s) [conformal]        = {v_1s_2s:.4f} Ha")
    check2_pass = v_same_1s > v_1s_2s
    print(f"    {'PASS' if check2_pass else 'FAIL'}: V_same({v_same_1s:.4f}) "
          f"{'>' if check2_pass else '<='} V_cross({v_1s_2s:.4f})")

    # --- Check 3: Monotonic decrease with shell separation ---
    print(f"\n  CHECK 3: V_ee monotonically decreasing with shell separation")
    pairs = [
        ((1, 0, 0), (1, 0, 0), True,  "(1s,1s) same-orbital"),
        ((1, 0, 0), (2, 0, 0), False, "(1s,2s)"),
        ((1, 0, 0), (3, 0, 0), False, "(1s,3s)"),
        ((2, 0, 0), (2, 0, 0), True,  "(2s,2s) same-orbital"),
        ((2, 0, 0), (3, 0, 0), False, "(2s,3s)"),
    ]
    prev_v = float('inf')
    all_monotonic = True
    for (s1, s2, same, label) in pairs:
        v = _vee(s1[0], s1[1], s1[2], s2[0], s2[1], s2[2], Z, same_orbital=same)
        d2 = 4.0 if same else _d2_phys(s1[0], s1[1], s1[2], s2[0], s2[1], s2[2], Z)
        mono = "ok" if v <= prev_v else "NOT MONOTONIC"
        if v > prev_v:
            all_monotonic = False
        print(f"    {label:25s}  d^2_phys={d2:>8.4f}  V_ee={v:>8.4f} Ha  [{mono}]")
        prev_v = v

    print(f"    {'PASS' if all_monotonic else 'FAIL'}: monotonicity "
          f"{'holds' if all_monotonic else 'VIOLATED'}")

    return check1_pass


# -----------------------------------------------------------------------
# Energy audits
# -----------------------------------------------------------------------

HE_EXACT_HA: float = -2.9037
LI_EXACT_HA: float = -7.4781


def _run_fci_audit(
    label: str,
    n_electrons: int,
    nuclear_charge: int,
    max_n: int,
    exact_energy: float,
    vee_method: str = 'chordal',
) -> dict:
    """Run FCI ground-state audit for one species."""
    from geovac.lattice_index import LatticeIndex

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
        "label": label, "Z": nuclear_charge, "n_el": n_electrons,
        "max_n": max_n, "kappa_ee": kappa_ee, "n_sd": idx.n_sd,
        "E_geovac": E, "E_exact": exact_energy,
        "abs_err": abs_err, "pct_err": pct_err,
    }


def run_audits() -> None:
    print(f"\n{'=' * 65}")
    print("ENERGY AUDIT — Conformal-weighted V_ee with kappa_ee = 5Z/2")
    print(f"{'=' * 65}")

    results = []
    r_he = _run_fci_audit("Helium", 2, 2, 4, HE_EXACT_HA)
    results.append(r_he)
    r_li = _run_fci_audit("Lithium", 3, 3, 4, LI_EXACT_HA)
    results.append(r_li)

    # Summary table (including previous results for comparison)
    print(f"\n{'=' * 65}")
    print("COMPARISON TABLE — All V_ee methods, max_n=4")
    print(f"{'=' * 65}")
    print(f"  {'Method':>20}  {'Species':>8}  {'E_GeoVac':>12}  {'E_exact':>10}  {'err %':>8}")
    print("  " + "-" * 65)

    # Previous results (hardcoded from earlier audits)
    prev = [
        ("chordal kappa=1",    "He",  None,        -2.9037, None),
        ("chordal kappa=1",    "Li",  -6.143094,   -7.4781, 17.85),
        ("chordal kappa=5Z/2", "He",  -2.875616,   -2.9037, 0.97),
        ("chordal kappa=5Z/2", "Li",  -0.412759,   -7.4781, 94.48),
        ("mulliken",           "Li",  -3.495767,   -7.4781, 53.25),
    ]
    for method, species, e_gv, e_ex, err in prev:
        if e_gv is not None:
            print(f"  {method:>20}  {species:>8}  {e_gv:>12.6f}  {e_ex:>10.4f}  {err:>7.2f}%")
        else:
            print(f"  {method:>20}  {species:>8}  {'N/A':>12}  {e_ex:>10.4f}  {'N/A':>8}")

    # New results
    for r in results:
        print(f"  {'conformal 5Z/2':>20}  {r['label']:>8}  {r['E_geovac']:>12.6f}  "
              f"{r['E_exact']:>10.4f}  {r['pct_err']:>7.2f}%")


def main() -> None:
    check1_pass = run_sanity_checks(Z=3.0)

    if not check1_pass:
        print("\n  [STOP] Check 1 FAILED. V_ee(1s,2s) not in [1, 5] Ha.")
        print("  NOT proceeding to energy audit.")
        return

    run_audits()
    print("\n[DONE]")


if __name__ == "__main__":
    main()
