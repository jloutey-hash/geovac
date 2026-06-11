#!/usr/bin/env python3
"""
debug_slater_vee.py
===================
Validate Slater F0 integrals and run He/Li energy audits with vee_method='slater'.

Step 1: Verify F0 normalization and known analytical values.
Step 2: Run He (Z=2, 2e) and Li (Z=3, 3e) FCI audits.

Known exact F0 values (in units of Z):
  F0(1s,1s) = 5Z/8   = 0.6250 Z
  F0(1s,2s) = 17Z/81  = 0.2099 Z
  F0(1s,2p) = 59Z/243 = 0.2428 Z  (independent of m)
  F0(2s,2s) = 77Z/512 = 0.1504 Z
  F0(2s,2p) = 83Z/512 = 0.1621 Z
  F0(2p,2p) = 501Z/2560 = 0.1957 Z  (same m_l pair)

Author: GeoVac Development Team, February 2026
"""

import os
import sys
import time
import warnings

import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.special import genlaguerre
from scipy.integrate import quad
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# -----------------------------------------------------------------------
# Step 1: Standalone F0 validation
# -----------------------------------------------------------------------

def _radial_wf(r: np.ndarray, n: int, l: int, Z: float) -> np.ndarray:
    """Un-normalised hydrogen-like radial wavefunction."""
    rho = 2.0 * Z * r / n
    L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    return rho ** l * np.exp(-rho / 2.0) * L_poly


def compute_f0(n1: int, l1: int, n2: int, l2: int, Z: float) -> float:
    """Compute F0(n1 l1, n2 l2) numerically."""
    r_max = 80.0 / Z

    # Normalisation
    norm1_sq, _ = quad(lambda r: _radial_wf(np.array(r), n1, l1, Z) ** 2 * r ** 2,
                       0, r_max, limit=200)
    norm2_sq, _ = quad(lambda r: _radial_wf(np.array(r), n2, l2, Z) ** 2 * r ** 2,
                       0, r_max, limit=200)
    c1 = 1.0 / np.sqrt(norm1_sq)
    c2 = 1.0 / np.sqrt(norm2_sq)

    def R1(r: float) -> float:
        return float(c1 * _radial_wf(np.array(r), n1, l1, Z))

    def R2(r: float) -> float:
        return float(c2 * _radial_wf(np.array(r), n2, l2, Z))

    def inner(r1: float) -> float:
        if r1 < 1e-30:
            return 0.0
        p1, _ = quad(lambda r2: R2(r2) ** 2 * r2 ** 2, 0, r1, limit=100)
        p2, _ = quad(lambda r2: R2(r2) ** 2 * r2, r1, r_max, limit=100)
        return p1 / r1 + p2

    val, _ = quad(lambda r1: R1(r1) ** 2 * inner(r1) * r1 ** 2, 0, r_max, limit=200)
    return val


def run_f0_validation() -> bool:
    """Validate F0 against known analytical values."""
    print("=" * 65)
    print("STEP 1 — Slater F0 integral validation")
    print("=" * 65)

    Z = 1.0  # Compute in units of Z (scale later)
    known = [
        ((1, 0, 1, 0), 5.0 / 8.0,    "5/8"),
        ((1, 0, 2, 0), 17.0 / 81.0,   "17/81"),
        ((1, 0, 2, 1), 59.0 / 243.0,  "59/243"),
        ((2, 0, 2, 0), 77.0 / 512.0,  "77/512"),
        ((2, 0, 2, 1), 83.0 / 512.0,  "83/512"),
    ]

    all_pass = True
    print(f"\n  {'Pair':>16}  {'F0 (num)':>10}  {'F0 (exact)':>12}  {'rel err':>10}  {'status':>6}")
    print("  " + "-" * 60)

    for (n1, l1, n2, l2), exact, label in known:
        t0 = time.perf_counter()
        numerical = compute_f0(n1, l1, n2, l2, Z)
        dt = time.perf_counter() - t0
        rel_err = abs(numerical - exact) / exact
        ok = rel_err < 0.01  # 1% tolerance for numerical integration
        if not ok:
            all_pass = False
        pair_str = f"({n1}{['s','p','d','f'][l1]},{n2}{['s','p','d','f'][l2]})"
        print(f"  {pair_str:>16}  {numerical:>10.6f}  {exact:>10.6f} ({label:>6})  "
              f"{rel_err:>9.4%}  {'PASS' if ok else 'FAIL'}  [{dt:.2f}s]")

    print(f"\n  Overall: {'ALL PASS' if all_pass else 'SOME FAILED'}")
    return all_pass


# -----------------------------------------------------------------------
# Step 2: Energy audits
# -----------------------------------------------------------------------

HE_EXACT_HA: float = -2.9037
LI_EXACT_HA: float = -7.4781


def _run_fci_audit(
    label: str,
    n_electrons: int,
    nuclear_charge: int,
    max_n: int,
    exact_energy: float,
    vee_method: str = 'slater',
) -> dict:
    """Run FCI ground-state audit for one species."""
    from geovac.lattice_index import LatticeIndex

    print(f"\n  {label}: Z={nuclear_charge}, n_el={n_electrons}, "
          f"max_n={max_n}, vee={vee_method}")

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

    # Print a few V_ee matrix elements for diagnostics
    vee = idx._vee_matrix
    states = idx.lattice.states
    print(f"    V_ee diagnostics:")
    for i in range(min(3, len(states))):
        ni, li, mi = states[i]
        for j in range(i, min(3, len(states))):
            nj, lj, mj = states[j]
            si = f"({ni}{['s','p','d','f'][li]})"
            sj = f"({nj}{['s','p','d','f'][lj]})"
            print(f"      V_ee{si},{sj} = {vee[i, j]:.6f} Ha")

    return {
        "label": label, "Z": nuclear_charge, "n_el": n_electrons,
        "max_n": max_n, "n_sd": idx.n_sd,
        "E_geovac": E, "E_exact": exact_energy,
        "abs_err": abs_err, "pct_err": pct_err,
    }


def run_audits() -> None:
    print(f"\n{'=' * 65}")
    print("STEP 2 — Energy audit with exact Slater F0 integrals")
    print(f"{'=' * 65}")

    results = []
    r_he = _run_fci_audit("Helium", 2, 2, 4, HE_EXACT_HA)
    results.append(r_he)
    r_li = _run_fci_audit("Lithium", 3, 3, 4, LI_EXACT_HA)
    results.append(r_li)

    # Summary table
    print(f"\n{'=' * 65}")
    print("COMPARISON TABLE")
    print(f"{'=' * 65}")
    print(f"  {'Method':>20}  {'Species':>8}  {'E_GeoVac':>12}  {'E_exact':>10}  {'err %':>8}")
    print("  " + "-" * 65)

    # Previous results for comparison
    prev = [
        ("chordal kappa=5Z/2", "He", -2.875616, -2.9037, 0.97),
        ("chordal kappa=5Z/2", "Li", -0.412759, -7.4781, 94.48),
        ("mulliken",           "Li", -3.495767, -7.4781, 53.25),
    ]
    for method, species, e_gv, e_ex, err in prev:
        print(f"  {method:>20}  {species:>8}  {e_gv:>12.6f}  {e_ex:>10.4f}  {err:>7.2f}%")

    for r in results:
        print(f"  {'slater F0':>20}  {r['label']:>8}  {r['E_geovac']:>12.6f}  "
              f"{r['E_exact']:>10.4f}  {r['pct_err']:>7.2f}%")


def main() -> None:
    f0_ok = run_f0_validation()

    if not f0_ok:
        print("\n  [STOP] F0 validation FAILED. Not proceeding to energy audit.")
        return

    run_audits()
    print("\n[DONE]")


if __name__ == "__main__":
    main()
