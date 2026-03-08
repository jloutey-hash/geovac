#!/usr/bin/env python3
"""
debug_slater_full.py
====================
Validate full Slater integrals (F^k + G^k) and run He/Li energy audits
with vee_method='slater_full'.

Step 1: Validate Wigner 3j symbols against known values.
Step 2: Validate c^k angular coefficients.
Step 3: Validate two-electron integrals <ab|cd> for He (1s,1s).
Step 4: Run He and Li FCI audits with slater_full.

Author: GeoVac Development Team, March 2026
"""

import os
import sys
import time
import warnings

import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# -----------------------------------------------------------------------
# Step 1: Wigner 3j validation
# -----------------------------------------------------------------------

def validate_wigner3j() -> bool:
    """Validate Wigner 3j symbols against known analytical values."""
    from geovac.lattice_index import LatticeIndex

    print("=" * 65)
    print("STEP 1 — Wigner 3j validation")
    print("=" * 65)

    # Known values: (j1,j2,j3,m1,m2,m3) -> value
    known = [
        # (0 0 0; 0 0 0) = 1
        ((0, 0, 0, 0, 0, 0), 1.0),
        # (1 1 0; 0 0 0) = (-1)^1 / sqrt(3) = -1/sqrt(3)
        ((1, 1, 0, 0, 0, 0), -1.0 / np.sqrt(3)),
        # (1 1 0; 1 -1 0) = (-1)^(1-1+0) / sqrt(3) = 1/sqrt(3)
        ((1, 1, 0, 1, -1, 0), 1.0 / np.sqrt(3)),
        # (1 1 2; 0 0 0) = sqrt(2/15)
        ((1, 1, 2, 0, 0, 0), np.sqrt(2.0 / 15.0)),
        # Triangle violation: (1 1 3; 0 0 0) = 0
        ((1, 1, 3, 0, 0, 0), 0.0),
        # m-rule violation: (1 1 0; 1 1 0) = 0 (m1+m2+m3 != 0)
        ((1, 1, 0, 1, 1, 0), 0.0),
    ]

    all_pass = True
    for args, exact in known:
        val = LatticeIndex._wigner3j(*args)
        err = abs(val - exact)
        ok = err < 1e-10
        if not ok:
            all_pass = False
        print(f"  3j{args} = {val:+.8f}  exact={exact:+.8f}  "
              f"err={err:.2e}  {'PASS' if ok else 'FAIL'}")

    print(f"\n  Overall: {'ALL PASS' if all_pass else 'SOME FAILED'}")
    return all_pass


# -----------------------------------------------------------------------
# Step 2: c^k coefficient validation
# -----------------------------------------------------------------------

def validate_ck() -> bool:
    """Validate angular coupling coefficients c^k."""
    from geovac.lattice_index import LatticeIndex

    print(f"\n{'=' * 65}")
    print("STEP 2 — c^k angular coefficient validation")
    print("=" * 65)

    # Create a dummy instance just to access the method
    idx = LatticeIndex.__new__(LatticeIndex)

    # c^0(0,0,0,0) = (-1)^0 * sqrt(1*1) * 3j(0,0,0;0,0,0) * 3j(0,0,0;0,0,0)
    #              = 1 * 1 * 1 * 1 = 1
    known = [
        # c^0(l=0,m=0,l'=0,m'=0) = 1
        ((0, 0, 0, 0, 0), 1.0),
        # c^0(l=1,m=0,l'=1,m'=0) should be 1 (same orbital, k=0)
        ((1, 0, 1, 0, 0), 1.0),
        # c^1(l=0,m=0,l'=1,m'=0): parity (0+1+1)=even
        # = (-1)^0 * sqrt(1*3) * 3j(0,1,1;0,0,0) * 3j(0,1,1;0,0,0)
        # 3j(0,1,1;0,0,0) = (-1)^(0-0+0)/sqrt(2*1+1) ... need to check
    ]

    all_pass = True
    for args, exact in known:
        val = idx._ck_coefficient(*args)
        err = abs(val - exact)
        ok = err < 1e-8
        if not ok:
            all_pass = False
        la, ma, lc, mc, k = args
        print(f"  c^{k}(l={la},m={ma},l'={lc},m'={mc}) = {val:+.8f}  "
              f"exact={exact:+.8f}  err={err:.2e}  {'PASS' if ok else 'FAIL'}")

    print(f"\n  Overall: {'ALL PASS' if all_pass else 'SOME FAILED'}")
    return all_pass


# -----------------------------------------------------------------------
# Step 3: Two-electron integral validation (He 1s1s)
# -----------------------------------------------------------------------

def validate_he_eri() -> bool:
    """
    For He with max_n=1 (only 1s orbital), the only ERI is:
    <1s 1s|1s 1s> = F0(1s,1s) = 5Z/8

    With Z=2: F0 = 5*2/8 = 1.25 Ha
    """
    from geovac.lattice_index import LatticeIndex

    print(f"\n{'=' * 65}")
    print("STEP 3 — He two-electron integral validation")
    print("=" * 65)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        idx = LatticeIndex(
            n_electrons=2,
            max_n=2,  # 1s, 2s, 2p states
            nuclear_charge=2,
            vee_method='slater_full',
        )

    # Check ERI values
    print(f"\n  ERI table has {len(idx._eri)} entries")
    states = idx.lattice.states

    # Print some ERI values
    for key, val in sorted(idx._eri.items())[:20]:
        a, b, c, d = key
        sa = f"({states[a][0]}{['s','p','d','f'][states[a][1]]},m={states[a][2]})"
        sb = f"({states[b][0]}{['s','p','d','f'][states[b][1]]},m={states[b][2]})"
        sc = f"({states[c][0]}{['s','p','d','f'][states[c][1]]},m={states[c][2]})"
        sd_s = f"({states[d][0]}{['s','p','d','f'][states[d][1]]},m={states[d][2]})"
        print(f"  <{sa}{sb}|{sc}{sd_s}> = {val:.6f}")

    # Key check: <1s,1s|1s,1s> = F0(1s,1s) = 5Z/8
    Z = 2
    exact_f0_1s1s = 5.0 * Z / 8.0
    eri_1s1s = idx._get_eri(0, 0, 0, 0)  # spatial index 0 = 1s
    err = abs(eri_1s1s - exact_f0_1s1s)
    ok = err < 0.01
    print(f"\n  <1s,1s|1s,1s> = {eri_1s1s:.6f}  exact={exact_f0_1s1s:.6f}  "
          f"err={err:.4f}  {'PASS' if ok else 'FAIL'}")

    return ok


# -----------------------------------------------------------------------
# Step 4: Energy audits
# -----------------------------------------------------------------------

HE_EXACT_HA: float = -2.9037
LI_EXACT_HA: float = -7.4781


def run_audit(
    label: str,
    n_electrons: int,
    nuclear_charge: int,
    max_n: int,
    exact_energy: float,
    vee_method: str = 'slater_full',
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

    if hasattr(idx, '_eri'):
        print(f"    ERI entries = {len(idx._eri)}")

    return {
        "label": label, "Z": nuclear_charge, "n_el": n_electrons,
        "max_n": max_n, "n_sd": idx.n_sd,
        "E_geovac": E, "E_exact": exact_energy,
        "abs_err": abs_err, "pct_err": pct_err,
    }


def run_audits() -> None:
    print(f"\n{'=' * 65}")
    print("STEP 4 — Energy audits with full Slater integrals")
    print(f"{'=' * 65}")

    results = []

    # He with slater_full (small basis first)
    r_he_small = run_audit("He full (n=2)", 2, 2, 2, HE_EXACT_HA, 'slater_full')
    results.append(r_he_small)

    # He with slater_full
    r_he = run_audit("He full (n=4)", 2, 2, 4, HE_EXACT_HA, 'slater_full')
    results.append(r_he)

    # He with slater (baseline comparison)
    r_he_base = run_audit("He (slater F0)", 2, 2, 4, HE_EXACT_HA, 'slater')
    results.append(r_he_base)

    # Li with slater_full (small basis first)
    r_li_small = run_audit("Li full (n=2)", 3, 3, 2, LI_EXACT_HA, 'slater_full')
    results.append(r_li_small)

    # Li with slater_full
    r_li = run_audit("Li full (n=3)", 3, 3, 3, LI_EXACT_HA, 'slater_full')
    results.append(r_li)

    # Li with slater (baseline)
    r_li_base = run_audit("Li (slater F0)", 3, 3, 4, LI_EXACT_HA, 'slater')
    results.append(r_li_base)

    # Summary
    print(f"\n{'=' * 65}")
    print("COMPARISON TABLE")
    print(f"{'=' * 65}")
    print(f"  {'Method':>25}  {'Species':>8}  {'E_GeoVac':>12}  {'E_exact':>10}  {'err %':>8}")
    print("  " + "-" * 70)

    for r in results:
        print(f"  {r['label']:>25}  {'He' if r['Z']==2 else 'Li':>8}  "
              f"{r['E_geovac']:>12.6f}  {r['E_exact']:>10.4f}  {r['pct_err']:>7.2f}%")


def main() -> None:
    ok1 = validate_wigner3j()
    if not ok1:
        print("\n  [STOP] Wigner 3j validation FAILED.")
        return

    ok2 = validate_ck()

    ok3 = validate_he_eri()

    run_audits()
    print("\n[DONE]")


if __name__ == "__main__":
    main()
