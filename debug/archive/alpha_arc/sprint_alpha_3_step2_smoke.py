"""Sprint alpha-PES Step 2 smoke test: backward-compat of multi_zeta_basis=False.

Verify that:
  1. build_balanced_hamiltonian(nah_spec(), R=3.5, multi_zeta_basis=False)
     gives bit-identical result to the legacy call without the kwarg.
  2. With multi_zeta_basis=True, the run succeeds and h1 is different
     (numerically, with substantial cross_vne magnitude change).
"""

import numpy as np

from geovac.molecular_spec import nah_spec
from geovac.balanced_coupled import build_balanced_hamiltonian


def main():
    spec = nah_spec(max_n=2)
    R = 3.5

    # Run 1: baseline (no kwarg)
    r0 = build_balanced_hamiltonian(spec, R=R, verbose=False)

    # Run 2: explicit multi_zeta_basis=False
    r1 = build_balanced_hamiltonian(
        spec, R=R, verbose=False, multi_zeta_basis=False,
    )

    # Run 3: multi_zeta_basis=True
    r2 = build_balanced_hamiltonian(
        spec, R=R, verbose=True, multi_zeta_basis=True,
    )

    # Backward-compat checks
    h1_diff_01 = np.max(np.abs(r0['h1'] - r1['h1']))
    eri_diff_01 = np.max(np.abs(r0['eri'] - r1['eri']))
    print(f"\n--- Backward-compat (baseline vs explicit False) ---")
    print(f"  max |h1_0 - h1_1|  = {h1_diff_01:.3e}")
    print(f"  max |eri_0 - eri_1| = {eri_diff_01:.3e}")
    assert h1_diff_01 < 1e-14, "Bit-exact backward-compat FAILED for h1"
    assert eri_diff_01 < 1e-14, "Bit-exact backward-compat FAILED for eri"
    print("  >>> bit-exact backward-compat PASS")

    # Multi-zeta vs baseline checks
    h1_diff_02 = np.max(np.abs(r0['h1'] - r2['h1']))
    cross_vne_diff = np.max(np.abs(r0['h1_cross_vne'] - r2['h1_cross_vne']))
    print(f"\n--- Multi-zeta vs baseline diff ---")
    print(f"  max |h1_0 - h1_2|     = {h1_diff_02:.4f} Ha")
    print(f"  max |h1_cross_vne diff| = {cross_vne_diff:.4f} Ha")
    print(f"  multi_zeta_enabled = {r2['multi_zeta_basis_enabled']}")
    print(f"  multi_zeta_diagnostics = {r2['multi_zeta_diagnostics']}")
    assert h1_diff_02 > 0.05, (
        f"Multi-zeta produced too small a diff: {h1_diff_02:.4f} Ha — "
        "may indicate wiring didn't dispatch."
    )
    print("  >>> multi_zeta dispatch active PASS")

    print("\nAll smoke checks PASS.")


if __name__ == '__main__':
    main()
