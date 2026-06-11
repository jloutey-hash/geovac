"""Diagnostic: try computing propagation number at n_max=3, N_t=5 separately
with verbose output to identify the failure mode.
"""

import sys
import time
import traceback
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from geovac.operator_system_lorentzian import LorentzianTruncatedOperatorSystem


def main():
    print("=" * 70)
    print("Propagation-number diagnostic at (n_max=3, N_t=5)")
    print("=" * 70)
    import numpy as np
    T = 2.0 * np.pi
    t_build = time.time()
    op_sys = LorentzianTruncatedOperatorSystem(n_max=3, N_t=5, T_max=T / 2.0)
    print(f"Build time: {time.time() - t_build:.2f} s")
    print(f"dim_K = {op_sys.dim_K}")
    print(f"dim(O) = {op_sys.dim}")
    print(f"num multipliers = {len(op_sys.multiplier_matrices)}")
    print(f"achievable envelope dim = {op_sys.achievable_envelope_dim}")

    print("\nAttempting compute_propagation_number(max_k=3, envelope='achievable')...")
    t_prop = time.time()
    try:
        prop_val, dim_seq = op_sys.compute_propagation_number(
            max_k=3, envelope="achievable", verbose=True,
        )
        print(f"  Prop: {prop_val}")
        print(f"  Dim seq: {dim_seq}")
        print(f"  Time: {time.time() - t_prop:.2f} s")
    except Exception as exc:
        print(f"  EXCEPTION: {type(exc).__name__}: {exc}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
