"""Quick smoke test of xcwg_wilson_loop_scaling at n_max=2 only.

Confirms:
  - Pilot S(L=4) ~ 0.248 reproduced
  - Pilot S(L=6) ~ 0.436 reproduced
  - New L=8 value computed
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

# Import the run_n_max function
from xcwg_wilson_loop_scaling import run_n_max

if __name__ == "__main__":
    result = run_n_max(
        2,
        target_lengths=[4, 6, 8],
        caps={4: 1_000_000, 6: 1_000_000, 8: 1_000_000},
    )
    print()
    print("==== smoke test n_max=2 ====")
    for L, info in result["wilson_action_by_length"].items():
        print(f"  L={L}: count={info['count']}, <S>={info['mean']:.6f}, std={info['std']:.6f}")
    print(f"  scaling: alpha={result['scaling_fit']['alpha']:.4f}, "
          f"R^2={result['scaling_fit']['r_squared']:.4f}, "
          f"verdict={result['scaling_fit']['verdict']}")
    # Pilot values for reference:
    print()
    print("Pilot reference (from xcwg_observables_pilot.json, sample-30):")
    print("  L=4: <S>=0.24769 (sample 30 of 44)")
    print("  L=6: <S>=0.43611 (sample 30 of 144)")
