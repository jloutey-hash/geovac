r"""
debug/compute_q5p_tc2c_higher_cutoff.py

Sprint Q5'-Tannakian-Closure TC-2c:\ finite-cutoff reconstruction at
higher $n_{\max}$.

Extends TC-2a (abelian factor at $n_{\max} = 2$, dim 15) and TC-2b
($SL_2$ factor on PW panel, dim 3) to $n_{\max} \in \{3, 4\}$. The
abelian-factor argument is identical at every cutoff (same Jordan-block
witness pattern, same naturality $+$ unit structure), so the only
question is whether the bit-exact $\Q$-linear computation scales.

Expected results:

| $n_{\max}$ | $N(n_{\max})$ | $3 N(n_{\max})$ | predicted dim Aut$^\otimes_{\mathrm{abelian}}$ |
|:----------:|:-------------:|:---------------:|:---------------------------------------------:|
| 1          | 2             | 6               | 6 (verified in TC-2a test)                    |
| 2          | 5             | 15              | 15 (verified in TC-2a)                        |
| 3          | 9             | 27              | 27 (TC-2c primary)                            |
| 4          | 14            | 42              | 42 (TC-2c stretch)                            |

Combined with the $SL_2$ factor on the PW panel (dim 3, $n_{\max}$-independent
since it acts on the orthogonal $j_{\max}$ axis), the combined dim at
each cutoff is $3 N(n_{\max}) + 3$.

Strategy
--------
Re-uses the TC-2a witness-variety pipeline:
``build_witness_panel``, ``parameterize_eta``,
``collect_naturality_constraints``, ``solve_linear_system``,
``extract_solution_structure``, ``verify_phi_recovery`` — all
parameterised by ``n_max``.

References
----------
- Sprint TC-2a memo (`debug/sprint_q5p_tc2a_aut_equality_memo.md`).
- Sprint TC-2b memo (`debug/sprint_q5p_tc2b_sl2_aut_equality_memo.md`).
- Deligne--Milne 1982 Theorems 2.3 and 2.11.
"""

from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path


def _load_tc2a_driver():
    """Import the TC-2a driver helpers."""
    debug_dir = Path(__file__).parent
    spec = importlib.util.spec_from_file_location(
        "tc2a_driver",
        debug_dir / "compute_q5p_tc2a_aut_equality.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


tc2a = _load_tc2a_driver()


def run_at_cutoff(n_max):
    """Bit-exact TC-2a pipeline at the given cutoff."""
    t_start = time.time()
    gens, panel = tc2a.build_witness_panel(n_max)
    n_gens = len(gens)
    predicted_dim = n_gens  # 3 N(n_max)

    etas, eta_T_mat, all_symbols = tc2a.parameterize_eta(gens)
    constraints, n_morphisms_per_pair = tc2a.collect_naturality_constraints(
        gens, panel, etas, eta_T_mat, n_max
    )
    n_total_morphisms = sum(n_morphisms_per_pair.values())
    tc2a.add_unit_normalization(constraints, eta_T_mat)

    result = tc2a.solve_linear_system(constraints, all_symbols)
    dim_solved = result["dim_solution_variety"]

    # Phi recovery check
    if dim_solved == predicted_dim:
        sol_struct = tc2a.extract_solution_structure(
            constraints, all_symbols, gens, etas, eta_T_mat
        )
        phi_check = tc2a.verify_phi_recovery(gens, sol_struct)
        n_free_vars = sol_struct.get("n_free_vars")
        eta_T_solved = sol_struct.get("eta_T_solved")
        phi_match = phi_check["all_match"]
    else:
        n_free_vars = None
        eta_T_solved = None
        phi_match = None

    elapsed = time.time() - t_start
    return {
        "n_max": n_max,
        "n_gens": n_gens,
        "predicted_dim_abelian": predicted_dim,
        "predicted_dim_combined": predicted_dim + 3,
        "computed_dim": dim_solved,
        "n_total_morphisms": n_total_morphisms,
        "n_constraints": result["n_constraints_kept"],
        "rank_A": result["rank_A"],
        "rank_aug": result["rank_aug"],
        "consistent": result["consistent"],
        "n_free_vars": n_free_vars,
        "eta_T_solved": eta_T_solved,
        "phi_recovery_all_match": phi_match,
        "abelian_match": dim_solved == predicted_dim,
        "wall_time_seconds": float(elapsed),
    }


def main():
    print("=" * 72)
    print("Sprint Q5'-Tannakian-Closure TC-2c")
    print("Higher-cutoff reconstruction at n_max in {3, 4}")
    print("=" * 72)

    results = {}

    for n_max in (3, 4):
        print(f"\nRunning at n_max = {n_max}...")
        r = run_at_cutoff(n_max)
        results[n_max] = r
        print(f"  n_gens (= 3 N({n_max})):  {r['n_gens']}")
        print(f"  morphisms in panel:       {r['n_total_morphisms']}")
        print(f"  constraints:              {r['n_constraints']}")
        print(f"  A shape:                  ({r['n_constraints']}, {3 * r['n_gens'] + 1})")
        print(f"  rank_A:                   {r['rank_A']}")
        print(f"  consistent:               {r['consistent']}")
        print(f"  computed dim:             {r['computed_dim']}")
        print(f"  predicted dim (abelian):  {r['predicted_dim_abelian']}")
        print(f"  abelian_match:            {r['abelian_match']}")
        if r["phi_recovery_all_match"] is not None:
            print(f"  Phi recovery:             {r['phi_recovery_all_match']}")
            print(f"  eta_T solved:             {r['eta_T_solved']}")
        print(f"  Combined dim (+SL_2):     {r['predicted_dim_combined']}")
        print(f"  Wall time:                {r['wall_time_seconds']:.2f} s")

    overall = all(r["abelian_match"] for r in results.values())
    verdict = "POSITIVE" if overall else "FAIL"

    print("\n" + "=" * 72)
    print(f"Overall TC-2c verdict: {verdict}")
    print("=" * 72)

    # Summary table
    print("\nReconstruction dim across cutoffs (combining TC-2a, 2b, 2c):")
    print(f"  {'n_max':>6}  {'3 N':>5}  {'abelian':>8}  {'+ SL_2':>8}  {'wall(s)':>8}")
    print(f"  {'-'*6}  {'-'*5}  {'-'*8}  {'-'*8}  {'-'*8}")
    print(f"  {1:>6}  {6:>5}  {6:>8}  {6+3:>8}  {'(test)':>8}")
    print(f"  {2:>6}  {15:>5}  {15:>8}  {15+3:>8}  {'(TC-2a)':>8}")
    for n_max, r in results.items():
        print(f"  {n_max:>6}  {r['n_gens']:>5}  {r['computed_dim']:>8}  {r['predicted_dim_combined']:>8}  {r['wall_time_seconds']:>8.2f}")

    out_data = {
        "sprint": "Q5-prime-Tannakian-Closure-TC-2c",
        "results_by_cutoff": {str(k): v for k, v in results.items()},
        "overall_verdict": verdict,
        "summary": {
            "n_max_tested": list(results.keys()),
            "abelian_match_all": overall,
        },
    }
    data_path = Path(__file__).parent / "data" / "sprint_q5p_tc2c_higher_cutoff.json"
    data_path.parent.mkdir(parents=True, exist_ok=True)
    with open(data_path, "w") as f:
        json.dump(out_data, f, indent=2)
    print(f"\nData saved to {data_path}")

    return out_data


if __name__ == "__main__":
    main()
