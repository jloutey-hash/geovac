r"""
Sprint Q5'-Tannakian-Closure TC-1c --- rigidity on
$\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$.

Goal
----
Verify the rigidity axioms bit-exact:\ for every finite-dim rep $V$ in
$\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$, the dual rep
$V^\vee = \mathrm{Hom}_\mathbb{Q}(V, \mathbb{Q})$ with the contragredient
Hopf action $X_g^{V^\vee} = -(X_g^V)^T$ admits evaluation
$\mathrm{ev}_V : V^\vee \otimes V \to \mathbf{1}$ and coevaluation
$\mathrm{coev}_V : \mathbf{1} \to V \otimes V^\vee$ satisfying the two
snake (zigzag) identities.

Closes the THIRD of four sprint-scale Tannakian-closure prerequisites
named by PS-4 (v3.69.0). TC-1a (abelian, v3.70.0) and TC-1b (symmetric
monoidal, v3.71.0) supply the substrate;\ TC-1d (fiber functor) is the
remaining sprint-scale stone before the multi-year TC-1e first stone.

Five rigidity panels at $n_{\max} \in \{2, 3\}$
-----------------------------------------------
For each of 4 reps in the panel (T1 unit, T2 trivial-2-dim, J2 Jordan-2,
J3 Jordan-3) at each cutoff:

1. **Dual rep contragredient action.** $X_g^{V^\vee} = -(X_g^V)^T$ matches
   `dual_rep(V).X(g)` bit-exact at every non-zero generator, with
   nilpotency and pairwise commutativity preserved. 1 check per rep.

2. **Evaluation intertwines.** $\mathrm{ev}_V$ is a valid rep morphism:\
   $\mathrm{ev}_V \cdot X_g^{V^\vee \otimes V} = 0$ at every non-zero generator
   (right side is zero because $\mathbf{1}$ has zero action). 1 check per rep.

3. **Coevaluation intertwines.** $\mathrm{coev}_V$ is a valid rep morphism:\
   $X_g^{V \otimes V^\vee} \cdot \mathrm{coev}_V = 0$ at every non-zero generator.
   1 check per rep.

4. **First snake identity.** $(\mathrm{ev}_V \otimes \mathrm{id}_{V^\vee}) \circ
   (\mathrm{id}_{V^\vee} \otimes \mathrm{coev}_V) = \mathrm{id}_{V^\vee}$
   bit-exact as a $\dim V \times \dim V$ matrix identity. 1 check per rep.

5. **Second snake identity.** $(\mathrm{id}_V \otimes \mathrm{ev}_V) \circ
   (\mathrm{coev}_V \otimes \mathrm{id}_V) = \mathrm{id}_V$
   bit-exact as a $\dim V \times \dim V$ matrix identity. 1 check per rep.

5 checks $\times$ 4 reps $\times$ 2 cutoffs $= 40$ bit-exact identities.

Plus per-cutoff structural identities:

- **Double dual** $(V^\vee)^\vee = V$ as strict equality of endo data
  (because $S^2 = \mathrm{id}$ on the abelian primitive Hopf algebra),
  for each of the 4 reps. 4 checks $\times$ 2 cutoffs $= 8$.
- **Unit self-dual** $\mathbf{1}^\vee = \mathbf{1}$. 1 check $\times$ 2 cutoffs $= 2$.

Total: $40 + 8 + 2 = 50$ bit-exact zero residuals expected.

Output
------
- ``debug/data/sprint_q5p_tc1c_rigidity.json``
- Console summary at end.

Discipline
----------
Bit-exact ``sympy.Rational`` / ``sympy.Integer`` throughout. No floats. No PSLQ.

References
----------
- TC-1a memo ``debug/sprint_q5p_tc1a_abelian_memo.md`` (abelian substrate).
- TC-1b memo ``debug/sprint_q5p_tc1b_monoidal_memo.md`` (symmetric monoidal).
- v3.61.0 Track A memo ``debug/sprint_q5p_stage2_hopf_memo.md``
  (abelian primitive Hopf $\mathcal{H}_{\mathrm{GV}}$, antipode $S(x_g) = -x_g$).
- Deligne, P.; Milne, J. S. ``Tannakian categories'' (1982), \S 1
  (rigidity in tensor categories).
- Mac Lane, S. *Categories for the Working Mathematician* (1998), Ch. VII
  (monoidal categories), Ch. VII.7 (duality and rigidity).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

from sympy import Integer, Matrix

from geovac.tannakian import (
    FinDimRep,
    coevaluation_morphism,
    dual_rep,
    evaluation_morphism,
    trivial_rep,
    unit_object,
    verify_coevaluation_intertwines,
    verify_double_dual_iso,
    verify_dual_action,
    verify_evaluation_intertwines,
    verify_snake_identity_first,
    verify_snake_identity_second,
    verify_unit_self_dual,
)


# =====================================================================
# Panel substrate
# =====================================================================


def build_panel(n_max: int) -> Dict[str, FinDimRep]:
    T1 = unit_object(n_max)
    T1.label = "T1_unit"
    T2 = trivial_rep(n_max, dim=2)
    T2.label = "T2"
    J2 = FinDimRep(
        n_max=n_max, dim=2,
        endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])},
        label="J2",
    )
    J3 = FinDimRep(
        n_max=n_max, dim=3,
        endos={(1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])},
        label="J3",
    )
    return {"T1": T1, "T2": T2, "J2": J2, "J3": J3}


# =====================================================================
# Axiom panels
# =====================================================================


REP_NAMES = ["T1", "T2", "J2", "J3"]


def run_dual_action_panel(panel: Dict[str, FinDimRep]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for name in REP_NAMES:
        v = verify_dual_action(panel[name])
        v["rep"] = name
        results.append(v)
    return results


def run_evaluation_panel(panel: Dict[str, FinDimRep]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for name in REP_NAMES:
        v = verify_evaluation_intertwines(panel[name])
        v["rep"] = name
        results.append(v)
    return results


def run_coevaluation_panel(panel: Dict[str, FinDimRep]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for name in REP_NAMES:
        v = verify_coevaluation_intertwines(panel[name])
        v["rep"] = name
        results.append(v)
    return results


def run_snake_first_panel(panel: Dict[str, FinDimRep]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for name in REP_NAMES:
        v = verify_snake_identity_first(panel[name])
        v["rep"] = name
        results.append(v)
    return results


def run_snake_second_panel(panel: Dict[str, FinDimRep]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for name in REP_NAMES:
        v = verify_snake_identity_second(panel[name])
        v["rep"] = name
        results.append(v)
    return results


def run_double_dual_panel(panel: Dict[str, FinDimRep]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for name in REP_NAMES:
        v = verify_double_dual_iso(panel[name])
        v["rep"] = name
        results.append(v)
    return results


def run_unit_self_dual(n_max: int) -> Dict[str, Any]:
    return verify_unit_self_dual(n_max)


def count_zero_residuals(
    dual_results: List[Dict[str, Any]],
    ev_results: List[Dict[str, Any]],
    coev_results: List[Dict[str, Any]],
    snake1_results: List[Dict[str, Any]],
    snake2_results: List[Dict[str, Any]],
    double_dual_results: List[Dict[str, Any]],
    unit_self_dual: Dict[str, Any],
) -> Tuple[int, int]:
    total = 0
    zero = 0
    # Five axiom panels (dual action / ev / coev / snake1 / snake2) at 1 per rep
    for panel in (
        dual_results, ev_results, coev_results,
        snake1_results, snake2_results,
    ):
        for r in panel:
            total += 1
            if r["bit_exact"]:
                zero += 1
    # Double dual: 1 per rep
    for r in double_dual_results:
        total += 1
        if r["bit_exact"]:
            zero += 1
    # Unit self-dual: 1 check
    total += 1
    if unit_self_dual["bit_exact"]:
        zero += 1
    return zero, total


# =====================================================================
# Main driver
# =====================================================================


def main() -> None:
    print("=" * 72)
    print("Sprint Q5'-Tannakian-Closure  TC-1c")
    print("Rigidity on Rep_fin(H_GV)")
    print("=" * 72)

    t_global = time.time()
    cutoffs = [2, 3]
    grand_zero = 0
    grand_total = 0
    persist: Dict[str, Any] = {
        "sprint": "Q5'-Tannakian-Closure TC-1c",
        "date_iso": "2026-06-06",
        "cutoffs": cutoffs,
        "results_by_cutoff": {},
    }

    for n_max in cutoffs:
        print(f"\n--- Cutoff n_max = {n_max} ---")
        t_c = time.time()
        panel = build_panel(n_max)
        dual_results = run_dual_action_panel(panel)
        ev_results = run_evaluation_panel(panel)
        coev_results = run_coevaluation_panel(panel)
        snake1_results = run_snake_first_panel(panel)
        snake2_results = run_snake_second_panel(panel)
        double_dual_results = run_double_dual_panel(panel)
        unit_self_dual = run_unit_self_dual(n_max)
        print(f"  [1] Dual action:           {sum(r['bit_exact'] for r in dual_results)}/{len(dual_results)} reps pass")
        print(f"  [2] Evaluation morphism:    {sum(r['bit_exact'] for r in ev_results)}/{len(ev_results)} reps pass")
        print(f"  [3] Coevaluation morphism:  {sum(r['bit_exact'] for r in coev_results)}/{len(coev_results)} reps pass")
        print(f"  [4] First snake identity:   {sum(r['bit_exact'] for r in snake1_results)}/{len(snake1_results)} reps pass")
        print(f"  [5] Second snake identity:  {sum(r['bit_exact'] for r in snake2_results)}/{len(snake2_results)} reps pass")
        print(f"  [+] Double dual (V^vv = V): {sum(r['bit_exact'] for r in double_dual_results)}/{len(double_dual_results)} reps pass")
        print(f"  [+] Unit self-dual:         {'PASS' if unit_self_dual['bit_exact'] else 'FAIL'}")
        zero, total = count_zero_residuals(
            dual_results, ev_results, coev_results,
            snake1_results, snake2_results,
            double_dual_results, unit_self_dual,
        )
        print(f"  --> {zero} / {total} bit-exact zero residuals ({(time.time() - t_c):.2f}s)")
        grand_zero += zero
        grand_total += total
        persist["results_by_cutoff"][str(n_max)] = {
            "dual_action": dual_results,
            "evaluation": ev_results,
            "coevaluation": coev_results,
            "snake_first": snake1_results,
            "snake_second": snake2_results,
            "double_dual": double_dual_results,
            "unit_self_dual": unit_self_dual,
            "zero_residuals": zero,
            "total_identities": total,
        }

    elapsed = time.time() - t_global
    print("\n" + "=" * 72)
    print(f"Total bit-exact zero residuals : {grand_zero} / {grand_total}")
    print(f"Total wall time                : {elapsed:.2f}s")
    print("=" * 72)

    persist["totals"] = {
        "zero_residuals": grand_zero,
        "total_identities": grand_total,
    }
    persist["wall_time_seconds"] = elapsed

    out_path = Path("debug/data/sprint_q5p_tc1c_rigidity.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(persist, f, indent=2, default=str)
    print(f"\nData written to: {out_path}")


if __name__ == "__main__":
    main()
