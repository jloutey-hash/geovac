r"""
Sprint Q5'-Tannakian-Closure TC-1d --- fiber functor $\omega:
\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}) \to
\mathrm{Vec}_\mathbb{Q}$.

Verifies the four Deligne--Milne 1982 fiber-functor properties bit-exact:
unit preservation, $\otimes$-preservation, exactness (kernel + cokernel
+ direct sum preservation), faithfulness.

Closes the THIRD of four sprint-scale Tannakian-closure prerequisites
named by PS-4 (v3.69.0).

Per cutoff $n_{\max} \in \{2, 3\}$:
1. $\omega(\mathbf{1}) = \mathbb{Q}$:                1
2. $\otimes$-preservation:                           4 pairs
3. Kernel preservation:                              4 morphisms
4. Cokernel preservation:                            4 morphisms
5. Direct sum preservation:                          3 pairs
6. Faithfulness:                                     4 tests

Per cutoff: 1 + 4 + 4 + 4 + 3 + 4 = 20. Two cutoffs: 40.

Built from `geovac/tannakian.py` (TC-1a/b/c/d infrastructure).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict

from sympy import Matrix, eye

from geovac.tannakian import (
    FinDimRep,
    RepMorphism,
    fiber_functor_object,
    trivial_rep,
    unit_object,
    verify_omega_faithful,
    verify_omega_preserves_cokernel,
    verify_omega_preserves_direct_sum,
    verify_omega_preserves_kernel,
    verify_omega_tensor_preservation,
    verify_omega_unit,
)


def build_panel(n_max: int) -> Dict[str, Any]:
    T1 = unit_object(n_max)
    T1.label = "T1"
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
    f1 = RepMorphism(T1, J2, Matrix([[1], [0]]), label="f1")
    f2 = RepMorphism(J2, T1, Matrix([[0, 1]]), label="f2")
    f3 = RepMorphism(J2, J2, eye(2), label="f3")
    f4 = RepMorphism(T1, T2, Matrix([[1], [0]]), label="f4")
    return {
        "T1": T1, "T2": T2, "J2": J2, "J3": J3,
        "f1": f1, "f2": f2, "f3": f3, "f4": f4,
    }


def main() -> None:
    print("=" * 72)
    print("Sprint Q5'-Tannakian-Closure  TC-1d")
    print("Fiber functor omega: Rep_fin(H_GV) -> Vec_Q")
    print("=" * 72)

    t_global = time.time()
    cutoffs = [2, 3]
    grand_zero = 0
    grand_total = 0
    persist: Dict[str, Any] = {
        "sprint": "Q5'-Tannakian-Closure TC-1d",
        "date_iso": "2026-06-06",
        "cutoffs": cutoffs,
        "results_by_cutoff": {},
    }

    for n_max in cutoffs:
        print(f"\n--- Cutoff n_max = {n_max} ---")
        panel = build_panel(n_max)
        t_c = time.time()

        v_unit = verify_omega_unit(n_max)
        unit_ok = v_unit["bit_exact"]

        tensor_pairs = [("T1", "J2"), ("J2", "J2"), ("J2", "J3"), ("J3", "J2")]
        tensor_results = []
        for n1, n2 in tensor_pairs:
            v = verify_omega_tensor_preservation(panel[n1], panel[n2])
            v["pair"] = [n1, n2]
            tensor_results.append(v)

        morph_names = ["f1", "f2", "f3", "f4"]
        ker_results = []
        for name in morph_names:
            v = verify_omega_preserves_kernel(panel[name])
            v["morphism"] = name
            ker_results.append(v)

        coker_results = []
        for name in morph_names:
            v = verify_omega_preserves_cokernel(panel[name])
            v["morphism"] = name
            coker_results.append(v)

        ds_pairs = [("T1", "T1"), ("T1", "J2"), ("J2", "J3")]
        ds_results = []
        for n1, n2 in ds_pairs:
            v = verify_omega_preserves_direct_sum(panel[n1], panel[n2])
            v["pair"] = [n1, n2]
            ds_results.append(v)

        T = trivial_rep(n_max, dim=2)
        id_T = RepMorphism(T, T, eye(2), label="id_T")
        zero_T = RepMorphism(T, T, Matrix.zeros(2, 2), label="zero_T")
        faithful_results = [
            verify_omega_faithful(id_T, zero_T),
            verify_omega_faithful(id_T, id_T),
            verify_omega_faithful(panel["f1"], panel["f1"]),
            verify_omega_faithful(panel["f2"], panel["f2"]),
        ]

        zero = (
            (1 if unit_ok else 0)
            + sum(1 for r in tensor_results if r["bit_exact"])
            + sum(1 for r in ker_results if r["bit_exact"])
            + sum(1 for r in coker_results if r["bit_exact"])
            + sum(1 for r in ds_results if r["bit_exact"])
            + sum(1 for r in faithful_results if r["bit_exact"])
        )
        total = (
            1 + len(tensor_results) + len(ker_results)
            + len(coker_results) + len(ds_results)
            + len(faithful_results)
        )

        print(f"  [1] omega(1) = Q:                {'OK' if unit_ok else 'FAIL'}")
        print(f"  [2] otimes-preservation:         {sum(1 for r in tensor_results if r['bit_exact'])}/{len(tensor_results)}")
        print(f"  [3] Kernel preservation:         {sum(1 for r in ker_results if r['bit_exact'])}/{len(ker_results)}")
        print(f"  [4] Cokernel preservation:       {sum(1 for r in coker_results if r['bit_exact'])}/{len(coker_results)}")
        print(f"  [5] Direct sum preservation:     {sum(1 for r in ds_results if r['bit_exact'])}/{len(ds_results)}")
        print(f"  [6] Faithfulness:                {sum(1 for r in faithful_results if r['bit_exact'])}/{len(faithful_results)}")
        print(f"  --> {zero} / {total} bit-exact zero residuals ({(time.time() - t_c):.2f}s)")

        grand_zero += zero
        grand_total += total
        persist["results_by_cutoff"][str(n_max)] = {
            "unit": v_unit,
            "tensor": tensor_results,
            "kernel": ker_results,
            "cokernel": coker_results,
            "direct_sum": ds_results,
            "faithful": faithful_results,
            "zero_residuals": zero,
            "total_identities": total,
        }

    elapsed = time.time() - t_global
    print("\n" + "=" * 72)
    print(f"Total bit-exact zero residuals : {grand_zero} / {grand_total}")
    print(f"Total wall time                : {elapsed:.2f}s")
    print("=" * 72)

    persist["totals"] = {"zero_residuals": grand_zero, "total_identities": grand_total}
    persist["wall_time_seconds"] = elapsed

    out_path = Path("debug/data/sprint_q5p_tc1d_fiber_functor.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(persist, f, indent=2, default=str)
    print(f"\nData written to: {out_path}")


if __name__ == "__main__":
    main()
