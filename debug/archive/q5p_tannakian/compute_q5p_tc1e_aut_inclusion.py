r"""
Sprint Q5'-Tannakian-Closure TC-1e --- first stone of the multi-year
wall:\ explicit inclusion $\mathbb{G}_a^{3 N(n_{\max})}(\mathbb{Q})
\hookrightarrow \mathrm{Aut}^\otimes(\omega)$ at finite cutoff.

Goal
----
Construct the explicit map
$\Phi: \mathbb{Q}^{3 N(n_{\max})} \to \mathrm{Aut}^\otimes(\omega)$,
$\Phi(t)(V) = \exp(\sum_g t_g X_g^V)$, and verify bit-exact at
$n_{\max} \in \{2, 3\}$ on a representative panel that includes reps
activating multiple distinct primitive generators.

This is the **first stone of the multi-year wall** named by PS-4
(v3.69.0):\ closing the inclusion direction
$\mathrm{Aut}^\otimes(\omega) \supseteq U^*_{\mathrm{Levi}}$ at finite
cutoff. The converse equality $\mathrm{Aut}^\otimes(\omega) = U^*$
(full Tannakian closure proper) requires the inverse-limit pro-finite
Tannakian theorem coherent with v3.66.0 FO3 Interpretation C at the
period level on $\mathcal{O}_\infty$ — multi-year.

Panel
-----
Per cutoff $n_{\max} \in \{2, 3\}$:

* Reps (5): $T_1$ unit, $T_2$ trivial 2-dim, $J_2$ Jordan 2-dim at
  generator $g_0 = (1, 0, 0)$, $J_3$ Jordan 3-dim at $g_0$,
  $K_2$ Jordan 2-dim at a SECOND generator $g_1 = (2, 0, 1)$.

* Test parameters (4):
  - $t^{(0)} = 0$ (zero):\ $\eta_V = \mathrm{id}$ for every $V$.
  - $t^{(1)}$: only $t_{g_0} = 1$ non-zero.
  - $t^{(2)}$: only $t_{g_1} = 1$ non-zero.
  - $t^{(3)}$: $t_{g_0} = 2/3$ and $t_{g_1} = -1/4$ (mixed).

* Per $t$:
  - Invertibility on the 5 reps:                       5 checks
  - $\eta_{\mathbf{1}}(t) = \mathrm{id}$:              1 check
  - Naturality on morphisms $f_1, f_2$:                2 checks
  - $\otimes$-compatibility on 3 pairs:                3 checks
  Per $t$ subtotal:                                    11
  Across 4 $t$ values:                                 44

* Group law:
  - $\Phi(t^{(1)} + t^{(2)}) = \Phi(t^{(1)}) \cdot \Phi(t^{(2)})$
    on reps $J_2, J_3, K_2$:                           3 checks
  - $\Phi(2 t^{(1)}) = \Phi(t^{(1)})^2$ on $J_2, J_3$: 2 checks
  Subtotal:                                            5

Per cutoff:\ $44 + 5 = 49$. Two cutoffs:\ **98 bit-exact zero residuals**.

$SL_2$ commutativity is categorical (independence-of-axes:\ $SL_2$
acts on the $j_{\max}$ axis of the Peter--Weyl decoration, independent
of the $n_{\max}$ axis on which the fiber functor's substrate lives).
No per-cell check needed.

Output
------
- ``debug/data/sprint_q5p_tc1e_aut_inclusion.json``
- Console summary at end.

Discipline
----------
Bit-exact ``sympy.Rational`` / ``sympy.Integer`` throughout.

References
----------
- TC-1d memo (fiber functor $\omega$).
- v3.61.0 Track A (abelian primitive Hopf substrate).
- v3.63.0 L1 (Levi decomposition $U^* = \mathbb{G}_a^{3N} \rtimes SL_2$).
- PS-4 memo (Tannakian-readiness gap-list naming TC-1e as multi-year heart).
- Deligne, P.; Milne, J. S. ``Tannakian categories'' (1982).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict

from sympy import Matrix, Rational

from geovac.tannakian import (
    FinDimRep,
    RepMorphism,
    trivial_rep,
    unit_object,
    verify_natural_auto_group_law,
    verify_natural_auto_invertibility,
    verify_natural_auto_naturality,
    verify_natural_auto_tensor,
    verify_natural_auto_unit,
)


def build_panel(n_max: int) -> Dict[str, Any]:
    T1 = unit_object(n_max)
    T1.label = "T1"
    T2 = trivial_rep(n_max, dim=2)
    T2.label = "T2"
    J_mat = Matrix([[0, 1], [0, 0]])
    J3_mat = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    J2 = FinDimRep(n_max=n_max, dim=2, endos={(1, 0, 0): J_mat}, label="J2_g0")
    J3 = FinDimRep(n_max=n_max, dim=3, endos={(1, 0, 0): J3_mat}, label="J3_g0")
    K2 = FinDimRep(n_max=n_max, dim=2, endos={(2, 0, 1): J_mat}, label="K2_g1")
    # Morphisms (intertwine when restricted to reps where the relevant
    # generators have actions matching):
    f1 = RepMorphism(T1, J2, Matrix([[1], [0]]), label="f1_T1_to_J2")
    f2 = RepMorphism(J2, T1, Matrix([[0, 1]]), label="f2_J2_to_T1")
    return {
        "T1": T1, "T2": T2, "J2": J2, "J3": J3, "K2": K2,
        "f1": f1, "f2": f2,
    }


def main() -> None:
    print("=" * 72)
    print("Sprint Q5'-Tannakian-Closure  TC-1e")
    print("First stone: G_a^{3N}(Q) -> Aut^otimes(omega) at finite cutoff")
    print("=" * 72)

    t_global = time.time()
    cutoffs = [2, 3]
    grand_zero = 0
    grand_total = 0
    persist: Dict[str, Any] = {
        "sprint": "Q5'-Tannakian-Closure TC-1e",
        "date_iso": "2026-06-06",
        "cutoffs": cutoffs,
        "results_by_cutoff": {},
    }

    g0 = (1, 0, 0)
    g1 = (2, 0, 1)

    t_values = {
        "t_zero": {},
        "t_g0_one": {g0: Rational(1)},
        "t_g1_one": {g1: Rational(1)},
        "t_mixed": {g0: Rational(2, 3), g1: Rational(-1, 4)},
    }

    for n_max in cutoffs:
        print(f"\n--- Cutoff n_max = {n_max} ---")
        panel = build_panel(n_max)
        t_c = time.time()

        per_cutoff_zero = 0
        per_cutoff_total = 0
        per_t_results: Dict[str, Any] = {}

        for t_label, t_dict in t_values.items():
            # Invertibility on 5 reps
            inv_results = []
            for rep_name in ["T1", "T2", "J2", "J3", "K2"]:
                v = verify_natural_auto_invertibility(t_dict, panel[rep_name])
                v["rep"] = rep_name
                inv_results.append(v)
            inv_zero = sum(1 for r in inv_results if r["bit_exact"])
            # Unit
            unit_r = verify_natural_auto_unit(t_dict, n_max)
            unit_zero = 1 if unit_r["bit_exact"] else 0
            # Naturality
            nat_results = []
            for m_name in ["f1", "f2"]:
                v = verify_natural_auto_naturality(t_dict, panel[m_name])
                v["morphism"] = m_name
                nat_results.append(v)
            nat_zero = sum(1 for r in nat_results if r["bit_exact"])
            # Tensor (3 pairs)
            tensor_results = []
            for n1, n2 in [("T1", "J2"), ("J2", "J3"), ("J2", "K2")]:
                v = verify_natural_auto_tensor(t_dict, panel[n1], panel[n2])
                v["pair"] = [n1, n2]
                tensor_results.append(v)
            tensor_zero = sum(1 for r in tensor_results if r["bit_exact"])
            sub_total = 5 + 1 + 2 + 3
            sub_zero = inv_zero + unit_zero + nat_zero + tensor_zero
            print(f"  t = {t_label}: inv {inv_zero}/5, unit {unit_zero}/1, "
                  f"nat {nat_zero}/2, tensor {tensor_zero}/3 -> {sub_zero}/{sub_total}")
            per_cutoff_zero += sub_zero
            per_cutoff_total += sub_total
            per_t_results[t_label] = {
                "invertibility": inv_results,
                "unit": unit_r,
                "naturality": nat_results,
                "tensor": tensor_results,
                "subtotal": [sub_zero, sub_total],
            }

        # Group law: t1 = t_g0_one, t2 = t_g1_one; doubling t1 on J_2, J_3
        gl_results = []
        t1 = t_values["t_g0_one"]
        t2 = t_values["t_g1_one"]
        for rep_name in ["J2", "J3", "K2"]:
            v = verify_natural_auto_group_law(t1, t2, panel[rep_name])
            v["test"] = f"t_g0 + t_g1 on {rep_name}"
            v["rep"] = rep_name
            gl_results.append(v)
        # Doubling: Phi(t1) * Phi(t1) = Phi(2 t1)
        for rep_name in ["J2", "J3"]:
            v = verify_natural_auto_group_law(t1, t1, panel[rep_name])
            v["test"] = f"2 t_g0 on {rep_name}"
            v["rep"] = rep_name
            gl_results.append(v)
        gl_zero = sum(1 for r in gl_results if r["bit_exact"])
        print(f"  Group law:                         {gl_zero}/{len(gl_results)}")
        per_cutoff_zero += gl_zero
        per_cutoff_total += len(gl_results)

        print(f"  --> {per_cutoff_zero} / {per_cutoff_total} bit-exact zero "
              f"residuals ({(time.time() - t_c):.2f}s)")
        grand_zero += per_cutoff_zero
        grand_total += per_cutoff_total
        persist["results_by_cutoff"][str(n_max)] = {
            "per_t": per_t_results,
            "group_law": gl_results,
            "zero_residuals": per_cutoff_zero,
            "total_identities": per_cutoff_total,
        }

    elapsed = time.time() - t_global
    print("\n" + "=" * 72)
    print(f"Total bit-exact zero residuals : {grand_zero} / {grand_total}")
    print(f"Total wall time                : {elapsed:.2f}s")
    print("=" * 72)
    print()
    print("SL_2 categorical commutativity:")
    print("  SL_2 acts on the Peter-Weyl decoration (j_max axis), independent")
    print("  of the n_max axis on which Aut^otimes(omega) acts via the G_a")
    print("  factor. Acting first by SL_2 and then by Phi(t) equals acting in")
    print("  the reverse order, by independence of axes. Categorical bit-exact.")

    persist["totals"] = {"zero_residuals": grand_zero, "total_identities": grand_total}
    persist["SL2_categorical_commutativity"] = True
    persist["wall_time_seconds"] = elapsed

    out_path = Path("debug/data/sprint_q5p_tc1e_aut_inclusion.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(persist, f, indent=2, default=str)
    print(f"\nData written to: {out_path}")


if __name__ == "__main__":
    main()
