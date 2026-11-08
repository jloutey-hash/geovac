r"""
Sprint Q5'-Tannakian-Closure TC-1b --- symmetric monoidal structure on
$\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$.

Goal
----
Verify the symmetric monoidal axioms bit-exact:\ the tensor product
$\otimes_\mathbb{Q}$ with the diagonal Hopf action determined by the
v3.61.0 Track A abelian primitive coproduct $\Delta(x) = x \otimes 1 +
1 \otimes x$ makes $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$
a symmetric monoidal category. Closes the SECOND of four sprint-scale
Tannakian-closure prerequisites named by PS-4 (v3.69.0).

Six axiom panels at $n_{\max} \in \{2, 3\}$
-------------------------------------------
1. **Tensor diagonal action.** Independent reconstruction of
   $X_g^{M \otimes N} = X_g^M \otimes I_N + I_M \otimes X_g^N$ compared
   against `tensor_rep`'s output. Bit-exact at every non-zero generator,
   plus nilpotency check for each resulting endomorphism. 4 pairs $\times$
   $\le 2$ checks $=$ 8 per cutoff.

2. **Tensor functoriality.** For pairs of composable morphism pairs
   $(f, g)$ and $(f', g')$,
   $(f' \otimes g') \circ (f \otimes g) = (f' \circ f) \otimes (g' \circ g)$
   bit-exact. 3 test compositions per cutoff.

3. **Unitor naturality + intertwining.** Verify $\lambda_M$ and $\rho_M$
   are valid rep morphisms (intertwine the action) for $M \in \{T_1, J_2, J_3\}$.
   3 reps $\times$ 2 unitors $=$ 6 per cutoff.

4. **Associator intertwines.** For triples $(M, N, P) \in \{(T_1, J_2, J_3),
   (J_2, T_1, J_2)\}$, verify $\alpha_{M, N, P}$ intertwines the diagonal
   action bit-exact. 2 triples per cutoff.

5. **Braiding intertwines AND symmetric.** For pairs $(M, N) \in \{(T_1, J_2),
   (J_2, J_3), (J_2, J_2)\}$, verify $\sigma_{M, N}$ intertwines and
   $\sigma_{N, M} \circ \sigma_{M, N} = \mathrm{id}_{M \otimes N}$.
   3 pairs $\times$ 2 conditions $=$ 6 per cutoff.

6. **Coherence diagrams.** Pentagon on $(T_1, J_2, J_3, J_2)$ + triangle on
   $(J_2, J_3)$ + hexagon on $(J_2, J_3, T_1)$. 3 per cutoff.

Per cutoff: $8 + 3 + 6 + 2 + 6 + 3 = 28$. Two cutoffs: $56$ bit-exact zero
residuals expected.

Output
------
- ``debug/data/sprint_q5p_tc1b_monoidal.json``
- Console summary at end.

Discipline
----------
Bit-exact ``sympy.Rational`` / ``sympy.Integer`` throughout. No floats. No PSLQ.

References
----------
- TC-1a memo ``debug/sprint_q5p_tc1a_abelian_memo.md`` (substrate).
- v3.61.0 Track A memo ``debug/sprint_q5p_stage2_hopf_memo.md``
  (abelian primitive coproduct $\Delta(x) = x \otimes 1 + 1 \otimes x$).
- Deligne, P.; Milne, J. S. ``Tannakian categories'' (1982).
- Mac Lane, S. *Categories for the Working Mathematician* (1998), Ch. VII
  (monoidal categories), Ch. XI (symmetric monoidal coherence).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

from sympy import Integer, Matrix, eye as sp_eye

from geovac.tannakian import (
    FinDimRep,
    RepMorphism,
    associator,
    braiding,
    compose,
    left_unitor,
    right_unitor,
    tensor_morphism,
    tensor_rep,
    trivial_rep,
    unit_object,
    verify_associator_intertwines,
    verify_braiding_intertwines,
    verify_braiding_symmetric,
    verify_hexagon_coherence,
    verify_pentagon_coherence,
    verify_tensor_diagonal_action,
    verify_tensor_functoriality,
    verify_triangle_coherence,
    verify_unitor_intertwines,
)


# =====================================================================
# Panel substrate
# =====================================================================


def build_panel(n_max: int) -> Dict[str, Any]:
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
    # Morphisms for functoriality tests
    f1 = RepMorphism(T1, J2, Matrix([[1], [0]]), label="f1")
    f2 = RepMorphism(J2, T1, Matrix([[0, 1]]), label="f2")
    # f2 ∘ f1 = (0,1) * (1, 0)^T = 0 → zero morphism T1 → T1
    g1 = RepMorphism(T1, T2, Matrix([[1], [0]]), label="g1")
    g2 = RepMorphism(T2, T1, Matrix([[1, 0]]), label="g2")
    # g2 ∘ g1 = 1 → identity on T1
    id_J2 = RepMorphism(J2, J2, sp_eye(2), label="id_J2")
    id_J3 = RepMorphism(J3, J3, sp_eye(3), label="id_J3")
    return {
        "T1": T1, "T2": T2, "J2": J2, "J3": J3,
        "f1": f1, "f2": f2, "g1": g1, "g2": g2,
        "id_J2": id_J2, "id_J3": id_J3,
    }


# =====================================================================
# Axiom panels
# =====================================================================


def run_tensor_diagonal_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    pairs = [("T1", "J2"), ("J2", "J2"), ("J2", "J3"), ("J3", "J2")]
    results: List[Dict[str, Any]] = []
    for n1, n2 in pairs:
        M = panel[n1]
        N = panel[n2]
        v = verify_tensor_diagonal_action(M, N)
        v["pair"] = [n1, n2]
        results.append(v)
    return results


def run_tensor_functoriality_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    # Three composition test cases.
    # Case 1: f2 ∘ f1: T1 -> J2 -> T1, g2 ∘ g1: T1 -> T2 -> T1
    f1 = panel["f1"]
    f2 = panel["f2"]
    g1 = panel["g1"]
    g2 = panel["g2"]
    id_J2 = panel["id_J2"]
    id_J3 = panel["id_J3"]
    # Case 2: id_J2 ∘ id_J2 = id_J2, id_J3 ∘ id_J3 = id_J3 (trivial)
    # Case 3: f1 ∘ id_T1, id_J2 ∘ f1 (composition with identity)
    # Use these three:
    tests = [
        (f1, g1, f2, g2),  # (T1->J2,T1->T2) and (J2->T1,T2->T1)
        (id_J2, id_J3, id_J2, id_J3),  # identities
        (f1, id_J3, id_J2, id_J3),  # f1 ⊗ id, id ⊗ id
    ]
    results: List[Dict[str, Any]] = []
    for (f, g, f_prime, g_prime) in tests:
        v = verify_tensor_functoriality(f, g, f_prime, g_prime)
        results.append(v)
    return results


def run_unitor_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    for name in ["T1", "J2", "J3"]:
        M = panel[name]
        v = verify_unitor_intertwines(M)
        v["rep"] = name
        results.append(v)
    return results


def run_associator_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    triples = [("T1", "J2", "J3"), ("J2", "T1", "J2")]
    results: List[Dict[str, Any]] = []
    for n1, n2, n3 in triples:
        v = verify_associator_intertwines(panel[n1], panel[n2], panel[n3])
        v["triple"] = [n1, n2, n3]
        results.append(v)
    return results


def run_braiding_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    pairs = [("T1", "J2"), ("J2", "J3"), ("J2", "J2")]
    results: List[Dict[str, Any]] = []
    for n1, n2 in pairs:
        v_int = verify_braiding_intertwines(panel[n1], panel[n2])
        v_sym = verify_braiding_symmetric(panel[n1], panel[n2])
        results.append({
            "pair": [n1, n2],
            "intertwines": v_int["bit_exact"],
            "symmetric": v_sym["bit_exact"],
            "bit_exact": v_int["bit_exact"] and v_sym["bit_exact"],
        })
    return results


def run_coherence_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    # Pentagon
    v_p = verify_pentagon_coherence(
        panel["T1"], panel["J2"], panel["J3"], panel["J2"],
    )
    results.append({"diagram": "pentagon", **v_p})
    # Triangle
    v_t = verify_triangle_coherence(panel["J2"], panel["J3"])
    results.append({"diagram": "triangle", **v_t})
    # Hexagon
    v_h = verify_hexagon_coherence(panel["J2"], panel["J3"], panel["T1"])
    results.append({"diagram": "hexagon", **v_h})
    return results


def count_zero_residuals(
    diag_results: List[Dict[str, Any]],
    fun_results: List[Dict[str, Any]],
    unit_results: List[Dict[str, Any]],
    assoc_results: List[Dict[str, Any]],
    braid_results: List[Dict[str, Any]],
    coh_results: List[Dict[str, Any]],
) -> Tuple[int, int]:
    total = 0
    zero = 0
    # Tensor diagonal: 2 checks per pair (each pair has up to 1 mismatch reported + nilpotency)
    for r in diag_results:
        # Use 2 checks per pair as accounting (one for matrix match, one for nilpotency).
        total += 2
        if r["bit_exact"]:
            zero += 2
        elif len(r["mismatches"]) == 0:
            zero += 1
        # else: at least the nilpotency check could pass; conservatively count it
    # Functoriality: 1 per test
    for r in fun_results:
        total += 1
        if r["bit_exact"]:
            zero += 1
    # Unitor: 2 checks per rep (left + right)
    for r in unit_results:
        total += 2
        if r["bit_exact"]:
            zero += 2
        else:
            if r["left_unitor_intertwines"]:
                zero += 1
            if r["right_unitor_intertwines"]:
                zero += 1
    # Associator: 1 per triple
    for r in assoc_results:
        total += 1
        if r["bit_exact"]:
            zero += 1
    # Braiding: 2 per pair (intertwines + symmetric)
    for r in braid_results:
        total += 2
        if r["intertwines"]:
            zero += 1
        if r["symmetric"]:
            zero += 1
    # Coherence: 1 per diagram
    for r in coh_results:
        total += 1
        if r["bit_exact"]:
            zero += 1
    return zero, total


# =====================================================================
# Main driver
# =====================================================================


def main() -> None:
    print("=" * 72)
    print("Sprint Q5'-Tannakian-Closure  TC-1b")
    print("Symmetric monoidal structure on Rep_fin(H_GV)")
    print("=" * 72)

    t_global = time.time()
    cutoffs = [2, 3]
    grand_zero = 0
    grand_total = 0
    persist: Dict[str, Any] = {
        "sprint": "Q5'-Tannakian-Closure TC-1b",
        "date_iso": "2026-06-06",
        "cutoffs": cutoffs,
        "results_by_cutoff": {},
    }

    for n_max in cutoffs:
        print(f"\n--- Cutoff n_max = {n_max} ---")
        t_c = time.time()
        panel = build_panel(n_max)
        diag = run_tensor_diagonal_panel(panel)
        fun = run_tensor_functoriality_panel(panel)
        unit_p = run_unitor_panel(panel)
        assoc = run_associator_panel(panel)
        braid = run_braiding_panel(panel)
        coh = run_coherence_panel(panel)
        print(f"  [1] Tensor diagonal action: {sum(r['bit_exact'] for r in diag)}/{len(diag)} pairs pass")
        print(f"  [2] Tensor functoriality:    {sum(r['bit_exact'] for r in fun)}/{len(fun)} tests pass")
        print(f"  [3] Unitor intertwines:       {sum(r['bit_exact'] for r in unit_p)}/{len(unit_p)} reps pass")
        print(f"  [4] Associator intertwines:   {sum(r['bit_exact'] for r in assoc)}/{len(assoc)} triples pass")
        print(f"  [5] Braiding intertwines+sym: {sum(r['bit_exact'] for r in braid)}/{len(braid)} pairs pass")
        print(f"  [6] Coherence diagrams:       {sum(r['bit_exact'] for r in coh)}/{len(coh)} diagrams pass")
        zero, total = count_zero_residuals(diag, fun, unit_p, assoc, braid, coh)
        print(f"  --> {zero} / {total} bit-exact zero residuals ({(time.time() - t_c):.2f}s)")
        grand_zero += zero
        grand_total += total
        persist["results_by_cutoff"][str(n_max)] = {
            "tensor_diagonal": diag,
            "functoriality": fun,
            "unitors": unit_p,
            "associator": assoc,
            "braiding": braid,
            "coherence": coh,
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

    out_path = Path("debug/data/sprint_q5p_tc1b_monoidal.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(persist, f, indent=2, default=str)
    print(f"\nData written to: {out_path}")


if __name__ == "__main__":
    main()
