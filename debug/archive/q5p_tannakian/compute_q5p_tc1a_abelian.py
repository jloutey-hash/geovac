r"""
Sprint Q5'-Tannakian-Closure TC-1a --- the rep category
$\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}}(n_{\max}))$ is
abelian, bit-exact at $n_{\max} \in \{2, 3\}$.

Goal
----
First sub-track of the Tannakian Reconstruction Foundation sprint (TC-1).
PS-4 (v3.69.0) named four Deligne--Milne 1982 prerequisites for
Tannakian closure;\ three of them (abelian / tensor / rigidity) are
sprint-scale bookkeeping on top of the v3.61.0 Track A abelian primitive
substrate. TC-1a closes the first:\ verify bit-exact at small reps and
small cutoffs that $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$
is abelian.

The rep category
----------------
For an abelian primitive Hopf algebra $\mathcal{H}_{\mathrm{GV}}(n_{\max})
= \mathrm{Sym}_{\mathbb{Q}}(V_{n_{\max}}) = \mathcal{O}(\mathbb{G}_a^{3 N(n_{\max})})$,
a finite-dim rational rep is a finite-dim $\mathbb{Q}$-vector space $M$
equipped with $3 N(n_{\max})$ pairwise commuting nilpotent endomorphisms
$\{X_g\}_g$ (one per primitive generator $g = (n, l, s)$). A morphism
$f: M \to M'$ is a $\mathbb{Q}$-linear map satisfying
$f \circ X_g^M = X_g^{M'} \circ f$ for every $g$.

Substrate panel
---------------
Five reps at each cutoff $n_{\max} \in \{2, 3\}$:

R0  = zero rep (dim 0)
T1  = trivial 1-dim rep (dim 1, all endos = 0)
T2  = trivial 2-dim rep (dim 2, all endos = 0)
J2  = "Jordan" 2-dim rep with $X_{(1, 0, 0)} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$,
        all other endos = 0
J3  = "Jordan" 3-dim rep with $X_{(1, 0, 0)} = $ canonical nilpotent
        $\begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$,
        all other endos = 0

Five morphisms:

f0  = 0 -> T1 (the unique zero morphism)
f1  = T1 -> J2 mapping $1 \mapsto e_1$ (the kernel of $X_{(1,0,0)}$)
f2  = J2 -> T1 mapping $e_2 \mapsto 1, e_1 \mapsto 0$ (the cokernel of $f_1$)
f3  = identity on J2
f4  = T1 -> T2 mapping $1 \mapsto e_1$ (test injection)

Bit-exact panel
---------------
For each cutoff $n_{\max} \in \{2, 3\}$, verify:

1. Zero object axioms ($R0$ has unique morphisms to/from every $R$):\ 5 reps $\times$ 2 directions = 10 identities per cutoff.
2. Direct sum universal property:\ for representative pairs of reps,
   $\pi_1 \circ \iota_1 = \mathrm{id}_{R_1}$, $\pi_2 \circ \iota_2 = \mathrm{id}_{R_2}$,
   $\pi_1 \circ \iota_2 = 0$, $\pi_2 \circ \iota_1 = 0$,
   $\iota_1 \circ \pi_1 + \iota_2 \circ \pi_2 = \mathrm{id}_S$.
   5 identities per pair $\times$ 3 pairs $= 15$ per cutoff.
3. Kernel universal property:\ for each test morphism, $f \circ \iota_{\ker} = 0$
   and $\iota_{\ker}$ is mono. 2 conditions $\times$ 4 morphisms = 8 per cutoff.
4. Cokernel universal property:\ $\pi_{\mathrm{coker}} \circ f = 0$ and
   $\pi_{\mathrm{coker}}$ is epi. 8 per cutoff.
5. Mono $=$ ker(coker):\ for each mono test morphism, kernel-of-cokernel
   recovers the original up to canonical isomorphism. 3 conditions $\times$
   2 test monos = 6 per cutoff.
6. Epi $=$ coker(ker):\ symmetric, 6 per cutoff.

Total per cutoff: $10 + 15 + 8 + 8 + 6 + 6 = 53$. Two cutoffs $\times 53$
$=$ 106 bit-exact identities expected.

Output
------
- ``debug/data/sprint_q5p_tc1a_abelian.json``
- Console summary at end.

Discipline
----------
Bit-exact ``sympy.Rational`` / ``sympy.Integer`` throughout. No floats.
No PSLQ. No transcendentals introduced.

References
----------
- ``geovac/tannakian.py`` (the substrate built this sprint).
- ``geovac/pro_system.py`` PS-2 + PS-3 (primitive_generators, n_primitive_generators).
- v3.61.0 Track A memo ``debug/sprint_q5p_stage2_hopf_memo.md``
  (abelian primitive Hopf substrate).
- v3.69.0 PS-4 memo ``debug/sprint_q5p_ps4_endo_rigidity_memo.md``
  (Tannakian-readiness gap-list).
- Deligne, P.; Milne, J. S. ``Tannakian categories'' (1982).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

from sympy import Integer, Matrix, Rational, eye as sp_eye

from geovac.tannakian import (
    FinDimRep,
    RepMorphism,
    cokernel,
    compose,
    direct_sum,
    kernel,
    trivial_rep,
    verify_cokernel_universal_property,
    verify_direct_sum_universal_property,
    verify_epi_eq_coker_ker,
    verify_kernel_universal_property,
    verify_mono_eq_ker_coker,
    verify_zero_object_axiom,
    zero_rep,
)


# =====================================================================
# Build the substrate panel at a given cutoff
# =====================================================================


def build_panel(n_max: int) -> Dict[str, Any]:
    """Build the five reps and five morphisms at cutoff `n_max`."""
    # Reps
    R0 = zero_rep(n_max)
    R0.label = "R0_zero"
    T1 = trivial_rep(n_max, dim=1)
    T1.label = "T1_trivial1"
    T2 = trivial_rep(n_max, dim=2)
    T2.label = "T2_trivial2"
    J2 = FinDimRep(
        n_max=n_max,
        dim=2,
        endos={(1, 0, 0): Matrix([[0, 1], [0, 0]])},
        label="J2_jordan",
    )
    J3 = FinDimRep(
        n_max=n_max,
        dim=3,
        endos={(1, 0, 0): Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])},
        label="J3_jordan",
    )
    # Morphisms
    f0 = RepMorphism(
        source=R0,
        target=T1,
        matrix=Matrix.zeros(1, 0),
        label="f0_zero",
        validate=False,
    )
    # f1: T1 -> J2 mapping 1 -> e_1 (a column vector [1, 0])
    f1 = RepMorphism(
        source=T1,
        target=J2,
        matrix=Matrix([[1], [0]]),
        label="f1_T1_to_J2",
    )
    # f2: J2 -> T1 mapping e_2 -> 1, e_1 -> 0 (a row vector [0, 1])
    f2 = RepMorphism(
        source=J2,
        target=T1,
        matrix=Matrix([[0, 1]]),
        label="f2_J2_to_T1",
    )
    # f3: identity on J2
    f3 = RepMorphism(
        source=J2,
        target=J2,
        matrix=sp_eye(2),
        label="f3_id_J2",
    )
    # f4: T1 -> T2 injection 1 -> e_1
    f4 = RepMorphism(
        source=T1,
        target=T2,
        matrix=Matrix([[1], [0]]),
        label="f4_T1_to_T2",
    )
    return {
        "reps": {"R0": R0, "T1": T1, "T2": T2, "J2": J2, "J3": J3},
        "morphisms": {"f0": f0, "f1": f1, "f2": f2, "f3": f3, "f4": f4},
    }


# =====================================================================
# Axiom panels
# =====================================================================


def run_zero_object_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Zero object axiom: for each rep R, the unique morphisms 0 -> R
    and R -> 0 exist (bit-exact construction).
    """
    results: List[Dict[str, Any]] = []
    for name, R in panel["reps"].items():
        v = verify_zero_object_axiom(R)
        v["rep"] = name
        results.append(v)
    return results


def run_direct_sum_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Direct sum universal property on a representative set of pairs."""
    pairs = [
        ("T1", "T1"),
        ("T1", "J2"),
        ("J2", "J3"),
    ]
    results: List[Dict[str, Any]] = []
    for n1, n2 in pairs:
        R1 = panel["reps"][n1]
        R2 = panel["reps"][n2]
        S, iota1, iota2, pi1, pi2 = direct_sum(R1, R2)
        v = verify_direct_sum_universal_property(R1, R2, S, iota1, iota2, pi1, pi2)
        v["pair"] = [n1, n2]
        results.append(v)
    return results


def run_kernel_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Kernel universal property on each non-zero test morphism."""
    morph_names = ["f1", "f2", "f3", "f4"]
    results: List[Dict[str, Any]] = []
    for name in morph_names:
        f = panel["morphisms"][name]
        K, iota = kernel(f)
        v = verify_kernel_universal_property(f, K, iota)
        v["morphism"] = name
        results.append(v)
    return results


def run_cokernel_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Cokernel universal property on each non-zero test morphism."""
    morph_names = ["f1", "f2", "f3", "f4"]
    results: List[Dict[str, Any]] = []
    for name in morph_names:
        f = panel["morphisms"][name]
        C, pi = cokernel(f)
        v = verify_cokernel_universal_property(f, C, pi)
        v["morphism"] = name
        results.append(v)
    return results


def run_mono_eq_ker_coker_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Mono = ker(coker) on injective test morphisms (f1, f4)."""
    morph_names = ["f1", "f4"]
    results: List[Dict[str, Any]] = []
    for name in morph_names:
        f = panel["morphisms"][name]
        v = verify_mono_eq_ker_coker(f)
        v["morphism"] = name
        results.append(v)
    return results


def run_epi_eq_coker_ker_panel(panel: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Epi = coker(ker) on surjective test morphisms (f2, f3)."""
    morph_names = ["f2", "f3"]
    results: List[Dict[str, Any]] = []
    for name in morph_names:
        f = panel["morphisms"][name]
        v = verify_epi_eq_coker_ker(f)
        v["morphism"] = name
        results.append(v)
    return results


# =====================================================================
# Main driver
# =====================================================================


def count_bit_exact_results(panels: Dict[str, List[Dict[str, Any]]]) -> Tuple[int, int]:
    """Count bit-exact axiom identities across all panels."""
    total = 0
    zero_residuals = 0
    for axiom_name, panel_results in panels.items():
        if axiom_name == "zero_object":
            # Per-rep: 2 directions (R0 -> R, R -> R0); each construction is direct.
            for r in panel_results:
                total += 2
                if r["bit_exact"]:
                    zero_residuals += 2
        elif axiom_name == "direct_sum":
            # 5 identities per pair
            for r in panel_results:
                total += 5
                if r["bit_exact"]:
                    zero_residuals += 5
        elif axiom_name == "kernel":
            # 2 conditions per morphism
            for r in panel_results:
                total += 2
                if r["bit_exact"]:
                    zero_residuals += 2
        elif axiom_name == "cokernel":
            for r in panel_results:
                total += 2
                if r["bit_exact"]:
                    zero_residuals += 2
        elif axiom_name == "mono_eq_ker_coker":
            # 3 conditions per mono test
            for r in panel_results:
                total += 3
                if r["bit_exact"]:
                    zero_residuals += 3
        elif axiom_name == "epi_eq_coker_ker":
            for r in panel_results:
                total += 3
                if r["bit_exact"]:
                    zero_residuals += 3
    return zero_residuals, total


def main() -> None:
    print("=" * 72)
    print("Sprint Q5'-Tannakian-Closure  TC-1a")
    print("Rep_fin(H_GV(n_max)) is abelian, bit-exact at n_max in {2, 3}")
    print("=" * 72)

    t_global = time.time()

    cutoffs = [2, 3]
    all_results: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}

    for n_max in cutoffs:
        print(f"\n--- Cutoff n_max = {n_max} ---")
        t_cutoff = time.time()
        panel = build_panel(n_max)
        results: Dict[str, List[Dict[str, Any]]] = {}

        results["zero_object"] = run_zero_object_panel(panel)
        print(f"  [1] Zero object axiom: 5 reps, "
              f"all bit-exact = "
              f"{all(r['bit_exact'] for r in results['zero_object'])}")

        results["direct_sum"] = run_direct_sum_panel(panel)
        print(f"  [2] Direct sum universal property: 3 pairs, "
              f"all bit-exact = "
              f"{all(r['bit_exact'] for r in results['direct_sum'])}")

        results["kernel"] = run_kernel_panel(panel)
        print(f"  [3] Kernel universal property: 4 morphisms, "
              f"all bit-exact = "
              f"{all(r['bit_exact'] for r in results['kernel'])}")

        results["cokernel"] = run_cokernel_panel(panel)
        print(f"  [4] Cokernel universal property: 4 morphisms, "
              f"all bit-exact = "
              f"{all(r['bit_exact'] for r in results['cokernel'])}")

        results["mono_eq_ker_coker"] = run_mono_eq_ker_coker_panel(panel)
        print(f"  [5] Mono = ker(coker): 2 morphisms, "
              f"all bit-exact = "
              f"{all(r['bit_exact'] for r in results['mono_eq_ker_coker'])}")

        results["epi_eq_coker_ker"] = run_epi_eq_coker_ker_panel(panel)
        print(f"  [6] Epi = coker(ker): 2 morphisms, "
              f"all bit-exact = "
              f"{all(r['bit_exact'] for r in results['epi_eq_coker_ker'])}")

        zero_resid, total = count_bit_exact_results(results)
        print(f"  --> {zero_resid} / {total} bit-exact zero residuals "
              f"({(time.time() - t_cutoff):.2f}s)")

        all_results[str(n_max)] = results

    elapsed_total = time.time() - t_global

    # Aggregate
    grand_zero = 0
    grand_total = 0
    for n_max_str, results in all_results.items():
        z, t = count_bit_exact_results(results)
        grand_zero += z
        grand_total += t

    print(f"\n{'=' * 72}")
    print(f"Total bit-exact zero residuals : {grand_zero} / {grand_total}")
    print(f"Total wall time                : {elapsed_total:.2f}s")
    print(f"{'=' * 72}")

    # Persist
    data: Dict[str, Any] = {
        "sprint": "Q5'-Tannakian-Closure TC-1a",
        "date_iso": "2026-06-06",
        "cutoffs": cutoffs,
        "results_by_cutoff": {
            str(n_max): {
                axiom: [
                    {k: (str(v) if isinstance(v, (list, tuple, dict, str, bool, int)) else str(v))
                     for k, v in r.items()}
                    for r in results
                ]
                for axiom, results in all_results[str(n_max)].items()
            }
            for n_max in cutoffs
        },
        "totals": {
            "zero_residuals": grand_zero,
            "total_identities": grand_total,
        },
        "wall_time_seconds": elapsed_total,
    }

    out_path = Path("debug/data/sprint_q5p_tc1a_abelian.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"\nData written to: {out_path}")


if __name__ == "__main__":
    main()
