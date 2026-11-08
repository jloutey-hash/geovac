r"""
Sprint Q5'-ProSystem-Lockdown PS-2 --- $U^*$-action compatibility with
the closed-form transition maps from PS-1 on the truncated
Camporesi--Higuchi pro-system, extended to $n_{\max} = 5$.

Goal
----
Verify that the Levi-decomposed cosmic-Galois group $U^*_{\mathrm{GeoVac}}
= \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ (v3.63.0 L1) acts compatibly
with PS-1's transition maps $P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k$
on every cocycle class in the Interpretation-C-closed family (v3.66.0
FO3:\ $\chi, \eta, F(s)$), across the entire $n_{\max} \le 5$ panel.

Three sub-verifications
-----------------------
1. **Hopf-hom cofiltered axiom.** Build $\Phi_{m, k}: \mathcal{H}_{\mathrm{GV}}(m)
   \to \mathcal{H}_{\mathrm{GV}}(k)$ as the block-diagonal three-Mellin-slot
   lift of PS-1's $P_{m, k}$, and verify $\Phi_{m, k} = \Phi_{n, k} \cdot
   \Phi_{m, n}$ bit-exact for every triple $(m, n, k)$ with $1 \le k < n < m
   \le 5$.  Structurally identical to PS-1's algebra-level cofiltered axiom
   replicated three times (one per Mellin slot $k \in \{0, 1, 2\}$);
   verified independently here.

2. **$\mathbb{G}_a$ generator compatibility.** For each of the $3 N(m)$
   translation generators $e_{(n, l), s}$ of $\mathbb{G}_a^{3 N(m)}$ and
   for every pair $(m, k)$ with $1 \le k < m \le 5$, verify the
   sector-locality rule:\ $\Phi_{m, k}$ sends $e_{(n, l), s}^{(m)} \mapsto
   e_{(n, l), s}^{(k)}$ if $n \le k$, and to zero otherwise.  Total
   generator checks across the panel:\ $\sum_{(m, k)} 3 N(m) = 435$.

3. **Class-level $U^*$-action compatibility on $\chi, \eta$.** Verify
   $P_{m, k}(U^* \cdot \psi^{(m)}) = U^* \cdot P_{m, k}(\psi^{(m)}) =
   \psi^{(k)}$ bit-exact for $\psi \in \{\chi, \eta\}$ across all
   10 pairs $1 \le k < m \le 5$.  Interpretation C (v3.66.0 FO3) fixes
   $U^*$ to act trivially on depth-0 classes; this routine is the
   pro-system-level bookkeeping that the triviality survives the
   transition.

The $SL_2$ factor (v3.63.0 L1) acts on the Peter--Weyl decoration which
is $n_{\max}$-independent;\ its commutativity with $P_{m, k}$ is a
categorical statement (independent factors commute structurally), not
something requiring a per-cell bit-exact check.

Continuum-side $U^*$-action
---------------------------
The non-trivial part of Interpretation C lives on the M3 (odd-$\zeta$)
slot of the continuum Mellin lift $F(s)$ (v3.66.0 FO2).  $F(s)$ is an
inverse-limit object, not a finite-cutoff cocycle class, and is
handled by PS-3 (inverse limit and extended $U^*$-action on the
limit object).  PS-2 verifies the finite-cutoff substrate;\ PS-3 lifts
to the continuum.

Bit-exact panel (all required to be zero residual)
---------------------------------------------------
- Hopf-hom cofiltered axiom:\ 10 triples in $\{1, \ldots, 5\}$, each at
  matrix dimension $3 N(k) \times 3 N(m)$.
- $\mathbb{G}_a$ generator compatibility:\ 435 generator checks total.
- Class-level $U^*$-action on $\chi$:\ 10 pairs, two boolean conditions
  per pair (lhs == rhs, lhs == $\psi^{(k)}$).
- Class-level $U^*$-action on $\eta$:\ 10 pairs, two boolean conditions.

Total bit-exact zero residuals expected: $10 + 435 + 20 + 20 = 485$.

Output
------
- ``debug/data/sprint_q5p_ps2_ustar_compatibility.json``
- Console summary at end.

Discipline
----------
Bit-exact ``sympy.Rational`` / ``sympy.Integer`` throughout. No floats.
No PSLQ.

References
----------
- PS-1 memo ``debug/sprint_q5p_ps1_transitions_memo.md``.
- v3.61.0 Track A memo ``debug/sprint_q5p_stage2_hopf_memo.md``
  (Hopf substrate, abelian primitive).
- v3.63.0 L1 memo ``debug/sprint_q5p_levi_synthesis_memo.md``
  (Levi decomposition).
- v3.66.0 FO3 memo ``debug/sprint_q5p_fo2_fo3_mt_period_memo.md``
  (Interpretation C of $U^*$-action).
- ``geovac/pro_system.py`` (PS-1 substrate + PS-2 additions).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

from sympy import Integer, Matrix, zeros as sp_zeros

from geovac.spectral_triple import FockSpectralTriple
from geovac.pro_system import (
    HopfTransition,
    MELLIN_SLOTS,
    N_sectors,
    TransitionMap,
    n_primitive_generators,
    primitive_generators,
    sectors_at_cutoff,
    verify_Ga_generator_compatibility,
    verify_class_action_compatibility,
    verify_hopf_cofiltered_axiom,
)


# =====================================================================
# Class extraction (re-using PS-1's helpers)
# =====================================================================


def _sector_idempotent(st: FockSpectralTriple, sector_idx: int) -> Matrix:
    N = st.dim_H
    M = sp_zeros(N, N)
    for i in range(N):
        if st._state_to_sector[i] == sector_idx:
            M[i, i] = Integer(1)
    return M


def jlo_class_per_sector(st: FockSpectralTriple) -> Dict[Tuple[int, int], int]:
    gamma = st.grading
    out: Dict[Tuple[int, int], int] = {}
    for s_idx, (n, l) in enumerate(st.sectors):
        e_s = _sector_idempotent(st, s_idx)
        out[(n, l)] = int((gamma * e_s).trace())
    return out


def cm_eta_class_per_sector(st: FockSpectralTriple) -> Dict[Tuple[int, int], int]:
    gamma = st.grading
    D = st.dirac_operator
    out: Dict[Tuple[int, int], int] = {}
    for s_idx, (n, l) in enumerate(st.sectors):
        e_s = _sector_idempotent(st, s_idx)
        out[(n, l)] = int((gamma * D * e_s).trace())
    return out


def class_vector(
    cls: Dict[Tuple[int, int], int],
    sectors: List[Tuple[int, int]],
) -> List[int]:
    return [cls[s] for s in sectors]


# =====================================================================
# Main driver
# =====================================================================


def main() -> None:
    print("=" * 72)
    print("Sprint Q5'-ProSystem-Lockdown  PS-2")
    print("U* compatibility with transitions, n_max = 5 extension")
    print("=" * 72)

    t_global = time.time()

    cutoffs = [1, 2, 3, 4, 5]
    triples = [(m, n, k) for m in cutoffs
               for n in range(2, m)
               for k in range(1, n)
               if k < n < m]
    pairs = [(m, k) for m in cutoffs for k in cutoffs if 1 <= k < m]

    # -----------------------------------------------------------------
    # [1] Construct spectral triples and extract chi, eta classes.
    # -----------------------------------------------------------------
    print("\n[1] Constructing FockSpectralTriple and class vectors:")
    triples_st: Dict[int, FockSpectralTriple] = {}
    chi_classes: Dict[int, Dict[Tuple[int, int], int]] = {}
    eta_classes: Dict[int, Dict[Tuple[int, int], int]] = {}
    for nm in cutoffs:
        t0 = time.time()
        st = FockSpectralTriple(n_max=nm)
        triples_st[nm] = st
        chi_classes[nm] = jlo_class_per_sector(st)
        eta_classes[nm] = cm_eta_class_per_sector(st)
        dt = time.time() - t0
        print(f"    n_max={nm}: dim_H={st.dim_H}, N_sectors={st.n_sectors},"
              f" 3*N={n_primitive_generators(nm)} generators ({dt:.2f}s)")

    # -----------------------------------------------------------------
    # [2] Hopf-hom cofiltered axiom at all triples.
    # -----------------------------------------------------------------
    print("\n[2] Hopf-hom cofiltered axiom Phi_{m,k} = Phi_{n,k} * Phi_{m,n}:")
    hopf_triple_results: List[Dict[str, Any]] = []
    hopf_all_bit_exact = True
    for (m, n, k) in triples:
        v = verify_hopf_cofiltered_axiom(m, n, k)
        hopf_triple_results.append({
            "m": v["m"], "n": v["n"], "k": v["k"],
            "dim_high": v["dim_high"], "dim_low": v["dim_low"],
            "bit_exact": v["bit_exact"],
            "residual_norm_squared": int(v["residual_norm_squared"]),
        })
        hopf_all_bit_exact = hopf_all_bit_exact and v["bit_exact"]
        status = "OK" if v["bit_exact"] else "FAIL"
        print(f"    Phi_({m},{n},{k}): dim {v['dim_low']}x{v['dim_high']} -> {status}")
    print(f"  Total triples: {len(triples)}; all bit-exact: {hopf_all_bit_exact}")

    # -----------------------------------------------------------------
    # [3] G_a generator compatibility at all pairs.
    # -----------------------------------------------------------------
    print("\n[3] G_a generator compatibility for each pair (m, k):")
    ga_pair_results: List[Dict[str, Any]] = []
    ga_all_bit_exact = True
    n_total_generators = 0
    for (m, k) in pairs:
        v = verify_Ga_generator_compatibility(m, k)
        n_total_generators += v["total_generators"]
        ga_pair_results.append({
            "m": v["m"], "k": v["k"],
            "total_generators": v["total_generators"],
            "n_survived": v["n_survived"],
            "n_killed": v["n_killed"],
            "expected_survived": v["expected_survived"],
            "expected_killed": v["expected_killed"],
            "all_bit_exact": v["all_bit_exact"],
        })
        ga_all_bit_exact = ga_all_bit_exact and v["all_bit_exact"]
        status = "OK" if v["all_bit_exact"] else "FAIL"
        print(f"    Pair ({m},{k}): {v['total_generators']} generators "
              f"({v['n_survived']} survive + {v['n_killed']} killed) -> {status}")
    print(f"  Total generator checks: {n_total_generators}; "
          f"all bit-exact: {ga_all_bit_exact}")

    # -----------------------------------------------------------------
    # [4] Class-level U*-action compatibility on chi and eta.
    # -----------------------------------------------------------------
    print("\n[4] Class-level U*-action compatibility on chi, eta:")
    class_pair_results: List[Dict[str, Any]] = []
    chi_all_bit_exact = True
    eta_all_bit_exact = True
    for (m, k) in pairs:
        sec_high = sectors_at_cutoff(m)
        sec_low = sectors_at_cutoff(k)
        psi_chi_high = class_vector(chi_classes[m], sec_high)
        psi_chi_low = class_vector(chi_classes[k], sec_low)
        psi_eta_high = class_vector(eta_classes[m], sec_high)
        psi_eta_low = class_vector(eta_classes[k], sec_low)
        v_chi = verify_class_action_compatibility(m, k, psi_chi_high, psi_chi_low)
        v_eta = verify_class_action_compatibility(m, k, psi_eta_high, psi_eta_low)
        class_pair_results.append({
            "m": m, "k": k,
            "chi_lhs_eq_rhs": v_chi["lhs_eq_rhs"],
            "chi_lhs_eq_psi_low": v_chi["lhs_eq_psi_low"],
            "chi_all_bit_exact": v_chi["all_bit_exact"],
            "eta_lhs_eq_rhs": v_eta["lhs_eq_rhs"],
            "eta_lhs_eq_psi_low": v_eta["lhs_eq_psi_low"],
            "eta_all_bit_exact": v_eta["all_bit_exact"],
        })
        chi_all_bit_exact = chi_all_bit_exact and v_chi["all_bit_exact"]
        eta_all_bit_exact = eta_all_bit_exact and v_eta["all_bit_exact"]
        chi_status = "OK" if v_chi["all_bit_exact"] else "FAIL"
        eta_status = "OK" if v_eta["all_bit_exact"] else "FAIL"
        print(f"    ({m},{k}): chi {chi_status}, eta {eta_status}")
    print(f"  Total pairs: {len(pairs)} per character; "
          f"chi: {chi_all_bit_exact}, eta: {eta_all_bit_exact}")

    # -----------------------------------------------------------------
    # [5] SL_2 commutativity --- categorical statement.
    # -----------------------------------------------------------------
    print("\n[5] SL_2 factor commutativity with P_{m,k}:")
    print("    SL_2 acts on the Peter--Weyl decoration (j_max axis,")
    print("    independent of n_max).  Acting first and then truncating")
    print("    or truncating first and then acting both yield the same")
    print("    result by the independence of axes.  Categorical bit-exact.")
    sl2_categorical = True

    # -----------------------------------------------------------------
    # [6] Summary
    # -----------------------------------------------------------------
    n_hopf_cofiltered = len(hopf_triple_results)
    n_ga_generators = n_total_generators
    n_class_chi = len(class_pair_results) * 2  # two boolean conditions per pair
    n_class_eta = len(class_pair_results) * 2
    n_total = n_hopf_cofiltered + n_ga_generators + n_class_chi + n_class_eta
    n_zero_residuals = (
        sum(1 for r in hopf_triple_results if r["bit_exact"])
        + sum(r["total_generators"] for r in ga_pair_results if r["all_bit_exact"])
        + sum(2 for r in class_pair_results if r["chi_all_bit_exact"])
        + sum(2 for r in class_pair_results if r["eta_all_bit_exact"])
    )

    print("\n[6] Bit-exact zero-residual summary:")
    print(f"    Hopf-hom cofiltered axiom triples : {n_hopf_cofiltered}")
    print(f"    G_a generator checks              : {n_ga_generators}")
    print(f"    Class-level chi compatibility     : {n_class_chi}")
    print(f"    Class-level eta compatibility     : {n_class_eta}")
    print(f"    Total bit-exact zero residuals    : {n_zero_residuals} / {n_total}")
    print()

    elapsed_total = time.time() - t_global
    print(f"Total wall time: {elapsed_total:.2f}s")

    # -----------------------------------------------------------------
    # Persist
    # -----------------------------------------------------------------
    data: Dict[str, Any] = {
        "sprint": "Q5'-ProSystem-Lockdown PS-2",
        "date_iso": "2026-06-06",
        "cutoffs": cutoffs,
        "n_triples": len(triples),
        "n_pairs": len(pairs),
        "mellin_slots": list(MELLIN_SLOTS),
        "N_sectors": {str(nm): N_sectors(nm) for nm in cutoffs},
        "n_primitive_generators": {str(nm): n_primitive_generators(nm) for nm in cutoffs},
        "hopf_cofiltered_results": hopf_triple_results,
        "hopf_cofiltered_all_bit_exact": hopf_all_bit_exact,
        "Ga_pair_results": ga_pair_results,
        "Ga_all_bit_exact": ga_all_bit_exact,
        "Ga_total_generator_checks": n_total_generators,
        "class_pair_results": class_pair_results,
        "class_chi_all_bit_exact": chi_all_bit_exact,
        "class_eta_all_bit_exact": eta_all_bit_exact,
        "SL2_categorical_commutativity": sl2_categorical,
        "bit_exact_summary": {
            "hopf_cofiltered_triples": n_hopf_cofiltered,
            "Ga_generator_checks": n_ga_generators,
            "class_chi_identities": n_class_chi,
            "class_eta_identities": n_class_eta,
            "total_zero_residuals": n_zero_residuals,
            "total_identities": n_total,
        },
        "wall_time_seconds": elapsed_total,
    }

    out_path = Path("debug/data/sprint_q5p_ps2_ustar_compatibility.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\nData written to: {out_path}")


if __name__ == "__main__":
    main()
