r"""
Sprint Q5'-ProSystem-Lockdown PS-1 — transition maps and cofiltered
axiom verification on the truncated Camporesi--Higuchi spectral triple
pro-system, extended to $n_{\max} = 5$.

Goal
----
Extend v3.60.0's pro-system functoriality result with three new ingredients:
1. The transition map $P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k$ as a
   closed-form algebra homomorphism (matrix object), not just a
   per-cell verification.
2. The cofiltered axiom $P_{m, k} = P_{n, k} \cdot P_{m, n}$ verified
   bit-exact for all triples $1 \le k < n < m$, $m \le 5$.
3. Extension of the per-sector closed-form panel and pull-back identity
   to $n_{\max} = 5$ (new cell, 6 new sectors, 20 total at level 5).

Together these promote the v3.60.0 pull-back compatibility result (verified
at consecutive cutoff pairs at the value level) to a closed-form inverse
system in the strict sense:\ the transitions are explicit
$N(k) \times N(m)$ 0/1 integer matrices, they commute as composed
homomorphisms by the cofiltered axiom, and the cocycle class data
pulls back through any chain bit-exactly.

This is PS-1 of the Pro-System-Lockdown sprint (sprint plan:\ four
sub-tracks PS-1/PS-2/PS-3/PS-4 toward locking down the pro-system in
preparation for Tannakian closure, multi-year).

Setup
-----
Cutoffs $n_{\max} \in \{1, 2, 3, 4, 5\}$. Sectors at cutoff $n_{\max}$
are $(n, l)$ with $1 \le n \le n_{\max}$, $0 \le l \le n$. Counts:
$N(1) = 2, N(2) = 5, N(3) = 9, N(4) = 14, N(5) = 20$.

Per-sector closed forms (v3.60.0):
- $\chi_{(n, l)} = +2$ if $l < n$, $-2n$ if $l = n$ (JLO HP$^{\mathrm{even}}$).
- $\eta_{(n, l)} = (2l + 1)(2n + 1)$ if $l < n$, $n(2n + 1)$ if $l = n$.

Pull-back $P^*_{m, k}[\psi^{(m)}] = [\psi^{(k)}]$ should hold bit-exact
at ALL pairs $1 \le k < m \le 5$, not just consecutive pairs, by
sector-locality.

Bit-exact panel (all required to be zero residual)
---------------------------------------------------
- Per-sector closed-form verifications:\ $5 \times 2 = 10$ verifications
  (one per cutoff per character), each covering all sectors at that
  cutoff. Total individual sector predictions:\ $\sum N(n_{\max}) =
  2 + 5 + 9 + 14 + 20 = 50$ per character, $100$ total.
- Cofiltered axiom triples $(m, n, k)$ with $1 \le k < n < m \le 5$:\
  $\binom{5}{3} = 10$ matrix-level identities.
- All-pairs pull-back $P^*_{m, k}[\psi^{(m)}] = [\psi^{(k)}]$ for
  $1 \le k < m \le 5$:\ $\binom{5}{2} = 10$ class-vector identities per
  character, $20$ total.

Total bit-exact zero residuals expected: $100 + 10 + 20 = 130$.

Output
------
- ``debug/data/sprint_q5p_ps1_transitions.json`` — bit-exact data.
- Console summary at end.

Discipline
----------
Bit-exact ``sympy.Rational`` / ``sympy.Integer`` throughout. No floats.
No PSLQ.

References
----------
- v3.60.0 Sprint Q5'-Stage1-Prosystem memo
  (``debug/sprint_q5p_prosystem_memo.md``).
- v3.66.0 Sprint Q5'-HP-Round3-FollowOns memo
  (``debug/sprint_q5p_hp_round3_followons_2026_06_06_memo.md``).
- ``geovac/pro_system.py`` (TransitionMap, cofiltered-axiom verifier).
- Paper 55 \S\ref{subsec:open_m2_m3}.
"""

from __future__ import annotations

import json
import time
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Tuple

from sympy import Integer, Matrix, Rational, zeros as sp_zeros

from geovac.spectral_triple import FockSpectralTriple
from geovac.pro_system import (
    TransitionMap,
    N_sectors,
    sectors_at_cutoff,
    verify_cofiltered_axiom,
)


# =====================================================================
# Helpers --- per-sector class extraction
# =====================================================================


def _sector_idempotent(st: FockSpectralTriple, sector_idx: int) -> Matrix:
    """Return e_s as dim_H x dim_H diagonal idempotent."""
    N = st.dim_H
    M = sp_zeros(N, N)
    for i in range(N):
        if st._state_to_sector[i] == sector_idx:
            M[i, i] = Integer(1)
    return M


def jlo_class_per_sector(
    st: FockSpectralTriple,
) -> Dict[Tuple[int, int], int]:
    """JLO HP^even leading vector chi_s = Tr(gamma e_s)."""
    gamma = st.grading
    out: Dict[Tuple[int, int], int] = {}
    for s_idx, (n, l) in enumerate(st.sectors):
        e_s = _sector_idempotent(st, s_idx)
        out[(n, l)] = int((gamma * e_s).trace())
    return out


def cm_eta_class_per_sector(
    st: FockSpectralTriple,
) -> Dict[Tuple[int, int], int]:
    """CM-eta residue class eta_s = Tr(gamma D e_s)."""
    gamma = st.grading
    D = st.dirac_operator
    out: Dict[Tuple[int, int], int] = {}
    for s_idx, (n, l) in enumerate(st.sectors):
        e_s = _sector_idempotent(st, s_idx)
        out[(n, l)] = int((gamma * D * e_s).trace())
    return out


def jlo_closed_form(n: int, l: int) -> int:
    """chi_{(n, l)} per-sector closed form."""
    return -2 * n if l == n else 2


def cm_eta_closed_form(n: int, l: int) -> int:
    """eta_{(n, l)} per-sector closed form."""
    return n * (2 * n + 1) if l == n else (2 * l + 1) * (2 * n + 1)


# =====================================================================
# Panel construction
# =====================================================================


def class_vector(
    cls: Dict[Tuple[int, int], int],
    sectors: List[Tuple[int, int]],
) -> List[int]:
    """Order class dict by canonical lex sector list."""
    return [cls[s] for s in sectors]


def main() -> None:
    print("=" * 72)
    print("Sprint Q5'-ProSystem-Lockdown  PS-1")
    print("Transition maps, cofiltered axiom, n_max = 5 extension")
    print("=" * 72)

    t_global = time.time()

    cutoffs = [1, 2, 3, 4, 5]
    triples: Dict[int, FockSpectralTriple] = {}
    n_actual: Dict[int, int] = {}
    n_closed: Dict[int, int] = {}
    dim_H: Dict[int, int] = {}
    jlo_classes: Dict[int, Dict[Tuple[int, int], int]] = {}
    eta_classes: Dict[int, Dict[Tuple[int, int], int]] = {}

    # -----------------------------------------------------------------
    # [1] Construct spectral triples at each cutoff, extract chi & eta.
    # -----------------------------------------------------------------
    print("\n[1] Constructing FockSpectralTriple at each cutoff:")
    for nm in cutoffs:
        t0 = time.time()
        st = FockSpectralTriple(n_max=nm)
        triples[nm] = st
        n_actual[nm] = st.n_sectors
        n_closed[nm] = N_sectors(nm)
        dim_H[nm] = st.dim_H
        jlo_classes[nm] = jlo_class_per_sector(st)
        eta_classes[nm] = cm_eta_class_per_sector(st)
        dt = time.time() - t0
        print(
            f"    n_max={nm}: dim_H={st.dim_H}, n_sectors={st.n_sectors}"
            f" (closed form {n_closed[nm]}), {dt:.2f}s"
        )

    n_sectors_match = all(n_actual[nm] == n_closed[nm] for nm in cutoffs)
    print(f"\n  Sector count closed form match (all cutoffs): {n_sectors_match}")

    # -----------------------------------------------------------------
    # [2] Per-sector closed-form verification at each cutoff.
    # -----------------------------------------------------------------
    print("\n[2] Per-sector closed-form verification:")
    closed_form_results: Dict[str, Any] = {"jlo": {}, "eta": {}}
    total_jlo_predictions = 0
    total_eta_predictions = 0
    jlo_all_match = True
    eta_all_match = True
    for nm in cutoffs:
        secs = sectors_at_cutoff(nm)
        jlo_mismatches = []
        eta_mismatches = []
        for (n, l) in secs:
            cf_jlo = jlo_closed_form(n, l)
            cf_eta = cm_eta_closed_form(n, l)
            obs_jlo = jlo_classes[nm][(n, l)]
            obs_eta = eta_classes[nm][(n, l)]
            total_jlo_predictions += 1
            total_eta_predictions += 1
            if obs_jlo != cf_jlo:
                jlo_mismatches.append({"sector": [n, l], "obs": obs_jlo, "cf": cf_jlo})
            if obs_eta != cf_eta:
                eta_mismatches.append({"sector": [n, l], "obs": obs_eta, "cf": cf_eta})
        closed_form_results["jlo"][nm] = {
            "sectors_tested": len(secs),
            "mismatches": jlo_mismatches,
            "all_match": len(jlo_mismatches) == 0,
        }
        closed_form_results["eta"][nm] = {
            "sectors_tested": len(secs),
            "mismatches": eta_mismatches,
            "all_match": len(eta_mismatches) == 0,
        }
        jlo_all_match = jlo_all_match and len(jlo_mismatches) == 0
        eta_all_match = eta_all_match and len(eta_mismatches) == 0
        print(
            f"    n_max={nm}: chi all-match={closed_form_results['jlo'][nm]['all_match']}, "
            f"eta all-match={closed_form_results['eta'][nm]['all_match']}"
        )
    print(
        f"  Total bit-exact predictions: chi={total_jlo_predictions}, "
        f"eta={total_eta_predictions} (all match: chi={jlo_all_match}, eta={eta_all_match})"
    )

    # -----------------------------------------------------------------
    # [3] Cofiltered axiom verification for all triples (m, n, k).
    # -----------------------------------------------------------------
    print("\n[3] Cofiltered axiom verification P_{m,k} = P_{n,k} * P_{m,n}:")
    triple_results: List[Dict[str, Any]] = []
    cofiltered_all = True
    for k, n, m in [(k, n, m)
                     for m in cutoffs
                     for n in cutoffs if 1 <= n < m
                     for k in cutoffs if 1 <= k < n]:
        if not (1 <= k < n < m <= 5):
            continue
        v = verify_cofiltered_axiom(m, n, k)
        triple_results.append({
            "m": v["m"],
            "n": v["n"],
            "k": v["k"],
            "dim_high": v["dim_high"],
            "dim_low": v["dim_low"],
            "bit_exact": v["bit_exact"],
            "residual_norm_squared": int(v["residual_norm_squared"]),
        })
        cofiltered_all = cofiltered_all and v["bit_exact"]
        status = "OK" if v["bit_exact"] else "FAIL"
        print(
            f"    ({m},{n},{k}): {P_dim_str(v)} -> {status} "
            f"(residual^2 = {int(v['residual_norm_squared'])})"
        )
    print(f"  Total triples tested: {len(triple_results)}; all bit-exact: {cofiltered_all}")

    # -----------------------------------------------------------------
    # [4] All-pairs pull-back identity on chi and eta classes.
    # -----------------------------------------------------------------
    print("\n[4] All-pairs pull-back P^*_{m,k}[psi^(m)] = [psi^(k)]:")
    pullback_results: List[Dict[str, Any]] = []
    pullback_all_chi = True
    pullback_all_eta = True
    for k, m in [(k, m) for m in cutoffs for k in cutoffs if 1 <= k < m]:
        P = TransitionMap(m, k)
        secs_low = sectors_at_cutoff(k)
        secs_high = sectors_at_cutoff(m)

        v_chi_high = class_vector(jlo_classes[m], secs_high)
        v_chi_low = class_vector(jlo_classes[k], secs_low)
        v_chi_pulled = P.apply_to_vector(v_chi_high)
        chi_match = v_chi_low == v_chi_pulled

        v_eta_high = class_vector(eta_classes[m], secs_high)
        v_eta_low = class_vector(eta_classes[k], secs_low)
        v_eta_pulled = P.apply_to_vector(v_eta_high)
        eta_match = v_eta_low == v_eta_pulled

        pullback_results.append({
            "m": m,
            "k": k,
            "chi_bit_exact": chi_match,
            "eta_bit_exact": eta_match,
        })
        pullback_all_chi = pullback_all_chi and chi_match
        pullback_all_eta = pullback_all_eta and eta_match
        print(f"    P^*_{{{m},{k}}}: chi {'OK' if chi_match else 'FAIL'}, "
              f"eta {'OK' if eta_match else 'FAIL'}")
    print(
        f"  Total pairs tested: {len(pullback_results)} per character; "
        f"chi all-OK: {pullback_all_chi}, eta all-OK: {pullback_all_eta}"
    )

    # -----------------------------------------------------------------
    # [5] Closed-form transition matrices --- explicit at n_max=5.
    # -----------------------------------------------------------------
    print("\n[5] Closed-form transition matrices at n_max = 5:")
    matrix_examples: Dict[str, Dict[str, Any]] = {}
    for (m, k) in [(5, 4), (5, 3), (5, 1)]:
        P = TransitionMap(m, k)
        M = P.matrix
        matrix_examples[f"P_{m}_{k}"] = {
            "n_high": m,
            "n_low": k,
            "shape": [M.rows, M.cols],
            "nnz": sum(1 for i in range(M.rows) for j in range(M.cols) if M[i, j] != 0),
            "trace_on_first_block": int(sum(M[i, i] for i in range(min(M.rows, M.cols)))),
        }
        print(
            f"    P_{{{m},{k}}}: shape {M.rows}x{M.cols}, "
            f"nnz = {matrix_examples[f'P_{m}_{k}']['nnz']}, "
            f"first-block trace = {matrix_examples[f'P_{m}_{k}']['trace_on_first_block']}"
        )

    # -----------------------------------------------------------------
    # [6] Bit-exact residual count summary.
    # -----------------------------------------------------------------
    n_per_sector = total_jlo_predictions + total_eta_predictions
    n_cofiltered = len(triple_results)
    n_pullback = 2 * len(pullback_results)
    n_total = n_per_sector + n_cofiltered + n_pullback
    n_zero_residuals = (
        n_per_sector * int(jlo_all_match and eta_all_match)
        + sum(1 for r in triple_results if r["bit_exact"])
        + sum(1 for r in pullback_results if r["chi_bit_exact"])
        + sum(1 for r in pullback_results if r["eta_bit_exact"])
    )
    # For closed-form: count each individual sector prediction, not the
    # per-cutoff "all_match" boolean.
    n_zero_residuals_correct = (
        (total_jlo_predictions if jlo_all_match else
            total_jlo_predictions - sum(len(r["mismatches"])
                                         for r in closed_form_results["jlo"].values()))
        + (total_eta_predictions if eta_all_match else
            total_eta_predictions - sum(len(r["mismatches"])
                                         for r in closed_form_results["eta"].values()))
        + sum(1 for r in triple_results if r["bit_exact"])
        + sum(1 for r in pullback_results if r["chi_bit_exact"])
        + sum(1 for r in pullback_results if r["eta_bit_exact"])
    )

    print("\n[6] Bit-exact zero-residual summary:")
    print(f"    Per-sector closed-form predictions: {n_per_sector} "
          f"(chi {total_jlo_predictions} + eta {total_eta_predictions})")
    print(f"    Cofiltered axiom triples           : {n_cofiltered}")
    print(f"    All-pairs pull-back identities     : {n_pullback}")
    print(f"    Total bit-exact zero residuals     : {n_zero_residuals_correct} / {n_total}")
    print()

    elapsed_total = time.time() - t_global
    print(f"Total wall time: {elapsed_total:.2f}s")

    # -----------------------------------------------------------------
    # Persist
    # -----------------------------------------------------------------
    data: Dict[str, Any] = {
        "sprint": "Q5'-ProSystem-Lockdown PS-1",
        "date_iso": "2026-06-06",
        "cutoffs": cutoffs,
        "dim_H": {str(nm): dim_H[nm] for nm in cutoffs},
        "N_sectors_actual": {str(nm): n_actual[nm] for nm in cutoffs},
        "N_sectors_closed_form": {str(nm): n_closed[nm] for nm in cutoffs},
        "sectors": {str(nm): [[n, l] for (n, l) in sectors_at_cutoff(nm)]
                    for nm in cutoffs},
        "jlo_class_vectors": {
            str(nm): [jlo_classes[nm][s] for s in sectors_at_cutoff(nm)]
            for nm in cutoffs
        },
        "eta_class_vectors": {
            str(nm): [eta_classes[nm][s] for s in sectors_at_cutoff(nm)]
            for nm in cutoffs
        },
        "closed_form_results": {
            "jlo": {str(nm): closed_form_results["jlo"][nm] for nm in cutoffs},
            "eta": {str(nm): closed_form_results["eta"][nm] for nm in cutoffs},
            "total_chi_predictions": total_jlo_predictions,
            "total_eta_predictions": total_eta_predictions,
            "chi_all_match": jlo_all_match,
            "eta_all_match": eta_all_match,
        },
        "cofiltered_axiom_results": triple_results,
        "cofiltered_all_bit_exact": cofiltered_all,
        "pullback_results": pullback_results,
        "pullback_all_chi_bit_exact": pullback_all_chi,
        "pullback_all_eta_bit_exact": pullback_all_eta,
        "transition_matrix_examples": matrix_examples,
        "bit_exact_summary": {
            "per_sector_predictions_total": n_per_sector,
            "cofiltered_axiom_triples": n_cofiltered,
            "pullback_identities": n_pullback,
            "total_zero_residuals": n_zero_residuals_correct,
            "total_identities": n_total,
        },
        "wall_time_seconds": elapsed_total,
    }

    out_path = Path("debug/data/sprint_q5p_ps1_transitions.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\nData written to: {out_path}")


def P_dim_str(v: Dict[str, Any]) -> str:
    return f"dim {v['dim_low']}x{v['dim_high']}"


if __name__ == "__main__":
    main()
