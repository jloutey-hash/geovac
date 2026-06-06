"""
Sprint Q5'-CH-1 driver — master Mellin engine sources on the truncated
Camporesi--Higuchi spectral triple.

Goal
----
Compute the master Mellin engine sources
    G_k(t) = Tr(D^k . exp(-t D^2))                        [trace]
    S_k(t) = Tr(gamma . D^k . exp(-t D^2))                [supertrace]
at k in {0, 1, 2} on the truncated Fock spectral triple at
n_max in {2, 3, 4}, then check whether the structure of the
small-t expansion aligns with the master Mellin partition M1/M2/M3
(Paper 32 §VIII, Theorem thm:pi_source_case_exhaustion and
Remark rem:master_mellin_domain).

Methodology
-----------
Rather than diagonalize D (slow, mixes chirality basis), reconstruct
G_k(t), S_k(t) from the moment sequences
    tr_j  = Tr(D^j)
    sup_j = Tr(gamma . D^j)
via the small-t Taylor expansion
    Tr(D^k . exp(-t D^2)) = sum_{j>=0} (-t)^j / j! . Tr(D^{k+2j}).

All matrices stay sympy.Rational. No floating point. No PSLQ at this
stage — at finite n_max, every coefficient is a rational and the M1/M2/M3
period content is structural, not numerical (periods enter only in the
Mellin transform and the n_max -> infinity limit).

Two D variants are computed for comparison:
    D_full   = Lambda + kappa . A   (full Dirac graph, kappa = -1/16)
    D_diag   = Lambda                (clean CH diagonal spectrum)
The D_diag variant is the "naive CH" Mellin engine, connecting directly
to the published M_k periods on the continuum sphere. The D_full variant
shows what the graph regulator contributes.

Output
------
JSON at debug/data/sprint_q5p_ch1_data.json containing for each n_max:
    spectrum, dim_H, moment sequences tr_j and sup_j,
    small-t coefficients of G_k(t), S_k(t) for k in {0, 1, 2}.

This is a scoping sprint: the verdict on the k-slot Tannakian-relevance
(per debug/sprint_q5p_k_slot_tannakian_memo.md) lives in the structural
pattern of the coefficients, not in any single numerical match.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List

import sympy as sp
from sympy import Integer, Matrix, Rational, eye as sp_eye, factorial

from geovac.spectral_triple import FockSpectralTriple


def _is_diagonal(M: Matrix) -> bool:
    """Cheap check for a diagonal sympy matrix."""
    N = M.shape[0]
    for i in range(N):
        for j in range(N):
            if i != j and M[i, j] != 0:
                return False
    return True


def trace_powers(M: Matrix, j_max: int) -> List[sp.Expr]:
    """Compute Tr(M^j) for j = 0, 1, ..., j_max iteratively.

    Fast path for diagonal matrices: Tr(M^j) = sum_i M_ii^j (no matrix mul).
    Returns list of exact sympy expressions (length j_max + 1).
    """
    N = M.shape[0]
    if _is_diagonal(M):
        diag = [M[i, i] for i in range(N)]
        return [sp.simplify(sum(d ** j for d in diag)) for j in range(j_max + 1)]
    powers: List[sp.Expr] = []
    M_j = sp_eye(N)
    for j in range(j_max + 1):
        powers.append(sp.simplify(M_j.trace()))
        if j < j_max:
            M_j = M_j * M
    return powers


def supertrace_powers(M: Matrix, gamma: Matrix, j_max: int) -> List[sp.Expr]:
    """Compute Tr(gamma . M^j) for j = 0, ..., j_max iteratively.

    Fast path when M is diagonal: (gamma . M^j)_ii = gamma_ii . M_ii^j.
    """
    N = M.shape[0]
    if _is_diagonal(M) and _is_diagonal(gamma):
        diag_M = [M[i, i] for i in range(N)]
        diag_g = [gamma[i, i] for i in range(N)]
        return [
            sp.simplify(sum(g * d ** j for g, d in zip(diag_g, diag_M)))
            for j in range(j_max + 1)
        ]
    powers: List[sp.Expr] = []
    M_j = sp_eye(N)
    for j in range(j_max + 1):
        powers.append(sp.simplify((gamma * M_j).trace()))
        if j < j_max:
            M_j = M_j * M
    return powers


def assemble_Gk_smallt(moments: List[sp.Expr], k: int, j_max: int) -> List[sp.Expr]:
    """Assemble the small-t Taylor coefficients of
        G_k(t) = Tr(D^k . exp(-t D^2)) = sum_{j>=0} (-t)^j / j! . Tr(D^{k+2j})

    Returns the list of coefficients [c_0, c_1, ...] where G_k(t) = sum c_j t^j,
    truncated at the largest j with k + 2j <= j_max.
    """
    coeffs: List[sp.Expr] = []
    j_eff_max = (j_max - k) // 2
    for j in range(j_eff_max + 1):
        idx = k + 2 * j
        sign = (Rational(-1) ** j) / sp.factorial(j)
        coeffs.append(sp.simplify(sign * moments[idx]))
    return coeffs


def compute_for_nmax(n_max: int, j_max: int) -> Dict:
    """Compute the master Mellin engine data at a given n_max."""
    t0 = time.time()
    st = FockSpectralTriple(n_max=n_max)
    D_full = st.dirac_operator
    Lambda = st.diagonal_part
    gamma = st.grading
    N = st.dim_H

    print(f"  n_max={n_max}: dim_H={N}, kappa={st._kappa}")

    # Spectrum on the diagonal Lambda (clean CH)
    spectrum_Lambda = []
    for i in range(N):
        spectrum_Lambda.append({
            "lambda": str(Lambda[i, i]),
            "chi": int(gamma[i, i]),
        })

    # Moment sequences
    print(f"    computing tr_Lambda_j for j=0..{j_max}...")
    tr_Lambda = trace_powers(Lambda, j_max)
    print(f"    computing sup_Lambda_j for j=0..{j_max}...")
    sup_Lambda = supertrace_powers(Lambda, gamma, j_max)

    # Full D moments (slower; cap j_max if needed)
    j_max_full = min(j_max, 6)  # full D matrix multiplication caps at j=6 for tractability
    print(f"    computing tr_D_j (full) for j=0..{j_max_full}...")
    tr_D_full = trace_powers(D_full, j_max_full)
    print(f"    computing sup_D_j (full) for j=0..{j_max_full}...")
    sup_D_full = supertrace_powers(D_full, gamma, j_max_full)

    # Assemble G_k(t) and S_k(t) small-t expansions
    G_Lambda = {k: assemble_Gk_smallt(tr_Lambda, k, j_max) for k in (0, 1, 2)}
    S_Lambda = {k: assemble_Gk_smallt(sup_Lambda, k, j_max) for k in (0, 1, 2)}
    G_full = {k: assemble_Gk_smallt(tr_D_full, k, j_max_full) for k in (0, 1, 2)}
    S_full = {k: assemble_Gk_smallt(sup_D_full, k, j_max_full) for k in (0, 1, 2)}

    elapsed = time.time() - t0
    print(f"    done in {elapsed:.1f}s")

    return {
        "n_max": n_max,
        "dim_H": N,
        "kappa": str(st._kappa),
        "spectrum_Lambda": spectrum_Lambda,
        "moments_Lambda": {
            "tr_Lambda_j": [str(x) for x in tr_Lambda],
            "sup_Lambda_j": [str(x) for x in sup_Lambda],
        },
        "moments_D_full": {
            "tr_D_j": [str(x) for x in tr_D_full],
            "sup_D_j": [str(x) for x in sup_D_full],
            "j_max_full": j_max_full,
        },
        "small_t_Gk_Lambda": {
            str(k): [str(c) for c in G_Lambda[k]] for k in (0, 1, 2)
        },
        "small_t_Sk_Lambda": {
            str(k): [str(c) for c in S_Lambda[k]] for k in (0, 1, 2)
        },
        "small_t_Gk_full": {
            str(k): [str(c) for c in G_full[k]] for k in (0, 1, 2)
        },
        "small_t_Sk_full": {
            str(k): [str(c) for c in S_full[k]] for k in (0, 1, 2)
        },
        "wall_seconds": elapsed,
    }


def summarize(output: Dict) -> None:
    """Print a headline summary."""
    print("\n" + "=" * 70)
    print("Sprint Q5'-CH-1 — Master Mellin Engine on Truncated CH Triple")
    print("=" * 70)

    for n_max in (2, 3, 4):
        key = f"n_max={n_max}"
        if key not in output:
            continue
        d = output[key]
        print(f"\n--- n_max = {n_max}, dim_H = {d['dim_H']} ---")

        print(f"\n  Chirality split of CH spectrum:")
        chi_pos = sum(1 for s in d["spectrum_Lambda"] if s["chi"] == 1)
        chi_neg = sum(1 for s in d["spectrum_Lambda"] if s["chi"] == -1)
        print(f"    chi=+1: {chi_pos}, chi=-1: {chi_neg}")

        print(f"\n  Lambda moments Tr(Lambda^j), j=0..6:")
        for j in range(7):
            print(f"    j={j}: {d['moments_Lambda']['tr_Lambda_j'][j]}")

        print(f"\n  Supertrace Lambda moments Tr(gamma . Lambda^j), j=0..6:")
        for j in range(7):
            print(f"    j={j}: {d['moments_Lambda']['sup_Lambda_j'][j]}")

        for k in (0, 1, 2):
            print(f"\n  G_{k}(t) on Lambda, small-t coefficients (k={k}):")
            for j, c in enumerate(d["small_t_Gk_Lambda"][str(k)][:5]):
                print(f"    t^{j}: {c}")
            print(f"  S_{k}(t) on Lambda, small-t coefficients (k={k}):")
            for j, c in enumerate(d["small_t_Sk_Lambda"][str(k)][:5]):
                print(f"    t^{j}: {c}")


def main() -> None:
    j_max = 10
    output: Dict = {"sprint": "Q5'-CH-1", "j_max": j_max}

    for n_max in (2, 3, 4):
        print(f"\n=== n_max = {n_max} ===")
        output[f"n_max={n_max}"] = compute_for_nmax(n_max, j_max)

    out_path = Path("debug/data/sprint_q5p_ch1_data.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print(f"\nOutput written: {out_path}")

    summarize(output)


if __name__ == "__main__":
    main()
