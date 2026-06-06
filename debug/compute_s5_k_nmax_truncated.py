"""
Sprint Q5'-S5-falsifier driver — master Mellin engine sources on the
truncated Bargmann-Segal Hardy-on-S^5 spectral triple.

Goal
----
Falsifier test for the chirality-parity selection rule established for
the Camporesi-Higuchi S^3 triple in Sprint Q5'-CH-1 (2026-06-05):

    On CH-S^3:
        Tr(D^k e^{-tD^2})       nonzero for k in {0,2}, identically 0 for k=1
        Tr(gamma D^k e^{-tD^2}) identically 0 for k in {0,2}, nonzero for k=1

The structural argument in CH-1 attributed this to two facts:
    (i)  CH spectrum {+-(n+1/2)} is +/-symmetric with equal multiplicity
    (ii) chirality alignment gamma_ii = chi_i, Lambda_ii = chi_i |lambda_i|

The S^5 Bargmann-Segal Hardy sector breaks both: it has
    (i)  spectrum {N + 3/2} STRICTLY POSITIVE (holomorphic sector only)
    (ii) no chirality grading at all (no spinor lift, per Paper 24 §V.4)

This sprint computes:
    Tr(L^k e^{-tL^2})         on the diagonal Euler operator L = N_hat + 3/2
    Tr(D^k e^{-tD^2})         on D = L + kappa A   (graph-regulated)
    Tr(chi L^k e^{-tL^2})     for trial Z/2 gradings chi (no canonical gamma)

The two trial gradings tested:
    chi_N(v) := (-1)^N_v       shell-parity
    chi_l(v) := (-1)^l_v       angular-parity

Both are natural Z/2 gradings on the Bargmann graph; neither has a
half-integer-weight content that would be the genuine Hardy-sector
analog of CH chirality (which is forbidden by Paper 24 §V.4).

Verdict gates
-------------
- BREAKS-AS-EXPECTED: Tr(L^k ...) nonzero for ALL k (no spectral
  cancellation), AND no trial chi produces the CH-1 parity rule. The
  rule was CH-specific.
- BREAKS-DIFFERENTLY: parity rule fails via some other mechanism.
- HOLDS-UNEXPECTEDLY: a natural trial chi reproduces the CH-1 rule.

Methodology
-----------
All matrices stay sympy.Rational. No floating point. Diagonal-fast-path
for L (the Euler operator IS diagonal). Full-D path uses the Bargmann
adjacency A from geovac.nuclear.bargmann_graph, with kappa = -1/16
matching the Paper 0 topological constant; we sample j_max=10 on the
diagonal and j_max=6 on the full D.
"""

from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Integer, Matrix, Rational, eye as sp_eye

from geovac.nuclear.bargmann_graph import build_bargmann_graph, BargmannGraph


KAPPA = Rational(-1, 16)  # Paper 0 universal topological constant


def fraction_to_rational(f: Fraction) -> Rational:
    """Convert fractions.Fraction to sympy.Rational without going through float."""
    return Rational(f.numerator, f.denominator)


def build_diag_Lambda(g: BargmannGraph) -> Matrix:
    """Bargmann diagonal Euler-shifted operator L = N_hat + 3/2 as sympy Matrix."""
    N = g.n_nodes
    L = sp.zeros(N, N)
    for i, d in enumerate(g.diagonal):
        L[i, i] = fraction_to_rational(d)
    return L


def build_adjacency(g: BargmannGraph) -> Matrix:
    """Bargmann adjacency (squared dipole matrix elements) as sympy Matrix.

    Note: the entries of `g.adjacency` are SQUARED matrix elements, hence rational.
    The actual dipole hopping matrix element is the (real) square root, which we
    use here as the off-diagonal piece of the "full D" perturbation. For the
    Mellin-engine moment computation we need entries of A^j to be exact rationals;
    the squared entries serve this purpose directly and are the natural choice
    for a discrete graph regulator. The structural verdict is INDEPENDENT of
    whether we use the squared or the square-root convention because the
    spectral-parity argument is about Tr(L^j) and L is diagonal — A enters only
    as a higher-order perturbation.
    """
    N = g.n_nodes
    A = sp.zeros(N, N)
    for (i, j), w in g.adjacency.items():
        ratw = fraction_to_rational(w)
        A[i, j] = ratw
        A[j, i] = ratw
    return A


def build_full_D(g: BargmannGraph) -> Matrix:
    """Full Dirac-graph regulator D = L + kappa A on the Bargmann graph.

    L is the diagonal Euler operator (matching the CH-1 driver's "Lambda").
    kappa A is the graph-regulator perturbation; kappa = -1/16 from Paper 0.
    """
    return build_diag_Lambda(g) + KAPPA * build_adjacency(g)


def _is_diagonal(M: Matrix) -> bool:
    N = M.shape[0]
    for i in range(N):
        for j in range(N):
            if i != j and M[i, j] != 0:
                return False
    return True


def trace_powers(M: Matrix, j_max: int) -> List[sp.Expr]:
    """Compute Tr(M^j) for j = 0..j_max iteratively. Fast path for diagonal."""
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


def chi_supertrace_powers_diagonal(L: Matrix, chi_diag: List[Rational],
                                   j_max: int) -> List[sp.Expr]:
    """Compute Tr(chi . L^j) when L is diagonal and chi is given as a list."""
    N = L.shape[0]
    diag_L = [L[i, i] for i in range(N)]
    return [
        sp.simplify(sum(c * d ** j for c, d in zip(chi_diag, diag_L)))
        for j in range(j_max + 1)
    ]


def chi_supertrace_powers_full(D: Matrix, chi_diag: List[Rational],
                               j_max: int) -> List[sp.Expr]:
    """Compute Tr(chi . D^j) for general D. chi is a diagonal grading.
    Returns Tr(diag(chi) . D^j) iteratively."""
    N = D.shape[0]
    chi = sp.zeros(N, N)
    for i, c in enumerate(chi_diag):
        chi[i, i] = c
    powers: List[sp.Expr] = []
    M_j = sp_eye(N)
    for j in range(j_max + 1):
        powers.append(sp.simplify((chi * M_j).trace()))
        if j < j_max:
            M_j = M_j * D
    return powers


def assemble_Gk_smallt(moments: List[sp.Expr], k: int, j_max: int) -> List[sp.Expr]:
    """Small-t coefficients of G_k(t) = Tr(M^k e^{-tM^2})
                              = sum_{j>=0} (-t)^j/j! Tr(M^{k+2j})."""
    coeffs: List[sp.Expr] = []
    j_eff_max = (j_max - k) // 2
    for j in range(j_eff_max + 1):
        idx = k + 2 * j
        sign = (Rational(-1) ** j) / sp.factorial(j)
        coeffs.append(sp.simplify(sign * moments[idx]))
    return coeffs


def build_chi_N(g: BargmannGraph) -> List[Rational]:
    """Trial grading chi_N(v) = (-1)^N_v (shell-parity)."""
    return [Rational((-1) ** N) for (N, l, m) in g.nodes]


def build_chi_l(g: BargmannGraph) -> List[Rational]:
    """Trial grading chi_l(v) = (-1)^l_v (angular-parity).
    Note: on shell N, l has same parity as N, so on every shell
    chi_l is constant — this collapses to chi_N up to a global sign."""
    return [Rational((-1) ** l) for (N, l, m) in g.nodes]


def trial_chi_balance(chi_diag: List[Rational]) -> Tuple[int, int]:
    """Count +/- entries of a trial grading."""
    pos = sum(1 for c in chi_diag if c == 1)
    neg = sum(1 for c in chi_diag if c == -1)
    return pos, neg


def compute_for_nmax(N_max: int, j_max: int) -> Dict:
    """Compute the S^5 master Mellin engine data at a given N_max."""
    t0 = time.time()
    g = build_bargmann_graph(N_max)
    L = build_diag_Lambda(g)
    A = build_adjacency(g)
    D_full = L + KAPPA * A
    N = g.n_nodes

    print(f"  N_max={N_max}: dim_H={N}, n_edges={len(g.adjacency)}, kappa={KAPPA}")

    spectrum_L = []
    for i, (Nv, lv, mv) in enumerate(g.nodes):
        spectrum_L.append({
            "node": [Nv, lv, mv],
            "lambda": str(L[i, i]),
        })

    chi_N_diag = build_chi_N(g)
    chi_l_diag = build_chi_l(g)
    posN, negN = trial_chi_balance(chi_N_diag)
    posl, negl = trial_chi_balance(chi_l_diag)
    print(f"    chi_N (shell-parity) balance: +{posN} / -{negN}")
    print(f"    chi_l (angular-parity) balance: +{posl} / -{negl}")

    # Moments on diagonal Euler operator L
    print(f"    computing Tr(L^j) for j=0..{j_max}...")
    tr_L = trace_powers(L, j_max)

    # Trial supertraces on L
    print(f"    computing Tr(chi_N . L^j)...")
    sup_chiN_L = chi_supertrace_powers_diagonal(L, chi_N_diag, j_max)
    print(f"    computing Tr(chi_l . L^j)...")
    sup_chil_L = chi_supertrace_powers_diagonal(L, chi_l_diag, j_max)

    # Full D moments (slower; cap j_max if needed)
    j_max_full = min(j_max, 6)
    print(f"    computing Tr(D^j) (full) for j=0..{j_max_full}...")
    tr_D_full = trace_powers(D_full, j_max_full)
    print(f"    computing Tr(chi_N . D^j) (full)...")
    sup_chiN_D = chi_supertrace_powers_full(D_full, chi_N_diag, j_max_full)
    print(f"    computing Tr(chi_l . D^j) (full)...")
    sup_chil_D = chi_supertrace_powers_full(D_full, chi_l_diag, j_max_full)

    # Assemble G_k, S_k small-t expansions
    G_L = {k: assemble_Gk_smallt(tr_L, k, j_max) for k in (0, 1, 2)}
    S_chiN_L = {k: assemble_Gk_smallt(sup_chiN_L, k, j_max) for k in (0, 1, 2)}
    S_chil_L = {k: assemble_Gk_smallt(sup_chil_L, k, j_max) for k in (0, 1, 2)}
    G_D_full = {k: assemble_Gk_smallt(tr_D_full, k, j_max_full) for k in (0, 1, 2)}
    S_chiN_D_full = {k: assemble_Gk_smallt(sup_chiN_D, k, j_max_full) for k in (0, 1, 2)}
    S_chil_D_full = {k: assemble_Gk_smallt(sup_chil_D, k, j_max_full) for k in (0, 1, 2)}

    elapsed = time.time() - t0
    print(f"    done in {elapsed:.1f}s")

    return {
        "N_max": N_max,
        "dim_H": N,
        "kappa": str(KAPPA),
        "n_edges": len(g.adjacency),
        "spectrum_L": spectrum_L,
        "chi_N_balance": [posN, negN],
        "chi_l_balance": [posl, negl],
        "moments_L": {
            "tr_L_j": [str(x) for x in tr_L],
            "sup_chiN_L_j": [str(x) for x in sup_chiN_L],
            "sup_chil_L_j": [str(x) for x in sup_chil_L],
        },
        "moments_D_full": {
            "tr_D_j": [str(x) for x in tr_D_full],
            "sup_chiN_D_j": [str(x) for x in sup_chiN_D],
            "sup_chil_D_j": [str(x) for x in sup_chil_D],
            "j_max_full": j_max_full,
        },
        "small_t_Gk_L": {
            str(k): [str(c) for c in G_L[k]] for k in (0, 1, 2)
        },
        "small_t_Sk_chiN_L": {
            str(k): [str(c) for c in S_chiN_L[k]] for k in (0, 1, 2)
        },
        "small_t_Sk_chil_L": {
            str(k): [str(c) for c in S_chil_L[k]] for k in (0, 1, 2)
        },
        "small_t_Gk_D_full": {
            str(k): [str(c) for c in G_D_full[k]] for k in (0, 1, 2)
        },
        "small_t_Sk_chiN_D_full": {
            str(k): [str(c) for c in S_chiN_D_full[k]] for k in (0, 1, 2)
        },
        "small_t_Sk_chil_D_full": {
            str(k): [str(c) for c in S_chil_D_full[k]] for k in (0, 1, 2)
        },
        "wall_seconds": elapsed,
    }


def check_parity_rule(data: Dict, n_max_key: str) -> Dict:
    """Check whether the CH-1 parity rule (Tr(L^j)=0 for j odd, Tr(chi L^j)=0
    for j even) holds at this N_max. Return verdict dict.

    CH-1 rule:
        Tr(Lambda^j) = 0 for j odd    (spectral +/-symmetry)
        Tr(gamma Lambda^j) = 0 for j even (chirality alignment)
    """
    d = data[n_max_key]
    tr_L = d["moments_L"]["tr_L_j"]
    sup_chiN_L = d["moments_L"]["sup_chiN_L_j"]
    sup_chil_L = d["moments_L"]["sup_chil_L_j"]

    # Tr(L^j) for j odd: should be 0 on CH, will NOT be 0 on Hardy (positive spec)
    spectral_parity_zeros = []
    for j in range(len(tr_L)):
        if j % 2 == 1:
            spectral_parity_zeros.append({
                "j": j,
                "Tr(L^j)": tr_L[j],
                "is_zero": tr_L[j] == "0",
            })

    # Tr(chi_N L^j) for j even: should be 0 on CH if chi_N aligns; check both
    chiN_even_zeros = []
    for j in range(len(sup_chiN_L)):
        if j % 2 == 0:
            chiN_even_zeros.append({
                "j": j,
                "Tr(chi_N L^j)": sup_chiN_L[j],
                "is_zero": sup_chiN_L[j] == "0",
            })

    chil_even_zeros = []
    for j in range(len(sup_chil_L)):
        if j % 2 == 0:
            chil_even_zeros.append({
                "j": j,
                "Tr(chi_l L^j)": sup_chil_L[j],
                "is_zero": sup_chil_L[j] == "0",
            })

    return {
        "spectral_parity_zeros (j odd)": spectral_parity_zeros,
        "chiN_chirality_zeros (j even)": chiN_even_zeros,
        "chil_chirality_zeros (j even)": chil_even_zeros,
    }


def summarize(output: Dict) -> None:
    """Print a headline summary."""
    print("\n" + "=" * 70)
    print("Sprint Q5'-S5-falsifier — Master Mellin Engine on Truncated S^5 Hardy")
    print("=" * 70)

    for N_max in (2, 3, 4):
        key = f"N_max={N_max}"
        if key not in output:
            continue
        d = output[key]
        print(f"\n--- N_max = {N_max}, dim_H = {d['dim_H']}, n_edges = {d['n_edges']} ---")
        print(f"  chi_N balance: +{d['chi_N_balance'][0]} / -{d['chi_N_balance'][1]}")
        print(f"  chi_l balance: +{d['chi_l_balance'][0]} / -{d['chi_l_balance'][1]}")

        print(f"\n  Bargmann L moments Tr(L^j), j=0..6:")
        for j in range(min(7, len(d['moments_L']['tr_L_j']))):
            print(f"    j={j}: {d['moments_L']['tr_L_j'][j]}")

        print(f"\n  Trial chi_N supertrace Tr(chi_N . L^j), j=0..6:")
        for j in range(min(7, len(d['moments_L']['sup_chiN_L_j']))):
            print(f"    j={j}: {d['moments_L']['sup_chiN_L_j'][j]}")

        print(f"\n  Trial chi_l supertrace Tr(chi_l . L^j), j=0..6:")
        for j in range(min(7, len(d['moments_L']['sup_chil_L_j']))):
            print(f"    j={j}: {d['moments_L']['sup_chil_L_j'][j]}")

        verdict = check_parity_rule(output, key)
        ch1_spec_holds = all(z["is_zero"] for z in verdict["spectral_parity_zeros (j odd)"])
        ch1_chiN_holds = all(z["is_zero"] for z in verdict["chiN_chirality_zeros (j even)"])
        ch1_chil_holds = all(z["is_zero"] for z in verdict["chil_chirality_zeros (j even)"])
        print(f"\n  CH-1 rule check at N_max={N_max}:")
        print(f"    Spectral +/- symmetry (Tr(L^j)=0 for j odd): {ch1_spec_holds}")
        print(f"    chi_N alignment (Tr(chi_N L^j)=0 for j even): {ch1_chiN_holds}")
        print(f"    chi_l alignment (Tr(chi_l L^j)=0 for j even): {ch1_chil_holds}")


def write_verdict(output: Dict) -> Dict:
    """Compute final verdict-against-gate."""
    verdicts = {}
    for N_max in (2, 3, 4):
        key = f"N_max={N_max}"
        if key not in output:
            continue
        v = check_parity_rule(output, key)
        ch1_spec = all(z["is_zero"] for z in v["spectral_parity_zeros (j odd)"])
        ch1_chiN = all(z["is_zero"] for z in v["chiN_chirality_zeros (j even)"])
        ch1_chil = all(z["is_zero"] for z in v["chil_chirality_zeros (j even)"])
        verdicts[key] = {
            "ch1_spectral_symmetry_holds": ch1_spec,
            "ch1_chiN_alignment_holds": ch1_chiN,
            "ch1_chil_alignment_holds": ch1_chil,
            "ch1_rule_holds": ch1_spec and (ch1_chiN or ch1_chil),
        }

    # Top-level verdict
    any_chi_holds = any(
        verdicts[k]["ch1_chiN_alignment_holds"] or verdicts[k]["ch1_chil_alignment_holds"]
        for k in verdicts
    )
    spec_holds = all(verdicts[k]["ch1_spectral_symmetry_holds"] for k in verdicts)

    if spec_holds and any_chi_holds:
        gate = "HOLDS-UNEXPECTEDLY"
    elif spec_holds and not any_chi_holds:
        gate = "BREAKS-DIFFERENTLY (spectral holds; chi does not)"
    elif not spec_holds and any_chi_holds:
        gate = "BREAKS-DIFFERENTLY (spectral fails; chi works)"
    else:
        gate = "BREAKS-AS-EXPECTED"

    return {"per_N_max": verdicts, "overall_gate": gate}


def main() -> None:
    j_max = 10
    output: Dict = {"sprint": "Q5'-S5-falsifier", "j_max": j_max}

    for N_max in (2, 3, 4):
        print(f"\n=== N_max = {N_max} ===")
        output[f"N_max={N_max}"] = compute_for_nmax(N_max, j_max)

    output["verdict"] = write_verdict(output)

    out_path = Path("debug/data/sprint_q5p_s5_falsifier_data.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
    print(f"\nOutput written: {out_path}")

    summarize(output)

    print("\n" + "=" * 70)
    print(f"OVERALL VERDICT: {output['verdict']['overall_gate']}")
    print("=" * 70)


if __name__ == "__main__":
    main()
