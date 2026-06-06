"""Part 2: explicit identification of non-primitive coproduct content and
comparison to v3.61.0 Track B closure drift.

Goal: given T1 · T2 = c · e_s in the enriched algebra, the multiplicativity
constraint on a putative primitive coproduct forces
  Delta(c · e_s) = c · (e_s ⊗ 1 + 1 ⊗ e_s)
to equal
  Delta(T1) · Delta(T2) = T1*T2 ⊗ 1 + T1 ⊗ T2 + T2 ⊗ T1 + 1 ⊗ T1*T2
                       = c · e_s ⊗ 1 + 1 ⊗ c · e_s + (T1 ⊗ T2 + T2 ⊗ T1).

The non-primitive content is exactly (T1 ⊗ T2 + T2 ⊗ T1) — its non-vanishing
forces non-primitive coproduct for the enriched substrate.

We compute this content bit-exactly and tag it against the v3.61.0 Track B
drift signature (kappa^2 = 1/256 timescale).
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Rational, Integer, Matrix, zeros as sp_zeros

from geovac.spectral_triple import FockSpectralTriple


def serialize_rational(x):
    if isinstance(x, sp.Rational):
        return f"{x.p}/{x.q}" if x.q != 1 else str(x.p)
    return str(x)


def matrix_nonzero(M: Matrix) -> bool:
    return any(M[i, j] != 0 for i in range(M.rows) for j in range(M.cols))


def main():
    print("=" * 70)
    print("Sprint Q5'-OffDiag-Dirac Part 2 — Non-primitive content extraction")
    print("=" * 70)

    triple = FockSpectralTriple(n_max=2)
    N = triple.dim_H

    # ------------------------------------------------------------------
    # Rebuild sector idempotents and transition generators
    # ------------------------------------------------------------------
    idem: Dict[Tuple[int, int], Matrix] = {}
    for sec in triple.sectors:
        e = sp_zeros(N, N)
        for i, k in enumerate(triple._state_to_sector):
            if triple.sectors[k] == sec:
                e[i, i] = Integer(1)
        idem[sec] = e

    D = triple.dirac_operator
    Lambda = triple.diagonal_part
    A_part = D - Lambda

    trans: Dict[Tuple[Tuple[int, int], Tuple[int, int]], Matrix] = {}
    for s_from in triple.sectors:
        for s_to in triple.sectors:
            if s_to == s_from:
                continue
            T = idem[s_to] * A_part * idem[s_from]
            if matrix_nonzero(T):
                trans[(s_from, s_to)] = T

    # ------------------------------------------------------------------
    # NON-PRIMITIVE CONTENT EXTRACTION
    # ------------------------------------------------------------------
    print("\n--- Non-primitive coproduct content ---")
    print("For each chain composition T2 * T1 = c * e_s, the non-primitive")
    print("content of putative primitive Delta is (T1 (x) T2 + T2 (x) T1).")
    print("We compute its eta-bivalued image: (eta x eta)(T1 (x) T2 + T2 (x) T1)")
    print("which is 2 * eta(T1) * eta(T2).")
    print()

    # Compute eta values for all transitions
    g = triple.grading
    eta_trans = {}
    for k, T in trans.items():
        eta_trans[k] = sp.trace(g * D * T)

    # Compute all chain compositions T2 · T1 in algebra and the
    # non-primitive content
    nonprim_panel: List[Dict] = []
    print(f"{'T1 -> T2':40} {'lands on':25} {'coeff':12} {'2*eta(T1)*eta(T2)':25}")
    print("-" * 110)
    for (s0, s1), T1 in trans.items():
        for (s2, s3), T2 in trans.items():
            if s1 != s2:
                continue  # nonchain
            P = T2 * T1  # T2 after T1
            if not matrix_nonzero(P):
                continue
            # Decompose: where does it land?
            # By matrix-algebra structure: ends in sector pair (s0, s3)
            # If s0 == s3: idempotent e_{s0}, with coefficient = P[i,i] for i in s0
            # If s0 != s3: transition T_{s0 -> s3}, with coefficient
            if s0 == s3:
                # Find sample diagonal entry
                idx = next(i for i in range(N) if idem[s0][i, i] != 0)
                coeff = P[idx, idx]
                land = f"e_{s0}"
            else:
                # Find sample nonzero entry pattern
                T_ref = trans.get((s0, s3))
                if T_ref is None:
                    # New transition not in our generator set (would be
                    # a 2-step path); flag
                    land = f"two-step T_{s0}->{s3} (NEW GENERATOR)"
                    coeff = None
                else:
                    # Find a reference entry of T_ref
                    found_ref = None
                    for i in range(N):
                        for j in range(N):
                            if T_ref[i, j] != 0:
                                found_ref = (i, j, T_ref[i, j])
                                break
                        if found_ref:
                            break
                    if found_ref:
                        ii, jj, refv = found_ref
                        coeff = sp.Rational(P[ii, jj]) / sp.Rational(refv)
                        land = f"T_{s0}->{s3}"
                    else:
                        land = "??"
                        coeff = None

            nonprim_value = 2 * eta_trans[(s0, s1)] * eta_trans[(s2, s3)]
            t1_str = f"T_{s0}->{s1}"
            t2_str = f"T_{s2}->{s3}"
            print(f"{t1_str + ', ' + t2_str:40} {land:25} {serialize_rational(coeff) if coeff is not None else 'NEW':12} {serialize_rational(nonprim_value):25}")

            nonprim_panel.append({
                "T1": t1_str,
                "T2": t2_str,
                "lands_in": land,
                "coeff": serialize_rational(coeff) if coeff is not None else "NEW_GENERATOR",
                "eta_T1": serialize_rational(eta_trans[(s0, s1)]),
                "eta_T2": serialize_rational(eta_trans[(s2, s3)]),
                "nonprimitive_eta_bivalue": serialize_rational(nonprim_value),
            })

    # ------------------------------------------------------------------
    # Check: do any chain compositions produce NEW transition generators
    # not in the basic kappa A enumeration? (i.e. two-step paths)
    # ------------------------------------------------------------------
    print()
    print("--- Two-step path check ---")
    two_step_paths = [r for r in nonprim_panel if "NEW" in r["lands_in"]]
    print(f"# of chain compositions landing on NEW (two-step) transitions: {len(two_step_paths)}")
    if two_step_paths:
        print("These would FORCE algebra-closure failure unless we include")
        print("two-step transitions in the basis.")
        for r in two_step_paths[:5]:
            print(f"  {r['T1']} · {r['T2']} -> {r['lands_in']}")
    else:
        print("All chain compositions stay within the (idempotents + single-step transitions) basis.")
        print("Algebra is closed under multiplication of the basic generators.")

    # ------------------------------------------------------------------
    # Comparison to Track B 1/65536 = 1/2^16 drift
    # ------------------------------------------------------------------
    print("\n--- Comparison to Track B drift signature ---")
    print("Track B: degree-3 closure residual = +-1/65536 = +-1/2^16 on (e_2, e_3) palindromes (n_max>=3)")
    print()
    print("The non-primitive content has the scale 2*eta(T1)*eta(T2).")
    print("eta(T) values are of form r/2^k with k in {6, 7} (i.e. 1/64 or 1/128).")
    print("So 2*eta(T1)*eta(T2) ~ 2 * (1/2^a) * (1/2^b) = 1/2^(a+b-1).")
    print()

    # Specifically, compute the (e_2, e_3) palindrome analog at the
    # Hopf-non-primitivity level
    # The (e_2, e_3) palindrome in Track B is the 4-tuple on sectors (2,0) (2,1)
    # In our enumeration: (2, 0) = sector index 2, (2, 1) = sector index 3.
    # The transitions between these are T_{(2,0)->(2,1)} and T_{(2,1)->(2,0)}.

    s2 = (2, 0)
    s3 = (2, 1)
    eta_T_23 = eta_trans.get((s2, s3))
    eta_T_32 = eta_trans.get((s3, s2))
    print(f"eta(T_{{{s2}->{s3}}}) = {serialize_rational(eta_T_23)}")
    print(f"eta(T_{{{s3}->{s2}}}) = {serialize_rational(eta_T_32)}")
    if eta_T_23 is not None and eta_T_32 is not None:
        nonprim_e2e3 = 2 * eta_T_23 * eta_T_32
        print(f"2 * eta(T_23) * eta(T_32) = {serialize_rational(nonprim_e2e3)}")
        # Compare to 1/65536
        target = Rational(1, 65536)
        ratio = nonprim_e2e3 / target if target != 0 else None
        print(f"Ratio to 1/65536 = {serialize_rational(ratio) if ratio is not None else 'undefined'}")
        print(f"1/65536 = {serialize_rational(Rational(1, 65536))} = 1/2^16")
        # Express nonprim_e2e3 in terms of 2^k
        if nonprim_e2e3 != 0:
            ab = abs(nonprim_e2e3)
            print(f"Non-prim content denom = {ab.q if isinstance(ab, sp.Rational) else 'symbolic'}")

    # ------------------------------------------------------------------
    # Lie algebra structure: which commutators are nonzero?
    # ------------------------------------------------------------------
    print("\n--- Commutator panel: [T_a, T_b] in the matrix algebra ---")
    print("Lie bracket of two transition generators.")
    print(f"{'T_a':25} {'T_b':25} {'[T_a, T_b]':15}")
    print("-" * 70)

    nz_commutators = []
    trans_keys = list(trans.keys())
    for i, ka in enumerate(trans_keys):
        for j, kb in enumerate(trans_keys):
            if i >= j:
                continue
            Ta, Tb = trans[ka], trans[kb]
            C = Ta * Tb - Tb * Ta
            if matrix_nonzero(C):
                ta_str = f"T_{ka[0]}->{ka[1]}"
                tb_str = f"T_{kb[0]}->{kb[1]}"
                nz_commutators.append((ta_str, tb_str))

    print(f"\n# of nonzero commutators [T_a, T_b]: {len(nz_commutators)} of {len(trans_keys)*(len(trans_keys)-1)//2}")
    if nz_commutators:
        print("First 10 nonzero commutators:")
        for (ta, tb) in nz_commutators[:10]:
            print(f"  [{ta}, {tb}] != 0")

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    results = {
        "n_max": 2,
        "nonprimitive_panel": nonprim_panel,
        "n_chain_compositions": len(nonprim_panel),
        "two_step_path_count": len(two_step_paths),
        "n_nonzero_commutators": len(nz_commutators),
        "nonzero_commutators_sample": nz_commutators[:20],
    }

    out_path = Path(__file__).parent / "data" / "sprint_q5p_offdiag_dirac_part2.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
