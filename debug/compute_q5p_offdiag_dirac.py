"""Sprint Q5'-OffDiag-Dirac — first scoping step of the multi-year cross-shell
off-diagonal Dirac perturbation substrate enrichment.

The v3.61.0 Track A finding established that the abelian primitive Hopf
candidate H_GV is forced by *sector-locality* of the cocycle classes chi_s
and eta_s (both depend only on the sector label (n, l), not on n_max).
Sector-locality is the structural condition that forces primitivity.

This sprint tests whether promoting the substrate to track the off-diagonal
kappa*A content explicitly (via transition operators e_{s' -> s} = e_s A e_{s'})
breaks sector-locality and yields non-primitive coproduct structure.

Driver discipline: bit-exact sympy.Rational throughout. No floats. No PSLQ.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Rational, Integer, Matrix, zeros as sp_zeros, eye as sp_eye

from geovac.spectral_triple import FockSpectralTriple


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def serialize_rational(x):
    if isinstance(x, sp.Rational):
        return f"{x.p}/{x.q}" if x.q != 1 else str(x.p)
    return str(x)


def serialize_matrix_entries(M: Matrix) -> List[List[str]]:
    return [[serialize_rational(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]


def matrix_nonzero(M: Matrix) -> bool:
    for i in range(M.rows):
        for j in range(M.cols):
            if M[i, j] != 0:
                return True
    return False


# ----------------------------------------------------------------------
# Step 1. Off-diagonal-tracking substrate construction
# ----------------------------------------------------------------------

def build_sector_idempotents(triple: FockSpectralTriple) -> Dict[Tuple[int, int], Matrix]:
    """Build diagonal sector projectors e_{(n, l)} on H."""
    N = triple.dim_H
    idem: Dict[Tuple[int, int], Matrix] = {}
    for sec in triple.sectors:
        e = sp_zeros(N, N)
        for i, k in enumerate(triple._state_to_sector):
            if triple.sectors[k] == sec:
                e[i, i] = Integer(1)
        idem[sec] = e
    return idem


def build_transition_generators(
    triple: FockSpectralTriple,
    idem: Dict[Tuple[int, int], Matrix],
) -> Dict[Tuple[Tuple[int, int], Tuple[int, int]], Matrix]:
    """Build transition generators e_{s' -> s} = e_s · A · e_{s'} where
    A is the (signed) off-diagonal Dirac perturbation kappa * A_graph.

    By construction:
      - e_{s' -> s} is zero if no edge connects sector s' to sector s.
      - Sector composition: e_{s1->s2} · e_{s3->s4} = delta_{s2, s3} e_{s1->s4}
        only via matrix multiplication on H (not the symbolic generator structure).

    We track only NON-ZERO transition generators (Gaunt-selected pairs).
    """
    D = triple.dirac_operator
    Lambda = triple.diagonal_part
    A_part = D - Lambda  # kappa * A_graph (the off-diagonal piece)

    trans: Dict[Tuple[Tuple[int, int], Tuple[int, int]], Matrix] = {}
    for s_from in triple.sectors:
        e_from = idem[s_from]
        for s_to in triple.sectors:
            if s_to == s_from:
                continue  # diagonal handled by idempotents
            e_to = idem[s_to]
            # e_{s' -> s} = e_to · A_part · e_from
            T = e_to * A_part * e_from
            if matrix_nonzero(T):
                trans[(s_from, s_to)] = T
    return trans


# ----------------------------------------------------------------------
# Step 3. Cocycle classes on the enriched substrate
# ----------------------------------------------------------------------

def chi_value(triple: FockSpectralTriple, M: Matrix) -> sp.Rational:
    """Compute chi(M) = Tr(gamma · M)."""
    g = triple.grading
    return sp.trace(g * M)


def eta_value(triple: FockSpectralTriple, M: Matrix) -> sp.Rational:
    """Compute eta(M) = Tr(gamma · D · M)."""
    g = triple.grading
    D = triple.dirac_operator
    return sp.trace(g * D * M)


def kappa_zero_check(triple: FockSpectralTriple, M: Matrix) -> sp.Rational:
    """Compute Tr(gamma · A_part · M) — the kappa*A-only contribution."""
    g = triple.grading
    Lambda = triple.diagonal_part
    D = triple.dirac_operator
    A_part = D - Lambda
    return sp.trace(g * A_part * M)


# ----------------------------------------------------------------------
# Step 4. Hopf candidate primitive-coproduct test on enriched substrate
# ----------------------------------------------------------------------

def matrix_product_in_algebra(
    M1: Matrix, M2: Matrix,
    idem: Dict[Tuple[int, int], Matrix],
    trans: Dict[Tuple[Tuple[int, int], Tuple[int, int]], Matrix],
) -> Dict[str, sp.Rational]:
    """Compute M1 * M2 in the enriched algebra and decompose it onto the
    basis {idempotents} ∪ {transitions}, returning the linear combination
    coefficients in a structured dict.

    Returns dict mapping generator-name -> coefficient (Rational).
    Residual not on the basis would indicate the algebra is non-closed
    (which would itself be a structural finding).
    """
    P = M1 * M2

    coeffs: Dict[str, sp.Rational] = {}

    # First, project onto idempotents: each e_s is rank-(dim_s) diagonal
    # with 1s on the sector states. We can extract the coefficient c_s by
    # picking one state's diagonal entry, since e_s is the unique generator
    # that has a nonzero entry on that state's diagonal.
    for s, e_s in idem.items():
        # Find the first nonzero diagonal index
        idx = None
        for i in range(e_s.rows):
            if e_s[i, i] != 0:
                idx = i
                break
        if idx is None:
            continue
        c = P[idx, idx]
        if c != 0:
            coeffs[f"e_{s}"] = c

    # Project onto transitions: e_{s'->s} has its nonzero entries at (i, j)
    # where state i is in sector s and state j is in sector s'.
    for (s_from, s_to), T in trans.items():
        # Find a nonzero entry index
        idx_pair = None
        for i in range(T.rows):
            for j in range(T.cols):
                if T[i, j] != 0:
                    idx_pair = (i, j, T[i, j])
                    break
            if idx_pair is not None:
                break
        if idx_pair is None:
            continue
        i_ref, j_ref, T_val = idx_pair
        P_val = P[i_ref, j_ref]
        if P_val == 0:
            continue
        # Coefficient on this transition basis element
        coeff = sp.Rational(P_val) / sp.Rational(T_val)
        if coeff != 0:
            coeffs[f"T_{s_from}->{s_to}"] = coeff

    return coeffs


def test_primitivity_forcing(
    triple: FockSpectralTriple,
    idem: Dict[Tuple[int, int], Matrix],
    trans: Dict[Tuple[Tuple[int, int], Tuple[int, int]], Matrix],
) -> Dict:
    """Test whether the matrix-algebra relations force non-primitive coproduct.

    Key principle: if x1 * x2 = x3 in the algebra (as matrix product), and
    we adopt PRIMITIVE coproduct Delta(x_g) = x_g ⊗ 1 + 1 ⊗ x_g on all
    generators including transitions, then multiplicativity gives

      Delta(x1 * x2) = Delta(x1) * Delta(x2)
                     = (x1 ⊗ 1 + 1 ⊗ x1) * (x2 ⊗ 1 + 1 ⊗ x2)
                     = x1 x2 ⊗ 1 + x1 ⊗ x2 + x2 ⊗ x1 + 1 ⊗ x1 x2.

    If x1 * x2 = x3 in the algebra, this becomes
      Delta(x3) = x3 ⊗ 1 + 1 ⊗ x3 + (x1 ⊗ x2 + x2 ⊗ x1)

    which is NOT primitive — it has non-primitive content (x1 ⊗ x2 + x2 ⊗ x1)
    iff x1 and x2 are non-trivially related to x3.

    So the test is: are there non-trivial multiplicative relations among
    the generators of the enriched algebra? If yes, primitive coproduct
    on generators CANNOT be extended to a multiplicative coproduct on the
    full algebra.
    """
    relations: List[Dict] = []

    # Sort generators for deterministic enumeration
    sector_keys = list(idem.keys())
    trans_keys = list(trans.keys())

    # (a) Idempotent × idempotent: e_s · e_t = delta_{s,t} e_s
    # These are EXISTING algebra relations; the diagonal part of the algebra
    # is commutative and each idempotent is a primitive generator (no
    # multiplicative relation to other idempotents that would induce
    # non-primitivity).
    for s in sector_keys:
        for t in sector_keys:
            if s == t:
                continue
            P = idem[s] * idem[t]
            if matrix_nonzero(P):
                relations.append({
                    "type": "idempotent_idempotent",
                    "factors": [f"e_{s}", f"e_{t}"],
                    "product": "nonzero (UNEXPECTED)",
                })

    # (b) Idempotent × transition: e_t · T_{s_from -> s_to}
    # Should equal T_{s_from -> s_to} if t = s_to, else zero
    for t in sector_keys:
        for (s_from, s_to), T in trans.items():
            P = idem[t] * T
            target_zero = (t != s_to)
            if target_zero:
                if matrix_nonzero(P):
                    relations.append({
                        "type": "idempotent_times_transition_should_be_zero",
                        "factors": [f"e_{t}", f"T_{s_from}->{s_to}"],
                        "expected": "zero", "result": "nonzero",
                    })
            else:
                # Should equal T
                if not (P - T).equals(sp_zeros(*P.shape)):
                    relations.append({
                        "type": "idempotent_times_transition_should_be_T",
                        "factors": [f"e_{t}", f"T_{s_from}->{s_to}"],
                        "expected": "T", "result": "different",
                    })

    # (c) Transition × Transition: composing transitions
    # T_{s2 -> s3} · T_{s0 -> s1} should be nonzero only if s1 == s2
    # (chain composition along the shell-transition graph).
    # This is THE LOAD-BEARING TEST for non-primitivity forcing.
    composition_relations: List[Dict] = []
    for (s0, s1), T1 in trans.items():
        for (s2, s3), T2 in trans.items():
            P = T2 * T1  # matrix product: T2 applied AFTER T1
            if matrix_nonzero(P):
                # If s1 != s2, this product should be zero by the
                # matrix-algebra structure (e_{s2} · e_{s1} = 0)
                if s1 != s2:
                    composition_relations.append({
                        "type": "transition_composition_nonchain",
                        "factors": [f"T_{s0}->{s1}", f"T_{s2}->{s3}"],
                        "expected": "zero (nonchain)", "result": "nonzero",
                    })
                else:
                    # CHAIN COMPOSITION: T_{s2 -> s3} · T_{s0 -> s1=s2}
                    # gives a nonzero product. The product lands in
                    # sector pair (s0, s3) — corresponds to either:
                    #   (i) a 2-step transition T_{s0 -> s3} (if s0 != s3),
                    #  (ii) a sector idempotent e_{s0} (if s0 == s3).
                    # This is the load-bearing structural relation.
                    decomp = matrix_product_in_algebra(T2, T1, idem, trans)
                    if decomp:
                        composition_relations.append({
                            "type": "transition_composition_chain",
                            "factors": [f"T_{s0}->{s1}", f"T_{s2}->{s3}"],
                            "intermediate_sector": str(s1),
                            "lands_in": list(decomp.keys()),
                            "decomp": {k: serialize_rational(v) for k, v in decomp.items()},
                        })

    return {
        "diagonal_relations_count": len(relations),
        "diagonal_relations": relations[:10],  # truncate for memo
        "composition_relations_count": len(composition_relations),
        "composition_relations": composition_relations,
    }


# ----------------------------------------------------------------------
# Main computation
# ----------------------------------------------------------------------

def run_nmax_2():
    print("=" * 70)
    print("Sprint Q5'-OffDiag-Dirac — n_max = 2")
    print("=" * 70)

    t0 = time.time()
    triple = FockSpectralTriple(n_max=2)
    print(f"dim H = {triple.dim_H}, n_sectors = {triple.n_sectors}")
    print(f"sectors = {triple.sectors}")
    print(f"kappa = {triple._kappa}")

    # ------------------------------------------------------------------
    # Step 1: Build off-diagonal-tracking substrate
    # ------------------------------------------------------------------
    idem = build_sector_idempotents(triple)
    trans = build_transition_generators(triple, idem)

    print(f"\n# of idempotent generators: {len(idem)}")
    print(f"# of transition generators: {len(trans)}")

    transition_list = []
    for (s_from, s_to), T in trans.items():
        nnz = sum(1 for i in range(T.rows) for j in range(T.cols) if T[i, j] != 0)
        transition_list.append({
            "from": str(s_from),
            "to": str(s_to),
            "nnz": nnz,
            "first_entry": serialize_rational(
                T[next((i for i in range(T.rows) if any(T[i, j] != 0 for j in range(T.cols))), 0),
                  next((j for j in range(T.cols) if any(T[i, j] != 0 for i in range(T.rows))), 0)]
            ) if nnz > 0 else "0",
        })

    print("\nTransition generators (sector_from -> sector_to, nnz, sample entry):")
    for t in transition_list:
        print(f"  {t['from']} -> {t['to']}: nnz={t['nnz']}, entry={t['first_entry']}")

    # ------------------------------------------------------------------
    # Step 3: Cocycle classes on enriched substrate
    # ------------------------------------------------------------------
    print("\n--- Cocycle class computation ---")

    cocycle_data: Dict[str, Dict[str, str]] = {}

    print("\n(a) Idempotents (recovery of v3.60.0 chi_s, eta_s)")
    for s, e in idem.items():
        chi_s = chi_value(triple, e)
        eta_s = eta_value(triple, e)
        kappa_eta = kappa_zero_check(triple, e)  # eta from kappa*A only
        print(f"  e_{s}: chi={chi_s}, eta={eta_s}, eta_kappa-only={kappa_eta}")
        cocycle_data[f"e_{s}"] = {
            "chi": serialize_rational(chi_s),
            "eta": serialize_rational(eta_s),
            "eta_kappa_only": serialize_rational(kappa_eta),
        }

    print("\n(b) Transitions (new — does off-diagonal content survive at class level?)")
    transition_cocycle_panel = {}
    for (s_from, s_to), T in trans.items():
        chi_T = chi_value(triple, T)
        eta_T = eta_value(triple, T)
        key = f"T_{s_from}->{s_to}"
        cocycle_data[key] = {
            "chi": serialize_rational(chi_T),
            "eta": serialize_rational(eta_T),
        }
        transition_cocycle_panel[key] = {
            "chi": chi_T,
            "eta": eta_T,
        }
        print(f"  T_{s_from}->{s_to}: chi={chi_T}, eta={eta_T}")

    # Count non-vanishing transition cocycle classes
    n_nonzero_chi = sum(1 for v in transition_cocycle_panel.values() if v["chi"] != 0)
    n_nonzero_eta = sum(1 for v in transition_cocycle_panel.values() if v["eta"] != 0)
    print(f"\n  Non-vanishing transition chi cocycles: {n_nonzero_chi} of {len(transition_cocycle_panel)}")
    print(f"  Non-vanishing transition eta cocycles: {n_nonzero_eta} of {len(transition_cocycle_panel)}")

    # ------------------------------------------------------------------
    # Step 4 + 5: Primitivity forcing test
    # ------------------------------------------------------------------
    print("\n--- Primitivity forcing test ---")
    primitivity_data = test_primitivity_forcing(triple, idem, trans)
    print(f"# of diagonal-relations problems: {primitivity_data['diagonal_relations_count']}")
    print(f"# of composition (chain) relations: {primitivity_data['composition_relations_count']}")

    # Detail the chain compositions: each one is a structural algebra relation
    # x1 * x2 = (linear combination) that FORCES non-primitive coproduct.
    chain_relations_summary = []
    for rel in primitivity_data["composition_relations"]:
        if rel.get("type") == "transition_composition_chain":
            chain_relations_summary.append(rel)

    print(f"\n# of non-trivial chain-composition relations: {len(chain_relations_summary)}")
    if chain_relations_summary:
        print("\nSample chain-composition relations (first 5):")
        for r in chain_relations_summary[:5]:
            print(f"  {r['factors'][0]} · {r['factors'][1]}")
            print(f"    intermediate sector: {r['intermediate_sector']}")
            print(f"    lands in: {r['lands_in']}")
            print(f"    decomp: {r['decomp']}")

    # ------------------------------------------------------------------
    # Step 6: Compare to Track B drift residual
    # ------------------------------------------------------------------
    print("\n--- Comparison to v3.61.0 Track B drift ---")
    print("Track B found a +- 1/65536 closure drift at degree 3 on (e_2, e_3)")
    print("palindromes (fixed point at n_max >= 3).")
    print("Track 1 found: kappa A contribution to Tr(gamma A · e_s) vanishes per-sector.")
    print()
    print("The off-diagonal Dirac matrix elements connecting shell n=1 to n=2 are:")
    Lambda = triple.diagonal_part
    A_part = triple.dirac_operator - Lambda
    print(f"  Sample A_part entries kappa values:")
    nonzero_count = 0
    sample_entries = []
    for i in range(A_part.rows):
        for j in range(A_part.cols):
            if A_part[i, j] != 0:
                lab_i = triple.labels[i]
                lab_j = triple.labels[j]
                s_i = triple._state_to_sector[i]
                s_j = triple._state_to_sector[j]
                if nonzero_count < 5:
                    print(f"    A[{i},{j}] = {A_part[i,j]} ({triple.sectors[s_i]} <- {triple.sectors[s_j]})")
                sample_entries.append({
                    "i": i, "j": j, "value": serialize_rational(A_part[i, j]),
                    "from_sector": str(triple.sectors[s_j]),
                    "to_sector": str(triple.sectors[s_i]),
                })
                nonzero_count += 1
    print(f"  Total off-diagonal nonzero entries in kappa*A: {nonzero_count}")

    # ------------------------------------------------------------------
    # Step 7: Lowest-order commutator in candidate Hopf
    # ------------------------------------------------------------------
    print("\n--- Lowest-order non-trivial commutator (matrix-algebra) ---")
    # Test commutator [T1, T2] on a sample chain pair
    sample_commutator = None
    if chain_relations_summary:
        first_chain = chain_relations_summary[0]
        s0_s1_str = first_chain['factors'][0].replace("T_", "")
        s2_s3_str = first_chain['factors'][1].replace("T_", "")
        # Parse "(1, 0)->(2, 0)"
        # We'll just compute [T1, T2] for the first two listed transitions
        trans_keys_list = list(trans.keys())
        if len(trans_keys_list) >= 2:
            k1 = trans_keys_list[0]
            k2 = trans_keys_list[1]
            T1 = trans[k1]
            T2 = trans[k2]
            comm = T1 * T2 - T2 * T1
            comm_nz = matrix_nonzero(comm)
            print(f"  [T_{k1[0]}->{k1[1]}, T_{k2[0]}->{k2[1]}] is {'nonzero' if comm_nz else 'zero'}")
            sample_commutator = {
                "pair": [f"T_{k1[0]}->{k1[1]}", f"T_{k2[0]}->{k2[1]}"],
                "nonzero": comm_nz,
            }

    # ------------------------------------------------------------------
    # Assemble results
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    results = {
        "sprint": "Q5'-OffDiag-Dirac",
        "date": "2026-06-05",
        "n_max": 2,
        "dim_H": triple.dim_H,
        "n_sectors": triple.n_sectors,
        "sectors": [str(s) for s in triple.sectors],
        "kappa": serialize_rational(triple._kappa),
        "n_idempotents": len(idem),
        "n_transitions": len(trans),
        "transitions": transition_list,
        "cocycle_classes": cocycle_data,
        "n_nonzero_chi_transition": n_nonzero_chi,
        "n_nonzero_eta_transition": n_nonzero_eta,
        "primitivity_test": {
            "n_diagonal_relations_problems": primitivity_data["diagonal_relations_count"],
            "n_chain_composition_relations": len(chain_relations_summary),
            "chain_relations_sample": chain_relations_summary[:10],
        },
        "off_diagonal_A_part": {
            "total_nonzero_entries": nonzero_count,
            "sample_entries": sample_entries[:10],
        },
        "sample_commutator": sample_commutator,
        "wall_time_s": elapsed,
    }
    return results


def run_nmax_check_at_higher_dim():
    """Quick check at n_max=3 for transition-generator count scaling."""
    print("\n" + "=" * 70)
    print("Sprint Q5'-OffDiag-Dirac — n_max = 3 (scaling sanity check)")
    print("=" * 70)

    triple = FockSpectralTriple(n_max=3)
    print(f"dim H = {triple.dim_H}, n_sectors = {triple.n_sectors}")
    idem = build_sector_idempotents(triple)
    trans = build_transition_generators(triple, idem)
    print(f"# of idempotents: {len(idem)}")
    print(f"# of transitions: {len(trans)}")

    # Quick scan: does any transition survive at the eta-class level?
    n_nz_chi = 0
    n_nz_eta = 0
    for T in trans.values():
        if chi_value(triple, T) != 0:
            n_nz_chi += 1
        if eta_value(triple, T) != 0:
            n_nz_eta += 1
    print(f"Non-vanishing transition chi: {n_nz_chi} of {len(trans)}")
    print(f"Non-vanishing transition eta: {n_nz_eta} of {len(trans)}")

    return {
        "n_max": 3,
        "dim_H": triple.dim_H,
        "n_sectors": triple.n_sectors,
        "n_idempotents": len(idem),
        "n_transitions": len(trans),
        "n_nonzero_chi_transition": n_nz_chi,
        "n_nonzero_eta_transition": n_nz_eta,
    }


def main():
    results_nmax2 = run_nmax_2()
    results_nmax3 = run_nmax_check_at_higher_dim()

    final = {
        "n_max_2": results_nmax2,
        "n_max_3_sanity": results_nmax3,
    }

    out_path = Path(__file__).parent / "data" / "sprint_q5p_offdiag_dirac.json"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(final, f, indent=2, default=str)

    print(f"\nWrote {out_path}")
    print(f"Total wall time: {results_nmax2['wall_time_s']:.2f} s")


if __name__ == "__main__":
    main()
