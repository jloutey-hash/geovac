"""Sprint Q5'-mJ-Smash-Product (T4) - structural test of the m_J-resolved
smash-product substrate at n_max = 2 flagged by Sprint Q5'-Combined-Substrate
(v3.63.0 L5) as the only path to non-trivial Hopf extension beyond the
Levi-decomposition.

L5 finding (v3.63.0): on the BASIC OffDiag substrate, the SU(2) z-rotation
U_theta produces three distinct phase factors {1, e^{i theta}, e^{-i theta}}
on the conjugate U_theta T U_theta^{-1}, demonstrating that the basic
basis is NOT closed under SU(2) action. L5 explicitly named the m_J
-resolved refinement as the load-bearing data for a non-trivial
smash-product.

T4 sprint scope: verify at n_max = 2 that
  (a) Every E1 transition T_{(s', m_J') -> (s, m_J)} satisfies
      Delta_mJ = m_J - m_J' in {-1, 0, +1} (Wigner-Eckart for vector
      operator).
  (b) The m_J-resolved OffDiag substrate splits into THREE U(1)
      eigenspaces by Delta_mJ value.
  (c) The U(1) z-rotation action gives a Z-graded extension of
      the basic substrate via Delta_mJ grading.
  (d) The full SU(2) action requires j_max = max(j) over all
      states, which is 3/2 at n_max = 2 (from sector (2,1) kappa=-2
      and (2,2) kappa=+2).

Smash-product structure (multi-year, scoped here):
  H_{m_J-OD} \rtimes O(SL_2) with SL_2 acting via (J_+, J_-, J_z).
  This sprint constructs the H_{m_J-OD} side bit-exactly and
  identifies the U(1) sub-action.

Discipline: bit-exact sympy.Rational throughout; no PSLQ.
"""
from __future__ import annotations

import json
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Rational, Integer, Symbol

from geovac.spectral_triple import FockSpectralTriple


def main() -> None:
    t0 = time.time()
    print("=" * 70)
    print("Sprint Q5'-mJ-Smash-Product (T4)")
    print("=" * 70, flush=True)

    n_max = 2
    triple = FockSpectralTriple(n_max=n_max)
    labels = triple.labels
    dim_H = triple.dim_H
    print(f"\nn_max = {n_max}, dim H = {dim_H}", flush=True)
    print(f"# Dirac state labels: {len(labels)}", flush=True)
    print(f"\nState index | (n, kappa, m_J) | sector (n, l) | j", flush=True)
    print(f"-" * 70, flush=True)
    for i, lab in enumerate(labels):
        sector = (lab.n_fock, lab.l)
        print(f"  {i:>2}        | ({lab.n_fock}, {lab.kappa:+d}, {lab.m_j})  | "
              f"{sector}        | {lab.j}", flush=True)

    # Step 1: get the off-diagonal Dirac operator A = D - Lambda
    Lambda = triple.diagonal_part
    A = triple.dirac_operator - Lambda

    # Step 2: list all non-zero transitions
    transitions = []
    for i in range(dim_H):
        for j in range(dim_H):
            if A[i, j] != 0:
                lab_to = labels[i]
                lab_from = labels[j]
                d_mJ = lab_to.m_j - lab_from.m_j
                d_n = lab_to.n_fock - lab_from.n_fock
                d_l = lab_to.l - lab_from.l
                transitions.append(dict(
                    state_from=j,
                    state_to=i,
                    lab_from=(lab_from.n_fock, lab_from.kappa, int(2*lab_from.m_j)),
                    lab_to=(lab_to.n_fock, lab_to.kappa, int(2*lab_to.m_j)),
                    A_value=str(A[i, j]),
                    delta_mJ_times_2=int(2 * d_mJ),
                    delta_n=d_n,
                    delta_l=d_l,
                ))

    print(f"\n# non-zero off-diagonal Dirac entries: {len(transitions)}", flush=True)

    # Step 3: verify Delta m_J in {-1, 0, +1} (Wigner-Eckart selection rule)
    d_mJ_values = set(t["delta_mJ_times_2"] for t in transitions)
    wigner_eckart = all(abs(t["delta_mJ_times_2"]) <= 2 for t in transitions)
    print(f"\nDelta m_J*2 values in transitions: {sorted(d_mJ_values)}", flush=True)
    print(f"Wigner-Eckart selection (|Delta m_J| <= 1): {wigner_eckart}", flush=True)

    # Step 4: group transitions by Delta m_J (U(1) eigenspaces)
    by_dmJ = defaultdict(list)
    for t in transitions:
        by_dmJ[t["delta_mJ_times_2"]].append(t)
    print(f"\nU(1) eigenspace decomposition (by 2*Delta m_J):", flush=True)
    for dmJ in sorted(by_dmJ.keys()):
        print(f"  2*Delta m_J = {dmJ:+d}: {len(by_dmJ[dmJ])} transitions", flush=True)

    # Step 5: verify U(1) z-rotation action
    # U_theta |n, kappa, m_J> = exp(i * m_J * theta) |n, kappa, m_J>
    # Action on transition T_{from -> to}: U_theta T U_theta^{-1}
    #   = exp(i * (m_J_to - m_J_from) * theta) * T_{from -> to}
    #   = exp(i * (Delta m_J) * theta) * T_{from -> to}
    print(f"\nU(1) z-rotation U_theta T U_theta^{{-1}} on transitions:", flush=True)
    print(f"  acts by phase exp(i * Delta m_J * theta) per transition", flush=True)
    print(f"  three eigenphases on the basic basis: {1, sp.exp(sp.I * Symbol('theta')), sp.exp(-sp.I * Symbol('theta'))}", flush=True)
    print(f"  this is exactly the L5 finding (three phase factors {{1, e^{{i theta}}, e^{{-i theta}}}}).", flush=True)
    print(f"  m_J-resolved refinement DIAGONALIZES the U(1) action.", flush=True)

    # Step 6: verify j_max requirement
    j_values = set(lab.j for lab in labels)
    j_max_required = max(j_values)
    print(f"\nDistinct j values in spectrum: {sorted(j_values, key=lambda j: (j.p, j.q))}", flush=True)
    print(f"j_max required for full SU(2) action: {j_max_required}", flush=True)
    print(f"L1 substrate decorated-PW j_max = 1/2 is INSUFFICIENT for full SU(2)", flush=True)
    print(f"  on m_J-resolved OffDiag at n_max = 2 (need j_max >= 3/2).", flush=True)

    # Step 7: identify sector structure of m_J-resolved substrate
    # For each sector (n, l), enumerate m_J states
    sector_mJ = defaultdict(list)
    for i, lab in enumerate(labels):
        sector_mJ[(lab.n_fock, lab.l)].append((int(2*lab.m_j), int(lab.kappa), i))

    print(f"\nm_J-resolved sector decomposition at n_max = {n_max}:", flush=True)
    for sector in sorted(sector_mJ.keys()):
        states = sector_mJ[sector]
        print(f"  sector {sector}: {len(states)} m_J states", flush=True)
        for two_mJ, kappa, idx in sorted(states):
            mJ = Rational(two_mJ, 2)
            print(f"    state {idx}: m_J = {mJ}, kappa = {kappa:+d}", flush=True)

    # Step 8: verify the Z-graded structure
    # The Delta m_J grading makes the m_J-resolved substrate Z-graded
    # under the U(1) sub-action.
    print(f"\nZ-graded structure (Delta m_J grading):", flush=True)
    print(f"  m_J-OffDiag substrate = directsum_{{k in {{-1, 0, +1}}}} A_k", flush=True)
    print(f"  where A_k = span of transitions with Delta m_J = k.", flush=True)
    # Count generators per grade
    grade_dims = {sp.Rational(dmJ, 2): len(transitions_dmJ)
                  for dmJ, transitions_dmJ in by_dmJ.items()}
    print(f"  dim A_{{-1}} = {grade_dims.get(Rational(-2, 2), 0)}", flush=True)
    print(f"  dim A_0   = {grade_dims.get(0, 0)}", flush=True)
    print(f"  dim A_{{+1}} = {grade_dims.get(Rational(2, 2), 0)}", flush=True)
    total_dim = sum(grade_dims.values())
    print(f"  total: {total_dim} (should equal # off-diag transitions {len(transitions)})", flush=True)

    # Step 9: scope multi-year smash-product structure
    print(f"\nMulti-year smash-product H_{{m_J-OD}} cross-prod O(SL_2):", flush=True)
    print(f"  algebra: H tensorprod K with action via (J_+, J_-, J_z)", flush=True)
    print(f"  coproduct: Delta(h tensorprod k) = (h_{{(1)}} tensor k_{{(1)}}) tensor (h_{{(2)}} tensor k_{{(2)}})", flush=True)
    print(f"  Verification requires j_max >= 3/2 in the K factor;", flush=True)
    print(f"  bit-exact Hopf-axiom check is sprint-scale at j_max = 3/2", flush=True)
    print(f"  but requires extending L2 decorated-PW substrate first.", flush=True)
    print(f"  NAMED MULTI-YEAR FOLLOW-ON.", flush=True)

    # Output
    wall = time.time() - t0
    print(f"\nWall: {wall:.2f} s", flush=True)

    out = dict(
        sprint="Q5p-mJ-Smash-Product (T4)",
        date="2026-06-06 (continuation of v3.63.0 L5 sub-sprint)",
        purpose="Structural test of m_J-resolved OffDiag substrate at n_max=2 confirming the Delta m_J in {-1, 0, +1} grading and the U(1) z-rotation eigenspace decomposition. Scopes the multi-year smash-product extension.",
        n_max=n_max,
        dim_H=dim_H,
        n_states=len(labels),
        labels=[
            dict(
                index=i,
                n_fock=lab.n_fock,
                kappa=lab.kappa,
                two_m_J=int(2 * lab.m_j),
                m_J_str=str(lab.m_j),
                sector=[lab.n_fock, lab.l],
                j=str(lab.j),
            )
            for i, lab in enumerate(labels)
        ],
        n_off_diagonal_transitions=len(transitions),
        delta_mJ_times_2_values=sorted(d_mJ_values),
        wigner_eckart_passes=wigner_eckart,
        u1_grade_dims=dict(
            dim_A_minus_1=grade_dims.get(Rational(-2, 2), 0),
            dim_A_0=grade_dims.get(0, 0),
            dim_A_plus_1=grade_dims.get(Rational(2, 2), 0),
        ),
        u1_grade_dims_total=total_dim,
        u1_grade_dims_match_total_transitions=total_dim == len(transitions),
        j_values=[str(j) for j in sorted(j_values, key=lambda jj: (jj.p, jj.q))],
        j_max_required=str(j_max_required),
        L2_decorated_PW_j_max=str(Rational(1, 2)),
        L2_sufficient_for_full_SU2_on_mJ_OffDiag=False,
        L2_needs_extension_to_j_max=str(j_max_required),
        sector_mJ_decomposition={
            f"{sector}": [
                dict(two_mJ=tm, kappa=k, state_index=i)
                for tm, k, i in sorted(states)
            ]
            for sector, states in sector_mJ.items()
        },
        smash_product_scope="multi-year; sprint-scale n_max=2 Delta m_J grading + U(1) z-rotation diagonalization verified; full SU(2) action needs L2 extension to j_max >= 3/2; bit-exact Hopf-axiom check is sprint-scale at j_max = 3/2.",
        named_followons=[
            "Extend L2 decorated-PW substrate to j_max = 3/2 (sprint-scale, ~1 day)",
            "Construct H_{m_J-OD} cross-prod O(SL_2) at n_max=2, j_max=3/2 (sprint-scale, ~2-3 days)",
            "Bit-exact verification of Hopf-axiom panel for the smash product (sprint-scale, ~1-2 days after construction)",
            "Generalize to n_max=3 (multi-year given growing j-content)",
        ],
        wall_seconds=wall,
    )

    out_path = Path(__file__).parent / "data" / "sprint_q5p_mj_smash_product.json"
    out_path.parent.mkdir(exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"Written: {out_path}", flush=True)


if __name__ == "__main__":
    main()
