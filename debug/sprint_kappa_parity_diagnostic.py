"""Sprint kappa-parity Z2 tapering diagnostic on relativistic Tier 2 chemistry.

Tests whether P_kappa = (-1)^{N_{kappa<0}} (Z-string on kappa<0 spinor
orbitals) commutes with the relativistic chemistry Hamiltonian
H_rel built by composed_qubit_relativistic.build_composed_hamiltonian_relativistic.

In the spinor JW encoding each DiracLabel(n, kappa, two_m_j) occupies a
SINGLE qubit (no separate alpha/beta spin doubling), so the kappa-parity
stabilizer is

    P_kappa^{block} = prod_{q in block, kappa_q < 0} Z_q

We check (a) global P_kappa over all spin-orbitals, (b) per-sub-block.

GO if commutator residual < 1e-10; PARTIAL if global-only; STOP otherwise.
"""

from __future__ import annotations

import time
from collections import defaultdict
from typing import Any, Dict, List, Tuple

from openfermion import QubitOperator
from openfermion.utils import commutator

from geovac.composed_qubit_relativistic import build_composed_hamiltonian_relativistic
from geovac.dirac_matrix_elements import DiracLabel
from geovac.molecular_spec import (
    lih_spec_relativistic,
    beh_spec_relativistic,
    cah_spec_relativistic,
)


def enumerate_relativistic_orbital_table(
    spec, alpha_num: float = 7.2973525693e-3
) -> Tuple[List[Tuple[Any, int, int, int]], int]:
    """Enumerate (sub_block_key, n_fock, kappa, two_m_j) for a relativistic spec.

    Mirrors ``build_composed_hamiltonian_relativistic`` Phase 1 enumeration.
    """
    from geovac.composed_qubit_relativistic import enumerate_dirac_labels

    table: List[Tuple[Any, int, int, int]] = []
    for blk_idx, blk in enumerate(spec.blocks):
        l_min = getattr(blk, 'l_min', 0)
        center_labels = enumerate_dirac_labels(blk.max_n, l_min=l_min)
        sb_key_center = (blk_idx, blk.label, 'center')
        for lab in center_labels:
            table.append((sb_key_center, lab.n_fock, lab.kappa, lab.two_m_j))
        if blk.has_h_partner:
            partner_max_n = blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            partner_labels = enumerate_dirac_labels(partner_max_n, l_min=0)
            sb_key_partner = (blk_idx, blk.label, 'partner')
            for lab in partner_labels:
                table.append((sb_key_partner, lab.n_fock, lab.kappa, lab.two_m_j))
    return table, len(table)


def build_kappa_parity_stabilizers(
    orbital_table: List[Tuple[Any, int, int, int]],
    mode: str = "global",
) -> List[QubitOperator]:
    """Build the kappa-parity Z-string stabilizers.

    In the spinor JW encoding, each DiracLabel is ONE qubit (no alpha/beta doubling).
    So ``P_kappa = ?_{q: kappa_q < 0} Z_q``.

    Parameters
    ----------
    orbital_table : list of (sb_key, n, kappa, two_m_j)
        Order MUST match the relativistic JW assembly order.
    mode : 'global' or 'per_block'
        'global': one P_kappa over ALL kappa<0 qubits.
        'per_block': one P_kappa per sub-block with at least one kappa<0 orbital.

    Returns
    -------
    list of QubitOperator
    """
    if mode == "global":
        neg_qubits = [i for i, (_, _, kappa, _) in enumerate(orbital_table) if kappa < 0]
        if not neg_qubits:
            return []
        op = QubitOperator(())
        for q in neg_qubits:
            op *= QubitOperator(((int(q), "Z"),))
        return [op]
    elif mode == "per_block":
        groups: Dict[Any, List[int]] = defaultdict(list)
        for i, (sb_key, _, kappa, _) in enumerate(orbital_table):
            if kappa < 0:
                groups[sb_key].append(i)
        ops: List[QubitOperator] = []
        for sb_key, qubits in groups.items():
            if not qubits:
                continue
            op = QubitOperator(())
            for q in qubits:
                op *= QubitOperator(((int(q), "Z"),))
            ops.append(op)
        return ops
    else:
        raise ValueError(f"unknown mode {mode!r}")


def commutator_residual(
    H: QubitOperator, P: QubitOperator
) -> float:
    """Return max |coefficient| of [H, P]."""
    c = commutator(H, P)
    if not c.terms:
        return 0.0
    return max(abs(v) for v in c.terms.values())


def run_diagnostic(
    spec_factory,
    spec_name: str,
    max_n: int = 2,
    include_pk: bool = True,
    include_breit: bool = False,
    alpha_num: float = 7.2973525693e-3,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Run kappa-parity commutation audit on a relativistic spec."""
    print(f"\n{'='*70}")
    print(f"  kappa-parity commutator audit: {spec_name}")
    print(f"  max_n={max_n}, PK={include_pk}, Breit={include_breit}")
    print(f"{'='*70}")

    t0 = time.perf_counter()
    spec = spec_factory(max_n=max_n, include_pk=include_pk)
    result = build_composed_hamiltonian_relativistic(
        spec, alpha_num=alpha_num, pk_in_hamiltonian=include_pk,
        verbose=False, include_breit=include_breit,
    )
    H = result["qubit_op"]
    Q = result["Q"]
    N_pauli = result["N_pauli"]
    t_build = time.perf_counter() - t0
    print(f"  Hamiltonian built: Q={Q}, N_pauli={N_pauli}, wall={t_build:.2f}s")

    orbital_table, Q_check = enumerate_relativistic_orbital_table(spec)
    assert Q_check == Q, f"orbital table length {Q_check} != Q {Q}"

    n_kappa_neg_total = sum(1 for (_, _, kappa, _) in orbital_table if kappa < 0)
    n_kappa_pos_total = Q - n_kappa_neg_total
    n_sub_blocks = len({sb for (sb, _, _, _) in orbital_table})
    print(f"  Orbital table: {n_kappa_neg_total} kappa<0 + {n_kappa_pos_total} kappa>0 "
          f"across {n_sub_blocks} sub-block(s)")

    # ---- Test 1: global P_kappa ----
    global_stabs = build_kappa_parity_stabilizers(orbital_table, mode="global")
    if not global_stabs:
        print("  [global] no kappa<0 orbitals; stabilizer trivial ? SKIP")
        residual_global = 0.0
    else:
        P_global = global_stabs[0]
        # Hermiticity sanity
        n_terms_P = len(P_global.terms)
        # All Pauli Z-string coefficients are real, so Hermitian by construction
        t1 = time.perf_counter()
        residual_global = commutator_residual(H, P_global)
        t_comm = time.perf_counter() - t1
        verdict = "PASS" if residual_global < 1e-10 else (
            "BORDERLINE" if residual_global < 1e-6 else "FAIL")
        print(f"  [global]  P_kappa on {n_kappa_neg_total} qubits, "
              f"|[H,P]|_max = {residual_global:.3e}  [{verdict}]  "
              f"(comm wall={t_comm:.2f}s)")

    # ---- Test 2: per-sub-block P_kappa ----
    per_block_stabs = build_kappa_parity_stabilizers(orbital_table, mode="per_block")
    print(f"  [per_block] {len(per_block_stabs)} sub-block stabilizer(s)")
    per_block_residuals: List[float] = []
    per_block_pass: List[bool] = []
    for i, P_i in enumerate(per_block_stabs):
        n_q = len(next(iter(P_i.terms)))
        res = commutator_residual(H, P_i)
        passed = res < 1e-10
        per_block_residuals.append(res)
        per_block_pass.append(passed)
        verdict = "PASS" if passed else (
            "BORDERLINE" if res < 1e-6 else "FAIL")
        print(f"    P_kappa^({i}): {n_q} qubits, |[H,P_i]|_max = {res:.3e}  [{verdict}]")

    # ---- Verdict ----
    if global_stabs and residual_global < 1e-10:
        if per_block_stabs and all(per_block_pass):
            verdict = "GO_PER_BLOCK"
            delta_Q = len(per_block_stabs)
        else:
            verdict = "GO_GLOBAL"
            delta_Q = 1
    elif global_stabs and residual_global < 1e-6:
        verdict = "BORDERLINE_GLOBAL"
        delta_Q = 0
    elif not global_stabs:
        verdict = "TRIVIAL_NO_KAPPA_NEG"
        delta_Q = 0
    else:
        verdict = "STOP"
        delta_Q = 0

    print(f"  ==> VERDICT: {verdict}  (?Q savings = {delta_Q})")

    return {
        "spec_name": spec_name,
        "max_n": max_n,
        "include_pk": include_pk,
        "include_breit": include_breit,
        "Q": Q,
        "N_pauli": N_pauli,
        "n_kappa_neg": n_kappa_neg_total,
        "n_kappa_pos": n_kappa_pos_total,
        "n_sub_blocks": n_sub_blocks,
        "residual_global": residual_global,
        "per_block_residuals": per_block_residuals,
        "n_per_block": len(per_block_stabs),
        "verdict": verdict,
        "delta_Q": delta_Q,
        "wall_build_s": t_build,
    }


def main() -> None:
    panel: List[Dict[str, Any]] = []

    # Primary test: LiH-rel max_n=2 with and without PK; with and without Breit
    panel.append(run_diagnostic(lih_spec_relativistic, "LiH_rel",
                                max_n=2, include_pk=True, include_breit=False))
    panel.append(run_diagnostic(lih_spec_relativistic, "LiH_rel",
                                max_n=2, include_pk=False, include_breit=False))
    panel.append(run_diagnostic(lih_spec_relativistic, "LiH_rel",
                                max_n=2, include_pk=True, include_breit=True))

    # Secondary: BeH-rel
    try:
        panel.append(run_diagnostic(beh_spec_relativistic, "BeH_rel",
                                    max_n=2, include_pk=True,
                                    include_breit=False))
    except Exception as e:
        print(f"\n[BeH_rel] FAILED to build: {type(e).__name__}: {e}")

    # Tertiary: CaH-rel (heavier Z, frozen core, larger Q)
    try:
        panel.append(run_diagnostic(cah_spec_relativistic, "CaH_rel",
                                    max_n=2, include_pk=True,
                                    include_breit=False))
    except Exception as e:
        print(f"\n[CaH_rel] FAILED to build: {type(e).__name__}: {e}")

    print("\n" + "="*70)
    print(" SUMMARY PANEL")
    print("="*70)
    print(f"  {'spec':<14}{'PK':>4}{'Brt':>5}{'Q':>5}{'N_P':>7}"
          f"{'kappa<0':>5}{'n_sb':>5}{'|[H,P]|_glob':>14}{'verdict':>22}")
    for row in panel:
        print(f"  {row['spec_name']:<14}"
              f"{'Y' if row['include_pk'] else 'N':>4}"
              f"{'Y' if row['include_breit'] else 'N':>5}"
              f"{row['Q']:>5}{row['N_pauli']:>7}"
              f"{row['n_kappa_neg']:>5}{row['n_sub_blocks']:>5}"
              f"{row['residual_global']:>14.3e}"
              f"{row['verdict']:>22}")

    import json, os
    os.makedirs("debug/data", exist_ok=True)
    out_path = "debug/data/kappa_parity_diagnostic.json"
    with open(out_path, "w") as f:
        json.dump(panel, f, indent=2)
    print(f"\n  Saved {out_path}")


if __name__ == "__main__":
    main()
