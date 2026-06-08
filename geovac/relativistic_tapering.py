"""Relativistic-builder Z2 tapering scaffold (Tier 2 / Dirac-on-S^3).

This module is the relativistic-spinor analog of ``geovac.z2_tapering``.
It targets the relativistic chemistry Hamiltonians produced by
``geovac.composed_qubit_relativistic.build_composed_hamiltonian_relativistic``
where each ``DiracLabel(n_fock, kappa, two_m_j)`` occupies a SINGLE
Jordan-Wigner qubit (no separate alpha/beta spin doubling because the
Dirac spinor's spin is already encoded in the sign of kappa).

Sprint kappa-parity-Z2 (2026-06-07) verdict: NEGATIVE-STRUCTURAL
---------------------------------------------------------------
We hypothesised that the kappa-parity stabilizer

    P_kappa = prod_{q: kappa_q < 0} Z_q

would commute with the relativistic chemistry Hamiltonian
(``j = l + 1/2`` vs ``j = l - 1/2`` chirality branch count parity).
Direct numerical audit on LiH_rel, BeH_rel, CaH_rel at max_n=2 with
PK on/off and Breit on/off returned commutator residuals 1.7e-2 to 1.3e-1
(7-9 orders of magnitude ABOVE the 1e-10 PASS gate), confirming a
clean structural negative.

**Mechanism (LiH_rel diagnostic, ``debug/sprint_kappa_parity_mechanism.py``):**
38% of the ERI tuples (1,008 of 2,676) carry an ODD shift in
N_{kappa<0}. The jj-coupled full-Gaunt angular coefficient X_k allows
processes that connect p_{3/2} (kappa = -2) to p_{1/2} (kappa = +1) at
the same l = 1: the parity selection (l_a + l_c + k) even holds
trivially for same-l pairs, so the Coulomb operator can transfer one
electron between kappa-branches. The Phillips-Kleinman barrier
contributes the same kind of odd-Delta term in h1 because PK is
diagonal in (l, m_j) but l-fixed orbitals come in both kappa-branches.

This is the relativistic analog of the chemistry-side "wrong by
structure" finding (CLAUDE.md §3 W1e period-class, Sprint H1 Yukawa
non-selection): the kappa-branch label is NOT a conserved quantum number
of the Dirac-Coulomb Hamiltonian. Only j and m_j are.

Public API
----------
The module ships an audit-driven gate (mirroring
``z2_tapering.apply_hopf_tapering``) but with ``drop_noncommuting=True``
by default: any stabilizer that fails the commutation audit is silently
dropped. With kappa-parity this means the kept stabilizer list is
empty and the tapered Hamiltonian equals the naive Hamiltonian (no
qubit savings).

This module is therefore retained as institutional scaffolding for
future relativistic-symmetry probes (e.g. m_j-parity within a fixed-j
sector, or simultaneous j/m_j conservation as a more elaborate
non-Pauli-string stabilizer). It documents the negative result via the
audit gate rather than guessing.

References
----------
- Diagnostic: ``debug/sprint_kappa_parity_diagnostic.py``,
  ``debug/sprint_kappa_parity_mechanism.py``,
  ``debug/sprint_ch_kappa_parity_z2_memo.md``.
- Companion module: ``geovac.z2_tapering`` (non-relativistic Hopf-U(1)
  ``m_l -> -m_l`` reflection).
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    from openfermion import QubitOperator, jordan_wigner  # noqa: F401
    from openfermion.transforms import taper_off_qubits
    from openfermion.utils import commutator
    _HAS_OPENFERMION = True
except ImportError:  # pragma: no cover
    _HAS_OPENFERMION = False


__all__ = [
    "enumerate_relativistic_orbital_table",
    "build_kappa_parity_stabilizers",
    "audit_kappa_parity_commutation",
    "relativistic_tapered_from_spec",
]


# ---------------------------------------------------------------------------
# Relativistic orbital enumeration (mirrors build_composed_hamiltonian_relativistic)
# ---------------------------------------------------------------------------

def enumerate_relativistic_orbital_table(
    spec: Any,
) -> List[Tuple[Any, int, int, int]]:
    """Enumerate ``(sub_block_key, n_fock, kappa, two_m_j)`` for a relativistic spec.

    Order matches Jordan-Wigner qubit ordering used by
    ``geovac.composed_qubit_relativistic.build_composed_hamiltonian_relativistic``.

    Parameters
    ----------
    spec : MolecularSpec
        Must have ``relativistic=True`` (or be a relativistic spec factory
        output).  Standard composed builder fields (``blocks``,
        ``has_h_partner``, ``max_n_partner``, ``l_min``) are used.

    Returns
    -------
    list of (sb_key, n_fock, kappa, two_m_j) tuples.
    """
    from geovac.composed_qubit_relativistic import enumerate_dirac_labels

    table: List[Tuple[Any, int, int, int]] = []
    for blk_idx, blk in enumerate(spec.blocks):
        l_min = getattr(blk, "l_min", 0)
        center_labels = enumerate_dirac_labels(blk.max_n, l_min=l_min)
        sb_key_center = (blk_idx, blk.label, "center")
        for lab in center_labels:
            table.append((sb_key_center, lab.n_fock, lab.kappa, lab.two_m_j))
        if blk.has_h_partner:
            partner_max_n = (
                blk.max_n_partner if blk.max_n_partner > 0 else blk.max_n
            )
            partner_labels = enumerate_dirac_labels(partner_max_n, l_min=0)
            sb_key_partner = (blk_idx, blk.label, "partner")
            for lab in partner_labels:
                table.append((sb_key_partner, lab.n_fock, lab.kappa, lab.two_m_j))
    return table


# ---------------------------------------------------------------------------
# kappa-parity stabilizer construction
# ---------------------------------------------------------------------------

def _z_string_from_indices(qubit_indices: Sequence[int]) -> "QubitOperator":
    op = QubitOperator(())
    for q in qubit_indices:
        op *= QubitOperator(((int(q), "Z"),))
    return op


def build_kappa_parity_stabilizers(
    orbital_table: Sequence[Tuple[Any, int, int, int]],
    mode: str = "per_block",
) -> List["QubitOperator"]:
    """Construct the kappa-parity Z-string stabilizers.

    For each sub-block (or globally), the stabilizer is

        P_kappa = prod_{q: kappa_q < 0} Z_q

    counting electrons in the ``j = l + 1/2`` (Dirac kappa < 0) branch.

    Parameters
    ----------
    orbital_table : sequence of ``(sb_key, n_fock, kappa, two_m_j)``
        From :func:`enumerate_relativistic_orbital_table`.
    mode : 'global' or 'per_block'
        'global': one combined ``P_kappa`` over all ``kappa<0`` qubits.
        'per_block': one ``P_kappa^(i)`` per sub-block.

    Returns
    -------
    list of QubitOperator
        Empty if no ``kappa<0`` orbital exists.
    """
    if not _HAS_OPENFERMION:
        raise ImportError(
            "openfermion is required for build_kappa_parity_stabilizers. "
            "Install with: pip install openfermion"
        )
    if mode not in ("global", "per_block"):
        raise ValueError(f"mode must be 'global' or 'per_block', got {mode!r}")

    if mode == "global":
        neg_qubits = [
            i for i, (_, _, kappa, _) in enumerate(orbital_table) if kappa < 0
        ]
        if not neg_qubits:
            return []
        return [_z_string_from_indices(neg_qubits)]

    groups: Dict[Any, List[int]] = defaultdict(list)
    for i, (sb_key, _, kappa, _) in enumerate(orbital_table):
        if kappa < 0:
            groups[sb_key].append(i)
    ops: List["QubitOperator"] = []
    for sb_key in groups:
        qubits = groups[sb_key]
        if not qubits:
            continue
        ops.append(_z_string_from_indices(qubits))
    return ops


# ---------------------------------------------------------------------------
# Audit gate
# ---------------------------------------------------------------------------

def _commutes(
    op_a: "QubitOperator", op_b: "QubitOperator", atol: float = 1e-10
) -> Tuple[bool, float]:
    c = commutator(op_a, op_b)
    max_coef = max((abs(v) for v in c.terms.values()), default=0.0)
    return (max_coef < atol), float(max_coef)


def audit_kappa_parity_commutation(
    H: "QubitOperator",
    stabilizers: Sequence["QubitOperator"],
    atol: float = 1e-10,
) -> Tuple[List[bool], List[float]]:
    """Return per-stabilizer (passed, residual) lists for ``[H, P_i]``."""
    if not _HAS_OPENFERMION:
        raise ImportError(
            "openfermion is required for audit_kappa_parity_commutation."
        )
    passed: List[bool] = []
    residuals: List[float] = []
    for P_i in stabilizers:
        ok, res = _commutes(H, P_i, atol=atol)
        passed.append(ok)
        residuals.append(res)
    return passed, residuals


# ---------------------------------------------------------------------------
# Full pipeline from a relativistic MolecularSpec
# ---------------------------------------------------------------------------

def relativistic_tapered_from_spec(
    spec: Any,
    mode: str = "per_block",
    alpha_num: float = 7.2973525693e-3,
    pk_in_hamiltonian: bool = True,
    include_breit: bool = False,
    sector_signs: Optional[Sequence[int]] = None,
    drop_noncommuting: bool = True,
    audit_tol: float = 1e-10,
    verbose: bool = False,
) -> Dict[str, Any]:
    """End-to-end kappa-parity tapering pipeline for a relativistic spec.

    Builds the relativistic chemistry Hamiltonian via
    ``build_composed_hamiltonian_relativistic``, audits the kappa-parity
    stabilizer(s), drops any that fail the commutation gate (default:
    drop), and applies ``openfermion.transforms.taper_off_qubits`` with
    the standard alpha/beta parity Z-strings + the surviving
    kappa-parity stabilizers.

    On the relativistic builder, the standard alpha/beta parity tapering
    of the non-rel pipeline does NOT apply (there is no spin doubling).
    Only the kappa-parity stabilizer(s) — if they survive the audit —
    contribute qubit savings.

    Sprint kappa-parity-Z2 (2026-06-07) verdict: NEGATIVE-STRUCTURAL.
    The audit drops 100% of stabilizers on every relativistic spec
    tested; ``delta_Q`` is 0 and the returned tapered operator is the
    naive operator. The function is retained as institutional scaffolding
    and for future probes on related symmetries.

    Returns
    -------
    dict with keys:
      qubit_op_naive : QubitOperator
      qubit_op_tapered : QubitOperator  (equal to naive if all stabilizers
                                          dropped by audit)
      Q_naive, Q_tapered : int
      delta_Q : int
      n_sub_blocks : int
      n_kappa_neg : int
      orbital_table : list of (sb_key, n_fock, kappa, two_m_j)
      stabilizers_built : list of QubitOperator (before audit)
      stabilizers_passed : list of bool
      stabilizers_residuals : list of float
      stabilizers_kept : list of QubitOperator (post-audit)
      audit_verdict : 'KEEP_<n>' or 'DROP_ALL'
    """
    if not _HAS_OPENFERMION:
        raise ImportError(
            "openfermion is required for relativistic_tapered_from_spec."
        )
    from geovac.composed_qubit_relativistic import (
        build_composed_hamiltonian_relativistic,
    )

    result = build_composed_hamiltonian_relativistic(
        spec,
        alpha_num=alpha_num,
        pk_in_hamiltonian=pk_in_hamiltonian,
        include_breit=include_breit,
        verbose=verbose,
    )
    qop_naive = result["qubit_op"]
    Q_naive = int(result["Q"])

    orbital_table = enumerate_relativistic_orbital_table(spec)
    if len(orbital_table) != Q_naive:
        raise RuntimeError(
            f"orbital_table length {len(orbital_table)} != Q={Q_naive}; "
            "spec enumeration is out of sync with the relativistic builder"
        )

    stabs_built = build_kappa_parity_stabilizers(orbital_table, mode=mode)
    passed, residuals = audit_kappa_parity_commutation(
        qop_naive, stabs_built, atol=audit_tol
    )
    if drop_noncommuting:
        stabs_kept = [P for P, ok in zip(stabs_built, passed) if ok]
    else:
        stabs_kept = list(stabs_built)

    n_sub_blocks = len({sb for (sb, _, _, _) in orbital_table})
    n_kappa_neg = sum(1 for (_, _, kappa, _) in orbital_table if kappa < 0)

    if not stabs_kept:
        qop_tapered = qop_naive
        Q_tapered = Q_naive
        audit_verdict = "DROP_ALL"
    else:
        n_stab = len(stabs_kept)
        if sector_signs is None:
            signs = [1] * n_stab
        else:
            signs = list(sector_signs)
            if len(signs) != n_stab:
                raise ValueError(
                    f"sector_signs length {len(signs)} != number of kept "
                    f"stabilizers {n_stab}"
                )
        signed_stabs = [s * sgn for s, sgn in zip(stabs_kept, signs)]
        qop_tapered = taper_off_qubits(qop_naive, signed_stabs)
        Q_tapered = 0
        for term in qop_tapered.terms:
            for q, _ in term:
                if q + 1 > Q_tapered:
                    Q_tapered = q + 1
        audit_verdict = f"KEEP_{n_stab}"

    return {
        "qubit_op_naive": qop_naive,
        "qubit_op_tapered": qop_tapered,
        "Q_naive": Q_naive,
        "Q_tapered": Q_tapered,
        "delta_Q": Q_naive - Q_tapered,
        "n_sub_blocks": n_sub_blocks,
        "n_kappa_neg": n_kappa_neg,
        "orbital_table": orbital_table,
        "stabilizers_built": stabs_built,
        "stabilizers_passed": passed,
        "stabilizers_residuals": residuals,
        "stabilizers_kept": stabs_kept,
        "audit_verdict": audit_verdict,
        "mode": mode,
    }
