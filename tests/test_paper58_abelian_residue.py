"""Backing tests for Paper 58 (Angular Sparsity Is an Atomic-Sector Property).

Covers the two observations MEASURED for the first time in the Paper-58 work,
which exist in no other backing:

  test_paper58_no_angle
      Obs. "The composed builder carries no angular geometry".  The paper's
      H2O discussion, and the claim that Corollary 1 is not well posed in the
      bent case, both rest on this.  If someone later adds angular geometry to
      the builder, this test SHOULD fail -- that is the point.

  test_paper58_swap_cost
      Obs. "Permutation symmetry costs qubits on this builder".  Enabling
      equivalent-atom-swap tapering REDUCES the qubit saving and inflates the
      Pauli count.  This also pins the documentation discrepancy: the
      docstring of extended_tapered_from_spec advertises "BeH2 +1"; the
      measured value is -1 (Hopf-only baseline) / -3 (Hopf + ell baseline).

Backfill still owed for Paper 58 (census, bra/ket certificate, m-rule theorem,
NaH ladder): those live in exploratory drivers and must be migrated
self-contained into tests/ before the paper is a finished artifact, per the
cite-permanent-records policy.  See the paper's backing table.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

pytest.importorskip("openfermion",
                    reason="tapering pipeline requires openfermion")


# ---------------------------------------------------------------------------
# Obs: the composed builder carries no angular geometry
# ---------------------------------------------------------------------------

def test_paper58_no_angle():
    """molecular_spec encodes no bond angle; H2O is radial-only.

    Three independent legs, so a single cosmetic edit cannot silently
    flip the result:
      (a) no angular vocabulary anywhere in the module source;
      (b) h2o_spec delegates to the generic hydride factory;
      (c) the built spec exposes a single R and no nuclei.
    """
    from geovac import molecular_spec as ms

    src = Path(ms.__file__).read_text(encoding="utf-8", errors="replace")

    # (a) No angular vocabulary. 'angle' as a whole word, degree symbols,
    #     or the water bond angle in either convention.
    angular_hits = re.findall(
        r"\bangle\b|\btheta\b|\bbent\b|104\.[0-9]|\bdeg\b", src, re.IGNORECASE)
    assert angular_hits == [], (
        f"molecular_spec.py now contains angular vocabulary {angular_hits!r}; "
        "Paper 58 Obs. 'no angular geometry' and the not-well-posed status of "
        "Corollary 1 in the bent case both need revisiting."
    )

    # (b) h2o_spec is the generic hydride factory, i.e. Z + one distance.
    h2o_src = src[src.index("def h2o_spec"):]
    h2o_body = h2o_src[:h2o_src.index("def hf_spec")]
    assert "hydride_spec(8" in h2o_body, (
        "h2o_spec no longer delegates to hydride_spec(8); its geometry "
        "parameterization has changed."
    )

    # (c) The built spec has one radial parameter and no nuclear positions.
    spec = ms.h2o_spec()
    assert spec.nuclei is None, (
        "h2o_spec now supplies nuclei; C2v may now be representable and "
        "Paper 58 Corollary 1 may have become testable."
    )
    assert isinstance(spec.R, float) and spec.R > 0
    # Five blocks (core, two bond pairs, two lone pairs), all sharing one R.
    assert len(spec.blocks) == 5, f"expected 5 H2O blocks, got {len(spec.blocks)}"


# ---------------------------------------------------------------------------
# Obs: permutation symmetry costs qubits
# ---------------------------------------------------------------------------

def _linear_symmetric_nuclei(Z_c: float, Z_o: float, R: float):
    """X-A-X on the z-axis.  Linear, so angle-free and unambiguous."""
    return [
        {"Z": float(Z_c), "position": (0.0, 0.0, 0.0), "label": "A"},
        {"Z": float(Z_o), "position": (0.0, 0.0, float(R)), "label": "X1"},
        {"Z": float(Z_o), "position": (0.0, 0.0, -float(R)), "label": "X2"},
    ]


def _n_pauli(qubit_op) -> int:
    return sum(1 for term in qubit_op.terms if term)


def _taper(spec, nuclei, *, ell: bool, swap: bool):
    from geovac.extended_tapering import extended_tapered_from_spec
    out = extended_tapered_from_spec(
        spec, use_hopf=True, use_ell_parity=ell,
        use_atom_swap=swap, use_inversion=False, nuclei=nuclei,
    )
    return out["delta_Q"], _n_pauli(out["qubit_op_tapered"])


@pytest.mark.slow
def test_paper58_swap_cost():
    """Atom-swap tapering reduces delta_Q and inflates Pauli count (BeH2).

    Guards the sign, not a fitted magnitude: swap must make delta_Q strictly
    WORSE and the Pauli count strictly larger, under both baselines.
    """
    from geovac import molecular_spec as ms

    spec = ms.beh2_spec()
    nuclei = _linear_symmetric_nuclei(4.0, 1.0, spec.R)

    # Sanity: the geometry really is swap-eligible, else the test is vacuous.
    from geovac.extended_tapering import (
        find_equivalent_atom_pairs, is_centrosymmetric,
    )
    assert is_centrosymmetric(nuclei), "BeH2 test geometry is not centrosymmetric"
    assert len(find_equivalent_atom_pairs(spec, nuclei)) >= 1, (
        "no equivalent atom pairs found; the swap leg would be vacuous"
    )

    for ell in (False, True):
        dq_base, pauli_base = _taper(spec, nuclei, ell=ell, swap=False)
        dq_swap, pauli_swap = _taper(spec, nuclei, ell=ell, swap=True)

        label = "hopf+ell" if ell else "hopf-only"
        assert dq_swap < dq_base, (
            f"[{label}] atom-swap did NOT cost qubits: "
            f"delta_Q {dq_base} -> {dq_swap}.  Paper 58 Obs. "
            "'permutation symmetry costs qubits' and the documentation "
            "discrepancy it reports would both need revisiting."
        )
        assert pauli_swap > pauli_base, (
            f"[{label}] atom-swap did NOT inflate Pauli count: "
            f"{pauli_base} -> {pauli_swap}"
        )

    # Pin the documentation discrepancy explicitly: the docstring claims +1.
    dq_base, _ = _taper(spec, nuclei, ell=False, swap=False)
    dq_swap, _ = _taper(spec, nuclei, ell=False, swap=True)
    assert dq_swap - dq_base < 0, (
        "extended_tapered_from_spec docstring advertises BeH2 atom-swap as "
        f"+1 qubit; measured {dq_swap - dq_base:+d} on the Hopf-only baseline."
    )


@pytest.mark.slow
def test_paper58_swap_null_control():
    """LiH has no equivalent atoms, so the swap flag must be a no-op.

    Without this, test_paper58_swap_cost could be passing because the swap
    path is broken in general rather than because permutation costs qubits.
    """
    from geovac import molecular_spec as ms

    spec = ms.lih_spec()
    nuclei = [
        {"Z": 3.0, "position": (0.0, 0.0, 0.0), "label": "Li"},
        {"Z": 1.0, "position": (0.0, 0.0, float(spec.R)), "label": "H"},
    ]
    from geovac.extended_tapering import find_equivalent_atom_pairs
    assert len(find_equivalent_atom_pairs(spec, nuclei)) == 0

    dq_base, pauli_base = _taper(spec, nuclei, ell=True, swap=False)
    dq_swap, pauli_swap = _taper(spec, nuclei, ell=True, swap=True)
    assert (dq_base, pauli_base) == (dq_swap, pauli_swap), (
        "swap flag changed the result for a molecule with no equivalent "
        f"atoms: delta_Q {dq_base}->{dq_swap}, Pauli {pauli_base}->{pauli_swap}"
    )
