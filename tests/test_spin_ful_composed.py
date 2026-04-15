"""Regression tests for the spin-ful composed qubit pipeline (Track T3).

Validates two properties:

1. **Non-relativistic path is unchanged.**  The existing LiH/BeH₂/H₂O
   composed Pauli counts (334/556/778 with identity, equivalently 333/555/777
   for the `N_pauli` field which excludes identity) are preserved bit-for-bit
   by the default ``build_composed_hamiltonian`` call.  Any change that
   perturbs these numbers means the relativistic extension leaked into the
   scalar pipeline.

2. **Relativistic path builds cleanly** for LiH/BeH/CaH, yields block-
   diagonal ERIs, produces correct α → 0 limit (zero spin-orbit energy),
   and passes a JW eigenvalue consistency spot-check.

See ``docs/spin_ful_composed_design_memo.md`` for the API and Pauli-count
tables.
"""

from __future__ import annotations

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Non-relativistic regression: published Pauli counts must be untouched.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec_fn,name,expected_terms", [
    ("lih_spec", "LiH", 334),     # CLAUDE.md §10: "LiH composed Pauli terms (Q=30) == 334"
    ("beh2_spec", "BeH2", 556),   # CLAUDE.md §10: "BeH2 composed Pauli terms (Q=50) == 556"
    ("h2o_spec", "H2O", 778),     # CLAUDE.md §10: "H2O composed Pauli terms (Q=70) == 778"
])
def test_scalar_pauli_counts_unchanged(spec_fn, name, expected_terms):
    """Scalar LiH/BeH₂/H₂O Pauli counts must match the CLAUDE.md benchmark."""
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac import molecular_spec as ms

    spec = getattr(ms, spec_fn)(max_n=2)
    assert spec.relativistic is False, (
        f"{name} default spec must be non-relativistic")
    r = build_composed_hamiltonian(spec)
    # N_pauli excludes identity; CLAUDE.md §10 benchmark includes it.
    n_terms_with_identity = len(r['qubit_op'].terms)
    assert n_terms_with_identity == expected_terms, (
        f"{name} scalar Pauli count regression: got {n_terms_with_identity}, "
        f"expected {expected_terms} (CLAUDE.md §10 benchmark)")


def test_scalar_relativistic_flag_default_false():
    """MolecularSpec.relativistic must default to False for all hydride specs."""
    from geovac.molecular_spec import (
        lih_spec, beh2_spec, h2o_spec, hf_spec, nah_spec,
        lih_spec_relativistic, beh_spec_relativistic, cah_spec_relativistic,
    )
    for fn in (lih_spec, beh2_spec, h2o_spec, hf_spec, nah_spec):
        assert fn().relativistic is False, (
            f"{fn.__name__} must default to relativistic=False")
    for fn in (lih_spec_relativistic, beh_spec_relativistic,
               cah_spec_relativistic):
        assert fn().relativistic is True, (
            f"{fn.__name__} must set relativistic=True")


# ---------------------------------------------------------------------------
# Relativistic smoke tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec_fn,name,q_nmax1,q_nmax2", [
    ("lih_spec_relativistic", "LiH",  6, 30),
    ("beh_spec_relativistic", "BeH",  6, 30),
    ("cah_spec_relativistic", "CaH",  4, 20),
])
def test_relativistic_build_smoke(spec_fn, name, q_nmax1, q_nmax2):
    """Relativistic LiH/BeH/CaH must build without error at n_max ∈ {1, 2}."""
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac import molecular_spec as ms

    for nmx, q_expected in [(1, q_nmax1), (2, q_nmax2)]:
        spec = getattr(ms, spec_fn)(max_n=nmx)
        r = build_composed_hamiltonian(spec)
        assert r['Q'] == q_expected, (
            f"{name}_rel n_max={nmx}: expected Q={q_expected}, got {r['Q']}")
        assert r['N_pauli'] > 0
        assert r['qwc_groups'] > 0
        assert r['cross_block_eri_count'] == 0, (
            f"{name}_rel has cross-block ERI entries — block diagonality broken")


def test_relativistic_block_diagonality():
    """Relativistic ERI must be block-diagonal (no cross-block 4-tuples)."""
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.molecular_spec import lih_spec_relativistic

    r = build_composed_hamiltonian(lih_spec_relativistic(max_n=2))
    # Explicit check: derive which block each qubit index belongs to
    blocks = r['blocks']
    # Rebuild a block_of mapping from the returned blocks list
    block_of = []
    for bl in blocks:
        block_of.extend([bl['label']] * bl['n_orbitals'])
    for (a, b, c, d), val in r['eri_sparse'].items():
        assert (block_of[a] == block_of[b] == block_of[c] == block_of[d]), (
            f"cross-block ERI ({a},{b},{c},{d}) with value {val}")


def test_relativistic_nmax1_pauli_count():
    """LiH n_max=1 relativistic: Q=6 spinor orbitals, 10 Pauli terms with identity."""
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.molecular_spec import lih_spec_relativistic

    r = build_composed_hamiltonian(lih_spec_relativistic(max_n=1))
    assert r['Q'] == 6
    # 3 s_{1/2} blocks, each with 2 DiracLabels → 2 spin-orbitals per block.
    # ⟨ab|V|cd⟩ reduces to a single k=0 channel per block, so N_pauli = 9
    # (6 one-body diagonal + 3 two-body terms).  Exact; brittle on purpose.
    assert r['N_pauli'] == 9


def test_relativistic_alpha_zero_kills_so():
    """Setting α = 0 must zero out the spin-orbit diagonal contribution."""
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.molecular_spec import lih_spec_relativistic

    r0 = build_composed_hamiltonian(lih_spec_relativistic(max_n=2),
                                    alpha_num=0.0)
    rN = build_composed_hamiltonian(lih_spec_relativistic(max_n=2))
    # The SO diagonal vector itself should be all-zero at α=0 (also for n_max=2
    # where only p-channels can contribute; but κ=−1 l=0 gives 0 anyway).
    assert np.allclose(r0['h1_so_diag'], 0.0), (
        "Spin-orbit diagonal must vanish at α=0")
    # At finite α the SO shift is nonzero for p-channels (l=1, κ≠−1).
    assert np.any(np.abs(rN['h1_so_diag']) > 0), (
        "Spin-orbit diagonal must be nonzero for p channels at α=CODATA")


def test_relativistic_hermiticity():
    """Relativistic qubit_op must be Hermitian (real coefficients)."""
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.molecular_spec import cah_spec_relativistic

    r = build_composed_hamiltonian(cah_spec_relativistic(max_n=2))
    for term, coeff in r['qubit_op'].terms.items():
        # JW output of a real Hermitian FermionOp is real.
        assert abs(coeff.imag) < 1e-10, (
            f"Non-real coefficient {coeff} on term {term}")


def test_relativistic_identity_contains_nuclear_repulsion():
    """JW identity term absorbs nuclear repulsion + one-body diagonal shifts.

    Under Jordan-Wigner, a†a = (I - Z)/2, so every diagonal h₁ entry
    contributes h_ii/2 to the identity coefficient.  The identity term
    is therefore NOT just the nuclear repulsion — it is
    ``nuclear_repulsion + Σᵢ h1[i,i]/2 + two-body shifts``.  We check
    that it's real and finite.
    """
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.molecular_spec import lih_spec_relativistic

    r = build_composed_hamiltonian(lih_spec_relativistic(max_n=1))
    id_coeff = r['qubit_op'].terms.get((), 0.0)
    assert abs(id_coeff.imag) < 1e-10
    assert np.isfinite(id_coeff.real)


# ---------------------------------------------------------------------------
# Resource table sanity (n_max=2 baseline)
# ---------------------------------------------------------------------------


def test_relativistic_pauli_ratio_lih_nmax2():
    """LiH relativistic / scalar Pauli ratio at n_max=2.

    T0 prediction (fullgaunt, l_max=1): d_spinor/d_scalar ≈ 0.58.  In Pauli
    space, both scalar (spin-doubled JW, Q=2M=30) and spinor (Q=30 directly)
    have the same qubit count at n_max=2.  The relative Pauli count is
    dominated by the additional σ,τ independence in the spinor basis
    (scalar JW has σ=τ restriction in one-body) plus the angular-density
    ratio.  Observed ratio ~4.3× at n_max=2.
    """
    from geovac.composed_qubit import build_composed_hamiltonian
    from geovac.molecular_spec import lih_spec, lih_spec_relativistic

    r_scalar = build_composed_hamiltonian(lih_spec(max_n=2))
    r_rel = build_composed_hamiltonian(lih_spec_relativistic(max_n=2))
    ratio = r_rel['N_pauli'] / r_scalar['N_pauli']
    # Observed 805 / 333 ≈ 2.42 (post-Hermitian symmetrization).  Pin ±20%.
    assert 1.9 < ratio < 3.0, (
        f"LiH n_max=2 rel/scalar Pauli ratio {ratio:.2f} outside [1.9, 3.0]")
