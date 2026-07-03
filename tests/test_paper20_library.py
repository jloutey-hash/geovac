"""Paper 20 / Paper 14 — molecule-library size regression pin.

The papers + CLAUDE.md headline the `hamiltonian()` library as **37 systems**
(35 composed molecules + the He and H2 reference atoms). A QA cert (2026-06-29)
found this count was asserted in prose but had NO test — `len(set(...))` was never
checked, so a registry add/drop could silently desync the headline. This pins it.
"""
import pytest

from geovac.ecosystem_export import _SYSTEM_REGISTRY


def test_library_has_37_distinct_systems() -> None:
    """The library headline (Papers 14/20, CLAUDE.md §1.5) is 37 distinct systems."""
    distinct = set(_SYSTEM_REGISTRY.values())
    assert len(distinct) == 37, (
        f"library has {len(distinct)} distinct systems, expected 37 "
        f"(35 composed molecules + He + H2). Registry desynced from the paper headline."
    )


def test_library_includes_he_and_h2_reference_atoms() -> None:
    """The '35 composed + He + H2' decomposition: the two reference atoms must be present."""
    canon = set(_SYSTEM_REGISTRY.values())
    for ref in ("He", "H2"):
        assert ref in canon, f"reference system {ref!r} missing from the library registry"


def test_composed_lih_one_norm_headline_pinned() -> None:
    """P20 abstract + tab:resources: composed LiH electronic-only 1-norm
    32.6 Ha (8th-cert backfill -- previously only one_norm > 0 was asserted)."""
    from geovac.ecosystem_export import hamiltonian
    H = hamiltonian('LiH')
    assert abs(H.one_norm - 32.6) < 0.1, H.one_norm


@pytest.mark.slow
def test_beh2_h2o_qpe_regime_one_norms_pinned() -> None:
    """P20 SBeH2-H2O + P14 tab:balanced cells (8th-cert backfill): the
    lambda convention here INCLUDES the identity coefficient (the QPE
    encoding cost); the BeH2 composed-with-PK 354.9 vintage came from the
    deprecated legacy builder path (live 373.4)."""
    from geovac.molecular_spec import beh2_spec, h2o_spec
    from geovac.balanced_coupled import build_balanced_hamiltonian
    from geovac.composed_qubit import build_composed_hamiltonian

    def lam_all(op):
        return sum(abs(c) for c in op.terms.values())

    bal_beh2 = build_balanced_hamiltonian(beh2_spec(), R=2.539)['qubit_op']
    assert len(bal_beh2.terms) == 2652
    assert abs(lam_all(bal_beh2) - 306.4) < 0.5, lam_all(bal_beh2)

    bal_h2o = build_balanced_hamiltonian(h2o_spec(), R=1.808)['qubit_op']
    assert len(bal_h2o.terms) == 5798
    assert abs(lam_all(bal_h2o) - 1511.1) < 1.5, lam_all(bal_h2o)

    comp_beh2 = build_composed_hamiltonian(beh2_spec())['qubit_op']
    assert abs(lam_all(comp_beh2) - 373.4) < 0.5, lam_all(comp_beh2)
    # the QPE-regime direction: balanced strictly cheaper than composed+PK
    assert lam_all(bal_beh2) < lam_all(comp_beh2)

    comp_h2o = build_composed_hamiltonian(h2o_spec())['qubit_op']
    assert abs(lam_all(comp_h2o) - 28055.2) < 5.0, lam_all(comp_h2o)
