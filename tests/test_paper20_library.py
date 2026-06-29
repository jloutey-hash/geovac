"""Paper 20 / Paper 14 — molecule-library size regression pin.

The papers + CLAUDE.md headline the `hamiltonian()` library as **37 systems**
(35 composed molecules + the He and H2 reference atoms). A QA cert (2026-06-29)
found this count was asserted in prose but had NO test — `len(set(...))` was never
checked, so a registry add/drop could silently desync the headline. This pins it.
"""
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
