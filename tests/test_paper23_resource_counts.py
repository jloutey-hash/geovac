"""Paper 23 (nuclear shell) resource-count backing tests.

Closes the NO-TEST gaps flagged in docs/claim_test_matrix.md (group4 pre-work,
2026-06-28): the deuteron / He-4 Pauli term counts and 1-norms were uncaught by any
test (loose 0<n<10000 guards only), which let the paper's 1-norm numbers drift from
the live code (deuteron 227->342 MeV; He-4 557/552->467/462 MeV) without detection.

These pin the reproducible code values that Paper 23 was corrected to (2026-06-28),
at the paper's hw=10 MeV, N_shells=2. The term-count structure (592 / 80 Z-only /
512 XY for the deuteron; 712 for He-4) was already correct in the paper; only the
1-norm magnitudes had drifted.
"""
from __future__ import annotations
import pytest

from geovac.nuclear.nuclear_hamiltonian import (
    build_deuteron_hamiltonian,
    build_he4_hamiltonian,
)


def _split(H_pauli):
    """(non-identity count, Z-only count, XY count, 1-norm non-identity)."""
    n_nonid = z_only = xy = 0
    l1 = 0.0
    for term, coeff in H_pauli.items():
        s = str(term)
        letters = set(s) - {'I'}
        if not letters:          # identity
            continue
        n_nonid += 1
        l1 += abs(coeff)
        if letters <= {'Z'}:
            z_only += 1
        elif letters & {'X', 'Y'}:
            xy += 1
    return n_nonid, z_only, xy, float(l1)


def test_deuteron_resource_counts():
    """Deuteron at hw=10: 592 non-I Pauli (80 Z-only + 512 XY), 1-norm ~342 MeV."""
    H = build_deuteron_hamiltonian(N_shells=2, hw=10.0)['H_pauli']
    n_nonid, z_only, xy, l1 = _split(H)
    assert n_nonid == 592, f"deuteron non-I Pauli: expected 592, got {n_nonid}"
    assert z_only == 80, f"deuteron Z-only: expected 80, got {z_only}"
    assert xy == 512, f"deuteron XY: expected 512, got {xy}"
    # Paper 23 Table I (corrected 2026-06-28): 1-norm (non-I) ~= 342.2 MeV
    assert abs(l1 - 342.2) < 1.0, f"deuteron 1-norm: expected ~342.2 MeV, got {l1:.2f}"


def test_he4_resource_counts_no_coulomb():
    """He-4 at hw=10, no Coulomb: 712 non-I Pauli, 1-norm ~467 MeV."""
    H = build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=False)['H_pauli']
    n_nonid, _, _, l1 = _split(H)
    assert n_nonid == 712, f"He-4 non-I Pauli: expected 712, got {n_nonid}"
    assert abs(l1 - 466.9) < 1.0, f"He-4 (no Coul) 1-norm: expected ~466.9 MeV, got {l1:.2f}"


def test_he4_resource_counts_with_coulomb():
    """He-4 at hw=10, with Coulomb: 712 non-I Pauli, 1-norm ~462 MeV."""
    H = build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=True)['H_pauli']
    n_nonid, _, _, l1 = _split(H)
    assert n_nonid == 712, f"He-4 non-I Pauli: expected 712, got {n_nonid}"
    assert abs(l1 - 462.4) < 1.0, f"He-4 (Coul) 1-norm: expected ~462.4 MeV, got {l1:.2f}"


def test_he4_coulomb_reduces_1norm():
    """Coulomb repulsion lowers the He-4 1-norm by a few MeV (structural check)."""
    l1_no = _split(build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=False)['H_pauli'])[3]
    l1_yes = _split(build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=True)['H_pauli'])[3]
    assert l1_yes < l1_no, "Coulomb should reduce the 1-norm"
    assert 0 < (l1_no - l1_yes) < 20, f"unexpected Coulomb 1-norm shift {l1_no - l1_yes:.2f}"
