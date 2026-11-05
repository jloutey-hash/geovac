"""Paper 23 (nuclear shell) resource-count backing tests.

Closes the NO-TEST gaps flagged in docs/claim_test_matrix.md (group4 pre-work,
2026-06-28): the deuteron / He-4 Pauli term counts and 1-norms were uncaught by any
test (loose 0<n<10000 guards only), which let the paper's 1-norm numbers drift from
the live code (deuteron 227->342 MeV; He-4 557/552->467/462 MeV) without detection.

These pin the reproducible code values at the paper's hw=10 MeV, N_shells=2.

The deuteron structure is 688 non-I (80 Z-only + 608 XY) and He-4 is 904; the
1-norms are 383.9 / 641.4 / 632.2 MeV (Minnesota V_S = -91.85 MeV).  Counts alone
cannot catch a wrong encoding (a real-coefficient XY/YX encoding had the same 688
terms and the wrong spectrum), so every Pauli Hamiltonian here is also checked
against its FCI matrix in the physical particle-number sector.
"""
from __future__ import annotations
from itertools import combinations

import numpy as np
import pytest

from geovac.nuclear.nuclear_hamiltonian import (
    build_deuteron_hamiltonian,
    build_he4_hamiltonian,
    pauli_sector_matrix,
)


def _sector_basis(n_p, n_n, N_p, N_n):
    """Occupation tuples with N_p protons and N_n neutrons (proton qubits first)."""
    out = []
    for pp in combinations(range(n_p), N_p):
        for nn in combinations(range(n_n), N_n):
            occ = [0] * (n_p + n_n)
            for a in pp:
                occ[a] = 1
            for b in nn:
                occ[n_p + b] = 1
            out.append(tuple(occ))
    return out


def _spectrum_error(data, N_p, N_n):
    """Max |eigenvalue difference| between the Pauli sum (sector) and the FCI matrix."""
    M = pauli_sector_matrix(
        data['H_pauli'], _sector_basis(data['Q_p'], data['Q_n'], N_p, N_n))
    assert np.abs(M - M.conj().T).max() < 1e-10, "sector matrix not Hermitian"
    H = data['H_matrix']
    H = H.toarray() if hasattr(H, 'toarray') else np.asarray(H)
    ev_pauli = np.linalg.eigvalsh(M)
    ev_fci = np.linalg.eigvalsh(H)
    assert ev_pauli.shape == ev_fci.shape
    return float(np.abs(ev_pauli - ev_fci).max())


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
    """Deuteron at hw=10: 688 non-I Pauli (80 Z-only + 608 XY), 1-norm ~383.9 MeV."""
    H = build_deuteron_hamiltonian(N_shells=2, hw=10.0)['H_pauli']
    n_nonid, z_only, xy, l1 = _split(H)
    assert n_nonid == 688, f"deuteron non-I Pauli: expected 688, got {n_nonid}"
    assert z_only == 80, f"deuteron Z-only: expected 80, got {z_only}"
    assert xy == 608, f"deuteron XY: expected 608, got {xy}"
    # Paper 23 Table I: 1-norm (non-I) ~= 383.9 MeV
    assert abs(l1 - 383.87) < 0.05, f"deuteron 1-norm: expected ~383.87 MeV, got {l1:.2f}"


def test_he4_resource_counts_no_coulomb():
    """He-4 at hw=10, no Coulomb: 904 non-I Pauli, 1-norm ~641.4 MeV."""
    H = build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=False)['H_pauli']
    n_nonid, _, _, l1 = _split(H)
    assert n_nonid == 904, f"He-4 non-I Pauli: expected 904, got {n_nonid}"
    assert abs(l1 - 641.39) < 0.05, f"He-4 (no Coul) 1-norm: expected ~641.39 MeV, got {l1:.2f}"


def test_he4_resource_counts_with_coulomb():
    """He-4 at hw=10, with Coulomb: 904 non-I Pauli, 1-norm ~632.2 MeV."""
    H = build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=True)['H_pauli']
    n_nonid, _, _, l1 = _split(H)
    assert n_nonid == 904, f"He-4 non-I Pauli: expected 904, got {n_nonid}"
    assert abs(l1 - 632.19) < 0.05, f"He-4 (Coul) 1-norm: expected ~632.19 MeV, got {l1:.2f}"


def test_deuteron_pauli_spectrum_matches_fci():
    """The encoded deuteron reproduces the 64-dim (1p, 1n) FCI spectrum."""
    err = _spectrum_error(build_deuteron_hamiltonian(N_shells=2, hw=10.0), 1, 1)
    assert err < 1e-9, f"deuteron Pauli vs FCI spectrum: max error {err:.3e} MeV"


@pytest.mark.parametrize("coulomb", [False, True])
def test_he4_pauli_spectrum_matches_fci(coulomb):
    """The encoded He-4 reproduces the 784-dim (2p, 2n) FCI spectrum."""
    data = build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=coulomb)
    err = _spectrum_error(data, 2, 2)
    assert err < 1e-9, f"He-4 Pauli vs FCI spectrum (Coulomb={coulomb}): max error {err:.3e} MeV"


def test_he4_coulomb_reduces_1norm():
    """Coulomb repulsion lowers the He-4 1-norm by a few MeV (structural check)."""
    l1_no = _split(build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=False)['H_pauli'])[3]
    l1_yes = _split(build_he4_hamiltonian(N_shells=2, hw=10.0, include_coulomb=True)['H_pauli'])[3]
    assert l1_yes < l1_no, "Coulomb should reduce the 1-norm"
    assert 0 < (l1_no - l1_yes) < 20, f"unexpected Coulomb 1-norm shift {l1_no - l1_yes:.2f}"


def test_paper23_composed_nuclear_electronic_counts():
    """Composed nuclear-electronic deuterium: Q=26, 710 non-I Pauli.

    Paper 23 sec "Resource counts and coefficient hierarchy" states the
    full-register count decomposes as nuclear + electronic + cross. This
    pins the decomposition AND its sum, so a future edit cannot leave a
    component larger than the total (which is exactly what happened when
    the 2026-08-22 moshinsky correction was propagated: the nuclear block
    went 592 -> 688 while the stated total stayed at 614).
    """
    from geovac.nuclear.nuclear_electronic import (
        build_deuterium_composed_hamiltonian,
        analyze_composed_hamiltonian,
    )
    res = build_deuterium_composed_hamiltonian()
    info = analyze_composed_hamiltonian(
        res['pauli_terms'], res['Q_nuc'], res['Q_elec'])

    assert info['Q_total'] == 26, f"Q_total = {info['Q_total']}, expected 26"
    assert info['Q_nuc'] == 16 and info['Q_elec'] == 10

    nuc = info['n_nuclear_only']
    ele = info['n_electronic_only']
    cross = info['n_cross_register']
    total = info['n_pauli_terms_non_identity']

    assert nuc == 688, f'nuclear block: expected 688, got {nuc}'
    assert ele == 10, f'electronic block: expected 10, got {ele}'
    assert cross == 12, f'cross-register hyperfine: expected 12, got {cross}'
    assert total == 710, f'composed non-I total: expected 710, got {total}'

    # The decomposition must actually close -- the invariant the stale
    # 614 violated.
    assert nuc + ele + cross == total, \
        f'decomposition does not close: {nuc} + {ele} + {cross} != {total}'
