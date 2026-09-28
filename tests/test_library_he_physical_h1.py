"""Guard: the pip-installable library's He qubit Hamiltonian is helium's
PHYSICAL Hamiltonian, not the retired hop-augmented model.

Until 2026-09-25 ``ecosystem_export._build_he`` built He with
``LatticeIndex(h1_method='hybrid')``: the exact hydrogenic diagonal
-Z^2/(2n^2) PLUS a kappa*(-A) = +1/16 one-body hop on every GeometricLattice
edge.  In the fixed k = Z hydrogenic basis T + V_ne is exactly diagonal, so
the hop is not part of the atom's Hamiltonian (it lowers He by ~48 mHa, more
than helium's correlation energy).  The library now uses h1_method='exact'.

What wrong answer each test rejects
-----------------------------------
* ``test_one_electron_block_is_hydrogenic`` reads the one-body operator back
  OUT OF THE QUBIT OPERATOR, independently of how it was built: in the N = 1
  sector every normal-ordered two-body term vanishes, so the N = 1 block of a
  JW-encoded Hamiltonian IS the spin-orbital h1.  It must be diagonal with
  spectrum {-Z^2/(2n^2)} of multiplicity 2 n^2 per shell.  The hop model
  gives off-diagonal entries of exactly 1/16 -> rejected.  A wrong Z, or a
  graph-Laplacian 'graph' h1, changes the spectrum -> rejected.
* ``test_library_he_resource_counts`` pins the physical Pauli count and
  1-norm (MEASURED 2026-09-25, debug/library_he_fix.py; the composed builder
  on the same He spec is an independent second route and agrees on the
  count exactly and on lambda to 5e-5 relative).  The hop model gives
  287 / 11.175 and 14078 / 74.207 -> rejected.
* ``test_vqe_benchmark_builder_matches_library`` rejects the two He builders
  drifting apart again (the benchmark builder carried the same hop).
* ``test_legacy_opt_in_is_the_hop_model`` shows the opt-in still reproduces
  the retired numbers, i.e. that the tests above discriminate a real
  difference rather than passing on any He operator.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from geovac.ecosystem_export import hamiltonian, _build_he

Z = 2

# (n_max) -> (Q, non-identity Pauli, non-identity 1-norm), physical h1.
PHYSICAL = {2: (10, 279, 10.800482), 3: (28, 14046, 72.582417)}
# Retired hop-augmented values (registry aliases / C21 RETIRED).
HOP_MODEL = {2: (287, 11.175482), 3: (14078, 74.207417)}


def _one_electron_block(qop, n_qubits: int) -> np.ndarray:
    """<e_p| P |e_q> summed over Pauli strings, e_q = single qubit q set.

    openfermion conventions: X|b> = |1-b>, Y|0> = i|1>, Y|1> = -i|0>,
    Z|b> = (-1)^b |b>.
    """
    H = np.zeros((n_qubits, n_qubits), dtype=complex)
    for term, coeff in qop.terms.items():
        flips = [k for k, p in term if p in ("X", "Y")]
        if len(flips) > 2:
            continue
        for q in range(n_qubits):
            bits = {q: 1}
            phase = 1.0 + 0j
            for k, p in term:
                b = bits.get(k, 0)
                if p == "Z":
                    phase *= -1 if b else 1
                elif p == "X":
                    bits[k] = 1 - b
                elif p == "Y":
                    phase *= 1j if b == 0 else -1j
                    bits[k] = 1 - b
            occ = [k for k, v in bits.items() if v]
            if len(occ) != 1:
                continue
            H[occ[0], q] += coeff * phase
    return H


def _expected_spectrum(n_max: int) -> np.ndarray:
    vals = []
    for n in range(1, n_max + 1):
        vals += [-Z * Z / (2.0 * n * n)] * (2 * n * n)  # 2 spins x n^2 states
    return np.sort(np.array(vals))


@pytest.mark.parametrize("n_max", [2, 3])
def test_one_electron_block_is_hydrogenic(n_max: int) -> None:
    ham = hamiltonian("He", max_n=n_max)
    Q = PHYSICAL[n_max][0]
    assert ham.n_qubits == Q
    h = _one_electron_block(ham._qubit_op, Q)
    assert np.abs(h.imag).max() < 1e-12
    h = h.real
    off = h - np.diag(np.diag(h))
    assert np.abs(off).max() < 1e-12, (
        f"library He one-body operator has off-diagonal entries "
        f"(max {np.abs(off).max():.6g}; the retired graph hop is 1/16)"
    )
    np.testing.assert_allclose(np.sort(np.diag(h)), _expected_spectrum(n_max),
                               atol=1e-12)


@pytest.mark.parametrize("n_max", [2, 3])
def test_library_he_resource_counts(n_max: int) -> None:
    Q, n_ref, lam_ref = PHYSICAL[n_max]
    ham = hamiltonian("He", max_n=n_max)
    terms = ham._qubit_op.terms
    n_pauli = sum(1 for k in terms if k)
    lam = float(sum(abs(c) for k, c in terms.items() if k))
    assert ham.metadata["h1_method"] == "exact"
    assert n_pauli == n_ref, f"n_max={n_max}: {n_pauli} != {n_ref}"
    assert lam == pytest.approx(lam_ref, abs=5e-6)


def test_vqe_benchmark_builder_matches_library() -> None:
    from geovac.vqe_benchmark import build_geovac_he
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, of_op, nq, _ = build_geovac_he(max_n=2)
    lib = hamiltonian("He", max_n=2)._qubit_op
    assert nq == 10
    assert set(of_op.terms) == set(lib.terms)
    for k, c in lib.terms.items():
        assert abs(of_op.terms[k] - c) < 1e-12


@pytest.mark.parametrize("n_max", [2])
def test_legacy_opt_in_is_the_hop_model(n_max: int) -> None:
    legacy = _build_he(max_n=n_max, include_graph_hop=True)
    terms = legacy._qubit_op.terms
    n_ref, lam_ref = HOP_MODEL[n_max]
    assert sum(1 for k in terms if k) == n_ref
    assert float(sum(abs(c) for k, c in terms.items() if k)) == \
        pytest.approx(lam_ref, abs=5e-6)
    h = _one_electron_block(legacy._qubit_op, PHYSICAL[n_max][0]).real
    off = np.abs(h - np.diag(np.diag(h)))
    assert np.isclose(off.max(), 1.0 / 16.0, atol=1e-12)
