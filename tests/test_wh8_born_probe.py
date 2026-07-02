"""WH8 Step-1 frozen falsifier -- the skeleton measure does NOT generate Born.

Pins the 2026-07-01 probe (debug/wh8_born_measure_probe.py, PI-directed):
the graph's only native measures are counting-type (rational, built from
discrete invariants), and they reproduce Born statistics exactly on the FLAT
locus and nowhere physical:

  A. TV(Born(flat), counting) = 0 exactly; TV = 1/3 exactly for amplitudes
     proportional to (1,2,3,4) on a 4-fold multiplet (rational identity).
  B. He graph-native CI ground state (max_n=3, slater_full, 0 params):
     ~91% weight on 1s^2; TV vs uniform-on-support ~ 0.86.
  C. The flat locus is dynamically unstable under the graph's own H
     (TV ~ t^2, exponent ~2).

Upgrades the Born rule's Class-1 (calibration-data) classification to a
TESTED NEGATIVE. WH8 (CLAUDE.md 1.7): the Born measure is the exchange
constant of the observation projection, unique given the projection lattice
(Gleason, dim >= 3). FALSIFIER DIRECTION: a skeleton-side rule
(state-independent, built from degeneracies / selection rules / rational
couplings) reproducing Born statistics on non-flat physical states would
break these pins and falsify the import classification.

Self-contained recompute (no debug/ import -- the durability rule).
"""
from __future__ import annotations

import warnings
from fractions import Fraction

import numpy as np
import pytest
from scipy.sparse.linalg import eigsh


def _tv(p: np.ndarray, q: np.ndarray) -> float:
    return 0.5 * float(np.sum(np.abs(p - q)))


def _born(vec: np.ndarray) -> np.ndarray:
    w = np.abs(vec) ** 2
    return w / w.sum()


def test_flat_locus_coincidence_is_exact_and_rational() -> None:
    """Leg A: counting = Born exactly on flat states; the generic deviation
    is the exact rational 1/3 for amplitudes (1,2,3,4) on d=4."""
    d = 4
    flat = np.ones(d) / np.sqrt(d)
    counting = np.ones(d) / d
    assert _tv(_born(flat), counting) == 0.0

    # exact rational arithmetic: born weights (1,4,9,16)/30 vs 1/4 each
    born_q = [Fraction(k * k, 30) for k in (1, 2, 3, 4)]
    tv_q = Fraction(1, 2) * sum(abs(b - Fraction(1, 4)) for b in born_q)
    assert tv_q == Fraction(1, 3)
    generic = np.array([1.0, 2.0, 3.0, 4.0])
    generic /= np.linalg.norm(generic)
    assert abs(_tv(_born(generic), counting) - 1.0 / 3.0) < 1e-12


def _he_ground_state():
    from geovac.lattice_index import LatticeIndex
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        li = LatticeIndex(n_electrons=2, max_n=3, nuclear_charge=2,
                          vee_method='slater_full')
        H = li.assemble_hamiltonian()
    evals, evecs = eigsh(H, k=1, which='SA')
    return li, H, float(evals[0]), evecs[:, 0]


def test_he_ground_state_born_diverges_from_all_counting_measures() -> None:
    """Leg B: the physical ground state is far from every skeleton
    counting-measure candidate (pinned 2026-07-01: TV_support ~ 0.858,
    top weight 0.908 on 1s^2)."""
    li, H, E0, psi = _he_ground_state()
    p = _born(psi)
    n_sd = len(p)
    assert n_sd == 378

    q_basis = np.ones(n_sd) / n_sd
    supp = p > 1e-14
    q_support = np.where(supp, 1.0, 0.0)
    q_support /= q_support.sum()

    tv_support = _tv(p, q_support)
    assert 0.80 < tv_support < 0.90, tv_support           # pinned 0.858
    assert _tv(p, q_basis) > 0.95                          # pinned 0.986

    # dominant configuration is 1s up, 1s down with ~91% weight
    top = int(np.argmax(p))
    labels = sorted(li.sp_states[o][:3] for o in li.sd_basis[top])
    assert labels == [(1, 0, 0), (1, 0, 0)], labels
    assert 0.85 < float(p[top]) < 0.95, float(p[top])      # pinned 0.908

    # sanity: the energy is the known max_n=3 value (graph-native CI series)
    assert abs(E0 - (-2.8404)) < 5e-3, E0


def test_flat_locus_dynamically_unstable_t_squared() -> None:
    """Leg C: a flat state on the GS support leaves the coincidence locus
    under the graph's own Hamiltonian at the perturbative t^2 rate."""
    _, H, _, psi = _he_ground_state()
    p = _born(psi)
    supp = np.where(p > 1e-14)[0]
    n_sd = len(p)
    psi0 = np.zeros(n_sd)
    psi0[supp] = 1.0 / np.sqrt(len(supp))
    counting = np.zeros(n_sd)
    counting[supp] = 1.0 / len(supp)

    E, V = np.linalg.eigh(H.toarray())
    c0 = V.T @ psi0

    def tv_at(t: float) -> float:
        psi_t = V @ (np.exp(-1j * E * t) * c0)
        return _tv(_born(psi_t), counting)

    assert tv_at(0.0) < 1e-12   # float dust from the V @ V.T round-trip
    tv1, tv2, tv5 = tv_at(0.01), tv_at(0.02), tv_at(0.5)
    assert tv5 > 3e-3                                       # pinned 6.1e-3
    slope = np.log(tv2 / tv1) / np.log(2.0)
    assert abs(slope - 2.0) < 0.1, slope                    # pinned 2.00
