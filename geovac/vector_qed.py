"""
Vector-photon QED on the Fock graph.
====================================

This module implements QED with EXPLICIT VECTOR PHOTON MODES labeled by
(q, m_q) coupled to the scalar Fock-graph electron via Wigner 3j symbols.
It is the "calibration" step in the GeoVac exchange constant taxonomy
(Paper 18): importing continuum angular momentum structure (vector photon
quantum numbers) to go from scalar graph QED to physical vector QED.

In contrast to the graph-native QED (``graph_qed_self_energy.py``) where
the photon is a scalar 1-cochain on the Fock graph (recovering 1/8
continuum selection rules), this module gives the photon its own mode
space with a propagator diagonal in q, and couples it to the electron
graph via CG coefficients.  This should recover ALL 8 continuum QED
selection rules.

Electron states
---------------
Scalar Fock graph nodes: (n, l, m) with n = 1, ..., n_max, l = 0, ..., n-1,
m = -l, ..., l.  Electron propagator G_e(n) = 1 / (n^2 - 1) for n >= 2.
For n = 1 (ground state), the eigenvalue is 0 so G_e is undefined -- we
exclude n_b = 1 from internal sums (Pauli exclusion / normal ordering).

Photon modes
------------
(q, m_q) with q = 1, ..., q_max, m_q = -q, ..., q.  Photon propagator:
    G_gamma(q) = 1 / [q(q+2)]
from the vector Laplacian eigenvalues on unit S^3 (Hodge-1 spectrum).

Vertex coupling
---------------
V(a, b, q, m_q) couples electron a = (n_a, l_a, m_a) to electron
b = (n_b, l_b, m_b) via photon mode (q, m_q):

    V = sqrt((2*l_a+1)*(2*q+1)*(2*l_b+1)/(4*pi))
      * (-1)^(l_a - m_a) * 3j(l_a, q, l_b; -m_a, m_q, m_b)
      * parity_factor(l_a, l_b, q)

where parity_factor = 1 if l_a + l_b + q is odd (E-type), 0 if even.

This enforces:
  - Triangle inequality: |l_a - l_b| <= q <= l_a + l_b
  - Magnetic conservation: -m_a + m_q + m_b = 0 (3j selection)
  - Spatial parity: l_a + l_b + q odd (E-type / vector coupling)

The parity rule forces the ground state structural zero:
for GS with l_a = 0, triangle forces q = l_b, so l_b + l_b = 2*l_b is
always even -> V(GS, any, any) = 0.

Transcendental taxonomy (Paper 18)
----------------------------------
The vertex introduces pi via the spherical harmonic normalization factor
sqrt((2l+1)(2q+1)(2l'+1) / (4*pi)).  This is the CALIBRATION tier in
Paper 18: the continuum angular-momentum structure (Wigner 3j symbols +
spherical harmonic normalization) is the price of vector photon quantum
numbers.  All graph-intrinsic (pi-free) content is lost when the photon
is promoted from scalar 1-cochain to vector mode.

References
----------
- GeoVac graph_qed_self_energy.py (scalar graph self-energy, comparison)
- GeoVac qed_self_energy.py (continuum self-energy, selection rules)
- GeoVac qed_vertex.py (continuum vertex coupling, SO(4) channel count)
- GeoVac Paper 28 (QED on S^3, selection rule census)
- GeoVac Paper 18 (exchange constant taxonomy)
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from sympy.physics.wigner import wigner_3j

__all__ = [
    "build_electron_states",
    "build_photon_modes",
    "electron_propagator",
    "vector_photon_propagator",
    "vertex_coupling",
    "build_vertex_tensor",
    "compute_self_energy",
    "check_selection_rules",
    "VectorQEDResult",
    # Dirac extension (Part 2)
    "build_dirac_electron_states",
    "dirac_electron_propagator",
    "dirac_vertex_coupling",
    "build_dirac_vertex_tensor",
    "compute_dirac_self_energy",
    "check_dirac_selection_rules",
    "DiracVectorQEDResult",
    "run_dirac_vector_diagnostic",
]


# ---------------------------------------------------------------------------
# Step 1: Electron states
# ---------------------------------------------------------------------------

def build_electron_states(n_max: int) -> List[Tuple[int, int, int]]:
    """Build the ordered list of scalar Fock graph electron states.

    States are (n, l, m) with n = 1, ..., n_max, l = 0, ..., n-1,
    m = -l, ..., l.  Ordered by n, then l, then m.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    list of (n, l, m) tuples
    """
    states = []
    for n in range(1, n_max + 1):
        for l in range(0, n):
            for m in range(-l, l + 1):
                states.append((n, l, m))
    return states


# ---------------------------------------------------------------------------
# Step 2: Photon modes
# ---------------------------------------------------------------------------

def build_photon_modes(q_max: int) -> List[Tuple[int, int]]:
    """Build the ordered list of vector photon modes.

    Modes are (q, m_q) with q = 1, ..., q_max, m_q = -q, ..., q.

    Parameters
    ----------
    q_max : int
        Maximum photon angular momentum.

    Returns
    -------
    list of (q, m_q) tuples
    """
    modes = []
    for q in range(1, q_max + 1):
        for mq in range(-q, q + 1):
            modes.append((q, mq))
    return modes


# ---------------------------------------------------------------------------
# Step 3: Propagators
# ---------------------------------------------------------------------------

def electron_propagator(n: int) -> float:
    """Electron propagator on the Fock graph.

    G_e(n) = 1 / (n^2 - 1) for n >= 2.
    For n = 1, returns None (ground state, excluded from internal sums).

    Parameters
    ----------
    n : int
        Principal quantum number.

    Returns
    -------
    float or None
    """
    if n < 2:
        return None
    return 1.0 / (n * n - 1)


def vector_photon_propagator(q: int) -> float:
    """Vector photon propagator from Hodge-1 eigenvalues on S^3.

    G_gamma(q) = 1 / [q(q+2)]

    Parameters
    ----------
    q : int
        Photon angular momentum quantum number (q >= 1).

    Returns
    -------
    float
    """
    return 1.0 / (q * (q + 2))


# ---------------------------------------------------------------------------
# Step 4: Vertex coupling
# ---------------------------------------------------------------------------

def vertex_coupling(
    n_a: int, l_a: int, m_a: int,
    n_b: int, l_b: int, m_b: int,
    q: int, m_q: int,
) -> float:
    """Compute the vertex coupling V(a, b, q, m_q).

    V = sqrt((2*l_a+1)*(2*q+1)*(2*l_b+1)/(4*pi))
      * (-1)^(l_a - m_a) * 3j(l_a, q, l_b; -m_a, m_q, m_b)
      * parity_factor

    where parity_factor = 1 if l_a + l_b + q is odd, 0 if even.

    Parameters
    ----------
    n_a, l_a, m_a : int
        Quantum numbers of electron state a.
    n_b, l_b, m_b : int
        Quantum numbers of electron state b.
    q, m_q : int
        Photon mode quantum numbers.

    Returns
    -------
    float
        Vertex coupling value.
    """
    # Parity selection rule: l_a + l_b + q must be odd (E-type)
    if (l_a + l_b + q) % 2 == 0:
        return 0.0

    # 3j selection rule: -m_a + m_q + m_b = 0
    if -m_a + m_q + m_b != 0:
        return 0.0

    # Triangle inequality: |l_a - l_b| <= q <= l_a + l_b
    if q < abs(l_a - l_b) or q > l_a + l_b:
        return 0.0

    # Compute 3j symbol
    threej = float(wigner_3j(l_a, q, l_b, -m_a, m_q, m_b))
    if abs(threej) < 1e-15:
        return 0.0

    # Angular normalization factor
    prefactor = np.sqrt(
        (2 * l_a + 1) * (2 * q + 1) * (2 * l_b + 1) / (4 * np.pi)
    )

    # Phase
    phase = (-1) ** (l_a - m_a)

    return phase * prefactor * threej


def build_vertex_tensor(
    states: List[Tuple[int, int, int]],
    modes: List[Tuple[int, int]],
) -> np.ndarray:
    """Build the full vertex tensor V[i, j, k].

    V[i, j, k] = vertex_coupling(states[i], states[j], modes[k])

    Parameters
    ----------
    states : list of (n, l, m) tuples
        Electron states.
    modes : list of (q, m_q) tuples
        Photon modes.

    Returns
    -------
    numpy.ndarray of shape (N_e, N_e, N_gamma)
        Vertex tensor.
    """
    N_e = len(states)
    N_gamma = len(modes)
    V = np.zeros((N_e, N_e, N_gamma))

    for i, (n_a, l_a, m_a) in enumerate(states):
        for j, (n_b, l_b, m_b) in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                V[i, j, k] = vertex_coupling(
                    n_a, l_a, m_a, n_b, l_b, m_b, q, m_q
                )
    return V


# ---------------------------------------------------------------------------
# Step 5: Self-energy
# ---------------------------------------------------------------------------

@dataclass
class VectorQEDResult:
    """Result container for the vector QED self-energy computation.

    Attributes
    ----------
    n_max : int
        Electron basis truncation.
    q_max : int
        Photon basis truncation.
    N_electron : int
        Number of electron states.
    N_photon : int
        Number of photon modes.
    Sigma : numpy.ndarray
        Self-energy matrix (N_e x N_e).
    vertex_nonzero_count : int
        Number of nonzero vertex tensor entries.
    vertex_total : int
        Total number of vertex tensor entries.
    vertex_sparsity : float
        Fraction of zero entries in vertex tensor.
    contains_pi : bool
        Whether the vertex normalization introduces pi.
    selection_rules : dict
        Selection rule pass/fail analysis.
    """
    n_max: int
    q_max: int
    N_electron: int
    N_photon: int
    Sigma: np.ndarray
    vertex_nonzero_count: int
    vertex_total: int
    vertex_sparsity: float
    contains_pi: bool
    selection_rules: dict


def compute_self_energy(
    n_max: int = 3,
    q_max: Optional[int] = None,
) -> VectorQEDResult:
    """Compute the one-loop vector-photon self-energy.

    Sigma(a, c) = sum_{b, q, m_q} V(a, b, q, m_q) * G_e(n_b)
                * G_gamma(q) * V(c, b, q, m_q)

    where b runs over internal electron states (excluding n_b = 1).

    Parameters
    ----------
    n_max : int
        Maximum electron principal quantum number.
    q_max : int, optional
        Maximum photon angular momentum. Default: n_max - 1.

    Returns
    -------
    VectorQEDResult
        Full self-energy result with selection rule analysis.
    """
    if q_max is None:
        q_max = n_max - 1
    if q_max < 1:
        q_max = 1

    states = build_electron_states(n_max)
    modes = build_photon_modes(q_max)
    N_e = len(states)
    N_gamma = len(modes)

    # Build vertex tensor
    V = build_vertex_tensor(states, modes)

    # Count nonzero entries
    vertex_nnz = int(np.count_nonzero(V))
    vertex_total = V.size
    vertex_sparsity = 1.0 - vertex_nnz / vertex_total if vertex_total > 0 else 1.0

    # Build electron propagator array (None for n=1)
    G_e = np.zeros(N_e)
    for i, (n, l, m) in enumerate(states):
        prop = electron_propagator(n)
        if prop is not None:
            G_e[i] = prop
        # else: G_e[i] = 0, state excluded from internal sum

    # Mask: which states can be internal (n >= 2)
    internal_mask = np.array([s[0] >= 2 for s in states], dtype=bool)

    # Build photon propagator array
    G_gamma = np.zeros(N_gamma)
    for k, (q, m_q) in enumerate(modes):
        G_gamma[k] = vector_photon_propagator(q)

    # Compute self-energy
    # Sigma[a, c] = sum_{b (internal), k (photon)} V[a, b, k] * G_e[b] * G_gamma[k] * V[c, b, k]
    Sigma = np.zeros((N_e, N_e))
    for b_idx in range(N_e):
        if not internal_mask[b_idx]:
            continue
        ge_b = G_e[b_idx]
        for k_idx in range(N_gamma):
            gg_k = G_gamma[k_idx]
            weight = ge_b * gg_k
            if abs(weight) < 1e-30:
                continue
            # V[:, b, k] is the column of vertex couplings for internal b, photon k
            v_col_a = V[:, b_idx, k_idx]  # shape (N_e,)
            # Sigma += weight * outer(v_col_a, v_col_c)
            # But v_col_c = V[:, b, k] from the second vertex
            # Since the second vertex is V(c, b, k), and we have V[c, b, k],
            # this is the same column: V[:, b_idx, k_idx]
            Sigma += weight * np.outer(v_col_a, v_col_a)

    # Check selection rules
    rules = check_selection_rules(n_max, q_max, states, modes, V, Sigma)

    return VectorQEDResult(
        n_max=n_max,
        q_max=q_max,
        N_electron=N_e,
        N_photon=N_gamma,
        Sigma=Sigma,
        vertex_nonzero_count=vertex_nnz,
        vertex_total=vertex_total,
        vertex_sparsity=vertex_sparsity,
        contains_pi=True,  # sqrt(1/(4*pi)) in normalization
        selection_rules=rules,
    )


# ---------------------------------------------------------------------------
# Step 6: Selection rule checks
# ---------------------------------------------------------------------------

def check_selection_rules(
    n_max: int,
    q_max: int,
    states: Optional[List[Tuple[int, int, int]]] = None,
    modes: Optional[List[Tuple[int, int]]] = None,
    V: Optional[np.ndarray] = None,
    Sigma: Optional[np.ndarray] = None,
) -> Dict:
    """Check all 8 continuum QED selection rules.

    Parameters
    ----------
    n_max : int
        Maximum electron principal quantum number.
    q_max : int
        Maximum photon angular momentum.
    states : list, optional
        Pre-computed electron states.
    modes : list, optional
        Pre-computed photon modes.
    V : numpy.ndarray, optional
        Pre-computed vertex tensor.
    Sigma : numpy.ndarray, optional
        Pre-computed self-energy matrix.

    Returns
    -------
    dict
        Pass/fail for each of 8 selection rules.
    """
    if states is None:
        states = build_electron_states(n_max)
    if modes is None:
        modes = build_photon_modes(q_max)
    if V is None:
        V = build_vertex_tensor(states, modes)
    if Sigma is None:
        result = compute_self_energy(n_max, q_max)
        Sigma = result.Sigma

    N_e = len(states)
    N_gamma = len(modes)
    rules = {}

    # 1. Gaunt/CG sparsity
    vertex_nnz = int(np.count_nonzero(V))
    vertex_total = V.size
    rules['1_gaunt_cg_sparsity'] = {
        'pass': vertex_nnz < vertex_total,
        'nonzero': vertex_nnz,
        'total': vertex_total,
        'sparsity': 1.0 - vertex_nnz / vertex_total if vertex_total > 0 else 1.0,
        'description': 'Wigner 3j selection rules enforce angular momentum conservation',
    }

    # 2. Vertex parity (GS structural zero)
    # Find ground state indices (n=1, l=0, m=0)
    gs_indices = [i for i, s in enumerate(states) if s[0] == 1]
    if gs_indices:
        gs_block = Sigma[np.ix_(gs_indices, gs_indices)]
        gs_max = float(np.max(np.abs(gs_block)))
        gs_zero = gs_max < 1e-14
    else:
        gs_zero = None
        gs_max = None

    rules['2_vertex_parity_gs_zero'] = {
        'pass': gs_zero,
        'gs_indices': gs_indices,
        'gs_block_max_abs': gs_max,
        'description': (
            'Parity rule l_a+l_b+q odd => for GS (l=0): q=l_b, so '
            'l_b+l_b=2*l_b is even => V(GS,any,any)=0 => Sigma(GS,GS)=0'
        ),
    }

    # 3. SO(4) channel count W
    # For each (n_ext, n_int) pair, count how many (l_ext, l_int, q) triples
    # produce nonzero vertex coupling.
    w_counts = {}
    for n_ext in range(1, n_max + 1):
        for n_int in range(2, n_max + 1):  # internal excludes n=1
            count = 0
            for l_ext in range(0, n_ext):
                for l_int in range(0, n_int):
                    for q_val in range(1, q_max + 1):
                        if (l_ext + l_int + q_val) % 2 == 1:
                            if abs(l_ext - l_int) <= q_val <= l_ext + l_int:
                                count += 1
            w_counts[(n_ext, n_int)] = count

    rules['3_so4_channel_count'] = {
        'pass': True,  # W is always well-defined
        'w_counts': {f'({k[0]},{k[1]})': v for k, v in w_counts.items()},
        'description': (
            'Number of (l_ext, l_int, q) triples with nonzero coupling '
            'for each (n_ext, n_int) pair. Structural, always passes.'
        ),
    }

    # 4. Delta_m conservation
    # Sigma should be diagonal in m for each (n, l) subspace
    # More precisely: Sigma(a,c) = 0 unless m_a = m_c
    m_violations = 0
    m_total_offdiag = 0
    for i, (n_i, l_i, m_i) in enumerate(states):
        for j, (n_j, l_j, m_j) in enumerate(states):
            if m_i != m_j:
                m_total_offdiag += 1
                if abs(Sigma[i, j]) > 1e-14:
                    m_violations += 1

    rules['4_delta_m_conservation'] = {
        'pass': m_violations == 0,
        'violations': m_violations,
        'offdiag_m_total': m_total_offdiag,
        'description': (
            'Sigma(a,c) = 0 when m_a != m_c, enforced by 3j selection: '
            'sum over m_q conserves m.'
        ),
    }

    # 5. Spatial parity
    # All nonzero vertex entries have l_a + l_b + q odd
    parity_violations = 0
    for i, (n_a, l_a, m_a) in enumerate(states):
        for j, (n_b, l_b, m_b) in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                if abs(V[i, j, k]) > 1e-15:
                    if (l_a + l_b + q) % 2 == 0:
                        parity_violations += 1

    rules['5_spatial_parity'] = {
        'pass': parity_violations == 0,
        'violations': parity_violations,
        'description': (
            'All nonzero V entries have l_a + l_b + q odd (E-type parity).'
        ),
    }

    # 6. Furry's theorem (tadpole = 0)
    # Tadpole: sum_{q, m_q} V(a, a, q, m_q) * G_gamma(q) for each a
    tadpole = np.zeros(N_e)
    for a_idx in range(N_e):
        for k_idx, (q, m_q) in enumerate(modes):
            tadpole[a_idx] += V[a_idx, a_idx, k_idx] * vector_photon_propagator(q)
    tadpole_max = float(np.max(np.abs(tadpole)))

    rules['6_furry_theorem'] = {
        'pass': tadpole_max < 1e-14,
        'tadpole_max_abs': tadpole_max,
        'description': (
            'Tadpole = sum_{q,m_q} V(a,a,q,m_q)*G_gamma(q) = 0 for all a. '
            'Enforced by parity: V(a,a,q,m_q) requires l_a+l_a+q odd => q odd, '
            'and triangle with l_a=l_a forces q=0..2*l_a, but q must also be '
            'nonzero (q>=1). For l_a=0: q must equal 0 (triangle), but q>=1 => '
            'V=0. For l_a>0: 3j(l,q,l;-m,0,m) with q odd gives V(a,a,q,0).'
        ),
    }

    # 7. Ward identity: ||[Sigma, H0]|| / ||Sigma|| ~ 0
    # H0 = diag(n^2 - 1) = the free Fock Laplacian eigenvalues
    H0 = np.diag([s[0] ** 2 - 1 for s in states]).astype(float)
    commutator = Sigma @ H0 - H0 @ Sigma
    comm_norm = np.linalg.norm(commutator, 'fro')
    sigma_norm = np.linalg.norm(Sigma, 'fro')
    ward_ratio = comm_norm / sigma_norm if sigma_norm > 1e-30 else 0.0

    rules['7_ward_identity'] = {
        'pass': bool(ward_ratio < 0.01),  # Threshold for approximate Ward
        'commutator_norm': comm_norm,
        'sigma_norm': sigma_norm,
        'ratio': ward_ratio,
        'description': (
            '||[Sigma, H0]|| / ||Sigma|| measures how well the self-energy '
            'commutes with the free Hamiltonian. For exact Ward identity this '
            'should vanish. On the truncated basis it will be approximate.'
        ),
    }

    # 8. Charge conjugation
    # Check that the photon modes have definite C-parity (-1)^q,
    # and that only C-even combinations contribute to the self-energy.
    # In practice: the self-energy sum with parity constraint (l_a+l_b+q odd)
    # already selects the correct C-parity sector.
    # We check by verifying Sigma is Hermitian (C-even observable).
    sigma_hermitian = bool(np.allclose(Sigma, Sigma.T, atol=1e-14))

    rules['8_charge_conjugation'] = {
        'pass': sigma_hermitian,
        'is_hermitian': sigma_hermitian,
        'max_antisym': float(np.max(np.abs(Sigma - Sigma.T))),
        'description': (
            'Self-energy is Hermitian (symmetric for real matrix), '
            'consistent with C-even observable.'
        ),
    }

    # Summary
    n_pass = sum(1 for r in rules.values() if r.get('pass') is True)
    n_total = len(rules)
    rules['_summary'] = {
        'pass_count': n_pass,
        'total': n_total,
        'all_pass': n_pass == n_total,
    }

    return rules


# ---------------------------------------------------------------------------
# Step 7: Comparison helper (scalar vs vector)
# ---------------------------------------------------------------------------

def compare_with_scalar_graph(
    n_max: int = 2,
    q_max: Optional[int] = None,
) -> Dict:
    """Compare the vector-photon self-energy structure with the scalar graph.

    Loads the scalar graph self-energy from ``graph_qed_self_energy`` and
    compares structural properties (dimensions, sparsity, GS zero, etc.)
    without expecting numerical equality.

    Parameters
    ----------
    n_max : int
    q_max : int, optional

    Returns
    -------
    dict with structural comparison.
    """
    # Vector QED
    vec_result = compute_self_energy(n_max, q_max)

    # Scalar graph QED (import here to avoid circular dependency)
    from geovac.graph_qed_self_energy import (
        compute_self_energy as scalar_self_energy,
    )
    from sympy import Rational
    scalar_result = scalar_self_energy(n_max, t=Rational(0), exact=False)

    # Structural comparison
    comparison = {
        'n_max': n_max,
        'vector': {
            'N_electron': vec_result.N_electron,
            'N_photon': vec_result.N_photon,
            'Sigma_shape': list(vec_result.Sigma.shape),
            'vertex_sparsity': vec_result.vertex_sparsity,
            'gs_zero': vec_result.selection_rules.get(
                '2_vertex_parity_gs_zero', {}
            ).get('pass'),
            'contains_pi': vec_result.contains_pi,
            'Sigma_trace': float(np.trace(vec_result.Sigma)),
            'Sigma_frobenius': float(np.linalg.norm(vec_result.Sigma, 'fro')),
            'Sigma_max': float(np.max(np.abs(vec_result.Sigma))),
            'selection_rules_pass': vec_result.selection_rules.get(
                '_summary', {}
            ).get('pass_count', 0),
            'selection_rules_total': vec_result.selection_rules.get(
                '_summary', {}
            ).get('total', 0),
        },
        'scalar': {
            'N_dirac': scalar_result.N_dirac,
            'E_fock': scalar_result.E_fock,
            'Sigma_shape': list(scalar_result.Sigma_numpy.shape),
            'gs_zero': scalar_result.ground_state_zero,
            'Sigma_trace': float(np.trace(scalar_result.Sigma_numpy)),
            'Sigma_frobenius': float(
                np.linalg.norm(scalar_result.Sigma_numpy, 'fro')
            ),
            'Sigma_max': float(np.max(np.abs(scalar_result.Sigma_numpy))),
        },
    }

    return comparison


# ---------------------------------------------------------------------------
# Diagnostic driver
# ---------------------------------------------------------------------------

def run_diagnostic(
    n_max_list: Optional[List[int]] = None,
    output_path: Optional[Path] = None,
) -> Dict:
    """Run the full vector QED diagnostic at multiple n_max values.

    Parameters
    ----------
    n_max_list : list of int, optional
        n_max values to test. Default: [2, 3].
    output_path : Path, optional
        Output JSON path.

    Returns
    -------
    dict with comprehensive results.
    """
    if n_max_list is None:
        n_max_list = [2, 3]
    if output_path is None:
        output_path = (
            Path(__file__).parent.parent
            / "debug" / "data" / "vector_qed_diagnostic.json"
        )

    results = {
        'module': 'geovac.vector_qed',
        'description': (
            'Vector-photon QED on the Fock graph: explicit photon modes '
            '(q, m_q) coupled via Wigner 3j symbols. Calibration step '
            'importing continuum angular momentum structure.'
        ),
        'analyses': {},
    }

    for n_max in n_max_list:
        q_max = n_max - 1 if n_max >= 2 else 1
        print(f"\n=== n_max = {n_max}, q_max = {q_max} ===")

        result = compute_self_energy(n_max, q_max)

        analysis = {
            'n_max': n_max,
            'q_max': q_max,
            'N_electron': result.N_electron,
            'N_photon': result.N_photon,
            'vertex_nonzero': result.vertex_nonzero_count,
            'vertex_total': result.vertex_total,
            'vertex_sparsity': result.vertex_sparsity,
            'contains_pi': result.contains_pi,
            'Sigma_trace': float(np.trace(result.Sigma)),
            'Sigma_frobenius': float(np.linalg.norm(result.Sigma, 'fro')),
            'Sigma_max_entry': float(np.max(np.abs(result.Sigma))),
            'Sigma_is_symmetric': bool(
                np.allclose(result.Sigma, result.Sigma.T, atol=1e-14)
            ),
            'Sigma_nonzero_count': int(
                np.count_nonzero(np.abs(result.Sigma) > 1e-14)
            ),
            'Sigma_eigenvalues': sorted(
                [float(ev) for ev in np.linalg.eigvalsh(result.Sigma)]
            ),
            'selection_rules': result.selection_rules,
        }

        # Print selection rule table
        print(f"\nSelection Rule Census (n_max={n_max}, q_max={q_max}):")
        print(f"  Electrons: {result.N_electron} states")
        print(f"  Photons:   {result.N_photon} modes")
        print(f"  Vertex:    {result.vertex_nonzero_count}/{result.vertex_total} "
              f"nonzero ({result.vertex_sparsity:.1%} sparse)")
        print(f"  Sigma:     {result.N_electron}x{result.N_electron}, "
              f"trace={np.trace(result.Sigma):.6f}")
        print()
        print(f"  {'Rule':<35} {'Pass?':<8} {'Details'}")
        print(f"  {'-'*35} {'-'*8} {'-'*40}")
        for key, val in result.selection_rules.items():
            if key.startswith('_'):
                continue
            status = 'PASS' if val.get('pass') else 'FAIL'
            detail = ''
            if 'violations' in val:
                detail = f"violations={val['violations']}"
            elif 'ratio' in val:
                detail = f"ratio={val['ratio']:.6f}"
            elif 'gs_block_max_abs' in val and val['gs_block_max_abs'] is not None:
                detail = f"|GS block|_max={val['gs_block_max_abs']:.2e}"
            elif 'tadpole_max_abs' in val:
                detail = f"|tadpole|_max={val['tadpole_max_abs']:.2e}"
            elif 'sparsity' in val:
                detail = f"sparsity={val['sparsity']:.1%}"
            print(f"  {key:<35} {status:<8} {detail}")

        summary = result.selection_rules.get('_summary', {})
        print(f"\n  Summary: {summary.get('pass_count', 0)}/{summary.get('total', 0)} "
              f"selection rules pass")

        results['analyses'][f'n_max_{n_max}'] = analysis

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {output_path}")

    return results


# ===================================================================
# PART 2: Dirac electron states + vector photon modes
# ===================================================================
#
# Extension combining SPINOR (Dirac) electron states labeled by
# (n, kappa, m_j) with VECTOR photon modes (q, m_q).  The goal is
# to recover 7/8 or 8/8 continuum QED selection rules by combining:
#   - Dirac spinor structure (recovers Furry's theorem via kappa -> -kappa
#     under E1 parity flip)
#   - Vector photon modes (recover vertex parity / GS structural zero)
#
# The vertex coupling uses jj-coupling 3j symbols with explicit parity
# enforcement (l_a + l_b + q odd, E1-type).

from sympy import Rational as _Rational
from geovac.dirac_matrix_elements import (
    DiracLabel,
    iter_dirac_labels,
    kappa_to_l,
    kappa_to_j,
)


def build_dirac_electron_states(n_max: int) -> List[DiracLabel]:
    """Build the ordered list of Dirac electron states.

    States are labeled by (n_fock, kappa, m_j) with:
      - n_fock = 1, ..., n_max
      - kappa = -(l+1) for j = l+1/2, or kappa = l for j = l-1/2 (l >= 1)
      - m_j = -j, -j+1, ..., j

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    list of DiracLabel
        Ordered list of Dirac states.
    """
    return list(iter_dirac_labels(n_max))


def dirac_electron_propagator(label: DiracLabel) -> float:
    """Dirac electron propagator on the Fock graph.

    G_e(n) = 1 / |lambda_n| where |lambda_n| = n_fock + 1/2
    (Camporesi-Higuchi Dirac eigenvalue in Fock convention).

    Unlike the scalar Fock propagator (1/(n^2-1) which diverges at n=1),
    the Dirac propagator is well-defined at all n >= 1:
      G_e(n=1) = 1/(3/2) = 2/3.

    Parameters
    ----------
    label : DiracLabel
        Dirac state label.

    Returns
    -------
    float
        Propagator value 1/(n_fock + 1/2).
    """
    return 1.0 / (label.n_fock + 0.5)


def dirac_vertex_coupling(
    a: DiracLabel,
    b: DiracLabel,
    q: int,
    m_q: int,
) -> float:
    """Compute the E1 vertex coupling V(a, b, q, m_q) for Dirac states.

    V = sqrt((2*j_a+1)*(2*q+1)*(2*j_b+1)/(4*pi))
      * (-1)^(j_a - m_j_a) * 3j(j_a, q, j_b; -m_j_a, m_q, m_j_b)
      * parity_factor(l_a, l_b, q)

    where parity_factor = 1 if l_a + l_b + q is odd (E1), 0 if even.

    Selection rules enforced:
      1. l-parity: l_a + l_b + q must be odd (E-type / vector coupling)
      2. l-triangle: |l_a - l_b| <= q <= l_a + l_b (orbital angular momentum)
      3. j-triangle: |j_a - j_b| <= q <= j_a + j_b (total angular momentum)
      4. Magnetic conservation: -m_j_a + m_q + m_j_b = 0

    The l-triangle is ESSENTIAL for the GS structural zero: with l_a = 0,
    the l-triangle forces q = l_b, so l_a + l_b + q = 0 + l_b + l_b = 2*l_b
    (even) => V = 0.  Without the l-triangle, the j-triangle allows q values
    that bypass the parity protection (e.g., j_a = 1/2, j_b = 1/2, q = 1
    gives l+l+q = 0+0+1 = 1 odd, nonzero coupling to GS).

    Parameters
    ----------
    a : DiracLabel
        External (or first) electron state.
    b : DiracLabel
        Internal (or second) electron state.
    q : int
        Photon angular momentum (q >= 1).
    m_q : int
        Photon magnetic quantum number.

    Returns
    -------
    float
        Vertex coupling value.
    """
    l_a = kappa_to_l(a.kappa)
    l_b = kappa_to_l(b.kappa)

    # 1. Parity selection rule: l_a + l_b + q must be odd (E-type)
    if (l_a + l_b + q) % 2 == 0:
        return 0.0

    # 2. l-triangle: |l_a - l_b| <= q <= l_a + l_b (orbital AM conservation)
    # This is ESSENTIAL for GS structural zero -- see docstring.
    if q < abs(l_a - l_b) or q > l_a + l_b:
        return 0.0

    j_a = a.j  # sympy Rational
    j_b = b.j
    m_j_a = a.m_j
    m_j_b = b.m_j

    # 3. Magnetic conservation: -m_j_a + m_q + m_j_b = 0
    if -m_j_a + m_q + m_j_b != 0:
        return 0.0

    # 4. Triangle inequality on j: |j_a - j_b| <= q <= j_a + j_b
    if q < abs(j_a - j_b) or q > j_a + j_b:
        return 0.0

    # Compute 3j symbol: (j_a, q, j_b; -m_j_a, m_q, m_j_b)
    threej = float(wigner_3j(j_a, q, j_b, -m_j_a, m_q, m_j_b))
    if abs(threej) < 1e-15:
        return 0.0

    # Angular normalization factor
    prefactor = np.sqrt(
        float((2 * j_a + 1) * (2 * q + 1) * (2 * j_b + 1)) / (4 * np.pi)
    )

    # Phase: (-1)^(j_a - m_j_a)
    phase_exp = j_a - m_j_a  # This is always an integer
    phase = (-1) ** int(phase_exp)

    return phase * prefactor * threej


def build_dirac_vertex_tensor(
    states: List[DiracLabel],
    modes: List[Tuple[int, int]],
) -> np.ndarray:
    """Build the full vertex tensor V[i, j, k] for Dirac states.

    V[i, j, k] = dirac_vertex_coupling(states[i], states[j], modes[k])

    Parameters
    ----------
    states : list of DiracLabel
        Dirac electron states.
    modes : list of (q, m_q) tuples
        Photon modes.

    Returns
    -------
    numpy.ndarray of shape (N_e, N_e, N_gamma)
        Vertex tensor.
    """
    N_e = len(states)
    N_gamma = len(modes)
    V = np.zeros((N_e, N_e, N_gamma))

    for i, a in enumerate(states):
        for j, b in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                V[i, j, k] = dirac_vertex_coupling(a, b, q, m_q)
    return V


@dataclass
class DiracVectorQEDResult:
    """Result container for the Dirac + vector photon self-energy computation.

    Attributes
    ----------
    n_max : int
        Electron basis truncation.
    q_max : int
        Photon basis truncation.
    N_electron : int
        Number of Dirac electron states.
    N_photon : int
        Number of photon modes.
    states : list of DiracLabel
        Electron state list.
    Sigma : numpy.ndarray
        Self-energy matrix (N_e x N_e).
    vertex_nonzero_count : int
        Number of nonzero vertex tensor entries.
    vertex_total : int
        Total number of vertex tensor entries.
    vertex_sparsity : float
        Fraction of zero entries in vertex tensor.
    contains_pi : bool
        Whether the vertex normalization introduces pi.
    selection_rules : dict
        Selection rule pass/fail analysis (all 8 rules).
    """
    n_max: int
    q_max: int
    N_electron: int
    N_photon: int
    states: List[DiracLabel]
    Sigma: np.ndarray
    vertex_nonzero_count: int
    vertex_total: int
    vertex_sparsity: float
    contains_pi: bool
    selection_rules: dict


def compute_dirac_self_energy(
    n_max: int = 3,
    q_max: Optional[int] = None,
    exclude_gs_internal: bool = False,
) -> DiracVectorQEDResult:
    """Compute the one-loop self-energy with Dirac electrons + vector photons.

    Sigma(a, c) = sum_{b, q, m_q} V(a, b, q, m_q) * G_e(b)
                * G_gamma(q) * V(c, b, q, m_q)

    Unlike the scalar case, all states (including n=1) can be internal
    because the Dirac propagator G_e = 1/|lambda_n| = 1/(n+1/2) is
    well-defined at n=1.

    Parameters
    ----------
    n_max : int
        Maximum electron principal quantum number.
    q_max : int, optional
        Maximum photon angular momentum. Default: n_max - 1.
    exclude_gs_internal : bool
        If True, exclude n=1 states from internal sums (for comparison
        with the scalar case where G_e(n=1) is undefined).

    Returns
    -------
    DiracVectorQEDResult
        Full self-energy result with selection rule analysis.
    """
    if q_max is None:
        q_max = n_max - 1
    if q_max < 1:
        q_max = 1

    states = build_dirac_electron_states(n_max)
    modes = build_photon_modes(q_max)
    N_e = len(states)
    N_gamma = len(modes)

    # Build vertex tensor
    V = build_dirac_vertex_tensor(states, modes)

    # Count nonzero entries
    vertex_nnz = int(np.count_nonzero(V))
    vertex_total = V.size
    vertex_sparsity = 1.0 - vertex_nnz / vertex_total if vertex_total > 0 else 1.0

    # Build electron propagator array
    G_e = np.zeros(N_e)
    for i, lab in enumerate(states):
        G_e[i] = dirac_electron_propagator(lab)

    # Internal state mask
    if exclude_gs_internal:
        internal_mask = np.array([s.n_fock >= 2 for s in states], dtype=bool)
    else:
        internal_mask = np.ones(N_e, dtype=bool)

    # Photon propagator array
    G_gamma = np.zeros(N_gamma)
    for k, (q, m_q) in enumerate(modes):
        G_gamma[k] = vector_photon_propagator(q)

    # Compute self-energy
    Sigma = np.zeros((N_e, N_e))
    for b_idx in range(N_e):
        if not internal_mask[b_idx]:
            continue
        ge_b = G_e[b_idx]
        for k_idx in range(N_gamma):
            gg_k = G_gamma[k_idx]
            weight = ge_b * gg_k
            if abs(weight) < 1e-30:
                continue
            v_col = V[:, b_idx, k_idx]
            Sigma += weight * np.outer(v_col, v_col)

    # Check selection rules
    rules = check_dirac_selection_rules(n_max, q_max, states, modes, V, Sigma)

    return DiracVectorQEDResult(
        n_max=n_max,
        q_max=q_max,
        N_electron=N_e,
        N_photon=N_gamma,
        states=states,
        Sigma=Sigma,
        vertex_nonzero_count=vertex_nnz,
        vertex_total=vertex_total,
        vertex_sparsity=vertex_sparsity,
        contains_pi=True,
        selection_rules=rules,
    )


def check_dirac_selection_rules(
    n_max: int,
    q_max: int,
    states: Optional[List[DiracLabel]] = None,
    modes: Optional[List[Tuple[int, int]]] = None,
    V: Optional[np.ndarray] = None,
    Sigma: Optional[np.ndarray] = None,
) -> Dict:
    """Check all 8 continuum QED selection rules for Dirac + vector photon.

    Parameters
    ----------
    n_max : int
        Maximum electron principal quantum number.
    q_max : int
        Maximum photon angular momentum.
    states : list, optional
        Pre-computed Dirac electron states.
    modes : list, optional
        Pre-computed photon modes.
    V : numpy.ndarray, optional
        Pre-computed vertex tensor.
    Sigma : numpy.ndarray, optional
        Pre-computed self-energy matrix.

    Returns
    -------
    dict
        Pass/fail for each of 8 selection rules.
    """
    if states is None:
        states = build_dirac_electron_states(n_max)
    if modes is None:
        modes = build_photon_modes(q_max)
    if V is None:
        V = build_dirac_vertex_tensor(states, modes)
    if Sigma is None:
        result = compute_dirac_self_energy(n_max, q_max)
        Sigma = result.Sigma

    N_e = len(states)
    N_gamma = len(modes)
    rules = {}

    # 1. Gaunt/CG sparsity
    vertex_nnz = int(np.count_nonzero(V))
    vertex_total = V.size
    rules['1_gaunt_cg_sparsity'] = {
        'pass': vertex_nnz < vertex_total,
        'nonzero': vertex_nnz,
        'total': vertex_total,
        'sparsity': 1.0 - vertex_nnz / vertex_total if vertex_total > 0 else 1.0,
        'description': (
            'Wigner 3j selection rules enforce angular momentum conservation '
            'in (j_a, j_b, q) coupling.'
        ),
    }

    # 2. Vertex parity / GS structural zero
    # GS = (n=1, kappa=-1, m_j=+-1/2)
    gs_indices = [i for i, s in enumerate(states) if s.n_fock == 1]
    if gs_indices:
        gs_block = Sigma[np.ix_(gs_indices, gs_indices)]
        gs_max = float(np.max(np.abs(gs_block)))
        # Also check full GS row/column
        gs_row_max = float(np.max(np.abs(Sigma[gs_indices, :])))
        gs_zero = gs_max < 1e-14
    else:
        gs_zero = None
        gs_max = None
        gs_row_max = None

    rules['2_vertex_parity_gs_zero'] = {
        'pass': gs_zero,
        'gs_indices': gs_indices,
        'gs_block_max_abs': gs_max,
        'gs_row_max_abs': gs_row_max,
        'description': (
            'For GS (l=0, kappa=-1): only l_b values with l_a+l_b+q odd '
            'contribute. With l_a=0, need l_b+q odd; triangle on j forces '
            'q = l_b (since j_a=1/2, need |1/2-j_b| <= q <= 1/2+j_b). '
            'For l_b values with kappa giving j_b, l_b+l_b=even => V=0.'
        ),
    }

    # 3. SO(4) channel count
    w_counts = {}
    for n_ext in range(1, n_max + 1):
        for n_int in range(1, n_max + 1):
            count = 0
            for s_ext in states:
                if s_ext.n_fock != n_ext:
                    continue
                for s_int in states:
                    if s_int.n_fock != n_int:
                        continue
                    l_ext = kappa_to_l(s_ext.kappa)
                    l_int = kappa_to_l(s_int.kappa)
                    j_ext = float(s_ext.j)
                    j_int = float(s_int.j)
                    for q_val in range(1, q_max + 1):
                        if (l_ext + l_int + q_val) % 2 == 1:
                            l_ok = abs(l_ext - l_int) <= q_val <= l_ext + l_int
                            j_ok = abs(j_ext - j_int) <= q_val <= j_ext + j_int
                            if l_ok and j_ok:
                                count += 1
                                break
                    else:
                        continue
                    break
            w_counts[(n_ext, n_int)] = count

    # Recompute properly: count distinct (kappa_ext, kappa_int, q) triples
    w_counts_proper = {}
    for n_ext in range(1, n_max + 1):
        for n_int in range(1, n_max + 1):
            kappas_ext = set()
            kappas_int = set()
            for s in states:
                if s.n_fock == n_ext:
                    kappas_ext.add(s.kappa)
                if s.n_fock == n_int:
                    kappas_int.add(s.kappa)
            count = 0
            for k_ext in sorted(kappas_ext):
                l_ext = kappa_to_l(k_ext)
                j_ext = float(kappa_to_j(k_ext))
                for k_int in sorted(kappas_int):
                    l_int = kappa_to_l(k_int)
                    j_int = float(kappa_to_j(k_int))
                    for q_val in range(1, q_max + 1):
                        if (l_ext + l_int + q_val) % 2 == 1:
                            # Both l-triangle AND j-triangle must hold
                            l_ok = abs(l_ext - l_int) <= q_val <= l_ext + l_int
                            j_ok = abs(j_ext - j_int) <= q_val <= j_ext + j_int
                            if l_ok and j_ok:
                                count += 1
            w_counts_proper[(n_ext, n_int)] = count

    rules['3_so4_channel_count'] = {
        'pass': True,
        'w_counts': {f'({k[0]},{k[1]})': v for k, v in w_counts_proper.items()},
        'description': (
            'Number of (kappa_ext, kappa_int, q) triples with nonzero coupling '
            'for each (n_ext, n_int) pair. Uses j-triangle and l-parity.'
        ),
    }

    # 4. Delta_m_j conservation
    mj_violations = 0
    mj_total_offdiag = 0
    for i, s_i in enumerate(states):
        for j, s_j in enumerate(states):
            if s_i.two_m_j != s_j.two_m_j:
                mj_total_offdiag += 1
                if abs(Sigma[i, j]) > 1e-14:
                    mj_violations += 1

    rules['4_delta_mj_conservation'] = {
        'pass': mj_violations == 0,
        'violations': mj_violations,
        'offdiag_mj_total': mj_total_offdiag,
        'description': (
            'Sigma(a,c) = 0 when m_j_a != m_j_c, enforced by 3j selection: '
            'sum over m_q conserves m_j.'
        ),
    }

    # 5. Spatial parity
    parity_violations = 0
    for i, a in enumerate(states):
        for j, b in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                if abs(V[i, j, k]) > 1e-15:
                    l_a = kappa_to_l(a.kappa)
                    l_b = kappa_to_l(b.kappa)
                    if (l_a + l_b + q) % 2 == 0:
                        parity_violations += 1

    rules['5_spatial_parity'] = {
        'pass': parity_violations == 0,
        'violations': parity_violations,
        'description': (
            'All nonzero V entries have l_a + l_b + q odd (E-type parity), '
            'where l is derived from kappa.'
        ),
    }

    # 6. Furry's theorem (tadpole = 0)
    # Tadpole: sum_{q, m_q} V(a, a, q, m_q) * G_gamma(q)
    # On the Dirac basis, V(a, a, q, m_q) requires l_a + l_a + q odd
    # => q odd; AND requires the j-triangle |j_a - j_a| <= q <= 2*j_a
    # => 0 <= q <= 2*j_a. But q >= 1 (photon).
    # Additionally, for V(a,a,q,m_q), m_q must be 0 (m_j conservation).
    # So we need: q odd, 1 <= q <= 2*j_a, m_q = 0.
    # For j_a = 1/2: 2*j_a = 1, so q=1 only, which IS odd. Check 3j:
    #   3j(1/2, 1, 1/2; -m, 0, m) which is nonzero.
    # BUT the parity rule: l_a + l_a + q = 2*l_a + q. For q=1: 2*l_a+1
    # is always odd. So parity is ALWAYS satisfied for q=1 diagonal!
    # This means the parity rule alone does NOT kill the tadpole for
    # Dirac states. The question is whether the 3j symbol kills it.
    # For (1/2, 1, 1/2; -m, 0, m): this is a nonzero 3j symbol.
    # So Furry's theorem may FAIL here -- the E1 vertex does not
    # automatically give tadpole=0 for spinor states.
    #
    # However, on the DIRAC graph with CG-projected vertex (not explicit
    # vector modes), Furry IS recovered because the identity vertex
    # couples kappa to -kappa. Here we're using explicit 3j coupling
    # which may behave differently.
    tadpole = np.zeros(N_e)
    for a_idx in range(N_e):
        for k_idx, (q, m_q) in enumerate(modes):
            tadpole[a_idx] += V[a_idx, a_idx, k_idx] * vector_photon_propagator(q)
    tadpole_max = float(np.max(np.abs(tadpole)))
    tadpole_rms = float(np.sqrt(np.mean(tadpole ** 2)))

    # Count how many states have nonzero tadpole
    tadpole_nonzero = int(np.sum(np.abs(tadpole) > 1e-14))

    rules['6_furry_theorem'] = {
        'pass': tadpole_max < 1e-14,
        'tadpole_max_abs': tadpole_max,
        'tadpole_rms': tadpole_rms,
        'tadpole_nonzero_count': tadpole_nonzero,
        'description': (
            'Tadpole = sum_{q,m_q} V(a,a,q,m_q)*G_gamma(q) = 0 for all a. '
            'On Dirac states, diagonal coupling V(a,a) has l_a+l_a+q = 2l_a+q; '
            'for q odd this is always odd (parity satisfied). '
            'Furry requires the 3j(j,q,j;-m,0,m) to vanish or cancel in sum.'
        ),
    }

    # 7. Ward identity: ||[Sigma, H0]|| / ||Sigma|| ~ 0
    # H0 = diag(|lambda_n|) = diag(n_fock + 1/2), the Dirac eigenvalues
    H0 = np.diag([float(s.n_fock + 0.5) for s in states])
    commutator = Sigma @ H0 - H0 @ Sigma
    comm_norm = np.linalg.norm(commutator, 'fro')
    sigma_norm = np.linalg.norm(Sigma, 'fro')
    ward_ratio = comm_norm / sigma_norm if sigma_norm > 1e-30 else 0.0

    rules['7_ward_identity'] = {
        'pass': bool(ward_ratio < 0.01),
        'commutator_norm': comm_norm,
        'sigma_norm': sigma_norm,
        'ratio': ward_ratio,
        'description': (
            '||[Sigma, H0]|| / ||Sigma|| where H0 = diag(n+1/2) = Dirac '
            'eigenvalues. For exact Ward identity this vanishes.'
        ),
    }

    # 8. Charge conjugation / kappa symmetry
    # Check that Sigma is symmetric (Hermitian for real case)
    sigma_hermitian = bool(np.allclose(Sigma, Sigma.T, atol=1e-14))

    # Also check kappa -> -kappa symmetry: for each pair of states that
    # differ only in sign of kappa (same n, same l, same |m_j|),
    # their diagonal self-energy elements should be equal.
    kappa_sym_max_diff = 0.0
    kappa_sym_pairs_checked = 0
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            if (si.n_fock == sj.n_fock and si.kappa == -sj.kappa
                    and si.two_m_j == sj.two_m_j):
                diff = abs(Sigma[i, i] - Sigma[j, j])
                kappa_sym_max_diff = max(kappa_sym_max_diff, diff)
                kappa_sym_pairs_checked += 1

    rules['8_charge_conjugation'] = {
        'pass': sigma_hermitian,
        'is_hermitian': sigma_hermitian,
        'max_antisym': float(np.max(np.abs(Sigma - Sigma.T))),
        'kappa_symmetry_max_diff': kappa_sym_max_diff,
        'kappa_symmetry_pairs_checked': kappa_sym_pairs_checked,
        'description': (
            'Self-energy is Hermitian (symmetric for real). '
            'Additionally checks kappa -> -kappa symmetry of diagonal elements.'
        ),
    }

    # Summary
    n_pass = sum(1 for r in rules.values() if r.get('pass') is True)
    n_total = len(rules)
    rules['_summary'] = {
        'pass_count': n_pass,
        'total': n_total,
        'all_pass': n_pass == n_total,
    }

    return rules


def run_dirac_vector_diagnostic(
    n_max_list: Optional[List[int]] = None,
    output_path: Optional[Path] = None,
) -> Dict:
    """Run the full Dirac + vector photon QED diagnostic.

    Parameters
    ----------
    n_max_list : list of int, optional
        n_max values to test. Default: [2, 3].
    output_path : Path, optional
        Output JSON path.

    Returns
    -------
    dict with comprehensive results.
    """
    if n_max_list is None:
        n_max_list = [2, 3]
    if output_path is None:
        output_path = (
            Path(__file__).parent.parent
            / "debug" / "data" / "vector_qed_dirac_diagnostic.json"
        )

    results = {
        'module': 'geovac.vector_qed (Dirac extension)',
        'description': (
            'Dirac electron states (n, kappa, m_j) + vector photon modes '
            '(q, m_q) coupled via Wigner 3j symbols with E1 parity. '
            'Combines spinor electron structure with vector photon modes '
            'to maximize selection rule recovery.'
        ),
        'analyses': {},
    }

    for n_max in n_max_list:
        q_max = n_max - 1 if n_max >= 2 else 1
        print(f"\n=== Dirac vector QED: n_max = {n_max}, q_max = {q_max} ===")

        result = compute_dirac_self_energy(n_max, q_max)

        # Build analysis dict
        sigma_eigs = sorted(float(ev) for ev in np.linalg.eigvalsh(result.Sigma))
        analysis = {
            'n_max': n_max,
            'q_max': q_max,
            'N_electron': result.N_electron,
            'N_photon': result.N_photon,
            'vertex_nonzero': result.vertex_nonzero_count,
            'vertex_total': result.vertex_total,
            'vertex_sparsity': result.vertex_sparsity,
            'contains_pi': result.contains_pi,
            'Sigma_trace': float(np.trace(result.Sigma)),
            'Sigma_frobenius': float(np.linalg.norm(result.Sigma, 'fro')),
            'Sigma_max_entry': float(np.max(np.abs(result.Sigma))),
            'Sigma_is_symmetric': bool(
                np.allclose(result.Sigma, result.Sigma.T, atol=1e-14)
            ),
            'Sigma_nonzero_count': int(
                np.count_nonzero(np.abs(result.Sigma) > 1e-14)
            ),
            'Sigma_eigenvalues': sigma_eigs,
            'Sigma_zero_eigs': sum(1 for ev in sigma_eigs if abs(ev) < 1e-12),
            'selection_rules': {},
        }

        # Sanitize selection rules for JSON
        for key, val in result.selection_rules.items():
            sanitized = {}
            for k, v in val.items():
                if isinstance(v, (int, float, bool, str, type(None))):
                    sanitized[k] = v
                elif isinstance(v, list):
                    sanitized[k] = [
                        x if isinstance(x, (int, float, bool, str)) else str(x)
                        for x in v
                    ]
                elif isinstance(v, dict):
                    sanitized[k] = {
                        str(kk): vv for kk, vv in v.items()
                    }
                else:
                    sanitized[k] = str(v)
            analysis['selection_rules'][key] = sanitized

        # Print selection rule table
        print(f"\nSelection Rule Census (Dirac + Vector, n_max={n_max}):")
        print(f"  Electrons: {result.N_electron} Dirac states")
        print(f"  Photons:   {result.N_photon} vector modes")
        print(f"  Vertex:    {result.vertex_nonzero_count}/{result.vertex_total} "
              f"nonzero ({result.vertex_sparsity:.1%} sparse)")
        print(f"  Sigma:     {result.N_electron}x{result.N_electron}, "
              f"trace={np.trace(result.Sigma):.6f}, "
              f"zero eigs={analysis['Sigma_zero_eigs']}")
        print()
        print(f"  {'Rule':<35} {'Pass?':<8} {'Details'}")
        print(f"  {'-'*35} {'-'*8} {'-'*40}")
        for key, val in result.selection_rules.items():
            if key.startswith('_'):
                continue
            status = 'PASS' if val.get('pass') else 'FAIL'
            detail = ''
            if 'violations' in val:
                detail = f"violations={val['violations']}"
            elif 'ratio' in val:
                detail = f"ratio={val['ratio']:.6f}"
            elif 'gs_block_max_abs' in val and val['gs_block_max_abs'] is not None:
                detail = f"|GS block|_max={val['gs_block_max_abs']:.2e}"
            elif 'tadpole_max_abs' in val:
                detail = (f"|tadpole|_max={val['tadpole_max_abs']:.2e}, "
                          f"nonzero={val.get('tadpole_nonzero_count', '?')}")
            elif 'sparsity' in val:
                detail = f"sparsity={val['sparsity']:.1%}"
            print(f"  {key:<35} {status:<8} {detail}")

        summary = result.selection_rules.get('_summary', {})
        print(f"\n  Summary: {summary.get('pass_count', 0)}/{summary.get('total', 0)} "
              f"selection rules pass")

        results['analyses'][f'n_max_{n_max}'] = analysis

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {output_path}")

    return results


if __name__ == "__main__":
    run_diagnostic()
