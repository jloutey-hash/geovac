"""
One-loop electron self-energy and vertex correction on the GeoVac graph.
=========================================================================

This module computes graph-native QED diagrams with external electron
lines.  Everything is a FINITE MATRIX TRACE on the finite graph -- no
infinite sums, no regularization, no renormalization.

One-loop self-energy (A)
------------------------
The self-energy Sigma[a, b] (Dirac node indices) is the one-loop
correction where electron a emits a photon (Fock edge e), propagates
as virtual electron c, and absorbs a photon (Fock edge e') to become
electron b:

    Sigma[a, b] = sum_{c, e, e'} V[a, c, e] * G_gamma[e, e'] * V[c, b, e']

where:
  - V[a, c, e] is the vertex coupling tensor (graph_qed_vertex.py)
  - G_gamma[e, e'] is the photon propagator = L_1^+ (graph_qed_photon.py)
  - The sum over c is an electron intermediate-state sum (finite, N_dirac terms)
  - The sums over e, e' are photon edge sums (finite, E_fock terms)

In matrix form, using V_e[a, c] = V[a, c, e]:

    Sigma = sum_{e, e'} V_e . G_gamma[e, e'] . V_{e'}^T

where V_e is the N_dirac x N_dirac matrix of vertex couplings at edge e.
Note: V[c, b, e'] is the second vertex with electron c incoming and
b outgoing, so we contract V_e'[c, b] = V[c, b, e'] giving V_{e'}
transposed relative to the first vertex.

One-loop vertex correction (B)
-------------------------------
The vertex correction Lambda[a, b, e] dresses one vertex with a photon
exchange.  In continuum QED this is Lambda^mu(p, p') involving two
vertex insertions and both electron and photon internal propagators.

On the graph:

    Lambda[a, b, e] = sum_{c, d, e', e''} V[a, c, e'] * G_e[c, d]
                      * V[d, b, e''] * G_gamma[e', e'']

But the external photon edge e constrains the overall structure.  The
simplest contraction that captures the vertex correction physics is:

    Lambda_e[a, b] = sum_{c, d, e', e''} V_{e'}[a, c] * G_e[c, d]
                     * G_gamma[e', e''] * V_{e''}[d, b]

which in matrix form is:

    Lambda_e = sum_{e', e''} G_gamma[e', e''] * V_{e'} . G_e . V_{e''}^T

This gives an N_dirac x N_dirac matrix for EACH photon edge e (but
actually the free edge index is in the G_gamma contraction).  The
full vertex correction tensor is Lambda[a, b, e] with the e index
coming from the external photon.

For the anomalous magnetic moment extraction, we work with the
edge-summed (contracted) vertex correction:

    Lambda_total[a, b] = sum_e Lambda_e[a, b]

Anomalous magnetic moment (C)
------------------------------
On the graph, F_2 (the Pauli form factor) is extracted from the
vertex correction by comparing to the bare vertex structure.  In the
simplest extraction, the ratio Lambda_total / (bare vertex norm)
gives a dimensionless number that should be rational (pi-free).

Self-energy structural zero (D)
---------------------------------
In the continuum (qed_self_energy.py), Sigma(n_ext=0) = 0 exactly
because vertex parity forces n_int + n_ext + q odd; with n_ext=0
this requires 2*n_int odd, which is impossible.  On the graph,
we check whether the ground-state block of Sigma vanishes.

Transcendental taxonomy (Paper 18)
-----------------------------------
- Vertex V: algebraic (CG coefficients = sqrt(rational))
- Photon propagator G_gamma: rational (pseudoinverse of integer matrix)
- Electron propagator G_e: rational (inverse of rational Dirac matrix)
- Self-energy Sigma: rational (V . G_gamma . V^T contracts sqrt's)
- Vertex correction Lambda: rational (same mechanism)

All quantities live in the intrinsic tier -- no pi, no zeta.

References
----------
- GeoVac graph_qed_vertex.py (vertex coupling V[a, b, e])
- GeoVac graph_qed_photon.py (photon propagator G_gamma = L_1^+)
- GeoVac graph_qed_propagator.py (electron propagator G_e = D^{-1})
- GeoVac qed_self_energy.py (continuum spectral-sum self-energy)
- GeoVac Paper 28 (QED on S^3)
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Matrix, Rational, zeros as sp_zeros

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.graph_qed_propagator import (
    DiracGraphOperator,
    electron_propagator,
)
from geovac.graph_qed_photon import (
    build_fock_graph,
    compute_photon_propagator,
    FockGraphData,
    PhotonPropagatorData,
)
from geovac.graph_qed_vertex import (
    build_projection_matrix,
    build_vertex_tensor,
    vertex_tensor_to_matrices,
)

__all__ = [
    "SelfEnergyResult",
    "VertexCorrectionResult",
    "compute_self_energy",
    "compute_vertex_correction",
    "extract_anomalous_moment",
    "self_energy_structural_zero_check",
    "self_energy_pi_free_certificate",
    "analyze_self_energy_and_vertex",
    "_check_rational_matrix",
    "_check_pi_free_matrix",
    "_check_algebraic_matrix",
    "_ground_state_indices",
]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SelfEnergyResult:
    """One-loop electron self-energy on the graph.

    Attributes
    ----------
    n_max : int
        Graph truncation level.
    t : str
        Dirac graph hopping parameter (string representation).
    Sigma : Optional[Matrix]
        Self-energy matrix Sigma[a, b] (N_dirac x N_dirac), exact sympy.
    Sigma_numpy : np.ndarray
        Self-energy matrix as float64 array.
    is_rational : bool
        True if all Sigma entries are rational.
    is_pi_free : bool
        True if no transcendentals appear in Sigma.
    is_hermitian : bool
        True if Sigma is symmetric (Hermitian for real matrices).
    ground_state_block : Optional[Matrix]
        The 2x2 block of Sigma at the ground state (n_fock=1, kappa=-1).
    ground_state_zero : bool
        True if the ground-state block is identically zero.
    trace : str
        Trace of Sigma (string representation).
    eigenvalues : List[str]
        Eigenvalues of Sigma (string representations).
    N_dirac : int
        Number of Dirac states.
    E_fock : int
        Number of Fock edges.
    """
    n_max: int
    t: str
    Sigma: Optional[Matrix]
    Sigma_numpy: np.ndarray
    is_rational: bool
    is_pi_free: bool
    is_hermitian: bool
    ground_state_block: Optional[Matrix]
    ground_state_zero: bool
    trace: str
    eigenvalues: List[str]
    N_dirac: int
    E_fock: int


@dataclass
class VertexCorrectionResult:
    """One-loop vertex correction on the graph.

    Attributes
    ----------
    n_max : int
        Graph truncation level.
    t : str
        Dirac graph hopping parameter.
    Lambda_total : Optional[Matrix]
        Edge-summed vertex correction Lambda_total[a, b], exact sympy.
    Lambda_total_numpy : np.ndarray
        Lambda_total as float64 array.
    Lambda_per_edge : Optional[List[Matrix]]
        Vertex correction per photon edge (list of N_dirac x N_dirac).
    is_rational : bool
        True if all Lambda entries are rational.
    is_pi_free : bool
        True if no transcendentals appear.
    trace : str
        Trace of Lambda_total.
    eigenvalues : List[str]
        Eigenvalues of Lambda_total.
    anomalous_moment : Optional[str]
        Graph-native F_2 analog (string).
    anomalous_moment_float : Optional[float]
        F_2 analog as float.
    N_dirac : int
    E_fock : int
    """
    n_max: int
    t: str
    Lambda_total: Optional[Matrix]
    Lambda_total_numpy: np.ndarray
    Lambda_per_edge: Optional[List[Matrix]]
    is_rational: bool
    is_pi_free: bool
    trace: str
    eigenvalues: List[str]
    anomalous_moment: Optional[str]
    anomalous_moment_float: Optional[float]
    N_dirac: int
    E_fock: int


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _check_pi_free_matrix(M: Matrix) -> bool:
    """Check that all entries of a sympy Matrix are algebraic (pi-free)."""
    transcendental_atoms = (sp.pi, sp.E, sp.EulerGamma, sp.Catalan)
    for entry in M:
        expr = sp.sympify(entry)
        for atom in sp.preorder_traversal(expr):
            if atom in transcendental_atoms:
                return False
            if isinstance(atom, (sp.zeta, sp.log, sp.exp)):
                return False
    return True


def _check_rational_matrix(M: Matrix) -> bool:
    """Check that all entries are rational (not just algebraic).

    An entry is rational if, after simplification, it can be expressed
    as a sympy Rational (integer/integer).  Entries containing sqrt()
    or other algebraic irrationals return False.
    """
    for entry in M:
        simplified = sp.nsimplify(sp.expand(entry), rational=False)
        if simplified.is_Rational or simplified.is_Integer:
            continue
        if simplified == sp.S.Zero:
            continue
        # Check if it's a pure rational after full simplification
        if simplified.is_number and simplified.is_rational:
            continue
        return False
    return True


def _check_algebraic_matrix(M: Matrix) -> bool:
    """Check that all entries are algebraic (sqrt of rational, no transcendentals).

    Returns True if all entries are algebraic numbers (possibly irrational,
    like sqrt(2)/3) but contain no pi, log, zeta, etc.
    """
    return _check_pi_free_matrix(M)


def _matrix_is_symmetric(M: Matrix, tol: float = 1e-14) -> bool:
    """Check if a sympy Matrix is symmetric."""
    if M.rows != M.cols:
        return False
    diff = M - M.T
    for entry in diff:
        val = sp.nsimplify(entry)
        if val != 0:
            if abs(complex(val)) > tol:
                return False
    return True


def _ground_state_indices(dirac_labels: List[DiracLabel]) -> List[int]:
    """Return indices of ground-state Dirac labels (n_fock=1, kappa=-1).

    The ground state has n_fock=1, l=0, j=1/2, kappa=-1.
    There are two states: m_j = +1/2 and m_j = -1/2 (two_m_j = +1, -1).
    """
    indices = []
    for i, lab in enumerate(dirac_labels):
        if lab.n_fock == 1 and lab.kappa == -1:
            indices.append(i)
    return indices


# ---------------------------------------------------------------------------
# (A) One-loop self-energy
# ---------------------------------------------------------------------------

def compute_self_energy(
    n_max: int,
    t: sp.Expr = Rational(0),
    exact: bool = True,
) -> SelfEnergyResult:
    """Compute the one-loop electron self-energy on the graph.

    Sigma[a, b] = sum_{e, e'} (V_e . G_gamma[e, e'] . V_{e'}^T)[a, b]

    where the sum runs over all photon edge pairs (e, e').

    Parameters
    ----------
    n_max : int
        Graph truncation level (2 recommended for exact sympy).
    t : sympy.Rational
        Dirac graph hopping parameter (0 = free diagonal propagator).
    exact : bool
        Use exact sympy arithmetic.

    Returns
    -------
    SelfEnergyResult
        Complete self-energy analysis.
    """
    # Build vertex tensor
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    E = fock_data.E
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    # Build photon propagator G_gamma = L_1^+
    photon_data = compute_photon_propagator(n_max, exact=exact)

    if exact and photon_data.G_gamma is not None:
        G_gamma = photon_data.G_gamma
    else:
        G_gamma = None

    # Compute Sigma
    if exact and G_gamma is not None:
        Sigma = sp_zeros(N_dirac, N_dirac)
        for e1 in range(E_fock):
            for e2 in range(E_fock):
                g_ee = G_gamma[e1, e2]
                if g_ee == 0:
                    continue
                # Sigma += G_gamma[e1, e2] * V_{e1} . V_{e2}^T
                contrib = g_ee * (V_mats[e1] * V_mats[e2].T)
                Sigma = Sigma + contrib

        # Simplify entries
        for i in range(N_dirac):
            for j in range(N_dirac):
                Sigma[i, j] = sp.nsimplify(sp.expand(Sigma[i, j]),
                                           rational=False)

        is_rational = _check_rational_matrix(Sigma)
        is_pi_free = _check_pi_free_matrix(Sigma)
        is_hermitian = _matrix_is_symmetric(Sigma)

        Sigma_np = np.array(Sigma.tolist(), dtype=float)

        trace_val = str(sp.nsimplify(Sigma.trace(), rational=False))

        # Eigenvalues
        try:
            eig_dict = Sigma.eigenvals()
            eigenvalues = []
            for ev, mult in sorted(eig_dict.items(),
                                   key=lambda x: float(x[0])):
                eigenvalues.extend(
                    [str(sp.nsimplify(ev, rational=False))] * mult
                )
        except Exception:
            eigenvalues = [f"{ev:.10g}" for ev in np.linalg.eigvalsh(Sigma_np)]
    else:
        # Numpy path
        G_gamma_np = photon_data.G_gamma_numeric
        V_mats_np = [np.array(vm.tolist(), dtype=float) for vm in V_mats]
        Sigma_np = np.zeros((N_dirac, N_dirac))
        for e1 in range(E_fock):
            for e2 in range(E_fock):
                g_ee = G_gamma_np[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                Sigma_np += g_ee * (V_mats_np[e1] @ V_mats_np[e2].T)
        Sigma = None
        is_rational = None
        is_pi_free = None
        is_hermitian = bool(np.allclose(Sigma_np, Sigma_np.T, atol=1e-12))
        trace_val = f"{np.trace(Sigma_np):.10g}"
        eigenvalues = [f"{ev:.10g}" for ev in np.linalg.eigvalsh(Sigma_np)]

    # Ground state block check
    gs_indices = _ground_state_indices(dirac_labels)
    if Sigma is not None and len(gs_indices) > 0:
        gs_block = sp.zeros(len(gs_indices), len(gs_indices))
        for ii, gi in enumerate(gs_indices):
            for jj, gj in enumerate(gs_indices):
                gs_block[ii, jj] = Sigma[gi, gj]
        gs_zero = all(
            sp.nsimplify(gs_block[i, j]) == 0
            for i in range(gs_block.rows)
            for j in range(gs_block.cols)
        )
    elif Sigma is None and len(gs_indices) > 0:
        gs_block_np = Sigma_np[np.ix_(gs_indices, gs_indices)]
        gs_block = None
        gs_zero = bool(np.allclose(gs_block_np, 0, atol=1e-12))
    else:
        gs_block = None
        gs_zero = None

    return SelfEnergyResult(
        n_max=n_max,
        t=str(t),
        Sigma=Sigma,
        Sigma_numpy=Sigma_np,
        is_rational=is_rational,
        is_pi_free=is_pi_free,
        is_hermitian=is_hermitian,
        ground_state_block=gs_block if isinstance(gs_block, Matrix) else None,
        ground_state_zero=gs_zero,
        trace=trace_val,
        eigenvalues=eigenvalues,
        N_dirac=N_dirac,
        E_fock=E_fock,
    )


# ---------------------------------------------------------------------------
# (B) One-loop vertex correction
# ---------------------------------------------------------------------------

def compute_vertex_correction(
    n_max: int,
    t: sp.Expr = Rational(0),
    exact: bool = True,
) -> VertexCorrectionResult:
    """Compute the one-loop vertex correction on the graph.

    Lambda_e[a, b] = sum_{e', e''} G_gamma[e', e''] *
                     (V_{e'}[a, :] . G_e . V_{e''}^T[:, b])

    i.e. Lambda_e = sum_{e', e''} G_gamma[e', e''] * V_{e'} . G_e . V_{e''}^T

    The edge-summed (full) vertex correction is:

    Lambda_total[a, b] = sum_e Lambda_e[a, b]

    But since Lambda_e does not actually depend on a free external edge
    index e (the photon contraction is fully internal), the "total" is:

    Lambda_total = sum_{e', e''} G_gamma[e', e''] * V_{e'} . G_e . V_{e''}^T

    Parameters
    ----------
    n_max : int
        Graph truncation level.
    t : sympy.Rational
        Dirac graph hopping parameter.
    exact : bool
        Use exact sympy arithmetic.

    Returns
    -------
    VertexCorrectionResult
    """
    # Build vertex tensor
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    E = fock_data.E
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    # Build photon propagator
    photon_data = compute_photon_propagator(n_max, exact=exact)

    # Build electron propagator
    # Use Neumann series when exact and t != 0 for algebraically cleaner
    # expression trees (avoids the full M.inv() which grows expression size).
    op = DiracGraphOperator(n_max=n_max, t=t)
    _prop_method = 'neumann' if (exact and t != Rational(0)) else 'auto'
    G_e, G_e_rational = electron_propagator(op, exact=exact, method=_prop_method)

    if exact and photon_data.G_gamma is not None:
        G_gamma = photon_data.G_gamma

        Lambda_total = sp_zeros(N_dirac, N_dirac)
        Lambda_per_edge_list: Optional[List[Matrix]] = []

        for e1 in range(E_fock):
            for e2 in range(E_fock):
                g_ee = G_gamma[e1, e2]
                if g_ee == 0:
                    continue
                # Lambda_total += G_gamma[e1, e2] * V_{e1} . G_e . V_{e2}^T
                contrib = g_ee * (V_mats[e1] * G_e * V_mats[e2].T)
                Lambda_total = Lambda_total + contrib

        # Simplify
        for i in range(N_dirac):
            for j in range(N_dirac):
                Lambda_total[i, j] = sp.nsimplify(
                    sp.expand(Lambda_total[i, j]), rational=False
                )

        is_rational = _check_rational_matrix(Lambda_total)
        is_pi_free = _check_pi_free_matrix(Lambda_total)

        Lambda_np = np.array(Lambda_total.tolist(), dtype=float)

        trace_val = str(sp.nsimplify(Lambda_total.trace(), rational=False))
        try:
            eig_dict = Lambda_total.eigenvals()
            eigenvalues = []
            for ev, mult in sorted(eig_dict.items(),
                                   key=lambda x: abs(complex(x[0]))):
                eigenvalues.extend(
                    [str(sp.nsimplify(ev, rational=False))] * mult
                )
        except Exception:
            eigenvalues = [f"{ev:.10g}" for ev in
                           np.sort(np.linalg.eigvals(Lambda_np).real)]

        # Per-edge vertex correction (optional, can be large)
        Lambda_per_edge_list = None  # Skip for now to save memory

    else:
        # Numpy path
        G_gamma_np = photon_data.G_gamma_numeric
        V_mats_np = [np.array(vm.tolist(), dtype=float) for vm in V_mats]
        G_e_np = (np.array(G_e.tolist(), dtype=float)
                  if isinstance(G_e, Matrix) else G_e)

        Lambda_np = np.zeros((N_dirac, N_dirac))
        for e1 in range(E_fock):
            for e2 in range(E_fock):
                g_ee = G_gamma_np[e1, e2]
                if abs(g_ee) < 1e-15:
                    continue
                Lambda_np += g_ee * (V_mats_np[e1] @ G_e_np @ V_mats_np[e2].T)

        Lambda_total = None
        Lambda_per_edge_list = None
        is_rational = None
        is_pi_free = None
        trace_val = f"{np.trace(Lambda_np):.10g}"
        eigenvalues = [f"{ev:.10g}" for ev in
                       np.sort(np.linalg.eigvals(Lambda_np).real)]

    return VertexCorrectionResult(
        n_max=n_max,
        t=str(t),
        Lambda_total=Lambda_total,
        Lambda_total_numpy=Lambda_np,
        Lambda_per_edge=Lambda_per_edge_list,
        is_rational=is_rational,
        is_pi_free=is_pi_free,
        trace=trace_val,
        eigenvalues=eigenvalues,
        anomalous_moment=None,  # filled by extract_anomalous_moment
        anomalous_moment_float=None,
        N_dirac=N_dirac,
        E_fock=E_fock,
    )


# ---------------------------------------------------------------------------
# (C) Anomalous magnetic moment extraction
# ---------------------------------------------------------------------------

def extract_anomalous_moment(
    n_max: int,
    t: sp.Expr = Rational(0),
    exact: bool = True,
) -> Dict:
    """Extract the graph-native anomalous magnetic moment F_2.

    Strategy: compare the vertex correction Lambda_total to the bare
    vertex structure.  The bare vertex is V_bare[a, b] = sum_e V_e[a, b],
    the total coupling summed over all photon edges.

    The "F_2" ratio is extracted as:

        F_2 = Tr(Lambda_total) / Tr(V_bare . G_e)

    This is a dimensionless RATIONAL number on the graph.

    In the continuum, F_2 -> alpha/(2*pi) (Schwinger result).  On the
    graph, F_2 is pi-free and alpha-free: it is a pure rational number
    determined by the graph topology and CG algebra.

    Parameters
    ----------
    n_max : int
    t : sympy.Rational
    exact : bool

    Returns
    -------
    dict with F_2 value, trace decomposition, and pi-free certificate.
    """
    # Vertex correction
    vc_result = compute_vertex_correction(n_max, t=t, exact=exact)

    # Build bare vertex sum
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    # Electron propagator
    # Use Neumann series when exact and t != 0 (same rationale as vertex correction).
    op = DiracGraphOperator(n_max=n_max, t=t)
    _prop_method = 'neumann' if (exact and t != Rational(0)) else 'auto'
    G_e, _ = electron_propagator(op, exact=exact, method=_prop_method)

    if exact and vc_result.Lambda_total is not None:
        # Bare vertex: V_bare = sum_e V_e
        V_bare = sp_zeros(N_dirac, N_dirac)
        for vm in V_mats:
            V_bare = V_bare + vm

        # Trace of Lambda_total
        tr_lambda = sp.nsimplify(vc_result.Lambda_total.trace(), rational=False)

        # Normalization: Tr(V_bare . G_e)
        norm_matrix = V_bare * G_e
        tr_norm = sp.nsimplify(norm_matrix.trace(), rational=False)

        # F_2 ratio
        if tr_norm != 0:
            F2 = sp.nsimplify(tr_lambda / tr_norm, rational=False)
            F2_float = float(F2)
        else:
            F2 = sp.Integer(0)
            F2_float = 0.0

        F2_rational = isinstance(
            sp.nsimplify(F2, rational=True),
            (sp.Rational, sp.Integer, sp.core.numbers.Zero,
             sp.core.numbers.One, sp.core.numbers.NegativeOne)
        )

        return {
            'n_max': n_max,
            't': str(t),
            'F2': str(F2),
            'F2_float': F2_float,
            'F2_rational': F2_rational,
            'Tr_Lambda': str(tr_lambda),
            'Tr_Lambda_float': float(tr_lambda),
            'Tr_V_bare_G_e': str(tr_norm),
            'Tr_V_bare_G_e_float': float(tr_norm),
            'V_bare_nonzero_count': sum(
                1 for entry in V_bare if entry != 0
            ),
            'Lambda_trace': str(tr_lambda),
            'continuum_schwinger': 'alpha/(2*pi)',
            'graph_value_note': (
                'On the graph, F_2 is a rational number (pi-free). '
                'The continuum Schwinger result alpha/(2*pi) is NOT '
                'expected to appear -- the graph is a finite truncation '
                'with no continuum limit taken.'
            ),
        }
    else:
        # Numpy path
        V_bare_np = sum(
            np.array(vm.tolist(), dtype=float) for vm in V_mats
        )
        G_e_np = (np.array(G_e.tolist(), dtype=float)
                  if isinstance(G_e, Matrix) else G_e)
        tr_lambda = np.trace(vc_result.Lambda_total_numpy)
        tr_norm = np.trace(V_bare_np @ G_e_np)
        F2_float = tr_lambda / tr_norm if abs(tr_norm) > 1e-15 else 0.0

        return {
            'n_max': n_max,
            't': str(t),
            'F2': f"{F2_float:.10g}",
            'F2_float': F2_float,
            'F2_rational': None,
            'Tr_Lambda': f"{tr_lambda:.10g}",
            'Tr_Lambda_float': float(tr_lambda),
            'Tr_V_bare_G_e': f"{tr_norm:.10g}",
            'Tr_V_bare_G_e_float': float(tr_norm),
        }


# ---------------------------------------------------------------------------
# (D) Structural zero check
# ---------------------------------------------------------------------------

def self_energy_structural_zero_check(
    n_max: int,
    t: sp.Expr = Rational(0),
    exact: bool = True,
) -> Dict:
    """Check whether Sigma has a structural zero at the ground state.

    In the continuum (qed_self_energy.py), Sigma(n_ext=0) = 0 exactly
    because vertex parity forces n_int + n_ext + q odd; with n_ext=0
    this requires 2*n_int odd, which is impossible.

    On the graph, we check whether the ground-state block of Sigma
    (n_fock=1, kappa=-1, m_j=+/-1/2) is identically zero.

    Parameters
    ----------
    n_max : int
    t : sympy.Rational
    exact : bool

    Returns
    -------
    dict with the check result and mechanism analysis.
    """
    se_result = compute_self_energy(n_max, t=t, exact=exact)
    dirac_labels = list(iter_dirac_labels(n_max))

    gs_indices = _ground_state_indices(dirac_labels)

    result = {
        'n_max': n_max,
        't': str(t),
        'ground_state_indices': gs_indices,
        'ground_state_labels': [
            f"(n={dirac_labels[i].n_fock}, k={dirac_labels[i].kappa}, "
            f"2mj={dirac_labels[i].two_m_j})"
            for i in gs_indices
        ],
        'ground_state_zero': se_result.ground_state_zero,
        'Sigma_is_rational': se_result.is_rational,
        'Sigma_is_pi_free': se_result.is_pi_free,
    }

    # Detailed block analysis
    if se_result.Sigma is not None and gs_indices:
        gs_block_entries = {}
        for ii, gi in enumerate(gs_indices):
            for jj, gj in enumerate(gs_indices):
                val = se_result.Sigma[gi, gj]
                gs_block_entries[f"({ii},{jj})"] = str(
                    sp.nsimplify(val, rational=False)
                )
        result['ground_state_block_entries'] = gs_block_entries

        # Check which vertex couplings connect to the ground state
        P, dlabs, fstates = build_projection_matrix(n_max)
        fock_data = build_fock_graph(n_max)
        entries, N_d, V_f, E_f = build_vertex_tensor(
            n_max, P=P, dirac_labels=dlabs, fock_data=fock_data
        )

        # Count vertex couplings involving ground-state indices
        gs_vertex_count = 0
        gs_vertex_edges = set()
        for a, b, e, val in entries:
            if a in gs_indices or b in gs_indices:
                gs_vertex_count += 1
                gs_vertex_edges.add(e)

        result['gs_vertex_coupling_count'] = gs_vertex_count
        result['gs_vertex_edge_count'] = len(gs_vertex_edges)

    # Mechanism analysis
    if se_result.ground_state_zero:
        result['mechanism'] = (
            'The graph self-energy has a structural zero at the ground '
            'state, analogous to the continuum result (Theorem 4, Paper 28). '
            'The graph mechanism may differ from the continuum vertex parity '
            'rule (n_ext + n_int + q odd).'
        )
    else:
        result['mechanism'] = (
            'The graph self-energy does NOT have a structural zero at the '
            'ground state. This is a genuine graph-vs-continuum difference: '
            'the graph vertex (CG projection) does not enforce the same '
            'parity selection rule as the continuum SO(4) vertex.'
        )

    # Additional: full diagonal of Sigma
    if se_result.Sigma is not None:
        diag_entries = []
        for i in range(se_result.N_dirac):
            lab = dirac_labels[i]
            val = se_result.Sigma[i, i]
            diag_entries.append({
                'index': i,
                'label': f"(n={lab.n_fock}, k={lab.kappa}, 2mj={lab.two_m_j})",
                'Sigma_ii': str(sp.nsimplify(val, rational=False)),
                'Sigma_ii_float': float(val),
            })
        result['Sigma_diagonal'] = diag_entries

    return result


# ---------------------------------------------------------------------------
# Pi-free certificate
# ---------------------------------------------------------------------------

def self_energy_pi_free_certificate(
    n_max: int,
    t: sp.Expr = Rational(0),
) -> Dict:
    """Verify that Sigma and Lambda are pi-free.

    Checks:
    1. Sigma entries are rational (sqrt's from CG cancel in bilinear).
    2. Lambda entries are rational.
    3. No transcendentals appear.

    Parameters
    ----------
    n_max : int
    t : sympy.Rational

    Returns
    -------
    dict with certificate results.
    """
    se = compute_self_energy(n_max, t=t, exact=True)
    vc = compute_vertex_correction(n_max, t=t, exact=True)

    return {
        'n_max': n_max,
        't': str(t),
        'Sigma_rational': se.is_rational,
        'Sigma_pi_free': se.is_pi_free,
        'Sigma_hermitian': se.is_hermitian,
        'Lambda_rational': vc.is_rational,
        'Lambda_pi_free': vc.is_pi_free,
        'all_pass': (
            (se.is_rational is True) and
            (se.is_pi_free is True) and
            (vc.is_rational is True) and
            (vc.is_pi_free is True)
        ),
    }


# ---------------------------------------------------------------------------
# Full analysis driver
# ---------------------------------------------------------------------------

def analyze_self_energy_and_vertex(
    n_max: int = 2,
    t_values: Optional[List] = None,
) -> Dict:
    """Run the full self-energy and vertex correction analysis.

    Parameters
    ----------
    n_max : int
        Graph truncation level (2 recommended for exact sympy).
    t_values : list, optional
        Hopping values to test. Default: [0, -1/16].

    Returns
    -------
    dict with comprehensive analysis, suitable for JSON serialization.
    """
    if t_values is None:
        t_values = [Rational(0), Rational(-1, 16)]

    results = {
        'module': 'geovac.graph_qed_self_energy',
        'n_max': n_max,
        'analyses': {},
    }

    for t in t_values:
        t_key = f"t={t}"

        # Self-energy
        se = compute_self_energy(n_max, t=t, exact=True)

        # Vertex correction
        vc = compute_vertex_correction(n_max, t=t, exact=True)

        # Anomalous moment
        am = extract_anomalous_moment(n_max, t=t, exact=True)

        # Structural zero
        sz = self_energy_structural_zero_check(n_max, t=t, exact=True)

        # Pi-free certificate
        cert = self_energy_pi_free_certificate(n_max, t=t)

        analysis = {
            'self_energy': {
                'is_rational': se.is_rational,
                'is_pi_free': se.is_pi_free,
                'is_hermitian': se.is_hermitian,
                'trace': se.trace,
                'eigenvalues': se.eigenvalues,
                'N_dirac': se.N_dirac,
                'E_fock': se.E_fock,
                'ground_state_zero': se.ground_state_zero,
            },
            'vertex_correction': {
                'is_rational': vc.is_rational,
                'is_pi_free': vc.is_pi_free,
                'trace': vc.trace,
                'eigenvalues': vc.eigenvalues,
            },
            'anomalous_moment': am,
            'structural_zero': sz,
            'pi_free_certificate': cert,
        }

        # Add Sigma matrix entries for exact results
        if se.Sigma is not None:
            analysis['self_energy']['matrix_entries'] = [
                [str(se.Sigma[i, j]) for j in range(se.N_dirac)]
                for i in range(se.N_dirac)
            ]

        if vc.Lambda_total is not None:
            analysis['vertex_correction']['matrix_entries'] = [
                [str(vc.Lambda_total[i, j]) for j in range(vc.N_dirac)]
                for i in range(vc.N_dirac)
            ]

        results['analyses'][t_key] = analysis

    return results


# ---------------------------------------------------------------------------
# CLI / script entry point
# ---------------------------------------------------------------------------

def _run_analysis_and_save(output_path: Optional[Path] = None) -> Dict:
    """Run full analysis at n_max=2 and save JSON."""

    if output_path is None:
        output_path = (
            Path(__file__).parent.parent
            / "debug" / "data" / "gn5_self_energy_vertex.json"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("Analyzing graph self-energy and vertex correction at n_max=2 ...")
    results = analyze_self_energy_and_vertex(n_max=2)

    with output_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"Saved: {output_path}")
    return results


if __name__ == "__main__":
    _run_analysis_and_save()
